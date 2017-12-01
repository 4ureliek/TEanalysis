package GAL::Parser::hapmap_genotypes;

use strict;
use vars qw($VERSION);
$VERSION = 0.2.0;

use base qw(GAL::Parser);
use GAL::Reader::DelimitedLine;

=head1 NAME

GAL::Parser::hapmap_genotype - <One line description of module's purpose here>

=head1 VERSION

This document describes GAL::Parser::hapmap_genotype version 0.2.0

=head1 SYNOPSIS

     use GAL::Parser::hapmap_genotype;

=for author to fill in:
     Brief code example(s) here showing commonest usage(s).
     This section will be as far as many users bother reading
     so make it as educational and exemplary as possible.

=head1 DESCRIPTION

=for author to fill in:
     Write a full description of the module and its features here.
     Use subsections (=head2, =head3) as appropriate.

=head1 METHODS

=cut

#-----------------------------------------------------------------------------

=head2 new

     Title   : new
     Usage   : GAL::Parser::hapmap_genotype->new();
     Function: Creates a GAL::Parser::hapmap_genotype object;
     Returns : A GAL::Parser::hapmap_genotype object
     Args    :

=cut

sub new {
	my ($class, @args) = @_;
	my $self = $class->SUPER::new(@args);
	return $self;
}

#-----------------------------------------------------------------------------

sub _initialize_args {
	my ($self, @args) = @_;

	######################################################################
	# This block of code handels class attributes.  Use the
	# @valid_attributes below to define the valid attributes for
	# this class.  You must have identically named get/set methods
	# for each attribute.  Leave the rest of this block alone!
	######################################################################
	my $args = $self->SUPER::_initialize_args(@args);
	my @valid_attributes = qw(file fh header_count); # Set valid class attributes here
	$self->set_attributes($args, @valid_attributes);
	######################################################################
}

#-----------------------------------------------------------------------------

sub header_count {
    my ($self, $header_count) = @_;
    
    if ($header_count) {
	$self->{header_count} = $header_count;
    }
    return $self->{header_count};
}

#-----------------------------------------------------------------------------

=head2 parse_record

 Title   : parse_record
 Usage   : $a = $self->parse_record();
 Function: Parse the data from a record.
 Returns : A hash ref needed by Feature.pm to create a Feature object
 Args    : A hash ref of fields that this sub can understand (In this case GFF3).

=cut

sub parse_record {
    my ($self, $data) = @_;
    
    # rsID alleles chrom pos strand assembly# center protLSID
    # assayLSID panelLSID QCcode NA18484 NA18485 NA18486
    
    # Col1: refSNP rs# identifier at the time of release (NB might merge 
    #       with another rs# in the future)
    # Col2: SNP alleles according to dbSNP
    # Col3: chromosome that SNP maps to 
    # Col4: chromosome position of SNP, in basepairs on reference sequence
    # Col5: strand of reference sequence that SNP maps to
    # Col6: version of reference sequence assembly
    # Col7: Center where genotyping was done
    # Col8: HapMap genotype center that produced the genotypes
    # Col9: LSID for HapMap protocol used for genotyping
    # Col10: LSID for HapMap assay used for genotyping
    # Col11: LSID for panel of individuals genotyped
    # Col12: QC-code, currently 'QC+' for all entries (for future use)
    # Col13: and on: observed genotypes of samples, one per column, sample 
    # identifiers in column headers (Coriell catalog numbers, example: 
    # NA10847). Duplicate samples have .dup suffix.
    
    my %record;
    @record{qw(rsID alleles chrom pos strand assembly center
	       protLSID assayLSID panelLSID QCcode)} = splice(@{$data}, 0, 11);
    
    my $seqid      = $record{chrom};
    my @protocol_data = split /:/, $record{protLSID};
    my $source     = join ':', ($record{center}, $protocol_data[4]);
    my $type       = 'SNV';
    my $start      = $record{pos};
    my $end        = $start;
    my $feature_id = join ':', ($source, $seqid, $start);
    my $score      = '.';
    my $strand     = $record{strand};
    my $phase      = '.';
    
    my $reference_seq ||= uc $self->fasta->seq($seqid, $start, $end);
    $reference_seq = $self->revcomp($reference_seq) if $strand eq '-';
    
    my %variant_seqs;
    my ($header) = $self->reader->headers;
    my @header_data = split /\s/, $header;
    splice(@header_data, 0, 11);
    @variant_seqs{@header_data} = @{$data};
    
    my @features;
    for my $individual (sort keys %variant_seqs) {
	my $this_feature_id = join ':', ($individual, $feature_id);
	my %variant_seqs_hash = map {$_ => 1 } split //, $variant_seqs{$individual};
	my @variant_seqs = keys %variant_seqs_hash;

	# Make sure that at least one variant sequence is not N
	next unless grep {$_ ne 'N'} @variant_seqs;
	# Make sure that at least one variant sequence differs from
        # the reference sequence.
	next unless grep {$_ ne $reference_seq} @variant_seqs;
	
	my $zygosity = scalar @variant_seqs > 1 ? 'heterozygous' : 'homozygous';
	
	my $attributes = {Reference_seq => [$reference_seq],
			  Variant_seq   => \@variant_seqs,
			  Zygosity      => [$zygosity],
			  ID            => [$this_feature_id],
			  individual    => [$individual],
		      };
	
	my $feature = {feature_id => $this_feature_id,
		       seqid      => $seqid,
		       source     => $source,
		       type       => $type,
		       start      => $start,
		       end        => $end,
		       score      => $score,
		       strand     => $strand,
		       phase      => $phase,
		       attributes => $attributes,
		   };
	push @features, $feature;
    }
    return scalar @features ? \@features : undef;
}

#-----------------------------------------------------------------------------

sub reader {
    my $self = shift;
    
    if (! $self->{reader}) {
	# The column names for the incoming data
	# my @field_names = qw(seqid source type start end score strand phase
	#	 	     attributes);
	my $header_count = $self->header_count;
	$header_count = defined $header_count ? $header_count : 1;
	my $reader = GAL::Reader::DelimitedLine->new(field_separator => qr/\s/,
						     header_count    => $header_count);
	$self->{reader} = $reader;
    }
    return $self->{reader};
}

#-----------------------------------------------------------------------------

=head1 DIAGNOSTICS

=for author to fill in:
     List every single error and warning message that the module can
     generate (even the ones that will "never happen"), with a full
     explanation of each problem, one or more likely causes, and any
     suggested remedies.

=over

=item C<< Error message here, perhaps with %s placeholders >>

[Description of error here]

=item C<< Another error message here >>

[Description of error here]

[Et cetera, et cetera]

=back

=head1 CONFIGURATION AND ENVIRONMENT

<GAL::Parser::hapmap_genotype> requires no configuration files or environment variables.

=head1 DEPENDENCIES

None.

=head1 INCOMPATIBILITIES

None reported.

=head1 BUGS AND LIMITATIONS

No bugs have been reported.

Please report any bugs or feature requests to:
barry.moore@genetics.utah.edu

=head1 AUTHOR

Barry Moore <barry.moore@genetics.utah.edu>

=head1 LICENCE AND COPYRIGHT

Copyright (c) 2010-2014, Barry Moore <barry.moore@genetics.utah.edu>.  All rights reserved.

    This module is free software; you can redistribute it and/or
    modify it under the same terms as Perl itself (See LICENSE).

=head1 DISCLAIMER OF WARRANTY

BECAUSE THIS SOFTWARE IS LICENSED FREE OF CHARGE, THERE IS NO WARRANTY
FOR THE SOFTWARE, TO THE EXTENT PERMITTED BY APPLICABLE LAW. EXCEPT WHEN
OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS AND/OR OTHER PARTIES
PROVIDE THE SOFTWARE "AS IS" WITHOUT WARRANTY OF ANY KIND, EITHER
EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE
ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE SOFTWARE IS WITH
YOU. SHOULD THE SOFTWARE PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL
NECESSARY SERVICING, REPAIR, OR CORRECTION.

IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING
WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MAY MODIFY AND/OR
REDISTRIBUTE THE SOFTWARE AS PERMITTED BY THE ABOVE LICENCE, BE
LIABLE TO YOU FOR DAMAGES, INCLUDING ANY GENERAL, SPECIAL, INCIDENTAL,
OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OR INABILITY TO USE
THE SOFTWARE (INCLUDING BUT NOT LIMITED TO LOSS OF DATA OR DATA BEING
RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR THIRD PARTIES OR A
FAILURE OF THE SOFTWARE TO OPERATE WITH ANY OTHER SOFTWARE), EVEN IF
SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE POSSIBILITY OF
SUCH DAMAGES.

=cut

1;
