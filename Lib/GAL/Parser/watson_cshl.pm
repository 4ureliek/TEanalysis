package GAL::Parser::watson_cshl;

use strict;
use vars qw($VERSION);

$VERSION = 0.2.0;
use base qw(GAL::Parser);
use GAL::Reader::DelimitedLine;

=head1 NAME

GAL::Parser::watson_cshl - Parse SNP files from James Watson's genome (CSHL version)

=head1 VERSION

This document describes GAL::Parser::watson_cshl version 0.2.0

=head1 SYNOPSIS

    use GAL::Parser::watson_cshl
    my $parser = GAL::Parser::watson_cshl->new(file => 'watson_cshl.txt');
    while (my $feature_hash = $parser->next_feature_hash) {
	print $parser->to_gff3($feature_hash) . "\n";
    }

=head1 DESCRIPTION

L<GAL::Parser::watson_cshl> parses SNP files from James Watson's
genome (CSHL version).

=head1 Constructor

New L<GAL::Parser::watson_cshl> objects are created by the class method
new.  Arguments should be passed to the constructor as a list (or
reference) of key value pairs.  All attributes of the
L<GAL::Parser::watson_cshl> object can be set in the call to new. An
simple example of object creation would look like this:

    my $parser = GAL::Parser::watson_cshl->new(file => 'watson_cshl.txt');

The constructor recognizes the following parameters which will set the
appropriate attributes:

=over

=item * C<< file => feature_file.txt >>

This optional parameter provides the filename for the file containing
the data to be parsed. While this parameter is optional either it, or
the following fh parameter must be set.

=item * C<< fh => feature_file.txt >>

This optional parameter provides a filehandle to read data from. While
this parameter is optional either it, or the following fh parameter
must be set.

=back

=cut

#-----------------------------------------------------------------------------

=head2 new

     Title   : new
     Usage   : GAL::Parser::watson_cshl->new();
     Function: Creates a GAL::Parser::watson_cshl object;
     Returns : A GAL::Parser::watson_cshl object
     Args    : See the attributes described above.

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
	my @valid_attributes = qw(); # Set valid class attributes here.
	$self->set_attributes($args, @valid_attributes);
	######################################################################
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
	my ($self, $record) = @_;

	# id chromosome coordinate reference_seq variant_seq
	# match_status rsid alternate_seq variant_count
	# alternate_seq_count total_coverage zygosity

	my $id         = $record->{id};
	my $seqid      = $record->{chromosome};
	my $source     = 'CSHL';

	my $type       = 'SNV';
	my $start      = $record->{coordinate};
	my $end        = $record->{coordinate};
	my $score      = '.';
	my $strand     = '+';
	my $phase      = '.';

	my $reference_seq = $record->{reference_seq};
	my @variant_seqs;
	push @variant_seqs, $record->{variant_seq};

	my @variant_reads;
	push @variant_reads, $record->{variant_count};

	if ($record->{alternate_seq} ne '.') {
		push @variant_seqs, $record->{alternate_seq};
		push @variant_reads,   $record->{alternate_seq_count};
	}

	# If we have reference_reads then push that seq to the variants
	my $reference_reads = $record->{total_coverage} - $record->{variant_count}  -
	  $record->{alternate_count};
	if ($reference_reads > 0) {
		push @variant_reads, $reference_reads;
		push @variant_seqs, $reference_seq;
	}

	my $total_reads = $record->{total_coverage};

	my $zygosity = scalar @variant_seqs > 1 ? 'heterozygous' : 'homozygous';
	my $their_zygosity = $record->{zygosity} eq 'het' ? 'heterozygous' : undef;

	my $intersected_snp = $record->{rs} =~ /rs\d+/ ? 'dbSNP:' . $record->{rs} : undef;

	my $attributes = {Reference_seq => [$reference_seq],
			  Variant_seq   => \@variant_seqs,
			  ID            => [$id],
			  Variant_reads => \@variant_reads,
			  Total_reads   => [$total_reads],
			  Zygosity      => [$zygosity],
			 };

	push @{$attributes->{Intersected_feature}}, $intersected_snp
	  if $intersected_snp;

	my $feature_data = {feature_id => $id,
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

	return $feature_data;
}

#-----------------------------------------------------------------------------

=head2 reader

 Title   : reader
 Usage   : $a = $self->reader
 Function: Return the reader object.
 Returns : A L<GAL::Reader::DelimitedLine> singleton.
 Args    : None

=cut

sub reader {
  my $self = shift;

  # The columns are:
  #
  # BCM_local_SNP_ID -- unique ID for referring to the SNPs ahead of
  # submission to dbSNP (we can talk about what and when to submit to
  # dbSNP).
  #
  # chromosome --  (self explanatory)
  #
  # coordinate -- (self explanatory)
  #
  # reference_allele -- plus strand reference base
  #
  # variant_allele -- plus strand variant base
  #
  # match_status -- a Y, N or "." if a dbSNP allele, Y if the variant
  # matches the dbSNP allele, or N if it doesn't; a "." if it's a novel
  # SNP.
  #
  # rs# -- the rsid if dbSNP, "novel" otherwise.
  #
  # alternate_allele -- usually a "." (surrogate for null). A, C, T or G
  # if a third allele is seen in the reads at the given position, it's
  # listed here.  I'm don't expect you to dis play 3d allele
  # information.
  #
  # variant_count -- number of reads in which variant allele was
  # seen. Can be 1 variants matching dbSNP alleles ("Y" in match_status
  # column), must be 2 for novel alleles, for dbSNP positions that don't
  # match the dbSNP alleles ("N" in match_status column) or for dbSNP
  # positions where there is an alternate allele.
  #
  # alternate_allele_count -- number of reads in which an
  # alternate_allele is seen. Generally these are seen in only one read
  # and are probably errors, and should not be mentioned. I n some rare
  # instances (134 times), both the variant allele and the alternate
  # allele are seen multiple times.
  #
  # total_coverage -- the total number of reads at a given SNP position.
  #
  # "genotype" -- "het" if the reference allele is seen at least
  # once. "." (null) if not. These are the sites that are confidently
  # heterozygotes. The others provisionally homozygote s, and in cases
  # where the coverage is deep enough probably they are.
  if (! $self->{reader}) {
    my @field_names = qw(id chromosome coordinate reference_seq
			 variant_seq match_status rsid alternate_seq
			 variant_count alternate_seq_count
			 total_coverage zygosity);
    my $reader = GAL::Reader::DelimitedLine->new(field_names => \@field_names);
    $self->{reader} = $reader;
  }
  return $self->{reader};
}

#-----------------------------------------------------------------------------

=head1 DIAGNOSTICS

L<GAL::Parser::watson_cshl> does not throw any warnings or errors.

=head1 CONFIGURATION AND ENVIRONMENT

L<GAL::Parser::watson_cshl> requires no configuration files or
environment variables.

=head1 DEPENDENCIES

L<GAL::Parser>
L<GAL::Reader::DelimitedLine>

=head1 INCOMPATIBILITIES

None reported.

=head1 BUGS AND LIMITATIONS

No bugs have been reported.

Please report any bugs or feature requests to:
barry.moore@genetics.utah.edu

=head1 AUTHOR

Barry Moore <barry.moore@genetics.utah.edu>

=head1 LICENCE AND COPYRIGHT

Copyright (c) 2010-2014, Barry Moore <barry.moore@genetics.utah.edu>.  All
rights reserved.

    This module is free software; you can redistribute it and/or
    modify it under the same terms as Perl itself (See LICENSE).

=head1 DISCLAIMER OF WARRANTY

BECAUSE THIS SOFTWARE IS LICENSED FREE OF CHARGE, THERE IS NO WARRANTY
FOR THE SOFTWARE, TO THE EXTENT PERMITTED BY APPLICABLE LAW. EXCEPT
WHEN OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS AND/OR OTHER
PARTIES PROVIDE THE SOFTWARE "AS IS" WITHOUT WARRANTY OF ANY KIND,
EITHER EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE. THE ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE
SOFTWARE IS WITH YOU. SHOULD THE SOFTWARE PROVE DEFECTIVE, YOU ASSUME
THE COST OF ALL NECESSARY SERVICING, REPAIR, OR CORRECTION.

IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING
WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MAY MODIFY AND/OR
REDISTRIBUTE THE SOFTWARE AS PERMITTED BY THE ABOVE LICENCE, BE LIABLE
TO YOU FOR DAMAGES, INCLUDING ANY GENERAL, SPECIAL, INCIDENTAL, OR
CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OR INABILITY TO USE THE
SOFTWARE (INCLUDING BUT NOT LIMITED TO LOSS OF DATA OR DATA BEING
RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR THIRD PARTIES OR A
FAILURE OF THE SOFTWARE TO OPERATE WITH ANY OTHER SOFTWARE), EVEN IF
SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH
DAMAGES.

=cut

1;
