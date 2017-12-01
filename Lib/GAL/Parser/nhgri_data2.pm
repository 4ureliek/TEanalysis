package GAL::Parser::nhgri_data2;

use strict;
use vars qw($VERSION);
$VERSION = 0.2.0;

use base qw(GAL::Parser);
use GAL::Reader::DelimitedLine;

=head1 NAME

GAL::Parser::nhgri_data2 - Parse NHGRI_DATA2 files

=head1 VERSION

This document describes GAL::Parser::nhgri_data2 version 0.2.0

=head1 SYNOPSIS

    my $parser = GAL::Parser::nhgri_data2->new(file => 'nhgri_data2.txt');

    while (my $feature_hash = $parser->next_feature_hash) {
	print $parser->to_gff3($feature_hash) . "\n";
    }

=head1 DESCRIPTION

L<GAL::Parser::nhgri_data2> provides a parser for NHGRI_DATA2 data.

=head1 Constructor

New L<GAL::Parser::nhgri_data2> objects are created by the class method
new.  Arguments should be passed to the constructor as a list (or
reference) of key value pairs.  All attributes of the
L<GAL::Parser::nhgri_data2> object can be set in the call to new. An
simple example of object creation would look like this:

    my $parser = GAL::Parser::nhgri_data2->new(file => 'nhgri_data2.txt');

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
     Usage   : GAL::Parser::nhgri_data2->new();
     Function: Creates a GAL::Parser::nhgri_data2 object;
     Returns : A GAL::Parser::nhgri_data2 object
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
	my @valid_attributes = qw(file fh); # Set valid class attributes here
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
	my ($self, $data) = @_;


	my %record;

	# We didn't get $data as a hash because we didn't give field names to the reader.
	# This allows us to get the data as an array and deal with multiple genotypes per line.

	# Chr LeftFlank RightFlank ref_allele muttype var_allele refseq transcript ref_aa var_aa aa_pos RS\# M87_4.NA M87_4.NA.score M87_4.NA.coverage
	# chrX 110529 110531 T SNP G NA NA NA NA 0 rs67516 229(G,T),rs6649863(G,T) GG 16 24 
	# chrX 314529 314531 G SNP A NA NA NA NA 0 -                               CG 17 15 

	@record{qw(Chr LeftFlank RightFlank ref_allele muttype var_allele refseq transcript ref_aa var_aa aa_pos RS)} = splice(@{$data}, 0, 12);

	# Ignore indels for now.
	return undef if $record{muttype} ne 'SNP';

	my $seqid      = $record{Chr};
	my $source     = 'NHGRI';
	my $type       = 'SNV';
	my $start      = $record{LeftFlank} + 1;
	my $end        = $start;
	my $feature_id = join ':', ($seqid, $source, $type, $start, $end);
	my $score      = '.';
	my $strand     = '+';
	my $phase      = '.';

	my $reference_seq = $record{ref_allele};
	my $alias_text = $record{RS} if $record{RS};
	$alias_text =~ s/\"//g;
	$alias_text =~ s/\(.*?\)//g;
	my @aliases = grep {$_ ne '-'} split /,/, $alias_text;

	my %individuals;
	my ($header) = $self->reader->headers;
	my @header_data = split /\s/, $header;
	splice(@header_data, 0, 12);
	while (@{$data}) {
	    my ($individual_id) = (shift @header_data, shift @header_data, shift @header_data);
	    $individuals{$individual_id} = [shift @{$data},   shift @{$data},   shift @{$data}];
	}

	my @features;
	for my $individual_id (sort keys %individuals) {

	    $source = $individual_id;

	    my $individual = $individuals{$individual_id};
	    my ($this_zygosity, $this_score, $coverage) = @{$individual};
	    $score = $this_score ? $this_score : $score;

	    # Ignore no-calls for now.
	    next if $this_zygosity eq 'NA';

	    my @variant_seqs = split '', $this_zygosity;
	    shift @variant_seqs if $variant_seqs[0] eq $variant_seqs[1];

	    next if (scalar @variant_seqs == 1 && $variant_seqs[0] eq $reference_seq);

	    my $zygosity = scalar @variant_seqs > 1 ? 'heterozygous' : 'homozygous';


	    my $attributes = {Reference_seq => [$reference_seq],
			      Variant_seq   => \@variant_seqs,
			      Zygosity      => [$zygosity],
			      ID            => [$feature_id],
			      Total_reads   => [$coverage],
			  };

	    push @{$attributes->{Alias}}, @aliases if @aliases;

	    my $feature_data = {feature_id => $feature_id,
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
	    push @features, $feature_data;
	}
	return scalar @features ? \@features : undef;
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

  if (! $self->{reader}) {
      
    # Change the field names in the next line to define how the hash keys
    # are named for each column in you input file.
    # If you don't pass field_names as an argument to the reader then
    # the reader will return an array rather than a hash.
    my @field_names = qw();
    my $reader = GAL::Reader::DelimitedLine->new(comment_pattern => qr/^\#/,
						 header_count    => 1,
						 );

    $self->{reader} = $reader;
  }
  return $self->{reader};
}

#-----------------------------------------------------------------------------

=head1 DIAGNOSTICS

L<GAL::Parser::nhgri_data2> does not throw any warnings or errors.

=head1 CONFIGURATION AND ENVIRONMENT

L<GAL::Parser::nhgri_data2> requires no configuration files or
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
