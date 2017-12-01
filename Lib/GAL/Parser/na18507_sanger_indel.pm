package GAL::Parser::na18507_sanger_indel;

use strict;
use vars qw($VERSION);

$VERSION = 0.2.0;
use base qw(GAL::Parser);
use GAL::Reader::DelimitedLine;

=head1 NAME

GAL::Parser::na18507_sanger_indel - Parse Sanger Center indel files
for NA18507

=head1 VERSION

This document describes GAL::Parser::na18507_sanger_indel version 0.2.0

=head1 SYNOPSIS

    my $parser = GAL::Parser::na18507_sanger_indel->new(file =>
             'na18507_sanger_indel.txt');

    while (my $feature_hash = $parser->next_feature_hash) {
	print $parser->to_gff3($feature_hash) . "\n";
    }

=head1 DESCRIPTION

L<GAL::Parser::na18507_sanger_indel> provides a parser for the Sanger
Center's indel data for NA18507.

=head1 Constructor

New L<GAL::Parser::na18507_sanger_indel> objects are created by the
class method new.  Arguments should be passed to the constructor as a
list (or reference) of key value pairs.  All attributes of the
L<GAL::Parser::na18507_sanger_indel> object can be set in the call to
new. An simple example of object creation would look like this:

    my $parser = GAL::Parser::na18507_sanger_indel->new(file =>
             'na18507_sanger_indel.txt');

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
     Usage   : GAL::Parser::na18507_sanger_indel->new();
     Function: Creates a GAL::Parser::na18507_sanger_indel object;
     Returns : A GAL::Parser::na18507_sanger_indel object
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
	my @valid_attributes = qw(); # Set valid class attributes here
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

	# 1       709353  */+ATAAT        16
	# 1       741363  */+AT   45
	# 1       757983  */-TTG  13
	# 1       768166  +CT/+CT 373
	# 1       770921  -C/*    723
	# 1       776933  */-AG   114

	# $self->fields([qw(chr pos ref_base con_base con_qual read_depth ave_hits_elsewhere)]);

	my @variant_seqs = grep {$_ ne '*'} split m|/|, $record->{seq_text};

	# Ignore records that are compound heterozygous.
	return undef if (scalar @variant_seqs == 2 && $variant_seqs[0] ne $variant_seqs[1]);

	my ($zygosity, $reference_seq, $type, $start, $end);

	$zygosity = scalar @variant_seqs == 1 ? 'heterozygous' : 'homozygous';

	for my $variant_seq (@variant_seqs) {
		if ($variant_seq =~ /^-/) {
			$type = 'nucleotide_deletion';
			($reference_seq = $variant_seq) =~ s/^-//;
			$variant_seq = '-';
			$start = $record->{location};
			$end   = $record->{location} + length($reference_seq) - 1;
		}
		else {
			$type = 'nucleotide_insertion';
			$variant_seq =~ s/^\+//;
			$reference_seq = '-';
			$start = $record->{location};
			$end   = $record->{location};
		}
	}

	my $seqid      = 'chr' . $record->{chr};
	my $source     = 'NA18507_Sanger';
	my $score      = $record->{score};
	my $strand     = '+';
	my $phase      = '.';
	my $id         = join ':', ($seqid, $source, $type, $start);

	my $attributes = {Reference_seq => [$reference_seq],
			  Variant_seq   => \@variant_seqs,
			  Zygosity      => [$zygosity],
			  ID            => [$id],
			 };

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

  if (! $self->{reader}) {
    my @field_names = qw(chr location seq_text score);
    my $reader = GAL::Reader::DelimitedLine->new(field_names => \@field_names);
    $self->{reader} = $reader;
  }
  return $self->{reader};
}

#-----------------------------------------------------------------------------

=head1 DIAGNOSTICS

L<GAL::Parser::na18507_sanger_indel> does not throw any warnings or errors.

=head1 CONFIGURATION AND ENVIRONMENT

L<GAL::Parser::na18507_sanger_indel> requires no configuration files or environment variables.

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
