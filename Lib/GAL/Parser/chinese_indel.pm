package GAL::Parser::chinese_indel;

use strict;
use vars qw($VERSION);


$VERSION = 0.2.0;
use base qw(GAL::Parser);
use GAL::Reader::DelimitedLine;

=head1 NAME

GAL::Parser::chinese_indel - Parse indel files from the first Chinese
genome

=head1 VERSION

This document describes GAL::Parser::chinese_indel version 0.2.0

=head1 SYNOPSIS

    my $parser = GAL::Parser::chinese_indel->new(file =>
             'chinese_indel.gff');

    while (my $feature_hash = $parser->next_feature_hash) {
	print $parser->to_gff3($feature_hash) . "\n";
    }

=head1 DESCRIPTION

L<GAL::Parser::chinese_indel> provides a parser for indel data from
the first Chinese genome to be published by Wang, et al. 2008
(http://www.ncbi.nlm.nih.gov/pubmed/18987735).

=head1 Constructor

New L<GAL::Parser::chinese_indel> objects are created by the class
method new.  Arguments should be passed to the constructor as a list
(or reference) of key value pairs.  All attributes of the
L<GAL::Parser::chinese_indel> object can be set in the call to new. An
simple example of object creation would look like this:

    my $parser = GAL::Parser::chinese_indel->new(file =>
             'chinese_indel.gff');

The constructor recognizes the following parameters which will set the
appropriate attributes:

=over 4

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
     Usage   : GAL::Parser::chinese_indel->new();
     Function: Creates a GAL::Parser::chinese_indel object;
     Returns : A GAL::Parser::chinese_indel object
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

	# $record is a hash reference that contains the keys assigned
	# in the $self->fields call in _initialize_args above

	# Fill in the first 8 columns for GFF3
	# See http://www.sequenceontology.org/resources/gff3.html for details.
	my $original_atts = $self->parse_attributes($record->{attributes});

	my $id         = $original_atts->{ID}[0];
	my $seqid      = $record->{seqid};
	my $source     = $record->{source};
	my $type       = $record->{type};
	my $start      = $record->{start} + 1; # Soap indels are space based
	my $end        = $record->{end};
	my $score      = $record->{score};
	my $strand     = $record->{strand};
	my $phase      = $record->{phase};

	# chr1 soap Indel   1409   14094 + . ID=YHIndel00001; status=novel; Type=+1; location=Transposons0034384:LINE/L1; Base=A
	# chr1 soap Indel   3824   38253 + . ID=YHIndel00002; status=novel; Type=-1; Base=C
	# chr1 soap Indel  29355  293553 + . ID=YHIndel00003; status=novel; Type=+1; location=Transposons0167541:LTR/MaLR; Base=C
	# chr1 soap Indel  39377  393784 + . ID=YHIndel00004; status=novel; Type=-1; location=Transposons0034394:LINE/L1; Base=G
	# chr1 soap Indel  53601  536043 + . ID=YHIndel00005; status=novel; Type=-3; Base=CTA
	# chr1 soap Indel 241491 2414923 + . ID=YHIndel00006; status=novel; Type=-1; Base=C
	# chr1 soap Indel 429836 4298384 + . ID=YHIndel00007; status=novel; Type=-2; Base=CC

	# Create the attributes hash

	my $indel_type = $original_atts->{Type}[0];
	my ($reference_seq, $variant_seq);
	if ($indel_type > 0) {
		$reference_seq = '-';
		$variant_seq   = $original_atts->{Base}[0];
		$type             = 'nucleotide_insertion';
		$start            = $record->{start} - 1;
		$end              = $record->{start} - 1;
	}
	elsif ($indel_type < 0) {
		$reference_seq = $original_atts->{Base}[0];
		$variant_seq   = '-';
		$type             = 'nucleotide_deletion';
		$start            = $record->{start} + 1;
		$end              = $record->{end};
	}

	my $attributes = {Reference_seq => [$reference_seq],
			  Variant_seq   => [$variant_seq],
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
 Returns : A GAL::Reader::DelimitedLine singleton.
 Args    : None

=cut

sub reader {
  my $self = shift;

  if (! $self->{reader}) {
    my @field_names = qw(seqid source type start end score strand phase attributes);
    my $reader = GAL::Reader::DelimitedLine->new(field_names => \@field_names);
    $self->{reader} = $reader;
  }
  return $self->{reader};
}

#-----------------------------------------------------------------------------

=head1 DIAGNOSTICS

L<GAL::Parser::chinese_indel> does not throw any warnings or errors.

=head1 CONFIGURATION AND ENVIRONMENT

L<GAL::Parser::chinese_indel> requires no configuration files or
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
