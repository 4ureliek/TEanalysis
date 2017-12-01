package GAL::Parser::venter_indel;

use strict;
use warnings;
use vars qw($VERSION);

$VERSION = 0.2.0;
use base qw(GAL::Parser);
use GAL::Reader::DelimitedLine;

=head1 NAME

GAL::Parser::venter_indel - Parse indel file from the Craig Venter
genome

=head1 VERSION

This document describes GAL::Parser::venter_indel version 0.2.0

=head1 SYNOPSIS

    use GAL::Parser::venter_indel;
    my $parser = GAL::Parser::venter_indel->new(file => 'venter_indel.gff');
    #my $parser = GAL::Parser::venter_indel->new(file => $venter_indel_gff3_file, fasta => $fasta_file);
    while (my $feature_hash = $parser->next_feature_hash) {
	print $parser->to_gff3($feature_hash) . "\n";
    }

=head1 DESCRIPTION

L<GAL::Parser::venter_indel> parses indel files from the Craig Venter
genome published by Levy, et al. 2007
(http://www.ncbi.nlm.nih.gov/pubmed/17803354)

indel is the insertion or deletion changes of sequences

The file mentioned is tab-delimited file of the following fileds: (chromosome; variant_id; variant_type; start; end; score; strand; phase; null; seq; zygosity)

=head1 Constructor

New L<GAL::Parser::venter_indel> objects are created by the class
method new.  Arguments should be passed to the constructor as a list
(or reference) of key value pairs.  All attributes of the
L<GAL::Parser::venter_indel> object can be set in the call to new. An
simple example of object creation would look like this:

    my $parser = GAL::Parser::venter_indel->new(file => 'venter_indel.gff');

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
     Usage   : GAL::Parser::venter_indel->new();
     Function: Creates a GAL::Parser::venter_indel object;
     Returns : A GAL::Parser::venter_indel object
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

	# 1 1104685014413 homozygous_indel 714051 714051 . + . . tccat Homozygous_Insertion
	# 1 1104685097444 homozygous_indel 747705 747740 . + . . CCTGGCCAGCAGATCCACCCTGTCTATACTACCTG Homozygous_Deletion
	# 1 1104685097445 homozygous_indel 751820 751820 . + . . T Homozygous_Insertion
	# 1 1104685097447 homozygous_indel 758024 758024 . + . . gtttt Homozygous_Insertion
	# 1 1104685097448 homozygous_indel 764762 764804 . + . . CACACACACCTGGACACACACACGTAGACACACACACCTAGA Homozygous_Deletion
	# 1 1104685097449 homozygous_indel 765122 765122 . + . . gaaa Homozygous_Insertion
	# 1 1104685097450 homozygous_indel 765666 765667 . + . . A Homozygous_Deletion
	# 1 1104685097451 homozygous_indel 768169 768169 . + . . CT Homozygous_Insertion
	# 1 1104685097452 homozygous_indel 778884 778933 . + . . AAACTGATGAACCCCGACCCTGATGAACGTGAGATGACCGCCGTGTGGT Homozygous_Deletion


	# $self->fields([qw(chromosome variant_id variant_type start end score
	# strand phase null seq zygosity)]);

	my $type = $record->{variant_type};
	my ($zygosity, $indel_type) = split /_/, $record->{zygosity};

	# Fill in the first 8 columns for GFF3
	my $id         = $record->{variant_id};
	my $seqid      = 'chr' . $record->{chromosome};
	my $source     = 'Levy_07';
	my $score      = '.';
	my $strand     = $record->{strand};
	my $phase      = '.';

	my ($start, $end);
	if ($indel_type eq 'Insertion') {
		$start = $record->{start};
		$end   = $record->{end};
	}
	else {
		$start = $record->{start} + 1;
		$end   = $record->{end};
	}

	my ($reference_seq, $variant_seq);
	if ($indel_type eq 'Deletion') {
	    $reference_seq = uc $record->{seq};
	    $variant_seq   = '-';
	    my $ref_fasta = uc $self->fasta->seq($seqid, $start, $end);
	    if ($ref_fasta ne $reference_seq) {
		$reference_seq = $self->revcomp($reference_seq);
	    }
	    if ($ref_fasta ne $reference_seq) {
		my $message = join '', ('The reference sequence given in the ' .
					'file does not match the reference '   .
					'fasta sequence'
					);
		my $code = 'reference_sequence_mismatch';
		my $message = "$ref_fasta | $reference_seq";
		$self->warn($code, $message);
	    }
	    $reference_seq = length($reference_seq) >= 50 ? '~' : $reference_seq;
	}
	elsif ($indel_type eq 'Insertion') {
	    $reference_seq = '-';
	    $variant_seq   = uc $record->{seq};
	    $variant_seq = $self->revcomp($variant_seq) if $strand == '-';
	    $strand = '+';
	}
	elsif ($type eq 'assembly_comparison_inversion') {
	    $reference_seq = uc $self->fasta->seq($seqid, $start, $end);
	    $reference_seq = length($reference_seq) >= 50 ? '~' : $reference_seq;
	    $variant_seq = $self->revcomp($reference_seq);
	    $variant_seq = length($variant_seq) >= 50 ? '~' : $variant_seq;
	    $strand = '+';
	}
	else {
	    my $message = "Unable to determine variant type";
	    my $code = "FATAL : unable_to_determine_variant_type : $type";
	    $self->throw($code, $message);
	}

	$zygosity = lc $zygosity;
	$zygosity ||= 'homozygous';

	$type = $type eq 'homozygous_indel' ? $indel_type : $type;
	my %type_map = ('Deletion'			=> 'deletion',
			'Insertion'			=> 'insertion',
			'assembly_comparison_inversion' => 'inversion',
			);
	$type = $type_map{$type};

	my $attributes = {Reference_seq => [$reference_seq],
			  Variant_seq   => [$variant_seq],
			  ID            => [$id],
			  Zygosity      => [$zygosity],
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
    my @field_names = qw(chromosome variant_id variant_type start end score
			  strand phase null seq zygosity);
    my $reader = GAL::Reader::DelimitedLine->new(field_names => \@field_names);
    $self->{reader} = $reader;
  }
  return $self->{reader};
}

#-----------------------------------------------------------------------------

=head1 DIAGNOSTICS

L<GAL::Parser::venter_indel> does not throw any warnings or errors.

=head1 CONFIGURATION AND ENVIRONMENT

L<GAL::Parser::venter_indel> requires no configuration files or
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

Copyright (c) 2012, Barry Moore <barry.moore@genetics.utah.edu>.  All
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
