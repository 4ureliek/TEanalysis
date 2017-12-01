package GAL::Parser::venter_snp;

use strict;
use warnings;
use vars qw($VERSION);

$VERSION = 0.2.0;
use base qw(GAL::Parser);
use GAL::Reader::DelimitedLine;

=head1 NAME

GAL::Parser::venter_snp - Parse SNP files from the Craig Venter genome

=head1 VERSION

This document describes GAL::Parser::venter_snp version 0.2.0

=head1 SYNOPSIS

    use GAL::Parser::venter_snp
    #my $parser = GAL::Parser::venter_snp->new(file => 'venter_snp.gff');
    my $parser = GAL::Parser::venter_snp->new(file => $venter_snp_gff3_file, fasta => $fasta_file);

    while (my $feature_hash = $parser->next_feature_hash) {
	print $parser->to_gff3($feature_hash) . "\n";
    }

=head1 DESCRIPTION

L<GAL::Parser::venter_snp> parses SNP files from the Craig Venter
genome published by Levy, et al. 2007
(http://www.ncbi.nlm.nih.gov/pubmed/17803354).

SNP stands for single nucleotide polymorphism

The file mentioned is tab-delimited file of the following fileds:(Chromosome; Variant identifier; Variant Type; Chromosome start position; Chromosome end position; Not used; strand; (alleles, RepeatMasker status, TandemRepeat status); post processing status)

=head1 Constructor

New L<GAL::Parser::venter_snp> objects are created by the class method
new.  Arguments should be passed to the constructor as a list (or
reference) of key value pairs.  All attributes of the
L<GAL::Parser::venter_snp> object can be set in the call to new. An
simple example of object creation would look like this:

    my $parser = GAL::Parser::venter_snp->new(file => 'venter_snp');

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
     Usage   : GAL::Parser::venter_snp->new();
     Function: Creates a GAL::Parser::venter_snp object;
     Returns : A GAL::Parser::venter_snp object
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

# 1.  Chromosome (NCBI 36)
#
# 2.  Variant identifier (unique to each variant)
#
# 3.  Variant Type
#         There are 7 different variant types in this file.  See  Fig. 4
#         in Levy, et al. for a description of the variant types.
#         "heterozygous_mixed_sequence_variant" in this file corresponds
#         to "complex" in the figure.
#
# 4.  Chromosome start position (NCBI 36, in space-based coordinates)
#
# 5.  Chromosome end position (NCBI 36, space-based coordinates)
#
# 6.  Not used.
#
# 7.  Orientation with respect to NCBI.  The sequence is give in column.
#     If the value in column 7 is "+" or "-", this means that the sequence
#     corresponds to the positive or negative strand of NCBI, respectively.
#     A "." is given for variants that had ambiguous mapping.
#
# 8.  This field is delimited by semicolons. In the first field, the
#     alleles of the variant are given (e.g. A/B).  For homozygous
#     variants, the first allele A matches NCBI 36, the second allele B
#     matches HuRef.  There can be more than 2 alleles, which may occur
#     when reads may pile up in repititive regions.
#
#     The second field "RMR" indicates RepeatMasker status.   RMR=1
#     indicates that the variant occurs in a region identified by
#     RepeatMasker; RMR=0 indicates that the variant does not.
#
#     The third field "TR" indicates TandemRepeat status. TR=1 indicates
#     that the variant occurs in a tandem repeat, as identified by
#     TandemRepeatFinder.  TR=0 indicates that the variant does not.
#
#
# 9.  This column indicates whether a variant was post-processed.
#     Method1 indicates the variant was kept in its original form and
#     not post-processed.  Method1_MSV_clean corresponds to a modified
#     variant call where heterozygous mixed sequence variants were
#     changed to indels.  Method2 indicates that the variant is composed of a
#     cluster of Method1 variants that were within 10 bp of each other
#     (see paper for procedure).  Method2_AmbiguousMapping indicates
#     that after the clustering, the position of the new variant could
#     not be easily placed on NCBI 36, and there may be some error in
#     the mapping position given.
#
# *Please note that the values provided in Table 4 of Levy, et al. are based on Method1
# variants.

# 1       1103675000194   homozygous_SNP  815416  815417  .       +       A/G;RMR=1;TR=0  Method1
# 1       1103675000195   heterozygous_SNP        817115  817116  .       +       A/C;RMR=1;TR=0  Method1

	# Headers use by parser
	# chromosome variant_id variant_type start end orientation seqs processing

	# Fill in the first 8 columns for GFF3
	my $id         = $record->{variant_id};
	my $seqid      = 'chr' . $record->{chromosome};
	my $source     = 'Celera';
	my $type       = ''; # Assigned below
	my $start      = ++$record->{start};
	my $end        = $record->{end};
	my $score      = '.';
	my $strand     = $record->{orientation};
	my $phase      = '.';


	my ($their_zygosity, $variant_type) =
	  $record->{variant_type} =~ /(.*?)_(.*)/;

	#return undef unless $variant_type eq 'SNP';

	my %type_map = (deletion               => 'nucleotide_deletion',
			insertion              => 'nucleotide_insertion',
			mixed_sequence_variant => 'sequence_alteration',
			MNP                    => 'MNP',
			SNP                    => 'SNV',
		       );

	$type = $type_map{$variant_type};

	my $reference_seq = uc $self->fasta->seq($seqid, $start, $end);
	$reference_seq = $self->revcomp($reference_seq) if $strand eq '-';

	my ($seq_text) = split /;/, $record->{seqs};
	my @variant_seqs = split m|/|, $seq_text;
	shift @variant_seqs if $their_zygosity eq 'homozygous';

	# if ($reference_seq ne $real_ref) {
	#     $self->warn('reference_seq_mismatch',
	#		('The reference sequence given for this record ' .
	#		'does not match the reference sequence in the ' .
	#		'given fasta files.'
	#		)
	#		);
	# }

	#unshift @variant_seqs, $reference_seq
	#  if $record->{variant_type} =~ /^heterozygous/;

	my $zygosity = scalar @variant_seqs > 1 ? 'heterozygous' : 'homozygous';

	# sort | uniq -c | sort -nr
	# 1624998 heterozygous_SNP
	# 1450860 homozygous_SNP
	#  128111 heterozygous_insertion
	#   92647 heterozygous_deletion
	#   22525 heterozygous_MNP
	#   21480 heterozygous_mixed_sequence_variant
	#   14838 homozygous_MNP

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
    my @field_names = qw(chromosome variant_id variant_type start end score
			  orientation seqs processing);
    my $reader = GAL::Reader::DelimitedLine->new(field_names => \@field_names);
    $self->{reader} = $reader;
  }
  return $self->{reader};
}

#-----------------------------------------------------------------------------

=head1 DIAGNOSTICS

=over 4

=item reference_seq_mismatch

The reference sequence given for this record does not match the
reference sequence in the given fasta files.

=back

=head1 CONFIGURATION AND ENVIRONMENT

L<GAL::Parser::venter_snp> requires no configuration files or
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
