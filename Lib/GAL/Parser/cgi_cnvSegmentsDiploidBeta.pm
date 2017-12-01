package GAL::Parser::cgi_cnvSegmentsDiploidBeta;

use strict;
use vars qw($VERSION);
$VERSION = 0.2.0;

use base qw(GAL::Parser);
use GAL::Reader::DelimitedLine;

=head1 NAME

GAL::Parser::cgi_cnvSegmentsDiploidBeta - Parse Complete Genomics
cnvSegmentsDiploidBeta files

=head1 VERSION

This document describes GAL::Parser::cgi_cnvSegmentsDiploidBeta version 0.2.0

=head1 SYNOPSIS

    my $parser = GAL::Parser::cgi_cnvSegmentsDiploidBeta->new(file => 'cgi_cnvSegmentsDiploidBeta.tsv');

    while (my $feature_hash = $parser->next_feature_hash) {
	print $parser->to_gff3($feature_hash) . "\n";
    }

=head1 DESCRIPTION

L<GAL::Parser::cgi_cnvSegmentsDiploidBeta> provides a parser for
cnvSegmentsDiploidBeta.tsv files.

=head1 Constructor

New L<GAL::Parser::cgi_cnvSegmentsDiploidBeta> objects are created by the class method
new.  Arguments should be passed to the constructor as a list (or
reference) of key value pairs.  All attributes of the
L<GAL::Parser::cgi_cnvSegmentsDiploidBeta> object can be set in the call to new. An
simple example of object creation would look like this:

    my $parser = GAL::Parser::cgi_cnvSegmentsDiploidBeta->new(file => 'cgi_cnvSegmentsDiploidBeta.txt');

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
     Usage   : GAL::Parser::cgi_cnvSegmentsDiploidBeta->new();
     Function: Creates a GAL::Parser::cgi_cnvSegmentsDiploidBeta object;
     Returns : A GAL::Parser::cgi_cnvSegmentsDiploidBeta object
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
	my ($self, $record) = @_;

	# $record is a hash reference that contains the keys assigned
	# in the $self->fields call in _initialize_args above

	#chr begin end avgNormalizedCvg relativeCvg calledPloidy calledCNVType ploidy Score CNV Type Score overlappingGene knownCNV repeats
	# chr1	521368 712000  58.2  1.11  N  hypervariable  0	0
	# chr1	712000 724000  51.5  0.98  2  =              45  45
	# chr1	17206000	17270000	78.6	1.50	3	+	41	48	CROCC;LOC100506687;RNU1-2	dgv.1:Variation_2301;dgv.1:Variation_3284;dgv.1:Variation_4210;dgv.1:Variation_4211;dgv.1:Variation_4212;dgv.2:Variation_6797;dgv.4:Variation_31610;dgv.5:Variation_34489;dgv.6:Variation_38949;dgv.6:Variation_38959;dgv.9:Variation_74226;dgv.9:Variation_74388;dgv.9:Variation_74390;dgv.9:Variation_84102	DNA:2;LINE:8;LTR:14;Low_complexity:1;SINE:29;SegDup:100;SelfChain:100;Simple_repeat:1;snRNA:1;tRNA:1
	# chr1	72766000	72812000	0.6	0.01	0	-	71	72		dgv.1:Variation_0383;dgv.1:Variation_1543;dgv.1:Variation_3293;dgv.1:Variation_4229;dgv.2:Variation_6832;dgv.3:Variation_22923;dgv.3:Variation_23231;dgv.4:Variation_30388;dgv.4:Variation_31627;dgv.5:Variation_34829;dgv.5:Variation_34830;dgv.5:Variation_34831;dgv.5:Variation_34832;dgv.5:Variation_34833;dgv.6:Variation_37834;dgv.7:Variation_43688;dgv.8:Variation_58405;dgv.9:Variation_60842;dgv.9:Variation_64030;dgv.9:Variation_74531;dgv.9:Variation_84304;dgv.9:Variation_97396	DNA:2;LINE:40;LTR:19;Low_complexity:1;SINE:8;SelfChain:1;Simple_repeat:2

	#------------------------------
	#  chr
	#  Index: 0
	#  Description:  Chromosome name.  The mitochondrion is
	#  		 represented as chrM, though this may be
	#  		 absent from CNV analyses. The pseudoautosomal
	#  		 regions within the sex chromosomes X and Y
	#  		 are reported at their coordinates on
	#  		 chromosome X
	#  Format: TEXT
	#  Example: chr1
	#-----------------------------
    	#  begin
	#  Index: 1
	#  Description: Beginning of segment.
	#  Format: INT
	#  Example: 72766000
	#-----------------------------
	#  end
	#  Index: 2
	#  Description: End of segment
	#  Format: INT
	#  Example: 72812000
	#-----------------------------
	#  avgNormalizedCvg
	#  Index: 3
	#  Description: Baseline-normalized average coverage over the
	#  		interval from begin to end. avgNormalizedCvg
	#  		is no-) in regions of the genome with very
	#  		high or very low GC content as well as in
	#  		regions with very low average coverage among
	#  		the baseline samples. See further description
	#  		in 'Detailed Ploidy and Coverage Information'
	#  Format: FLOAT
	#  Example:  58.2
	#-----------------------------
	#  relativeCvg
	#  Index: 4
	#  Description: avgNormalizedCvg divided by estimate of
	#  		diploid median average adjusted coverage. Value if
	#  		avgNormalizedCvg .Nis Nis
	#  Format: FLOAT
	#  Example: 7.13
	#-----------------------------
	#  calledPloidy
	#  Index: 5
	#  Description: Called ploidy for the segment. Typically an
	#  		integer in the range [0,1,...,MAX_PLOIDY]; 'N'
	#  		when calledCNVType is invariant,
	#  		hypervarialble or 'N'.
	#  Format: TEXT
	#  Example: 0, 1, N
	#-----------------------------
	#  calledCNVType
	#  Index: 6
	#  Description: Classification of called ploidy to one of six categories
	#               - (hyphen): a reduction in copy number relative to the nominal expectation (diploid for autosomes, sex-appropriate for sex chromosomes).
	#               = (equals): a match to the nominal expectation.
	#               + (plus): an increase relative to the nominal expectation.
	#               hypervariable: coverage not interpretable as a discrete ploidy due to high diversity of coverage levels in the sequenced genome and a set  genomes.
	#               N: whole genome coverage has no-;
	#  Format: TEXT
	#-----------------------------
	#  ploidyScore
	#  Index: 7
	#  Description: Phred-like confidence that the segment has the called ploidy.
	#  Format: INT
	#  Example: 20
	#-----------------------------
	#  CNVTypeScore
	#  Index: 8
	#  Description: Phred-like confidence that the calledCNVType is correct.
	#  Format: INT
	#  Example: 20
	#-----------------------------
	#  overlappingGene
	#  Index: 9
	#  Data Type: TEXT
	#  Description: Gene(s) overlapping called segment, with minimum overlap of a single base pair.
	#  Format: CROCC;LOC100506687;RNU1-2
	#-----------------------------
	#  knownCNV
	#  Index: 10

	#  Description: Known CNVs in the Database for Genomic
	#  Variants that overlap called segment. Overlap requires that
	#  the CNV segment in DGV covers at least 80% of Complete
	#  Genomics called CNV segment, allowing a single-window error
	#  in the boundary on each side of the called segment.
	#  Format: TEXT
	#  Example: dgv.1:Variation_3310;dgv.1:Variation_4247
	#-----------------------------
	#  repeats
	#  Index: 11
	#  Data Type: TEXT
	#  Description: Percent of called CNV segment that overlaps
	#  		with each category of genomic
	#  		repeats. Categories include: DNA, LINE,
	#  		Low_Complexity, SINE, Satellite, SegDup,
	#  		Self-chain, Simple_Repeats, scRNA, tRNA, and
	#  		snRNA. If the amount of overlap for a category
	#  		is less than 1%, category is not reported.
	#  Format: TEXT
	#  Example:
	#------------------------------

	# Fill in the first 8 columns for GFF3
	# See http://www.sequenceontology.org/gff3.html for details.

	return undef unless $record->{calledCNVType} =~ /^[+-]$/;

	#chr begin end avgNormalizedCvg relativeCvg calledPloidy calledCNVType ploidy Score CNV Type Score overlappingGene knownCNV repeats
	my $seqid      = $record->{chr};
	my $source     = 'cgiCNV';
	my $type       = ($record->{calledCNVType} eq '+' ? 'copy_number_gain' :
			  $record->{calledCNVType} eq '-' ? 'copy_number_loss' :
			  'copy_number_variation');
	my $start      = $record->{begin};
	my $end        = $record->{end};
	my $score      = $record->{Score};
	my $strand     = '+';
	my $phase      = '.';
	my $feature_id = "cgiCNV_$seqid:$start-$end";

	# Create the attributes hash
	# See http://www.sequenceontology.org/gvf.html

	# Assign the reference and variant sequences:
	# Reference_seq=A
	# The order of the Variant_seqs matters for indexing Variant_effect
	# - see below
	# Variant_seq=G,A
	my $variant_seqs = '~';

	# If you need to look up the reference sequence
	# Then be sure that you pass (fasta => '/path/to/fasta') as an
	# argument when instantiating the class.
	my $reference_seq = uc $self->fasta->seq($seqid, $start, $end);
	$reference_seq = '~' if length $reference_seq > 25;

	my $attributes = {Reference_seq       => [$reference_seq],
			  Variant_seq         => [$variant_seqs],
			  ID                  => [$feature_id],
			  cgi_relativeCvg     => [$record->{relativeCvg}],
			  cgi_calledPloidy    => [$record->{calledPloidy}],
			  cgi_calledCNVType   => [$record->{calledCNVType}],
			  avgNormalizedCvg    => [$record->{avgNormalizedCvg}],
			  cgi_ploidyScore     => [$record->{ploidyScore}],
			  cgi_CNVTypeScore    => [$record->{CNVTypeScore}],
			  cgi_overlappingGene => [$record->{overlappingGene}],
			  cgi_knownCNV        => [$record->{knownCNV}],
			  cgi_repeats         => [$record->{repeats}],
	};

	map {delete $attributes->{$_} unless length $attributes->{$_}[0]} keys %{$attributes};
        for my $values (values %{$attributes}) {
            map {s/;/,/g} @{$values};
	    }

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

	return $feature_data;
}

#-----------------------------------------------------------------------------

=head2 parse_metadata

 Title   : parse_metadata
 Usage   : $a = $self->parse_metadata
 Function: Call the reader to return metadata and then parse and store it.
 Returns : N/A
 Args    : A hash reference of parse_metadata

=cut

sub parse_metadata {
    my $self = shift;

    my $raw_data_array = $self->reader->metadata;

    my %raw_data;
    for my $line (@{$raw_data_array}) {
	my ($key, $value) = split /\s+/, $line, 2;
	$raw_data{$key} = $value;
    }

    #ASSEMBLY_ID GS000007756-ASM
    #FORMAT_VERSION 2.0
    #GENERATED_AT 2012-Feb-18 21:32:18.134149
    #GENERATED_BY CallCNVs
    #GENOME_REFERENCE NCBI build 37
    #MAX_PLOIDY 10
    #SAMPLE GS01003-DNA_A01
    #SOFTWARE_VERSION 2.0.2.20
    #TYPE CNV-SEGMENTS
    #WINDOW_SHIFT 2000
    #WINDOW_WIDTH 2000
    #GENE_ANNOTATIONS NCBI build 37.2
    #DGV_VERSION 9

    my %metadata;
    $metadata{'gvf-version'} = '1.06';
    $metadata{'feature-ontology'} = 'http://sourceforge.net/projects/song/files/SO_Feature_Annotation/sofa_2_4_1/sofa_2_4_1.obo/download';
    $metadata{'file-version'} = '1.0';
    $metadata{'individual-id'} = $raw_data{ASSEMBLY_ID};
    $metadata{'technology-platform-class'} = 'SRS';
    $metadata{'technology-platform-name'} = 'Complete_Genomics';
    $metadata{'technology-platform-version'} = '2.0.2.20';
    $metadata{'sequencing-scope'} = 'whole_genome';
    $metadata{'sequence-alignment'} = 'Complete Genomics analysis pipeline';
    $metadata{'variant-calling'} =
	('CNVs called with ' . $raw_data{GENERATED_BY} . ' ' .
	 'version ' . $raw_data{SOFTWARE_VERSION} . ' ' .
	 'with parameters' .
	 'WINDOW_SHIFT=' . $raw_data{WINDOW_SHIFT} . ', ' .
	 'WINDOW_WIDTH=' . $raw_data{WINDOW_WIDTH} . ', ' .
	 'GENE_ANNOTATIONS=' . $raw_data{GENE_ANNOTATIONS} . ', ' .
	 'DGV_VERSION=' . $raw_data{DGV_VERSION});
    $metadata{'sample-description'} = $raw_data{SAMPLE};
    $metadata{'genomic-source'} = 'somatic';
    $self->{metadata} = \%metadata;
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

    #chr begin end avgNormalizedCvg relativeCvg calledPloidy calledCNVType ploidyScore CNVTypeScore overlappingGene knownCNV repeats

    my @field_names = qw(chr begin end avgNormalizedCvg relativeCvg
    			 calledPloidy calledCNVType ploidyScore CNVTypeScore
    			 overlappingGene knownCNV repeat);

    my $reader = GAL::Reader::DelimitedLine->new(field_names     => \@field_names,
						 field_separator => "\t",
						 comment_pattern => qr/^(\#|>)/,
	);
    $self->{reader} = $reader;
  }
  return $self->{reader};
}

#-----------------------------------------------------------------------------

=head1 DIAGNOSTICS

L<GAL::Parser::cgi_cnvSegmentsDiploidBeta> does not throw any warnings or errors.

=head1 CONFIGURATION AND ENVIRONMENT

L<GAL::Parser::cgi_cnvSegmentsDiploidBeta> requires no configuration files or
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
