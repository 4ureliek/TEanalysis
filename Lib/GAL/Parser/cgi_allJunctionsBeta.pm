package GAL::Parser::cgi_allJunctionsBeta;

use strict;
use vars qw($VERSION);
$VERSION = '0.01';

use base qw(GAL::Parser);
use GAL::Reader::DelimitedLine;
use Storable qw(dclone);

=head1 NAME

GAL::Parser::cgi_allJunctionsBeta - Parse Complete Genomics
cnvSegmentsDiploidBeta files

=head1 VERSION

This document describes GAL::Parser::cgi_allJunctionsBeta version 0.01

=head1 SYNOPSIS

    my $parser = GAL::Parser::cgi_allJunctionsBeta->new(file => 'cgi_allJunctionsBeta.tsv');

    while (my $feature_hash = $parser->next_feature_hash) {
	print $parser->to_gff3($feature_hash) . "\n";
    }

=head1 DESCRIPTION

L<GAL::Parser::cgi_allJunctionsBeta> provides a parser for
cnvSegmentsDiploidBeta.tsv files.

=head1 Constructor

New L<GAL::Parser::cgi_allJunctionsBeta> objects are created by the class method
new.  Arguments should be passed to the constructor as a list (or
reference) of key value pairs.  All attributes of the
L<GAL::Parser::cgi_allJunctionsBeta> object can be set in the call to new. An
simple example of object creation would look like this:

    my $parser = GAL::Parser::cgi_allJunctionsBeta->new(file => 'cgi_allJunctionsBeta.txt');

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
     Usage   : GAL::Parser::cgi_allJunctionsBeta->new();
     Function: Creates a GAL::Parser::cgi_allJunctionsBeta object;
     Returns : A GAL::Parser::cgi_allJunctionsBeta object
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

	# 00 Id
	# 01 LeftChr
	# 02 LeftPosition
	# 03 LeftStrand
	# 04 LeftLength
	# 05 RightChr
	# 06 RightPosition
	# 07 RightStrand
	# 08 RightLength
	# 09 StrandConsistent
	# 10 Interchromosomal
	# 11 Distance
	# 12 DiscordantMatePairAlignments
	# 13 JunctionSequenceResolved
	# 14 TransitionSequence
	# 15 TransitionLength
	# 16 LeftRepeatClassification
	# 17 RightRepeatClassification
	# 18 LeftGenes
	# 19 RightGenes
	# 20 XRef
	# 21 DeletedTransposableElement
	# 22 KnownUnderrepresentedRepeat
	# 23 FrequencyInBaselineGenomeSet
	# 24 AssembledSequence
	# 25 EventId
	# 26 Type
	# 27 RelatedJunctions

	#ASSEMBLY_ID	GS000007756-ASM
	#SOFTWARE_VERSION	2.0.2.20
	#GENERATED_BY	cgatools
	#GENERATED_AT	2012-Feb-19 00:25:16.982799
	#FORMAT_VERSION	2.0
	#GENOME_REFERENCE	NCBI build 37
	#SAMPLE	GS01003-DNA_A01
	#TYPE	JUNCTIONS
	#DBSNP_BUILD	dbSNP build 132
	#GENE_ANNOTATIONS	NCBI build 37.2
	
	# >Id	LeftChr	LeftPosition	LeftStrand	LeftLength	RightChr	RightPosition	RightStrand	RightLength	StrandConsistent	Interchromosomal	Distance	DiscordantMatePairAlignments	JunctionSequenceResolved	TransitionSequence	TransitionLength	LeftRepeatClassification	RightRepeatClassification	LeftGenes	RightGenes	XRef	DeletedTransposableElement	KnownUnderrepresentedRepeat	FrequencyInBaselineGenomeSet	AssembledSequence	EventId	Type	RelatedJunctions
	# 1519	chr1	825766	-	149	chr1	5726936	-	800	Y	N	4901170	7	Y		0	AluYc:SINE:Alu;SegDup;Self chain	AluSx:SINE:Alu;SegDup;Self chain						0.92	accaaaatgactattctttctaccctcctaGTTCAGACATAGCCTGAGACTTTTTTTTTTTTGAGATGAAGTCTCACTCTGTCACTCAGGCTGGAGTGCAGTGGCATGGTCTCGGCTCATTGCAATCTCTACCTCCCGGGTTCAAGTGATTCTCCTGCCTCAGCTTCCAGAGTAGCTGGGGCTATAGGCACCGGCCACCACACCCAGCAAATTTCCGTATTATTAGTAGAGATGGGGTTTCACCATGTTGGTCAGACTGGTCTTGAACTCCTGACCTCAGGTGATCTGCCCGCCTCTGCCTCCCAAAATGCTGAGATTACAGATGTGAGCCACTgtgcccggccgcctgagacattttggacga	1152	complex	3198;3199;4309;4314;4316
	# 4316	chr1	826300	+	447	chr1	5727488	+	173	Y	N	4901188	9	Y		0	SegDup;Self chain	SegDup;Self chain						1.00	tttggttgccacaaaccatcgaaacaaagaTGCAGATCATTGATGTAAAATTACAGTTAGTTCCTTCCCACTCCTTTTCAGCTTCTCTTCGTTGCTATGACCCAGCGTCTCCTGTGTCAGTTTTCAGTCTATTGTCTCCCAGCTTCTAGTGCACCTTTCAATATGTGCACTGTGATAAACTGGGAAACACTGTTCAATATACCTTCTGGAAGTGAACATTCTGCAGGCATCTAGGTAGAGGATGGAGAGACTGCAGGGGGCAGGAGCTCTCTGGCTGGGCCTTGATCATGTTCAAGCCCCAACCACAGACCTAGGCGTGGTCCCTCAGCCACCTTGTAGCCTTGGCTTGCAACATCTCGACATGGAAACCAAAACGCTGCACGGCCCATGTGATATGAAAGTTCTTGAAAAGTTGCCCAGACCCcctcttttaccccttgtgcaacctgcacac	1152	complex	1519;3198;3199;4309;4314
	# 4314	chr1	826313	+	223	chr1	5728168	+	128	Y	N	4901855	5	Y	CTTTTGTAACCTGCACACAGTGACCTGTATTCTAGAGGGTCCACACAGAGCTGCCATTCCTTCTGCCAGACCCTGTGGGACTCAGGCATTCTGGAGGCTTCCTGCCCTACAAAGGCAGCCAGACTCCCGCCATGCATCCCTGCACCAGCGGCTCACGGCTAGCTCCCTCATCTGCATTCCAGTGGCTCACAGCCAGCTCCCTCCCCTGCATTCCAGTGGCTCATGGCCAGCTCCCTCACCTGCACCAGCGGCTCATGGCCAGCTCCCTCCCCT	273	SegDup;Self chain	SegDup;Self chain;Tandem period 100;Tandem period 34;Tandem period 66						1.00	gaggatggagagactgcagggggcaggagcTCTCTGGCTGGGCCTTGTTCAAGCCCCAACCACAGAGACCTAGGCATGGTCCCTCAGCCACCTTGTAGCCTTGGCTTGCAACATCTCGACATGGAAACCAAAACGCTGCACGGCCCATGTGATATGAAAATTCTTGAAAAGTTGCCCAGACCCCCTCTTGTACCCCTTTTGTAACCTGCACACAGTGACCTGTATTCTAGAGGGTCCACACAGAGCTGCCATTCCTTCTGCCAGACCCTGTGGGACTCAGGCATTCTGGAGGCTTCCTGCCCTACAAAGGCAGCCAGACTCCCGCCATGCATCCCTGCACCAGCGGCTCACGGCTAGCTCCCTCATCTGCATTCCAGTGGCTCACAGCCAGCTCCCTCCCCTGCATTCCAGTGGCTCATGGCCAGCTCCCTCACCTGCACCAGCGGCTCATGGCCAGCTCCCTCCCCTacactccagcggctcaccgccggctccctc	1152	complex	1519;3198;3199;4309;4316
	# 4318	chr1	869490	+	442	chr1	870330	+	317	Y	N	840	84	Y		0	Self chain;Tandem period 28	Self chain;Tandem period 28	NM_152486:+	NM_152486:+	rs70954822 (chr1:869484-870324)			0.67	gctcctaggtcactgcgcaggactgctccgTTACAGGTGGGCAGGGGAGGCTGCTCCGTTACAGGTGGGCAGGGGAGGCGGCTGCGTTACAGGTGTGCAGGGGAGGCGGCTGCGTTACAGGTGGGCAGGGGAGGCGGCTGCGTTACAGGTGGGCGGGGGAGGCGGCTGCGTTAcaggtgggcgggcggtgctgcaggaggact	3446	deletion	
	# 4317	chr1	964038	+	457	chr1	964525	+	288	Y	N	487	41	Y		0	Self chain;Tandem period 61	Self chain	NM_198576:+	NM_198576:+				0.35	gcattcgagcctgtgtgcggtgcatgggccAAGGGCGCCCACACCCACGCCACCCTTTCCGAAGGAACCGAGCCCCAGCCCCTCATGGGCCAAGGGCACCCACAGCCACGCCACCCTTTCCGAAGGAACCGAGCCCCAGCCCCTCGTGGGCCAAGGGCGCCCACAGCCACGCCACCCTCTCCCAAGGAACCGAGCCCCAGCCCCTCTGGGGCCTGCCAattgccagagagccccagtgctccacccac	3445	deletion	
	# 4315	chr1	1232627	+	342	chr1	1233165	+	84	Y	N	538	8	Y		0	Self chain;Tandem period 226;Tandem period 46	Self chain;Tandem period 226;Tandem period 46	NM_030649:-	NM_030649:-				0.52	agtgcacaccagcacacacccagcttgcagGAGCCACACACAGGCGTGTGCACGTGTGTGGGGCAGGGGCCATCCCCAGTGGCACGTGTGTGTGCACAGGCACGGGGCAGGGGCCATCCCCGGTGGCACATGTGTGCACGGGCTTGGGGCAGGGGCCATCCCCGGTGGCACGTGTGTGTGCACAGGCACGGGGCAGGGGCCATCCCCGGTGGCACGTGTGTGTGCACGGGcttggggcaggggcaccaagggccacacct	3444	deletion	
	# 4313	chr1	2911548	+	415	chr1	2911862	+	297	Y	N	314	5	Y		0	AluYg6:SINE:Alu	L1ME3A:LINE:L1				AluYg6 2.2% (chr1:2911548-2911856)		0.62	tgagggagaatgaaaacttatgtgtttaaaCCACTGTATATTTGACATTTGTCCCTCTCAAATGCTTTCAAACCGAAGCCTAATTAATACTCATGGCAACCAACCTGAATTCAGCTCTCTAATGTTAAACCATCCTTGCATTCCTCTGATAAATCCAATTTGGTCACGATATTTTGTTTTTAAATTTATATTTGCAGGTTTCACATCGCTGGATTCAATTTTTGGGTATGTTATTTATTACTTTTAAGAAATTAATGTTTGTACATGAATTTACCTAAATGTCTCCTTCTTATGTTGTCTTTTTCTGATCCTGAATATCAAAGTAAttttagcatcattaagtcatttgatggaga	3443	deletion	
	# 4312	chr1	2918978	+	317	chr1	2919453	+	343	Y	N	475	12	N		64	ERVL-B4-int:LTR:ERVL	ERVL-B4-int:LTR:ERVL						0.25		3442	deletion	

	return if $record->{Type} eq 'artifact';

	my %type_map = ('deletion'           		       => 'deletion',
			'tandem-duplication' 		       => 'tandem_duplication',
			'complex'            		       => 'complex_structural_alteration',
			'interchromosomal'   		       => 'translocation',
			'distal-duplication' 		       => 'duplication',
			'probable-inversion' 		       => 'inversion',
			'inversion'          		       => 'inversion',   
			'distal-duplication-by-mobile-element' => 'mobile_element_insertion',
	    );

	#my @positions = sort {$a <=> $b} ($record->{LeftPosition},
	#				  $record->{LeftPosition} + $record->{LeftLength},
	#				  $record->{RightPosition},
	#				  $record->{RightPosition} + $record->{RightLength},
	#);
	
	my $feature_id = 'cgiSvJunction_' . sprintf("%06d", $record->{Id});
	my $seqid      = $record->{LeftChr};
	my $source     = 'cgiJunction';
	my $type       = $type_map{$record->{Type}};
	my $start      = $record->{LeftPosition},
	my $end        = $record->{LeftPosition} + $record->{LeftLength},
	my $score      = '.';
	my $strand     = $record->{LeftStrand};
	my $phase      = '.';
	
	my $variant_seqs  = '~';
	my $reference_seq = '.';
	
	my $attributes = {ID                               => [$feature_id],
			  Reference_seq       		   => [$reference_seq],
			  Variant_seq         		   => [$variant_seqs],
			  cgi_StrandConsistent  	   => [$record->{StrandConsistent}],
			  cgi_Interchromosomal  	   => [$record->{Interchromosomal}],
			  cgi_Distance                     => [$record->{Distance}],
			  cgi_DiscordantMatePairAlignments => [$record->{DiscordantMatePairAlignments}],
			  cgi_JunctionSequenceResolved     => [$record->{JunctionSequenceResolved}],
			  cgi_TransitionSequence           => [$record->{TransitionSequence}],
			  cgi_TransitionLength             => [$record->{TransitionLength}],
			  cgi_LeftRepeatClassification     => [$record->{LeftRepeatClassification}],
			  cgi_RightRepeatClassification    => [$record->{RightRepeatClassification}],
			  cgi_LeftGenes                    => [$record->{LeftGenes}],
			  cgi_RightGenes                   => [$record->{RightGenes}],
			  cgi_XRef                         => [$record->{XRef}],
			  cgi_DeletedTransposableElement   => [$record->{DeletedTransposableElement}],
			  cgi_KnownUnderrepresentedRepeat  => [$record->{KnownUnderrepresentedRepeat}],
			  cgi_FrequencyInBaselineGenomeSet => [$record->{FrequencyInBaselineGenomeSet}],
			  cgi_AssembledSequence            => [$record->{AssembledSequence}],
			  cgi_EventId                      => [$record->{EventId}],
			  cgi_Type                         => [$record->{Type}],
			  cgi_RelatedJunctions             => [$record->{RelatedJunctions}],
	};
	
	map {delete $attributes->{$_} unless length $attributes->{$_}[0]} keys %{$attributes};
	for my $values (values %{$attributes}) {
	    map {s/;/,/g} @{$values};
	}

        my $featureL = {feature_id => $feature_id,
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

	my @features;
	if ($record->{RightChr}) {
	    $featureL->{feature_id} .= '_L';
	    $featureL->{attributes}{ID}[0] = $featureL->{feature_id};

	    my $featureR = dclone($featureL);
	    $featureR->{feature_id} = $feature_id . '_R';
	    $featureR->{seqid}      = $record->{RightChr};
	    $featureR->{start}      = $record->{RightPosition};
	    $featureR->{end}        = $record->{RightPosition} + $record->{RightLength};
	    $featureR->{strand}     = $record->{RightStrand};
	    $featureR->{attributes}{ID}[0] = $featureR->{feature_id};
	    push @features, $featureL, $featureR;
	}
	else {
	    push @features, $featureL;
	}
	return \@features;
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

    my @field_names = qw(Id LeftChr LeftPosition LeftStrand LeftLength
    			 RightChr RightPosition RightStrand RightLength StrandConsistent
    			 Interchromosomal Distance DiscordantMatePairAlignments
    			 JunctionSequenceResolved TransitionSequence TransitionLength
    			 LeftRepeatClassification RightRepeatClassification LeftGenes
    			 RightGenes XRef DeletedTransposableElement
    			 KnownUnderrepresentedRepeat FrequencyInBaselineGenomeSet
    			 AssembledSequence EventId Type RelatedJunctions);

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

L<GAL::Parser::cgi_allJunctionsBeta> does not throw any warnings or errors.

=head1 CONFIGURATION AND ENVIRONMENT

L<GAL::Parser::cgi_allJunctionsBeta> requires no configuration files or
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
    modify it under the same terms as Perl itself.

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
