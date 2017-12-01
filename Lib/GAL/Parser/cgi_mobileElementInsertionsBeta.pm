package GAL::Parser::cgi_mobileElementInsertionsBeta;

use strict;
use vars qw($VERSION);
$VERSION = '0.01';

use base qw(GAL::Parser);
use GAL::Reader::DelimitedLine;

=head1 NAME

GAL::Parser::cgi_mobileElementInsertionsBeta - Parse Complete Genomics
cnvSegmentsDiploidBeta files

=head1 VERSION

This document describes GAL::Parser::cgi_mobileElementInsertionsBeta version 0.01

=head1 SYNOPSIS

    my $parser = GAL::Parser::cgi_mobileElementInsertionsBeta->new(file => 'cgi_mobileElementInsertionsBeta.tsv');

    while (my $feature_hash = $parser->next_feature_hash) {
	print $parser->to_gff3($feature_hash) . "\n";
    }

=head1 DESCRIPTION

L<GAL::Parser::cgi_mobileElementInsertionsBeta> provides a parser for
cnvSegmentsDiploidBeta.tsv files.

=head1 Constructor

New L<GAL::Parser::cgi_mobileElementInsertionsBeta> objects are created by the class method
new.  Arguments should be passed to the constructor as a list (or
reference) of key value pairs.  All attributes of the
L<GAL::Parser::cgi_mobileElementInsertionsBeta> object can be set in the call to new. An
simple example of object creation would look like this:

    my $parser = GAL::Parser::cgi_mobileElementInsertionsBeta->new(file => 'cgi_mobileElementInsertionsBeta.txt');

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
     Usage   : GAL::Parser::cgi_mobileElementInsertionsBeta->new();
     Function: Creates a GAL::Parser::cgi_mobileElementInsertionsBeta object;
     Returns : A GAL::Parser::cgi_mobileElementInsertionsBeta object
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
	
	#ASSEMBLY_ID	GS000007756-ASM
	#SOFTWARE_VERSION	2.0.2.20
	#GENERATED_BY	ExportWorkflow&ReportEngine
	#GENERATED_AT	2012-02-19 00:24:27.408501
	#FORMAT_VERSION	2.0
	#GENOME_REFERENCE	NCBI build 37
	#SAMPLE	GS01003-DNA_A01
	#TYPE	MEI
	#GENE_ANNOTATIONS	NCBI build 37.2
	#MEI_1000G_ANNOTATIONS	INITIAL-DATA-RELEASE
	
	# >Chromosome	InsertRangeBegin	InsertRangeEnd	Strand	ElementType	ElementTypeScore	ElementSequenceBegin	ElementSequenceEnd	NextBestElementType	InsertionScore	InsertionDnbCount	InsertionLeftDnbCount	InsertionRightDnbCount	ReferenceDnbCount	GeneOverlap	XRef	FrequencyInBaseline	NovelEventCountForInsertionScore	KnownEventSensitivityForInsertionScore
	# chr1	825329	825775	-	AluY	0	8	279	AluYa1	354	20	19	1	N			1.000	841	0.990
	# chr1	1584407	1584563	-	AluYa4	0	0	221	AluYa5	563	21	17	4	3	CDK11B:-:INTRON		0.846	669	0.977
	# chr1	1647620	1647754	-	AluYa1	24	10	292	AluYa5	730	30	1	29	173	CDK11A:-:INTRON		1.000	591	0.970
	# chr1	1650934	1650964	-	AluSc	19	94	192	AluYa1	78	4	4	0	148	CDK11A:-:INTRON		0.981	1499	0.997
	# chr1	2295693	2295723	-	L1HS	122	5063	5217	L1PA13	311	21	21	0	176	MORN1:-:INTRON		0.173	897	0.992
	# chr1	3995420	3995492	-	AluY	0	148	259	AluYd3	47	4	0	4	130			1.000	1912	0.998
	# chr1	5726192	5726578	-	AluY	0	1	262	AluYa1	419	23	16	7	N			1.000	772	0.984
	# chr1	5729479	5729972	+	AluSp	70	54	201	AluSq	145	8	3	5	N			0.788	1196	0.997

	my $seqid      = $record->{Chromosome};
	my $source     = 'cgiMEI';
	my $type       = 'insertion';
	my $start      = $record->{InsertRangeBegin};
	my $end        = $record->{InsertRangeEnd};
	my $score      = '.';
	my $strand     = $record->{Strand};
	my $phase      = '.';
	my $feature_id = 'cgiMEI_' . sprintf("%06d", $self->counter('MEI'));

	my $variant_seqs = '~';
	my $reference_seq = '-';

	my $attributes = {Reference_seq       => [$reference_seq],
			  Variant_seq         => [$variant_seqs],
			  ID                  => [$feature_id],
			  cgi_ElementType            		     => [$record->{ElementType}],
			  cgi_ElementTypeScore       		     => [$record->{ElementTypeScore}],
			  cgi_ElementSequenceBegin   		     => [$record->{ElementSequenceBegin}],
			  cgi_ElementSequenceEnd     		     => [$record->{ElementSequenceEnd}],
			  cgi_NextBestElementType    		     => [$record->{NextBestElementType}],
			  cgi_InsertionScore         		     => [$record->{InsertionScore}],
			  cgi_InsertionDnbCount      		     => [$record->{InsertionDnbCount}],
			  cgi_InsertionLeftDnbCount  		     => [$record->{InsertionLeftDnbCount}],
			  cgi_InsertionRightDnbCount 		     => [$record->{InsertionRightDnbCount}],
			  cgi_ReferenceDnbCount      		     => [$record->{ReferenceDnbCount}],
			  cgi_GeneOverlap                            => [$record->{GeneOverlap}],
			  cgi_XRef                                   => [$record->{XRef}],
			  cgi_FrequencyInBaseline                    => [$record->{FrequencyInBaseline}],
			  cgi_NovelEventCountForInsertionScore       => [$record->{NovelEventCountForInsertionScore}],
			  cgi_KnownEventSensitivityForInsertionScore => [$record->{KnownEventSensitivityForInsertionScore}],
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

    my @field_names = qw(Chromosome InsertRangeBegin InsertRangeEnd
    			 Strand ElementType ElementTypeScore ElementSequenceBegin
    			 ElementSequenceEnd NextBestElementType InsertionScore
    			 InsertionDnbCount InsertionLeftDnbCount InsertionRightDnbCount
    			 ReferenceDnbCount GeneOverlap XRef FrequencyInBaseline
    			 NovelEventCountForInsertionScore
    			 KnownEventSensitivityForInsertionScore);


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

L<GAL::Parser::cgi_mobileElementInsertionsBeta> does not throw any warnings or errors.

=head1 CONFIGURATION AND ENVIRONMENT

L<GAL::Parser::cgi_mobileElementInsertionsBeta> requires no configuration files or
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
