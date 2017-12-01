package GAL::Parser::cgi_allSvEventsBeta;

use strict;
use vars qw($VERSION);
$VERSION = '0.01';

use base qw(GAL::Parser);
use GAL::Reader::DelimitedLine;
use Storable qw(dclone);

=head1 NAME

GAL::Parser::cgi_allSvEventsBeta - Parse Complete Genomics
cnvSegmentsDiploidBeta files

=head1 VERSION

This document describes GAL::Parser::cgi_allSvEventsBeta version 0.01

=head1 SYNOPSIS

    my $parser = GAL::Parser::cgi_allSvEventsBeta->new(file => 'cgi_allSvEventsBeta.tsv');

    while (my $feature_hash = $parser->next_feature_hash) {
	print $parser->to_gff3($feature_hash) . "\n";
    }

=head1 DESCRIPTION

L<GAL::Parser::cgi_allSvEventsBeta> provides a parser for
cnvSegmentsDiploidBeta.tsv files.

=head1 Constructor

New L<GAL::Parser::cgi_allSvEventsBeta> objects are created by the class method
new.  Arguments should be passed to the constructor as a list (or
reference) of key value pairs.  All attributes of the
L<GAL::Parser::cgi_allSvEventsBeta> object can be set in the call to new. An
simple example of object creation would look like this:

    my $parser = GAL::Parser::cgi_allSvEventsBeta->new(file => 'cgi_allSvEventsBeta.txt');

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
     Usage   : GAL::Parser::cgi_allSvEventsBeta->new();
     Function: Creates a GAL::Parser::cgi_allSvEventsBeta object;
     Returns : A GAL::Parser::cgi_allSvEventsBeta object
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

	# 00 EventId
	# 01 Type
	# 02 RelatedJunctionIds
	# 03 MatePairCounts
	# 04 FrequenciesInBaselineGenomeSet
	# 05 OriginRegionChr
	# 06 OriginRegionBegin
	# 07 OriginRegionEnd
	# 08 OriginRegionLength
	# 09 OriginRegionStrand
	# 10 DestinationRegionChr
	# 11 DestinationRegionBegin
	# 12 DestinationRegionEnd
	# 13 DestinationRegionLength
	# 14 DestinationRegionStrand
	# 15 DisruptedGenes
	# 16 ContainedGenes
	# 17 GeneFusions
	# 18 RelatedMobileElement
	# 19 MobileElementChr
	# 20 MobileElementBegin
	# 21 MobileElementEnd
	# 22 MobileElementStrand
	
	#GENERATED_BY	cgatools
	#GENERATED_AT	2012-Feb-19 00:25:18.457521
	#SOFTWARE_VERSION	2.0.2.20
	#FORMAT_VERSION	2.0
	#TYPE	SV-EVENTS
	
	# >EventId	Type	RelatedJunctionIds	MatePairCounts	FrequenciesInBaselineGenomeSet	OriginRegionChr	OriginRegionBegin	OriginRegionEnd	OriginRegionLength	OriginRegionStrand	DestinationRegionChr	DestinationRegionBegin	DestinationRegionEnd	DestinationRegionLength	DestinationRegionStrand	DisruptedGenes	ContainedGenes	GeneFusions	RelatedMobileElement	MobileElementChr	MobileElementBegin	MobileElementEnd	MobileElementStrand
	# 1	complex	1;110;3809	3;23;12	0.00;0.63;0.63											LOC100287484							
	# 2	complex	2;3;4;5;6;7;8;9;10;11;12;13;14;15;16;17;18;19;21;22;23;24;25;26;27;28;29;30;31;32;33;34;43;49;50;51;52;53;54;55;56;57;58;59;60;61;62;63;64;65;66;67;68;69;70;71;72;73;74;75;76;77;78;79;80;81;82;83;111;112;113;114;115;116;169;170;178;240;308;363;364;365;366;456;457;458;459;460;461;462;463;464;465;466;467;468;469;470;603;604;605;606;607;609;610;611;612;613;614;615;616;617;618;619;620;621;622;623;624;625;626;627;628;629;630;631;632;633;634;635;636;637;638;639;640;641;642;643;644;645;646;647;648;674;675;676;677;678;679;680;681;682;683;684;692;697;728;729;730;731;732;733;734;735;736;802;803;804;805;806;807;808;809;810;811;812;813;872;891;930;1172;1207;1208;1274;1275;1276;1277;1278;1279;1280;1281;1282;1283;1284;1285;1286;1287;1288;1521;1522;1523;1524;1525;1526;1527;1528;1529;1530;1531;1532;1533;1534;1535;1536;1537;1538;1539;1540;1541;1542;1543;1544;1548;1549;1601;1602;1603;1605;1606;1607;1608;1609;1610;1611;1612;1891;1892;1893;1894;1895;2008;2062;2064;2065;2066;2121;2624;3113;3120;3121;3122;3123;3124;3125;3700;3701;3702;3703;3704;3705;3706;3707;3708;3709;3710;3711;3713;3714;3723;3724;3725;3726;3727;3728;3729;3730;3731;3732;3733;3734;3735;3736;3737;3738;3739;3740;3741;3742;3743;3744;3745;3746;3747;3748;3749;3750;3751;3752;3753;3754;3755;3756;3757;3758;3759;3760;3761;3762;3764;3776;3803;3804;3805;3866;3871;3899;3900;4390	3;3;4;5;5;6;4;3;8;3;3;9;4;3;5;6;3;3;4;3;3;3;6;4;5;6;4;3;9;3;3;3;59;6;3;9;3;3;3;6;13;9;10;11;8;6;3;4;4;3;25;3;3;4;9;5;6;3;3;7;9;4;6;14;14;7;5;3;4;3;3;5;4;3;3;4;3;6;16;8;4;10;11;15;3;8;16;3;4;5;3;12;7;4;3;3;3;77;4;4;3;6;7;3;3;12;12;30;3;8;13;7;3;3;11;4;8;15;7;3;11;21;4;3;3;5;3;4;3;3;7;6;4;4;13;4;3;5;3;3;25;4;8;10;4;3;32;3;9;3;5;5;7;4;3;18;3;5;3;3;27;12;4;3;11;3;38;20;14;11;5;4;4;4;3;3;3;8;4;16;3;5;3;29;5;13;3;4;4;3;11;3;4;4;3;3;3;4;3;8;4;4;3;7;6;3;4;5;3;5;44;3;6;3;3;3;4;3;4;3;4;3;3;3;7;13;3;3;4;4;8;4;4;8;3;8;3;3;3;3;31;4;3;4;3;3;7;3;3;4;5;4;3;4;3;3;7;4;8;3;3;3;3;3;3;7;4;5;42;8;34;17;10;8;5;6;4;4;3;3;7;6;3;8;5;5;7;3;3;3;17;4;3;3;3;5;4;6;6;4;4;3;4;3;4;3;8;3;3;10;3;3;3;3;3;4;9;3	1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;0.96;0.96;0.96;1.00;1.00;1.00;1.00;1.00;1.00;1.00;0.75;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;0.11;0.22;0.11;0.22;0.59;0.59;0.98;0.98;0.33;0.77;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;0.52;0.52;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;0.89;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;0.04;0.65;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;0.70;0.07;0.98;0.06;0.30;0.30;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;0.26;0.26;1.00;1.00;1.00;0.26;0.26;1.00;1.00;1.00;1.00;1.00;1.00;0.98;0.98;0.98;0.98;0.98;1.00;0.19;0.96;0.96;0.96;0.48;0.92;0.41;0.74;0.74;0.74;0.74;0.37;0.37;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;0.00;0.00;0.00;0.00;0.11;0.15;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;0.19;0.63;0.93;0.93;0.93;0.65;0.07;0.96;0.96;0.15											DUX4C;IL1RAPL1;LOC100287759;LOC100288801;LOC727758							
	# 3	probable-inversion	20	4	0.00	chrY	28814196	58840770	30026574	-													
	# 4	interchromosomal	35	6	0.67	chrX	3195494	3195494	0	-	chrY	14172830	14172830	0	+								
	# 5	complex	36;1547;3712	5;23;14	0.00;0.00;0.00																		
	# 6	complex	37;86;87;158;824;1234;1256;1628;2007;3130;3137;3163;3765;3849;4322	18;5;36;160;21;6;6;89;12;11;63;24;9;45;7	1.00;1.00;1.00;1.00;0.15;0.10;0.63;1.00;1.00;0.46;1.00;1.00;0.52;0.98;0.88											PARN;STAG1		PARN/STAG1					
	# 7	distal-duplication	38;3715	7;12	0.54;0.52	chrX	131394197	131394278	81	-	chrX	131393575	131393591	16	+								
	# 8	artifact	39	3	0.00																		
	# 9	complex	40;100;101;108;149;171;173;179;267;410;411;414;417;424;425;437;438;439;440;441;444;452;659;685;686;687;688;689;690;691;693;694;718;719;721;724;738;825;952;955;956;983;1237;1289;1328;1329;1337;1338;1377;1545;1613;1614;1615;1616;1617;1623;1736;1739;1746;1918;1919;1920;1921;1922;1923;1924;1925;1926;1927;1928;1929;1930;1931;1932;1933;1934;1935;1936;1937;2216;2326;2329;2334;3132;3143;3208;3655;3717;3763;3801;3802;3843;3949;3950;4102;4211;4213;4214;4215;4216;4217;4218;4319;4320	4;1330;575;10;8;6;8;9;22;5;162;38;15;49;30;4;3;13;6;8;3;58;78;49;5;5;4;7;5;3;33;3;4;47;4;34;19;5;808;4;13;38;28;3;3;7;125;269;6;100;45;98;29;26;5;88;72;5;29;86;19;27;14;126;50;11;8;75;17;12;4;99;18;29;50;6;226;369;3;13;86;6;21;16;7;5;642;41;15;43;5;16;100;156;396;87;10;9;9;189;167;288;730;302	0.07;0.98;0.00;0.87;0.94;0.96;0.98;1.00;1.00;0.00;1.00;0.00;1.00;1.00;0.90;0.73;0.73;0.98;0.63;0.71;0.00;0.19;0.15;1.00;0.67;0.74;0.96;0.96;1.00;1.00;1.00;0.26;0.15;1.00;0.22;1.00;1.00;1.00;1.00;0.54;0.94;0.38;0.37;0.81;0.04;0.83;1.00;0.96;0.71;1.00;1.00;1.00;1.00;1.00;0.22;1.00;0.77;0.00;0.96;1.00;1.00;1.00;0.94;1.00;1.00;1.00;1.00;1.00;0.67;0.96;0.00;1.00;1.00;1.00;1.00;1.00;1.00;1.00;0.75;0.48;1.00;0.81;0.98;0.98;0.31;0.87;1.00;0.00;0.07;0.00;0.96;0.79;1.00;1.00;1.00;1.00;0.79;0.00;0.85;1.00;1.00;1.00;0.98;0.98											CDC27;FAM156B;LOC100288870;LOC284801;LOC401629;NCRNA00164;PMF1;SLIT1;TOP		NCRNA00164/LOC284801;TOP/CDC27;TSS-UPSTREAM[LOC284801]/NCRNA00164					
	# 10	artifact	41	3	0.00																		
	# 11	complex	42;1564;3719	8;15;11	0.48;0.48;0.48																		
	# 12	probable-inversion	44	7	0.62	chrX	896835	948129	51294	-													
	# 13	interchromosomal	45	37	0.94	chr20	11281583	11281583	0	-	chrX	81651234	81651234	0	+								

	return if $record->{Type} =~ /^(artifact|complex)$/;

	my %type_map = ('deletion'           		       => 'deletion',			    
			'tandem-duplication' 		       => 'tandem_duplication',		    
			'complex'            		       => 'complex_structural_alteration',  
			'interchromosomal'   		       => 'interchromosomal_breakpoint',    
			'distal-duplication' 		       => 'duplication',		    
			'probable-inversion' 		       => 'inversion',			    
			'inversion'          		       => 'inversion',			    
			'distal-duplication-by-mobile-element' => 'mobile_element_insertion',       
	    );

	my $feature_id = 'cgiSvEvent_' . sprintf("%06d", $record->{EventId});
	my $source     = 'cgiSvEvent';
	my $type       = $type_map{$record->{Type}};
	my $score      = '.';
	my $phase      = '.';
	
	my $variant_seqs  = '~';
	my $reference_seq = '.';
	
	my @positions;
	my ($orgChr, $orgBegin, $orgEnd, $orgStrand, $destChr, $destBegin, $destEnd, $destStrand) =
	    @{$record}{qw(OriginRegionChr
			  OriginRegionBegin         
			  OriginRegionEnd
			  OriginRegionStrand
			  DestinationRegionChr
			  DestinationRegionBegin
			  DestinationRegionEnd
			  DestinationRegionStrand
			  )};

	my $featureOrg  = {};
	my $featureHold = {};
	if ($orgChr && $destChr) {
	    # Double locus
	    if ($orgChr eq $destChr) {
		# Intrachromosomal
		if ($orgBegin < $destEnd && $orgEnd > $destBegin) {
		    # Loci overlap
		    my @positions =  sort {$a <=> $b} ($orgBegin, $orgEnd, $destBegin, $destEnd);
		    $featureOrg->{seqid} = $orgChr;
		    $featureOrg->{start} = shift @positions;
		    $featureOrg->{end}   = pop @positions;
		    $featureOrg->{strand} = $orgStrand;
		}
		else {
		    # Double locus - loci don't overlap
		    @{$featureOrg}{qw(seqid start end strand)}  = ($orgChr, $orgBegin, $orgEnd, $orgStrand);
		    @{$featureHold}{qw(seqid start end strand)} = ($destChr, $destBegin, $destEnd, $orgStrand);
		}
	    }
	    else {
		# Interchromosomal
		@{$featureOrg}{qw(seqid start end strand)}  = ($orgChr, $orgBegin, $orgEnd, $orgStrand);
		@{$featureHold}{qw(seqid start end strand)} = ($destChr, $destBegin, $destEnd, $orgStrand);
	    }
	}
	else {
	    # Single locus
	    @{$featureOrg}{qw(seqid start end strand)}  = ($orgChr, $orgBegin, $orgEnd, $orgStrand);
	}
	
	my $attributes = {ID                                 => [$feature_id],
			  Reference_seq       		     => [$reference_seq],
			  Variant_seq         		     => [$variant_seqs],
			  cgi_Type    		 	     => [$record->{Type}],
			  cgi_MatePairCounts     	     => [$record->{MatePairCounts}],
			  cgi_FrequenciesInBaselineGenomeSet => [$record->{FrequenciesInBaselineGenomeSet}],
			  cgi_DisruptedGenes 		     => [$record->{DisruptedGenes}],
			  cgi_ContainedGenes 		     => [$record->{ContainedGenes}],
			  cgi_GeneFusions                    => [$record->{GeneFusions}],
			  cgi_RelatedMobileElement 	     => [$record->{RelatedMobileElement}],
	};
	
	map {delete $attributes->{$_} unless length $attributes->{$_}[0]} keys %{$attributes};
        for my $values (values %{$attributes}) {
            map {s/;/,/g} @{$values};
	    }
	
        @{$featureOrg}{qw(feature_id    source   type   score   phase   attributes)} = 
	                 ($feature_id, $source, $type, $score, $phase, $attributes);

	my @features;

	if ($destChr) {
	    my $featureDest = dclone($featureOrg);
	    $featureOrg->{feature_id} = $feature_id . '_O';
	    $featureOrg->{attributes}{ID}[0] = $featureOrg->{feature_id};
	    $featureDest->{feature_id} = $feature_id . '_D';
	    $featureDest->{attributes}{ID}[0] = $featureDest->{feature_id};

	    @{$featureDest}{qw(seqid start end strand)} = ($destChr, $destBegin, $destEnd, $destStrand);
	    push @features, $featureDest;
	}
	push @features, $featureOrg;

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

    my @field_names = qw(EventId Type RelatedJunctionIds
    			 MatePairCounts FrequenciesInBaselineGenomeSet
    			 OriginRegionChr OriginRegionBegin
    			 OriginRegionEnd OriginRegionLength
    			 OriginRegionStrand DestinationRegionChr
    			 DestinationRegionBegin DestinationRegionEnd
    			 DestinationRegionLength
    			 DestinationRegionStrand DisruptedGenes
    			 ContainedGenes GeneFusions
    			 RelatedMobileElement MobileElementChr
    			 MobileElementBegin MobileElementEnd
    			 MobileElementStrand);


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

L<GAL::Parser::cgi_allSvEventsBeta> does not throw any warnings or errors.

=head1 CONFIGURATION AND ENVIRONMENT

L<GAL::Parser::cgi_allSvEventsBeta> requires no configuration files or
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
