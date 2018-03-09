#!/usr/bin/perl
#######################################################
# Author  :  Aurelie Kapusta (https://github.com/4ureliek), with the help of Edward Chuong
# email   :  4urelie.k@gmail.com  
# Purpose :  Writen to test enrichment of TEs in a set of simple features (ChIP-seq for example) 
#######################################################
use strict;
use warnings;
use Carp;
use Getopt::Long;
use Bio::SeqIO;
use Statistics::R; #required to get the Binomial Test p-values
use Data::Dumper;
use vars qw($BIN);
use Cwd 'abs_path';
BEGIN { 	
	$BIN = abs_path($0);
	$BIN =~ s/(.*)\/.*$/$1/;
	unshift(@INC, "$BIN/Lib");
}
use TEshuffle;

#-------------------------------------------------------------------------------
#------------------------------- DESCRIPTION -----------------------------------
#-------------------------------------------------------------------------------
#flush buffer
$| = 1;

my $VERSION = "4.4";
my $SCRIPTNAME = "TE-analysis_Shuffle_bed.pl";
my $CHANGELOG;
set_chlog();
sub set_chlog {
$CHANGELOG = "
#  - v1.0 = Mar 2016 
#           based on TE-analysis_Shuffle_v3+.pl, v3.3, 
#           but adapted to more general input files = bed file corresponding to any features to test.
#  - v2.0 = Mar 2016 
#           attempt of making this faster by removing any length info and allowing overlaps 
#           of shuffled features (since they are independent tests, it's OK)
#           Also added the possibility of several files to -e and -i
#  - v2.1 = Oct 2016
#           remove empty column of length info from output
#           get enrichment by age categories if age file provided
#           bug fix for total counts of hit features when upper levels (by class or family, by Rname was probably OK)
#           Changes in stats, bug fix; use R for the binomial test
#  - v3.0 = Oct 25 2016
#           TEshuffle.pm for subroutines shared with the shuffle_tr script
#  - v4.0 = Nov 28 2016
#           different choices to shuffle the TEs:
#             shufflebed = completely random positions, but same chromosome
#             shuffle inside the current TE positions, same chromosome
#             shuffle each TE, keeping its distance to a TSS, same chromosome - thanks to: Cedric Feschotte, Ed Chuong
#           make subfolders for each input file (for the shuffled outputs)
#  - v4.1 = Nov 30 2016
#           Minor bug fix: when random was 0, values were not reported 
#           (still interesting to see them in the obs, even if no stats possible)
#  - v4.2 = Dec 01 2017
#           Bug fix for when long TEs were shuffled to the position of small TEs that are
#             too close to the start of the genomic sequence (led to negative starts).
#             This is now checked for -s rm and -s tss:
#               for -s rm, the TE is shifted of as many bp as needed
#               for -s tss the start is simply changed to 1 (to avoid having a TE placed closer to a tss)
#           Also added the option to use -r file in -s rm as well (to check ends) and chift the TE if needed.
#  - v4.3 = Jan 23 2018
#           Avoid passing to subs variables that don't need to + use upercase syntax
#           Cleaned up a few things
#           Bug fix for cases when the 4th col of the bed file is not a unique ID => identifier is now a cat of cols 1 to 4
#  - v4.4 = Mar 08 2018
#           Minor changes in logs and check options
#           Was assuming a 5 columns bed file => fix
#           Change for -s tss: if the TE ends up out of the scaffold/chr, put on the other side and if still out, shift it to be inside
#           Also, skip if not in the annotation file
#           Make -r mandatory for all
\n";
	return 1;
}

my $USAGE;
set_usage();
sub set_usage {
	$USAGE = "
Synopsis (v$VERSION):

    perl $SCRIPTNAME -f features.bed [-o <nt>] -q features_to_shuffle [-n <nb>] -s shuffling_type -r <genome.sizes> [-b]
            [-a <annotations>] 
            [-e <genome.gaps>] [-d] [-i <include.range>] [-x] 
            [-w <bedtools_path>] [-l <if_nonTE>] [-t <filterTE>] [-c] [-g <TE.age.tab>] [-v] [-h]

    /!\\ REQUIRES: Bedtools, at least v18 (but I advise updating up to the last version; v26 has a better shuffling)
    /!\\ Previous outputs, if any, will be moved as *.previous [which means previous results are only saved once]

    Typically, for the 3 types, the mandatory arguments are:
    perl $SCRIPTNAME -f features.bed -q rm.out -r genome.sizes -s rm 
    perl $SCRIPTNAME -f features.bed -q rm.out -r genome.sizes -s tss -a annotations.gtf
    perl $SCRIPTNAME -f features.bed -q rm.out -r genome.sizes -s bed -e genome.gaps
	
   CITATION:
    - include the version of the script + link to the GitHub page (https://github.com/4ureliek/TEanalysis)
    - Cite Kapusta et al. (2013) PLoS Genetics (DOI: 10.1371/journal.pgen.1003470)
        (http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1003470)
      Also cite Trizzani, Kapusta and Brown (2018) BioRxiv if possible 
        (https://www.biorxiv.org/content/early/2018/02/20/268771 - https://doi.org/10.1101/268771)
    - for BEDtools, please cite 
      Quinlan AR and Hall IM (2010) Bioinformatics (DOI: 10.1093/bioinformatics/btq033)

   DESCRIPTION:
    Features provided in -s will be overlapped with -f file (which must be simple intervals in bed format), 
       without (no_boot) or with (boot) shuffling (on same chromosome)
       One feature may overlap with several repeats and all are considered.
       However, if there are several fragments of the same repeat in a feature, it will be counted
       only one time. Could be an issue for larger features, but otherwise can't normalize by total
       number of features.
       There are 3 options for shuffling, but all will shuffle on the same chromosome:
          - with bedtool shuffle (=> random position)
            [will require files for -r and -e that can be generated by this script, see below]
          - shuffle among the current TE positions (still random, but less)
          - shuffle the TEs keeping their distance to the closest TSS of an annotation file provided
            [will require ensembl gtf, or gencode gtf or gff3 annotations]

       Note that because TEs are often fragmented + there are inversions, the counts of TEs are likely inflated;
       this also means that when TEs are shuffled, there are more fragments than TEs. Some should be moved non independently, 
       or the input file should be corrected when possible to limit that issue 
       [not implemented in this script for now, but you may edit the RM file to merge TE fragments]

    Two-tailed permutation test and a binomial test are done on the counts of overlaps. 
       The results are in a .stats.txt file. Note that high bootstraps takes a lot of time. 
       Note that for low counts, expected and/or observed, stats likely don't mean much.


   MANDATORY ARGUMENTS:	
    -f,--feat     => (STRING) ChIPseq peaks, chromatin marks, etc, in bed format
                              /!\\ Script assumes no overlap between features
    -q,--query    => (STRING) Features to shuffle = TE file
                              For now, can only be the repeat masker .out or the .bed file generated by the TE-analysis_pipeline script
                              See -l and -t for filters, and -g for age data       
    -s,--shuffle  => (STRING) Shuffling type. Should be one of the following:
                               -s bed  => use bedtools shuffle, random position on the same chromosome
                               -s rm   => shuffle inside the current TE positions; still random, but less; also, with this
                                          the stats won't be affected with overall depletion or enrichment of TEs in the -f file
                               -s tss  => shuffle the TEs on the same chromosome keeping their distance
                                          to the closest TSS of an annotation file provided
                                          (TSS are shuffled and then assigned to a TE; 
                                          same TSS will be assigned multiple TEs if fewer TSS than TEs).
                                          Thanks to: Edward Chuong & Cedric Feschotte for the suggestion
    -r,--range    => (STRING) To know the maximum value in a given chromosome/scaffold. 
                              File should be: sequence_name \\t length
                              Can be files from UCSC, files *.chrom.sizes
                              If you don't have such file, use -b (--build) and provide the genome fasta file for -r  

   MANDATORY ARGUMENTS IF USING -s tss:
    -a,--annot    => (STRING) gtf or gff; annotations to load all unique TSS: will be used to set the distance
                              between each TE and the closest TSS that will be kept while randomization
                              Note that it requires transcript lines
                              
   MANDATORY ARGUMENTS IF USING -s bed:                           
    -e,--excl     => (STRING) This will be used as -excl for bedtools shuffle: \"coordinates in which features from -i should not be placed.\"
                              More than one file may be provided (comma separated), they will be concatenated 
                              (in a file = first-file-name.cat.bed).
                              By default, at least one file is required = assembly gaps, and it needs to be the first file
                              if not in bed format. Indeed, you may provide the UCSC gap file, with columns as:
                                  bin, chrom, chromStart, chromEnd, ix, n, size, type, bridge
                              it will be converted to a bed file. 
                              If you do nothave this file, you may provide the genome file in fasta format
                              and add the option -d (--dogaps), to generate a bed file corresponding to assembly gaps.
                              If you need to generate the <genome.gaps> file but you would also like to add more files to the -e option, 
                              just do a first run with no bootstraps (in this example the genome.range is also being generated):
                                 perl ~/bin/$SCRIPTNAME -f input.bed -q genome.out -s rm -r genome.fa -b -e genome.fa -d -n 0                                
                              Other files may correspond to regions of low mappability, for example for hg19:
                              http://www.broadinstitute.org/~anshul/projects/encode/rawdata/blacklists/hg19-blacklist-README.pdf
                              Notes: -> when the bed file is generated by this script, any N stretch > 50nt will be considered as a gap 
                                        (this can be changed in the load_gap subroutine)         
                                     -> 3% of the shuffled feature may overlap with these regions 
                                        (this can be changed in the shuffle subroutine).

   OPTIONAL ARGUMENTS IF USING -s bed:
    -d,--dogaps   => (BOOL)   See above; use this and provide the genome fasta file if no gap file (-g)
                              If several files in -e, then the genome needs to be the first one.
                              This step is not optimized, it will take a while (but will create the required file)                       
    -i,--incl     => (STRING) To use as -incl for bedtools shuffle: \"coordinates in which features from -i should be placed.\"
                              Bed of gff format. Could be intervals close to TSS for example.
                              More than one file (same format) may be provided (comma separated), 
                              they will be concatenated (in a file = first-file-name.cat.bed)
    -x,--x        => (BOOL)   to add the -noOverlapping option to the bedtools shuffle command line, 
                              and therefore NOT allow overlaps between the shuffled features.
                              This may create issues mostly if -i is used (space to shuffle may be too small to shuffle features)
                             
   OTHER OPTIONAL ARGUMENTS (for all three -s):                                  
    -b,--build    => (BOOL)   See above; use this and provide the genome fasta file if no range/lengths file (-r)
                              This step may take a while but will create the required file	
    -o,--overlap  => (INT)    Minimal length (in nt) of intersection in order to consider the TE included in the feature.
                              Default = 10 (to match the TEanalysis-pipeline.pl)
    -n,--nboot    => (STRING) number of bootsraps with shuffled -s file
                              Default = 100 for faster runs; use higher -n for good pvalues 
                              (-n 10000 is best for permutation test but this will take a while)
                              If set to 0, no bootstrap will be done
    -w,--where    => (STRING) if BEDtools are not in your path, provide path to BEDtools bin directory                             

   OPTIONAL ARGUMENTS FOR TE FILTERING (for all -s): 
    -l,--low      => (STRING) To set the behavior regarding non TE sequences: all, no_low, no_nonTE, none
                                 -l all = keep all non TE sequences (no filtering)
                                 -l no_low [default] = keep all besides low_complexity and simple_repeat
                                 -l no_nonTE = keep all except when class = nonTE
                                 -l none = everything is filtered out (nonTE, low_complexity, simple_repeat, snRNA, srpRNA, rRNA, tRNA/tRNA, satellite)
    -t,--te       => (STRING) <type,name>
                              run the script on only a subset of repeats. Not case sensitive.
                              The type can be: name, class or family and it will be EXACT MATCH unless -c is chosen as well
                              ex: -t name,nhAT1_ML => only fragments corresponding to the repeat named exactly nhAT1_ML will be looked at
                                  -t class,DNA => all repeats with class named exactly DNA (as in ...#DNA/hAT or ...#DNA/Tc1)
                                  -t family,hAT => all repeats with family named exactly hAT (so NOT ...#DNA/hAT-Charlie for example)
    -c,--contain  => (BOOL)   to check if the \"name\" determined with -filter is included in the value in Repeat Masker output, instead of exact match
                              ex: -t name,HERVK -c => all fragments containing HERVK in their name
                                  -t family,hAT -c => all repeats with family containing hAT (...#DNA/hAT, ...#DNA/hAT-Charlie, etc)
    -g,--group    => (STRING) provide a file with TE age: 
                                 Rname  Rclass  Rfam  Rclass/Rfam  %div(avg)  lineage  age_category
                              At least Rname and lineage are required (other columns can be \"na\"),
                              and age_category can be empty. But if age_category has values, it will 
                              be used as well. Typically:
                                  TE1  LTR  ERVL-MaLR  LTR/ERVL-MaLR  24.6  Eutheria  Ancient
                                  TE2  LTR  ERVL-MaLR  LTR/ERVL-MaLR   9.9  Primates  LineageSpe
  
   OPTIONAL ARGUMENTS (GENERAL): 
    -v,--version  => (BOOL)   print the version
    -h,--help     => (BOOL)   print this usage
\n";

	return 1;
}

#-------------------------------------------------------------------------------
#------------------------------ LOAD AND CHECK ---------------------------------
#-------------------------------------------------------------------------------
my ($IN,$SHUFFLE,$STYPE,$TSSFILE,$EXCLUDE,$DOGAPS,$BUILD,$DOBUILD,$F_REGEX,$ALLOW,$NOOVERLAPS,$V,$HELP);
my $INTERS = 10;
my $NBOOT = 10;
my $INCL = "na";
my $NONTE = "no_low";
my $FILTER = "na";
my $TEAGE = "na";
my $BEDTOOLS = "";
GetOptions(
	  'feat=s'		=> \$IN,
	  'query=s'     => \$SHUFFLE,
	  'shuffle=s'   => \$STYPE,
	  'overlap=s'   => \$INTERS,
	  'nboot=s'     => \$NBOOT,
	  'annot=s'     => \$TSSFILE,			 	  
	  'range=s'     => \$BUILD,
	  'build'       => \$DOBUILD,
	  'excl=s'		=> \$EXCLUDE,
	  'dogaps'      => \$DOGAPS,
	  'incl=s'		=> \$INCL,
	  'x'		    => \$NOOVERLAPS,
	  'low=s'		=> \$NONTE,
	  'te=s'		=> \$FILTER,
	  'contain'     => \$F_REGEX,
	  'group=s'     => \$TEAGE,
	  'where=s'     => \$BEDTOOLS,
	  'version'     => \$V,
	  'help'		=> \$HELP,);

#Check options, if files exist, etc
check_options();
sub check_options {
	die "\n --- $SCRIPTNAME version $VERSION\n\n" if ($V);
	die $USAGE if ($HELP);
	die "\n SOME MANDATORY ARGUMENTS MISSING, CHECK USAGE:\n$USAGE" if (! $IN || ! $SHUFFLE || ! $STYPE  || ! $BUILD);
	die "\n -f $IN is not a bed file?\n\n" unless ($IN =~ /\.bed$/);
	die "\n -f $IN does not exist?\n\n" if (! -e $IN);
	die "\n -q $SHUFFLE is not in a proper format? (not .out, .bed, .gff or .gff3)\n\n" unless ($SHUFFLE =~ /\.out$/ || $SHUFFLE =~ /\.bed$/ || $SHUFFLE =~ /\.gff$/ || $SHUFFLE =~ /\.gff3$/);
	die "\n -q $SHUFFLE does not exist?\n\n" if (! -e $SHUFFLE);
	die "\n -s $STYPE should be one of the following: bed, rm or tss\n\n" if ($STYPE ne "bed" && $STYPE ne "rm" && $STYPE ne "tss");	
	die "\n -r $BUILD does not exist?\n\n" if (! -e $BUILD);

	#deal with conditional mandatory stuff
	if ($STYPE eq "tss") {
		die "\n -s tss was set, but -a is missing?\n\n" if (! $TSSFILE);
		die "\n -a $TSSFILE does not exist?\n\n" if ($TSSFILE && ! -e $TSSFILE);	
	} elsif ($STYPE eq "bed") {
		die "\n -s bed was set, but -e is missing?\n\n" if (! $EXCLUDE);
		die "\n -e $EXCLUDE does not exist?\n\n" if ($EXCLUDE !~ /,/ && ! -e $EXCLUDE); #if several files, can't check existence here
		die "\n -i $INCL does not exist?\n\n" if ($INCL ne "na" && $INCL !~ /,/ && ! -e $INCL); #if several files, can't check existence here
	}

	#Now the rest
	die "\n -n $NBOOT but should be an integer\n\n" if ($NBOOT !~ /\d+/);
	die "\n -i $INTERS but should be an integer\n\n" if ($INTERS !~ /\d+/);
	die "\n -w $BEDTOOLS does not exist?\n\n" if ($BEDTOOLS ne "" && ! -e $BEDTOOLS);
	die "\n -t requires 2 values separated by a coma (-t <name,filter>; use -h to see the usage)\n\n" if ($FILTER ne "na" && $FILTER !~ /,/);
	die "\n -g $TEAGE does not exist?\n\n" if ($TEAGE ne "na" && ! -e $TEAGE);
	($DOGAPS)?($DOGAPS = "y"):($DOGAPS = "n");
	($DOBUILD)?($DOBUILD = "y"):($DOBUILD = "n");
	($F_REGEX)?($F_REGEX = "y"):($F_REGEX="n");
	$BEDTOOLS = $BEDTOOLS."/" if ($BEDTOOLS ne "" && substr($BEDTOOLS,-1,1) ne "/"); #put the / at the end of path if not there
	($NOOVERLAPS)?($NOOVERLAPS = "-noOverlapping"):($NOOVERLAPS = "");
}

#-------------------------------------------------------------------------------
#----------------------------------- MAIN --------------------------------------
#-------------------------------------------------------------------------------
#Quick log on files used
print STDERR "\n --- $SCRIPTNAME v$VERSION started, with:\n";
print STDERR "     input file = $IN\n";
print STDERR "     features to shuffle = $SHUFFLE\n";
print STDERR "     shuffling type = $STYPE\n";
print STDERR "     genome chromosome sizes = $BUILD\n";
my $BEDV = `$BEDTOOLS bedtools --version`;
chomp $BEDV;
print STDERR "     bedtools version = $BEDV\n";

#Prep outputs
print STDERR " --- prepping output directories and files\n";
my $DIR = $IN.".shuffle-".$STYPE.".".$NBOOT;
print STDERR "     output directory = $DIR\n";
my ($STATS,$OUT,$OUTB,$TEMP) = TEshuffle::prep_out("bed",$DIR,$NBOOT,$FILTER,$IN,$STYPE,$NONTE);

#Chosomosome sizes / Genome range
my ($OKSEQ,$BUILD_FILE);
print STDERR " --- loading build (genome range)\n";
($OKSEQ,$BUILD_FILE) = TEshuffle::load_build($BUILD,$DOBUILD);

#prep steps if shuffling type is bed
my $EXCL;
if ($STYPE eq "bed") {
	#Files to exclude for shuffling
	print STDERR " --- getting ranges to exclude in the shuffling of features from $EXCLUDE\n";
	my @exclude = ();
	if ($EXCLUDE =~ /,/) {
		($DOGAPS eq "y")?(print STDERR "     several files provided, -d chosen, genome file (fasta) should be the first one\n"):
						 (print STDERR "     several files provided, assembly gaps should be the first one\n");
		@exclude = split(",",$EXCLUDE) if ($EXCLUDE =~ /,/);
	} else {
		$exclude[0] = $EXCLUDE;
	}
	$exclude[0] = TEshuffle::load_gap($exclude[0],$DOGAPS);
	print STDERR "     concatenating files for -e\n" if ($EXCLUDE =~ /,/);
	 if ($EXCLUDE =~ /,/) {
	 	$EXCL = TEshuffle::concat_beds(\@exclude);
	 } else {
	 	$EXCL = $exclude[0];
	 }	

	#If relevant, files to include for shuffling
	if ($INCL ne "na" && $INCL =~ /,/) {
		print STDERR " --- concatenating $INCL files to one file\n";
		my @include = split(",",$INCL);
		$INCL = TEshuffle::concat_beds(\@include);
	}
}

#Load TEage if any
my $AGE = ();
if ($TEAGE ne "na") {
	print STDERR " --- Loading TE ages from $TEAGE\n";
	$AGE = TEshuffle::load_TEage($TEAGE,$V);
}

#Now features to shuffle (need to be after in case there was $OKSEQ loaded)
print STDERR " --- checking file in -s, print in .bed if not a .bed or gff file\n";
print STDERR "     filtering TEs based on filter ($FILTER) and non TE behavior ($NONTE)\n" unless ($FILTER eq "na");
print STDERR "     + getting genomic counts for each repeat\n";
print STDERR "     + load all TE positions in a hash (since -s is set to rm)\n" if ($STYPE eq "rm") ;
my ($TOSHUFF_F,$PARSEDRM,$RM,$RM_C) = TEshuffle::RMtobed($SHUFFLE,$OKSEQ,$FILTER,$F_REGEX,$NONTE,$AGE,"y",$STYPE); #Note: $RM and $RM_C are empty unless $STYPE eq rm

#prep steps if shuffling type is tss
my ($TSSBED,$CLOSEST,$ALLTSS);
if ($STYPE eq "tss")  {
	#sort TEs
	my $bedsort = $BEDTOOLS."bedtools sort";
	my $sorted = $TOSHUFF_F;
	$sorted =~ s/\.bed$/\.sorted\.bed/;
	if (! -e $sorted) {
		print STDERR " --- sorting features of $TOSHUFF_F\n";	
		print STDERR "     $bedsort -i $TOSHUFF_F > $sorted\n";	
		`$bedsort -i $TOSHUFF_F > $sorted` ;
	}	
	$TOSHUFF_F = $sorted;
	print STDERR " --- loading the tss from $TSSFILE\n";
	#print the tss in a bed file => use bedtools closest
	($TSSBED,$ALLTSS) = TEshuffle::load_and_print_tss($TSSFILE);	
	print STDERR " --- sorting features in the tss file\n" unless (-e "$TSSBED.bed");
	print STDERR "     $bedsort -i $TSSBED > $TSSBED.bed\n" unless (-e "$TSSBED.bed");
	`$bedsort -i $TSSBED > $TSSBED.bed` unless (-e "$TSSBED.bed");
	$TSSBED = $TSSBED.".bed";
	#get the closest tss if relevant
	my $tssclosest = $TOSHUFF_F.".closest.tss";
	my $closestbed = $BEDTOOLS."closestBed";
	print STDERR " --- getting closest tss for each feature in $TOSHUFF_F, with the command line below\n";
	print STDERR "     $closestbed -a $TOSHUFF_F -b $TSSBED -D b -t first > $tssclosest\n"; #I want only one entry per TE, therefore -t first
	`$closestbed -a $TOSHUFF_F -b $TSSBED -D b -t first > $tssclosest`;
	print STDERR " --- loading distance to TSS\n";
	$CLOSEST = TEshuffle::load_closest_tss($tssclosest);
}

#Get total number of features in input file (= counting number of lines with stuff in it)
print STDERR " --- getting number and length of input\n";
my %IN_FEAT = ();
get_features_info($IN);
print STDERR "     number of features = $IN_FEAT{'nb'}\n";

#Join -i file with -s
my $INTERSectBed = $BEDTOOLS."intersectBed";
print STDERR " --- intersecting with command lines:\n";
print STDERR "        $INTERSectBed -a $TOSHUFF_F -b $IN -wo > $DIR/no_boot.joined\n";
`$INTERSectBed -a $TOSHUFF_F -b $IN -wo > $DIR/no_boot.joined`;

#Process the joined files
print STDERR " --- checking intersections of $IN with features in $TOSHUFF_F (observed)\n";
my $OBS = ();
$OBS = check_for_overlap("$DIR/no_boot.joined","no_boot",$OUT,$OBS);

#Now bootstrap runs
print STDERR " --- running $NBOOT bootstraps now (to get significance of the overlaps)\n";
if ($STYPE eq "rm") {
	print STDERR "        with intersect command line similar to the one above,\n";
	print STDERR "        and TEs shuffled among ALL TE positions, keeping length info\n";
} elsif ($STYPE eq "tss") {
	print STDERR "        with intersect command line similar to the one above,\n";
	print STDERR "        and TEs shuffled keeping same distance to TSS (using $TSSBED)\n";
} elsif ($STYPE eq "bed")  {
	print STDERR "        with intersect command line similar to the one above, and shuffle command line:\n";
	if ($INCL eq "na") {
		print STDERR "        ".$BEDTOOLS."shuffleBed -i $TOSHUFF_F -excl $EXCL -f 2 $NOOVERLAPS -g $BUILD -chrom -maxTries 10000\n";
    } else {
    	print STDERR "        ".$BEDTOOLS."shuffleBed -incl $INCL -i $TOSHUFF_F -excl $EXCL -f 2 $NOOVERLAPS -g $BUILD -chrom -maxTries 10000\n";
	}
}

my $BOOTS = ();
if ($NBOOT > 0) {
	foreach (my $i = 1; $i <= $NBOOT; $i++) {
		print STDERR "     ..$i bootstraps done\n" if (($i == 10) || ($i == 100) || ($i == 1000) || (($i > 1000) && (substr($i/1000,-1,1) == 0)));	
		my $shuffled;
		$shuffled = TEshuffle::shuffle_tss($TOSHUFF_F,$TEMP,$i,$ALLTSS,$CLOSEST,$OKSEQ) if ($STYPE eq "tss");
		$shuffled = TEshuffle::shuffle_rm($TOSHUFF_F,$TEMP,$i,$RM,$RM_C,$OKSEQ) if ($STYPE eq "rm");	
		$shuffled = TEshuffle::shuffle_bed($TOSHUFF_F,$TEMP,$i,$EXCL,$INCL,$BUILD_FILE,$BEDTOOLS,$NOOVERLAPS) if ($STYPE eq "bed");
		`$INTERSectBed -a $shuffled -b $IN -wo > $TEMP/boot.$i.joined`;
		$BOOTS = check_for_overlap("$TEMP/boot.$i.joined","boot.".$i,$OUTB,$BOOTS);
		`cat $OUTB >> $OUTB.CAT.boot.txt` if (-e $OUTB);
		`rm -Rf $TEMP/boot.$i.joined $shuffled`; #these files are now not needed anymore, all is stored
	}
}

#Stats now
if ($NBOOT > 0) {
	print STDERR " --- getting and printing counts stats\n";
	print_stats();
	print STDERR "     Stats printed in: $STATS.txt\n";
	
	#save disk space
	print STDERR " --- saving disk space by deleting shuffled and joined files\n";
	`rm -Rf $TEMP`;	
}

#end
print STDERR " --- $SCRIPTNAME done\n";
print STDERR "\n";
exit;

#-------------------------------------------------------------------------------
#-------------------------------- SUBROUTINES ----------------------------------
#-------------------------------------------------------------------------------
sub get_features_info {
	my $nb = `grep -c -E "\\w" $IN`;
#	my $len = `more $file | awk '{SUM += (\$3-\$2)} END {print SUM}'`; #this assumes no overlaps, trusting user for now
	chomp($nb);
#	chomp($len);
	$IN_FEAT{'nb'} = $nb;
#	$IN_FEAT{'len'} = $len;				
	return 1;
}

#-------------------------------------------------------------------------------
sub check_for_overlap {
	#NB called in 2 situations, OBS and the bootstraps
	#   check_for_overlap("$TEMP/no_boot.joined","no_boot",$OUT,$INTERS,$AGE);
	#   $BOOTS = check_for_overlap("$TEMP/boot.$i.joined","boot.".$i,$OUTB,$INTERS,$AGE,$BOOTS);
	my ($file,$fileid,$out,$counts) = @_;
	my %check = ();
	open(my $fh, "<$file") or confess "\n   ERROR (sub check_for_overlap): could not open to read $file!\n";
	LINE: while(<$fh>){
		chomp(my $l = $_);
		#FYI:
		# chr1	4522383	4522590	1111;18.9;4.6;1.0;chr1;4522383;4522590;(190949381);-;B3;SINE/B2;(0);216;1;1923	.	-	chr1	4496315	4529218	[ID] [score] [strand] [overlap_len]
		my @l = split(/\s+/,$l);	
		next LINE unless ($l[-1] >= $INTERS);
		my @rm = split(";",$l[3]);
		my $Rnam = $rm[9];
		$Rnam =~ s/\//_/; #making sure no / in repeat name
		my ($Rcla,$Rfam) = TEshuffle::get_Rclass_Rfam($Rnam,$rm[10]);
		my $id = $l[9]."#".$l[6]."#".$l[7]."#".$l[8];
		$id = $id."#".$l[11] if ($l[11]);
		
		#Increment in the data structure, but only if feature not already counted for this category
		unless ($check{$id}{'tot'}) {
			$counts->{$fileid}{'tot'}{'tot'}{'tot'}{'tot'}++;
		}
		unless ($check{$id}{$Rcla}) {
			$counts->{$fileid}{$Rcla}{'tot'}{'tot'}{'tot'}++;	
		}
		unless ($check{$id}{$Rcla.$Rfam}) {
			$counts->{$fileid}{$Rcla}{$Rfam}{'tot'}{'tot'}++;
		}
		unless ($check{$id}{$Rcla.$Rfam.$Rnam}) {
			$counts->{$fileid}{$Rcla}{$Rfam}{$Rnam}{'tot'}++;	
		}
						
		#Need to check if a feature is counted several times in the upper classes
		$check{$id}{'tot'}=1;
		$check{$id}{$Rcla}=1;
		$check{$id}{$Rcla.$Rfam}=1;
		$check{$id}{$Rcla.$Rfam.$Rnam}=1;
		
		#Now, age categories if any
		if ($TEAGE ne "na") {
			if ($AGE->{$Rnam}) {
				my $cat1 = $AGE->{$Rnam}[4];
				my $cat2 = $AGE->{$Rnam}[5];
				unless ($check{$id}{'age'}) { #easier to load tot hit with these keys for the print_out sub
					$counts->{$fileid}{'age'}{'cat.1'}{'tot'}{'tot'}++; 
					$counts->{$fileid}{'age'}{'cat.2'}{'tot'}{'tot'}++;
				}
				unless ($check{$id}{$cat1}) {
					$counts->{$fileid}{'age'}{'cat.1'}{$cat1}{'tot'}++;
				}
				if ($cat2 && ! $check{$id}{$cat2}) {
					$counts->{$fileid}{'age'}{'cat.2'}{$cat2}{'tot'}++;
				}
				$check{$id}{'age'}=1;
				$check{$id}{$cat1}=1;
				$check{$id}{$cat2}=1;		
			} else {
				print STDERR "        WARN: $Rnam not in $TEAGE?\n" unless ($check{$Rnam});
				$check{$Rnam} = 1;
			}
		}
	}
	close ($fh);		
	#Print only if observed; the hash for bootstraps is incrementing
	if ($fileid eq "no_boot") {
		print STDERR "     printing details in $out\n";	
		print_out($counts,$fileid,$out);
	}
	return $counts;
}

#-------------------------------------------------------------------------------
sub print_out {
	my ($counts,$fileid,$out) = @_;	
	foreach my $Rclass (keys %{$counts->{$fileid}}) {
		print_out_sub($fileid,$Rclass,"tot","tot",$counts,$out.".Rclass") if ($Rclass ne "age");
		foreach my $Rfam (keys %{$counts->{$fileid}{$Rclass}}) {
			print_out_sub($fileid,$Rclass,$Rfam,"tot",$counts,$out.".Rfam") if ($Rclass ne "age");			
			foreach my $Rname (keys %{$counts->{$fileid}{$Rclass}{$Rfam}}) {					
				print_out_sub($fileid,$Rclass,$Rfam,$Rname,$counts,$out.".Rname") if ($Rclass ne "age");
				print_out_sub($fileid,$Rclass,$Rfam,$Rname,$counts,$out.".age1") if (($Rclass eq "age") && ($Rfam eq "cat.1"));				
				print_out_sub($fileid,$Rclass,$Rfam,$Rname,$counts,$out.".age2") if (($Rclass eq "age") && ($Rfam eq "cat.2"));
			}
		}
	}
    return 1;
}

#-------------------------------------------------------------------------------
sub print_out_sub {
	my ($fileid,$key1,$key2,$key3,$counts,$out) = @_;
	my $tothit = $counts->{$fileid}{'tot'}{'tot'}{'tot'}{'tot'};
	my $hit = $counts->{$fileid}{$key1}{$key2}{$key3}{'tot'};	
	my $unhit = $IN_FEAT{'nb'} - $hit;
#	my $len = $counts->{$fileid}{$key1}{$key2}{$key3}{'len'}{'tot'};
	open (my $fh, ">>", $out) or confess "ERROR (sub print_out_sub): can't open to write $out $!\n";
				#fileid, class, fam, name, hits, total features loaded, unhit feat, total feat hit all categories
	print $fh "$fileid\t$key1\t$key2\t$key3\t$hit\t$IN_FEAT{'nb'}\t$unhit\t$tothit\n";
	close $fh;
    return 1;
}

#-------------------------------------------------------------------------------
sub print_stats {
	#get the bootstrapped avg values, sd, agregate all values
	my $exp = get_stats_data();
	$exp = TEshuffle::binomial_test_R($exp,"bed");
	
	#now print; permutation test + binomial test with avg lengths
	my $midval = $NBOOT/2;
	open (my $fh, ">", $STATS.".txt") or confess "ERROR (sub print_stats): can't open to write $STATS.txt $!\n";	
	print $fh "#Script $SCRIPTNAME, v$VERSION\n";
	print $fh "#Aggregated results + stats\n";
	print $fh "#Features in input file (counts):\n\t$IN_FEAT{'nb'}\n";
	print $fh "#With $NBOOT bootstraps for exp (expected); sd = standard deviation; nb = number; len = length; avg = average\n";
	print $fh "#Two tests are made (permutation and binomial) to assess how significant the difference between observed and random, so two pvalues are given\n";
	print $fh "#For the two tailed permutation test:\n";
	print $fh "#if rank is < $midval and pvalue is not \"ns\", there are significantly fewer observed values than expected \n";
	print $fh "#if rank is > $midval and pvalue is not \"ns\", there are significantly higher observed values than expected \n";
	print $fh "#The binomial test is done with binom.test from R, two sided\n";
	
	print $fh "\n#Level_(tot_means_all)\t#\t#\t#COUNTS\t#\t#\t#\t#\t#\t#\t#\t#\t#\t#\t#\t#\t#\n";
	print $fh "#Rclass\tRfam\tRname\tobs_hits\t%_obs_(%of_features)\tobs_tot_hits\tnb_of_trials(nb_of_TE_in_genome)\texp_avg_hits\texp_sd\t%_exp_(%of_features)\texp_tot_hits(avg)\tobs_rank_in_exp\t2-tailed_permutation-test_pvalue(obs.vs.exp)\tsignificance\tbinomal_test_proba\tbinomial_test_95%_confidence_interval\tbinomial_test_pval\tsignificance\n\n";
	foreach my $Rclass (keys %{$exp}) { #loop on all the repeat classes; if not in the obs then it will be 0 for obs values			
		foreach my $Rfam (keys %{$exp->{$Rclass}}) {			
			foreach my $Rname (keys %{$exp->{$Rclass}{$Rfam}}) {
				#observed
				my ($obsnb,$obsper) = (0,0);
				$obsnb = $OBS->{'no_boot'}{$Rclass}{$Rfam}{$Rname}{'tot'} if ($OBS->{'no_boot'}{$Rclass}{$Rfam}{$Rname}{'tot'});
				$obsper = $obsnb/$IN_FEAT{'nb'}*100 unless ($obsnb == 0);
				#expected
				my $expper = 0;
				my $expavg = $exp->{$Rclass}{$Rfam}{$Rname}{'avg'};	
				$expper = $expavg/$IN_FEAT{'nb'}*100 unless ($expavg == 0);
				#stats
				my $pval_nb = $exp->{$Rclass}{$Rfam}{$Rname}{'pval'};		
				$pval_nb = "na" if ($expavg == 0 && $obsnb == 0);									
				#Now print stuff
				print $fh "$Rclass\t$Rfam\t$Rname\t";
				print $fh "$obsnb\t$obsper\t$OBS->{'no_boot'}{'tot'}{'tot'}{'tot'}{'tot'}\t"; 
				print $fh "$PARSEDRM->{$Rclass}{$Rfam}{$Rname}\t"; 
				print $fh "$expavg\t$exp->{$Rclass}{$Rfam}{$Rname}{'sd'}\t$expper\t$exp->{'tot'}{'tot'}{'tot'}{'avg'}\t";			
				my $sign = TEshuffle::get_sign($pval_nb);				
				print $fh "$exp->{$Rclass}{$Rfam}{$Rname}{'rank'}\t$pval_nb\t$sign\t";								
				#Binomial
				my $b_pval = $exp->{$Rclass}{$Rfam}{$Rname}{'binom_pval'};
				$b_pval = "na" unless ($b_pval =~ /\d/);
				$sign = TEshuffle::get_sign($b_pval);
				print $fh "$exp->{$Rclass}{$Rfam}{$Rname}{'binom_prob'}\t$exp->{$Rclass}{$Rfam}{$Rname}{'binom_conf'}\t$exp->{$Rclass}{$Rfam}{$Rname}{'binom_pval'}\t$sign\n";	
			}
		}
	}
close $fh;
    return 1;
}

#-------------------------------------------------------------------------------
sub get_stats_data {
	my $exp = initialize_exp(); #0 values for all the ones seen in obs => so that even if not seen in exp, will be there

	#agregate data
	my ($nb_c,$nb_f,$nb_r,$nb_a1,$nb_a2) = ();
	foreach my $round (keys %{$BOOTS}) {
		foreach my $Rclass (keys %{$BOOTS->{$round}}) {
			push(@{$nb_c->{$Rclass}{'tot'}{'tot'}},$BOOTS->{$round}{$Rclass}{'tot'}{'tot'}{'tot'}) if ($Rclass ne "age");	
			foreach my $Rfam (keys %{$BOOTS->{$round}{$Rclass}}) {
				push(@{$nb_f->{$Rclass}{$Rfam}{'tot'}},$BOOTS->{$round}{$Rclass}{$Rfam}{'tot'}{'tot'}) if ($Rclass ne "age");		
				foreach my $Rname (keys %{$BOOTS->{$round}{$Rclass}{$Rfam}}) {
					push(@{$nb_r->{$Rclass}{$Rfam}{$Rname}},$BOOTS->{$round}{$Rclass}{$Rfam}{$Rname}{'tot'}) if ($Rclass ne "age");	
					push(@{$nb_a1->{$Rclass}{$Rfam}{$Rname}},$BOOTS->{$round}{$Rclass}{$Rfam}{$Rname}{'tot'}) if (($Rclass eq "age") && ($Rfam eq "cat.1"));
					push(@{$nb_a2->{$Rclass}{$Rfam}{$Rname}},$BOOTS->{$round}{$Rclass}{$Rfam}{$Rname}{'tot'}) if (($Rclass eq "age") && ($Rfam eq "cat.2"));
				}
			}
		}		
	}
	
	#get avg, sd and p values now => load in new hash, that does not have the fileID
	foreach my $round (keys %{$BOOTS}) {
		foreach my $Rclass (keys %{$BOOTS->{$round}}) {
			$exp = get_stats_data_details($Rclass,"tot","tot",$nb_c->{$Rclass}{'tot'}{'tot'},$exp) if ($Rclass ne "age");	
			foreach my $Rfam (keys %{$BOOTS->{$round}{$Rclass}}) {
				$exp = get_stats_data_details($Rclass,$Rfam,"tot",$nb_f->{$Rclass}{$Rfam}{'tot'},$exp) if ($Rclass ne "age");	
				foreach my $Rname (keys %{$BOOTS->{$round}{$Rclass}{$Rfam}}) {
					$exp = get_stats_data_details($Rclass,$Rfam,$Rname,$nb_r->{$Rclass}{$Rfam}{$Rname},$exp) if ($Rclass ne "age");
					$exp = get_stats_data_details($Rclass,$Rfam,$Rname,$nb_a1->{$Rclass}{$Rfam}{$Rname},$exp) if (($Rclass eq "age") && ($Rfam eq "cat.1"));
					$exp = get_stats_data_details($Rclass,$Rfam,$Rname,$nb_a2->{$Rclass}{$Rfam}{$Rname},$exp) if (($Rclass eq "age") && ($Rfam eq "cat.2"));
				}
			}
		}		
	}	
	return($exp);
}

#-------------------------------------------------------------------------------
sub initialize_exp {
	my $exp = ();
	#obs:
	#$counts->{$fileid}{'age'}{'cat.1'}{'tot'}{'tot'}
	foreach my $f (keys %{$OBS}) {
		foreach my $k1 (keys %{$OBS->{$f}}) {
			foreach my $k2 (keys %{$OBS->{$f}{$k1}}) {
				foreach my $k3 (keys %{$OBS->{$f}{$k1}{$k2}}) {
					$exp->{$k1}{$k2}{$k3}{'avg'}=0;
					$exp->{$k1}{$k2}{$k3}{'sd'}=0;
					$exp->{$k1}{$k2}{$k3}{'p'}=0; 
					$exp->{$k1}{$k2}{$k3}{'rank'}="na";
					$exp->{$k1}{$k2}{$k3}{'pval'}="na";
				}
			}
		}
	}
	return $exp;
}

#-------------------------------------------------------------------------------	
sub get_stats_data_details {
	my ($key1,$key2,$key3,$agg_data,$exp) = @_;	
	#get average and sd of the expected
	($exp->{$key1}{$key2}{$key3}{'avg'},$exp->{$key1}{$key2}{$key3}{'sd'}) = TEshuffle::get_avg_and_sd($agg_data);
	
	my $obs = $OBS->{'no_boot'}{$key1}{$key2}{$key3}{'tot'};
#	print STDERR "FYI: no observed value for {$key1}{$key2}{$key3}{'tot'}\n" unless ($obs);
	$obs = 0 unless ($obs);	
	
	#Get the rank of the observed value in the list of expected + pvalue for the permutation test
	my $rank = 1; #pvalue can't be 0, so I have to start there - that does mean there will be a rank nboot+1
	my @data = sort {$a <=> $b} @{$agg_data};
	EXP: foreach my $exp (@data) {
		last EXP if ($exp > $obs);
		$rank++;
	}	
	$exp->{$key1}{$key2}{$key3}{'rank'}=$rank;
	if ($rank <= $NBOOT/2) {
		$exp->{$key1}{$key2}{$key3}{'pval'}=$rank/$NBOOT*2;
	} else {
		$exp->{$key1}{$key2}{$key3}{'pval'}=($NBOOT+2-$rank)/$NBOOT*2; #+2 so it is symetrical (around nboot+1)
	}
	
	#Binomial test
	#get all the values needed for binomial test in R => do them all at once
	my $n = $PARSEDRM->{$key1}{$key2}{$key3} if ($PARSEDRM->{$key1}{$key2}{$key3});
	$n = 0 unless ($n);
#	print STDERR "        WARN: no value for total number (from RM output), for {$key1}{$key2}{$key3}? => no binomial test\n" if ($n == 0);
	my $p = 0;		
#	print STDERR "        WARN: n must be > x but n = $n and x = $obs for {$key1}{$key2}{$key3} => no binomial test\n" if ($n < $obs);	
	$p=$exp->{$key1}{$key2}{$key3}{'avg'}/$n unless ($n == 0); #should not happen, but could
	$exp->{$key1}{$key2}{$key3}{'p'} = $p;
	$exp->{$key1}{$key2}{$key3}{'n'} = $n;
	$exp->{$key1}{$key2}{$key3}{'x'} = $obs;
	return($exp);
}



