#!/usr/bin/perl
#######################################################
# Author  :  Zev Kronenberg (https://github.com/zeeev), for the v1.0
#            Modifications by Aurelie Kapusta (https://github.com/4ureliek) after v2.0
# email   :  4urelie.k@gmail.com  
# Purpose :  Originally writen to generate data (observed vs expected) shown in Figure 5, 
#            Kapusta et al. 2013 PLoS Genetics
#            (http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1003470)
#            But highly modified since then to be more useful
#######################################################
use strict;
use warnings;
use Carp;
use Getopt::Long;
use Bio::SeqIO;

use vars qw($BIN);
use Cwd 'abs_path';
BEGIN { 	
	$BIN = abs_path($0);
	$BIN =~ s/(.*)\/.*$/$1/;
	unshift(@INC, "$BIN/Lib");
}

use Statistics::R; #required to get the Binomial Test p-values for the TE stuff
use GAL::Annotation; #if issues, there is an alternative subroutine not using this module, see usage
use TEshuffle;
#use Data::Dumper;

#-----------------------------------------------------------------------------
#------------------------------- DESCRIPTION ---------------------------------
#-----------------------------------------------------------------------------
#flush buffer
$| = 1;

my $VERSION = "6.5";
my $SCRIPTNAME = "TE-analysis_Shuffle.pl";
my $CHANGELOG;
set_chlog();
sub set_chlog {
	$CHANGELOG = "
#	- v1.0 = 2012
#	- v2.0 = Jan 11 2016
#		Too many changes to list here, but concept remains the same
#	- v3.0 = Jan 26 2016
#       Change to bedtools instead of using the Set::IntervalTree perl module. 
#       This means there will be many files printed, but will be much faster
#       However this also means that there will not be a verification that transcripts
#          are located in an acceptable range (not in the -excl file). 
#          Assumes they are OK. Should be though.
#	- v3.1 = Jan 29 2016
#       Delete intermediate files, for space
#	- v3.2 = Feb 03 2016
#       Bug fix in rank (and therefore p value), was inverted (1000 instead of 1)
#       Delete the temp folders
#	- v3.3 = Feb 08/10 2016
#       Count transcript hits as well, as a category
#       Two-tailed test
#       Bug fix for -m 1 (would return only 0s for observed values)
#       Print more stuff in the stats.txt file so that no need of R
#	- v4.0 = Mar 28 2016
#       Allow skipping one of the inputs [make -p OR -l non mandatory]
#       Get results by repeats, like TE-analysis_Shuffle_bed_v2+.pl
#       + basically integrate all the improvements made there:
#         = Add binomial test as well
#         = Filter on some repeats
#         = possibility of several files to -e and -i
#	- v4.1 = Mar 31 2016
#       Few bug fix
#	- v4.2 = Apr 5 2016
#       Bug fix in stats by repeats (use of -f)
#       Correct rank for permutation when last rank (pvalue can't be 0)
#	- v5.0 = Oct 25 2016
#       Bug fix in stats for permutations 
#       Use R for binomial test
#       Get enrichment by age categories if age file provided
#       TEshuffle.pm for subroutines shared with the shuffle_bed script
#	- v6.0 = Mar 09 2017 [debugged]
#      different choices to shuffle the TEs:
#         - shufflebed = completely random positions, but same chromosome
#         - shuffle inside the current TE positions, same chromosome
#         - shuffle each TE, keeping its distance to a TSS, same chromosome - thanks to: Cedric Feschotte, Ed Chuong
#      make subfolders for each input file (for the shuffled outputs)
#      Also report observed values even if 0 in expected (interesting to see them in the obs, even if no stats possible)
#      Needed to round \$x for the binomial test in R (closest integer)
#	- v6.1 = Apr 04 2017
#      the 'transcript' cat is in fact genes => fix the total counts to get proper %
#	- v6.2 = Dec 01 2017
#      Bug fix for when long TEs were shuffled to the position of small TEs that are
#         too close to the start of the genomic sequence (led to negative starts).
#         This is now checked for -s rm and -s tss:
#            for -s rm, the TE is shifted of as many bp as needed
#            for -s tss the start is simply changed to 1 (to avoid having a TE placed closer to a tss)
#      Also added the option to use -r file in -s rm as well (to check ends) and shift the TE if needed.
#	- v6.3 = Mar 08 2018
#      Change for -s tss: if the TE ends up out of the scaffold/chr, put on the other side and if still out, 
#         shift it to be inside
#      Also, skip if not in the annotation file
#      Make -r mandatory for all
#	- v6.4 = Mar 08 2019
#      Option to keep the expected values, so the distributions can be plotted, and standardized,
#         to compare observed values (make it default)
#      Add bedtools version in log
#      Minor cosmetic stuff
#	- v6.5 = Mar 26 2019
#      Add a column with the correct count of features, not just the total count of exons
#      Bug fix - TSS_polyA was likely underestimated because could be overwritten...
\n";
	return 1;
}

my $USAGE;
set_usage();
sub set_usage {
	$USAGE = "
Synopsis (v$VERSION):

    perl $SCRIPTNAME -l lncRNA.gff [-p prot.gff] [-o <nt>] [-m <nb>] -q features_to_shuffle [-n <nb>] 
            -s shuffling_type -r <genome.sizes> [-b] 
            [-a <annotations>] [-f]
            [-e <genome.gaps>] [-d] [-i <include.range>] [-x] 
            [-w <bedtools_path>] [-u <no_low>] [-t <type,name>] [-c] [-g <TE.age.tab>] [-v] [-h]

    /!\\ REQUIRES
            - Bedtools, v2.25+
            - GAL::Annotation version later than Jan 2016 [update of is_coding]
              see https://github.com/The-Sequence-Ontology/GAL
              If issues with it, set the --just option to load transcripts without GAL
    /!\\ Previous outputs, if any, will be moved as *.previous (which only saves results once)

    Typically, for the 3 types, the mandatory arguments are as follow (one of -l or -p is required):
    perl $SCRIPTNAME -l|-p my_data.gff -q rm.out -r genome.sizes -s rm 
    perl $SCRIPTNAME -l|-p my_data.gff -q rm.out -r genome.sizes -s tss -a annotations.gtf
    perl $SCRIPTNAME -l|-p my_data.gff -q rm.out -r genome.sizes -s bed -r genome.range -e genome.gaps

	Note that -r is advised for -s rm (but won't affect -s tss)

  CITATIONS:
    - Include the version of the script + link to the GitHub page (https://github.com/4ureliek/TEanalysis)
    - Cite Kapusta et al. (2013) PLoS Genetics (DOI: 10.1371/journal.pgen.1003470)
      (http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1003470)
    - For BEDtools, Quinlan AR and Hall IM (2010) Bioinformatics (DOI: 10.1093/bioinformatics/btq033)

  SOME HISTORY:
    Originally writen by Zev Kronenberg, PhD (https://github.com/zeeev) to generate 
    data (observed vs expected) shown in Figure 5 of Kapusta et al. 2013 PLoS Genetics
    But highly modified since then by Aurelie Kapusta.
    Thanks to Edward Chuong, PhD, for the suggestion of the shuffling '-s tss', 
    and thanks to E. Chuong, Cedric Feschotte PhD and Javier Hernandez PhD, for helpful 
    discussions during the development of this script.

  DESCRIPTION:
    Features provided in -s will be overlapped with -p and/or -l files, without (no_boot) 
    or with (boot) shuffling (on same chromosome). One of -p or -l is mandatory. H
    aving both in the same run means that they are intersected with the same TE files, 
    which may be better for comparisons, but does not seem necessary with high bootstraps.
       
    A random transcript per gene is selected: use -m to do several repetitions of no_boot
    
    For each bootstrap (-n) with shuffling features in -s, transcripts are randomly selected as well
    Note that high bootstraps takes a lot of time.
    Shuffling is done by default with allowing overlaps between shuffled features,
    because it is faster and OK when over representation of specific repeats are considered.  
    Note that because TEs are often fragmented + there are inversions, the counts 
    for the exonized TEs is likely inflated; this also means that when TEs are shuffled, 
    there are more fragments than TEs. Some should be moved non independently, 
    or the input file should be corrected when possible to limit that issue 
    [not implemented in this script for now]
    
    Note that one exon may have several types of overlaps (e.g. \"SPL\" and \"exonized\"),
    but each exon is counted only one time for each category (important for \"exonized\").
    Similarly for TEs, each hit is counted unless it's the same repeat name / family / class (depending on the level)
   
    If you need to generate the <genome.gaps> file but you would also like to add more files to the -e option, 
    just do a first run with no bootstraps (in this example the genome.range is also being generated):
    perl ~/bin/$SCRIPTNAME -l input.gtf -q genome.out -r genome.fa -b -e genome.fa -d -n 0

    Two-tailed permutation test is done on the counts of overlaps for categories
    and the results are in a *.stats.cat.txt file
    If -f is used then stats are also made on each repeat, with two-tailed 
    permutation and binomial tests and the results are in a *.stats.TE.txt file.
    Note that the output *.stats.cat.txt is basically included in the output *.stats.TE.txt,
    with values of tot tot tot in the columns Rclass, Rfam and Rname
    The use of -f will take longer but requires fewer bootsraps, 
    because binomial test is more sensitive.  
  
  MANDATORY ARGUMENTS:
    -p,--prot     => (STRING) protein coding gff3 file; one of -p or -l is mandatory
    -l,--lnc      => (STRING) lncRNAs gff3 file; one of -p or -l is mandatory
    -q,--query    => (STRING) Features to shuffle = TE file
                              Repeat masker .out or the .bed file generated by the TE-analysis_pipeline                            
    -s,--shuffle  => (STRING) Shuffling type. Should be one of the following:
                               -t bed  => use bedtools shuffle, random position on the same chromosome
                               -t rm   => shuffle inside the current TE positions; still random, but less
                               -t tss  => shuffle the TEs on the same chromosome keeping their distance
                                          to the closest TSS of an annotation file provided
                                          (TSS are shuffled and then assigned to a TE; 
                                          same TSS will be assigned multiple TEs if fewer TSS than TEs)
                                          thanks to: Cedric Feschotte and Edward Chuong for the idea
     -r,--range    => (STRING) To know the maximum value in a given chromosome/scaffold. 
                              File should be: Name \\t length
                              Can be files from UCSC, files *.chrom.sizes
                              If you don't have such file, use -b (--build) and provide the genome fasta file for -r
                                                                                    
   MANDATORY ARGUMENTS IF USING -s tss:
    -a,--annot    => (STRING) gtf or gff; annotations to load all unique TSS: will be used to set the distance
                              between each TE and the closest TSS that will be kept while randomization
                              Note that it requires transcript lines
    -j,--just     => (BOOL)   GAL::Annotation seems to not really work on gtf files (won't load transcripts)
                              When that happens, the script will switch to alternative way of loading
                              transcripts. Or set this option directly to skip using GAL.

   MANDATORY ARGUMENTS IF USING -s bed:                               
    -e,--excl     => (STRING) This will be used as -excl for bedtools shuffle: 
                              \"coordinates in which features from -i should not be placed.\"
                              More than one file may be provided (comma separated), they will be concatenated 
                              (in a file = first-file-name.cat.bed).
                              By default, at least one file is required = assembly gaps, and it needs to be the first file
                              if not in bed format. Indeed, you may provide the UCSC gap file, with columns as:
                                  bin, chrom, chromStart, chromEnd, ix, n, size, type, bridge
                              it will be converted to a bed file. 
                              If you do not have this file, you may provide the genome file in fasta format
                              and add the option -d (--dogaps), to generate a bed file corresponding to assembly gaps.
                              If you need to generate the <genome.gaps> file but you would also 
                              like to add more files to the -e option, just do a first run with 
                              no bootstraps (in this example the genome.range is also being generated):
                                 perl ~/bin/$SCRIPTNAME -l lncRNA.gff -p prot.gff -q genome.out -s rm -r genome.fa -b -e genome.fa -d -n 0                                
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
       
   OTHER OPTIONAL ARGUMENTS (for all -s):
    -b,--build    => (BOOL)   See above; use this and provide the genome fasta file if no range/lengths file (-r)
                              This step may take a while but will create the required file	
    -o,--overlap  => (INT)    Minimal length (in nt) of intersection in order to consider the TE included in the feature.
                              Default = 10 (to match the TEanalysis-pipeline.pl)
    -m,--more     => (INT)    Even in the no_boot, a random transcript is picked. Set this number to do repetitions for no_boot.
                              Default = 1 (still need it done 1 time; set this to 0 is equivalent to 1)
                              For binomial test, the observed value will be the average, rounded to closest integer
    -n,--nboot    => (STRING) number of bootsraps with shuffled -s file
                              Default = 100 for faster runs; use higher -n for good pvalues 
                              (-n 10000 is best for permutation test but this will take a while)
                              If set to 0, no bootstrap will be done
    -f,--full     => (BOOL)   Use -f to also do stats for each repeat separately (separated output, with binomial test as well)
                              Results will be in a file *.stats.TE.txt
                              Note that the output *.stats.cat.txt is basically included in the output *.stats.TE.txt,
                              with values of tot tot tot in the columns Rclass, Rfam and Rname
                              This will take longer but requires fewer bootsraps, because binomial test is more sensitive
    -w,--where    => (STRING) if BEDtools are not in your path, provide path to BEDtools bin directory                             

   OPTIONAL ARGUMENTS FOR TE FILTERING (for all -s): 
    -k,--keep     => (STRING) To set which intermediate files to keep / print.
                                 -k all  = print the values to plot distribution + keep all intermediate files
                                 -k dist [default] = print the values to plot distribution
                                 -k none = everything is deleted
    -u,--u        => (STRING) To set the behavior regarding non TE sequences: all, no_low, no_nonTE, none
                                 -u all = keep all non TE sequences (no filtering)
                                 -u no_low [default] = keep all besides low_complexity and simple_repeat
                                 -u no_nonTE = keep all except when class = nonTE
                                 -u none = everything is filtered out 
                                    (nonTE, low_complexity, simple_repeat, snRNA, srpRNA, rRNA, tRNA/tRNA, satellite)
    -t,--te       => (STRING) <type,name>
                              run the script on only a subset of repeats. Not case sensitive.
                              The type can be: name, class or family and it will be EXACT MATCH unless -c is chosen as well
                              ex: -t name,nhAT1_ML => only fragments corresponding to the repeat named exactly nhAT1_ML will be looked at
                                  -t class,DNA => all repeats with class named exactly DNA (as in ...#DNA/hAT or ...#DNA/Tc1)
                                  -t family,hAT => all repeats with family named exactly hAT (so NOT ...#DNA/hAT-Charlie for example)
    -c,--contain  => (BOOL)   to check if the \"name\" determined with -filter is included in 
                              the value in Repeat Masker output, instead of exact match
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

#-----------------------------------------------------------------------------
#------------------------------ LOAD AND CHECK -------------------------------
#-----------------------------------------------------------------------------
my ($SHUFFLE,$STYPE,$TSSFILE,$JUST,$FULL,$EXCLUDE,$DOGAPS,$BUILD,$DOBUILD,$F_REGEX,$ALLOW,$NOOVERLAPS,$V,$HELP);
my ($PROT,$LINC) = ("n","n");
my $INTERS = 10;
my $MORE = 0;
my $NBOOT = 10;
my $INCL = "na";
my $NONTE = "no_low";
my $FILTER = "na";
my $TEAGE = "na";
my $BEDTOOLS = "";
my $KEEP = "dist";
my $CATOUT = "y"; #removed from options, not really relevant to ask for choice
my $OPT_SUCCESS = GetOptions(
			 	  'prot=s'		=> \$PROT,
			 	  'lnc=s'		=> \$LINC,
			 	  'more=s'      => \$MORE,
			 	  'query=s'     => \$SHUFFLE,
			 	  'shuffle=s'   => \$STYPE,
			 	  'annot=s'     => \$TSSFILE,
			 	  'just'        => \$JUST,				 	  
			 	  'overlap=s'   => \$INTERS,
			 	  'nboot=s'     => \$NBOOT,
			 	  'full'        => \$FULL,
			 	  'range=s'     => \$BUILD,
			 	  'build'       => \$DOBUILD,
			 	  'excl=s'		=> \$EXCLUDE,
			 	  'dogaps'      => \$DOGAPS,
			 	  'incl=s'		=> \$INCL,
			 	  'x'		    => \$NOOVERLAPS,
			 	  'u=s'		    => \$NONTE,
			 	  'te=s'		=> \$FILTER,
			 	  'contain'     => \$F_REGEX,
			 	  'group=s'     => \$TEAGE,
			 	  'where=s'     => \$BEDTOOLS,
			 	  'keep'        => \$KEEP,
			 	  'version'     => \$V,
			 	  'help'		=> \$HELP,);
			 	  
#Check options, if files exist, etc
die "\n --- $SCRIPTNAME version $VERSION\n\n" if $V;
die $USAGE if ($HELP);
die "\n SOME MANDATORY ARGUMENTS MISSING, CHECK USAGE:\n$USAGE" if ((! $LINC && ! $PROT) || ! $SHUFFLE || ! $STYPE || ! $BUILD);
die "\n One of -l or -p needs to be provided\n\n" if ($PROT eq "n" && $LINC eq "n" );
die "\n -t is required\n\n" if (! $STYPE);
die "\n -q is required\n\n" if (! $SHUFFLE);
die "\n -r $BUILD does not exist?\n\n" if (! -e $BUILD);

die "\n -p $PROT does not exist?\n\n"  if ($PROT ne "n" && ! -e $PROT);
die "\n -p $PROT is not a gff file?\n\n" unless ($PROT eq "n" || $PROT =~ /\.gff$/ || $PROT =~ /\.gff3$/);
die "\n -l $LINC does not exist?\n\n"  if ($LINC ne "n" && ! -e $LINC);
die "\n -l $LINC is not a gff file?\n\n" unless ($LINC eq "n" || $LINC =~ /\.gff$/ || $LINC =~ /\.gff3$/);

die "\n -q $SHUFFLE is not in a proper format? (not .out, .bed, .gff or .gff3)\n\n" unless ($SHUFFLE =~ /\.out$/ || $SHUFFLE =~ /\.bed$/ || $SHUFFLE =~ /\.gff$/ || $SHUFFLE =~ /\.gff3$/);
die "\n -q $SHUFFLE does not exist?\n\n" if (! -e $SHUFFLE);
die "\n -s $STYPE should be one of the following: bed, rm or tss\n\n" if ($STYPE ne "bed" && $STYPE ne "rm" && $STYPE ne "tss");
#deal with conditional mandatory stuff
die "\n -s tss was set, but -a is missing?\n\n" if ($STYPE eq "tss" && ! $TSSFILE);
die "\n -a $TSSFILE does not exist?\n\n" if ($TSSFILE && ! -e $TSSFILE);
if ($STYPE eq "bed") {
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

($FULL)?($FULL = "y"):($FULL = "n");
($DOGAPS)?($DOGAPS = "y"):($DOGAPS = "n");
($DOBUILD)?($DOBUILD = "y"):($DOBUILD = "n");
($F_REGEX)?($F_REGEX = "y"):($F_REGEX="n");
$BEDTOOLS = $BEDTOOLS."/" if ($BEDTOOLS ne "" && substr($BEDTOOLS,-1,1) ne "/"); #put the / at the end of path if not there
($NOOVERLAPS)?($NOOVERLAPS = "-noOverlapping"):($NOOVERLAPS = "");
$MORE = 1 if ($MORE == 0); #1 rep if set to 0, same thing here

#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------
#Prep steps
print STDERR "\n --- $SCRIPTNAME v$VERSION started, with:\n";
print STDERR "     input lncRNA file = $LINC\n" if ($LINC ne "n");
print STDERR "     input mRNA file = $PROT\n" if ($PROT ne "n");
print STDERR "     features to shuffle = $SHUFFLE\n";
print STDERR "     shuffling type = $STYPE\n";
my $BEDV = $BEDTOOLS."bedtools --version";
my $BEDVER = `$BEDV`;
chomp $BEDVER;
print STDERR "     bedtools version = $BEDVER\n";

#Outputs
print STDERR " --- prepping output directories and files\n";
my $INPUT;
($LINC)?($INPUT = $LINC):($INPUT = $PROT);
my $DIR = $INPUT.".shuffle-".$STYPE.".".$NBOOT;
print STDERR "     output directory = $DIR\n";
my ($STATS,$DISTRIB,$OUTL,$OUTLB,$TEMP_L,$OUTP,$OUTPB,$TEMP_B,$TEMP) = TEshuffle::prep_out("tr",$DIR,$NBOOT,$FILTER,$INPUT,$STYPE,$NONTE,$KEEP,$LINC,$PROT,$SHUFFLE);	

#Chosomosome sizes / Genome range
print STDERR " --- loading build (genome range)\n";
my ($OKSEQ,$BUILD_FILE) = TEshuffle::load_build($BUILD,$DOBUILD);

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
	($EXCLUDE =~ /,/)?($EXCL = TEshuffle::concat_beds(\@exclude)):($EXCL = $exclude[0]);

	#If relevant, files to include for shuffling
	if (($INCL ne "na") && ($INCL =~ /,/)) {
		print STDERR " --- concatenating $INCL files to one file\n";
		my @include = split(",",$INCL);
		$INCL = TEshuffle::concat_beds(\@include);
	}
}

#Load TEage if any
print STDERR " --- Loading TE ages from $TEAGE\n" unless ($TEAGE eq "na");
my $AGE = ();
$AGE = TEshuffle::load_TEage($TEAGE,$V) unless ($TEAGE eq "na");

#Now features to shuffle (need to be after in case there was $OKSEQ loaded)
print STDERR " --- checking file in -s, print in .bed if not a .bed or gff file\n";
print STDERR "     filtering TEs based on filter ($FILTER) and non TE behavior ($NONTE)\n" unless ($FILTER eq "na");
print STDERR "     + getting genomic counts for each repeat\n";
print STDERR "     + load all TE positions in a hash (since $STYPE is set to rm)\n" if ($STYPE eq "rm") ;
my ($TOSHUFF_FILE,$PARSEDRM,$RM,$RM_C) = TEshuffle::RMtobed($SHUFFLE,$OKSEQ,$FILTER,$F_REGEX,$NONTE,$AGE,"y",$STYPE); #Note: $RM and $RM_C are empty unless $STYPE eq rm

#prep steps if shuffling type is tss
my ($TSSBED,$CLOSEST,$ALLTSS);
if ($STYPE eq "tss")  {
	#sort TEs
	my $bedsort = $BEDTOOLS."bedtools sort";
	my $sorted = $TOSHUFF_FILE;
	$sorted =~ s/\.bed$/\.sorted\.bed/;
	print STDERR " --- sorting features of $TOSHUFF_FILE\n" unless (-e $sorted);	
	print STDERR "     $bedsort -i $TOSHUFF_FILE > $sorted\n" unless (-e $sorted);	
	`$bedsort -i $TOSHUFF_FILE > $sorted` unless (-e $sorted);
	$TOSHUFF_FILE = $sorted;
	print STDERR " --- loading the tss from $TSSFILE\n";
	#print the tss in a bed file => use bedtools closest
	($TSSBED,$ALLTSS) = TEshuffle::load_and_print_tss($TSSFILE);	
	print STDERR " --- sorting features in the tss file\n" unless (-e "$TSSBED.bed");
	print STDERR "     $bedsort -i $TSSBED > $TSSBED.bed\n" unless (-e "$TSSBED.bed");
	`$bedsort -i $TSSBED > $TSSBED.bed` unless (-e "$TSSBED.bed");
	$TSSBED = $TSSBED.".bed";
	#get the closest tss if relevant
	my $tssclosest = $TOSHUFF_FILE.".closest-tss.".TEshuffle::filename($TSSFILE);
	my $CLOSESTBed = $BEDTOOLS."closestBed";
	print STDERR " --- getting closest tss for each feature in $TOSHUFF_FILE\n";
	if (-e $tssclosest) {
		print STDERR "     $tssclosest exists, skipping\n";
	} else {
		print STDERR "     with the command line below\n";
		print STDERR "     $CLOSESTBed -a $TOSHUFF_FILE -b $TSSBED -D b -t first > $tssclosest\n"; #I want only one entry per TE, therefore -t first
		`$CLOSESTBed -a $TOSHUFF_FILE -b $TSSBED -D b -t first > $tssclosest`;
	}
	print STDERR " --- loading distance to TSS\n";
	$CLOSEST = TEshuffle::load_closest_tss($tssclosest);
}

#Load the gff file(s)
print STDERR " --- Load gene IDs / transcript IDs for:\n";
my %WHICHGENE = ();
my $L_TR = (); #trinfos
my $P_TR = (); #trinfos
my %FLAGGAL = (0 => 0, 1 => 0);
print STDERR "  -> $LINC\n" unless ($LINC eq "n");
($L_TR) = read_gff_gal($LINC,0) if ($LINC ne "n" && ! $JUST);
($L_TR) = load_gene_tr($LINC,0) if ($LINC ne "n" && ($JUST || $FLAGGAL{0} == 1)); #if GAL::Annotation is a problem
print STDERR "  -> $PROT\n" unless ($PROT eq "n");
($P_TR) = read_gff_gal($PROT,1) if ($PROT ne "n" && ! $JUST);
($P_TR) = load_gene_tr($PROT,1) if ($PROT ne "n" && ($JUST || $FLAGGAL{1} == 1)); #if GAL::Annotation is a problem

#Join -p and/or -l files
my $INTERSECTBED = $BEDTOOLS."intersectBed";
print STDERR " --- Intersect with command lines:\n";
print STDERR "      $INTERSECTBED -a $TOSHUFF_FILE -b $LINC -wo > $TEMP_L/no_boot.joined\n" unless ($LINC eq "n");
system "$INTERSECTBED -a $TOSHUFF_FILE -b $LINC -wo > $TEMP_L/no_boot.joined" unless ($LINC eq "n");
print STDERR "      $INTERSECTBED -a $TOSHUFF_FILE -b $PROT -wo > $TEMP_B/no_boot.joined\n" unless ($PROT eq "n");
system "$INTERSECTBED -a $TOSHUFF_FILE -b $PROT -wo > $TEMP_B/no_boot.joined" unless ($PROT eq "n");

#Process the joined files with -m X repeats
print STDERR " --- Check intersection(s) with features in $TOSHUFF_FILE (observed)\n";
print STDERR "     (if -m set, there will be several rounds of random transcript selection)\n";
my $NO_BOOT = ();
my $NO_BOOTS_TOT_EXONS = (); #will contain all the count info of the dataset for all categories
for(my $j = 1; $j <= $MORE; $j++) {
	print STDERR "     ..$j rounds done\n" if ($j == 10 || $j == 100 || $j == 1000 || ($j > 1000 && substr($j/1000,-1,1) == 0));	
	($NO_BOOT,$NO_BOOTS_TOT_EXONS) = 
		check_for_featured_overlap("$TEMP_L/no_boot.joined",$L_TR,"no_boot.".$j,'transcript',$OUTL,$NO_BOOT,$NO_BOOTS_TOT_EXONS) 
		unless ($LINC eq "n");
	($NO_BOOT,$NO_BOOTS_TOT_EXONS) = 
		check_for_featured_overlap("$TEMP_B/no_boot.joined",$P_TR,"no_boot.".$j,'mRNA',$OUTP,$NO_BOOT,$NO_BOOTS_TOT_EXONS) 
		unless ($PROT eq "n");
	`cat $OUTL >> $CATOUT.no-boot.txt` if ($CATOUT && -e $OUTL);
	`cat $OUTP >> $CATOUT.no-boot.txt` if ($CATOUT && -e $OUTP);
}

#Now bootstrap runs
print STDERR " --- Run $NBOOT bootstraps now (to get significance of the overlaps)\n";
my $BOOTS = ();
my $BOOTS_TOT_EXONS = (); #will contain all the count info of the dataset for all categories
if ($NBOOT > 0) {
	foreach (my $i = 1; $i <= $NBOOT; $i++) {
		print STDERR "     ..$i bootstraps done\n" if (($i == 10) || ($i == 100) || ($i == 1000) || (($i > 1000) && (substr($i/1000,-1,1) == 0)));	
		my $SHUFFLED;
		$SHUFFLED = TEshuffle::shuffle_tss($TOSHUFF_FILE,$TEMP,$i,$ALLTSS,$CLOSEST,$OKSEQ) if ($STYPE eq "tss");
		$SHUFFLED = TEshuffle::shuffle_rm($TOSHUFF_FILE,$TEMP,$i,$RM,$RM_C,$OKSEQ) if ($STYPE eq "rm");	
		$SHUFFLED = TEshuffle::shuffle_bed($TOSHUFF_FILE,$TEMP,$i,$EXCL,$INCL,$BUILD_FILE,$BEDTOOLS,$NOOVERLAPS) if ($STYPE eq "bed");
 		system "      $INTERSECTBED -a $SHUFFLED -b $LINC -wo > $TEMP_L/boot.$i.joined" unless ($LINC eq "n");
 		system "      $INTERSECTBED -a $SHUFFLED -b $PROT -wo > $TEMP_B/boot.$i.joined" unless ($PROT eq "n");	
		($BOOTS,$BOOTS_TOT_EXONS) = 
			check_for_featured_overlap("$TEMP_L/boot.$i.joined",$L_TR,"boot.".$i,'transcript',$OUTLB,$BOOTS,$BOOTS_TOT_EXONS) 
			unless ($LINC eq "n");
		($BOOTS,$BOOTS_TOT_EXONS) = 
			check_for_featured_overlap("$TEMP_B/boot.$i.joined",$P_TR,"boot.".$i,'mRNA',$OUTPB,$BOOTS,$BOOTS_TOT_EXONS) 
			unless ($PROT eq "n");		
		`cat $OUTLB >> $CATOUT.boot.txt` if ($CATOUT && -e $OUTLB);
		`cat $OUTPB >> $CATOUT.boot.txt` if ($CATOUT && -e $OUTPB);
		`rm -Rf $SHUFFLED` unless ($KEEP eq "all"); #these files are now not needed anymore, all is stored
		`rm -Rf $TEMP_L/boot.$i.joined` unless ($KEEP eq "all" || $LINC eq "n");
		`rm -Rf $TEMP_B/boot.$i.joined` unless ($KEEP eq "all" || $PROT eq "n");
	}
}
`rm -Rf $TEMP_L` unless ($KEEP eq "all" || $LINC eq "n");
`rm -Rf $TEMP_B` unless ($KEEP eq "all" || $PROT eq "n");
`rm -Rf $TEMP` unless ($KEEP eq "all");


#Gather all results and print outputs
print STDERR " --- Get and print stats\n" if ($NBOOT > 0);
if ($NBOOT > 0) {
	#get the boot and no_boot total_exons values, avg and sd
	print STDERR "     Get number of exons (total and hit)\n";
	my $no_boot_exons = get_exon_data($NO_BOOTS_TOT_EXONS);
	my $boot_exons = get_exon_data($BOOTS_TOT_EXONS);		
	#now print
	print_cat_data($no_boot_exons,$boot_exons);
	print_rep_data($no_boot_exons,$boot_exons) if ($FULL eq "y");
}

print STDERR " --- $SCRIPTNAME done\n";
print STDERR "     Stats for categeories printed in: $STATS.cat.txt\n" if ($NBOOT > 0);
print STDERR "     Stats for TEs printed in: $STATS.TE.txt\n" if ($NBOOT > 0 && $FULL eq "y");
print STDERR "\n";
exit;


#-----------------------------------------------------------------------------
#-------------------------------- SUBROUTINES --------------------------------
#-----------------------------------------------------------------------------
sub read_gff_gal {
	my ($gff3_file,$coding) = @_;
    my %trinfo;
	my $gene_count=0;
	#load annotations through GAL::Annotation
    my $annotation = GAL::Annotation->new($gff3_file);
    my $features = $annotation->features;
    my $genes = $features->search({type => 'gene'});
    print STDERR "     GAL::Annotation has finished loading, now going through it\n";
	my $type = "transcript";
	GENE: while (my $gene = $genes->next) {
		if($coding eq 1){
			next GENE unless $gene->is_coding; #function updated Jan 2016 by Barry Moore to return true if any child is mRNA or has CDS exons
		}
		my $gene_id = $gene->feature_id;
		my $seqid   = $gene->seqid;
		next GENE if ($STYPE eq "bed" && ! $OKSEQ->{$seqid}); #if not in build of stuff OK to shuffle on, remove here as well; only relevant for -s bed though
		my @tr = $gene->transcripts;
		TRANSCRIPT: foreach my $tr (@tr) {
			my $tr_id = $tr->feature_id;
			my $tr_strand = $tr->strand;
			if ($tr_strand !~ /\+|-/) {
				print STDERR "     Warning: transcript strand for $tr_id is undetermined ($tr_strand)\n";
				next TRANSCRIPT;
			}
			
			#Check if transcript is coding or not
			$type = "mRNA" if ($tr->has_CDS);
			my @exons = sort { $a->start <=> $b->start } $tr->exons;
# 			if ($type eq 'mRNA') {
# 				@exons = sort { $a->start <=> $b->start } $transcript->CDSs;
# 			} else {
# 				@exons = sort { $a->start <=> $b->start } $transcript->exons;
# 			}
			#Now get info of number of exons in this transcript
			$trinfo{$gene_id}{$type}{$tr_id}{'st'}=$tr->start;
			$trinfo{$gene_id}{$type}{$tr_id}{'en'}=$tr->end;
			$trinfo{$gene_id}{$type}{$tr_id}{'nb'}=scalar(@exons);
			$WHICHGENE{$tr_id}=$gene_id;
		}
		$gene_count++;		
	}
    print STDERR "        total genes loaded (type=$type): $gene_count\n";	
    $FLAGGAL{$coding} = 1 if ($gene_count == 0);
    print STDERR "        WARN: gene count = 0, transcript info will (try to) be loaded without GAL::Annotation\n" if ($gene_count == 0);
	return (\%trinfo);
}

#-----------------------------------------------------------------------------
sub load_gene_tr {
	#Not as solid as using GAL::Annotation, but it is an alternative, in case issues with GAL
	my ($file,$coding) = @_;
    print STDERR "     Loading transcripts / genes relationships without using GAL::Annotation\n";
    my $gene_count=0;
	my %trinfo = ();
	my ($gid,$trid,$type);
	open(my $fh, "<$file") or confess "\n   ERROR (sub load_gene_tr): could not open to read $file!\n";
	LINE: while(<$fh>) {
		chomp(my $l = $_);
		next LINE if (substr($l,0,1) eq "#");
		my @l = split('\s+',$l);
		next LINE if (($STYPE eq "bed") && (! $OKSEQ->{$l[0]})); #if not in build of stuff OK to shuffle on, remove here as well; only relevant for -s bed though
		my $id = $l[8];
		$id = $1 if $id =~ /^ID=(.+?);/;
		if ($l[2] eq "gene") {
			$gid = $id;
			$gene_count++;
		} elsif ($l[2] eq "transcript") {
			$trid = $id;
			$type = "transcript";
			$type = "mRNA" if ($l =~ /protein_coding/ || $coding == 1); #if coding, it should be in the "gene_type" or the "transcript_type", but if not assume when coding is set to 1, it's coding
			$trinfo{$gid}{$type}{$trid}{'st'} = $l[3];
			$trinfo{$gid}{$type}{$trid}{'en'} = $l[4];
			$WHICHGENE{$trid}=$gid;
		} else {	
			$trinfo{$gid}{$type}{$trid}{'nb'}++ if ($l[2] eq "exon"); #count exons; includes UTRs for pc genes
# 			$trinfo{$gid}{$type}{$trid}{'nb'}++ if (($type eq "mRNA") && ($l[2] eq "CDS")); #count number of coding exons only
#			$trinfo{$gid}{$type}{$trid}{'nb'}++ if (($type eq "transcript") && ($l[2] eq "exon")); #count number of exons
		}
	}
	close ($fh);
	print STDERR "        total genes loaded (type=$type): $gene_count\n";	
	return (\%trinfo); #looping through keys will get transcripts => put in an array for each gene later
}

#-----------------------------------------------------------------------------
sub check_for_featured_overlap {
	my ($file,$trinfo,$fileid,$type,$out,$counts,$total_exons) = @_;
	my %chosen_tr = ();
	my %check = ();
	my %checkTE = ();
	
	#now loop
	open(my $fh, "<$file") or confess "\n   ERROR (sub check_for_featured_overlap): could not open to read $file!\n";
	LINE: while(<$fh>){
		chomp(my $l = $_);
		next LINE if (substr($l,0,1) eq "#");
		my @l = split(/\s+/,$l);

#FYI:
# chr1	4522383	4522590	1111;18.9;4.6;1.0;chr1;4522383;4522590;(190949381);-;B3;SINE/B2;(0);216;1;1923	.	-	chr1	Cufflinks	gene	4496315	4529218	.	+	.	ID=XLOC_000001;Name=uc007aez.1;
# chr1	4522383	4522590	1111;18.9;4.6;1.0;chr1;4522383;4522590;(190949381);-;B3;SINE/B2;(0);216;1;1923	.	-	chr1	Cufflinks	transcript	4496316	4523815	.	+	.	ID=TCONS_00000002;Parent=XLOC_000001;
		
# 		if ($l[8] eq "transcript") {
# 			#TO DO: count intron hits when transcript hit but not exon hit, using a flag; for now it does not matter
# 		} elsif {
		if ($l[8] eq "exon") {
			my $tridf = $l[14];
			my $trid = $tridf;
			$trid = $1 if $trid =~ /Parent=(.+?);/;
			next LINE unless (defined $WHICHGENE{$trid}); #checked for non coding when coding are looked at
			
			#get a random tr for this gene, but only the first time this gene is met, and keep which tr is chosen			
			my $gid = $WHICHGENE{$trid};			
			$chosen_tr{$gid} = random_tr($trinfo,$gid,$type) unless (defined $chosen_tr{$gid});
			my $chosen = $chosen_tr{$gid};		
			next LINE if ($trid ne $chosen); #skip if current transcript is not the chosen one
			
			my $ilen = $l[-1]; #last value of the line is intersection length
			next LINE if ($ilen < $INTERS);

			#now check what category of overlap this exon is;
			my $cat = overlap_category(\@l,$trinfo,$gid,$type,$trid);
		
			#now increment in the data structure		
			#since only one transcript per gene, there should be no worry here about unique counts, 1 exon can only be counted one time in a category; 
			#however unique exon hits count need a check, and there could be TE overlaps fucking things up, so better safe than sorry	
			unless (defined $check{$tridf}{$cat}) { 
				($counts->{$cat}{$type}{$fileid}{'tot'}{'tot'}{'tot'}{'nr'})?
				($counts->{$cat}{$type}{$fileid}{'tot'}{'tot'}{'tot'}{'nr'}++):($counts->{$cat}{$type}{$fileid}{'tot'}{'tot'}{'tot'}{'nr'}=1);
			}
			$check{$tridf}{$cat}=1;
			unless (defined $check{$gid}{'hit'}) { #counting each tr hit only one time => equivalent to a number of genes, not transcripts
#			unless (defined $check{$chosen}{'hit'}) { #what I had before
				($counts->{'transcript'}{$type}{$fileid}{'tot'}{'tot'}{'tot'}{'nr'})?
				($counts->{'transcript'}{$type}{$fileid}{'tot'}{'tot'}{'tot'}{'nr'}++):($counts->{'transcript'}{$type}{$fileid}{'tot'}{'tot'}{'tot'}{'nr'}=1);
				#duplicate, but it's easier that way:
				($counts->{'transcript'}{$type}{$fileid}{'tot'}{'tot'}{'tot'}{'tot'})?
				($counts->{'transcript'}{$type}{$fileid}{'tot'}{'tot'}{'tot'}{'tot'}++):($counts->{'transcript'}{$type}{$fileid}{'tot'}{'tot'}{'tot'}{'tot'}=1); 
			}	
			$check{$gid}{'hit'}=1;
#			$check{$chosen}{'hit'}=1;
			
			#Do the repeats stuff if relevant
			unless ($FULL eq "n") {
				my @l = split(/\s+/,$l);	
				next LINE unless ($ilen >= $INTERS);
				my @rm = split(";",$l[3]);
				my $Rnam = $rm[9];
				my ($Rcla,$Rfam) = TEshuffle::get_Rclass_Rfam($Rnam,$rm[10]);
				#Increment in the data structure, but only if relevant = avoid counting hits several times
				unless ($checkTE{$tridf}{$cat}{$type}{'tot'}) {
					($counts->{$cat}{$type}{$fileid}{'tot'}{'tot'}{'tot'}{'tot'})?
					($counts->{$cat}{$type}{$fileid}{'tot'}{'tot'}{'tot'}{'tot'}++):($counts->{$cat}{$type}{$fileid}{'tot'}{'tot'}{'tot'}{'tot'}=1);
				}	
				unless ($checkTE{$tridf}{$cat}{$type}{$Rcla}) {
					($counts->{$cat}{$type}{$fileid}{$Rcla}{'tot'}{'tot'}{'tot'})?
					($counts->{$cat}{$type}{$fileid}{$Rcla}{'tot'}{'tot'}{'tot'}++):($counts->{$cat}{$type}{$fileid}{$Rcla}{'tot'}{'tot'}{'tot'}=1);			
				}
				unless ($checkTE{$tridf}{$cat}{$type}{$Rcla.$Rfam}) {
					($counts->{$cat}{$type}{$fileid}{$Rcla}{$Rfam}{'tot'}{'tot'})?
					($counts->{$cat}{$type}{$fileid}{$Rcla}{$Rfam}{'tot'}{'tot'}++):($counts->{$cat}{$type}{$fileid}{$Rcla}{$Rfam}{'tot'}{'tot'}=1);
				}
				unless ($checkTE{$tridf}{$cat}{$type}{$Rcla.$Rfam.$Rnam}) {
					($counts->{$cat}{$type}{$fileid}{$Rcla}{$Rfam}{$Rnam}{'tot'})?
					($counts->{$cat}{$type}{$fileid}{$Rcla}{$Rfam}{$Rnam}{'tot'}++):($counts->{$cat}{$type}{$fileid}{$Rcla}{$Rfam}{$Rnam}{'tot'}=1);	
				}				
				#Need to check if a feature is counted several times in the upper classes
				$checkTE{$tridf}{$cat}{$type}{'tot'}=1;
				$checkTE{$tridf}{$cat}{$type}{$Rcla}=1;
				$checkTE{$tridf}{$cat}{$type}{$Rcla.$Rfam}=1;
				$checkTE{$tridf}{$cat}{$type}{$Rcla.$Rfam.$Rnam}=1;
				
				#Age categories if any; only increment per exon & age category, not per TE
				if ($AGE->{$Rnam}) {
					unless ($checkTE{$tridf}{'age'}) { #easier to load tot hit with these keys for the print_out sub
						($counts->{$cat}{$type}{$fileid}{'age'}{'cat.1'}{'tot'}{'tot'})?
						($counts->{$cat}{$type}{$fileid}{'age'}{'cat.1'}{'tot'}{'tot'}++):($counts->{$cat}{$type}{$fileid}{'age'}{'cat.1'}{'tot'}{'tot'}=1); 
					}
					unless ($checkTE{$tridf}{$AGE->{$Rnam}[4]}) {
						($counts->{$cat}{$type}{$fileid}{'age'}{'cat.1'}{$AGE->{$Rnam}[4]}{'tot'})?
						($counts->{$cat}{$type}{$fileid}{'age'}{'cat.1'}{$AGE->{$Rnam}[4]}{'tot'}++):($counts->{$cat}{$type}{$fileid}{'age'}{'cat.1'}{$AGE->{$Rnam}[4]}{'tot'}=1);
					}
					if (($AGE->{$Rnam}[5]) && (! $checkTE{$tridf}{$AGE->{$Rnam}[5]})) {
						($counts->{$cat}{$type}{$fileid}{'age'}{'cat.2'}{$AGE->{$Rnam}[5]}{'tot'})?
						($counts->{$cat}{$type}{$fileid}{'age'}{'cat.2'}{$AGE->{$Rnam}[5]}{'tot'}++):($counts->{$cat}{$type}{$fileid}{'age'}{'cat.2'}{$AGE->{$Rnam}[5]}{'tot'}=1);
					}
					$checkTE{$tridf}{'age'}=1;
					$checkTE{$tridf}{$AGE->{$Rnam}[4]}=1;
					$checkTE{$tridf}{$AGE->{$Rnam}[5]}=1;
					#tot cat.2 is the same as cat.1 since that's just a key thing.
					$counts->{$cat}{$type}{$fileid}{'age'}{'cat.2'}{'tot'}{'tot'}=$counts->{$cat}{$type}{$fileid}{'age'}{'cat.1'}{'tot'}{'tot'};
				}	
			}	
		}
	}
	close ($fh);
	
	#Get the counts of all features in the set
	$total_exons = load_feat_counts($total_exons,$trinfo,\%chosen_tr,$type,$fileid);
	
	return ($counts,$total_exons);
}

#-----------------------------------------------------------------------------
sub random_tr {
	my ($trinfo,$gene_id,$type) = @_;
	my @trid = keys (%{$trinfo->{$gene_id}{$type}});
	my $r = int(rand(scalar(@trid)));
	return ($trid[$r]);
}

#-----------------------------------------------------------------------------
sub load_feat_counts {
	#fileid contains the run ID
	my ($total_exons,$trinfo,$chosen_tr,$type,$fileid) = @_;
	foreach my $gid (keys %{$trinfo}) {
		my $trid;
		if (! $chosen_tr->{$gid}) {
			$trid = random_tr($trinfo,$gid,$type); 
		} else {
			#this gene was already encountered as overlapping => use that
			$trid =  $chosen_tr->{$gid};
		}	
		#Extract tr infos
		my ($Trstart,$Trend,$Trex) = ($trinfo->{$gid}{$type}{$trid}{'st'},$trinfo->{$gid}{$type}{$trid}{'en'},$trinfo->{$gid}{$type}{$trid}{'nb'});
		
		#Not get the various features
		$total_exons->{$type}{'transcript'}{$fileid}++;
		
		#Now get the counts for the rest
		$total_exons->{$type}{'TSS_polyA'}{$fileid}++;
		$total_exons->{$type}{'TSS'}{$fileid}++;
		$total_exons->{$type}{'polyA'}{$fileid}++;
		$total_exons->{$type}{'exonized'}{$fileid}+=$Trex;
		if ($Trex > 1) {
			#if at elast 2 exons, there will be first and last stuff
			$total_exons->{$type}{'TSS_5SPL'}{$fileid}++;
			$total_exons->{$type}{'3SPL_polyA'}{$fileid}++;
			$total_exons->{$type}{'5SPL'}{$fileid}+=$Trex-1;#minus the last exon
			$total_exons->{$type}{'3SPL'}{$fileid}+=$Trex-1;#minus the first exon
		}
		if ($Trex > 2) {
			#if at least 3, there will be middle exons
			$total_exons->{$type}{'3SPL_exon_5SPL'}{$fileid}+=$Trex-2; #minus the first and last exons
		}
	}	
	return ($total_exons);
}

#-----------------------------------------------------------------------------
sub overlap_category {
	my ($l,$infos,$gid,$type,$trid) = @_;
	my ($Trstart,$Trend,$Trex) = ($infos->{$gid}{$type}{$trid}{'st'},$infos->{$gid}{$type}{$trid}{'en'},$infos->{$gid}{$type}{$trid}{'nb'});
	
	#FYI, structure of $l
    #chr1	4522383	4522590	1111;18.9;4.6;1.0;chr1;4522383;4522590;(190949381);-;B3;SINE/B2;(0);216;1;1923	.	-	chr1	Cufflinks	transcript	4496316	4523815	.	+	.	ID=TCONS_00000002;Parent=XLOC_000001;
	my ($Gst,$Gen) = ($l->[1],$l->[2]);  #TE coordinates
	my ($st,$en) = ($l->[9],$l->[10]);   #exon coordinates
	my $strand = $l->[12];	
	my $cat = "exonized"; #the default
	
	#Check the TSS_polyA with Tr corrdinates first, indep of strand. 
	#Could also be below with overhang both sides, but cleaner to double check with transcript coordinates	
	#TR:       |========|--------|========|--------|========|          #strand does not matter here
    #TE:  [=======================================================]
    #    Gst                                                     Gen
	if ($Gst<$Trstart && $Gen>$Trend) {
		return("TSS_polyA");
	}
		
	#Now the rest; easiest is to set what are the exons
	my $ExType = "MIDDLE";
	$ExType = "FIRST"  if (($strand eq "+" && $st == $Trstart) || ($strand eq "-" && $en == $Trend));
	$ExType = "LAST"   if (($strand eq "+" && $en == $Trend)   || ($strand eq "-" && $st == $Trstart));
	$ExType = "SINGLE" if ($st == $Trstart && $en == $Trend);
	
	if ($Gst < $st) {
		if ($Gen > $en) { # overhang TE start AND end side
			if ($ExType eq "FIRST") {
				#If + strand:
				#          st       en
				#TR:       |========|--------|========|--------|========|    
				#TE:  [==================================]                   #the TE could be overlapping more than this exon, doesn't matter
				return("TSS_5SPL");		
			} elsif ($ExType eq "LAST") {
				#If + strand:
				#                                              st       en
				#TR:       |========|--------|========|--------|========|    
				#TE:                            [==================================] 
				return("3SPL_polyA");
			} else {
				#                            st       en
				#TR:       |========|--------|========|--------|========|
				#TE:  [==================================]        
				return("3SPL_exon_5SPL");          
			}			
		} else {  #overhang TE start side only
			#          st       en
			#TR:       |========|--------|========|--------|========|
			#TE:                     [========]
			($strand eq "+")?($cat = "3SPL"):($cat = "5SPL");
			#          st       en
			#TR:       |========|--------|========|--------|========|
			#TE:  [========]
			$cat = "TSS"   if ($strand eq "+" && ($Trex == 1 || $ExType eq "FIRST"));
			$cat = "polyA" if ($strand eq "-" && ($Trex == 1 || $ExType eq "LAST"));
			return ($cat);
		}
	} elsif ($Gen > $en) { # => overhang only end side
		#                            st       en
		#TR:       |========|--------|========|--------|========|
		#TE:                              [========]
		($strand eq "+")?($cat = "5SPL"):($cat = "3SPL");
		#                                             st       en
		#TR:       |========|--------|========|--------|========|
		#TE:                                               [========]
		$cat = "polyA" if ($strand eq "+" && ($Trex == 1 || $ExType eq "LAST"));
		$cat = "TSS"   if ($strand eq "-" && ($Trex == 1 || $ExType eq "FIRST"));
		return ($cat);
	}
	
	#Gst was > st AND Gen was < en
	#                            st       en
	#TR:       |========|--------|========|--------|========|
	#TE:                           [====]
	#                             Gst  Gen
	return ($cat);
}

#-----------------------------------------------------------------------------
sub print_cat_data {
	my ($no_boot_exons,$boot_exons) = @_;
	#get the no_boot values, avg and sd
	print STDERR "     Get data for each category of overlap\n";
	my $obs = ();
	$obs = get_cat_data($NO_BOOT,0,"na",$obs);
	my $exp = initialize_cat_exp($obs); #0 values for all the ones seen in obs => so that even if not seen in exp, will be there
	$exp = get_cat_data($BOOTS,$NBOOT,$obs,$exp);
	
	my $midval = $NBOOT/2;
	open (my $fh, ">", $STATS.".cat.txt") or confess "ERROR (sub print_stats): can't open to write $STATS.cat.txt $!\n";	
	print $fh "#Script $SCRIPTNAME, v$VERSION\n";
	print $fh "#Aggregated results + stats\n";
	print $fh "#With $MORE repetitions for obs (observed) and $NBOOT bootstraps for exp (expected); sd = standard deviation; nb = number; len = length; avg = average\n";
	print $fh "#Two tests are made (permutation and binomial) to assess how significant the difference between observed and random, so two pvalues are given\n";
	print $fh "#For the two tailed permutation test:\n";
	print $fh "#if rank is < $midval and pvalue is not \"ns\", there are significantly fewer observed values than expected \n";
	print $fh "#if rank is > $midval and pvalue is not \"ns\", there are significantly higher observed values than expected \n";
	print $fh "#The binomial test is done with binom.test from R, two sided\n";
	print $fh "#The category \"gene\" corresponds to the hit of at least one feature of any mature transcript of that gene\n";
	print $fh "#For all categories besides \"gene\", counts are of exons\n";
	print $fh "\n#trancript_type\tcagtegory_id\toverlap_category\tobs_mean\tobs_sd\t%_obs\tobs_tot\tobs_tot_sd\texp_mean\texp_sd\t%_exp\texp_tot\texp_tot_sd\t";
	print $fh "obs_rank_in_exp\t2-tailed_permutation-test_pvalue(obs.vs.exp)\tsignificance\n\n";
	my %o = ('TSS_polyA'=>0,
			 'TSS'=>1,
			 'TSS_5SPL'=>2,	
			 '5SPL'=>3,
			 '3SPL'=>4,
			 '3SPL_exon_5SPL'=>5,
			 'exonized'=>6,
			 '3SPL_polyA'=>7,
			 'polyA'=>8,
			 'transcript'=>9
			 );
	
	foreach my $cat (keys %{$obs}) {
		foreach my $type (keys %{$obs->{$cat}}) {
			my $pval = $exp->{$cat}{$type}{'pval'};
			my $obsper = 0;
			$obsper = $obs->{$cat}{$type}{'avg'}/$no_boot_exons->{$type}{$cat}{'avg'}*100 unless ($no_boot_exons->{$type}{$cat}{'avg'} == 0);		
			my $expper = 0;
			$expper = $exp->{$cat}{$type}{'avg'}/$boot_exons->{$type}{$cat}{'avg'}*100 unless ($boot_exons->{$type}{$cat}{'avg'} == 0);
			my $sign = TEshuffle::get_sign($pval);
			print $fh "$type\t$o{$cat}\t$cat\t$obs->{$cat}{$type}{'avg'}\t$obs->{$cat}{$type}{'sd'}\t$obsper\t$no_boot_exons->{$type}{$cat}{'avg'}\t$no_boot_exons->{$type}{$cat}{'sd'}\t";
			print $fh "$exp->{$cat}{$type}{'avg'}\t$exp->{$cat}{$type}{'sd'}\t$expper\t$boot_exons->{$type}{$cat}{'avg'}\t$boot_exons->{$type}{$cat}{'sd'}\t$exp->{$cat}{$type}{'rank'}\t$pval\t$sign\n";		
		}
	}
	close $fh;
	return 1;
}

#-----------------------------------------------------------------------------
sub print_rep_data {
	my ($no_boot_exons,$boot_exons) = @_;
	print STDERR "     Get data for each repeat, family and class (total and per category)\n";
	my $te_obs = ();
	$te_obs = get_te_data($NO_BOOT,0,"na",$te_obs);	
	my $te_exp = initialize_te_exp($te_obs); #0 values for all the ones seen in obs => so that even if not seen in exp, will be there
	$te_exp = get_te_data($BOOTS,$NBOOT,$te_obs,$te_exp);
	$te_exp = TEshuffle::binomial_test_R($te_exp,"tr");
	my $midval = $NBOOT/2;
	open (my $fh, ">", $STATS.".TE.txt") or confess "ERROR (sub print_stats): can't open to write $STATS.TEs.txt $!\n";	
	print $fh "#Script $SCRIPTNAME, v$VERSION\n";
	print $fh "#Aggregated results + stats\n";
	print $fh "#With $MORE repetitions for obs (observed) and $NBOOT bootstraps for exp (expected)\n";
	print $fh "sd = standard deviation; nb = number; avg = average\n";
	print $fh "#Two tests are made (permutation and binomial) to assess how significant the difference between observed and random, so two pvalues are given\n";
	print $fh "#For the two tailed permutation test:\n";
	print $fh "#if rank is < $midval and pvalue is not \"ns\", there are significantly fewer observed values than expected \n";
	print $fh "#if rank is > $midval and pvalue is not \"ns\", there are significantly higher observed values than expected \n";
	print $fh "#The binomial test is done with binom.test from R, two sided\n";

	print $fh "\n#\t#\tLevel_(tot_means_all)\t#\t#\tCOUNTS\t#\t#\t#\t#\t#\t#\t#\t#\t#\t#\t#\t#\t#\t#\t#\n";
	print $fh "#Type\tCategory\tRclass\tRfam\tRname\t";
	print $fh "obs_nb_of_hits\tobs_nb_sd\t%_obs_nb_(%of_features)\tobs_tot_nb_of_hits\tobs_tot_sd\t";
	print $fh "nb_of_trials(nb_of_TE_in_genome)\t";
	print $fh "exp_nb_of_hits\texp_nb_sd\t%_exp_nb_(%of_features)\texp_tot_nb_of_hits\texp_tot_sd\t";
	print $fh "obs_rank_in_exp\t2-tailed_permutation-test_pvalue(obs.vs.exp)\tsignificance\tbinomal_test_proba\tbinomial_test_95%_confidence_interval\t_binomial_test_pval\n\n";

	foreach my $cat (keys %{$te_exp}) {
		foreach my $type (keys %{$te_exp->{$cat}}) {		
			foreach my $Rclass (keys %{$te_exp->{$cat}{$type}}) { 		
				foreach my $Rfam (keys %{$te_exp->{$cat}{$type}{$Rclass}}) {
					foreach my $Rname (keys %{$te_exp->{$cat}{$type}{$Rclass}{$Rfam}}) {
# 						print STDERR "obs value = $te_obs->{$cat}{$type}{$Rclass}{$Rfam}{$Rname}{'avg'}\n";
						#observed
						my ($te_obsnb,$te_obssd,$te_obsper) = (0,0,0);			
						$te_obsnb = $te_obs->{$cat}{$type}{$Rclass}{$Rfam}{$Rname}{'avg'} if ($te_obs->{$cat}{$type}{$Rclass}{$Rfam}{$Rname}{'avg'});
						$te_obssd = $te_obs->{$cat}{$type}{$Rclass}{$Rfam}{$Rname}{'sd'} if ($te_obs->{$cat}{$type}{$Rclass}{$Rfam}{$Rname}{'sd'});
						$te_obsper = $te_obsnb/$no_boot_exons->{$type}{'avg'}*100 unless ($te_obsnb == 0);
						$te_obs->{$cat}{$type}{'tot'}{'tot'}{'tot'}{'avg'} = 0 unless ($te_obs->{$cat}{$type}{'tot'}{'tot'}{'tot'}{'avg'});						
						#expected
						my $te_expper = 0;
						my $te_expavg = $te_exp->{$cat}{$type}{$Rclass}{$Rfam}{$Rname}{'avg'};	
						$te_expper = $te_expavg/$boot_exons->{$type}{'avg'}*100 unless ($te_expavg == 0);
						#stats
						my $pval = $te_exp->{$cat}{$type}{$Rclass}{$Rfam}{$Rname}{'pval'};		
						$pval = "na" if (($te_expavg == 0) && ($te_obsnb == 0));									
						#Now print stuff
						print $fh "$type\t$cat\t$Rclass\t$Rfam\t$Rname\t";
						print $fh "$te_obsnb\t$te_obssd\t$te_obsper\t$no_boot_exons->{$type}{'avg'}\t$no_boot_exons->{$type}{'sd'}\t"; 
						print $fh "$PARSEDRM->{$Rclass}{$Rfam}{$Rname}\t"; 
						print $fh "$te_expavg\t$te_exp->{$cat}{$type}{$Rclass}{$Rfam}{$Rname}{'sd'}\t$te_expper\t$boot_exons->{$type}{'avg'}\t$boot_exons->{$type}{'sd'}\t";			
						my $sign = TEshuffle::get_sign($pval);				
						print $fh "$te_exp->{$cat}{$type}{$Rclass}{$Rfam}{$Rname}{'rank'}\t$pval\t$sign\t";
						#Binomial
				        $sign = "nd"; #reinitialize
						$sign = TEshuffle::get_sign($te_exp->{$cat}{$type}{$Rclass}{$Rfam}{$Rname}{'binom_pval'});
						print $fh "$te_exp->{$cat}{$type}{$Rclass}{$Rfam}{$Rname}{'binom_prob'}\t$te_exp->{$cat}{$type}{$Rclass}{$Rfam}{$Rname}{'binom_conf'}";
						print $fh "\t$te_exp->{$cat}{$type}{$Rclass}{$Rfam}{$Rname}{'binom_pval'}\t$sign\n";	
					}
				}		
			}
		}
	}
	close $fh;
    return 1;
}

#-----------------------------------------------------------------------------
sub get_exon_data {
	my $tot_ex = shift;	
	my %exons = ();
	foreach my $type (keys %{$tot_ex}) {
		my @data = ();
		foreach my $cat (keys %{$tot_ex->{$type}}) {
			foreach my $round (keys %{$tot_ex->{$type}{$cat}}) {
				push(@data,$tot_ex->{$type}{$cat}{$round});	
			}
			#get average and standard deviation from @data
			($exons{$type}{$cat}{'avg'},$exons{$type}{$cat}{'sd'}) = TEshuffle::get_avg_and_sd(\@data);
		}	
	}
	return(\%exons);
}

#-----------------------------------------------------------------------------
sub initialize_cat_exp {
	my $obs = shift;
	my $exp = ();
	#obs:
	#$cat_data->{$cat}{$type}{'pval'}
	foreach my $k1 (keys %{$obs}) {
		foreach my $k2 (keys %{$obs->{$k1}}) {
			foreach my $k3 (keys %{$obs->{$k1}{$k2}}) {
				$exp->{$k1}{$k2}{$k3}{'avg'}=0;
				$exp->{$k1}{$k2}{$k3}{'sd'}=0;
				$exp->{$k1}{$k2}{$k3}{'rank'}="na";
				$exp->{$k1}{$k2}{$k3}{'pval'}="na";
			}
		}
	}
	return $exp;
}

#-----------------------------------------------------------------------------
sub get_cat_data {
	my ($all_data,$n,$obs,$cat_data) = @_;
	foreach my $cat (keys %{$all_data}) {
		foreach my $type (keys %{$all_data->{$cat}}) {
			my @data = ();
			foreach my $round (keys %{$all_data->{$cat}{$type}}) {
				my $hits = $all_data->{$cat}{$type}{$round}{'tot'}{'tot'}{'tot'}{'nr'};
				$hits = 0 unless ($hits);
				push(@data,$hits);	
			}
			#get average and standard deviation from @data
			($cat_data->{$cat}{$type}{'avg'},$cat_data->{$cat}{$type}{'sd'}) = TEshuffle::get_avg_and_sd(\@data);
				
			#Now get he rank of the observed value in the list of expected => get a p value
			unless ($n == 0) {
				($obs->{$cat}{$type}{'avg'},$obs->{$cat}{$type}{'sd'}) = (0,"na") unless ($obs->{$cat}{$type}{'avg'});
			
				my $rank = 1; #pvalue can't be 0, so I have to start there - that does mean there will be a rank nboot+1
				my @data = sort {$a <=> $b} @data;
				TEshuffle::print_distrib(\@data,$DISTRIB,$type,$cat) if ($KEEP ne "none");
				EXP: foreach my $exp (@data) {
					last EXP if ($exp > $obs->{$cat}{$type}{'avg'});
					$rank++ if ($exp < $obs->{$cat}{$type}{'avg'});
				}
				$cat_data->{$cat}{$type}{'rank'}=$rank;
				if ($rank <= $NBOOT/2) {
					$cat_data->{$cat}{$type}{'pval'}=$rank/$NBOOT*2; #*2 because 2 tailed
				} else {
					$cat_data->{$cat}{$type}{'pval'}=($NBOOT+2-$rank)/$NBOOT*2;  #+2 so it is symetrical (around nboot+1)
				}			
			}
		}
	}
	return($cat_data);
}

#-----------------------------------------------------------------------------
sub initialize_te_exp {
	my $obs = shift;
	my $exp = ();
	#obs:
	#$te_data->{$cat}{$type}{$key1}{$key2}{$key3}{'rank'};
	foreach my $c (keys %{$obs}) {
		foreach my $t (keys %{$obs->{$c}}) {
			foreach my $k1 (keys %{$obs->{$c}{$t}}) {
				foreach my $k2 (keys %{$obs->{$c}{$t}{$k1}}) {
					foreach my $k3 (keys %{$obs->{$c}{$t}{$k1}{$k2}}) {
						$exp->{$c}{$t}{$k1}{$k2}{$k3}{'avg'}=0;
						$exp->{$c}{$t}{$k1}{$k2}{$k3}{'sd'}=0;
						$exp->{$c}{$t}{$k1}{$k2}{$k3}{'rank'}="na";
						$exp->{$c}{$t}{$k1}{$k2}{$k3}{'pval'}="na";
					}
				}
			}
		}	
	}		
	return $exp;
}

#-----------------------------------------------------------------------------
sub get_te_data {
	my ($counts,$nboot,$te_obs,$te_data) = @_;
	# $counts->{$cat}{$type}{$fileid}{'tot'}{'tot'}{'tot'}{'tot'}
	# $counts->{$cat}{$type}{$fileid}{$Rcla}{'tot'}{'tot'}{'tot'}
	# $counts->{$cat}{$type}{$fileid}{$Rcla}{$Rfam}{'tot'}{'tot'}
	# $counts->{$cat}{$type}{$fileid}{$Rcla}{$Rfam}{$Rnam}{'tot'}

	#agregate data [Since I can easily have the categories details here, put them so I can compare]
	my ($nb_c,$nb_f,$nb_r,$nb_a1,$nb_a2) = ();
	foreach my $cat (keys %{$counts}) {
		foreach my $type (keys %{$counts->{$cat}}) {
			foreach my $round (keys %{$counts->{$cat}{$type}}) {
				push(@{$nb_c->{$cat}{$type}{'tot'}{'tot'}{'tot'}},$counts->{$cat}{$type}{$round}{'tot'}{'tot'}{'tot'}{'tot'});	
				foreach my $Rclass (keys %{$counts->{$cat}{$type}{$round}}) {
					push(@{$nb_c->{$cat}{$type}{$Rclass}{'tot'}{'tot'}},$counts->{$cat}{$type}{$round}{$Rclass}{'tot'}{'tot'}{'tot'}) if ($Rclass ne "age");		
					foreach my $Rfam (keys %{$counts->{$cat}{$type}{$round}{$Rclass}}) {
						push(@{$nb_f->{$cat}{$type}{$Rclass}{$Rfam}{'tot'}},$counts->{$cat}{$type}{$round}{$Rclass}{$Rfam}{'tot'}{'tot'}) if ($Rclass ne "age");		
						foreach my $Rname (keys %{$counts->{$cat}{$type}{$round}{$Rclass}{$Rfam}}) {
							my $count = $counts->{$cat}{$type}{$round}{$Rclass}{$Rfam}{$Rname}{'tot'};
							$count = 0 if (! $count); #Might happen, if NO_BOOT has no overlap
							push(@{$nb_r->{$cat}{$type}{$Rclass}{$Rfam}{$Rname}},$count);
							push(@{$nb_a1->{$cat}{$type}{$Rclass}{$Rfam}{$Rname}},$count) if ($Rclass eq "age" && $Rfam eq "cat.1");								
							push(@{$nb_a2->{$cat}{$type}{$Rclass}{$Rfam}{$Rname}},$count) if ($Rclass eq "age" && $Rfam eq "cat.2");							
						}
					}
				}
			}
		}	
	}
	
	#get avg, sd and p values now => load in new hash, that does not have the fileID
	foreach my $cat (keys %{$counts}) {
		foreach my $type (keys %{$counts->{$cat}}) {
			foreach my $round (keys %{$counts->{$cat}{$type}}) {
				foreach my $Rclass (keys %{$counts->{$cat}{$type}{$round}}) {
					$te_data = get_te_data_details($nboot,$cat,$type,$Rclass,"tot","tot",$nb_c->{$cat}{$type}{$Rclass}{'tot'}{'tot'},$te_data,$te_obs) 
						if ($Rclass ne "age");	
					foreach my $Rfam (keys %{$counts->{$cat}{$type}{$round}{$Rclass}}) {
						$te_data = get_te_data_details($nboot,$cat,$type,$Rclass,$Rfam,"tot",$nb_f->{$cat}{$type}{$Rclass}{$Rfam}{'tot'},$te_data,$te_obs) 
							if ($Rclass ne "age");	
						foreach my $Rname (keys %{$counts->{$cat}{$type}{$round}{$Rclass}{$Rfam}}) {
							my $nb = $nb_r->{$cat}{$type}{$Rclass}{$Rfam}{$Rname};							
							$te_data = get_te_data_details($nboot,$cat,$type,$Rclass,$Rfam,$Rname,$nb,$te_data,$te_obs);
							if ($Rclass eq "age" && $Rfam eq "cat.1") {
								$nb = $nb_a1->{$cat}{$type}{$Rclass}{$Rfam}{$Rname};
								$te_data = get_te_data_details($nboot,$cat,$type,$Rclass,$Rfam,$Rname,$nb,$te_data,$te_obs);
							}
							if ($Rclass eq "age" && $Rfam eq "cat.2") {
								$nb = $nb_a2->{$cat}{$type}{$Rclass}{$Rfam}{$Rname};
								$te_data = get_te_data_details($nboot,$cat,$type,$Rclass,$Rfam,$Rname,$nb,$te_data,$te_obs);
							}
						}
					}
				}		
			}
		}
	}		
	$counts = (); #empty this
	return($te_data);
}

#-----------------------------------------------------------------------------	
sub get_te_data_details {
	my ($nboot,$cat,$type,$key1,$key2,$key3,$agg_data,$te_data,$te_obs) = @_;	
	#get average and sd of the expected	
	($te_data->{$cat}{$type}{$key1}{$key2}{$key3}{'avg'},$te_data->{$cat}{$type}{$key1}{$key2}{$key3}{'sd'}) = TEshuffle::get_avg_and_sd($agg_data);
	
#	print STDERR "{$cat}{$type}{$key1}{$key2}{$key3}\n   AVG  = $te_data->{$cat}{$type}{$key1}{$key2}{$key3}{'avg'}; SD = $te_data->{$cat}{$type}{$key1}{$key2}{$key3}{'sd'}\n";
		
	if ($nboot > 0) {	
		my $observed = $te_obs->{$cat}{$type}{$key1}{$key2}{$key3}{'avg'};		
#		print STDERR "FYI: no observed value for {$cat}{'no_boot'}{$key1}{$key2}{$key3}{'tot'}\n" unless ($observed);
		$observed = 0 unless ($observed);				
		#Get the rank of the observed value in the list of expected + pvalue for the permutation test
		my $rank = 1; #pvalue can't be 0
		my @data = sort {$a <=> $b} @{$agg_data} if ($agg_data->[1]);
		EXP: foreach my $exp (@data) {
			last EXP if ($exp > $observed);
			$rank++;
		}
		$te_data->{$cat}{$type}{$key1}{$key2}{$key3}{'rank'}=$rank;
		if ($rank <= $NBOOT/2) {
			$te_data->{$cat}{$type}{$key1}{$key2}{$key3}{'pval'}=$rank/$NBOOT*2; #*2 because 2 tailed
		} else {
			$te_data->{$cat}{$type}{$key1}{$key2}{$key3}{'pval'}=($NBOOT+2-$rank)/$NBOOT*2;  #+2 so it is symetrical (around nboot+1)
		}
				
		#Binomial test
		#get all the values needed for binomial test in R => do them all at once
		#N, number of trials - here, total number of TEs in the genome
		my $n = $PARSEDRM->{$key1}{$key2}{$key3} if ($PARSEDRM->{$key1}{$key2}{$key3});
		$n = 0 unless ($n);
		my $p = 0;		
		$p=$te_data->{$cat}{$type}{$key1}{$key2}{$key3}{'avg'}/$n unless ($n == 0); #should not happen, but could
		$te_data->{$cat}{$type}{$key1}{$key2}{$key3}{'p'} = $p;
		$te_data->{$cat}{$type}{$key1}{$key2}{$key3}{'n'} = $n;
		#x cannot be non integer - round to closest int
		my $rounded = sprintf "%.0f", $observed;
		$te_data->{$cat}{$type}{$key1}{$key2}{$key3}{'x'} = $rounded;
	}
	return($te_data);
}


