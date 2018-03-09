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
use Statistics::R; #required to get the Binomial Test p-values

use vars qw($BIN);
use Cwd 'abs_path';
BEGIN { 	
	$BIN = abs_path($0);
	$BIN =~ s/(.*)\/.*$/$1/;
	unshift(@INC, "$BIN/Lib");
}
use GAL::Annotation; #if issues, there is an alternative subroutine not using this module, see usage
use TEshuffle;
use Data::Dumper;

#-----------------------------------------------------------------------------
#------------------------------- DESCRIPTION ---------------------------------
#-----------------------------------------------------------------------------
#flush buffer
$| = 1;

my $version = "6.3";
my $scriptname = "TE-analysis_Shuffle.pl";
my $changelog = "
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
#      Also added the option to use -r file in -s rm as well (to check ends) and chift the TE if needed.
#	- v6.3 = Mar 08 2018
#           Change for -s tss: if the TE ends up out of the scaffold/chr, put on the other side and if still out, shift it to be inside
#           Also, skip if not in the annotation file
#           Make -r mandatory for all
\n";

my $usage = "
Synopsis (v$version):

    perl $scriptname -l lncRNA.gff -p prot.gff [-o <nt>] [-m <nb>] -q features_to_shuffle [-n <nb>] 
            -s shuffling_type -r <genome.sizes> [-b] 
            [-a <annotations>] [-f]
            [-e <genome.gaps>] [-d] [-i <include.range>] [-x] 
            [-w <bedtools_path>] [-u <no_low>] [-t <type,name>] [-c] [-g <TE.age.tab>] [-v] [-h]

    /!\\ REQUIRES - Bedtools, v2.25+
	              - GAL::Annotation version later than Jan 2016 [update of is_coding]
	                see https://github.com/The-Sequence-Ontology/GAL
	                If issues with it, set the --just option                
    /!\\ Previous outputs, if any, will be moved as *.previous [which only saves results once]

    Typically, for the 3 types, the mandatory arguments are:
    perl $scriptname -l lncRNA.gff -p prot.gff -q rm.out -r genome.sizes -s rm 
    perl $scriptname -l lncRNA.gff -p prot.gff -q rm.out -r genome.sizes -s tss -a annotations.gtf
    perl $scriptname -l lncRNA.gff -p prot.gff -q rm.out -r genome.sizes -s bed -r genome.range -e genome.gaps

	Note that -r is advised for -s rm (but won't affect -s tss)

  CITATIONS:
    - Cite Kapusta et al. (2013) PLoS Genetics (DOI: 10.1371/journal.pgen.1003470)
      (http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1003470)
      but should also include the version and link GitHub (https://github.com/4ureliek/TEanalysis)
    - for BEDtools, Quinlan AR and Hall IM (2010) Bioinformatics (DOI: 10.1093/bioinformatics/btq033)

  SOME HISTORY:
    Originally writen by Zev Kronenberg (https://github.com/zeeev) to generate 
    data (observed vs expected) shown in Figure 5 of Kapusta et al. 2013 PLoS Genetics
    But highly modified since then by A. Kapusta
    Thanks to Edward Chuong, PhD, for the suggestion of the shuffling '-t tss', 
    and thanks to E. Chuong and Cedric Feschotte for helpful discussions during the development of this script

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
    or the input file should be corrected when possible to limit that issue [not implemented in this script for now]
    
    Note that one exon may have several types of overlaps (e.g. \"SPL\" and \"exonized\"),
    but each exon is counted only one time for each category (important for \"exonized\").
    Similarly for TEs, each hit is counted unless it's the same repeat name / family / class (depending on the level)
   
    If you need to generate the <genome.gaps> file but you would also like to add more files to the -e option, 
    just do a first run with no bootstraps (in this example the genome.range is also being generated):
    perl ~/bin/$scriptname -l input.gtf -q genome.out -r genome.fa -b -e genome.fa -d -n 0

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
                              When that happens, the script will die and switch to alternative way of loading
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
                                 perl ~/bin/$scriptname -l lncRNA.gff -p prot.gff -q genome.out -s rm -r genome.fa -b -e genome.fa -d -n 0                                
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

#-----------------------------------------------------------------------------
#------------------------------ LOAD AND CHECK -------------------------------
#-----------------------------------------------------------------------------
my ($shuffle,$stype,$tssfile,$just,$full,$exclude,$dogaps,$build,$dobuild,$f_regexp,$allow,$nooverlaps,$v,$help);
my ($prot,$linc) = ("n","n");
my $inters = 10;
my $more = 0;
my $nboot = 10;
my $incl = "na";
my $nonTE = "no_low";
my $filter = "na";
my $TEage = "na";
my $bedtools = "";
my $catout = "y"; #removed from options, not really relevant to ask for choice
my $opt_success = GetOptions(
			 	  'prot=s'		=> \$prot,
			 	  'lnc=s'		=> \$linc,
			 	  'more=s'      => \$more,
			 	  'query=s'     => \$shuffle,
			 	  'shuffle=s'   => \$stype,
			 	  'annot=s'     => \$tssfile,
			 	  'just'        => \$just,				 	  
			 	  'overlap=s'   => \$inters,
			 	  'nboot=s'     => \$nboot,
			 	  'full'        => \$full,
			 	  'range=s'     => \$build,
			 	  'build'       => \$dobuild,
			 	  'excl=s'		=> \$exclude,
			 	  'dogaps'      => \$dogaps,
			 	  'incl=s'		=> \$incl,
			 	  'x'		    => \$nooverlaps,
			 	  'u=s'		    => \$nonTE,
			 	  'te=s'		=> \$filter,
			 	  'contain'     => \$f_regexp,
			 	  'group=s'     => \$TEage,
			 	  'where=s'     => \$bedtools,
			 	  'version'     => \$v,
			 	  'help'		=> \$help,);
			 	  
#Check options, if files exist, etc
die "\n --- $scriptname version $version\n\n" if $v;
die $usage if ($help);
die "\n SOME MANDATORY ARGUMENTS MISSING, CHECK USAGE:\n$usage" if ((! $linc && ! $prot) || ! $shuffle || ! $stype || ! $build);
die "\n One of -l or -p needs to be provided\n\n" if (($prot eq "n") && ($linc eq "n" ));
die "\n -t is required\n\n" if ( ! $stype);
die "\n -q is required\n\n" if ( ! $shuffle);
die "\n -r $build does not exist?\n\n" if (! -e $build);

die "\n -p $prot does not exist?\n\n"  if (($prot ne "n") && (! -e $prot));
die "\n -p $prot is not a gff file?\n\n" unless (($prot eq "n") || ($prot =~ /\.gff$/) || ($prot =~ /\.gff3$/));
die "\n -l $linc does not exist?\n\n"  if (($linc ne "n") && (! -e $linc));
die "\n -l $linc is not a gff file?\n\n" unless (($linc eq "n") || ($linc =~ /\.gff$/) || ($linc =~ /\.gff3$/));

die "\n -q $shuffle is not in a proper format? (not .out, .bed, .gff or .gff3)\n\n" unless (($shuffle =~ /\.out$/) || ($shuffle =~ /\.bed$/) || ($shuffle =~ /\.gff$/) || ($shuffle =~ /\.gff3$/));
die "\n -q $shuffle does not exist?\n\n" if (! -e $shuffle);
die "\n -s $stype should be one of the following: bed, rm or tss\n\n" if (($stype ne "bed") && ($stype ne "rm") && ($stype ne "tss"));
#deal with conditional mandatory stuff
die "\n -s tss was set, but -a is missing?\n\n" if (($stype eq "tss") && (! $tssfile));
die "\n -a $tssfile does not exist?\n\n" if (($tssfile) && (! -e $tssfile));
if ($stype eq "bed") {
	die "\n -s bed was set, but -e is missing?\n\n" if (! $exclude);
	die "\n -e $exclude does not exist?\n\n" if (($exclude !~ /,/) && (! -e $exclude)); #if several files, can't check existence here
	die "\n -i $incl does not exist?\n\n" if (($incl ne "na") && ($incl !~ /,/) && (! -e $incl)); #if several files, can't check existence here
}
#Now the rest
die "\n -n $nboot but should be an integer\n\n" if ($nboot !~ /\d+/);
die "\n -i $inters but should be an integer\n\n" if ($inters !~ /\d+/);
die "\n -w $bedtools does not exist?\n\n" if (($bedtools ne "") && (! -e $bedtools));
die "\n -t requires 2 values separated by a coma (-t <name,filter>; use -h to see the usage)\n\n" if (($filter ne "na") && ($filter !~ /,/));
die "\n -g $TEage does not exist?\n\n" if (($TEage ne "na") && (! -e $TEage));

($full)?($full = "y"):($full = "n");
($dogaps)?($dogaps = "y"):($dogaps = "n");
($dobuild)?($dobuild = "y"):($dobuild = "n");
($f_regexp)?($f_regexp = "y"):($f_regexp="n");
$bedtools = $bedtools."/" if (($bedtools ne "") && (substr($bedtools,-1,1) ne "/")); #put the / at the end of path if not there
($nooverlaps)?($nooverlaps = "-noOverlapping"):($nooverlaps = "");
$more = 1 if ($more == 0); #1 rep if set to 0, same thing here


#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------
#Prep steps
print STDERR "\n --- $scriptname v$version started, with:\n";
print STDERR "     input lncRNA file = $linc\n" if ($linc ne "n");
print STDERR "     input mRNA file = $prot\n" if ($prot ne "n");
print STDERR "     features to shuffle = $shuffle\n";
print STDERR "     shuffling type = $stype\n";

#Outputs
print STDERR " --- prepping output directories and files\n";
my ($input);
($linc)?($input = $linc):($input = $prot);
my $dir = $input.".shuffle-".$stype.".".$nboot;
print STDERR "     output directory = $dir\n";
my ($stats,$outl,$outlb,$temp_l,$outp,$outpb,$temp_p,$temp) = TEshuffle::prep_out("tr",$dir,$nboot,$filter,$input,$stype,$nonTE,$linc,$prot,$shuffle);

#Chosomosome sizes / Genome range
print STDERR " --- loading build (genome range)\n";
my ($okseq,$build_file) = TEshuffle::load_build($build,$dobuild);

#prep steps if shuffling type is bed
my $excl;
if ($stype eq "bed") {
	#Files to exclude for shuffling
	print STDERR " --- getting ranges to exclude in the shuffling of features from $exclude\n";
	my @exclude = ();
	if ($exclude =~ /,/) {
		($dogaps eq "y")?(print STDERR "     several files provided, -d chosen, genome file (fasta) should be the first one\n"):
						 (print STDERR "     several files provided, assembly gaps should be the first one\n");
		@exclude = split(",",$exclude) if ($exclude =~ /,/);
	} else {
		$exclude[0] = $exclude;
	}
	$exclude[0] = TEshuffle::load_gap($exclude[0],$dogaps);
	print STDERR "     concatenating files for -e\n" if ($exclude =~ /,/);
	($exclude =~ /,/)?($excl = TEshuffle::concat_beds(\@exclude)):($excl = $exclude[0]);

	#If relevant, files to include for shuffling
	if (($incl ne "na") && ($incl =~ /,/)) {
		print STDERR " --- concatenating $incl files to one file\n";
		my @include = split(",",$incl);
		$incl = TEshuffle::concat_beds(\@include);
	}
}

#Load TEage if any
print STDERR " --- Loading TE ages from $TEage\n" unless ($TEage eq "na");
my $age = ();
$age = TEshuffle::load_TEage($TEage,$v) unless ($TEage eq "na");

#Now features to shuffle (need to be after in case there was $okseq loaded)
print STDERR " --- checking file in -s, print in .bed if not a .bed or gff file\n";
print STDERR "     filtering TEs based on filter ($filter) and non TE behavior ($nonTE)\n" unless ($filter eq "na");
print STDERR "     + getting genomic counts for each repeat\n";
print STDERR "     + load all TE positions in a hash (since $stype is set to rm)\n" if ($stype eq "rm") ;
my ($toshuff_file,$parsedRM,$rm,$rm_c) = TEshuffle::RMtobed($shuffle,$okseq,$filter,$f_regexp,$nonTE,$age,"y",$stype); #Note: $rm and $rm_c are empty unless $stype eq rm

#prep steps if shuffling type is tss
my ($tssbed,$closest,$alltss);
if ($stype eq "tss")  {
	#sort TEs
	my $bedsort = $bedtools."bedtools sort";
	my $sorted = $toshuff_file;
	$sorted =~ s/\.bed$/\.sorted\.bed/;
	print STDERR " --- sorting features of $toshuff_file\n" unless (-e $sorted);	
	print STDERR "     $bedsort -i $toshuff_file > $sorted\n" unless (-e $sorted);	
	`$bedsort -i $toshuff_file > $sorted` unless (-e $sorted);
	$toshuff_file = $sorted;
	print STDERR " --- loading the tss from $tssfile\n";
	#print the tss in a bed file => use bedtools closest
	($tssbed,$alltss) = TEshuffle::load_and_print_tss($tssfile);	
	print STDERR " --- sorting features in the tss file\n" unless (-e "$tssbed.bed");
	print STDERR "     $bedsort -i $tssbed > $tssbed.bed\n" unless (-e "$tssbed.bed");
	`$bedsort -i $tssbed > $tssbed.bed` unless (-e "$tssbed.bed");
	$tssbed = $tssbed.".bed";
	#get the closest tss if relevant
	my $tssclosest = $toshuff_file.".closest.tss";
	my $closestBed = $bedtools."closestBed";
	print STDERR " --- getting closest tss for each feature in $toshuff_file, with the command line below\n";
	print STDERR "     $closestBed -a $toshuff_file -b $tssbed -D b -first > $tssclosest\n"; #I want only one entry per TE, therefore -t first
	`$closestBed -a $toshuff_file -b $tssbed -D b -t first > $tssclosest`;
	print STDERR " --- loading distance to TSS\n";
	$closest = TEshuffle::load_closest_tss($tssclosest);
}

#Load the gff file(s)
print STDERR " --- Load gene IDs / transcript IDs for:\n";
my $whichgene = ();
my $l_tr = ();
my $p_tr = ();
my $flagGAL = 0;
print STDERR "  -> $linc\n" unless ($linc eq "n");
($l_tr,$whichgene,$flagGAL) = read_gff($linc,$okseq,0,$whichgene,$stype) if (($linc ne "n") && (! $just));
($l_tr,$whichgene) = load_gene_tr($linc,$okseq,0,$whichgene,$stype) if (($linc ne "n") && (($just) || ($flagGAL == 1))); #if GAL::Annotation is a problem
print STDERR "  -> $prot\n" unless ($prot eq "n");
($p_tr,$whichgene,$flagGAL) = read_gff($prot,$okseq,1,$whichgene,$stype) if (($prot ne "n") && (! $just));
($p_tr,$whichgene) = load_gene_tr($prot,$okseq,1,$whichgene,$stype) if (($prot ne "n") && (($just) || ($flagGAL == 1))); #if GAL::Annotation is a problem

#Join -p and/or -l files
my $intersectBed = $bedtools."intersectBed";
print STDERR " --- Intersect with command lines:\n";
print STDERR "      $intersectBed -a $toshuff_file -b $linc -wo > $temp_l/no_boot.joined\n" unless ($linc eq "n");
system "$intersectBed -a $toshuff_file -b $linc -wo > $temp_l/no_boot.joined" unless ($linc eq "n");
print STDERR "      $intersectBed -a $toshuff_file -b $prot -wo > $temp_p/no_boot.joined\n" unless ($prot eq "n");
system "$intersectBed -a $toshuff_file -b $prot -wo > $temp_p/no_boot.joined" unless ($prot eq "n");

#Process the joined files with -m X repeats
print STDERR " --- Check intersection(s) with features in $toshuff_file (observed)\n";
print STDERR "     (if -m set, there will be several rounds of random transcript selection)\n";
my $no_boot = ();
my $no_boot_tot_exons = ();
for(my $j = 1; $j <= $more; $j++) {
	print STDERR "     ..$j rounds done\n" if (($j == 10) || ($j == 100) || ($j == 1000) || (($j > 1000) && (substr($j/1000,-1,1) == 0)));	
	($no_boot,$no_boot_tot_exons) = 
		check_for_featured_overlap("$temp_l/no_boot.joined",$l_tr,"no_boot.".$j,'transcript',$outl,$inters,$no_boot,$no_boot_tot_exons,$whichgene,$full) 
		unless ($linc eq "n");
	($no_boot,$no_boot_tot_exons) = 
		check_for_featured_overlap("$temp_p/no_boot.joined",$p_tr,"no_boot.".$j,'mRNA',$outp,$inters,$no_boot,$no_boot_tot_exons,$whichgene,$full) 
		unless ($prot eq "n");
	`cat $outl >> $catout.no-boot.txt` if (($catout) && (-e $outl));
	`cat $outp >> $catout.no-boot.txt` if (($catout) && (-e $outp));
}

#Now bootstrap runs
print STDERR " --- Run $nboot bootstraps now (to get significance of the overlaps)\n";
my $boots = ();
my $boots_tot_exons = ();
if ($nboot > 0) {
	foreach (my $i = 1; $i <= $nboot; $i++) {
		print STDERR "     ..$i bootstraps done\n" if (($i == 10) || ($i == 100) || ($i == 1000) || (($i > 1000) && (substr($i/1000,-1,1) == 0)));	
		my $shuffled;
		$shuffled = TEshuffle::shuffle_tss($toshuff_file,$temp,$i,$alltss,$closest,$okseq) if ($stype eq "tss");
		$shuffled = TEshuffle::shuffle_rm($toshuff_file,$temp,$i,$rm,$rm_c,$okseq) if ($stype eq "rm");	
		$shuffled = TEshuffle::shuffle_bed($toshuff_file,$temp,$i,$excl,$incl,$build_file,$bedtools,$nooverlaps) if ($stype eq "bed");
 		system "      $intersectBed -a $shuffled -b $linc -wo > $temp_l/boot.$i.joined" unless ($linc eq "n");
 		system "      $intersectBed -a $shuffled -b $prot -wo > $temp_p/boot.$i.joined" unless ($prot eq "n");	
		($boots,$boots_tot_exons) = 
			check_for_featured_overlap("$temp_l/boot.$i.joined",$l_tr,"boot.".$i,'transcript',$outlb,$inters,$boots,$boots_tot_exons,$whichgene,$full) 
			unless ($linc eq "n");
		($boots,$boots_tot_exons) = 
			check_for_featured_overlap("$temp_p/boot.$i.joined",$p_tr,"boot.".$i,'mRNA',$outpb,$inters,$boots,$boots_tot_exons,$whichgene,$full) 
			unless ($prot eq "n");		
		`cat $outlb >> $catout.boot.txt` if (($catout) && (-e $outlb));
		`cat $outpb >> $catout.boot.txt` if (($catout) && (-e $outpb));
		`rm -Rf $shuffled`; #these files are now not needed anymore, all is stored
		`rm -Rf $temp_l/boot.$i.joined` unless ($linc eq "n");
		`rm -Rf $temp_p/boot.$i.joined` unless ($prot eq "n");
	}
}
`rm -Rf $temp_l` unless ($linc eq "n");
`rm -Rf $temp_p` unless ($prot eq "n");

#Stats now
print STDERR " --- Get and print stats\n" if ($nboot > 0);
print_stats($stats,$no_boot,$more,$no_boot_tot_exons,$boots,$nboot,$boots_tot_exons,$parsedRM,$age,$full,$scriptname,$version) if ($nboot > 0);

#end
print STDERR " --- $scriptname done\n";
print STDERR "     Stats for categeories printed in: $stats.cat.txt\n" if ($nboot > 0);
print STDERR "     Stats for TEs printed in: $stats.TE.txt\n" if (($nboot > 0) && ($full eq "y"));
print STDERR "\n";
exit;


#-----------------------------------------------------------------------------
#-------------------------------- SUBROUTINES --------------------------------
#-----------------------------------------------------------------------------
sub load_gene_tr {
	#Not as solid as using GAL::Annotation, but it is an alternative
	my ($file,$okseq,$coding,$whichgene,$stype) = @_;
    print STDERR "     Loading transcripts / genes relationships without using GAL::Annotation\n";
	my %tr = ();
	my ($gid,$trid,$type);
	open(my $fh, "<$file") or confess "\n   ERROR (sub load_gene_tr): could not open to read $file!\n";
	LINE: while(<$fh>) {
		chomp(my $l = $_);
		next LINE if (substr($l,0,1) eq "#");
		my @l = split('\s+',$l);
		next LINE if (($stype eq "bed") && (! $okseq->{$l[0]})); #if not in build of stuff OK to shuffle on, remove here as well; only relevant for -s bed though
		my $id = $l[8];
		$id = $1 if $id =~ /^ID=(.+?);/;
		if ($l[2] eq "gene") {
			$gid = $id;
		} elsif ($l[2] eq "transcript") {
			$trid = $id;
			$type = "transcript";
			$type = "mRNA" if (($l =~ /protein_coding/) || ($coding == 1)); #if coding, it should be in the "gene_type" or the "transcript_type", but if not assume when coding is set to 1, it's coding
			$tr{$gid}{$type}{$trid}{'st'} = $l[3];
			$tr{$gid}{$type}{$trid}{'en'} = $l[4];
			$whichgene->{$trid}=$gid;
		} else {	
			$tr{$gid}{$type}{$trid}{'nb'}++ if ($l[2] eq "exon"); #count exons; includes UTRs for pc genes
# 			$tr{$gid}{$type}{$trid}{'nb'}++ if (($type eq "mRNA") && ($l[2] eq "CDS")); #count number of coding exons only
#			$tr{$gid}{$type}{$trid}{'nb'}++ if (($type eq "transcript") && ($l[2] eq "exon")); #count number of exons
		}
	}
	close ($fh);
	return (\%tr,$whichgene); #looping through keys will get transcripts => put in an array for each gene later
}

#-----------------------------------------------------------------------------
sub read_gff {
	my ($gff3_file,$okseq,$coding,$whichgene,$stype) = @_;
    my %trinfo;
	my $gene_count=0;
	my $flagGAL = 0;
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
		next GENE if (($stype eq "bed") && (! $okseq->{$seqid})); #if not in build of stuff OK to shuffle on, remove here as well; only relevant for -s bed though
		my @tr = $gene->transcripts;
		TRANSCRIPT: foreach my $tr (@tr) {
			my $tr_id = $tr->feature_id;
			my $tr_strand = $tr->strand;
			if($tr_strand !~ /\+|-/){
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
			$whichgene->{$tr_id}=$gene_id;
		}
		$gene_count++;		
	}
    print STDERR "        total genes loaded (type=$type): $gene_count\n";	
    $flagGAL = 1 if ($gene_count == 0);
    print STDERR "        WARN: gene count = 0, transcript info will (try to) be loaded without GAL::Annotation\n" if ($gene_count == 0);
	return (\%trinfo,$whichgene,$flagGAL);
}

#-----------------------------------------------------------------------------
sub check_for_featured_overlap {
	my ($file,$infos,$fileid,$type,$out,$inters,$counts,$total_exons,$whichgene,$full) = @_;
	my %chosen_tr = ();
	my %checkgenes = ();
	my %check = ();
	my %checkTE = ();

	#initialize exon count for this run
	$total_exons->{$type}{$fileid}{'tot'}=0;
	$total_exons->{$type}{$fileid}{'hit'}=0;	
	
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
			my $trid = $l[14];
			$trid = $1 if $trid =~ /Parent=(.+?);/;
			next LINE unless (defined $whichgene->{$trid}); #checked for non coding when coding are looked at
			
			#get a random tr for this gene, but only the first time this gene is met, and keep which tr is chosen			
			my $gid = $whichgene->{$trid};			
			$chosen_tr{$gid} = random_tr($infos,$gid,$type) unless (defined $chosen_tr{$gid}); #%infos contain infos about the transcript => start, end, number of exons in it
			my $chosen = $chosen_tr{$gid};
			$total_exons->{$type}{$fileid}{'tot'}+=$infos->{$gid}{$type}{$chosen}{'nb'} 
				unless (defined $checkgenes{$gid}); #increment with exon numbers of the transcript chosen for this gene
			$checkgenes{$gid}=1; #store that exons of a random transcript for this gene have been counted			
			next LINE if ($trid ne $chosen); #skip if not randomly chosen transcript that is hit	
			
			my $ilen = $l[-1]; #last value of the line is intersection length
			next LINE unless ($ilen >= $inters);

			#now check what category of overlap this exon is;
			my $cat = overlap_category(\@l,$infos,$gid,$type,$trid);
		
			#now increment in the data structure		
			#since only one transcript per gene, there should be no worry here about unique counts, 1 exon can only be counted one time in a category; 
			#however unique exon hits count need a check, and there could be TE overlaps fucking things up, so better safe than sorry	
			unless (defined $check{$l[14]}{$cat}) { 
				($counts->{$cat}{$type}{$fileid}{'tot'}{'tot'}{'tot'}{'nr'})?
				($counts->{$cat}{$type}{$fileid}{'tot'}{'tot'}{'tot'}{'nr'}++):($counts->{$cat}{$type}{$fileid}{'tot'}{'tot'}{'tot'}{'nr'}=1);
			}
			$check{$l[14]}{$cat}=1;
			$total_exons->{$type}{$fileid}{'hit'}++ unless (defined $check{$l[14]}{'hit'});
			$check{$l[14]}{'hit'}=1;
			unless (defined $check{$gid}{'hit'}) { #counting each tr hit only one time => equivalent to a number of genes, not transcripts
#			unless (defined $check{$chosen}{'hit'}) {
				($counts->{'transcript'}{$type}{$fileid}{'tot'}{'tot'}{'tot'}{'nr'})?
				($counts->{'transcript'}{$type}{$fileid}{'tot'}{'tot'}{'tot'}{'nr'}++):($counts->{'transcript'}{$type}{$fileid}{'tot'}{'tot'}{'tot'}{'nr'}=1);
				#duplicate, but it's easier that way:
				($counts->{'transcript'}{$type}{$fileid}{'tot'}{'tot'}{'tot'}{'tot'})?
				($counts->{'transcript'}{$type}{$fileid}{'tot'}{'tot'}{'tot'}{'tot'}++):($counts->{'transcript'}{$type}{$fileid}{'tot'}{'tot'}{'tot'}{'tot'}=1); 
			}	
			$check{$gid}{'hit'}=1;
#			$check{$chosen}{'hit'}=1;
			
			#Do the repeats stuff if relevant
			unless ($full eq "n") {
				my @l = split(/\s+/,$l);	
				next LINE unless ($ilen >= $inters);
				my @rm = split(";",$l[3]);
				my $Rnam = $rm[9];
				my ($Rcla,$Rfam) = TEshuffle::get_Rclass_Rfam($Rnam,$rm[10]);
				#Increment in the data structure, but only if relevant = avoid counting hits several times
				unless ($checkTE{$l[14]}{$cat}{$type}{'tot'}) {
					($counts->{$cat}{$type}{$fileid}{'tot'}{'tot'}{'tot'}{'tot'})?
					($counts->{$cat}{$type}{$fileid}{'tot'}{'tot'}{'tot'}{'tot'}++):($counts->{$cat}{$type}{$fileid}{'tot'}{'tot'}{'tot'}{'tot'}=1);
				}	
				unless ($checkTE{$l[14]}{$cat}{$type}{$Rcla}) {
					($counts->{$cat}{$type}{$fileid}{$Rcla}{'tot'}{'tot'}{'tot'})?
					($counts->{$cat}{$type}{$fileid}{$Rcla}{'tot'}{'tot'}{'tot'}++):($counts->{$cat}{$type}{$fileid}{$Rcla}{'tot'}{'tot'}{'tot'}=1);			
				}
				unless ($checkTE{$l[14]}{$cat}{$type}{$Rcla.$Rfam}) {
					($counts->{$cat}{$type}{$fileid}{$Rcla}{$Rfam}{'tot'}{'tot'})?
					($counts->{$cat}{$type}{$fileid}{$Rcla}{$Rfam}{'tot'}{'tot'}++):($counts->{$cat}{$type}{$fileid}{$Rcla}{$Rfam}{'tot'}{'tot'}=1);
				}
				unless ($checkTE{$l[14]}{$cat}{$type}{$Rcla.$Rfam.$Rnam}) {
					($counts->{$cat}{$type}{$fileid}{$Rcla}{$Rfam}{$Rnam}{'tot'})?
					($counts->{$cat}{$type}{$fileid}{$Rcla}{$Rfam}{$Rnam}{'tot'}++):($counts->{$cat}{$type}{$fileid}{$Rcla}{$Rfam}{$Rnam}{'tot'}=1);	
				}				
				#Need to check if a feature is counted several times in the upper classes
				$checkTE{$l[14]}{$cat}{$type}{'tot'}=1;
				$checkTE{$l[14]}{$cat}{$type}{$Rcla}=1;
				$checkTE{$l[14]}{$cat}{$type}{$Rcla.$Rfam}=1;
				$checkTE{$l[14]}{$cat}{$type}{$Rcla.$Rfam.$Rnam}=1;
				
				#Age categories if any; only increment per exon & age category, not per TE
				if ($age->{$Rnam}) {
					unless ($checkTE{$l[14]}{'age'}) { #easier to load tot hit with these keys for the print_out sub
						($counts->{$cat}{$type}{$fileid}{'age'}{'cat.1'}{'tot'}{'tot'})?
						($counts->{$cat}{$type}{$fileid}{'age'}{'cat.1'}{'tot'}{'tot'}++):($counts->{$cat}{$type}{$fileid}{'age'}{'cat.1'}{'tot'}{'tot'}=1); 
					}
					unless ($checkTE{$l[14]}{$age->{$Rnam}[4]}) {
						($counts->{$cat}{$type}{$fileid}{'age'}{'cat.1'}{$age->{$Rnam}[4]}{'tot'})?
						($counts->{$cat}{$type}{$fileid}{'age'}{'cat.1'}{$age->{$Rnam}[4]}{'tot'}++):($counts->{$cat}{$type}{$fileid}{'age'}{'cat.1'}{$age->{$Rnam}[4]}{'tot'}=1);
					}
					if (($age->{$Rnam}[5]) && (! $checkTE{$l[14]}{$age->{$Rnam}[5]})) {
						($counts->{$cat}{$type}{$fileid}{'age'}{'cat.2'}{$age->{$Rnam}[5]}{'tot'})?
						($counts->{$cat}{$type}{$fileid}{'age'}{'cat.2'}{$age->{$Rnam}[5]}{'tot'}++):($counts->{$cat}{$type}{$fileid}{'age'}{'cat.2'}{$age->{$Rnam}[5]}{'tot'}=1);
					}
					$checkTE{$l[14]}{'age'}=1;
					$checkTE{$l[14]}{$age->{$Rnam}[4]}=1;
					$checkTE{$l[14]}{$age->{$Rnam}[5]}=1;
				}
				$counts->{$cat}{$type}{$fileid}{'age'}{'cat.2'}{'tot'}{'tot'}=$counts->{$cat}{$type}{$fileid}{'age'}{'cat.1'}{'tot'}{'tot'};
			}	
		}
	}
	close ($fh);
		
	#Add to total exons the number of exons in a random transcript of each gene that was not already counted		
	GENE: foreach my $gene (keys %{$infos}) {
		next GENE if (defined $checkgenes{$gene});
		my $chosen = random_tr($infos,$gene,$type);
		$total_exons->{$type}{$fileid}{'tot'}+=$infos->{$gene}{$type}{$chosen}{'nb'};
	}
	return ($counts,$total_exons);
}

#-----------------------------------------------------------------------------
sub random_tr {
	my ($infos,$gene_id,$type) = @_;
	my @trid = keys (%{$infos->{$gene_id}{$type}});
	my $r = int(rand(scalar(@trid)));
	return ($trid[$r]);
}

#-----------------------------------------------------------------------------
sub overlap_category {
	my ($l,$infos,$gid,$type,$trid) = @_;
	my ($Trstart,$Trend,$Trex) = ($infos->{$gid}{$type}{$trid}{'st'},$infos->{$gid}{$type}{$trid}{'en'},$infos->{$gid}{$type}{$trid}{'nb'});
	my ($Gstart,$Gend,$st,$en,$strand) = ($l->[1],$l->[2],$l->[9],$l->[10],$l->[12]);	
	my $cat = "exonized"; #the default

	#Check the TSS_polyA with Tr corrdinates first, indep of strand. Could also be below with overhang both sides, but cleaner to double check with transcript coordinates
	$cat = "TSS_polyA" if (($Gstart<$Trstart) && ($Gend>$Trend));
	
	#Now the rest; easiest is to set what are the exons
	my $ExType = "MIDDLE";
	$ExType = "FIRST" if ((($strand eq "+") && ($st == $Trstart)) || (($strand eq "-") && ($en == $Trend)));
	$ExType = "LAST" if ((($strand eq "+") && ($en == $Trend)) || (($strand eq "-") && ($st == $Trstart)));
	
	if ($Gstart < $st) {
		if ($Gend > $en) { # overhang TE start AND end side
			$cat = "3SPL_exon_5SPL";
			$cat = "TSS_5SPL" if ($ExType eq "FIRST");
			$cat = "3SPL_polyA" if ($ExType eq "LAST");
		} else {  #overhang TE start side only
			($strand eq "+")?($cat = "3SPL"):($cat = "5SPL");
			$cat = "TSS" if (($strand eq "+") && (($Trex == 1) || ($ExType eq "FIRST")));
			$cat = "polyA" if (($strand eq "-") && (($Trex == 1) || ($ExType eq "LAST")));			
		}
	} elsif ($Gend > $en) { # => overhang only end side
		($strand eq "+")?($cat = "5SPL"):($cat = "3SPL");
		$cat = "polyA" if (($strand eq "+") && (($Trex == 1) || ($ExType eq "LAST")));
		$cat = "TSS" if (($strand eq "-") && (($Trex == 1) || ($ExType eq "FIRST")));		
	}
	return ($cat);
}

#-----------------------------------------------------------------------------
sub print_stats {
	my ($out,$no_boot,$more,$no_boot_tot_ex,$boots,$nboot,$boot_tot_ex,$parsedRM,$age,$full,$scriptname,$version) = @_;

	#get the boot and no_boot total_exons values, avg and sd
	print STDERR "     Get number of exons (total and hit)\n";
	my $no_boot_exons = get_exon_data($no_boot_tot_ex);
	my $boot_exons = get_exon_data($boot_tot_ex);	
			
	#now print
	print_cat_data($no_boot,$no_boot_exons,$more,$boots,$boot_exons,$nboot,$out,$scriptname,$version);
	print_rep_data($no_boot,$no_boot_exons,$more,$boots,$boot_exons,$nboot,$out,$parsedRM,$age,$scriptname,$version) if ($full eq "y");
	return 1;
}

#-----------------------------------------------------------------------------
sub print_cat_data {
	my ($no_boot,$no_boot_exons,$more,$boots,$boot_exons,$nboot,$out,$scriptname,$version) = @_;
	#get the no_boot values, avg and sd
	print STDERR "     Get data for each category of overlap\n";
	my $obs = ();
	$obs = get_cat_data($no_boot,0,"na",$out,$obs);
	my $exp = initialize_cat_exp($obs); #0 values for all the ones seen in obs => so that even if not seen in exp, will be there
	$exp = get_cat_data($boots,$nboot,$obs,$out,$exp);
	
	my $midval = $nboot/2;
	open (my $fh, ">", $out.".cat.txt") or confess "ERROR (sub print_stats): can't open to write $out.cat.txt $!\n";	
	print $fh "#Script $scriptname, v$version\n";
	print $fh "#Aggregated results + stats\n";
	print $fh "#With $more repetitions for obs (observed) and $nboot bootstraps for exp (expected); sd = standard deviation; nb = number; len = length; avg = average\n";
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
			$obsper = $obs->{$cat}{$type}{'avg'}/$no_boot_exons->{$type}{'avg'}*100 unless ($no_boot_exons->{$type}{'avg'} == 0);		
			my $expper = 0;
			$expper = $exp->{$cat}{$type}{'avg'}/$boot_exons->{$type}{'avg'}*100 unless ($boot_exons->{$type}{'avg'} == 0);
			my $sign = TEshuffle::get_sign($pval);
			print $fh "$type\t$o{$cat}\t$cat\t$obs->{$cat}{$type}{'avg'}\t$obs->{$cat}{$type}{'sd'}\t$obsper\t$no_boot_exons->{$type}{'avg'}\t$no_boot_exons->{$type}{'sd'}\t";
			print $fh "$exp->{$cat}{$type}{'avg'}\t$exp->{$cat}{$type}{'sd'}\t$expper\t$boot_exons->{$type}{'avg'}\t$boot_exons->{$type}{'sd'}\t$exp->{$cat}{$type}{'rank'}\t$pval\t$sign\n";		
		}
	}
	close $fh;
	return 1;
}

#-----------------------------------------------------------------------------
sub print_rep_data {
	my ($no_boot,$no_boot_exons,$more,$boots,$boot_exons,$nboot,$out,$parsedRM,$age,$scriptname,$version) = @_;
	print STDERR "     Get data for each repeat, family and class (total and per category)\n";
	my $te_obs = ();
	$te_obs = get_te_data($no_boot,0,"na",$parsedRM,$te_obs);	
	my $te_exp = initialize_te_exp($te_obs); #0 values for all the ones seen in obs => so that even if not seen in exp, will be there
	$te_exp = get_te_data($boots,$nboot,$te_obs,$parsedRM,$te_exp);
	$te_exp = TEshuffle::binomial_test_R($te_exp,"tr");
	my $midval = $nboot/2;
	open (my $fh, ">", $out.".TE.txt") or confess "ERROR (sub print_stats): can't open to write $out.TEs.txt $!\n";	
	print $fh "#Script $scriptname, v$version\n";
	print $fh "#Aggregated results + stats\n";
	print $fh "#With $more repetitions for obs (observed) and $nboot bootstraps for exp (expected)\n";
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
						print $fh "$parsedRM->{$Rclass}{$Rfam}{$Rname}\t"; 
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
		foreach my $round (keys %{$tot_ex->{$type}}) {
			push(@data,$tot_ex->{$type}{$round}{'tot'});	
		}
		#get average and standard deviation from @data
		($exons{$type}{'avg'},$exons{$type}{'sd'}) = TEshuffle::get_avg_and_sd(\@data);
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
	my ($all_data,$n,$obs,$out,$cat_data) = @_;
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
				EXP: foreach my $exp (@data) {
					last EXP if ($exp > $obs->{$cat}{$type}{'avg'});
					$rank++ if ($exp < $obs->{$cat}{$type}{'avg'});
				}
				$cat_data->{$cat}{$type}{'rank'}=$rank;
				if ($rank <= $nboot/2) {
					$cat_data->{$cat}{$type}{'pval'}=$rank/$nboot*2; #*2 because 2 tailed
				} else {
					$cat_data->{$cat}{$type}{'pval'}=($nboot+2-$rank)/$nboot*2;  #+2 so it is symetrical (around nboot+1)
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
	my ($counts,$nboot,$te_obs,$parsedRM,$te_data) = @_;
	# $counts->{$cat}{$type}{$fileid}{'tot'}{'tot'}{'tot'}{'tot'}
	# $counts->{$cat}{$type}{$fileid}{$Rcla}{'tot'}{'tot'}{'tot'}
	# $counts->{$cat}{$type}{$fileid}{$Rcla}{$Rfam}{'tot'}{'tot'}
	# $counts->{$cat}{$type}{$fileid}{$Rcla}{$Rfam}{$Rnam}{'tot'}

	#agregate data [I could easily have the categories details here, so put them; so I can compare]
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
							push(@{$nb_r->{$cat}{$type}{$Rclass}{$Rfam}{$Rname}},$counts->{$cat}{$type}{$round}{$Rclass}{$Rfam}{$Rname}{'tot'});
							push(@{$nb_a1->{$cat}{$type}{$Rclass}{$Rfam}{$Rname}},$counts->{$cat}{$type}{$round}{$Rclass}{$Rfam}{$Rname}{'tot'}) 
								if (($Rclass eq "age") && ($Rfam eq "cat.1"));
							push(@{$nb_a2->{$cat}{$type}{$Rclass}{$Rfam}{$Rname}},$counts->{$cat}{$type}{$round}{$Rclass}{$Rfam}{$Rname}{'tot'}) 
								if (($Rclass eq "age") && ($Rfam eq "cat.2"));
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
					$te_data = get_te_data_details($cat,$type,$Rclass,"tot","tot",$nb_c->{$cat}{$type}{$Rclass}{'tot'}{'tot'},$te_data,$te_obs,$nboot,$parsedRM) 
						if ($Rclass ne "age");	
					foreach my $Rfam (keys %{$counts->{$cat}{$type}{$round}{$Rclass}}) {
						$te_data = get_te_data_details($cat,$type,$Rclass,$Rfam,"tot",$nb_f->{$cat}{$type}{$Rclass}{$Rfam}{'tot'},$te_data,$te_obs,$nboot,$parsedRM) 
							if ($Rclass ne "age");	
						foreach my $Rname (keys %{$counts->{$cat}{$type}{$round}{$Rclass}{$Rfam}}) {
							$te_data = get_te_data_details($cat,$type,$Rclass,$Rfam,$Rname,$nb_r->{$cat}{$type}{$Rclass}{$Rfam}{$Rname},$te_data,$te_obs,$nboot,$parsedRM);
							$te_data = get_te_data_details($cat,$type,$Rclass,$Rfam,$Rname,$nb_a1->{$cat}{$type}{$Rclass}{$Rfam}{$Rname},$te_data,$te_obs,$nboot,$parsedRM) 
								if (($Rclass eq "age") && ($Rfam eq "cat.1"));
							$te_data = get_te_data_details($cat,$type,$Rclass,$Rfam,$Rname,$nb_a2->{$cat}{$type}{$Rclass}{$Rfam}{$Rname},$te_data,$te_obs,$nboot,$parsedRM) 
								if (($Rclass eq "age") && ($Rfam eq "cat.2"));
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
	my ($cat,$type,$key1,$key2,$key3,$agg_data,$te_data,$te_obs,$nboot,$parsedRM) = @_;	
	#get average and sd of the expected	
	($te_data->{$cat}{$type}{$key1}{$key2}{$key3}{'avg'},$te_data->{$cat}{$type}{$key1}{$key2}{$key3}{'sd'}) = TEshuffle::get_avg_and_sd($agg_data);		
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
		if ($rank <= $nboot/2) {
			$te_data->{$cat}{$type}{$key1}{$key2}{$key3}{'pval'}=$rank/$nboot*2; #*2 because 2 tailed
		} else {
			$te_data->{$cat}{$type}{$key1}{$key2}{$key3}{'pval'}=($nboot+2-$rank)/$nboot*2;  #+2 so it is symetrical (around nboot+1)
		}
				
		#Binomial test
		#get all the values needed for binomial test in R => do them all at once
		my $n = $parsedRM->{$key1}{$key2}{$key3} if ($parsedRM->{$key1}{$key2}{$key3});
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


