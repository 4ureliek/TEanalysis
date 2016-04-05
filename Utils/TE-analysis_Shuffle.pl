#!/usr/bin/perl
#######################################################
# Author  :  Zev Kronenberg (https://github.com/zeeev), for the v1.0
#            Modifications by Aurelie Kapusta (https://github.com/4ureliek) after v2.0
# email   :  4urelie.k@gmail.com  
# Purpose :  Writen to generate data (observed vs expected) shown in Figure 5, 
#            Kapusta et al. 2013 PLoS Genetics
#            (http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1003470)
#            But highly modified to be more useful
#######################################################
use strict;
use warnings;
use Carp;
use Getopt::Long;
use Bio::SeqIO;
use GAL::Annotation; #if issues, there is an alternative subroutine not using this module, see usage
use Math::CDF;
#use Data::Dumper;

#-----------------------------------------------------------------------------
#------------------------------- DESCRIPTION ---------------------------------
#-----------------------------------------------------------------------------
#flush buffer
$| = 1;

my $version = "4.2";
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
\n";

# TO DO:
#  - allow any bed or gff file for -s

my $usage = "
Synopsis (v$version):

    perl $scriptname -l lncRNA.gff -p prot.gff [-o <X>] [-m <X>] -s features_to_shuffle [-n <X>] [-f] -e exclude.range [-d] -r genome.range [-b] [-i include.range] [-a] [-l no_low] [-t <type,name>] [-c] [-w <bedtools_path>] [-v] [-h]

    /!\\ REQUIRES - Bedtools, v2.25+
	              - GAL::Annotation version later than Jan 2016 [update of is_coding]
	                see https://github.com/The-Sequence-Ontology/GAL
	                If issues, open the script in a text editor and comment line 15,	                
	                as well as the lines with read_gff subroutine, and uncomment lines with subroutine load_gene_tr
    /!\\ Previous outputs, if any, will be moved as *.previous [which only saves results once]
  
  CITATION:
    - You may cite Kapusta et al. (2013) PLoS Genetics (DOI: 10.1371/journal.pgen.1003470)
      but should also include the GitHub link to this script
    - for BEDtools, Quinlan AR and Hall IM (2010) Bioinformatics (DOI: 10.1093/bioinformatics/btq033)

  DESCRIPTION:
    - Features provided in -s will be overlapped with -p and/or -l files,
       without (no_boot) or with (boot) shuffling (on same chromosome)
       One of -p or -l is mandatory. Having both in the same run means that they are 
       intersected with the same TE files, which may be better for comparisons, 
       but does not seem necessary with high bootstraps.
       
    - A random transcript per gene is selected: use -m to do several repetitions of no_boot
    
    - For each bootstrap (-n) with shuffling features in -s, transcripts are randomly selected as well
       Note that high bootstraps takes a lot of time.
       Shuffling is done by default with allowing overlaps between shuffled features, bceause it is
       faster and OK when over representation of specific repeats are considered.  
       Note that because TEs are often fragmented + there are inversions, the counts for the exonized TEs is likely inflated;
       this also means that when TEs are shuffled, there are more fragments than TEs. Some should be moved non independently, 
       or the input file should be corrected when possible to limit that issue [not implemented in this script for now]
                  
    - Note that one exon may have several types of overlaps (e.g. \"SPL\" and \"exonized\"),
       but each exon is counted only one time for each category (important for \"exonized\").
       However for the results per repeat, each hit is counted, unless it's the same TE
   
    - If you need to generate the <genome.gaps> file but you would also like to add more files to the -e option, 
       just do a first run with no bootstraps (in this example the genome.range is also being generated):
       perl ~/bin/$scriptname -f input.bed -s genome.out -r genome.fa -b -e genome.fa -d -n 0

    - Two-tailed permutation test is done on the counts of overlaps for categories
      and the results are in a *.stats.cat.txt file
      If -f is used then stats are also made on each repeat, with two-tailed 
      permutation and binomial tests and the results are in a *.stats.TE.txt file.
      Note that the output *.stats.cat.txt is basically included in the output *.stats.TE.txt,
      with values of tot tot tot in the columns Rclass, Rfam and Rname
      The use of -f will take longer but requires fewer bootsraps, 
      because binomial test is more sensitive
  
  MANDATORY ARGUMENTS:
    -p,--prot     => (STRING) protein coding gff3 file; one of -p or -l is mandatory
    -l,--lnc      => (STRING) lncRNAs gff3 file; one of -p or -l is mandatory
    -s,--shuffle  => (STRING) Features to shuffle = TE file
                              Repeat masker .out or the .bed file generated by the TE-analysis_pipeline                            
    -r,--range    => (STRING) To know the maximum value in a given chromosome/scaffold. 
                              File should be: Name \\t length
                              Can be files from UCSC, files *.chrom.sizes
                              If you don't have such file, use -b (--build) and provide the genome fasta file for -r                              
    -e,--excl     => (STRING) This will be used as -excl for bedtools shuffle: \"coordinates in which features from -i should not be placed.\"
                              More than one file may be provided (comma separated), they will be concatenated 
                              (in a file = first-file-name.cat.bed).
                              By default, at least one file is required = assembly gaps, and it needs to be the first file
                              if not in bed format. Indeed, you may provide the UCSC gap file, with columns as:
                                  bin, chrom, chromStart, chromEnd, ix, n, size, type, bridge
                              it will be converted to a bed file. Additionally, you may provide the genome file in fasta format
                              and add the option -d (--dogaps), to generate a bed file corresponding to assembly gaps.
                              Other files may correspond to regions of low mappability, for example for hg19:
                              http://www.broadinstitute.org/~anshul/projects/encode/rawdata/blacklists/hg19-blacklist-README.pdf
                              Notes: -> when the bed file is generated by this script, any N stretch > 50nt will be considered as a gap 
                                        (this can be changed in the load_gap subroutine)         
                                     -> 3% of the shuffled feature may overlap with these regions 
                                        (this can be changed in the shuffle subroutine).	
  OPTIONAL ARGUMENTS:
    -o,--overlap  => (INT)    Minimal length (in nt) of intersection in order to consider the TE included in the feature.
                              Default = 10 (to match the TEanalysis-pipeline.pl)
    -m,--more     => (INT)    Even in the no_boot, a random transcript is picked. Set this number to do repetitions for no_boot.
                              Default = 1 (still need it done 1 time; set this to 0 is equivalent to 1)
    -n,--nboot    => (STRING) number of bootsraps with shuffled -s file
                              Default = 100 for faster runs; use higher -n for good pvalues 
                              (-n 10000 is best for permutation test but this will take a while)
                              If set to 0, no bootstrap will be done
    -f,--full     => (BOOL)   Use -f to also do stats for each repeat separately (separated output, with binomial test as well)
                              Results will be in a file *.stats.TE.txt
                              Note that the output *.stats.cat.txt is basically included in the output *.stats.TE.txt,
                              with values of tot tot tot in the columns Rclass, Rfam and Rname
                              This will take longer but requires fewer bootsraps, because binomial test is more sensitive
    -b,--build    => (BOOL)   See above; use this and provide the genome fasta file if no range/lengths file (-r)
                              This step may take a while but will create the required file	
    -d,--dogaps   => (BOOL)   See above; use this and provide the genome fasta file if no gap file (-g)
                              If several files in -e, then the genome needs to be the first one.
                              This step is not optimized, it will take a while (but will create the required file)                       
    -i,--incl     => (STRING) To use as -incl for bedtools shuffle: \"coordinates in which features from -i should be placed.\"
                              Bed of gff format. Could be intervals close to transcripts for example.
                              More than one file (same format) may be provided (comma separated), 
                              they will be concatenated (in a file = first-file-name.cat.bed)
    -a,--add      => (BOOL)   to add the -noOverlapping option to the bedtools shuffle command line, 
                              and therefore NOT allow overlaps between the shuffled features.
                              This may create issues mostly if -i is used (space to shuffle into smaller than features to shuffle)
    -u,--u        => (STRING) To set the behavior regarding non TE sequences: all, no_low, no_nonTE, none
                                 -t all = keep all non TE sequences (no filtering)
                                 -t no_low [default] = keep all besides low_complexity and simple_repeat
                                 -t no_nonTE = keep all except when class = nonTE
                                 -t none = everything is filtered out (nonTE, low_complexity, simple_repeat, snRNA, srpRNA, rRNA, tRNA/tRNA, satellite)
    -t,--te       => (STRING) <type,name>
                              run the script on only a subset of repeats. Not case sensitive.
                              The type can be: name, class or family and it will be EXACT MATCH unless -c is chosen as well
                              ex: -a name,nhAT1_ML => only fragments corresponding to the repeat named exactly nhAT1_ML will be looked at
                                  -a class,DNA => all repeats with class named exactly DNA (as in ...#DNA/hAT or ...#DNA/Tc1)
                                  -a family,hAT => all repeats with family named exactly hAT (so NOT ...#DNA/hAT-Charlie for example)
    -c,--contain  => (BOOL)   to check if the \"name\" determined with -filter is included in the value in Repeat Masker output, instead of exact match
                              ex: -a name,HERVK -c => all fragments containing HERVK in their name
                                  -a family,hAT -c => all repeats with family containing hAT (...#DNA/hAT, ...#DNA/hAT-Charlie, etc)
    -w,--where    => (STRING) if BEDtools are not in your path, provide path to BEDtools bin directory
    -v,--version  => (BOOL)   print the version
    -h,--help     => (BOOL)   print this usage
\n";

#-----------------------------------------------------------------------------
#------------------------------ LOAD AND CHECK -------------------------------
#-----------------------------------------------------------------------------

my ($shuffle,$full,$exclude,$dogaps,$build,$dobuild,$f_regexp,$allow,$nooverlaps,$v,$help);
my ($prot,$linc) = ("n","n");
my $inters = 10;
my $more = 0;
my $nboot = 10;
my $incl = "na";
my $nonTE = "no_low";
my $filter = "na";
my $bedtools = "";
my $catout = "y"; #removed from options, not really relevant to ask for choice
my $opt_success = GetOptions(
			 	  'prot=s'		=> \$prot,
			 	  'lnc=s'		=> \$linc,
			 	  'more=s'      => \$more,
			 	  'shuffle=s'   => \$shuffle,
			 	  'overlap=s'   => \$inters,
			 	  'nboot=s'     => \$nboot,
			 	  'full'        => \$full,
			 	  'range=s'     => \$build,
			 	  'build'       => \$dobuild,
			 	  'excl=s'		=> \$exclude,
			 	  'dogaps'      => \$dogaps,
			 	  'incl=s'		=> \$incl,
			 	  'add'		    => \$nooverlaps,
			 	  'u=s'		    => \$nonTE,
			 	  'te=s'		=> \$filter,
			 	  'contain'     => \$f_regexp,
			 	  'where=s'     => \$bedtools,
			 	  'version'     => \$v,
			 	  'help'		=> \$help,);
			 	  
#Check options, if files exist, etc
die "\n --- $scriptname version $version\n\n" if $v;
die $usage if $help || ! $opt_success;
die $usage unless (($shuffle) && ($exclude) && (($prot ne "n") || ($linc ne "n" )));
die "\n -p $prot does not exist?\n\n"  if (($prot ne "n") && (! -e $prot));
die "\n -p $prot is not a gff file?\n\n" unless (($prot eq "n") || ($prot =~ /\.gff$/) || ($prot =~ /\.gff3$/));
die "\n -l $linc does not exist?\n\n"  if (($linc ne "n") && (! -e $linc));
die "\n -l $linc is not a gff file?\n\n" unless (($linc eq "n") || ($linc =~ /\.gff$/) || ($linc =~ /\.gff3$/));
die "\n -s $shuffle does not exist?\n\n"  if (! -e $shuffle);
die "\n -s $shuffle is not in a proper format (not .out, .bed, .gff or .gff)?\n\n" unless (($shuffle =~ /\.out$/) || ($shuffle =~ /\.bed$/) || ($shuffle =~ /\.gff$/) || ($shuffle =~ /\.gff3$/));
die "\n -i $incl does not exist?\n\n"  if (($incl ne "na") && ($exclude !~ /,/) && (! -e $incl));
die "\n -r $build does not exist?\n\n"  if (! -e $build);
die "\n -e $exclude does not exist?\n\n"  if (($exclude !~ /,/) && (! -e $exclude));
die "\n -n $nboot but should be an integer\n\n" if ($nboot !~ /\d+/);
die "\n -i $inters but should be an integer\n\n" if ($inters !~ /\d+/);
die "\n -m $more but should be an integer\n\n" if ($more !~ /\d+/);
die "\n -w $bedtools does not exist?\n\n"  if (($bedtools ne "") && (! -e $bedtools));
die "\n -t requires 2 values separated by a coma (-t <name,filter>; use -h to see the usage)\n\n" if (($filter ne "na") && ($filter !~ /,/));
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
print STDERR "\n --- $scriptname v$version\n";

#Genome range
print STDERR " --- loading build (genome range)\n";
my ($okseq,$build_file) = load_build($build,$dobuild);

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
$exclude[0] = load_gap($exclude[0],$dogaps);
print STDERR "     concatenating files for -e\n" if ($exclude =~ /,/);
my $excl;
($exclude =~ /,/)?($excl = concat_beds(\@exclude)):($excl = $exclude[0]);

#If relevant, files to include for shuffling
if (($incl ne "na") && ($incl =~ /,/)) {
	print STDERR " --- concatenating $incl files to one file\n";
	my @include = split(",",$incl);
	$incl = concat_beds(\@include);
}

#Now features to shuffle
print STDERR " --- checking file in -s, print in .bed if not a .bed or gff file\n";
print STDERR "     filtering TEs based on filter ($filter) and non TE behavior ($nonTE)\n" unless ($filter eq "na");
print STDERR "     + getting genomic counts for each repeat\n" unless ($full eq "n");
my ($toshuff_file,$parsedRM) = RMtobed($shuffle,$okseq,$filter,$f_regexp,$nonTE,$full);

#Load the gff files
print STDERR " --- Load gene IDs / transcript IDs for for:\n";
my $whichgene = ();
my $l_tr = ();
my $p_tr = ();
print STDERR "  -> $linc\n" unless ($linc eq "n");
($l_tr,$whichgene) = read_gff($linc, $okseq, 0, $whichgene) unless ($linc eq "n"); #comment if GAL::Annotation is a problem
# ($l_tr,$whichgene) = load_gene_tr($linc,$okseq,$whichgene) unless ($linc eq "n"); #uncomment if GAL::Annotation is a problem
print STDERR "  -> $prot\n" unless ($prot eq "n");
($p_tr,$whichgene) = read_gff($prot, $okseq, 1, $whichgene) unless ($prot eq "n"); #comment if GAL::Annotation is a problem
# ($p_tr,$whichgene) = load_gene_tr($prot,$okseq,$whichgene) unless ($prot eq "n"); #uncomment if GAL::Annotation is a problem

#Outputs
my ($stats,$input);
($linc)?($input = $linc):($input = $prot);
my ($f_type,$f_name) = split(",",$filter) unless ($filter eq "na");	
($filter eq "na")?($stats = "$input.nonTE-$nonTE.$nboot.boot.stats"):($stats = "$input.nonTE-$nonTE.$f_name.$nboot.boot.stats");		
my ($outl,$outlb,$temp_l,$outp,$outpb,$temp_p) = ("$linc.no_boot","$linc.boot","$linc.temp","$prot.no_boot","$prot.boot","$prot.temp");
my $temp_s = "$toshuff_file.temp";
cleanup_out($outl,$outp,$outlb,$outpb,$catout,$stats,$temp_l,$temp_p,$temp_s,$nboot,$linc,$prot);

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
	($no_boot,$no_boot_tot_exons) = check_for_featured_overlap("$temp_l/no_boot.joined",$l_tr,"no_boot.".$j,'transcript',$outl,$inters,$no_boot,$no_boot_tot_exons,$whichgene,$full) unless ($linc eq "n");
	($no_boot,$no_boot_tot_exons) = check_for_featured_overlap("$temp_p/no_boot.joined",$p_tr,"no_boot.".$j,'mRNA',$outp,$inters,$no_boot,$no_boot_tot_exons,$whichgene,$full) unless ($prot eq "n");
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
# 		print STDERR "  -- BOOTSTRAP $i\n";	
# 		print STDERR "     Shuffle features with command line:\n";
		my $shuffled = shuffle($toshuff_file,$temp_s,$i,$excl,$incl,$build_file,$bedtools,$nooverlaps);
# 		print STDERR " --- Intersect with command lines:\n";
# 		print STDERR "      $intersectBed -a $shuffled -b $linc -wo > $temp_l/boot.$i.joined\n" unless ($linc eq "n");
 		system "      $intersectBed -a $shuffled -b $linc -wo > $temp_l/boot.$i.joined" unless ($linc eq "n");
# 		print STDERR "      $intersectBed -a $shuffled -b $prot -wo > $temp_p/boot.$i.joined\n" unless ($prot eq "n");
 		system "      $intersectBed -a $shuffled -b $prot -wo > $temp_p/boot.$i.joined" unless ($prot eq "n");
# 		print STDERR " --- Check intersections of $linc and $prot with features in $shuffled (expected)\n";
		($boots,$boots_tot_exons) = check_for_featured_overlap("$temp_l/boot.$i.joined",$l_tr,"boot.".$i,'transcript',$outlb,$inters,$boots,$boots_tot_exons,$whichgene,$full) unless ($linc eq "n");
		($boots,$boots_tot_exons) = check_for_featured_overlap("$temp_p/boot.$i.joined",$p_tr,"boot.".$i,'mRNA',$outpb,$inters,$boots,$boots_tot_exons,$whichgene,$full) unless ($prot eq "n");;
		`cat $outlb >> $catout.boot.txt` if (($catout) && (-e $outlb));
		`cat $outpb >> $catout.boot.txt` if (($catout) && (-e $outpb));
		`rm -Rf $shuffled`; #these files are now not needed anymore, all is stored
		`rm -Rf $temp_l/boot.$i.joined` unless ($linc eq "n");
		`rm -Rf $temp_p/boot.$i.joined` unless ($prot eq "n");
	}
}
`rm -Rf $temp_s`; #these folders are not needed anymore
`rm -Rf $temp_l` unless ($linc eq "n");
`rm -Rf $temp_p` unless ($prot eq "n");

#Stats now
print STDERR " --- Get and print stats\n" if ($nboot > 0);
print_stats($stats,$no_boot,$more,$no_boot_tot_exons,$boots,$nboot,$boots_tot_exons,$parsedRM,$full,$scriptname,$version) if ($nboot > 0);

#end
print STDERR " --- $scriptname done\n";
print STDERR "     Stats for categeories printed in: $stats.cat.txt\n" if ($nboot > 0);
print STDERR "     Stats for TEs printed in: $stats.TE.txt\n" if (($nboot > 0) && ($full eq "y"));
print STDERR "\n";
exit;


#-----------------------------------------------------------------------------
#-------------------------------- SUBROUTINES --------------------------------
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Get build if needed and get chromosomes included in it
# my ($okseq,$build_file) = load_build($build,$dobuild);
#-----------------------------------------------------------------------------
sub load_build{ 
    my ($file,$dobuild) = @_;
    my %build = ();
    my $build_file;
    if ($dobuild eq "n") {
    	$build_file = $file;
		open (my $in, '<', $file) or confess "ERROR (sub load_build): can't open to read $file $!\n";
		GAP_LINE: while(<$in>){
			chomp;
			my @l = split /\t/, $_;
			$build{$l[0]} = $l[1]; 
		}
		close $in;
    } else {
		print STDERR "     printing chrom. range file for $file\n";
    	$build_file = $file.".build.tab";
    	open (my $lenfh, ">", $build_file) or confess "ERROR (sub load_build): can't open to write $build_file $!\n";
		my $fa = Bio::SeqIO->new(-file => $file, -format => "fasta") or confess "ERROR (sub parse_build): can not open as fasta object $file $!\n";
		my %replen = ();
		while( my $seq = $fa->next_seq() ) {
			my $name = $seq->display_id;
			my $len = $seq->length;
			$build{$name} = $len;
			print $lenfh "$name\t$len\n";
		}
		undef $fa;
		close $lenfh;    
    }
    return (\%build,$build_file);
}

#-----------------------------------------------------------------------------
# loading assembly gaps
# my $excl = load_gap($gaps,$dogaps);
#-----------------------------------------------------------------------------
sub load_gap {
    my ($file,$dogaps) = @_;
    my $bed;  
	if ($dogaps eq "y") { 
	    #gaps not in file yet, get them from fasta file in bed format
		$bed = print_gap_bed($file);
	} else {
		#file with gaps, check if bed already, or if bed already made
		($file =~ /\.bed$/)?($bed = $file):($bed = $file.".bed");
		if ($file !~ /\.bed$/) {
			if (-e $bed) {
				print STDERR "     $bed exists, skipping printing $file in bed format\n";
			} else {
				print STDERR "     printing bed file corresponding to assembly gaps in $file (limiting to features > 50nt)\n";
				open (my $in, '<', $file) or confess "ERROR (sub gap_file): can't open to read $file $!\n";
				open (my $out, '>', $bed) or confess "ERROR (sub gap_file): can't open to write $bed $!\n";
				LINE: while(<$in>){
					chomp(my $l = $_);
					next LINE if (substr($l,0,1) eq "#");
					my @l = split(/\t/,$l);			
					#bin, chrom, chromStart, chromEnd, ix, n, size, type, bridge)
					print $out "$l[1]\t$l[2]\t$l[3]\n" unless ($l[3]-$l[2] < 50);
				}
				close $in;
				close $out;
			}	
		}
	}
    return ($bed);
}

#-----------------------------------------------------------------------------
# printing assembly gaps if needed, in bed format
# $file = print_gap($file) if ($ifbuild eq "n");
#-----------------------------------------------------------------------------
sub print_gap_bed {
	my $file = shift;
	print STDERR "     loading assembly gaps from $file (stretches of Ns > 50nt)\n";
	my $fasta = Bio::SeqIO->new(-file => $file, -format => "fasta") or confess "ERROR (sub print_gap_bed): failed to create SeqIO object from $file $!\n";
	my $gaps = "$file.gaps.bed";
	open my $gapfh, ">", $gaps or confess "ERROR (sub print_gap): could not open to write $gaps $!\n";
	while( my $seqfa = $fasta->next_seq() ) {
		my $head = $seqfa->display_id;
		my $seq = $seqfa->seq;
		$seq = uc($seq); #uppercases just in case
		my @seq = split("",$seq);
		my $Nstart = 0;
		my $Nlen = 0;
		NUCL: for (my $i = 0; $i <= $#seq; $i++){
			#go to next nt if this is not a gap to extend or first nt after a gap
			my $nt = $seq[$i];
			next NUCL if (($nt ne "N") && ($Nlen == 0));
			if ($nt eq "N") {
				$Nstart = $i+1 if ($Nlen == 0); #do not increment if extending gap
				$Nlen++; #will be 1 if was 0
			} else { #end of a gap, print and reintialize (already skipped if length is < 10)
				my $Nend = $Nlen + $Nstart - 1;
				print $gapfh "$head\t$Nstart\t$Nend\n" unless ($Nlen < 10);
				$Nstart = 0;
				$Nlen = 0;
			}
		}
	}
	undef $fasta;
	close $gapfh;
	return ($gaps);
}

#-----------------------------------------------------------------------------
# Concatenate files with cat
# my $excl = concat_beds(\@exclude) if ($exclude =~ /,/);
# $incl = concat_beds(\@include);
#-----------------------------------------------------------------------------
sub concat_beds {
	my $files = shift;
	my $concat = $files->[0].".cat.bed";
	my $list = "";
	foreach my $file (@{$files}) {
		#system "sed -i -e '$a\' $file"; #add newline if not one, otherwise it will create issues with cat, but problem here with the $, perl is looking for a variable
		my $check = `tail -1 $file | wc -l`;
		chomp($check);
		`echo >> $file` if ($check == 0);
		print STDERR "     $file has no newline at the end (tail -1 $file | wc -l returned $check)\n" if ($check == 0);
		$list = $list." $file";
	}
	`cat $list > $concat`;
	return ($concat);
}

#-----------------------------------------------------------------------------
# Convert RMoutput .out file to bed if needed
# my ($toshuff_file,$parsedRM) = RMtobed($shuffle,$okseq,$filter,$f_regexp,$nonTE,$full);
#-----------------------------------------------------------------------------
sub RMtobed {
	my ($file,$okseq,$filter,$f_regexp,$full) = @_;
	my $parsed = ();
	my ($f_type,$f_name) = split(",",$filter) unless ($filter eq "na");	
	my $ok = $file;
	$ok = $1 if ($file =~ /(.*)\.out$/);	
	$ok = $1 if ($file =~ /(.*)\.bed$/);	
	($filter eq "na")?($ok = $ok.".nonTE-$nonTE.bed"):($ok = $ok.".$f_name.bed");		
	if (-e $ok) { #if has been filtered same way, OK, just get dictionary
		print STDERR "     -> $ok exists, just getting dictionary from it\n" unless ($full eq "n");
		$parsed = getparsedRM($ok,$parsed,"file") unless ($full eq "n");
		return ($ok,$parsed);
	}
	print STDERR "     -> $file is in bed format, but $ok does not exist; generating it...\n" if ($file =~ /(.*)\.bed$/);
	#now it means RM.out, or that the proper bed file does not exist => make it a bed file + load the parsing hash	
	open(my $fh, "<$file") or confess "\n   ERROR (sub RMtobed): could not open to read $file $!\n";
	open(my $bed_fh, ">$ok") or confess "\n   ERROR (sub RMtobed): could not open to write $ok $!\n";
	LINE: while(<$fh>) {
		chomp(my $l = $_);
		my @l = ();
		my ($Rname,$classfam);
		if ($file =~ /(.*)\.out$/) { #if .out	
			$l =~ s/^\s+//;
			next LINE if (($l =~ /^[Ss]core|^SW|^#/) || ($l !~ /\w/));
			$l = $1 if ($l =~ /^(.*)\*$/); #remove the star
			@l = split('\s+',$l);			
		} else { #if .bed file
			my @all = split('\s+',$l);
			@l = split(";",$all[3]);
		}	
		next LINE unless (defined $okseq->{$l[4]}); #if not in build of stuff OK to shuffle on, remove here as well
		$l[8] =~ s/C/-/; #correct strand to match convention
		($Rname,$classfam) = ($l[9],$l[10]);	
		my ($Rclass,$Rfam) = get_Rclass_Rfam($Rname,$classfam);
		
		#now filter non TE or not, based on -nonTE value
		next LINE if (($nonTE eq "none") && ($Rclass =~ /nonTE|snRNA|rRNA|tRNA|snoRNA|scRNA|srpRNA|[Ll]ow_complexity|[Ss]imple_repeat|[Ss]atellite/));
		next LINE if (($nonTE eq "no_nonTE") && ($Rclass =~ /nonTE/));
		next LINE if (($nonTE eq "no_low") && ($Rclass =~ /[Ll]ow_complexity|[Ss]imple_repeat/));

		#filter out stuff from -filter if relevant
		if ($filter ne "na") {				
			if ($f_regexp eq "y") {
				#check if what's det as the filter is included in the names
				next LINE unless ((($f_type eq "name") && ($Rname =~ /$f_name/))
							   || (($f_type eq "class") && ($Rclass =~ /$f_name/)) 
							   || (($f_type eq "family") && ($Rfam =~ /$f_name/)));
			} else {
				next LINE unless ((($f_type eq "name") && ($f_name eq $Rname))
							   || (($f_type eq "class") && ($f_name eq $Rclass)) 
							   || (($f_type eq "family") && ($f_name eq $Rfam)));
			}	
		}
		
		#create unique ID = the RMout line
		my ($chr,$start,$end,$strand) = ($l[4],$l[5],$l[6],$l[8]);
		my $ID = $l[0];
		for (my $i=1; $i<=$#l;$i++) {
			$ID = $ID.";".$l[$i];
		}
		$ID =~ s/\s//; #should not need this, but just in case
		
# 		#correct the start if base is 0
# 		$start=$start+1 if ($base == 0);
		
		#print in bed format
		print $bed_fh "$chr\t$start\t$end\t$ID\t.\t$strand\n"; #with ID being the whole line => easy to acces to RMoutput original info
		
		#get dictionary
		$parsed = getparsedRM(\@l,$parsed,"line") unless ($full eq "n");
	}
	close ($fh);
	close ($bed_fh);
	return ($ok,$parsed);
}

#-----------------------------------------------------------------------------
# Get infos from RMout
# $parsed = getparsedRM($file,$parsed,"file");
# $parsed = getparsedRM(\@l,$parsed,"line");
#-----------------------------------------------------------------------------
sub getparsedRM {
	my ($info,$parsed,$type) = @_;
	if ($type eq "line") { #meaning it's read from the .out while being converted in .bed
		$parsed = getparsedRMline($parsed,$info);
	} else { #meaning it's already a bed file => open and read it
		open(my $fh, "<$info") or confess "\n   ERROR (sub getparsedRM): could not open to read $info!\n";
		LINE: while(<$fh>) {
			chomp(my $l = $_);
			$l =~ s/^\s+//;
			next LINE if (($l =~ /^[Ss]core|^SW|^#/) || ($l !~ /\w/));
			my @l = split('\s+',$l);
			my @RMline = split(";",$l[3]);
			$parsed = getparsedRMline($parsed,\@RMline);
		}
		close $fh;	
	}
	return($parsed);
}

#-----------------------------------------------------------------------------
# Get counts and length for each TE
# $parsed = getparsedRMline($parsed,$info)
# $parsed = getparsedRMline($parsed,\@RMline)
#-----------------------------------------------------------------------------
sub getparsedRMline {
	my ($parsed,$l) = @_;
	my ($Rname,$classfam) = ($l->[9],$l->[10]);
	my ($Rclass,$Rfam) = get_Rclass_Rfam($Rname,$classfam);
	#now feed dictionary
	($parsed->{'tot'}{'tot'}{'tot'})?($parsed->{'tot'}{'tot'}{'tot'}++):($parsed->{'tot'}{'tot'}{'tot'}=1);
	($parsed->{$Rclass}{'tot'}{'tot'})?($parsed->{$Rclass}{'tot'}{'tot'}++):($parsed->{$Rclass}{'tot'}{'tot'}=1);
	($parsed->{$Rclass}{$Rfam}{'tot'})?($parsed->{$Rclass}{$Rfam}{'tot'}++):($parsed->{$Rclass}{$Rfam}{'tot'}=1);
	($parsed->{$Rclass}{$Rfam}{$Rname})?($parsed->{$Rclass}{$Rfam}{$Rname}++):($parsed->{$Rclass}{$Rfam}{$Rname}=1);
	return $parsed;
}

#-----------------------------------------------------------------------------
# Load geneIDs{type}{trIDs}
# my $l_tr = load_gene_tr($lincs,$okseq,$whichgene);
# my $p_tr = load_gene_tr($prot,$okseq,$whichgene);
#-----------------------------------------------------------------------------
sub load_gene_tr {
	#Not as solid as using GAL::Annotation, but it is an alternative
	my ($file,$okseq,$whichgene) = @_;
	my %tr = ();
	my ($gid,$trid,$type);
	open(my $fh, "<$file") or confess "\n   ERROR (sub load_gene_tr): could not open to read $file!\n";
	LINE: while(<$fh>) {
		chomp(my $l = $_);
		next LINE if (substr($l,0,1) eq "#");
		my @l = split('\s+',$l);
		next LINE unless (defined $okseq->{$l[0]}); #if not in build of stuff OK to shuffle on, remove here as well
		my $id = $l[8];
		$id = $1 if $id =~ /^ID=(.+?);/;
		if ($l[2] eq "gene") {
			$gid = $id;
		} elsif ($l[2] eq "transcript") {
			$trid = $id;
			$type = "transcript";
			$type = "mRNA" if ($l =~ /protein_coding/); #if coding, it should be in the "gene_type" or the "transcript_type"
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
# Load geneIDs{type}{trIDs} with GAL::Annotation
# ($l_tr,$whichgene) = read_gff($linc, $okseq, 0, $whichgene);
# ($p_tr,$whichgene) = read_gff($linc, $okseq, 0, $whichgene);
#-----------------------------------------------------------------------------
sub read_gff {
	my ($gff3_file,$okseq,$coding,$whichgene) = @_;
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
		next GENE unless (defined $okseq->{$seqid});
		
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
	return (\%trinfo,$whichgene);
}

#-----------------------------------------------------------------------------
# Cleanup outputs
# cleanup_out($outl,$outp,$outlb,$outpb,$catout,$stats,$temp_l,$temp_p,$temp_s,$nboot,$linc,$prot);
#-----------------------------------------------------------------------------
sub cleanup_out {
	my ($outl,$outp,$outlb,$outpb,$catout,$stats,$temp_l,$temp_p,$temp_s,$nboot,$linc,$prot) = @_;
	`mv $outl $outl.previous` if (-e $outl);
	`mv $outp $outp.previous` if (-e $outp);
	`mv $outlb $outlb.previous` if (-e $outlb);
	`mv $outpb $outpb.previous` if (-e $outpb);
	`mv $stats.cat.txt $stats.cat.txt.previous` if (-e $stats.".cat.txt");
	`mv $stats.TE.txt $stats.TE.txt.previous` if (-e $stats.".TE.txt");
	`mv $stats.details.txt $stats.details.txt.previous` if (-e $stats.".details.txt");
	`rm -Rf $temp_l` if (-e $temp_l);
	`rm -Rf $temp_p` if (-e $temp_p);
	`rm -Rf $temp_s` if (-e $temp_s);
	`mkdir $temp_l` unless ($linc eq "n");
	`mkdir $temp_p` unless ($prot eq "n");
	`mkdir $temp_s` if ($nboot > 0);
	return 1;
}	

#-----------------------------------------------------------------------------
# Shuffle, with bedtools
# my $shuffled = shuffle($toshuff_file,$temp_s,$i,$excl,$incl,$build_file,$bedtools,$nooverlaps);
#-----------------------------------------------------------------------------
sub shuffle {
	my ($toshuff_file,$temp_s,$nb,$excl,$incl,$build,$bedtools,$nooverlaps) = @_;
	my $out = $temp_s."/shufffled".$nb;
	my $bed_shuffle = $bedtools."shuffleBed";
#  	($incl eq "na")?(print STDERR "      $bed_shuffle -i $toshuff_file -excl $excl -f 2 $shuff_allow_overlaps -g $build -chrom -maxTries 10000 > $out\n"):
#  	                (print STDERR "      $bed_shuffle -incl $incl -i $toshuff_file -excl $excl -f 2 $shuff_allow_overlaps -g $build -chrom -maxTries 10000 > $out\n");
	($incl eq "na")?(system "$bed_shuffle -i $toshuff_file -excl $excl -f 2 $nooverlaps -g $build -chrom -maxTries 10000 > $out"):
	                (system "$bed_shuffle -incl $incl -i $toshuff_file -excl $excl -f 2 $nooverlaps -g $build -chrom -maxTries 10000 > $out");	
	return ($out);

}

#-----------------------------------------------------------------------------
# Check overlap with TEs and count
# Keep only one transcript per gene -> requires to load transcript IDs per gene
#  	($no_boot,$no_boot_tot_exons) = check_for_featured_overlap($linc,$l_tr,"no_boot.".$j,'transcript',$outl,$inters,$no_boot,$no_boot_tot_exons,$whichgene,$full);
#	($no_boot,$no_boot_tot_exons) = check_for_featured_overlap($prot,$p_tr,"no_boot.".$j,'mRNA',$outp,$inters,$no_boot,$no_boot_tot_exons,$whichgene,$full);
# 	($boots,$boots_tot_exons) = check_for_featured_overlap($linc,$l_tr,"boot.".$i,'transcript',$outlb,$inters,$boots,$boots_tot_exons,$whichgene,$full);
# 	($boots,$boots_tot_exons) = check_for_featured_overlap($prot,$p_tr,"boot.".$i,'mRNA',$outpb,$inters,$boots,$boots_tot_exons,$whichgene,$full);
#-----------------------------------------------------------------------------
sub check_for_featured_overlap {
	my ($file,$infos,$fileid,$type,$out,$inters,$counts,$total_exons,$whichgene,$full) = @_;
	my %chosen_tr = ();
	my %checkgenes = ();
	my %check = ();

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
			$total_exons->{$type}{$fileid}{'tot'}+=$infos->{$gid}{$type}{$chosen}{'nb'} unless (defined $checkgenes{$gid}); #increment with exon numbers of the transcript chosen for this gene
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
				($counts->{$cat}{$type}{$fileid}{'tot'}{'tot'}{'tot'}{'nr'})?($counts->{$cat}{$type}{$fileid}{'tot'}{'tot'}{'tot'}{'nr'}++):($counts->{$cat}{$type}{$fileid}{'tot'}{'tot'}{'tot'}{'nr'}=1);
			}
			$check{$l[14]}{$cat}=1;
			$total_exons->{$type}{$fileid}{'hit'}++ unless (defined $check{$l[14]}{'hit'});
			$check{$l[14]}{'hit'}=1;
			unless (defined $check{$chosen}{'hit'}) { #counting each tr hit only one time
				($counts->{'transcript'}{$type}{$fileid}{'tot'}{'tot'}{'tot'}{'nr'})?($counts->{'transcript'}{$type}{$fileid}{'tot'}{'tot'}{'tot'}{'nr'}++):($counts->{'transcript'}{$type}{$fileid}{'tot'}{'tot'}{'tot'}{'nr'}=1);
				($counts->{'transcript'}{$type}{$fileid}{'tot'}{'tot'}{'tot'}{'tot'})?($counts->{'transcript'}{$type}{$fileid}{'tot'}{'tot'}{'tot'}{'tot'}++):($counts->{'transcript'}{$type}{$fileid}{'tot'}{'tot'}{'tot'}{'tot'}=1); #duplicate, but it's easier that way
			}	
			$check{$chosen}{'hit'}=1;
			
			#Do the repeats stuff if relevant
			unless ($full eq "n") {
				my @l = split(/\s+/,$l);	
				next LINE unless ($ilen >= $inters);
				my @rm = split(";",$l[3]);
				my $Rnam = $rm[9];
				my ($Rcla,$Rfam) = get_Rclass_Rfam($Rnam,$rm[10]);
#				($counts->{'tot'}{$type}{$fileid}{'tot'}{'tot'}{'tot'}{'tot'})?($counts->{'tot'}{$type}{$fileid}{'tot'}{'tot'}{'tot'}{'tot'}++):($counts->{'tot'}{$type}{$fileid}{'tot'}{'tot'}{'tot'}{'tot'}=1);
				($counts->{$cat}{$type}{$fileid}{'tot'}{'tot'}{'tot'}{'tot'})?($counts->{$cat}{$type}{$fileid}{'tot'}{'tot'}{'tot'}{'tot'}++):($counts->{$cat}{$type}{$fileid}{'tot'}{'tot'}{'tot'}{'tot'}=1);
				($counts->{$cat}{$type}{$fileid}{$Rcla}{'tot'}{'tot'}{'tot'})?($counts->{$cat}{$type}{$fileid}{$Rcla}{'tot'}{'tot'}{'tot'}++):($counts->{$cat}{$type}{$fileid}{$Rcla}{'tot'}{'tot'}{'tot'}=1);
				($counts->{$cat}{$type}{$fileid}{$Rcla}{$Rfam}{'tot'}{'tot'})?($counts->{$cat}{$type}{$fileid}{$Rcla}{$Rfam}{'tot'}{'tot'}++):($counts->{$cat}{$type}{$fileid}{$Rcla}{$Rfam}{'tot'}{'tot'}=1);
				($counts->{$cat}{$type}{$fileid}{$Rcla}{$Rfam}{$Rnam}{'tot'})?($counts->{$cat}{$type}{$fileid}{$Rcla}{$Rfam}{$Rnam}{'tot'}++):($counts->{$cat}{$type}{$fileid}{$Rcla}{$Rfam}{$Rnam}{'tot'}=1);	
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
# Get random transcript
# $chosen_tr = random_tr($infos,$gene_id,$type);
#-----------------------------------------------------------------------------
sub random_tr {
	my ($infos,$gene_id,$type) = @_;
	my @trid = keys (%{$infos->{$gene_id}{$type}});
	my $r = int(rand(scalar(@trid)));
	return ($trid[$r]);
}

#----------------------------------------------------------------------------
# Get class fam
# my ($Rclass,$Rfam) = get_Rclass_Rfam($Rname,$classfam);
#----------------------------------------------------------------------------
sub get_Rclass_Rfam {
	my($Rname,$classfam) = @_;
	my ($Rclass,$Rfam);
	if ($classfam =~ /\//) {
		my $incaseof;
		($Rclass,$Rfam,$incaseof) = split(/\//, $classfam); #should not happen but in case there was a / in family
		$Rfam = $Rfam."_".$incaseof if ($incaseof);
	} else {
		$Rfam = $classfam;
		$Rfam=~ s/^(.*)\..*$/$1/;
		$Rclass = $classfam;
		$Rclass =~ s/^.*\.(.*)$/$1/;
	}
	$Rfam = $Rfam."--int" if (($Rfam =~ /ERV/) && (($Rname =~ /[-_][iI]nt/) || ($Rname =~ /[-_]I$/)));
	return ($Rclass,$Rfam);
}

#-----------------------------------------------------------------------------
# Get the overlap category
# my $cat = overlap_category($l,$infos,$gid,$type,$trid);
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
# Print Stats (permutation test + binomial test)
# print_stats($stats,$no_boot,$more,$no_boot_tot_exons,$boots,$nboot,$boots_tot_exons,$parsedRM,$full,$scriptname,$version) if ($nboot > 0);
#-----------------------------------------------------------------------------
sub print_stats {
	my ($out,$no_boot,$more,$no_boot_tot_ex,$boots,$nboot,$boot_tot_ex,$parsedRM,$full,$scriptname,$version) = @_;

	#get the boot and no_boot total_exons values, avg and sd
	print STDERR "     Get number of exons (total and hit)\n";
	my $no_boot_exons = get_exon_data($no_boot_tot_ex);
	my $boot_exons = get_exon_data($boot_tot_ex);	
			
	#now print cat data
	print_cat_data($no_boot,$no_boot_exons,$more,$boots,$boot_exons,$nboot,$out,$scriptname,$version);
	print_rep_data($no_boot,$no_boot_exons,$more,$boots,$boot_exons,$nboot,$out,$parsedRM,$scriptname,$version) if ($full eq "y");
	return 1;
}
sub print_cat_data {
	my ($no_boot,$no_boot_exons,$more,$boots,$boot_exons,$nboot,$out,$scriptname,$version) = @_;
	#get the no_boot values, avg and sd
	print STDERR "     Get data for each category of overlap\n";
	my $obs = get_cat_data($no_boot,0,"na",$out);
	my $exp = get_cat_data($boots,$nboot,$obs,$out);
	
	my $midval = $nboot/2;
	open (my $fh, ">", $out.".cat.txt") or confess "ERROR (sub print_stats): can't open to write $out.cat.txt $!\n";	
	print $fh "#Script $scriptname, v$version\n";
	print $fh "#Aggregated results + stats\n";
	print $fh "#With $more repetitions for obs (observed) and $nboot bootstraps for exp (expected); sd = standard deviation\n";
	print $fh "#The p value shows how significant the difference between observed and random is (two tailed permutation test)\n";
	print $fh "#If rank is < $midval and pvalue is not \"na\", there are significantly fewer observed values than expected \n";
	print $fh "#If rank is > $midval and pvalue is not \"na\", there are significantly higher observed values than expected \n";
	print $fh "#The category \"transcript\" corresponds to at least one feature hit per transcript\n";
	print $fh "#For all categories besides \"transcript\", counts are of exons\n";
	print $fh "\n#trancript_type\tcagtegory_id\toverlap_category\tobs_mean\tobs_sd\t%_obs\tobs_tot\tobs_tot_sd\texp_mean\texp_sd\t%_exp\texp_tot\texp_tot_sd\tobs_rank_in_exp\t2-tailed_permutation-test_pvalue(obs.vs.exp)\tsignificance\n\n";
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
			my $sign = get_sign($pval);
			print $fh "$type\t$o{$cat}\t$cat\t$obs->{$cat}{$type}{'avg'}\t$obs->{$cat}{$type}{'sd'}\t$obsper\t$no_boot_exons->{$type}{'avg'}\t$no_boot_exons->{$type}{'sd'}\t$exp->{$cat}{$type}{'avg'}\t$exp->{$cat}{$type}{'sd'}\t$expper\t$boot_exons->{$type}{'avg'}\t$boot_exons->{$type}{'sd'}\t$exp->{$cat}{$type}{'rank'}\t$pval\t$sign\n";		
		}
	}
	close $fh;
	return 1;
}
sub print_rep_data {
	my ($no_boot,$no_boot_exons,$more,$boots,$boot_exons,$nboot,$out,$parsedRM,$scriptname,$version) = @_;
	print STDERR "     Get data for each repeat, family and class (total and per category)\n";
	my $te_obs = get_te_data($no_boot,0,"na",$out,$parsedRM);	
	my $te_exp = get_te_data($boots,$nboot,$te_obs,$out,$parsedRM);	
	my $midval = $nboot/2;
	open (my $fh, ">", $out.".TE.txt") or confess "ERROR (sub print_stats): can't open to write $out.TEs.txt $!\n";	
	print $fh "#Script $scriptname, v$version\n";
	print $fh "#Aggregated results + stats\n";
	print $fh "#With $more repetitions for obs (observed) and $nboot bootstraps for exp (expected)\n";
	print $fh "sd = standard deviation; nb = number; avg = average\n";
	print $fh "#Two tests are made (permutation and binomial) to assess how significant the difference between observed and random, so two pvalues are given\n";
	print $fh "#The binomial test is not 2 sided, so it only tests for enrichment in observed values\n";
	print $fh "#A pvalue cannot be 0, so when the probability is 1 (meaning enrichment of this repeat in the observed set), pvalue is set to na, but the significance is ***\n";
	print $fh "#For the two tailed permutation test:\n";
	print $fh "#if rank is < $midval and pvalue is not \"ns\", there are significantly fewer observed values than expected \n";
	print $fh "#if rank is > $midval and pvalue is not \"ns\", there are significantly higher observed values than expected \n";
	print $fh "\n#\t#\tLevel_(tot_means_all)\t#\t#\tCOUNTS\t#\t#\t#\t#\t#\t#\t#\t#\t#\t#\t#\t#\t#\t#\t#\n";
	print $fh "#Type\tCategory\tRclass\tRfam\tRname\t";
	print $fh "obs_nb_of_hits\tobs_nb_sd\t%_obs_nb_(%of_features)\tobs_tot_nb_of_hits\tobs_tot_sd\t";
	print $fh "exp_nb_of_hits\texp_nb_sd\t%_exp_nb_(%of_features)\texp_tot_nb_of_hits\texp_tot_sd\t";
	print $fh "obs_rank_in_exp\t2-tailed_permutation-test_pvalue(obs.vs.exp)\tsignificance\tnb_binomal_test_proba(obs.vs.exp)\tnb_binomal_test_pvalue(1-proba)\tsignificance\n\n";

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
						my $pval_nb = $te_exp->{$cat}{$type}{$Rclass}{$Rfam}{$Rname}{'pval'};		
						$pval_nb = "na" if (($te_expavg == 0) && ($te_obsnb == 0));									
						#Now print stuff
						print $fh "$type\t$cat\t$Rclass\t$Rfam\t$Rname\t";
						print $fh "$te_obsnb\t$te_obssd\t$te_obsper\t$no_boot_exons->{$type}{'avg'}\t$no_boot_exons->{$type}{'sd'}\t"; 
						print $fh "$te_expavg\t$te_exp->{$cat}{$type}{$Rclass}{$Rfam}{$Rname}{'sd'}\t$te_expper\t$boot_exons->{$type}{'avg'}\t$boot_exons->{$type}{'sd'}\t";			
						my $sign = get_sign($pval_nb);				
						print $fh "$te_exp->{$cat}{$type}{$Rclass}{$Rfam}{$Rname}{'rank'}\t$pval_nb\t$sign\t";
						$te_exp->{$cat}{$type}{$Rclass}{$Rfam}{$Rname}{'binom_proba'} = "na" unless ($te_exp->{$cat}{$type}{$Rclass}{$Rfam}{$Rname}{'binom_proba'});
						$te_exp->{$cat}{$type}{$Rclass}{$Rfam}{$Rname}{'binom_pval'} = "na" unless ($te_exp->{$cat}{$type}{$Rclass}{$Rfam}{$Rname}{'binom_pval'});
						$sign = get_sign($te_exp->{$cat}{$type}{$Rclass}{$Rfam}{$Rname}{'binom_pval'});
						$sign = "***" if (($te_exp->{$cat}{$type}{$Rclass}{$Rfam}{$Rname}{'binom_proba'} ne "na") && ($te_exp->{$cat}{$type}{$Rclass}{$Rfam}{$Rname}{'binom_proba'} == 1));
						print $fh "$te_exp->{$cat}{$type}{$Rclass}{$Rfam}{$Rname}{'binom_proba'}\t$te_exp->{$cat}{$type}{$Rclass}{$Rfam}{$Rname}{'binom_pval'}\t$sign\n";	
					}
				}		
			}
		}
	}
	close $fh;
    return 1;
}

#-----------------------------------------------------------------------------
# Get counts from exon_tot hash
# my $no_boot_exons = get_exon_data($no_boot_tot_ex);
# my $boot_exons = get_exon_data($boot_tot_ex);	
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
		($exons{$type}{'avg'},$exons{$type}{'sd'}) = get_avg_and_sd(\@data);
	}
	return(\%exons);
}

#-----------------------------------------------------------------------------
# Get data for each category + for repeats; separate them just for clarity
# my $obs = get_cat_data($no_boot,0,"na",$out);
# my $exp = get_cat_data($boots,$nboot,$obs,$out);
#-----------------------------------------------------------------------------
sub get_cat_data {
# $counts->{$cat}{$type}{$fileid}{'tot'}{'tot'}{'tot'}{'nr'}
	my ($all_data,$n,$obs,$out) = @_;
	my %cat_data = ();
	foreach my $cat (keys %{$all_data}) {
		foreach my $type (keys %{$all_data->{$cat}}) {
			my @data = ();
			foreach my $round (keys %{$all_data->{$cat}{$type}}) {
				my $hits = $all_data->{$cat}{$type}{$round}{'tot'}{'tot'}{'tot'}{'nr'};
				$hits = 0 unless ($hits);
				push(@data,$hits);	
			}
			#get average and standard deviation from @data			
			($cat_data{$cat}{$type}{'avg'},$cat_data{$cat}{$type}{'sd'}) = get_avg_and_sd(\@data);
			
			#Now get he rank of the observed value in the list of expected => get a p value
			unless ($n == 0) {
				($obs->{$cat}{$type}{'avg'},$obs->{$cat}{$type}{'sd'}) = (0,"na") unless ($obs->{$cat}{$type}{'avg'});
				my $rank = 1; #pvalue can't be 0
				@data = sort {$a <=> $b} @data;
				EXP: foreach my $exp (@data) {
					last EXP if ($exp > $obs->{$cat}{$type}{'avg'});
					$rank++ if ($exp < $obs->{$cat}{$type}{'avg'});
				}
				$cat_data{$cat}{$type}{'rank'}=$rank;
				if ($rank <= $nboot/2) {
					$cat_data{$cat}{$type}{'pval'}=$rank/$nboot*2; #*2 because 2 tailed
				} else {
					$cat_data{$cat}{$type}{'pval'}=(abs($nboot-$rank+1))/$nboot*2; #+1 because pvalue can't be 0
				}			
			}
		}
	}
	return(\%cat_data);
}

#-----------------------------------------------------------------------------
# Get the stats values for each repeats now -> includes binomial
# my $te_obs = get_te_data($no_boot,0,"na",$out,$parsedRM);
# my $te_exp = get_te_data($boots,$nboot,$te_obs,$out,$parsedRM);
#-----------------------------------------------------------------------------
sub get_te_data {
	my ($counts,$n,$te_obs,$out,$parsedRM) = @_;
	# $counts->{'tot'}{$type}{$fileid}{'tot'}{'tot'}{'tot'}{'tot'}
	# $counts->{$cat}{$type}{$fileid}{'tot'}{'tot'}{'tot'}{'tot'}
	# $counts->{$cat}{$type}{$fileid}{$Rcla}{'tot'}{'tot'}{'tot'}
	# $counts->{$cat}{$type}{$fileid}{$Rcla}{$Rfam}{'tot'}{'tot'}
	# $counts->{$cat}{$type}{$fileid}{$Rcla}{$Rfam}{$Rnam}{'tot'}
	my $te_data = ();

	#agregate data [I could easily have the categories details here, so put them; so I can compare]
	my ($nb_c,$len_c,$nb_f,$len_f,$nb_r,$len_r) = ();
	foreach my $cat (keys %{$counts}) {
		foreach my $type (keys %{$counts->{$cat}}) {
			foreach my $round (keys %{$counts->{$cat}{$type}}) {
				push(@{$nb_c->{$cat}{$type}{'tot'}{'tot'}{'tot'}},$counts->{$cat}{$type}{$round}{'tot'}{'tot'}{'tot'}{'tot'});	
				foreach my $Rclass (keys %{$counts->{$cat}{$type}{$round}}) {
					push(@{$nb_c->{$cat}{$type}{$Rclass}{'tot'}{'tot'}},$counts->{$cat}{$type}{$round}{$Rclass}{'tot'}{'tot'}{'tot'});	
					foreach my $Rfam (keys %{$counts->{$cat}{$type}{$round}{$Rclass}}) {
						push(@{$nb_f->{$cat}{$type}{$Rclass}{$Rfam}{'tot'}},$counts->{$cat}{$type}{$round}{$Rclass}{$Rfam}{'tot'}{'tot'});	
						foreach my $Rname (keys %{$counts->{$cat}{$type}{$round}{$Rclass}{$Rfam}}) {
							push(@{$nb_r->{$cat}{$type}{$Rclass}{$Rfam}{$Rname}},$counts->{$cat}{$type}{$round}{$Rclass}{$Rfam}{$Rname}{'tot'});	
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
					$te_data = get_te_data_details($cat,$type,$Rclass,"tot","tot",$nb_c->{$cat}{$type}{$Rclass}{'tot'}{'tot'},$te_data,$te_obs,$nboot,$parsedRM,$n);
					foreach my $Rfam (keys %{$counts->{$cat}{$type}{$round}{$Rclass}}) {
						$te_data = get_te_data_details($cat,$type,$Rclass,$Rfam,"tot",$nb_f->{$cat}{$type}{$Rclass}{$Rfam}{'tot'},$te_data,$te_obs,$nboot,$parsedRM,$n);
						foreach my $Rname (keys %{$counts->{$cat}{$type}{$round}{$Rclass}{$Rfam}}) {
							$te_data = get_te_data_details($cat,$type,$Rclass,$Rfam,$Rname,$nb_r->{$cat}{$type}{$Rclass}{$Rfam}{$Rname},$te_data,$te_obs,$nboot,$parsedRM,$n);
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
# sub get_te_data_details
# called by get_tes_data, to get average, sd, rank and p value for all the lists
# and stats when boots
#-----------------------------------------------------------------------------	
sub get_te_data_details {
	my ($cat,$type,$key1,$key2,$key3,$agg_data,$te_data,$te_obs,$nboot,$parsedRM,$n) = @_;	
	#get average and sd of the expected
	($te_data->{$cat}{$type}{$key1}{$key2}{$key3}{'avg'},$te_data->{$cat}{$type}{$key1}{$key2}{$key3}{'sd'}) = get_avg_and_sd($agg_data);
	$te_data->{$cat}{$type}{$key1}{$key2}{$key3}{'avg'} = 0 unless ($te_data->{$cat}{$type}{$key1}{$key2}{$key3}{'avg'});
	$te_data->{$cat}{$type}{$key1}{$key2}{$key3}{'sd'} = "na" unless ($te_data->{$cat}{$type}{$key1}{$key2}{$key3}{'sd'});
	
	if ($n > 0) {
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
			$te_data->{$cat}{$type}{$key1}{$key2}{$key3}{'pval'}=(abs($nboot-$rank+1))/$nboot*2; #+1 because pvalue can't be 0
		}
		
		#Binomial test
		# http://search.cpan.org/~callahan/Math-CDF-0.1/CDF.pm
		# $prob = pbinom($x, $n, $p)
		# Generates cumulative probabilities from the Binomial distribution. 
		# This is the probability of having $x or fewer successes in $n trials when each trial has a $p probability of success. 
		# $x should be in the range of [0,$n], $n should be in the range (0,Inf) and $p should be in [0,1].			
		my $n = $parsedRM->{$key1}{$key2}{$key3};
		my $x = $observed;
		my $p = $te_data->{$cat}{$type}{$key1}{$key2}{$key3}{'avg'}/$n;
		$te_data->{$cat}{$type}{$key1}{$key2}{$key3}{'binom_proba'} = Math::CDF::pbinom($x, $n, $p);
		$te_data->{$cat}{$type}{$key1}{$key2}{$key3}{'binom_pval'} = 1 - $te_data->{$cat}{$type}{$key1}{$key2}{$key3}{'binom_proba'};	
	}
	return($te_data);
}

#-----------------------------------------------------------------------------
# sub get_avg_and_sd, get average and standard deviation from a list
#-----------------------------------------------------------------------------	
sub get_avg_and_sd{
	my($data) = @_;
#   confess "ERROR (sub print_stats/get_avg_and_sd): empty data array $!\n" if (! @$data);
    return ($data->[0],"na") unless ($data->[1]); #returning 0,0 was an issue if -m =1
    
	my $total = 0;
	foreach (@$data) {
		$total += $_;
	}
	my $avg = $total / @$data;

	my $sqtotal = 0;
	foreach(@$data) {
		$sqtotal += ($avg-$_) ** 2;
	}
	my $sd = ($sqtotal / (@$data-1)) ** 0.5;
	return ($avg,$sd);
}

#-----------------------------------------------------------------------------
# sub get_sign, just get the *** correspondance from pvalue
#-----------------------------------------------------------------------------	
sub get_sign {
	my $pval = shift;
	my $sign;
	if ($pval eq "na") {	
		$sign = "na";
	} else {
		$sign = "ns" if ($pval > 0.05);
		$sign = "*" if ($pval <= 0.05);
		$sign = "**" if ($pval <= 0.01);
		$sign = "***" if ($pval <= 0.001);
	}
	return($sign);
}


