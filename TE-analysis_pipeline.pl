#!/usr/bin/perl -w
#######################################################
# Author  :  Aurelie Kapusta
# version :  see below, $version
# email   :  4urelie.k@gmail.com
# Purpose :  Pipeline to analyse TE composition in features (exons of transcripts, coding or non coding, transcription factor binding sites, ChIP-seq data, etc)
#            See documentation for more details
######################################################
BEGIN{
   #what to do on: kill -s ALRM <pid> so I can check where it is if it stalls
   $SIG{ALRM}  = sub {print STDERR "SIGALRM received\n"; print STDERR Carp::longmess; print "\n";};
   #what to do on ^C
   $SIG{INT}  = sub {print STDERR "SIGINT received\n"; print STDERR "\n\n".Carp::longmess; exit;};
   #add a folder in INC
   #unshift(@INC, "~/bin/BioPerl-1.6.901");
}

#load modules
use warnings;
use strict;
use Carp;
use Getopt::Long;
use Data::Dumper;

#keep STDOUT and STDERR from buffering
select((select(STDERR), $|=1)[0]); #make STDERR buffer flush immediately
select((select(STDOUT), $|=1)[0]); #make STDOUT buffer flush immediately

my $version = "4.16";
my $changelog;
set_changelog();
sub set_changelog { 
	$changelog = "
#	v1.0 = Mar 2013
#      [...]
#	v4.0 = Nov 2014
#          - put all in one big pipeline script (use of subroutines), instead of 4 different scripts
#          - correct few bugs that still remained
#          - more flexibility in the input file	+ in required files => more dynamic pipeline depending on what users want to do
#	v4.1 = 5 Jan 2015
#          - bug fix
#          - up and dw: wrong regions were kept after subtraction (if statement changed to unless)
#                       and corrections for the cuts was wrong as well (a cc stupid error)
#          - Now --int added in family name when family = ERV and name of TE has I or int at the end.
#	v4.2 = 11 Feb 2015
#          - bug fix
#              = -nonTE option was not taking argument
#              = some file names stuff
#              = -dir option fixed
#          - added number of different chromosomes/scaffolds in the TE-ratios output
#            Indeed, when are on 2 different sequences, less likely to be 2 pieces of the same TE
#	v4.3 = 12 Feb 2015
#          - bug fix
#              = -addcol (were added in ExSt file but not in TrInfos files)
#	v4.4 = 17 Feb 2015
#          - bug fix
#              = -parse on added columns, was doing the opposite (excluding)...
#	v4.5 = 02 Mar 2015
#          - bug fix
#              -f bed => did not have the TEov filtering
#              TEov: if there is a % in the argument then it will filter on % of the feature coverage and not nt amount
#	v4.6 = 10 Feb 2016
#          - bug fix
#             for -addcol with no -myf (when input = gtf file)
#	v4.7 = 20 Sep 2016
#          - update to read properly recent Gencode gff3
#   v4.8 = 04 Oct 2016
#          - update usage and help
#   v4.9 = 09 Jan 2017
#          - bug fix in the intersection calculation. Fix it, but also change to use the -wo instead
#   v4.10 = 25 Jan 2017
#          - Add count of features from the bed file in summary file
#   v4.11 = 06 Feb 2017
#          - Check on \$add to avoid error \"Use of uninitialized value in concatenation (.) or string....\"
#   v4.12 = 27 Mar 2017
#          - Bug fix \$add, still errors
#          - Bug fix \$listtojoin, was erased when no subtract
#          - Bug fix amounts
#   v4.13 = 08 Jun 2018
#          - Bug fix die errors if -RMparsed file not provided or if a repeat is missing from it
#          - Few cosmetic changes
#   v4.14 = 20 Jul 2018
#          - deal with differences between old and new parseRM.pl outputs
#          - Few cosmetic changes
#   v4.15 = 24 & 26 Jul 2018
#          - small bug fix to avoid dying at the TEratio printing step when a repeat is not in the parsedRM file
#   v4.16 = 15 Aug 2018
#          - fixed that -parse option thing, since the format is more strict now (8 columns input file).
#          - bug fix introduced with v4.14 for parsedRM file loading, if old parsedRM format the masked length was 0

# TO DO: 
#  - check what is used in TEinfoRMP, remove useless stuff 
#  - global vars in uc, and no need to pass them to subs / cleaning writing!
#  - Do a utils script to integrate data from several runs as summary or TEratios, like the Coverage one, but with mosaic plots
#  - When -parse, previous files have to be deleted or it won't actually filter, it's annoying. Solve that.
#  - Fix intron TrInfos - use exon coordinates to check and not introns, to see if SPL overlap. However, it's not really more informative than the exon output. 
\n";
	return 1;
}

my $usage;
set_usage();
sub set_usage {
	$usage = "\nUsage [$version]: 
    perl TE-analysis_pipeline_v4+.pl -i <inputfile> [-dir] -f <format> [-myf <col_details.txt>] 
    -RMout <RepeatMasker.out> [-RMparsed RM.out.parsed] [-base <RM base>] [-TE <TE.tab>] [-TEage] [-nonTE <X>] 
    [-fa <genome.fa>] [-subtract <what-to-subtract>] [-subid <name>] [-noselfsub] [-bedtools <path/to/bins>] 
    [-addcol <col1,col2,etc>] [-filter <col,filter>] [-parse <col_nb,filter>] [-cut <X,X,X>]
    [-v] [-clean] [-chlog] [-h] [-help]
   
   SYNOPSIS
    Type -help for detailed explanations + on how to read the outputs.
    Pipeline to analyse TE composition in features (exons of transcripts, coding or non coding, 
       transcription factor binding sites, ChIP-seq data, etc)
    
   REQUIREMENTS
    BEDtools is required for this pipeline to run (Quinlan AR and Hall IM, 2010. Bioinformatics)
    Available at: https://github.com/arq5x/bedtools2
    
   CITATION - please put the GitHub link in Methods, and cite:
      - For this pipeline using -f gtf, cite Kapusta et al. (2013) PLoS Genetics
      - For this pipeline using -f bed, cite Lynch et al. (2015) Cell Reports
      - For the use of BEDtools, Quinlan AR and Hall IM (2010) Bioinformatics
      
   DEBUGGING
     First thing to do = check that your input files are not encoded as Classic Mac (CR).
     Also, double check the usage and the doc, just in case.
     Then, if you really think your input files are OK, shoot me an email or open an issue on GitHub with the errors, 
     your command line and sample files reproducing the error if possible.
   
   DETAILS OF OPTIONS (MD = mandatory and OPT = optional):
    -i         => MD  - (STRING) input file, see more details below (-f usage)
                         in .bed format for -f bed; in gtf, gff3, or any tabulated file for -f gtf
    -dir       => OPT - (BOOL) add this if -i corresponds to a folder containing files to work with (need to contain ONLY these files!!)
    -f         => MD  - (STRING) This sets the type of analysis (kind of related to the format of the input file)
                         chose between -f gtf (complex) or -f bed (simple)
                         Check -help for more information.
                         gtf (default) = complex analysis (for transcripts - TSS, exons, splicing sites, etc)
                               gtf and gff files will work, but only if all the info is in it (see -help)
                         bed = simple analysis (for TF binding sites, ChIP-seq data etc)
                               5 columns bedfile is required: chr start end unique_ID score/. strand                        
    -myf       => OPT - (STRING) if input file is not formated as a gtf/gff3 or bed: use this option with a text file to set the column numbers.
                         To generate an example of the text file to provide, type this option alone (with the path/name of the file to create)
                         When you run the pipeline, you still need to provide -f to determine the type of analysis.
    -RMout     => MD  - (STRING) repeat masker output file .out
                         Even if you already have the .out.bed file, put .out in command line
                         Obviously requires to be for the same assembly file/version than the input file (ex. mm9 or mm10, hg19 or hg38, etc)
    -TEov      => OPT - (INT) minimal length (in nt) of intersection in order to consider the TE included in the feature.
                         Default = 10
                        (STRING) put a % after a number (FLOAT) to filter out when less than X% of the feature overlaps with the TE
                         Ex: -TEov 80% will skip the line if less than 80% of the feature overlaps with the TE
    -RMparsed  => OPT - (STRING) repeat masker output parsed with parseRM.pl script (with or without the -lib option)
                         Typically: <RMout.out>.parseRM.all-repeats.tab
                         See documentation for more details on how to get this file.
                         If not provided, over represented TE families won't be determined
    -base      => OPT - (INT) 0 or 1. Typically, if the file is from Repeat Masker website chose 1, if from UCSC chose 0.
                         If base 0 is chosen, 1 will be added to start of TE coordinates
                         Default = 1
    -TE        => OPT - (STRING) file with TE information, tab delimited. Minimum columns = Rname Rclass Rfam Rclass/Rfam
                         Class and family info are required to parse correctly the files, and for each element
                         they will be extracted from files set in options, with this priority: -TE > -RMparsed > -RMout
                         This file here is mostly to allow providing TE age info, or to modify (some) TE classes.
                         If the -TEage flag is set, then 7 columns are needed: Rname Rclass Rfam Rclass/Rfam %div AGE Ancient/LineageSpe
                         (6th will be used for age but 7th will be added in the output as well)        
    -TEage     => OPT - (BOOL) add this flag if the file in -TE contains age info
                         If this is not set, no age parsing will be available (\"na\" will be put in some of the corresponding columns)
    -nonTE     => OPT - (STRING) to set the behavior regarding non TE sequences: all, no_low, no_nonTE, none
                         Class info will be extracted from files set in options, with this priority: -TE > -RMparsed > -RMout
                            all = keep all non TE sequences (no filtering)
                            no_low = keep all besides low_complexity and simple_repeat
                            no_nonTE = keep all except when class = nonTE
                            none (default) = everything is filtered out (nonTE, low_complexity, simple_repeat, snRNA, srpRNA, rRNA, tRNA/tRNA, satellite)
    -fa        => OPT - (STRING) corresponding genome file (fasta). 
                         Required only when -f gtf, to make sure that surrounding intergenic sequences won't be outside of sequences
                         If not provided, intergenic regions won't be looked at 
    -subtract  => OPT - (STRING) gtf or bed file to subtract from introns and intergenetic regions. Typically, Ensembl gene annotation.
                         Only relevant if -f gtf is used
    -subid     => OPT - (STRING) \"ID\" to show in output files after subtraction of file set in -subtract (for introns and intergenic regions)
                         Default = sub
    -noselfsub => OPT - (BOOL) chose this option to avoid subtracting the input file from itself (from introns and intergenetic regions)
    -bedtools  => OPT - (STRING) if BEDtools are not in your path, provide path to BEDtools bin directory
    -addcol    => OPT - (STRING) add columns (no limit) of the input file to the outputs, separated by a coma and NO SPACE (ex: -addcol 10,11)
                         Note that first column of the file = 0
                         These added columns can be used for -parse
    -filter    => OPT - (STRING) to filter on a specific column, to keep lines where <filter> is found in the column <col>
                         Can't be used on added columns. Lines not matching it won't be printed at all in any of the files.
                         You can:
                         (i)  use a number corresponding to the -myf file column numbers, to the column number of a bed file, 
                              or to columns 0 to 7 of a gtf file
                              For example, use -filter 0,chr1 to analyze only chromosome 1
                         (ii) use the identifier of the feature if you want to filter a gtf file or a custom file
                              For example, use -filter gene_type,lincRNA or -filter transcript_type,lincRNA to keep only lincRNAs.
                              If input file is a gtf file (with \"transcript\" lines), then you can also use a special filter 
                              = gene_type,intragenic or transcript_type,intragenic to look at everything that is not lincRNA. 
                              This will exclude all lincRNAs but keep all the non coding stuff, unless they are < 200nt [should exclude all small RNAs]
    -parse     => OPT - (STRING) to filter on a specific column before joining with TEs, to keep lines where <filter> is found in the column <col>
                         Can be used ONLY on added columns. Lines not matching it won't be printed at all in any of the files.
                         Typically, this allows to quickly check if some subsets differ (ex: -addcol 9,10,11 -parse 10,liver) 
                         if column 10 has the tissue of maximum expression in the original file.
    -cut       => OPT - (STRING) to set size of intergenic regions analyzed (downstream and upstream), in nt
                         Default = 10000,5000,1000
    -v         => OPT - (BOOL) chose this to make the script talk to you
                         print the version if only option
    -clean     => OPT - (BOOL) to use alone with -i, to delete any previous outputs generated by this script for this input file 
                         If -i was a directory, add the -dir flag
                         (RMout.bed won't be deleted -> delete it manually if needed)
    -chlog     => OPT - (BOOL) to print the change log between versions
    -h         => OPT - (BOOL) to print this usage
    -help      => OPT - (BOOL) to print a more detailed doc for this pipeline
    
";
	return 1;
}

my $longhelp;
set_longhelp();
sub set_longhelp {
	$longhelp = "\nSome documentation for TE-analysis_pipeline.pl [$version]
    Author       :  Aurelie Kapusta 
	Last update  :  04 Oct 2016

	--------------------------------
	   INPUT FILES
	--------------------------------
	With -f bed (typically for TF binding sites, ChIPseq data, etc)
	  -i input file needs to be in the 5 columns bed format:
	     chr start end unique_ID score/. strand
	  If you only have a 3 columns format, just type:
         sed 's/$/   .       ./' peaks.bed > peaks.mod.bed
        (white spaces are tabs, in the command line you can get them by pressing ctrl+v and then the tab key)
    
    With -f gtf (typically for gene or transcripts annotations, coding or not):
       -i input file needs to be in gtf, gff, gff3 or any tabulated format containing the required information
       If a non-standard file is used (e.g. tabulated file from transcript assemblies), type:
          perl TE-analysis_pipeline_v4+.pl -myf myfile.conf
          This will create a .conf file that should be edited to set the columns of the various required info (such as gene_ID, transcript_ID, etc)
          Note that the same column can refer to different features (such as id or name, for genes or transcripts). No need to make 4 columns with the same info in it in the input file!
    
    TE info:
      -RMout should be the repeat masker output of the SAME assembly, in .out format
         The pipeline will create a bed file from it (quite long step if the file is big), 
         with nonTE annotations filtered out (or not, depends on your choice of -nonTE)    
         During that conversion, TE class and family are updated using the files from 1) -TE and 2) -RMparsed if provided.
         Many pre-masked assemblies can be found at http://www.repeatmasker.org/genomicDatasets/RMGenomicDatasets.html
         
      -RMparsed is the repeat masker output file of -RMout but parsed with the parseRM_simple.pl script (with or without the -lib option)
         This script can be found at: https://github.com/4ureliek/Parsing-RepeatMasker-Outputs
         (you can run it with or without -lib option will work, column numbers will be corrected based on column headers)
         Typically, the file is <RMout.out>.parseRM.all-repeats.tab
         If not provided, over represented TE families won't be determined
         Some pre-parsed RM outputs can be found in the Data directory of this pipeline

	NOTE THAT: This pipeline expects that all repeat names are different, there will be issues if they are not
         (this may happen for user's de novo libraries used for RM annotation; also, in older libraries there were 
         repeats with same name but different class/fam - if this happens then the first occurence will be the one in the hash with TE info;
         you can also correct that by using -TE). Note that repeat names will be matched in lower cases, so
         if some repeat names are different just thanks to the case that will be a problem.
         When the ERV can be identified as internal (-int or -I associated) then family of the TE is renamed with --int

	--------------------------------
	   FILTERING
	--------------------------------
	-filter
	-addcol and -parse
	
	This pipeline expects that info in the added columns are by transcripts. 
	If they correspond to exons, then the Exon Structure (ExSt) file and TrInfos output files will have it wrong. 
	However the -parse will work correctly if you use -clean beforehand (or manually delete the ExSt file) 

	Note that by default, overlaps < 9nt are filtered out (-TEov option)
	
	--------------------------------
	   FILES CREATED DURING THE RUN
	--------------------------------
	- creates \"Exon Structure\" file, with a line per exon with some info added like exon type (FIRST, LAST, etc), transcript coordinates and mature transcript length
	- extract features:
		exons (split non coding and coding -> split in CDS and UTRs)
		introns
		XXkb up and down (depending what is set - ex. 10, 5 and 1)
	- for introns and up+dw, <what-to-subtract.gtf or bed> will be subtracted from the sets => outputs are <file>.subtract.<name>.bed
	- then all resulting files are intersected with TEs (RMout.out)

	This script relies on Transcript information, so gene info is not mandatory. 
	However, it is better for analysis - note that GTF from the UCSC table browser use the same ID for gene and transcript

	When CDS, in ExSt file, exon nb will be for CDS exons only
	ExSt also contains features that would be filtered out. 
	It is basically generated once per input file, so you should delete it if the input file has been modified.

	UTRs won't be in exonSt file, just as the original \"exon\"
	UTR output - ALL UTRs => some can be undetermined for 5' or 3' so it might be relevant to check that output

	--------------------------------
	   READ THE OUTPUTS
	--------------------------------
	Typically, outputs can be used to generate figures and supp tables as in Kapusta et al. (2013) PLoS Genetics. 
	Supp tables (~TrInfos files) are a bit different because were generated by an older version of this pipeline.	
	See also: piRNA paper, Vinny's paper [for the bed format]
	Relevant outputs (and not intermediary files) are as follow:
	
	## ExSt file
	----------------
	   [ONLY for -f gtf => transcript analysis]
	   File used internally, summarizing all info about each exons (but can be useful, to get the whole transcript structure)
	
	## TrInfos files 
	----------------
	   [ONLY for -f gtf => transcript analysis]
	   Very useful table with all info of TE content per transcripts. One line per TE overlapping with a feature (exon, up and down...)
	   Note: Forget the intron TrInfos file... \"SPL\" may reveal overlap with a boundary of an intron piece after subtraction. 
	         Need to fix that, but it's not major. Not used to count features in the _Summary file anyway
	
	## _Summary.tab
	----------------
	   [ONLY for -f gtf => transcript analysis]
	   Summarizes number of features (TSS, polyA etc) overlapping with TEs       
	
	## CAT.tab
	----------------
	   Contains info of amount and class of different TE super families. Good quick way to look at global TE composition
	   % are to get proportions. They are not % regarding total amount of nt in the data (divide by the value at #total_length(nt) to get that)
	   When -f gtf (complex analysis) then individual files are concatenated in *concat.CAT.tab
	
	## CAT-class.tab
	----------------
	   Same as CAT.tab but by TE classes
	   This file also contains the total % of TEs in the data.
	   When -f gtf (complex analysis) then individual files are concatenated in *concat.CAT-class.tab
	
	## AGE.tab
	----------------
	   Contains info of amounts, but by age and not by class
	   When -f gtf (complex analysis) then individual files are concatenated in *concat.AGE.tab
	
	## TE-ratios.tab
	----------------
	# To check if there is any over representation of TEs
	# These do not have stats in them; significance should be tested using the \"nrCounts\" columns
	# Column details are as follow:
		Rname = repeat name from the repeat masker output
		Rclass = class
		Rfam = family
	IN SET:
		Len_masked = total length in nt covered by this repeat in the set
		%_masked = (Len_masked) / (total amount of TE in the set in nt) *100
		Counts = count of all fragments
		nrCounts = number of fragments corrected using the Repeat Masker Interrupted Repeats track (eg if there is one deletion or one insertion in a TE, 2 fragments but 1 corrected fragment)
		%nrCounts = (nrCounts) / (total amount of nrCounts in the set) *100
		NbChrs = number of different chromosomes / scaffolds where these counts are located.
		         Indeed, when 2 fragments are on 2 different sequences, less likely to be 2 pieces of the same TE... 
		         (nrCount is not always accurate and needs to be double checked mostly for low numbers, by looking at the TEjoin file.
		         Indeed, lowering numbers may affect significance)
	IN GENOME:
		Len_masked = total length in nt covered by this repeat in the genome
		%_masked = (Len_masked) / (total amount, in the genome, that is masked by all TEs that are also overlapping the peaks) *100
		Counts = count of all fragments
		nrCounts = number of fragments corrected using the Repeat Masker Interrupted Repeats track	
		%nrCounts = (nrCounts) / (total amount of nrCounts in the genome for these TEs) *100	
	AGE_INFORMATION: [if provided]
		Lineage = lineage, mostly from http://www.repeatmasker.org/cgi-bin/ViewRepeat?id=XXX (where XXX = repeat name), but also some personal checks
		Ancient/LineageSpe = Eutherian shared / lineage specific after the split
		avg_%div = pondered average of all % divergence from repeat masker (all fragments)	
	RATIOS:
		Len = ratio that tells you if a TE is potentially over represented
		nb = ratio based on nrCounts = used for stats but better to relie on RATION Len for figures if any
";
	return 1;
}	
		
################################################################################
# Get arguments/options, check some of them
################################################################################
my ($ft,$myf,$addcol,$fa,$bedtools) = ("gtf","na","na","na","na");
my ($subtract,$subid,$selfsub) = ("na","sub","yes"); 
my ($RMbase,$TEage,$nonTE,$TEov) = (1,"na","none",10);
my ($filter,$parse) = ("all,all","all,all");
my $cut = "10000,5000,1000";
my ($INPUT,$dir,$RMout,$RMPARSED,$TE,$clean,$chlog,$h,$help,$v);
GetOptions ('i=s'        => \$INPUT, 
            'dir'        => \$dir, 
            'f=s'        => \$ft, 
            'myf=s'      => \$myf, 
            'RMout=s'    => \$RMout, 
            'TEov=s'     => \$TEov, 
            'base=s'     => \$RMbase, 
            'RMparsed=s' => \$RMPARSED, 
            'TE=s'       => \$TE, 
            'TEage'      => \$TEage, 
            'nonTE=s'    => \$nonTE, 
            'fa=s'       => \$fa, 
            'subtract=s' => \$subtract, 
            'subid=s'    => \$subid, 
            'noselfsub'  => \$selfsub, 
            'bedtools=s' => \$bedtools, 
            'filter=s'   => \$filter, 
            'parse=s'    => \$parse, 
            'addcol=s'   => \$addcol, 
            'cut=s'      => \$cut, 
            'clean'      => \$clean, 
            'chlog'      => \$chlog, 
            'h'          => \$h, 
            'help'       => \$help, 
            'v'          => \$v);

#check step to see if mandatory arguments are provided + if help/changelog
die "\n version $version\n\n" if (! $INPUT && ! $RMout && ! $h && ! $help && ! $chlog && $v);
die $changelog if ($chlog);
die $longhelp if ($help);
die $usage if ($myf eq "na" && ! $clean && (! $INPUT || ! $RMout || $h));

#print sample file for column informations for custom format, if -myf set without mandatory arguments
myformat($myf) if (! $INPUT && ! $RMout && $myf ne "na");

die "\n   ERROR (main): input file $INPUT does not exist?\n\n" unless (-e $INPUT); 
die "\n   ERROR (main): RMout $RMout does not exist?\n\n" unless (-e $RMout); 

#get the list of input files + clean if needed
my @listin = ();
my %realins = ();
if ($dir) {	
	my $core;
	$INPUT = $1 if ($INPUT =~ /^(.*)\/$/); #avoid the / at the end of path
	my @list = `ls $INPUT`;	
	foreach my $in (@list) {
		chomp ($in);
		my $full = $INPUT."/".$in;
		push(@listin,$full);
		#need to get core names in case cleaning required
		if ($clean) {
			$core = $1 if ($in =~ /^(.*)\.TEjoin\./);
			$realins{$INPUT."/".$core}=1;
		}
	}
} else {
	push(@listin,$INPUT);
	$realins{$INPUT}=1 if ($clean);	
}
clean_out(\%realins) if ($clean);

#check some more options
my ($fcol,$fname) = split (",",$filter);
my ($pcol,$pname) = split (",",$parse);
#check -f option values
die "\n   ERROR (main): check -f option (use -h if you need to see the usage)\n\n" if ($ft ne "gtf" && $ft ne "bed"); #even if not set it will be "gtf"
#check filter and parse options
die "\n   ERROR (main): -filter requires 2 values separated by a coma (-filter <col,filter>; use -h if you need to see the usage)\n\n" if ($filter ne "all,all" && $filter !~ /,/);
die "\n   ERROR (main): -filter col should be numeric is -myf is chosen (use -h if you need to see the usage)\n\n" if ($filter ne "all,all" && $fcol !~ /\d/ && $myf ne "na");
die "\n   ERROR (main): -filter col should be numeric < 7 or non digit (use -h if you need to see the usage)\n\n" if ($filter ne "all,all" && ($fcol !~ /[0-7]/ || $fcol =~ /\D/) && $myf ne "na");
die "\n   ERROR (main): -filter intragenic can't be chosen if -f bed is chosen (use -h if you need to see the usage)\n\n" if (($fname eq "intragenic") && ($myf ne "na") && ($ft ne "gtf"));
die "\n   ERROR (main): -parse requires 2 values separated by a coma (use -h if you need to see the usage)\n\n" if (($parse ne "all,all") && ($parse !~ /,/));
die "\n   ERROR (main): use of -parse ($parse) without any columns added with -addcol (use -h if you need to see the usage)\n\n" if (($addcol eq "na") && ($parse ne "all,all"));
#if relevant, check that $addcol and $cut are numerical and , only
die "\n   ERROR (main): check -addcol option (use -h if you need to see the usage)\n\n" if ($addcol ne "na" && $addcol !~ /[,\d]/);
die "\n   ERROR (main): check -cut option (use -h if you need to see the usage)\n\n" if ($cut !~ /^[0-9,]+$/);
#if relevant, check extension of the file to subtract
die "\n   ERROR (main): check -subtract, file is not .gtf or .bed - if it is, please add the correct extension\n\n" if ($subtract ne "na" && $subtract !~ /.*\.bed$/ && $subtract !~ /.*\.gtf$/);

#Now start everything => log in STDERR if -v basically
if ($v) {
	print STDERR "\n ------------------------------------------------------------------------------------------------\n";
	print STDERR   " --- Script for TE analysis started (v$version)\n";
	print STDERR   "      - Input file = $INPUT\n" unless ($dir);
	print STDERR "        - Input files are located in $INPUT\n" if ($dir);
	print STDERR   "        -> Format / Analyze type = $ft\n";
	print STDERR   "        -> Custom columns will be extracted from $myf\n" if ($myf ne "na");
	print STDERR   "        -> Filtering of input file = $filter\n" if ($filter ne "all,all");	
	print STDERR   "        -> These columns will be added in the output files = $addcol\n" if ($addcol ne "na");
	print STDERR   "        -> Filtering on these added columns = $parse\n" if ($parse ne "all,all");
	print STDERR   "      - Repeat Masker output file .out = $RMout\n";
	print STDERR   "        -> base = $RMbase\n";
    print STDERR   "        -> minimal length of intersection = $TEov nt\n" if (substr($TEov,-1) ne "%");
    print STDERR   "        -> minimal overlap of feature and TE = $TEov of the feature\n" if (substr($TEov,-1) eq "%");
	print STDERR   "        -> parsed file = $RMPARSED (over/under represention of TE families will be determined)\n" if ($RMPARSED);
	print STDERR   "        -> TE file = $TE\n" if ($TE);
	print STDERR   "        -> lineage/age info will be added in outputs + TE amounts will be parsed by age\n" if ($TEage ne "na");
	print STDERR   "        -> nonTE repeats will be filtered as -nonTE $nonTE (see usage or full documentation)\n";
	if ($ft eq "gtf") {
		if ($fa ne "na") {
			print STDERR   "      - Genome file (fasta) = $fa\n";
			print STDERR   "        -> up and downstream regions will be analyzed\n";
			print STDERR   "        -> of lengths = $cut nt\n";
		} elsif ($fa eq "na") { 
			print STDERR   "      - No genome file (-fa option) provided\n";
			print STDERR   "        -> up and downstream regions won't be analyzed\n";
		}		
		print STDERR   "      - $subtract will be subtracted from introns and intergenic regions\n" if ($subtract ne "na");
		($selfsub eq "yes")?(print STDERR "        -> option -noselfsub not chosen  => $INPUT will be subtracted from introns and surrounding regions as well\n"):(print STDERR "        -> option -noselfsub chosen  => $INPUT won't be subtracted from introns and surrounding regions\n");
		print STDERR   "        -> subtract ID = $subid (will be in file names)\n" if ($subtract ne "na");	
	}
	$bedtools = $1 if (($bedtools ne "na") && ($bedtools =~ /^(.*)\/$/)); #avoid the / at the end of pathif ($bedtools ne "na");
	print STDERR   "        -> BEDtools path = $bedtools\n" if ($bedtools ne "na");
    print STDERR " ------------------------------------------------------------------------------------------------\n";
}


#Get TE infos if relevant
my $TEinfo = ();
my $TEinfoRMP = ();
my $TE_RMP = ();
if ($TE) {
	$TEinfo = get_TEs_infos($TE,$v);
} else {
	($TEinfo->{"na"} = "na");
}
if ($RMPARSED) {
	($TEinfoRMP,$TE_RMP) = get_TEs_infos($RMPARSED,$v);
} else {
	$TEinfoRMP->{"na"} = "na",$TE_RMP->{"na"} = "na";
}

#Then write TEs in bed format if needed => do some filtering out (nonTE stuff)
$RMout = RMtobed($RMout,$RMbase,$TEinfoRMP,$TEinfo,$nonTE,$v);

#Convert to bed the file to be subtracted unless it was given in bed format
$subtract = gtftobed($subtract,$v) if ($subtract ne "na");

#Get lengths from genome file if relevant
my $lengths = get_lengths_noBio($fa,$v) if ($fa ne "na");

#Get what column will be what - variable depending of $ft value
my $set = set($ft,$myf,$v);


#extract all relevant info from the input file; if -f gtf => make Exon Structure file + get up/dw etc
#####################################################
foreach my $in (@listin) {
	print "\n --- Dealing with $in " if (($v) && ($dir));
	my $addn = "$fname.$pname";
	my $listtojoin = ();
	my $listoffiles = ();
	my $name = $in;
	$name = $1 if ($in =~ /(.*)\.gtf|gff|gff3$/);
	my $ExSt = "$name.ExSt.$addn.tab";
	my $TrInfos = ();
	my $feat = ();
	my $ExInfo = ();
	if ($ft eq "gtf") {
		my ($pc,$big);
		#Load input file and get all info + filtering value
		my $TrMatLen = ();
		$TrMatLen->{"na"} = "na"; #initialize a value to avoid messing with arguments for get_tr_infos
		$TrMatLen = get_TrMatLen($in,$v) if ($fname eq "intragenic"); #only possible if $myf ne "na" (option checked at the beginning of the script)
		($pc,$big,$TrInfos,$listtojoin) = get_tr_infos($in,$set,$name,$filter,$addcol,$parse,$addn,$TrMatLen,$v);		
		push(@{$listoffiles},"$name.nc.$addn.introns") if (-e "$name.nc.$addn.introns.bed"); 
		push(@{$listoffiles},"$name.pc.$addn.introns") if (-e "$name.pc.$addn.introns.bed"); 
		push(@{$listoffiles},"$name.pc.$addn.CDS.introns") if (-e "$name.pc.$addn.CDS.introns.bed"); 

		#print/load exon structure file => -filter and -parse will be done there
		unless (-e $ExSt) {
			#BIG TABLE STRUCTURE: ($chr,$type,$feat,$currstart,$currend,$strand,$gene_id,$tr_id,$gene_name,$tr_name,$gene_type,$tr_type,$colstoadd);	
			my @big = @{$big};
			@big = sort {
				# by chr
				($a -> [0] cmp $b -> [0]) 
				# by gene
				|| ($a -> [6] cmp $b -> [6])
				# by transcript 
				|| ($a -> [7] cmp $b -> [7])
				# by coordinates
				|| ($a -> [3] <=> $b -> [3])
				|| ($a -> [4] <=> $b -> [4])
			} @big;	
			#Now load the big table and print ExSt file
			($feat,$ExInfo,$TrInfos) = print_ExSt($ExSt,$TrInfos,\@big,$pc,$addcol,$v);		
			#Done with @big => undef it since it's big, no need to use the memory
			undef @big;	
		} else {
			print STDERR   "\n --- Exon Structure file seems to exist\n" if ($v);
			($feat,$ExInfo,$TrInfos) = load_ExSt($ExSt,$TrInfos,$v);
		}		

		#Get UTRs, only if there were coding stuff
		if (-e "$name.pc.$addn.CDS.bed") {
			get_UTRs($bedtools,$name,$pc,$TrInfos,$addn,$v);
			push(@{$listtojoin},("$name.pc.$addn.UTRs","$name.pc.$addn.5UTRs","$name.pc.$addn.3UTRs"));
		}
	
		#Get up and down regions, unless no genome file
		if ($fa ne "na") {
			my @cut = split(",",$cut);
			my $len = $cut[0];
			get_up_dw($name,$len,$TrInfos,$lengths,$fa,$addn,$v); #unless ((-e "$name.nc.$len.up.bed") && (-e "$name.nc.$len.dw.bed"));
			push(@{$listoffiles},("$name.nc.$addn.up","$name.nc.$addn.dw"));
			push(@{$listoffiles},("$name.pc.$addn.up","$name.pc.$addn.dw")) if (-e "$name.pc.$addn.CDS.bed");
		}
	
		#Now subtract stuff if relevant
		if (($subtract ne "na") || ($selfsub eq "yes")) {
			$listtojoin = subtract($listtojoin,$listoffiles,$subtract,$subid,$name,$in,$selfsub,$bedtools,$TrInfos,$addn,$cut,$v);
		} else {
			push(@{$listtojoin},@{$listoffiles});
		}	
	} else {
		#Bed file - just join it
		$name = $1 if ($in =~ /(.*)\.bed/);
		push(@{$listtojoin},$name);
	}
	
	#Join files with TEs [longest step]
	#####################################################
	my $listtoanalyse = TE_join($listtojoin,$RMout,$bedtools,$v); #first list that has the extension of the files in it

	#Parse the joined files with TEs
	#####################################################
	my $feat_nb = `wc -l $in`;
	$feat_nb =~ s/^\s*([0-9]+?)\s+.*$/$1/;
	$feat->{'all'}{'all'}{$in}=$feat_nb; #save number of features by file for bed ft
	my $totlen = get_amounts($listtojoin,$name,$addn,$TrInfos,$cut,$ft,$v); #using files before joining
	my $countTE = parse_join($listtoanalyse,$TrInfos,$ExInfo,$TEinfoRMP,$TE_RMP,$TEinfo,$TEage,$addcol,$cut,$totlen,$TEov,$ft,$v);

	#Summary file
	summary($in,$countTE,$feat,$addn,$ft,$v);	
}
#30 subs later - Done!
print STDERR   "\n --- Script for TE analysis is done\n" if ($v);
print STDERR " ------------------------------------------------------------------------------------------------\n\n" if ($v);
exit;




##########################################################################################################
# SUBROUTINES
##########################################################################################################
#----------------------------------------------------------------------------
# clean previous outputs generated by this script
# clean_out($realins);
#----------------------------------------------------------------------------
sub clean_out {
	my $reals = shift;
	print STDERR "\n --- Deleting previous outputfiles (except RMout.bed and <fa>.fa.lengths):\n";
	foreach my $in (keys %{$reals}) { 
		print STDERR "    For inputfile = $in\n";

		$in = $1 if ($in =~ /(.*)\.gtf/);
		$in = $1 if ($in =~ /(.*)\.bed/);
		
		print STDERR "      - $in.*TEjoin.*\n";
		`rm -Rf $in.*TEjoin.*`;
		print STDERR "      - $in.*.AGE.tab\n";
		`rm -Rf $in.*.AGE.tab`;
		print STDERR "      - $in.*.CAT.tab\n";
		`rm -Rf $in.*.CAT.tab`;
		print STDERR "      - $in.*.CAT-class.tab\n";
		`rm -Rf $in.*.CAT-class.tab`;
		print STDERR "      - $in.*.TEs-ratios.tab\n";
		`rm -Rf $in.*.TEs-ratios.tab`;
		print STDERR "      - $in.*.amounts.txt\n";
		`rm -Rf $in.*.amounts.txt`;	
		print STDERR "      - $in.*.concat.*.tab\n";
		`rm -Rf $in.*.concat.*.tab`;
		print STDERR "      - $in.*_Summary.tab\n";
		`rm -Rf $in.*_Summary.tab`;	
		print STDERR "      - $in.pc.*\n";
		`rm -Rf $in.pc.*`;
		print STDERR "      - $in.nc.*\n";
		`rm -Rf $in.nc.*`;
		print STDERR "      - $in.ExSt.*.tab\n";
		`rm -Rf $in.ExSt.*.tab`;	
		print STDERR "      - $in.nc-pc*\n";
		`rm -Rf $in.nc-pc*`;
		print STDERR "\n";
	}
	exit;	
}

#----------------------------------------------------------------------------
# subroutine to print a sample file with the column details and exit
# myformat($myf) if ((! $in) && (! $RMout) && ($myf);
# called by main
#----------------------------------------------------------------------------
sub myformat {
	my $myf = shift;
	open(my $myf_fh, ">$myf") or confess "\n   ERROR (sub myformat): could not open to write $myf $!\n";
	print $myf_fh "###################################
# Configuration file for the script TE-analysis_gtf_pipeline_ak.pl
# Allow use of a custom input file -> set columns
###################################
# Numbers correspond to column numbers; note that first column = 0
# If -bed is chosen: type, feat, gene_id, gene_name and tr_name won't be considered, and tr_id will correspond to the unique ID.
# Spaces don't matter, they will be removed anyway
# Don't remove the # between numbers and comments
# If you run on a server, you can edit this file with emacs, vim, etc
###################################

chr             = 0   #chromosome or scaffold name
type            = 1   #havana, cufflink...
feat            = 2   #feature, e.g exon, CDS, start,stop, mRNA or gene... Note that mRNA and gene lines will be ignored necause not required for the pipeline to run
start           = 3   #start of the feature
end             = 4   #end of the feature
strand          = 5   #strand of the feature
gene_id         = 6   #gene ID
gene_name       = 6   #gene name; if none just use same as gene_id
transcript_id   = 7   #transcript_ID (use this as the unique ID if -f bed is chosen)
transcript_name = 7   #transcript name; if none, just use same as tr_id
gene_type       = 8   #gene biotype (e.g. protein_coding, lincRNA, snRNA, processed_transcript, etc) - required to split coding and non coding elements.
transcript_type = 8   #can be same as gene biotype
\n\n";
	close ($myf_fh);
	exit;	
}

#----------------------------------------------------------------------------
# get TE infos; for -TE or for -RMparsed
# ($TE)?($TEinfo = get_TEs_infos($TE,$v)):($TEinfo->{"na"} = "na");
# ($RMPARSED)?(($TEinfoRMP,$TE_RMP) = get_TEs_infos($RMPARSED,$v)):($TEinfoRMP->{"na"} = "na",$TE_RMP->{"na"} = "na");
# called by main
#----------------------------------------------------------------------------
sub get_TEs_infos {
	my ($in,$v) = @_;
	print STDERR " --- Loading TE info from $in\n" if ($v);
	my %TEs = ();
	my %TE_RMP = ();
	open(my $in_fh, "<", $in) or confess "\n   ERROR (sub get_TEs_infos): could not open to read $in $!\n";
	my $i = 0;
	my $r = 0;
	my $ifn = "y";
	LINE: while(<$in_fh>) {
		chomp (my $line = $_);
		if ($i == 0 && $in =~ /.out.parseRM.*.tab$/ && $line =~ /Rfullname/ && $line =~ /MED_LEN_MASKED/) {
			$r = 1 if ($line =~ /Rlen/);
			$ifn = "n";
		}	
		$i++;
		next LINE if ($line !~ /\w/ || substr($line,0,5) eq "Rname" || substr($line,0,1) eq "#");		
		my @TEs = split('\t', $line); 
		
		#Deal with name, class, fam:
		my ($Rn,$Rc,$Rf) = ($TEs[0],$TEs[1],$TEs[2]);		
		#make sure no / in family
		$Rf =~ s/\//_/;
		#add the -int if not there
		$Rf = $Rf."--int" if ($Rf =~ /ERV/ && ($Rn =~ /[-_][iI]nt/ || $Rn =~ /[-_]I$/));
		my $Rcf = $Rc."/".$Rf;
		
		if ($in =~ /.parseRM.*.tab$/) {
			#Now the rest of the columns will differ:	
			#columns in new parseRM.pl:
			# 3		4			5								6				7			8
			# Rlen 	FRG_NB_all	FRG_NB_Reconstructed_repeats	LEN_MASKED_NR	AVG_%DIV	MED_%DIV	
			# 9			10			11			12			13				14
			# AVG_%DEL	MED_%DEL	AVG_%INS	MED_%INS	AVG_LEN_MASKED	%_GENOME
			# 15			16					17
			# LEN_OVERLAP	%_OVERLAP_(GENOME)	%_OVERLAP_(LEN_MASKED)
		
			#columns in old parseRM.pl (Rlen not necessarily there):
			# 5			6					7			8			9
			# FRG_NB	FRG_NB_StartToEnd	NR_FRG_NB	AVG_%DIV	MED_%DIV	
			# 10			11		12			13			14			15				16				17
			# AVG_%DEL	MED_%DEL	AVG_%INS	MED_%INS	LEN_MASKED	AVG_LEN_MASKED	MED_LEN_MASKED	%_GENOME	
			# 18			19					20
			# LEN_OVERLAP	%_OVERLAP_(GENOME)	%_OVERLAP_(LEN_MASKED)
			
			my ($frg,$frgnr,$ad,$len,$pgm,$leno,$pgmo);
			if ($ifn eq "y") {
				($frg,$frgnr,$ad) = ($TEs[4],$TEs[5],$TEs[7]);
				($len,$pgm) = ($TEs[6],$TEs[14]);
				if ($TEs[15] && $TEs[16]) {
					($leno,$pgmo) = ($TEs[15],$TEs[16]);
				} else {
					($leno,$pgmo) = (0,0);
				}
			} else {
				($frg,$frgnr,$ad) = ($TEs[5+$r],$TEs[7+$r],$TEs[8+$r]);
				($len,$pgm) = ($TEs[14+$r],$TEs[17+$r]);
				if ($TEs[18+$r] && $TEs[19+$r]) {
					($leno,$pgmo) = ($TEs[18+$r],$TEs[19+$r]);
				} else {
					($leno,$pgmo) = (0,0);
				}
			}

			if ($pgm eq "nd") {
				print STDERR "     /!\\ The % of the genome covered by each repeat is missing in $in\n";
				print STDERR "          Please rerun parseRM.pl with the -f option\n";
				print STDERR "     ... exiting\n\n";
				exit;
			}
			
			#edit the list	
			#class famm class/fam frg nr_frg avg%div len_masked %genome_masked len_overlap %genome_overlap(for this repeat) 
			@TEs = ($Rc,$Rf,$Rcf,$frg,$frgnr,$ad,$len,$pgm,$leno,$pgmo);			
					
			#load TE_RMP hash		
			$TE_RMP{'l'}{lc($Rn)}+=($len-$leno);
			$TE_RMP{'p'}{lc($Rn)}+=($pgm-$pgmo);
			$TE_RMP{'cnr'}{lc($Rn)}+=$frgnr;
			$TE_RMP{'ctot'}{lc($Rn)}+=$frg;

		}
		#now load TE list
		$TEs{lc($Rn)} = \@TEs;	#??? The only thing I use from this is the %div?? Double check, and store that only if yes!
	}	
	close $in_fh;
	if ($in =~ /.parseRM.*.tab$/) {
		return (\%TEs,\%TE_RMP);		
	} else {
		return \%TEs;
	}
}

#----------------------------------------------------------------------------
# Convert RMoutput .out file to bed + filter nonTE + update class / fam if relevant
# ($RMout) = RMtobed($RMout,$RMbase,$TEinfoRMP,$TEinfo,$nonTE,$v);
# called by main
#----------------------------------------------------------------------------
sub RMtobed {
	my ($RMout,$base,$TEinfoRMP,$TEinfo,$nonTE,$v) = @_;
	my $bed = $RMout;
	$bed =~ s/(.*)\.out$/$1/;
	if (! $TEinfo->{"na"} || ! $TEinfoRMP->{"na"}) {
		$bed = $bed.".class.nonTE-$nonTE.bed"; 
	} else {
		$bed = $bed.".nonTE-$nonTE.bed";	
	}
	print STDERR " --- Converting $RMout to bed\n" if ($v);
	unless (-e $bed) {
		print STDERR "     Filtering nonTE repeats based on -nonTE: $nonTE\n" if ($v);
		print STDERR "     Updating TE class and family based on -TE\n" if ($v && ! $TEinfo->{"na"});
		print STDERR "     Updating TE class and family (if relevant) based on -RMparsed\n" if (($v) && (! $TEinfoRMP->{"na"}) && ($TEinfo->{"na"}));
		print STDERR "     Updating TE class and family (if relevant) based on 1) -TE and 2) -RMparsed\n" if (($v) && (! $TEinfo->{"na"}) && (! $TEinfoRMP->{"na"}));
		open(my $fh, "<$RMout") or confess "\n   ERROR (sub RMtobed): could not open to read $RMout!\n";
		open(my $bed_fh, ">$bed") or confess "\n   ERROR (sub RMtobed): could not open to write $bed!\n";
		LINE: while(<$fh>) {
			chomp(my $l = $_);
			$l =~ s/^\s+//;
			next LINE if (($l =~ /^[Ss]core|^SW|^#/) || ($l !~ /\w/));
			$l = $1 if ($l =~ /^(.*)\*$/); #remove the star
			my @l = split('\s+',$l);
			$l[8] =~ s/C/-/; #correct strand
			
			#now get the TE info hash, unless already defined for this element.
			#in older libraries there were repeats with same name but different classfam - well if it happens then the first occurence will be the one in the hash
			my ($Rclass,$Rfam,$Rclassfam) = get_Rclass_Rfam($l[10],$l[9]);
			my $rname = lc($l[9]);
						
			#update class and family if relevant (will have the additional --int)
			if ($TEinfo->{$rname}) {
				($Rclass,$Rfam,$Rclassfam) = ($TEinfo->{$rname}[0],$TEinfo->{$rname}[1],$TEinfo->{$rname}[2]);	
			} elsif ($TEinfoRMP->{$rname}) { #use this only if not in $TEinfo
				($Rclass,$Rfam,$Rclassfam) = ($TEinfoRMP->{$rname}[0],$TEinfoRMP->{$rname}[1],$TEinfoRMP->{$rname}[2]);
			}		
			
			#now filter non TE or not, based on -nonTE value
			next LINE if (($nonTE eq "none") && ($Rclass =~ /nonTE|snRNA|rRNA|tRNA|snoRNA|scRNA|srpRNA|[Ll]ow_complexity|[Ss]imple_repeat|[Ss]atellite|ARTEFACT/));
			next LINE if (($nonTE eq "no_nonTE") && ($Rclass =~ /nonTE/));
			next LINE if (($nonTE eq "no_low") && ($Rclass =~ /[Ll]ow_complexity|[Ss]imple_repeat/));
			
			#now create unique ID + Rclassfam will be updated (in cased changed)
			my ($chr,$start,$end,$strand) = ($l[4],$l[5],$l[6],$l[8]);
			my $ID = $l[0];
			for (my $i=1; $i<=9;$i++) {
				$ID = $ID.";".$l[$i];
			}
			$ID = $ID.";".$Rclassfam;
			for (my $i=11; $i<=$#l;$i++) {
				$ID = $ID.";".$l[$i];
			}
			$ID =~ s/\s//; #should not need this since it would come from error in TEinfo input file for ex, but somebody could still have the issue
			
			#correct the start if base is 0
			$start=$start+1 if ($base == 0);
			
			#now print 
# 			$chr =~ s/gi\|.+\|gb/gb/ if ($gb == 1); #For Rachel
			print $bed_fh "$chr\t$start\t$end\t$ID\t.\t$strand\n"; #with ID being the whole line => easy to acces to RMoutput
		}
		close ($fh);
		close ($bed_fh);
		print STDERR "     => $bed\n" if ($v);
	} else {
		print STDERR "       $bed exists, skipping (delete if -TE or -parsedRM were changed)\n" if ($v);
	}
	return ($bed);
}

#----------------------------------------------------------------------------
# Get Rclassfam from RMout
# my ($Rclass,$Rfam,$Rclassfam) = get_Rclass_Rfam($classfam,$Rname);
# called by RMtobed
#----------------------------------------------------------------------------
sub get_Rclass_Rfam {
	my($classfam,$Rname) = @_;
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
	my $Rclassfam = "$Rclass/$Rfam";
	return ($Rclass,$Rfam,$Rclassfam);
}

#----------------------------------------------------------------------------
# Convert gtf file to bed
# gtftobed($subtract,$v) if (($subtract) && ($subtract !~ /.*\.bed$/));
# called by main
#----------------------------------------------------------------------------
sub gtftobed {
	my ($gtf,$v) = @_;
	my $bed;
	if ($gtf =~ /(.*)\.gtf/) {
		$bed = $1.".bed";
	} elsif ($gtf =~ /.*\.bed/){
		$bed = $gtf;
		return ($bed);
	}	
	print STDERR "\n --- Converting $gtf to bed\n" if ($v);
	unless (-e $bed) {	
		open(my $gtf_fh, "<$gtf") or confess "\n   ERROR (sub gtftobed): could not open to read $gtf!\n";
		open(my $bed_fh, ">$bed") or confess "\n   ERROR (sub gtftobed): could not open to write $bed!\n";
		LINE: while(<$gtf_fh>) {
			chomp(my $line = $_);
			$line =~ s/\s+"/"/g;
			$line =~ s/";\s+/";/g;
			next LINE if (substr($line,0,1) eq "#");
			next LINE if (substr($line,0,2) eq "MT"); #get rid of mitochondrial DNA
			my @line = split (/\s+/,$line);		
			my ($chr,$feat,$start,$end,$strand) = ($line[0],$line[2],$line[3],$line[4],$line[6]);
			if ($feat eq "exon") {
				#correct for the fact that there are only numbers for chromosomes in ensembl gtf (not gencode for ex)
				$chr ="chr".$chr if ($chr =~ /^\d+$/); 
				#Get strand +/- (Note some stuff might be a * => no strand info, e.g. single exons lncRNAs...)
				$strand = "+" if ($strand eq "1");
				$strand = "-" if ($strand eq "-1");
				#now, get gene and transcript ids + gene and tr biotypes - different depending on the gtf type
				my ($tr_id);
				my @info = split(';',$line[8]);
				for (my $i = 0; $i < $#info; $i++) {
					my ($id,$value);
					($info[$i] =~ /"/)?(($id,$value) = split(/"/,$info[$i])):(($id,$value) = split(/\s+/,$info[$i]));
					$tr_id = $value if ($id eq "transcript_id");
				}
				my $ID = $tr_id."#".$chr.":".$start."-".$end;
				print $bed_fh "$chr\t$start\t$end\t$ID\t.\t$strand\n" unless ($end - $start <1);
			}
		}
		close ($gtf_fh);
		close ($bed_fh);
	} else {
		print STDERR "       $bed exists, skipping\n" if ($v);
	}	
	return ($bed);
}

#----------------------------------------------------------------------------
# Get lengths of all sequences and store that by sequence ID. Note that if some are not unique, it just replaces by last length.
# This sub does not use BioPerl - avoid having to index the genome
# my $lengths = get_lengths_noBio($fa,$v) if ($fa);
#----------------------------------------------------------------------------
sub get_lengths_noBio {
	my ($fa,$v) = @_;
	
	print STDERR "\n --- Getting sequences lengths for $fa\n" if ($v);
	my %len = ();	
	my $lengthfile = "$fa.lengths";	
	
	if (-e $lengthfile) {
		print STDERR "     -> lengths have been previously calculated ($lengthfile exists) => extracting\n" if ($v);
		#extract lengths now
		open (my $lengths_fh, "<", $lengthfile) or confess "   \nERROR (sub get_lengths): could not open $lengthfile $!\n";
		while (<$lengths_fh>) {
			chomp (my $line = $_);
			my ($id,$len) = split(/\s+/,$line);
			$len{$id}=$len;
		}	
		close ($lengths_fh);
	} else {
		#looping through fasta file
		my $id = "";
		my $l = 0;
		my $c = 0;
		open (my $fa_fh, "<", $fa) or confess "   \nERROR (sub get_lengths): could not open $fa $!\n";
		open (my $len_fh, ">", $lengthfile) or warn "   \nERROR (sub get_lengths): could not create $lengthfile, but lengths will be calculated $!\n";
		while (<$fa_fh>) {
			chomp (my $line = $_);
			if (substr($line,0,1) eq ">") {
				#first get and print unless first header
				unless ($c == 0) {
					print $len_fh "$id\t$l\n";
					$len{$id}=$l;
				}
				$c=1;
				#store header and reinitialize length
				my @id = split (/\s+/,$line);
				$id = $id[0];
				$id =~ s/>//;
				$l = 0;
			} else {
				#get length; could be more than one line so increment
				$l+=length($line);
			}
		}
		#get and print len last sequence
		print $len_fh "$id\t$l\n";
		$len{$id}=$l;
		
		close ($fa_fh);
		close ($len_fh);
	}
	return (\%len);
}

#----------------------------------------------------------------------------
# subroutine to set format
# my $set = set($ft,$myf,$v);
# called by main
#----------------------------------------------------------------------------
sub set {
	my ($ft,$myf,$v) = @_;
	my %set = ();
	if ($myf ne "na") { #then load custom file
	print STDERR "\n --- Reading custom column information from $myf\n" if ($v);
		open(my $myf_fh, "<", $myf) or confess "\n   ERROR (sub set): could not open $myf!\n";
		while(<$myf_fh>) {
			chomp (my $line = $_);
			if (($line !~ /^#/) && ($line =~ /=/)) { #deal only with lines that have a "=" and avoid commented full lines
				$line =~ s/\s+//g; #remove spaces
				my @line = split(/[=#]/,$line);
				$set{$line[0]}=$line[1]; #some values may be empty, if -bed
				confess "\n   ERROR (sub set): $line[1] is not a number? Check the file from -myf for formatting issues\n" if ($line[1] !~ /\d/);
			}	
		}	
		close $myf_fh;		
	} else {
		if ($ft eq "gtf"){
		#1       ensembl exon    183114  183240  .       +       .       gene_id "ENSG00000279928"; gene_version "1"; transcript_id "ENST00000624431"; transcript_version "1"; exon_number "2"; gene_name "FO538757.3"; gene_source "ensembl"; gene_biotype "protein_coding"; transcript_name "FO538757.3-201"; transcript_source "ensembl"; transcript_biotype "protein_coding"; exon_id "ENSE00003759581"; exon_version "1";
			$set{"chr"} = 0;
			$set{"type"} = 1;
			$set{"feat"} = 2;
			$set{"start"} = 3;
			$set{"end"} = 4;
			$set{"strand"} = 6;
			$set{"info"} = 8;
		} elsif ($ft eq "bed") {
			$set{"chr"} = 0;
			$set{"start"} = 1;
			$set{"end"} = 2;
			$set{"transcript_id"} = 3;
			$set{"strand"} = 5;
		}
	}
	return (\%set);
}

#----------------------------------------------------------------------------
# subroutine to extract gtf info line
# my $TrMatLen = get_TrMatLen($in,$v) if (($fname eq "intragenic") && ($myf ne "na") && ($ft eq "gtf"));
# called by main
#----------------------------------------------------------------------------
sub get_TrMatLen {
	my ($in,$v)	= @_;
	my $TrMatLen = ();
	print STDERR "\n --- Getting mature length of transcripts to allow intragenic filtering\n" if ($v);
	open(my $in_fh, "<", $in) or confess "\n   ERROR (sub get_tr_infos): could not open to read $in!\n";
	LINE: while(<$in_fh>) {
		chomp (my $line = $_);
		next LINE if (($line !~ /\w/) && ($line =~ /^MT/) && (substr($line,0,1) eq "#")); #non blank + get rid of mitochondrial DNA or commented lines
		$line =~ s/\s+"/"/g;
		$line =~ s/";\s+/";/g;
		my @l = split (/\s+/,$line);	
		my ($feat,$st,$en,$info) = ($l[2],$l[3],$l[4],$l[8]);
		next LINE unless (($feat eq "transcript") || ($feat eq "mRNA"));
		my @info = split(';',$info);
		my ($tr_id,$gene_id);
		for (my $i = 0; $i < $#info; $i++) {
			my ($id,$value);
			($info[$i] =~ /"/)?(($id,$value) = split(/"/,$info[$i])):(($id,$value) = split(/\s+/,$info[$i]));
			$tr_id = $value if ($id eq "transcript_id");
			$gene_id = $value if ($id eq "gene_id");
		}
		my $tr_id_full = $gene_id."#".$tr_id; 
		$TrMatLen->{$tr_id_full} = $en - $st +1;
	}	
	return ($TrMatLen);
}		

#----------------------------------------------------------------------------
# subroutine to extract gtf info line
# get_gtf_info($line,$set);
# called by get_tr_infos
#----------------------------------------------------------------------------
sub get_gtf_values {
	my ($line,$set,$filter,$addcol)	= @_;
	my %set = %{$set};
	$line =~ s/\s+"/"/g;
	$line =~ s/"; /";/g;
	my @line = split (/\t/,$line);		
	
	#Get easy values
	my ($chr,$type,$feat,$start,$end,$strand) = ($line["$set{'chr'}"],$line["$set{'type'}"],$line["$set{'feat'}"],$line["$set{'start'}"],$line["$set{'end'}"],$line["$set{'strand'}"]);

	#Return if it is mRNA, transcript or gene
	return ("na",$chr,$type,$feat) if (($feat eq "gene") || ($feat eq "transcript") || ($feat eq "mRNA"));
	
	# Deal with columns to add
	my @colstoadd = ();
	unless ($addcol eq "na") {
		my @addcols = split (",",$addcol);	
		foreach my $col (@addcols) {
			push (@colstoadd,$line[$col]);
		}
	}
		
	#correct for the fact that there are only numbers for chromosomes in ensembl gtf (not gencode for ex)
	$chr = "chr".$chr if (($chr =~ /^\d+$/) || ($chr =~ /^\w$/)); #one letter = X, Y 
	
	#Get strand +/- (Note some stuff might be a * => no strand info, e.g. single exons lncRNAs...)
	$line[$set{"strand"}] = "+" if ($strand eq "1");
	$line[$set{"strand"}] = "-" if ($strand eq "-1");
	
	#Deal with filtering
	my ($fcol,$fname) = split (",",$filter);
	my $filtering = "na";
	$filtering = $line[$fcol] if ($fcol =~ /\d/); #if -myf or if <7 then just column number, it was checked before
	
	#now, get gene and transcript ids + gene and tr biotypes - different if real gtf or not. Deal with filtering and columns to add as well.
	my ($gene_id,$tr_id,$gene_name,$tr_name,$gene_type,$tr_type);
	my %v = ('gene_id' => "na",
             'transcript_id' => "na",
             'gene_name' => "na",
             'transcript_name' => "na",
             'gene_type' => "na",
             'transcript_type' => "na");
	if ($set{"info"}) { #if that's defined, means it is gtf format and not the -myf
		my $info = $line[$set{"info"}];
		my @info = split(';',$info);
		for (my $i = 0; $i <= $#info; $i++) {
			my ($id,$value);
			if ($info[$i] =~ /"/) {
				($id,$value) = split(/"/,$info[$i]);
			} elsif ($info[$i] =~ /=/) {
				($id,$value) = split(/=/,$info[$i]);
			} else {
				($id,$value) = split(/\s+/,$info[$i]);
			}
			$id =~ s/_biotype$/_type/; #not same appelation in gencode and ensembl
			$v{$id}=$value;
		}
		($gene_id,$tr_id,$gene_name,$tr_name,$gene_type,$tr_type) = ($v{"gene_id"},$v{"transcript_id"},$v{"gene_name"},$v{"transcript_name"},$v{"gene_type"},$v{"transcript_type"});
		
		#get filtering
		if (($fcol !~ /[0-7]/) && ($fcol ne "na")) {
			FIL: foreach my $key (keys %v) {
				if ($key eq $fcol) {
					$filtering = $v{$key};
					last FIL;
				}		
			}
		}	
	} else { #it is a custom file => just load the values
		($gene_id,$tr_id,$gene_name,$tr_name,$gene_type,$tr_type) = ($line[$set{"gene_id"}],$line[$set{"transcript_id"}],$line[$set{"gene_name"}],$line[$set{"transcript_name"}],$line[$set{"gene_type"}],$line[$set{"transcript_type"}]);
	}	

#NB:	
	#ENSEMBL
	#gene_id "ENSG00000279457"; gene_version "1"; transcript_id "ENST00000623834"; transcript_version "1"; exon_number "6"; gene_name "FO538757.2"; gene_source "ensembl"; gene_biotype "protein_coding"; transcript_name "FO538757.2-201"; transcript_source "ensembl"; transcript_biotype "protein_coding"; exon_id "ENSE00003756131"; exon_version "1";
	
	#GENCODE	
	#gene_id "ENSG00000223972.5"; transcript_id "ENSG00000223972.5"; gene_type "transcribed_unprocessed_pseudogene"; gene_status "KNOWN"; gene_name "DDX11L1"; transcript_type "transcribed_unprocessed_pseudogene"; transcript_status "KNOWN"; transcript_name "DDX11L1"; level 2; havana_gene "OTTHUMG00000000961.2";

	#GENCODE v21
	#ID=exon:ENST00000473358.1:1;Parent=ENST00000473358.1;gene_id=ENSG00000243485.4;transcript_id=ENST00000473358.1; gene_type=lincRNA;gene_status=KNOWN;gene_name=MIR1302-2;transcript_type=lincRNA;transcript_status=KNOWN;transcript_name=MIR1302-2-001;exon_number=1;exon_id=ENSE00001947070.1;level=2;transcript_support_level=5;tag=not_best_in_genome_evidence,dotter_confirmed,basic;havana_gene=OTTHUMG00000000959.2;havana_transcript=OTTHUMT00000002840.1

	return ($filtering,$chr,$type,$feat,$start,$end,$strand,$gene_id,$tr_id,$gene_name,$tr_name,$gene_type,$tr_type,\@colstoadd);
}		

#----------------------------------------------------------------------------
# subroutine to extract all relevant info from the input file, if -f gtf
# ($pc,$big,$TrInfos,$listtojoin) = get_tr_infos($in,$set,$name,$filter,$addcol,$parse,$addn,$TrMatLen,$v);
# called by main
#----------------------------------------------------------------------------
sub get_tr_infos {
	my ($in,$set,$name,$filter,$addcol,$parse,$addn,$TrMatLen,$v) = @_;	
	my $ExSt = "$name.ExSt.$addn.tab";
	my %pc = ();
	my %TrInfos = ();
	my @big = ();
	my @listtojoin = ();
	print STDERR "\n --- Reading $in and extracting informations (getting exons, introns, transcripts info)...\n" if ($v);
	print STDERR "     Also filtering based on -filter ($filter)\n" if ($filter ne "all,all" && $v);
	print STDERR "     Also filtering based on -parse ($parse)\n" if ($parse ne "all,all" && $v);
	my %ifstuff = ();
	($ifstuff{'pc'},$ifstuff{'nci'},$ifstuff{'pci'},$ifstuff{'cdsi'}) = (0,0,0,0);
	my $prevtr;
	my %prevend = ();
	open(my $ncex_fh, ">$name.nc.$addn.exons.bed") or confess "\n   ERROR (sub get_tr_infos): could not open to write $name.nc.exons.$addn.bed!\n";
	open(my $ncint_fh, ">$name.nc.$addn.introns.bed") or confess "\n   ERROR (sub get_tr_infos): could not open to write $name.nc.$addn.introns.bed!\n";
	push(@listtojoin,"$name.nc.$addn.exons");
	open(my $codex_fh, ">$name.pc.$addn.exons.bed") or confess "\n   ERROR (sub get_tr_infos): could not open to write $name.pc.$addn.exons.bed!\n";
	open(my $cds_fh, ">$name.pc.$addn.CDS.bed") or confess "\n   ERROR (sub get_tr_infos): could not open to write $name.pc.$addn.CDS.bed!\n";	
	open(my $codint_fh, ">$name.pc.$addn.introns.bed") or confess "\n   ERROR (sub get_tr_infos): could not open to write $name.pc.$addn.introns.bed!\n";
	open(my $cdsint_fh, ">$name.pc.$addn.CDS.introns.bed") or confess "\n   ERROR (sub get_tr_infos): could not open to write $name.pc.$addn.CDS.introns.bed!\n";
	open(my $in_fh, "<", $in) or confess "\n   ERROR (sub get_tr_infos): could not open to read $in!\n";
	LINE: while(<$in_fh>) {
		chomp (my $line = $_);
		next LINE if (($line !~ /\w/) || ($line =~ /^MT/) || (substr($line,0,1) eq "#")); #non blank + get rid of mitochondrial DNA or commented lines
		my ($filtering,$chr,$type,$feat,$currstart,$currend,$strand,$gene_id,$tr_id,$gene_name,$tr_name,$gene_type,$tr_type,$colstoadd) = get_gtf_values($line,$set,$filter,$addcol);
		next LINE if (($feat eq "gene") || ($feat eq "transcript") || ($feat eq "mRNA"));
		$tr_type = $gene_type unless ($tr_type);
				
		#sometimes same transcript ID for different genes... So need to merge them
		my $tr_id_full = $gene_id."#".$tr_id; 
		my $currID = $tr_id_full."#".$chr.":".$currstart."-".$currend;
		
		#Filter out if relevant: -filter
		my ($fcol,$fname) = split (",",$filter);
		#1) deal with the special filter = intragenic: if non coding and not lncRNA and not < 200 nt => keep
		next LINE if ($TrMatLen->{$tr_id_full} && $filtering ne "protein_coding" && $fname eq "intragenic" && $filtering ne "lincRNA" && $TrMatLen->{$tr_id_full} <= 200);
		#2) simply filter out if it does not match
		next LINE if (($filter ne "all,all") && ($fname ne $filtering));
				
		#Filter out if relevant: -parse	
		if ($parse ne "all,all" && $addcol ne "na") {
			my ($pcol,$pname) = split (",",$parse);
			my $skip = 1;
			my @line = split('\t',$line);
			if ($line[$pcol] ne $pname) {
				next LINE;
			}
		}
				
		#OK, carry on now
		#store the whole "line", with the added columns
		my @forbig = ($chr,$type,$feat,$currstart,$currend,$strand,$gene_id,$tr_id,$gene_name,$tr_name,$gene_type,$tr_type);
		if ($addcol ne "na") {
			foreach my $col (@{$colstoadd}) {
				push(@forbig,$col);
			}
		}	
				
		#BELOW 2 commented unless statements - to filter some stuff out
		#1) not analyze stuff that are not on assembled chromosomes for mouse or human
		#unless ((($gtf =~ /mm10/) || ($gtf =~ /hg19/)) && ($chr !~ /chr[0-9XY][0-9]*/)) { 
		#2) not analyze stuff that are removed from Myotis_7x.fasta
		#unless (($gtf =~ /Myoluc2/) && ($chr =~ /chrAA/)) {
							
		#initialize to avoid problem later to deal with CDS and UTRs
		$pc{$tr_id_full}{'ATG'} = "nd" unless $pc{$tr_id_full}{'ATG'};
		$pc{$tr_id_full}{'STOP'} = "nd" unless $pc{$tr_id_full}{'STOP'};

		#all exons (will contain UTRs as well) + do introns
		my $intronend = $currstart-1 if ($prevtr);
		my $start = $currstart;
		my $end = $currend;
		my $exlen = $end - $start + 1;
		my $matlen = $exlen;
		my $nbofex = 1;
		if ($feat eq "exon"){ #Sometimes UTRs are defined, but not always if I remember well... So I'll stick with generating them with bedtools using full exons.
							  #TO DO = extract them all and compare results, to see if that's true
			#Values to get intron, only if not first line basically
			my $intronstart = $prevend{'exon'}+1 if ($prevtr);			
			if ($tr_type =~ /protein_coding/) {	
				$ifstuff{'pc'} = 1;
				print $codex_fh "$chr\t$currstart\t$currend\t$currID\t.\t$strand\n";
				#UTR introns + first coding exon introns basically
				if (($prevtr) && ($prevtr eq $tr_id_full) && ($intronend - $intronstart > 1)) {
					print $codint_fh "$chr\t$intronstart\t$intronend\t$currID\t.\t$strand\n";
					$ifstuff{'pci'} = 1;
				}	
			} else { #non coding stuff
				print $ncex_fh "$chr\t$currstart\t$currend\t$currID\t.\t$strand\n";
				if (($prevtr) && ($prevtr eq $tr_id_full) && ($intronend - $intronstart > 1)) {
					print $ncint_fh "$chr\t$intronstart\t$intronend\t$currID\t.\t$strand\n";
					$ifstuff{'nci'} = 1;
				}
			}
			$prevend{'exon'} = $currend; #Reinitialize "prev" values for introns
		
			#store all lines to loop on that later and get ExSt (unless the file exists)
			push (@big, \@forbig) unless ((-e $ExSt) || ($currend - $currstart <1));
			
			#get transcript start, end, mature len and tot nb of exons in a hash table 
			#the way it's done, it won't matter if input file is not ordered			
			if ($TrInfos{$tr_id_full}) {
				my ($trstart,$trend,$trnbofex,$trmatlen) = ($TrInfos{$tr_id_full}->[1],$TrInfos{$tr_id_full}->[2],$TrInfos{$tr_id_full}->[4],$TrInfos{$tr_id_full}->[5]);
				$start = $trstart if ($trstart <= $currstart); #=> replace start
				$end = $trend if ($trend >= $currend); #=> replace end
				$matlen = $trmatlen+$exlen; #get total mature length
				$nbofex = $trnbofex + 1; #get total number of exons
			}
			#initialize or replace values
			my @forTrI = ($chr,$start,$end,$strand,$nbofex,$matlen,$gene_type,$tr_type,$gene_name,$tr_name);
			$TrInfos{$tr_id_full} = \@forTrI;
	
		#coding exons
		} elsif ($feat eq "CDS") {
			#print bed file
			print $cds_fh "$chr\t$currstart\t$currend\t$currID\t.\t$strand\n";
			#store all lines to loop on that later and get ExSt
			push (@big, \@forbig) unless ((-e $ExSt) || ($currend - $currstart <1));				
			#get introns
			my $intronstart = $prevend{'cds'}+1 if ($prevend{'cds'});
			if (($prevend{'cds'}) && ($prevtr eq $tr_id_full) && ($intronend - $intronstart > 1)) {
				print $cdsint_fh "$chr\t$intronstart\t$intronend\t$currID\t.\t$strand\n";
				$ifstuff{'cdsi'} = 1;
			}
			$prevend{'cds'} = $currend; #Reinitialize "prev" values for introns
			#get CDS start and end
			if ($pc{$tr_id_full}{'INFO'}) {
				my ($CDSstart,$CDSend,$CDSex,$CDSlen) = ($pc{$tr_id_full}{'INFO'}->[1],$pc{$tr_id_full}{'INFO'}->[2],$pc{$tr_id_full}{'INFO'}->[4],$pc{$tr_id_full}{'INFO'}->[5]);
				$start = $CDSstart if ($CDSstart <= $currstart); #=> replace start
				$end = $CDSend if ($CDSend >= $currend); #=> replace end
				$matlen = $CDSlen+$exlen; #get total length of CDS exons
				$nbofex = $CDSex + 1; #get total number of CDS exons
			}
			#initialize or replace values
			my @forCDSI = ($chr,$start,$end,$strand,$nbofex,$matlen,$gene_type,$tr_type,$gene_name,$tr_name);
			$pc{$tr_id_full}{'INFO'} = \@forCDSI;;
	
		#start/stop
		} elsif ($feat eq "start_codon") {
			($strand eq "+")?($pc{$tr_id_full}{'ATG'} = $currstart):($pc{$tr_id_full}{'ATG'} = $currend);
		} elsif ($feat eq "stop_codon") {
			($strand eq "+")?($pc{$tr_id_full}{'STOP'} = $currstart-1):($pc{$tr_id_full}{'STOP'} = $currend+1);
		}
		#Reinitialize "prev" values for introns
		$prevtr = $tr_id_full;	
		
		#} #to be commented if non chr stuff are to be analyzed too
		#} #to be commented if AAPE scaffolds from Myotis are to be analyzed too
	}
	close($ncex_fh);
	close($codex_fh);
	close($cds_fh);
	close($in_fh);
	close($ncint_fh);
	close($codint_fh);
	close($cdsint_fh);
	if ($ifstuff{'pc'} == 1) {
		push(@listtojoin,"$name.pc.$addn.CDS");
		unless ($ifstuff{'pci'} == 1) {
			unlink "$name.pc.$addn.introns.bed";
			print STDERR "         NB: no coding introns (only single exon stuff)\n" if ($v);
		}
		unless ($ifstuff{'cdsi'} == 1) {
			unlink "$name.pc.$addn.CDS.introns.bed";
			print STDERR "         NB: no CDS introns (only single exon stuff)\n" if ($v);
		}	
	} else {`rm -Rf $name.pc.$addn.*`;}
	
	unless ($ifstuff{'nci'} == 1) {
		unlink "$name.nc.$addn.introns.bed";
		print STDERR "         NB: no non-coding introns (only single exon stuff)\n" if ($v);
	}
		
	#TrInfos :
	#  0    1      2    3       4       5       6          7        8          9        10    11      ...
	#($chr,$start,$end,$strand,$nbofex,$matlen,$gene_type,$tr_type,$gene_name,$tr_name, ADDED COLUMNS IF ANY
				
	return (\%pc,\@big,\%TrInfos,\@listtojoin);
}

#----------------------------------------------------------------------------
# subroutine to get Exson Structure file, if -f gtf
# ($feat,$ExInfo,$TrInfos) = print_ExSt($ExSt,$TrInfos,\@big,$pc,$v);
# called by main
#----------------------------------------------------------------------------
sub print_ExSt {
	my ($ExSt,$TrInfos,$big,$pc,$addcol,$v) = @_;
	my %set = %{$set};
	my %TrInfos = %{$TrInfos};
	my @big = @{$big};
	my %ExCount = ();
	my %CDSCount = ();
	my ($featC,$ExInfos);
	my @l = (); #exon structure line
	print STDERR "\n --- Printing Exon Structure file (ExSt)...\n" if ($v);
	open(my $exst_fh, ">", $ExSt) or confess "\n   ERROR (sub get_ExSt): could not open to write $ExSt!\n\n";
	print $exst_fh "#If feature = CDS, then values \"tr\" are about CDS part only (start, end, length, mature length, number of exons). \n#chr\ttype\tfeature\tExon_start\tExon_end\tstrand\texon_nb\tExon_type\ttr_id\ttr_name\ttr_biotype\ttr_chr\ttr_start\ttr_end\ttr_len\ttr_mat_len\ttr_tot_nb_exons\tgene_id\tgene_name\tgene_biotype\n\n";
	for (my $i = 0; $i <= $#big; $i++){
		my ($chr,$type,$feat,$start,$end,$strand,$gene_id,$tr_id,$gene_name,$tr_name,$gene_type,$tr_type) = ($big[$i]->[0],$big[$i]->[1],$big[$i]->[2],$big[$i]->[3],$big[$i]->[4],$big[$i]->[5],$big[$i]->[6],$big[$i]->[7],$big[$i]->[8],$big[$i]->[9],$big[$i]->[10],$big[$i]->[11]);

		#sometimes same transcript ID for different genes... So need to merge them
		my $tr_id_full = $gene_id."#".$tr_id; 
						
		#pull Tr infos -> ExSt
		if ($feat eq "exon"){	
			my ($trchr,$trstart,$trend,$trtotex,$trmatlen) = ($TrInfos{$tr_id_full}->[0],$TrInfos{$tr_id_full}->[1],$TrInfos{$tr_id_full}->[2],$TrInfos{$tr_id_full}->[4],$TrInfos{$tr_id_full}->[5]);		
			my $trlen = $trend - $trstart + 1;
		
			#get exon type; table is sorted so it's easy	
			($ExCount{$tr_id_full})?($ExCount{$tr_id_full}++):($ExCount{$tr_id_full}=1); #will increment		
			my $ExType;
			if ($trtotex==1) {
				$ExType = "SINGLE";
			} elsif ((($ExCount{$tr_id_full} == 1) && ($strand eq "+")) || (($ExCount{$tr_id_full} == $trtotex) && ($strand eq "-"))) {
				$ExType = "FIRST";
			} elsif ((($ExCount{$tr_id_full} == 1) && ($strand eq "-")) || (($ExCount{$tr_id_full} == $trtotex) && ($strand eq "+"))) {
				$ExType = "LAST";
			} else {
				$ExType = "MIDDLE";
			}		
			print $exst_fh "$chr\t$type\t$feat\t$start\t$end\t$strand\t$ExCount{$tr_id_full}\t$ExType\t$tr_id\t$tr_name\t$tr_type\t$trchr\t$trstart\t$trend\t$trlen\t$trmatlen\t$trtotex\t$gene_id\t$gene_name\t$gene_type\t";
			@l = ($chr,$type,$feat,$start,$end,$strand,$ExCount{$tr_id_full},$ExType,$tr_id,$tr_name,$tr_type,$trchr,$trstart,$trend,$trlen,$trmatlen,$trtotex,$gene_id,$gene_name,$gene_type);
		} elsif ($feat eq "CDS") {
			my $CDStype = "indet";
			($CDSCount{$tr_id_full})?($CDSCount{$tr_id_full}++):($CDSCount{$tr_id_full}=1); #will increment	
			my ($CDSchr,$CDSstart,$CDSend,$CDStotex,$CDSmatlen) = ($pc->{$tr_id_full}{'INFO'}->[0],$pc->{$tr_id_full}{'INFO'}->[1],$pc->{$tr_id_full}{'INFO'}->[2],$pc->{$tr_id_full}{'INFO'}->[4],$pc->{$tr_id_full}{'INFO'}->[5]);		
			my $CDSlen = $CDSend - $CDSstart + 1;
			
			#deal with cases that don't have start or stop codons annotated but still have UTRs...
			$pc = correct_pc_hash($pc,$strand,$tr_id_full);
			$CDStype = "FIRST"  if ((($strand eq "+") && ($start == $pc->{$tr_id_full}{'ATG'}))  || (($strand eq "-") && ($end   == $pc->{$tr_id_full}{'ATG'})));
			$CDStype = "LAST"   if ((($strand eq "+") && ($end   == $pc->{$tr_id_full}{'STOP'})) || (($strand eq "-") && ($start == $pc->{$tr_id_full}{'STOP'})));	
			$CDStype = "MIDDLE" if ((($strand eq "+") && (($start > $pc->{$tr_id_full}{'ATG'}) && ($end < $pc->{$tr_id_full}{'STOP'}))) || (($strand eq "-") && (($start > $pc->{$tr_id_full}{'STOP'}) && ($end < $pc->{$tr_id_full}{'ATG'}))));
			$CDStype = "SINGLE" if ((($strand eq "+") && ($start == $pc->{$tr_id_full}{'ATG'}) && ($end == $pc->{$tr_id_full}{'STOP'})) || (($strand eq "-") && ($start == $pc->{$tr_id_full}{'STOP'}) && ($end == $pc->{$tr_id_full}{'ATG'})));
			print $exst_fh "$chr\t$type\t$feat\t$start\t$end\t$strand\t$CDSCount{$tr_id_full}\t$CDStype\t$tr_id\t$tr_name\t$tr_type\t$CDSchr\t$CDSstart\t$CDSend\t$CDSlen\t$CDSmatlen\t$CDStotex\t$gene_id\t$gene_name\t$gene_type\t";
			@l = ($chr,$type,$feat,$start,$end,$strand,$CDSCount{$tr_id_full},$CDStype,$tr_id,$tr_name,$tr_type,$CDSchr,$CDSstart,$CDSend,$CDSlen,$CDSmatlen,$CDStotex,$gene_id,$gene_name,$gene_type);
		}
		
		#Print extra columns if relevant
		if ($addcol ne "na") {
			my @cols = split(",",$addcol);
			for (my $add = 12; $add <= $#cols+12; $add++){
				print $exst_fh "$big[$i]->[$add]\t";
				push(@l,$big[$i]->[$add]);		
			}
		}
		print $exst_fh "\n";		
		
		#Now do the "loading"
		($featC,$ExInfos,$TrInfos) = load_ExSt_line(\@l,$featC,$ExInfos,$TrInfos);	
	}	
	close $exst_fh;
	return ($featC,$ExInfos,$TrInfos);
}

#----------------------------------------------------------------------------
# Load Exon Structure line into the hash
# ($feat,$ExInfos,$TrInfos) = load_ExSt_line(\@l,$feat,$ExInfos,$TrInfos)
# called by print_ExSt and load_ExSt
#----------------------------------------------------------------------------
sub load_ExSt_line {	
	my ($l,$featC,$ExInfos,$TrInfos) = @_;
	my @l = @{$l};
#ExSt (@l)
#  0     1      2      3      4     5        6                       7        8      9          10        11      12        13      14      15         16        17        18          19
#"$chr\t$type\t$feat\t$start\t$end\t$strand\t$ExCount{$tr_id_full}\t$ExType\t$tr_id\t$tr_name\t$tr_type\t$trchr\t$trstart\t$trend\t$trlen\t$trmatlen\t$trtotex\t$gene_id\t$gene_name\t$gene_type\t";
#"$chr\t$type\t$feat\t$start\t$end\t$strand\t$CDSCount{$tr_id_full}\t$CDStype\t$tr_id\t$tr_name\t$tr_type\t$CDSchr\t$CDSstart\t$CDSend\t$CDSlen\t$CDSmatlen\t$CDStotex\t$gene_id\t$gene_name\t$gene_type\t"
	
	my ($chr,$st,$en,$strand,$ExCount,$ExType,$tr_id,$Trstart,$Trend,$tr_type,$gene_id,$feat)     
	= ($l[0],$l[3],$l[4],$l[5],$l[6],$l[7],$l[8],$l[12],$l[13],$l[10],$l[17],$l[2]);
	my $ExID = $gene_id."#".$tr_id."#".$chr.":".$st."-".$en; # gene_id + tr-id + chr:start-end of exon
	my $Exlen = $en - $st + 1;	
	$ExInfos->{$ExID}="$ExCount\t$ExType";
	
	#add additional columns in the TrInfos hash
	# 0     1      2    3       4       5       6          7        8          9
	#($chr,$start,$end,$strand,$nbofex,$matlen,$gene_type,$tr_type,$gene_name,$tr_name)
	my $tr_id_full = $gene_id."#".$tr_id;
	for (my $i = 20, my $j = 10; $i <= $#l; $i++, $j++) {
		$TrInfos->{$tr_id_full}[$j] = $l[$i];
	}
	
	#get TSS and polyA, for tot and nr counts [separate coding and non coding], but only for exons
	if ($feat eq "exon") {
		my $TSS;
		my $polyA;	
		if ($l[5] eq "-") {
			$TSS = $chr.":".$Trend;
			$polyA = $chr.":".$Trstart;
		} else { #If + or if no strand info - TO DO = deal with * (no strand info); for now is keeps start = start, but should globally be ignored for TSS and polyA counts
			$TSS = $chr.":".$Trstart;
			$polyA = $chr.":".$Trend;
		}		
		$featC->{'TSS'}{$tr_type}{$TSS}=1; #keys = nr TSS [for each transcript type]
		$featC->{'polyA'}{$tr_type}{$polyA}=1; #keys = nr polyA [for each transcript type]		
		$featC->{'tr'}{$tr_type}{$tr_id_full}=1; #keys = nb of transcripts [for each transcript type]
		
		#get splicing sites when relevant
		my ($SPL5,$SPL3) = ("na","na");
		if ($ExType eq "MIDDLE"){
			($strand eq "-")?(($SPL5,$SPL3) = ($chr.":".$st,$chr.":".$en)):(($SPL5,$SPL3) = ($chr.":".$en,$chr.":".$st));
		} elsif ($ExType eq "FIRST") {
			($strand eq "-")?($SPL3 = $chr.":".$en):($SPL5 = $chr.":".$en);
		} elsif ($ExType eq "LAST") {
			($strand eq "-")?($SPL5 = $chr.":".$st):($SPL3 = $chr.":".$st);
		}
		$featC->{'3SPL'}{$tr_type}{$SPL3}=1 if ($SPL3 ne "na");
		$featC->{'5SPL'}{$tr_type}{$SPL5}=1 if ($SPL5 ne "na");
	}
	return ($featC,$ExInfos,$TrInfos);
}

#----------------------------------------------------------------------------
# Load Exon Structure
# ($feat,$ExInfo,$TrInfos) = load_ExSt($ExSt,$TrInfos,$v);
# called by main
#----------------------------------------------------------------------------
sub load_ExSt {	
	my ($ExSt,$TrInfos,$v) = @_;
	print STDERR "     -> Loading some required info from $ExSt file\n" if ($v);
	my @ExSt = ();
	my $feat = ();
	my $ExInfos = ();
	open(my $in_fh, "<", $ExSt) or confess "\n   ERROR (sub load_ExSt): could not open to read $ExSt $!\n";
	LINE: while(<$in_fh>) {
		chomp (my $l = $_);
		next LINE if (($l !~ /\w/) || (substr($l,0,1) eq "#")); #non blank + non comment
		my @l = split('\s+', $l);
		#Load lines
		($feat,$ExInfos,$TrInfos) = load_ExSt_line(\@l,$feat,$ExInfos,$TrInfos)
	}
	close $in_fh;
	return ($feat,$ExInfos,$TrInfos);
}

#----------------------------------------------------------------------------
# subroutine to deal with cases of pc genes that don't have start or stop codons annotated but still have UTRs
# $pc = correct_pc_hash($pc,$strand,$tr_id_full);
# called by print_ExSt and get_UTRs
#----------------------------------------------------------------------------
sub correct_pc_hash {
	my ($pc,$strand,$tr_id_full) = @_;
	my ($CDSstart,$CDSend) = ($pc->{$tr_id_full}{'INFO'}->[1],$pc->{$tr_id_full}{'INFO'}->[2]);	
	if ($strand eq "+") {
		$pc->{$tr_id_full}{'ATG'} = $CDSstart if ($pc->{$tr_id_full}{'ATG'} eq "nd");
		$pc->{$tr_id_full}{'STOP'} = $CDSend if ($pc->{$tr_id_full}{'STOP'} eq "nd");
	} elsif ($strand eq "-") {
		$pc->{$tr_id_full}{'ATG'}  = $CDSend if ($pc->{$tr_id_full}{'ATG'} eq "nd");
		$pc->{$tr_id_full}{'STOP'} = $CDSstart if ($pc->{$tr_id_full}{'STOP'} eq "nd");
	}
	return $pc;
}	

#----------------------------------------------------------------------------
# subroutine to get UTRs, if -f gtf
# get_UTRs($bedtools,$name,$pc,$TrInfos,$addn,$v);
# called by main
#----------------------------------------------------------------------------
sub get_UTRs {
	my ($bedtools,$name,$pc,$TrInfos,$addn,$v) = @_;
	print STDERR "\n --- Getting UTRs...\n" if ($v);
	my $utrtemp = "$name.pc.$addn.UTRs.temp.bed";
	#			subtract b from a
	($bedtools eq "na")?(system "subtractBed -a $name.pc.$addn.exons.bed -b $name.pc.$addn.CDS.bed > $utrtemp"):(system "$bedtools/subtractBed -a $name.pc.$addn.exons.bed -b $name.pc.$addn.CDS.bed > $utrtemp");
	# => IDs kept are exon IDs => can be used to retrieve infos about this exon => that's why no need to add UTRs in ExonSt file.

	open (my $utrtemp_fh, "<",$utrtemp) or confess "\n   ERROR (sub get_UTRs): could not open to read $utrtemp!\n";
	open (my $utr_fh, ">$name.pc.$addn.UTRs.bed") or confess "\n   ERROR (sub get_UTRs): could not open to write $name.pc.$addn.UTRs.bed!\n";
	open (my $utr5_fh, ">$name.pc.$addn.5UTRs.bed") or confess "\n   ERROR (sub get_UTRs): could not open to write $name.pc.$addn.5UTRs.bed!\n";
	open (my $utr3_fh, ">$name.pc.$addn.3UTRs.bed") or confess "\n   ERROR (sub get_UTRs): could not open to write $name.pc.$addn.3UTRs.bed!\n";
	while (<$utrtemp_fh>) {
		chomp(my $line = $_);
		my ($chr,$start,$end,$ID,$x,$strand) = split('\s+',$line);
		my ($gene_id,$tr_id) = split ("#",$ID);
		my $tr_id_full = $gene_id."#".$tr_id; 
		my @infos = split ('\t',$TrInfos->{$tr_id_full});
		#now, check if 5'UTR or 3'UTR
		my $utr = "indet";
		#deal with cases that don't have start or stop codons annotated but still have UTRs
		$pc = correct_pc_hash($pc,$strand,$tr_id_full);
		$utr="5UTR" if ((($end <= $pc->{$tr_id_full}{'ATG'}) && ($strand eq "+")) || (($start >= $pc->{$tr_id_full}{'ATG'}) && ($strand eq "-")));
		$utr="3UTR" if ((($start >= $pc->{$tr_id_full}{'STOP'}) && ($strand eq "+")) || (($end <= $pc->{$tr_id_full}{'STOP'}) && ($strand eq "-")));
		$ID = $ID."#".$utr;
		print $utr_fh "$chr\t$start\t$end\t$ID\t.\t$strand\n";
		print $utr5_fh "$chr\t$start\t$end\t$ID\t.\t$strand\n" if ($utr eq "5UTR");
		print $utr3_fh "$chr\t$start\t$end\t$ID\t.\t$strand\n" if ($utr eq "3UTR");
	}
	close $utrtemp_fh;
	close $utr_fh;
	close $utr5_fh;
	close $utr3_fh;
	#unlink "$name.pc.$addn.UTRs.temp.bed";
	return;
}

#----------------------------------------------------------------------------
# subroutine to get upstream and downstream regions
# get_up_dw($name,$len,$TrInfos,$lengths,$fa,$addn,$v);
# called by main
#----------------------------------------------------------------------------
sub get_up_dw {
	my ($name,$len,$TrInfos,$lengths,$fa,$addn,$v) = @_;
	my %TrInfos = %{$TrInfos};
	my %lengths = %{$lengths};
	print STDERR "\n --- Getting upstream and downstream regions...\n" if ($v);

	#Notes:
	#If +
	#up-newstart = start-len
	#up-newend = start-1
	#dw-newstart = end+1
	#dw-newend = end+len

	#If -
	#up-newstart = end+1
	#up-newend = end+len
	#dw-newstart = start-len
	#dw-newend = start-1
	
	open(my $pcup_fh, ">$name.pc.$addn.up.bed") or confess "\n   ERROR (sub get_up_dw): could not open to write $name.pc.$addn.up.bed!\n" if (-e "$name.pc.$addn.CDS.bed");
	open(my $pcdw_fh, ">$name.pc.$addn.dw.bed") or confess "\n   ERROR (sub get_up_dw): could not open to write $name.pc.$addn.dw.bed!\n" if (-e "$name.pc.$addn.CDS.bed");
	open(my $up_fh, ">$name.nc.$addn.up.bed") or confess "\n   ERROR (sub get_up_dw): could not open to write $name.nc.$addn.up.bed!\n";
	open(my $dw_fh, ">$name.nc.$addn.dw.bed") or confess "\n   ERROR (sub get_up_dw): could not open to write $name.nc.$addn.dw.bed!\n";
	my %printwarn;
	KEY: foreach my $key (sort keys %TrInfos) {
		my ($chr,$trstart,$trend,$strand,$tr_type) = ($TrInfos{$key}->[0],$TrInfos{$key}->[1],$TrInfos{$key}->[2],$TrInfos{$key}->[3],$TrInfos{$key}->[7]);
		
		#make sure I get the chr if it is a GL random thing
		my $chrcorr = $chr;
		if ($chr =~ /GL/) {
			foreach my $key (sort keys %lengths) {
				$chrcorr = $key if ($chr =~ /^chr.*_$key/);
			}
		}
		unless ($lengths{$chrcorr}) {
			warn "\t\t$chr not found in $fa.lengths.txt file - skip\n" unless ($printwarn{$chrcorr});
			$printwarn{$chrcorr} = 1;
			next KEY;
		}
				
		my ($upstart,$dwstart,$upend,$dwend);
		if ($strand eq "+") {
			($trstart - $len < 0)?($upstart = 1):($upstart = $trstart - $len);
			($trstart - 1 == 0)?($upend = 1):($upend = $trstart - 1);
			$dwstart = $trend + 1;
			($trend + $len >= $lengths{$chrcorr})?($dwend = $lengths{$chrcorr}):($dwend = $trend + $len);	
		} else {
			$upstart = $trend + 1;
			($trend + $len >= $lengths{$chrcorr})?($upend = $lengths{$chrcorr}):($upend = $trend + $len);
			($trstart - $len < 0)?($dwstart = 1):($dwstart = $trstart - $len);
			($trstart - 1 == 0)?($dwend = 1):($dwend = $trstart);
		}
		if ($tr_type =~ /protein_coding/) {
			print $pcup_fh "$chr\t$upstart\t$upend\t$key\t.\t$strand\n" unless ($upend - $upstart < 1);
			print $pcdw_fh "$chr\t$dwstart\t$dwend\t$key\t.\t$strand\n" unless ($dwend - $dwstart < 1);
		} else {
			print $up_fh "$chr\t$upstart\t$upend\t$key\t.\t$strand\n" unless ($upend - $upstart < 1);
			print $dw_fh "$chr\t$dwstart\t$dwend\t$key\t.\t$strand\n" unless ($dwend - $dwstart < 1);
		}
	}
	close $pcup_fh if (-e "$name.$addn.pc.CDS");
	close $pcdw_fh if (-e "$name.$addn.pc.CDS");
	close $up_fh;
	close $dw_fh;
	return;
}

#----------------------------------------------------------------------------
# subroutine to subtract stuff [requires BEDtools]
# $listtojoin = subtract($listtojoin,$listoffiles,$subtract,$subid,$name,$in,$selfsub,$bedtools,$TrInfos,$addn,$cut,$v) if (($subtract) || ($selfsub eq "yes"));
# called by main
#----------------------------------------------------------------------------
sub subtract {
	my ($listtojoin,$listoffiles,$tosub,$subid,$name,$in,$selfsub,$bedtools,$TrInfos,$addn,$v) = @_;
	my @listtojoin = @{$listtojoin};
	my @listoffiles = @{$listoffiles};
	my @listtojointemp = ();
	#First get what needs to be subtracted in one file + log it, depending on the cases
	my $tosubfull;
	print STDERR "\n --- Subtracting annotated stuff from up & down regions and introns:\n" if ($v);
	if (($tosub ne "na") && ($selfsub ne "yes")) { #meaning just subtract what is defined by option -subtract
		print STDERR "     => $tosub\n" if ($v);
		$tosubfull = $tosub;
	}
	if (($tosub eq "na") && ($selfsub eq "yes")) { #meaning no file def with -subtract but should subtract input file
		if (-e "$name.pc.$addn.exons.bed") {
			$tosubfull = "$name.nc-pc.$addn.cat.bed";			
			unless (-e $tosubfull) {
				print STDERR "     => $tosubfull\n     (NB: command line executed = cat $name.nc.$addn.exons.bed $name.pc.$addn.exons.bed)\n" if ($v);
				system "cat $name.nc.$addn.exons.bed $name.pc.$addn.exons.bed > $tosubfull";
			} else {
				print STDERR "     => $tosubfull\n     (NB: $name.nc.$addn.exons.bed $name.pc.$addn.exons.bed were already concatenated)\n" if ($v);
			}
		} else {
			$tosubfull = "$name.nc.$addn.exons.bed";
			print STDERR "     => $tosubfull (no pc exons)\n" if ($v);
		}	
	}
	if (($tosub ne "na") && ($selfsub eq "yes")) { #meaning both files to subtract
		if (-e "$name.pc.$addn.exons.bed") {
			$tosubfull = "$name.nc-pc.$addn.$subid.cat.bed";
			unless (-e $tosubfull) {
				print STDERR "     => $tosubfull\n     (NB: command line executed = cat $tosub $name.nc.$addn.exons.bed $name.pc.$addn.exons.bed)\n" if ($v);
				system "cat $tosub $name.nc.$addn.exons.bed $name.pc.$addn.exons.bed > $tosubfull";
			} else {
				print STDERR "     => $tosubfull\n     (NB: $tosub and $name.nc.$addn.exons.bed and $name.pc.$addn.exons.bed were already concatenated)\n" if ($v);
			}
		} else {
			$tosubfull = "$name.nc.$addn.$subid.cat.bed";
			unless (-e $tosubfull) {
				print STDERR "     => $tosubfull\n     (NB: no pc exons; command line executed = cat $tosub $name.nc.$addn.exons.bed)\n" if ($v);
				system "cat $tosub $name.nc.$addn.exons.bed > $tosubfull";
			} else {
				print STDERR "     => $tosubfull\n     (NB: $tosub and $name.nc.$addn.exons.bed were already concatenated)\n" if ($v);
			}
		}
	}
	
	#Now actually subtract from files if relevant
	foreach my $f (@listoffiles) {
		print STDERR "      - from $f.bed\n" if ($v);
		my $out = "$f.$subid.bed";
		#			subtract b from a
		unless (-e $out) {
			($bedtools eq "na")?(system "subtractBed -a $f.bed -b $tosubfull > $out"):(system "$bedtools/subtractBed -a $f.bed -b $tosubfull > $out");
			print STDERR "           -> $out\n" if ($v);
		} else {
			print STDERR "         $out already exists, skipping...\n" if ($v);
		}
		push(@listtojointemp,"$f.$subid");
	}	
		
	#correct up/dw to keep only adjacent fragments, before joining
	print STDERR "\n --- Correcting coordinates to keep only fragments adjacent to the transcripts (for up and dw regions, if extracted)\n" if ($v);
	SUBTRACTED: foreach my $f (@listtojointemp) {
		if (($f =~ /\.up\./) || ($f =~ /\.dw\./)) {
			print STDERR "      - for $f.bed\n" if ($v);
			push(@listtojoin,"$f.corr");
			if (-e  "$f.corr.bed") {
				print STDERR "         $f.corr.bed exists, skipping\n" if ($v);
				next SUBTRACTED;
			} 
			print STDERR "           -> $f.corr.bed\n" if ($v);
			open(my $fh, "<", "$f.bed") or confess "\n   ERROR (sub subtract): could not open to read $f.bed $!\n"; 
			open(my $corrfh, ">", "$f.corr.bed") or confess "\n   ERROR (sub subtract): could not open to write $f.corr.bed $!\n"; 
			
			LINE: while(<$fh>) {
				chomp (my $l = $_);
				next LINE if ($l !~ /\w/);
				my ($chr,$st,$en,$id,$x,$strand) = split('\t',$l);
				#filter out lines where region is not a fragment adjacent to the transcript
				my ($gene_id,$tr_id,$coords) = split("#",$id);
				my $tr_id_full = "$gene_id#$tr_id";
				my ($tstart,$tend) = ($TrInfos->{$tr_id_full}[1],$TrInfos->{$tr_id_full}[2]); 					
				# if +, when upstream   it is en = tstart-1; when downstream it is st = tend+1
				# if -, when downstream it is en = tstart-1; when upstream   it is st = tend+1;
				next LINE unless (($en == $tstart-1) || ($st == $tend+1)); 
				print $corrfh "$l\n";
			}
			close $fh;
			close $corrfh;			
		} else {
			push(@listtojoin,$f);
		}
	}
	return (\@listtojoin);	
}

#----------------------------------------------------------------------------
# Join files with TEs
# $listtoanalyse = TE_join($listtojoin,$RMout,$bedtools,$v);
# called by main
#----------------------------------------------------------------------------
sub TE_join {
	my ($listtojoin,$RMoutbed,$bedtools,$v) = @_;
	my @listtojoin = @{$listtojoin};
	my @listtoanalyse = ();
	print STDERR "\n --- Joining files with $RMout:\n" if ($v);
	for (my $i=0;$i<=$#listtojoin;$i++) {
		print STDERR "      - $listtojoin[$i].bed\n" if ($v);
		my $out = "$listtojoin[$i].TEjoin.bed";
		push(@listtoanalyse,$out);
		unless (-e $out) {
			($bedtools eq "na")?(system "intersectBed -a $RMoutbed -b $listtojoin[$i].bed -wo > $out"):(system "$bedtools/intersectBed -a $RMoutbed -b $listtojoin[$i].bed -wo > $out");
		} else {
			print STDERR "          $out already exists, skipping...\n" if ($v);
		}
	}	
	return (\@listtoanalyse);
}

#----------------------------------------------------------------------------
# get amounts of files before join
# my $totlen = get_amounts($listtojoin,$name,$addn,$TrInfos,$cut,$ft,$v);
# called by main
#----------------------------------------------------------------------------
sub get_amounts {
	my ($list,$name,$addn,$TrInfos,$cut,$ft,$v) = @_;	
	print STDERR "\n --- Getting genomic coverage of all subsets before joining:\n" if ($v);
	my $amount = "$name.$addn.amounts.txt";
	`rm -Rf $amount`;
	print STDERR "      => $amount\n" if ($v);
	my @cut = split(",",$cut);
	my $genlen = ();
	foreach my $file (@{$list}) {
		print STDERR "         - $file.bed\n" if ($v);
		#load the bed file into a table
		if ($ft eq "bed" || ($ft eq "gtf" && ($file =~ /\.exons/ || $file =~ /\.CDS/ || $file =~ /\.UTR/ || $file =~ /\.introns/))) {
			$genlen = get_amounts_sub($amount,$file,$genlen,"all");
		} elsif ($ft eq "gtf" && ($file =~ /\.up\./ || $file =~ /\.dw\./)) {
			print STDERR "         => " if ($v);
			foreach my $c (@cut) {
				print STDERR "$c " if ($v);
				$genlen = get_amounts_sub($amount,$file,$genlen,$c);
			}
			print STDERR "nt \n" if ($v);
		}
	}
	return($genlen);
}

#----------------------------------------------------------------------------
# get amounts of files before join - subroutine
# $genlen = get_amounts_sub($file,\@bed,$genlen,$c);
# $genlen = get_amounts_sub($file,\@bed,$genlen,"all");
# called by get_amounts
#----------------------------------------------------------------------------
sub get_amounts_sub {
	my ($amount,$file,$genlen,$c) = @_;
	#load the file in a table
	open(my $fh, "<", "$file.bed") or confess "\n   ERROR (sub get_amounts_sub): could not open to read $file.bed $!\n";	
	my @bed = ();
	LINE: while(<$fh>) {
		chomp (my $l = $_);
		next LINE if ($l !~ /\w/);
		my @l = split('\t',$l); 	
		push(@bed, \@l);
	}
	close($fh);
	#sort the bed file
	@bed = sort {($a -> [0] cmp $b -> [0]) || ($a -> [1] <=> $b -> [1]) || ($a -> [2] <=> $b -> [2])} @bed;
	
	#now loop to get amount
	my $totlen = 0;
	my $overlap = 0;
	my ($chrU,$startU,$endU);
	for (my $i = 0; $i <= $#bed; $i++){	
		my ($chr,$start,$end,$strand) = ($bed[$i]->[0],$bed[$i]->[1],$bed[$i]->[2],$bed[$i]->[5]); 
		my $len = $end - $start + 1;
		#correct start and end for $c when relevant
		if (($c ne "all") && ($len>$c)) {
			($start = $end - $c + 1) if ((($strand eq "+") && ($file =~ /\.up\./)) || (($strand eq "-") && ($file =~ /\.dw\./)));
			($end = $start + $c - 1) if ((($strand eq "+") && ($file =~ /\.dw\./)) || (($strand eq "-") && ($file =~ /\.up\./)));
		}
		if ($i==0) {
			$totlen = $end - $start + 1; #Just get length of first line, initialize totlen
			($chrU,$startU,$endU) = ($chr,$start,$end);
		} else {
			$totlen = ($endU - $startU + 1) if ($i==1); #first round = second line => need to get length of first line
			$totlen += ($end - $start + 1); #add current length	
			#get overlap:
			$overlap += ($endU - $start + 1) if (($chr eq $chrU) && ($start >= $startU) && ($start <= $endU));
			($chrU,$startU,$endU) = ($chr,$start,$end);
		}	
	}
	$genlen->{$file}{$c} = $totlen - $overlap;	
	open(my $amount_fh, ">>", $amount) or confess "\n   ERROR (sub get_amounts_sub): could not open to write $amount $!\n";
	print $amount_fh "$file\t$c\t$genlen->{$file}{$c}\t(total length = $totlen; overlaps between features = $overlap)\n";
	close $amount_fh;
	undef (@bed);
	return ($genlen);
}

#----------------------------------------------------------------------------
# get amounts of files before join - subroutine
# my $countTE = parse_join($listtoanalyse,$TrInfos,$ExInfo,$TEinfoRMP,$TE_RMP,$TEinfo,$TEage,$addcol,$cut,$totlen,$TEov,$ft,$v);
# called by main
#----------------------------------------------------------------------------
sub parse_join {
	my ($listtoanalyse,$TrInfos,$ExInfo,$TEinfoRMP,$TE_RMP,$TEinfo,$TEage,$addcol,$cut,$totlen,$TEov,$ft,$v) = @_;
	print STDERR "\n --- Parsing joined files\n" if ($v);	
	my $countTE;
	FILE: foreach my $file (@{$listtoanalyse}) {
	# I. load the join file into a table and sort it
		my $join = load_join($file,$TEov,$v);
		unless ($join->[0][0]) {
			print STDERR "         -> No TE intersected, skip analysis\n" if ($v);	
			next FILE;
		}
		#sort the table:    by Gname                     by Gstart                    by Gend                     by $st                          by $en  
		@{$join} = sort { ($a -> [4] cmp $b -> [4]) || ($a -> [5] <=> $b -> [5]) || ($a -> [6] <=> $b -> [6]) || ($a -> [18] <=> $b -> [18]) || ($a -> [19] <=> $b -> [19])} @{$join};

	# II. loop to get TrInfos and other values+ print output files	
		print STDERR "         -> Extracting parsing information + printing outputs:\n" if ($v);
		my @cut = split(",",$cut);
		my @addcol = ();
		($addcol =~ /,/)?(@addcol = split (",",$addcol)):(push(@addcol,$addcol));
		my $add = $#addcol;
		$add = 0 if ($addcol eq "na");
		$file =~ s/\.bed$//;
		
		#0) Deal with simple bed file
		if ($ft eq "bed") {
			my $TEcat = ();	
			($countTE,$TEcat) = loop_join_simple($join,$TEinfoRMP,$TEinfo,$TEage,$countTE,$TEov,$v);
			print_OUT($file,$TEcat,$TEinfo,$TEinfoRMP,$TE_RMP,$TEage,$totlen,"all",$v);		
		#1) Deal with all exon/intron stuff
		} elsif (($file =~ /\.exons/) || ($file =~ /\.CDS/) || ($file =~ /\.UTR/) || ($file =~ /\.introns\./)) {		
			my $TrI = "$file.TrInfos.tab";		
			my $TrI_type = "E";
			$TrI_type = "C" if ($file =~ /\.CDS/);
			my $TrI_header =  "#Gene_ID\tGene_name\tGene_biotype\tTr_ID\tTr_name\tTr_biotype\tTrChr\tTrStart\tTrEnd\tTrStrand\tTrlen\tTrMatureLen\tExonNb\tExonType\tExChr\tExStart\tExEnd\tExLen\tlen_thisTE-Exon\t%thisTE-Exon\t%thisTE-MatureTr\tlen_allTEs-MatureTr(nr)\t%allTEs-MatureTr(nr)\tCategory\tTEname\tTEclass\tTEfam\tTEchr\tTEstart\tTEend\tTEstrand\tRMblock\tTEage(1)\tTEage(2)\tTE_avg_%div\n\n";
			if ($file =~ /\.introns\./) {
				$TrI_type = "I";
				$TrI_header = "#Gene_ID\tGene_name\tGene_biotype\tTr_ID\tTr_name\tTr_biotype\tTrChr\tTrStart\tTrEnd\tTrStrand\tTrlen\tTrMatureLen\tIntronChr\tIntronStart\tIntronEnd\tIntronLen\tlen_thisTE-Intron\t%thisTE-Intron\tCategory\tTEname\tTEclass\tTEfam\tTEchr\tTEstart\tTEend\tTEstrand\tRMblock\tTEage(1)\tTEage(2)\tTE_avg_%div\n\n";
			}
			my $TEcat = ();	
			($countTE,$TEcat) = loop_join($TrI,$TrI_type,$TrI_header,$TrInfos,$ExInfo,$add,$join,$TEinfoRMP,$TEinfo,$TEage,$countTE,$TEov,"all",$v);
			print_OUT($file,$TEcat,$TEinfo,$TEinfoRMP,$TE_RMP,$TEage,$totlen,"all",$v);

		#2) Deal with up and down
		} elsif (($file =~ /\.up\./) || ($file =~ /\.dw\./)) {
			foreach my $c (@cut) {
				my $TrI = "$file.$c.TrInfos.tab";
				my $TrI_type = "up/dw";
				$TrI_type = "up" if ($file =~ /\.up\./);
				$TrI_type = "dw" if ($file =~ /\.dw\./);
				my $TrI_header = "#Gene_ID\tGene_name\tGene_biotype\tTr_ID\tTr_name\tTr_biotype\tTrChr\tTrStart\tTrEnd\tTrStrand\tTrlen\tTrMatureLen\tup/dw\tup/dw_cat\tregionChr\tregionStart\tregionEnd\tregion_real_Len\tlen_thisTE-region\t%thisTE-region\tCategory\tTEname\tTEclass\tTEfam\tTEchr\tTEstart\tTEend\tTEstrand\tRMblock\tTEage(1)\tTEage(2)\tTE_avg_%div\n\n";
				my $TEcat = ();
				($countTE,$TEcat) = loop_join($TrI,$TrI_type,$TrI_header,$TrInfos,$ExInfo,$add,$join,$TEinfoRMP,$TEinfo,$TEage,$countTE,$TEov,$c,$v);
				print_OUT($file,$TEcat,$TEinfo,$TEinfoRMP,$TE_RMP,$TEage,$totlen,$c,$v);
			}
		}	
	}
	return ($countTE);
}

#----------------------------------------------------------------------------
# Load the join file in a table + filter if intersection is < 10nt
# my $join = load_join($file,$TEov,$v);
# called by parse_join
# TO DO: implement a filter on TE name, class or family.
#----------------------------------------------------------------------------
sub load_join {
	my ($in,$TEov,$v) = @_;
	print STDERR "      - loading $in in a table\n" if ($v);
	my @join = ();
	open(my $fh, "<", $in) or confess "\n   ERROR (sub parse_join/load_join): could not open to read $in $!\n";
	LINE: while (<$fh>) {
		chomp(my $l = $_);
		my @l = split('\s+', $l);
		my $Ilen = $l[-1]; #last value of the line is intersection length, -wo from Bedtools
		
		#filter if < TEov	
		if (substr($TEov,-1) eq "%") {
			my ($st,$en) = ($l[18],$l[19]);
			my $Plen = $Ilen/($en-$st+1)*100;
			my $TEovP = substr($TEov,0,-1); #should remove the last %
			next LINE if ($Plen<$TEovP);
		} elsif (substr($TEov,-1) ne "%") {
			next LINE if ($Ilen<$TEov);
		}
		
		#Correct @l
		my $startBASE1 = $l[1];
		my @RMout = split(';', $l[3]);
		$RMout[5] = $startBASE1; #correction of the start in case it was base 0
		@l = @l[4..$#l];
		unshift(@l,@RMout);
		
#deprecated now that -wo in intersectBed has the intersection length; if uncommented, the "#filter if < TEov" need to be moved after this
# 		my ($Istart,$Iend,$Ilen);
# 		my ($Gst,$Gen,$st,$en) = ($l[5],$l[6],$l[18],$l[19]);
# 		#Get intersection coordinates
# 		($Gst > $st)?($Istart = $Gst):($Istart = $st);
# 		($Gen < $en)?($Iend = $Gen):($Iend = $en);
# 		$Ilen = $Iend-$Istart+1; 
		
		#push in table
		push @join, \@l;
	}
	close ($fh);
	return (\@join);
}

#----------------------------------------------------------------------------
# print transcript infos
# ($countTE,$TEcat) = loop_join($TrI,$TrI_type,$TrI_header,$TrInfos,$ExInfo,$add,$join,$TEinfoRMP,$TEinfo,$TEage,$countTE,$TEov,"all",$v);
# ($countTE,$TEcat) = loop_join($TrI,$TrI_type,$TrI_header,$TrInfos,$ExInfo,$add,$join,$TEinfoRMP,$TEinfo,$TEage,$countTE,$TEov,$c,$v);
# called by parse_join
#----------------------------------------------------------------------------			
sub loop_join {		
	my ($TrI,$TrI_type,$TrI_header,$TrInfos,$ExInfo,$add,$join,$TEinfoRMP,$TEinfo,$TEage,$countTE,$TEov,$c,$v) = @_;
	my @join = @{$join};
	my %trinfosline = ();
	my $TEcat = (); 
	my $p = ();
	my %countchr = ();
	print STDERR "            - $TrI\n" if ($v);
	LINE: for (my $i = 0; $i <= $#join; $i++){
		#current line:
		my ($Gname,$Gstart,$Gend,$Rstrand,$Rname,$RclassFam,$block) = ($join->[$i][4],$join->[$i][5],$join->[$i][6],$join->[$i][8],$join->[$i][9],$join->[$i][10],$join->[$i][14]);
		my $Glen = $Gend - $Gstart + 1;
		my ($Rclass, $Rfam) = split ("/",$RclassFam);
		my ($st,$en,$ID,$strand) = ($join->[$i][18],$join->[$i][19],$join->[$i][20],$join->[$i][22]);
		my ($gene_id,$tr_id,$co) = split("#",$ID);
		my $Tr_ID_full = $gene_id."#".$tr_id;
		my $tr_type = $TrInfos->{$Tr_ID_full}[7];
		
		#correct start and end for $c when relevant	(up and dw)	
		my $len = $en - $st + 1;
		if (($c ne "all") && ($len>$c)) {		
			($st = $en - $c + 1) if ((($strand eq "+") && ($TrI_type eq "up")) || (($strand eq "-") && ($TrI_type eq "dw")));
			($en = $st + $c - 1) if ((($strand eq "+") && ($TrI_type eq "dw")) || (($strand eq "-") && ($TrI_type eq "up")));
		}
		#=> skip if there is no TE intersection anymore, or if new intersection is < $TEov nt
		next LINE if (($Gend - $st +1 <$TEov) || ($en - $Gstart +1 <$TEov));
		
		#get TE contribution to feature + count stuff
		my ($exnb,$ExType);
		($TrI_type =~ /E|C/)?(($exnb,$ExType) = split(/\t/,$ExInfo->{$ID})):($ExType = "na");
		my ($cat,$TElenEx,$TEperEx,$TEperTr);
		($cat,$countTE,$TElenEx,$TEperEx,$TEperTr) = get_TE_contrib($TrI_type,$TrInfos,$ExType,$countTE,$Gname,$Gstart,$Gend,$st,$en,$strand,$Tr_ID_full,$tr_type);
	
		#get TE lineage info if relevant
		my ($age1,$age2,$div) = get_TE_age_info($TEinfoRMP,$TEinfo,$TEage,$Rname);
	
		#save TrInfo line stuff for when it will be printed; reset the ID as well to take in account corrections (+get real st and en for introns)
		my $line_id = $i."##".$Tr_ID_full."##".$Gname."##".$st."##".$en."##".$ID;
		($TrI_type eq "E")?($trinfosline{$line_id}->[0] = "$TElenEx\t$TEperEx\t$TEperTr"):($trinfosline{$line_id}->[0] = "$TElenEx\t$TEperEx"); 
		$trinfosline{$line_id}->[1] = "$cat\t$Rname\t$Rclass\t$Rfam\t$Gname\t$Gstart\t$Gend\t$Rstrand\t$block\t$age1\t$age2\t$div";
		
		#count number of instances of TEs (using block to least overestimate counts)
		my $Rfullname = $Rname."#".$RclassFam;
		my $RfullID = $Rfullname."-#-".$block.":".$Gname;
		($TEcat->{'c_nr'}{$RfullID})? ($TEcat->{'c_nr'}{$RfullID}++):($TEcat->{'c_nr'}{$RfullID}=1); #all fragments will be values

		#intersection coordinates
		my ($Istart,$Iend);
		($Gstart > $st)?($Istart = $Gstart):($Istart = $st);
		($Gend < $en)?($Iend = $Gend):($Iend = $en);
		if ($i == 0) {
			($TEcat->{'tr'}{$Tr_ID_full},$TEcat->{'lf'}{$Rfullname},$TEcat->{'l'}{$RclassFam},$TEcat->{'lc'}{$Rclass}) = ($TElenEx,$TElenEx,$TElenEx,$TElenEx);
			($p->{'tr'}{$Tr_ID_full}{'e'},$p->{'lf'}{$Rfullname}{'e'},$p->{'l'}{$RclassFam}{'e'},$p->{'lc'}{$Rclass}{'e'}) = ($Iend,$Iend,$Iend,$Iend);
			if ($TEage ne "na") {
				($TEcat->{'a1'}{$age1},$TEcat->{'a2'}{$age2}) = ($TElenEx,$TElenEx);
				($p->{'a1'}{$age1}{'e'},$p->{'a2'}{$age2}{'e'}) = ($Iend,$Iend);
			}							
			next LINE;
		}		
		if ($Gname eq $join->[$i-1][4]) { #it's ordered by chr so this applies to all checks
			($p,$TEcat) = get_TEamounts($p,$TEcat,$Istart,$Iend,"tr",$Tr_ID_full);
			($p,$TEcat) = get_TEamounts($p,$TEcat,$Istart,$Iend,"lf",$Rfullname);
			($p,$TEcat) = get_TEamounts($p,$TEcat,$Istart,$Iend,"l",$RclassFam);
			($p,$TEcat) = get_TEamounts($p,$TEcat,$Istart,$Iend,"lc",$Rclass);
			($p,$TEcat) = get_TEamounts($p,$TEcat,$Istart,$Iend,"a1",$age1) if ($TEage ne "na");
			($p,$TEcat) = get_TEamounts($p,$TEcat,$Istart,$Iend,"a2",$age2) if ($TEage ne "na");
		} else { # $Gname ne $join->[$i-1][4], can't overlap
			$TEcat->{'tr'}{$Tr_ID_full}+=$TElenEx;
			$TEcat->{'lf'}{$Rfullname}+=$TElenEx;
			$TEcat->{'l'}{$RclassFam}+=$TElenEx;
			$TEcat->{'lc'}{$Rclass}+=$TElenEx;
			($p->{'tr'}{$Tr_ID_full}{'e'},$p->{'lf'}{$Rfullname}{'e'},$p->{'l'}{$RclassFam}{'e'},$p->{'lc'}{$Rclass}{'e'}) = ($Iend,$Iend,$Iend,$Iend);
			if ($TEage ne "na") {
				$TEcat->{'a1'}{$age1}+=$TElenEx;
				$TEcat->{'a2'}{$age2}+=$TElenEx;
				($p->{'a1'}{$age1}{'e'},$p->{'a2'}{$age2}{'e'}) = ($Iend,$Iend);
			}						
		}
		#get nb of chrs
		($TEcat->{'chr'}{$Rfullname})?($TEcat->{'chr'}{$Rfullname}+=1):($TEcat->{'chr'}{$Rfullname}=1) unless ($countchr{$Rfullname}{$Gname});
		$countchr{$Rfullname}{$Gname}=1;
	}
	
	#Now print
	open(my $tr_fh, ">", $TrI) or confess "\n   ERROR (sub parse_join/loop_join): could not open to write $TrI $!\n";
	print $tr_fh $TrI_header; 					
	#TrInfos hash:
	#  0    1      2    3       4       5       6          7        8          9        10    11      ...
	#($chr,$start,$end,$strand,$nbofex,$matlen,$gene_type,$tr_type,$gene_name,$tr_name, ADDED COLUMNS IF ANY
	
	foreach my $key (sort keys %trinfosline) {
		my ($i,$Tr_ID_full,$chr,$currstart,$currend,$ID) = split("##",$key); #$i."##".$Tr_ID_full."##".$Gname."##".$st."##".$end."##".$ID;
		my ($gene_id,$tr_id,$coords) = split("#",$ID);
		my ($Trchr,$Trstart,$Trend,$Trstrand,$nbofex,$Trmatlen,$gene_type,$tr_type,$gene_name,$tr_name) = ($TrInfos->{$Tr_ID_full}[0],$TrInfos->{$Tr_ID_full}[1],$TrInfos->{$Tr_ID_full}[2],$TrInfos->{$Tr_ID_full}[3],$TrInfos->{$Tr_ID_full}[4],$TrInfos->{$Tr_ID_full}[5],$TrInfos->{$Tr_ID_full}[6],$TrInfos->{$Tr_ID_full}[7],$TrInfos->{$Tr_ID_full}[8],$TrInfos->{$Tr_ID_full}[9]);			
		my $Trlen = $Trend - $Trstart +1;
		my $ExLen = $currend - $currstart +1;
		print $tr_fh "$gene_id\t$gene_name\t$gene_type\t$tr_id\t$tr_name\t$tr_type\t$Trchr\t$Trstart\t$Trend\t$Trstrand\t$Trlen\t$Trmatlen\t";
		print $tr_fh "$ExInfo->{$ID}\t" if (($TrI_type eq "E") || ($TrI_type eq "C"));
		print $tr_fh "up\t$c\t" if ($TrI_type eq "up");
		print $tr_fh "dw\t$c\t" if ($TrI_type eq "dw");
		print $tr_fh "$chr\t$currstart\t$currend\t$ExLen\t$trinfosline{$key}->[0]\t";
		if (($TrI_type eq "E") || ($TrI_type eq "C")) {			
			my $perTE_matTr = $TEcat->{'tr'}{$Tr_ID_full}/$Trmatlen*100;
			print $tr_fh "$TEcat->{'tr'}{$Tr_ID_full}\t$perTE_matTr\t";
		}
		print $tr_fh "$trinfosline{$key}->[1]\t";

		#print added columns if relevant
		unless ($add == 0) {
			for (my $a=10; $a<=$add+10; $a++) {
				print $tr_fh "$TrInfos->{$Tr_ID_full}[$a]\t";
			}
		}	
		print $tr_fh "\n";
	}
	close($tr_fh);
	return($countTE,$TEcat);		
}	

#----------------------------------------------------------------------------
# Get the TE contribution ($cat) + count stuff hit by TEs (TSS, polyA, transcripts)
# ($cat,$countTE,$TElenEx,$TEperEx,$TEperTr) = get_TE_contrib($TrI_type,$TrInfos,$ExType,$countTE,$Gname,$Gstart,$Gend,$st,$en,$strand,$Tr_ID_full,$tr_type);
# called by parse_join/loop_join
# TO DO: it works, but code should be shorten here with all these if statements stuff.
#----------------------------------------------------------------------------
sub get_TE_contrib {			
	my ($TrI_type,$TrInfos,$ExType,$countTE,$Gname,$Gstart,$Gend,$st,$en,$strand,$Tr_ID_full,$tr_type) = @_;
	#TrInfos:
		#  0    1      2    3       4       5       6          7        8          9        10    11      ...
		#($chr,$start,$end,$strand,$nbofex,$matlen,$gene_type,$tr_type,$gene_name,$tr_name, ADDED COLUMNS IF ANY

	#Get feature length
	my $Exlen = $en - $st +1;
	
	#get percentage in this exon + in this transcript, made of this TE DNA
	my ($Istart,$Iend,$Ilen);
	($Gstart > $st)?($Istart = $Gstart):($Istart = $st);
	($Gend < $en)?($Iend = $Gend):($Iend = $en);
	$Ilen = $Iend-$Istart+1;
	my $TEperEx = $Ilen/$Exlen * 100;
	my $TEperTr = $Ilen/$TrInfos->{$Tr_ID_full}[5] * 100;
	
	#=> check for TSS, PolyA, SPL, BOTH-SPL - else = exonized
	my $cat = "indet";
	if ($TrI_type eq "E"){	
		#Initialize SPL sites (will be used only when the coords are actually SPL sites)
		my ($SPL5,$SPL3);
		($strand eq "+")?(($SPL5,$SPL3) = ($Gname.":".$en,$Gname.":".$st)):(($SPL5,$SPL3) = ($Gname.":".$st,$Gname.":".$en));
		#get TSS and polyA (correct start and end of transcript based on strand)
		my ($Trstart,$Trend);
		($strand eq "+")?(($Trstart,$Trend) = ($TrInfos->{$Tr_ID_full}[1],$TrInfos->{$Tr_ID_full}[2])):(($Trstart,$Trend) = ($TrInfos->{$Tr_ID_full}[2],$TrInfos->{$Tr_ID_full}[1]));
		my $TSS = $Gname.":".$Trstart;
		my $polyA = $Gname.":".$Trend;
	
		$cat = "exonized";
		if ($Gstart < $st) { #overhang TE start side
			if ($Gend > $en) { # and end side
				if ($ExType eq "SINGLE") {			
					$cat = "TSS_polyA";
					$cat = "ATG_STOP" if ($TrI_type eq "C");
				} elsif ($ExType eq "FIRST") {
					$cat = "TSS_5SPL";
					$cat = "ATG_5SPL" if ($TrI_type eq "C");
				} elsif ($ExType eq "LAST") {
					$cat = "3SPL_polyA";
					$cat = "3SPL_STOP" if ($TrI_type eq "C");
				} else {
					$cat = "3SPL_exon_5SPL";
					$cat = "indet" if ($ExType eq "indet"); #some UTRs
				}
			} else { #not end side
				#strand -
				if ($strand eq "-") { 
					if (($ExType eq "SINGLE") || ($ExType eq "LAST")) {
						$cat = "polyA";
						$cat = "STOP" if ($TrI_type eq "C");
					} else {
						$cat = "5SPL";
					}
				}
				#strand +
				else {
					if (($ExType eq "SINGLE") || ($ExType eq "FIRST")) {
						$cat = "TSS";
						$cat = "ATG" if ($TrI_type eq "C");
					} else {
						$cat = "3SPL";
					}
				}
			}
		} elsif ($Gend > $en) { # => overhang only end side
			#strand -
			if ($strand eq "-") { 
				if (($ExType eq "SINGLE") || ($ExType eq "FIRST")) {
					$cat = "TSS";
					$cat = "ATG" if ($TrI_type eq "C");;
				} else {
					$cat = "3SPL";
				}
			}
			#strand +
			else {
				if (($ExType eq "SINGLE") || ($ExType eq "LAST")) {
					$cat = "polyA";
					$cat = "STOP" if ($TrI_type eq "C");;
				} else {
					$cat = "5SPL";
				}
			}
		}
		#get total counts - not for introns or up and down
		$countTE->{'TSS'}{$tr_type}{$TSS}=1 if ($cat =~ /TSS/);
		$countTE->{'polyA'}{$tr_type}{$polyA}=1 if ($cat =~ /polyA/);
		$countTE->{'3SPL'}{$tr_type}{$SPL3}=1 if ($cat =~ /3SPL/);
		$countTE->{'5SPL'}{$tr_type}{$SPL5}=1 if ($cat =~ /5SPL/);
		$countTE->{'tr'}{$tr_type}{$Tr_ID_full}=1;
	} elsif ($TrI_type eq "I") {
		$cat = "intronized";
		if (($Gstart < $st) &&  ($Gend > $en)) { #overhang both sides
			$cat = "5SPL_3SPL";
		} elsif (($Gstart < $st) || ($Gend > $en)) { # overhang only one side
			$cat = "SPL";
		}
	} elsif (($TrI_type eq "up") || ($TrI_type eq "dw")) {
		$countTE->{'tr-up'}{$tr_type}{$Tr_ID_full}++ if ($TrI_type eq "up");
		$countTE->{'tr-dw'}{$tr_type}{$Tr_ID_full}++ if ($TrI_type eq "dw");
		$cat = "included";
		if ($Gstart < $st) { # overhang TE start side
			if ($Gend > $en) { # and end side
				$cat = "TSS_Full" if ($TrI_type eq "up");
				$cat = "PolyA_Full" if ($TrI_type eq "dw");
			} else { # overhang only start side
				$cat = "partial";
				$cat = "PolyA" if (($strand eq "+") && ($TrI_type eq "dw"));
				$cat = "TSS" if (($strand eq "-") && ($TrI_type eq "up"));
			} 
		} elsif ($Gend > $en) { # overhang only end side
			$cat = "partial";
			$cat = "TSS" if (($strand eq "+") && ($TrI_type eq "up"));
			$cat = "PolyA" if (($strand eq "-") && ($TrI_type eq "dw"));
		}
	}	
	return ($cat,$countTE,$Ilen,$TEperEx,$TEperTr);
}

#----------------------------------------------------------------------------
# Loop join simple - get values
# my ($countTE,$TEcat) = loop_join_simple($join,$TEinfoRMP,$TEinfo,$TEage,$countTE,$TEov,$v);
# called by parse_join
#----------------------------------------------------------------------------			
sub loop_join_simple {		
	my ($join,$TEinfoRMP,$TEinfo,$TEage,$countTE,$TEov,$v) = @_;
	my @join = @{$join};
	my $TEcat = (); 
	my $p = ();
	my %countchr = ();
#NB structure of countTE = $countTE->{'tr'}{$tr_type}{$Tr_ID_full}=1;
	LINE: for (my $i = 0; $i <= $#join; $i++){
		#current line:
		my ($Gname,$Gstart,$Gend,$Rstrand,$Rname,$RclassFam,$block,$Ilen) = ($join->[$i][4],$join->[$i][5],$join->[$i][6],$join->[$i][8],$join->[$i][9],$join->[$i][10],$join->[$i][14],$join->[$i][-1]); #last column is now overlap
		my $Glen = $Gend - $Gstart + 1;
		my ($Rclass, $Rfam) = split ("/",$RclassFam);
		my ($st,$en,$ID,$strand) = ($join->[$i][18],$join->[$i][19],$join->[$i][20],$join->[$i][22]);
		my $Exlen = $en - $st +1; 

		#get the intersection start and end
		my ($Istart,$Iend);
		($Gstart > $st)?($Istart = $Gstart):($Istart = $st);
		($Gend < $en)?($Iend = $Gend):($Iend = $en);
		
		#also count number of features that have an overlap
		$countTE->{'all'}{'all'}{$ID}=1;
		
		#get TE lineage info if relevant
		my ($age1,$age2,$div) = get_TE_age_info($TEinfoRMP,$TEinfo,$TEage,$Rname);
		
		#count number of instances of TEs (using block to least overestimate counts)		
		my $Rfullname = $Rname."#".$RclassFam;
		my $RfullID = $Rfullname."-#-".$block.":".$Gname;
		($TEcat->{'c_nr'}{$RfullID})? ($TEcat->{'c_nr'}{$RfullID}++):($TEcat->{'c_nr'}{$RfullID}=1); #all fragments will be values
						
		if ($i == 0) {
			($TEcat->{'lf'}{$Rfullname},$TEcat->{'l'}{$RclassFam},$TEcat->{'lc'}{$Rclass}) = ($Ilen,$Ilen,$Ilen);
			($p->{'lf'}{$Rfullname}{'e'},$p->{'l'}{$RclassFam}{'e'},$p->{'lc'}{$Rclass}{'e'}) = ($Iend,$Iend,$Iend);
			if ($TEage ne "na") {
				($TEcat->{'a1'}{$age1},$TEcat->{'a2'}{$age2}) = ($Ilen,$Ilen);
				($p->{'a1'}{$age1}{'e'},$p->{'a2'}{$age2}{'e'}) = ($Iend,$Iend);
			}							
			next LINE;
		}		
		if ($Gname eq $join->[$i-1][4]) { #it's ordered by chr so this applies to all checks
			($p,$TEcat) = get_TEamounts($p,$TEcat,$Istart,$Iend,"lf",$Rfullname);
			($p,$TEcat) = get_TEamounts($p,$TEcat,$Istart,$Iend,"l",$RclassFam);
			($p,$TEcat) = get_TEamounts($p,$TEcat,$Istart,$Iend,"lc",$Rclass);
			($p,$TEcat) = get_TEamounts($p,$TEcat,$Istart,$Iend,"a1",$age1) if ($TEage ne "na");
			($p,$TEcat) = get_TEamounts($p,$TEcat,$Istart,$Iend,"a2",$age2) if ($TEage ne "na");
		} else { # $Gname ne $join->[$i-1][4], can't overlap
			$TEcat->{'lf'}{$Rfullname}+=$Ilen;
			$TEcat->{'l'}{$RclassFam}+=$Ilen;
			$TEcat->{'lc'}{$Rclass}+=$Ilen;
			($p->{'lf'}{$Rfullname}{'e'},$p->{'l'}{$RclassFam}{'e'},$p->{'lc'}{$Rclass}{'e'}) = ($Iend,$Iend,$Iend);
			if ($TEage ne "na") {
				$TEcat->{'a1'}{$age1}+=$Ilen;
				$TEcat->{'a2'}{$age2}+=$Ilen;
				($p->{'a1'}{$age1}{'e'},$p->{'a2'}{$age2}{'e'}) = ($Iend,$Iend);
			}						
		}
		#get nb of chrs
		($TEcat->{'chr'}{$Rfullname})?($TEcat->{'chr'}{$Rfullname}+=1):($TEcat->{'chr'}{$Rfullname}=1) unless ($countchr{$Rfullname}{$Gname});
		$countchr{$Rfullname}{$Gname}=1;
	}
	return($countTE,$TEcat);		
}	

#----------------------------------------------------------------------------
# get TEage info for this line in @join
# my ($age1,$age2,$div) = get_TE_age_info($TEinfoRMP,$TEinfo,$TEage,$Rname);
# called by parse_join/loop_join and parse_join/print_CAT
#----------------------------------------------------------------------------			
sub get_TE_age_info {		
	my ($TEinfoRMP,$TEinfo,$TEage,$Rname) = @_;
	my ($age1,$age2,$div) = ("na","na","na");
	if (($TEage ne "na") && (! $TEinfo->{"na"}) && ($TEinfo->{lc($Rname)}[4])) {
		($div,$age1,$age2) = ($TEinfo->{lc($Rname)}[3],$TEinfo->{lc($Rname)}[4],$TEinfo->{lc($Rname)}[5]);
	} elsif ($TEinfoRMP->{lc($Rname)}[5]) {
		#       0     1   2         3   4      5 
		#@TEs = class fam class/fam frg nr_frg avg%div len_masked %genome_masked len_overlap %genome_overlap(for this repeat) 			
		$div = $TEinfoRMP->{lc($Rname)}[5];
	}
	return ($age1,$age2,$div);
}	

#----------------------------------------------------------------------------
# get genomic TE amount in various subsets (transcripts, classfam, age)
# ($p,$TEcat) = get_TEamounts($p,$TEcat,$Istart,$Iend,"lf",$Rfullname);
# ($p,$TEcat) = get_TEamounts($p,$TEcat,$Istart,$Iend,"l",$RclassFam);
# ($p,$TEcat) = get_TEamounts($p,$TEcat,$Istart,$Iend,"lc",$Rclass);
# ($p,$TEcat) = get_TEamounts($p,$TEcat,$Istart,$Iend,"a1",$age1) if ($TEage ne "na");
# ($p,$TEcat) = get_TEamounts($p,$TEcat,$Istart,$Iend,"a2",$age2) if ($TEage ne "na");
# called by parse_join/loop_join
#----------------------------------------------------------------------------	
sub get_TEamounts {			
	my ($p,$hash,$Istart,$Iend,$d,$id) = @_;		
	my $pIend = $p->{$d}{$id}{'e'};	
	my $Ilen = $Iend - $Istart +1;
	if ($hash->{$d}{$id}) {		
		if ($Istart < $pIend) { #if intersection coordinates overlap with previous ones; could be the same one	
			if ($Iend > $pIend) {
				$Ilen = $Iend - $pIend +1; #TO DO: check +1??
				$hash->{$d}{$id}+=$Ilen;
			} 
		} else {
			$hash->{$d}{$id}+=$Ilen;
		}
	} else { #hash value is new, set it
		$hash->{$d}{$id}=$Ilen;
	}
	#(re)set values
	$p->{$d}{$id}{'e'} = $Iend;
	return ($p,$hash);
}

#----------------------------------------------------------------------------
# print category output file
# print_OUT($file,$TEcat,$TEinfo,$TEinfoRMP,$TE_RMP,$TEage,$totlen,$c,$v);
# called by parse_join
#----------------------------------------------------------------------------
sub print_OUT {
	my ($file,$TEcat,$TEinfo,$TEinfoRMP,$TE_RMP,$TEage,$totlen,$c,$v) = @_;
	
	#I. TEs-ratios.out
	#I.   Get values
# 	$TEcat->{'c_nr'}{$Rfullname."--".$block}
# 	$TEcat->{'c_all'}{$Rfullname."--".$Gstart."-".$Gend}
	my %TEcat_counts = ();
	my $totnrTEs = 0;
	foreach my $id (sort keys %{$TEcat->{'c_nr'}}) {
		my ($Rfullname,$block) = split("-#-",$id);
		my $nb = $TEcat->{'c_nr'}{$id};
		($TEcat_counts{'c_all'}{$Rfullname})?($TEcat_counts{'c_all'}{$Rfullname}+=$nb):($TEcat_counts{'c_all'}{$Rfullname}=$nb);
		($TEcat_counts{'c_nr'}{$Rfullname})?($TEcat_counts{'c_nr'}{$Rfullname}++):($TEcat_counts{'c_nr'}{$Rfullname}=1);
		$totnrTEs+=1;
		my ($Rname,$Rclassfam) = split("#",$Rfullname);
		my ($Rclass,$Rfam) = split("/",$Rclassfam);
		$TEcat->{'ccf'}{$Rclassfam}+=1;
		$TEcat->{'cc'}{$Rclass}+=1;
	}
	my ($totGnrTEs,$totGlenTEs,$totlenTEs) = (0,0,0);
	foreach my $id (sort keys %{$TEcat->{'lf'}}) {
		$totlenTEs+=$TEcat->{'lf'}{$id};
		my ($Rname,$Rclass,$Rfam) = split(/#|\//,$id);	
		if ($RMPARSED) {
			if (defined $TE_RMP->{'l'}{lc($Rname)}) {
                $totGlenTEs+=$TE_RMP->{'l'}{lc($Rname)};
                $totGnrTEs+=$TE_RMP->{'cnr'}{lc($Rname)};
            } else {
				print STDERR "            WARN: $Rname#$Rclass/$Rfam is in the RMout but not in the parsedRM file\n";
				print STDERR "                  (/!\\ this will generate \"0\" or \"na\" values in the TE ratio file)\n";
            }
		}              
	}		
	$file =~ s/\.TEjoin$//;
	
	#I. print for TEs-ratios.out
	my $all = "$file.$c.TEs-ratios.tab";
	print STDERR "            - $all\n" if ($v);
	open(my $allfh, ">", $all) or confess "\n   ERROR (sub parse_join/print_OUT): could not open to write $all $!\n";
	print $allfh "#\t#\t#\tIN_DATA_DSET\t#\t#\t#\t#\t#\tIN_GENOME\t#\t#\t#\t#\tAGE_INFORMATION\t#\t#\tRATIOS\t#\n"; 
	my $onehead = "Len_masked\t%_masked\tCounts\tnrCounts\t%nrCounts";
	print $allfh "#Rname\tRclass\tRfam\t$onehead\tNbChrs\t$onehead\tLineage\tAncient/LineageSpe\tavg_%div\tLen\tnrCounts\n\n"; 					
	foreach my $id (sort keys %{$TEcat->{'lf'}}) {
		my ($Rname,$Rclass,$Rfam) = split(/#|\//,$id);
		my ($a1,$a2,$div) = get_TE_age_info($TEinfoRMP,$TEinfo,$TEage,$Rname);
		#get in set info
		my $in_set_per_len = $TEcat->{'lf'}{$id} / $totlenTEs *100; 
		my $in_set_per_nrc = $TEcat_counts{'c_nr'}{$id} / $totnrTEs *100; 
		my $in_set_chrs = "nd";
		$in_set_chrs = $TEcat->{'chr'}{$id} if ($TEcat->{'chr'}{$id});
		my $in_set = "$TEcat->{'lf'}{$id}\t$in_set_per_len\t$TEcat_counts{'c_all'}{$id}\t$TEcat_counts{'c_nr'}{$id}\t$in_set_per_nrc\t$in_set_chrs";
		#get in genome and ratios inf relevant
		my $in_genome = "na\tna\tna\tna\tna";
		my ($r_len,$r_nrc) = ("na","na");
		if ($RMPARSED && $TE_RMP->{'l'}{lc($Rname)} && $totGlenTEs != 0) {
			my $in_G_per_len = $TE_RMP->{'l'}{lc($Rname)} / $totGlenTEs *100;
			my $in_G_per_nrc = $TE_RMP->{'cnr'}{lc($Rname)} / $totGnrTEs *100;
			$in_genome = "$TE_RMP->{'l'}{lc($Rname)}\t$in_G_per_len\t$TE_RMP->{'ctot'}{lc($Rname)}\t$TE_RMP->{'cnr'}{lc($Rname)}\t$in_G_per_nrc";
			$r_len = $in_set_per_len / $in_G_per_len;
			$r_nrc = $in_set_per_nrc / $in_G_per_nrc;
		}			
		print $allfh "$Rname\t$Rclass\t$Rfam\t$in_set\t$in_genome\t$a1\t$a2\t$div\t$r_len\t$r_nrc\n";
	}
	close $allfh;
	
	
	#II. CAT.out
	#    Will use: $totlen->{$file}{$c}, $TEcat
# 	$TEcat->{'tr'}{$Tr_ID_full},$TEcat->{'lf'}{$Rfullname},$TEcat->{'l'}{$RclassFam},$TEcat->{'lc'}{$Rclass}),$TEcat->{'a1'}{$age1},$TEcat->{'a2'}{$age2}
	my %Cat = ();
	my %class = ();
	#Initialize to 0 all values
	my @cols = ("Length","nrCounts","%len","%nrCounts");
	foreach my $k (@cols) {
		for(my $i=0;$i<=15;$i++) {
			$Cat{$k}->[$i] = 0;
		}
		for(my $i=0;$i<=5;$i++) {
			$class{$k}->[$i] = 0;
		}
	}		
	# CLASSES:
	#    0   1    2    3   4      5      
	#    LTR LINE SINE DNA OTHERS NONTE
	# ALL CAT:
	#    0          1        2          3              4                  5                  6                   7                   8
	#    ERV/LTRint ERV/LTRs LTR/Others nonLTR/LINE/L1 nonLTR/LINE/L2,3,4 nonLTR/LINE/Others nonLTR/SINE/B1,ALUs nonLTR/SINE/B2,MIRs nonLTR/SINE/B3,B4,others 
	#    9       10        11            12                          13         14     15
	#    DNA/hAT DNA/TcMar DNA/Helitrons DNA/Tlr,Mavericks,Polintons DNA/Others Others nonTE
	my %tot=();
	foreach my $id (sort keys %{$TEcat->{'l'}}) { #keys are class/fam		
		my $n = get_cat_nb($id,"cat");
		$Cat{'Length'}->[$n] += $TEcat->{'l'}{$id};
		$Cat{'nrCounts'}->[$n] += $TEcat->{'ccf'}{$id};
		($tot{'cat'}{'l'})?($tot{'cat'}{'l'}+=$TEcat->{'l'}{$id}):($tot{'cat'}{'l'}=$TEcat->{'l'}{$id});
		($tot{'cat'}{'c'})?($tot{'cat'}{'c'}+=$TEcat->{'ccf'}{$id}):($tot{'cat'}{'c'}=$TEcat->{'ccf'}{$id});
	}
	foreach my $id (sort keys %{$TEcat->{'l'}}) { #keys are class/fam		
		my $n = get_cat_nb($id,"cat");
		($Cat{'Length'}->[$n]==0)? ($Cat{'%len'}->[$n]=0):($Cat{'%len'}->[$n]=$Cat{'Length'}->[$n]/$tot{'cat'}{'l'}*100);
		($Cat{'nrCounts'}->[$n]==0)? ($Cat{'%nrCounts'}->[$n]=0):($Cat{'%nrCounts'}->[$n]=$Cat{'nrCounts'}->[$n]/$tot{'cat'}{'c'}*100);
	}
	foreach my $id (sort keys %{$TEcat->{'lc'}}) { #keys are class	
		my $n = get_cat_nb($id,"class");
		$class{'Length'}->[$n] += $TEcat->{'lc'}{$id};
		$class{'nrCounts'}->[$n] += $TEcat->{'cc'}{$id};
		($tot{'class'}{'l'})?($tot{'class'}{'l'}+=$TEcat->{'lc'}{$id}):($tot{'class'}{'l'}=$TEcat->{'lc'}{$id});
		($tot{'class'}{'c'})?($tot{'class'}{'c'}+=$TEcat->{'cc'}{$id}):($tot{'class'}{'c'}=$TEcat->{'cc'}{$id});
	}
	foreach my $id (sort keys %{$TEcat->{'lc'}}) { #keys are class	
		my $n = get_cat_nb($id,"class");
		($class{'Length'}->[$n]==0)? ($class{'%len'}->[$n]=0):($class{'%len'}->[$n]=$class{'Length'}->[$n]/$tot{'class'}{'l'}*100);
		($class{'nrCounts'}->[$n]==0)? ($class{'%nrCounts'}->[$n]=0):($class{'%nrCounts'}->[$n]=$class{'nrCounts'}->[$n]/$tot{'class'}{'c'}*100);
	}
	
	my $cat = $file.".$c.CAT.tab";
	print STDERR "            - $cat\n" if ($v);
	open(my $catfh, ">", $cat) or confess "\n   ERROR (sub parse_join/print_OUT): could not open to write $cat $!\n";
	print $catfh "#file:\t$file\ncut:\t$c\n#total_length(nt):\t$totlen->{$file}{$c}\n";
	print $catfh "Rclass\tERV-int/LTR\tERV-rest/LTR\tLTR/Others\tnonLTR/LINE/L1\tnonLTR/LINE/L2,3,4\tnonLTR/LINE/Others\tnonLTR/SINE/B1,ALUs\tnonLTR/SINE/B2,MIRs\tnonLTR/SINE/B3,B4,others\tDNA/hAT\tDNA/TcMar\tDNA/Helitrons\tDNA/Tlr,Mavericks,Polintons\tDNA/Others\tOthers\tnonTE\n";
	foreach my $k (@cols) {
		print $catfh "$k\t";
		for(my $i=0;$i<=15;$i++) {
			print $catfh "$Cat{$k}->[$i]\t"
		}
		print $catfh "\n";
	}	
	print $catfh "\n\n";
	close $catfh;
	
	my $class = $file.".$c.CAT-class.tab";
	print STDERR "            - $class\n" if ($v);
	open(my $classfh, ">", $class) or confess "\n   ERROR (sub parse_join/print_OUT): could not open to write $class $!\n";
	print $classfh "#file:\t$file\ncut:\t$c\n#total_length(nt):\t$totlen->{$file}{$c}\n\n#Rclass\tLTR\tLINE\tSINE\tDNA\tOthers\tnonTE\tTOTAL\t%tot\n";
	foreach my $k (@cols) {
		print $classfh "$k\t";
		my ($tot,$per) = (0,0);
		for(my $i=0;$i<=5;$i++) {
			print $classfh "$class{$k}->[$i]\t";
			$tot+=$class{$k}->[$i];
		}
		if ($k eq "Length") {
			my $per = 0;
			$per = $tot / $totlen->{$file}{$c} *100 if ($totlen->{$file}{$c} > 0);
			print $classfh "$tot\t$per";
		}
		print $classfh "\n";
		
	}	
	print $classfh "\n\n";
	close $classfh;

	#III. AGE.out if relevant
	#($TEcat->{'a1'}{$age1},$TEcat->{'a2'}{$age2})
	if ($TEage ne "na") {
		#get values for the CAT-AGE.out
		my $catage = $file.".$c.AGE.tab";
		print STDERR "            - $catage\n" if ($v);
		open(my $catagefh, ">", $catage) or confess "\n   ERROR (sub parse_join/print_OUT): could not open to write $catage $!\n";
		print $catagefh "#file:\t$file\n#total_length(nt):\t$totlen->{$file}{$c}\n\n";
		foreach my $id (sort keys %{$TEcat->{'a1'}}) { #keys are lineages
			print $catagefh "$id\t$TEcat->{'a1'}{$id}\n";
		}
		print $catagefh "\n";
		foreach my $id (sort keys %{$TEcat->{'a2'}}) { #keys are Ancient / Lineage Spe
			print $catagefh "$id\t$TEcat->{'a2'}{$id}\n";
		}
		print $catagefh "\n\n";
		close $catagefh;
	}			
	return;
}

#----------------------------------------------------------------------------
# get the number of category
# my $n = get_cat_nb($id,"cat");
# called by parse_join/print_OUT
#----------------------------------------------------------------------------
sub get_cat_nb {
	my ($id,$type) = @_;
	my ($RClass,$RFam) = split("/",$id);
	my $n;
	if ($type eq "cat") {
	# ALL CAT:
	#    0          1        2          3              4                  5                  6                   7                   8
	#    ERV/LTRint ERV/LTRs LTR/Others nonLTR/LINE/L1 nonLTR/LINE/L2,3,4 nonLTR/LINE/Others nonLTR/SINE/B1,ALUs nonLTR/SINE/B2,MIRs nonLTR/SINE/B3,B4,others 
	#    9       10        11            12                          13         14     15
	#    DNA/hAT DNA/TcMar DNA/Helitrons DNA/Tlr,Mavericks,Polintons DNA/Others Others nonTE

		$n = 14; #default coordinate in the all cat table
		# SINES -> Alu, MIR
		if ($RClass eq "SINE") {
			$n = 8;
			if (($RFam eq "ALU") || ($RFam eq "Alu")) {
				$n = 6;
			} elsif (($RFam eq "B2") || ($RFam eq "MIR")) {
				$n = 7;
			}
		}
		# LTR -> int or non int
		elsif (($RClass eq "LTR") || ($RClass eq "ERV")) {
			$n = 2;
			if ($RFam =~ /ERV/) {
				if ($RFam =~ /--int$/) {
					$n = 0;
				} else {
					$n = 1;
				} 
			}	
		} 
		# LINEs -> LINE1, LINE2 RTE CR1, others
		elsif ($RClass eq "LINE") {
			$n = 5;
			if ($RFam eq "L1") {
				$n = 3;
			} elsif (($RFam eq "L2") || ($RFam eq "CR1") || ($RFam eq "RTE")) {
				$n = 4;
			} 
		} 
		# DNA -> hAT, TcMar, others
		elsif ($RClass eq "DNA") {
			$n = 13;
			if ($RFam =~ /hAT/) {
				$n = 9;
			} elsif ($RFam =~ /TcMar/) {
				$n = 10;
			} elsif ($RFam =~ /Helitron/) {
				$n = 11;
			} elsif ($RFam =~ /[Tt]lr|[Mm]averick|[Pp]olinton/) {
				$n = 12;
			} 		
		} elsif ($RClass eq "RC"){ #Helitrons as well
			$n = 11;
		# non TE stuff
		} elsif ($RClass eq "nonTE"){
			$n = 15;
		}
	} elsif ($type eq "class") {
		# CLASSES:
		#    0   1    2    3   4      5      
		#    LTR LINE SINE DNA OTHERS NONTE
		$n = 4; #default coordinate in the all cat table
		if ($RClass eq "SINE") {
			$n = 2;
		} elsif ($RClass eq "LINE") {
			$n = 1;
		} elsif (($RClass eq "LTR") || ($RClass eq "ERV")) {
			$n = 0;
		} elsif (($RClass eq "DNA") || ($RClass eq "RC")) {
			$n = 3;
		} elsif ($RClass eq "nonTE"){
			$n = 5;
		}
	}	
	return $n;
}

#----------------------------------------------------------------------------
# Print summary file
# summary($in,$countTE,$feat,$addn,$ft,$v);
# called by main
#----------------------------------------------------------------------------
sub summary {
	my ($in,$countTE,$feat,$addn,$ft,$v) = @_;
	$in =~ s/(.*)\.gtf/$1/;
	print STDERR "\n --- Printing summary file and concatenating other outputs\n" if ($v);
	my $s = $in.".".$addn."_Summary.tab";
	print STDERR "      - $s\n" if ($v);	
	open(my $fh, ">", $s) or confess "\n   ERROR (sub summary): could not open to write $s $!\n";	
	print $fh "#SUMMARY NUMBERS FOR $in\n";

	#get numbers of TSS polyA - for each transcript type
# 	$feat->{'TSS'}{$tr_type}{$TSS}=1; #keys = nr TSS [for each transcript type]
# 	$feat->{'polyA'}{$tr_type}{$polyA}=1; #keys = nr polyA [for each transcript type]
# 	$feat->{'tr'}{$tr_type}{$tr_id_full}=1 #keys = nb of transcripts [for each transcript type]
# 	$featC->{'3SPL'}{$tr_type}{$SPL3}=1 if ($SPL3 ne "na");
# 	$featC->{'5SPL'}{$tr_type}{$SPL5}=1 if ($SPL5 ne "na");
#	$feat->{'all'}{'all'}{$in}=1 #WHEN BED FORMAT
	my %tot = ();
	foreach my $id (sort keys %{$feat}) {
		foreach my $tr_type (sort keys %{$feat->{$id}}) { 
			foreach my $id2 (sort keys %{$feat->{$id}{$tr_type}}) { 
				$tot{$id}{$tr_type}+=$feat->{$id}{$tr_type}{$id2}; #when bed format it will simply count total number of features
			}
		}
	}

# 	$countTE->{'TSS'}{$tr_type}{$TSS}=1 if ($cat =~ /TSS/);
# 	$countTE->{'polyA'}{$tr_type}{$polyA}=1 if ($cat =~ /polyA/);
# 	$countTE->{'3SPL'}{$tr_type}{$SPL3}=1 if ($cat =~ /3SPL/);
# 	$countTE->{'5SPL'}{$tr_type}{$SPL5}=1 if ($cat =~ /5SPL/);
# 	$countTE->{'tr'}{$tr_type}{$Tr_ID_full}=1;
	
	print $fh "#Counts of TSS and polyA in TEs, and transcripts with at least one TE fragment (of minimum 10nt) in the features (e.g. exons, intergenic regions)\n\n";
	print $fh "#\t(nr) means that unique TSS, splicing sites or polyA are considered (if 10 transcripts share a TSS the count will be 1).\n";
	print $fh "#\ttype\tin_TE(nr)\tTotal_In_Set(nr)\t%(nr)\n";
	foreach my $id (sort keys %{$countTE}) { #category
		foreach my $type (sort keys %{$countTE->{$id}}) { #type of transcripts
			my $nb = 0;
			foreach my $id2 (sort keys %{$countTE->{$id}{$type}}) { 
			#keys are a given TSS or a given splicing site, etc
				$nb++;
			}
			my ($tot,$per) = (0,0);
			(exists $tot{$id}{$type})?($tot = $tot{$id}{$type}):($tot = $tot{'tr'}{$type});
			$per = $nb / $tot *100 unless (! $tot || $tot == 0);
			print $fh "$id\t$type\t$nb\t$tot\t$per\n";
		}		
	}	
	close $fh;	

	#cat other outputs for easier access when $ft eq gtf
	if ($ft eq "gtf") {
		`rm -Rf $in.$addn.concat.CAT.tab` if (-e "$in.$addn.concat.CAT.tab");
		`cat $in.*.CAT.tab > $in.$addn.concat.CAT.tab`;
		`rm -Rf $in.$addn.concat.CAT-class.tab` if (-e "$in.$addn.concat.CAT.tab");
		`cat $in.*.CAT-class.tab > $in.$addn.concat.CAT-class.tab`;
		`rm -Rf $in.$addn.concat.AGE.tab` if (-e "$in.$addn.concat.AGE.tab");
		`cat $in.*.AGE.tab > $in.$addn.concat.AGE.tab` if ($TEage ne "na");
	}
	return;
}


