TEanalysis pipeline
=====
version 4.5

USAGE:
Last update  :  Mars 2015

	perl TE-analysis_pipeline_v4+.pl -i <inputfile> [-dir] -f <format> [-myf <col_details.txt>] -RMout <RepeatMasker.out> [-RMparsed RM.out.parsed] [-base <RM base>] [-TE <TE.tab>] [-TEage] [-nonTE <X>] [-fa <genome.fa>] [-subtract <what-to-subtract>] [-subid <name>] [-noselfsub] [-bedtools <path/to/bins>] [-addcol <col1,col2,etc>] [-filter <col,filter>] [-parse <col_nb,filter>] [-cut <X,X,X>] [-v] [-clean] [-chlog] [-h] [-help]
   
   SYNOPSIS
    Type -help for detailed explanations + on how to read the outputs.
    Pipeline to analyse TE composition in features (exons of transcripts, coding or non coding, transcription factor binding sites, ChIP-seq data, etc)
    
   REQUIREMENTS
    BEDtools is required for this pipeline to run (Quinlan AR and Hall IM, 2010. Bioinformatics)
    Available at: https://github.com/arq5x/bedtools2
    
   CITATION - please cite in the Methods:
      - For this pipeline using -f gtf, cite Kapusta et al. (2013) PLoS Genetics
      - For this pipeline using -f bed, cite Lynch et al. (2015) Cell Reports
      - For the use of BEDtools, Quinlan AR and Hall IM (2010) Bioinformatics
      
   DEBUGGING
     First thing to do = check that your input files are not encoded as Classic Mac (CR).
     Also, double check the usage and the doc, just in case.
     Then, if you really think your input files are OK, shoot me an email with the errors, your command line and sample files reproducing the error if possible.
   
   DETAILS OF OPTIONS (MD = mandatory and OPT = optional):
    -i         => MD  - (STRING) input file
    -dir       => OPT - (BOOL) add this if -i corresponds to a folder containing files to work with (need to contain ONLY these files!!)
    -f         => MD  - (STRING) This means more type of analysis than format of the input file.
                         chose between -f gtf or -f bed (check -help for more information)
                         gtf (default) = complex analysis (for transcripts - TSS, exons, splicing sites, etc)
                         bed = simple analysis (for TF binding sites, ChIP-seq data etc)
                               5 columns bedfile is required: chr start end unique_ID score/. strand                        
    -myf       => OPT - (STRING) if input file is not formated as a gtf/gff3 or bed: use this option with a text file to set the column numbers.
                         You still need to provide -f to determine the type of analysis.
                         To generate an example of the text file to provide, type this option alone (with the path/name of the file to create)
    -RMout     => MD  - (STRING) repeat masker output file .out
                         Even if you already have the .out.bed file, put .out in command line
                         Obviously requires to be for the same assembly file/version than the input file (ex. mm9 or mm10, hg19 or hg38, etc)
    -TEov      => OPT - (INT) minimal length (in nt) of intersection in order to consider the TE included in the feature.
                         Default = 10
                        (STRING) put a % after a number (FLOAT) to filter out when less than X% of the feature overlaps with the TE
                         Ex: -TEov 80% will skip the line if less than 80% of the feature overlaps with the TE
    -RMparsed  => OPT - (STRING) repeat masker output parsed with parseRM.pl script (with or without the -lib option)
                         Typically, the file is <RMout.out>.parseRM.all-repeats.tab (see documentation for more details on how to get this file)
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
                         (i)  use a number corresponding to the -myf file column numbers, to the column number of a bed file, or to columns 0 to 7 of a gtf file
                              For example, use -filter 0,chr1 to analyze only chromosome 1
                         (ii) use the identifier of the feature if you want to filter a gtf file or a custom file
                              For example, use -filter gene_type,lincRNA or -filter transcript_type,lincRNA to keep only lincRNAs.
                              If input file is a gtf file (with \"transcript\" lines), then you can also use a special filter 
                              = gene_type,intragenic or transcript_type,intragenic to look at everything that is not lincRNA. 
                              This will exclude all lincRNAs but keep all the non coding stuff, unless they are < 200nt [should exclude all small RNAs]
    -parse     => OPT - (STRING) to filter on a specific column before joining with TEs, to keep lines where <filter> is found in the column <col>
                         Can be used ONLY on added columns. Lines not matching it won't be printed at all in any of the files.
                         Typically, this allows to quickly check if some subsets differ (ex: -addcol 10 -parse 10,liver)
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
    
    

DOCUMENTATION
    Author       :  Aurelie Kapusta 
	Last update  :  11 Feb 2015

	- creates \"Exon Structure\" file, with a line per exon with some info added like exon type (FIRST, LAST, etc), transcript coordinates and mature transcript length
	- extract features:
		exons (split non coding and coding -> split in CDS and UTRs)
		introns
		XXkb up and down (depending what is set - ex. 10, 5 and 1)
	- for introns and up+dw, <what-to-subtract.gtf or bed> will be subtracted from the sets => outputs are <file>.subtract.<name>.bed
	- then all resulting files are intersected with TEs (RMout.out)

	This script relies on Transcript information, so gene info is not mandatory. However, it is better for analysis - note that GTF from the UCSC table browser use the same ID for gene and transcript


	When CDS, in ExSt file, exon nb will be for CDS exons only
	ExSt also contains stuff that would be filtered out. It is basically generated once per input file. Delete it if input file has been modified.

	UTRs won't be in exonSt file, just as the original \"exon\"
	UTR output - ALL UTRs => some can be undet for 5' or 3' so it might be relevant to check that output

	--------------------------------
	   FILTERING
	--------------------------------
	-filter
	-addcol and -parse
	
	This pipeline expects that info in the added columns are by transcripts. 
	If they correspond to exons, then the Exon Structure (ExSt) file and TrInfos output files will have it wrong. 
	However the -parse will work correctly if you use -clean beforehand (or manually delete the ExSt file) 

	+ filter out if intersection is < 9nt

	--------------------------------
	   REPEAT MASKER AND TE FILES
	--------------------------------
	When RM.out is converted to bed (quite long step if the file is big), nonTE sequences are filtered out and class and family are updated using 1) -TE and 2) -RMparsed if provided.
	
	This pipeline expects that all repeat names are different. 
	In older libraries there were repeats with same name but different classfam - well if it happens then the first occurence will be the one in the hash with TE info.
	You can also correct that by using -TE or by editing the file provided in -parseRM.

	To get -parseRM file, you can use this script: https://github.com/4ureliek/Parsing-RepeatMasker-Outputs/blob/master/parseRM.pl
	(with or without -lib option will work, column numbers will be corrected based on column headers)

	repeat names will be matched in lower cases. If some repeat names are different just thanks to the case that will be a problem.

	When the ERV can be identified as internal (-int or -I associated) then family of the TE is renamed with --int

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
	   File used internally summarizing all info about each exons
	
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

