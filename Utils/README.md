TE-analysis_Coverage
=====
version 3.2
Last update  :  Mar 2016

    perl TE-analysis_Coverage.pl -in <file.tab> -type <X>
                    [-RM <X>] [-lib <repeats.fa>] [-big] [-TEs <TEclass>] [-cons] [-force <X>] [-R] [-Rfile]
	                [-filter <type,filter>] [-contain] [-skip <X>] [-min_frg <X>] [-min_len <X>] [-cat <X>] 
	                [-v] [-h] [-chlog]
	
	DESCRIPTION:
    This script will output the coverage of each repeat consensus in the input file
    Use -filter to to restrict this to some repeats
    This can be in the genome (RM.out) or after intersection with some features (TE-analysis_pipeline: https://github.com/4ureliek/TEanalysis)
    Output files can be directly used to plot the coverage in R (all repeats in one file with -big), use -Rfile to get command lines
	
    /!\\ Because of indels, convert coordinates from genome to consensus is not super precise. 
    This is why even features with 1 nt (TSS, etc) may span more than 1 nt.
    You can use -force Rstart or -force Rend to convert both coordinates from the start or the end in consensus respectively

	MANDATORY ARGUMENT:	
     -in      =>   (STRING) input file, see below for the 3 possibilities:
     -type    =>   (STRING) 3 options:
                               -type RMout   if the input file is a repeat masker .out output file
                               -type TEjoin  if the input file is a TEjoin.bed file (from TE-analysis_pipeline, v4+)
                               -type TrInfo if the input file is a exons.TEjoin.TrInfos.tab file (from TE-analysis_pipeline, v4+)
                                             /!\\ you also need to provide the original repeat masker .out file used for the analysis using -RM
     
	OPTIONAL ARGUMENTS 
     -RM      => (STRING)   If -type TrInfo, then use this to provide the repeat masker output file .out (needed for coordinates in consensus)
     -lib     => (STRING)   Provide the library file (fasta) to obtain consensus lengths (otherwise calculated from the .out info)
                            Requries Bio::SeqIO
     -big     =>   (BOOL)   add an output: a big tabulated file with all repeats in it, in addition to a file per repeat
                            Requires Array::Transpose::Ragged
                            Required for -R or -Rfile
     -TEs     => (STRING)   Optional input file to correct class and/or family of (some) TEs
                            with TE information as follow, 3 first columns mandatory: 
                            Rname \\t Rclass \\t Rfam \\t Rclass/Rfam \\t etc (only first 3 columns will be used)
     -cons    =>   (BOOL)   Print an output file that has the consensus lengths for all repeats
                            This can be important because of the merging of repeats between .align and .out of Repeat Masker
                            (coordinates in consensus could be wrong; use this option to check how many instances,
                            and if it could have a great impact skip the inconsistencies with -skip in filtering options)
     -force   =>   (STRING) Not relevant for -type RMout.
                            use -force Rstart to convert both coordinates from the start in consensus 
                            use -force Rend to convert both coordinates from the end in consensus
                            Note this can create issues (values outside of the consensus length)
     -R       =>   (BOOL)   To directly get the images of the plots (pdf). -big is required
                            Requires Statistics::R
     -Rfile   =>   (BOOL)   To print a file with R command lines to plot the coverage graphs. -big is required                            
                            Behavior will be different if bigtable or not

	OPTIONAL FILTERING
     -filter  => (STRING)   To run the script on only a subset of repeats. Not case sensitive.
                            The type can be: name, class or family and it will be EXACT MATCH unless -contain is chosen as well
                              ex: name,nhAT1_ML => only fragments corresponding to the repeat named exactly nhAT1_ML will be looked at
                                   class,DNA => all repeats with class named exactly DNA (as in ...#DNA/hAT or ...#DNA/Tc1)
                                   family,hAT => all repeats with family named exactly hAT (so NOT ...#DNA/hAT-Charlie for example)
     -contain =>   (BOOL)   To check if the \"name\" determined with -filter is included in the value in Repeat Masker output, instead of exact match
                            Note that presence of [ or ] in names will create errors
                               ex: name,HERVK => all fragments containing HERVK in their name
                                   family,hAT => all repeats with family containing hAT (...#DNA/hAT, ...#DNA/hAT-Charlie, etc)
     -skip    =>    (INT)   To skip all instances of consensus length issues when > Xnt, that can arise from the merging of repeats 
                            between .align and .out of Repeat Masker (coordinates in consensus could be wrong)
                            Typically, -skip 1
     -min_frg =>    (INT)   Filter for a minimum number of fragments (e.g. coverage not so useful if less than 3 fragments for example)
                            Typically, -min_frg 3   
     -min_len =>    (INT)   Filter for a minimum length of the repeat consensus (in nt)
                            Can be useful if lots of unclasified short repeats
                            Typically, -min_len 80         
     -cat     => (STRING)   Relevant only if -TrInfo. Will filter on the category column to plot some features.
                            If there are several transcripts with the same feature it won't be counted several times.
                            5 possibilities: 
                               -cat exonized [default] to plot all exonized pieces of the TE. 
                                              -> This includes the \"exonized\" category as well as all partial overlaps from other categories
                               -cat TSS      to plot all TSS located in the TE (i.e. for the 3 categories: TSS, TSS_5SPL, TSS_polyA)
                               -cat polyA    to plot all polyA sites located in the TE (i.e. for the 3 categories: polyA, 3SPL_polyA, TSS_polyA)
                               -cat 5SPL     to plot all 5' SPL sites located in the TE (i.e. for the 3 categories: 5SPL, 3SPL_exon_5SPL, TSS_5SPL)
                               -cat 3SPL     to plot all 3' SPL sites located in the TE (i.e. for the 3 categories: 3SPL, 3SPL_exon_5SPL, 3SPL_polyA)

	OTHER OPTIONAL ARGUMENTS
     -v       =>   (BOOL)   verbose mode, makes the script talk to you (print in STDERR)
     -v       =>   (BOOL)   print version if only option
     -h,-help =>   (BOOL)   print this usage
     -chlog   =>   (BOOL)   print changelog


TE-analysis_get_big_tables
=====
version 1.0
Last update  :  Apr 2016

    perl TE-analysis_get_big_tables.pl -d directory [-chlog] [-h] [-help]
   
	SYNOPSIS
     Util for TE-analysis_pipeline.pl (https://github.com/4ureliek/TEanalysis)
     Use to merge results of several runs: will print tables that will 
     facilitate plotting / compararisons 
    
	REQUIREMENTS
     - R installed on the computer / server
     - perl module Statistics::R
    
	CITATION
     - For the use of this script, you may cite Kapusta et al. (2013) PLoS Genetics (DOI: 10.1371/journal.pgen.1003470)
       but also include the GitHub link to this script

	MANDATORY ARGUMENTS:
     -d,--dir     => (STRING) directory with a bunch of the outputs from TE-analysis_pipeline.pl
                              Will be processed:
                                *_Summary.tab
                                *.concat.CAT-class.tab
                                *.concat.CAT.tab
                                *.concat.AGE.tab (if any)

	OPTIONAL ARGUMENTS:
     -c,--chlog   => (BOOL)   print log of changes
     -v,--version => (BOOL)   print the version
     -h,--help    => (BOOL)   print this usage


TE-analysis_Shuffle
=====
version 4.0
Last update  :  Mar 31 2016

    perl TE-analysis_Shuffle.pl -l lncRNA.gff -p prot.gff [-o <X>] [-m <X>] -s features_to_shuffle [-n <X>] 
                               [-f] -e exclude.range [-d] -r genome.range [-b] [-i include.range] [-a] 
                               [-l no_low] [-t <type,name>] [-c] [-w <bedtools_path>] [-v] [-h]

    /!\\ REQUIRES - Bedtools, v2.25+
	              - GAL::Annotation version later than Jan 2016 [update of is_coding]
	                see https://github.com/The-Sequence-Ontology/GAL
	                If issues, open the script in a text editor and comment line 15,	                
	                as well as the lines with read_gff subroutine, and uncomment lines with subroutine load_gene_tr
    /!\\ Previous outputs, if any, will be moved as *.previous [which only saves results once]
  
	CITATION:
    - For the use of this script, you may cite Kapusta et al. (2013) PLoS Genetics (DOI: 10.1371/journal.pgen.1003470)
      but also include the GitHub link to this script
    - for BEDtools, Quinlan AR and Hall IM (2010) Bioinformatics (DOI: 10.1093/bioinformatics/btq033)

	DESCRIPTION:
    Features provided in -s will be overlapped with -p and/or -l files,
    without (no_boot) or with (boot) shuffling (on same chromosome)
    One of -p or -l is mandatory. Having both in the same run means that they are 
    intersected with the same TE files, which may be better for comparisons, 
    but does not seem necessary with high bootstraps.
       
    A random transcript per gene is selected: use -m to do several repetitions of no_boot
    
    For each bootstrap (-n) with shuffling features in -s, transcripts are randomly selected as well
    Note that high bootstraps takes a lot of time.
    Shuffling is done by default with allowing overlaps between shuffled features, bceause it is
    faster and OK when over representation of specific repeats are considered.  
    Note that because TEs are often fragmented + there are inversions, the counts for the exonized TEs is likely inflated;
    this also means that when TEs are shuffled, there are more fragments than TEs. Some should be moved non independently, 
    or the input file should be corrected when possible to limit that issue [not implemented in this script for now]
    
    Note that one exon may have several types of overlaps (e.g. "SPL" and "exonized"),
    but each exon is counted only one time for each category (important for "exonized").
    However for the results per repeat, each hit is counted, unless it's the same TE
   
    If you need to generate the <genome.gaps> file but you would also like to add more files to the -e option, 
    just do a first run with no bootstraps (in this example the genome.range is also being generated):
    perl ~/bin/$scriptname -f input.bed -s genome.out -r genome.fa -b -e genome.fa -d -n 0

    Two-tailed permutation test is done on the counts of overlaps for categories
    and the results are in a .stats.cat.txt file.
    If -f is used then stats are also made on each repeat, with two-tailed 
    permutation and binomial tests and the results are in a .stats.TE.txt file.     
  
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

    
TE-analysis_Shuffle_bed
=====
version 2.0
Last update  :  Mar 21 2016

    perl TE-analysis_Shuffle_bed -f features.bed [-o <nt>] -s features_to_shuffle [-n <nb>] 
                                 -r <genome.range> [-b] -e <genome.gaps> [-d] [-i <include.range>] 
                                [-a] [-l <if_nonTE>] [-t <filterTE>] [-c] [-w <bedtools_path>] [-v] [-h]

    /!\\ REQUIRES: Bedtools, versionr ecent enough to include the Shuffle
    /!\\ Previous outputs, if any, will be moved as *.previous [which means previous results are only saved once]

	DESCRIPTION:
    Features provided in -s will be overlapped with -i file (which must be simple intervals in bed format), 
       without (no_boot) or with (boot) shuffling (on same chromosome)
       One feature may overlap with several repeats and all are considered.
       Note that because TEs are often fragmented + there are inversions, the counts of TEs are likely inflated;
       this also means that when TEs are shuffled, there are more fragments than TEs. Some should be moved non independently, 
       or the input file should be corrected when possible to limit that issue [not implemented in this script for now]
       
    If you need to generate the <genome.gaps> file but you would also like to add more files to the -e option, 
       just do a first run with no bootstraps (in this example the genome.range is also being generated):
       perl ~/bin/$scriptname -f input.bed -s genome.out -r genome.fa -b -e genome.fa -d -n 0

    Two-tailed permutation test ans a binomial test are done on the counts of overlaps. 
       The results are in a .stats.txt file. Note that high bootstraps takes a lot of time. 
  
	MANDATORY ARGUMENTS:	
    -f,--feat     => (STRING) ChIPseq peaks, chromatin marks, etc, in bed format
                              /!\\ Script assumes no overlap between peaks
    -s,--shuffle  => (STRING) Features to shuffle = TE file
                              For now, can only be the repeat masker .out or the .bed file generated by the TE-analysis_pipeline script
                              No filtering on TEs in this script                           
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
    -n,--nboot    => (STRING) number of bootsraps with shuffled -s file
                              Default = 100 for faster runs; use higher -n for good pvalues 
                              (-n 10000 is best for permutation test but this will take a while)
                              If set to 0, no bootstrap will be done
    -b,--build    => (BOOL)   See above; use this and provide the genome fasta file if no range/lengths file (-r)
                              This step may take a while but will create the required file	
    -d,--dogaps   => (BOOL)   See above; use this and provide the genome fasta file if no gap file (-g)
                              If several files in -e, then the genome needs to be the first one.
                              This step is not optimized, it will take a while (but will create the required file)                       
    -i,--incl     => (STRING) To use as -incl for bedtools shuffle: \"coordinates in which features from -i should be placed.\"
                              Bed of gff format. Could be intervals close to TSS for example.
                              More than one file (same format) may be provided (comma separated), 
                              they will be concatenated (in a file = first-file-name.cat.bed)
    -a,--add      => (BOOL)   to add the -noOverlapping option to the bedtools shuffle command line, 
                              and therefore NOT allow overlaps between the shuffled features.
                              This may create issues mostly if -i is used (space to shuffle into smaller than features to shuffle)
    -l,--low      => (STRING) To set the behavior regarding non TE sequences: all, no_low, no_nonTE, none
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