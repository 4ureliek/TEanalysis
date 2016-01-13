TE-analysis_Coverage
=====
version 3.0
Last update  :  Nov 2015

	perl TE-analysis_Coverage.pl -in <file.tab> -type <X>
                    [-RM <X>] [-lib <repeats.fa>] [-big] [-TEs <TEclass>] [-cons] [-force <X>] [-R] [-Rfile]
	                [-filter <type,filter>] [-contain] [-skip <X>] [-min_frg <X>] [-min_len <X>] [-cat <X>] 
	                [-v] [-h] [-chlog]
	
	This script will output the coverage of each repeat consensus in the input file
	Use -filter to to restrict this to some repeats
	This can be in the genome (RM.out) or after intersection with some features (TE-analysis_pipeline: https://github.com/4ureliek/TEanalysis)
	Output files can be directly used to plot the coverage in R (all repeats in one file with -big), use -Rfile to get command lines
	
	/!\ Because of indels, convert coordinates from genome to consensus is not super precise. 
    This is why even features with 1 nt (TSS, etc) may span more than 1 nt.
    You can use -force Rstart or -force Rend to convert both coordinates from the start or the end in consensus respectively

    MANDATORY ARGUMENT:	
     -in      =>   (STRING) input file, see below for the 3 possibilities:
     -type    =>   (STRING) 3 options:
                               -type RMout   if the input file is a repeat masker .out output file
                               -type TEjoin  if the input file is a TEjoin.bed file (from TE-analysis_pipeline, v4+)
                               -type TrInfo if the input file is a exons.TEjoin.TrInfos.tab file (from TE-analysis_pipeline, v4+)
                                             /!\ you also need to provide the original repeat masker .out file used for the analysis using -RM
     
    OPTIONAL ARGUMENTS 
     -RM      => (STRING)   If -type TrInfo, then use this to provide the repeat masker output file .out (needed for coordinates in consensus)
     -lib     => (STRING)   Provide the library file (fasta) to obtain consensus lengths (otherwise calculated from the .out info)
                            Requries Bio::SeqIO
     -big     =>   (BOOL)   add an output: a big tabulated file with all repeats in it, in addition to a file per repeat
                            Requires Array::Transpose::Ragged
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
     -R       =>   (BOOL)   To directly get the images of the plots (pdf)
                            Requires Statistics::R
     -Rfile   =>   (BOOL)   To print a file with R command lines to plot the coverage graphs
                            Behavior will be different if bigtable or not

    OPTIONAL FILTERING
     -filter  => (STRING)   To run the script on only a subset of repeats. Not case sensitive.
                            The type can be: name, class or family and it will be EXACT MATCH unless -contain is chosen as well
                              ex: name,nhAT1_ML => only fragments corresponding to the repeat named exactly nhAT1_ML will be looked at
                                   class,DNA => all repeats with class named exactly DNA (as in ...#DNA/hAT or ...#DNA/Tc1)
                                   family,hAT => all repeats with family named exactly hAT (so NOT ...#DNA/hAT-Charlie for example)
     -contain =>   (BOOL)   To check if the "name" determined with -filter is included in the value in Repeat Masker output, instead of exact match
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
                                              -> This includes the "exonized" category as well as all partial overlaps from other categories
                               -cat TSS      to plot all TSS located in the TE (i.e. for the 3 categories: TSS, TSS_5SPL, TSS_polyA)
                               -cat polyA    to plot all polyA sites located in the TE (i.e. for the 3 categories: polyA, 3SPL_polyA, TSS_polyA)
                               -cat 5SPL     to plot all 5' SPL sites located in the TE (i.e. for the 3 categories: 5SPL, 3SPL_exon_5SPL, TSS_5SPL)
                               -cat 3SPL     to plot all 3' SPL sites located in the TE (i.e. for the 3 categories: 3SPL, 3SPL_exon_5SPL, 3SPL_polyA)

    OTHER OPTIONAL ARGUMENTS
     -v       =>   (BOOL)   verbose mode, makes the script talk to you (print in STDERR)
     -v       =>   (BOOL)   print version if only option
     -h,-help =>   (BOOL)   print this usage
     -chlog   =>   (BOOL)   print changelog




TE-analysis_Shuffle
=====
version 2.0
Last update  :  Jan 13 2016

	$scriptname -p prot.gff -l linc.gff -s features_to_shuffle [-n <X>] -r genome.build [-b] -g genome.gaps [-d] [-c <X>] [-t] [-v] [-h]

		/!\\ REQUIRES Set::IntervalTree version 0.02 and won't work with v0.10
					  GAL::Annotation version later than Jan 2016 [update of is_coding]
		/!\\ Previous outputs will be moved as *.previous

	CITATION:
		- Kapusta et al. (2013) PLoS Genetics (DOI: 10.1371/journal.pgen.1003470)
		
	Description:
		Features provided in -s will be overlapped with -p and -l files, without (no_boot) or with (boot) shuffling
		Use -m to do several repetitions of no_boot (one random transcript selected for each round)
		Each bootstrap (-n) a new random transcript will be selected AND features in -s are reshuffled   
		Note that one exon may have several types of overlaps (e.g. \"SPL\" and \"exonized\"), 
		but if several TEs are exonized then the exon is counted only one time for the Exonized category.
		Output files can be processed in R (use -print to get an example of command lines with the STDERR)
		and columns are as follow:
		transcript_type\tboot_and/or_round_value\toverlap_category\tnb_overlaps\tnb_uniq_exons_in_this_category\ttotal_nb_exons_loaded\tunhit_exons_in_this_category
  
  	MANDATORY ARGUMENT:	
    -p,--prot     => (STRING) protein coding gff3 file
    -l,--linc     => (STRING) lncRNAs gff3 file
    -s,--shuffle  => (STRING) Features to shuffle = TE file, in gff3 format
                              The Repeat Masker Utils contains a script to do so, written by R. Hubley
    -r,--range    => (STRING) To know the maximum value in a given chromosome/scaffold. 
                              File should be: Name \\t length
                              From UCSC, files *.chrom.sizes
                              If you don't have such file, use -b (--build) and provide the genome fasta file for -r
    -g,--gap      => (STRING) \"gap file\" bin, chrom, chromStart, chromEnd, ix, n, size, type, bridge
                              Corresponds to the UCSC assembly gap file, or you can get it with this script:
                              https://github.com/4ureliek/DelGet/blob/master/Utilities/fasta_get_gaps.pl
                              If you don't have such file, use -d (--dogaps) and provide the genome fasta file for -g
                              (the one from UCSC may contain gaps of lenght = 1 and it will create issues)
	
 	 OPTIONAL ARGUMENTS:
    -i,--inter    => (INT)    Minimal length (in nt) of intersection in order to consider the TE included in the feature.
                              Default = 10 (to match the TEanalysis-pipeline)
    -n,--nboot    => (STRING) number of bootsraps
                              Default = 100
                              If set to 1, no bootstrap will be done
    -m,--more     => (INT)    Even in the no_boot, a random transcript is picked. Set this number to do repetitions for no_boot.
                              Default = 0 (this would mean X times more bootstraps)
    -b,--build    => (BOOL)   See above; use this and provide the genome fasta file if no range/lengths file (-r)
                              This step may take a while but will create the required file						
    -d,--dogaps   => (BOOL)   See above; use this and provide the genome fasta file if no gap file (-g)
                              This step is not optimized, it will take a while (but will create the required file)
    -c,--cat      => (STRING) Concatenate the outputs from -p and -l, boot and no boot. Must provide core filename
                              Typically, -c linc.prots.cat
    -t,--text     => (BOOL)   [NOT THERE YET] To get examples of command lines to use in R to process the outputs
    -v,--version  => (BOOL)   print the version
    -h,--help     => (BOOL)   print this usage