#!/usr/bin/perl -w
#######################################################
# Author  :  Aurelie Kapusta
# version :  see below / see changelog
# email   :  4urelie.k@gmail.com  
# github  :  https://github.com/4ureliek?tab=repositories
#######################################################
#load modules
use warnings;
use strict;
use Carp;
use Getopt::Long;

my $version = "1.0";
my $scriptname = "TE-analysis_get_big_tables.pl";
my $changelog = "
# Change log for $scriptname:
#	v1.0 = 1 Apr 2016
\n";

my $usage = "\nUsage [$version]: 
    perl $scriptname -d directory [-chlog] [-h] [-help]
   
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
\n";

#-----------------------------------------------------------------------------
#------------------------------ LOAD AND CHECK -------------------------------
#-----------------------------------------------------------------------------

my ($dir,$chlog,$v,$help);
my $opt_success = GetOptions(
			 	  'dir=s'		=> \$dir,
			 	  'chlog'       => \$chlog,
			 	  'version'     => \$v,
			 	  'help'		=> \$help,);
			 	  
#Check options, if files exist, etc
die "\n --- $scriptname version $version\n\n" if $v;
die $changelog if $chlog;
die $usage if $help || ! $opt_success;
die $usage unless ($dir);
die "\n -d $dir does not exist?\n\n"  if (! -e $dir);


#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------
print STDERR "\n --- $scriptname v$version\n";

#Loading any summary files
print STDERR " --- Loading *._Summary.tab files\n";
my $summary = load_summary($dir);

#Loading any CAT-class files
print STDERR " --- Loading *.concat.CAT-class.tab files\n";
my @files = `ls $dir | grep "concat.CAT-class.tab"`;
my $class = load_cat($dir,\@files);

#Loading any CAT files
print STDERR " --- Loading *.concat.CAT.tab files\n";
@files = `ls $dir | grep "concat.CAT.tab"`;
my $cat = load_cat($dir,\@files); #same sub works, even if it will fo the totlen 2 times

#Loading any AGE files
print STDERR " --- Loading *.concat.AGE.tab files (if any)\n";
my ($age,$totlen) = load_age($dir);

#Print data now
print STDERR " --- Print data in tables\n";
my $out = print_plots($dir,$summary,$class,$cat,$age,$totlen);

print STDERR " --- $scriptname done\n";
print STDERR "      => $out\n\n";

exit;

#-----------------------------------------------------------------------------
#-------------------------------- SUBROUTINES --------------------------------
#-----------------------------------------------------------------------------
#----------------------------------------------------------------------------
# get a filename from a full path
# my $name = filename($filename);
#----------------------------------------------------------------------------
sub filename {
	my($file) = shift;
	$file =~ s/.*\/(.*)$/$1/;
	return $file;
}

#-----------------------------------------------------------------------------
# Load all the summary files
# my $summary = load_summary($dir);
#-----------------------------------------------------------------------------
sub load_summary{ 
    my ($dir) = @_;   
    my @files = `ls $dir | grep "_Summary.tab"`;
    my $summary = ();
    foreach my $file (@files) {
    	chomp($file);
    	open(my $fh, "<$dir/$file") or confess "\n   ERROR (sub load_summary): could not open to read $dir/$file $!\n"; 
    	LINE: while (<$fh>) {
    		chomp(my $l = $_);
    		my @l = split(/\s+/,$l);
    		my $fname = filename($file); #in the summary it's missing the parse info if any
    		$fname =~ s/_Summary.tab//;
    		next LINE if ((substr($l,0,1) eq "#") || ($l !~ /\w/));   		
			$summary->{$l[0]}{$l[1]}{$fname}{'in_TEs'}=$l[2];
			$summary->{$l[0]}{$l[1]}{$fname}{'in_set'}=$l[3];
			$summary->{$l[0]}{$l[1]}{$fname}{'nr_per'}=$l[4];
    	}
    	close($fh);
    }
    return ($summary);
    
#STRUCTURE:

# #SUMMARY NUMBERS FOR UCSC_lncRNAs_with-DE.mm10				
# #Counts of TSS and polyA in TEs, and transcripts with at least one TE fragment (of minimum 10nt) in the features (e.g. exons, intergenic regions)				
# 				
# #	(nr) means that unique TSS, splicing sites or polyA are considered (if 10 transcripts share a TSS the count will be 1).			
# # 	type	in_TE(nr)	Total_In_Set(nr)	%(nr)
# 3SPL	na		687			3207				21.42188962
# 5SPL	na		532			3226				16.49101054
# TSS	na		139			1367				10.16825165
# polyA	na		426			1521				28.00788955
# tr	na		1609		1929				83.41109383
# tr-dw	na		679			1929				35.19958528
# tr-up	na		777			1929				40.27993779
}

#-----------------------------------------------------------------------------
# Load all the .concat.CAT-class.tab AND all the .concat.CAT.tab files
# my $class = load_cat($dir,\@files);
# my $cat = load_cat($dir,\@files);
#-----------------------------------------------------------------------------
sub load_cat{ 
    my ($dir,$files) = @_;   
    my $class = ();
    my ($fname,$type);    
    my @classes;
    foreach my $file (@{$files}) {
    	chomp($file);    	
    	open(my $fh, "<$dir/$file") or confess "\n   ERROR (sub load_cat): could not open to read $dir/$file $!\n"; 
    	LINE: while (<$fh>) {
    		chomp(my $l = $_);
    		next LINE if ($l !~ /\w/);  
    		my @l = split(/\s+/,$l);    		
    		#First line of a new file:
    		($fname,$type) = get_name_type($l[1]) if ($l[0] eq "#file:"); 	
    		next LINE if ($l[0] eq "#file:");
    		#Get the full type for the up and dw, but also add 0s
    		if ($l[0] eq "cut:") {
    			my $type2 = $l[1];
#    			$type2 = sprintf("%06d",$l[1]) unless ($type2 eq "all"); #not that great
				$type = $type."_".$type2;
				next LINE;
			}
#			#Get total length of this set
#    		$totlen->{$fname}{$type} = $l[1] if ($l[0] eq "#total_length(nt):");
			next LINE if ($l[0] eq "#total_length(nt):");
			#Now get values for all classes:
			@classes = @l if ($l[0] eq "#Rclass");
			next LINE if ($l[0] eq "#Rclass");						
			#all other lines = data to put in hash
			CLASS: for (my $c = 1; $c < $#classes; $c++) {
				next CLASS if ($classes[$c] eq "TOTAL");
				$class->{$l[0]}{$fname}{$type}{$classes[$c]}=$l[$c];
			}
    	}
    	close($fh);
    }
    return ($class);
    
#STRUCTURE of the CAT_class (CAT is very similar): 

# #file:	UCSC_lncRNAs_with-DE.mm10.nc.all.all.exons							
# cut:	all							
# #total_length(nt):	2945120							
# 								
# #Rclass	LTR	LINE	SINE	DNA	Others	nonTE	TOTAL	%tot
# Length	377539	193221	287768	34838	5383	0	898749	30.51654941
# nrCounts	1245	671	2346	233	37	0		
# %len	42.00716774	21.49888345	32.01872825	3.876276914	0.598943643	0		
# %nrCounts	27.47131509	14.80582524	51.76522507	5.141218005	0.816416593	0		
}

#-----------------------------------------------------------------------------
# Get the type of input, if exons, introns, up or dw
# my ($fname,$type) = get_name_type($l[1]) if ($l[0] eq "#file:");
#-----------------------------------------------------------------------------
sub get_name_type { 
	my $name = shift;    
    my ($fname,$type) = ("nd","nd");
    ($fname,$type) = ($1,"02.exons") if ($name =~ /^(.*)\.exons$/);  
    ($fname,$type) = ($1,"04.introns") if ($name =~ /^(.*)\.introns\..+?$/);  
    ($fname,$type) = ($1,"01.up") if ($name =~ /^(.*)\.up\..+?\.corr$/);  
    ($fname,$type) = ($1,"03.dw") if ($name =~ /^(.*)\.dw\..+?\.corr$/);   
    ($fname,$type) = ($1,"02.CDS") if ($name =~ /^(.*)\.CDS$/);  
    ($fname,$type) = ($1,"02.UTR") if ($name =~ /^(.*)\.UTRs$/);      
	return ($fname,$type);
}   

#-----------------------------------------------------------------------------
# Load all the .concat.AGE.tab
# my ($age,$totlen) = load_age($dir);
#-----------------------------------------------------------------------------
sub load_age { 
    my ($dir) = @_;   
    my @files = `ls $dir | grep "concat.AGE.tab"`;
    return "na" unless ($files[0]);
    my $age = ();
    my ($fname,$type);
    my $totlen = ();
    foreach my $file (@files) {
    	chomp($file);
    	open(my $fh, "<$dir/$file") or confess "\n   ERROR (sub load_age): could not open to read $dir/$file $!\n"; 
    	LINE: while (<$fh>) {
    		chomp(my $l = $_);
    		my @l = split(/\s+/,$l);
    		next LINE if ($l !~ /\w/); 
    		($fname,$type) = get_name_type($l[1]) if ($l[0] eq "#file:");
    		$totlen->{$fname}{$type} = $l[1] if ($l[0] eq "#total_length(nt):");
    		next LINE if (substr($l,0,1) eq "#");   
    		$l[0] = "01_Mus/rat" if ($l[0] eq "01_Mus/Rat");		
			$age->{$fname}{$type}{$l[0]}=$l[1];
    	}
    	close($fh);
    }
    return ($age,$totlen);

# STRUCTURE:

# #file:	UCSC_lncRNAs_with-DE.mm10.nc.all.all.exons
# #total_length(nt):	2945120
# 	
# 01_Mus	217968
# 01_Mus/Rat	15575
# 01_Mus/rat	13069
# 02_Muridae	226860
# 03_Euarchontoglires/Primates?	2094
# 03_Primates	1387
# 03_Rodentia	261485
# 04_Euarchontoglires	2750
# 04_Glires	182
# 04_Primates/Glires	2226
# 05_Boreotheria	11047
# 05_Euarchontoglires	2504
# 06_Boreotheria/Eutheria?	275
# 07_Eutheria	90293
# 07_Eutheria+Metatheria	3650
# 08_Eutheria/Theria?	3340
# 09_Theria	1548
# 10_Theria/Mammalia	2789
# 11_Mammalia	18365
# 12_Amniota	1409
# 14_Tetrapoda	676
# 	
# Ancient	132882
# LineageSpe	774285
}   

#-----------------------------------------------------------------------------
# Now print stuff in a way that will make it easy to plot
# print_plots($dir,$summary,$class,$cat,$age,$totlen);
#-----------------------------------------------------------------------------
sub print_plots {
	my ($dir,$summary,$class,$cat,$age,$totlen) = @_;
	$dir = $1 if $dir =~ /(.*)\/$/;
	my $out = $dir.".tab";
	open(my $fh, ">",$out) or confess "\n   ERROR (sub RMtobed): could not open to write $out $!\n";
	
	#I. SUMMARY
	print $fh "#I. Summary:\n";
	print $fh "#File\tinput_type(if_any)\tCategory\tsubset\t%(nr)\tin_TE(nr)\tTotal_In_Set(nr)\n";
	foreach my $cat (keys %{$summary}) {
		foreach my $in_type (keys %{$summary->{$cat}}) {
			foreach my $file (keys %{$summary->{$cat}{$in_type}}) {
				my ($fname,$subset) = get_subset($file);
				print $fh "$fname\t$in_type\t$cat\t$subset\t$summary->{$cat}{$in_type}{$file}{'nr_per'}\t$summary->{$cat}{$in_type}{$file}{'in_TEs'}\t$summary->{$cat}{$in_type}{$file}{'in_set'}\n";
			}	
		}	
	}
	
	#II. CAT CLASSES
	print $fh "\n#II. Repeat classes:\n";
	print $fh "#File\tcategory\tsubset\ttype\t";
	CATPREP: foreach my $cat (keys %{$class}) {
		foreach my $file (keys %{$class->{$cat}}) {
			foreach my $type (sort keys %{$class->{$cat}{$file}}) {
				foreach my $Rclass (keys %{$class->{$cat}{$file}{$type}}) {
					print $fh "$Rclass\t";
				}
				last CATPREP;
			}
		}
	}
	print $fh "\n";	
	foreach my $cat (keys %{$class}) {
		foreach my $file (keys %{$class->{$cat}}) {
		
#			print STDERR "Looping in class hash: $file\n";
		
			foreach my $type (sort keys %{$class->{$cat}{$file}}) {	
				my ($fname,$subset) = get_subset($file);
				print $fh "$fname\t$cat\t$subset\t$type\t";	
				foreach my $Rclass (keys %{$class->{$cat}{$file}{$type}}) {
					print $fh "$class->{$cat}{$file}{$type}{$Rclass}\t";
				}
				print $fh "\n";
			}	
		}	
	}
	
 	#III. AGE - use the totlen to get % instead of just raw amounts
	close $fh if ($age eq "na");
	return ($out) if ($age eq "na");
	print $fh "\n#III. Age (in % of set size):\nFile\tsubset\ttype\t";	
	$age = populate_down3($age); #I need to get 0 values for the empty ones
	AGEPREP: foreach my $file (keys %{$age}) {		
		foreach my $type (sort keys %{$age->{$file}}) {
			foreach my $agecat (sort keys %{$age->{$file}{$type}}) {
				print $fh "$agecat\t";
			}
			last AGEPREP;
		}	
	}	
	print $fh "\n";
	foreach my $file (sort keys %{$age}) {	
		TYPE: foreach my $type (keys %{$age->{$file}}) {
			my ($fname,$subset) = get_subset($file);
			next TYPE if (($fname =~ /\.nc$/) && (($type eq "02.CDS") || ($type eq "02.UTR")));
			next TYPE if (($fname =~ /\.pc$/) && ($type eq "02.exons"));
			next TYPE if (($fname =~ /\.pc\.all/) && ($type ne "04.introns") && ($subset eq "all.CDS"));
			next TYPE if (($fname =~ /\.pc$/) && (($subset =~ /ncRNA/) || ($subset =~ /genic/)));
			print $fh "$fname\t$subset\t$type\t";
			foreach my $agecat (sort keys %{$age->{$file}{$type}}) {				
				my $per = "na";
				$per = $age->{$file}{$type}{$agecat}/$totlen->{$file}{$type} if ($totlen->{$file}{$type});
				print $fh "$per\t";
			}
			print $fh "\n";
		}		
	}
	print $fh "\n";	

	
# 	AGEPREP: foreach my $agecat (sort keys %{$age}) {
# 		foreach my $file (keys %{$age->{$agecat}}) {
# 			foreach my $type (sort keys %{$age->{$agecat}{$file}}) {
# 				print $fh "$file\t";
# 			}			
# 		}
# 		last AGEPREP;		
# 	}
# 	print $fh "\n#\t";
# 	AGEPREP2: foreach my $agecat (sort keys %{$age}) {
# 		foreach my $file (keys %{$age->{$agecat}}) {
# 			foreach my $type (sort keys %{$age->{$agecat}{$file}}) {
# 				print $fh "$type\t";
# 				$age->{$agecat}{$file}{$type} = 0 unless ($age->{$agecat}{$file}{$type});
# 			}
# 		}
# 		last AGEPREP2;
# 	}	
# 	print $fh "\n";
# 	AGECAT2: foreach my $agecat (sort keys %{$age}) {		
# 		next AGECAT2 if (($agecat ne "LineageSpe") && ($agecat ne "Ancient"));
# 		print $fh "$agecat\t";
# 		foreach my $file (keys %{$age->{$agecat}}) {
# 			foreach my $type (sort keys %{$age->{$agecat}{$file}}) {
# 				my $per = $age->{$agecat}{$file}{$type}/$totlen->{$file}{$type};
# 				print $fh "$age->{$agecat}{$file}{$type}\t";
# 				print $fh "$per\t";
# 			}
# 		}
# 		print $fh "\n";
# 	}	
# 	AGECAT: foreach my $agecat (sort keys %{$age}) {
# 		next AGECAT if (($agecat eq "LineageSpe") || ($agecat eq "Ancient"));
# 		print $fh "$agecat\t";
# 		foreach my $file (keys %{$age->{$agecat}}) {
# 			foreach my $type (sort keys %{$age->{$agecat}{$file}}) {
# 				my $per = $age->{$agecat}{$file}{$type}/$totlen->{$file}{$type};
# 				print $fh "$age->{$agecat}{$file}{$type}\t";
# 				print $fh "$per\t";
# 			}
# 		}
# 		print $fh "\n";
# 	}		
	close $fh;
	return ($out);
}		
	
#-----------------------------------------------------------------------------
# populate down so that there are 0 values for empty ones
# $age = populate_down3($age);
#-----------------------------------------------------------------------------
sub populate_down3 {
	my $hash = shift; 
	#get lists 
	my (@keys1,@keys2,@keys3) = ();	
	foreach my $key1 (keys %{$hash}) {
		push(@keys1,$key1);
		foreach my $key2 (keys %{$hash->{$key1}}) {
			push(@keys2,$key2);
			foreach my $key3 (keys %{$hash->{$key1}{$key2}}) {
				push(@keys3,$key3);
			}
		}
	}
	#Now fill up
	foreach my $key1 (@keys1) {
		KEY2: foreach my $key2 (@keys2) {
			foreach my $key3 (@keys3) {
				$hash->{$key1}{$key2}{$key3}=0 unless ($hash->{$key1}{$key2}{$key3});
			}
		}
	}				
	return ($hash);
}   

#-----------------------------------------------------------------------------
# get subset = filter, parsed, from TE analysis pipeline
# my ($fname,$subset) = get_subset($file);
#-----------------------------------------------------------------------------
sub get_subset {
	my $file = shift;
	my @file = split(/\./,$file);
	my $parsed = pop(@file);
	my $filter = pop(@file);
	my $subset = $filter.".".$parsed;	
	my $fname = join('.',@file);
	return ($fname,$subset);
}   

				

