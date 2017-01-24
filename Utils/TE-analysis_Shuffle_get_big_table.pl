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
my $scriptname = "TE-analysis_Shuffle_get_big_table.pl";
my $changelog = "
# Change log for $scriptname:
#	v1.0 = 29 Nov 2016
\n";

my $usage = "\nUsage [$version]: 
    perl $scriptname -d directory [-b] [-chlog] [-h] [-help]
   
    SYNOPSIS
     Util for TE-analysis_Shuffle_XX.pl (https://github.com/4ureliek/TEanalysis/tree/master/Utils)
     Use to merge results of several runs (useful to compare conditions)
    
    CITATION
     - the GitHub link to this script

    MANDATORY ARGUMENTS:
     -d,--dir   => (STRING) directory with the subdirectories containing
                            a bunch of the *.stats.tab.txt outputs TE-analysis_Shuffle_XX.pl
	
    OPTIONAL ARGUMENTS:
     -c,--chlog => (BOOL)   print log of changes
     -v,--v     => (BOOL)   print the version
     -h,--help  => (BOOL)   print this usage
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
$dir =~ s/\/$// if ($dir =~ /\/$/); #remove the / at the end of directory if any

#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------
print STDERR "\n --- $scriptname v$version\n";

#Loading any summary files
print STDERR " --- Loading *.stats.tab.txt files\n";
my ($stats,$infos) = load_stats($dir);

#Print data now
print STDERR " --- Print data in tables\n";
my $out = print_tables($dir,$stats,$infos);

print STDERR " --- $scriptname done\n";
print STDERR "      => $out\n\n";

exit;

#-----------------------------------------------------------------------------
#-------------------------------- SUBROUTINES --------------------------------
#-----------------------------------------------------------------------------
#----------------------------------------------------------------------------
# get a filename from a full path
#----------------------------------------------------------------------------
sub filename {
	my($file) = shift;
	$file =~ s/.*\/(.*)$/$1/;
	return $file;
}

#-----------------------------------------------------------------------------
# Load all the stats files
#-----------------------------------------------------------------------------
sub load_stats{ 
    my ($dir) = @_;   
    my @statfiles = `ls $dir/*/*.stats.tab.txt`;
    my $stats = ();
    my $infos = ();
    my @files = ();
    my @ids = ();
    foreach my $file (@statfiles) {
    	chomp($file);
    	open(my $fh, "<$file") or confess "\n   ERROR (sub load_stats): could not open to read $file $!\n"; 
    	LINE: while (<$fh>) {
    		chomp(my $l = $_);
    		next LINE if ((substr($l,0,1) eq "#") || ($l !~ /\w/)); 
    		my @l = split(/\t/,$l);    		
    		my $fname = filename($file);
    		$fname =~ s/\.stats\.tab\.txt$//;
    		$infos->{$fname}{'nb'} = $l[1] if ($l =~ /^\t/);
    		next LINE if ($l =~ /^\t/);
    		#Finish loading infos, and removing from the columns
     		$infos->{$fname}{'exp_tot'}=splice(@l,10,1);
     		$infos->{$fname}{'obs_tot'}=splice(@l,5,1);     			
    		#Now keep evrything by repeat and filename
    		my @name = splice(@l,0,3);   
    		my $id = join("\t",@name);
    		push(@ids,$id);
    		my $data = join("\t",@l);
    		$stats->{$id}{$fname} = $data;
    		push(@files,$fname);
    	}
    	close($fh);
    }
    #data is loaded, but need to fill up for all files
    $stats = fillup_stats($stats,\@files,\@ids);   
    return ($stats,$infos);
}

#-----------------------------------------------------------------------------
# need to fill up data frame
#-----------------------------------------------------------------------------
sub fillup_stats {
	my ($stats,$files,$ids) = @_;	
	foreach my $id (@{$ids}) {
		foreach my $fname (@{$files}) {
			$stats->{$id}{$fname} = "0\t0\t0\t0\t0\t0\tna\tna\tna\tna\tna\tna\tna" unless ($stats->{$id}{$fname});
		}
	}
	return $stats;
}   

#-----------------------------------------------------------------------------
# Now print stuff in a way that will make it easy to compare things
#-----------------------------------------------------------------------------
sub print_tables {
	my ($dir,$stats,$infos) = @_;
	#out file
	$dir = $1 if $dir =~ /(.*)\/$/;
	$dir = "_" if ($dir eq "."); #if current dir, would make a ..tab file
	my $out = $dir.".tab";
	prep_out($out,$infos);
	
	#now print data
	open(my $fh, ">>",$out) or confess "\n   ERROR (sub print_tables): could not open to write $out $!\n";	
	foreach my $id (keys %{$stats}) {
		print $fh "$id\t";
		foreach my $fname (keys %{$stats->{$id}}) {
			print STDERR "$id - $fname => data = $stats->{$id}{$fname}\n" if $id =~ /RLTR11D/;
			print $fh "$stats->{$id}{$fname}\t";
		}
		print $fh "\n";
	}
	close $fh;
	return ($out);
}		

#-----------------------------------------------------------------------------
# prep output file
#-----------------------------------------------------------------------------
sub prep_out {
	my ($out,$infos) = @_;	
	#prep headers
	my $head1 = "COUNTS\t#\t#\t#\t#\t#\tPERMUTATION_TEST\t#\t#\tBINOMIAL_TEST\t#\t#\t#\t";
	my $head2 = "obs_hits\t%_obs_(%of_features)\tnb_of_trials(nb_of_TE_in_genome)\texp_avg_hits\texp_sd\t%_exp_(%of_features)\tobs_rank_in_exp\t2-tailed_permutation-test_pvalue(obs.vs.exp)\tsignificance\tbinomal_test_proba\tbinomial_test_95%_confidence_interval\tbinomial_test_pval\tsignificance\t";
	
	#Now print all
	open(my $fh, ">",$out) or confess "\n   ERROR (sub prep_out): could not open to write $out $!\n";
	print $fh "#\t#\t#\t";
	foreach my $fname (keys %{$infos}) {	
		print $fh "$fname\tnb_features=\t$infos->{$fname}{'nb'}\t";
		print $fh "Obs_tot_hit=\t$infos->{$fname}{'obs_tot'}\tExp_tot_hits(avg)=\t$infos->{$fname}{'exp_tot'}\t";
		print $fh "#\t#\t#\t#\t#\t#\t";		
	}
	print $fh "\n#Level_(tot_means_all)\t#\t#\t";	
	foreach my $fname (keys %{$infos}) {
		print $fh $head1;
	}
	print $fh "\n#Rclass\tRfam\tRname\t";
	foreach my $fname (keys %{$infos}) {
		print $fh $head2;
	}
	print $fh "\n";	
	close $fh;
	return 1;
}   

				

