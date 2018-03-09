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

my $VERSION = "1.1";
my $SCRIPTNAME = "TE-analysis_Shuffle_get_big_table.pl";
my $CHANGELOG = "
# Change log for $SCRIPTNAME:
#	v1.0 = 29 Nov 2016
#	v1.1 = 25 Jan 2018
#          Make more outputs
#          upercase norm
\n";

my $USAGE = "\nUsage [$VERSION]: 
    perl $SCRIPTNAME -d directory -n bootstraps [-c] [-h] [-help]
   
    SYNOPSIS
     Util for TE-analysis_Shuffle_XX.pl (https://github.com/4ureliek/TEanalysis/tree/master/Utils)
     Use to merge results of several runs (useful to compare conditions)
    
    CITATION
     - the GitHub link to this script

    MANDATORY ARGUMENTS:
     -d,--dir   => (STRING) directory with the subdirectories containing
                            a bunch of the *.stats.tab.txt outputs TE-analysis_Shuffle_XX.pl
     -n,--nboot => (INT)    number of bootstraps
	
    OPTIONAL ARGUMENTS:
     -m,--more  => (BOOL)   print more outputs: a matrix of just the pvalues, 
                            and if the TE is enriched or depleted
                            as well as a file with odds ratios (observed/expected average)
     -c,--chlog => (BOOL)   print log of changes
     -v,--v     => (BOOL)   print the version
     -h,--help  => (BOOL)   print this usage
\n";

#-----------------------------------------------------------------------------
#------------------------------ LOAD AND CHECK -------------------------------
#-----------------------------------------------------------------------------
my ($DIR,$NB,$MORE,$CHLOG,$V,$HELP);
my $opt_success = GetOptions(
			 	  'dir=s'		=> \$DIR,
			 	  'nboot=s'		=> \$NB,
			 	  'more'		=> \$MORE,
			 	  'chlog'       => \$CHLOG,
			 	  'version'     => \$V,
			 	  'help'		=> \$HELP,);
			 	  
#Check options, if files exist, etc
die "\n --- $SCRIPTNAME version $VERSION\n\n" if $V;
die $CHANGELOG if $CHLOG;
die $USAGE if $HELP || ! $opt_success;
die $USAGE unless ($DIR);
die $USAGE unless ($NB);
die "\n -d $DIR does not exist?\n\n"  if (! -e $DIR);
$DIR =~ s/\/$// if ($DIR =~ /\/$/); #remove the / at the end of directory if any

#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------
print STDERR "\n --- $SCRIPTNAME v$VERSION\n";

#Loading any summary files
print STDERR " --- Loading *.stats.tab.txt files\n";
my %STATS = ();
my %INFOS = ();
load_stats();

#Print data now
print STDERR " --- Print data in tables\n";
my $out = print_tables();

print STDERR " --- $SCRIPTNAME done\n";
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
#    my @statfiles = `ls $DIR/*/*.stats.tab.txt`;
	my @statfiles = `ls $DIR/*.stats.tab.txt`;
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
    		$INFOS{$fname}{'nb'} = $l[1] if ($l =~ /^\t/);
    		next LINE if ($l =~ /^\t/);
    		#Finish loading infos, and removing from the columns
     		$INFOS{$fname}{'exp_tot'}=splice(@l,10,1);
     		$INFOS{$fname}{'obs_tot'}=splice(@l,5,1);     			
    		#Now keep evrything by repeat and filename
    		my @name = splice(@l,0,3);   
    		my $id = join("\t",@name);
    		push(@ids,$id);
    		my $data = join("\t",@l);
    		$STATS{$id}{$fname} = $data;
    		push(@files,$fname);
    	}
    	close($fh);
    }
    #data is loaded, but need to fill up for all files
    fillup_stats(\@files,\@ids);   
    return 1;
}

#-----------------------------------------------------------------------------
# need to fill up data frame
#-----------------------------------------------------------------------------
sub fillup_stats {
	my ($files,$ids) = @_;	
	foreach my $id (@{$ids}) {
		foreach my $fname (@{$files}) {
			$STATS{$id}{$fname} = "0\t0\t0\t0\t0\t0\tna\tna\tna\tna\tna\tna\tna" unless ($STATS{$id}{$fname});
		}
	}
	return 1;
}   

#-----------------------------------------------------------------------------
# Now print stuff in a way that will make it easy to compare things
#-----------------------------------------------------------------------------
sub print_tables {
	#out file
	$DIR = $1 if $DIR =~ /(.*)\/$/;
	$DIR = "_" if ($DIR eq "."); #if current dir, would make a ..tab file
	my $out = $DIR.".tab";	
	prep_out($out,"full");
	open(my $fh, ">>",$out) or confess "\n   ERROR (sub print_tables): could not open to write $out $!\n";

	my ($outb,$outp,$outo,$outs,$fhb,$fhp,$fho,$fhs);
	if ($MORE) {
		$outb = $DIR.".binomial-pval.tab";
		$outp = $DIR.".permutation-pval.tab";
		$outo = $DIR.".odd_ratios.tab";
		$outs = $DIR.".simple.tab";
		prep_out($outb,"pval");
		prep_out($outp,"pval");
		prep_out($outo,"or");
		prep_out($outs,"simple");
		open($fhb, ">>",$outb) or confess "\n   ERROR (sub print_tables): could not open to write $outb $!\n";	
		open($fhp, ">>",$outp) or confess "\n   ERROR (sub print_tables): could not open to write $outp $!\n";
		open($fho, ">>",$outo) or confess "\n   ERROR (sub print_tables): could not open to write $outo $!\n";
		open($fhs, ">>",$outs) or confess "\n   ERROR (sub print_tables): could not open to write $outs $!\n";	
	}	
	foreach my $id (sort keys %STATS) {
		print $fh "$id\t";
		print $fhb "$id\t" if ($MORE);
		print $fhp "$id\t" if ($MORE);
		print $fho "$id\t" if ($MORE);
		print $fhs "$id\t" if ($MORE);
		foreach my $fname (sort keys %{$STATS{$id}}) {
			print $fh "$STATS{$id}{$fname}\t";
			if ($MORE) {
				my @data = split(/\t/,$STATS{$id}{$fname});
				my $rank = $data[6];
				my $pval_p = $data[7];
				my $pval_b = $data[11];
				if ($rank <= $NB/2) {
					print $fhb "$pval_b\tDEPLETED\t";
					print $fhp "$pval_p\tDEPLETED\t";
				} else { #elsif ($rank > $NB/2) {
					print $fhb "$pval_b\tENRICHED\t";
					print $fhp "$pval_p\tENRICHED\t";
				}
				my $or = $data[0]/$data[3]; #obs_hits/exp_avg_hits
				print $fho "$or\t";
				print $fhs "$data[0]\t$rank\t$data[8]\t$data[12]\t";
			}	
		}
		print $fh "\n";
		print $fhb "\n" if ($MORE);
		print $fhp "\n" if ($MORE);
		print $fho "\n" if ($MORE);
		print $fhs "\n" if ($MORE);
	}
	close $fh;
	close $fhb if ($MORE);
	close $fhp if ($MORE);
	close $fho if ($MORE);
	return ($out);
}		

#-----------------------------------------------------------------------------
# prep output file
#-----------------------------------------------------------------------------
sub prep_out {
	my ($out,$type) = @_;	
	#prep headers
	my $head1 = "COUNTS\t#\t#\t#\t#\t#\tPERMUTATION_TEST\t#\t#\tBINOMIAL_TEST\t#\t#\t#\t";
	my $head2 = "obs_hits\t%_obs_(%of_features)\tnb_of_trials(nb_of_TE_in_genome)\texp_avg_hits\texp_sd\t%_exp_(%of_features)\tobs_rank_in_exp\t2-tailed_permutation-test_pvalue(obs.vs.exp)\tsignificance\tbinomal_test_proba\tbinomial_test_95%_confidence_interval\tbinomial_test_pval\tsignificance\t";
	my $head2_b = "obs_hits\tobs_rank_in_exp\tperm_sign\tbinom_sign\t";
	
	#Now print all
	open(my $fh, ">",$out) or confess "\n   ERROR (sub prep_out): could not open to write $out $!\n";
	print $fh "#\t#\t#\t";
	foreach my $fname (sort keys %INFOS) {
		print $fh "$fname\t" if ($type eq "or");
		print $fh "$fname\t#\t" if ($type eq "pval");
		print $fh "$fname\t#\t#\t#\t" if ($type eq "simple");
		print $fh "$fname\tnb_features=\t$INFOS{$fname}{'nb'}\t" if ($type eq "full");
		print $fh "Obs_tot_hit=\t$INFOS{$fname}{'obs_tot'}\tExp_tot_hits(avg)=\t$INFOS{$fname}{'exp_tot'}\t" if ($type eq "full");
		print $fh "#\t#\t#\t#\t#\t#\t" if ($type eq "full");
	}
	if ($type eq "full") {
		print $fh "\n#Level_(tot_means_all)\t#\t#\t";
		foreach my $fname (sort keys %INFOS) {
			print $fh $head1;
		}
	}
	if ($type eq "full" || $type eq "simple") {	
		print $fh "\n#Rclass\tRfam\tRname\t";
		foreach my $fname (sort keys %INFOS) {
			if ($type eq "full") {
				print $fh $head2;
			} else {
				print $fh $head2_b;
			}			
		}
	}	
	print $fh "\n";	
	close $fh;
	return 1;
}   

				

