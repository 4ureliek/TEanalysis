#!/usr/bin/perl -w
#######################################################
# SUBROUTINES
# FOR SCRIPTS = TE-analysis_Shuffle_tr.pl & TE-analysis_Shuffle_bed.pl
# => defined as TEshuffle package
######################################################
package TEshuffle;
use Carp;

#-----------------------------------------------------------------------------
# Get build if needed and get chromosomes included in it
# my ($okseq,$build_file) = TEshuffle::load_build($build,$dobuild);
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
# my $excl = TEshuffle::load_gap($gaps,$dogaps);
#-----------------------------------------------------------------------------
sub load_gap {
    my ($file,$dogaps) = @_;
    my $bed;  
	if ($dogaps eq "y") { 
	    #gaps not in file yet, get them from fasta file in bed format
		$bed = print_gap_bed($file);
	} else {
		print STDERR "     printing bed file corresponding to assembly gaps in $file (limiting to features > 50nt)\n";
		#file with gaps, check if bed already
		($file =~ /\.bed$/)?($bed = $file):($bed = $file.".bed");
		if ($file !~ /\.bed$/) {
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
    return ($bed);
}

#-----------------------------------------------------------------------------
# printing assembly gaps if needed, in bed format
# $file = TEshuffle::print_gap($file) if ($ifbuild eq "n");
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
# $incl = TEshuffle::concat_beds(\@include);
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
# Load TE age file if any
# $age = TEshuffle::load_TEage($TEage,$v) unless ($TEage eq "na");
#-----------------------------------------------------------------------------
sub load_TEage {
	my $in = shift;
	my %TEs = ();
	open(my $in_fh, "<", $in) or confess "\n   ERROR (sub load_TEage): could not open to read $in $!\n";
	LINE: while(<$in_fh>) {
		chomp (my $line = $_);
		next LINE if (($line !~ /\w/) || (substr($line,0,5) eq "Rname") || (substr($line,0,1) eq "#"));
		my @TEs = split('\t', $line); 		
		my $Rname = shift @TEs; #Rname not in values now, will be the key
		$Rname =~ s/\//_/; #making sure no / in repeat name
		$TEs[1] =~ s/\//_/; #making sure no / in family
		$TEs[1] = $TEs[1]."--int" if (($TEs[1] =~ /ERV/) && (($Rname =~ /[-_][iI]nt/) || ($Rname =~ /[-_]I$/)));	
		$TEs[2] = $TEs[0]."/".$TEs[1]; #correct class/fam
		$TEs{$Rname} = \@TEs;
	}	
	close $in_fh;
	return (\%TEs);
}

#-----------------------------------------------------------------------------
# Convert RMoutput .out file to bed if needed
# my ($toshuff_file,$parsedRM) = TEshuffle::RMtobed($shuffle,$okseq,$filter,$f_regexp,$nonTE,$age,"y");
# my ($toshuff_file,$parsedRM) = TEshuffle::RMtobed($shuffle,$okseq,$filter,$f_regexp,$nonTE,$age,$full);
#-----------------------------------------------------------------------------
sub RMtobed {
	my ($file,$okseq,$filter,$f_regexp,$nonTE,$age,$full) = @_;
	my $parsed = ();
	my ($f_type,$f_name) = split(",",$filter) unless ($filter eq "na");	
	my $ok = $file;
	$ok = $1 if ($file =~ /(.*)\.out$/);	
	$ok = $1 if ($file =~ /(.*)\.bed$/);	
	($filter eq "na")?($ok = $ok.".nonTE-$nonTE.bed"):($ok = $ok.".$f_name.bed");		
	if (-e $ok) { #if has been filtered same way, OK, just get dictionary
		print STDERR "     -> $ok exists, just getting dictionary from it\n";
		$parsed = getparsedRM($ok,$parsed,"file",$age);
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
		$parsed = getparsedRM(\@l,$parsed,"line",$age) unless ($full eq "n");
	}
	close ($fh);
	close ($bed_fh);
	return ($ok,$parsed);
}

#-----------------------------------------------------------------------------
# Get infos from RMout
# $parsed = getparsedRM($file,$parsed,"file",$age);
# $parsed = getparsedRM(\@l,$parsed,"line",$age);
#-----------------------------------------------------------------------------
sub getparsedRM {
	my ($info,$parsed,$type,$age) = @_;
	if ($type eq "line") { #meaning it's read from the .out while being converted in .bed
		$parsed = getparsedRMline($parsed,$info,$age);
	} else { #meaning it's already a bed file => open and read it
		open(my $fh, "<$info") or confess "\n   ERROR (sub getparsedRM): could not open to read $info!\n";
		LINE: while(<$fh>) {
			chomp(my $l = $_);
			$l =~ s/^\s+//;
			next LINE if (($l =~ /^[Ss]core|^SW|^#/) || ($l !~ /\w/));
			my @l = split('\s+',$l);
			my @RMline = split(";",$l[3]);
			$parsed = getparsedRMline($parsed,\@RMline,$age);
		}
		close $fh;	
	}
	return($parsed);
}

#-----------------------------------------------------------------------------
# Get counts and length for each TE
# $parsed = getparsedRMline($parsed,$info,$age)
# $parsed = getparsedRMline($parsed,\@RMline,$age)
#-----------------------------------------------------------------------------
sub getparsedRMline {
	my ($parsed,$l,$age) = @_;
	my ($Rname,$classfam) = ($l->[9],$l->[10]);
	my ($Rclass,$Rfam) = get_Rclass_Rfam($Rname,$classfam);
	#now feed dictionary
	($parsed->{'tot'}{'tot'}{'tot'})?($parsed->{'tot'}{'tot'}{'tot'}++):($parsed->{'tot'}{'tot'}{'tot'}=1);
	($parsed->{$Rclass}{'tot'}{'tot'})?($parsed->{$Rclass}{'tot'}{'tot'}++):($parsed->{$Rclass}{'tot'}{'tot'}=1);
	($parsed->{$Rclass}{$Rfam}{'tot'})?($parsed->{$Rclass}{$Rfam}{'tot'}++):($parsed->{$Rclass}{$Rfam}{'tot'}=1);
	($parsed->{$Rclass}{$Rfam}{$Rname})?($parsed->{$Rclass}{$Rfam}{$Rname}++):($parsed->{$Rclass}{$Rfam}{$Rname}=1);
	#Also feed age if relevant
	if ($age->{$Rname}) {
		($parsed->{'age'}{'cat.1'}{'tot'})?($parsed->{'age'}{'cat.1'}{'tot'}++):($parsed->{'age'}{'cat.1'}{'tot'}=1);
		($parsed->{'age'}{'cat.1'}{$age->{$Rname}[4]})?($parsed->{'age'}{'cat.1'}{$age->{$Rname}[4]}++):($parsed->{'age'}{'cat.1'}{$age->{$Rname}[4]}=1);
		($parsed->{'age'}{'cat.2'}{$age->{$Rname}[5]})?($parsed->{'age'}{'cat.2'}{$age->{$Rname}[5]}++):($parsed->{'age'}{'cat.2'}{$age->{$Rname}[5]}=1);	
	}	
	$parsed->{'age'}{'cat.2'}{'tot'}=$parsed->{'age'}{'cat.1'}{'tot'};
	return $parsed;
}

#-----------------------------------------------------------------------------
# Shuffle, with bedtools
# my $shuffled = TEshuffle::shuffle($toshuff_file,$temp_s,$i,$excl,$incl,$build_file,$bedtools,$nooverlaps);
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

#----------------------------------------------------------------------------
# Get class fam
# my ($Rclass,$Rfam) = get_Rclass_Rfam($Rname,$classfam);
# my ($Rcla,$Rfam) = TEshuffle::get_Rclass_Rfam($Rnam,$rm[10]);
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
# sub get_avg_and_sd
#-----------------------------------------------------------------------------	
sub get_avg_and_sd{
	my($data) = @_;
    warn "WARN: (sub print_stats/get_avg_and_sd): empty data array $!\n" if (not @$data);
    return ($data->[0],"na") if (@$data == 1); #returning 0,0 was an issue if -m =1
        
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
# sub get_sign
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

#-----------------------------------------------------------------------------
# sub binomial_test_R
# $exp = binomial_test_R($exp);
#-----------------------------------------------------------------------------
sub binomial_test_R {
	my ($exp,$type) = @_;
	#Start R bridge
	my $R = Statistics::R->new() ;
	$R->startR ;
	
	#2 types of hash structures
	if ($type eq "bed") {
		foreach my $k1 (keys %{$exp}) { 	
			foreach my $k2 (keys %{$exp->{$k1}}) {			
				KEY: foreach my $k3 (keys %{$exp->{$k1}{$k2}}) {
					my ($x,$n,$p) = ($exp->{$k1}{$k2}{$k3}{'x'},$exp->{$k1}{$k2}{$k3}{'n'},$exp->{$k1}{$k2}{$k3}{'p'});
					#set values
					($exp->{$k1}{$k2}{$k3}{'binom_prob'},$exp->{$k1}{$k2}{$k3}{'binom_conf'},$exp->{$k1}{$k2}{$k3}{'binom_pval'}) = ("na","na","na");
					print STDERR "        WARN: no value for total number (from parsed RM table), for {$k1}{$k2}{$k3}? => no binomial test\n" if ($n == 0);
					next KEY if ($n == 0);								
					print STDERR "        WARN: expected overlaps (avg) > total number (from parsed RM table), for {$k1}{$k2}{$k3}, likely due to multiple peaks overlapping with the repeat => no binomial test\n" if ($p > 1);				
					next KEY if ($p > 1);								
					$R->send(qq`d<-binom.test($x,$n,$p,alternative="two.sided")`);				
					my $data = $R->get('d');		
					my @data = @{$data};
					DATA: for (my $i=0; $i<$#data; $i++) {		
						$exp->{$k1}{$k2}{$k3}{'binom_pval'} = $data->[$i+2] if ($data->[$i] eq "p-value");
						$exp->{$k1}{$k2}{$k3}{'binom_conf'} = "$data->[$i+2];$data->[$i+3]" if ($data->[$i] eq "confidence");
					}
					$exp->{$k1}{$k2}{$k3}{'binom_prob'} = $data->[-1];
				}
			}
		}
	} elsif ($type eq "tr") {
		foreach my $cat (keys %{$exp}) {
			foreach my $t (keys %{$exp->{$cat}}) { 	
				foreach my $k1 (keys %{$exp->{$cat}{$t}}) { 	
					foreach my $k2 (keys %{$exp->{$cat}{$t}{$k1}}) {			
						KEY: foreach my $k3 (keys %{$exp->{$cat}{$t}{$k1}{$k2}}) {
							my ($x,$n,$p) = ($exp->{$cat}{$t}{$k1}{$k2}{$k3}{'x'},$exp->{$cat}{$t}{$k1}{$k2}{$k3}{'n'},$exp->{$cat}{$t}{$k1}{$k2}{$k3}{'p'});
							#set values
							($exp->{$cat}{$t}{$k1}{$k2}{$k3}{'binom_prob'},$exp->{$cat}{$t}{$k1}{$k2}{$k3}{'binom_conf'},$exp->{$cat}{$t}{$k1}{$k2}{$k3}{'binom_pval'}) = ("na","na","na");
							print STDERR "        WARN: no value for total number (from parsed RM table), for {$k1}{$k2}{$k3}? => no binomial test\n" if ($n == 0);
							next KEY if ($n == 0);								
							print STDERR "        WARN: expected overlaps (avg) > total number (from parsed RM table), for {$k1}{$k2}{$k3}, likely due to multiple peaks overlapping with the repeat => no binomial test\n" if ($p > 1);				
							next KEY if ($p > 1);								
							$R->send(qq`d<-binom.test($x,$n,$p,alternative="two.sided")`);				
							my $data = $R->get('d');		
							my @data = @{$data};
							DATA: for (my $i=0; $i<$#data; $i++) {		
								$exp->{$cat}{$t}{$k1}{$k2}{$k3}{'binom_pval'} = $data->[$i+2] if ($data->[$i] eq "p-value");
								$exp->{$cat}{$t}{$k1}{$k2}{$k3}{'binom_conf'} = "$data->[$i+2];$data->[$i+3]" if ($data->[$i] eq "confidence");
							}
							$exp->{$cat}{$t}{$k1}{$k2}{$k3}{'binom_prob'} = $data->[-1];
						}
					}
				}
			}
		}		
	}	
	$R->stopR() ;
	return($exp);
}


#ensure last returned value is true (to load as require in scripts)
1;