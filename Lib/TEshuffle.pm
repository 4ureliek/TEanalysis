#!/usr/bin/perl -w
#######################################################
# SUBROUTINES
# FOR SCRIPTS = TE-analysis_Shuffle_tr.pl & TE-analysis_Shuffle_bed.pl
# => defined as TEshuffle package
######################################################
package TEshuffle;
#use Data::Dumper;
use strict;
use warnings;
use Carp;

#----------------------------------------------------------------------------
# get a filename from a full path
# my $genone_name = DelGet::filename($genone);
#----------------------------------------------------------------------------
sub filename {
	my($name) = shift;
	$name =~ s/.*\/(.*)$/$1/;
	return $name;
}

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
# my ($toshuff_file,$parsedRM,$rm,$rm_c) = TEshuffle::RMtobed($shuffle,$okseq,$filter,$f_regexp,$nonTE,$age,"y",$stype);
# my ($toshuff_file,$parsedRM) = TEshuffle::RMtobed($shuffle,$okseq,$filter,$f_regexp,$nonTE,$age,$full);
#-----------------------------------------------------------------------------
sub RMtobed {
	my ($file,$okseq,$filter,$f_regexp,$nonTE,$age,$full,$stype) = @_;
	my $parsed = ();
	my $rm = ();
	my $rm_c = ();
	my ($f_type,$f_name) = split(",",$filter) unless ($filter eq "na");	
	my $ok = $file;
	$ok = $1 if ($file =~ /(.*)\.out$/);	
	$ok = $1 if ($file =~ /(.*)\.bed$/);
	if ($filter eq "na") {
		$ok = $ok.".nonTE-$nonTE" unless ($ok =~ /nonTE-$nonTE/);
	} else {
		$ok = $ok.".$f_name";
	}	
	$ok = $ok.".bed";
	if (-e $ok) { #if has been filtered same way, OK, just get dictionary
		print STDERR "     -> $ok exists, just getting dictionary from it\n";
		($parsed,$rm,$rm_c) = getparsedRM($ok,$parsed,"file",$age,$rm,$rm_c,$stype);
		return ($ok,$parsed,$rm,$rm_c);
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
		next LINE if (($stype eq "bed") && (! $okseq->{$l[4]})); #if not in build of stuff OK to shuffle on when stype = bed, remove here as well
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
		($parsed,$rm,$rm_c) = getparsedRM(\@l,$parsed,"line",$age,$rm,$rm_c,$stype) unless ($full eq "n");
	
	}
	close ($fh);
	close ($bed_fh);	
	return ($ok,$parsed,$rm,$rm_c);
}

#-----------------------------------------------------------------------------
# Get infos from RMout
# ($parsed,$rm,$rm_c) = getparsedRM($ok,$parsed,"file",$age,$rm,$rm_c,$stype);
# ($parsed,$rm,$rm_c) = getparsedRM(\@l,$parsed,"line",$age,$rm,$rm_c,$stype);
#-----------------------------------------------------------------------------
sub getparsedRM {
	my ($info,$parsed,$type,$age,$rm,$rm_c,$stype) = @_;
	if ($type eq "line") { #meaning it's read from the .out while being converted in .bed
		($parsed,$rm,$rm_c) = getparsedRMline($parsed,$info,$age,$rm,$rm_c,$stype);
	} else { #meaning it's already a bed file => open and read it
		open(my $fh, "<$info") or confess "\n   ERROR (sub getparsedRM): could not open to read $info!\n";
		LINE: while(<$fh>) {
			chomp(my $l = $_);
			$l =~ s/^\s+//;
			next LINE if (($l =~ /^[Ss]core|^SW|^#/) || ($l !~ /\w/));
			my @l = split('\s+',$l);
			my @RMline = split(";",$l[3]);
			($parsed,$rm,$rm_c) = getparsedRMline($parsed,\@RMline,$age,$rm,$rm_c,$stype);
		}
		close $fh;	
	}
	return ($parsed,$rm,$rm_c);
}

#-----------------------------------------------------------------------------
# Get counts and length for each TE
# ($parsed,$rm,$rm_c) = getparsedRMline($parsed,$info,$age,$rm,$rm_c,$stype)
# ($parsed,$rm,$rm_c) = getparsedRMline($parsed,\@RMline,$age,$rm,$rm_c,$stype)
#-----------------------------------------------------------------------------
sub getparsedRMline {
	my ($parsed,$l,$age,$rm,$rm_c,$stype) = @_;
	my ($Rname,$classfam) = ($l->[9],$l->[10]);
	$Rname =~ s/\//_/; #making sure no / in repeat name
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
		#cat2 tot = cat1 tot:
		$parsed->{'age'}{'cat.2'}{'tot'}=$parsed->{'age'}{'cat.1'}{'tot'};
	}
	#Load arrays to shuffle the TEs inside RMout, if $stype is rm
	if ($stype eq "rm") {
		my $middle = $l->[5] + int(($l->[6] - $l->[5]) / 2 + 0.5);
		push(@{$rm_c->{$l->[4]}},$middle); #=> I can shuffle this, and then associate these coordinates with the other one => reconstruct a bed file, sort and then and intersect	
		push(@{$rm->{$l->[4]}},$l); #to keep all infos
	}
	return ($parsed,$rm,$rm_c);
}

#-----------------------------------------------------------------------------
# Load TSS from annotation gtf or gff
# ($tss,$alltss) = load_and_print_tss($tssfile);
#-----------------------------------------------------------------------------
sub load_and_print_tss {
	my $file = shift;
	my $tssbed = $file;
	$tssbed = $1 if ($file =~ /^(.*)\.(gtf|gff|gff3)$/);
	$tssbed = $tssbed.".tss";
	print STDERR "     $tssbed exists, skipping printing and load only\n" if (-e $tssbed);	
	my $alltss = load_all_tss($tssbed) if (-e $tssbed);	
	return ($tssbed,$alltss) if (-e $tssbed);	
	my $i = 0;
	my %alltss = ();
	open (my $in, '<', $file) or confess "ERROR (sub print_tss): can't open to read $file $!\n";
	open (my $out, '>', $tssbed) or confess "ERROR (sub print_tss): can't open to write $tssbed $!\n";
	LINE: while (<$in>) {
		chomp(my $l = $_);
		next LINE if (substr($l,0,1) eq "#");
		my @l = split(/\t/,$l);		 
		next LINE unless ($l[2] =~ /transcript/); #could be "processed transcript"...?
		my $trst;
		($l[6] eq "+")?($trst=$l[3]):($trst=$l[4]);
		my $st=$trst-1;
		my $en=$trst+1;
		my $tssid = "tss_".$i;
		$l[0] = "chr".$l[0] if ($l[0]=~/^[0-9]|^[MXYZW]$|^MT$/); #ensembl gtf only has chr number
		print $out "$l[0]\t$st\t$en\t$tssid\t$l[6]\n"; #may as well keep strand info
		$i++;		
		#Load the TSS in hash as well
		my @tss = ($trst,$l[6]);
		push(@{$alltss{$l[0]}},\@tss);
	}
	close $in;
	close $out;
	return ($tssbed,\%alltss);
}

#-----------------------------------------------------------------------------
# Load all TSS
# my $alltss = load_all_tss($tssbed);
#-----------------------------------------------------------------------------
sub load_all_tss {
	my $tssbed = shift;
	my %alltss = ();
	open(my $fh, "<$tssbed") or confess "\n   ERROR (sub load_all_tss): could not open to read $tssbed $!\n";
	LINE: while(<$fh>) {
		chomp(my $l = $_);
		my @l = split('\s+',$l);
		#chr1	3073252	3073254	tss_0	+
		my @tss = ($l[1]+1,$l[4]);
		push(@{$alltss{$l[0]}},\@tss);
	}
	return (\%alltss);
}

#-----------------------------------------------------------------------------
# Load distance to closes
# $closest = TEshuffle::load_closest_tss($tssclosest);
#-----------------------------------------------------------------------------
sub load_closest_tss {
	my $tssclosest = shift;
	my %closest = ();
	open(my $fh, "<$tssclosest") or confess "\n   ERROR (sub load_closest_tss): could not open to read $tssclosest $!\n";
	LINE: while(<$fh>) {
		chomp(my $l = $_);
		my @l = split('\s+',$l);
		#chr19	3000005	3000195	470;23.0;4.5;2.7;chr19;3000005;3000195;(58431371);+;IAPLTR3-int;LTR/ERVK;5531;5745;(1076);2062961	.	+	chr19	3137089	3137091	tss_116236	-	136895
		$closest{$l[0]}{$l[3]}=$l[-1]; #key = chr, then RM line, value = distance to TSS	
	}
	return (\%closest);
}

#-----------------------------------------------------------------------------
# Cleanup previous outputs
# my ($stats,$out,$outb,$temp) =                                TEshuffle::prep_out("bed",$dir,$nboot,$filter,$input,$stype,$nonTE);
# my ($stats,$outl,$outlb,$temp_l,$outp,$outpb,$temp_p,$temp) = TEshuffle::prep_out("tr", $dir,$nboot,$filter,$input,$stype,$nonTE,$linc,$prot,$shuffle);
#-----------------------------------------------------------------------------
sub prep_out {
	my ($type,$dir,$nboot,$filter,$input,$stype,$nonTE,$linc,$prot,$shuffle) = @_;	
	#common steps
	my ($f_type,$f_name) = split(",",$filter) unless ($filter eq "na");		
	my $inputname = filename($input);
	my $stats;
	($filter eq "na")?($stats = "$dir/$inputname.shuffle-$stype.nonTE-$nonTE.boot-$nboot.stats.tab"):
                      ($stats = "$dir/$inputname.shuffle-$stype.nonTE-$nonTE.boot-$nboot.$f_name.stats.tab");		
	`rm -Rf $dir.previous` if ((-e $dir) && (-e $dir.".previous"));
	`mv $dir $dir.previous` if (-e $dir);
	`mkdir $dir`;
	#bed type - fewer directories
	if ($type eq "bed") {
		my ($out,$outb,$temp) = ("$dir/$inputname.no_boot","$dir/$inputname.boot","$dir/$inputname.temp");
		`mkdir $temp` if ($nboot > 0);
		return ($stats,$out,$outb,$temp);
	} 
	#tr type - more directories
	if ($type eq "tr") {
		my $lincname = filename($linc);
		my $protname = filename($prot);
		my ($outl,$outlb,$temp_l) = ("$dir/$lincname.no_boot","$dir/$lincname.boot","$dir/$lincname.temp");
		my ($outp,$outpb,$temp_p) = ("$dir/$protname.no_boot","$dir/$protname.boot","$dir/$protname.temp");
		`mkdir $temp_l` unless ($linc eq "n"); #also used for the -m 
		`mkdir $temp_p` unless ($prot eq "n"); #also used for the -m 
		my $shufname = filename($shuffle);
		my $temp = "$dir/$shufname.temp";
		`mkdir $temp` if ($nboot > 0);
		return ($stats,$outl,$outlb,$temp_l,$outp,$outpb,$temp_p,$temp);
	}
	return 1;
}	

#-----------------------------------------------------------------------------
# Shuffle, but keep distance to closest TSS
# $shuffled = TEshuffle::shuffle_tss($toshuff_file,$temp,$i,$alltss,$closest) if ($stype eq "tss");
#-----------------------------------------------------------------------------
sub shuffle_tss {
	my ($toshuff_file,$temp,$nb,$alltss,$closest) = @_;
	my $out = $temp."/shufffled".$nb;	
	open (my $fh, '>', $out) or confess "\n   ERROR (sub shuffle_tss): can't open to write $out $!\n";
	CHR: foreach my $chr (sort keys %{$closest}) {
		warn "     WARN: $chr not found in the -a file\n" if (! $alltss->{$chr});
		next CHR if (! $alltss->{$chr});
		fisher_yates_shuffle($alltss->{$chr}); # permutes @array in place; one array per chr => keep chr distribution
		my $i = 0;
		foreach my $te (keys %{$closest->{$chr}}) {
			$i = 0 if (! $alltss->{$chr}[$i]); #basically, if end of the TSS for that chromosome, restart
			my $tss = $alltss->{$chr}[$i];	
			#get start and end, depends on strands
			my @te = split(";",$te);
			my $len = $te[6] - $te[5];
			my ($st,$en);
			my $dist = $closest->{$chr}{$te};
			if ($tss->[1] eq "+") {
				($dist > 0)?($st = $tss->[0] + $dist):($en = $tss->[0] - abs($dist)); #or + $dist
				($dist > 0)?($en = $st + $len):($st = $en - $len);		
			} else { 
				($dist > 0)?($en = $tss->[0] - $dist):($st = $tss->[0] + abs($dist)); #or - $dist
				($dist > 0)?($st = $en - $len):($en = $st + $len);		
			}
			$st = 1 if ($st <0); #added v4.2
			print $fh "$chr\t$st\t$en\t$te\t.\t$te[8]\n"; #shuffled TE, same chromosome, with its distance to TSS maintained
			$i++;
		}
	}
	close ($fh);
	return ($out);
}

#-----------------------------------------------------------------------------
# Shuffle, but positions instead => keek current TE distribution
# $shuffled = TEshuffle::shuffle_rm($toshuff_file,$temp,$i,$rm,$rm_c,$okseq) if ($stype eq "rm");	
#-----------------------------------------------------------------------------
sub shuffle_rm {
	my ($toshuff_file,$temp,$nb,$rm,$rm_c,$okseq) = @_;
	my $out = $temp."/shufffled".$nb;
	open (my $fh, '>', $out) or confess "\n   ERROR (sub shuffle_rm): can't open to write $out $!\n";
	CHR: foreach my $chr (sort keys %{$rm_c}) {
		warn "     WARN: $chr not found in the -a file\n" if (! $rm_c->{$chr});
		next CHR if (! $rm_c->{$chr});
		fisher_yates_shuffle($rm_c->{$chr}) if ($rm_c->{$chr}[1]); # permutes @array in place; one array per chr => keep chr distribution
		my $i = 0;
		foreach my $te (@{$rm->{$chr}}) {			
			my $l = join(";",@{$te}); #=> reconstitute the bed id
			my $middle = $rm_c->{$chr}[$i];			
			my $midlen = int(($te->[6] - $te->[5]) / 2 + 0.5);
			my $st = $middle - $midlen;
			my $en = $middle + $midlen;
			if ($st < 0) {
				my $shift = -$st; #put back in positive
				$st = $st+$shift+1;
				$en = $en+$shift+1;	
			}
			if ($okseq->{$chr} && $en > $okseq->{$chr}) {
				my $shift = $en - $okseq->{$chr};
				$st = $st-$shift;
				$en = $okseq->{$chr};
			}
			print $fh "$chr\t$st\t$en\t$l\t.\t$te->[8]\n"; 
			$i++;
		}	
	}
	close $fh;
	return ($out);
}

#-----------------------------------------------------------------------------
# fisher_yates_shuffle
#-----------------------------------------------------------------------------
sub fisher_yates_shuffle { #http://www.perlmonks.org/?node_id=1869
    my $array = shift;
    my $i = @$array;
    while ( --$i ) {
        my $j = int rand( $i+1 );
        @$array[$i,$j] = @$array[$j,$i];
    }
}

#-----------------------------------------------------------------------------
# Shuffle, with bedtools
# my $shuffled = TEshuffle::shuffle_bed($toshuff_file,$temp_s,$i,$excl,$incl,$build_file,$bedtools,$nooverlaps);
#-----------------------------------------------------------------------------
sub shuffle_bed {
	my ($toshuff_file,$temp_s,$nb,$excl,$incl,$build,$bedtools,$nooverlaps) = @_;
	my $out = $temp_s."/shufffled".$nb;
	my $bed_shuffle = $bedtools."shuffleBed";
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
	VALA: foreach (@$data) {
		my $val = $_;
		print STDERR "UNDEF in sub get_avg_and_sd\n" unless ($_);
#		print Dumper($data) unless ($_);
#		next VALA unless ($_);
		$total += $val;
	}
	my $avg = $total / @$data;

	my $sqtotal = 0;
	VALM: foreach(@$data) {
#		next VALM unless ($_);
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
# $exp = TEshuffle::binomial_test_R($exp,"bed");
# $te_exp = TEshuffle::binomial_test_R($te_exp,"tr_rep");
#-----------------------------------------------------------------------------
sub binomial_test_R {
	my ($exp,$type) = @_;
	#Start R bridge
	my $R = Statistics::R->new() ;
	$R->startR;
	
	#2 types of hash structures
	if ($type eq "bed") {
		foreach my $k1 (keys %{$exp}) { 	
			foreach my $k2 (keys %{$exp->{$k1}}) {			
				KEY: foreach my $k3 (keys %{$exp->{$k1}{$k2}}) {
					my $x = $exp->{$k1}{$k2}{$k3}{'x'};
					my $n = $exp->{$k1}{$k2}{$k3}{'n'};
					my $p = $exp->{$k1}{$k2}{$k3}{'p'};
					#set values
					$exp->{$k1}{$k2}{$k3}{'binom_prob'} = "na"; 
					$exp->{$k1}{$k2}{$k3}{'binom_conf'} = "na";
					$exp->{$k1}{$k2}{$k3}{'binom_pval'} = "na";					
					#check now
					my $skip = check_binom_values($k1,$k2,$k3,$x,$p,$n);
					next KEY if ($skip eq "y");
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
							my $x = $exp->{$cat}{$t}{$k1}{$k2}{$k3}{'x'};
							my $n = $exp->{$cat}{$t}{$k1}{$k2}{$k3}{'n'};
							my $p = $exp->{$cat}{$t}{$k1}{$k2}{$k3}{'p'};
							#set values
							$exp->{$cat}{$t}{$k1}{$k2}{$k3}{'binom_prob'} = "na"; 
							$exp->{$cat}{$t}{$k1}{$k2}{$k3}{'binom_conf'} = "na";
							$exp->{$cat}{$t}{$k1}{$k2}{$k3}{'binom_pval'} = "na";
							#check now
							my $skip = check_binom_values($k1,$k2,$k3,$x,$p,$n);
							next KEY if ($skip eq "y");
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

#-----------------------------------------------------------------------------
# check the x, n and p values & print warnings
#-----------------------------------------------------------------------------
sub check_binom_values {
	my ($k1,$k2,$k3,$x,$p,$n) = @_;
	my $skip = "n";
	return("y") if ($p == 0);
	if (! $n || $n == 0) {
		print STDERR "        WARN: no value for total number (from parsed RM table), for {$k1}{$k2}{$k3}? => no binomial test\n";
		$skip = "y";
	}
	if ($p > 1) {
		print STDERR "        WARN: expected overlaps (avg) > total number (from parsed RM table) for {$k1}{$k2}{$k3}, likely due to multiple peaks overlapping with the repeat => no binomial test\n";
		$skip = "y";
	}
	if ($n < $x) {
		print STDERR "        WARN: we should have n > x but n = $n and x = $x for {$k1}{$k2}{$k3}, likely due to multiple peaks overlapping with the repeat => no binomial test\n";
		$skip = "y";
	}
	return $skip;
}

#ensure last returned value is true (to load as require in scripts)
1;