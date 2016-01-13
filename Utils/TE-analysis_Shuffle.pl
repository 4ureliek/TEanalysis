#!/usr/bin/perl
#######################################################
# Author  :  Zev Kronenberg (https://github.com/zeeev), for the v1.0
#            Modifications by Aurelie Kapusta (https://github.com/4ureliek) after v2.0
# email   :  4urelie.k@gmail.com  
# Purpose :  Writen to generate data (observed vs expected) shown in Figure 5, 
#            Kapusta et al. 2013 PLoS Genetics
#            (http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1003470)
#######################################################
use strict;
use warnings;
use Carp;
use Getopt::Long;
use Set::IntervalTree;
use GAL::Annotation;
use POSIX;
use Set::IntSpan::Fast;
use Bio::SeqIO;
use List::Util 'shuffle';

#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------
#flush buffer
$| = 1;

my $version = "2.0";
my $scriptname = "TE-analysis_Shuffle.pl";
my $changelog = "
#	- v1.0 = 2012
#	- v2.0 = Jan 11 2016
#		Too many changes to list here, but concept remains the same
\n";
my $usage = "
Synopsis (v$version):

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
    -t,--text     => (BOOL)   To get examples of command lines to use in R to process the outputs
    -v,--version  => (BOOL)   print the version
    -h,--help     => (BOOL)   print this usage
\n";

my ($prot,$lincs,$gap_file,$shuffle,$build,$ifbuild,$ifgaps,$catout,$rprint,$v,$help);
my $more = 0;
my $nboot = 100;
my $inters = 10;
my $opt_success = GetOptions(
			 	  'prot=s'		=> \$prot,
			 	  'lincs=s'		=> \$lincs,
			 	  'shuffle=s'   => \$shuffle,
			 	  'inter=s'     => \$inters,
			 	  'nboot=s'     => \$nboot,
			 	  'more=s'      => \$more,
			 	  'range=s'		=> \$build,
			 	  'build'       => \$ifbuild,
			 	  'gap=s'		=> \$gap_file,
			 	  'dogaps'      => \$ifgaps,
			 	  'cat=s'		=> \$catout,
			 	  'text'        => \$rprint,
			 	  'version'     => \$v,
			 	  'help'		=> \$help,);

#Check options, files
die "\n --- $scriptname version $version\n\n" if $v;
die $usage if $help || ! $opt_success;
die $usage unless $prot && $gap_file && $shuffle;
die "\n -p $prot does not exist?\n\n"  if (! -e $prot);
die "\n -l $lincs does not exist?\n\n"  if (! -e $lincs);
die "\n -s $shuffle does not exist?\n\n"  if (! -e $shuffle);
die "\n -r $build does not exist?\n\n"  if (! -e $build);
die "\n -g $gap_file does not exist?\n\n"  if (! -e $gap_file);
die "\n -n $nboot but should be an integer\n\n" if ($nboot !~ /\d+/);
die "\n -i $inters but should be an integer\n\n" if ($inters !~ /\d+/);
die "\n -m $more but should be an integer\n\n" if ($more !~ /\d+/);
($ifbuild)?($ifbuild = "n"):($ifbuild = "y");
($ifgaps)?($ifgaps = "n"):($ifgaps = "y");

#Main stuff now
print STDERR "\n --- $scriptname v$version\n";
print STDERR " --- loading build from $build...\n";
my $build_r = parse_build($build,$ifbuild);

print STDERR " --- loading gaps from $gap_file\n";
my $gaps = load_gap($gap_file,$ifgaps);

print STDERR  " --- Getting acceptable chromosome ranges...\n";
my ($okay_genomic_ranges, $n_ranges, $okay_genomic_ranges_tree)  = exclude_left($gaps, $build_r);

print STDERR " --- loading lncRNAs from $lincs\n";
my ($linc_features,$ln_exons) = read_gff($lincs, $okay_genomic_ranges_tree, 'transcript', 0);

print STDERR " --- loading coding reference from $prot\n";
my  ($ref_seq_coding,$pc_exons) = read_gff($prot, $okay_genomic_ranges_tree, 'mRNA', 1);

print STDERR " --- loading features to shuffle from $shuffle\n";
my $features = feature_to_shuffle($shuffle, $okay_genomic_ranges_tree);

$more = 1 if ($more == 0); #1 rep if set to 0, same thing here
#clean up outputs
my $outl = "$lincs.no_boot";
my $outp = "$prot.no_boot";
my $outlb = "$lincs.boot";
my $outpb = "$prot.boot";
#move to "previous" outputs
`mv $outl $outl.previous` if (-e $outl);
`mv $outp $outp.previous` if (-e $outp);
`mv $outlb $outlb.previous` if (-e $outlb);
`mv $outpb $outpb.previous` if (-e $outpb);
`mv $catout.no-boot.txt $catout.no-boot.txt.previous` if (($catout) && (-e "$catout.no-boot.txt"));
`mv $catout.boot.txt $catout.boot.txt.previous` if (($catout) && (-e "$catout.boot.txt"));
#Now loop
foreach (my $j = 1; $j <= $more; $j++) {
	print STDERR " --- getting overlaps (observed)\n";
	print STDERR "  -- ROUND $j\n" unless ($more == 1);		
	print STDERR "   - for $lincs\n";
	check_for_featured_overlap($features, $linc_features, $ln_exons, 'no_boot.'.$j, 'lncRNA', $outl, $inters);
	print STDERR "   - for $prot\n";
	check_for_featured_overlap($features, $ref_seq_coding, $pc_exons, 'no_boot.'.$j, 'protein', $outp, $inters);
	`cat $outl >> $catout.no-boot.txt` if (($catout) && (-e $outl));
	`cat $outp >> $catout.no-boot.txt` if (($catout) && (-e $outp));

	if ($nboot > 0) {
		print STDERR "  -> getting overlaps (with $nboot bootstraps => \"expected\") for both $lincs and $prot\n";
		for (my $i = 1; $i <= $nboot; $i++) {
			print STDERR "     => bootstrap $i\n";
			my ($linc_features,$ln_exons) = read_gff($lincs, $okay_genomic_ranges_tree, 'transcript', 0);
			my ($ref_seq_coding,$pc_exons) = read_gff($prot, $okay_genomic_ranges_tree, 'mRNA', 1);
			
			my $s_features = shuffle_features_no_overlap($features, $okay_genomic_ranges, $n_ranges);		
			check_for_featured_overlap($s_features, $linc_features, $ln_exons, 'boot.'.$j.".".$i, 'lncRNA', $outlb, $inters);
			check_for_featured_overlap($s_features, $ref_seq_coding, $pc_exons, 'boot.'.$j.".".$i, 'mRNA', $outpb, $inters);
			`cat $outlb >> $catout.boot.txt` if (($catout) && (-e $outlb));
			`cat $outpb >> $catout.boot.txt` if (($catout) && (-e $outpb));
		}
	}
}
my $cmdlines = $lincs;
($cmdlines =~ /\//)?($cmdlines =~ s/(.*)\/.*$/$1\/TEanalysis_Shuffle.R.example.txt/):($cmdlines = "TEanalysis_Shuffle.R.example.txt");
print_Rcmdlines($cmdlines,$scriptname,$version) if ($rprint);
print STDERR " --- $scriptname done\n";
print STDERR "     R command lines in $cmdlines\n" if ($rprint);
print STDERR "\n";
exit;

#-----------------------------------------------------------------------------
#-------------------------------- SUBROUTINES --------------------------------
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# In order to add a gap from the beginning of the genome to the first exome
# regions you need to parse the build so that you know from 0 -> start and  
# 999 -> end of each seqid.  data is $BUILD{seqid} = length;
# my $build_r = parse_build($build);
#-----------------------------------------------------------------------------
sub parse_build{
    my ($file,$ifbuild) = @_;
    my %build;
    my %chrs;
    if ($ifbuild eq "y") {
		open (my $in, '<', $file) or confess "ERROR (sub parse_build): can't open to read $file $!\n";
		GAP_LINE: while(<$in>){
			chomp;
			my @l = split /\t/, $_;
			$chrs{$l[0]} = $l[1]; 
		}
		close $in;
    } else {
    	my $lengths = $file.".lengths.tab";
    	open (my $lenfh, ">", $lengths) or confess "ERROR (sub parse_build): can't open to write $lengths $!\n";
		my $fa = Bio::SeqIO->new(-file => $file, -format => "fasta") or confess "ERROR (sub parse_build): can't open to read $file $!\n";
		my %replen = ();
		while( my $seq = $fa->next_seq() ) {
			my $name = $seq->display_id;
			my $len = $seq->length;
			$chrs{$name} = $len;
			print $lenfh "$name\t$len\n";
		}
		undef $fa;
		close $lenfh;    
    }
    #get the range
    while(my ($feat, $val) = each %chrs){
		my $set = Set::IntSpan::Fast->new();
		$set->add_range(1, $val); #add_range( $from, $to ) <= the whole chromosome for now
		$build{$feat} = $set;
    }
    return \%build;
}

#-----------------------------------------------------------------------------
# loading assembly gaps
# my $gaps = load_gap($gap_file);
#-----------------------------------------------------------------------------
sub load_gap {
    my ($file,$ifbuild) = @_;
    my %gaps;  
    open (my $in, '<', $file) or confess "ERROR (sub load_gap): can't open to read $file $!\n";
	$file = print_gap($file) if ($ifbuild eq "n");
	GAP_LINE: while(<$in>){
		chomp;
		next GAP_LINE if $_ =~ /^\#/;
		next GAP_LINE if ! defined $_;
		my ($bin, $chrom, $chromStart, $chromEnd, $ix, $n, $size, $type, $bridge) = split /\t/, $_;     
		my $set = Set::IntSpan::Fast->new() if ! defined $gaps{$chrom};
		$gaps{$chrom} = $set if ! defined $gaps{$chrom};
		$gaps{$chrom}->add_range($chromStart + 1, $chromEnd + 1);
	}
    close $in;
    return (\%gaps);
}

#-----------------------------------------------------------------------------
# printing assembly gaps if needed
# $file = print_gap($file) if ($ifbuild eq "n");
#-----------------------------------------------------------------------------
sub print_gap {
	my $file = shift;
	my $fasta = Bio::SeqIO->new(-file => $file, -format => "fasta") or confess "ERROR (sub print_gap): failed to create SeqIO object from $file $!\n";
	my $gaps = "$file.gaps.tab";
	open my $gapfh, ">", $gaps or confess "ERROR (sub print_gap): could not open to write $gaps $!\n";
	print $gapfh "#bin\tchrom\tchromStart\tchromEnd\tix\tn\tsize\ttype\tbridge\n";
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
			} else { #end of a gap, print and reintialize (already skipped if length is 0)
				my $Nend = $Nlen + $Nstart - 1;
				print $gapfh ".\t$head\t$Nstart\t$Nend\t.\t.\t$Nlen\t.\t.\n";
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
# Getting acceptable chromosome ranges
# my ($okay_genomic_ranges, $n_ranges, $okay_genomic_ranges_tree)  = exclude_left($gaps, $build_r);
#-----------------------------------------------------------------------------
sub exclude_left {
    my ($remove, $total) = @_;
    my %okay_ranges;
    my %n_ranges;
    my %okay_range_tree;
    while(my($chr, $set) = each %{$total}){
		if(defined $remove->{$chr}){
			my $new = $set->diff($remove->{$chr}); #Return a set containing all the elements that are in this set but not the supplied set
			my @ranges = split /\,/, $new->as_string();
			RANGE: foreach my $range (@ranges){
				my @start_end = split /-/, $range;
				next RANGE if ! defined $start_end[1]; #Somehow that happens
				$okay_range_tree{$chr} = Set::IntervalTree->new if ! defined $okay_range_tree{$chr};
				$okay_range_tree{$chr}->insert($range, $start_end[0], $start_end[1]); #PerlIntervalTree::insert(SV *value, long low, long high), used afterwards to check all features
			}
			$okay_ranges{$chr} = \@ranges; #used for the sfuffle of the TEs
			$n_ranges{$chr} = scalar @ranges; 
		}
    }
    return (\%okay_ranges, \%n_ranges, \%okay_range_tree);
}

#-----------------------------------------------------------------------------
# Load gff3 files
# my ($linc_features,$ln_exons) = read_gff($lincs, $okay_genomic_ranges_tree, 'transcript', 0);
# my ($ref_seq_coding,$pc_exons) = read_gff($prot, $okay_genomic_ranges_tree, 'mRNA', 1);
#-----------------------------------------------------------------------------
sub read_gff {
	my ($gff3_file, $chr_pos_tree, $type, $coding) = @_;
    my %SHIT;
	
	#load annotations through GAL::Annotation
    my $annotation = GAL::Annotation->new($gff3_file);
    my $features = $annotation->features;
    my $genes = $features->search({type => 'gene'});
    print STDERR "     GAL::Annotation has finished loading\n";

	GENE: while (my $gene = $genes->next) {
		if($coding eq 1){
			next GENE unless $gene->is_coding; #function updated Jan 2016 by Barry Moore to return true if any child is mRNA or has CDS exons
		}
		my $gene_id = $gene->feature_id;
		my $seqid   = $gene->seqid;
		my $begin   = $gene->start;
		my $end     = $gene->end;

		if (! defined $chr_pos_tree->{$seqid}){
			print STDERR "     Warning: $seqid is not found in build\n";
			next GENE;
		}
		
		my @range = @{$chr_pos_tree->{$seqid}->fetch($begin, $end)};
		if (! defined $range[0]){
			print STDERR "     Warning: $begin-$end isn't in an acceptable genomic region\n";
			next GENE;
		}
		
		my @tr = $gene->transcripts;		
		my @shuffled = shuffle(@tr); #instead of random number, shuffle, so can just go through if not coding
		my $tr_id = $shuffled[0]->feature_id;
		#if coding, need to check that the one picked at random is not non coding
		my $n_tr = scalar @tr;
		my $got_one = 0;
		if (($coding eq 1) && (! $shuffled[0]->has_CDS)) {
			TRCODING: for (my $j = 1; $j < $n_tr; $j++) { #could do a shift, more elegant, but I'd rather avoid an "until" loop even if in theory there will be at least one coding transcript
				$tr_id = $shuffled[$j]->feature_id;
				$got_one = 1 if ($shuffled[$j]->has_CDS);
				last TRCODING if ($shuffled[$j]->has_CDS);
			}
			print STDERR "     Warning: $gene_id did not have coding transcripts (???)\n" if (($coding eq 1) && ($got_one == 0));
		} 		
		$SHIT{$gene_id} = $tr_id;		
	}

    print STDERR "     loaded one random transcript per gene (type=$type)\n";
    my %MASTER_TREE;
    my $total_exons_in_data;
    my $total_transcripts;
    my $total_exons;
	
	FEAT: foreach my $gene (keys %SHIT){
		my $transcript_id = $SHIT{$gene};
		my ($transcript) = $features->search({feature_id => [$transcript_id]});
		my $seqid  = $transcript->seqid;
		my $transcript_strand = $transcript->strand;
		if($transcript_strand !~ /\+|-/){
			print STDERR "     Warning: $transcript_strand is a bullshit strand ($transcript_id)\n";
			next FEAT;
		}

		#feature coordinates; load them on genomic                                                                                                      
		my $transcript_start = $transcript->start;
		my $transcript_end   = $transcript->end;
		#($transcript_start, $transcript_end) = ($transcript_end, $transcript_start) if $transcript_strand eq '-'; #remove this otherwise ranges in Interval::Tree don't work
		my @range = @{$chr_pos_tree->{$seqid}->fetch($transcript_start, $transcript_end)};
		#print STDERR " -> $transcript_id: ".scalar(@range)." intervals found in @range.\n";
		if (! defined $range[0]){
			print STDERR "     Warning: $transcript_start-$transcript_end isn't in an acceptable genomic region ($transcript_id)\n";
			next FEAT;
		}
		my $range = $range[0];

		$MASTER_TREE{$seqid}{$range} = Set::IntervalTree->new if ! defined $MASTER_TREE{$seqid}{$range};
		my %transcript_info = ('type'       => $type,
							   'name'       => $transcript->feature_id,
							   'start'      => $transcript_start,
							   'strand'     => $transcript_strand,
							   'end'        => $transcript_end,
							   'exon_info'  => undef,
							   'exons'      => undef,
							   'exon_count' => undef);

		my $exon_tree = Set::IntervalTree->new;
		my @exons;
		if ($transcript->type eq 'mRNA') {
			@exons = sort { $a->start <=> $b->start } $transcript->CDSs;
			print STDERR "exon array @exons\n";
			if(! defined $exons[0]){
				print STDERR "     Warning: no CDS for $transcript_id";
				next FEAT;
			}
		} else {
			@exons = sort { $a->start <=> $b->start } $transcript->exons;
		}
		my $exon_count = 0;
		my $exon_total = scalar @exons;
		foreach my $exon(@exons){
			$exon_count++;
			$total_exons_in_data++;
			my $exon_name  = $exon->feature_id;
			my $exon_start = $exon->start;
			my $exon_end   = $exon->end;
			#($exon_start, $exon_end) = ($exon_end, $exon_start) if $exon->strand eq '-'; #remove this otherwise ranges in Interval::Tree don't work
			my %exon_info = ('name' => $exon_name,
							 'start' => $exon_start,
							 'end' => $exon_end,
							 'exon_count' => $exon_count,
							 'total_exons' => $exon_total);

			$exon_tree->insert(\%exon_info, $exon_start, $exon_end);
			$transcript_info{exon_info}{$exon_count} = \%exon_info;
		}

		$transcript_info{exons} = $exon_tree;
		$transcript_info{exon_count} = $exon_count;
		#Insert a range into the interval tree and associate it with a perl scalar; the first value needs to be the low of the range
		$MASTER_TREE{$seqid}{$range}->insert(\%transcript_info, $transcript_start, $transcript_end);
		$total_transcripts++;
	}
    print STDERR "        total transcripts (type=$type): $total_transcripts\n";	
    print STDERR "        total exons in transcripts (type=$type): $total_exons_in_data\n";
    return (\%MASTER_TREE,$total_exons_in_data);
}

#-----------------------------------------------------------------------------
# Load the features to shuffle = not shuffled at that point
# my $features = feature_to_shuffle($shuffle, $okay_genomic_ranges_tree);
#-----------------------------------------------------------------------------
sub feature_to_shuffle {
	my ($file, $okay) = @_;   
	my %features;
    open (my $IN, '<', $file) or confess "ERROR (sub feature_to_shuffle): Can't open to read $file$!\n";
	GFF3_LINE: while(<$IN>){   
		chomp;
		next GFF3_LINE if $_ =~ /^\#/;
		next GFF3_LINE if ! defined $_;
		last GFF3_LINE if $_ =~ /^\>/; # ?  
		my ($seqid, $source, $type, $start, $end, $score, $strand, $phase, $attributes) = split /\t/, $_;
		next GFF3_LINE if ! defined $seqid;
		my @attributes = split /=|;/, $attributes;
		next GFF3_LINE if ! defined $okay->{$seqid};
		my ($range) = @{$okay->{$seqid}->fetch($start, $end)};
		if (! defined $range){
			print STDERR "     Warning: $start-$end isn't in an acceptable genomic region\n";
			next GFF3_LINE;
		}
		$features{$attributes[1]."#".$seqid."#".$start."#".$end} = [$start, $end, $range, $seqid, $strand]; #need to make the key unique
	}
    return \%features;
}

#-----------------------------------------------------------------------------
# Shuffle the features
# my $s_features = shuffle_features_no_overlap($features, $okay_genomic_ranges, $n_ranges);	
#-----------------------------------------------------------------------------
sub shuffle_features_no_overlap {
    my ($features_to_shuffle, $okay_ranges, $n_ranges) = @_;
    my %new_feature_pos;
    my %new_feature_overlap_check;   
    my $n_features = scalar (keys %{$features_to_shuffle});
	FEATURE: while(my ($feature_id, $values) = each %{$features_to_shuffle}){     
		my ($start, $end, $range, $seqid, $strand) = @{$values};
		next FEATURE if ! defined $okay_ranges->{$seqid};
		my $flag = 0;
		my $new_start = 0;
		my $new_end   = 0;
		my $random_range;
		my $range_info;
		while($flag == 0){
			#Get new start and new end
			$random_range = int(rand($n_ranges->{$seqid}));
			$range_info   = $okay_ranges->{$seqid}[$random_range];
			my ($range_min, $range_max) = split /-/, $range_info;
			next FEATURE if (! $range_max); #Seems like sometimes this was undef... But I could not track this error, $okay_ranges when dumped shows min and max for all
			$range_max = $range_max - ($end - $start);
			my $new_start = int($range_min + rand($range_max - $range_min));
			my $new_end   = $new_start + ($end - $start); 
			#Now save
			$new_feature_overlap_check{$seqid}{$range_info} = Set::IntervalTree->new if (! defined $new_feature_overlap_check{$seqid}{$random_range});
			my @overlap = @{$new_feature_overlap_check{$seqid}{$range_info}->fetch($new_start, $new_end)};
			$flag = 1 if (! defined $overlap[0]); #make sure no overlap with the randosized feature
		}
		$new_feature_overlap_check{$seqid}{$range_info}->insert(\$feature_id, $new_start, $new_end);
		$new_feature_pos{$feature_id} = [$new_start, $new_end, $range_info, $seqid];
	}
    return \%new_feature_pos;
}

#-----------------------------------------------------------------------------
# Check overlap with TEs and count
# check_for_featured_overlap($features, $linc_features, $ln_exons, 'no_boot', 'lncRNA', $outl, $inters);
# check_for_featured_overlap($s_features, $linc_features, $ln_exons, $i, 'lncRNA', $outl, $inters);
#-----------------------------------------------------------------------------
sub check_for_featured_overlap { #This is checked looping through TEs
	my ($tes, $feat, $n_exons, $boot_val, $type, $out, $inters) = @_;
	my %type_overlap = ('TSS_polyA'      => 0,
						'TSS_SPL5'       => 0,
						'TSS'            => 0,
						'SPL5'           => 0,
						'SPL3'           => 0,
						'SPL3_exon_SPL5' => 0,
						'Exonized'       => 0,
						'SPL3_polyA'     => 0,
						'polyA'          => 0);
    my %hit_exons;
	my %hit_exons_cat;
    my %n_transcripts;
    my %intron; #to check if there are intron hits only
	my %exonized;
	while(my($feature_id, $value) = each %{$tes}){ #loaded TEs with their coordinates and the OK range they are in
		my ($start, $end, $random_range, $seqid) = @{$value};
		if(defined $feat->{$seqid}{$random_range}){
			#fetch any feature overlapping the TE coordinates => fetch returns an arrayref of perl objects whose ranges overlap the specified range.
			my @overlapping_features =  @{$feat->{$seqid}{$random_range}->fetch($start, $end)}; 
			FEATURE_HIT: foreach my $transcript (@overlapping_features){
				my $exon_results   = $transcript->{'exons'};
				my @exons          = @{$exon_results->fetch($start, $end)}; #Get all exons of that transcript
				$intron{$transcript->{'name'}} = 0 if (! $exons[1]); #if only one exon, can't have intron only hit
				my $strand = $transcript->{strand};
				$n_transcripts{$transcript->{'name'}} = 1; #count unique transcript hits [note: could be intron only]
				my $TEov = 0;

				#if the feature (here, transcript) is contained within the repetitive element.                                                                     
				if(is_contained($start, $end, $transcript->{start}, $transcript->{end})){
					$TEov = $transcript->{end}-$transcript->{start}+1;
					next FEATURE_HIT if ($TEov < $inters); #Won't happen in default but user could up it
					$type_overlap{TSS_polyA}++;
					$intron{$transcript->{'name'}} = 0; #doesn't matter if no introns, it has an exon hit
					foreach my $exon_id (@exons){ 
						$hit_exons{TSS_polyA}{$exon_id->{name}} = 1; #count each exon of the transcript as hit
						$hit_exons_cat{TSS_polyA}{$exon_id->{name}} = 1; #"label" all exons of the transcript with TSS_polyA
					}
					next FEATURE_HIT;
				}
				
				#now if NOT TSS_polyA, check the rest [could be condensed but it's easier to read that way]			
				#positive strand rules           
				if($strand eq '+'){
				#if a feature hangs off the start.                                                                                     
					if($start < $transcript->{start}){
						$TEov = $transcript->{start}-$start+1;
						$intron{$transcript->{'name'}} = 0 unless ($TEov < $inters);
						if ($end < $transcript->{exon_info}{1}{end}) {
							$type_overlap{TSS}++ unless ($TEov < $inters);
							$hit_exons{TSS}{$exons[0]->{name}} = 1 unless ($TEov < $inters);
							next FEATURE_HIT;
						} elsif ($end > $transcript->{exon_info}{1}{end}) {
							$type_overlap{TSS_SPL5}++ unless ($TEov < $inters);
							$hit_exons{TSS_SPL5}{$exons[0]->{name}} = 1 unless ($TEov < $inters);
						}
					}
					#if a feature hangs off the end.                                                                                            
					if($end > $transcript->{end}){
						$TEov = $end-$transcript->{end}+1;
						$intron{$transcript->{'name'}} = 0 unless ($TEov < $inters);
						if ($start > $transcript->{exon_info}{$transcript->{exon_count}}{start}) {
							$type_overlap{polyA}++ unless ($TEov < $inters);
							$hit_exons{polyA}{$exons[-1]->{name}} = 1 unless ($TEov < $inters);
							next FEATURE_HIT;
						} elsif ($start < $transcript->{exon_info}{$transcript->{exon_count}}{start}) {	
							$type_overlap{SPL3_polyA}++ unless ($TEov < $inters);
							$hit_exons{SPL3_polyA}{$exons[-1]->{name}} = 1 unless ($TEov < $inters);
						}				
					}
				}
				#negative strand rules: start of transcript is $transcript->{end}                                                                                               
				else{
					#if a feature hangs off the start.                                                                                          
					if($end > $transcript->{end}){
						$TEov = $end-$transcript->{end}+1;
						$intron{$transcript->{'name'}} = 0 unless ($TEov < $inters);
						if ($start > $transcript->{exon_info}{1}{end}) {
							$type_overlap{TSS}++ unless ($TEov < $inters);
							$hit_exons{TSS}{$exons[-1]->{name}} = 1 unless ($TEov < $inters);
							next FEATURE_HIT;
						} elsif ($start < $transcript->{exon_info}{1}{end}) {
							$type_overlap{TSS_SPL5}++ unless ($TEov < $inters);
							$hit_exons{TSS_SPL5}{$exons[-1]->{name}} = 1 unless ($TEov < $inters);
						}
					}
					#if a feature hangs off the end.                                                                                            
					if($start < $transcript->{start}){
						$TEov = $transcript->{start}-$start+1;
						$intron{$transcript->{'name'}} = 0 unless ($TEov < $inters);
						if ($end < $transcript->{exon_info}{$transcript->{exon_count}}{start}) {
							$type_overlap{polyA}++ unless ($TEov < $inters);
							$hit_exons{polyA}{$exons[0]->{name}} = 1 unless ($TEov < $inters);
							next FEATURE_HIT;
						} elsif ($end > $transcript->{exon_info}{$transcript->{exon_count}}{start}) {
							$type_overlap{SPL3_polyA}++ unless ($TEov < $inters);
							$hit_exons{SPL3_polyA}{$exons[0]->{name}} = 1 unless ($TEov < $inters);
						}				
					}
				}
				
				#Not check if when not TSS or polyA it also gives some other SPL
				EXON: foreach my $exon (@exons){
					foreach my $ovtype (keys %hit_exons) {
						next EXON if (defined $hit_exons{$ovtype}{$exon->{name}}); #this feature overlap with this exon was already taken in account
					}
					my $Trst =  $transcript->{exon_info}{$transcript->{exon_count}}{start};
					my $Tren =  $transcript->{exon_info}{$transcript->{exon_count}}{end};
					#does the exon contain the feature.                                                                                           
					if(is_contained($exon->{start}, $exon->{end}, $start, $end)){
						$TEov = $exon->{end}-$exon->{start}+1;	
						$intron{$transcript->{'name'}} = 0 unless ($TEov < $inters);				
						next EXON if ($TEov < $inters); #with the default should not happen (Repeat Masker won't mask 10nt only) but user may set it up higher
						$type_overlap{Exonized}++ unless ($exonized{$exon->{name}}); 
						$exonized{$exon->{name}} = 1;
						$hit_exons{Exonized}{$exon->{name}} = 1;
						#not ++ if already one TE exonized, otherwise it fucks up exon counts (each exon can only be counted one time for a given feature)
						next EXON;
					}
					#does the feature contain the exon   
					if(is_contained($start, $end, $exon->{start}, $exon->{end})){
						$TEov = $end-$start+1;
						$intron{$transcript->{'name'}} = 0 unless ($TEov < $inters);
						next EXON if ($TEov < $inters); #with the default should not happen (exons > 10 nt usually...) but user may set it up higher
						$type_overlap{SPL3_exon_SPL5}++;
						$hit_exons{SPL3_exon_SPL5}{$exon->{name}} = 1;
						next EXON;
					} else {
						if ($strand eq '+') {
							if ($start > $Trst) {
								$TEov = $start-$Trst+1;
								$type_overlap{SPL3}++ unless ($TEov < $inters);
								$hit_exons{SPL3}{$exon->{name}} = 1 unless ($TEov < $inters);
								$intron{$transcript->{'name'}} = 0 unless ($TEov < $inters);
							} elsif ($end < $Tren) {
								$TEov = $Tren-$end+1;
								$type_overlap{SPL5}++ unless ($TEov < $inters);
								$hit_exons{SPL5}{$exon->{name}} = 1 unless ($TEov < $inters);
								$intron{$transcript->{'name'}} = 0 unless ($TEov < $inters);
							}
						} else {
							if ($start > $Trst) {
								$TEov = $start-$Trst+1;
								$type_overlap{SPL5}++ unless ($TEov < $inters);
								$hit_exons{SPL5}{$exon->{name}} = 1 unless ($TEov < $inters);
								$intron{$transcript->{'name'}} = 0 unless ($TEov < $inters);
							} elsif ($end < $Tren) {
								$TEov = $Tren-$end+1;
								$type_overlap{SPL3}++ unless ($TEov < $inters);
								$hit_exons{SPL3}{$exon->{name}} = 1 unless ($TEov < $inters);
								$intron{$transcript->{'name'}} = 0 unless ($TEov < $inters);
							}
						}				
					}
				}
			}
		}
	}
	
	open my $outfh, ">>", $out or confess "ERROR (sub check_for_featured_overlap): can't open to write $out $!\n";
	PRINT_RESULTS: while(my($key, $value) = each %type_overlap){
		my $n_uniq_exons = scalar keys %{$hit_exons{$key}};
		my $unhit_exons = $n_exons - $n_uniq_exons; #same number for each run, but may change between runs
		print $outfh "$type\t$boot_val\t$key\t$value\t$n_uniq_exons\t$n_exons\t$unhit_exons\n";
	}
	close $outfh;
    
    my $n_introns = 0;
    IFINTRON: foreach my $trname (keys %n_transcripts) {   	
    	next IFINTRON if (defined $intron{$trname});
    	$n_introns++;
    }    
    my $n = keys %n_transcripts;
	print STDERR "        type=$type had $n unique transcripts hit(s) and $n_introns intron only hit(s)\n";
	`rm -Rf $out` if (($n == 0) && ($n_introns == 0)); 
}

#-----------------------------------------------------------------------------
# does the left contain the right
# called by check_for_featured_overlap
#-----------------------------------------------------------------------------
sub is_contained {
    my ($c_start, $c_end, $start, $end) = @_;
    my $results = 0;
    return $results if $start < $c_start;
    return $results if $end   > $c_end;   
    return 1;
}

#-----------------------------------------------------------------------------
# Print R command lines (examples)
# print_Rcmdlines($cmdlines,$scriptname,$version) if ($rprint);
#-----------------------------------------------------------------------------
sub print_Rcmdlines {
    my ($file,$scriptname,$version) = @_;
	
	open my $fh, ">", $cmdlines or confess "ERROR (sub print_Rcmdlines): can't open to write $file $!\n";
	print $fh "#Example of R command lines to process the output files
#From the script $scriptname, v$version
#Prior to running these, concatenate the outputs of -l and -p if -c not used
#Run these command lines for both the boot and no_boot files (X=boot or X=no_boot in the read.csv below)

setwd(\"/Users/my/Documents/ProjectName\")
library(ggplot2)
dat<-read.csv(\"file.cat.X\", sep=\"\\t\", header=FALSE)
#dat
new.dat<-aggregate(dat\$V4, by=list(dat\$V1,dat\$V3),mean)
new.dat\$sd<-aggregate(dat\$V4, by=list(dat\$V1,dat\$V3),sd)\$x
new.dat\$mean_exon<-aggregate(dat\$V5, by=list(dat\$V1,dat\$V3),mean)\$x
new.dat\$total<-aggregate(dat\$V6, by=list(dat\$V1,dat\$V3),mean)\$x
new.dat\$unhit<-aggregate(dat\$V7, by=list(dat\$V1,dat\$V3),mean)\$x
#new.dat
dat2<-dat[grep(dat\$V3, pattern=\"TSS_PolyA\", invert=TRUE),]
ggplot(dat2, aes(x=V3, y=V4))+geom_boxplot()+facet_grid(.~V1)
new.dat";

	close $fh;

    return 1;
}






# #-----------------------------------------------------------------------------
# sub fisher_yates_shuffle {
#     my $array = shift;
#     my $i;
#     for ($i = @$array; --$i; ) {
#         my $j = int rand ($i+1);
#         next if $i == $j;
#         @$array[$i,$j] = @$array[$j,$i];
#     }
# }

# #-----------------------------------------------------------------------------
# sub print_te {
#     my $features = shift;
# #    $features{$attributes[1]} = [$start, $end, $range, $seqid, $strand];
#     while (my($id, $data) = each %{$features}){
# 		my ($start, $end, $range, $seqid) = @{$data};
# 		print join "\t", ($seqid,$start,$end,$id);
# 		print "\n";
#     }
# }


