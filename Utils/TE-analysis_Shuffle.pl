#!/usr/bin/perl

#Author: Zev Kronenberg
#Date: 2012
#Writen to generate data shown in Figure 5 in Kapusta et al. 2013 PLoS Genetics (http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1003470)

use strict;
use warnings;
use Getopt::Long;
use Set::IntervalTree;
use GAL::Annotation;
use POSIX;
use Data::Dumper;
use Set::IntSpan::Fast;
#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------

#flush that fucking buffer
$| = 1;


my $usage = "

Synopsis:

BOOTSTRAP_GENOMIC_COORDINATES -p perminant_gff -l linc.gff -b features_to_shuffle -n 100 -r built.txt -g hg19_repeat_mask

Description:

No really, how the hell do you use this thing!

";


my ($help);
my $perminant;
my $lincs;
my $gap_file;
my $boot;
my $nboot = 100;
my $build;
my $opt_success = GetOptions('help'          => \$help,
			     'perminant=s'   => \$perminant,
			     'boot=s'        => \$boot,
			     'range=s'       => \$build,
			     'lincs=s'       => \$lincs,
			     'gap=s'         => \$gap_file);

die $usage if $help || ! $opt_success;


die $usage unless $perminant && $gap_file && $boot;

my $build_r         = parse_build($build);
print STDERR "loaded build...\n";
my $gaps                 = load_gap($gap_file);
print STDERR "loaded gaps...\n";
my ($okay_genomic_ranges, $n_gaps, $okay_genomic_ranges_tree)  = exclude_left($gaps, $build_r);
print STDERR  "found acceptable chromosome ranges...\n";
my $link_features   = read_gff($lincs, $okay_genomic_ranges_tree, 'transcript', 0);
print STDERR "loaded lincRNAs\n";
my $ref_seq_coding  = read_gff($perminant, $okay_genomic_ranges_tree, 'mRNA', 1);
print STDERR "loaded ref_seq...\n";
my $s_features      = feature_to_shuffle($boot, $okay_genomic_ranges_tree);
print STDERR "loaded features to shuffle...\n";

check_for_featured_overlap($s_features,  $ref_seq_coding, 'no_boot', 'protein');
check_for_featured_overlap($s_features,  $link_features,  'no_boot', 'lncRNA');



#-----------------------------------------------------------------------------
#-------------------------------- SUBROUTINES --------------------------------
#-----------------------------------------------------------------------------
sub load_gap {
    my $file = shift;
    my %gaps;
    
    open (my $IN, '<', $file) or die "Can't open $file for reading\n$!\n";
  GAP_LINE: while(<$IN>){
      chomp;
      my ($bin, $chrom, $chromStart, $chromEnd, $ix, $n, $size, $type, $bridge) = split /\t/, $_;
      
      next GAP_LINE if $_ =~ /^\#/;
      next GAP_LINE if ! defined $_;
      my $set = Set::IntSpan::Fast->new() if ! defined $gaps{$chrom};
      $gaps{$chrom} = $set if ! defined $gaps{$chrom};
      $gaps{$chrom}->add_range($chromStart + 1, $chromEnd + 1);
  }
    close $IN;
    return (\%gaps);
}

#-----------------------------------------------------------------------------
# In order to add a gap from the beginning of the genome to the first exome
# regions you need to parse the build so that you know from 0 -> start and  
# 999 -> end of each seqid.  data is $BUILD{seqid} = length;

sub parse_build{
    my $build_dat = shift;
    my %chrs;
    my %build;
    if (defined $build_dat){
        open (my $IN, '<', $build_dat) or die "Can't open $build_dat for reading\n$!\n";
      GAP_LINE: while(<$IN>){
          chomp;
          my @l = split /\t/, $_;
          $chrs{$l[0]} = $l[2];
      }
    }else{
	%chrs = (chr1  =>  249250621,
		 chr2  =>  243199373,
		 chr3  =>  198022430,
		 chr4  =>  191154276,
		 chr5  =>  180915260,
		 chr6  =>  171115067,
		 chr7  =>  159138663,
		 chr8  =>  146364022,
		 chr9  =>  141213431,
		 chr10 =>  135534747,
		 chr11 =>  135006516,
		 chr12 =>  133851895,
		 chr13 =>  115169878,
		 chr14 =>  107349540,
		 chr15 =>  102531392,
		 chr16 =>  90354753,
		 chr17 =>  81195210,
		 chr18 =>  78077248,
		 chr19 =>  59128983,
		 chr20 =>  63025520,
		 chr21 =>  48129895,
		 chr22 =>  51304566,
		 chrX  =>  155270560,
		 chrY  =>  59373566,
		 chrM  =>  16571,
	    );
    }
    while(my ($feat, $val) = each %chrs){
	my $set = Set::IntSpan::Fast->new();
	$set->add_range(1, $val);
	$build{$feat} = $set;
    }
    return \%build;
}
#-----------------------------------------------------------------------------
sub exclude_left {
    my ($remove, $total) = @_;
    my %new_data;
    my %n_ranges;
    my %okay_range_tree;

    while(my($chr, $set) = each %{$total}){
	if(defined $remove->{$chr}){
	    my $new = $set->diff($remove->{$chr});
	    my @ranges = split /\,/, $new->as_string();
	    RANGE: foreach my $range (@ranges){
		my @start_end = split /-/, $range;
		next RANGE if ! defined $start_end[1];
		$okay_range_tree{$chr} = Set::IntervalTree->new if ! defined $okay_range_tree{$chr};
		$okay_range_tree{$chr}->insert($range, $start_end[0], $start_end[1]);
	    }
	    $new_data{$chr} = \@ranges;
	    $n_ranges{$chr} = scalar @ranges; 
	}
    }
    return (\%new_data, \%n_ranges, \%okay_range_tree);
}
#-----------------------------------------------------------------------------
sub feature_to_shuffle {
    
    my ($file, $okay) = @_;
    
    my %features;
    open (my $IN, '<', $file) or die "Can't open $file for reading\n$!\n";
  GFF3_LINE: while(<$IN>){
    
      chomp;
      next GFF3_LINE if $_ =~ /^\#/;
      next GFF3_LINE if ! defined $_;
      last GFF3_LINE if $_ =~ /^\>/;
      
      my ($seqid, $soruce, $type, $start, $end, $score, $strand, $phase, $attributes) = split /\t/, $_;
      next GFF3_LINE if ! defined $seqid;
      my @attributes = split /=|;/, $attributes;
      next GFF3_LINE if! defined $okay->{$seqid};
      my ($range) = @{$okay->{$seqid}->fetch($start, $end)};
      if (! defined $range){
	  print STDERR "Warning: not valid range:: $_\n";
	  next GFF3_LINE;
      }
      $features{$attributes[1]} = [$start, $end, $range, $seqid, $strand];
  }
    return \%features;
}

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
	  $random_range = int(rand($n_ranges->{$seqid}));
	  $range_info   = $okay_ranges->{$seqid}[$random_range];
	  ($new_start, $new_end) = new_start_stop($seqid, $start, $end, $okay_ranges, $random_range);
	  $new_feature_overlap_check{$seqid}{$range_info} = Set::IntervalTree->new if ! defined $new_feature_overlap_check{$seqid}{$random_range};
	  my @overlap = @{$new_feature_overlap_check{$seqid}{$range_info}->fetch($new_start, $new_end)};
	  $flag = 1 if ! defined $overlap[0];
      }

      $new_feature_overlap_check{$seqid}{$range_info}->insert(\$feature_id, $new_start, $new_end);
      $new_feature_pos{$feature_id} = [$new_start, $new_end, $range_info, $seqid];
  }
    return \%new_feature_pos;
}
#-----------------------------------------------------------------------------
sub new_start_stop {
    
    my ($seqid, $start, $end, $okay_ranges, $random_range) = @_;
    my $length = $end - $start;
    
    my $range = $okay_ranges->{$seqid}[$random_range];
    my ($range_min, $range_max) = split /-/, $range;
    $range_max = $range_max - $length;

    my $new_start = int($range_min + rand($range_max - $range_min));
    my $new_end   = $new_start + $length; 

    return ($new_start, $new_end);
}

#-----------------------------------------------------------------------------
sub check_for_featured_overlap {

    my ($shuffled, $features,  $boot_val, $type) = @_;
    my %type_overlap;
    my %hit_exons;
    my %n_transcripts;
    my %intron;

  REP: while(my($feature_id, $value) = each %{$shuffled}){
      my ($start, $end, $random_range, $seqid) = @{$value};
      $start = $start + 10;
      $end   = $end   - 10;
      if(defined $features->{$seqid}{$random_range}){
	 
          my @overlapping_features =  @{$features->{$seqid}{$random_range}->fetch($start, $end)};
        FEATURE_HIT: foreach my $transcript (@overlapping_features){
	    my $intron_flag = 0;
	    
	    my $exon_results   = $transcript->{'exons'};
            my @exons          = @{$exon_results->fetch($start, $end)};
            my $n_hits         = scalar @exons;
            my $strand = $transcript->{strand};
	    $n_transcripts{$transcript->{'name'}} = 1;
            #if the feature is contained within the repetative element.                                                                     
            if(is_contained($start, $end, $transcript->{start}, $transcript->{end})){
		$type_overlap{TSS_PolyA}++;
                foreach my $exon_id (@exons){
                    $hit_exons{$exon_id->{name}} = 1;
		}
                next FEATURE_HIT;
            }

            #positive strand rules                                                                                                          
            if($strand eq '+'){
                #if a feature hangs off the start.                                                                                          
                if($start < $transcript->{start}){
                    $type_overlap{TSS}++      if $end < $transcript->{exon_info}{1}{end};
                    $type_overlap{TSS_SPL}++  if $end > $transcript->{exon_info}{1}{end};
		    $intron_flag = 1;
		}
                #if a feature hangs off the end.                                                                                            
                if($end > $transcript->{end}){
                    $type_overlap{PolyA}++      if $start > $transcript->{exon_info}{$transcript->{exon_count}}{start};
                    $type_overlap{PolyA_SPL}++  if $start < $transcript->{exon_info}{$transcript->{exon_count}}{start};
		    $intron_flag = 1;
		}
            }
            #negative strand rules                                                                                                          
            else{
                #if a feature hangs off the start.                                                                                          
                if($end > $transcript->{start}){
                    $type_overlap{TSS}++      if $start > $transcript->{exon_info}{1}{end};
                    $type_overlap{TSS_SPL}++  if $start < $transcript->{exon_info}{1}{end};
		    $intron_flag = 1; 
		}
                #if a feature hangs off the end.                                                                                            
                if($start < $transcript->{end}){
                    $type_overlap{PolyA}++      if $end < $transcript->{exon_info}{$transcript->{exon_count}}{start};
                    $type_overlap{PolyA_SPL}++  if $end > $transcript->{exon_info}{$transcript->{exon_count}}{start};
		    $intron_flag = 1;
		}
            }

	   $intron{$transcript->{'name'}} = 1 if ! defined $exons[0] && $intron_flag == 0;
          EXON: foreach my $exon (@exons){
              #count the exon                                                                                                               
	      $hit_exons{$exon->{name}} = 1;
              #does the exon contain the feature.                                                                                           
              if(is_contained($exon->{start}, $exon->{end}, $start, $end)){
                  $type_overlap{Exonized}++;
		  next FEATURE_HIT;
              }
              #does the feature contain the exon.              
	      if(is_contained($start, $end, $exon->{start}, $exon->{end})){
		  $type_overlap{Both_SPL}++;
		  next EXON;
              }
              else{
                  $type_overlap{SPL}++;
              }
          }
	}
      }
  }

    my $n_uniq_exons = scalar keys %hit_exons;
  PRINT_RESULTS: while(my($key, $value) = each %type_overlap){
      print "$type\t$boot_val\t$key\t$value\t$n_uniq_exons\n";
  }
    
    my $n_introns = keys %intron;
    my $n = keys %n_transcripts;
	print STDERR "$type had $n uniq transcripts hit\n";
	print STDERR "$type had $n_introns intron hits\n";
}




#-----------------------------------------------------------------------------
# does the left contain the right;
sub is_contained {

    my ($c_start, $c_end, $start, $end) = @_;
    my $results = 0;

    return $results if $start < $c_start;
    return $results if $end   > $c_end;
    
    return 1;
}
#-----------------------------------------------------------------------------
sub fisher_yates_shuffle {
    my $array = shift;
    my $i;
    for ($i = @$array; --$i; ) {
        my $j = int rand ($i+1);
        next if $i == $j;
        @$array[$i,$j] = @$array[$j,$i];
    }
}

#-----------------------------------------------------------------------------
sub print_te {

    my $features = shift;

#    $features{$attributes[1]} = [$start, $end, $range, $seqid, $strand];

    while (my($id, $data) = each %{$features}){
	my ($start, $end, $range, $seqid) = @{$data};
	print join "\t", ($seqid,$start,$end,$id);
	print "\n";
    }
}

#-----------------------------------------------------------------------------

sub read_gff {

    my ($gff3_file, $chr_pos_tree, $type, $coding) = @_;
    my %SHIT;

    my $annotation = GAL::Annotation->new($gff3_file);
    my $features = $annotation->features;
    my $genes = $features->search({type => 'gene'});
    print STDERR "GAL annotation has finished working\n";

  GENE: while (my $gene = $genes->next) {
      if($coding eq 1){
          next GENE unless $gene->is_coding;
      }

      my $gene_id = $gene->feature_id;
      my $seqid   = $gene->seqid;
      my $begin   = $gene->start;
      my $end     = $gene->end;

      if (! defined $chr_pos_tree->{$seqid}){
          print STDERR "Warning: $seqid is not found in build\n";
          next GENE;
      }
      my @range = @{$chr_pos_tree->{$seqid}->fetch($begin, $end)};
      if (! defined $range[0]){
          print STDERR "Warning: $begin-$end isn't in an acceptable genomic region\n";
          next GENE;
      }

      my @mrnas = $gene->transcripts;
      my $n_mrnas = scalar @mrnas;
      my $random_index = int(rand($n_mrnas));

      my $mrna_id = $mrnas[$random_index]->feature_id;
      $SHIT{$gene_id} = $mrna_id;

  }

    print STDERR "loaded random transcripts\n";
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
          print STDERR "Warning: $transcript_strand is a bullshit strand\n";
          next FEAT;
      }

      #feature coordinates                                                                                                            
      my $transcript_start = $transcript->start;
      my $transcript_end   = $transcript->end;

      ($transcript_start, $transcript_end) = ($transcript_end, $transcript_start)
          if $transcript_strand eq '-';

      my @range = @{$chr_pos_tree->{$seqid}->fetch($transcript_start, $transcript_end)};
      if (! defined $range[0]){
          print STDERR "Warning: $transcript_start-$transcript_end isn't in an acceptable genomic region\n";
          next FEAT;
      }
      my $range = $range[0];

      $MASTER_TREE{$seqid}{$range} = Set::IntervalTree->new if ! defined $MASTER_TREE{$seqid}{$range};
      my %transcript_info = ('type'       => 'mRNA',
                             'name'       => $transcript->feature_id,
                             'start'      => $transcript_start,
                             'strand'     => $transcript_strand,
                             'end'        => $transcript_end,
                             'exon_info'  => undef,
                             'exons'      => undef,
                             'exon_count' => undef
          );

      my $exon_tree = Set::IntervalTree->new;
      my @exons;
      if ($transcript->type eq 'mRNA') {
          @exons = sort { $a->start <=> $b->start } $transcript->CDSs;
          if(! defined $exons[0]){
              print STDERR "Warning: no CDS for $transcript_id";
              next FEAT;
          }
      }
      else {
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
          ($exon_start, $exon_end) = ($exon_end, $exon_start)
              if $exon->strand eq '-';
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
      $MASTER_TREE{$seqid}{$range}->insert(\%transcript_info, $transcript_start, $transcript_end);
      $total_transcripts++;
  }

    print STDERR "total exons in $type:  $total_exons_in_data\n";
    print STDERR "total transcripts of $type: $total_transcripts\n";
    return \%MASTER_TREE;

}

#----------------------------------------------------------------------------- 
