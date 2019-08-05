#!/usr/bin/perl
use strict;

use Test::More;

BEGIN {
	use lib '../../';
	#TEST 1
	use_ok('GAL::Schema::Result::Feature::mrna');
	use_ok('GAL::Annotation');
}

my $path = $0;
$path =~ s/[^\/]+$//;
$path ||= '.';
chdir($path);

# TEST 2
my $object = GAL::Schema::Result::Feature::mrna->new();
isa_ok($object, 'GAL::Schema::Result::Feature::mrna');

ok(my $schema = GAL::Annotation->new('data/dmel-4-r5.46.genes.gff',
				     'data/dmel-4-chromosome-r5.46.fasta'),
   '$schema = GAL::Annotation->new("dmel.gff", "dmel.fasta"');

ok(my $features = $schema->features, '$features = $schema->features');

#--------------------------------------------------------------------------------
# Plus strand
#--------------------------------------------------------------------------------
ok(my ($mRNA_plus) = $features->search({feature_id => 'FBtr0089157'}),
   '$features->search({feature_id => "FBtr0089157"})');

ok(my $mrna_seq = $mRNA_plus->seq, '$mrna_seq = $mRNA->seq');

ok(my $CDSs = $mRNA_plus->CDSs, '$mRNA_plus->CDSs');
ok(my $first_CDS = $CDSs->first, '$first_CDS = $CDSs->first');
ok($first_CDS->type eq 'CDS', '$first_CDSs->type');
ok(substr($first_CDS->seq, 0, 3) eq 'ATG', 'substr($first_CDS->seq, 0, 3) eq "ATG"');

ok(my ($last_CDS) = sort({$b->start <=> $a->start} $mRNA_plus->CDSs), '$last_CDS');
ok(substr($last_CDS->seq, -3, 3) =~ /^(TAG|TGA|TAA)$/, 'substr($last_CDS->seq, -3, 3) =~ /^(TAG|TGA|TAA)$/');

ok(my $gseq = $mRNA_plus->CDS_seq_genomic, '$mRNA_plus->CDS_seq_genomic');
ok($gseq =~ /^[ATGCN]+$/i, '$gseq =~ /^[ATGCN]+$/i');

ok($mRNA_plus->protein_seq =~ /^[A-Z]+$/, '$mRNA_plus->protein_seq =~ /^[A-Z]+$/');

ok(my $translation_start = $mRNA_plus->translation_start,'$mRNA_plus->translation_start');
ok(my $fp_UTRs = $mRNA_plus->five_prime_UTRs, '$mRNA_plus->five_prime_UTRs');
while (my $fp_UTR = $fp_UTRs->next) {
  ok($fp_UTR->my_end < $translation_start, '$fp_UTR->end < $translation_start');
}

ok(my $translation_end  = $mRNA_plus->translation_end, '$mRNA_plus->translation_end');
ok(my $tp_UTRs = $mRNA_plus->three_prime_UTRs, '$mRNA_plus->three_prime_UTRs');
while (my $tp_UTR = $tp_UTRs->next) {
  ok($tp_UTR->my_start > $translation_end, '$tp_UTR->start < $translation_end');
}

ok(! ($mRNA_plus->map2CDS(93055))[0],        '! ($mRNA_plus->map2CDS(93055))[0]');
ok(($mRNA_plus->map2CDS(93056))[0] == 1,     '($mRNA_plus->map2CDS(93056))[0] == 1');
ok(($mRNA_plus->map2CDS(93057))[0] == 2,     '($mRNA_plus->map2CDS(93057))[0] == 2');
ok(($mRNA_plus->map2CDS(93058))[0] == 3,     '($mRNA_plus->map2CDS(93058))[0] == 3');
ok(($mRNA_plus->map2CDS(93059))[0] == 4,     '($mRNA_plus->map2CDS(93059))[0] == 4');
ok(($mRNA_plus->map2CDS(93060))[0] == 5,     '($mRNA_plus->map2CDS(93060))[0] == 5');
ok(! ($mRNA_plus->map2CDS(103263))[0],       '! ($mRNA_plus->map2CDS(93055))[0]');
ok(($mRNA_plus->map2CDS(103264))[0] == 310,  '($mRNA_plus->map2CDS(93056))[0] == 1');
ok(($mRNA_plus->map2protein(93056))[0] == 1, '($mRNA_plus->map2protein(93056))[0] ==1');
ok(($mRNA_plus->map2protein(93057))[0] == 1, '($mRNA_plus->map2protein(93057))[0] ==1');
ok(($mRNA_plus->map2protein(93058))[0] == 1, '($mRNA_plus->map2protein(93058))[0] ==1');
ok(($mRNA_plus->map2protein(93059))[0] == 2, '($mRNA_plus->map2protein(93059))[0] ==2');
ok(($mRNA_plus->map2protein(93060))[0] == 2, '($mRNA_plus->map2protein(93060))[0] ==2');
ok(($mRNA_plus->map2protein(93061))[0] == 2, '($mRNA_plus->map2protein(93061))[0] ==2');
ok($mRNA_plus->CDS_start == 93056,           '$mRNA_plus->CDS_start');
ok($mRNA_plus->CDS_end == 129118,            '$mRNA_plus->CDS_end');
ok($mRNA_plus->CDS_length == 2256,           '$mRNA_plus->CDS_length');
ok($mRNA_plus->protein_length == 751,        '$mRNA_plus->protein_length');

ok(! $mRNA_plus->phase_at_location(93055),   '$mRNA_plus->phase_at_location');
ok($mRNA_plus->phase_at_location(93056) == 0, '$mRNA_plus->phase_at_location');
ok($mRNA_plus->phase_at_location(93057) == 2, '$mRNA_plus->phase_at_location');
ok($mRNA_plus->phase_at_location(93058) == 1, '$mRNA_plus->phase_at_location');
ok($mRNA_plus->phase_at_location(93059) == 0, '$mRNA_plus->phase_at_location');

ok(! $mRNA_plus->frame_at_location(93055),    '$mRNA_plus->frame_at_location');
ok($mRNA_plus->frame_at_location(93056) == 0, '$mRNA_plus->frame_at_location');
ok($mRNA_plus->frame_at_location(93057) == 1, '$mRNA_plus->frame_at_location');
ok($mRNA_plus->frame_at_location(93058) == 2, '$mRNA_plus->frame_at_location');
ok($mRNA_plus->frame_at_location(93059) == 0, '$mRNA_plus->frame_at_location');

ok(! $mRNA_plus->codon_at_location(93055),        '$mRNA_plus->codon_at_location');
ok($mRNA_plus->codon_at_location(93056) eq 'ATG', '$mRNA_plus->codon_at_location');
ok($mRNA_plus->codon_at_location(93057) eq 'ATG', '$mRNA_plus->codon_at_location');
ok($mRNA_plus->codon_at_location(93058) eq 'ATG', '$mRNA_plus->codon_at_location');
ok($mRNA_plus->codon_at_location(93059) eq 'CCT', '$mRNA_plus->codon_at_location');

#--------------------------------------------------------------------------------
# Minus strand
#--------------------------------------------------------------------------------

ok(my ($mRNA_minus) = $features->search({feature_id => 'FBtr0304048'}),
   '$features->search({feature_id => "FBtr0304048"})');

ok(my $mrna_seq = $mRNA_minus->seq, '$mrna_seq = $mRNA->seq');

ok(my $CDSs = $mRNA_minus->CDSs, '$mRNA_minus->CDSs');
ok(my $first_CDS = $CDSs->first, '$first_CDS = $CDSs->first');
ok($first_CDS->type eq 'CDS', '$first_CDSs->type');
ok(substr($first_CDS->seq, 0, 3) eq 'ATG', 'substr($first_CDS->seq, 0, 3) eq "ATG"');

ok(my ($last_CDS) = sort({$a->start <=> $b->start} $mRNA_minus->CDSs), '$first_CDS = $CDSs->last');
ok(substr($last_CDS->seq, -3, 3) =~ /^(TAG|TGA|TAA)$/, 'substr($last_CDS->seq, -3, 3) =~ /^(TAG|TGA|TAA)$/');

ok(my $gseq = $mRNA_minus->CDS_seq_genomic, '$mRNA_minus->CDS_seq_genomic');
ok($gseq =~ /^[ATGCN]+$/i, '$gseq =~ /^[ATGCN]+$/i');
ok($mRNA_minus->protein_seq =~ /^[A-Z]+$/, '$mRNA_minus->protein_seq =~ /^[A-Z]+$/');

ok(my $translation_start = $mRNA_minus->translation_start,'$mRNA_minus->translation_start');
ok(my $fp_UTRs = $mRNA_minus->five_prime_UTRs, '$mRNA_minus->five_prime_UTRs');
while (my $fp_UTR = $fp_UTRs->next) {
  ok($fp_UTR->my_end > $translation_start, '$fp_UTR->my_end > $translation_start');
}

ok(my $translation_end  = $mRNA_minus->translation_end, '$mRNA_minus->translation_end');
ok(my $tp_UTRs = $mRNA_minus->three_prime_UTRs, '$mRNA_minus->three_prime_UTRs');
while (my $tp_UTR = $tp_UTRs->next) {
  ok($tp_UTR->my_start < $translation_end, '$tp_UTR->my_start < $translation_end');
}

#ok(! ($mRNA_minus->map2CDS(93055))[0],        '! ($mRNA_minus->map2CDS(93055))[0]');
#ok(($mRNA_minus->map2CDS(93056))[0] == 1,     '($mRNA_minus->map2CDS(93056))[0] == 1');
#ok(($mRNA_minus->map2CDS(93057))[0] == 2,     '($mRNA_minus->map2CDS(93057))[0] == 2');
#ok(($mRNA_minus->map2CDS(93058))[0] == 3,     '($mRNA_minus->map2CDS(93058))[0] == 3');
#ok(($mRNA_minus->map2CDS(93059))[0] == 4,     '($mRNA_minus->map2CDS(93059))[0] == 4');
#ok(($mRNA_minus->map2CDS(93060))[0] == 5,     '($mRNA_minus->map2CDS(93060))[0] == 5');
#ok(! ($mRNA_minus->map2CDS(103263))[0],       '! ($mRNA_minus->map2CDS(93055))[0]');
#ok(($mRNA_minus->map2CDS(103264))[0],         '($mRNA_minus->map2CDS(93056))[0] == 1');
#ok(($mRNA_minus->map2protein(93056))[0] == 1, '($mRNA_minus->map2protein(93056))[0] ==1');
#ok(($mRNA_minus->map2protein(93057))[0] == 1, '($mRNA_minus->map2protein(93057))[0] ==1');
#ok(($mRNA_minus->map2protein(93058))[0] == 1, '($mRNA_minus->map2protein(93058))[0] ==1');
#ok(($mRNA_minus->map2protein(93059))[0] == 2, '($mRNA_minus->map2protein(93059))[0] ==2');
#ok(($mRNA_minus->map2protein(93060))[0] == 2, '($mRNA_minus->map2protein(93060))[0] ==2');
#ok(($mRNA_minus->map2protein(93061))[0] == 2, '($mRNA_minus->map2protein(93061))[0] ==2');
#ok($mRNA_minus->CDS_start == 93056,           '$mRNA_minus->CDS_start');
#ok($mRNA_minus->CDS_end == 129118,            '$mRNA_minus->CDS_end');
#ok($mRNA_minus->CDS_length == 2256,           '$mRNA_minus->CDS_length');
#ok($mRNA_minus->protein_length == 751,        '$mRNA_minus->protein_length');

#ok($mRNA_minus->phase_at_location(93056), '$mRNA_minus->phase_at_location');
#ok($mRNA_minus->phase_at_location(93057), '$mRNA_minus->phase_at_location');
#ok($mRNA_minus->phase_at_location(93058), '$mRNA_minus->phase_at_location');
#ok($mRNA_minus->phase_at_location(93059), '$mRNA_minus->phase_at_location');
#
#ok($mRNA_minus->frame_at_location(93056), '$mRNA_minus->frame_at_location');
#ok($mRNA_minus->frame_at_location(93057), '$mRNA_minus->frame_at_location');
#ok($mRNA_minus->frame_at_location(93058), '$mRNA_minus->frame_at_location');
#ok($mRNA_minus->frame_at_location(93059), '$mRNA_minus->frame_at_location');
#
#ok($mRNA_minus->codon_at_location(93056), '$mRNA_minus->codon_at_location');
#ok($mRNA_minus->codon_at_location(93057), '$mRNA_minus->codon_at_location');
#ok($mRNA_minus->codon_at_location(93058), '$mRNA_minus->codon_at_location');
#ok($mRNA_minus->codon_at_location(93059), '$mRNA_minus->codon_at_location');

#--------------------------------------------------------------------------------
# Infer UTRs
#--------------------------------------------------------------------------------

system('rm data/dmel-4-r5.46.genes.sqlite') if
  -e 'data/dmel-4-r5.46.genes.sqlite';

ok(my $schema_genes = GAL::Annotation->new('data/dmel-4-r5.46.genes.gff',
					    'data/dmel-4-chromosome-r5.46.fasta'),
   '$schema = GAL::Annotation->new("dmel.gff", "dmel.fasta"');

ok(my $features_genes = $schema_genes->features,
   '$features = $schema->features');

ok(my $mrnas1 = $features_genes->search({type => 'mRNA'}),
   '$features_genes->search({type => \'mRNA\'});');

my %annotated_UTRs;
while (my $mrna = $mrnas1->next) {
  $mrna->infer_five_prime_UTR();
  $mrna->infer_three_prime_UTR();

  for my $fpUTR ($mrna->five_prime_UTRs) {
    my $key = join ':', $fpUTR->get_values(qw(seqid type start end));
    $annotated_UTRs{$key}++;
  }
  for my $tpUTR ($mrna->three_prime_UTRs) {
    my $key = join ':', $tpUTR->get_values(qw(seqid type start end));
    $annotated_UTRs{$key}++;
  }
}
my $UTRs_annotated = join "\n", sort keys %annotated_UTRs;

system('rm data/dmel-4-r5.46.no_UTR.sqlite') if
  -e 'data/dmel-4-r5.46.no_UTR.sqlite';

ok(my $schema_no_UTR = GAL::Annotation->new('data/dmel-4-r5.46.no_UTR.gff',
					    'data/dmel-4-chromosome-r5.46.fasta'),
   '$schema = GAL::Annotation->new("dmel.gff", "dmel.fasta"');

ok(my $features_no_UTR = $schema_no_UTR->features,
   '$features = $schema->features');

ok(my $mrnas2 = $features_no_UTR->search({type => 'mRNA'}),
   '$features_no_UTR->search({type => \'mRNA\'});');

my %infered_UTRs;
while (my $mrna = $mrnas2->next) {
  $mrna->infer_five_prime_UTR();
  $mrna->infer_three_prime_UTR();

  for my $fpUTR ($mrna->five_prime_UTRs) {
    my $key = join ':', $fpUTR->get_values(qw(seqid type start end));
    $infered_UTRs{$key}++;
  }
  for my $tpUTR ($mrna->three_prime_UTRs) {
    my $key = join ':', $tpUTR->get_values(qw(seqid type start end));
    $infered_UTRs{$key}++;
  }
}
my $UTRs_infered = join "\n", sort keys %infered_UTRs;

ok($UTRs_annotated eq $UTRs_infered,
   'Inferred UTRs equal to annoated UTRs');
#compare to above

done_testing();

################################################################################
################################# Ways to Test #################################
################################################################################

__END__



=head3
# Various other ways to say "ok"
ok($this eq $that, $test_name);

is  ($this, $that,    $test_name);
isnt($this, $that,    $test_name);

# Rather than print STDERR "# here's what went wrong\n"
diag("here's what went wrong");

like  ($this, qr/that/, $test_name);
unlike($this, qr/that/, $test_name);

cmp_ok($this, '==', $that, $test_name);

is_deeply($complex_structure1, $complex_structure2, $test_name);

can_ok($module, @methods);
isa_ok($object, $class);

pass($test_name);
fail($test_name);

BAIL_OUT($why);
=cut
