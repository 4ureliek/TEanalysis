#!/usr/bin/perl
use strict;

use Test::More;
use List::MoreUtils qw(uniq);

BEGIN {
	use lib '../../';
	use_ok('GAL::Schema::Result::Feature');
	use_ok('GAL::Annotation');
}

my $path = $0;
$path =~ s/[^\/]+$//;
$path ||= '.';
chdir($path);

my $object = GAL::Schema::Result::Feature->new();
isa_ok($object, 'GAL::Schema::Result::Feature');

system('rm data/dmel-4-r5.46.genes.sqlite') if
  -e 'data/dmel-4-r5.46.genes.sqlite';

ok(my $schema = GAL::Annotation->new('data/dmel-4-r5.46.genes.gff',
				     'data/dmel-4-chromosome-r5.46.fasta'),
   '$schema = GAL::Annotation->new("dmel.gff", "dmel.fasta"');

ok(my $features = $schema->features, '$features = $schema->features');

#--------------------------------------------------------------------------------
# Test $feature->seqid
# Test $feature->source
# Test $feature->type
# Test $feature->start
# Test $feature->end
# Test $feature->score
# Test $feature->strand
# Test $feature->phase
# Test $feature->attributes_hash
# Test $feature->attribute_value
# Test $feature->seq
# Test $feature->feature_seq
# Test $feature->gc_content
# Test $feature->genomic_seq
# Test $feature->length
# Test $feature->my_start
# Test $feature->my_end
# Test $feature->genomic_length
# Test $feature->get_feature_bins
# Test $feature->to_gff3
# Test $feature->locus
# Test $feature->get_values
#--------------------------------------------------------------------------------

my $count;
while (my $feature = $features->next) {
  ok(my $seqid  = $feature->seqid,  '$feature->seqid');
  ok(my $source = $feature->source, '$feature->source');
  ok(my $type   = $feature->type,   '$feature->type');
  ok(my $start  = $feature->start,  '$feature->start');
  ok(my $end    = $feature->end,    '$feature->end');
  ok(my $score  = $feature->score,  '$feature->score');
  ok(my $strand = $feature->strand, '$feature->strand');
  my $phase  = $feature->phase;
  ok(defined $phase, 'defined $phase');
  ok(my $atts   = $feature->attributes_hash,
     '$feature->attributes_hash');
  ok(my ($id1) = @{$atts->{ID}},
     '@{$atts->{ID}}');
  ok(my ($id2) = $feature->attribute_value('ID'),
     '$feature->attribute_value("Parent")');
  ok($id1 eq $id2, '$id1 eq $id2,');
  ok(my $seq1 = $feature->seq,
     '$feature->seq');
  ok(my $seq2 = $feature->feature_seq,
     '$feature->feature_seq');
  ok($seq1 eq $seq2, '$seq1 eq $seq2');
  ok(my $gc = $feature->gc_content,
     '$feature->gc_content');
  ok($gc >= 0 && $gc <= 1, '$gc >= 0 && $gc <= 1');
  ok(my $genomic_seq = $feature->genomic_seq,
     '$feature->genomic_seq');
  if ($strand eq '-') {
    ok($seq1 ne $genomic_seq, '$seq1 eq reverse $genomic_seq');
  }
  else {
    ok($seq1 eq $genomic_seq, '$seq1 eq reverse $genomic_seq');
  }
  ok(my $length = $feature->length,
     '$feature->length');
  ok($length =~ /^\d+$/, '$length =~ /^\d+$/');
  ok(my $my_start = $feature->my_start,
     '$feature->my_start');
  if ($strand eq '-') {
    ok($my_start eq $end, '$my_start eq $end');
  }
  else {
    ok($my_start eq $start, '$my_start eq $start');
  }
  ok(my $my_end = $feature->my_end,
     '$feature->my_end');
  if ($strand eq '-') {
    ok($my_end eq $start, '$my_end eq $start');
  }
  else {
    ok($my_end eq $end, '$my_end eq $end');
  }
  ok(my $genomic_length = $feature->genomic_length,
     '$feature->genomic_length');
  ok($genomic_length =~ /^\d+$/, '$genomic_length =~ /^\d+$/');
  ok(my $feature_bins = $feature->get_feature_bins,
     '$feature->get_feature_bins');
  ok(ref $feature_bins eq 'ARRAY', 'ref $feature_bins eq "ARRAY"');
  ok(my $to_gff3 = $feature->to_gff3,
     '$feature->to_gff3');
  ok(my $locus = $feature->locus,
     '$feature->locus');
  ok(my @values = $feature->get_values(qw(seqid source type start end score strand phase)),
     '$feature->get_values');
  my $j1 = join ',', @values;
  my $j2 = join ',', ($seqid, $source, $type, $start, $end, $score, $strand, $phase);
  ok($j1 eq $j2, '$j1 eq $j2');

  ok(my @values = $feature->get_values(qw(phase strand score end start type source seqid)),
     '$feature->get_values');
  my $j3 = join ',', @values;
  my $j4 = join ',', ($phase, $strand, $score, $end, $start, $type, $source, $seqid);
  ok($j3 eq $j4, '$j3 eq $j4');
  last if $count++ > 100;
}

#--------------------------------------------------------------------------------
# Test $feature->annotation
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# Test $feature->throw
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# Test $feature->warn
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# Test $feature->inflate_result
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# Test $feature->to_gff3_recursive
# Test $feature->get_recursive_children
#--------------------------------------------------------------------------------

ok(my $genes = $schema->features->search({type => 'gene'}),
   '$features = $schema->features->search({type => "gene"})');

$count = 0;
my @lines;
while (my $gene = $genes->next) {
  ok(push @lines, $gene->to_gff3_recursive,      '$gene->to_gff3_recursive');
  my @children;
  ok(! $gene->get_recursive_children(\@children), '$gene->get_recursive_children');
  for my $child (@children) {
    ok($child->type, '$child->type');
  }
  last if $count++ > 10;
}

my $count = scalar grep {$_ !~ /recursive/} @lines;
my $uniq_count = scalar grep {$_ !~ /recursive/} uniq(@lines);
ok($count == $uniq_count, 'to_gff3_recursive lines are uniq');

#--------------------------------------------------------------------------------
# Test $feature->get_transcript_types
#--------------------------------------------------------------------------------

# To get a list of all of the subs and throws:
# Select an empty line and then: C-u M-| grep -nP '^sub ' ../Schema::Result::Feature.pm
# Select an empty line and then: C-u M-| grep -C2 -P '\>throw(' ../Schema::Result::Feature.pm

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
