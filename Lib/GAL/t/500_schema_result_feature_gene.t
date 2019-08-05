#!/usr/bin/perl
use strict;

use Test::More;

BEGIN {
	use lib '../../';
	use_ok('GAL::Schema::Result::Feature::gene');
	use_ok('GAL::Annotation');
}

my $path = $0;
$path =~ s/[^\/]+$//;
$path ||= '.';
chdir($path);

my $object = GAL::Schema::Result::Feature::gene->new();
isa_ok($object, 'GAL::Schema::Result::Feature::gene');

system('rm data/dmel-4-r5.46.genes.sqlite') if
  -e 'data/dmel-4-r5.46.genes.sqlite';

ok(my $schema = GAL::Annotation->new('data/dmel-4-r5.46.genes.gff',
                                     'data/dmel-4-chromosome-r5.46.fasta'),
   '$schema = GAL::Annotation->new("dmel.gff", "dmel.fasta"');


ok(my $features = $schema->features, '$features = $schema->features');

ok(my $genes = $features->search({type => 'gene'}),
   '$genes->search({type => "gene"})');

#-----------------------------------------------------------------------------------------
# Test $gene->transcripts
#-----------------------------------------------------------------------------------------

while (my $gene = $genes->next) {
  ok(my $transcripts = $gene->transcripts, '$transcripts = $gene->transcripts');
  while (my $transcript = $transcripts->next) {
    ok($transcript->start =~ /^\d+$/, '$transcript->start =~ /^\d+$/');
  }
}

#-----------------------------------------------------------------------------------------
# Test $gene->infer_introns
#-----------------------------------------------------------------------------------------

my $count;
$genes->reset;
while (my $gene = $genes->next) {
  ok(! $gene->infer_introns, '$gene->infer_introns');
  last if ++$count > 10;
}

#-----------------------------------------------------------------------------------------
# Test $gene->splice_complexity
#-----------------------------------------------------------------------------------------

my %correct_spclx = (FBgn0025740 => 0.0236879432624113,
		     FBgn0004859 => 0.0812490391918168,
		     FBgn0264617 => 0.282275592262552,
		     FBgn0017545 => 0.125650529356979,
		     FBgn0085432 => 0.275248541992092,
		     FBgn0011747 => 0.0132822280969565,
		     FBgn0052000 => 0.318241826073237,
		     FBgn0053978 => 0.645084857326208,
		     FBgn0039889 => 0.0156105699049662,
		     FBgn0039890 => 0.0775155459744926,
		     FBgn0051998 => 0.0811987276075674,
);

$count = 0;
$genes->reset;
while (my $gene = $genes->next) {
  next unless $gene->transcripts->count > 1;
  ok(my $spclx = $gene->splice_complexity, '$gene->splice_complexity');
  my $feature_id = $gene->feature_id;
  if (exists $correct_spclx{$feature_id}) {
    ok($correct_spclx{$feature_id} eq $spclx, '$correct_spclx{$feature_id} == $spclx');
  }
  last if ++$count > 10;
}

#-----------------------------------------------------------------------------------------
# Test $gene->mRNAs
#-----------------------------------------------------------------------------------------

while (my $gene = $genes->next) {
  ok(my $mRNAs = $gene->mRNAs, '$mRNAs = $gene->mRNAs');
  while (my $mRNA = $mRNAs->next) {
    ok($mRNA->start =~ /^\d+$/, '$mRNA->start =~ /^\d+$/');
  }
}

#-----------------------------------------------------------------------------------------
# Test $gene->is_coding
#-----------------------------------------------------------------------------------------

my %correct_coding = (FBgn0040037 => 1,
		      FBgn0052011 => 0,
		      FBgn0052010 => 0,
		      FBgn0052009 => 0,
		      FBgn0025740 => 1,
		      FBgn0004859 => 1,
		      FBgn0264617 => 0,
		      FBgn0017545 => 1,
		      FBgn0263851 => 0,
		      FBgn0085432 => 1,
		     );

$genes->reset;
while (my $gene = $genes->next) {
  my $is_coding = $gene->is_coding;
  ok($is_coding =~ /^[0-1]$/, '$is_coding =~ /^[0-1]$/');
  my $feature_id = $gene->feature_id;
  if (exists $correct_coding{$feature_id}) {
    ok($correct_coding{$feature_id} eq $is_coding, '$correct_coding{$feature_id} == $is_coding');
  }
  last if ++$count > 20;
}

#-----------------------------------------------------------------------------------------
# Test $gene->get_transcript_types
#-----------------------------------------------------------------------------------------

$genes->reset;
while (my $gene = $genes->next) {
  ok($gene->get_transcript_types, '$gene->get_transcript_types');
}

# To get a list of all of the subs and throws:
# Select an empty line and then: C-u M-| grep -nP '^sub ' ../Schema::Result::Feature::gene.pm
# Select an empty line and then: C-u M-| grep -C2 -P '\>throw(' ../Schema::Result::Feature::gene.pm

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
