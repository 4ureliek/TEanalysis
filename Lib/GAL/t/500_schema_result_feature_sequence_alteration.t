#!/usr/bin/perl
use strict;

use Test::More;

BEGIN {
	use lib '../../';
	use_ok('GAL::Schema::Result::Feature::sequence_alteration');
	use_ok('GAL::Annotation');
}

my $path = $0;
$path =~ s/[^\/]+$//;
$path ||= '.';
chdir($path);

my $object = GAL::Schema::Result::Feature::sequence_alteration->new();
isa_ok($object, 'GAL::Schema::Result::Feature::sequence_alteration');

system('rm data/10Gen_Chinese_SNV_partial.sqlite') if
  -e 'data/10Gen_Chinese_SNV_partial.sqlite';

ok(my $schema = GAL::Annotation->new('data/10Gen_Chinese_SNV_partial.gvf'),
   '$schema = GAL::Annotation->new("data/10Gen_Chinese_SNV_partial.gvf"');

ok(my $seq_alts = $schema->features, '$features = $schema->features');

ok(my $snvs = $seq_alts->search({type => 'SNV'}),
   '$snvs->search({type => "SNVs"})');

#--------------------------------------------------------------------------------
# Test $seq_alt->variant_seq
#--------------------------------------------------------------------------------

while (my $snv = $snvs->next) {
  ok(my @alts = $snv->variant_seq, '$snv->variant_seq');
  #print join ',', @alts;
  #print "\n";
}

#--------------------------------------------------------------------------------
# Test $seq_alt->variant_seq_no_ref
#--------------------------------------------------------------------------------

$snvs->reset;
while (my $snv = $snvs->next) {
  ok(my @alts = $snv->variant_seq_no_ref, '$snv->variant_seq_no_ref');
  #print join ',', @alts;
  #print "\n";
}

#--------------------------------------------------------------------------------
# Test $seq_alt->reference_seq
#--------------------------------------------------------------------------------

$snvs->reset;
while (my $snv = $snvs->next) {
  ok(my @alts = $snv->reference_seq, '$snv->reference_seq');
  #print join ',', @alts;
  #print "\n";
}

#--------------------------------------------------------------------------------
# Test $seq_alt->variant_reads
#--------------------------------------------------------------------------------

$snvs->reset;
while (my $snv = $snvs->next) {
  ok(my @alts = $snv->variant_reads, '$snv->variant_reads');
  #print join ',', @alts;
  #print "\n";
}

#--------------------------------------------------------------------------------
# Test $seq_alt->total_reads
#--------------------------------------------------------------------------------

$snvs->reset;
while (my $snv = $snvs->next) {
  ok(my @alts = $snv->total_reads, '$snv->total_reads');
  #print join ',', @alts;
  #print "\n";
}

#--------------------------------------------------------------------------------
# Test $seq_alt->genotype
#--------------------------------------------------------------------------------

$snvs->reset;
while (my $snv = $snvs->next) {
  ok(my @alts = $snv->genotype, '$snv->genotype');
  #print join ',', @alts;
  #print "\n";
}

##--------------------------------------------------------------------------------
## Test $seq_alt->variant_effect
##--------------------------------------------------------------------------------
#
#$snvs->reset;
#while (my $snv = $snvs->next) {
#  ok(my @alts = $snv->variant_effect, '$snv->variant_effect');
#  #print join ',', @alts;
#  #print "\n";
#}
#
##--------------------------------------------------------------------------------
## Test $seq_alt->variant_copy_number
##--------------------------------------------------------------------------------
#
#$snvs->reset;
#while (my $snv = $snvs->next) {
#  ok(my @alts = $snv->variant_copy_number, '$snv->variant_copy_number');
#  #print join ',', @alts;
#  #print "\n";
#}
#
##--------------------------------------------------------------------------------
## Test $seq_alt->reference_copy_number
##--------------------------------------------------------------------------------
#
#$snvs->reset;
#while (my $snv = $snvs->next) {
#  ok(my @alts = $snv->reference_copy_number, '$snv->reference_copy_number');
#  #print join ',', @alts;
#  #print "\n";
#}

# To get a list of all of the subs and throws:
# Select an empty line and then: C-u M-| grep -nP '^sub ' ../Schema::Result::Feature::sequence_alteration.pm
# Select an empty line and then: C-u M-| grep -C2 -P '\>throw(' ../Schema::Result::Feature::sequence_alteration.pm

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
