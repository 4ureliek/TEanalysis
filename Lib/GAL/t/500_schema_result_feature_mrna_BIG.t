#!/usr/bin/perl
use strict;

use Test::More;

my $data_path = $ENV{GAL_BIG_TEST};

if( ! $data_path) {
  plan skip_all => 'Big data tests require a data path: export GAL_BIG_TEST="/path/to/GAL/big/data" to run';
}

eval "require Time::HiRes; ; Time::HiRes->import( qw(tv_interval gettimeofday) );";
if ($@) {
  plan skip_all => 'Time::HiRes required for test';
}

use lib '../../';
use_ok('GAL::Annotation');
use_ok('GAL::Schema::Result::Feature::mrna');

my $path = $0;
$path =~ s/[^\/]+$//;
$path ||= '.';
chdir($path);

my $gff3_file  = $data_path . '/gff/refGene_hg19.gff3';
my $fasta_file = $data_path . '/fasta/chrAll.fa';

my $annotation = GAL::Annotation->new($gff3_file,
				      $fasta_file);

my $features = $annotation->features;

# Do the search and get an interator for all matching features
my $mrnas = $features->search({type => 'mRNA'});

# Iterate over the features

my $t0 = [gettimeofday()];
my $count = 1;
while (my $mrna = $mrnas->next) {
  my $length = $mrna->length;
  my $elapsed = tv_interval($t0, [gettimeofday()]);
  my $avg = sprintf("%.4f", $elapsed/$count++);
  print "$count: $avg\r";
  #print "$length\n";
  #print '';
}

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
