#!/usr/bin/perl
use strict;

use Test::More 'no_plan'; # tests => 10;

BEGIN {
	use lib '../../';
	#TEST 1
	use_ok('GAL::List::Numeric');
}

my $path = $0;
$path =~ s/[^\/]+$//;
$path ||= '.';
chdir($path);

my @rand_list = map {int(rand(1000)) + 1} (1 .. 1000);

# TEST 
my $list_numeric = GAL::List::Numeric->new(list => \@rand_list);
isa_ok($list_numeric, 'GAL::List::Numeric');

# TEST 
ok($list_numeric->count == 1000, '$list->count');

# TEST 
ok($list_numeric->stats, '$list_numeric->stats');

# TEST 
ok($list_numeric->count, '$list_numeric->count');

# TEST 
ok($list_numeric->mean, '$list_numeric->mean');

# TEST 
ok($list_numeric->sum, '$list_numeric->sum');

# TEST 
ok($list_numeric->variance, '$list_numeric->variance');

# TEST 
ok($list_numeric->standard_deviation, '$list_numeric->standard_deviation');

# TEST 
ok($list_numeric->min, '$list_numeric->min');

# TEST 
ok($list_numeric->mindex, '$list_numeric->min_dex');

# TEST 
ok($list_numeric->max, '$list_numeric->max');

# TEST 
ok($list_numeric->maxdex, '$list_numeric->maxdex');

# TEST 
ok($list_numeric->sample_range, '$list_numeric->sample_range');

# TEST 
ok($list_numeric->percentile(25), '$list_numeric->percentile(25)');

# TEST 
ok($list_numeric->median, '$list_numeric->median');

# TEST 
ok($list_numeric->harmonic_mean, '$list_numeric->harmonic_mean');

# TEST 
ok($list_numeric->geometric_mean, '$list_numeric->geometric_mean');

# TEST 
ok($list_numeric->mode, '$list_numeric->mode');

# TEST 
ok($list_numeric->trimmed_mean(0.1,0.1), '$list_numeric->trimmed_mean(0.1,0.1)');

# TEST 
ok($list_numeric->frequency_distribution(10), '$list_numeric->frequency_distribution(10)');

# TEST 
ok($list_numeric->least_squares_fit, '$list_numeric->least_squares_fit');

# TEST 
ok($list_numeric->bins, '$list_numeric->bins');

# TEST 
ok($list_numeric->bin_range, '$list_numeric->bin_range');

# TEST 
ok($list_numeric->fd, '$list_numeric->fd');

# TEST 
ok($list_numeric->cfd, '$list_numeric->cfd');

# TEST 
ok($list_numeric->relative_fd, '$list_numeric->relative_fd');

# TEST 
ok($list_numeric->relative_cfd, '$list_numeric->relative_cfd');


################################################################################
################################# Ways to Test #################################
################################################################################

__END__

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
