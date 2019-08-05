#!/usr/bin/perl
use strict;

use Test::More 'no_plan'; # tests => 10;

BEGIN {
	use lib '../../';
	#TEST 1
	use_ok('GAL::List::Numeric');
	#TEST 2
	use_ok('GAL::List::Categorical');
}

my $path = $0;
$path =~ s/[^\/]+$//;
$path ||= '.';
chdir($path);

# TEST 2
my $list_categorical = GAL::List::Categorical->new(list => [qw(red red red blue blue
							       green yellow orange orange
							       purple purple purple purple)]);
isa_ok($list_categorical, 'GAL::List::Categorical');

# TEST 
ok($list_categorical->count == 13, '$list_categorical->count');

# TEST 
my $list_categorical_ref = $list_categorical->list;
ok(ref $list_categorical_ref eq 'ARRAY', '$x = $list_categorical->list');

# TEST 
my @list_ary = $list_categorical->list;
ok(scalar @list_ary == 13, '@x = $list_categorical->list');

# TEST 
my $category_counts = $list_categorical->category_counts;
ok($category_counts->{blue} == 2, '$list_categorical->category_counts');

# TEST 
ok($list_categorical->count_uniq == 6, '$list_categorical->count_uniq');

# TEST 
ok($list_categorical->maxstr eq 'yellow', '$list_categorical->maxstr');

# TEST 
ok($list_categorical->minstr eq 'blue', '$list_categorical->minstr');

# TEST 
my $list_categorical_text = join ' ', $list_categorical->list;
my $list_categorical_shuffle = join ' ', $list_categorical->shuffle;
ok($list_categorical_text ne $list_categorical_shuffle, '$list_categorical->shuffle');

# TEST 
my $list_categorical_uniq = join ' ', sort $list_categorical->uniq;
my $list_categorical_uniq_test = 'blue green orange purple red yellow';
ok($list_categorical_uniq eq $list_categorical_uniq_test, '$list_categorical->uniq');

# TEST 
ok($list_categorical->random_pick(1), '$list_categorical->random_pick');

# TEST 
my $list_numeric = GAL::List::Numeric->new(list => [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]);
isa_ok($list_numeric, 'GAL::List::Numeric');

# TEST 
ok($list_numeric->max == 9, '$list_numeric->max');

# TEST 
ok($list_numeric->min == 0, '$list_numeric->min');

# TEST 
ok($list_numeric->sum == 45, '$list_numeric->sum');


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
