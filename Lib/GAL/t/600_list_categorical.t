#!/usr/bin/perl
use strict;

use Test::More 'no_plan'; # tests => 10;

BEGIN {
	use lib '../../';
	#TEST 1
	use_ok('GAL::List::Categorical');
}

my $path = $0;
$path =~ s/[^\/]+$//;
$path ||= '.';
chdir($path);

my $list = GAL::List::Categorical->new(list => [qw(red red red blue
						   blue green yellow
						   orange orange
						   purple purple
						   purple purple)]);
# TEST 2
#isaok($list, 'GAL::List::Categorical');

# TEST 3
ok($list->count == 13, '$list->count');

# TEST 4
ok($list->cardinality == 6, '$list->cardinality');

# TEST 5
ok($list->count_uniq == 6, '$list->count_uniq');

# TEST 6
ok($list->maxstr eq 'yellow', '$list->maxstr');

# TEST 7
ok($list->minstr eq 'blue' , '$list->minstr');

# TEST 8
my @random = $list->shuffle;
ok(scalar @random == 13, '$list->shuffle');

# TEST 9
my @uniq = $list->uniq;
ok(scalar @uniq == 6, '$list->uniq');

# TEST 10
my ($random_element) = $list->random_pick(1);
ok($random_element =~ /^(red|blue|green|yellow|orange|purple)$/,
   '$list->random_pick');

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
