#!/usr/bin/perl
use strict;

use Test::More 'no_plan'; # tests => 10;

BEGIN {
	use lib '../../../';
	#TEST 1
	use_ok('GAL::Parser::annovar_exonic');
}

my $path = $0;
$path =~ s/[^\/]+$//;
$path ||= '.';
chdir($path);

my $parser = GAL::Parser::annovar_exonic->new(file => '../data/NA12878.annovar_input_short.exonic_variant_function');

# TEST 2
isa_ok($parser, 'GAL::Parser::annovar_exonic');

while (my $feature = $parser->next_feature_hash) {
	print $parser->to_gff3($feature) . "\n";
}

# TEST 3
ok($parser->get_features, '$parser->get_features');



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
