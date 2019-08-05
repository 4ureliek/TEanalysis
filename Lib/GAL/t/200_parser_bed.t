#!/usr/bin/perl
use strict;

use Test::More;

BEGIN {
	use lib '../../';
	use_ok('GAL::Parser::bed');
}

my $path = $0;
$path =~ s/[^\/]+$//;
$path ||= '.';
chdir($path);

my $parser = GAL::Parser::bed->new(file => 'data/bed_test_genes.bed');

isa_ok($parser, 'GAL::Parser::bed');

ok(my $record = $parser->next_record, '$parser->next_record');

ok($parser->parse_record($record), '$parser->parse_record');

while (my $feature = $parser->next_feature_hash) {
  ok($feature, 'feature parses');
}

done_testing();

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
