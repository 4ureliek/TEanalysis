#!/usr/bin/perl
use strict;

use Test::More tests => 11;

BEGIN {
	use lib '../../';
	# TEST 1
	use_ok('GAL::Parser');
}

my $path = $0;
$path =~ s/[^\/]+$//;
$path ||= '.';
chdir($path);

my $parser = GAL::Parser->new(file => 'data/soap_snp.gff');

# TEST 2
isa_ok($parser, 'GAL::Parser');

# TEST
ok(! $parser->annotation, '$parser->annotation');

# TEST
ok($parser->file, '$parser->file');

# TEST
ok($parser->fh, '$parser->fh');

# TEST
ok($parser->reader, '$parser->reader');

# TEST
ok($parser->next_record, '$parser->next_record');

# TEST
ok(my $feat_hash = $parser->next_feature_hash, '$parser->next_feature_hash');

# TEST
ok($parser->to_gff3($feat_hash), '$parser->to_gff3');

# TEST
ok($parser->parse_record, '$parser->parse_record');

# TEST 11
my $attribute_text = 'ID=12345; name=Gene1; Parent=6789,9876;';
ok(my $att_hash = $parser->parse_attributes($attribute_text), '$parser->parse_attributes');

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
