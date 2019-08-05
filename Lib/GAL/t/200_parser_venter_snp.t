#!/usr/bin/perl
use strict;

use Test::More;

BEGIN {
	use lib '../../';
	use_ok('GAL::Parser::venter_snp');
}

my $path = $0;
$path =~ s/[^\/]+$//;
$path ||= '.';
chdir($path);

my $parser = GAL::Parser::venter_snp->new(file  => 'data/venter_snp.gff',
					  fasta => 'data/fasta_hg18/chr22_hg18.fasta');

isa_ok($parser, 'GAL::Parser::venter_snp');

ok(my $record = $parser->next_record, '$parser->next_record');

ok($parser->parse_record($record), '$parser->parse_record');

while (my $variant = $parser->next_feature_hash) {
  ok($variant, 'variant parses');
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
