#!/usr/bin/perl
use strict;

use Test::More;

BEGIN {
	use lib '../../';
	use_ok('GAL::Schema::Result::Relationship');
}

my $path = $0;
$path =~ s/[^\/]+$//;
$path ||= '.';
chdir($path);

my $object = GAL::Schema::Result::Relationship->new();
isa_ok($object, 'GAL::Schema::Result::Relationship');

# To get a list of all of the subs and throws:
# Select an empty line and then: C-u M-| grep -nP '^sub ' ../Schema::Result::Relationship.pm
# Select an empty line and then: C-u M-| grep -C2 -P '\>throw(' ../Schema::Result::Relationship.pm

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
