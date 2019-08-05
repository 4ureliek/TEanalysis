#!/usr/bin/perl
use strict;

use Test::More tests => 8;

BEGIN {
	use lib '../../';
	#TEST 1
	use_ok('GAL::Storage');
}

my $path = $0;
$path =~ s/[^\/]+$//;
$path ||= '.';
chdir($path);

# TEST 2
my $storage = GAL::Storage->new();
isa_ok($storage, 'GAL::Storage');

# TEST 3
ok($storage->dsn(), '$storage->dsn()');

# TEST 4
ok($storage->scheme() eq 'DBI', '$storage->scheme()');

# TEST 5
ok($storage->user('fred') eq 'fred', '$storage->user(\'fred\')');

# TEST 6
ok($storage->password('12345') == 12345, '$storage->password(\'12345\')');

# TEST 7
ok($storage->database() =~ /gal_database_\d{8}_[A-z0-9]/, '$storage->database()');

# TEST 8
ok($storage->driver() eq 'SQLite', '$storage->driver()');

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
