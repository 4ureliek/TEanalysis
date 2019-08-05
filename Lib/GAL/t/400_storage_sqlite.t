#!/usr/bin/perl
use strict;

use Test::More tests => 8;

BEGIN {
	use lib '../../';
	#TEST 1
	use_ok('GAL::Storage::SQLite');
}

my $path = $0;
$path =~ s/[^\/]+$//;
$path ||= '.';
chdir($path);

# TEST 2
my $storage = GAL::Storage::SQLite->new(dsn => 'DBI:SQLite:data/test_storage.sqlite');
isa_ok($storage, 'GAL::Storage::SQLite');

# TEST 3
ok(! $storage->drop_database, '$storage->drop_database');

# TEST 4
ok($storage->scheme, '$storage->scheme');

# TEST 5
ok($storage->driver, '$storage->driver');

# TEST 6
ok($storage->database, '$storage->database');

# TEST 7
ok($storage->dbh, '$storage->dbh');

# TEST 8
ok($storage->index_database, '$storage->index_database');


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
