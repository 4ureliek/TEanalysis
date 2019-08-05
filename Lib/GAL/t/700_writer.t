#!/usr/bin/perl
use strict;

use Test::More tests => 3;

BEGIN {
	use lib '../../';
	#TEST 1
	use_ok('GAL::Writer');
}

my $path = $0;
$path =~ s/[^\/]+$//;
$path ||= '.';
chdir($path);

# TEST 2
my $writer = GAL::Writer->new();
isa_ok($writer, 'GAL::Writer');

can_ok($writer, qw(file
		   fh
		   data_source
		   get_metadata_text
		   get_comment_text
		   get_feature_text
		   get_features_text
		   write_metadata
		   write_comments
		   write_feature
		   write_features
		   write_file
		 ));

################################################################################
################################# Ways to Test #################################
################################################################################

__END__

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
