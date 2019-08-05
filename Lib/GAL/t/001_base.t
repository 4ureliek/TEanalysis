#!/usr/bin/perl
use strict;

use Test::More;

BEGIN {
	use lib '../../';
	# TEST 1
	use_ok('GAL::Base');
}

my $path = $0;
$path =~ s/[^\/]+$//;
$path ||= '.';
chdir($path);

# TEST 2
my $base = GAL::Base->new();
isa_ok($base, 'GAL::Base');

# TEST
# ok($base->throw('error_code', 'Test throw'), '$base->throw');

# TEST 3
ok($base->warn('warn_code', 'Test warn'), '$base->warn');
ok($base->info('info_code', 'Test info'), '$base->info');

# TEST 4
my $wrap_text = 'x' x 100;
$wrap_text = $base->wrap_text($wrap_text, 50);
my $wrap_test = 'x' x 50;
$wrap_test .= "\n";
$wrap_test .= 'x' x 50;
ok($wrap_text eq $wrap_test, '$base->wrap_text');

# TEST 5
my $trim_text = ' Hello World ';
$trim_text = $base->trim_whitespace($trim_text);
my $trim_test = 'Hello World';
ok($trim_text = $trim_test, '$base->trim_whitespace');

# TEST 6-9
my %args_hash  = (test => 1);
my @args_array = (test => 2);
my $args_href  = {test => 3};
my $args_aref  = [test => 4];

my $p_args_hash  = $base->prepare_args(%args_hash);
my $p_args_array = $base->prepare_args(@args_array);
my $p_args_href  = $base->prepare_args($args_href);
my $p_args_aref  = $base->prepare_args($args_aref);

ok($p_args_hash->{test}  == 1, '$base->prepare_args(%args_hash)');
ok($p_args_array->{test} == 2, '$base->prepare_args(@args_array)');
ok($p_args_href->{test}  == 3, '$base->prepare_args($args_href)');
ok($p_args_aref->{test}  == 4, '$base->prepare_args($args_aref)');

# TEST
# ok($base->set_attributes, '$base->set_attributes');

# TEST 10
my @nts;
for my $code (qw(A C G T U M R W S Y K V H D B N)) {
	push @nts, $base->expand_iupac_nt_codes($code);
}
my $nts_string = join '', @nts;
is($nts_string,
   'ACGTTACAGATCGCTGTACGACTAGTCGTACGT',
   '$base->expand_iupac_nt_codes($code)');

# TEST 11
ok($base->load_module('GAL::Annotation'), '$base->load_module');

# TEST 12
my $seq = $base->revcomp('AAATTTCCCGGG');
is($seq, 'CCCGGGAAATTT', '$base->revcomp');

# TEST 13
my %feature = (seqid => 'chr1',
	       start => 1,
	       end  => 10000
	      );
my @bins = $base->get_feature_bins(\%feature);
my $bin_text = join ' ', @bins;
is($bin_text, 'chr1:1:0 chr1:5:0 chr1:4:0 chr1:3:0 chr1:2:0', '$base->get_feature_bins');

# TEST 14
my $pep = $base->translate('AATGCCCGACCATTAGCCC', 1, 15);
my $test_pep = 'MPDH*';
is($pep, $test_pep, '$base->translate');

# TEST 15
my $genetic_code = $base->genetic_code;
is($genetic_code->{ATG}, 'M', '$base->genetic_code');

# TEST 16
ok($base->time_stamp =~ /^\d{8}$/, '$base->time_stamp');

# TEST 17
my $rand_text = $base->random_string;
ok($rand_text =~ /^[A-z0-9]{8}$/, '$base->random_string');

# TEST 18-19
ok($base->float_lt(0.000123, 0.000124), '$base->float_lt');
ok(! $base->float_lt(0.000124, 0.000123), '$base->float_lt');

# TEST 20-22
ok($base->float_le(0.000123, 0.000124), '$base->float_le');
ok($base->float_le(0.000123, 0.000123), '$base->float_le');
ok(! $base->float_le(0.000124, 0.000123), '$base->float_le');

# TEST 23-24
ok($base->float_gt(0.000124, 0.000123), '$base->float_gt');
ok(! $base->float_gt(0.000123, 0.000124), '$base->float_gt');

# TEST 25-27
ok($base->float_gt(0.000124, 0.000123), '$base->float_gt');
ok($base->float_ge(0.000123, 0.000123), '$base->float_ge');
ok(! $base->float_gt(0.000123, 0.000124), '$base->float_gt');

ok(my @transcripts = $base->get_transcript_types(), 'my @transcripts = $base->get_transcript_types()');
ok(grep {m/transcript/} @transcripts, 'grep {m/transcript//} @transcripts');
done_testing;

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
