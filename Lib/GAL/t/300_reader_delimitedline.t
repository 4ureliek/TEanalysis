#!/usr/bin/perl
use strict;

use Test::More tests => 25;

BEGIN {
	use lib '../../';
	#TEST 1
	use_ok('GAL::Reader::DelimitedLine');
}

my $path = $0;
$path =~ s/[^\/]+$//;
$path ||= '.';
chdir($path);

my $reader = GAL::Reader::DelimitedLine->new();
isa_ok($reader, 'GAL::Reader::DelimitedLine');

can_ok($reader, qw(
		    new
		    _initialize_args
		    field_names
		    field_separator
		    end_of_data
		    comment_pattern
		    metadata_pattern
		    header_count
		    next_record
		    headers
		    current_line
		    comments
		    metadata
		 ));

ok($reader->file('./data/delimitedline_test.txt'), '$reader->file');

ok($reader->field_names(qw(seqid source type start end score strand phase
			   attributes)), '$reader->field_names');

ok($reader->end_of_data(qr/^\#\#FASTA/), '$reader->end_of_data(qr//)');

ok($reader->comment_pattern(qr/^\#[^\#]/),  '$reader->comment_pattern(qr/^\#[^\#]/');

ok($reader->comment_parser(sub {my ($r, $p, $l) = @_;
				$l =~ s/^\#\s*//g;
				$l =~ s/\s/_/g;
				return $l}),
  '$reader->comment_parser(sub{})');

ok($reader->metadata_pattern(qr/^\#\#/), '$reader->metadata_pattern(qr//, sub{}');

ok($reader->metadata_parser(   sub {my ($r, $p, $l) = @_;
					       $l =~ s/^\#\#//;
					       my %h = split(/\s+/, $l, 2);
					       return \%h
					     }),
  '$reader->metadata_parser(sub{})');

ok($reader->header_count(1), '$reader->header_count = 1');

ok($reader->next_record, '$reader->next_record');

ok($reader->current_line =~ /^chr22\t/, '$reader->current_line');

ok(grep({$_->{metadata} eq 'value'} $reader->metadata), 'Recovered in-file metadata');
ok(grep({/^\#This_is_an_in-file_comment$/} $reader->comments), 'Recovered in-file comments');
ok(grep({/^seqid\tsource\ttype\tstart\tend\tscore\tstrand\tphase\tattributes$/}
	$reader->headers),
   'Recovered in-file headers');

ok($reader->comments('This is a comment'), '$reader->comments');

ok((my @comments = $reader->comments), '$reader->comments');

ok(grep({/^This is a comment$/} @comments), 'Recovered the given comment');

ok($reader->metadata('key value'), '$reader->metadata($metadata)');

ok(my @metadata = $reader->metadata, '$reader->metadata');

ok(grep({$_->{key} eq 'value'} @metadata), 'Recovered the given metadata');

ok($reader->headers('This is a header'), '$reader->headers');

ok(my @headers = $reader->headers, '$reader->headers');

ok(grep({/^This is a header$/} @headers), 'Recovered the given header');

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
