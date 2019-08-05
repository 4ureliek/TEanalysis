#!/usr/bin/perl

use strict;
use warnings;

use lib '../';

# This script uses DBI directly to access the data.  It can process 4
# million features in 2 minutes.

use DBI;

my $dbh = DBI->connect('DBI:mysql:database=pg_chinese_snp');
my $statement = 'select f.feature_id, a.value from feature f, attribute a where f.feature_id = a.feature_id and a.key = "reference_allele"';
my $sth   = $dbh->prepare($statement);
my $rv = $sth->execute;
while (my @row = $sth->fetchrow_array) {
	print join "\t", @row;
	print "\n";
}
