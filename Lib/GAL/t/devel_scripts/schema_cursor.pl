#!/usr/bin/perl

use strict;
use warnings;

use lib '../';

# This script uses a resultset cursor to access the data.  It can
# process 4 million features in 1.3 hours.

use GAL::Schema;
my $schema = GAL::Schema->connect('DBI:mysql:database=pg_chinese_snp');

my $feature_rs    = $schema->resultset('Feature');
my $attributes_rs = $schema->resultset('Attribute');

my $feature_cursor = $feature_rs->cursor;
while (my @feature_columns = $feature_cursor->next) {
	my $id = $feature_columns[0];
	my $ref_allele_rs = $attributes_rs->search({feature_id => $id,
						    key    => 'reference_allele'});
	my $ref_allele_cursor = $ref_allele_rs->cursor;
	my @ref_allele_columns = $ref_allele_cursor->next;
	my $ref_allele = $ref_allele_columns[3];
	print "$id\t$ref_allele\n";
}
