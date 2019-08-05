#!/usr/bin/perl

use strict;
use warnings;

use lib '../';

# This script uses full DBIx::Class object inflation.  It takes about
# 2.6 hours to process 4 million features.

use GAL::Schema;
my $schema = GAL::Schema->connect('DBI:mysql:database=pg_chinese_snp');

my $feature_rs    = $schema->resultset('Feature');

while (my $feature = $feature_rs->next) {
	my $id         = $feature->id;
	my $ref_allele = $feature->attributes->find({key => 'reference_allele'})->value;
	print "$id\t$ref_allele\n";
}
