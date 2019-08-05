#!/usr/bin/perl

use strict;
use warnings;

use lib '../';
use DBIx::Class::ResultClass::HashRefInflator;

# This script uses DBIx::Class::ResultClass::HashRefInflator to skip
# object creation.  It can process 4 million features in about 1.6
# hours

use GAL::Schema;
my $schema = GAL::Schema->connect('DBI:mysql:database=pg_chinese_snp');

my $feature_rs    = $schema->resultset('Feature');
my $attributes_rs = $schema->resultset('Attribute');

$feature_rs->result_class('DBIx::Class::ResultClass::HashRefInflator');

while (my $feature = $feature_rs->next) {
	my $id         = $feature->{feature_id};

	my $ref_allele_rs = $attributes_rs->search({feature_id => $id,
						    key    => 'reference_allele'});
	$ref_allele_rs->result_class('DBIx::Class::ResultClass::HashRefInflator');
	my $ref_allele_hash = $ref_allele_rs->first;
	my $ref_allele = $ref_allele_hash->{value};

	print "$id\t$ref_allele\n";
}
