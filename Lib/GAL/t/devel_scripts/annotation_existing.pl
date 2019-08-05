#!/usr/bin/perl
use strict;
use warnings;

use GAL::Annotation;

my $annotation = GAL::Annotation->new(dsn => 'DBI:mysql:pg_chinese_snp');

my $feature_rs = $annotation->schema->resultset('Feature');

while (my $feature = $feature_rs->next) {
        my $id         = $feature->id;
        my $ref_allele = $feature->attributes->find({att_key => 'reference_allele'})->att_value;
        print "$id\t$ref_allele\n";
}

