#!/usr/bin/perl
use strict;
use warnings;

use GAL::Annotation;

my $annotation = GAL::Annotation->new(dsn     => 'DBI:SQLite:database=gal_annotation_test',
				      parser  => 'gff3',
				      storage => 'SQLite',
				     );

$annotation->load_file('./data/soap_snp.gff');

while (my $feature = $annotation->next_feature) {
        my $id         = $feature->id;
        print "$id\n";
}
