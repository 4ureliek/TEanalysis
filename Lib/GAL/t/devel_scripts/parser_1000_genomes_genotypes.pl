#!/usr/bin/perl
use strict;

use lib '../../';
use GAL::Parser::1000_genomes_genotypes;

my $parser = GAL::Parser::1000_genomes_genotypes->new(file => '../data/CEU.low_coverage.2010_07.genotypes.short.vcf');

while (my $feature_hash = $parser->next_feature_hash) {
	print $parser->to_gff3($feature_hash) . "\n";
}
