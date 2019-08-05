#!/usr/bin/perl
use strict;

use lib '../../';
use GAL::Parser::samtools_pileup;

my $parser = GAL::Parser::samtools_pileup->new(file => 'data/samtools_pileup.txt');

while (my $feature_hash = $parser->next_feature_hash) {
	print $parser->to_gff3($feature_hash) . "\n";
}
