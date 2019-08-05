#!/usr/bin/perl
use strict;

use lib '../../';
use GAL::Parser::smith_ant_snps;

my $parser = GAL::Parser::smith_ant_snps->new(file => '../data/smith_ant_snps.gff',
					      fasta => '../data/fasta_ant/');

while (my $feature_hash = $parser->next_feature_hash) {
	print $parser->to_gff3($feature_hash) . "\n";
}
