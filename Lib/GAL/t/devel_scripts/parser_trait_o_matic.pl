#!/usr/bin/perl
use strict;

use lib '../../';
use GAL::Parser::trait_o_matic;

my $parser = GAL::Parser::trait_o_matic->new(file => 'data/trait_o_matic.gff');

while (my $feature_hash = $parser->next_feature_hash) {
	print $parser->to_gff3($feature_hash) . "\n";
}
