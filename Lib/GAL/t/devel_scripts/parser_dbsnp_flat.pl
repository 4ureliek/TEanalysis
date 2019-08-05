#!/usr/bin/perl
use strict;
use warnings;

use lib '../../';
use GAL::Parser::dbsnp_flat;

my $parser = GAL::Parser::dbsnp_flat->new(file     => 'data/dbsnp_flat.txt',
					  assembly => 'reference',
					 );

while (my $feature = $parser->next_feature_hash) {
	print $parser->to_gff3($feature);
	print '';
}
