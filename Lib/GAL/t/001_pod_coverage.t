#!/usr/bin/perl

use strict;
use warnings;

use Test::More;
#use Test::Pod::Coverage;

plan(skip_all => '001_pod_coverage.t under developement');

# Select an empty line and use C-u M-| to update the list of tests;
# find ../../ -name '*.pm' | sort | perl -ane 'chomp;s/^[\.|\/]*//;s/\.pm$//;s/\//::/g;print "pod_coverage_ok(\"$_\", \"$_ POD is covered.\");\n"' | wc

# Use M-! to count the number of tests
#find ../../ -name '*.pm' | sort | perl -ane 'chomp;s/^[\.|\/]*//;s/\.pm$//;s/\//::/g;print "pod_coverage_ok(\"$_\", \"$_ POD is covered.\");\n"' | wc

pod_coverage_ok("GAL::Annotation", "GAL::Annotation POD is covered");
pod_coverage_ok("GAL::Base", "GAL::Base POD is covered");
pod_coverage_ok("GAL::Feature", "GAL::Feature POD is covered");
pod_coverage_ok("GAL::Feature::sequence_alteration", "GAL::Feature::sequence_alteration POD is covered");
pod_coverage_ok("GAL::Feature::sequence_feature", "GAL::Feature::sequence_feature POD is covered");
pod_coverage_ok("GAL::HOLD::FeatureFactory", "GAL::HOLD::FeatureFactory POD is covered");
pod_coverage_ok("GAL::HOLD::Schema", "GAL::HOLD::Schema POD is covered");
pod_coverage_ok("GAL::Index", "GAL::Index POD is covered");
pod_coverage_ok("GAL::List", "GAL::List POD is covered");
pod_coverage_ok("GAL::List::Categorical", "GAL::List::Categorical POD is covered");
pod_coverage_ok("GAL::List::Numeric", "GAL::List::Numeric POD is covered");
pod_coverage_ok("GAL::List::Span", "GAL::List::Span POD is covered");
pod_coverage_ok("GAL::Parser", "GAL::Parser POD is covered");
pod_coverage_ok("GAL::Parser::basic_snp", "GAL::Parser::basic_snp POD is covered");
pod_coverage_ok("GAL::Parser::complete_genomics", "GAL::Parser::complete_genomics POD is covered");
pod_coverage_ok("GAL::Parser::complete_genomics_new", "GAL::Parser::complete_genomics_new POD is covered");
pod_coverage_ok("GAL::Parser::cosmic", "GAL::Parser::cosmic POD is covered");
pod_coverage_ok("GAL::Parser::dbsnp_flat", "GAL::Parser::dbsnp_flat POD is covered");
pod_coverage_ok("GAL::Parser::gff3", "GAL::Parser::gff3 POD is covered");
pod_coverage_ok("GAL::Parser::illumina_indel", "GAL::Parser::illumina_indel POD is covered");
pod_coverage_ok("GAL::Parser::illumina_snp", "GAL::Parser::illumina_snp POD is covered");
pod_coverage_ok("GAL::Parser::maq_cns2snp", "GAL::Parser::maq_cns2snp POD is covered");
pod_coverage_ok("GAL::Parser::na18507_sanger_indel", "GAL::Parser::na18507_sanger_indel POD is covered");
pod_coverage_ok("GAL::Parser::quake", "GAL::Parser::quake POD is covered");
pod_coverage_ok("GAL::Parser::samtools_pileup", "GAL::Parser::samtools_pileup POD is covered");
pod_coverage_ok("GAL::Parser::soap_indel", "GAL::Parser::soap_indel POD is covered");
pod_coverage_ok("GAL::Parser::soap_snp", "GAL::Parser::soap_snp POD is covered");
pod_coverage_ok("GAL::Parser::solid", "GAL::Parser::solid POD is covered");
pod_coverage_ok("GAL::Parser::template", "GAL::Parser::template POD is covered");
pod_coverage_ok("GAL::Parser::template_sequence_alteration", "GAL::Parser::template_sequence_alteration POD is covered");
pod_coverage_ok("GAL::Parser::trait_o_matic", "GAL::Parser::trait_o_matic POD is covered");
pod_coverage_ok("GAL::Parser::venter_indel", "GAL::Parser::venter_indel POD is covered");
pod_coverage_ok("GAL::Parser::venter_snp", "GAL::Parser::venter_snp POD is covered");
pod_coverage_ok("GAL::Parser::watson_cshl", "GAL::Parser::watson_cshl POD is covered");
pod_coverage_ok("GAL::Reader", "GAL::Reader POD is covered");
pod_coverage_ok("GAL::Reader::RecordParser", "GAL::Reader::RecordParser POD is covered");
pod_coverage_ok("GAL::Reader::DelimitedLine", "GAL::Reader::DelimitedLine POD is covered");
pod_coverage_ok("GAL::Reader::Template", "GAL::Reader::Template POD is covered");
pod_coverage_ok("GAL::Schema", "GAL::Schema POD is covered");
pod_coverage_ok("GAL::Schema::Result::Attribute", "GAL::Schema::Result::Attribute POD is covered");
pod_coverage_ok("GAL::Schema::Result::Feature", "GAL::Schema::Result::Feature POD is covered");
pod_coverage_ok("GAL::Schema::Result::Feature::cds", "GAL::Schema::Result::Feature::cds POD is covered");
pod_coverage_ok("GAL::Schema::Result::Feature::exon", "GAL::Schema::Result::Feature::exon POD is covered");
pod_coverage_ok("GAL::Schema::Result::Feature::five_prime_utr", "GAL::Schema::Result::Feature::five_prime_utr POD is covered");
pod_coverage_ok("GAL::Schema::Result::Feature::gene", "GAL::Schema::Result::Feature::gene POD is covered");
pod_coverage_ok("GAL::Schema::Result::Feature::intron", "GAL::Schema::Result::Feature::intron POD is covered");
pod_coverage_ok("GAL::Schema::Result::Feature::mrna", "GAL::Schema::Result::Feature::mrna POD is covered");
pod_coverage_ok("GAL::Schema::Result::Feature::protein", "GAL::Schema::Result::Feature::protein POD is covered");
pod_coverage_ok("GAL::Schema::Result::Feature::sequence_alteration", "GAL::Schema::Result::Feature::sequence_alteration POD is covered");
pod_coverage_ok("GAL::Schema::Result::Feature::sequence_feature", "GAL::Schema::Result::Feature::sequence_feature POD is covered");
pod_coverage_ok("GAL::Schema::Result::Feature::template", "GAL::Schema::Result::Feature::template POD is covered");
pod_coverage_ok("GAL::Schema::Result::Feature::three_prime_utr", "GAL::Schema::Result::Feature::three_prime_utr POD is covered");
pod_coverage_ok("GAL::Schema::Result::Feature::transcript", "GAL::Schema::Result::Feature::transcript POD is covered");
pod_coverage_ok("GAL::Schema::Result::Relationship", "GAL::Schema::Result::Relationship POD is covered");
pod_coverage_ok("GAL::SchemaAnnotation", "GAL::SchemaAnnotation POD is covered");
pod_coverage_ok("GAL::Storage", "GAL::Storage POD is covered");
pod_coverage_ok("GAL::Storage::SQLite", "GAL::Storage::SQLite POD is covered");
pod_coverage_ok("GAL::Storage::mysql", "GAL::Storage::mysql POD is covered");
