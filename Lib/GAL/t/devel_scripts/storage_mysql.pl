#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use GAL::Annotation;
use IO::Prompt;

#-----------------------------------------------------------------------------
#----------------------------------- MAIN ------------------------------------
#-----------------------------------------------------------------------------

my $usage = "

Synopsis:

gal_02_load_annotations

Description:

Options:

  create     Drop (if exists) and load the database.  Default
	     will just open it.
	     to add to  the existing database if it exists.
  format     The format of the input file (a subclass of GAL::Parser).
  dsn        The DBI connection string (e.g. dbi:mysql:database).
	     Currently only expected to work with MySQL.
  user       The database user.  Defaults first to \$ENV{\${driver}_user}
	     where \$driver is the driver parsed from the dsn. ['']
  password   The database password.   Defaults first to
	     \$ENV{\${driver}_password} where \$driver is the driver
	     parsed from the dsn.
  prompt     Prompt for the password on the command line.

";

my ($help, $create, $format, $dsn, $user, $password, $prompt, $temp_dir);
my $opt_success = GetOptions('help'         => \$help,
			     'create'      => \$create,
			     'format=s'     => \$format,
			     'dsn=s'        => \$dsn,
			     'user=s'       => \$user,
			     'password=s'   => \$password,
			     'prompt'       => sub {$password = prompt("Password: ", -tty,  -echo => '*')->{value}},
			      );

die $usage if $help || ! $opt_success;

my @files = @ARGV;
die $usage unless scalar @files;
die $usage if grep {! -r $_}  @files;

my %storage_args = (class	  => 'mysql',
		    dsn	  => $dsn,
		    user	  => $user,
		    password	  => $password,
		   );

my $annotations = GAL::Annotation->new(storage => \%storage_args,
					    );

if ($create) {
  $annotations->storage->drop_database;

  $annotations->load_file(file   => \@files,
			  format => $format,
			 );
}

my $schema = $annotations->schema;

my $genes = $schema->resultset('Feature')->search({ type => 'gene' });

while (my $gene = $genes->next) {
  print $gene->feature_id . "\n";
  my $mRNAs = $gene->children->search({type => 'mRNA'});
  while (my $mRNA = $mRNAs->next) {
    print "\t" . $mRNA->feature_id . "\n";
    my $CDSs = $mRNA->children->search({type => 'CDS'});
    while (my $CDS = $CDSs->next) {
      print "\t\t";
      print join "\t", ($CDS->feature_id,
			$CDS->end - $CDS->start,
		       );
      print "\n";
    }
  }
}


# my $exons = $gene_annotations->set(type => 'exons');
#
# $exons->cluster_and_create_genes;
#
# my $mRNAs->$gene_annotations->set(type => 'mRNA');
# $mRNAs->infer_features(infer_type => 'CDS');
#
# my $CDSs->$gene_annotations->set(type => 'CDS');
#
# my $cds_footprint = $CDSs->footprint;
#
# my $snp_annotations = GAL::Annotation->new(dsn  => 'dbi:mysql:NA19240_snps');
#
# $snp_annotations->load_file(file    => 'NA19240_snps.gff3',
#			      format  => 'solid',
#			     );
#
# my @genes = $gene_annotations->set(seqid  => 'chr1',
#				     type   => 'gene');
#
# for my $gene (@genes) {
#
#   for my $mRNA ($gene->children(type => 'mRNA')) {
#
#     my $CDSs = $mRNA->children->(type => 'CDS');
#     my $CDS  = $CDSs->longest;
#
#     my $snps = $snp_annotations->overlaps(feature => $CDS,
#					  type    => 'SNP');
#
#     while (my $snp = $snps->next) {
#
#       my $ref_allele = $snp->seq;
#
#       # Only dealing with one variant_allele;
#       my ($variant_allele) = $snp->attributes(variant_allele);
#       my $codons = $gene_annotations->set(type   => 'codon',
#					  parent => $CDS);
#
#       my $codon = $codons->overlaps(feature => $snp,
#				    parent  => $CDS,
#				    type    => 'codon');
#
#       my $codon_position  = $codon->my_coordinate($snp->start);
#       my $ref_codon_seq   = $codon->seq;
#       my $variant_codon_seq = $ref_codon_seq;
#       substr($variant_codon_seq, $codon_position - 1, $snp->length, $variant_allele);
#
#       $snp->variant_effect(feature => $CDS);
#     }
#   }
# }

#-----------------------------------------------------------------------------
#-------------------------------- SUBROUTINES --------------------------------
#-----------------------------------------------------------------------------
