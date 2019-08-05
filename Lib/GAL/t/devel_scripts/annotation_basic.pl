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

annotation_basic   --mode --format gff3                      \
                   --dsn 'dbi:SQLite:dbname=feature.sqlite'  \
                   genome_annotation.gff3

Description:

A test script for parsing, storing and manipulating an annotation set.

Options:

  parser       The format of the input file (a subclass of GAL::Parser).
  storage      The storage engine to use (a subclass of GAL::Storage).
  fasta        The path to fasta files associated with the annotation.
  mode         The mode for loading files (overwrite, append, if_empty).
  dsn          The DBI connection string (e.g. dbi:mysql:database) for the
	       features
  user         The database user.  Defaults first to \$ENV{\${driver}_user}
	       where \$driver is the driver parsed from the dsn. ['']
  password     The database password.   Defaults first to
	       \$ENV{\${driver}_password} where \$driver is the driver
	       parsed from the dsn.
  prompt       Prompt for the password on the command line.

";

my ($help, $parser, $storage, $access, $mode, $format, $fasta, $dsn, $user, $password, $prompt);
my $opt_success = GetOptions('help'       => \$help,
			     'parser=s'   => \$parser,
			     'storage=s'  => \$storage,
			     'fasta=s'    => \$fasta,
			     'mode=s'     => \$mode,
			     'dsn=s'      => \$dsn,
			     'user=s'     => \$user,
			     'password=s' => \$password,
			     'prompt'     => sub {$password = prompt("Password: ", -tty,  -echo => '*')->{value}},
			      );

die $usage if $help || ! $opt_success;

my $feature_file = shift;

my $storage_args = {class    => $storage,
		    dsn	     => $dsn,
		    user     => $user,
		    password => $password,
		   };

my $parser_args = {class => $parser,
		  };

my $fasta_args = {path => $fasta};

my $feature_store = GAL::Annotation->new(storage => $storage_args,
					 parser  => $parser_args,
					 fasta   => $fasta_args,
					);

if ($feature_file) {
  $feature_store->load_files(files => $feature_file,
			     mode  => $mode,
			    );
    die;
}

my $features = $feature_store->schema;

my $mrnas = $features->resultset('Feature')->search({type => 'mRNA'});
while (my $mrna = $mrnas->next) {
  my $mrna_id = $mrna->feature_id;
  my $CDSs = $mrna->CDSs;
  while (my $CDS = $CDSs->next) {
    my $bins = $CDS->feature_bins;
    my $start = $CDS->start;
    my $end   = $CDS->end;
    print join "\t", ($mrna_id,
		      $CDS->feature_id,
		      $start,
		      $end,
		     );
    print "\n";
    print '';
  }
}

#-----------------------------------------------------------------------------
#-------------------------------- SUBROUTINES --------------------------------
#-----------------------------------------------------------------------------
