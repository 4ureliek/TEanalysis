=pod

=head1 Summary

The Genome Annotation Library (GAL) is a object-oriented Perl
programming library (and a growing collection of associated
scripts) that is designed to provide a simple and intuitive
interface for accomplishing Sequence Ontology empowered analyses of
genomic annotation data.

=head1 Synopsis

A very simple GAL script that infers introns for mRNAs is seen here:

  use GAL::Annotation;
  my $features = GAL::Annotation->new($gff_file, $fasta_file)->features;

  # Search for mRNAs and get an iterator;
  my $mrnas = $features->search({type => 'mRNA'});
  print STDERR 'Processing ' . $mrnas->count . ' mRNA features';
  # Iterate mRNAs
  while (my $mrna = $mrnas->next) {
    # Infer (if not existing) and return iterator
    my $introns = $mrna->introns;
    while (my $intron = $introns->next) {
      # Use the intron object for access to data
      my $id    = $intron->feature_id;
      my ($start, $end) = $intron->get_values(qw(start end));
      print $intron->fasta if $intron->gc_content > 50;
      print STDERR join ', ', ($feature_id,
			       $intron->locus,
			       $intron->gc_content
			      );
    }
  }

=head1 Description

GAL provides a collection of object classes for sequence features
who's inheritance structure is based on the is_a releationships
between sequence features as described in the Sequence Ontology
(e.g. an mRNA is_a type of transcript and can thus inherit it's
methods).  In addition feature classes provide methods which will
retrieve iterators for child features as described by the SO part_of
relationships and captured by the GFF3 Parent attribute. For
example, transcript objects can return iterators for the exons which
are 'part_of' them.  mRNA objects inherit this behavior to return
exon iterators from the transcript class and in addition can return
an iterator for the CDS features that are 'part_of' them.

=head1 Storage

GAL uses an SQLite relational database backend to store parsed
sequence features.  This allows flexible searching of sequence
features and their attributes.  Results of these searches are
returned as iterators, that, in addition to iterating over the
returned features, provide useful methods for summarizing details of
their lists.  This flexibility in searching and random access to the
data comes at the cost of additional time needed to load, index and
query the database.  Sometimes, all that you need is the ability to
parse and iterate over the entire file one feature at a time.  In
these cases the parser classes can be used directly, as described
below, and features will be rapidly parsed and returned (in a simple
hash sturcture rather than feature objects).

=head1 Class Structure

GAL uses DBIx::Class to provide sequence feature objects from the
relational database backend.

=head1 Parsers

While GAL is specifically designed for parsing GFF3 as the input
format for seuqence features, an abstract Parser class allows
straight-forward subclassing such that any input file, where the
contents can be converted into genome locatable sequence features,
can can be parsed if an appropriate parser subclass is written.  In
addition to the default ability to parse GFF3, several other parser
classes are provided with GAL, and new classes are easily added.

=head1 Using GAL

As described above you can use the full functionality of GAL by
allowing it to load features into a relational database or you can use
light-weight parser classes directly to rapid sequential access to the
data in hash structures.

In the first case a GAL::Annotation object is used to access feature
ojects:

  use GAL::Annotation;
  my $features = GAL::Annotation->new($gff_file, $fasta_file)->features;
  while (my $feature = $features->next) {
      print $feature->feature_id . "\n";
      print $feature->seq . "\n";
  }

In the second case a GAL::Parser object (or subclass) is used to
access features as hash data structures:

  use GAL::Parser;
  my $parser = GAL::Parser->new($gff_file);
  while (my $feature = $parser->next_feature_hash) {
      # Not $feature is a hash reference not an object!
      print $feature->{feature_id} . "\n";
      # But this doesn't work $feature->{seq} - no OO sugar.
  }

=head1 GAL Classes

The following index can be used to explore the documentation in
the GAL library.

L<GAL::Annotation> - The primary class used to initiate a GAL-based
annotation analysis.  Start here is you're exploring the library for
the first time and want access to the full OO feature classes and DB
backend.  Provided a GFF3 file objects of this class will load the
features into a database, index the database and can return an
iterator object based to the features in the data base.

L<GAL::Base> - A base class providing shared utility methods for all
objects library-wide.  This class is not instantiated directly by the
user.

L<GAL::Feature> - Experimental and not currently in use.

L<GAL::Parser> - Base class for all L<GAL::Parser> subclasses.  By
default it parsers GFF3/GVF format.  L<GAL::Parser> is used internally
by the GAL::Annotation object to parse feature files and can also be
used directly to access features as hash references
($parser->next_feature_hash).

L<GAL::Parser::gff3> - A parser class used to parse GFF3 file format.
This is the default parser for the library.

L<GAL::Parser::gff3_fast> - A fast GFF3 parser class that makes the
assumption that the incoming GFF3 file is well formated and does no
validation.  It may break in unfreindly ways if the data is not well
formated.

L<GAL::Parser::ncbi_blast_tab> - A parser for NCBI Blast+ output that
is formated with the

L<GAL::Parser::VCFv4_0> - A parser for VCF v4.0 files

L<GAL::Parser::basic_snp> - A parser for SNVs in basic tab-delimited
format.

L<GAL::Parser::hgmd_indel> - A parser for indel flat files from HGMD.

L<GAL::Parser::hgmd_snv.html> - A parser for SNV flat files from HGMD.

L<GAL::Parser::maq_cns2snp.html> - A parser for MAQ cns2snp output
files.

L<GAL::Parser::samtools_pileup.html> - A parser for samtools pileup
output.

L<GAL::Parser::ucsc_gene_table.html> - A parser for UCSC refGene (and
other) tab-delimited files.

L<GAL::Parser::ucsc_kg_table.html> - A parser for UCSC knownGene
tab-delimited files.

L<GAL::Parser::template.html> - A template to use as a starting point
for writing new GAL::Parser subclasses for sequence feature formats.

L<GAL::Storage> - Provides backend storage in an SQLite database.

L<GAL::Storage::SQLite> - Provides backend storage in an SQLite
database.

L<GAL::Schema> - A class that provides basic configuration to
DBIx::Class.

L<GAL::Schema::Result::Feature> - A DBIx::Class subclass that
describes the structure of the feature table in the database.  The
class also impliments methods that are inherited by
GAL::Schema::Result::Feature subclasses.

L<GAL::Schema::Result::Feature::sequence_feature> - A basic
sequence_feature class that provides methods for all subclasses.  The
methods are actually implimented in GAL::Schema::Result::Feature, but
inherited here so that all sequence_feature types for which no
subclass has been written will inherit all basic sequence_feature
behavior from this class.

L<GAL::Schema::Result::Feature::gene> - A class to provide gene
specific methods.

L<GAL::Schema::Result::Feature::transcript> - A class to provide
transcript specific methods.

L<GAL::Schema::Result::Feature::mrna> - A class to provide mRNA
specific methods.

L<GAL::Schema::Result::Feature::exon> - A class to provide exon
specific methods.

L<GAL::Schema::Result::Feature::intron> - A class to provide intron
specific methods.

L<GAL::Schema::Result::Feature::cds> - A class to provide CDS specific
methods.

L<GAL::Schema::Result::Feature::protein> - A class to provide protein
specific methods.

L<GAL::Schema::Result::Feature::five_prime_utr> - A class to provide
five_prime_UTR specific methods.

L<GAL::Schema::Result::Feature::three_prime_utr> - A class to provide
three_prime_UTR specific methods.

L<GAL::Schema::Result::Feature::sequence_alteration> - A class to
provide sequence_alteration specific methods.

L<GAL::Schema::Result::Feature::template> - A template class for use
in writing new GAL::Schema::Result::Feature subclasses.

L<GAL::Schema::Result::Attribute>

L<GAL::Schema::Result::Relationship>

L<GAL::Interval> - Classes that provide interval (sequence ranges)
aggregation and analysis functions for GAL.  Most of these are
currently in developement and thus not may not be functional or
stable.

L<GAL::Interval::Span>

L<GAL::List> - List aggregation and analysis methods

L<GAL::List::Categorical> - List aggregation and analysis functions
for categorical lists

L<GAL::List::Numeric> - List aggregation and analysis functions for
numerical lists
