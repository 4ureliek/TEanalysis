package GAL::Annotation;

use strict;
use warnings;

use base qw(GAL::Base);
use GAL::Schema;
use Scalar::Util qw(weaken);

use vars qw($VERSION);
$VERSION = 0.2.1;

=head1 NAME

L<GAL::Annotation> - The Genome Annotation Library

=head1 VERSION

This document describes GAL::Annotation version 0.2.1

=head1 SYNOPSIS

    use GAL::Annotation;

    my $annotation = GAL::Annotation->new('annotations.gff3',
					  'genome.fa');

    # Get a smart iterator object.
    my $features = $annotation->features;
    my $feature_count = $features->count;
    my @types = $features->types->uniq;

    # Do a search and get iterator of features matching the
    # search criteria.
    my $mrnas = $features->search({type => 'mRNA'});

    # Iterate over the features.
    while (my $mrna = $mrnas->next) {
      # Get the feature ID
      my $mrna_id = $mrna->feature_id;
      my $protein_seq = $mrna->protein_seq;
      # Get all the exons for this mRNA
      my $exons = $mrna->exons;
      while (my $exon = $exons->next) {
	my $length = $exon->length;
	my $gc_content = $exon->gc_content;
      }
      # Introns don't exist in the dataset, so GAL
      # will infer them on the fly.
      my $introns = $mrna->introns;
      while (my $intron = $introns->next) {
        my $seq = $intron->seq;
      }
    }

=head1 DESCRIPTION

L<The Genome Annotation Library
(GAL)|http://www.sequenceontology.org/software/GAL.html> is a
collection of object oriented L<Perl|http://www.perl.org> modules
developed by the L<Sequence Ontology|http://www.sequenceontology.org>
that make working with genome annotations simple and intuitive.

GAL was designed to work with sequence features stored in GFF3 files
(or sequence alterations stored in GVF files) and GAL provides many
parsers for converting sequence features in other formats to GFF3 and
GVF.

Sequence features are represented as objects in GAL based on the
relationships in the L<Sequence
Ontology|http://www.sequenceontology.org/browser/obob.cgi>.  This
allows the objects to inherit the appropriate behavior from their
parent types based on the is_a relationships within the ontology and
it also allows the programmer to easily traverse the part_of
relationships while writing code.

=head1 INHERITS FROM

L<GAL::Base>

=head1 INHERITED BY

None

=head1 USES

L<GAL::Schema>

=head1 CONSTRUCTOR

New Annotation objects are created by the class method new.  Arguments
should be passed to the constructor as a list (or reference) of key
value pairs.  All attributes of the Annotation object can be set in
the call to new, but reasonable defaults will be used where ever
possilbe to keep object creation simple.  An simple example of object
creation would look like this:

    my $feat_store = GAL::Annotation->new($gff_file);
    my $feat_store = GAL::Annotation->new($gff_file, $fasta_file);

The resulting object would use a GFF3 parser and SQLite storage by
default The first example would not have access to feature sequence,
the second one would.

A more complex object creation might look like this:

    my $feat_store = GAL::Annotation->new(parser  => {class => gff3},
					  storage => {class => mysql,
						      dsn   => 'dbi:SQLite:database'
						      user  => 'me',
						      password => 'secret'
					  fasta   =>  '/path/to/fasta/files/'
					  );

The constructor recognizes the following parameters which will set the
appropriate attributes:

=over 4

=item * C<< parser => parser_subclass [gff3] >>

This optional parameter defines which parser subclass to instantiate.
This parameter will default to parsing GFF3 if not provided.  See L<GAL::Parser>
for a complete list of available parser classes.

=item * C<< storage => storage_subclass [SQLite] >>

This optional parameter defines which storage subclass to instantiate.
Currently available storage classes are SQLite (the default) and
mysql.

=item * C<< fasta => '/path/to/fasta/files/ >>

This optional parameter defines a path to a fasta file or a collection
of fasta files that correspond the annotated features.  The IDs (first
contiguous non-whitespace charachters) of the fasta headers must
correspond to the sequence IDs (seqids) in the annotated features.
The fasta parameter is optional, but if the fasta attribute is not set
then the features will not have access to their sequence.  Access to
the sequence in provided by Bio::DB::Fasta.

=back

=cut

#-----------------------------------------------------------------------------
#-------------------------------- Constructor --------------------------------
#-----------------------------------------------------------------------------

=head2 new

     Title   : new
     Usage   : GAL::Annotation->new($gff3_file, $fasta_file);
     Function: Creates a GAL::Annotation object;
     Returns : A GAL::Annotation object
     Args    : A list (ref) of key value pairs for the attributes specified above.

=cut

sub new {
    my ($class, @args) = @_;

    my ($args, $feature_file) = _check_for_simple_args(@args);

    my $self = $class->SUPER::new($args);
    if (! defined $feature_file) {
      $self->throw('missing_feature_file',
		   'GAL::Annotation requires at least a feature file');
    }
    $self->load_files($feature_file) if $feature_file;
    return $self;
}

#-----------------------------------------------------------------------------

sub _initialize_args {
    my ($self, @args) = @_;

    ######################################################################
    # This block of code handles class attributes.  Use the
    # @valid_attributes below to define the valid attributes for
    # this class.  You must have identically named get/set methods
    # for each attribute.  Leave the rest of this block alone!
    ######################################################################
    my $args = $self->SUPER::_initialize_args(@args);
    # Set valid class attributes here
    my @valid_attributes = qw(fasta parser storage);
    $self->set_attributes($args, @valid_attributes);
    ######################################################################
}

#-----------------------------------------------------------------------------

sub _check_for_simple_args {

    my @args = @_;

    return {} unless @args;

    # If one or two files are passed in we'll assume they are a GFF3
    # and Fasta file.
    my ($args, $feature_file, $fasta_file);
    # If one or two args
    if (@args <= 2) {
	if (defined $args[0] && # and args are defined and not
	    ! ref $args[0]   && # and not a reference
	    $args[0] !~ /^(parser|storage|fasta)$/) { # and not an option
	    ($feature_file, $fasta_file) = @args;
	    my $sqlite_file;
	    ($sqlite_file = $feature_file) =~ s/\.[^\.]+$//;
	    $sqlite_file .= '.sqlite';
	    $args = ({parser  => {class => 'gff3'},
		      storage => {class    => 'SQLite',
				  database => $sqlite_file}
		  });
	    $args->{file}  = $feature_file;
	    $args->{fasta} = $fasta_file if $fasta_file;
	}
	else {
	    my %args;
	    %args = @args if defined $args[1];
	    $args = \%args;
	}
	return ($args, $feature_file);
    }
    return {};
}

#-----------------------------------------------------------------------------
#-------------------------------- Attributes ---------------------------------
#-----------------------------------------------------------------------------

=head1  ATTRIBUTES

All attributes can be supplied as parameters to the GAL::Annotation
constructor as a list (or referenece) of key value pairs.

=cut

=head2 parser

 Title   : parser
 Usage   : $parser = $self->parser(class => gff3);
 Function: Create or return a parser object.
 Returns : A GAL::Parser::subclass object.
 Args    : (class => gal_parser_subclass_name)
	   See L<GAL::Parser> and its subclasses for more arguments.
 Notes   : The parser object is created as a singleton, but it
	   can be changed by passing new arguments to a call to
	   parser.

=cut

sub parser {
  my ($self, @args) = @_;

  if (! $self->{parser} || @args) {
    my $class = 'GAL::Parser';
    $self->load_module($class);
    my $parser = $class->new(@args);
    my $weak_self = $self;
    weaken $weak_self;
    $parser->annotation($weak_self);
    $self->{parser} = $parser;
  }
  return $self->{parser};
}

#-----------------------------------------------------------------------------

=head2 storage

  Title   : storage
  Usage   : $storage = $self->storage(class => 'SQLite');
  Function: Create or return a storage object.
  Returns : A GAL::Storage::subclass object.
  Args    : (class => gal_storage_subclass)
	    See GAL::Storage and its subclasses for more arguments.
  Notes   : The storage object is created as a singleton and can not be
	    destroyed or recreated after being created.

=cut

sub storage {
  my ($self, @args) = @_;

  if (! $self->{storage}) {
    my $class = 'GAL::Storage';
    $self->load_module($class);
    my $storage = $class->new(@args);
    my $weak_self = $self;
    weaken $weak_self;
    $storage->annotation($weak_self);
    $self->{storage} = $storage;
  }
  return $self->{storage};
}

#-----------------------------------------------------------------------------

=head2 fasta

The fasta attribute is provided by L<GAL::Base>, see that module for more
details.

=cut

#-----------------------------------------------------------------------------
#---------------------------------- Metohds ----------------------------------
#-----------------------------------------------------------------------------

=head1 Methods

=head2 features

  Title   : features
  Usage   : $self->features();
  Function: Return an iterator L<GAL::Schema::Result::Feature> object.
            This object is a subclass of L<DBIx::Class::ResultSet> for
            all features.
  Returns : A L<GAL::Schema::Result::Feature> object
  Args    : N/A

=cut

sub features {
  my ($self, $query) = shift;
  my $features = $self->schema->resultset('Feature');
  return $features;
}

#-----------------------------------------------------------------------------

=head2 schema

  Title   : schema
  Usage   : $self->schema();
  Function: Create and/or return the L<DBIx::Class::Schema> object
  Returns : L<DBIx::Class::Schema> object.
  Args    : N/A - Arguments are provided by the <GAL::Storage> object.

=cut

sub schema {
    my $self = shift;

    # Create the schema if it doesn't exist
    if (! $self->{schema}) {
	# We should be able to reuse the Annotation dbh here!
	# my $schema = GAL::Schema->connect([sub {$self->storage->dbh}]);
	my $schema = GAL::Schema->connect($self->storage->dsn,
					  $self->storage->user,
					  $self->storage->password,
					  );
	# If we're using SQLite we should be using a larger cache_size;
	# but why doesn't this seem to help?  Do this in the Storage object
	# $schema->storage->dbh_do(sub {
	#			     my ($storage, $dbh) = @_;
	#			     $dbh->do('PRAGMA cache_size = 800000');
	#			   },
	#			  );
	my $weak_self = $self;
	weaken $weak_self;
	$schema->annotation($weak_self); # See GAL::SchemaAnnotation
	$self->{schema} = $schema;
    }
    return $self->{schema};
}

#-----------------------------------------------------------------------------

=head2 load_files

 Title   : load_files
 Usage   : $a = $self->load_files(qw(file1.gff3 file2.gff3));
 Function: Parse and store all of the features in a file. If a single
	   file is given as an argument and if there are gff[3] and
	   sqlite versions of that files base name then time stamps
	   are compared and the database is only (re)loaded if the
	   GFF3 file is newer.
 Returns : N/A
 Args    : A list of files.
 Notes   : Default

=cut

sub load_files {

  my ($self, $files) = @_;

  # If a single file
  if (! ref $files) {
    my ($file_ext) = $files =~ /\.([^\.]+)$/;
    my $file_base;
    ($file_base = $files) =~ s/\.[^\.]+$//;
    my $sqlite_file = $file_base . '.sqlite';
    # If an SQLite versions of the file exists
    if (-e $files && -e $sqlite_file) {
      # If the SQLite file is newer
      if ((stat($sqlite_file))[9] > (stat($files))[9]) {
	# Then load that SQLite file
	$self->storage->database($sqlite_file);
	return;
      }
    }
  }
  $self->storage->load_files($files);
}

#-----------------------------------------------------------------------------

=head1 DIAGNOSTICS

=over

=item C<< missing_feature_file >>

C<GAL::Annotation> must have at least a feature file to run, but you
did not provide one.

missing_feature_file

=head1 CONFIGURATION AND ENVIRONMENT

<GAL::Annotation> requires no configuration files or environment variables.

=head1 DEPENDENCIES

Modules in GAL/lib use the following modules:

Bio::DB::Fasta
Carp
DBD::SQLite
DBI
List::Util
Scalar::Util
Set::IntSpan::Fast
Statistics::Descriptive
Text::RecordParser

Some scripts in GAL/bin and/or GAL/lib/GAL/t use the following modules:

Data::Dumper
FileHandle
Getopt::Long
IO::Prompt
List::MoreUtils
TAP::Harness
Test::More
Test::Pod::Coverage
URI::Escape
XML::LibXML::Reader

=head1 INCOMPATIBILITIES

None reported.

=head1 BUGS AND LIMITATIONS

I'm sure there are plenty of bugs right now - please let me know if you find one.

Please report any bugs or feature requests to:
barry.moore@genetics.utah.edu

=head1 AUTHOR

Barry Moore <barry.moore@genetics.utah.edu>

=head1 LICENCE AND COPYRIGHT

Copyright (c) 2010-2014, Barry Moore <barry.moore@genetics.utah.edu>.  All rights reserved.

    This module is free software; you can redistribute it and/or
    modify it under the same terms as Perl itself (See LICENSE).

=head1 DISCLAIMER OF WARRANTY

BECAUSE THIS SOFTWARE IS LICENSED FREE OF CHARGE, THERE IS NO WARRANTY
FOR THE SOFTWARE, TO THE EXTENT PERMITTED BY APPLICABLE LAW. EXCEPT WHEN
OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS AND/OR OTHER PARTIES
PROVIDE THE SOFTWARE "AS IS" WITHOUT WARRANTY OF ANY KIND, EITHER
EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE
ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE SOFTWARE IS WITH
YOU. SHOULD THE SOFTWARE PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL
NECESSARY SERVICING, REPAIR, OR CORRECTION.

IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING
WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MAY MODIFY AND/OR
REDISTRIBUTE THE SOFTWARE AS PERMITTED BY THE ABOVE LICENCE, BE
LIABLE TO YOU FOR DAMAGES, INCLUDING ANY GENERAL, SPECIAL, INCIDENTAL,
OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OR INABILITY TO USE
THE SOFTWARE (INCLUDING BUT NOT LIMITED TO LOSS OF DATA OR DATA BEING
RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR THIRD PARTIES OR A
FAILURE OF THE SOFTWARE TO OPERATE WITH ANY OTHER SOFTWARE), EVEN IF
SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE POSSIBILITY OF
SUCH DAMAGES.

=cut

1;
