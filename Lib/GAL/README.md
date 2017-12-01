GAL::Annotation - Genome Annotation Library

![GAL](https://github.com/The-Sequence-Ontology/GAL/blob/master/gal-image.jpg)

#VERSION

This document describes GAL::Annotation version 0.01

#SYNOPSIS

```perl 
use GAL::Annotation;

# Assuming defaults (GFF3 parser and SQLite storage)
my $annot = GAL::Annotation->new(qw(file.gff file.fasta);
my $features = $annot->features;

# Either way, once you have features - get to work.
my $mrnas = $features->search( {type => 'mRNA'} );
while (my $mrna = $mrnas->next) {
    print $mrna->feature_id . "\n";
    # Introns don't exist in the GFF3 file so infer them
    my $intons = $mrna->introns;
    while (my $intron = $introns->next) {
        print join "\n", (
            $intron->start,
            $intron->end,
            $intron->seq,
        );
    }
}    

```

#DESCRIPTION

The Genome Annotation Library (GAL) is a collection of modules that
strive to make working with genome annotations simple, intuitive and
fast. Users of GAL first create an annotation object which in turn
will contain Parser, Storage and Schema objects. The parser allows
features to be loaded into GAL's storage from a variety of
formats. The storage object specifies how the features should be
stored, and the schema object provides flexible query and iterator
functions over the features.  In addtion, Index objects (not yet
implimented) provide additional key/value mapped look up tables, and
List objects provide aggregation and analysis functionality for lists
of feature attributes.

A wide variety of parsers are available to convert sequence features
from various formats, and new parsers are easy to write. See
GAL::Parser for more details. Currently SQLite and MySQL storage
options are available (a fast RAM storage engine is on the TODO
list). Schema objects are provided by DBIx::Class and a familiarity
with that package is necessary to fully understand how to query and
iterate over feature objects.

#CONSTRUCTOR

New Annotation objects are created by the class method `new`. Arguments
should be passed to the constructor as a list (or reference) of key
value pairs. All attributes of the Annotation object can be set in the
call to new, but reasonable defaults will be used where ever possilbe
to keep object creation simple. An simple example of object creation
would look like this:

```perl
my $feat_store = GAL::Annotation->new($gff_file);
my $feat_store = GAL::Annotation->new($gff_file, $fasta_file);
```

The resulting object would use a GFF3 parser and SQLite storage by
default The first example would not have access to feature sequence,
the second one would.

A more complex object creation might look like this:

```perl
my $feat_store = GAL::Annotation->new(
    parser  => { class => gff3 },
    storage => {
        class    => mysql,
        dsn      => 'dbi:mysql:database' user => 'me',
        password => 'secret' fasta => '/path/to/fasta/files/'
    }
);
```

The constructor recognizes the following parameters which will set the
appropriate attributes:

* `parser => parser_subclass [gff3]`

    This optional parameter defines which parser subclass to
    instantiate. This parameter will default to gff3 if not provided.
    See GAL::Parser for a complete list of available parser classes.

* `storage => storage_subclass [SQLite]`

    This optional parameter defines which storage subclass to
    instantiate. Currently available storage classes are SQLite (the
    default) and mysql.

* `fasta => '/path/to/fasta/files/`

    This optional parameter defines a path to a fasta file or a
    collection of fasta files that correspond the annotated features.
    The IDs (first contiguous non-whitespace charachters) of the fasta
    headers must correspond to the sequence IDs (seqids) in the
    annotated features. The fasta parameter is optional, but if the
    fasta attribute is not set then the features will not have access
    to their sequence. Access to the sequence in provided by
    Bio::DB::Fasta.

### new ###

       Title   : new
       Usage   : GAL::Annotation->new();
       Function: Creates a GAL::Annotation object;
       Returns : A GAL::Annotation object
       Args    : A list of key value pairs for the attributes specified above.

# ATTRIBUTES #

All attributes can be supplied as parameters to the GAL::Annotation
constructor as a list (or referenece) of key value pairs.

### parser ###

     Title   : parser
     Usage   : $parser = $self->parser();
     Function: Create or return a parser object.
     Returns : A GAL::Parser::subclass object.
     Args    : (class => gal_parser_subclass)
               See GAL::Parser and its subclasses for more arguments.
     Notes   : The parser object is created as a singleton, but it
               can be changed by passing new arguments to a call to
               parser.

### storage ###

    Title   : storage
    Usage   : $storage = $self->storage();
    Function: Create or return a storage object.
    Returns : A GAL::Storage::subclass object.
    Args    : (class => gal_storage_subclass)
              See GAL::Storage and its subclasses for more arguments.
    Notes   : The storage object is created as a singleton and can not be
              destroyed or recreated after being created.

### fasta ###

    The fasta attribute is provided by GAL::Base, see that module for
    more details.

# METHODS #

### features ###

    Title   : features
    Usage   : $self->features();
    Function: Return a GAL::Schema::Result::Feature object (a
              DBIx::Class::ResultSet for all features).
    Returns : A GAL::Schema::Result::Feature object
    Args    : N/A

### schema ###

    Title   : schema
    Usage   : $self->schema();
    Function: Create and/or return the DBIx::Class::Schema object
    Returns : DBIx::Class::Schema object.
    Args    : N/A - Arguments are provided by the GAL::Storage object.

### load_files ###

    Title   : load_files
    Usage   : $a = $self->load_files();
    Function: Parse and store all of the features in a file. If a single
              file is given as an argument and if there are gff[3] and
              sqlite versions of that files base name then time stamps
              are compared and the database is only (re)loaded if the
              GFF3 file is newer.
    Returns : N/A
    Args    : A list of files.
    Notes   : Default

# DIAGNOSTICS #

<GAL::Annotation> currently does not throw any warnings or errors,
but most other modules in the library do, and details of those
errors can be found in those modules.

# CONFIGURATION AND ENVIRONMENT #

<GAL::Annotation> requires no configuration files or environment
variables.

# DEPENDENCIES #

Modules in GAL/lib use the following modules:

* Bio::DB::Fasta 
* Carp 
* DBD::SQLite 
* DBI 
* List::Util 
* Scalar::Util
* Set::IntSpan::Fast Statistics::Descriptive 
* Text::RecordParser

Some script in GAL/bin and/or GAL/lib/GAL/t use the following
modules:

* Data::Dumper 
* FileHandle Getopt::Long 
* IO::Prompt 
* List::MoreUtils
* TAP::Harness 
* Test::More 
* Test::Pod::Coverage 
* URI::Escape
* XML::LibXML::Reader

# INCOMPATIBILITIES #

None reported.

# BUGS AND LIMITATIONS #

I'm sure there are plenty of bugs right now - please let me know if
you find one.

Please report any bugs or feature requests to:
<barry.moore@genetics.utah.edu>

# AUTHOR #

Barry Moore <barry.moore@genetics.utah.edu>

# LICENCE AND COPYRIGHT #

Copyright (c) 2010-2015, Barry Moore <barry.moore@genetics.utah.edu>. All
rights reserved.

This module is free software; you can redistribute it and/or
modify it under the same terms as Perl itself.

# DISCLAIMER OF WARRANTY #

BECAUSE THIS SOFTWARE IS LICENSED FREE OF CHARGE, THERE IS NO WARRANTY
FOR THE SOFTWARE, TO THE EXTENT PERMITTED BY APPLICABLE LAW. EXCEPT
WHEN OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS AND/OR OTHER
PARTIES PROVIDE THE SOFTWARE "AS IS" WITHOUT WARRANTY OF ANY KIND,
EITHER EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE. THE ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE
SOFTWARE IS WITH YOU. SHOULD THE SOFTWARE PROVE DEFECTIVE, YOU ASSUME
THE COST OF ALL NECESSARY SERVICING, REPAIR, OR CORRECTION.

IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING
WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MAY MODIFY AND/OR
REDISTRIBUTE THE SOFTWARE AS PERMITTED BY THE ABOVE LICENCE, BE LIABLE
TO YOU FOR DAMAGES, INCLUDING ANY GENERAL, SPECIAL, INCIDENTAL, OR
CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OR INABILITY TO USE THE
SOFTWARE (INCLUDING BUT NOT LIMITED TO LOSS OF DATA OR DATA BEING
RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR THIRD PARTIES OR A
FAILURE OF THE SOFTWARE TO OPERATE WITH ANY OTHER SOFTWARE), EVEN IF
SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH
DAMAGES.

