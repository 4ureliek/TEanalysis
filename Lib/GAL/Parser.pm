package GAL::Parser;

use strict;
use vars qw($VERSION);


$VERSION = 0.2.0;
use base qw(GAL::Base);

=head1 NAME

GAL::Parser - Parser objects for the Genome Annotation Library

=head1 VERSION

This document describes GAL::Parser version 0.2.0

=head1 SYNOPSIS

    use GAL::Parser;

    my $parser = GAL::Parser::gff3->new(file => 'data/feature.gff3');

    while (my $feature_hash = $parser->next_feature_hash) {
	print $parser->to_gff3($feature_hash) . "\n";
    }


=head1 DESCRIPTION

<GAL::Parser> is not intended to be instantiated by itself, but rather
functions as a base class for <GAL::Parser> subclasses.  It provides a
variety of attributes and methods that are generally applicable to all
parsers.  While parsers are intended for use from within
<GAL::Annotation> objects. Parsers can be instantiated seperately from
the rest of the GAL library and there are many use cases for the
parsers as stand alone objects.  Anytime you just need fast access to
iterate over all features in a flat file and are happy to have hashes
of those features you should just use the parser directly without the
<GAL::Annotation> object.

=head1 CONSTRUCTOR

New GAL::Parser::subclass objects are created by the class method
new.  Arguments should be passed to the constructor as a list (or
reference) of key value pairs.  All attributes of the Parser object
can be set in the call to new. An simple example
of object creation would look like this:

    my $parser = GAL::Parser::gff3->new(file => 'data/feature.gff3');

The constructor recognizes the following parameters which will set the
appropriate attributes:

=over 4

=item * C<< file => feature_file.txt >>

This optional parameter provides the filename for the file containing
the data to be parsed. While this parameter is optional either it, or
the following fh parameter must be set.

=item * C<< fh => feature_file.txt >>

This optional parameter provides a filehandle to read data from. While
this parameter is optional either it, or the following fh parameter
must be set.

=item * C<< annotation => $gal_annotation_object >>

This parameter is not intended for public use as a setter, but it is
available for use as a getter.  When a parser is instantiated via
<GAL::Annotation>, a weakened copy of the <GAL::Annotation> object is
stored in the parser.

=back

=cut

#-----------------------------------------------------------------------------
#-------------------------------- Constructor --------------------------------
#-----------------------------------------------------------------------------

=head2 new

     Title   : new
     Usage   : GAL::Parser::subclass->new();
     Function: Creates a GAL::Parser object;
     Returns : A GAL::Parser object
     Args    : (file => $file)
	       (fh   => $FH)

=cut

sub new {
	my ($class, @args) = @_;
	my $self = $class->SUPER::new(@args);
	return $self;
}

#-----------------------------------------------------------------------------

sub _initialize_args {
	my ($self, @args) = @_;

	######################################################################
	# This block of code handels class attributes.  Use the
	# @valid_attributes below to define the valid attributes for
	# this class.  You must have identically named get/set methods
	# for each attribute.  Leave the rest of this block alone!
	######################################################################
	my $args = $self->SUPER::_initialize_args(@args);
	# Set valid class attributes here
	my @valid_attributes = qw(annotation class fh file);
	$self->set_attributes($args, @valid_attributes);
	######################################################################
	return $args;
}

#-----------------------------------------------------------------------------
#-------------------------------- Attributes ---------------------------------
#-----------------------------------------------------------------------------

=head1  ATTRIBUTES

All attributes can be supplied as parameters to the constructor as a
list (or referenece) of key value pairs.

=head2 annotation

 Title   : annotation
 Usage   : $a = $self->annotation()
 Function: Get/Set the value of annotation.
 Returns : The value of annotation.
 Args    : A value to set annotation to.

=cut

sub annotation {
  my ($self, $annotation) = @_;
  $self->{annotation} = $annotation if $annotation;
  return $self->{annotation};
}

# #-----------------------------------------------------------------------------
#
# =head2 class
#
#  Title   : class
#  Usage   : $a = $self->class()
#  Function: Get/Set the value of class.
#  Returns : The value of class.
#  Args    : A value to set class to.
#
# =cut
#
# sub class {
#   my ($self, $class) = @_;
#
#   if ($class) {
#       $class =~ s/GAL::Parser:://;
#       $class = 'GAL::Parser::' . $class;
#       $self->{class} = $class;
#   }
#   $self->{class} ||= 'GAL::Parser::gff3';
#   return $self->{class};
# }
#
# #-----------------------------------------------------------------------------

=head2 file

 Title   : file
 Usage   : $a = $self->file();
 Function: Get/Set the value of file.
 Returns : The value of file.
 Args    : A value to set file to.

=cut

sub file {
	my ($self, $file) = @_;
	return $self->reader->file($file);
}

#-----------------------------------------------------------------------------

=head2 fh

 Title   : fh
 Usage   : $a = $self->fh();
 Function: Get/Set the filehandle.
 Returns : The filehandle.
 Args    : A filehandle.

=cut

sub fh {
  my ($self, $fh) = @_;
  return $self->reader->fh($fh);
}

#-----------------------------------------------------------------------------

=head2 reader

 Title   : reader
 Usage   : $a = $self->reader();
 Function: Get/Set the reader.
 Returns : The value of reader.
 Args    : A value to set reader to.

=cut

sub reader {
  my $self = shift;

  if (! $self->{reader}) {
    my @field_names = qw(seqid source type start end score strand phase
			 attributes);
    $self->load_module('GAL::Reader::DelimitedLine');
    my $reader = GAL::Reader::DelimitedLine->new(field_names => \@field_names);
    $self->{reader} = $reader;
  }
  return $self->{reader};
}

#-----------------------------------------------------------------------------
#---------------------------------- Methods ----------------------------------
#-----------------------------------------------------------------------------

=head1 METHODS

=head2 next_record

 Title   : next_record
 Usage   : $a = $self->next_record();
 Function: Return the next record from the reader in whatever format
	   that reader specifies.
 Returns : The next record from the reader.
 Args    : N/A

=cut

sub next_record {
	shift->reader->next_record;
}

#-----------------------------------------------------------------------------

=head2 metadata

 Title   : metadata
 Usage   : $a = $self->metadata
 Function: Return the metadata as parsed by the parse_metadata
 Returns : A hash (reference)
 Args    : A hash reference of metadata

=cut

sub metadata {
    my ($self, $metadata) = @_;

    if ($metadata) {
	$self->{metadata} = $metadata;
    }
    if (! $self->{metadata}) {
	$self->parse_metadata;
    }
    return $self->{metadata};
}

#-----------------------------------------------------------------------------

=head2 counter

 Title   : counter
 Usage   : $a = $self->counter('name', $start);
	   $a = $self->counter('name');
 Function: Return a counter associated with the given name.  A counter with
	   that name will only created once during the life of the
	   object and will return incrementing values each time called
	   thereafter.  A start value may optionally be given.  This
	   can be given each time, but will only be used the first
	   time given and discarded thereafter.
 Returns : An incrementing value for the counter
 Args    : A name for the counter and optionally a integer each for a start and
	   step value.

=cut

sub counter {
  my ($self, $name, $start) = @_;

  if (exists $self->{counters}{$name}) {
    return $self->{counters}{$name}++;
  }

  $start ||= 1;

  $self->{counters}{$name} = $start + 1;
  return $start;
}

#-----------------------------------------------------------------------------

=head2 next_feature_hash

 Title   : next_feature_hash
 Usage   : $a = $self->next_feature_hash;
 Function: Return the next record from the parser as a 'feature hash'.
 Returns : A hash or hash reference.
 Args    : N/A

The feature hash has the following format:

%feature = (feature_id => $feature_id,
	    seqid      => $seqid,
	    source     => $source,
	    type       => $type,
	    start      => $start,
	    end        => $end,
	    score      => $score,
	    strand     => $strand,
	    phase      => $phase,
	    attributes => {ID => $feature_id,
			   Parent => $parent_id,
			  }
	   );

This hash follows the format layed out by the GFF3 specification
(http://www.sequenceontology.org/resources/gff3.html).  Please see
that document for details on constrants for each value and attribute.

=cut

sub next_feature_hash {
  my $self = shift;

  my $feature;

  # If a previous record has returned multiple features then
  # grab them off the stack first instead of reading a new one
  # from the file.
  if (ref $self->{_feature_stack} eq 'ARRAY' &&
      scalar @{$self->{_feature_stack}} > 0) {
    $feature = shift @{$self->{_feature_stack}};
    return wantarray ? %{$feature} : $feature;
  }

  # Allow parse_record to return undef to ignore a record, but
  # still keep parsing the file.
  until ($feature) {
    # Get the next record from the file.
    my $record = $self->next_record;
    return undef if ! defined $record;

    # Parser the record - probably overridden by a subclass.
    $feature = $self->parse_record($record);

    # Allow parsers to return more than one feature.
    # This allows the parser to expand a single record into
    # multiple features.
    if (ref $feature eq 'ARRAY') {
      my $this_feature = shift @{$feature};
      push @{$self->{_feature_stack}}, @{$feature};
      $feature = $this_feature;
    }
  }

  $self->validate_feature($feature);

  return wantarray ? %{$feature} : $feature;
}

#-----------------------------------------------------------------------------

=head2  all_feature_hashes

 Title   : all_feature_hashes
 Alias   : all
 Usage   : $self->all_feature_hashes();
 Function: Returns an array(ref) of all (remaining) features (as hashes) 
           available via the current filehandle.
 Returns : An array (or array reference) 
 Args    : N/A

=cut

sub all {shift->all_feature_hashes()};

sub all_feature_hashes {

    my $self = shift;

    my @features;
    while (my $feature = $self->next_feature_hash) {
	push @features, $feature if $feature;
    }
    return wantarray ? @features : \@features;
}

#-----------------------------------------------------------------------------

=head2  to_gff3

 Title   : to_gff3
 Usage   : $self->to_gff3($feature_hash)
 Function: Returns a string of GFF3 formatted text for a given feature hash
 Returns : A string in GFF3 format.
 Args    : A feature hash reference in the form returned by next_feature_hash

=cut

sub to_gff3 {
  my ($self, $feature) = @_;

  my %ATTRB_ORDER = (
		     ID                   =>  1,
		     Name                 =>  2,
		     Alias                =>  3,
		     Parent               =>  3,
		     Target               =>  4,
		     Gap                  =>  5,
		     Derives_from         =>  6,
		     Note                 =>  7,
		     Dbxref               =>  8,
		     Ontology_term        =>  9,
		     Variant_seq          =>  10,
		     Reference_seq        =>  11,
		     Variant_reads        =>  12,
		     Total_reads          =>  13,
		     Genotype             =>  14,
		     Variant_effect       =>  15,
		     Variant_copy_number  =>  16,
		    );

  my $attribute_text;
  for my $key (sort {($ATTRB_ORDER{$a} || 99) <=> ($ATTRB_ORDER{$b} || 99) ||
		       $a cmp $b}
	       keys %{$feature->{attributes}}) {
    my $value_text = join ',', @{$feature->{attributes}{$key}};
    $attribute_text .= "$key=$value_text;";
  }

  my $gff3_text = join "\t", ($feature->{seqid},
			      $feature->{source},
			      $feature->{type},
			      $feature->{start},
			      $feature->{end},
			      $feature->{score},
			      $feature->{strand},
			      $feature->{phase},
			      $attribute_text,
			     );

  return $gff3_text;
}

#-----------------------------------------------------------------------------

=head2 parse_metadata

 Title   : parse_metadata
 Usage   : $a = $self->parse_metadata
 Function: Call the reader to return metadata and then parse and store it.
 Returns : N/A
 Args    : A hash reference of parse_metadata

=cut

sub parse_metadata {
  my $self = shift;

  my $raw_data = $self->reader->metadata;
  # Subclasses can override and do something here.
  $self->{metadata} = $raw_data;
}

#-----------------------------------------------------------------------------

=head2 parse_record

 Title   : parse_record
 Usage   : $a = $self->parse_record();
 Function: Parse the data from a record.
 Returns : Feature data as a hash (or reference);
 Args    : A data structure of feature data that this method (probably
	   overridden by a subclass) understands.

=cut

sub parse_record {
  my ($self, $record) = @_;

  my $attributes = $self->parse_attributes($record->{attributes});

  my %feature = (seqid      => $record->{seqid},
		 source     => $record->{source},
		 type       => $record->{type},
		 start      => $record->{start},
		 end        => $record->{end},
		 score      => $record->{score},
		 strand     => $record->{strand},
		 phase      => $record->{phase},
		 attributes => $attributes,
		);

  if (exists $attributes->{ID}         &&
      ref $attributes->{ID} eq 'ARRAY' &&
      $attributes->{ID}[0]) {
    $feature{feature_id} = $attributes->{ID}[0];
  }
  else {
    $feature{feature_id} =  $self->create_unique_id(\%feature);
    $attributes->{ID}[0] = $feature{feature_id}
  }

  return wantarray ? %feature : \%feature;
}

#-----------------------------------------------------------------------------

=head2 parse_attributes

 Title   : parse_attributes
 Usage   : $a = $self->parse_attributes($attrb_text);
 Function: Parse the attributes from a GFF3 column 9 formatted string of text.
 Returns : A hash (or reference) of attribute key value pairs.
 Args    : A GFF3 column 9 formated string of text.

=cut

sub parse_attributes {
	my ($self, $attrb_text) = @_;

	my @attrbs = split /\s*;\s*/, $attrb_text;
	my %attrb_hash;
	for my $attrb (@attrbs) {
		my ($tag, $value_text) = split /=/, $attrb;
		my @values = split /,/, $value_text;
		push @{$attrb_hash{$tag}}, @values;
	}
	return wantarray ? %attrb_hash : \%attrb_hash;
}

#-----------------------------------------------------------------------------

=head2 create_unique_id

 Title   : create_unique_id
 Usage   : $feature_id = $self->create_unique_id(\%feature_hash, 8);
 Function: Creates a uniq_id for the feature of the form gene_0000001
	   where the prefix (gene) is the feature type and the numerical
	   sufix is an incrementing count padded with zeros to the length
	   specified or 8 by default.
 Returns : A unique ID
 Args    : A feature hash and an integer (12 max) for padding the
	   numerical portion of the ID.

=cut

sub create_unique_id {
  my ($self, $feature, $pad) = @_;

  $pad ||= 8;
  my $type = $feature->{type};

  my $id = "${type}_" . sprintf("%0${pad}s", ++$self->{id_counter}{$type});

  return $id;
}

#-----------------------------------------------------------------------------

=head2 validate_feature

 Title   : validate_feature
 Usage   : $self->validate_feature(%feature_hash);
 Function: Provides very basic syntactic validation for features.
 Returns : 1 on success - dies on errors
 Args    : A feature hash as a hash reference.

=cut

sub validate_feature {

   my ($self, $feature) = @_;

   $self->throw('missing_feature_id', $self->to_gff3($feature))
     unless exists $feature->{feature_id} && $feature->{feature_id};;
   $self->throw('missing_seqid', $self->to_gff3($feature))
     unless exists $feature->{seqid} || $feature->{seqid};
   $self->throw('missing_source', $self->to_gff3($feature))
     unless exists $feature->{source} && $feature->{source};
   $self->throw('missing_type', $self->to_gff3($feature))
     unless exists $feature->{type} && $feature->{type};
   $self->throw('missing_start', $self->to_gff3($feature))
     unless exists $feature->{start};
   $self->throw('invalid_start', $self->to_gff3($feature))
     unless $feature->{start} =~ /^\d+$/;
   $self->throw('missing_end', $self->to_gff3($feature))
     unless exists $feature->{end};
   $self->throw('invalid_end', $self->to_gff3($feature))
     unless $feature->{end} =~ /^\d+$/;
   $self->throw('missing_score', $self->to_gff3($feature))
     unless exists $feature->{score} && defined $feature->{score};
   $self->throw('missing_strand', $self->to_gff3($feature))
     unless exists $feature->{strand};
   $self->throw('invalid_strand', $self->to_gff3($feature))
     unless $feature->{strand} =~ /^[\.\+\-]$/;
   $self->throw('missing_phase', $self->to_gff3($feature))
     unless exists $feature->{phase} && defined $feature->{phase};
   $self->throw('invalid_phase', $self->to_gff3($feature))
     unless $feature->{phase} =~ /^[\.012]/;
   $self->throw('missing_attributes', $self->to_gff3($feature))
     unless exists $feature->{attributes} && $feature->{attributes};
   $self->throw('missing_attribute_ID', $self->to_gff3($feature))
     unless exists $feature->{attributes}{ID};
   $self->throw('invalid_attribute_ID', $self->to_gff3($feature))
     unless ref $feature->{attributes}{ID} eq 'ARRAY';
   $self->throw('missing_attribute_ID', $self->to_gff3($feature))
     unless $feature->{attributes}{ID}[0];
   return 1;
}

#-----------------------------------------------------------------------------

=head1 DIAGNOSTICS

<GAL::Parser> currently does not throw any warnings or errors, but
subclasses may, and details of those errors can be found in those
modules.

=head1 CONFIGURATION AND ENVIRONMENT

<GAL::Parser> requires no configuration files or environment variables.

=head1 DEPENDENCIES

<GAL::Reader>

<GAL::Reader> and subclasses of <GAL::Parser> have other dependencies.

=head1 INCOMPATIBILITIES

None reported.

=head1 BUGS AND LIMITATIONS

No bugs have been reported.

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
