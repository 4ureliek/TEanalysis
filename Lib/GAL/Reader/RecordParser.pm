package GAL::Reader::RecordParser;

use strict;
use vars qw($VERSION);


$VERSION = 0.2.0;
use base qw(GAL::Reader);
use Text::RecordParser;

=head1 NAME

GAL::Reader::RecordParser - Record Parsing using Text::RecordParser

=head1 VERSION

This document describes GAL::Reader::RecordParser version 0.2.0

=head1 SYNOPSIS

    use GAL::Reader::RecrodParser
    $reader = GAL::Reader::RecrodParser->new(file => 'annotation_file.txt',
					     record_separator => "\n",
					     field_separator  => "\t",
					     bind_fields => [qw(seqid source
								type start end
								score strand
								phase attrb
								)
							    ],
					    );
    $reader->next_record, '$reader->next_record');

=head1 DESCRIPTION

<GAL::Reader::RecordParser> provides flexible record reading via
Text::RecordParser.  It is not intended for general library use, but
rather as a GAL::Reader subclass for developers of GAL::Parser
subclasses.  There is however no reason why it couldn't also be used
as a stand alone module for other purposes.

=head1 CONSTRUCTOR

New GAL::Reader::RecordParser objects are created by the class method new.
Arguments should be passed to the constructor as a list (or reference)
of key value pairs.  All attributes of the Reader object can be set in
the call to new. An simple example of object creation would look like
this:

    $reader = GAL::Reader::RecrodParser->new(file => 'annotation_file.txt',
					     record_separator => "\n",
					     field_separator  => "\t",
					     bind_fields => [qw(seqid source
								type start end
								score strand
								phase attrb
								)
							    ],
					    );

The constructor recognizes the following parameters which will set the
appropriate attributes:

=over 4

=item * C<< record_separator >>

This optional parameter defines the pattern by which records are
separated.  The default is a new line.

=item * C<< field_separator >>

This optional parameter defines the pattern by which fields are
separated.  The default is a comma.

=item * C<< comment >>

This optional parameter defines the pattern by which comment lines
are identified.

=item * C<< bind_fields => [qw(seqid source type start end)] >>

This attribute provides an orderd list that describes the field names
of the columns in the delimited file.  If this attribute is not set then
Text::RecordParser will automatically assume that the first line of text
in the file are headers.

=cut

#-----------------------------------------------------------------------------

=head2 new

     Title   : new
     Usage   : GAL::Reader::RecordParser->new();
     Function: Creates a GAL::Reader::RecordParser object;
     Returns : A GAL::Reader::RecordParser object
     Args    :

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
	my @valid_attributes = qw(file fh record_separator field_separator
				  comment_delimiter fields);
	$self->set_attributes($args, @valid_attributes);
	######################################################################
	return $args;
}

#-----------------------------------------------------------------------------

=head2 record_separator

 Title   : record_separator
 Usage   : $a = $self->record_separator("\n");
 Function: Gets/set the pattern to use as the record separator.
 Returns : The pattern to use as the record separator.
 Args    : The pattern to use as the record separator.

=cut

sub record_separator {
  my ($self, $record_separator) = @_;
  $self->{record_separator} = $record_separator if $record_separator;
  return $self->{record_separator};
}

#-----------------------------------------------------------------------------

=head2 field_separator

 Title   : field_separator
 Usage   : $a = $self->field_separator("\t");
 Function: Gets/set the pattern to use as the field separator.
 Returns : The pattern to use as the field separator.
 Args    : The pattern to use as the field separator.

=cut

sub field_separator {
  my ($self, $field_separator) = @_;
  $self->{field_separator} = $field_separator if $field_separator;
  return $self->{field_separator};
}

#-----------------------------------------------------------------------------

=head2 comment

 Title   : comment
 Usage   : $a = $self->comment(qr/^#/);
 Function: Takes a regex to apply to a record to see if it looks like a
	   comment to skip.
 Returns : The stored regular expression
 Args    : A regex to apply to a record to see if it looks like a
	   comment to skip.

=cut

sub comment {
  my ($self, $comment) = @_;
  $self->{comment} = $comment if $comment;
  return $self->{comment};
}

#-----------------------------------------------------------------------------

=head2 bind_fields

 Title   : bind_fields
 Usage   : $a = $self->bind_fields();
 Function: Takes an array of field names to use as the key values when
	   a hash is returned from C<next_record>.
 Returns : The array reference of field names used as key values for hashes
	   returned by C<next_record>.
 Args    : An array of field names to use as the key values for hashes
	   returned from C<next_record>.

=cut

sub bind_fields {
  my ($self, $bind_fields) = @_;
  $self->{bind_fields} = $bind_fields if $bind_fields;
  return wantarray ? @{$self->{bind_fields}} : $self->{bind_fields};
}

#-----------------------------------------------------------------------------

=head2 _external_reader

 Title   : _external_reader
 Usage   : $a = $self->_external_reader();
 Function: Get the external_reader.
 Returns : A Text::RecordParser object as a singleton.
 Args    : None

=cut

sub _external_reader {
  my $self = shift;
  if (! $self->{_external_reader}) {
    my $_external_reader;
    my $reader_args = {fh               => $self->fh,
		       record_separator => $self->record_separator,
		       field_separator  => $self->field_separator,
		       comment          => $self->comment,
		      };

    $_external_reader = Text::RecordParser->new($reader_args);
    $_external_reader->bind_fields($self->bind_fields);
    $self->{_external_reader} = $_external_reader;
  }
  return $self->{_external_reader};
}

#-----------------------------------------------------------------------------

=head2 next_record

 Title   : next_record
 Usage   : $a = $self->next_record();
 Function: Return the next record from the _external_reader
 Returns : The next record from the _external_reader.
 Args    : N/A

=cut

sub next_record {
	my $self = shift;
	return $self->_external_reader->fetchrow_hashref;
}

#-----------------------------------------------------------------------------

=head1 DIAGNOSTICS

This module does not throw any errors or warning messages.

=head1 CONFIGURATION AND ENVIRONMENT

<GAL::Reader::RecordParser> requires no configuration files or environment variables.

=head1 DEPENDENCIES

<GAL::Reader>
<Text::RecordParser>

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
