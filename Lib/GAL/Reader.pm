package GAL::Reader;

use strict;
use vars qw($VERSION);

$VERSION = 0.2.0;
use base qw(GAL::Base);

=head1 NAME

GAL::Reader - Reader objects for the Genome Annotation Library

=head1 VERSION

This document describes GAL::Reader version 0.2.0

=head1 SYNOPSIS

  use GAL::Reader::DelimitedLine;
  $reader = GAL::Reader::DelimitedLine->new(field_names => \@field_names);

=head1 DESCRIPTION

<GAL::Reader>, via it's subclasses, provides file reading access for a
variety of file formats.  The reader objects don't parse the
information in the files themselves, but rather provide standard
access to broad categories of formats such as tab-delimited,
multi-line record, XML files and others.  <GAL::Reader> should not be
instantiated on it's own, but rather acts as a base class for
functionality common to all readers.

=head1 CONSTRUCTOR

New GAL::Reader::subclass objects are created by the class method new.
Arguments should be passed to the constructor as a list (or reference)
of key value pairs.  All attributes of the Reader object can be set in
the call to new. An simple example of object creation would look like
this:

  $reader = GAL::Reader::DelimitedLine->new(field_names => \@field_names);

The constructor recognizes the following parameters which will set the
appropriate attributes:

=over 4

=item * C<< file => feature_file.txt >>

This optional parameter defines what file to parse. While this
parameter is optional either it, or the following fh parameter must be
set before the first call to next_record.

=item * C<< fh => $FH >>

This optional parameter provides a file handle to parse. While this
parameter is optional, either it or the previous must be set before
the first call to next_record.

=back

Other attributes are used by subclasses of <GAL::Reader>.  See those modules
for details.

=cut

#-----------------------------------------------------------------------------
#-------------------------------- Constructor --------------------------------
#-----------------------------------------------------------------------------

=head2 new

     Title   : new
     Usage   : GAL::Reader::DelimitedLine->new();
     Function: Creates a GAL::Reader::DelimitedLine object;
     Returns : A GAL::Reader::DelimitedLine object
     Args    : See <GAL::Reader::DelimitedLine> and other subclasses.

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
	my @valid_attributes = qw(file fh);
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

=head2 file

 Title   : file
 Usage   : $a = $self->file();
 Function: Get/Set the value of file.
 Returns : The value of file.
 Args    : The filename of a readable file.

=cut

sub file {
  my ($self, $value) = @_;
  $self->{file} = $value if defined $value;
  return $self->{file};
}

#-----------------------------------------------------------------------------

=head2 fh

 Title   : fh
 Usage   : $a = $self->fh();
 Function: Get/Set the filehandle.  Once the file handle is set, the same
	   file handle is returned until another file handle is passed in.
 Returns : A filehandle
 Args    : A filehandle

=cut

sub fh {
  my ($self, $fh) = @_;
  $self->{fh} = $fh if defined $fh;
  if (! $self->{fh}) {
    my $file = $self->file;
    if ($file) {
      # TODO $self->open_file($file, 'read');
      $self->throw('file_doesnt_exist', $file) if ! -e $file;
      $self->throw('unreadable_file', $file)   if ! -r $file;
      open($fh, '<', $file) or $self->throw('failed_to_open_file', $file);
      $self->{fh} = $fh;
    }
  }
  return $self->{fh};
}

#-----------------------------------------------------------------------------
#---------------------------------- Methods ----------------------------------
#-----------------------------------------------------------------------------

=head1 METHODS

=head2 _external_reader

 Title   : _external_reader
 Usage   : $a = $self->_external_reader();
 Function: Get/Set the _external_reader if one is used.  This allows, for
	   example, Text::RecordParser or XML::Twig to be easily added as
	   an external reader by subclasses.
 Returns : The value of _external_reader.
 Args    : A value to set _external_reader to.

=cut

sub _external_reader {
	my ($self, $reader) = @_;
	$self->{_external_reader} = $reader if $reader;
	return $self->{_external_reader};

}

#-----------------------------------------------------------------------------

=head2 next_record

 Title   : next_record
 Usage   : $a = $self->next_record();
 Function: Return the next record from the external_reader.  Note that this
	   method must be overridden by a sublass of GAL::Reader.
 Returns : The next record from the external_reader.
 Args    : N/A

=cut

sub next_record {
	my $self = shift;
	my $err_code = 'method_must_be_overridden : GAL::Reader::next_record';
	my $caller = ref $self;
	my $err_msg  = ("The method GAL::Reader::next_record must be "       .
			"overridden by subclasses of GAL::Reader.  Send an " .
			"angry e-mail to the author of $caller!");
	$self->throw($err_code, $err_msg);
}

#-----------------------------------------------------------------------------

=head1 DIAGNOSTICS

=over

=item C<< file_doesnt_exist >>

GAL::Reader tried to open a file that did not exist.  Please check you
file path and filename.

=item C<< cant_read_file >>

GAL::Reader tried to open a file that exists, but could not be read.
Please check your file permissions.

=item C<< failed_to_open_file >>

GAL::Reader tried to open a file, but failed.

=item C<< method_must_be_overridden >>

GAL::Parser::next_record must be overridden by a subclasses of
<GAL::Parser>, but this was not done.  Please contact the author of
the subclass that you were using.

=back

=head1 CONFIGURATION AND ENVIRONMENT

<GAL::Reader> requires no configuration files or environment variables.

=head1 DEPENDENCIES

<GAL::Base>

Subclasses of GAL::Reader may have additional dependencies.

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
