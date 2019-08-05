package GAL::Writer;

use strict;
use vars qw($VERSION);

$VERSION = 0.2.0;
use base qw(GAL::Base);

=head1 NAME

GAL::Writer - Writer objects for the Genome Annotation Library

=head1 VERSION

This document describes GAL::Writer version 0.2.0

=head1 SYNOPSIS

  use GAL::Writer::GFF3;
  $writer = GAL::Writer::GFF3->new(field_names => \@field_names);

=head1 DESCRIPTION

<GAL::Writer>, via it's subclasses, provides file writing access for a
variety of file formats.  <GAL::Writer> should not be instantiated on
it's own, but rather acts as a base class for functionality common to
all writers.

=head1 CONSTRUCTOR

New <GAL::Writer> objects are instatiated via its' subclasses.
Arguments should be passed to the constructor as a list (or reference)
of key value pairs.  All attributes of the Writer object can be set in
the call to new or later as method calls. An simple example of object
creation would look like this:

  $writer = GAL::Writer::GFF3->new(file => $filename);

The constructor recognizes the following parameters which will set the
appropriate attributes:

=over 4

=item * C<< file => file.txt >>

This optional parameter defines the path/name of the output file to
write to.  If neither this attribute or the file attributes are set
the writer defaults to writing to STDOUT.

=item * C<< fh => $FH >>

This optional parameter provides a file handle open for writing. If
neither this attribute or the file attributes are set the writer
defaults to writing to STDOUT.

=back

Other attributes are used by subclasses of <GAL::Writer>.  See those modules
for details.

=cut

#-----------------------------------------------------------------------------
#-------------------------------- Constructor --------------------------------
#-----------------------------------------------------------------------------

=head2 new

     Title   : new
     Usage   : GAL::Writer::GFF3->new();
     Function: Creates a GAL::Writer::GFF3 object;
     Returns : A GAL::Writer::GFF3 object
     Args    : See <GAL::Writer::GFF3> and other subclasses.

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
 Args    : The filename of a writeable file.

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
 Returns : A filehandle open for writing.
 Args    : A filehandle open for writing.

=cut

sub fh {
  my ($self, $fh) = @_;
  $self->{fh} = $fh if defined $fh;
  if (! $self->{fh}) {
    my $file = $self->file;
    if ($file) {
      # TODO $self->open_file($file, 'write');
      open($fh, '>', $file) or
	$self->throw('failed_to_open_file_for_writing', $file);
      $self->{fh} = $fh;
    }
  }
  return $self->{fh};
}

#-----------------------------------------------------------------------------

=head2 data_source

 Title   : data_source
 Usage   : $a = $self->data_source();
 Function: Get/Set the data source that will provide the data for the
	   writer to write.  The data source can be either a
	   GAL::Parser object or a GAL::Schema::ResultSet object.
 Returns : A filehandle open for writing.
 Args    : A filehandle open for writing.

=cut

sub data_source {
  my ($self, $data_source) = @_;
  $self->{data_source} = $data_source if defined $data_source;
  if (! $self->{data_source}) {
      $self->throw('no_data_source_provided');
  }
  return $self->{data_source};
}

#-----------------------------------------------------------------------------
#---------------------------------- Methods ----------------------------------
#-----------------------------------------------------------------------------

=head1 METHODS

=head2 get_metadata_text

 Title   : get_metadata_text
 Usage   : $a = $self->get_metadata_text();
 Function: Get the metadata structured appropriately as text.
 Returns : The metadata as a string of text.
 Args    : N/A

=cut

sub get_metadata_text {
  my ($self) = @_;

	my $err_code = 'method_must_be_overridden : GAL::Writer::get_metadata_text';
	my $caller = ref $self;
	my $err_msg  = ("The method GAL::Writer::metadata_text must be "       .
			"overridden by subclasses of GAL::Writer.  Send an " .
			"angry e-mail to the author of $caller!");
	$self->throw($err_code, $err_msg);
}

#--------------------------------------------------------------------------------

=head2 get_comment_text

 Title   : get_comment_text
 Usage   : $a = $self->get_comment_text();
 Function: Get the comment structured appropriately as text.
 Returns : The comment as a string of text.
 Args    : N/A

=cut

sub get_comment_text {
  my ($self) = @_;

	my $err_code = 'method_must_be_overridden : GAL::Writer::get_comment_text';
	my $caller = ref $self;
	my $err_msg  = ("The method GAL::Writer::comment_text must be "       .
			"overridden by subclasses of GAL::Writer.  Send an " .
			"angry e-mail to the author of $caller!");
	$self->throw($err_code, $err_msg);
}

#--------------------------------------------------------------------------------

=head2 get_feature_text

 Title   : get_feature_text
 Usage   : $a = $self->get_feature_text($feature);
 Function: Get the feature structured appropriately as text.
 Returns : The feature as a string of text.
 Args    : N/A

=cut

sub get_feature_text {
  my ($self) = @_;

	my $err_code = 'method_must_be_overridden : GAL::Writer::get_feature_text';
	my $caller = ref $self;
	my $err_msg  = ("The method GAL::Writer::feature_text must be "       .
			"overridden by subclasses of GAL::Writer.  Send an " .
			"angry e-mail to the author of $caller!");
	$self->throw($err_code, $err_msg);
}

#--------------------------------------------------------------------------------

=head2 get_features_text

 Title   : get_features_text
 Usage   : $a = $self->get_features_text();
 Function: Get the features structured appropriately as text.
 Returns : The all of the features as a string of text.
 Args    : N/A

=cut

sub get_features_text {
  my ($self) = @_;

	my $err_code = 'method_must_be_overridden : GAL::Writer::get_features_text';
	my $caller = ref $self;
	my $err_msg  = ("The method GAL::Writer::features_text must be "       .
			"overridden by subclasses of GAL::Writer.  Send an " .
			"angry e-mail to the author of $caller!");
	$self->throw($err_code, $err_msg);
}


#--------------------------------------------------------------------------------

=head2 write_metadata

 Title   : write_metadata
 Usage   : $a = $self->write_metadata();
 Function: Write the metadata.
 Returns : N/A
 Args    : N/A

=cut

sub write_metadata {
  my ($self) = @_;

	my $err_code = 'method_must_be_overridden : GAL::Writer::write_metadata';
	my $caller = ref $self;
	my $err_msg  = ("The method GAL::Writer::metadata must be "       .
			"overridden by subclasses of GAL::Writer.  Send an " .
			"angry e-mail to the author of $caller!");
	$self->throw($err_code, $err_msg);
}

#--------------------------------------------------------------------------------

=head2 write_comments

 Title   : write_comments
 Usage   : $a = $self->write_comments();
 Function: Write the comments.
 Returns : N/A
 Args    : N/A

=cut

sub write_comments {
  my ($self) = @_;

	my $err_code = 'method_must_be_overridden : GAL::Writer::write_comments';
	my $caller = ref $self;
	my $err_msg  = ("The method GAL::Writer::comments must be "       .
			"overridden by subclasses of GAL::Writer.  Send an " .
			"angry e-mail to the author of $caller!");
	$self->throw($err_code, $err_msg);
}

#--------------------------------------------------------------------------------

=head2 write_feature

 Title   : write_feature
 Usage   : $a = $self->write_feature($feature);
 Function: Write the given feature.
 Returns : A GAL::Schema::Result::Feature object or a feature hash.
 Args    : N/A

=cut

sub write_feature {
  my ($self) = @_;

	my $err_code = 'method_must_be_overridden : GAL::Writer::write_feature';
	my $caller = ref $self;
	my $err_msg  = ("The method GAL::Writer::feature must be "       .
			"overridden by subclasses of GAL::Writer.  Send an " .
			"angry e-mail to the author of $caller!");
	$self->throw($err_code, $err_msg);
}

#--------------------------------------------------------------------------------

=head2 write_features

 Title   : write_features
 Usage   : $a = $self->write_features();
 Function: Write all the features.
 Returns : N/A
 Args    : N/A

=cut

sub write_features {
  my ($self) = @_;

	my $err_code = 'method_must_be_overridden : GAL::Writer::write_features';
	my $caller = ref $self;
	my $err_msg  = ("The method GAL::Writer::features must be "       .
			"overridden by subclasses of GAL::Writer.  Send an " .
			"angry e-mail to the author of $caller!");
	$self->throw($err_code, $err_msg);
}

#--------------------------------------------------------------------------------

=head2 write_file

 Title   : write_file
 Usage   : $a = $self->write_file();
 Function: Write all of the metadata, comments, headers and features.
 Returns : N/A
 Args    : N/A

=cut

sub write_file {
  my ($self) = @_;

	my $err_code = 'method_must_be_overridden : GAL::Writer::write_file';
	my $caller = ref $self;
	my $err_msg  = ("The method GAL::Writer::file must be "       .
			"overridden by subclasses of GAL::Writer.  Send an " .
			"angry e-mail to the author of $caller!");
	$self->throw($err_code, $err_msg);
}

#--------------------------------------------------------------------------------

=head1 DIAGNOSTICS

=over

=item C<< failed_to_open_file >>

GAL::Writer tried to open a file, but failed.

=item C<< method_must_be_overridden >>

GAL::Parser::next_record must be overridden by a subclasses of
<GAL::Parser>, but this was not done.  Please contact the author of
the subclass that you were using.

=back

=head1 CONFIGURATION AND ENVIRONMENT

<GAL::Writer> requires no configuration files or environment variables.

=head1 DEPENDENCIES

<GAL::Base>

Subclasses of GAL::Writer may have additional dependencies.

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
