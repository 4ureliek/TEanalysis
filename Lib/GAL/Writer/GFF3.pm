package GAL::Writer::GFF3;

use strict;
use vars qw($VERSION);

$VERSION = 0.2.0;
use base qw(GAL::Writer);

=head1 NAME

GAL::Writer::GFF3 - Writer::GFF3 objects for the Genome Annotation Library

=head1 VERSION

This document describes GAL::Writer::GFF3 version 0.2.0

=head1 SYNOPSIS

  use GAL::Writer::GFF3;
  $writer = GAL::Writer::GFF3->new(field_names => \@field_names);

=head1 DESCRIPTION

<GAL::Writer::GFF3>, via it's subclasses, provides file writing access for a
variety of file formats.  <GAL::Writer::GFF3> should not be instantiated on
it's own, but rather acts as a base class for functionality common to
all writers.

=head1 CONSTRUCTOR

New <GAL::Writer::GFF3> objects are instatiated via its' subclasses.
Arguments should be passed to the constructor as a list (or reference)
of key value pairs.  All attributes of the Writer::GFF3 object can be set in
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

Other attributes are used by subclasses of <GAL::Writer::GFF3>.  See those modules
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

  my $err_code = 'method_not_fully_implimented';
  my $caller = ref $self;
  my $err_msg = ("The method GAL::Writer::GFF3::get_metadata_text " .
		 "is not fully implimented yet.  Send an "          .
		 "angry e-mail to the author of $caller!");

  $self->warn($err_code, $err_msg);

  return '##gff-version 3';

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

  my $err_code = 'method_not_fully_implimented';
  my $caller = ref $self;
  my $err_msg = ("The method GAL::Writer::GFF3::get_comment_text " .
		 "is not fully implimented yet.  Send an "         .
		 "angry e-mail to the author of $caller!");
  $self->warn($err_code, $err_msg);
  return;
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
  my ($self, $feature) = @_;

  my $attribute_text;
  for my $key (sort _sort_attribute_keys keys %{$feature->{attributes}}) {
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
  my $feature
}

#--------------------------------------------------------------------------------

sub _sort_attribute_keys {

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

  return ((($ATTRB_ORDER{$a} || 99) <=>
	   ($ATTRB_ORDER{$b} || 99)) ||
	  $a cmp $b);

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

	my $err_code = 'method_must_be_overridden : GAL::Writer::GFF3::get_features_text';
	my $caller = ref $self;
	my $err_msg  = ("The method GAL::Writer::GFF3::features_text must be "       .
			"overridden by subclasses of GAL::Writer::GFF3.  Send an " .
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

  # Continue here.  Add a datasource method
  #my $datasource = $self->datasource;
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

	my $err_code = 'method_must_be_overridden : GAL::Writer::GFF3::write_comments';
	my $caller = ref $self;
	my $err_msg  = ("The method GAL::Writer::GFF3::comments must be "       .
			"overridden by subclasses of GAL::Writer::GFF3.  Send an " .
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
  my ($self, $feature) = @_;
  print $self->get_feature_text($feature);
  print "\n";
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

	my $err_code = 'method_must_be_overridden : GAL::Writer::GFF3::write_features';
	my $caller = ref $self;
	my $err_msg  = ("The method GAL::Writer::GFF3::features must be "       .
			"overridden by subclasses of GAL::Writer::GFF3.  Send an " .
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

	my $err_code = 'method_must_be_overridden : GAL::Writer::GFF3::write_file';
	my $caller = ref $self;
	my $err_msg  = ("The method GAL::Writer::GFF3::file must be "       .
			"overridden by subclasses of GAL::Writer::GFF3.  Send an " .
			"angry e-mail to the author of $caller!");
	$self->throw($err_code, $err_msg);
}

#--------------------------------------------------------------------------------

=head1 DIAGNOSTICS

=over

=item C<< failed_to_open_file >>

GAL::Writer::GFF3 tried to open a file, but failed.

=item C<< method_must_be_overridden >>

GAL::Parser::next_record must be overridden by a subclasses of
<GAL::Parser>, but this was not done.  Please contact the author of
the subclass that you were using.

=back

=head1 CONFIGURATION AND ENVIRONMENT

<GAL::Writer::GFF3> requires no configuration files or environment variables.

=head1 DEPENDENCIES

<GAL::Base>

Subclasses of GAL::Writer::GFF3 may have additional dependencies.

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
