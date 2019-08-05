package GAL::Reader::Template;

use strict;
use vars qw($VERSION);


$VERSION = 0.2.0;
use base qw(GAL::Reader);

=head1 NAME

GAL::Reader::Template - <One line description of module's purpose here>

=head1 VERSION

This document describes GAL::Reader::Template version 0.2.0

=head1 SYNOPSIS

    use GAL::Reader::Template;
    $reader = GAL::Reader::Template->new(file => 'annotation_file.gff');
    $reader->next_record, '$reader->next_record');

=head1 DESCRIPTION

<GAL::Reader::Template> provides << BLAH BLAH BLAH >> reading
capability to GAL.  It is not intended for general library use, but
rather as a GAL::Reader subclass for developers of GAL::Parser
subclasses.  There is however no reason why it couldn't also be used
as a stand alone module for other purposes.

=head1 CONSTRUCTOR

New GAL::Reader::Template objects are created by the class method new.
Arguments should be passed to the constructor as a list (or reference)
of key value pairs.  All attributes of the Reader object can be set in
the call to new. An simple example of object creation would look like
this:

    $reader = GAL::Reader::Template->new(file => 'annotation_file.gff');

The constructor recognizes the following parameters which will set the
appropriate attributes:

=over 4

=item * C<< some_attribute => $some_value >>

This optional attribute provides << BLAH BLAH BLAH >>


The following attributes are inhereted from GAL::Reader:


=item * C<< file => feature_file.txt >>

This optional parameter defines what file to parse. While this
parameter is optional either it, or the following fh parameter must be
set before the first call to next_record.

=item * C<< fh => $FH >>

This optional parameter provides a file handle to parse. While this
parameter is optional, either it or the previous must be set before
the first call to next_record.

=back

=cut

#-----------------------------------------------------------------------------
#------------------------------- Constructor ---------------------------------
#-----------------------------------------------------------------------------

=head2 new

     Title   : new
     Usage   : GAL::Reader::Template->new();
     Function: Creates a GAL::Reader::Template object;
     Returns : A GAL::Reader::Template object
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
	my @valid_attributes = qw();
	$self->set_attributes($args, @valid_attributes);
	######################################################################
	return $args;
}

#-----------------------------------------------------------------------------
#------------------------------- Attributes ----------------------------------
#-----------------------------------------------------------------------------
#
# =head2 attribute
#
#  Title   : attribute
#  Usage   : $a = $self->attribute()
#  Function: Get/Set the value of attribute.
#  Returns : The value of attribute.
#  Args    : A value to set attribute to.
#
# =cut
#
# sub attribute {
#   my ($self, $attribute) = @_;
#
#   if ($attribute) {
#   $self->{attribute} ||= 'GAL::Parser::gff3';
#   return $self->{attribute};
# }
#

#-----------------------------------------------------------------------------
#---------------------------------- Methods ----------------------------------
#-----------------------------------------------------------------------------

=head2 next_record

 Title   : next_record
 Usage   : $a = $self->next_record();
 Function: Return the next record from the reader
 Returns : The next record from the reader.
 Args    : N/A

=cut

sub next_record {
	my $self = shift;
	my $line;
	my $fh = $self->fh;
	while (my $line = <$fh>) {
	  last unless $line =~ /\s*\#/;
	}
	my @fields = split /\t/, $line;
	return wantarray ? @fields : \@fields;
}

#-----------------------------------------------------------------------------

=head1 DIAGNOSTICS

=for author to fill in:
     List every single error and warning message that the module can
     generate (even the ones that will "never happen"), with a full
     explanation of each problem, one or more likely causes, and any
     suggested remedies.

=over

=item C<< Error message here, perhaps with %s placeholders >>

[Description of error here]

=item C<< Another error message here >>

[Description of error here]

[Et cetera, et cetera]

=back

=head1 CONFIGURATION AND ENVIRONMENT

<GAL::Reader::Template> requires no configuration files or environment variables.

=head1 DEPENDENCIES

None.

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
