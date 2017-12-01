package GAL::List::Categorical;

use strict;
use vars qw($VERSION);
use List::Util;
use List::MoreUtils;

# List::UtilsBy qw(sort_by nsort_by rev_sort_by rev_nsort_by max_by
#		 max_by min_by min_by uniq_by partition_by count_by
#		 zip_by unzip_by extract_by weighted_shuffle_by
#		 bundle_by);

$VERSION = 0.2.0;
use base qw(GAL::List);

=head1 NAME

GAL::List::Categorical - Provide functions for categorical lists

=head1 VERSION

This document describes GAL::List::Categorical version 0.2.0

=head1 SYNOPSIS

    use GAL::List::Categorical;
    my $list = GAL::List::Categorical->new(list => [qw(red blue green)];
    my $count    = $list_catg->count;
    $list_ref    = $list_catg->list;
    @list_ary    = $list_catg->list;
    $class       = $list_catg->class;
    $cardinality = $list_catg->cardinality
    $count_uniq  = $list_catg->count_uniq;
    $max_str     = $list_catg->maxstr;
    $min_str     = $list_catg->minstr;
    @shff_list   = $list_catg->shuffle;
    @uniq_list   = $list_catg->uniq;
    $item        = $list_catg->random_pick;

    # Create a list by passing an iterator that has a 'next' method.
    # The values of the list will be accessed for each object via the
    # method passed by the method attribute.

    my $obj_list = GAL::List::Categorical->new(iterator => \@objs
					       method   => 'type');

    # Create a list by calling the given method on each object in
    # list.

    my $obj_list = GAL::List::Categorical->new(list   => \@objs
					       method => 'type');

    # Create a list by getting the value of the given key from each
    # hash in the list.  This is an alias for the method attribute -
    # GAL::List::Categorical will check each item in the list to see
    # if it is an object or a hash reference and do the right thing.

    my $hash_list = GAL::List::Categorical->new(list => \@objs
						key  => 'type');

=head1 DESCRIPTION

<GAL::List::Categorical> provides methods to work with list of
categorical data.

=head1 CONSTRUCTOR

To construct a GAL::List::Categorical simply pass a list as an array
reference to the list attribute.  See the method and key attributes for
additional ways to pass a list.

    my $list = GAL::List::Categorical->new(list => [qw(red blue green)];

=cut

#-----------------------------------------------------------------------------
#-------------------------------- Constructor --------------------------------
#-----------------------------------------------------------------------------

=head2 new

     Title   : new
     Usage   : GAL::List::Categorical->new()
     Function: Creates a GAL::List::Categorical object;
     Returns : A GAL::List::Categorical object
     Args    : list => \@list
	       But see the method, key, and iterator attributes for other
	       ways to construct a list.

=cut

sub new {
	my ($class, @args) = @_;
	my $self = $class->SUPER::new(@args);
	$self->_parse_list;
	return $self;
}

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
}

#-----------------------------------------------------------------------------
#--------------------------------- Attributes --------------------------------
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
#------------------------------------ Methods --------------------------------
#-----------------------------------------------------------------------------

=head2 count_table

 Title   : count_table
 Usage   : $a = $self->count_table()
 Function: Returns a string with a text based table of the counts of
           categories.
 Returns : A string
 Args    : A scalar containing the name for the list.

=cut

sub count_table {
  my ($self, $header) = @_;

  my $table = Text::Table->new('|', $header, '|', 'Count', '|');
  my @data;
  my %counts = $self->category_counts();
  for my $key (sort keys %counts) {
    push @data, ['|', $key, '|', $counts{$key}, '|'];
  }
  $table->load(@data);
  my $title_rule  = $table->rule('=', '=');
  my $body_rule   = $table->rule('-', '+');
  my $bottom_rule = $table->rule('-', '-');
  my $text_table;
  $text_table .= $title_rule;
  $text_table .= $table->title;
  $text_table .= $title_rule;
  $text_table .= join $body_rule, $table->body;
  $text_table .= $bottom_rule;
  return $text_table;
}

#-----------------------------------------------------------------------------

=head2 count

 Title   : count
 Usage   : $count = $self->count()
 Function: Return the count of elements in the list.
 Returns : A integer
 Args    : N/A

=cut

#-----------------------------------------------------------------------------

=head2 cardinality

 Title   : cardinality
 Usage   : $a = $self->cardinality()
 Function: Returns the number of uniq elements in the list
 Returns : An integer
 Args    : N/A

=cut

#-----------------------------------------------------------------------------

=head2 count_uniq

 Title   : count_uniq
 Usage   : $a = $self->count_uniq()
 Function: An alias for cardinality - returns the number of uniq elements in
	   the list
 Returns : An integer
 Args    :

=cut

#-----------------------------------------------------------------------------

=head2 category_counts

 Title   : category_counts
 Usage   : $a = $self->category_counts()
 Function: Return a hash(reference) with uniq elements as keys and element
	   counts as values.
 Returns : A hash(reference)
 Args    : N/A

=cut

sub category_counts {
    my $self = shift;
    my %hash = $self->_hash_list;
    return wantarray ? %hash : \%hash;
}

#-----------------------------------------------------------------------------

=head2 maxstr

 Title   : maxstr
 Usage   : $a = $self->maxstr()
 Function: Return the max value of the list sorted alphanumerically
 Returns : Scalar string
 Args    : N/A

=cut

sub maxstr {
  my $self = shift;
  my @list = $self->list;
  return List::Util::maxstr(@list);
}

#-----------------------------------------------------------------------------

=head2 minstr

 Title   : minstr
 Usage   : $a = $self->minstr()
 Function: Return the min value of the list sorted alphanumerically
 Returns : Scalar string
 Args    : N/A

=cut

sub minstr {
  return List::Util::minstr(shift->list);
}

#-----------------------------------------------------------------------------

=head2 shuffle

 Title   : shuffle
 Usage   : $a = $self->shuffle()
 Function: Returns the elements of the list in random order
 Returns : An array/array reference
 Args    : N/A

=cut

#-----------------------------------------------------------------------------

=head2 uniq

 Title   : uniq
 Usage   : $a = $self->uniq()
 Function: Return the uniq elements in the list
 Returns : An array/reference
 Args    :

=cut

#-----------------------------------------------------------------------------

=head2 random_pick

 Title   : random_pick
 Usage   : $a = $self->random_pick()
 Function: Return and element from the list at random.
 Returns : A scalar element from the list
 Args    : N/A

=cut

#-----------------------------------------------------------------------------

=head1 DIAGNOSTICS

<GAL::List::Categorical> currently does not throw any warnings or
errors.

=head1 CONFIGURATION AND ENVIRONMENT

<GAL::List::Categorical> requires no configuration files or environment variables.

=head1 DEPENDENCIES

<BAL::Base>

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
