package GAL::List;

use strict;
use vars qw($VERSION);

$VERSION = 0.2.0;
use base qw(GAL::Base);
use List::Util;
use List::MoreUtils;
use Scalar::Util qw(blessed);
use Text::Table;

=head1 NAME

C<GAL::List> - List aggregation and analysis functions for GAL

=head1 VERSION

This document describes C<GAL::List> version 0.2.0

=head1 SYNOPSIS

    use GAL::List::Categorical;
    my $list = GAL::List::Categorical->new(list => [qw(red red red blue blue
						       green yellow orange orange
						       purple purple purple purple)]);
    my $count    = $list->count;
    $list_ref    = $list->list;
    @list_ary    = $list->list;
    $catg_counts = $list->cardinality;
    $count_uniq  = $list->count_uniq;
    $max_str     = $list->maxstr;
    $min_str     = $list->minstr;
    @shff_list   = $list->shuffle;
    @uniq_list   = $list->uniq;
    $item        = $list->random_pick;

    use GAL::List::Numric;
    $list = GAL::List::Numeric->new(list => [0, 1, 2, 3, 4, 5, 6,
                                             7, 8, 9]);
    $max = $list->max;
    $min = $list->min;
    $sum = $list->sum;

=head1 DESCRIPTION

C<GAL::List> serves as a base class for the modules below it and
provides basic list summarization details.  It is not intended to be
used on it's own.  You should use it's subclasses instead.

=head1 CONSTRUCTOR

To construct a C<GAL::List> subclass simply pass it an appropriate list.

  my $colors = GAL::List::Categorical->new(list => [qw(red red red
                                                       blue green
						       yellow orange
						       orange purple
						       purple)]);

  my $numbers = GAL::List::Numeric->new(list => [0, 1, 2, 3, 4, 5, 6,
                                                 7, 8, 9]);


=cut

#-----------------------------------------------------------------------------
#-------------------------------- Constructor --------------------------------
#-----------------------------------------------------------------------------

=head2 new

     Title   : new
     Usage   : GAL::List->new(list => $array_ref)
     Function: Creates a C<GAL::List> object;
     Returns : A C<GAL::List> object
     Args    : Any valid attributes as described below.

=cut

sub new {
	my ($class, @args) = @_;
	my $self = $class->SUPER::new(@args);
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
	my @valid_attributes = qw(list iterator method);
	$self->set_attributes($args, @valid_attributes);
	######################################################################
}

#-----------------------------------------------------------------------------
#-------------------------------- Attributes ---------------------------------
#-----------------------------------------------------------------------------

=head1  ATTRIBUTES

All attributes can be supplied as parameters to the constructor as a
list (or referenece) of key value pairs.

=head2 list

 Title   : list
 Usage   : $ref   = $self->list();
           @array = $self->list();
 Function: Get/Set the contents of the objects list.
 Returns : The elements of the list as an array or a array reference.
 Args    : One of the following:
               - An array reference of scalar values: The values will
                 be added to the list.
               - An array reference of hashes: A key name must also be
                 supplied to L<GAL::List::key|/key> and this key will be used
                 to retrieve data from each hash reference.
               - An array reference of objects: A method name must also be
                 supplied to L<GAL::List::method|/method> and this method will be used
                 to retrieve data from each hash reference.
=cut

sub list {
  my ($self, $list) = @_;
  if ($list) {
    my $err_msg = ('GAL::List requires an array reference as the ' .
		   'first argument, but you gave a  ' . ref $list
		  );
    my $err_code = 'list_or_reference_required : ' . ref $list;
    $self->warn($err_code, $err_msg) unless ref $list eq 'ARRAY';
    $self->{list} = $list;
  }
  $self->{list} ||= [];
  return wantarray ? @{$self->{list}} : $self->{list};
}

#-----------------------------------------------------------------------------

=head2 method

 Title   : method
 Usage   : my $list = GAL::List(iterator => $features,
                                method   => 'length');
 Function: When an iterator object is passed to the C<GAL::List>
 	   constructor or object (see the iterator attribute below),
 	   this attribute sets the method used to retrieve the
 	   appropriate value from the iterator.
 Returns : The value of method.
 Args    : A value to set method to.

=cut

sub method {
  my ($self, $method) = @_;
  $self->{method} = $method if $method;
  return $self->{method};
}

#-----------------------------------------------------------------------------

=head2 key

 Title   : key
 Usage   : my $list = GAL::List(list => $hash_ref,
                                key  => score)
 Function: When a hash reference is passed to the list attribute, this
 	   attribute sets the key used to retrieve the appropriate
 	   value from each hash.
 Returns : The value of key.
 Args    : A value to set key to.

=cut

sub key {
  my ($self, $key) = @_;
  $self->{key} = $key if $key;
  return $self->{key};
}

#-----------------------------------------------------------------------------

=head2 iterator

 Title   : iterator
 Usage   : $a = $self->iterator()
 Function: Get/Set the value of iterator.
 Returns : The value of iterator.
 Args    : A value to set iterator to.

=cut

 sub iterator {
   my ($self, $iterator) = @_;
   $self->{iterator} = $iterator if $iterator;
   return $self->{iterator};
 }

#-----------------------------------------------------------------------------
#---------------------------------- Methods ----------------------------------
#-----------------------------------------------------------------------------

=head1 METHODS

=cut

sub _parse_list {
  my $self = shift;

  my $list     = $self->list();
  my $iterator = $self->iterator();
  my $method   = $self->method();
  my $key      = $self->method();

 VALUE:
  for my $value (@{$list}) {
    if (! ref $value) {
      next VALUE;
    }
    elsif (my $class = Scalar::Util::blessed $value) {
      if (! $method) {
	$self->throw('no_method_given_for_objects', "CLASS=$class; METHOD=$method");
      }
      if (! defined $method || ! $value->can($method)) {
	$self->warn('method_does_not_exist', "CLASS=$class; METHOD=$method");
	next VALUE;
      }
      $value = $value->$method;
    }
    elsif (ref $value eq 'HASH') {
      if (! defined $key || ! exists $value->{$key}) {
	$self->warn('key_does_not_exist', ("KEY=$key; "));
	next VALUE;
      }
      $value = $value->{$key};
    }
  }

  my @xtra_values;
  if ($iterator) {
    if (my $class = Scalar::Util::blessed $iterator) {
      if (! defined $method || ! $iterator->can($method)) {
	$self->warn('method_does_not_exist', "CLASS=$class; METHOD=$method");
      }
      else {
      ITR_LOOP:
	while (my $obj = $iterator->next) {
	  my $value = $obj->method;
	  push @xtra_values, $value;
	}
      }
    }
  }

  push @{$self->{list}}, @xtra_values if @xtra_values;
}

#-----------------------------------------------------------------------------

sub _hash_list {
  my $self = shift;

  if (! exists $self->{hash_list}) {
    map {$self->{hash_list}{$_}++} $self->list;
  }

  return wantarray ? %{$self->{hash_list}} : $self->{hash_list};

}

#-----------------------------------------------------------------------------

=head2 count

 Title   : count
 Usage   : $list = $self->count
 Function: Return the number of elements in the list
 Returns : Integer
 Args    : N/A

=cut

sub count {
  my @list = shift->list;
  return scalar @list;
}

#-----------------------------------------------------------------------------

=head2 cardinality

 Title   : cardinality
 Usage   : $list = $self->cardinality
 Function: Return the number of uniq elements in the list
 Returns : Integer
 Args    : N/A

=cut

sub cardinality {
  my $self = shift;
  my @uniq = $self->uniq;
  return scalar @uniq;
}

#-----------------------------------------------------------------------------

=head2 count_uniq

 Title   : count_uniq
 Usage   : $list = $self->count_uniq()
 Function: An alias for cardinality
 Returns : Integer
 Args    : N/A

=cut

sub count_uniq {
    return shift->cardinality;
}

#-----------------------------------------------------------------------------

=head2 category_counts

 Title   : category_counts
 Usage   : $list = $self->category_counts()
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

=head2 shuffle

 Title   : shuffle
 Usage   : $list = $self->shuffle()
 Function: Returns the elements of the list in random order
 Returns : An array/array reference
 Args    : N/A

=cut

sub shuffle {
  return wantarray ? List::Util::shuffle(shift->list) :
    [List::Util::shuffle(shift->list)];
}

#-----------------------------------------------------------------------------

=head2 in_place_shuffle

 Title   : in_place_shuffle
 Usage   : $list = $self->in_place_shuffle()
 Function: Shuffles the elements of a list in place
 Returns : N/A
 Args    : N/A

=cut

sub in_place_shuffle {
  my $self = shift;
  $self->list([$self->shuffle]);
}

#-----------------------------------------------------------------------------

=head2 uniq

 Title   : uniq
 Usage   : $list = $self->uniq()
 Function: Returns uniq values from the list.
 Returns : An array (ref) of unique values from the list.
 Args    : N/A

=cut

sub uniq {
  my $self = shift;
  my @uniq = List::MoreUtils::uniq($self->list);
  return wantarray ? @uniq : \@uniq;
}

#-----------------------------------------------------------------------------

=head2 random_pick

 Title   : random_pick
 Usage   : $list = $self->random_pick(8)
 Function: Return a random selection from the list of a given length.
 Returns : An array (ref) of randomly selected values from the list.
 Args    : An integer for the size of the returned list. Default is one.

=cut

sub random_pick {
    my ($self, $count) = @_;

    my @random_list;
    for (1 .. $count) {
      my $random = int(rand($self->count));
      push @random_list, $self->{list}[$random];
    }
    return wantarray ? @random_list : \@random_list;
}

#-----------------------------------------------------------------------------

=head1 DIAGNOSTICS

=over

=item C<< list_or_reference_required >>

L<GAL::List::list|/list> require an array or a reference to any array
be passed as an argument, but you have passed something else.

Keep in mind that several of C<GAL::List>'s methods are provided by
List::Util, and errors not found here may be thrown by that module.

=item C<< no_method_given_for_iterator' >>

C<GAL::List> can accept an iterator object that provides a next
method.  For each object returned by calling next C<GAL::List> will
call the method provided to L<GAL::List::method|/method>.  A method name
must be provided to L<GAL::List::method|/method> with an iterator is
provided to L<GAL::List::iterator|/iterator>.

=item C<< method_does_not_exist >>

An iterator was provided to L<GAL::List::iterator|/iterator> and a
method was provided to L<GAL::List::method|/method>, but the objects
returned by the iterator do not impliment the method provided to
L<GAL::List::method|/method>.

=item C<< key_does_not_exist' >>

If a hash reference is provided as the value for L<GAL::List::list|/list> then
a key must be provided to L<GAL::List::key|/key>. C<GAL::List> will then use the
key to retrieve it's value from the hash.  A hash was provided to
L<GAL::List::list|/list>, but no key was given.

=back

=head1 CONFIGURATION AND ENVIRONMENT

C<GAL::List> requires no configuration files or environment variables.

=head1 DEPENDENCIES

L<GAL::Base>
L<List::Util>

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
