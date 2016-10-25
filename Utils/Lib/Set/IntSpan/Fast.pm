package Set::IntSpan::Fast;

use warnings;
use strict;
use Carp;
use Data::Types qw(is_int);
use List::Util qw(min max);

=head1 NAME

Set::IntSpan::Fast - Fast handling of sets containing integer spans.

=head1 VERSION

This document describes Set::IntSpan::Fast version 1.15

=cut

BEGIN {
  our $VERSION = '1.15';
  our @ISA;
  eval "use Set::IntSpan::Fast::XS ()";
  if ( $@ ) {
    if ( $@ =~ /^Can't\s+locate/ ) {
      eval "use Set::IntSpan::Fast::PP ()";
      die $@ if $@;
      @ISA = qw( Set::IntSpan::Fast::PP );
    }
    else {
      die $@;
    }
  }
  else {
    @ISA = qw( Set::IntSpan::Fast::XS );
  }
}

1;
__END__

=head1 SYNOPSIS

    use Set::IntSpan::Fast;
    
    my $set = Set::IntSpan::Fast->new();
    $set->add(1, 3, 5, 7, 9);
    $set->add_range(100, 1_000_000);
    print $set->as_string(), "\n";    # prints 1,3,5,7,9,100-1000000

=head1 DESCRIPTION

C<Set::IntSpan::Fast> represents sets of integers. It is optimised for
sets that contain contiguous runs of values.

    1-1000, 2000-10000      # Efficiently handled

Sets that don't have this characteristic may still be represented but
some of the performance and storage space advantages will be lost.
Consider using bit vectors if your set does not typically contain
clusters of values.

Sets may be infinite - assuming you're prepared to accept that infinity
is actually no more than a fairly large integer. Specifically the
constants C<Set::IntSpan::Fast::NEGATIVE_INFINITY> and
C<Set::IntSpan::Fast::POSITIVE_INFINITY> are defined to be -(2^31-1) and
(2^31-2) respectively. To create an infinite set invert an empty one:

    my $inf = Set::IntSpan::Fast->new()->complement();

Sets need only be bounded in one direction - for example this is the set
of all positive integers (assuming you accept the slightly feeble
definition of infinity we're using):

    my $pos_int = Set::IntSpan::Fast->new();
    $pos_int->add_range(1, $pos_int->POSITIVE_INFINITY);

=head2 Set representation

The internal representation used is extremely simple: a set is
represented as a list of integers. Integers in even numbered positions
(0, 2, 4 etc) represent the start of a run of numbers while those in odd
numbered positions represent the ends of runs. As an example the set (1,
3-7, 9, 11, 12) would be represented internally as (1, 2, 3, 8, 11, 13).

=head2 Comparision with Set::IntSpan

The C<Set::IntSpan> module represents sets of integers as a number of
inclusive ranges, for example '1-10,19-23,45-48'. Because many of its
operations involve linear searches of the list of ranges its overall
performance tends to be proportional to the number of distinct ranges.
This is fine for small sets but suffers compared to other possible set
representations (bit vectors, hash keys) when the number of ranges
grows large.

This module also represents sets as ranges of values but stores those
ranges in order and uses a binary search for many internal operations
so that overall performance tends towards O log N where N is the number
of ranges.

=head2 C<Set::IntSpan::Fast::XS>

If L<Set::IntSpan::Fast::XS> is installed it will automatically be used
when a new C<Set::IntSpan::Fast> is created. There is no need to change
any code; the XS module is automatically detected and loaded.

If you have a C compiler consider installing C<Set::IntSpan::Fast::XS>
for even better performance.

=head1 INTERFACE

=head2 C<new>

Create a new set. Any arguments will be processed by a call to
C<add_from_string>:

    my $set = Set::IntSpan::Fast->new( '1, 3, 5, 10-100' );

Because C<add_from_string> handles multiple arguments this will work:

    my @nums = ( 1, 2, 3, 4, 5 );
    my $set = Set::IntSpan::Fast->new( @nums );

Bear in mind though that this validates each element of the array is it
would if you called C<add_from_string> so for large sets it will be
slightly more efficient to create an empty set and then call C<add>.

=head2 C<invert>

Complement the set. Because our notion of infinity is actually
disappointingly finite inverting a finite set results in another finite
set. For example inverting the empty set makes it contain all the
integers between C<NEGATIVE_INFINITY> and C<POSITIVE_INFINITY> inclusive.

As noted above C<NEGATIVE_INFINITY> and C<POSITIVE_INFINITY> are actually just
big integers.

=head2 C<copy>

Return an identical copy of the set.

    my $new_set = $set->copy();

=cut

sub copy {
    my $self  = shift;
    my $class = ref $self;
    my $copy  = $class->new();
    # This might not work for subclasses - in which case they should
    # override copy
    @$copy = @$self;
    return $copy;
}

=head2 C<add( $number ... )>

Add the specified integers to the set. Any number of arguments may be
specified in any order. All arguments must be integers between
C<Set::IntSpan::NEGATIVE_INFINITY> and C<Set::IntSpan::POSITIVE_INFINITY>
inclusive.

=head2 C<remove( $number ... )>

Remove the specified integers from the set. It is not an error to remove
non-members. Any number of arguments may be specified.

=head2 C<add_range( $from, $to )>

Add the inclusive range of integers to the set. Multiple ranges may be
specified:

    $set->add_range(1, 10, 20, 22, 15, 17);

Each pair of arguments constitute a range. The second argument in each
pair must be greater than or equal to the first.

=head2 C<add_from_string( $string )>

Add items to a set from a string representation, of the same form as
C<as_string>. Multiple strings may be supplied:

    $set->add_from_string( '1-10, 30-40', '100-200' );

is equivalent to

    $set->add_from_string( '1-10, 30-40, 100-200' );

By default items are separated by ',' and ranges delimited by '-'. You
may select different punctuation like this:

    $set->add_from_string( 
        { sep => ';', range => ':' },
        '1;3;5;7:11;19:27'
    );
    
When supplying an options hash in this way the C<sep> and C<range>
option may be either a regular expression or a literal string.

    $set->add_from_string( 
        { sep => qr/:+/, range => qr/[.]+/ },
        '1::3::5:7...11:19..27'
    );

Any embedded whitespace in the string will be ignored.

=head2 C<remove_range( $from, $to )>

Remove the inclusive range of integers from the set. Multiple ranges may
be specified:

    $set->remove_range(1, 10, 20, 22, 15, 17);

Each pair of arguments constitute a range. The second argument in each
pair must be greater than or equal to the first.

=head2 C<remove_from_string( $string )>

Remove items to a set from a string representation, of the same form as
C<as_string>. As with C<add_from_string> the punctuation characters may
be specified.

=head2 C<merge( $set ... )>

Merge the members of the supplied sets into this set. Any number of sets
may be supplied as arguments.

=head2 C<empty>

Empty a set.

=head2 Operators

=head3 C<complement>

Returns a new set that is the complement of this set. See the comments
about our definition of infinity above.

=head3 C<union( $set ... )>

Return a new set that is the union of this set and all of the supplied
sets.

    $un = $set->union( $other_set );

=head3 C<intersection( $set )>

Return a new set that is the intersection of this set and all the supplied
sets.

    $in = $set->intersection( $other_set );

=head3 C<xor( $set )>

Return a new set that contains all of the members that are in this set
or the supplied set but not both. Can actually handle more than two sets
in which case it returns a set that contains all the members that are in
some of the sets but not all of the sets.

=head3 C<diff( $set )>

Return a set containing all the elements that are in this set but not the
supplied set.

=head2 Tests

=head3 C<is_empty>

Return true if the set is empty.

=head3 C<contains( $number )>

Return true if the specified number is contained in the set.

=head3 C<contains_any($number, $number, $number ...)>

Return true if the set contains any of the specified numbers.

=head3 C<contains_all($number, $number, $number ...)>

Return true if the set contains all of the specified numbers.

=head3 C<contains_all_range( $low, $high )>

Return true if all the numbers in the range C<$low> to C<$high> (inclusive)
are in the set.

=head3 C<cardinality( [ $clip_lo, $clip_hi ] )>

Returns the number of members in the set. If a clipping range is supplied
return the count of members that fall within that inclusive range.

=head3 C<superset( $set )>

Returns true if this set is a superset of the supplied set. A set is
always a superset of itself, or in other words

    $set->superset( $set )
    
returns true.

=head3 C<subset( $set )>

Returns true if this set is a subset of the supplied set. A set is
always a subset of itself, or in other words

    $set->subset( $set )
    
returns true.

=head3 C<equals( $set )>

Returns true if this set is identical to the supplied set.

=head2 Getting set contents

=head3 C<as_array>

Return an array containing all the members of the set in ascending order.

=head3 C<as_string>

Return a string representation of the set.

    my $set = Set::IntSpan::Fast->new();
    $set->add(1, 3, 5, 7, 9);
    $set->add_range(100, 1_000_000);
    print $set->as_string(), "\n";    # prints 1,3,5,7,9,100-1000000

You may optionally supply a hash containing C<sep> and C<range> options:

    print $set->as_string({ sep => ';', range => '*' ), "\n";
        # prints 1;3;5;7;9;100*1000000

=head3 C<iterate_runs( [ $clip_lo, $clip_hi ] )>

Returns an iterator that returns each run of integers in the set in
ascending order. To iterate all the members of the set do something
like this:

    my $iter = $set->iterate_runs();
    while (my ( $from, $to ) = $iter->()) {
        for my $member ($from .. $to) {
            print "$member\n";
        }
    }

If a clipping range is specified only those members that fall within
the range will be returned.

=head2 Constants

The constants C<NEGATIVE_INFINITY> and C<POSITIVE_INFINITY> are exposed. As
noted above these are infinitely smaller than infinity but they're the
best we've got. They're not exported into the caller's namespace so if you
want to use them you'll have to use their fully qualified names:

    $set->add_range(1, Set::IntSpan::Fast::POSITIVE_INFINITY);

=head1 DIAGNOSTICS


=head3 C<< Range list must have an even number of elements >>

The lists of ranges passed to C<add_range> and C<remove_range> consist
of a number of pairs of integers each of which specify the start and end
of a range.

=head3 C<< Range limits must be integers >>

You may only add integers to sets.

=head3 C<< Range limits must be in ascending order >>

When specifying a range in a call to C<add_range> or C<remove_range> the
range bounds must be in ascending order. Multiple ranges don't need to
be in any particular order.

=head3 C<< Value out of range >>

Sets may only contain values in the range C<NEGATIVE_INFINITY> to
C<POSITIVE_INFINITY> inclusive.

=head3 C<< That's very kind of you - but I expect you meant complement() >>

The method that complements a set is called C<complement>.

=head3 C<< I need two sets to compare >>

C<superset> and C<subset> need two sets to compare. They may be called
either as a function:

    $ss = Set::IntSpan::Fast::subset( $s1, $s2 )
    
or as a method:

    $ss = $s1->subset( $s2 );

=head3 C<< Invalid Range String >>

The range string must only contain a comma separated list of ranges, with a hyphen used as the range limit separator. e.g. "1,5,8-12,15-29".

=head1 CONFIGURATION AND ENVIRONMENT

=for author to fill in: A full explanation of any configuration
system(s) used by the module, including the names and locations of any
configuration files, and the meaning of any environment variables or
properties that can be set. These descriptions must also include details
of any configuration language used.

Set::IntSpan::Fast requires no configuration files or environment
variables.

=head1 DEPENDENCIES

    Data::Types
    List::Util

=head1 INCOMPATIBILITIES

Although this module was conceived as a replacement for C<Set::IntSpan>
it isn't a drop-in replacement.

=head1 BUGS AND LIMITATIONS

No bugs have been reported.

Please report any bugs or feature requests to 
C<bug-set-intspan-fast@rt.cpan.org>, or through the web interface at
L<http://rt.cpan.org>.

=head1 AUTHOR

Andy Armstrong C<< <andy@hexten.net> >>

=head1 CREDITS

K. J. Cheetham L<< http://www.shadowcatsystems.co.uk/ >> for
add_from_string, remove_from_string. I butchered his code so any
errors are mine.

=head1 LICENCE AND COPYRIGHT

Copyright (c) 2006-2008, Andy Armstrong C<< <andy@hexten.net> >>. All
rights reserved.

This module is free software; you can redistribute it and/or modify it
under the same terms as Perl itself. See L<perlartistic>.

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
LIABLE TO YOU FOR DAMAGES, INCLUDING ANY GENERAL, SPECIAL,
INCIDENTAL, OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OR
INABILITY TO USE THE SOFTWARE (INCLUDING BUT NOT LIMITED TO LOSS OF
DATA OR DATA BEING RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR
THIRD PARTIES OR A FAILURE OF THE SOFTWARE TO OPERATE WITH ANY OTHER
SOFTWARE), EVEN IF SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE
POSSIBILITY OF SUCH DAMAGES.
