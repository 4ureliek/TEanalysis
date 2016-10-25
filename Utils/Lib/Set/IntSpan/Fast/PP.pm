package Set::IntSpan::Fast::PP;

use warnings;
use strict;
use Carp;
use Data::Types qw(is_int);
use List::Util qw(min max);

=head1 NAME

Set::IntSpan::Fast::PP - Pure Perl implementation.

=head1 VERSION

This document describes Set::IntSpan::Fast::PP version 1.15

=cut

our $VERSION = '1.15';

use constant POSITIVE_INFINITY => 2**31 - 2;
use constant NEGATIVE_INFINITY => -2**31 + 100;

sub new {
  my $class = shift;
  my $self = bless [], $class;
  $self->add_from_string( @_ ) if @_;
  return $self;
}

sub invert {
  my $self = shift;

  if ( $self->is_empty() ) {

    # Empty set
    @$self = ( NEGATIVE_INFINITY, POSITIVE_INFINITY );
  }
  else {

    # Either add or remove infinity from each end. The net
    # effect is always an even number of additions and deletions
    if ( $self->[0] == NEGATIVE_INFINITY ) {
      shift @{$self};
    }
    else {
      unshift @{$self}, NEGATIVE_INFINITY;
    }

    if ( $self->[-1] == POSITIVE_INFINITY ) {
      pop @{$self};
    }
    else {
      push @{$self}, POSITIVE_INFINITY;
    }
  }
}

sub copy {
  my $self  = shift;
  my $class = ref $self;
  my $copy  = $class->new();
  # This might not work for subclasses - in which case they should
  # override copy
  @$copy = @$self;
  return $copy;
}

sub empty { @{ $_[0] } = () }

sub add {
  my $self = shift;
  $self->add_range( $self->_list_to_ranges( @_ ) );
}

sub remove {
  my $self = shift;
  $self->remove_range( $self->_list_to_ranges( @_ ) );
}

sub add_range {
  my $self = shift;

  $self->_iterate_ranges(
    @_,
    sub {
      my ( $from, $to ) = @_;

      my $fpos = $self->_find_pos( $from );
      my $tpos = $self->_find_pos( $to + 1, $fpos );

      $from = $self->[ --$fpos ] if ( $fpos & 1 );
      $to   = $self->[ $tpos++ ] if ( $tpos & 1 );

      splice @$self, $fpos, $tpos - $fpos, ( $from, $to );
    }
  );
}

sub add_from_string {
  my $self = shift;

  my $ctl          = {};
  my $match_number = qr/\s* (-?\d+) \s*/x;
  my $match_single = qr/^ $match_number $/x;
  my $match_range;

  my @to_add = ();

  # Iterate args. Default punctuation spec prepended.
  for my $el ( { sep => qr/,/, range => qr/-/, }, @_ ) {

    # Allow parsing options to be set.
    if ( 'HASH' eq ref $el ) {
      %$ctl = ( %$ctl, %$el );
      for ( values %$ctl ) {
        $_ = quotemeta( $_ ) unless ref $_ eq 'Regexp';
      }
      $match_range = qr/^ $match_number $ctl->{range} $match_number $/x;
    }
    else {
      for my $part ( split $ctl->{sep}, $el ) {
        if ( my ( $start, $end ) = ( $part =~ $match_range ) ) {
          push @to_add, $start, $end;
        }
        elsif ( my ( $el ) = ( $part =~ $match_single ) ) {
          push @to_add, $el, $el;
        }
        else {
          croak "Invalid range string"
           unless $part =~ $match_single;
        }
      }
    }
  }

  $self->add_range( @to_add );
}

sub remove_range {
  my $self = shift;

  $self->invert();
  $self->add_range( @_ );
  $self->invert();
}

sub remove_from_string {
  my $self = shift;

  $self->invert();
  $self->add_from_string( @_ );
  $self->invert();
}

sub merge {
  my $self = shift;

  for my $other ( @_ ) {
    my $iter = $other->iterate_runs();
    while ( my ( $from, $to ) = $iter->() ) {
      $self->add_range( $from, $to );
    }
  }
}

sub compliment {
  croak "That's very kind of you - but I expect you meant complement()";
}

sub complement {
  my $new = shift->copy();
  $new->invert();
  return $new;
}

sub union {
  my $new = shift->copy;
  $new->merge( @_ );
  return $new;
}

sub intersection {
  my $self  = shift;
  my $class = ref $self;
  my $new   = $class->new();
  $new->merge( map { $_->complement() } $self, @_ );
  $new->invert();
  return $new;
}

sub xor {
  my $self = shift;
  return $self->union( @_ )
   ->intersection( $self->intersection( @_ )->complement() );
}

sub diff {
  my $self  = shift;
  my $other = shift;
  return $self->intersection( $other->union( @_ )->complement() );
}

sub is_empty {
  my $self = shift;
  return @$self == 0;
}

*contains = *contains_all;

sub contains_any {
  my $self = shift;

  for my $i ( @_ ) {
    my $pos = $self->_find_pos( $i + 1 );
    return 1 if $pos & 1;
  }

  return;
}

sub contains_all {
  my $self = shift;

  for my $i ( @_ ) {
    my $pos = $self->_find_pos( $i + 1 );
    return unless $pos & 1;
  }

  return 1;
}

sub contains_all_range {
  my ( $self, $lo, $hi ) = @_;

  croak "Range limits must be in ascending order" if $lo > $hi;

  my $pos = $self->_find_pos( $lo + 1 );
  return ( $pos & 1 ) && $hi < $self->[$pos];
}

sub cardinality {
  my $self = shift;

  my $card = 0;
  my $iter = $self->iterate_runs( @_ );
  while ( my ( $from, $to ) = $iter->() ) {
    $card += $to - $from + 1;
  }

  return $card;
}

sub superset {
  my $other = pop;
  return $other->subset( reverse( @_ ) );
}

sub subset {
  my $self = shift;
  my $other = shift || croak "I need two sets to compare";
  return $self->equals( $self->intersection( $other ) );
}

sub equals {
  return unless @_;

  # Array of array refs
  my @edges = @_;
  my $medge = scalar( @edges ) - 1;

  POS: for ( my $pos = 0;; $pos++ ) {
    my $v = $edges[0]->[$pos];
    if ( defined( $v ) ) {
      for ( @edges[ 1 .. $medge ] ) {
        my $vv = $_->[$pos];
        return unless defined( $vv ) && $vv == $v;
      }
    }
    else {
      for ( @edges[ 1 .. $medge ] ) {
        return if defined $_->[$pos];
      }
    }

    last POS unless defined( $v );
  }

  return 1;
}

sub as_array {
  my $self = shift;
  my @ar   = ();
  my $iter = $self->iterate_runs();
  while ( my ( $from, $to ) = $iter->() ) {
    push @ar, ( $from .. $to );
  }

  return @ar;
}

sub as_string {
  my $self = shift;
  my $ctl = { sep => ',', range => '-' };
  %$ctl = ( %$ctl, %{ $_[0] } ) if @_;
  my $iter = $self->iterate_runs();
  my @runs = ();
  while ( my ( $from, $to ) = $iter->() ) {
    push @runs,
     $from == $to ? $from : join( $ctl->{range}, $from, $to );
  }
  return join( $ctl->{sep}, @runs );
}

sub iterate_runs {
  my $self = shift;

  if ( @_ ) {

    # Clipped iterator
    my ( $clip_lo, $clip_hi ) = @_;

    my $pos = $self->_find_pos( $clip_lo ) & ~1;
    my $limit = ( $self->_find_pos( $clip_hi + 1, $pos ) + 1 ) & ~1;

    return sub {
      TRY: {
        return if $pos >= $limit;

        my @r = ( $self->[$pos], $self->[ $pos + 1 ] - 1 );
        $pos += 2;

        # Catch some edge cases
        redo TRY if $r[1] < $clip_lo;
        return   if $r[0] > $clip_hi;

        # Clip to range
        $r[0] = $clip_lo if $r[0] < $clip_lo;
        $r[1] = $clip_hi if $r[1] > $clip_hi;

        return @r;
      }
    };
  }
  else {

    # Unclipped iterator
    my $pos   = 0;
    my $limit = scalar( @$self );

    return sub {
      return if $pos >= $limit;
      my @r = ( $self->[$pos], $self->[ $pos + 1 ] - 1 );
      $pos += 2;
      return @r;
    };
  }

}

sub _list_to_ranges {
  my $self   = shift;
  my @list   = sort { $a <=> $b } @_;
  my @ranges = ();
  my $count  = scalar( @list );
  my $pos    = 0;
  while ( $pos < $count ) {
    my $end = $pos + 1;
    $end++ while $end < $count && $list[$end] <= $list[ $end - 1 ] + 1;
    push @ranges, ( $list[$pos], $list[ $end - 1 ] );
    $pos = $end;
  }

  return @ranges;
}

# Return the index of the first element >= the supplied value. If the
# supplied value is larger than any element in the list the returned
# value will be equal to the size of the list.
sub _find_pos {
  my $self = shift;
  my $val  = shift;
  my $low  = shift || 0;

  my $high = scalar( @$self );

  while ( $low < $high ) {
    my $mid = int( ( $low + $high ) / 2 );
    if ( $val < $self->[$mid] ) {
      $high = $mid;
    }
    elsif ( $val > $self->[$mid] ) {
      $low = $mid + 1;
    }
    else {
      return $mid;
    }
  }

  return $low;
}

sub _iterate_ranges {
  my $self = shift;
  my $cb   = pop;

  my $count = scalar( @_ );

  croak "Range list must have an even number of elements"
   if ( $count % 2 ) != 0;

  for ( my $p = 0; $p < $count; $p += 2 ) {
    my ( $from, $to ) = ( $_[$p], $_[ $p + 1 ] );
    croak "Range limits must be integers"
     unless is_int( $from ) && is_int( $to );
    croak "Range limits must be in ascending order"
     unless $from <= $to;
    croak "Value out of range"
     unless $from >= NEGATIVE_INFINITY && $to <= POSITIVE_INFINITY;

    # Internally we store inclusive/exclusive ranges to
    # simplify comparisons, hence '$to + 1'
    $cb->( $from, $to + 1 );
  }
}

1;
__END__

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
