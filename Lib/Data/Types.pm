package Data::Types;

use strict;
require Exporter;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

$VERSION = '0.09';

@ISA = qw(Exporter);

@EXPORT_OK = qw(is_whole to_whole is_count to_count is_int to_int is_real
                to_real is_decimal to_decimal is_float to_float is_string
                to_string );

@EXPORT = qw();

%EXPORT_TAGS = (
    all     => \@EXPORT_OK,
    whole   => [qw(is_whole to_whole)],
    count   => [qw(is_count to_count)],
    int     => [qw(is_int to_int)],
    decimal => [qw(is_decimal to_decimal)],
    real    => [qw(is_real to_real)],
    float   => [qw(is_float to_float)],
    string  => [qw(is_string to_string)],
    is      => [qw(is_whole is_int is_real is_decimal is_float is_string)],
    to      => [qw(to_whole to_int to_real to_decimal to_float to_string)],
);

use constant DEF_PRECISION => 5;

################################################################################
# FUNCTIONS                                                                    #
################################################################################

sub is_whole ($) {
    return unless defined $_[0];
    return unless $_[0] =~ /^\d+$/;
    return 1;
}

sub to_whole ($) {
    return unless defined $_[0];
    my ($num) = $_[0] =~ /([+-]?(?:\d+(?:\.\d*)?|\.\d+))/;
    return unless defined $num && $num >= 0;
    sprintf "%.0f", $num;
}

sub is_count ($) {
    return unless $_[0];
    return unless $_[0] =~ /^\d+$/;
    return 1;
}

sub to_count ($) {
    return unless $_[0];
    my ($num) = $_[0] =~ /([+-]?(?:\d+(?:\.\d*)?|\.\d+))/;
    return unless $num && $num > .5;
    sprintf "%.0f", $num;
}

sub is_int ($) {
    return unless defined $_[0] && $_[0] ne '';
    return unless $_[0] =~ /^[+-]?\d+$/;
    return 1;
}

sub to_int ($) {
    return unless defined $_[0] && $_[0] ne '';
    my ($num) = $_[0] =~ /([+-]?(?:\d+(?:\.\d*)?|\.\d+))/;
    return unless defined $num;
    sprintf "%.0f", $num;
}

sub is_decimal ($) {
    return unless defined $_[0] && $_[0] ne '';
    return unless $_[0] =~ /^[+-]?(?:\d+(?:\.\d*)?|\.\d+)$/;
    return 1;
}

sub to_decimal ($;$) {
    return unless defined $_[0] && $_[0] ne '';
    my ($num) = $_[0] =~ /([+-]?(?:\d+(?:\.\d*)?|\.\d+))/;
    return unless defined $num;
    $_[1] ||= DEF_PRECISION;
    sprintf "%.$_[1]f", $num;
}

#sub is_real ($) {
#    return unless defined $_[0] && $_[0] ne '';
#    return unless $_[0] =~ /^[+-]?\d*\.?\d*$/;
#    return 1;
#}

#sub to_real ($) {
#    return unless defined $_[0] && $_[0] ne '';
#    sprintf "%f", $_[0] =~ /([+-]?\d*\.?\d*)/;
#}

# These may need to be separated in the future, in order to identify non-decimal
# real numbers.
*is_real = *is_decimal;
*to_real = *to_decimal;

sub is_float ($) {
    return unless defined $_[0] && $_[0] ne '';
    return unless $_[0] =~ /^([+-]?)(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?$/;
    return 1;
}

sub to_float ($;$) {
    return unless defined $_[0] && $_[0] ne '';
    my ($num) = $_[0] =~ /(([+-]?)(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?)/;
    return unless defined $num;
    my $type = $num =~ /e|E/ ? 'e' : 'f';
    $_[1] ||= DEF_PRECISION;
    sprintf "%.$_[1]$type", $num;
#    sprintf "%g", $_[0] =~ /(([+-]?)(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?)/;
}

sub is_string ($) { defined $_[0] && ! ref $_[0] }

sub to_string ($;$) {
    return unless defined $_[0];
    return $_[1] ? substr("$_[0]", 0, $_[1]) : "$_[0]";
}

1;
__END__

=head1 NAME

Data::Types - Validate and convert data types.

=head1 SYNOPSIS

  use Data::Types qw(:all);

  my $whole = 4.5;
  $whole = to_whole($whole) unless is_whole($whole);

  my $int = 1.2;
  $int = to_int($int) unless is_int($int);

  my $decimal = '1.2foo';
  $decimal = to_decimal($decimal) unless is_decimal($decimal);

  my $real = '1.2foo';
  $real = to_real($real) unless is_real($real);

  my $float = '1.2foo';
  $float = to_float($float) unless is_float($float);

  my $string = [];
  $string = to_string($string) unless is_string($string);

=head1 DESCRIPTION

This module exports a number of functions that are useful for validating and
converting data types. It is intended for use in applications where data types
are more important than they typically are in Perl -- e.g., database
applications.

=head1 EXPORT

No functions are exported by default, though each function may be exported
explicitly (see L<"Functions">, below, for a list of functions available for
export). The following export tags are supported:

=over 4

=item :whole

Exports is_whole() and to_whole().

=item :count

Exports is_count() and to_count().

=item :int

Exports is_int() and to_int().

=item :decimal

Exports is_decimal() and to_decimal().

=item :real

Exports is_real() and to_real().

=item :float

Exports is_float() and to_float().

=item :string

Exports is_string() and to_string().

=item :is

Exports all validation functions: is_whole(), is_int(), is_real(), is_decimal(),
is_float(), and is_string().

=item :to

Exports all conversion functions: to_whole(), to_int(), to_real(), to_decimal(),
to_float(), and to_string().

=item :all

Exports all functions.

=back

=head1 FUNCTIONS

=head2 is_whole

  my $bool = is_whole($val);

Returns true if $val is a whole number (including 0), and false if it is not.
The regular expression used to test the wholeness of $val is C</^\d+$/>.

  my $bool = is_whole(1); # Returns true.
  $bool = is_whole(-1);   # Returns false.
  $bool = is_whole(0);    # Returns true.

=head2 to_whole

  my $whole = to_whole($val);

Converts $val to a whole number and returns it. Numbers will be rounded to the
nearest whole. If $val is a mixture of numbers and letters, to_whole() will
extract the first decimal number it finds and convert that number to a whole
number.

  my $whole = to_whole(10);     # Returns 10.
  $whole = to_whole(0);         # Returns 0.
  $whole = to_whole(.22);       # Returns 0.
  $whole = to_whole(-2);        # Returns undef.
  $whole = to_whole('foo3.56'); # Returns 4.
  $whole = to_whole('foo');     # Returns undef.

=head2 is_count

  my $bool = is_count($val);

Returns true if $val is a counting number (1, 2, 3, ...), and false if it is
not. The regular expression used to test whether $val is a counting number is
C</^\d+$/>.

  my $bool = is_count(1); # Returns true.
  $bool = is_count(-1);   # Returns false.
  $bool = is_count(0);    # Returns false.

=head2 to_count

  my $count = to_count($val);

Converts $val to a counting number and returns it. Numbers will be rounded to
the nearest counting number. Note that since 0 (zero) is not considered a
counting number by this module, it will not be returned. If $val is a mixture
of numbers and letters, to_count() will extract the first decimal number it
finds and convert that number to a counting number.

  my $count = to_count(10);     # Returns 10.
  $count = to_count(0);         # Returns undef.
  $count = to_count(.22);       # Returns undef (rounded down to 0).
  $count = to_count(-2);        # Returns undef.
  $count = to_count('foo3.56'); # Returns 4.
  $count = to_count('foo');     # Returns undef.

=head2 is_int

  my $bool = is_int($val);

Returns true if $val is an integer, and false if it is not. Numbers may be
preceded by a plus or minus sign. The regular expression used to test for an
integer in $val is C</^[+-]?\d+$/>.

  my $bool = is_int(0); # Returns true.
  $bool = is_int(22);   # Returns true.
  $bool = is_int(-22);  # Returns false.
  $bool = is_int(3.2);  # Returns false.

=head2 to_int

  my $int = to_int($val);

Converts $val to an integer. If $val is a decimal number, it will be rounded to
the nearest integer. If $val is a mixture of numbers and letters, to_int() will
extract the first decimal number it finds and convert that number to an integer.

  my $int = to_int(10.5);  # Returns 10.
  $int = to_int(10.51);    # Returns 11.
  $int = to_int(-0.22);    # Returns 0.
  $int = to_int(-6.51);    # Returns 7.
  $int = to_int('foo');    # Returns undef.

=head2 is_decimal

  my $bool = is_decimal($val);

Returns true if $val is a decimal number, and false if it is not. Numbers may be
preceded by a plus or minus sign. The regular expression used to test $val is
C</^[+-]?(?:\d+(?:\.\d*)?|\.\d+)$/>.

  my $bool = is_decimal(10)    # Returns true.
  $bool = is_decimal(10.8)     # Returns true.
  $bool = is_decimal(-33.48)   # Returns true.
  $bool = is_decimal((1.23e99) # Returns false.

=head2 to_decimal

  my $dec = to_decimal($val);
  $dec = to_decimal($val, $precision);

Converts $val to a decimal number. The optional second argument sets the
precision of the number. The default precision is 5. If $val is a mixture of
numbers and letters, to_decimal() will extract the first decimal number it
finds.

  my $dec = to_decimal(0);         # Returns 0.00000.
  $dec = to_decimal(10.5);         # Returns 10.5.
  $dec = to_decimal(10.500009);    # Returns 10.50001.
  $dec = to_decimal(10.500009, 7); # Returns 10.5000090.
  $dec = to_decimal('foo10.3')     # Returns 10.30000.
  $dec = to_decimal('foo-4.9')     # Returns -4.90000.
  $dec = to_decimal('foo')         # Returns undef.

=head2 is_real

  my $bool = is_real($val);

Returns true if $val is a real number, and false if it is not.

B<Note:> This function is currently equivalent to is_decimal(), since this
module cannot identify non-decimal real numbers (e.g., irrational numbers). This
implementation may change in the future.

=head2 to_real

  my $real = to_real($val);
  $real = to_real($val, $precision);

Converts $val to a real number.

B<Note:> Currently, this function is the equivalent of to_decimal(), since this
module cannot identify non-decimal real numbers (e.g., irrational numbers). This
implementation may change in the future.

=head2 is_float

  my $bool = is_real($val);

Returns true if $val is a float, and false if it is not. The regular expression
used to test $val is C</^([+-]?)(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?$/>.

  my $bool = is_real(30);   # Returns true.
  $bool = is_real(1.23e99); # Returns true.
  $bool = is_real('foo');   # Returns false.

=head2 to_float

  my $dec = to_float($val);
  $dec = to_float($val, $precision);

Converts $val to a float. The optional second argument sets the precision of the
number. The default precision is 5. If $val is a mixture of numbers and letters,
to_float() will extract the first float it finds.

  my $float = to_float(1.23);          # Returns 1.23000.
  $float = to_float(1.23e99);          # Returns 1.23000e+99.
  $float = to_float(1.23e99, 1);       # Returns 1.2e+99.
  $float = to_float('foo-1.23');       # Returns -1.23000.
  $float = to_float('ick_1.23e99foo'); # Returns 1.23000e+99.

=head2 is_string

  my $bool = is_string($val);

Returns true if $val is a string, and false if it is not. All defined
non-references are considered strings.

  my $bool = is_string('foo'); # Returns true.
  $bool = is_string(20001);    # Returns true.
  $bool = is_string([]);       # Returns false.
  $bool = is_string(undef);    # Returns false.

=head2 to_string

  my $string = to_string($val);
  $string = to_string($val, $length);

Converts $val into a string. If $val is a reference, the string value of the
reference will be returned. Such a value may be a memory address, or some other
value, if the stringification operator has been overridden for the object stored
in $val. If the optional second argument $length is passed, to_string() will
truncate the string to that length. If $length is 0 (zero), it will not limit
the length of the return string. If $val is undefined, to_string() will return
undef.

  my $string = to_string('foo');   # Returns 'foo'.
  $string = to_string([]);         # Returns 'ARRAY(0x101bec14)'.
  $string = to_string(undef);      # Returns undef.
  $string = to_string('hello', 4); # Returns 'hell'.

=head1 SUPPORT

This module is stored in an open L<GitHub
repository|http://github.com/theory/data-types/>. Feel free to fork and
contribute!

Please file bug reports via L<GitHub
Issues|http://github.com/theory/data-types/issues/> or by sending mail to
L<bug-Data-Types.cpan.org|mailto:bug-Data-Types.cpan.org>.

Patches against Class::Meta are welcome. Please send bug reports to
<bug-data-types@rt.cpan.org>.

=head1 AUTHOR

David E. Wheeler <david@justatheory.com>

=head1 SEE ALSO

L<perlfaq4|perlfaq4/"How do I determine whether a scalar is anumber/whole/integer/float?">
lists the most of the regular expressions used to identify the different numeric
types used in this module.

L<String::Checker|String::Checker> also does some data type validation.

L<String::Scanf|String::Scanf> reimplements the C C<sscanf()> function in
perl, and also does data type validation and conversion.

L<Regexp::Common|Regexp::Common> contains many useful common regular expressions
(surprise!), including some that can be used to identify data types.

Arthur Bergman's L<types|types> pragma, offers compile-time data types for
Perl 5.8.0. The data types include int, float, and string. I highly recommend
using this prgrma for fast, static data types.

=head1 COPYRIGHT AND LICENSE

Copyright (c) 2002-2011, David E. Wheeler. Some Rights Reserved.

This module is free software; you can redistribute it and/or modify it under the
same terms as Perl itself.

=cut
