package GAL::Parser::gff3_fast;

use strict;
use vars qw($VERSION);


$VERSION = 0.2.0;
use base qw(GAL::Parser);
use GAL::Reader::DelimitedLine;

=head1 NAME

GAL::Parser::gff3_fast - Parse GFF3_FAST files

=head1 VERSION

This document describes GAL::Parser::gff3_fast version 0.2.0

=head1 SYNOPSIS

    use GAL::Parser::gff3_fast;
    my $parser = GAL::Parser::gff3_fast->new(file => 'feature.gff3_fast');

    while (my $feature_hash = $parser->next_feature_hash) {
	print $parser->to_gff3_fast($feature_hash) . "\n";
    }

=head1 DESCRIPTION

L<GAL::Parser::gff3_fast> provides GFF3_FAST parsing ability for the GAL library.

=head1 Constructor

New L<GAL::Parser::gff3_fast> objects are created by the class method new.
Arguments should be passed to the constructor as a list (or reference)
of key value pairs.  All attributes of the L<GAL::Parser::gff3_fast> object
can be set in the call to new. An simple example of object creation
would look like this:

    my $parser = GAL::Parser::gff3_fast->new(file => 'feature.gff3_fast');

The constructor recognizes the following parameters which will set the
appropriate attributes:

=over

=item * C<< file => feature_file.txt >>

This optional parameter provides the filename for the file containing
the data to be parsed. While this parameter is optional either it, or
the following fh parameter must be set.

=item * C<< fh => feature_file.txt >>

This optional parameter provides a filehandle to read data from. While
this parameter is optional either it, or the following fh parameter
must be set.

=back

=cut

#-----------------------------------------------------------------------------

=head2 new

     Title   : new
     Usage   : GAL::Parser::gff3_fast->new();
     Function: Creates a GAL::Parser::gff3_fast object;
     Returns : A GAL::Parser::gff3_fast object
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
  my @valid_attributes = qw(file fh); # Set valid class attributes here
  $self->set_attributes($args, @valid_attributes);
  ######################################################################
}

#-----------------------------------------------------------------------------

sub next_feature_hash {

  my $self = shift;

  my $fh = $self->fh;
  my %f;
  my $line;
 LINE:
  while (! $line) {
    $line = <$fh>;
    return undef if ! defined $line;
    chomp $line;
    undef $line if $line =~ /^(#.*|\s*)$/;
  }
  @f{qw(seqid source type start end score strand phase attributes)} = split /\t/, $line;
  my %a = split /=|;/, $f{attributes};
  map {$_ = [split /,/, $_]} values %a;
  @f{qw(feature_id attributes)} = ($a{ID}[0], \%a);
  return wantarray ? %f : \%f;
}

#-----------------------------------------------------------------------------

sub file {

  my ($self, $file) = @_;

  if ($file) {
    $self->throw('unreadable_file', $file) unless -r $file;
    $self->{file} = $file if $file;
    $self->{fh} = undef;
  }
  return $self->{file};
}

#-----------------------------------------------------------------------------

sub fh {

  my ($self, $fh) = @_;

  $self->{fh} = $fh if $fh;
  if (! $self->{fh} && $self->{file}) {
    open(my $FH, '<', $self->{file}) or $self->throw('cant_open_file_for_reading', $self->{file});
    $self->{fh} = $FH;
  }
  return $self->{fh};
}


=head1 DIAGNOSTICS

L<GAL::Parser::gff3_fast> does not throw any warnings or errors.

=head1 CONFIGURATION AND ENVIRONMENT

L<GAL::Parser::gff3_fast> requires no configuration files or environment variables.

=head1 DEPENDENCIES

L<GAL::Parser>
L<GAL::Reader::DelimitedLine>

=head1 INCOMPATIBILITIES

None reported.

=head1 BUGS AND LIMITATIONS

No bugs have been reported.

Please report any bugs or feature requests to:
barry.moore@genetics.utah.edu

=head1 AUTHOR

Barry Moore <barry.moore@genetics.utah.edu>

=head1 LICENCE AND COPYRIGHT

Copyright (c) 2010-2014, Barry Moore <barry.moore@genetics.utah.edu>.  All
rights reserved.

    This module is free software; you can redistribute it and/or
    modify it under the same terms as Perl itself (See LICENSE).

=head1 DISCLAIMER OF WARRANTY

BECAUSE THIS SOFTWARE IS LICENSED FREE OF CHARGE, THERE IS NO WARRANTY
FOR THE SOFTWARE, TO THE EXTENT PERMITTED BY APPLICABLE LAW. EXCEPT
WHEN OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS AND/OR OTHER
PARTIES PROVIDE THE SOFTWARE "AS IS" WITHOUT WARRANTY OF ANY KIND,
EITHER EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE. THE ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE
SOFTWARE IS WITH YOU. SHOULD THE SOFTWARE PROVE DEFECTIVE, YOU ASSUME
THE COST OF ALL NECESSARY SERVICING, REPAIR, OR CORRECTION.

IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING
WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MAY MODIFY AND/OR
REDISTRIBUTE THE SOFTWARE AS PERMITTED BY THE ABOVE LICENCE, BE LIABLE
TO YOU FOR DAMAGES, INCLUDING ANY GENERAL, SPECIAL, INCIDENTAL, OR
CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OR INABILITY TO USE THE
SOFTWARE (INCLUDING BUT NOT LIMITED TO LOSS OF DATA OR DATA BEING
RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR THIRD PARTIES OR A
FAILURE OF THE SOFTWARE TO OPERATE WITH ANY OTHER SOFTWARE), EVEN IF
SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH
DAMAGES.

=cut

1;
