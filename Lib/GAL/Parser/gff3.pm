package GAL::Parser::gff3;

use strict;
use vars qw($VERSION);


$VERSION = 0.2.0;
use base qw(GAL::Parser);
use GAL::Reader::DelimitedLine;

=head1 NAME

GAL::Parser::gff3 - Parse GFF3 files

=head1 VERSION

This document describes GAL::Parser::gff3 version 0.2.0

=head1 SYNOPSIS

    use GAL::Parser::gff3;
    my $parser = GAL::Parser::gff3->new(file => 'feature.gff3');

    while (my $feature_hash = $parser->next_feature_hash) {
	print $parser->to_gff3($feature_hash) . "\n";
    }

=head1 DESCRIPTION

L<GAL::Parser::gff3> provides GFF3 parsing ability for the GAL library.

=head1 Constructor

New L<GAL::Parser::gff3> objects are created by the class method new.
Arguments should be passed to the constructor as a list (or reference)
of key value pairs.  All attributes of the L<GAL::Parser::gff3> object
can be set in the call to new. An simple example of object creation
would look like this:

    my $parser = GAL::Parser::gff3->new(file => 'feature.gff3');

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
     Usage   : GAL::Parser::gff3->new();
     Function: Creates a GAL::Parser::gff3 object;
     Returns : A GAL::Parser::gff3 object
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

=head2 parse_record

 Title   : parse_record
 Usage   : $a = $self->parse_record();
 Function: Parse the data from a record.
 Returns : A hash (or reference) of feature data.
 Args    : A hash reference of feature data from the reader

=cut

sub parse_record {
	my ($self, $record) = @_;

	my $attributes = $self->parse_attributes($record->{attributes});

	my %feature_hash = (seqid      => $record->{seqid},
			    source     => $record->{source},
			    type       => $record->{type},
			    start      => $record->{start},
			    end        => $record->{end},
			    score      => $record->{score},
			    strand     => $record->{strand},
			    phase      => $record->{phase},
			    attributes => $attributes,
			   );

	if (exists $attributes->{ID}         &&
	    ref $attributes->{ID} eq 'ARRAY' &&
	    $attributes->{ID}[0]) {
	  $feature_hash{feature_id} = $attributes->{ID}[0];
	}
	else {
	  $feature_hash{feature_id} =  $self->create_unique_id(\%feature_hash);
	  $attributes->{ID}[0] = $feature_hash{feature_id};
	}

	return wantarray ? %feature_hash : \%feature_hash;
}

#-----------------------------------------------------------------------------

=head2 parse_attributes

 Title   : parse_attributes
 Usage   : $a = $self->parse_attributes($attrb_text);
 Function: Parse the attributes from column 9 in a GFF3 file.
 Returns : A hash (or reference) of tag/value(s) pairs.  Values are always
           array references even if they are only a single value.
           For example %attributes = (tag => [$value]);
 Args    : The text from column 9 in a GFF3 file.

=cut

sub parse_attributes {
	my ($self, $attrb_text) = @_;

	my @attrbs = split /\s*;\s*/, $attrb_text;
	my %attrb_hash;
	for my $attrb (@attrbs) {
		my ($tag, $value_text) = split /=/, $attrb;
		my @values = split /,/, $value_text;
		push @{$attrb_hash{$tag}}, @values;
	}
	return wantarray ? %attrb_hash : \%attrb_hash;
}

#-----------------------------------------------------------------------------

=head2 reader

 Title   : reader
 Usage   : $a = $self->reader
 Function: Return the reader object.
 Returns : A GAL::Reader::DelimitedLine singleton.
 Args    : None

=cut

sub reader {
  my $self = shift;

  if (! $self->{reader}) {
    my @field_names = qw(seqid source type start end score strand phase
			 attributes);
    my $reader = GAL::Reader::DelimitedLine->new(field_separator => "\t",
						 field_names     => \@field_names,
						 end_of_data     => qr/\#\#FASTA/);
    $self->{reader} = $reader;
  }
  return $self->{reader};
}

#-----------------------------------------------------------------------------

=head1 DIAGNOSTICS

L<GAL::Parser::gff3> does not throw any warnings or errors.

=head1 CONFIGURATION AND ENVIRONMENT

L<GAL::Parser::gff3> requires no configuration files or environment variables.

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
