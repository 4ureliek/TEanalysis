package GAL::Parser::cuffdiff_gene_exp;

use strict;
use warnings;
use vars qw($VERSION);


$VERSION = 0.2.0;
use base qw(GAL::Parser);
use GAL::Reader::DelimitedLine;

=head1 NAME

GAL::Parser::cuffdiff_gene_exp - Parse CUFFDIFF_GENE_EXP files

=head1 VERSION

This document describes GAL::Parser::cuffdiff_gene_exp version 0.2.0

=head1 SYNOPSIS

    my $parser = GAL::Parser::cuffdiff_gene_exp->new(file => 'gene_exp.diff');

    while (my $feature_hash = $parser->next_feature_hash) {
	print $parser->to_gff3($feature_hash) . "\n";
    }

=head1 DESCRIPTION

<GAL::Parser::cuffdiff_gene_exp> provides a parser for cuffdiff
gene_exp.diff files.

=head1 Constructor

New <GAL::Parser::cuffdiff_gene_exp> objects are created by the class method
new.  Arguments should be passed to the constructor as a list (or
reference) of key value pairs.  All attributes of the
<GAL::Parser::cuffdiff_gene_exp> object can be set in the call to new. An
simple example of object creation would look like this:

    my $parser = GAL::Parser::cuffdiff_gene_exp->new(file => 'gene_exp.diff');

The constructor recognizes the following parameters which will set the
appropriate attributes:

=over

=item * C<< file => feature_file.txt >>

This optional parameter provides the filename for the file containing
the data to be parsed. While this parameter is optional either it, or
the following fh parameter must be set.

=item * C<< fh => $fh >>

This optional parameter provides a filehandle to read data from. While
this parameter is optional either it, or the file option should be sued
=back

=cut

#-----------------------------------------------------------------------------

=head2 new

     Title   : new
     Usage   : GAL::Parser::cuffdiff_gene_exp->new();
     Function: Creates a GAL::Parser::cuffdiff_gene_exp object;
     Returns : A GAL::Parser::cuffdiff_gene_exp object
     Args    : See the attributes described above.

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
	my @valid_attributes = qw(); # Set valid class attributes here
	$self->set_attributes($args, @valid_attributes);
	######################################################################
}

#-----------------------------------------------------------------------------

=head2 parse_record

 Title   : parse_record
 Usage   : $a = $self->parse_record();
 Function: Parse the data from a record.
 Returns : A hash ref needed by Feature.pm to create a Feature object
 Args    : A hash ref of fields that this sub can understand (In this case GFF3).

=cut

sub parse_record {
	my ($self, $record) = @_;

	# $record is a hash reference that contains the keys assigned
	# in @field_names array in the reader method below.

	# Fill in the first 8 columns for GFF3
	# See http://www.sequenceontology.org/resources/gff3.html for details.

	# test_id gene_id gene locus sample_1 sample_2 status value_1 value_2 log2_fc test_stat p_value q_value significant
	# XLOC_000002 XLOC_000002 MIR1302-2 1:30365-30503 CTRL CAPER NOTEST 0 0 0 0 1 1 no


	my ($seqid, $start, $end) = split /:|-/, $record->{locus};
	my $id         = $record->{gene};
	my $source     = 'cuffdiff';
	my $type       = 'gene';
	my $score      = $record->{test_stat};
	my $strand     = '.';
	my $phase      = '.';

	# Create the attribute hash reference.  Note that all values
	# are array references - even those that could only ever have
	# one value.  This is for consistency in the interface.
	# Suggested keys include (from the GFF3 spec), but are not
	# limited to: ID, Name, Alias, Parent, Target, Gap,
	# Derives_from, Note, Dbxref and Ontology_term. Note that
	# attribute names are case sensitive. "Parent" is not the same
	# as "parent". All attributes that begin with an uppercase
	# letter are reserved for later use. Attributes that begin
	# with a lowercase letter can be used freely by applications.

	my $attributes = {ID              =>  [$id],
			  cd_gene_id      =>  [$record->{gene_id}],
			  cd_gene         =>  [$record->{gene}],
			  cd_sample1      =>  [$record->{sample_1}],
			  cd_sample2      =>  [$record->{sample_2}],
			  cd_status       =>  [$record->{status}],
			  cd_value1       =>  [$record->{value_1}],
			  cd_value2       =>  [$record->{value_2}],
			  cd_test_stat    =>  [$record->{test_stat}],
			  cd_p_value      =>  [$record->{p_value}],
			  cd_q_value      =>  [$record->{q_value}],
			  cd_significant  =>  [$record->{significant}],
			 };

	my $feature_data = {feature_id => $id,
			    seqid      => $seqid,
			    source     => $source,
			    type       => $type,
			    start      => $start,
			    end        => $end,
			    score      => $score,
			    strand     => $strand,
			    phase      => $phase,
			    attributes => $attributes,
			   };

	return $feature_data;
}

#-----------------------------------------------------------------------------

=head2 reader

 Title   : reader
 Usage   : $a = $self->reader
 Function: Return the reader object.
 Returns : A L<GAL::Reader::DelimitedLine> singleton.
 Args    : None

=cut

sub reader {
  my $self = shift;

  if (! $self->{reader}) {
      my @field_names = qw(test_id gene_id gene locus sample_1 sample_2
                           status value_1 value_2 log2_fc test_stat p_value
                           q_value significant);

    my $reader = GAL::Reader::DelimitedLine->new(field_names  => \@field_names,
						 header_count => 1);
    $self->{reader} = $reader;
  }
  return $self->{reader};
}

#-----------------------------------------------------------------------------

=head1 DIAGNOSTICS

L<GAL::Parser::cuffdiff_gene_exp> does not throw any warnings or errors.

=head1 CONFIGURATION AND ENVIRONMENT

L<GAL::Parser::cuffdiff_gene_exp> requires no configuration files or environment variables.

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
