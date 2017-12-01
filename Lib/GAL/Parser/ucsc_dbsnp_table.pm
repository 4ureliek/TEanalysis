package GAL::Parser::ucsc_dbsnp_table;

use strict;
use vars qw($VERSION);
$VERSION = 0.2.0;

use base qw(GAL::Parser);
use GAL::Reader::DelimitedLine;

=head1 NAME

GAL::Parser::ucsc_dbsnp_table - Parse UCSC_DBSNP_TABLE files

=head1 VERSION

This document describes GAL::Parser::ucsc_dbsnp_table version 0.2.0

=head1 SYNOPSIS

    my $parser = GAL::Parser::ucsc_dbsnp_table->new(file => 'ucsc_dbsnp_table.txt');

    while (my $feature_hash = $parser->next_feature_hash) {
	print $parser->to_gff3($feature_hash) . "\n";
    }

=head1 DESCRIPTION

L<GAL::Parser::ucsc_dbsnp_table> provides a parser for UCSC_DBSNP_TABLE data.

=head1 Constructor

New L<GAL::Parser::ucsc_dbsnp_table> objects are created by the class method
new.  Arguments should be passed to the constructor as a list (or
reference) of key value pairs.  All attributes of the
L<GAL::Parser::ucsc_dbsnp_table> object can be set in the call to new. An
simple example of object creation would look like this:

    my $parser = GAL::Parser::ucsc_dbsnp_table->new(file => 'ucsc_dbsnp_table.txt');

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
     Usage   : GAL::Parser::ucsc_dbsnp_table->new();
     Function: Creates a GAL::Parser::ucsc_dbsnp_table object;
     Returns : A GAL::Parser::ucsc_dbsnp_table object
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
	my @valid_attributes = qw(file fh); # Set valid class attributes here
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

	# bin chrom chromStart chromEnd name score strand refNCBI refUCSC
	# observed molType class valid avHet avHetSE func locType weight

	my $seqid      = $record->{chrom};
	my $source     = ucsc_dbsnp;
	my $start      = $record->{start};
	my $end        = $record->{end};
	my $score      = $record->{score};
	my $strand     = $record->{strand};
	my $phase      = $record->{phase};

	my %type_map = ('unknown'	 => 'sequence_alteration',
			'single' 	 => 'SNV',
			'in-del' 	 => 'indel',
			'het'    	 => 'sequence_alteration',
			'microsatellite' => 'microsatellite',
			'named'          => 'sequence_alteration',
			'mixed'          => 'sequence_alteration',
			'mnp'            => 'MNP',
			'insertion'      => 'insertion',
			'deletion'       => 'deletion',
			);


	my $type       = 'SNV';
	my $feature_id = join ':', ($seqid, $source, $type, $start);


	my $reference_seq = $record->{refUCSC};
	my $ncbi_ref = $record->{refNCBI};
	if ($reference_seq ne $ncbi_ref) {
	    $self->warn('ncbi_ref_differs_from_ucsc_ref', "NCBI $ncbi_ref UCSC $reference_seq");
	}

	my @variant_seqs = split /\//, $record->{observed};

	if ($strand eq '-') {
	    @variant_seqs = map {$self->revcomp($_)} @variant_seqs;
	    $strand eq '+';
	}
	@variant_seqs = grep {$_ ne $reference_seq} @variant_seqs;


	my $attributes = {Reference_seq => [$reference_seq],
			  Variant_seq   => \@variant_seqs,
			  ID            => [$feature_id],
			 };

	my $feature_data = {feature_id => $feature_id,
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

  # bin chrom chromStart chromEnd name score strand refNCBI refUCSC
  # observed molType class valid avHet avHetSE func locType weight

  if (! $self->{reader}) {
    my @field_names = qw(bin chrom chromStart chromEnd name score strand
			 refNCBI refUCSC observed molType class valid
			 avHet avHetSE func locType weight);
    my $reader = GAL::Reader::DelimitedLine->new(field_names => \@field_names);
    $self->{reader} = $reader;
  }
  return $self->{reader};
}

#-----------------------------------------------------------------------------

=head1 DIAGNOSTICS

L<GAL::Parser::ucsc_dbsnp_table> does not throw any warnings or errors.

=head1 CONFIGURATION AND ENVIRONMENT

L<GAL::Parser::ucsc_dbsnp_table> requires no configuration files or
environment variables.

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
