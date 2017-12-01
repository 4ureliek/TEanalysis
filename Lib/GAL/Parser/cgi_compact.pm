package GAL::Parser::cgi_compact;

use strict;
use vars qw($VERSION);

$VERSION = 0.2.0;
use base qw(GAL::Parser);
use GAL::Reader::DelimitedLine;

=head1 NAME

GAL::Parser::cgi_compact - Parse Complete Genomics compact files

=head1 VERSION

This document describes GAL::Parser::cgi_compact version 0.2.0

=head1 SYNOPSIS

    my $parser = GAL::Parser::cgi_compact->new(file =>
          'CGI-Variations-Compact.csv');

    while (my $feature_hash = $parser->next_feature_hash) {
	print $parser->to_gff3($feature_hash) . "\n";
    }

=head1 DESCRIPTION

L<GAL::Parser::cgi_compact> provides a parser for the
CGI-Variations-Compact.csv format which Complete Genomics provided
with the sequence data for the Feb 5, 2009 sequencing of the anonymous
CEPH (cell-line) genome NA07022.

=head1 Constructor

New L<GAL::Parser::cgi_compact> objects are created by the class
method new.  Arguments should be passed to the constructor as a list
(or reference) of key value pairs.  All attributes of the
L<GAL::Parser::cgi_compact> object can be set in the call to new. An
simple example of object creation would look like this:

    my $parser = GAL::Parser::cgi_compact->new(file =>
          'CGI-Variations-Compact.csv');

The constructor recognizes the following parameters which will set the
appropriate attributes:

The following attributes are inhereted from L<GAL::Parser>.

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
     Usage   : GAL::Parser::cgi_compact->new();
     Function: Creates a GAL::Parser::cgi_compact object;
     Returns : A GAL::Parser::cgi_compact object
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
	my @valid_attributes = qw(); # Set attributes here.
	$self->set_attributes($args, @valid_attributes);
	######################################################################

	# Set the column headers from your incoming data file here
	# These will become the keys in your $record hash reference below.
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

	return undef unless $record->{locus} =~ /^\d+$/;

	# locus,contig,begin,end,vartype1,vartype2,reference,seq1,seq2,totalScore
	# 6,chr1,31843,31844,snp,snp,A,G,G,235
	# 21,chr1,36532,36533,snp,snp,A,G,G,36
	# 23,chr1,36970,36971,snp,snp,G,C,C,109
	# 24,chr1,37154,37155,snp,snp,T,G,G,181
	# 25,chr1,37354,37355,=,snp,C,C,G,73
	# 26,chr1,37623,37624,snp,snp,T,C,C,29
	# 27,chr1,38033,38034,=,snp,A,A,G,54

	# $self->fields([qw(locus contig begin end vartype1 vartype2 reference seq1 seq2 totalScore)]);
	# Fill in the first 8 columns for GFF3
	# See http://www.sequenceontology.org/resources/gff3.html for details.
	my $id         = sprintf 'CG_%09d', $record->{locus};
	my $seqid      = $record->{contig};
	my $source     = 'CGI';

	my %types = map {$_, 1} ($record->{vartype1}, $record->{vartype2});
	my $has_ref_seq;
	$has_ref_seq++ if $types{'='};
	delete $types{'='};

	my ($type) = scalar keys %types == 1 ? keys %types : '';

	my %type_map = (snp		    => 'SNV',
			ins		    => 'nucleotide_insertion',
			del		    => 'nucleotide_deletion',
			inv		    => 'inversion',
		       );

	$type = $type_map{$type} || 'sequence_alteration';

	my $start      = $record->{begin} + 1;
	my $end        = $record->{end};
	my $score      = $record->{totalScore};
	my $strand     = '+';
	my $phase      = '.';

	my $reference_seq = $record->{reference} || '-';
	my %variant_hash  = map {$_ => 1} ($record->{seq1}, $record->{seq2});
	$variant_hash{$reference_seq}++ if $has_ref_seq;
	my @variant_seqs = map {$_ ||= '-'} keys %variant_hash;

	my $zygosity = scalar @variant_seqs > 1 ? 'heterozygous' : 'homozygous';

	my $attributes = {Reference_seq => [$reference_seq],
			  Variant_seq   => \@variant_seqs,
			  Zygosity         => [$zygosity],
			  ID               => [$id],
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
 Returns : A GAL::Reader::DelimitedLine singleton.
 Args    : None

=cut

sub reader {
  my $self = shift;

  if (! $self->{reader}) {
    my @field_names = qw(locus contig begin end vartype1 vartype2 reference
			 seq1 seq2 totalScore);
    my $reader = GAL::Reader::DelimitedLine->new(field_separator   => ',',
						 field_names       => \@field_names,
						 comment_pattern => qr/^[^\d]/,
						);
    $self->{reader} = $reader;
  }
  return $self->{reader};
}

#-----------------------------------------------------------------------------

=head1 DIAGNOSTICS

<GAL::Parser::cgi_compact> does not throw any warnings or errors.

=head1 CONFIGURATION AND ENVIRONMENT

L<GAL::Parser::cgi_compact> requires no configuration files or
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
