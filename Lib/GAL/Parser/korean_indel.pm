package GAL::Parser::korean_indel;

use strict;
use vars qw($VERSION);


$VERSION = 0.2.0;
use base qw(GAL::Parser);
use GAL::Reader::DelimitedLine;

=head1 NAME

GAL::Parser::korean_indel - Parse indel files from the first Korean
genome

=head1 VERSION

This document describes GAL::Parser::korean_indel version 0.2.0

=head1 SYNOPSIS

    my $parser = GAL::Parser::korean_indel->new(file => 'korean_indel.gff');

    while (my $feature_hash = $parser->next_feature_hash) {
	print $parser->to_gff3($feature_hash) . "\n";
    }

=head1 DESCRIPTION

L<GAL::Parser::korean_indel> provides a parser for indel files from the
first Korean genome published by Ahn, et al. 2009
(http://www.ncbi.nlm.nih.gov/pubmed/19470904).

=head1 Constructor

New L<GAL::Parser::korean_indel> objects are created by the class
method new.  Arguments should be passed to the constructor as a list
(or reference) of key value pairs.  All attributes of the
L<GAL::Parser::korean_indel> object can be set in the call to
new. An simple example of object creation would look like this:

    my $parser = GAL::Parser::korean_indel->new(file => 'korean_indel.gff');

The constructor recognizes the following parameters which will set the
appropriate attributes:

=over 4

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
     Usage   : GAL::Parser::korean_indel->new();
     Function: Creates a GAL::Parser::korean_indel object;
     Returns : A GAL::Parser::korean_indel object
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
	my @valid_attributes = qw(); # Set valid class attributes here.
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

	# Create the attributes hash

	# N/A chr1 713662 13 -2:AG    AGGGAGAGAGAAAGGAAGAGACGATGAGAGAC AGAGAGAAGGAGAGAGAAAGTACAAAAGAACG	HET Non_genic Other
	# N/A chr1 714130 49  5:GAATG TGGAACGCACTCGAATGGAATGGAACGGACAT GAATGGAATGGAATGGAACGGACACGAATGGA	HET Non_genic Other
	# N/A chr1 715647 78 -2:AG    TAATGGAATGGACTTGAATGGAATAGAATGGA AGAGACTCGAATGGAATGGAATGCAATGGAAT	HET Non_genic Other
	# N/A chr1 780560 13  2:AT    TACGGGTGTATCTGTGTATTGTGTATGCACAC ACGAGCATATGTGTACATGAATTTGTATTGCA	HET Non_genic Other
	# N/A chr1 780622 27 -2:TA    CACATGTGTTTAATGCGAACACGTGTCATGTG TATGTGTTCACATGCATGTGTGTCTGTGTACT	HET Non_genic Other
	# N/A chr1 794457 9  -1:A     TGCTGTGACAAAAAAGCAGGGAAAGGGAATTT AAAAAAAAAAAAGCAAACAACAACAACAAAAA	HET Non_genic Other
	# N/A chr1 805541 10 -1:G     GTTAGCTGTGTTTTTTGTTGTTGTTGTTTTTT GGGGTTTTTTTTGTATAACATTATGTTAAGGT	HET Non_genic Other
	# N/A chr1 806031 30 -1:T     CAGTGTAGCCATCTGGTCCAGGCTTTTCTTTG TTGCTGGGTTTTTTATTACTGATGCAATCTTC	HET Non_genic Other
	# N/A chr1 806276 26 -1:T     TTATTTTTGAGTTTGGTAATTTGAGTATTCCC TTTTTTTCTTAGTCAATCTAGATAAAATTTTG	HET Non_genic Other
	# N/A chr1 806790 32  1:C     TGTATCAACATTTGTTGTGTTCTCATAAACTT TGTAATACATGGAGATTTCTGGTCCACATATG	HET Non_genic Other

	# $self->fields([qw(transcript_id chromosome location total_reads seq context1
	#                   context2 genotype gene_name gene_part)]);

	my ($seq_size, $seq) =  split /:/, $record->{seq};

	my $seqid      = $record->{chromosome};
	my $source     = 'Illumina_GA';
	my $type       = $seq_size < 0 ? 'nucleotide_deletion' : 'nucleotide_insertion';
	my $start      = $seq_size < 0 ? $record->{location} : $record->{location} - 1;
	my $end        = $seq_size < 0 ? $record->{location} - $seq_size - 1 : $record->{location} - 1;
	my $score      = '.';
	my $strand     = '+';
	my $phase      = '.';

	my $id = join ":", ($source, $type, $seqid, $start);

	my $reference_seq = $seq_size < 0 ? $seq : '-';
	my @variant_seqs;
	if ($record->{zygosity} eq 'HET') {
		push @variant_seqs, ($seq_size < 0 ? $seq : '-')
	}
	else {
		push @variant_seqs, ($seq_size < 0 ? '-'     : $seq);
	}

	my $total_reads = $record->{total_reads};

	my $zygosity = $record->{zygosity} eq 'HET' ? 'heterozygous' : 'homozygous';

	my $intersected_gene;
	$intersected_gene = $record->{gene_name} ne 'Non_genic' ? 'gene:HGNC:' . $record->{gene_name} : undef;

	#perl -lane 'print $F[9] unless $F[9] eq "Other"' KOREF-solexa-indel-X30_d3D50E20.gff | sort | uniq -c | sort -nr
	#  127516 Intron
	#     319 3UTR
	#      49 CDS
	#      27 5UTR

	my %type_map = ('Intron' => 'intron',
			'3UTR'   => 'three_prime_UTR',
			'CDS'    => 'CDS',
			'5UTR'   => 'five_prime_UTR',
			'Other'  =>  undef,
		       );

	my $intersected_gene_part;
	$intersected_gene_part = $record->{gene_part};

	my @intersected_features;

	push @intersected_features, $intersected_gene      if $intersected_gene;
	push @intersected_features, $intersected_gene_part if $intersected_gene_part;

	my $attributes = {Reference_seq => [$reference_seq],
			  Variant_seq   => \@variant_seqs,
			  Total_reads   => [$total_reads],
			  Zygosity      => [$zygosity],
			  ID            => [$id],
			 };

	$attributes->{Intersected_feature} = \@intersected_features if scalar @intersected_features;

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
    my @field_names = qw(transcript_id chromosome location total_reads seq context1
			  context2 zygosity gene_name gene_part);
    my $reader = GAL::Reader::DelimitedLine->new(field_names     => \@field_names);
    $self->{reader} = $reader;
  }
  return $self->{reader};
}

#-----------------------------------------------------------------------------

=head1 DIAGNOSTICS

L<GAL::Parser::korean_indel> does not throw any warnings or errors.

=head1 CONFIGURATION AND ENVIRONMENT

L<GAL::Parser::korean_indel> requires no configuration files or
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
