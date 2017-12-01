package GAL::Parser::hgmd_indel;

use strict;
use vars qw($VERSION);
$VERSION = 0.2.0;

use base qw(GAL::Parser);
use GAL::Reader::DelimitedLine;

=head1 NAME

GAL::Parser::hgmd_indel - Parse HGMD_INDEL files

=head1 VERSION

This document describes GAL::Parser::hgmd_indel version 0.2.0

=head1 SYNOPSIS

    my $parser = GAL::Parser::hgmd_indel->new(file => 'hgmd_indel.txt');

    while (my $feature_hash = $parser->next_feature_hash) {
	print $parser->to_gff3($feature_hash) . "\n";
    }

=head1 DESCRIPTION

L<GAL::Parser::hgmd_indel> provides a parser for HGMD_INDEL data.

=head1 Constructor

New L<GAL::Parser::hgmd_indel> objects are created by the class method
new.  Arguments should be passed to the constructor as a list (or
reference) of key value pairs.  All attributes of the
L<GAL::Parser::hgmd_indel> object can be set in the call to new. An
simple example of object creation would look like this:

    my $parser = GAL::Parser::hgmd_indel->new(file => 'hgmd_indel.txt');

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
     Usage   : GAL::Parser::hgmd_indel->new();
     Function: Creates a GAL::Parser::hgmd_indel object;
     Returns : A GAL::Parser::hgmd_indel object
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

  # Skip records if we don't have genomic coordinates or Variant_seq info
  if (lc $record->{genomic_coordinates_hg19} eq 'null') {
    $self->warn('missing_genomic_coordinates', $self->reader->current_line);
    return undef;
  }
  if (lc $record->{hgvs} eq 'null') {
    $self->warn('missing_variant_reference_seqs', $self->reader->current_line);
    return undef;
  }

  my ($seqid, $start_end_strand) = split /:/,
    $record->{genomic_coords_hg19};
  my ($start_end, $strand) = split /\s/, $start_end_strand;
  my ($start, $end) = split /-/, $start_end;
  $end ||= $start;

  my $feature_id = 'HGMD_INDEL_' . sprintf("%06d", $self->counter('INDEL'));
  my $source     = 'HGMD';
  my $score      = '.';
  my $strand     = '+';
  my $phase      = '.';

  my ($variant_seq, $reference_seq);

  my $hgvs = $record->{hgvs};
  my ($rna_id, $hgvs_detail) = split /:/, $hgvs;
  my ($hgvs_locus, $hgvs_detail) =
    $hgvs_detail =~ /([0-9_\-+*\?])((del|dup|ins)([A-Z]+))/;

  my $type;
  my ($del, $ins) = @{$record}{qw(deleted_sequence inserted_sequence)};
  if ($del && $ins) {
    $type = length $del == length $ins ? 'MNP' : 'complex_substitution';
    $reference_seq = $del;
    $variant_seq   = $ins;
  }
  elsif ($del && ! $ins) {
    $variant_seq = '-';
    $reference_seq  = $del;
    $type = 'nucleotide_deletion';
  }
  elsif (! $del && $ins) {
    $variant_seq = $ins;
    $reference_seq  = '-';
    $type = 'nucleotide_insertion';
    $end = $start;
  }
  else {
    $self->warn('no_variant_or_reference_seq', $self->reader->current_line);
    return undef;
  }

  my @dbxrefs;
  push @dbxrefs, ('dbSNP:' . $record->{dbsnp}) if $record->{dbsnp} ne 'null';
  push @dbxrefs, ('PMID:'  . $record->{pmid})  if $record->{pmid}  ne 'null';
  push @dbxrefs, ('ENTREZ:'  . $record->{pmid})
    if $record->{pmid}  =~ /NO ID|ABST|HGOL|LSDB/;

  my @variant_effects;
  push @variant_effects, ('gene_variant 0 gene ' . $record->{gene});
  push @variant_effects, ('transcript_variant 0 transcript ' .
			     $record->{acc_num});
  if ($record->{codon}) {
    push @variant_effects, ('coding_sequence_variant 0 mRNA ' .
			    $record->{acc_num});
  }

  my @seq_context;
  if ($record->{sequence_context_hg19} ne 'null') {
    @seq_context = split /\[.*?\]/, $record->{sequence_context_hg19};
  }

  my $attributes = {Reference_seq => [$reference_seq],
		    Variant_seq   => [$variant_seq],
		    ID            => [$feature_id],
		    hgmd_type     => [$record->{type}],
		    hgmd_class    => [$record->{tag}],
		   };

    $attributes->{Variant_effect} = \@variant_effects if @variant_effects;
    $attributes->{Sequence_context} = \@seq_context if @seq_context;
    $attributes->{Alias} = [$record->{hgvs}] if ($record->{hgvs} ne 'null');
    $attributes->{hgmd_disease} = [$record->{disease}]
      if ($record->{disease} ne 'null');

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

  if (! $self->{reader}) {

    # Change the field names in the next line to define how the hash keys
    # are named for each column in you input file.
    # If you don't pass field_names as an argument to the reader then
    # the reader will return an array rather than a hash.
    my @field_names = qw(type
			 tag
			 gene
			 disease
			 hgvs
			 genomic_coords_hg18
			 genomic_coords_hg19
			 sequence_context_hg18
			 sequence_context_hg19
			 hgmd_sequence_context
			 deleted_sequence
			 inserted_sequence
			 codon
			 author
			 journal
			 vol
			 page
			 year
			 pmid
			 acc_num
			 dbsnp
			 entrezid
		       );

    my $reader = GAL::Reader::DelimitedLine->new(field_names     => \@field_names,
						 field_separator => "\t",
						 comment_pattern => "^\/\/",
						 header_count    => 1);
    $self->{reader} = $reader;
  }
  return $self->{reader};
}

#-----------------------------------------------------------------------------

=head1 DIAGNOSTICS

L<GAL::Parser::hgmd_indel> does not throw any warnings or errors.

=head1 CONFIGURATION AND ENVIRONMENT

L<GAL::Parser::hgmd_indel> requires no configuration files or
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
