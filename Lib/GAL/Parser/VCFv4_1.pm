package GAL::Parser::VCFv4_1;

use strict;
use vars qw($VERSION);
$VERSION = 0.2.0;

use base qw(GAL::Parser);
use GAL::Reader::DelimitedLine;

=head1 NAME

GAL::Parser::VCFv4_1 - A GAL Parser for VCF version 4.1 files

=head1 VERSION

This document describes GAL::Parser::VCFv4_1 version 0.2.0

=head1 SYNOPSIS

     use GAL::Parser::VCFv4_1;

=for author to fill in:
     Brief code example(s) here showing commonest usage(s).
     This section will be as far as many users bother reading
     so make it as educational and exemplary as possible.

=head1 DESCRIPTION

=for author to fill in:
     Write a full description of the module and its features here.
     Use subsections (=head2, =head3) as appropriate.

=head1 METHODS

=cut

#-----------------------------------------------------------------------------

=head2 new

     Title   : new
     Usage   : GAL::Parser::VCFv4_1->new();
     Function: Creates a GAL::Parser::VCFv4_1 object;
     Returns : A GAL::Parser::VCFv4_1 object
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
	my @valid_attributes = qw(file fh header_count); # Set valid class attributes here
	$self->set_attributes($args, @valid_attributes);
	######################################################################
}

#-----------------------------------------------------------------------------

sub header_count {
    my ($self, $header_count) = @_;

    if ($header_count) {
	$self->{header_count} = $header_count;
    }
    return $self->{header_count};
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
  my ($self, $data) = @_;

  my %record;

  # We didn't get $data as a hash because we didn't give field names to the reader.
  # This allows us to get the data as an array and deal with multiple genotypes per line.
  @record{qw(chrom pos id ref alt qual filter info format)} = splice(@{$data}, 0, 9);

  my $seqid         = $self->get_seqid(\%record);
  my $source        = $self->get_source(\%record);
  my $type          = $self->get_type(\%record);
  my $start         = $self->get_start(\%record);
  my $end           = $self->get_end(\%record);;
  my $score         = $self->get_score(\%record);
  my $strand        = '+';
  my $phase         = '.';

  my $vcf_info_attributes   = $self->parse_vcf_info(\%record);
  my $vcf_format_attributes = $self->parse_vcf_format(\%record, $data);

  my $attributes = $self->build_attributes(\%record, $vcf_info_attributes,
					   $vcf_format_attributes);

  my $feature = {seqid      => $seqid,
		 source     => $source,
		 type       => $type,
		 start      => $start,
		 end        => $end,
		 score      => $score,
		 strand     => $strand,
		 phase      => $phase,
		 attributes => $attributes,
		};

  $self->validate_feature_id($feature);

  return $feature;
}

#-----------------------------------------------------------------------------

sub build_attributes {

  my ($self, $record, $vcf_info_attributes, $vcf_format_attributes) = @_;

  #map {$_ = 'vcf' . $_} keys %{$vcf_info_attributes};

  my %attributes = (%{$vcf_info_attributes}, %{$vcf_format_attributes});
  $attributes{Reference_seq} = [$record->{ref}];
  $attributes{Variant_seq} = [split /,/, $record->{alt}];

  return wantarray ? %attributes : \%attributes;

}

#-----------------------------------------------------------------------------

sub get_seqid {

  my ($self, $record) = @_;

  my $seqid = $record->{chrom};

  return $seqid;

}

#-----------------------------------------------------------------------------

sub get_source {

  my ($self, $record) = @_;

  return '.';
}

#-----------------------------------------------------------------------------

sub parse_vcf_info {

  my ($self, $record) = @_;

  my %info;
  my @pairs = split /;/, $record->{info};
  for my $pair (@pairs) {
    my ($key, $value) = split /=/, $pair;
    $value = defined $value ? $value : '';
    ($key, $value) = $self->parse_vcf_info_key_value($record, $key, $value);
    push @{$info{$key}}, @{$value};
  }
  return wantarray ? %info : \%info;
}

#-----------------------------------------------------------------------------

sub parse_vcf_info_key_value {

  my ($self, $key, $value) = @_;

  $key = 'vcf_' . $key;
  my @values = split /,/, $value;

  return $key, \@values;
}

#-----------------------------------------------------------------------------

sub get_type {
  my ($self, $record) = @_;
  return 'SNV';
}

#-----------------------------------------------------------------------------

sub get_start {
  my ($self, $record) = @_;
  return $record->{pos};
}

#-----------------------------------------------------------------------------

sub get_end {
  my ($self, $record) = @_;
  return $record->{pos};
}

#-----------------------------------------------------------------------------

sub get_score {
  my ($self, $record) = @_;
  return $record->{qual};
}

#-----------------------------------------------------------------------------

sub get_feature_id {
  my ($self, $record) = @_;
  my $feature_id = $record->{id};
  return $feature_id;
}

#-----------------------------------------------------------------------------

sub validate_feature_id {
  my ($self, $feature) = @_;

  my $feature_id = $feature->{attributes}{ID};
  if (! $feature_id ||
      $feature_id eq '.') {
    $feature_id = join ':', @{$feature}{qw(seqid type start end)};
    $feature->{attributes}{ID} = [$feature_id];
    $feature->{feature_id} = $feature_id;
  }
}

#-----------------------------------------------------------------------------

sub get_reference_seq {
  my ($self, $record) = @_;
  return uc $record->{ref};
}

#-----------------------------------------------------------------------------

sub get_variant_seqs {

  my ($self, $record) = @_;

  my @variant_seqs = split /,/, $record->{alt};
  map {$_ = uc $_} @variant_seqs;
  return wantarray ? @variant_seqs : \@variant_seqs;
}

#-----------------------------------------------------------------------------

sub _parse_vcf_metadata {

  my ($reader, $line) = @_;

  ##fileformat=VCFv4_1
  ##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
  ##INFO=<ID=HM2,Number=0,Type=Flag,Description="HapMap2 membership">
  ##INFO=<ID=HM3,Number=0,Type=Flag,Description="HapMap3 membership">
  ##reference=human_b36_both.fasta
  ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
  ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
  ##FORMAT=<ID=CB,Number=1,Type=String,Description="Called by S(Sanger), M(UMich), B(BI)">
  ##rsIDs=dbSNP b129 mapped to NCBI 36.3, August 10, 2009
  ##INFO=<ID=AC,Number=.,Type=Integer,Description="Allele count in genotypes">
  ##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">

  my ($key, $values_text) = split /=/, $line;

  my @values = split ',', $values_text;

  my %data;

  if ($key eq 'INFO') {

    # ##INFO=<ID=ID,Number=number,Type=type,Description=”description”>
    #
    # Possible Types for INFO fields are: Integer, Float, Flag, Character, and String.
    #
    # The Number entry is an Integer that describes the number of values
    # that can be included with the INFO field. For example, if the INFO
    # field contains a single number, then this value should be 1; if the
    # INFO field describes a pair of numbers, then this value should be 2
    # and so on. If the field has one value per alternate allele then this
    # value should be 'A'; if the field has one value for each possible
    # genotype (more relevant to the FORMAT tags) then this value should be
    # 'G'.  If the number of possible values varies, is unknown, or is
    # unbounded, then this value should be '.'. The 'Flag' type indicates
    # that the INFO field does not contain a Value entry, and hence the
    # Number should be 0 in this case. The Description value must be
    # surrounded by double-quotes. Double-quote character can be escaped
    # with backslash (\") and backslash as \\.

    %data = ($key => \@values);
  }
  elsif ($key eq 'FILTER') {
    # FILTERs that have been applied to the data should be described as follows:
    #
    # ##FILTER=<ID=ID,Description=”description”>
    #
    # Likewise, Genotype fields specified in the FORMAT field should be described as follows:
    #
    %data = ($key => \@values);
  }
  elsif ($key eq 'FORMAT') {

    # ##FORMAT=<ID=ID,Number=number,Type=type,Description=”description”>
    #
    # Possible Types for FORMAT fields are: Integer, Float, Character,
    # and String (this field is otherwise defined precisely as the
    # INFO field).
    %data = ($key => \@values);
  }
  elsif ($key eq 'ALT') {

    # Symbolic alternate alleles for imprecise structural variants:
    #
    # ##ALT=<ID=type,Description=description>
    #
    # The ID field indicates the type of structural variant, and can
    # be a colon-separated list of types and subtypes. ID values are
    # case sensitive strings and may not contain whitespace or angle
    # brackets. The first level type must be one of the following:
    #
    #     DEL Deletion relative to the reference
    #     INS Insertion of novel sequence relative to the reference
    #     DUP Region of elevated copy number relative to the reference
    #     INV Inversion of reference sequence
    #     CNV Copy number variable region (may be both deletion and duplication)
    #
    # The CNV category should not be used when a more specific category can be applied. Reserved subtypes include:
    #
    #     DUP:TANDEM Tandem duplication
    #     DEL:ME Deletion of mobile element relative to the reference
    #     INS:ME Insertion of a mobile element relative to the reference
    #
    # In addition, it is highly recommended (but not required) that
    # the header include tags describing the reference and contigs
    # backing the data contained in the file.  These tags are based on
    # the SQ field from the SAM spec; all tags are optional (see the
    # VCF example above).
    #
    # Breakpoint assemblies for structural variations may use an external file:
    #
    %data = ($key => \@values);
  }
  elsif ($key eq 'assembly') {
    # ##assembly=url
    #
    # The URL field specifies the location of a fasta file containing
    # breakpoint assemblies referenced in the VCF records for
    # structural variants via the BKPTID INFO key.
    #
    %data = ($key => \@values);
  }
  elsif ($key eq 'contig') {

    # As with chromosomal sequences it is highly recommended (but not
    # required) that the header include tags describing the contigs
    # referred to in the VCF file. This furthermore allows these
    # contigs to come from different files. The format is identical to
    # that of a reference sequence, but with an additional URL tag to
    # indicate where that sequence can be found. For example:.
    #
    # ##contig=<ID=ctg1,URL=ftp://somewhere.org/assembly.fa,...>
    #
    %data = ($key => \@values);
  }
  elsif ($key eq 'SAMPLE') {
    # As explained below it is possible to define sample to genome mappings:
    #
    # ##SAMPLE=<ID=S_ID,Genomes=G1_ID;G2_ID; ...;GK_ID,Mixture=N1;N2; ...;NK,Description=S1;S2; ...; SK >
    #
    %data = ($key => \@values);
  }
  elsif ($key eq 'PEDIGREE' || $key eq 'pedigreeDB') {

    # As well as derivation relationships between genomes using the
    # following syntax:
    #
    # ##PEDIGREE=<Name_0=G0-ID,Name_1=G1-ID,...,Name_N=GN-ID>
    #
    # or a link to a database:
    #
    # ##pedigreeDB=<url>
    %data = ($key => \@values);
  }
  else {
    %data = ($key => \@values);
  }
  return wantarray ? %data : \%data;
}

#-----------------------------------------------------------------------------

sub _parse_vcf_header {
  my ($reader, $line) = @_;

  #CHROM POS ID REF ALT QUAL FILTER INFO FORMAT SAMPLE1

  my @samples = split /\s+/, $line;
  my @headers = splice(@samples, 0, 9);
  my %header = (header  => \@headers,
		samples => \@samples,
	       );

  return wantarray ? %header : \%header;
}

#-----------------------------------------------------------------------------

sub parse_vcf_format {

  my ($self, $record, $data) = @_;



  my %format_attributes = (Individual_format => [$record->{format}],
			   Individual_data   => $data,
			  );

  return \%format_attributes;

  # Keep this code here for use in the GVF writer;
  #my @features;
  #for my $individual_id (sort keys %individual_data) {
  #
  #  my @individual_values = split /:/, $individual_data{$individual_id};
  #  my %individual_hash;
  #  map {$individual_hash{$_} = shift @individual_values} @format_order;
  #
  #  my ($var1, $phased, $var2) =
  #    split //, $individual_hash{GT};
  #  my $this_ref = $reference_seq;
  #
  #  my %variant_seq_hash = map {$_ => 1} @all_seqs[$var1, $var2];
  #  my @variant_seqs = keys %variant_seq_hash;
  #
  #  my @these_all_seqs = ($this_ref, @variant_seqs);
  #  my @lengths = map {length($_)} @these_all_seqs;
  #
  #  my $this_start = $start;
  #  my $this_end   = $end;
  #  #INDEL
  #  if (grep {$_ != 1} @lengths) {
  #    my ($min_length) = sort {$a <=> $b} @lengths;
  #    map {substr($_, 0, $min_length, '')} @these_all_seqs;
  #    map {$_ ||= '-'} @these_all_seqs;
  #    my $start_adjust = length($reference_seq) - length($these_all_seqs[0]);
  #    $this_start += $start_adjust;
  #    my $end_adjust = length($these_all_seqs[0]) - 1;
  #    $this_end = $this_start + $end_adjust;
  #
  #    $this_ref = shift @these_all_seqs;
  #    @variant_seqs = @these_all_seqs;
  #
  #    if ($this_ref eq '-') {
  #	$type = 'insertion';
  #    }
  #    elsif (grep {$_ eq '-'} @variant_seqs) {
  #	$type = 'deletion';
  #    }
  #    else {
  #	$type = 'complex_substitution';
  #	  }
  #  }
  #
  #  my $read_count = $individual_hash{DP} || $info{DP}[0] || undef;
  #
  #  # Make sure that at least one variant sequence differs from
  #  # the reference sequence.
  #  next unless grep {$_ ne $this_ref} @variant_seqs;
  #
  #  # Format Column Keys
  #  # GT : genotype, encoded as alleles values separated by either of
  #  # ”/” or “|”, e.g. The allele values are 0 for the reference
  #  # allele (what is in the reference sequence), 1 for the first
  #  # allele listed in ALT, 2 for the second allele list in ALT and
  #  # so on. For diploid calls examples could be 0/1 or 1|0 etc. For
  #  # haploid calls, e.g. on Y, male X, mitochondrion, only one
  #  # allele value should be given. All samples must have GT call
  #  # information; if a call cannot be made for a sample at a given
  #  # locus, ”.” must be specified for each missing allele in the GT
  #  # field (for example ./. for a diploid). The meanings of the
  #  # separators are:
  #  #       / : genotype unphased
  #  #       | : genotype phased
  #  # DP : read depth at this position for this sample (Integer)
  #  # FT : sample genotype filter indicating if this genotype was
  #  # “called” (similar in concept to the FILTER field). Again, use
  #  # PASS : to indicate that all filters have been passed, a
  #  # semi-colon separated list of codes for filters that fail, or
  #  # ”.” to indicate that filters have not been applied. These
  #  # values should be described in the meta-information in the same
  #  # way as FILTERs (Alphanumeric String)
  #  # GL : three floating point log10-scaled likelihoods for
  #  # AA,AB,BB genotypes where A=ref and B=alt; not applicable if
  #  # site is not biallelic. For example: GT:GL
  #  # 0/1:-323.03,-99.29,-802.53 (Numeric)
  #  # GQ : genotype quality, encoded as a phred quality
  #  # -10log_10p(genotype call is wrong) (Numeric)
  #  # HQ : haplotype qualities, two phred qualities comma separated
  #  # (Numeric)
  #
  #  my $source = $individual_id;
  #  my $feature_id = join ':', ($seqid, $source, $start);
}


##-----------------------------------------------------------------------------
#
#sub reader {
#  my ($self, $record) = @_;
#
#
#  # INFO Column Keys
#  # VFC 4.1
#
#  # INFO additional information: (String, no white-space, semi-colons,
#  # or equals-signs permitted; commas are permitted only as delimiters
#  # for lists of values) INFO fields are encoded as a
#  # semicolon-separated series of short keys with optional values in
#  # the format: <key>=<data>[,data]. Arbitrary keys are permitted,
#  # although the following sub-fields are reserved (albeit optional):
#  #  AA : ancestral allele
#  #  AC : allele count in genotypes, for each ALT allele, in the same order as listed
#  #  AF : allele frequency for each ALT allele in the same order as listed: use this when estimated from primary data, not called genotypes
#  #  AN : total number of alleles in called genotypes
#  #  BQ : RMS base quality at this position
#  #  CIGAR : cigar string describing how to align an alternate allele to the reference allele
#  #  DB : dbSNP membership
#  #  DP : combined depth across samples, e.g. DP=154
#  #  END : end position of the variant described in this record (for use with symbolic alleles)
#  #  H2 : membership in hapmap2
#  #  H3 : membership in hapmap3
#  #  MQ : RMS mapping quality, e.g. MQ=52
#  #  MQ0 : Number of MAPQ == 0 reads covering this record
#  #  NS : Number of samples with data
#  #  SB : strand bias at this position
#  #  SOMATIC : indicates that the record is a somatic mutation, for cancer genomics
#  #  VALIDATED : validated by follow-up experiment
#  #  1000G : membership in 1000 Genomes
#
#  my @pairs = split /;/, $record->{info};
#  my %info_attributes;
#  for my $pair (@pairs) {
#    my ($key, $value) = split /=/, $pair;
#    my @values = split /,/, $value;
#    # if ($key eq 'AA') {
#    #
#    # }
#    my $gvf_key = 'vcf' . $key;
#    $info_attributes{$gvf_key => \@values};
#  }
#  return wantarray ? %info_attributes : \%info_attributes;
#}

#-----------------------------------------------------------------------------

sub reader {
  my $self = shift;

  if (! $self->{reader}) {
    # We don't give field names here so that the record will be returned
    # as an array rather than a hash.  We splice out the first
    # 9 columns above and the rest in genotype data.
    # my @field_names = qw(chrom pos id ref alt qual filter info format);
    my $reader =
      GAL::Reader::DelimitedLine->new(metadata_pattern => qr/^\#\#/,
				      metadata_parser  => \&_parse_vcf_metadata,
				      header_count     => 1,
				      header_parser    => \&_parse_vcf_header,
				     );
    $self->{reader} = $reader;
  }
  return $self->{reader};
}

#-----------------------------------------------------------------------------

=head1 DIAGNOSTICS

=for author to fill in:
     List every single error and warning message that the module can
     generate (even the ones that will "never happen"), with a full
     explanation of each problem, one or more likely causes, and any
     suggested remedies.

=over

=item C<< Error message here, perhaps with %s placeholders >>

[Description of error here]

=item C<< Another error message here >>

[Description of error here]

[Et cetera, et cetera]

=back

=head1 CONFIGURATION AND ENVIRONMENT

<GAL::Parser::VCFv4_1> requires no configuration files or environment variables.

=head1 DEPENDENCIES

None.

=head1 INCOMPATIBILITIES

None reported.

=head1 BUGS AND LIMITATIONS

No bugs have been reported.

Please report any bugs or feature requests to:
barry.moore@genetics.utah.edu

=head1 AUTHOR

Barry Moore <barry.moore@genetics.utah.edu>

=head1 LICENCE AND COPYRIGHT

Copyright (c) 2010-2014, Barry Moore <barry.moore@genetics.utah.edu>.  All rights reserved.

    This module is free software; you can redistribute it and/or
    modify it under the same terms as Perl itself (See LICENSE).

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
LIABLE TO YOU FOR DAMAGES, INCLUDING ANY GENERAL, SPECIAL, INCIDENTAL,
OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OR INABILITY TO USE
THE SOFTWARE (INCLUDING BUT NOT LIMITED TO LOSS OF DATA OR DATA BEING
RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR THIRD PARTIES OR A
FAILURE OF THE SOFTWARE TO OPERATE WITH ANY OTHER SOFTWARE), EVEN IF
SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE POSSIBILITY OF
SUCH DAMAGES.

=cut

1;
