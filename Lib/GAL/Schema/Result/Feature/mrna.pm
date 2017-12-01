package GAL::Schema::Result::Feature::mrna;

use strict;
use warnings;
use base qw(GAL::Schema::Result::Feature::transcript);

=head1 NAME

GAL::Schema::Result::Feature::mrna - A mRNA object for the GAL Library 

=head1 VERSION

This document describes GAL::Schema::Result::Feature::mrna version 0.2.0

=head1 SYNOPSIS

    use GAL::Annotation;
    my $feat_store = GAL::Annotation->new(storage => $feat_store_args,
					  parser  => $parser_args,
					  fasta   => $fasta_args,
					 );

    $feat_store->load_files(files => $feature_file,
			    mode  => 'overwrite',
			    );

    my $features = $feat_store->schema->resultset('Feature');

    my $mrnas = $features->search({type => 'mRNA'});
	
=head1 DESCRIPTION

<GAL::Schema::Result::Feature::mrna> provides a <GAL::Schema::Result::Feature>
subclass for mRNA specific behavior. mRNA, messenger RNA is the intermediate 
molecule between DNA and protein. It includes UTR and coding sequences, but not
introns.

=head1 METHODS

=cut

#-----------------------------------------------------------------------------

=head2 CDSs

 Title   : CDSs
 Usage   : $CDSs = $self->CDSs
 Function: Get the mRNA's coding sequences sorted in the order of the mRNA's strand.
 Returns : A DBIx::Class::Result object loaded up with CDS's
 Args    : None

=cut

sub CDSs {

  my $self = shift;

  if (wantarray) {
      return grep {$_->type eq 'CDS'} $self->children->all;
  }
  else {
      my $sort_order = ($self->strand eq '-' ?
			{'-desc' => 'end'}   :
			{'-asc'  => 'start'});

      my $CDSs = $self->children({type => 'CDS'},
				  {order_by => $sort_order,
				   distinct => 1});
      return $CDSs;
  }
}

#-----------------------------------------------------------------------------

=head2 CDS_seq_genomic

 Title   : CDS_seq_genomic
 Usage   : $seq = $self->CDS_seq_genomic
 Function: Return the genomic sequence (on the plus strand of the genome) of
	   the full CDS for this mRNA.
 Returns : A DNA sequence
 Args    : None

=cut

sub CDS_seq_genomic {

  my $self = shift;

  my $CDS_seq_genomic = '';
  for my $CDS (sort {$a->start <=> $b->start} $self->CDSs) {
    my $this_seq = $CDS->genomic_seq;
    $CDS_seq_genomic .= $this_seq if $this_seq;
  }
  if (! $CDS_seq_genomic) {
      $self->warn('no_CDS_seq_available', $self->to_gff3);
  }
  return $CDS_seq_genomic;
}

#-----------------------------------------------------------------------------

=head2 CDS_seq

 Title   : CDS_seq
 Usage   : $seq = $self->CDS_seq
 Function: Return the genomic sequence (on the coding direction) 
	of the full CDS for this mRNA.
 Returns : A DNA sequence
 Args    : None

=cut

sub CDS_seq {

  my $self = shift;

  my $CDS_seq = $self->CDS_seq_genomic;
  if ($self->strand eq '-') {
    $CDS_seq = $self->annotation->revcomp($CDS_seq);
  }
  return $CDS_seq;
}

#-----------------------------------------------------------------------------

=head2 protein_seq

 Title   : protein_seq
 Usage   : $seq = $self->protein_seq
 Function: Return the protein sequence for this mRNA.
 Returns : Protein sequence
 Args    : None

=cut

sub protein_seq {

  my $self = shift;

  my $CDS_seq = $self->CDS_seq;
  my $protein_seq = $self->annotation->translate($CDS_seq);
  $protein_seq ||= '';
  if (! $protein_seq) {
      $self->warn('no_protein_sequence_available', $self->to_gff3);
  }
  $protein_seq =~ s/\*$// if $protein_seq;
  return $protein_seq;
}

#-----------------------------------------------------------------------------

=head2 map2CDS

 Title   : map2CDS
 Usage   : @CDS_coordinates = $self->map2CDS(@genomic_coordinates);
	   ($CDS_coordinate) = $self->map2CDS($genomic_coordinate);
 Function: Transform genomic coordinates to CDS coordinates.
 Returns : An array(ref) of integers or an empty list if the genomic
	   coordinates given do not overlap the CDS of this mRNA.
 Args    : An array of integers

=cut

sub map2CDS {

  my ($self, @coordinates) = @_;

  my ($CDS_start) = $self->genome2me($self->CDS_start);
  my ($CDS_end)   = $self->genome2me($self->CDS_end);

  my @CDS_coordinates = $self->genome2me(@coordinates);
  for my $coordinate (@CDS_coordinates) {
    if ($coordinate && $coordinate >= $CDS_start && $coordinate <= $CDS_end) {
      $coordinate = $coordinate - $CDS_start + 1;
    }
    else {
      $coordinate = undef;
    }
  }

  return wantarray ? @CDS_coordinates : \@CDS_coordinates;
}

#-----------------------------------------------------------------------------

=head2 map2protein

 Title   : map2protein
 Usage   : @protein_coordinates = $self->map2protein(@genomic_coordinates);
 Function: Transform genomic coordinates to protein sequence coordinates.
	   Note that 3 genomic coordinates will return the same protein
	   coordinate if they fall in the same codon.
 Returns : An array(ref) of integers or an empty list if the genomic
	   coordinates given do not overlap the protein (CDS) of this mRNA.
 Args    : An array of integers

=cut

sub map2protein {

  my ($self, @coordinates) = @_;

  my @protein_coordinates = $self->map2CDS(@coordinates);
  map {$_ = int(($_ - 1) / 3) + 1 if $_} @protein_coordinates;

  return wantarray ? @protein_coordinates : \@protein_coordinates;
}

#-----------------------------------------------------------------------------

=head2 CDS_start

 Title   : CDS_start
 Usage   : $start = $self->CDS_start
 Function: Returns the genomic coordinate of the start of this CDS.
 Returns : An integer
 Args    : None

=cut

sub CDS_start {

  my $self = shift;
  my $strand = $self->strand;
  my @CDSs = sort {$a->start <=> $b->start} $self->CDSs;
  my $CDS_start = $strand eq '-' ? $CDSs[-1]->end : $CDSs[0]->start;
  return $CDS_start;
}

#-----------------------------------------------------------------------------

=head2 CDS_end

 Title   : CDS_end
 Usage   : $end = $self->CDS_end
 Function: Returns the genomic coordinate of the end of this CDS.
 Returns : An integer
 Args    : None

=cut

sub CDS_end {

  my $self = shift;
  my $strand = $self->strand;
  my @CDSs = sort {$a->start <=> $b->start} $self->CDSs;
  my $CDS_end = $strand eq '-' ? $CDSs[0]->start : $CDSs[-1]->end;
  return $CDS_end;
}

#-----------------------------------------------------------------------------

=head2 CDS_length

 Title   : CDS_length
 Usage   : $length = $self->CDS_length
 Function: Return the length of the CDS for this mRNA.
 Returns : An integer
 Args    : None

=cut

sub CDS_length {

  my $self = shift;
  my $CDS_length;
  map {$CDS_length += $_->genomic_length} $self->CDSs;
  return $CDS_length;
}

#-----------------------------------------------------------------------------

=head2 protein_length

 Title   : protein_length
 Usage   : $length = $self->protein_length
 Function: Return the length of the protein for this mRNA. If the CDS length
	   is not divisible by three this method will truncate the integral
	   portion of the protein length.
 Returns : An integer
 Args    : None

=cut

sub protein_length {
  my $self = shift;
  return length $self->protein_seq;
}

#-----------------------------------------------------------------------------

=head2 phase_at_location

 Title   : phase_at_location
 Usage   : $phase = $self->phase_at_location($genomic_coordinate)
 Function: Return the phase (how many nts 3' does the next codon begin) for
	   a given genomic coordinate or undef if the coordinate does not overlap
	   the CDS of this mRNA.
 Returns : An integer
 Args    : An integer (genomic coordinate).

=cut

sub phase_at_location {

  my ($self, $location) = @_;

  my %mod2phase = (1 => 0,
		   2 => 2,
		   0 => 1,
		  );
  my ($CDS_location) = $self->map2CDS($location);
  return undef unless $CDS_location;
  my $modulus = $CDS_location % 3;
  return $mod2phase{$modulus};
}

#-----------------------------------------------------------------------------

=head2 frame_at_location

 Title   : frame_at_location
 Usage   : $frame = $self->frame_at_location($genomic_coordinate);
 Function: Return the reading frame (what position within a codon) for a given
	   genomic coordinate or undef if the coordinate does not overlap
	   the CDS of this mRNA.
 Returns : An integer
 Args    : An integer

=cut

sub frame_at_location {

  my ($self, $location) = @_;
  my %mod2frame = (1 => 0,
		   2 => 1,
		   0 => 2,
		  );
  my ($CDS_location) = $self->map2CDS($location);
  return undef unless $CDS_location;
  my $modulus = $CDS_location % 3;
  my $frame = $mod2frame{$modulus};
  return $frame;
}

#-----------------------------------------------------------------------------

=head2 codon_at_location

 Title   : codon_at_location
 Usage   : $codon = $self->codon_at_location($genomic_coordinate);
	   ($codon, $frame) = $self->codon_at_location
 Function: Return the codon (and frame in list context) that overlaps a given
	   genomic coordinate.  Returns undef in the genomic coordinate does
	   not overlap a complete codon.
 Returns : An scalar string and in list context also returns an integer.
 Args    : An integer

=cut

sub codon_at_location {

  my ($self, $location) = @_;

  my $CDS_sequence = $self->CDS_seq;
  my ($CDS_location) = $self->map2CDS($location);
  return undef unless $CDS_location;
  my $frame = $self->frame_at_location($location);
  my $codon_start = $CDS_location - $frame;
  $codon_start--;
  return undef if (($codon_start + 3) > length $CDS_sequence);
  my $codon = substr($CDS_sequence, $codon_start, 3);
  return wantarray ? ($codon, $frame) : $codon;
}

#-----------------------------------------------------------------------------

=head2 translation_start

 Title   : translation_start
 Usage   : $tans_start = $mrna->translation_start;
 Function: Get the coordinate of the first nt of the start codon.
 Returns : An integer
 Args    : None

=cut

sub translation_start {
  my $self = shift;

  my $strand = $self->strand;
  my $sort_sub = ($strand eq '-'            ?
		  sub {$b->start <=> $a->start} :
		  sub {$a->start <=> $b->start});

  my ($first_CDS) = sort $sort_sub $self->CDSs;
  return $first_CDS->my_start;
}

#-----------------------------------------------------------------------------

=head2 translation_end

 Title   : translation_end
 Usage   : $tans_end = $mrna->translation_end;
 Function: Get the coordinate of the last nt of the stop codon.
 Returns : An integer
 Args    : None

=cut

sub translation_end {
  my $self = shift;

  my $strand = $self->strand;
  my $sort_sub = ($strand eq '-'            ?
		  sub {$a->start <=> $b->start} :
		  sub {$b->start <=> $a->start});

  my ($last_CDS) = sort $sort_sub $self->CDSs;
  return $last_CDS->my_end;
}

#-----------------------------------------------------------------------------

=head2 five_prime_UTRs

 Title   : five_prime_UTRs
 Usage   : $five_prime_UTRs = $self->five_prime_UTRs
 Function: Get the features five_prime_UTRs
 Returns : A DBIx::Class::Result object loaded up with five_prime_UTRs
 Args    : None

=cut

sub five_prime_UTRs {
  my $self = shift;

  if (wantarray) {
      return grep {$_->type eq 'five_prime_UTR'} $self->children->all;
  }
  else {
      my $sort_order = ($self->strand eq '-' ?
			{'-desc' => 'end'}   :
			{'-asc'  => 'start'});

      my $five_prime_UTRs = $self->children({type => 'five_prime_UTR'},
					    {order_by => $sort_order,
					     distinct => 1});
      return $five_prime_UTRs;
  }
}

#-----------------------------------------------------------------------------

=head2 three_prime_UTRs

 Title   : three_prime_UTRs
 Usage   : $three_prime_UTRs = $self->three_prime_UTRs
 Function: Get the features three_prime_UTRs
 Returns : A DBIx::Class::Result object loaded up with three_prime_UTRs
 Args    : None

=cut

sub three_prime_UTRs {
  my $self = shift;

  if (wantarray) {
      return grep {$_->type eq 'three_prime_UTR'} $self->children->all;
  }
  else {
    my $sort_order = ($self->strand eq '-' ?
		      {'-desc' => 'end'}   :
		      {'-asc'  => 'start'});

    my $three_prime_UTRs = $self->children({type => 'three_prime_UTR'},
					   {order_by => $sort_order,
					    distinct => 1});
    return $three_prime_UTRs;
  }
}

#-----------------------------------------------------------------------------

=head2 infer_five_prime_UTR

 Title   : infer_five_prime_UTR
 Usage   : $five_prime_UTR = $self->infer_five_prime_UTR
 Function: Infer five_prime_UTR for the features
 Returns : A DBIx::Class::Result object loaded up with five_prime_UTR
 Args    : None

=cut

sub infer_five_prime_UTR {
  my ($self, $five_prime_UTR) = @_;

  $five_prime_UTR ||= $self->children({type => 'five_prime_UTR'},
				      {order_by => { -asc => 'start' },
				       distinct => 1});


  my $strand = $self->strand;
  my $sort_sub = ($strand eq '-'            ?
		  sub {$b->start <=> $a->start} :
		  sub {$a->start <=> $b->start});

  my $translation_start = $self->translation_start;

  my @coordinates;
  for my $exon (sort $sort_sub $self->exons) {
    if ($strand eq '-') {
      last if $exon->end <= $translation_start;
      unshift @coordinates, [$exon->start, $exon->end];
    }
    else {
      last if $exon->start >= $translation_start;
      push @coordinates, [$exon->start, $exon->end];
    }
  }

  my %template = (seqid  => $self->seqid,
		  source => $self->source,
		  type   => 'five_prime_UTR',
		  score  => '.',
		  strand => $strand,
		  phase  => '.',
 );

  my $parent_id = $self->feature_id;

  my $count = 1;
  my @five_prime_UTR;
  for my $coordinate_pair (@coordinates) {
    my ($start, $end) = @{$coordinate_pair};
    if ($strand eq '-') {
      $start = $translation_start + 1
	if $start <= $translation_start;
    }
    else {
      $end = $translation_start - 1
	if $end >= $translation_start;
    }
    my $feature_id = $parent_id . ':five_prime_UTR:';
    $feature_id .= sprintf("%03d", $count++);
    my @attrbs = ({att_key    => 'ID',
		   att_value  => $feature_id},
		  {att_key    => 'Parent',
		   att_value  => $parent_id},
		 );
    my @rels = {parent => $parent_id,
		child  => $feature_id};

    my %five_prime_UTR = %template;
    @five_prime_UTR{qw(feature_id start end attributes my_parents)} =
      ($feature_id, $start, $end, \@attrbs, \@rels);
    my $five_prime_UTR = $five_prime_UTR->find_or_create(\%five_prime_UTR);
    bless $five_prime_UTR, 'GAL::Schema::Result::Feature::five_prime_UTR';
    push @five_prime_UTR, $five_prime_UTR;
  }
  $five_prime_UTR->set_cache(\@five_prime_UTR);
  print '';
}

#-----------------------------------------------------------------------------

=head2 infer_three_prime_UTR

 Title   : infer_three_prime_UTR
 Usage   : $three_prime_UTR = $self->infer_three_prime_UTR
 Function: Infer three_prime_UTR for the features
 Returns : A DBIx::Class::Result object loaded up with three_prime_UTR
 Args    : None

=cut

sub infer_three_prime_UTR {
  my ($self, $three_prime_UTRs) = @_;

  $three_prime_UTRs ||= $self->children({type => 'three_prime_UTR'},
				       {order_by => { -asc => 'start' },
					distinct => 1});


  my $strand = $self->strand;
  my $sort_sub = ($strand eq '-'            ?
		  sub {$a->start <=> $b->start} :
		  sub {$b->start <=> $a->start});

  my $translation_end = $self->translation_end;

  my @coordinates;
  for my $exon (sort $sort_sub $self->exons) {
    if ($strand eq '-') {
      last if $exon->start >= $translation_end;
      push @coordinates, [$exon->start, $exon->end];
    }
    else {
      last if $exon->end <= $translation_end;
      unshift @coordinates, [$exon->start, $exon->end];
    }
  }

  my %template = (seqid  => $self->seqid,
		  source => $self->source,
		  type   => 'three_prime_UTR',
		  score  => '.',
		  strand => $strand,
		  phase  => '.',
 );

  my $parent_id = $self->feature_id;

  my $count = 1;
  my @three_prime_UTRs;
  for my $coordinate_pair (@coordinates) {
    my ($start, $end) = @{$coordinate_pair};
    if ($strand eq '-') {
      $end = $translation_end - 1
	if $end >= $translation_end;
    }
    else {
      $start = $translation_end + 1
	if $start <= $translation_end;
    }
    my $feature_id = $parent_id . ':three_prime_UTR:';
    $feature_id .= sprintf("%03d", $count++);
    my @attrbs = ({att_key    => 'ID',
		   att_value  => $feature_id},
		  {att_key    => 'Parent',
		   att_value  => $parent_id},
		 );
    my @rels = {parent => $parent_id,
		child  => $feature_id};

    my %three_prime_UTR = %template;
    @three_prime_UTR{qw(feature_id start end attributes my_parents)} =
      ($feature_id, $start, $end, \@attrbs, \@rels);
    my $three_prime_UTR = $three_prime_UTRs->find_or_create(\%three_prime_UTR);
    bless $three_prime_UTR, 'GAL::Schema::Result::Feature::three_prime_UTR';
    push @three_prime_UTRs, $three_prime_UTR;
  }
  $three_prime_UTRs->set_cache(\@three_prime_UTRs) if @three_prime_UTRs;
  print '';
}

#-----------------------------------------------------------------------------

=head2 five_prime_UTR_seq_genomic

 Title   : five_prime_UTR_seq_genomic
 Usage   : $seq = $self->five_prime_UTR_seq_genomic
 Function: Get the transcripts spliced five_prime_UTR genomic sequence (not
   reverse complimented for minus strand features.
 Returns : A text string of the five_prime_UTR spliced genomic sequence.
 Args    : None

=cut

sub five_prime_UTR_seq_genomic {
  my $self = shift;
  my $five_prime_UTR_seq_genomic;
  map {$five_prime_UTR_seq_genomic .= $_->genomic_seq} $self->five_prime_UTRs->all;
  return $five_prime_UTR_seq_genomic;
}

#-----------------------------------------------------------------------------

=head2 five_prime_UTR_seq

 Title   : five_prime_UTR_seq
 Usage   : $seq = $self->five_prime_UTR_seq
 Function: Get the transcripts spliced five_prime_UTR sequence reverse
   complimented for minus strand features.
 Returns : A text string of the five_prime_UTR spliced sequence.
 Args    : None

=cut

sub five_prime_UTR_seq {
  my $self = shift;
  my $five_prime_UTR_seq = $self->five_prime_UTR_seq_genomic;
  if ($self->strand eq '-') {
    $five_prime_UTR_seq =
      $self->annotation->revcomp($five_prime_UTR_seq);
  }
  return $five_prime_UTR_seq;
}

#-----------------------------------------------------------------------------

=head2 three_prime_UTR_seq_genomic

 Title   : three_prime_UTR_seq_genomic
 Usage   : $seq = $self->three_prime_UTR_seq_genomic
 Function: Get the transcripts spliced three_prime_UTR genomic sequence (not
   reverse complimented for minus strand features.
 Returns : A text string of the three_prime_UTR spliced genomic sequence.
 Args    : None

=cut

sub three_prime_UTR_seq_genomic {
  my $self = shift;
  my $three_prime_UTR_seq_genomic;
  map {$three_prime_UTR_seq_genomic .= $_->genomic_seq} $self->three_prime_UTRs->all;
  return $three_prime_UTR_seq_genomic;
}

#-----------------------------------------------------------------------------

=head2 three_prime_UTR_seq

 Title   : three_prime_UTR_seq
 Usage   : $seq = $self->three_prime_UTR_seq
 Function: Get the transcripts spliced three_prime_UTR sequence reverse
   complimented for minus strand features.
 Returns : A text string of the three_prime_UTR spliced sequence.
 Args    : None

=cut

sub three_prime_UTR_seq {
  my $self = shift;
  my $three_prime_UTR_seq = $self->three_prime_UTR_seq_genomic;
  if ($self->strand eq '-') {
    $three_prime_UTR_seq =
      $self->annotation->revcomp($three_prime_UTR_seq);
  }
  return $three_prime_UTR_seq;
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

<GAL::Schema::Result::Feature::mrna> requires no configuration files or environment variables.

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
