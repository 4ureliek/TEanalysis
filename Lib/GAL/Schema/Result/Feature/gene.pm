package GAL::Schema::Result::Feature::gene;

use strict;
use warnings;
use base qw(GAL::Schema::Result::Feature::sequence_feature);

=head1 NAME

GAL::Schema::Result::Feature::gene -  A gene object for the GAL Library

=head1 VERSION

This document describes GAL::Schema::Result::Feature::gene version 0.2.0

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

    my $genes = $features->search({type => 'gene'});
    while (my $gene = $genes->next) {
      my $mrnas = $gene->mRNAs;
      while (my $mrna = $mrnas->next) {
	my $id    = $mrna->feature_id;
	my $start = $mrna->start;
	my $end   = $mrna->end;
      }
    }

=head1 DESCRIPTION

<GAL::Schema::Result::Feature::gene> provides a <GAL::Schema::Result::Feature>
subclass for gene specific behavior.A gene is comprised all of the sequence 
elements necessary to encode a function transcript, including regulatory regions, 
transcribed regions and/or other functional sequence regions.


Gene= "The functional and physical unit of heredity passed from parent to offspring. Genes are pieces of DNA, and most genes contain the information for making a specific protein." (UMLS-NLM, 2014)
Through this class many methods offers information about the selected gene such as: returning the transcripts and their types; infer introns within the gene; alternative splicing; and may check if the gene has mRNA or not 
=head1 METHODS

=cut

#-----------------------------------------------------------------------------

=head2 transcripts

 Title   : transcripts
 Usage   : $transcripts = $self->transcripts
 Function: This method will get the gene's transcript sorted by start with the longest
	   transcripts first when start locations are equal.
 Returns : A DBIx::Class::Result object loaded up with transcripts.
 Args    : None

=cut

sub transcripts {

  my $self = shift;

  #TODO: GAL::lib::GAL::Schema::Result::Feature::gene::transcripts
  #TODO: should use SO directly.

  my $sort_order;
  if ($self->strand eq '-') {
    $sort_order = [{ -desc => 'end' },
		   { -asc =>  'start'},
		  ];
  }
  else {
    $sort_order = [{ -asc  => 'start' },
		   { -desc =>  'end'},
		  ];
  }

  my $transcript_types = $self->get_transcript_types;
  my $transcripts = $self->children({type => $transcript_types},
				    {order_by => $sort_order,
				    distinct => 1});

  return wantarray ? $transcripts->all : $transcripts;

}

#-----------------------------------------------------------------------------

=head2 infer_introns

 Title   : infer_introns
 Usage   : $self->infer_introns
 Function: Infer introns for all transcripts belonging to the gene.
 Returns : N/A
 Args    : None

=cut

sub infer_introns {

  my $self = shift;
  for my $transcript ($self->transcripts) {
    $transcript->infer_introns;
  }
}

#-----------------------------------------------------------------------------

=head2 splice_complexity

 Title   : splice_complexity
 Usage   : $splice_complexity = $self->splice_complexity
 Function: Get the transcript splice_complexity (PMID: 19236712)
 Returns : A numerical value
 Args    : None

=cut

sub splice_complexity {

  my $self = shift;

  my @transcripts = $self->transcripts;
  return undef unless scalar @transcripts > 1;

  my @AEDs;
  my %seen_pairs;
  for my $transcriptA (@transcripts) {
    my $idA = $transcriptA->feature_id;
    for my $transcriptB (@transcripts) {
      my $idB = $transcriptB->feature_id;
      my $pair = join '-', sort($idA, $idB);
      next if ($idA eq $idB || exists $seen_pairs{$pair});
      $seen_pairs{$pair}++;
      my $AED = $transcriptA->AED($transcriptB);
      push @AEDs, $AED;
    }
  }

  my $total_AED;
  map {$total_AED += $_} @AEDs;
  return $total_AED / scalar @AEDs;
}

#-----------------------------------------------------------------------------

=head2 mRNAs

 Title   : mRNAs
 Usage   : $mRNAs = $self->mRNAs
 Function: Get the genes mRNA features
 Returns : A DBIx::Class::Result object loaded up with mRNA features.
 Args    : None

=cut

sub mRNAs {
  my $self = shift;

  my $sort_order;
  if ($self->strand eq '-') {
      $sort_order = [{ -desc => 'end' },
		     { -asc =>  'start'},
	  ];
  }
  else {
      $sort_order = [{ -asc  => 'start' },
		     { -desc =>  'end'},
	  ];
  }


  #TODO: should use SO directly.
  my $mRNAs = $self->children({type => 'mRNA'},
			      {order_by => $sort_order,
			       distinct => 1});
  return wantarray ? $mRNAs->all : $mRNAs;
}

#-----------------------------------------------------------------------------
# This probably shouldn't be here, but leaving it as comment for now 10/1/12 BM

# sub feature_id {
#
#     my $self = shift;
#     return $self->id;
#
# }
#-----------------------------------------------------------------------------

=head2 is_coding

 Title   : is_coding
 Usage   : $is_coding = $self->is_coding
 Function: Return true if this gene has mRNA - this function will be
           deprecated in favor of has_mRNA
 Returns : 1 or undef
 Args    : None

=cut

sub is_coding {
    my $self = shift;
    return 1 if ($self->has_mRNA || $self->has_CDS);
}

#-----------------------------------------------------------------------------

=head2 has_mRNA

 Title   : has_mRNA
 Usage   : $has_mRNA = $self->has_mRNA
 Function: Return true if this gene has any mRNA
 Returns : 1 or undef
 Args    : None

=cut

sub has_mRNA {

  my $self = shift;
  return 1 if grep {$_->type eq 'mRNA'} $self->children(undef, {distinct => 1})->all;
}

#-----------------------------------------------------------------------------

=head2 has_CDS

 Title   : has_CDS
 Usage   : $has_CDS = $self->has_CDS
 Function: Return true if this gene has and CDS children
 Returns : 1 or undef
 Args    : None

=cut

sub has_CDS {

  my $self = shift;
  my @all_children;
  $self->get_recursive_children(\@all_children);
  return 1 if grep {$_->type eq 'CDS'} @all_children;
}

#-----------------------------------------------------------------------------

##############################################################################
## Consolidate three copies of this method 1) Here 2) Feature/gene.pm and
## 3) Base.pm
##############################################################################

=head2 get_transcript_types

 Title   : get_transcript_types
 Usage   : $transcript_types = $self->get_transcript_types
 Function: Get a list of all features that are types of transcripts.
 Returns : An array or ref
 Args    : None

=cut

sub get_transcript_types {

  my @transcripts = qw(EST
		       RNase_MRP_RNA
		       RNase_P_RNA
		       SRP_RNA
		       SRP_RNA_primary_transcript
		       Y_RNA
		       aberrant_processed_transcript
		       alternatively_spliced_transcript
		       antisense_RNA
		       antisense_primary_transcript
		       capped_mRNA
		       capped_primary_transcript
		       class_II_RNA
		       class_I_RNA
		       consensus_mRNA
		       dicistronic_mRNA
		       dicistronic_primary_transcript
		       dicistronic_transcript
		       edited_mRNA
		       edited_transcript
		       edited_transcript_by_A_to_I_substitution
		       enhancerRNA
		       enzymatic_RNA
		       exemplar_mRNA
		       guide_RNA
		       lnc_RNA
		       mRNA
		       mRNA_recoded_by_codon_redefinition
		       mRNA_recoded_by_translational_bypass
		       mRNA_region
		       mRNA_with_frameshift
		       mRNA_with_minus_1_frameshift
		       mRNA_with_minus_2_frameshift
		       mRNA_with_plus_1_frameshift
		       mRNA_with_plus_2_frameshift
		       mature_transcript
		       mature_transcript_region
		       miRNA_primary_transcript
		       mini_exon_donor_RNA
		       monocistronic_mRNA
		       monocistronic_primary_transcript
		       monocistronic_transcript
		       ncRNA
		       nc_primary_transcript
		       piRNA
		       polyadenylated_mRNA
		       polycistronic_mRNA
		       polycistronic_primary_transcript
		       polycistronic_transcript
		       pre_edited_mRNA
		       primary_transcript
		       primary_transcript_region
		       processed_transcript
		       protein_coding_primary_transcript
		       pseudogenic_transcript
		       rRNA
		       rRNA_cleavage_RNA
		       rRNA_primary_transcript
		       rasiRNA
		       recoded_mRNA
		       regional_centromere_outer_repeat_transcript
		       riboswitch
		       ribozyme
		       scRNA
		       scRNA_primary_transcript
		       siRNA
		       small_regulatory_ncRNA
		       snRNA
		       snRNA_primary_transcript
		       snoRNA
		       snoRNA_primary_transcript
		       spliced_leader_RNA
		       stRNA
		       tRNA
		       tRNA_primary_transcript
		       tasiRNA
		       tasiRNA_primary_transcript
		       telomerase_RNA
		       tmRNA_primary_transcript
		       trans_spliced_mRNA
		       trans_spliced_transcript
		       transcript
		       transcript_bound_by_nucleic_acid
		       transcript_bound_by_protein
		       transcript_region
		       transcript_with_translational_frameshift
		       vault_RNA
		     );
  return wantarray ? @transcripts : \@transcripts;
}

#----------------------------------------------------------------------------

=head1 DIAGNOSTICS

This module does not throw any error or warning messages.

=head1 CONFIGURATION AND ENVIRONMENT

<GAL::Schema::Result::Feature::gene> requires no configuration files or environment variables.

=head1 DEPENDENCIES

<GAL::Schema::Result::Feature::sequence_feature>

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
