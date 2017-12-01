package GAL::Schema::Result::Feature;

use strict;
use warnings;
use base qw/DBIx::Class/;
use Set::IntSpan::Fast;
use Storable qw(dclone);
use List::MoreUtils qw(uniq);

=head1 NAME

GAL::Schema::Result::Feature - Base class for all sequence features

=head1 VERSION

This document describes GAL::Schema::Result::Feature version 0.2.0

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
    while (my $mrna = $mrnas->next) {
      my $cdss = $mrna->CDSs;
      while (my $cds = $cdss->next) {
	my $cds   = $cds->feature_id;
	my $start = $cds->start;
	my $end   = $cds->end;
      }
    }

=head1 DESCRIPTION

<GAL::Schema::Result::Feature> is the base class for all sequence
feature classes in the GAL's <DBIx::Class> based feature retrival
system.  You don't instantiate objects from it or any of it's
subclasses yourself, <DBIx::Class> does that for you, hence there is
no constructor and there are no attributes available, but you will
want to have a look at the methods as they provide the base
functionality for all sequence feature subclasses.  to all subclasses
as well.  The first 8 methods provide access to the features primary
data as defined by the first 8 columns of GFF3
(http://www.sequenceontology.org/).

=head1 METHODS

=cut

# Map feature types to a parent that has an subclass.

my %FEATURE_MAP;

map {$FEATURE_MAP{$_} = 'transcript'}
  GAL::Schema::Result::Feature::get_transcript_types();

map {$FEATURE_MAP{$_} = lc $_}  qw(gene
				     transcript
				     mRNA
				     exon
				     intron
				     three_prime_UTR
				     five_prime_UTR
				   );

# Eventually get this from SO at runtime.
my @sequence_alterations = qw(copy_number_variation deletion indel
			      insertion duplication tandem_duplication
			      transgenic_insertion inversion substitution
			      MNP SNV SNP point_mutation transition
			      purine_transition A_to_G_transition
			      G_to_A_transition pyrimidine_transition
			      C_to_T_transition T_to_C_transition transversion
			      purine_to_pyrimidine_transversion
			      A_to_C_transversion A_to_T_transversion
			      G_to_C_transversion G_to_T_transversion
			      pyrimidine_to_purine_transversion
			      C_to_A_transversion C_to_G_transversion
			      T_to_A_transversion T_to_G_transversion
			      complex_substitution sequence_length_variation
			      simple_sequence_length_variation translocation);
map {$FEATURE_MAP{$_} = 'sequence_alteration'} @sequence_alterations;

# The ! to make it sort before any letters
my %ATT_ORDER = (ID              => '!01',
		 Variant_seq	 => '!02',
		 Reference_seq   => '!03',
		 Variant_effect  => '!04',
		 Genotype	 => '!05',
		 Variant_codon   => '!06',
		 Reference_codon => '!07',
		 Variant_aa	 => '!08',
		 Reference_aa	 => '!09',
		 Name		 => '!10',
		 Alias	   	 => '!11',
		 Parent	   	 => '!12',
		 Target	   	 => '!13',
		 Gap		 => '!14',
		 Derives_from	 => '!15',
		 Dbxref	         => '!16',
		 Ontology_term   => '!17',
		 Is_circular     => '!18',
    );

__PACKAGE__->load_components(qw/Core/);
__PACKAGE__->table('feature');
__PACKAGE__->add_columns(qw(feature_idx feature_id seqid source type
                            start end score strand phase bin));

__PACKAGE__->set_primary_key('feature_idx');
__PACKAGE__->has_many(attributes   => 'GAL::Schema::Result::Attribute',
                      'feature_idx');
__PACKAGE__->has_many(my_parents  => 'GAL::Schema::Result::Relationship',
                      {'foreign.child'   => 'self.feature_id'});
__PACKAGE__->has_many(my_children => 'GAL::Schema::Result::Relationship',
                      {'foreign.parent'  => 'self.feature_id'});
__PACKAGE__->many_to_many(parents  => 'my_parents', 'your_parents');
__PACKAGE__->many_to_many(children => 'my_children', 'your_children');

#-----------------------------------------------------------------------------

=head2 seqid

 Title   : seqid
 Usage   : $seqid = $self->seqid
 Function: Get the features sequence ID (i.e. what chromosome/contig it's on).
 Returns : The value of the sequence ID as text.
 Args    : None

=cut

#-----------------------------------------------------------------------------

=head2 source

 Title   : source
 Usage   : $source = $self->source
 Function: Get the features source.
 Returns : The value of the source as text.
 Args    : None

=cut

#-----------------------------------------------------------------------------

=head2 type

 Title   : type
 Usage   : $type = $self->type
 Function: Get the features type - a value constrained to be a child of the
	   SO's sequence feature term.
 Returns : The value of the type as text.
 Args    : None

=cut

#-----------------------------------------------------------------------------

=head2 start

 Title   : start
 Usage   : $start = $self->start
 Function: Get the features start on the sequence described by seqid.
 Returns : The value of the start an integer.
 Args    : None

=cut

#-----------------------------------------------------------------------------

=head2 end

 Title   : end
 Usage   : $end = $self->end
 Function: Get the features end on the sequence described by seqid.
 Returns : The value of the end as an integer.
 Args    : None

=cut

#-----------------------------------------------------------------------------

=head2 score

 Title   : score
 Usage   : $score = $self->score
 Function: Get the features score.
 Returns : The value of the score as text.
 Args    : None

=cut

#-----------------------------------------------------------------------------

=head2 strand

 Title   : strand
 Usage   : $strand = $self->strand
 Function: Get the features strand.
 Returns : The value of the strand as text.
 Args    : None

=cut

#-----------------------------------------------------------------------------

=head2 phase

 Title   : phase
 Usage   : $phase = $self->phase
 Function: Get the features phase.
 Returns : The value of the phase as text.
 Args    : None

=cut

#-----------------------------------------------------------------------------

=head2 attributes_hash

 Title   : attributes_hash
 Usage   : $self->attributes_hash
 Function: Return the attributes as a hash (or reference) with all values as
	   array references.  For consistency, even those values, such as ID,
	   that can have only one value are still returned as array
	   references.
 Returns : A hash or hash reference of attribute key/value pairs.
 Args    : None.

=cut

sub attributes_hash {
  my $self = shift;

  if (! $self->{my_attributes}) {
    map {push @{$self->{my_attributes}{$_->att_key}}, $_->att_value}
	   $self->attributes->all;
  }
  return wantarray ? %{$self->{my_attributes}} :
      dclone($self->{my_attributes});
}

#-----------------------------------------------------------------------------

=head2 attribute_exists

 Title   : attribute_exists
 Usage   : $self->attribute_exists($tag)
 Function: Return the value(s) of a particular attribute as an array or
	   reference.  Note that for consistency, values are always returned
	   as arrays (or reference) even in cases where only a single value
	   could exist such as ID.
 Returns : An array or reference of values.
 Args    : None

=cut

sub attribute_exists {
  my ($self, $tag) = @_;

  my $attributes = $self->attributes_hash;
  return 1 if exists $attributes->{$tag};
}

#-----------------------------------------------------------------------------

=head2 attribute_value

 Title   : attribute_value
 Usage   : $self->attribute_value($tag)
 Function: Return the value(s) of a particular attribute as an array or
	   reference.  Note that for consistency, values are always returned
	   as arrays (or reference) even in cases where only a single value
	   could exist such as ID.
 Returns : An array or reference of values.
 Args    : None

=cut

sub attribute_value {
  my ($self, $tag) = @_;

  my $attributes = $self->attributes_hash;
  my $values = [];
  if (exists $attributes->{$tag}) {
    $values = $attributes->{$tag};
  }
  return wantarray ? @{$values} : $values;
}

#-----------------------------------------------------------------------------

=head2 as_fasta

 Title   : as_fasta
 Usage   : $self->as_fasta
 Function: Returns a string of the feature formatted as a fasta
           sequence.
 Returns : A text string of the features sequence formated as
           fasta.
 Args    : None

=head2 seq

An alias for feature_seq

=cut

sub as_fasta {
  my $self = shift;
  my $header = '>' . $self->feature_id . ' ';
  $header .= 'locus='     . $self->locus  . ';';
  $header .= 'strand='    . $self->strand . ';';
  $header .= 'type='      . $self->type   . ';';
  $header .= 'parent='    . join(',', ($self->attribute_value('Parent'))) . ';' 
    if $self->attribute_exists('Parent');
  $header .= 'sib_count=' . join(':', ($self->sibling_part_of)) . ';'
    if $self->count_siblings;

  my $seq = $self->seq;
  $seq = $self->annotation->wrap_text($seq);
  my $fasta = "$header\n$seq";
  return $fasta;
}

#-----------------------------------------------------------------------------

=head2 feature_seq

 Title   : feature_seq
 Usage   : $self->feature_seq
 Function: Returns the features sequence as text 5' to 3' in the context of
	   the feature, so features on the reverse strand are reverse
	   complimented.
 Returns : A text string of the features sequence.
 Args    : None

=head2 seq

An alias for feature_seq

=cut

sub seq {return shift->feature_seq}

sub feature_seq {
  my $self = shift;

  my $seq = $self->annotation->fasta->seq($self->seqid,
					  $self->start,
					  $self->end
					 );
  $self->warn('feature_returned_no_sequence',
	      "Be sure you passed a fasta file and that this seqid exists" .
	      "If all else fails, delete the fasta sequence .index file "  .
	      "and try again.")
    if ! $seq;

  $seq = $self->annotation->revcomp($seq) if
    $self->strand eq '-';

  return $seq;
}

#-----------------------------------------------------------------------------

=head2 gc_content

 Title   : gc_content
 Usage   : $self->gc_content
 Function: Returns the ratio of the count of G or C nts over the sequence
	   length.
 Returns : A fractional value from 0 to 1
 Args    : None

=cut

sub gc_content {
  my $self = shift;

  my $seq = $self->seq;
  my $length = length $seq;
  my $gc = $seq =~ tr/GC/GC/;
  my $gc_ratio = $gc / $length;
  return $gc_ratio;
}

#-----------------------------------------------------------------------------

=head2 genomic_seq

 Title   : genomic_seq
 Usage   : $self->genomic_seq
 Function: Returns the features sequence as text 5' to 3' in the context of
	   the genome, so features on the reverse strand are not reverse
	   complimented.
 Returns : A text string of the features sequence on the genome.
 Args    : None

=cut

sub genomic_seq {
  my $self = shift;

  my $seq = $self->annotation->fasta->seq($self->seqid,
					  $self->start,
					  $self->end
					 );
  return $seq;
}

#-----------------------------------------------------------------------------

=head2 length

 Title   : length
 Usage   : $self->length
 Function: Returns the features length on the genome.  If a feature is
	   not contiguous on the genome (i.e. a transcript) then this method
	   will be (or should be) overridden by that subclass to provide the
	   'spliced' length of that feature.
 Returns : An integer.
 Args    : None

=cut

sub length {
  my $self = shift;

  return ($self->end - $self->start) + 1;
}

#-----------------------------------------------------------------------------

=head2 my_start

 Title   : my_start
 Usage   : $self->my_start
 Function: Returns the begin 5' coordinate on the feature's strand
	    (i.e. this is the same as calling start for a + strand
	    feature and end for a - strand feature.
 Returns : An integer.
 Args    : None

=cut

sub my_start {
  my $self = shift;
  return $self->strand eq '-' ? $self->end : $self->start;
}

#-----------------------------------------------------------------------------

=head2 my_end

 Title   : my_end
 Usage   : $self->my_end
 Function: Returns the begin 5' coordinate on the feature's strand
	    (i.e. this is the same as calling end for a + strand
	    feature and start for a - strand feature.
 Returns : An integer.
 Args    : None

=cut

sub my_end {
  my $self = shift;
  return $self->strand eq '-' ? $self->start : $self->end;
}

#-----------------------------------------------------------------------------

=head2 genomic_length

 Title   : genomic_length
 Usage   : $self->genomic_length
 Function: Returns the features genomic_length on the genome.  If a feature is
	   not contiguous on the genome (i.e. a transcript) then this method
	   will still provide the genomic length of (end - start + 1)
 Returns : An integer.
 Args    : None

=cut

sub genomic_length {
  my $self = shift;

  return ($self->end - $self->start) + 1;
}

#-----------------------------------------------------------------------------

=head2 annotation

 Title   : annotation
 Usage   : $self->annotation
 Function: Each feature has a weakened reference of the GAL::Annotation
	   object that created it and this method provides access to that
	   object.
 Returns : A GAL::Annotation object.
 Args    : None

=cut

sub annotation {
  return shift->result_source->schema->annotation;
}

#-----------------------------------------------------------------------------

=head2 throw

 Title   : throw
 Usage   : $self->throw
 Function: This is a convinience function that forwards throwings on to the
	   L<GAL::Annotation> object so that L<GAL::Base> can handle them.
	   You could get the same behaviour by doing
	   C<$self->annotation->throw(message => 'throwing'>.
 Returns : Nothing
 Args    : The following list: ($code, $message)

=cut

sub throw {
  return shift->annotation->throw(@_);
}

#-----------------------------------------------------------------------------

=head2 warn

 Title   : warn
 Usage   : $self->warn
 Function: This is a convinience function that forwards warnings on to the
	   L<GAL::Annotation> object so that L<GAL::Base> can handle them.
	   You could get the same behaviour by doing
	   C<$self->annotation->warn(message => 'warning'>.
 Returns : Nothing
 Args    : The following list: ($code, $message)

=cut

sub warn {
  return shift->annotation->warn(@_);
}

#-----------------------------------------------------------------------------

# This is a private method, but is not named as one due to contraints on
# the method name imposed by DBIx::Class.

sub inflate_result {
	my $self = shift;
	my $feature = $self->next::method(@_);

	my $subtype = $FEATURE_MAP{$feature->type} || 'sequence_feature';

	if ($subtype) {
		my $subclass = __PACKAGE__ . '::' . $subtype;
		$self->ensure_class_loaded($subclass);
		bless $feature, $subclass;
	}

	return $feature;
}

#-----------------------------------------------------------------------------

sub get_feature_bins {
    my $self = shift;
    $self->annotation->get_feature_bins($self);
}

#-----------------------------------------------------------------------------

=head2 as_gff3_recursive

 Title   : as_gff3_recursive
 Alias   : to_gff3_recursive
 Usage   : $self->as_gff3_recursive
 Function: Return a list of the string representation of this feature
	   and all of its part_of children in GFF3 format.
 Returns : An array or ref
 Args    : N/A

=cut

sub to_gff3_recursive {shift->as_gff3_recursive(@_)}

sub as_gff3_recursive {

  my $self = shift;

  my @gff_lines;
  my @features;

  $self->get_recursive_children(\@features);

  unshift @features, $self;

  @features = sort {($a->seqid cmp $b->seqid ||
		     $a->start <=> $b->start ||
		     $b->end   <=> $a->end)
		    } @features;

  for my $feature (@features) {
    push @gff_lines, $feature->as_gff3;
  }
  @gff_lines = uniq @gff_lines;
  return wantarray ? @gff_lines : \@gff_lines;
}

#-----------------------------------------------------------------------------

=head2 as_gff3

 Title   : as_gff3
 Alias   : to_gff3
 Usage   : $self->as_gff3
 Function: Return a string representation of this feature in GFF3 format.
 Returns : A string of GFF3 text
 Args    : N/A

=cut

sub to_gff3 {shift->as_gff3(@_)}

sub as_gff3 {

  my $self = shift;

  my $seqid      = $self->seqid;
  my $source     = $self->source;
  my $type       = $self->type;
  my $start      = $self->start;
  my $end        = $self->end;
  my $score      = $self->score;
  my $strand     = $self->strand;
  my $phase      = $self->phase;
  my $attributes = $self->attributes_hash;

  # Convert attribute values to comma sep list
  # ->attributes_hash uses dclone to prevent
  # mapping the actual attributes hash.
  map {$_ = join ',', @{$_}} values %{$attributes};

  my $attrb_text;

  for my $key (sort {($ATT_ORDER{$a} || $a) cmp ($ATT_ORDER{b} || $b)}
	       keys %{$attributes}) {
    $attrb_text .= "$key=" . $attributes->{$key} . ';';
  }

  $attrb_text ||= '';

  my $gff3_text = join "\t", ($seqid, $source, $type, $start, $end,
			      $score, $strand, $phase, $attrb_text);

  return $gff3_text;
}

#-----------------------------------------------------------------------------

=head2 get_recursive_children

 Title   : get_recursive_children
 Usage   : $self->get_recursive_children(\@children);
 Function: Collect all part_of children of an object recursively adding them to
           the provided array reference.
 Returns : Fills the supplied array reference.
 Args    : An array reference.

=cut

sub get_recursive_children {

  my ($self, $list) = @_;

  my $children = $self->children();

  while (my $child = $children->next) {
    push @{$list}, $child;
    $child->get_recursive_children($list);
  }
  return;
}

#-----------------------------------------------------------------------------

=head2 locus

 Title   : locus
 Usage   : $locus = $feature->locus;
 Function: Get the locus of a feature formatted as seqid:start-end
 Returns : A string formatted as seqid:start-end.
 Args    : None

=cut

sub locus {
  return sprintf('%s:%u-%u', shift->get_values(qw(seqid start end)));
}

#-----------------------------------------------------------------------------

=head2 siblings

 Title   : siblings
 Usage   : $sibs = $self->siblings;
 Function: Get an interator of all sibling features (from the first parent if
           there are multiple).
 Returns : A GAL::Schema::ResultSet object or a list of features in
           list context
 Args    : N/A

=cut

sub siblings {

  my $self = shift;

  my $parents = $self->parents;
  return $parents unless $parents->count;
  my $sibs = $parents->first->children->search(type => $self->type);
  return wantarray ? $sibs->all : $sibs;
}

#-----------------------------------------------------------------------------

=head2 sibling_part_of

 Title   : sibling_part_of
 Usage   : ($part, $of) = $self->sibling_part_of;
 Function: Return an array (ref) of the location in the list of sibs
           and the number of sibs that a given feature has.
 Returns : A GAL::Schema::ResultSet object or a list of features in
           list context
 Args    : N/A

=cut

sub sibling_part_of {

  my $self = shift;

  my $sibs  = $self->siblings;
  my $sibling_part_of = $sibs->count;
  my %part;
  my $part_count = 0;
  while (my $sib = $sibs->next) {
    $part_count++;
    last if $sib->feature_id eq $self->feature_id;
  }

  return wantarray ? ($part_count, $sibling_part_of) :
    [$part_count, $sibling_part_of];
}

#-----------------------------------------------------------------------------

=head2 count_parents

 Title   : count_parents
 Usage   : @values = $self->count_parents();
 Function: Count the number of parents.
 Returns : An integer as a scalar.
 Args    : N/A

=cut

sub count_parents {

  my $self = shift;

  return $self->parents->count;

}

#-----------------------------------------------------------------------------

=head2 count_siblings

 Title   : count_siblings
 Usage   : @values = $self->count_siblings();
 Function: Count the number of siblings.
 Returns : An integer as a scalar.
 Args    : N/A

=cut

sub count_siblings {

  my $self = shift;

  return $self->siblings->count;

}

#-----------------------------------------------------------------------------

=head2 count_children

 Title   : count_children
 Usage   : @values = $self->count_children();
 Function: Get the number of 1st degree children
 Returns : An array or ref
 Args    : N/A

=cut

sub count_children {

  my $self = shift;

  return $self->children->count;

}

#-----------------------------------------------------------------------------

=head2 get_children_IDs

 Title   : get_children_IDs
 Usage   : @IDs = $self->get_children_IDs();
 Function: Get a list of IDs of first degree children
 Returns : An array or ref
 Args    : N/A

=cut

sub get_children_IDs {

  my $self = shift;

  my @ids;
  my $children = $self->children;
  while (my $child = $children->next) {
    push @ids, $child->feature_id;
  }
  return wantarray ? @ids : \@ids;
}

#-----------------------------------------------------------------------------

=head2 get_values

 Title   : get_values
 Usage   : @values = $self->get_values(qw(seqid start end)
 Function: Get a list of values as an array(ref).
 Returns : An array or ref
 Args : A list of fields to retrieve (seqid source type start end
        score strand phase).

=cut

sub get_values {

  my ($self, @fields) = @_;

  my %allowed = map {$_ => 1}
    qw(seqid source type start end score strand phase);

  my @values;
  for my $field (@fields) {
    $self->throw('invalid_field_requested', $field)
      unless $allowed{$field};
    push @values, $self->$field;
  }
  return wantarray ? @values : \@values;
}

#-----------------------------------------------------------------------------

=head2 as_hash

 Title   : as_hash
 Usage   : $feature_hash = $self->as_hash
           %feature_hash = $self->as_hash;
 Function: Get feature values as a hash(ref)
 Returns : A hash or ref
 Args    : NA 

=cut

sub as_hash {

  my $self = shift;

  my %feature_hash;
  @feature_hash{qw(seqid source type start end score strand phase)} = 
      $self->get_values(qw(seqid source type start end score strand phase));

  my $attributes = $self->attributes_hash;

  $feature_hash{attributes} = $attributes;

  return wantarray ? %feature_hash : \%feature_hash;
}

#-----------------------------------------------------------------------------

=head2 overlaps

 Title   : overlaps
 Usage   : $rv = $feature->overlaps($another_feature);
 Function: Returns true if the two features overlap (share any
           nucleotide in common).
 Returns : 1 if the two features overlap, otherwise undef
 Args    : A GAL::Schema::Result::Feature object

=cut

sub overlaps {
    my ($self, $feature) = @_;
    return 1 if $feature->start <= $self->end && $feature->end >= $self->start;
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

#-----------------------------------------------------------------------------

=head1 DIAGNOSTICS


=over 4

L<GAL::Schema::Result::Feature> does not throw any warnings or errors
at this time.

=back

=head1 CONFIGURATION AND ENVIRONMENT

<GAL::Schema::Result::Feature> requires no configuration files or
environment variables.

=head1 DEPENDENCIES

<DBIx::Class>

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
