package GAL::Schema::Result::Feature::sequence_alteration;

use strict;
use warnings;
use base qw(GAL::Schema::Result::Feature::sequence_feature);

=head1 NAME

GAL::Schema::Result::Feature::sequence_alteration - A
sequence_alteration object for the GAL library

=head1 VERSION

This document describes
GAL::Schema::Result::Feature::sequence_alteration version 0.2.0

=head1 SYNOPSIS

    use GAL::Annotation
    my $var_store = GAL::Annotation->new(storage => $var_store_args,
					 parser  => $parser_args,
					 fasta   => $fasta_args,
					);

      $var_store->load_files(files => $variant_file,
			     mode  => 'overwrite',
			    );
    my $variants = $var_store->schema->resultset('Feature');
    my $snvs = $variants->search(\%where);

    while (my @snv = $snvs->next) {
      my $snv_id = $snv->feature_id;
      my $start  = $snv->start;
      my $reference_seq = $snv->reference_seq;
      my @variant_seqs  = $snv->variant_seqs;
    }

=head1 DESCRIPTION

<GAL::Schema::Result::Feature::sequence_alteration> provides a
<GAL::Schema::Result::Feature> subclass for sequence_alteration
specific behavior.  See documentation for
<GAL::Schema::Result::Feature> for more methods and details.

=head1 METHODS

=cut

#-----------------------------------------------------------------------------

=head2 variant_seq

 Title   : variant_seq
 Usage   : $variant_seq = $self->variant_seq
 Function: Get the value(s) for the features Variant_seq attribute.
 Returns : An array or reference of sequences.
 Args    : None

=cut

sub variant_seq {
  my $self = shift;
  my $variant_seqs = $self->attribute_value('Variant_seq');
  return wantarray ? @{$variant_seqs} : $variant_seqs;
}

#-----------------------------------------------------------------------------

=head2 variant_seq_no_ref

 Title   : variant_seq_no_ref
 Usage   : $variant_seq_no_ref = $self->variant_seq_no_ref
 Function: Get the value(s) for the features Variant_seq attribute, but
           remove the reference sequence if it is present.
 Returns : An array or reference of sequences.
 Args    : None

=cut

sub variant_seq_no_ref {
  my $self = shift;
  my $reference_seq = $self->reference_seq;
  my @variant_seqs_no_ref = grep {$_ ne $reference_seq}
    $self->attribute_value('Variant_seq');
  return wantarray ? @variant_seqs_no_ref : \@variant_seqs_no_ref;
}

#-----------------------------------------------------------------------------

=head2 reference_seq

 Title   : reference_seq
 Usage   : $reference_seq = $self->reference_seq
 Function: Get the value(s) for the features Reference_seq attribute.
 Returns : An array or reference of sequences.
 Args    : None

=cut

sub reference_seq {
  my $self = shift;
  my $reference_seq = $self->attribute_value('Reference_seq');
  return $reference_seq->[0];
}

#-----------------------------------------------------------------------------

=head2 Variant_reads

 Title   : Variant_reads
 Usage   : $Variant_reads = $self->Variant_reads
 Function: Get the value(s) for the features Variant_reads attribute.
 Returns : An array or reference of sequences.
 Args    : None

=cut

sub variant_reads {
  my $self = shift;
  my $variant_reads = $self->attribute_value('Variant_reads');
  return wantarray ? @{$variant_reads} : $variant_reads;
}

#-----------------------------------------------------------------------------

=head2 total_reads

 Title   : total_reads
 Usage   : $total_reads = $self->total_reads
 Function: Get the value(s) for the features Total_reads attribute.
 Returns : An array or reference of sequences.
 Args    : None

=cut

sub total_reads {
  my $self = shift;
  my $total_reads = $self->attribute_value('Total_reads');
  return wantarray ? @{$total_reads} : $total_reads;
}

#-----------------------------------------------------------------------------

=head2 genotype

 Title   : genotype
 Usage   : $genotype = $self->genotype
 Function: Get the value(s) for the features Genotype attribute.
 Returns : An array or reference of sequences.
 Args    : None

=cut

sub genotype {
  my $self = shift;
  my $genotype = $self->attribute_value('Genotype');
  return $genotype->[0];
}

#-----------------------------------------------------------------------------

=head2 variant_effect

 Title   : variant_effect
 Usage   : $variant_effect = $self->variant_effect
 Function: Get the value(s) for the features Variant_effect attribute.
 Returns : An array or reference of sequences.
 Args    : None

variant_effect describes the effect of that sequence variation.

=cut

sub variant_effect {
  my $self = shift;
  my $variant_effect = $self->attribute_value('Variant_effect');
  return wantarray ? @{$variant_effect} : $variant_effect;
}

#-----------------------------------------------------------------------------

=head2 variant_copy_number

 Title   : variant_copy_number
 Usage   : $variant_copy_number = $self->variant_copy_number
 Function: Get the value(s) for the features Variant_copy_number attribute.
 Returns : An array or reference of sequences.
 Args    : None

"A copy number variation (CNV) is when the number of copies of a particular gene varies from one individual to the next. Following the completion of the Human Genome Project, it became apparent that the genome experiences gains and losses of genetic material." (NHGRI, 2014)

=cut

sub variant_copy_number {
  my $self = shift;
  my $variant_copy_number = $self->attribute_value('Variant_copy_number');
  return wantarray ? @{$variant_copy_number} : $variant_copy_number;
}

#-----------------------------------------------------------------------------

=head2 reference_copy_number

 Title   : reference_copy_number
 Usage   : $reference_copy_number = $self->reference_copy_number
 Function: Get the value(s) for the features Reference_copy_number attribute.
 Returns : An array or reference of sequences.
 Args    : None

=cut

sub reference_copy_number {
  my $self = shift;
  my $reference_copy_number = $self->attribute_value('Reference_copy_number');
  return $reference_copy_number->[0];
}

#-----------------------------------------------------------------------------

=head1 DIAGNOSTICS

<GAL::Schema::Result::Feature::sequence_alteration> throws no warnings
or errors.

=head1 CONFIGURATION AND ENVIRONMENT

<GAL::Schema::Result::Feature::sequence_alteration> requires no
configuration files or environment variables.

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
