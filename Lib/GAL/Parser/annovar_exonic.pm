package GAL::Parser::annovar_exonic;

use strict;
use vars qw($VERSION);
$VERSION = 0.2.0;

use base qw(GAL::Parser::VCFv4_1);
use GAL::Reader::DelimitedLine;

=head1 NAME

GAL::Parser::annovar_exonic - Parse ANNOVAR_EXONIC files

=head1 VERSION

This document describes GAL::Parser::annovar_exonic version 0.2.0

=head1 SYNOPSIS

    my $parser = GAL::Parser::annovar_exonic->new(file => 'annovar_exonic.txt');

    while (my $feature_hash = $parser->next_feature_hash) {
	print $parser->to_gff3($feature_hash) . "\n";
    }

=head1 DESCRIPTION

L<GAL::Parser::annovar_exonic> provides a parser for ANNOVAR_EXONIC data.

=head1 Constructor

New L<GAL::Parser::annovar_exonic> objects are created by the class method
new.  Arguments should be passed to the constructor as a list (or
reference) of key value pairs.  All attributes of the
L<GAL::Parser::annovar_exonic> object can be set in the call to new. An
simple example of object creation would look like this:

    my $parser = GAL::Parser::annovar_exonic->new(file => 'annovar_exonic.txt');

The constructor recognizes the following parameters which will set the
appropriate attributes:

=item * C<< file => feature_file.txt >>

This optional parameter provides the filename for the file containing
the data to be parsed. While this parameter is optional either it, or
the following fh parameter must be set.

=item * C<< fh => feature_file.txt >>

This optional parameter provides a filehandle to read data from. While
this parameter is optional either it, or the following fh parameter
must be set.

=cut

#-----------------------------------------------------------------------------

sub parse_vcf_info_key_value {
  
  my ($self, $record, $key, $value) = @_;
  my %feature_map = (Transcript => 'transcript');
  my %allele_map;
  my @alleles = split /,/, $record->{alt};       

    my $reference_seq = $record->{ref_seq};
    my $variant_seq = $record->{var_seq};

    my $effect_type = $record->{effect_type};
    my ($effect, $type);
    if ($effect_type eq 'frameshift insertion') {
	$effect = 'frameshift_elongation';
	$type   = 'insertion';
    }
    elsif ($effect_type eq 'frameshift deletion') {
	$effect = 'frameshift_truncation';
	$type   = 'deletion';
    }
   elsif ($effect_type eq 'frameshift block substitution') { 
	  $effect = 'frameshift_variant';
	  $type = 'complex_substitution';
    }
   elsif ($effect_type eq 'stopgain SNV') {
          $effect = 'stop_gained';
	  $type = 'SNV';
    }
   elsif ($effect_type eq 'stoploss SNV') {
 	  $effect = 'stop_lost';
	  $type = 'SNV';
    }
   elsif ($effect_type eq 'nonframeshift insertion') {
	  $effect = 'inframe_insertion';
	  $type = 'insertion';
    }
   elsif ($effect_type eq 'nonframeshift deletion') { 
	  $effect = 'inframe_deletion';
	  $type = 'deletion';
    }
   elsif ($effect_type eq 'nonframeshift block substitution') {
 	  $effect = 'inframe_variant';
	if (length($reference_seq) == length($variant_seq)) {
		$type = 'MNP';
         }
	else {
		$type = 'complex_substitution';
	     }
    }
   elsif ($effect_type eq 'nonsynonymous SNV') {
  	  $effect = 'missense_variant';
	  $type = 'SNV';
    }
   elsif ($effect_type eq 'synonymous SNV') {
 	  $effect = 'synonymous_variant';
	  $type = 'SNV';
    }
  elsif ($effect_type eq 'unknown') {
	 $effect = 'sequence_variant';
	 $type = 'sequence_alteration';
    }
  else {
	$self->warn('unknown_annovar_effect_type', $effect_type)
   }

#    $effect = exists $effect_map{$effect} ? $effect_map{$effect} : 'sequence_variant';

    my $type_name  = uc substr($type, 0, 3);
    my $feature_id = sprintf("$type_name%06d", $self->counter($type));
			     
    # Create the attributes hash
    # See http://www.sequenceontology.org/gvf.html
    
    # Assign the reference and variant sequences:
    # Reference_seq=A
    # The order of the Variant_seqs matters for indexing Variant_effect
    # - see below
    # Variant_seq=G,A
    my $zygosity = $record->{zygosity};
    my @variant_seqs;
    if ($zygosity eq 'het') {
	@variant_seqs = ($reference_seq, $variant_seq);
	$zygosity = 'heterozygous';
    }
    elsif ($zygosity eq 'hom') {
	@variant_seqs = ($variant_seq);
	$zygosity = 'homozygous';
    }
    elsif ($zygosity eq 'unknown') {
	@variant_seqs = ($variant_seq);    
	$zygosity = 'homozygous';
	$self->warn('unknown_zygosity_defaults_to_homozygous', "$type $seqid:$start-$end");
    }
    else {
	$self->throw('unknown_zygosity', "$seqid:$start-$end has zygosity=$zygosity");
    }
    
    # For sequence_alteration features the suggested keys include:
    # ID, Reference_seq, Variant_seq, Variant_reads Total_reads,
    # Zygosity, Intersected_feature, Variant_effect, Copy_number
    
    my $variant_effect = $record->{effect};

    my $attributes = {Reference_seq  => [$reference_seq],
		      Variant_seq    => \@variant_seqs,
		      Variant_effect => [$effect],
		      Zygosity       => [$zygosity],
		      ID             => [$feature_id],
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
 Returns : A GAL::Reader::DelimitedLine singleton.
 Args    : None

=cut

sub reader {
  my $self = shift;

  if (! $self->{reader}) {
      
      # Example of Annovar exonic_variation output
      # line_number  line27
      # effect       frameshift
      # type         insertion
      # feature_data SAMD11:NM_152486:exon8:c.831_832insCCCT:p.E277fs,
      # chrom        1
      # start        866511
      # end          866511
      # ref_seq      -
      # var_seq      CCCT
      # zygosity     hom
      # score        2553.96
      # unknown_1    62
      # unknown_2    41.19

    my @field_names = qw(line_number effect_type feature_data chrom start end ref_seq var_seq zygosity score unknown_1 unknown_2);
    my $reader = GAL::Reader::DelimitedLine->new(field_names     => \@field_names,
						 field_separator => "\t");
    $self->{reader} = $reader;
  }
  return $self->{reader};
}

#-----------------------------------------------------------------------------

=head1 DIAGNOSTICS

L<GAL::Parser::annovar_exonic> does not throw any warnings or errors.

=head1 CONFIGURATION AND ENVIRONMENT

L<GAL::Parser::annovar_exonic> requires no configuration files or
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

Copyright (c) 2012, Barry Moore <barry.moore@genetics.utah.edu>.  All
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
