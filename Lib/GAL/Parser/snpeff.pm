package GAL::Parser::snpeff;

use strict;
use vars qw($VERSION);
$VERSION = 0.2.0;

use base qw(GAL::Parser::VCFv4_1);
use GAL::Reader::DelimitedLine;

=head1 NAME

GAL::Parser::snpeff - <One line description of module's purpose here>

=head1 VERSION

This document describes GAL::Parser::snpeff version 0.2.0

=head1 SYNOPSIS

     use GAL::Parser::snpeff;

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

=cut

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

=cut

#-----------------------------------------------------------------------------

sub parse_vcf_info_key_value {

  my ($self, $record, $key, $value) = @_;

  my %effects;
  my @all_effects;
  my %effect_map;
  my @var_effects;
#  my $ensembl_ids = [];
  my @ensembl_ids;

  if ($key eq 'EFF') {
    my @annotations = split /,/, $value;
    foreach my $annotation (@annotations) {
    #  my $effect =~ /EFF=(\w+)/g;
    #  push (@var_effects, $1);
    #  my @var_effects =~ s/\)$//; 
      my @values = split /\||\(/, $annotation;
      push (@var_effects, $values[0]);
      push (@ensembl_ids, $values[9]);
    }  

    $key = 'Variant_effect';

    %effect_map = (
      'INTERGENIC' => {
        effect => ['intergenic_variant'],
        feature => 'intergenic_region',
      },
      'UPSTREAM' => {
        effect => ['upstream_gene_variant'],
        feature => 'gene',
      },
      'UTR_5_PRIME' => {
        effect => ['5_prime_UTR_variant'],
        feature => 'five_prime_UTR',
      },
      'UTR_5_DELETED' => {
        effect => ['5_prime_UTR_truncation'],
        feature => 'five_prime_UTR',
      },
      'START_GAINED' => {
        effect => ['5_prime_UTR_premature_start_codon_gain_variant'],
        feature => 'five_prime_UTR',
      },
      'SPLICE_SITE_ACCEPTOR' => {
        effect => ['splice_acceptor_variant'],
        feature => 'splice_acceptor',
      }, 
      'SPLICE_SITE_DONOR' => {
        effect => ['splice_donor_variant'],
        feature => 'splice_donor',
      },
      'START_LOST' => {
        effect => ['start_lost'],
        feature => 'mRNA',
      },
      'SYNONYMOUS_START' => {
        effect => ['start_retained_variant'],
        feature => 'mRNA',
      },
      'CDS' => {
        effect => ['coding_sequence_variant'],
        feature => 'mRNA',
      },
      'GENE' => {
        effect => ['gene_variant'],
        feature => 'gene',
      },
      'TRANSCRIPT' => {
        effect => ['transcript_variant'],
        feature => 'transcript',
      },
      'EXON' => {
        effect => ['exon_variant'],
        feature => 'exon',
      },
      'EXON_DELETED' => {
        effect => ['exon_loss_variant'],
        feature => 'exon',
      },
      'NON_SYNONYMOUS_CODING' => {
        effect => ['missense_variant'],
        feature => 'mRNA',
      },
      'SYNONYMOUS_CODING' => {
        effect => ['synonymous_variant'],
        feature => 'mRNA',
      },
      'FRAME_SHIFT' => {
        effect => ['frameshift_variant'],
        feature => 'mRNA',
      },
      'CODON_CHANGE' => {
        effect => ['coding_sequence_variant'],
        feature => 'mRNA',
      },
      'CODON_INSERTION' => {
        effect => ['inframe_insertion'],
        feature => 'mRNA',
      },
      'CODON_CHANGE_PLUS_CODON_INSERTION' => {
        effect => ['disruptive_inframe_insertion'],
        feature => 'mRNA',
      },
      'CODON_DELETION' => {
        effect => ['inframe_deletion'],
        feature => 'mRNA',
      },
      'CODON_CHANGE_PLUS_CODON_DELETION' => {
        effect => ['disruptive_inframe_deletion'],
        feature => 'mRNA',
      },
      'STOP_GAINED' => {
        effect => ['stop_gained'],
        feature => 'mRNA',
      },
      'SYNONYMOUS_STOP' => {
        effect => ['stop_retained_variant'],
        feature => 'mRNA',
      },
      'STOP_LOST' => {
        effect => ['stop_lost'],
        feature => 'mRNA',
      },
      'INTRON' => {
        effect => ['intron_variant'],
        feature => 'intron',
      },
      'UTR_3_PRIME' => {
        effect => ['3_prime_UTR_variant'],
        feature => 'three_prime_UTR',
      },
      'UTR_3_DELETED' => {
        effect => ['3_prime_UTR_truncation'],
        feature => 'three_prime_UTR',
      },
      'DOWNSTREAM' => {
        effect => ['downstream_gene_variant'],
        feature => 'gene',
      }, 
      'INTRON_CONSERVED' => {
        effect => ['conserved_intron_variant'],
        feature => 'intron',
      },
      'INTERGENIC_CONSERVED' => {
        effect => ['conserved_intergenic_variant'],
        feature => 'intergenic_region',
      },
      'INTRAGENIC' => {
        effect => ['intragenic_variant'],
        feature => 'gene',
      }
    );

    my $so_var_effects;
    my $so_feature_type;

    if (scalar @ensembl_ids > 1) {
      if (scalar @ensembl_ids == scalar @var_effects) {
        for my $idx (0..$#ensembl_ids) {
          my $var_effect = $var_effects[$idx];
          my $ensembl_id = $ensembl_ids[$idx];
          if (exists $effect_map{$var_effect}) {
            $so_var_effects = $effect_map{$var_effect}{effect};
            $so_feature_type = $effect_map{$var_effect}{feature};
            for my $so_var_effect (@{$so_var_effects}) {
              push @{$effects{$so_var_effect}{$so_feature_type}}, $ensembl_id;
            }
          }
          else {
            $self->warn('unknown_effect_type', "$var_effect cannot be mapped to Sequence Ontology");
            push @{$effects{$var_effect}{'sequence_feature'}}, $ensembl_id;
          }
        }
      }
      else {
        my $locus = join':', (@{$record}{qw(CHROM POS ID)});
        $self->throw('effect_id_list_mismatch', $locus, "Different number of variant effects compared to accession IDs");
      }
    }
    else {
      for my $so_var_effect (@{$so_var_effects}) {
        push @{$effects{$so_var_effect}{$so_feature_type}}, @ensembl_ids;
      }
    }

=cut

    for my $effect (keys %effects) {
      for my $so_var_effect (keys %{$effects{$effect}}) {
        for my $so_feature_type (keys %{$effects{$so_var_effect}{$so_feature_type}}) {
          my $this_effect = "$so_var_effects . $so_feature_type ";
          $this_effect .= join ' ', @{$effects{$so_var_effect}{$so_feature_type}};
          push @all_effects, $this_effect;
        }
      }
    }

=cut

    for my $effect (keys %effects) {
      for my $so_feature_type (keys %{$effects{$effect}}) {
        my $this_effect = "$effect . $so_feature_type ";
          $this_effect .= join ' ', @{$effects{$effect}{$so_feature_type}};
          push @all_effects, $this_effect;
        }
      }
    
    $value = \@all_effects;
  }
  else {
    $key = 'vcf_' . $key;
    my @values = split /,/, $value;
    $value = \@values;
  }
  return $key, $value;
}

#--------------------------------------------------------------------------------- 

=cut

		my $variant_effect = join ' ', ($sequence_variant,
						$effect_hash{Genotype} - 1,
						$sequence_feature);
		push @variant_effects, $variant_effect;
	    }
		my $variant_effect   = join ' ', ($sequence_variant, $effect_hash{Genotype} - 1, $sequence_feature, $effect_hash{Gene_name});
		push @variant_effects, $variant_effect;
	    }
		my $variant_effect   = join ' ', ($sequence_variant, $effect_hash{Genotype} - 1, $sequence_feature, $effect_hash{Transcript_ID});
		push @variant_effects, $variant_effect;
	    }
		my $variant_effect   = join ' ', ($sequence_variant, $effect_hash{Genotype} - 1, $sequence_feature, $effect_hash{Transcript_ID});
		push @variant_effects, $variant_effect;
		
	    }
    my @all_seqs = ($reference_seq, @variant_seqs);
    my @lengths = map {length($_)} @all_seqs;

    #INDEL
    if (grep {$_ != 1} @lengths) {
	my ($min_length) = sort {$a <=> $b} @lengths;
	map {substr($_, 0, $min_length, '')} @all_seqs;
	map {$_ ||= '-'} @all_seqs;
	my $start_adjust = length($reference_seq) - length($all_seqs[0]);
	$start += $start_adjust;
	my $end_adjust = length($all_seqs[0]) - 1;
	$end = $start + $end_adjust;
	
	$reference_seq = shift @all_seqs;
	@variant_seqs = @all_seqs;
	
	if ($reference_seq eq '-') {
	    $type = 'insertion';
	}
	elsif (grep {$_ eq '-'} @variant_seqs) {
	    $type = 'deletion';
	}
	else {
	    $type = 'complex_substitution';
	}
    }
    
    my $read_count = $info{DP}[0] || undef;
    
    # Make sure that at least one variant sequence differs from
    # the reference sequence.
    next unless grep {$_ ne $reference_seq} @variant_seqs;
    
    # Format Column Keys
    # GT : genotype, encoded as alleles values separated by either of
    # ”/” or “|”, e.g. The allele values are 0 for the reference
    # allele (what is in the reference sequence), 1 for the first
    # allele listed in ALT, 2 for the second allele list in ALT and
    # so on. For diploid calls examples could be 0/1 or 1|0 etc. For
    # haploid calls, e.g. on Y, male X, mitochondrion, only one
    # allele value should be given. All samples must have GT call
    # information; if a call cannot be made for a sample at a given
    # locus, ”.” must be specified for each missing allele in the GT
    # field (for example ./. for a diploid). The meanings of the
    # separators are:
    #       / : genotype unphased
    #       | : genotype phased
    # DP : read depth at this position for this sample (Integer)
    # FT : sample genotype filter indicating if this genotype was
    # “called” (similar in concept to the FILTER field). Again, use
    # PASS : to indicate that all filters have been passed, a
    # semi-colon separated list of codes for filters that fail, or
    # ”.” to indicate that filters have not been applied. These
    # values should be described in the meta-information in the same
    # way as FILTERs (Alphanumeric String)
    # GL : three floating point log10-scaled likelihoods for
    # AA,AB,BB genotypes where A=ref and B=alt; not applicable if
    # site is not biallelic. For example: GT:GL
    # 0/1:-323.03,-99.29,-802.53 (Numeric)
    # GQ : genotype quality, encoded as a phred quality
    # -10log_10p(genotype call is wrong) (Numeric)
    # HQ : haplotype qualities, two phred qualities comma separated
    # (Numeric)
    
    my $feature_id = join ':', ($seqid, $source, $start);
    
    my $attributes = {Reference_seq => [$reference_seq],
		      Variant_seq   => \@variant_seqs,
		      ID            => [$feature_id],
    };
    $attributes->{Total_reads} = [$read_count] if $read_count;
    $attributes->{Variant_effect} = \@variant_effects
	if scalar @variant_effects;
    
    my $feature = {feature_id => $feature_id,
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
     return $feature;
}

#-----------------------------------------------------------------------------

sub reader {
    my $self = shift;


    if (! $self->{reader}) {
	# We don't give field names here so that the record will be returned
	# as an array rather than a hash.  We splice out the first
	# 9 columns above and the rest in genotype data.
	# my @field_names = qw(chrom pos id ref alt qual filter info format);
	my $reader = GAL::Reader::DelimitedLine->new(comment_pattern => qr/^\#\#/,
						     header_count    => 1
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

<GAL::Parser::snpeff> requires no configuration files or environment variables.

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
