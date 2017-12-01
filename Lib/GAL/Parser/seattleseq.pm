package GAL::Parser::seattleseq;

use strict;
use vars qw($VERSION);
$VERSION = 0.2.0;

use base qw(GAL::Parser::VCFv4_1);
use GAL::Reader::DelimitedLine;

=head1 NAME

GAL::Parser::SeattleSeq - <One line description of module's purpose here>

=head1 VERSION

This document describes GAL::Parser::SeattleSeq version 0.2.0

=head1 SYNOPSIS

     use GAL::Parser::SeattleSeq;

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

sub parse_vcf_info {

  my ($self, $record) = @_;

  my %info;
  my @pairs = split /;/, $record->{info};
  for my $pair (@pairs) {
    my ($key, $value) = split /=/, $pair;
    $value = defined $value ? $value : '';
    ($key, $value) = $self->parse_vcf_info_key_value($record, $key, $value);
    #my %order = (GM => 'A', FG => 'B');
    #for my $key (sort { $order{$a} cmp $order{$b} } keys %order) {
    #  next;
    #} 
    push @{$info{$key}}, @{$value};
  }
  return wantarray ? %info : \%info;
}

#-----------------------------------------------------------------------------

sub parse_vcf_info_key_value {
  
  my ($self, $record, $key, $value) = @_;

  my %effects;
  my @all_effects;
  my %effect_map;  
  my @var_effects;
  my $gm_ids =[];

#  my @alleles = split /,/, $record->{alt};
 
  if ($key eq 'FG') {
    #$gm_ids = $record->{info}{GM};
    $key = 'Variant_effect';
  
    @var_effects = split /,/, $value;

    %effect_map = (
      'frameshift-near-splice' 	=> {
        effect => [qw(frameshift_variant splice_region_variant)],
        feature => 'mRNA',
      },
      'frameshift' => {
        effect => ['frameshift_variant'],
	feature => 'mRNA',
      },
      'codingComplex-near-splice' => {
	effect => [qw(coding_sequence_variant splice_region_variant)],
	feature => 'mRNA',
      },	
      'codingComplex' => {
	effect => ['coding_sequence_variant'],
	feature => 'mRNA',
      },	
      'coding-near-splice' => {
	effect => ['splice_region_variant'],
	feature => 'mRNA',
      },
      'coding' => {
	effect => ['coding_sequence_variant'],
	feature => 'mRNA',
      },
      'splice-5' => {
	effect => ['splice_donor_variant'],
	feature => 'mRNA',
      },
      'splice-3' => {
	effect => ['splice_acceptor_variant'],
	feature => 'mRNA',
      },
      'utr-5' => {
	effect => ['5_prime_UTR_variant'],
	feature => 'mRNA',
      },
      'utr-3' => {
	effect => ['3_prime_UTR_variant'],
	feature => 'mRNA',
      },
      'intron' => {
	effect => ['intron_variant'],
	feature => '',
      },
      'near-gene-5' => {
	effect => ['2KB_upstream_variant'],
	feature => '',
      },
      'near-gene-3' => {
	effect => ['500B_downstream_variant'],
	feature => '',
      },
      'intergenic' => {
	effect => ['intergenic_variant'],
	feature => '',
      },
      'missense' => {
       effect => ['missense_variant'],
       feature => '',
      },
      'stop-gained' => {
       effect => ['stop_gained'],
       feature => '',
      },
      'stop-lost' => {
       effect => ['stop_lost'],
       feature => '',
      },
      'coding-notMod3' => {
       effect => ['frameshift_variant'],
       feature => 'mRNA',
      },
      'coding-notMod3-near-splice' => {
       effect => ['frameshift_variant'],
       feature => 'mRNA',
      },
      'coding-synonymous' => {
       effect => ['synonymous_variant'],
       feature => 'mRNA',
      },
      'coding-synonymous-near-splice' => {
       effect => [qw(synonymous_variant splice_region_variant)], 
       feature => 'mRNA',
      },
      'missense-near-splice' => {
       effect => [qw(missense_variant splice_region_variant)],
       feature => 'mRNA',
      },
      'stop-gained-near-splice' => {
       effect => [qw(stop_lost splice_region_variant)],	
       feature => 'mRNA',
      },	
      'stop-lost-near-splice' => {
       effect => [qw(stop_lost splice_region_variant)],	
       feature => 'mRNA',
      }
    );
    
    my $so_var_effects;
    my $so_feature_type;

    if (scalar @{$gm_ids} > 1) {
      if (scalar @{$gm_ids} == scalar @var_effects) { 
        for my $idx (0..$#{$gm_ids}) {
          my $var_effect = $var_effects[$idx];
          my $gm_id = $gm_ids->[$idx];
          if (exists $effect_map{$var_effect}) { 
            $so_var_effects = $effect_map{$var_effect}{effect};
            $so_feature_type = $effect_map{$var_effect}{feature};
            for my $so_var_effect (@{$so_var_effects}) {
	      push @{$effects{$so_var_effect}{$so_feature_type}}, $gm_id;
            }
          } 
          else {
            $self->warn('unknown_effect_type', "$var_effect cannot be mapped to Sequence Ontology"); 
	    push @{$effects{$var_effect}{'sequence_feature'}}, $gm_id;
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
        push @{$effects{$so_var_effect}{$so_feature_type}}, $gm_ids;
      }
    }

    for my $effect (keys %effects) {
      for my $so_var_effect (keys %{$effects{$effect}}) {
        for my $so_feature_type (keys %{$effects{$effect}{$so_var_effect}}) {
          my $this_effect = "$so_var_effect . $so_feature_type ";
          $this_effect .= join ' ', @{$effects{$effect}{$so_var_effect}{$so_feature_type}};
          push @all_effects, $this_effect;
        }
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
   my @effect_list;
    for my $effect (keys %variant_effects) {
	my $effect_text = join ' ',
	($effect, '.', 'transcript',
	 @{$variant_effects{$effect}});
	push @effect_list, $effect_text;
    }

    my @format_order = split /:/, $record{format};

    my %individual_data;
    my ($header) = $self->reader->headers;
    my @header_data = split /\s+/, $header;
    splice(@header_data, 0, 9);
    @individual_data{@header_data} = @{$data};

    my @features;
    for my $individual_id (sort keys %individual_data) {

	my @individual_values = split /:/, $individual_data{$individual_id};
	my %individual_hash;
	map {$individual_hash{$_} = shift @individual_values} @format_order;

	my ($var1, $phased, $var2) =
	    split //, $individual_hash{GT};
	my $this_ref = $reference_seq;

	my %variant_seq_hash = map {$_ => 1} @all_seqs[$var1, $var2];
	my @variant_seqs = keys %variant_seq_hash;

	my @these_all_seqs = ($this_ref, @variant_seqs);
	my @lengths = map {length($_)} @these_all_seqs;

	my $this_start = $start;
	my $this_end   = $end;
	#INDEL
	if (grep {$_ != 1} @lengths) {
	  my ($min_length) = sort {$a <=> $b} @lengths;
	  map {substr($_, 0, $min_length, '')} @these_all_seqs;
	  map {$_ ||= '-'} @these_all_seqs;
	  my $start_adjust = length($reference_seq) - length($these_all_seqs[0]);
	  $this_start += $start_adjust;
	  my $end_adjust = length($these_all_seqs[0]) - 1;
	  $this_end = $this_start + $end_adjust;

	  $this_ref = shift @these_all_seqs;
	  @variant_seqs = @these_all_seqs;

	  if ($this_ref eq '-') {
	    $type = 'insertion';
	  }
	  elsif (grep {$_ eq '-'} @variant_seqs) {
	    $type = 'deletion';
	  }
	  else {
	    $type = 'complex_substitution';
	  }
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

<GAL::Parser::SeattleSeq> requires no configuration files or environment variables.

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

Copyright (c) 2012, Barry Moore <barry.moore@genetics.utah.edu>.  All rights reserved.

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
