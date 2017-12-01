package GAL::Parser::vep;

use strict;
use vars qw($VERSION);
$VERSION = 0.2.0;

use base qw(GAL::Parser::VCFv4_1);
use GAL::Reader::DelimitedLine;

=head1 NAME

GAL::Parser::vep - <One line description of module's purpose here>

=head1 VERSION

This document describes GAL::Parser::vep version 0.2.0

=head1 SYNOPSIS

     use GAL::Parser::vep;

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


sub parse_vcf_info_key_value {

  my ($self, $record, $key, $value) = @_;

  my %feature_map = (Transcript => 'transcript');

  my %allele_map;

  my @alleles = split /,/, $record->{alt};

  if ($key eq 'CSQ') {
    my @values = split /\|/, $value;
    $key   = 'Variant_effect';

    my %effects;
    my @all_effects;

    my ($allele, $gene, $feature_id_txt, $feature_type, $conseq_txt,
	$cDNA_pos, $cds_pos, $protein_pos, $amino_acids, $codons,
	$existing) = splice(@values, 0, 10);

    my @conseqs = split /&/, $conseq_txt;

    for my $conseq (@conseqs) {
      $feature_type = (exists $feature_map{$feature_type} ?
		       $feature_map{$feature_type}        :
		       $feature_type);

      my @feature_ids = split /&/, $feature_id_txt;
      push @{$effects{$conseq}{$allele}{$feature_type}}, @feature_ids
	if scalar @feature_ids;
    }
    for my $effect (keys %effects) {
      for my $allele (keys %{$effects{$effect}}) {
	for my $feature_type (keys %{$effects{$effect}{$allele}}) {
	  my $this_effect = "$effect $allele $feature_type ";
	  $this_effect .= join ' ', @{$effects{$effect}{$allele}{$feature_type}};
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

<GAL::Parser::vep> requires no configuration files or environment variables.

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
