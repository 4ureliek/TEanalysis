package GAL::Parser::dbsnp_flat;

use strict;
use vars qw($VERSION);


$VERSION = 0.2.0;
use base qw(GAL::Parser);

=head1 NAME

GAL::Parser::template - <One line description of module's purpose here>

=head1 VERSION

This document describes GAL::Parser::template version 0.2.0

=head1 SYNOPSIS

     use GAL::Parser::template;

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
     Usage   : GAL::Parser::template->new();
     Function: Creates a GAL::Parser::template object;
     Returns : A GAL::Parser::template object
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
	my @valid_attributes = qw(assembly); # Set valid class attributes here
	$self->set_attributes($args, @valid_attributes);
	######################################################################

	# Set the column headers from your incoming data file here
	# These will become the keys in your $record hash reference below.
	$self->record_separator("\n\n");
	$self->field_separator(undef);
	$self->fields([qw(data)]);

}

#-----------------------------------------------------------------------------

=head2 assembly

 Title   : assembly
 Usage   : $assembly = $self->assembly();
 Function: Get/set the assembly.  For these dbSNP records, each variant is 
           mapped to multiple assemblies - here you specify which assembly
           you want to use.
 Returns : A text string representing the assembly (i.e. reference HuRef
           etc. for human).
 Args    : A text string representing the assembly (i.e. reference HuRef
           etc. for human).

=cut

sub assembly {
	my ($self, $assembly) = @_;

	$self->{assembly} = $assembly if $assembly;

	return $self->{assembly};
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

	$self->throw(message => 'This module needs testing an verification before use.  In particular Variant_seq may be on the wrong strand');

	my ($record, $errors) = $self->_parse_dbsnp_flat($record->{data});

	return undef unless $record;

	my $id         = $record->{id};
	my $source     = 'dbSNP';
	my $type       = $record->{type};
	my $start      = $record->{chr_start};
	my $end        = $record->{chr_end};
	my $score      = '.';
	my $strand     = '+';
	my $phase      = '.';

	my @error_codes;
	for my $error (@{$errors}) {
		push @error_codes, $error->[0];
		$self->warn(@{$error}[0,1]);
	}

	my @ctgs = @{$record->{ctgs}};

	my @features;
	for my $ctg (@ctgs) {

		my @seqs = @{$record->{seqs}};

		my @ctg_error_codes;
		for my $error (@{$ctg->{errors}}) {
			push @ctg_error_codes, $error->[0];
			my $error_message = join "\t", ('ERROR',
							$id,
							$error->[0],
							$error->[1],
						       );
			  $self->warn(message => $error_message);
		}


		my $attributes = {ID            => $id,
				  Variant_seq   => \@seqs,
				  Errors        => [@error_codes, @ctg_error_codes],
				 };

		delete $attributes->{Errors} unless scalar @{$attributes->{Errors}};

		my $feature = {feature_id => $id,
			       seqid      => 'chr' . $ctg->{chr},
			       source     => $source,
			       type       => $type,
			       start      => $ctg->{chr_start},
			       end        => $ctg->{chr_end},
			       score      => $score,
			       strand     => $strand,
			       phase      => $phase,
			       attributes => $attributes,
			      };

		push @features, $feature;
	}

	return \@features;
}

#-----------------------------------------------------------------------------

=head2 _parse_dbsnp_flat

 Title   : _parse_dbsnp_flat
 Usage   : $record = $self->_parse_dbsnp_flat($record->{data});
 Function: Parses one record from a dbSNP ASN1 flatfile format.
 Returns : Returns a hashref in a format more typical of GAL/Parser subclasses
 Args    : The data from one record in a dbSNP ASN1 flat file

=cut

sub _parse_dbsnp_flat {
	my ($self, $text_data) = @_;


	#Skip junk that's not really an rs record.
	return undef if $text_data !~ /^rs\d+\s+\|\s+/;

	my %record;
	my @errors;

	my @lines = split /\n/, $text_data;

	my $selected_assembly = $self->assembly;

      LINE:
	for my $line (@lines) {
		my ($class, @data) = split /\s+\|\s+/, $line;
		next LINE if $class =~ /^ss\d+$/;

		if ($class eq 'CTG') {

			my ($assembly, $chr, $chr_start, $ctg, $ctg_start, $ctg_end, $loctype, $orient) = @data;
			$assembly =~ s/assembly=//    or push @errors, ['CTG_invalid_assemly_tag', $line];
			next LINE unless $assembly eq $selected_assembly;

			my @ctg_errors;
			$chr         =~ s/chr=//       or push @errors, ['CTG_invalid_chr_tag', $line];
			$chr_start   =~ s/chr-pos=//   or push @errors, ['CTG_invalid_chr_pos_tag', $line];
			$ctg_start   =~ s/ctg-start=// or push @errors, ['CTG_invalid_ctg-start_tag', $line];
			$ctg_end     =~ s/ctg-end=//   or push @errors, ['CTG_invalid_ctg-end_tag', $line];
			$loctype     =~ s/loctype=//   or push @errors, ['CTG_invalid_loctype_tag', $line];
			$orient      =~ s/orient=//    or push @errors, ['CTG_invalid_orient_tag', $line];

			my $chr_end;
			if ($chr_start !~ /^[0-9]+$/) {
				$chr_start = '.';
				push @ctg_errors, ['CTG_invalid_chr-pos_value', $line];
			}
			elsif ($ctg_start =~ /^[0-9]+$/ &&
			       $ctg_end   =~ /^[0-9]+$/
			      ) {
				$chr_end = $chr_start + abs($ctg_start - $ctg_end);
			}
			else {
				push @ctg_errors, ['CTG_invalid_ctg-start/end_value', $line];
				$chr_end = '.';
			}

			my %ctg = (assembly  => $assembly,
				   chr       => $chr,
				   chr_start => $chr_start,
				   chr_end   => $chr_end,
				   ctg       => $ctg,
				   ctg_start => $ctg_start,
				   ctg_end   => $ctg_end,
				   loctype   => $loctype,
				   orient    => $orient,
				   errors    => \@ctg_errors,
				  );

			push @{$record{ctgs}}, \%ctg;

			# CTG | assembly=Celera | chr=1 | chr-pos=19193787 | NW_927841.1 | ctg-start=3670230 | ctg-end=3670230 | loctype=2 | orient=+
			# CTG | assembly=HuRef | chr=1 | chr-pos=19115539 | NW_001838573.1 | ctg-start=2491270 | ctg-end=2491270 | loctype=2 | orient=+
			# CTG | assembly=reference | chr=1 | chr-pos=20742048 | NT_004610.18 | ctg-start=3693803 | ctg-end=3693803 | loctype=2 | orient=+

			# BUT NOTE multiple locations on the reference assembly!!!!!
			# CTG | assembly=reference | chr=1 | chr-pos=1327197 | NT_004350.18 | ctg-start=815966 | ctg-end=815966 | loctype=2 | orient=+
			# CTG | assembly=reference | chr=1 | chr-pos=? | NT_113871.1 | ctg-start=111608 | ctg-end=111608 | loctype=2 | orient=+


		}
		elsif ($class eq 'LOC') {

			# I'm going to skip the LOC info for now - we do a better job of this with our own code.

			# LOC | KCNAB2 | locus_id=8514 | fxn-class=near-gene-3 | mrna_acc=NM_003636.2
			# LOC | KCNAB2 | locus_id=8514 | fxn-class=near-gene-3 | mrna_acc=NM_172130.1
			# LOC | TMED5 | locus_id=50999 | fxn-class=utr-3 | mrna_acc=NM_016040.3
			# LOC | TMEM51 | locus_id=55092 | fxn-class=utr-3 | mrna_acc=NM_018022.2
			# LOC | TMEM51 | locus_id=55092 | fxn-class=utr-3 | mrna_acc=NM_001136216.1
			# LOC | TMEM51 | locus_id=55092 | fxn-class=utr-3 | mrna_acc=NM_001136217.1
			# LOC | TMEM51 | locus_id=55092 | fxn-class=utr-3 | mrna_acc=NM_001136218.1
			# LOC | ATP2B4 | locus_id=493 | fxn-class=utr-3 | mrna_acc=NM_001001396.1

		}
		elsif ($class eq 'VAL') {

			# I'm going to skip the VAL info for now.

			# VAL | validated=YES | min_prob=? | max_prob=? | notwithdrawn | byFrequency
			# VAL | validated=NO | min_prob=? | max_prob=? | notwithdrawn
			# VAL | validated=YES | min_prob=? | max_prob=? | notwithdrawn | byHapMap
			# VAL | validated=YES | min_prob=? | max_prob=? | notwithdrawn | byCluster | by2Hit2Allele
			# VAL | validated=YES | min_prob=95 | max_prob=99 | notwithdrawn | byCluster | byFrequency | byOtherPop | by2Hit2Allele | byHapMap
		}
		elsif ($class eq 'SNP') {

			my ($seq_text, $het, $se_het) = @data;
			$seq_text =~ s/alleles=// or push @errors, ['SNP_invalid_seqs_tag', $line];
			$seq_text =~ s/\'//g;
			$het         =~ s/het=//         or push @errors, ['SNP_invalid_het_tag', $line];
			$se_het      =~ s/se\(het\)=//   or push @errors, ['SNP_invalid_se(het)_tag', $line];

			$record{het}     = $het;
			$record{se_het}  = $se_het;

			# SNP | alleles='-/T' | het=0.18 | se(het)=0.24
			# SNP | alleles='-/T' | het=0.18 | se(het)=0.24
			# SNP | alleles='C/T' | het=? | se(het)=?
			# SNP | alleles='C/T' | het=0.2 | se(het)=0.2443
			# SNP | alleles='A/G' | het=0.39 | se(het)=0.2098
			# SNP | alleles='A/G' | het=0.3 | se(het)=0.2458
			# SNP | alleles='A/G' | het=0.01 | se(het)=0.0839

			# alleles='(LARGEDELETION)/-' Here the sequence deleted is sequence found on the reference
			# alleles='(LARGEINSERTION)-' Here the sequence inserted is not sequence found on the reference
			# alleles='(TG)17/18/19/21/22/23/24/G/T' TG is repeated 17, 18, 19 etc. times as well as a G and a T as variants.
			# alleles='(HETEROZYGOUS)/C/T' Just remove the HETEROZYGOUS
			# alleles='(147 BP INSERTION)/-'
			# alleles='(1110 BP DELETION)/-'

			my @seqs  = split /\//, $seq_text;
			if ($seq_text =~ /[^ATGC-]/) {
				if ($seq_text =~ /\(LARGEDELETION\)/) {
					unshift @seqs,  '-';
					$record{type} = 'deletion';
				}
				if ($seq_text =~ /\(LARGEINSERTION\)/) {
					unshift @seqs, '~';
					$record{type} = 'insertion';
				}



				# FIX THIS !!!!!!!!!!!!!!!!

				if ($seq_text =~ /\(\d+ BP DELETION\)/) {
					unshift @seqs,  '-';
					$record{type} = 'deletion';
				}
				if (/\(([ATGC]+)\)(\d+)/) {
					my $repeat_seq = $1;
					my $count = $2;
					unshift @seqs, $count;
					map {$_ = $repeat_seq x $count if /^\d+$/} @seqs;
				}
				@seqs = grep {/^[ATGC-]+$/} @seqs;
			}

			$record{seqs} = \@seqs;
		}
		elsif ($class eq 'SEQ') {

			# I'm going to skip the SEQ info for now.

			# SEQ | 47271455 | source-db=remap | seq-pos=3818 | orient=-
			# SEQ | 8922276 | source-db=remap | seq-pos=1762 | orient=+
			# SEQ | 209977050 | source-db=remap | seq-pos=1763 | orient=+
			# SEQ | 209977051 | source-db=remap | seq-pos=1867 | orient=+
			# SEQ | 209977053 | source-db=remap | seq-pos=1794 | orient=+
			# SEQ | 209977055 | source-db=remap | seq-pos=1836 | orient=+

		}
		elsif ($class =~ /^rs\d+$/) {

			# $type needs to have the option of being set somewhere else.
			my ($id, $organism, $tax_id, $this_type, $genotype, $submit_link, $update);
			($id, $organism, $tax_id, $this_type, $genotype, $submit_link, $update) = ($class, @data);

			# Here we allow some other line to set the type and preserve that.
			$record{type} ||= $this_type;

			my %type_map = ('in-del'                       => 'sequence_alteration',
					'microsatellite'               => 'microsatellite',
					'mixed'                        => 'sequence_alteration',
					'multinucleotide-polymorphism' => 'MNP',
					'named-locus'                  => 'sequence_alteration',
					'snp'                          => 'SNV',
					'no-variation'                 => 'no_variation', # Add no_variation to SO
					'heterozygous'                 => 'sequence_alteration',
				       );

			# If we can't map it (maybe because some other line set it) we leave it alone.
			$record{type} = $type_map{$record{type}} || $record{type};

			$genotype    =~ s/genotype=//      or push @errors, ['rs_invalid_genotype_tag', $line];
			$submit_link =~ s/submitterlink=// or push @errors, ['rs_invalid_submitterlink_tag', $line];
			$update      =~ s/updated //       or push @errors, ['rs_invalid_update_tag', $line];

			$record{id}          = $id;
			$record{organism}    = $organism;
			$record{tax_id}      = $tax_id;
			$record{genotype}    = $genotype;
			$record{submit_link} = $submit_link;
			$record{update}      = $update;

			# rs242 | human | 9606 | in-del | genotype=NO | submitterlink=YES | updated 2008-03-28 13:52

		}
		else {
			push @errors, ['Uknown_line_type', $line];
		}
	}

	# If we didn't pick up any contigs the user probably failed to
	# specify one or misspelled it.
	unless (ref $record{ctgs} && scalar @{$record{ctgs}}) {
		$self->warn(message => ("Error\tNo_CTGs_tags_parsed\t",
					$record{id}));
		}

	return \%record, \@errors;
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

<GAL::Parser::template> requires no configuration files or environment variables.

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
