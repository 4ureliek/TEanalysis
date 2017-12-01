package GAL::Parser::trait_o_matic;

use strict;
use vars qw($VERSION);

$VERSION = 0.2.0;
use base qw(GAL::Parser);
use GAL::Reader::DelimitedLine;

=head1 NAME

GAL::Parser::trait_o_matic - Parse TRAIT-O-MATIC SNP files

=head1 VERSION

This document describes GAL::Parser::trait_o_matic version 0.2.0

=head1 SYNOPSIS

    my $parser = GAL::Parser::trait_o_matic->new(file =>
	     'trait_o_matic.gff');

    while (my $feature_hash = $parser->next_feature_hash) {
	print $parser->to_gff3($feature_hash) . "\n";
    }

=head1 DESCRIPTION

L<GAL::Parser::trait_o_matic> provides a parser for TRAIT-O-MATIC SNP
data (http://snp.med.harvard.edu/).

=head1 Constructor

New L<GAL::Parser::trait_o_matic> objects are created by the class
method new.  Arguments should be passed to the constructor as a list
(or reference) of key value pairs.  All attributes of the
L<GAL::Parser::trait_o_matic> object can be set in the call to new. An
simple example of object creation would look like this:

    my $parser = GAL::Parser::trait_o_matic->new(file => 'trait_o_matic.gff');

The constructor recognizes the following parameters which will set the
appropriate attributes:

=over

=item * C<< file => feature_file.txt >>

This optional parameter provides the filename for the file containing
the data to be parsed. While this parameter is optional either it, or
the following fh parameter must be set.

=item * C<< fh => feature_file.txt >>

This optional parameter provides a filehandle to read data from. While
this parameter is optional either it, or the following fh parameter
must be set.

=back

=cut

#-----------------------------------------------------------------------------

=head2 new

     Title   : new
     Usage   : GAL::Parser::trait_o_matic->new();
     Function: Creates a GAL::Parser::trait_o_matic object;
     Returns : A GAL::Parser::trait_o_matic object
     Args    : See the attributes described above.

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
	my @valid_attributes = qw(); # Set valid class attributes here
	$self->set_attributes($args, @valid_attributes);
	######################################################################
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

	# ==> /home/bmoore/Projects/10-05-17_Trait-o-Matic_Genomes/PGP01_NA20431_Church.gff <==
	# chr1	SIQ2	SNP	2266060	2266060	.	+	.	alleles A/C;ref_allele T
	# ==> /home/bmoore/Projects/10-05-17_Trait-o-Matic_Genomes/PGP02_NA21070_Halamka.gff <==
	# chr1	maq	SNP	3751339	3751339	.	+	.	alleles T;ref_allele C;read_depth 255
	# ==> /home/bmoore/Projects/10-05-17_Trait-o-Matic_Genomes/PGP03_NA21660_Dyson.gff <==
	# chr1	maq	SNP	1708759	1708759	.	+	.	alleles A;ref_allele T;read_depth 8
	# ==> /home/bmoore/Projects/10-05-17_Trait-o-Matic_Genomes/PGP04_NA21677_Angrist.gff <==
	# chr21	.	SNP	9719804	9719804	255	+	.	alleles C;ref_allele T;read_depth 135
	# ==> /home/bmoore/Projects/10-05-17_Trait-o-Matic_Genomes/PGP05_21687_Maxey.gff <==
	# chr1	maq	SNP	3751339	3751339	.	+	.	alleles T;ref_allele C;read_depth 79
	# ==> /home/bmoore/Projects/10-05-17_Trait-o-Matic_Genomes/PGP06_NA21730_Pinker.gff <==
	# chr1	maq	SNP	2431201	2431201	.	+	.	alleles C/G;ref_allele G;read_depth 22
	# ==> /home/bmoore/Projects/10-05-17_Trait-o-Matic_Genomes/PGP07_NA21731_Batchelder.gff <==
	# chr1	maq	SNP	3751339	3751339	.	+	.	alleles C/T;ref_allele C;read_depth 78
	# ==> /home/bmoore/Projects/10-05-17_Trait-o-Matic_Genomes/PGP08_NA21781_Lapidus.gff <==
	# chr1	maq	SNP	1708759	1708759	.	+	.	alleles A;ref_allele T;read_depth 23
	# ==> /home/bmoore/Projects/10-05-17_Trait-o-Matic_Genomes/PGP09_NA21833_Gill.gff <==
	# chr10	MAQ	SNP	10000011	10000011	.	+	.	alleles C
	# ==> /home/bmoore/Projects/10-05-17_Trait-o-Matic_Genomes/PGP10_NA21846_Sherley.gff <==
	# chr1	maq	SNP	1946645	1946645	.	+	.	alleles A/C;ref_allele A;read_depth 34
	# ==> /home/bmoore/Projects/10-05-17_Trait-o-Matic_Genomes/PGP11_NA_GatesSr.gff <==
	# chr10	-	-	53571	53571	.	+	.	alleles C
	# ==> /home/bmoore/Projects/10-05-17_Trait-o-Matic_Genomes/PGP12_NA_GatesJr.gff <==
	# chr10	-	-	60592	60592	.	+	.	alleles T
	# ==> /home/bmoore/Projects/10-05-17_Trait-o-Matic_Genomes/PGP_NA07022_Anonymous.gff <==
	# chr10	SIQ2	SNP	10000005	10000005	.	+	.	alleles G;ref_allele T
	# ==> /home/bmoore/Projects/10-05-17_Trait-o-Matic_Genomes/PGP_NA12156_Anonymous.gff <==
	# chr1	MAQ	SNP	861078	861078	.	+	.	alleles C;ref allele C
	# ==> /home/bmoore/Projects/10-05-17_Trait-o-Matic_Genomes/PGP_NA12878_Anonymous.gff <==
	# chr1	MAQ	SNP	861078	861078	.	+	.	alleles C;ref allele C
	# ==> /home/bmoore/Projects/10-05-17_Trait-o-Matic_Genomes/PGP_NA18507_Anonymous.gff <==
	# chr1	NA18507	SNP	18455	18455	.	+	.	alleles T/G;id 000000001
	# ==> /home/bmoore/Projects/10-05-17_Trait-o-Matic_Genomes/PGP_NA18517_Anonymous.gff <==
	# chr1	MAQ	SNP	59374	59374	.	+	.	alleles G;ref allele A
	# ==> /home/bmoore/Projects/10-05-17_Trait-o-Matic_Genomes/PGP_NA18555_Anonymous.gff <==
	# chr1	MAQ	SNP	59374	59374	.	+	.	alleles G;ref allele A
	# ==> /home/bmoore/Projects/10-05-17_Trait-o-Matic_Genomes/PGP_NA18956_Anonymous.gff <==
	# chr1	MAQ	SNP	59374	59374	.	+	.	alleles G;ref allele A
	# ==> /home/bmoore/Projects/10-05-17_Trait-o-Matic_Genomes/PGP_NA19129_Anonymous.gff <==
	# chr1	MAQ	SNP	59374	59374	.	+	.	alleles G;ref allele A
	# ==> /home/bmoore/Projects/10-05-17_Trait-o-Matic_Genomes/PGP_NA19240_Anonymous.gff <==
	# chr1	asm	SNP	28095	28095	.	+	.	alleles G;ref_allele A;db_xref dbsnp:rs806727;call_score 61;ID 819
	# ==> /home/bmoore/Projects/10-05-17_Trait-o-Matic_Genomes/PGP_NA_-Ai.gff <==
	# chr1	TRGFF	SNP	798785	798785	.	+	.	 alleles A;ref_allele G
	# ==> /home/bmoore/Projects/10-05-17_Trait-o-Matic_Genomes/PGP_NA_Anonymous_Chinese_Male.gff <==
	# chr1	SoapSNP	SNP	4793	4793	25.0	+	.	ID YHSNP0128643;alleles A/G;counts 48/26;ref_allele A;status novel
	# ==> /home/bmoore/Projects/10-05-17_Trait-o-Matic_Genomes/PGP_NA_Anonymous_Korean_Male.gff <==
	# chr10	GSNAP	SNP	100000568	100000568	.	+	.	alleles T/A;ref_allele A;read_depth 30
	# ==> /home/bmoore/Projects/10-05-17_Trait-o-Matic_Genomes/PGP_NA_D-kgao.gff <==
	# chr1	TRGFF	SNP	850871	850871	.	+	.	 alleles C;ref_allele G
	# ==> /home/bmoore/Projects/10-05-17_Trait-o-Matic_Genomes/PGP_NA_Gaq-o.gff <==
	# chr1	TRGFF	SNP	516656	516656	.	+	.	 alleles G;ref_allele A
	# ==> /home/bmoore/Projects/10-05-17_Trait-o-Matic_Genomes/PGP_NA_-Gubi.gff <==
	# chr1	TRGFF	SNP	257	257	.	+	.	 alleles A/C;ref_allele A
	# ==> /home/bmoore/Projects/10-05-17_Trait-o-Matic_Genomes/PGP_NA_Kim.gff <==
	# chr10	maq	SNP	56397	56397	.	+	.	alleles C/T;ref_allele C
	# ==> /home/bmoore/Projects/10-05-17_Trait-o-Matic_Genomes/PGP_NA_Quake.gff <==
	# chr1	HELICOS	SNP	237448542	237448542	.	+	.	alleles A/A;ref_allele G
	# ==> /home/bmoore/Projects/10-05-17_Trait-o-Matic_Genomes/PGP_NA_Tutu.gff <==
	# chr1	TRGFF	SNP	6120	6120	.	+	.	 alleles G/C;ref_allele G
	# ==> /home/bmoore/Projects/10-05-17_Trait-o-Matic_Genomes/PGP_NA_Venter.gff <==
	# chr10	CV	snp	10000005	10000005	.	+	.	alleles T/G;db_xref dbsnp:rs4747840;id 1103649859189;method Method1;ref_allele T;rmr 0;tr 1
	# ==> /home/bmoore/Projects/10-05-17_Trait-o-Matic_Genomes/PGP_NA_Watson.gff <==
	# chr1	JW	snp	41921	41921	.	+	.	alleles G/C;SNP BJW-1117373;counts 2/2;ref_allele G

	my $feature_id = join ':', @{$record}{qw(seqid source type start end)};
	my ($seqid, $source, $type, $start, $end, $score, $strand, $phase, $attributes_text) =
	    @{$record}{qw(seqid source type start end score strand phase attributes)};

	# Sources seen in the trait-o-matic files
	# -
	# .
	# asm
	# CV
	# GSNAP
	# HELICOS
	# JW
	# maq
	# MAQ
	# NA18507
	# SIQ2
	# SoapSNP
	# TRGFF

	$source =~ s/^asm$/ASM/;
	$source =~ s/^maq$/MAQ/;
	$source =~ s/^(-|\.|CV|JW|NA18507)$/PGP/;

	# -
	# snp
	# SNP
	$type = 'SNV';

	my %their_attributes;
	my @pairs = split ';', $attributes_text;
	for my $pair (@pairs) {
	    $pair =~ s/^\s*//;
	    $pair =~ s/\s*$//;
	    my ($key, $value_text) = split /\s/, $pair;
	    my @values = split /\//, $value_text;
	    push @{$their_attributes{$key}}, @values;
	}

	# Here are all of the attriute keys they use
	# alleles    = Variant_seq
	# call_score = Don't use
	# counts     = Variant_reads
	# db_xref    = Dbxref
	# id	     = ID
	# ID	     = ID
	# location   = Don't use
	# method     = Don't use
	# read_depth = Total_reads
	# ref	     = Reference_seq
	# ref_allele = Reference_seq
	# rmr	     = Don't use
	# SNP	     = ID
	# status     = Don't use
	# tr	     = Don't use

	my @variant_seqs = @{$their_attributes{alleles}};
	my ($reference_seq, @variant_reads, @dbxrefs, $total_reads);
	if (exists $their_attributes{ref}) {
	    $reference_seq = $their_attributes{ref}[0];
	}
	elsif (exists $their_attributes{ref_allele}) {
	    $reference_seq  = $their_attributes{ref_allele}[0];
	}
	if (exists $their_attributes{counts}) {
	    @variant_reads = $their_attributes{counts};
	}
	if (exists $their_attributes{dbxref}) {
	    @dbxrefs = @{$their_attributes{dbxref}};
	}
	if (exists $their_attributes{id}) {
	    $feature_id = $their_attributes{id}[0];
	}
	elsif (exists $their_attributes{ID}) {
	    $feature_id = $their_attributes{ID}[0];
	}
	if (exists $their_attributes{SNP}) {
	    $feature_id = $their_attributes{SNP}[0];
	}
	if (exists $their_attributes{read_depth}) {
	    $total_reads = $their_attributes{read_depth}[0];
	}

	$reference_seq ||= uc $self->fasta->seq($seqid, $start, $end);
	$reference_seq = $self->revcomp($reference_seq) if $strand eq '-';

	my $genotype = scalar @variant_seqs > 1 ? 'heterozygous' : 'homozygous';

	my %attributes = (ID            => [$feature_id],
			  Variant_seq   => \@variant_seqs,
			  Reference_seq => [$reference_seq],
			  );

	$attributes{Variant_reads} = \@variant_reads if scalar @variant_reads;
	$attributes{Total_reads}   = [$total_reads]  if defined $total_reads;
	$attributes{Dbxref}        = \@dbxrefs       if scalar @dbxrefs;

	my $feature_data = {feature_id => $feature_id,
			    seqid      => $seqid,
			    source     => $source,
			    type       => $type,
			    start      => $start,
			    end        => $end,
			    score      => $score,
			    strand     => $strand,
			    phase      => $phase,
			    attributes => \%attributes,
			   };
	return $feature_data;
}

#-----------------------------------------------------------------------------

=head2 reader

 Title   : reader
 Usage   : $a = $self->reader
 Function: Return the reader object.
 Returns : A L<GAL::Reader::DelimitedLine> singleton.
 Args    : None

=cut

sub reader {
  my $self = shift;

  if (! $self->{reader}) {
    my @field_names = qw(seqid source type start end score strand phase attributes);
    my $reader = GAL::Reader::DelimitedLine->new(field_names => \@field_names);
    $self->{reader} = $reader;
  }
  return $self->{reader};
}

#-----------------------------------------------------------------------------

=head1 DIAGNOSTICS

L<GAL::Parser::trait_o_matic> does not throw any warnings or errors.

=head1 CONFIGURATION AND ENVIRONMENT

L<GAL::Parser::trait_o_matic> requires no configuration files or environment variables.

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

Copyright (c) 2010-2014, Barry Moore <barry.moore@genetics.utah.edu>.  All
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
