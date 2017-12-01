package GAL::Parser::samtools_pileup;

use strict;
use vars qw($VERSION);


$VERSION = 0.2.0;
use base qw(GAL::Parser);
use GAL::Reader::DelimitedLine;

=head1 NAME

GAL::Parser::samtools_pileup - Parse SAMTOOLS pileup files

=head1 VERSION

This document describes GAL::Parser::samtools_pileup version 0.2.0

=head1 SYNOPSIS

    use GAL::Parser::samtools_pileup;
    my $parser = GAL::Parser::samtools_pileup->new(file => 'samtools_pileup.txt');
    while (my $feature_hash = $parser->next_feature_hash) {
	print $parser->to_gff3($feature_hash) . "\n";
    }

=head1 DESCRIPTION

L<GAL::Parser::samtools_pileup> provides a parser for SAMTOOLS pileup
files (http://samtools.sourceforge.net/pileup.shtml).

=head1 Constructor

New L<GAL::Parser::samtools_pileup> objects are created by the class
method new.  Arguments should be passed to the constructor as a list
(or reference) of key value pairs.  All attributes of the
L<GAL::Parser::samtools_pileup> object can be set in the call to
new. An simple example of object creation would look like this:

    my $parser = GAL::Parser::samtools_pileup->new(file => 'samtools_pileup.txt');

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
     Usage   : GAL::Parser::samtools_pileup->new();
     Function: Creates a GAL::Parser::samtools_pileup object;
     Returns : A GAL::Parser::samtools_pileup object
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


	# This parser was written to convert variant output created
	# with SAMtools like this:

	# From a sorted BAM alignment, raw SNP and indel calls are acquired by:
	#
	# 1. samtools pileup -vcf ref.fa aln.bam > raw.pileup
	# 
	# samtools pileup -vcf ref.fa aln.bam > raw.pileup The resultant output
	# should be further filtered by:
	# 
	# 1. samtools.pl varFilter raw.pileup | awk '$6>=20' > final.pileup  
	# 
	# samtools.pl varFilter raw.pileup | awk '$6>=20' > final.pileup to rule
	# out error-prone variant calls caused by factors not considered in the
	# statistical model.

	# With output that looks like this

	# chr1  10109  a  W  79  79   60  31  .$...t,.T.t,,tt.,tt..,,ttT,,,t,t  %%%%BB%%%ACC@=BB..%9AB<9>A'B9%7
	# chr1	10177  a  C  18  33   60  3   c,C     		      		53@       
	# chr1	13116  T  G  11  39   60  6   g.GgG,  		      		B39A>1
	# chr1	13118  A  G  8   39   60  6   g.GgG,  		      		B?;AA(
	# chr1	14464  A  W  22  41   60  5   t,Tt.   		      		AA@9A
	# chr1	14653  C  T  28  36   60  4   .TtT    		      		%>%>            
	# chr1	14907  A  R  85  147  60  18  ,,,g,.gGGGGggG.g^~G^~g  		ACB5BB7>@BB3>@<:A%
	# chr1	14930  A  R  27  99   60  19  g,.gGGGGggG.gGg.ggg?    		BB<AB@<5@B6BA=A4=1
	# chr1	15208  G  A  6   27   60  4   ,aa,    		      		>AB%
	# chr1	15211  T  G  39  39   60  4   gggg    		      		>@B@            

	# Het indels
	# chrX  305488   *  +G/+G  7  35  35  3  +G    *  2  1  0  0  0
	# chrX  719509   *  -C/-C  4  77  93 73  -C    *  3  0  0  0  0
	# chrX  787510   *  +G/+G 50 119  40  4  +G    *  4  0  0  0  0
	        
	# Homo indels
	# chrX  108090   *  +CT/*   5   98  40  3  +CT   *   3  0  0  0  0
	# chrX  136584   *  */-c    3   23 218  6  *    -c   5  1  0  0  0
	# chrX  273736   *  */-cca  29  29  26  3  *  -cca   2  1  0  0  0

	# With columns described as this:

	# http://samtools.sourceforge.net/samtools.shtml  See both main pileup
	# docs and the details for the -c option.
	# 1.  chromosome name
	# 2.  coordinate
	# 3.  reference base
	# 4.  consensus base
	# 5.  Phred-scaled consensus quality
	# 6.  SNP quality (i.e. the Phred-scaled probability of the consensus
	#     being identical to the reference)
	# 7.  root mean square (RMS) mapping quality of the reads covering
	#     the site
	# 8.  read bases
	# 9.  read qualities
	# 10. alignment mapping qualities
	
	# Mapped to record keys:

	# 1.  seqid
	# 2.  start
	# 3.  reference_seq
	# 4.  variant_seq
	# 5.  consensus_phred_qual
	# 6.  snv_phred_qual
	# 7.  rms
	# 8.  variant_reads
	# 9.  read_qual
	# 10. aln_map_qual

	my $seqid      = $record->{seqid};
	my $source     = 'SAMtools';
	my $type       = 'SNV';
	my $start      = $record->{start};
	my $id         = join ':', ($seqid, $source, $type, $start);
	my $end;
	my $score      = $record->{snv_phred_qual};
	my $strand     = '+';
	my $phase      = '.';

	# Create the attributes hash
	# See http://www.sequenceontology.org/gvf.html

	my $reference_seq = uc $record->{reference_seq};
	my @variant_seqs;
	if ($reference_seq eq '*') {
	    my %var_indel_hash;
	    @var_indel_hash{(split m|/|, $record->{variant_seq})} = ();
	    @variant_seqs = keys %var_indel_hash;
	    map {uc} @variant_seqs;
	    if (grep {/^\+/} @variant_seqs) {
		map {s/\+//} @variant_seqs;
		$type = 'insertion';
		$reference_seq = '-';
		map {$_ = $_ eq '*' ? $reference_seq : $_} @variant_seqs;
	    }
	    else {
		map {s/\-//} @variant_seqs;
		$type = 'deletion';
		$start = $strand eq '+' ? $start + 1 : $start - 1;
		my $length = 0;
		map {$length = length($_) > $length ? length($_) : $length} @variant_seqs;
		$end = $start + $length - 1;
		$reference_seq = uc $self->fasta->seq($seqid, $start, $end);
		map {$_ = $_ eq '*' ? $reference_seq : '-'} @variant_seqs;
	    }
	}
	else {
	    @variant_seqs = $self->expand_iupac_nt_codes($record->{variant_seq});
	}
	$end ||= $start;

	my $total_reads = $record->{variant_reads};

	my $zygosity = scalar @variant_seqs == 1 ? 'homozygous' : 'heterozygous';

	# Create the attribute hash reference.  Note that all values
	# are array references - even those that could only ever have
	# one value.  This is for consistency in the interface to
	# Features.pm and it's subclasses.  Suggested keys include
	# (from the GFF3 spec), but are not limited to: ID, Name,
	# Alias, Parent, Target, Gap, Derives_from, Note, Dbxref and
	# Ontology_term. Note that attribute names are case
	# sensitive. "Parent" is not the same as "parent". All
	# attributes that begin with an uppercase letter are reserved
	# for later use. Attributes that begin with a lowercase letter
	# can be used freely by applications.

	# For sequence_alteration features the suggested keys include:
	# ID, Reference_seq, Variant_seq, Variant_reads Total_reads,
	# Zygosity, Variant_effect, Copy_number

	my $attributes = {Reference_seq => [$reference_seq],
			  Variant_seq   => \@variant_seqs,
			  Total_reads   => [$total_reads],
			  Zygosity      => [$zygosity],
			  ID            => [$id],
			 };

	my $feature_data = {feature_id => $id,
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
 Returns : A L<GAL::Reader::DelimitedLine> singleton.
 Args    : None

=cut

sub reader {
  my $self = shift;

  if (! $self->{reader}) {
    my @field_names = qw(seqid start reference_seq variant_seq
			  consensus_phred_qual snv_phred_qual rms
			  variant_reads read_qual aln_map_qual);
    my $reader = GAL::Reader::DelimitedLine->new(field_names => \@field_names);
    $self->{reader} = $reader;
  }
  return $self->{reader};
}

#-----------------------------------------------------------------------------

=head1 DIAGNOSTICS

L<GAL::Parser::samtools_pileup> does not throw any warnings or errors.

=head1 CONFIGURATION AND ENVIRONMENT

L<GAL::Parser::samtools_pileup> requires no configuration files or
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
