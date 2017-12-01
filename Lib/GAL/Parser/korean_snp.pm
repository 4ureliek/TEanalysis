package GAL::Parser::korean_snp;

use strict;
use vars qw($VERSION);


$VERSION = 0.2.0;
use base qw(GAL::Parser);
use GAL::Reader::DelimitedLine;

=head1 NAME

GAL::Parser::korean_snp - Parse SNV files from the first Korean Genome

=head1 VERSION

This document describes GAL::Parser::korean_snp version 0.2.0

=head1 SYNOPSIS

    my $parser = GAL::Parser::korean_snp->new(file => 'korean_snp.gff');

    while (my $feature_hash = $parser->next_feature_hash) {
	print $parser->to_gff3($feature_hash) . "\n";
    }

=head1 DESCRIPTION

L<GAL::Parser::korean_snp> provides a parser for SNV files from the
first Korean genome published by Ahn, et al. 2009
(http://www.ncbi.nlm.nih.gov/pubmed/19470904).

=head1 Constructor

New L<GAL::Parser::korean_snp> objects are created by the class method new.
Arguments should be passed to the constructor as a list (or reference)
of key value pairs.  All attributes of the L<GAL::Parser::korean_snp> object
can be set in the call to new. An simple example of object creation
would look like this:

    my $parser = GAL::Parser::korean_snp->new(file => 'korean_snp.gff');

The constructor recognizes the following parameters which will set the
appropriate attributes:

=over 4

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
     Usage   : GAL::Parser::korean_snp->new();
     Function: Creates a GAL::Parser::korean_snp object;
     Returns : A GAL::Parser::korean_snp object
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

	my $id         = $record->{id};
	my $seqid      = $record->{chromosome};
	my $source     = 'KOBIC';
	my $type       = 'SNV';
	my $start      = $record->{location};
	my $end        = $record->{location};
	my $score      = '.';
	my $strand     = '+';
	my $phase      = '.';

	# Het with ref: get var_alleles and remove ref.  ref_reads and var_reads OK
	# chr10   56397   C       CT      rs12262442      28      C/T     17      11
	# chr10   61776   T       CT      rs61838967      15      T/C     7       8
	# chr10   65803   T       CT      KOREFSNP1       27      T/C     19      8
	# chr10   68106   C       AC      KOREFSNP2       43      C/A     22      21
	# chr10   84136   C       CT      rs4607995       24      C/T     10      13
	# chr10   84238   A       AT      rs10904041      22      A/T     5       16

	# Het but not ref: get var_alleles.  assign var_reads to correct var and calculate alter_var_reads from total - var_reads
	# chr10   12625631        A       GT      rs2815636       42      A/G     0       21
	# chr10   13864035        A       CT      rs5025431       27      A/T     0       15
	# chr10   14292681        G       AC      rs11528656      29      G/A     0       18
	# chr10   14771944        C       AG      rs3107794       29      C/G     0       15
	# chr10   15075637        A       CG      rs9731518       29      A/G     4       16

	# Homozygous get var_alleles and use only one.  ref_reads and var_reads OK
	# chr10   168434  T       GG      rs7089889       20      T/G     0       20
	# chr10   173151  T       CC      rs7476951       19      T/C     0       19
	# chr10   175171  G       TT      rs7898275       25      G/T     0       25
	# chr10   175358  C       TT      rs7910845       26      C/T     0       26

	# $self->fields([qw(chromosome location ref_seq var_seqs id total_reads read_seqs ref_reads var_reads)]);

	my $reference_seq = $record->{ref_seq};
	my %variant_seqs  = map {$_, 1} split //, $record->{var_seqs};
	my @variant_seqs = keys %variant_seqs;

	my $total_reads = $record->{total_reads};

	# chr10   56397   C       CT      rs12262442      28      C/T     17      11
	my @read_seqs = split m|/|, $record->{read_seqs};
	my %read_counts = ($read_seqs[0] => $record->{ref_reads},
			   $read_seqs[1] => $record->{var_reads},
			   );

	my @variant_reads = map {$read_counts{$_} || $total_reads - $record->{var_reads}} @variant_seqs;

	my $zygosity = scalar @variant_seqs > 1 ? 'heterozygous' : 'homozygous';

	my $attributes = {Reference_seq => [$reference_seq],
			  Variant_seq   => \@variant_seqs,
			  Variant_reads => \@variant_reads,
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
 Returns : A GAL::Reader::DelimitedLine singleton.
 Args    : None

=cut

sub reader {
  my $self = shift;

  if (! $self->{reader}) {
    my @field_names = qw(chromosome location ref_seq var_seqs id
			  total_reads read_seqs ref_reads var_reads);
    my $reader = GAL::Reader::DelimitedLine->new(field_names => \@field_names);
    $self->{reader} = $reader;
  }
  return $self->{reader};
}

#-----------------------------------------------------------------------------

=head1 DIAGNOSTICS

L<GAL::Parser::korean_snp> does not throw any warnings or errors.

=head1 CONFIGURATION AND ENVIRONMENT

L<GAL::Parser::korean_snp> requires no configuration files or environment variables.

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
