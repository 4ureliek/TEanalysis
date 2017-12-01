package GAL::Parser::chinese_snp;

use strict;
use vars qw($VERSION);

$VERSION = 0.2.0;
use base qw(GAL::Parser);
use GAL::Reader::DelimitedLine;

=head1 NAME

GAL::Parser::chinese_snp - Parse SNV files from the first Chinese
genome

=head1 VERSION

This document describes GAL::Parser::chinese_snp version 0.2.0

=head1 SYNOPSIS

    my $parser = GAL::Parser::chinese_snp->new(file =>
             'chinese_snp.gff');

    while (my $feature_hash = $parser->next_feature_hash) {
	print $parser->to_gff3($feature_hash) . "\n";
    }

=head1 DESCRIPTION

L<GAL::Parser::chinese_snp> provides a parser for SNV files from the
first Chinese genome published by Wang, et al. 2008
(http://www.ncbi.nlm.nih.gov/pubmed/18987735).

=head1 Constructor

New L<GAL::Parser::chinese_snp> objects are created by the class
method new.  Arguments should be passed to the constructor as a list
(or reference) of key value pairs.  All attributes of the
L<GAL::Parser::chinese_snp> object can be set in the call to new. An
simple example of object creation would look like this:

    my $parser = GAL::Parser::chinese_snp->new(file =>
             'chinese_snp.gff');

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
     Usage   : GAL::Parser::chinese_snp->new();
     Function: Creates a GAL::Parser::chinese_snp object;
     Returns : A GAL::Parser::chinese_snp object
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

	my $original_atts = $self->parse_attributes($record->{attributes});

	my $seqid      = $record->{seqid};
	my $source     = $record->{source};
	my $type       = $record->{type};
	my $start      = $record->{start};
	my $end        = $record->{end};
	my $score      = $record->{score};
	my $strand     = $record->{strand};
	my $phase      = '.';
	my $alias_xref = $original_atts->{ID}[0];

	$type = $type eq 'SNP' ? 'SNV' : $type;

	my $id = join ':', ($seqid, $source, $type, $start);

	# chr1 SoapSnp SNP SNP 4793 4793 25 + . ID=YHSNP0128643; status=novel; ref=A; allele=A/G; support1=48; support2=26;
	# chr1SoapSNPSNP6434643448+.ID=YHSNP0128644; status=novel; ref=G; allele=A/G; support1=10; support2=11;
	# chr1SoapSNPSNP938969389651+.ID=rs4287120; status=dbSNP; ref=T; allele=C/T; support1=5; support2=4; location=MSTB1:LTR/MaLR;
	# chr1SoapSNPSNP22570722570743+.ID=rs6603780; status=dbSNP; ref=C; allele=C/G; support1=23; support2=12;
	# chr1SoapSNPSNP22583922583931+.ID=rs6422503; status=dbSNP; ref=C; allele=A/C; support1=13; support2=5; location=L1P2:LINE/L1;
	# chr1SoapSNPSNP52684952684976+.ID=YHSNP0128645; status=novel; ref=G; allele=G/T; support1=14; support2=12; location=L1MD3:LINE/L1;
	# chr1SoapSNPSNP55473155473130+.ID=rs1832728; status=dbSNP; ref=T; allele=C/T; support1=37; support2=12; location=Mitochondrial:Mt-tRNA;
	# chr1SoapSNPSNP55535355535328+.ID=rs7349153; status=dbSNP; ref=T; allele=C/T; support1=37; support2=9;
	# chr1SoapSNPSNP55537155537122+.ID=rs9283150; status=dbSNP; ref=G; allele=A/G; support1=46; support2=27;
	# chr1SoapSNPSNP55677955677945+.ID=rs3949348; status=dbSNP; ref=A; allele=A/G; support1=37; support2=13;
	# chr1    SoapSNP SNP     774913  774913  74      +       .       ID=rs2905062; status=dbSNP; ref=G; allele=A/A; support1=26; location=MSTD:LTR/MaLR;
	# chr1    SoapSNP SNP     775852  775852  93      +       .       ID=rs2980300; status=dbSNP; ref=T; allele=C/C; support1=29;
	# chr1    SoapSNP SNP     777262  777262  43      +       .       ID=rs2905055; status=dbSNP; ref=G; allele=T/T; support1=12;

	my $reference_seq = $original_atts->{ref}[0];
	my @variant_seqs  = split m|/|, $original_atts->{allele}[0];

	shift @variant_seqs if $variant_seqs[0] eq $variant_seqs[1];

	my $support1 = (ref $original_atts->{support1} eq 'ARRAY' ?
			$original_atts->{support1}[0]             :
			0
		       );
	my $support2 = (ref $original_atts->{support2} eq 'ARRAY' ?
			$original_atts->{support2}[0]             :
			0
		       );

	my @variant_reads = ($support1, $support2);

	my $total_reads = $support1 + $support2;

	my $zygosity = scalar @variant_seqs > 1 ? 'heterozygous' : 'homozygous';

	my $attributes = {ID            => [$id],
			  Reference_seq => [$reference_seq],
			  Variant_seq   => \@variant_seqs,
			  Variant_reads => \@variant_reads,
			  Total_reads   => [$total_reads],
			  Zygosity      => [$zygosity],
			 };

	if ($alias_xref =~ /^rs\d+$/) {
	    push @{$attributes->{Dbxref}}, 'dbSNP:' . $alias_xref;
	}
	else {
	    push @{$attributes->{Alias}}, $alias_xref;
	}

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
    my @field_names = qw(seqid source type start end score strand phase attributes);
    my $reader = GAL::Reader::DelimitedLine->new(field_names => \@field_names);
    $self->{reader} = $reader;
  }
  return $self->{reader};
}

#-----------------------------------------------------------------------------

=head1 DIAGNOSTICS

L<GAL::Parser::chinese_snp> does not throw any warnings or errors.

=head1 CONFIGURATION AND ENVIRONMENT

L<GAL::Parser::chinese_snp> requires no configuration files or
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
