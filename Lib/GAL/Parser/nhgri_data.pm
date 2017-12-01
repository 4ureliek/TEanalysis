package GAL::Parser::nhgri_data;

use strict;
use vars qw($VERSION);
$VERSION = 0.2.0;

use base qw(GAL::Parser);
use GAL::Reader::DelimitedLine;

=head1 NAME

GAL::Parser::nhgri_data - Parse NHGRI_DATA files

=head1 VERSION

This document describes GAL::Parser::nhgri_data version 0.2.0

=head1 SYNOPSIS

    my $parser = GAL::Parser::nhgri_data->new(file => 'nhgri_data.txt');

    while (my $feature_hash = $parser->next_feature_hash) {
	print $parser->to_gff3($feature_hash) . "\n";
    }

=head1 DESCRIPTION

L<GAL::Parser::nhgri_data> provides a parser for NHGRI_DATA data.

=head1 Constructor

New L<GAL::Parser::nhgri_data> objects are created by the class method
new.  Arguments should be passed to the constructor as a list (or
reference) of key value pairs.  All attributes of the
L<GAL::Parser::nhgri_data> object can be set in the call to new. An
simple example of object creation would look like this:

    my $parser = GAL::Parser::nhgri_data->new(file => 'nhgri_data.txt');

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
     Usage   : GAL::Parser::nhgri_data->new();
     Function: Creates a GAL::Parser::nhgri_data object;
     Returns : A GAL::Parser::nhgri_data object
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
	my @valid_attributes = qw(file fh); # Set valid class attributes here
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

	# I modified the original file a bit to have only one individual per file.
	# That makes the last column the genotype for that individual.

	# Chr  LeftFlank RightFlankref_allele var_allele type refseq ref_aa var_aa aa_pos RS# genotype
	# chrX 139125960 139125962 C T NC 	NA 	    NA NA 0 "rs396467(C,T)"   TT
	# chrX 139126277 139126279 T C NC 	NA 	    NA NA 0 "rs411420(C,T)"   CC
	# chrX 139623329 139623331 G A Intronic RP1-177G6.2 NA NA 0 "rs41304542(A,G)" AG
	# chrX 139624427 139624429 T C 3'UTR 	RP1-177G6.2 NA NA 0 "rs11095830(C,T)" CC
	# chrX 139690983 139690985 C T 3'UTR 	AK054921    NA NA 0 "rs5907055(C,T)"  TT

	my $seqid      = $record->{Chr};
	my $source     = 'NHGRI';
	my $type       = 'SNV';
	my $start      = $record->{LeftFlank} + 1;
	my $end        = $start;
	my $feature_id = join ':', ($seqid, $source, $type, $start, $end);
	my $score      = '.';
	my $strand     = '+';
	my $phase      = '.';

	my $reference_seq = $record->{ref_allele};
	
	# Ignore no-calls for now.
	return undef if $record->{genotype} eq 'NA';
	# Ignore indels for now.
	return undef if length $record->{genotype} > 2;

	my @variant_seqs = split '', $record->{genotype};
	shift @variant_seqs if $variant_seqs[0] eq $variant_seqs[1];

	return undef if (scalar @variant_seqs == 1 && $variant_seqs[0] eq $reference_seq);

	my $zygosity = scalar @variant_seqs > 1 ? 'heterozygous' : 'homozygous';

	my $alias = $record->{RS} if $record->{RS};

	my $attributes = {Reference_seq => [$reference_seq],
			  Variant_seq   => \@variant_seqs,
			  Zygosity      => [$zygosity],
			  ID            => [$feature_id],
			 };

	push @{$attributes->{Alias}}, $alias if $alias;

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
 Returns : A L<GAL::Reader::DelimitedLine> singleton.
 Args    : None

=cut

sub reader {
  my $self = shift;

  if (! $self->{reader}) {
      
    # Change the field names in the next line to define how the hash keys
    # are named for each column in you input file.
    # If you don't pass field_names as an argument to the reader then
    # the reader will return an array rather than a hash.
    my @field_names = qw(Chr  LeftFlank RightFlank ref_allele var_allele type refseq ref_aa var_aa aa_pos RS genotype);
    my $reader = GAL::Reader::DelimitedLine->new(field_names     => \@field_names,
						 field_separator => "\t");
    $self->{reader} = $reader;
  }
  return $self->{reader};
}

#-----------------------------------------------------------------------------

=head1 DIAGNOSTICS

L<GAL::Parser::nhgri_data> does not throw any warnings or errors.

=head1 CONFIGURATION AND ENVIRONMENT

L<GAL::Parser::nhgri_data> requires no configuration files or
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
