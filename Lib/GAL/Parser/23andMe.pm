package GAL::Parser::23andMe;

use strict;
use vars qw($VERSION);

$VERSION = 0.2.0;
use base qw(GAL::Parser);
use GAL::Reader::DelimitedLine;

=head1 NAME

GAL::Parser::23andMe - Parse 23andMe genotype files

=head1 VERSION

This document describes GAL::Parser::23andMe version 0.2.0

=head1 SYNOPSIS

    my $parser = GAL::Parser::23andMe->new(file => 'genotypes.txt');

    while (my $feature_hash = $parser->next_feature_hash) {
	print $parser->to_23andMe($feature_hash) . "\n";
    }

=head1 DESCRIPTION

L<GAL::Parser::23andMe> provides a parser for 23andMe files.

=head1 Constructor

New L<GAL::Parser::23andMe> objects are created by the class method
new.  Arguments should be passed to the constructor as a list (or
reference) of key value pairs.  All attributes of the Parser object
can be set in the call to new. An simple example of object creation
would look like this:

    my $parser = GAL::Parser::23andMe->new(file => 'snp.txt');

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
     Usage   : GAL::Parser::23andMe->new();
     Function: Creates a GAL::Parser::23andMe object;
     Returns : A GAL::Parser::23andMe object
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
 Args    : A hash ref of fields that this sub can understand.

=cut

sub parse_record {
	my ($self, $record) = @_;



	# rsid chromosome position genotype
	# rs3094315  1 742429  AA
	# rs12562034 1 758311  GG
	# rs3934834  1 995669  CC
	# rs9442372  1 1008567 AA
	# rs3737728  1 1011278 AA
	
	my $seqid      = 'chr' . $record->{chromosome};
	my $source     = '23andMe';
	my $type       = 'SNV';
	my $start      = $record->{position};
	my $end        = $record->{position};
	my $score      = '.';
	my $strand     = '+';
	my $phase      = '.';
	my $id         = $record->{rsid};

        my $reference_seq ||= uc $self->fasta->seq($seqid, $start, $end);
	my %variant_hash  = map {$_ => 1} split '', $record->{genotype};
	my @variant_seqs = keys %variant_hash;

	return undef if @variant_seqs == 1 && $variant_seqs[0] eq $reference_seq;
	
	my $zygosity = scalar @variant_seqs > 1 ? 'heterozygous' : 'homozygous';

	my $attributes = {Variant_seq   => \@variant_seqs,
			  Reference_seq => [$reference_seq],
			  ID            => [$id],
			  Zygosity      => [$zygosity],
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
	my @field_names = qw(rsid chromosome position genotype);
	my $reader = GAL::Reader::DelimitedLine->new(field_names => \@field_names);
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

<GAL::Parser::23andMe> requires no configuration files or environment variables.

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
