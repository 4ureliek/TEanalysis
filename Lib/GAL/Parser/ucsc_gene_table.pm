package GAL::Parser::ucsc_gene_table;

use strict;
use vars qw($VERSION);
$VERSION = 0.2.0;

use base qw(GAL::Parser);
use GAL::Reader::DelimitedLine;

=head1 NAME

GAL::Parser::ucsc_gene_table - Parse UCSC gene table files.

    knownGene.txt
    refGene.txt
    ccdsGene.txt
    wgEncodeGencodeManualV4.txt
    vegaGene.txt
    ensGene.txt

=head1 VERSION

This document describes GAL::Parser::ucsc_gene_table version 0.2.0

=head1 SYNOPSIS

    my $parser = GAL::Parser::ucsc_gene_table->new(file => 'ucsc_gene_table.txt');

    while (my $feature_hash = $parser->next_feature_hash) {
	print $parser->to_gff3($feature_hash) . "\n";
    }

=head1 DESCRIPTION

L<GAL::Parser::ucsc_gene_table> provides a parser for UCSC_GENE_TABLE data.

=head1 Constructor

New L<GAL::Parser::ucsc_gene_table> objects are created by the class method
new.  Arguments should be passed to the constructor as a list (or
reference) of key value pairs.  All attributes of the
L<GAL::Parser::ucsc_gene_table> object can be set in the call to new. An
simple example of object creation would look like this:

    my $parser = GAL::Parser::ucsc_gene_table->new(file => 'ucsc_gene_table.txt');

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
     Usage   : GAL::Parser::ucsc_gene_table->new();
     Function: Creates a GAL::Parser::ucsc_gene_table object;
     Returns : A GAL::Parser::ucsc_gene_table object
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

	# bin name chrom strand txStart txEnd cdsStart cdsEnd
	# exonCount exonStarts exonEnds score name2 cdsStartStat
	# cdsEndStat exonFrames

	my $feature_id = $record->{name};
	my $seqid      = $record->{chrom};
	my $source     = 'UCSC';
	my $type       = 'mRNA';
	my $start      = $record->{txStart};
	my $end        = $record->{txEnd};
	my $score      = $record->{score};
	my $strand     = $record->{strand};
	my $phase      = '.';

	my $attributes = {ID => [$feature_id]};

	$attributes->{Dbxref} = [$record->{name2}] if $record->{name2};

	my $transcript = {feature_id => $feature_id,
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

	my @features = ($transcript);

	my @exon_starts = split ',', $record->{exonStarts};
	my @exon_ends   = split ',', $record->{exonEnds};

	if (@exon_starts != scalar @exon_ends) {
	    $self->warn('mis_matched_exon_start_ends',
			"Mismatched exon starts and ends $name",
			);
	}

	my @exon_pairs;
	for my $i (0 .. scalar @exon_starts - 1) {
	    my ($start, $end) = ($exon_starts[$i], $exon_ends[$i]);
	    if ($start > $end) {
		$self->warn('negative_length_exon',
			    "Negative length exon ($name, $start, $end)"
			    );
	    }
	    if ($start == $end) {
		$self->warn('zero_length_exon',
			    "Zero length exon : $name, $start, $end"
			    );
	    }
	    push @exon_pairs, [$start, $end];
	}

	if ($self->any_overlaps(@exon_pairs)) {
	    $self->warn('exons_overlap',
			"Exons overlap $name",
			);
	}

	my @cds_pairs;
	for my $pair (@exon_pairs) {
	    last if $cdsEnd - $cdsStart < 3; #Dont allow a CDS < 3nt long.
	    my ($start, $end) = @{$pair};
	    next if $end   < $cdsStart;
	    last if $start > $cdsEnd;
	    $start = $cdsStart if ($start < $cdsStart && $end > $cdsStart);
	    $end   = $cdsEnd   if ($start < $cdsEnd && $end > $cdsEnd);
	    push @cds_pairs, [$start, $end];
	}

	my $exons = $self->build_child_features(parent      => $transcript,
						type        => 'exon',
						coordinates => \@exon_pairs
						);
	my $CDSs = $self->build_child_features(parent       => $transcript,
					       type         => 'CDS',
					       coordinates => \@cds_pairs
					       );

	push @features, (@{$exons}, @{$CDSs});
	return $feature;
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
    my @field_names = qw(bin name chrom strand txStart txEnd cdsStart
			 cdsEnd exonCount exonStarts exonEnds score name2
			 cdsStartStat cdsEndStat exonFrames);
    my $reader = GAL::Reader::DelimitedLine->new(field_names => \@field_names);
    $self->{reader} = $reader;
  }
  return $self->{reader};
}

#-----------------------------------------------------------------------------

sub any_overlaps {

    my ($self, $pairs) = @_;

    for my $i (0 .. scalar @{@pairs} - 1) {
	my $pair_i = $pairs->[$i];
	for my $j (0 .. scalar @pairs - 1) {
	    next if $i == $j;
	    my $pair_j = $pairs->[$j];
	    # Return 1 unless these two don't overlap
	    return 1 unless ($pair_i->[1] < $pair_j->[0] or
			     $pair_i->[0] > $pair_j->[1]);
	}
    }
    # We never overlaped so return false
    return 0;
}

#-----------------------------------------------------------------------------

sub build_child_features {

    my %args = @_;

    my $parent = $args{parent};
    my $type   = $args{type};
    my $coords = $args{coordinates};
    my $parent_id = $parent->{feature_id};

    my @features;
    my $count;
    for my $pair (@{$coords}) {
	my ($start, $end) = @{$pair};
	my %feature = %{$parent};
	my $attributes = {gene_id       => [$parent_id],
			  transcript_id => [$parent_id],
		      };

	@feature{qw(type start end attributes)} =
	    ($type, $start, $end, $attributes);

	push @features, \%feature;
    }
}
return \@features;
}

#-----------------------------------------------------------------------------

=head1 DIAGNOSTICS

L<GAL::Parser::ucsc_gene_table> does not throw any warnings or errors.

=head1 CONFIGURATION AND ENVIRONMENT

L<GAL::Parser::ucsc_gene_table> requires no configuration files or
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
