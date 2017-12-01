package GAL::Schema::Result::Feature::five_prime_utr;

use strict;
use warnings;
use base qw(GAL::Schema::Result::Feature::sequence_feature);

=head1 NAME

GAL::Schema::Result::Feature::five_prime_utr - A five_prime_UTR object
for the GAL Library

=head1 VERSION

This document describes GAL::Schema::Result::Feature::five_prime_utr
version 0.2.0

=head1 SYNOPSIS

    use GAL::Annotation;
    my $feat_store = GAL::Annotation->new(storage => $feat_store_args,
					  parser  => $parser_args,
					  fasta   => $fasta_args,
					 );

    $feat_store->load_files(files => $feature_file,
			    mode  => 'overwrite',
			    );

    my $features = $feat_store->schema->resultset('Feature');

    my $mrnas = $features->search({type => 'mRNA'});
    while (my $mrna = $mrnas->next) {
      my $five_prime_utrs = $mrna->five_prime_UTRs;
      while (my $five_prime_utr = $five_prime_utrs->next) {
	my $utr_id  = $five_prime_utr->feature_id;
	my $start   = $five_prime_utr->start;
	my $end     = $five_prime_utr->end;
      }
    }

=head1 DESCRIPTION

<GAL::Schema::Result::Feature::five_prime_utr> provides a
<GAL::Schema::Result::Feature> subclass for five_prime_utr specific
behavior.

=head1 METHODS

There are currently no five_prime_utr specific methods implemented for
<GAL::Schema::Result::Feature::five_prime_utr>.  See
<GAL::Schema::Result::Feature> for inherited methods.

=cut

#-----------------------------------------------------------------------------



#-----------------------------------------------------------------------------

=head1 DIAGNOSTICS

<GAL::Schema::Result::Feature::five_prime_utr> currently throws no
warnings or errors.

=head1 CONFIGURATION AND ENVIRONMENT

<GAL::Schema::Result::Feature::five_prime_utr> requires no
configuration files or environment variables.

=head1 DEPENDENCIES

<GAL::Schema::Result::Feature::sequence_feature>

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
