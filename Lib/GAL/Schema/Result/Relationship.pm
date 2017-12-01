package GAL::Schema::Result::Relationship;
use base qw/DBIx::Class/;

=head1 NAME

GAL::Schema::Result::Relationship - Access to relationships for <GAL::Schema::Result::Feature> objects.

=head1 VERSION

This document describes GAL::Schema::Result::Relationship version 0.2.0

=head1 SYNOPSIS

    my $exons = $self->children->search({type => 'exon'});
    my $transcripts = $self->parents->search({type => 'mRNA'});

=head1 DESCRIPTION

The GAL::Schema::Result::Relationship class is not intended for public
used by GAL::Schema::Result::Feature to find parents and children of
features.

=cut

#-----------------------------------------------------------------------------

__PACKAGE__->load_components(qw/ Core /);
__PACKAGE__->table('relationship');
__PACKAGE__->add_columns(qw/ parent child /);
__PACKAGE__->set_primary_key(qw /parent child /);
__PACKAGE__->belongs_to('your_parents'  => 'GAL::Schema::Result::Feature', {'foreign.feature_id' => 'self.parent'});
__PACKAGE__->belongs_to('your_children' => 'GAL::Schema::Result::Feature', {'foreign.feature_id' => 'self.child'});

#-----------------------------------------------------------------------------

=head1 METHODS

=head2 parent

 Title   : parent
 Usage   : $parent = $self->parent
 Function: Get the features parent.
 Returns : The value of the parent as text.
 Args    : None

=head2 child

 Title   : child
 Usage   : $child = $self->child
 Function: Get the features child.
 Returns : The value of the child as text.
 Args    : None

=head2 relationship

 Title   : relationship
 Usage   : $relationship = $self->relationship
 Function: Get the features relationship.
 Returns : The value of the relationship as text.
 Args    : None
 Note    : The relationship column in the GAL schema - and
           hence this method - are not currently in use.
           They are here for future developement.

=head1 DIAGNOSTICS

<GAL::Schema::Result::Relationship> currently throws no errors or warning messages.

=head1 CONFIGURATION AND ENVIRONMENT

<GAL::Schema::Result::Relationship> requires no configuration files or environment variables.

=head1 DEPENDENCIES

<DBIx::Class>

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
