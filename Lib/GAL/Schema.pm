package GAL::Schema;
use base qw/DBIx::Class::Schema/;

=head1 NAME

GAL::Schema - DBIx::Class functionality for the Genome Annotation Library

=head1 VERSION

This document describes GAL::Schema version 0.2.0

=head1 SYNOPSIS

    use GAL::Schema;
    my $schema = GAL::Schema->connect($dsn,
                                      $user,
                                      $password,
                                     );

=head1 DESCRIPTION

<GAL::Schema> provides access to all of the query and iteration goodness of
the excellent DBIx::Class package.  To understand how to tap into all of this
you will want to be familiar with <DBI>, <DBIx::Class>, and to a lesser extent
<SQL::Abstract>.  All of those projects are well documented on CPAN.

This particular module doesn't really do much execpt provide a base class
the modules below it, and load them up.  One additional function provided here
is to load GAL::SchemaAnnotation which provides a single method 'annotation'
so that each feature constructed by DBIx::Class can have access to a weakened
copy of the GAL::Annotation object that is managing it.

=cut

#-----------------------------------------------------------------------------


__PACKAGE__->load_namespaces();
__PACKAGE__->load_components(qw/ +GAL::SchemaAnnotation /);

#-----------------------------------------------------------------------------

=head1 DIAGNOSTICS

This package does not throw any errors or warnings itself.  If you are
getting errors that seem like they are coming from here, they are probably
being thrown by <DBIx::Class>.

=head1 CONFIGURATION AND ENVIRONMENT

<GAL::Schema> requires no configuration files or environment variables.

=head1 DEPENDENCIES

<DBIx::Class::Schema>

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
