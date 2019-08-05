package GAL::Storage;

use strict;
use vars qw($VERSION);


$VERSION = 0.2.0;
use base qw(GAL::Base);
use DBI;

=head1 NAME

GAL::Storage - Base class for GAL storage engines

=head1 VERSION

This document describes GAL::Storage version 0.2.0

=head1 SYNOPSIS

my $feat_store_args = {class    => $storage,
		       dsn      => $feature_dsn,
		       user     => $user,
		       password => $password,
		   };

my $parser_args = {class => $parser,
		  };

my $fasta_args = {path => $fasta};

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
  my $mrna_id = $mrna->feature_id;
}

=head1 DESCRIPTION

The L<GAL::Storage> class provides a base class for the storage
classes available to GAL.  It is not intended to be used directly -
you should use one of it's subclasses instead.  It provides many of
the general attributes used by the sublcasses and several general
methods.

=head1 CONSTRUCTOR

New GAL::Storage::subclass objects are created by the class method
new.  Arguments should be passed to the constructor as a list (or
reference) of key value pairs.  All attributes of the Storage object
can be set in the call to new. An simple example of object creation
would look like this:

    my $parser = GAL::Storage::SQLite->new(dsn => 'dbi:SQLite:db_name);

The constructor recognizes the following parameters which will set the
appropriate attributes:

=over 4

=item * C<annotation>

This is a read only attribute that provides access to a weakened
version of the L<GAL::Annotation> object that created this storage

=item * C<dsn'>

dsn => 'dbi:SQLite:db_name

This optional parameter defines the data source name of the database
to open.  By default Storage will use and SQLite database with a
random filename, but see the comment for the database attribute below.

=item * C<scheme>

This is a read only parameter that is set to 'DBI';

=item * C<driver >

This is a read only parameter that defines the database driver
(i.e. what RDMS you are going to use).  It is set by each storage
subclass.

=item * C<database>

database => 'db_name'

This optional parameter defines the database name.  You don't need to
specify both the database name and the dsn as they both contain the
database name. Since the driver and the scheme are set by the class
you could give either the dsn or the database name it it will work.

=item * C<user => 'user_name'>

This optional parameter defines the user name for connecting to the
database.

=item * C<password => 'password'>

This optional parameter defines the password for connecting to the
database.

=back

=head1 METHODS

=cut

#-----------------------------------------------------------------------------
#-------------------------------- Constructor --------------------------------
#-----------------------------------------------------------------------------

=head2 new

     Title   : new
     Usage   : GAL::Storage->new()
     Function: Creates a Storage object;
     Returns : A Storage object
     Args    : See L</attributes> above

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
	# Set valid class attributes here
	my @valid_attributes = qw(annotation class dsn database user password);
	$self->set_attributes($args, @valid_attributes);
	######################################################################
}

#-----------------------------------------------------------------------------
#--------------------------------- Attributes --------------------------------
#-----------------------------------------------------------------------------

=head2 annotation

 Title   : annotation
 Usage   : $annotation = $self->annotation()
 Function: Get a weakened copy of the L<GAL::Annotation> object that owns
	   this storage object.
 Returns : A L<GAL::Annotation> object.
 Args    : None

=cut

sub annotation {
  my ($self, $annotation) = @_;
  $self->{annotation} = $annotation if $annotation;
  return $self->{annotation};
}

#-----------------------------------------------------------------------------

=head2 dsn

 Title   : dsn
 Usage   : $dsn = $self->dsn('dbi:SQLite:db_name')
 Function: Get/Set the value of the Database Source Name (DSN).
 Returns : A text value for the DSN
 Args    : A text value of the form scheme:driver:db_name.

=cut

sub dsn {
	my ($self, $dsn) = @_;

	if (exists  $self->{dsn} &&
	    defined $self->{dsn} &&
	    ! defined $dsn) {
	  return $self->{dsn};
	}

	# If we don't have a dsn, then create one
	$dsn ||= '';
	my ($scheme, $driver, $attr_string, $attr_hash, $driver_dsn)
	  = DBI->parse_dsn($dsn);
	$scheme = $self->scheme;
	my $self_driver = $self->driver;
	$self->warn('conflicting_db_drivers',
		    "You gave: $driver, using $self_driver instead\n"
		   )
	  if $driver && $self_driver ne $driver;
	$driver = $self->driver($driver);

	my ($database, %attributes);
	if ($driver_dsn && $driver_dsn =~ /=/) {
	  my %attributes = map {split /=/} split /;/, $driver_dsn;
	  $database = $attributes{dbname} || $attributes{database};
	}
	else {
	  $database = $driver_dsn;
	}

	$database = $self->database($database);
	my $attribute_txt = '';
	map {$attribute_txt .= $_ . '=' . $attributes{$_} . ';'} keys %attributes;

	$self->{dsn} = join ':', ($scheme, $driver, $database, $attribute_txt);
	$self->{dsn} =~ s/:$//;

	return $self->{dsn};
}

#-----------------------------------------------------------------------------

=head2 scheme

 Title   : scheme
 Usage   : $scheme = $self->scheme
 Function: Get the value of scheme.
 Returns : The value 'DBI'.
 Args    : None

=cut

sub scheme {
  return 'DBI';
}

#-----------------------------------------------------------------------------

=head2 driver

 Title   : driver
 Usage   : $driver = $self->driver
 Function: Get  the value of driver.
 Returns : The value 'SQLite'
 Args    : None

=cut

sub driver {
  return 'SQLite';
}

#-----------------------------------------------------------------------------

=head2 database

 Title   : database
 Usage   : $db = $self->database('db_name');
 Function: Get/Set the value of the database name.
 Returns : The value of the database name.
 Args    : A text value to set the database name to.

=cut

sub database {
      my ($self, $database) = @_;

      $self->{database} = $database if $database;
      $self->{database} ||= join '_', ('gal_database',
				       $self->time_stamp,
				       $self->random_string(8)
				      );
      return $self->{database};
}

#-----------------------------------------------------------------------------

=head2 user

 Title   : user
 Usage   : $a = $self->user()
 Function: Get/Set the value of the database user.
 Returns : The value of user.
 Args    : A value to set user to.

=cut

sub user {
	my ($self, $user) = @_;
	$self->{user} = $user if $user;
	return $self->{user};
}

#-----------------------------------------------------------------------------

=head2 password

 Title   : password
 Usage   : $a = $self->password()
 Function: Get/Set the value of the database password.
 Returns : The value of password.
 Args    : A value to set password to.

=cut

sub password {
	my ($self, $password) = @_;
	$self->{password} = $password if $password;
	return $self->{password};
}

#-----------------------------------------------------------------------------
#---------------------------------- Methods ----------------------------------
#-----------------------------------------------------------------------------

sub _load_schema {
  my ($self, $dbh) = @_;

  $self->throw('subclass_must_override_this_method',
	       ('Method _load_schema must be implimented by ' .
		'subclass'),
	      );
}

#-----------------------------------------------------------------------------

sub load_files {
  my $self = shift;
  $self->throw('subclass_must_override_this_method',
	       'Method load files must be implimented by sublcass',
	      );
}

#-----------------------------------------------------------------------------

=head2 add_features_to_buffer

 Title   : add_features_to_buffer
 Usage   : $self->add_features_to_buffer(@features)
 Function: Add feature hashes to a buffer for loading into storage.
 Returns : N/A
 Args    : A reference to an array of feature hashes.
 Alias   : features2buffer

=head2 features2buffer

An alias for add_features_to_buffer.

=cut

sub features2buffer {shift->add_features_to_buffer(@_)}

sub add_features_to_buffer {

	my ($self, $features) = @_;

	$features = ref $features eq 'ARRAY' ? $features : [$features];

	$self->{_feature_buffer} ||= [];

	#my $max_feature_buffer = $self->config('MAX_FEATURE_BUFFER')
	my $max_feature_buffer = 10_000;
	if (scalar @{$self->{_feature_buffer}} + scalar @{$features} >= $max_feature_buffer) {
	  push @{$self->{_feature_buffer}}, @{$features};
	  $self->flush_buffer;
	}
	else {
	  push @{$self->{_feature_buffer}}, @{$features};
	}
}

#-----------------------------------------------------------------------------

=head2 flush_buffer

 Title   : flush_buffer
 Usage   : $self->flush_buffer()
 Function: Get/Set value of flush_buffer.
 Returns : Value of flush_buffer.
 Args    : Value to set add_feature to.

=cut

sub flush_buffer {

	my $self = shift;

	$self->{_feature_buffer} ||= [];
	if (scalar @{$self->{_feature_buffer}}) {
		$self->add_features($self->{_feature_buffer});
		$self->{_feature_buffer} = undef;
	}
}

#-----------------------------------------------------------------------------

=head2 prepare_features

 Title   : prepare_features
 Usage   : $self->prepare_features()
 Function: Normalizes feature hashes produced by the parsers and seperates
	   the attributes and relationships for bulk insert into the database;
 Returns : An array reference of feature data and an array reference of
	   attribute data.
 Args    : A feature hash or array of feature hashes.

=cut

sub prepare_features {

  my ($self, $feature_hashes) = @_;

  $feature_hashes = ref $feature_hashes eq 'ARRAY' ? $feature_hashes :
    [$feature_hashes];

  my (@features, @attributes, @relationships);

  my $feature_idx;
  for my $feature_hash (@{$feature_hashes}) {
    my $feature_id = $feature_hash->{feature_id};
    my ($bin) = $self->get_feature_bins($feature_hash);
    my $attributes = $feature_hash->{attributes};
    my @parents = ref $attributes->{Parent} eq 'ARRAY' ? @{$attributes->{Parent}} : ();
    my @feature_data = (++$feature_idx, @{$feature_hash}{qw(feature_id seqid source
							    type start end score strand
							    phase)},
			$bin);
    push @features, \@feature_data;

    for my $key (keys %{$attributes}) {
      my @values = @{$attributes->{$key}};
      for my $value (@values) {
	#push @attributes, [undef, $feature_id, $key, $value];
	push @attributes, [$feature_idx, $key, $value];
      }
    }

    for my $parent (@parents) {
      my @relationship_data = ($parent, $feature_id);
      push @relationships, \@relationship_data;
    }
  }
  return (\@features, \@attributes, \@relationships);
}

#-----------------------------------------------------------------------------

sub add_features {
  my $self = shift;
  $self->throw('subclass_must_override_this_method : add_features',
	       ('Method must be implimented by subclass : ' .
	        'add_features'),
	      );
}

#-----------------------------------------------------------------------------

sub create_database {
  my $self = shift;
  $self->throw('subclass_must_override_this_method',
	       'Method must be implimented by subclass',
	      );
}

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#
# The methods below a turned off because for now we're using only DBIx::Class
# for object retrival, but I plan to write direct access to the data soon to
# avoid the overhead imposed by DBIx::Class
#
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#
#=head2 get_children
#
# Title   : get_children
# Usage   : $self->get_children()
# Function: Get/Set value of get_children.
# Returns : Value of get_children.
# Args    : Value to set get_children to.
#
#=cut
#
#sub get_children {
#	my $self = shift;
#	$self->throw('error_code',
#                    'Method must be implimented by subclass'
#		    );
#}
#
##-----------------------------------------------------------------------------
#
#=head2 get_children_recursively
#
# Title   : get_children_recursively
# Usage   : $self->get_children_recursively()
# Function: Get/Set value of get_children_recursively.
# Returns : Value of get_children_recursively.
# Args    : Value to set get_children_recursively to.
#
#=cut
#
#sub get_children_recursively {
#  my $self = shift;
#  $self->throw('error_code',
#                'Method must be implimented by subclass'
#	      );
#}
#
##-----------------------------------------------------------------------------
#
#=head2 get_parents
#
# Title   : get_parents
# Usage   : $self->get_parents()
# Function: Get/Set value of get_parents.
# Returns : Value of get_parents.
# Args    : Value to set get_parents to.
#
#=cut
#
#sub get_parents {
#  my $self = shift;
#  $self->throw('error_code',
#                'Method must be implimented by subclass'
#	      );
#}
#
##-----------------------------------------------------------------------------
#
#=head2 get_parents_recursively
#
# Title   : get_parents_recursively
# Usage   : $self->get_parents_recursively()
# Function: Get/Set value of get_parents_recursively.
# Returns : Value of get_parents_recursively.
# Args    : Value to set get_parents_recursively to.
#
#=cut
#
#sub get_parents_recursively {
#  my $self = shift;
#  $self->throw('error_code',
#                'Method must be implimented by subclass'
#	      );
#}
#
##-----------------------------------------------------------------------------
#
#=head2 get_all_features
#
# Title   : get_all_features
# Usage   : $self->get_all_features()
# Function: Get/Set value of get_all_features.
# Returns : Value of get_all_features.
# Args    : Value to set get_all_features to.
#
#=cut
#
#sub get_all_features {
#  my $self = shift;
#  $self->throw('error_code',
#                'Method must be implimented by subclass'
#	      );
#}
#
##-----------------------------------------------------------------------------
#
#=head2 get_features_by_type
#
# Title   : get_features_by_type
# Usage   : $self->get_features_by_type()
# Function: Get/Set value of get_features_by_type.
# Returns : Value of get_features_by_type.
# Args    : Value to set get_features_by_type to.
#
#=cut
#
#sub get_features_by_type {
#  my $self = shift;
#  $self->throw('error_code',
#                'Method must be implimented by subclass'
#	      );
#}
#
##-----------------------------------------------------------------------------
#
#=head2 get_recursive_features_by_type
#
# Title   : get_recursive_features_by_type
# Usage   : $self->get_recursive_features_by_type()
# Function: Get/Set value of get_recursive_features_by_type.
# Returns : Value of get_recursive_features_by_type.
# Args    : Value to set get_recursive_features_by_type to.
#
#=cut
#
#sub get_recursive_features_by_type {
#  my $self = shift;
#  $self->throw('error_code',
#                'Method must be implimented by subclass'
#	      );
#}
#
##-----------------------------------------------------------------------------
#
#=head2 get_feature_by_id
#
# Title   : get_feature_by_id
# Usage   : $self->get_feature_by_id()
# Function: Get/Set value of get_feature_by_id.
# Returns : Value of get_feature_by_id.
# Args    : Value to set get_feature_by_id to.
#
#=cut
#
#sub get_feature_by_id {
#  my $self = shift;
#  $self->throw('error_code',
#                'Method must be implimented by subclass'
#	      );
#}
#
##-----------------------------------------------------------------------------
#
#=head2 filter_features
#
# Title   : filter_features
# Usage   : $self->filter_features()
# Function: Get/Set value of filter_features.
# Returns : Value of filter_features.
# Args    : Value to set filter_features to.
#
#=cut
#
#sub filter_features {
#  my $self = shift;
#  $self->throw('error_code',
#                'Method must be implimented by subclass'
#		);
#}
#
##-----------------------------------------------------------------------------
#
#=head2 foo
#
# Title   : foo
# Usage   : $a = $self->foo()
# Function: Get/Set the value of foo.
# Returns : The value of foo.
# Args    : A value to set foo to.
#
#=cut
#
#sub foo {
#	my ($self, $value) = @_;
#	$self->{foo} = $value if defined $value;
#	return $self->{foo};
#}
#
##-----------------------------------------------------------------------------

=head1 DIAGNOSTICS

=over

=item C<subclass_must_override_this_method>

<GAL::Storage> subclasses are required to override a number of methods.  You
have called one of those methods and the author of the subclass has not
written a method to override it - you should complain bitterly to the author
of the subclass.

=back

=head1 CONFIGURATION AND ENVIRONMENT

<GAL::Storage> requires no configuration files or environment variables.

=head1 DEPENDENCIES

L<GAL::Base>
L<DBI>

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
