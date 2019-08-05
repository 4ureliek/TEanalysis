package GAL::Storage::SQLite;

use strict;
use vars qw($VERSION);

$VERSION = 0.2.0;
use base qw(GAL::Storage);

=head1 NAME

GAL::Storage::SQLite - SQLite feature storage for GAL

=head1 VERSION

This document describes GAL::Storage::SQLite version 0.2.0

=head1 SYNOPSIS

    use GAL::Storage::SQLite;
    my $storage = GAL::Storage::SQLite->new(dsn => 'dbi:SQLite:db_name');

=head1 DESCRIPTION

The L<GAL::Storage::SQLite> class provides SQLite based storage to
GAL.

=head1 CONSTRUCTOR

New GAL::Storage::SQLite objects are created by the class method new.
Arguments should be passed to the constructor as a list (or reference)
of key value pairs.  All attributes of the Storage object can be set
in the call to new. An simple example of object creation would look
like this:

    my $parser = GAL::Storage::SQLite->new(dsn => 'dbi:SQLite:db_name);

The constructor recognizes the following parameters which will set the
appropriate attributes:

=over 4

=item * C<scheme>

This is a read only parameter that is set to 'DBI';

=item * C<driver >

This is a read only parameter that is set to 'SQLite';

=item * C<database>

database => 'db_name'

This optional parameter defines the database name.  You don't need to
specify both the database name and the dsn as they both contain the
database name. Since the driver and the scheme are set by the class
you could give either the dsn or the database name and it will work.

B<The following attributes are inherited from> L<GAL::Storage>

=item * C<annotation>

This is a read only attribute that provides access to a weakened
version of the L<GAL::Annotation> object that created this storage

=item * C<dsn => 'dbi:SQLite:db_name'>

dsn => 'dbi:SQLite:db_name

This optional parameter defines the data source name of the database
to open.  By default Storage will use and SQLite database with a
random filename, but see the comment for the database attribute below.

=item * C<user => 'user_name'>

This optional parameter defines the user name for connecting to the
database.

=item * C<password => 'password'>

This optional parameter defines the password for connecting to the
database.

=back

=head1 CONSTRUCTOR

=cut

#-----------------------------------------------------------------------------

=head2 new

     Title   : new
     Usage   : GAL::Storage::SQLite->new();
     Function: Creates a Storage object;
     Returns : A L<GAL::Storage::SQLite> object
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
	# Set valid class attributes here
	my @valid_attributes = qw();
	$self->set_attributes($args, @valid_attributes);
	######################################################################
}

=head1 Attributes

=cut

#-----------------------------------------------------------------------------

=head2 mode

     Title   : mode
     Usage   : GAL::Storage::SQLite->mode('append');
     Function: Get/set the value of mode. This defines
               whether the database is overwritten or
               appended.
     Returns : The value of mode
     Args    : The valude append or the default overwrite.

=cut

sub mode {

  my ($self, $mode) = @_;

  $self->{mode} = $mode if $mode;
  $self->{mode} ||= 'overwrite';
  return $self->{mode};
}

#-----------------------------------------------------------------------------

=head2 scheme

 Title   : scheme
 Usage   : $a = $self->scheme();
 Function: Return the RDMS scheme.
 Returns : DBI
 Args    : N/A

=cut

sub scheme {'dbi'}

#-----------------------------------------------------------------------------

=head2 driver

 Title   : driver
 Usage   : $a = $self->driver();
 Function: Get the value of driver.
 Returns : SQLite
 Args    : N/A

=cut

sub driver {'SQLite'}

#-----------------------------------------------------------------------------

=head2 dbh

 Title   : dbh
 Usage   : $a = $self->dbh();
 Function: Get/Set the database handle
 Returns : The a DBI::SQLite object
 Args    : A DBI::SQLite object

=cut

sub dbh {
  my $self = shift;
  $self->{dbh} ||= $self->_open_or_create_database;
  return $self->{dbh};
}

#-----------------------------------------------------------------------------

=head2 database

 Title   : database
 Usage   : $a = $self->database();
 Function: Get/Set the database name.
 Returns : A database name.
 Args    : A database name.

=cut

sub database {
  my ($self, $database) = @_;

  if (! $database && ! $self->{database}) {
    $database = $self->random_string . '.sqlite';
  }
  $self->{database} = $database if $database;
  return $self->{database};
}

#-----------------------------------------------------------------------------

# Title   : _load_schema
# Usage   : $a = $self->_load_schema();
# Function: Load the database schema
# Returns : N/A
# Args    : N/A

sub _load_schema {
  my ($self, $dbh) = @_;

  $dbh->do("DROP TABLE IF EXISTS feature");
  $dbh->do("DROP TABLE IF EXISTS attribute");
  $dbh->do("DROP TABLE IF EXISTS relationship");
  $dbh->do("CREATE TABLE feature ("    .
	   "feature_idx INTEGER PRIMARY KEY AUTOINCREMENT, " .
	   "feature_id 	TEXT, "    .
	   "seqid      	TEXT, "    .
	   "source     	TEXT, "    .
	   "type       	TEXT, "    .
	   "start      	INTEGER, " .
	   "end        	INTEGER, " .
	   "score      	TEXT, "    .
	   "strand     	TEXT, "    .
	   "phase      	TEXT,"     .
	   "bin        	TEXT)"
	  );
  $dbh->do("CREATE TABLE attribute ("  .
	   "feature_idx  TEXT, " .
	   "att_key      TEXT, " .
	   "att_value    TEXT)"
	  );
  $dbh->do("CREATE TABLE relationship ("  .
	   "parent TEXT, " .
	   "child  TEXT, " .
	   "UNIQUE(parent, child) " .
	   "ON CONFLICT REPLACE)"
	  );
}

#-----------------------------------------------------------------------------

=head2 drop_database

  Title   : drop_database
  Usage   : $self->drop_database();
  Function: Drop the database;
  Returns : N/A
  Args    : N/A

=cut

sub drop_database {

  my $self = shift;

  my $database = $self->database;
  if (-e $database) {
    $self->warn('dropping_database', $database);
    `rm $database`;
  }
}

#-----------------------------------------------------------------------------

# Title   : _open_or_create_database
# Usage   : $self->_open_or_create_database();
# Function: Create the database if it doesnt exists.
# Returns : A DBD::SQLite object
# Args    : Nothing

sub _open_or_create_database {
  my $self = shift;

  my $dbh;
  my $dsn = $self->dsn;
  my $database = $self->database;

  if (! -e $database) {
    $dbh = DBI->connect($self->dsn);
    # This one is supposed to speed up write significantly, but doesn't
    # $dbh->do('PRAGMA default_synchronous = OFF');
    # $dbh->do('PRAGMA synchronous = 0');
    $self->_load_schema($dbh);
    # http://search.cpan.org/~adamk/DBD-SQLite-1.29/lib/DBD/SQLite.pm#Transactions
    # $dbh->{AutoCommit} = 0;
  }
  else {
    $self->warn('using_existing_database',
		("A database by the name $database already " .
		 "existed, so I'll use it as is."),
	       );
    $dbh = DBI->connect($self->dsn);
    # $dbh = DBI->connect("dbi:SQLite:dbname=:memory");
    # $dbh->sqlite_backup_from_file($self->database);
  }
  return $dbh
}

#-----------------------------------------------------------------------------

=head2 index_database

 Title   : index_database
 Usage   : $self->index_database();
 Function: Create indeces on the database.
 Returns : Nothing
 Args    : None

=cut

sub index_database {

  my $self = shift;
  my $dbh  = $self->dbh;

  $self->info('indexing_database', $self->database);
  # Create feature indeces
  $self->info('indexing_feature_id', $self->database);
  $dbh->do("CREATE INDEX feat_id_index ON feature (feature_id)");
  $self->info('indexing_locus', $self->database);
  $dbh->do("CREATE INDEX feat_locus_index ON feature (bin, seqid, start, end)");
  $self->info('indexing_feature_type', $self->database);
  $dbh->do("CREATE INDEX feat_type_index ON feature (type)");

  # Create attribute indeces
  $self->info('indexing_feature_attributes', $self->database);
  $dbh->do("CREATE INDEX att_feature_idx_index ON attribute (feature_idx, att_key, att_value)");

  # Create relationship indeces
  $self->info('indexing_relationships', $self->database);
  $dbh->do("CREATE INDEX rel_parent_child_index ON relationship (parent, child)");
  $dbh->do("CREATE INDEX rel_child_parent_index ON relationship (child, parent)");

  # Analyze
  $self->info('analyzing_database', $self->database);
  $dbh->do("ANALYZE");

  $self->info('finished_indexing_database', $self->database);
  return 1;
}

#-----------------------------------------------------------------------------

=head2 load_files

 Title   : load_files
 Usage   : $self->load_files();
 Function: Load a file(s) into the database
 Returns : Nothing
 Args    : An scalar string or array reference containing the name(s) of
	   files to load.

=cut

# TODO: Allow GAL::Storage::subclass::load_files to set the parser.

sub load_files {

  my ($self, $files) = @_;
  $files = ref $files eq 'ARRAY' ? $files : [$files];

  my $parser = $self->annotation->parser;

  $self->drop_database if $self->mode eq 'overwrite';

  for my $file (@{$files}) {
    $parser->file($file);
    $self->info('loading_database', $file);
    my $counter = 0;
    while (my $feature = $parser->next_feature_hash) {
      $self->add_features_to_buffer($feature);
      if (++$counter % 10000 == 0) {
	$self->info('loading_database', "Loaded $counter features");
      }
    }
    $self->flush_buffer;
    $self->info('finished_loading_database', $file);
  }
  $self->index_database;
}

#-----------------------------------------------------------------------------

=head2 add_features

 Title   : add_features
 Usage   : $self->add_features();
 Function: Add features to the database.
 Returns : Nothing
 Args    : An array reference of feature hashes.

=cut

sub add_features {
  my ($self, $features) = @_;

  my ($feat_rows, $att_rows, $rel_rows) = $self->prepare_features($features);
  my $dbh = $self->dbh;

  # "feature_idx INT, " .
  # "feature_id  VARCHAR(255), " .
  # "seqid       VARCHAR(255), " .
  # "source      VARCHAR(255), " .
  # "type        VARCHAR(255), " .
  # "start       INT, "          .
  # "end         INT, "          .
  # "score       VARCHAR(255), " .
  # "strand      VARCHAR(1), "   .
  # "phase       VARCHAR(1),"    .
  # "bin         VARCHAR(15))"

  my $feat_stmt = "INSERT INTO feature VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)";
  my $feat_sth  = $dbh->prepare($feat_stmt);

  my %feat_idx_map;

  # http://search.cpan.org/~adamk/DBD-SQLite-1.29/lib/DBD/SQLite.pm#Transactions
  $dbh->begin_work;
  for my $feat_row (@{$feat_rows}) {
    my $old_feat_idx = $feat_row->[0];
    $feat_row->[0] = undef;
    my $rv = $feat_sth->execute(@{$feat_row});
    if ($rv) {
      my ($auto_inc) = $dbh->selectrow_array('SELECT last_insert_rowid()');
      $feat_idx_map{$old_feat_idx} = $auto_inc;
    }
    else {
      my $warn_message = ('The following data failed to insert into the ' .
			  'feature table : ');
      my $warn_code = 'bad_feature_table_insert :';
      my $data = join ', ', @{$feat_row};
      $warn_message .= $data;
      $self->warn($warn_code,
		  $warn_message);
    }
  }

  # Removed # "attribute_id INT NOT NULL AUTO_INCREMENT, " .
  # "feature_idx VARCHAR(255), " .
  # "att_key     VARCHAR(255), "    .
  # "att_value   TEXT, "  .

  #my $att_stmt = "INSERT INTO attribute VALUES (?, ?, ?, ?)";
  my $att_stmt = "INSERT INTO attribute VALUES (?, ?, ?)";
  my $att_sth = $dbh->prepare($att_stmt);

  for my $att_row (@{$att_rows}) {
    next unless exists $feat_idx_map{$att_row->[0]};
    $att_row->[0] = $feat_idx_map{$att_row->[0]};
    my $rv = $att_sth->execute(@{$att_row});
    unless ($rv) {
      my $warn_message = ('The following data failed to insert into the ' .
			  'attribute table : ');
      my $warn_code = 'bad_attribute_table_insert :';
      my $data = join ', ', @{$att_row};
      $warn_message .= $data;
      $warn_code    .= $data;
      $self->warn($warn_code,
		  $warn_message);
    }
  }

  # "parent   VARCHAR(255), " .
  # "child    VARCHAR(255), "    .
  # REMOVED # "relationship VARCHAR(255)) "

  my $rel_stmt = "INSERT INTO relationship VALUES (?, ?)";
  my $rel_sth = $dbh->prepare($rel_stmt);

  for my $rel_row (@{$rel_rows}) {
    my $rv = $rel_sth->execute(@{$rel_row});
    unless ($rv) {
      my $warn_message = ('The following data failed to insert into the ' .
			  'relationship table : ');
      my $warn_code = 'bad_relationship_table_insert :';
      my $data = join ', ', @{$rel_row};
      $warn_message .= $data;
      $warn_code    .= $data;
      $self->warn($warn_code,
		  $warn_message);
    }
  }
  $dbh->commit;
}

#-----------------------------------------------------------------------------
#
#=head2 get_children
#
# Title   : get_children
# Usage   : $self->get_children();
# Function: Get/Set value of get_children.
# Returns : Value of get_children.
# Args    : Value to set get_children to.
#
#=cut
#
#sub get_children {
#	my $self = shift;
#	$self->not_implemented('get_children');
#}
#
##-----------------------------------------------------------------------------
#
#=head2 get_children_recursively
#
# Title   : get_children_recursively
# Usage   : $self->get_children_recursively();
# Function: Get/Set value of get_children_recursively.
# Returns : Value of get_children_recursively.
# Args    : Value to set get_children_recursively to.
#
#=cut
#
#sub get_children_recursively {
#  my $self = shift;
#  $self->not_implemented('get_children_recursively');
#}
#
##-----------------------------------------------------------------------------
#
#=head2 get_parents
#
# Title   : get_parents
# Usage   : $self->get_parents();
# Function: Get/Set value of get_parents.
# Returns : Value of get_parents.
# Args    : Value to set get_parents to.
#
#=cut
#
#sub get_parents {
#  my $self = shift;
#  $self->not_implemented('get_parents');
#}
#
##-----------------------------------------------------------------------------
#
#=head2 get_parents_recursively
#
# Title   : get_parents_recursively
# Usage   : $self->get_parents_recursively();
# Function: Get/Set value of get_parents_recursively.
# Returns : Value of get_parents_recursively.
# Args    : Value to set get_parents_recursively to.
#
#=cut
#
#sub get_parents_recursively {
#  my $self = shift;
#  $self->not_implemented('get_parents_recursively');
#}
#
##-----------------------------------------------------------------------------
#
#=head2 get_all_features
#
# Title   : get_all_features
# Usage   : $self->get_all_features();
# Function: Get/Set value of get_all_features.
# Returns : Value of get_all_features.
# Args    : Value to set get_all_features to.
#
#=cut
#
#sub get_all_features {
#  my $self = shift;
#  $self->not_implemented('get_all_features');
#}
#
##-----------------------------------------------------------------------------
#
#=head2 get_features_by_type
#
# Title   : get_features_by_type
# Usage   : $self->get_features_by_type();
# Function: Get/Set value of get_features_by_type.
# Returns : Value of get_features_by_type.
# Args    : Value to set get_features_by_type to.
#
#=cut
#
#sub get_features_by_type {
#  my $self = shift;
#  $self->not_implemented('get_features_by_type');
#}
#
##-----------------------------------------------------------------------------
#
#=head2 get_recursive_features_by_type
#
# Title   : get_recursive_features_by_type
# Usage   : $self->get_recursive_features_by_type();
# Function: Get/Set value of get_recursive_features_by_type.
# Returns : Value of get_recursive_features_by_type.
# Args    : Value to set get_recursive_features_by_type to.
#
#=cut
#
#sub get_recursive_features_by_type {
#  my $self = shift;
#  $self->not_implemented('get_recursive_features_by_type');
#}
#
##-----------------------------------------------------------------------------
#
#=head2 get_feature_by_id
#
# Title   : get_feature_by_id
# Usage   : $self->get_feature_by_id();
# Function: Get/Set value of get_feature_by_id.
# Returns : Value of get_feature_by_id.
# Args    : Value to set get_feature_by_id to.
#
#=cut
#
#sub get_feature_by_id {
#  my $self = shift;
#  $self->not_implemented('get_feature_by_id');
#}
#
##-----------------------------------------------------------------------------
#
#=head2 filter_features
#
# Title   : filter_features
# Usage   : $self->filter_features();
# Function: Get/Set value of filter_features.
# Returns : Value of filter_features.
# Args    : Value to set filter_features to.
#
#=cut
#
#sub filter_features {
#  my $self = shift;
#  $self->not_implemented('filter_features');
#}
#
##-----------------------------------------------------------------------------
#
#=head2 foo
#
# Title   : foo
# Usage   : $a = $self->foo();
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

=item C<using_existing_database>

A database by the same name specified already existed and you haven't
asked for it to be overwritten, so it will be used and possible
appended to.  If you want to overwrite this database you should use
the attribute C<mode => overwrite> when you create the storage object.

=item C<bad_feature_table_insert>

An error occured while trying to insert data into the feature table.
The offending data is given as a comma separated list.  Make sure that
all of the data is compatible with the database schema.  This may be a
problem in the parser you are using, or more likely the incoming data
is not in the correct format.

=item C<bad_attribute_table_insert>

Same problem as C<bad_feature_table_insert> above, but for the
attribute table.

=item C<bad_relationship_table_insert>

Same problem as C<bad_feature_table_insert> above, but for the
attribute table.

=back

=head1 CONFIGURATION AND ENVIRONMENT

<GAL::Storage::SQLite> requires no configuration files or environment
variables.

=head1 DEPENDENCIES

L<GAL::Storage>
L<DBI>
L<DBD::SQLite>

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
