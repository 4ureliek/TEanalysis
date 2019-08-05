package GAL::Run;

use strict;
use vars qw($VERSION);


$VERSION = 0.2.0;
use base qw(GAL::Base);

=head1 NAME

GAL::Run - <One line description of module's purpose here>

=head1 VERSION

This document describes GAL::Run version 1.0.4

=head1 SYNOPSIS

     use GAL::Run;

=for author to fill in:
     Brief code example(s) here showing commonest usage(s).
     This section will be as far as many users bother reading
     so make it as educational and exemplary as possible.

=head1 DESCRIPTION

=for author to fill in:
     Write a full description of the module and its features here.
     Use subsections (=head2, =head3) as appropriate.

=head1 METHODS

=cut

#-----------------------------------------------------------------------------
#                                 Constructor
#-----------------------------------------------------------------------------

=head2 new

     Title   : new
     Usage   : GAL::Run->new();
     Function: Creates a Run object;
     Returns : A Run object
     Args    :

=cut

sub new {
	my $class = shift;
	my $self = {};
	bless $self, $class;
	$self->_initialize_args(@_);
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
	# Set valid class attributes here
	my @valid_attributes = qw(path command stdout stderr );
	$self->set_attributes(\@args, \@valid_attributes);
	######################################################################
}

#-----------------------------------------------------------------------------
#                                 Attributes
#-----------------------------------------------------------------------------

=head2 path

 Title   : path
 Usage   : $a = $self->path();
 Function: Get/Set the value of path.
 Returns : Value of path.
 Args    : A value to set for path.

=cut

sub path {
    my ($self, $value) = @_;

    $self->{path} = $value if defined $value;
    return $self->{path};
}

#-----------------------------------------------------------------------------

=head2 command

 Title   : command
 Usage   : $a = $self->command();
 Function: Get/Set the value of command.
 Returns : Value of command.
 Args    : A value to set for command.

=cut

sub command {
    my ($self, $value) = @_;

    $self->{command} = $value if defined $value;
    return $self->{command};
}

#-----------------------------------------------------------------------------

=head2 stdout

 Title   : stdout
 Usage   : $a = $self->stdout();
 Function: Get/Set the path of the file where STDOUT is redirected.
 Returns : The path of the file where STDOUT is redirected.
 Args    : The path of the file where STDOUT is redirected.

=cut

sub stdout {
    my ($self, $value) = @_;

    $self->{stdout} = $value if defined $value;
    $self->{stdout} ||= 'gal_run_stdout_' . $self->random_string(8);

    return $self->{stdout};
}

#-----------------------------------------------------------------------------

=head2 stderr

 Title   : stderr
 Usage   : $a = $self->stderr();
 Function: Get/Set the path of the file where STDERR is redirected.
 Returns : The path of the file where STDERR is redirected.
 Args    : The path of the file where STDERR is redirected.

=cut

sub stderr {
    my ($self, $value) = @_;

    $self->{stderr} = $value if defined $value;
    $self->{stderr} ||= 'gal_run_stderr_' . $self->random_string(8);
    return $self->{stderr};
}

#-----------------------------------------------------------------------------
#                                    Methods
#-----------------------------------------------------------------------------

=head2 run

 Title   : run
 Usage   : $a = $self->run();
 Function: Run a system command
 Returns : The return value of the command run
 Args    : A hash/ref of the following:
	     cl_args => A scalar or arrayref of
			command line arguments.
	     path    => The path to the command
			executable.
	     command => The command to run.
	     stdout  => The name of a file
			for STDOUT output.
	     stderr  => The name of a file
			for STDERR output.

=cut

sub run {

  my ($self, @args_array) = @_;

  $self->clean_up;

  my %args = ref $args_array[0] ? %{$args_array[0]} : @args_array;

  my $path    = $self->path($args{path});
  my $command = $self->command($args{command});
  my $stdout  = $self->stdout($args{stdout});
  my $stderr  = $self->stderr($args{stderr});

  my @cl_args = ref $args{cl_args} ? @{$args{cl_args}} : $args{cl_args};

  $command = $path . '/' . $command;
  $command =~ s|/+|/|g;
  $command .= join ' ', ('', @cl_args);
  $command .= " > $stdout 2> $stderr";

  $self->debug('running_command', $command);
  my $rv = system($command);

  $self->return_value($rv);
  return $rv;
}

#-----------------------------------------------------------------------------

=head2 get_stdout

 Title   : get_stdout
 Usage   : $a = $self->get_stdout();
 Function: Return STDOUT from the command as a scalar.
 Returns : A scalar containing the output to STDOUT.
 Args    : N/A

=cut

sub get_stdout {
  my $self = shift;
  if (! $self->{stdout_contents}) {
    my $fh = $self->stdout_fh();
    $self->{stdout_contents} = join '', (<$fh>);
  }
  return $self->{stdout_contents};
}

#-----------------------------------------------------------------------------

=head2 get_stderr

 Title   : get_stderr
 Usage   : $a = $self->get_stderr();
 Function: Return STDERR from the command as a scalar.
 Returns : A scalar containing the errput to STDERR.
 Args    : N/A

=cut

sub get_stderr {
  my $self = shift;
  if (! $self->{stderr_contents}) {
    my $fh = $self->stderr_fh();
    $self->{stderr_contents} = join '', (<$fh>);
  }
  return $self->{stderr_contents};
}

#-----------------------------------------------------------------------------

=head2 stdout_fh

 Title   : stdout_fh
 Usage   : $a = $self->stdout_fh();
 Function: Retuns a filehandle to the contents of STDOUT
 Returns : A filehandle
 Args    : N/A

=cut

sub stdout_fh {
  my $self = shift;

  if (! $self->{stdout_fh}) {
    my $file = $self->{stdout};
    open(my $IN, '<', $file) or
      $self->throw('cant_open_file_for_reading', $file);
    $self->{stdout_fh} = $IN;
  }
  return $self->{stdout_fh};
}

#-----------------------------------------------------------------------------

=head2 stderr_fh

 Title   : stderr_fh
 Usage   : $a = $self->stderr_fh();
 Function: Retuns a filehandle to the contents of STDERR
 Returns : A filehandle
 Args    : N/A

=cut

sub stderr_fh {
  my $self = shift;

  if (! $self->{stderr_fh}) {
    my $file = $self->{stderr};
    open(my $IN, '<', $file) or
      $self->throw('cant_open_file_for_reading', $file);
    $self->{stderr_fh} = $IN;
  }
  return $self->{stderr_fh};
}

#-----------------------------------------------------------------------------

=head2

 Title   : return_value
 Usage   : $a = $self->return_value();
 Function: Returns the return value of the last executed command
 Returns : A scalar
 Args    : The return value to set

=cut

sub return_value {
    my ($self, $value) = @_;

    $self->{return_value} = $value if defined $value;
    return $self->{return_value};
}

#-----------------------------------------------------------------------------

=head2

 Title   : clean_up
 Usage   : $a = $self->clean_up();
 Function: Deletes temp files for STDOUT and STDERR from most recent run.
           Will also delete any additional files passed.
 Returns : N/A
 Args    : An array of filenames to delete.

=cut

sub clean_up {
    my ($self, @files) = @_;

    system("rm -f $self->{stdout}") if -e $self->{stdout};
    system("rm -f $self->{stderr}") if -e $self->{stderr};
    delete @{$self}{qw(return_value stderr_contents stdout_contents
		       stdout_fh stderr_fh)};
    for my $file (@files) {
	next unless -e $file;
	system("rm -f $file");
    }
}

#-----------------------------------------------------------------------------

=head2

 Title   : method
 Usage   : $a = $self->method();
 Function:
 Returns :
 Args    :

=cut

sub method {
    my $self = shift;
    $self->throw(message => 'This method is not yet implimented');
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

<GAL::Run> requires no configuration files or environment variables.

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
