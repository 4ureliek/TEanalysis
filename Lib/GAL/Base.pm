package GAL::Base;

use strict;
use vars qw($VERSION);
use Carp qw(croak cluck);
use Bio::DB::Fasta;
use Scalar::Util qw(blessed);

$VERSION = 0.2.1;

=head1 NAME

L<GAL::Base> - Base class for the Genome Annotation Library

=head1 VERSION

This document describes L<GAL::Base> version 0.2.1

=head1 SYNOPSIS

    use base qw(GAL::Base);

=head1 DESCRIPTION

L<GAL::Base> serves as a base class for all of the other classes in
GAL.  It is not intended to be instantiated directly, but rather to be
used with the 'use base' pragma by the other modules.  L<GAL::Base>
provides object instantiation, argument preparation and attribute
setting functions for other classes during object construction.  In
addition it provides a wide range of utility functions that are
expected to be applicable throughout the library.

=head1 INHERITS FROM

None

=head1 INHERITED BY

=over 4

=item L<Annotation.pm>

=item L<Base.pm>

=item L<Index.pm>

=item L<Interval.pm>

=item L<List.pm>

=over 4

=item L<List::Categorical.pm>

=item L<List::Numeric.pm>

=back

=item L<Manual.pm>

=item L<Parser.pm>

=over 4

=item L<Parser::23andMe.pm>

=item L<Parser::basic_snp.pm>

=item L<Parser::bed.pm>

=item L<Parser::cgi_complete.pm>

=item L<Parser::dbsnp_flat.pm>

=item L<Parser::gff3.pm>

=item L<Parser::gff3_fast.pm>

=item L<Parser::hapmap_genotypes.pm>

=item L<Parser::hgmd_indel.pm>

=item L<Parser::hgmd_snv.pm>

=item L<Parser::ncbi_blast_tab.pm>

=item L<Parser::samtools_pileup.pm>

=item L<Parser::soap_snp.pm>

=item L<Parser::solid_snp.pm>

=item L<Parser::template.pm>

=item L<Parser::template_sequence_alteration.pm>

=item L<Parser::ucsc_gene_table.pm>

=item L<Parser::vaast_smp.pm>

=item L<Parser::VCFv4_0.pm>

=back

=item L<Reader.pm>

=over 4

=item L<Reader::DelimitedLine.pm>

=item L<Reader::RecordParser.pm>

=item L<Reader::Template.pm>

=back

=item L<Run.pm>

=item L<Schema.pm>

=over 4

=item L<Schema::Result::Attribute.pm>

=item L<Schema::Result::Feature::cds.pm>

=item L<Schema::Result::Feature::exon.pm>

=item L<Schema::Result::Feature::five_prime_utr.pm>

=item L<Schema::Result::Feature::gene.pm>

=item L<Schema::Result::Feature::intron.pm>

=item L<Schema::Result::Feature::mrna.pm>

=item L<Schema::Result::Feature::protein.pm>

=item L<Schema::Result::Feature::sequence_alteration.pm>

=item L<Schema::Result::Feature::sequence_feature.pm>

=item L<Schema::Result::Feature::template.pm>

=item L<Schema::Result::Feature::three_prime_utr.pm>

=item L<Schema::Result::Feature::transcript.pm>

=item L<Schema::Result::Feature.pm>

=item L<Schema::Result::Relationship.pm>

=item L<Schema::ResultSet::Feature.pm>

=back

=item L<SchemaAnnotation.pm>

=item L<Storage.pm>

=over 4

=item L<Storage::mysql.pm>

=item L<Storage::RAM.pm>

=item L<Storage::SQLite.pm>

=back

=back

=head1 USES

=over

=item * L<Carp>

=item * L<Bio::DB::Fasta>

=item * L<Scalar::Util>

=back

=head1 CONSTRUCTOR

L<GAL::Base> is not intended to by instantiated on it's own.  It does
however, handle object creation for the rest of the library.  Each
class in GAL calls:

    my $self = $class->SUPER::new(@args);

This means that GAL::Base - at the bottom of the inheritance chain
does the actual object creation.  It creates the new object based on
the calling class, and then checks to see if there is a 'class'
argument provided.  If there is it calls the class attribute, which
L<GAL::Base> provides, and reblesses the object using the value to the
class argument as a subclass of the current calling class.  This
allows this kind of an idiom for any class that has subclasses:

    my $parser = GAL::Parser->(class => 'gff3');

That invocation would return a L<GAL::Parser::gff3> object instead of a
L<GAL::Parser> object.

L<GAL::Base> provides the following attributes to all other classes in
the GAL.  Each attribute is described full below in the
L</"Attributes"> section.

=over 4

=item * C<< fasta => '/path/to/fasta/files/' >>

This optional parameter defines a path to fasta files that are (or
will be) indexed and available to all objects in the library.  The
object created is an instance of L<Bio::DB::Fasta>.

=item * C<< class => 'subclass_name' >>

This optional parameter will cause the current object to be reblessed
as the given subclass of the calling object.

=cut

#-----------------------------------------------------------------------------
#--------------------------------- Constructor -------------------------------
#-----------------------------------------------------------------------------

=head2 new

     Title   : new
     Usage   : GAL::SomeClass->new();
     Function: Creates on object of the calling class
     Returns : An object of the calling class
     Args    : See the attributes described above.

=cut

sub new {
	my ($class, @args) = @_;
	$class = ref($class) || $class;
	my $self = bless {}, $class;

	my $args = $self->prepare_args(@args);

	# Call attribute class first so that objects get reblessed
	# into the appropriate subclass before they set their
	# attributes.
	if ($args->{class}) {
	  $self = $self->class($args->{class});
	  delete $args->{class};
	}

	$self->_initialize_args($args);
	return $self;
}

#-----------------------------------------------------------------------------

sub _initialize_args {
	my ($self, $args) = @_;

	######################################################################
	# This block of code handels class attributes.  Use the
	# @valid_attributes below to define the valid attributes for
	# this class.  You must have identically named get/set methods
	# for each attribute.  Leave the rest of this block alone!
	######################################################################
	my @valid_attributes = qw(fasta class); # Set valid class attributes here
	$self->set_attributes($args, @valid_attributes);
	######################################################################
	return $args;
}

#-----------------------------------------------------------------------------
#-------------------------------- Attributes ---------------------------------
#-----------------------------------------------------------------------------

=head1  ATTRIBUTES

All attributes can be supplied as parameters to the constructor as a
list (or referenece) of key value pairs.

=cut

=head2 fasta

  Title   : fasta
  Usage   : $fasta = $self->fasta($fasta_dir);
  Function: Provides a Bio::DB::Fasta object
  Returns : A Bio::DB::Fasta object
  Args    : A fasta file or directory of fasta files.

=cut

 sub fasta {
   my ($self, $fasta_path) = @_;

   if ($fasta_path) {
     # $fasta_path ||= $self->config('default_fasta_path');
     $self->info('fasta_database_loading_indexing', "Loading (and possibly indexing) $fasta_path");
     my $fasta = Bio::DB::Fasta->new($fasta_path);
     $self->{fasta} = $fasta;
   }
   if (! $self->{fasta}) {
     $self->throw('missing_fasta_datasource', 'Your code has requested ' .
		 'access to fasta sequence data, but no fasta ' .
		 'file has been provided');
   }
   return $self->{fasta};
 }

#-----------------------------------------------------------------------------

=head2 class

  Title   : class
  Usage   : $class = $self->class($subclass_name);
  Function: Reblesses an object into the appropriate subclass.
  Returns : A object blessed into the given subclass
  Args    : A class name of a subclass of the calling object.

=cut

 sub class {
   my ($self, $subclass) = @_;

   my $class = ref($self);
   $subclass =~ s/\.pm$//;
   $subclass =~ s/$class//;
   $subclass = $class . "::$subclass";
   $self->load_module($subclass);
   bless $self, $subclass;
   return $self;
 }

#-----------------------------------------------------------------------------
#------------------------------------ Methods --------------------------------
#-----------------------------------------------------------------------------

=head1 METHODS

=head2 verbosity

 Title   : verbosity
 Usage   : $base->verbosity($level);
 Function: Set the level of verbosity written to STDERR by the code.
 Returns : None
 Args    : Arguments can be either the words debug, info, unique, warn,
           fatal or their numerical equivalents as given below.

           debug  | 1: Print all FATAL, WARN, INFO, and DEBUG messages.  Produces
		       a lot of output.
	   info   | 2: Print all FATAL, WARN, and INFO messages. This is the
		       default.
	   unique | 3: Print only the first occurence of each error/info code.
	   warn   | 4: Don't print INFO messages.
	   fatal  | 5: Don't print INFO or WARN messages. Still dies with
		       message on FATAL errors.

=cut

sub verbosity {
  my ($self, $verbosity) = @_;

  if ($verbosity) {
    $verbosity = lc $verbosity;
    $verbosity = ($verbosity =~ /^d/ ? 1 :
		  $verbosity =~ /^i/ ? 2 :
		  $verbosity =~ /^u/ ? 3 :
		  $verbosity =~ /^w/ ? 4 :
		  $verbosity =~ /^f/ ? 5 :
		  $verbosity);
    if (! grep {$verbosity eq $_} qw(1 2 3 4 5)) {
      $self->warn('invalid_verbosity_level',
		  "$verbosity - setting verbosity level to info");
      $verbosity = 2;
    }
    $self->{verbosity} = $verbosity;
  }
  $self->{verbosity} ||= '2';
  return $self->{verbosity};
}

#-----------------------------------------------------------------------------

=head2 handle_message

 Title   : handle_message
 Usage   : $base->handle_message($level, $code, $message);
 Function: Handle a message and print to STDERR accordingly.
 Returns : None
 Args    : level  : FATAL, WARN, INFO, DEBUG
	   code   : $info_code # single_word_code_for_info
	   message: $info_msg  # Free text description of info

          Note that the 'code' above should not contain whitespace,
          should generally by lowercase unless uppercase is needed for
          clarity, and should be standardized throughout the library
          such that a code like 'file_not_found' is used consistently
          rather than 'file_not_found' in one place and
          'file_does_not_exist' in another.  All codes should be
          documented in the module and in Manual.pm

          The 'message' above can be arbitrary text and should provide
          details that will help explain this instance of the error,
          such as "Expected an integer but got the text 'ABCDEFG'
          instead".

=cut

sub handle_message {
	my ($self, $level, $code, $message) = @_;

	my $caller = ref($self);
	$level ||= 'UNKNOWN';
	$message ||= $caller;
	if (! $code) {
	  $code = 'unspecified_code';
	  $message = join '', ("GAL::Base is handling a message " .
				"for $caller without an error code.  "  .
				"Complain to the author of $caller!");
	}
	chomp $message;
	$message .= "\n";
	my $verbosity = $self->verbosity;

	if ($level eq 'FATAL') {
	  $message = join ' : ', ('FATAL', $code, $message);
	  croak $message;
	}
	elsif ($level eq 'WARN') {
	  return if ($verbosity > 4);
	  return if ($verbosity == 3 && $self->{unique_codes}{$code}++);
	  $message = join ' : ', ('WARN', $code, $message);
	  print STDERR $message;
	}
	elsif ($level eq 'INFO') {
	  return if ($verbosity > 3);
	  return if ($verbosity == 3 && $self->{unique_codes}{$code}++);
	  $message = join ' : ', ('INFO', $code, $message);
	  print STDERR $message;
	}
	elsif ($level eq 'DEBUG') {
	  return if ($verbosity != 1);
	  my $sub = (caller(4))[3];
	  #my $sub = $data[3];
	  $message = join ' : ', ('DEBUG', $code, "($caller\:\:$sub) $message");
	  print STDERR $message;
	  print '';
	}
	else {
	  $message = join '', ("GAL::Base is handling a message " .
			       "for $caller without an error level.  "  .
			       "Complain to the author of $caller!\n");
	  chomp $message;
	  $message = join ' : ', ('UNKNOWN', $code, $message);
	  croak $message;
	}
}

#-----------------------------------------------------------------------------

=head2 throw

 Title   : throw
 Usage   : $base->throw($error_code, $error_message);
 Function: Throw an error - print an error message to STDERR and die.
 Returns : None
 Args    : A text string for the error code and a text string for the
           error message. See the documentation for the method
           handle_message in this module for more details.

=cut

sub throw {
	shift->handle_message('FATAL', @_);
}

#-----------------------------------------------------------------------------

=head2 warn

 Title   : warn
 Usage   : $base->warn($warning_code, $warning_message);
 Function: Send a warning message to STDERR.
 Returns : None
 Args    : A text string for the error code and a text string for the
           error message. See the documentation for the method
           handle_message in this module for more details.

=cut

sub warn {
	shift->handle_message('WARN', @_);
}

#-----------------------------------------------------------------------------

=head2 info

 Title   : info
 Usage   : $base->info($info_code, $info_message);
 Function: Send a INFO message.
 Returns : None
 Args    : A text string for the info code and a text string for the
           info message. See the documentation for the method
           handle_message in this module for more details.

=cut

sub info {
	shift->handle_message('INFO', @_);
      }

#-----------------------------------------------------------------------------

=head2 debug

 Title   : debug
 Usage   : $base->debug($debug_code, $debug_message);
 Function: Send a DEBUG message.
 Returns : None
 Args    : A text string for the debug code and a text string for the
           debug message. See the documentation for the method
           handle_message in this module for more details.

=cut

sub debug {
	shift->handle_message('DEBUG', @_);
      }

#-----------------------------------------------------------------------------

=head2 wrap_text

 Title   : wrap_text
 Usage   : $text = $self->wrap_text($text, 50);
 Function: Wrap text to the specified column width.  Default width is 50
           characters.
 Returns : Returns wrapped text as a string.
 Args    : A string of text and an optional integer value for the wrapped
           width.

=cut

sub wrap_text {
	my ($self, $text, $cols) = @_;
	$cols ||= 50;
	$text =~ s/(.{0,$cols})/$1\n/g;
	$text =~ s/\n+$//;
	return $text;
}
#-----------------------------------------------------------------------------

=head2 trim_whitespace

 Title   : trim_whitespace
 Usage   : $trimmed_text = $self->trim_whitespace($text);
 Function: Trim leading and trailing whitespace from text;
 Returns : Returns trimmed text as a string.
 Args    : A text string.

=cut

sub trim_whitespace {
	my ($self, $text) = @_;
	$text =~ s/^\s+//;
	$text =~ s/\s+$//;
	return $text;
}

#-----------------------------------------------------------------------------

=head2 prepare_args

 Title   : prepare_args
 Usage   : $args = $self->prepare_args(@_);
 Function: Take a list of key value pairs that may be structured as an
           array, a hash or and array or hash reference and return
	   them as a hash or hash reference depending on calling
	   context.
 Returns : Hash or hash reference.
 Args    : An array, hash or reference to either.

=cut

sub prepare_args {

	my ($self, @args) = @_;

	my %args_hash;

	if (! $args[0]) {
	  # If no args are passed, don't do anything just return an empty
	  # hash(ref).
	}
	elsif (scalar @args == 1 && ref $args[0] eq 'ARRAY') {
		%args_hash = @{$args[0]};
	}
	elsif (scalar @args == 1 && ref $args[0] eq 'HASH') {
		%args_hash = %{$args[0]};
	}
	elsif (scalar @args % 2 == 0) {
		%args_hash = @args;
	}
	else {
		my $class = ref($self);
		my $err_code = 'invalid_arguments_to_prepare_args';
		my $err_msg  = ("Bad arguments passed to $class. A list "   .
				"of key value pairs or a reference to "     .
				"such a list was expected, But we got:\n"   .
				join ' ', @args);
		$self->throw($err_code, $err_msg);
	}

	return wantarray ? %args_hash : \%args_hash;
}

#-----------------------------------------------------------------------------

=head2 set_attributes

 Title   : set_attributes
 Usage   : $base->set_attributes($args, @valid_attributes);
 Function: Take a hash reference of arguments and a list (or reference) of
           valid attribute names and call the methods to set those
	   attribute values.
 Returns : None
 Args    : A hash reference of arguments and an array or array reference of
	   valid attributes names.

=cut

sub set_attributes {

	my $self = shift;
	my $args = shift;

	# Allow @valid_attributes to be passed as array or arrayref.
	my @valid_attributes = ref($_[0]) eq 'ARRAY' ? @{$_[0]} : @_;

	my $package = __PACKAGE__;
	my $caller = ref($self);

	$args = $self->prepare_args($args);

	for my $attribute (@valid_attributes) {
		next unless exists $args->{$attribute};
		if (exists $self->{$attribute}) {
			my $package = __PACKAGE__;
			my $caller = ref($self);
			my $warning_message =
			  ("$package is about to reset the value of $attribute " .
			   "on behalf of a $caller object.  This may be "   .
			   "a bad idea.");
			$self->warn('resetting_attribute_values', $warning_message);
		}
		$self->$attribute($args->{$attribute});
		delete $args->{$attribute};
	}
}

#-----------------------------------------------------------------------------

=head2 expand_iupac_nt_codes

 Title   : expand_iupac_nt_codes
 Usage   : @nucleotides = $self->expand_iupac_nt_codes('W');
 Function: Expands an IUPAC ambiguity codes to an array of nucleotides
 Returns : An array or array ref of nucleotides (ATGC-);
 Args    : An IUPAC Nucleotide ambiguity code (ACGTUMRWSYKVHDBNX-) or an
           array (or ref) of such.

=cut

sub expand_iupac_nt_codes {
	my ($self, @codes) = @_;

	my %iupac_code_map = ('A' => ['A'],
			      'C' => ['C'],
			      'G' => ['G'],
			      'T' => ['T'],
			      'U' => ['T'],
			      'M' => ['A', 'C'],
			      'R' => ['A', 'G'],
			      'W' => ['A', 'T'],
			      'S' => ['C', 'G'],
			      'Y' => ['C', 'T'],
			      'K' => ['G', 'T'],
			      'V' => ['A', 'C', 'G'],
			      'H' => ['A', 'C', 'T'],
			      'D' => ['A', 'G', 'T'],
			      'B' => ['C', 'G', 'T'],
			      'N' => ['A', 'C', 'G', 'T'],
			      'X' => ['A', 'C', 'G', 'T'],
			      '-' => ['-'],
			     );

	my @nts;
	for my $code (@codes) {
	  my $nts = $iupac_code_map{$code};
	  $self->throw('invalid_ipuac_nucleotide_code', $code)
	    unless $nts;
	  push @nts, @{$nts};
	}


	return wantarray ? @nts : \@nts;
}

#-----------------------------------------------------------------------------

=head2 load_module

 Title   : load_module
 Usage   : $base->load_module(Some::Module);
 Function: Do runtime loading (require) of a module/class.
 Returns : 1 on success - throws exception on failure
 Args    : A valid module name.

=cut

sub load_module {

	my ($self, $module_name) = @_;
	eval "require $module_name";
	if ($@) {
	  my $self_class = ref $self;
	  my $err_code = 'failed_to_load_module';
	  my $err_msg  = "$module_name in $self_class\n$@";
	  $self->throw($err_code, $err_msg);
	}
	return 1;
}

#-----------------------------------------------------------------------------

=head2 revcomp

 Title   : revcomp
 Usage   : $base->revcomp($sequence);
 Function: Get the reverse compliment of a nucleotide sequence
 Returns : The reverse complimented sequence
 Args    : A nucleotide sequence (ACGTRYMKSWHBVDNX).  Input sequence is case
           insensitive and case is maintained.

=cut

sub revcomp {

  my ($self, $sequence) = @_;

  my $revcomp_seq = reverse $sequence;
  $revcomp_seq =~ tr/acgtrymkswhbvdnxACGTRYMKSWHBVDNX/tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX/;
  return $revcomp_seq;
}

#-----------------------------------------------------------------------------

=head2 get_feature_bins

 Title   : get_feature_bins
 Usage   : $base->get_feature_bins($feature);
 Function: Get the genome bins for a given range.  GAL uses a feature binning
           scheme to speed range queries to a database.  This binning
           scheme is similar to (but not identical to) the scheme
           described in the UCSC Genome Browser publication
           (http://www.ncbi.nlm.nih.gov/pubmed/12045153) and on the
           UCSC Genome Browser wiki
           (http://genomewiki.ucsc.edu/index.php/Bin_indexing_system).
           Breifly, for a given genomic range this method will
           calculate the bins on the genome of sizes 512 MB, 64 MB, 8
           MB, 1 MB, and 128 KB that completely contains the feature.
           These bins are returned as a list sorted smallest bin to
           largest so that the first value returned is the smallest
           bin that completely contains the feature.
 Returns : An array of bins that the given range is compeltely contained
           within sorted smallest bin to largest so that the first
           value return is smallest possible bin in which the feature
           is fully contained.
 Args    : A hash reference with key values seqid, start, end (i.e. a
	   feature hash) or an object with methods seqid, start, end.

=cut

sub get_feature_bins {

    my ($self, $feature) = @_;

    my ($seqid, $start, $end);
    if (ref $feature eq 'HASH') {
      ($seqid, $start, $end) = @{$feature}{qw(seqid start end)};
    }
    elsif (blessed $feature && $feature->can('seqid') &&
	   $feature->can('start') && $feature->can('end')) {
      ($seqid, $start, $end) = ($feature->seqid, $feature->start,
				$feature->end)
    }
    else {
      my $data = ref $feature || $feature;
      $self->throw('invalid_arguments_to_get_feature_bins', $data);
    }
    my @feature_bins;
    my $count;
    my $single_bin;
    for my $bin_size (128_000, 1_000_000, 8_000_000, 64_000_000,
		      512_000_000) {
      $count++;
      my $start_bin = int($start/$bin_size);
      my $end_bin   = int($end/$bin_size);
      my @these_bins = map {$_ = join ':', ($seqid, $count, $_)} ($start_bin .. $end_bin);
	if (! $single_bin && scalar @these_bins == 1) {
	    $single_bin = shift @these_bins;
	}
	unshift @feature_bins, @these_bins;
    }
    unshift @feature_bins, $single_bin;
    return wantarray ? @feature_bins : \@feature_bins;
}

#-----------------------------------------------------------------------------

=head2 translate

 Title   : translate
 Usage   : $base->translate($sequence, $offset, $length);
 Function: Translate a nucleotide sequence to an amino acid sequence
 Returns : An amino acid sequence
 Args    : The sequence as a scalar, an integer offset from which to begin
           translation (default 0), and an integer length to translate
           (default is full length of input sequence).

=cut

sub translate {
  my ($self, $sequence, $offset, $length) = @_;

  my $genetic_code = $self->genetic_code;

  $offset ||= 0;
  $sequence =~ s/\s+|\d+//g;
  $length ||= length($sequence);

  my $polypeptide;
  for (my $i = (0 + $offset); $i < $length; $i += 3) {
    my $codon = uc substr($sequence, $i, 3);
    my $aa = $genetic_code->{$codon};
    $polypeptide .= $aa;
  }
  return $polypeptide;
}

#-----------------------------------------------------------------------------

=head2 genetic_code

 Title   : genetic_code
 Usage   : $base->genetic_code;
 Function: Returns a hash reference of the genetic code.  Currently only
           supports the standard genetic code without ambiguity codes.
 Returns : A hash reference of the genetic code
 Args    : None

=cut

sub genetic_code {
  my $self = shift;

  #TODO: Need to add ambiguity codes to codons.

  return {AAA => 'K',
	  AAC => 'N',
	  AAG => 'K',
	  AAT => 'N',
	  ACA => 'T',
	  ACC => 'T',
	  ACG => 'T',
	  ACT => 'T',
	  AGA => 'R',
	  AGC => 'S',
	  AGG => 'R',
	  AGT => 'S',
	  ATA => 'I',
	  ATC => 'I',
	  ATG => 'M',
	  ATT => 'I',
	  CAA => 'Q',
	  CAC => 'H',
	  CAG => 'Q',
	  CAT => 'H',
	  CCA => 'P',
	  CCC => 'P',
	  CCG => 'P',
	  CCT => 'P',
	  CGA => 'R',
	  CGC => 'R',
	  CGG => 'R',
	  CGT => 'R',
	  CTA => 'L',
	  CTC => 'L',
	  CTG => 'L',
	  CTT => 'L',
	  GAA => 'E',
	  GAC => 'D',
	  GAG => 'E',
	  GAT => 'D',
	  GCA => 'A',
	  GCC => 'A',
	  GCG => 'A',
	  GCT => 'A',
	  GGA => 'G',
	  GGC => 'G',
	  GGG => 'G',
	  GGT => 'G',
	  GTA => 'V',
	  GTC => 'V',
	  GTG => 'V',
	  GTT => 'V',
	  TAA => '*',
	  TAC => 'Y',
	  TAG => '*',
	  TAT => 'Y',
	  TCA => 'S',
	  TCC => 'S',
	  TCG => 'S',
	  TCT => 'S',
	  TGA => '*',
	  TGC => 'C',
	  TGG => 'W',
	  TGT => 'C',
	  TTA => 'L',
	  TTC => 'F',
	  TTG => 'L',
	  TTT => 'F',
	 };
}

#-----------------------------------------------------------------------------

=head2 amino_acid_data

 Title   : amino_acid_data
 Usage   : $base->amino_acid_data($aa, $value);
 Function: Returns data about an amino acid - either a specific value or a
           hash reference of all data for that amino acid.
 Returns : A string representing a data point about an given amino acid or a
           hash reference with all data points about the given amino acid.
 Args    : 1) An amino acid in either:
	      a) single-letter code
	      b) three-letter
	   2) Optionally a code for the value to return.  Any of:
	      a) name : The full name of the amino acid.
	      b) one_letter : The one-letter code for the amino acid.
	      c) three_letter : The three-letter code for the amino acid.
	      d) polarity : The polarity of the amino acid.
	      e) charge : The charge of the amino acid's side chain (at pH 7.4).
	      f) hydropathy : The hydropathy index of the amino acid.
	      g) weight : The molecular weight of the amino acid in g/mol.
	      h) size : Returns small or large.
	      i) h_bond : return 1 or 0 indicating if the amino acid can form hydrogen bonds
	      j) aromaticity : Returns aromatic, aliphatic or undef.

=cut

sub amino_acid_data {

  my ($self, $aa, $datum) = @_;

  my %aa321 = (Ala => 'A',
	       Arg => 'R',
	       Asn => 'N',
	       Asp => 'D',
	       Cys => 'C',
	       Glu => 'E',
	       Gln => 'Q',
	       Gly => 'G',
	       His => 'H',
	       Ile => 'I',
	       Leu => 'L',
	       Lys => 'K',
	       Met => 'M',
	       Phe => 'F',
	       Pro => 'P',
	       Ser => 'S',
	       Thr => 'T',
	       Trp => 'W',
	       Tyr => 'Y',
	       Val => 'V',
	       Sec => 'U',
	       Pyl => 'O',
	      );

  if (length($aa) == 3) {
    $aa = $aa321{$aa}
  }

  my %aa_data = (A => {name         => 'Alanine',
		       one_letter   => 'A',
		       three_letter => 'Ala',
		       polarity     => 'nonpolar',
		       charge       => 'neutral',
		       hydropathy   => 1.8,
		       weight       => 71.09,
		       size         => 'small',
		       h_bond       => 0,
		       aromaticity  => 'aliphatic',
		      },
		 R => {name         => 'Arginine',
		       one_letter   => 'R',
		       three_letter => 'Arg',
		       polarity     => 'polar',
		       charge       => 'positive',
		       hydropathy   => -4.5,
		       weight       => 156.19,
		       size         => 'large',
		       h_bond       => 1,
		       aromaticity  => undef,
		      },
		 N => {name         => 'Asparagine',
		       one_letter   => 'N',
		       three_letter => 'Asn',
		       polarity     => 'polar',
		       charge       => 'neutral',
		       hydropathy   => -3.5,
		       weight       => 114.11,
		       size         => 'small',
		       h_bond       => 1,
		       aromaticity  => undef,
		      },
		 D => {name         => 'Aspartic acid',
		       one_letter   => 'D',
		       three_letter => 'Asp',
		       polarity     => 'polar',
		       charge       => 'negative',
		       hydropathy   => -3.5,
		       weight       => 115.09,
		       size         => 'small',
		       h_bond       => 1,
		       aromaticity  => undef,
		      },
		 C => {name         => 'Cysteine',
		       one_letter   => 'C',
		       three_letter => 'Cys',
		       polarity     => 'polar',
		       charge       => 'neutral',
		       hydropathy   => 2.5,
		       weight       => 103.15,
		       size         => 'small',
		       h_bond       => 1,
		       aromaticity  => undef,
		      },
		 E => {name         => 'Glutamic acid',
		       one_letter   => 'E',
		       three_letter => 'Glu',
		       polarity     => 'polar',
		       charge       => 'negative',
		       hydropath    => -3.5,
		       weight       => 129.12,
		       size         => 'large',
		       h_bond       => 1,
		       aromaticity  => undef,
		     },
		 Q => {name         => 'Glutamine',
		       one_letter   => 'Q',
		       three_letter => 'Gln',
		       polarity     => 'polar',
		       charge       => 'neutral',
		       hydropathy   => -3.5,
		       weight       => 128.14,
		       size         => 'large',
		       h_bond       => 1,
		       aromaticity  => undef,
		      },
		 G => {name         => 'Glycine',
		       one_letter   => 'G',
		       three_letter => 'Gly',
		       polarity     => 'nonpolar',
		       charge       => 'neutral',
		       hydropathy   => -0.4,
		       weight       => 57.05,
		       size         => 'small',
		       h_bond       => 0,
		       aromaticity  => 'aliphatic',
		      },
		 H => {name         => 'Histidine',
		       one_letter   => 'H',
		       three_letter => 'His',
		       polarity     => 'polar',
		       charge       => 'positive',
		       hydropathy   => -3.2,
		       weight       => 137.14,
		       size         => 'large',
		       h_bond       => 1,
		       aromaticity  => 'aromatic',
		      },
		 I => {name         => 'Isoleucine',
		       one_letter   => 'I',
		       three_letter => 'Ile',
		       polarity     => 'nonpolar',
		       charge       => 'neutral',
		       hydropathy   => 4.5,
		       weight       => 113.16,
		       size         => 'large',
		       h_bond       => 0,
		       aromaticity  => 'aliphatic',
		      },
		 L => {name         => 'Leucine',
		       one_letter   => 'L',
		       three_letter => 'Leu',
		       polarity     => 'nonpolar',
		       charge       => 'neutral',
		       hydropathy   => 3.8,
		       weight       => 113.16,
		       size         => 'large',
		       h_bond       => 0,
		       aromaticity  => 'aliphatic',
		      },
		 K => {name         => 'Lysine',
		       one_letter   => 'K',
		       three_letter => 'Lys',
		       polarity     => 'polar',
		       charge       => 'positive',
		       hydropathy   => -3.9,
		       weight       => 128.17,
		       size         => 'large',
		       h_bond       => 1,
		       aromaticity  => undef,
		      },
		 M => {name         => 'Methionine',
		       one_letter   => 'M',
		       three_letter => 'Met',
		       polarity     => 'nonpolar',
		       charge       => 'neutral',
		       hydropathy   => 1.9,
		       weight       => 131.19,
		       size         => 'large',
		       h_bond       => 0,
		       aromaticity  => undef,
		      },
		 F => {name         => 'Phenylalanine',
		       one_letter   => 'F',
		       three_letter => 'Phe',
		       polarity     => 'nonpolar',
		       charge       => 'neutral',
		       hydropathy   => 2.8,
		       weight       => 147.18,
		       size         => 'large',
		       h_bond       => 0,
		       aromaticity  => 'aromatic',
		      },
		 P => {name         => 'Proline',
		       one_letter   => 'P',
		       three_letter => 'Pro',
		       polarity     => 'nonpolar',
		       charge       => 'neutral',
		       hydropathy   => -1.6,
		       weight       => 97.12,
		       size         => 'small',
		       h_bond       => 0,
		       aromaticity  => 'aliphatic',
		      },
		 S => {name         => 'Serine',
		       one_letter   => 'S',
		       three_letter => 'Ser',
		       polarity     => 'polar',
		       charge       => 'neutral',
		       hydropathy   => -0.8,
		       weight       => 87.08,
		       size         => 'small',
		       h_bond       => 1,
		       aromaticity  => undef,
		      },
		 T => {name         => 'Threonine',
		       one_letter   => 'T',
		       three_letter => 'Thr',
		       polarity     => 'polar',
		       charge       => 'neutral',
		       hydropathy   => -0.7,
		       size         => 'small',
		       h_bond       => 1,
		       aromaticity  => undef,
		      },
		 W => {name         => 'Tryptophan',
		       one_letter   => 'W',
		       three_letter => 'Trp',
		       polarity     => 'nonpolar',
		       charge       => 'neutral',
		       hydropathy   => -0.9,
		       weight       => 186.21,
		       size         => 'large',
		       h_bond       => 1,
		       aromaticity  => 'aromatic',
		      },
		 Y => {name         => 'Tyrosine',
		       one_letter   => 'Y',
		       three_letter => 'Tyr',
		       polarity     => 'polar',
		       charge       => 'neutral',
		       hydropathy   => -1.3,
		       weight       => 163.18,
		       size         => 'large',
		       h_bond       => 1,
		       aromaticity  => 'aromatic',
		      },
		 V => {name         => 'Valine',
		       one_letter   => 'V',
		       three_letter => 'Val',
		       polarity     => 'nonpolar',
		       charge       => 'neutral',
		       hydropathy   => 4.2,
		       weight       => 99.14,
		       size         => 'small',
		       h_bond       => 0,
		       aromaticity  => 'aliphatic',
		      },
		 U => {name         => 'Selenocysteine',
		       one_letter   => 'U',
		       three_letter => 'Sec',
		       polarity     => undef,
		       charge       => undef,
		       hydropathy   => undef,
		       weight       => undef,
		       size         => undef,
		       h_bond       => undef,
		       aromaticity  => undef,
		      },
		 O => {name         => 'Pyrrolysine',
		       one_letter   => 'O',
		       three_letter => 'Pyl',
		       polarity     => undef,
		       charge       => undef,
		       hydropathy   => undef,
		       weight       => undef,
		       size         => undef,
		       h_bond       => undef,
		       aromaticity  => undef,
		      }
		);


  if (defined $datum) {
    if (exists $aa_data{$aa}{$datum}) {
      return $aa_data{$aa}{$datum};
    }
    else {
      $self->throw('invalid_aa_datum_code',
		   "$datum passed to GAL::Base::amino_acid_data");
    }
  }
  else {
      return $aa_data{$aa};
  }
}

#-----------------------------------------------------------------------------

=head2 time_stamp

 Title   : time_stamp
 Usage   : $base->time_stamp;
 Function: Returns a YYYYMMDD time_stamp
 Returns : A YYYYMMDD time_stamp
 Args    : None

=cut

sub time_stamp {

  my $self = shift;

  my ($sec,$min,$hour,$mday,$mon,$year,$wday,
      $yday,$isdst) = localtime(time);
  my $time_stamp = sprintf("%02d%02d%02d", $year + 1900,
			   $mon + 1, $mday);
  return $time_stamp;
}

#-----------------------------------------------------------------------------

=head2 random_string

 Title   : random_string
 Usage   : $base->random_string(8);
 Function: Returns a random alphanumeric string
 Returns : A random alphanumeric string of a given length.  Default length is
           8 charachters.
 Args    : The length of the string to be returned [8]

=cut

sub random_string {
  my ($self, $length) = @_;
  $length ||= 8;
  my $random_string = join "", map { unpack "H*", chr(rand(256)) } (1 .. $length);
  return substr($random_string, 0, $length);
}

#-----------------------------------------------------------------------------

=head2 float_lt

 Title   : float_lt
 Usage   : $base->float_lt(0.0000123, 0.0000124, 7);
 Function: Return true if the first number given is less than (<) the
	   second number at a given level of accuracy.  Default accuracy
           length is 6 decimal places.
 Returns : 1 if the first number is less than the second, otherwise 0.
 Args    : The two values to compare and optionally a integer value for
	   the accuracy.  Accuracy defaults to 6 decimal places.

=cut

sub float_lt {
	my ($self, $A, $B, $accuracy) = @_;

	$accuracy ||= 6;
	$A = sprintf("%.${accuracy}f", $A);
	$B = sprintf("%.${accuracy}f", $B);
	$A =~ s/\.//;
	$B =~ s/\.//;
	return $A < $B;
}

#-----------------------------------------------------------------------------

=head2 float_le

 Title   : float_le
 Usage   : $base->float_le(0.0000123, 0.0000124, 7);
 Function: Return true if the first number given is less than or equal to
	   (<=) the second number at a given level of accuracy.  Default
           accuracy length is 6 decimal places.
 Returns : 1 if the first number is less than or equal to the second,
	   otherwise 0.
 Args    : The two values to compare and optionally a integer value for
	   the accuracy.  Accuracy defaults to 6 decimal places.

=cut

sub float_le {
	my ($self, $A, $B, $accuracy) = @_;

	$accuracy ||= 6;
	$A = sprintf("%.${accuracy}f", $A);
	$B = sprintf("%.${accuracy}f", $B);
	$A =~ s/\.//;
	$B =~ s/\.//;
	return $A <= $B;
}

#-----------------------------------------------------------------------------

=head2 float_gt

 Title   : float_gt
 Usage   : $base->float_gt(0.0000123, 0.0000124, 7);
 Function: Return true if the first number given is greater than (>) the
	   second number at a given level of accuracy.  Default accuracy
           length is 6 decimal places.
 Returns : 1 if the first number is greater than the second, otherwise 0
 Args    : The two values to compare and optionally a integer value for
	   the accuracy.  Accuracy defaults to 6 decimal places.

=cut

sub float_gt {
	my ($self, $A, $B, $accuracy) = @_;

	$accuracy ||= 6;
	$A = sprintf("%.${accuracy}f", $A);
	$B = sprintf("%.${accuracy}f", $B);
	$A =~ s/\.//;
	$B =~ s/\.//;
	return $A > $B;
}

#-----------------------------------------------------------------------------

=head2 float_ge

 Title   : float_ge
 Usage   : $base->float_ge(0.0000123, 0.0000124, 7);
 Function: Return true if the first number given is greater than or equal
	   to (>=) the second number at a given level of accuracy.  Default
           accuracy length is 6 decimal places.
 Returns : 1 if the first number is greater than or equal to the second,
	   otherwise 0.
 Args    : The two values to compare and optionally a integer value for
	   the accuracy.  Accuracy defaults to 6 decimal places.

=cut

sub float_ge {
	my ($self, $A, $B, $accuracy) = @_;

	$accuracy ||= 6;
	$A = sprintf("%.${accuracy}f", $A);
	$B = sprintf("%.${accuracy}f", $B);
	$A =~ s/\.//;
	$B =~ s/\.//;
	return $A >= $B;
}

#-----------------------------------------------------------------------------

=head2 open_file

 Title   : open_file
 Usage   : $base->open_file($file);
 Function: Open a given file for reading and return a filehandle.
 Returns : A filehandle
 Args    : A file path/name

=cut

sub open_file {

  my ($self, $file) = @_;

  if (! defined $file) {
    $self->throw('file_does_not_exist', $file);
  }
  if (! -e $file) {
    $self->throw('file_does_not_exist', $file);
  }
  if (! -r $file) {
    $self->throw('cant_read_file', $file);
  }
  open(my $FH, '<', $file) || $self->throw('cant_open_file_for_reading', $file);

  return $FH;
}

#-----------------------------------------------------------------------------

##############################################################################
## Consolidate three copies of this method 1) Here 2) Feature/gene.pm and
## 3) Base.pm
##############################################################################

=head2 get_transcript_types

 Title   : get_transcript_types
 Usage   : $transcript_types = $self->get_transcript_types
 Function: Get a list of all features that are types of transcripts.
 Returns : An array or ref of feature names for transcript and it's children.
 Args    : None

=cut

sub get_transcript_types {

    my $self = shift;

    my @transcripts = qw(EST
			 RNase_MRP_RNA
			 RNase_P_RNA
			 SRP_RNA
			 SRP_RNA_primary_transcript
			 Y_RNA
			 aberrant_processed_transcript
			 alternatively_spliced_transcript
			 antisense_RNA
			 antisense_primary_transcript
			 capped_mRNA
			 capped_primary_transcript
			 class_II_RNA
			 class_I_RNA
			 consensus_mRNA
			 dicistronic_mRNA
			 dicistronic_primary_transcript
			 dicistronic_transcript
			 edited_mRNA
			 edited_transcript
			 edited_transcript_by_A_to_I_substitution
			 enhancerRNA
			 enzymatic_RNA
			 exemplar_mRNA
			 guide_RNA
			 lnc_RNA
			 mRNA
			 mRNA_recoded_by_codon_redefinition
			 mRNA_recoded_by_translational_bypass
			 mRNA_region
			 mRNA_with_frameshift
			 mRNA_with_minus_1_frameshift
			 mRNA_with_minus_2_frameshift
			 mRNA_with_plus_1_frameshift
			 mRNA_with_plus_2_frameshift
			 mature_transcript
			 mature_transcript_region
			 miRNA_primary_transcript
			 mini_exon_donor_RNA
			 monocistronic_mRNA
			 monocistronic_primary_transcript
			 monocistronic_transcript
			 ncRNA
			 nc_primary_transcript
			 piRNA
			 polyadenylated_mRNA
			 polycistronic_mRNA
			 polycistronic_primary_transcript
			 polycistronic_transcript
			 pre_edited_mRNA
			 primary_transcript
			 primary_transcript_region
			 processed_transcript
			 protein_coding_primary_transcript
			 pseudogenic_transcript
			 rRNA
			 rRNA_cleavage_RNA
			 rRNA_primary_transcript
			 rasiRNA
			 recoded_mRNA
			 regional_centromere_outer_repeat_transcript
			 riboswitch
			 ribozyme
			 scRNA
			 scRNA_primary_transcript
			 siRNA
			 small_regulatory_ncRNA
			 snRNA
			 snRNA_primary_transcript
			 snoRNA
			 snoRNA_primary_transcript
			 spliced_leader_RNA
			 stRNA
			 tRNA
			 tRNA_primary_transcript
			 tasiRNA
			 tasiRNA_primary_transcript
			 telomerase_RNA
			 tmRNA_primary_transcript
			 trans_spliced_mRNA
			 trans_spliced_transcript
			 transcript
			 transcript_bound_by_nucleic_acid
			 transcript_bound_by_protein
			 transcript_region
			 transcript_with_translational_frameshift
			 vault_RNA
  );
  return wantarray ? @transcripts : \@transcripts;
}

#-----------------------------------------------------------------------------

#=head2  validate_feature
#
# Title   : validate_feature
# Usage   : $self->validate_feature($feature_hash)
# Function: Validates a feature hash.
# Returns : NA
# Args    : A feature hash reference in the form returned by next_feature_hash
#
#=cut
#
#sub validate_feature {
#  my ($self, $feature) = @_;
#
#  my ($feature_id, $seqid, $source, $type, $start, $end, $score,
#      $strand, $phase, $attributes) =
#	@{$feature}{qw(feature_id seqid source type start end
#		       score strand phase attributes)};
#
#  my $feature_text = $self->to_gff3($feature);
#
#  # Validate seqid
#  unless ($self->is_valid_seqid) {
#    my $error_code = 'invalid_characters_in_seqid_column';
#    $self->throw($error_code, "($seqid) $feature_text");
#  }
#  # Validate source - No validation
#  # Validate type
#  unless ($self->is_valid_sequence_feature($type)) {
#    my $error_code = 'invalid_type_value';
#    $self->warn($error_code, "($type) $feature_text");
#  }
#  # Validate start
#  unless ($self->is_integer($start)) {
#    my $error_code = 'invalid_start_value';
#    $self->throw($error_code, "($start) $feature_text");
#  }
#  # Validate end
#  unless ($self->is_integer($end)) {
#    my $error_code = 'invalid_end_value';
#    $self->throw($error_code, "($end) $feature_text");
#  }
#  # Validate end > start
#  unless ($start <= $end) {
#    my $error_code = 'invalid_start_gt_end_value';
#    $self->throw($error_code, "($start > $end) $feature_text");
#  }
#  # Validate score
#  unless ($self->is_real_number($score)) {
#    my $error_code = 'invalid_score_value';
#    $self->warn($error_code, "($score) $feature_text");
#  }
#  # Validate strand
#  unless ($strand =~ /^[.+-]/) {
#    my $error_code = 'invalid_strand_value';
#    $self->throw($error_code, "($strand) $feature_text");
#  }
#  # Validate phase
#  unless ($phase =~ /^[012.]/) {
#    my $error_code = 'invalid_phase_value';
#    $self->warn($error_code, "($phase) $feature_text");
#  }
#  # Validate attributes
#  $self->validate_attributes($attributes);
#
#  # feature_id
#  $feature->{feature_id} ||= $feature->{attributes}{ID}[0];
#  if ($feature->{feature_id} ne $feature->{attributes}{ID}[0]) {
#    my $fid = $feature->{feature_id};
#    my $aid = $feature->{attributes}{ID}[0];
#    my $error_code = 'feature_id_ne_attribute_ID';
#    $self->throw($error_code, "($fid ne $aid) $feature_text");
#  }
#}
#
##-----------------------------------------------------------------------------
#
#=head2  is_valid_gff3_attribute_key
#
# Title   : is_valid_gff3_attribute_key
# Usage   : $self->is_valid_gff3_attribute_key($key)
# Function: Validates a ...
# Returns : 1 if valid undef otherwise.
# Args    : An attribute key
#
#=cut
#
#sub validate_attributes {
#
#  $self->{valid_gff3_attributes_keys} ||= {};
#
#  return 1 if exists
#
#}
#
##-----------------------------------------------------------------------------
#
#=head2  is_valid_gvf_attribute_key
#
# Title   : is_valid_gvf_attribute_key
# Usage   : $self->is_valid_gvf_attribute_key($key)
# Function: Validates a ...
# Returns : 1 if valid undef otherwise.
# Args    : An attribute key
#
#=cut
#
#sub validate_attributes {
#
#
#}
#
##-----------------------------------------------------------------------------
#
#=head2  validate_attributes
#
# Title   : validate_attributes
# Usage   : $self->validate_attributes($feature_hash)
# Function: Validates a ...
# Returns : 1 if valid undef otherwise.
# Args    : A string of text.
#
#=cut
#
#sub validate_attributes {
#
#  my $attributes = shift;
#
#  if (! defined $attributes) {
#    my $error_code = 'missing_attributes';
#    $self->warn($error_code, $feature_text);
#  }
#  elsif (ref $attributes ne 'HASH') {
#    my $error_code = 'attributes_must_be_hash';
#    $self->warn($error_code, $feature_text);
#  }
#  else {
#    for my $key (keys %{$attributes}) {
#      unless ($self->is_valid_gff3_attribute_key($key) ||
#	      $self->is_valid_gvf_attribute_key($key)) {
#	my $error_code = 'invalid_attribute_key';
#	$self->throw($error_code, "($key) $feature_text");
#      }
#      unless (defined $attributes->{$key}) {
#	my $error_code = 'missing_attribute_value';
#	$self->throw($error_code, "($key) $feature_text");
#      }
#      unless (ref $attributes->{$key} eq 'ARRAY') {
#	my $error_code = 'attribute_value_must_be_array';
#	$self->throw($error_code, "($key) $feature_text");
#      }
#      my @values = @{$attributes->{$key}};
#      # Validate GFF3 Attributes
#      # ID
#      if ($key eq 'ID') {
#	if (scalar @values > 1) {
#	  my $error_code = 'invalid_ID_attribute_multiple_values';
#	  my $all_ids = join ' ', @values;
#	  $self->throw($error_code, "($all_ids) $feature_text");
#	}
#      }
#      # Name
#      elsif ($key eq 'Name') {
#	# No validation on Name values
#      }
#      # Alias
#      elsif ($key eq 'Alias') {
#	# No validation on Alias values
#      }
#      # Parent
#      elsif ($key eq 'Parent') {
#	# No validation on Alias values
#      }
#      # Target
#      elsif ($key eq 'Target') {
#      }
#      # Gap
#      elsif ($key eq 'Gap') {
#      }
#      # Derives_from
#      elsif ($key eq 'Derives_from') {
#      }
#      # Note
#      elsif ($key eq 'Note') {
#      }
#      # Dbxref
#      elsif ($key eq 'Dbxref') {
#      }
#      # Ontology_term
#      elsif ($key eq 'Ontology_term') {
#	# TODO: Validate Ontology_term
#      }
#      # Is_circular
#      elsif ($key eq 'Is_circular') {
#      }
#      # Validate GVF attributes
#      # Variant_seq
#      elsif ($key eq 'Variant_seq') {
#      }
#      # Reference_seq
#      elsif ($key eq 'Reference_seq') {
#      }
#      # Variant_reads
#      elsif ($key eq 'Variant_reads') {
#      }
#      # Total_reads
#      elsif ($key eq 'Total_reads') {
#      }
#      # Zygosity
#      elsif ($key eq 'Zygosity') {
#      }
#      # Variant_freq
#      elsif ($key eq 'Variant_freq') {
#      }
#      # Variant_effect
#      elsif ($key eq 'Variant_effect') {
#      }
#      # Start_range
#      elsif ($key eq 'Start_range') {
#      }
#      # End_range
#      elsif ($key eq 'End_range') {
#      }
#      # Phased
#      elsif ($key eq 'Phased') {
#      }
#      # Genotype
#      elsif ($key eq 'Genotype') {
#      }
#      # Individual
#      elsif ($key eq 'Individual') {
#      }
#      # Variant_copy_number
#      elsif ($key eq 'Variant_copy_number') {
#      }
#      # Reference_copy_number
#      elsif ($key eq 'Reference_copy_number') {
#      }
#      # Variant_codon
#      elsif ($key eq 'Variant_codon') {
#      }
#      # Reference_codon
#      elsif ($key eq 'Reference_codon') {
#      }
#      # Variant_aa
#      elsif ($key eq 'Variant_aa') {
#      }
#      # Reference_aa
#      elsif ($key eq 'Reference_aa') {
#      }
#      # Breakpoint_detail
#      elsif ($key eq 'Breakpoint_detail') {
#      }
#      # Sequence_context
#      elsif ($key eq 'Sequence_context') {
#      }
#      # All others
#      else {
#	if ($key =~ /^[A-Z]/) {
#	  my $error_code = 'non_standard_attribute_must_begin_with_lowercase';
#	  $self->throw($error_code, "($key) $feature_text");
#	}
#      }
#    }
#
##-----------------------------------------------------------------------------
#
#=head2  is_integer
#
# Title   : is_integer
# Usage   : $self->is_integer($feature_hash)
# Function: Validates a ...
# Returns : 1 if valid undef otherwise.
# Args    : A string of text.
#
#=cut
#
#sub is_integer {
#
#
#}
#
##-----------------------------------------------------------------------------
#
#=head2  is_real_number
#
# Title   : is_real_number
# Usage   : $self->is_real_number($feature_hash)
# Function: Validates a ...
# Returns : 1 if valid undef otherwise.
# Args    : A string of text.
#
#=cut
#
#sub is_real_number {
#
#
#}
#
##-----------------------------------------------------------------------------
#
#=head2  is_valid_aa_sequence
#
# Title   : is_valid_aa_sequence
# Usage   : $self->is_valid_aa_sequence($feature_hash)
# Function: Validates a ...
# Returns : 1 if valid undef otherwise.
# Args    : A string of text.
#
#=cut
#
#sub is_valid_aa_sequence {
#
#
#}
#
##-----------------------------------------------------------------------------
#
#=head2  is_valid_iupac_nt
#
# Title   : is_valid_iupac_nt
# Usage   : $self->is_valid_iupac_nt($feature_hash)
# Function: Validates a ...
# Returns : 1 if valid undef otherwise.
# Args    : A string of text.
#
#=cut
#
#sub is_valid_iupac_nt {
#
#
#}
#
##-----------------------------------------------------------------------------
#
#=head2  is_valid_seqid
#
# Title   : is_valid_seqid
# Usage   : $self->is_valid_seqid($feature_hash)
# Function: Validates a ...
# Returns : 1 if valid undef otherwise.
# Args    : A string of text.
#
#=cut
#
#sub is_valid_seqid {
#
#
#}
#
##-----------------------------------------------------------------------------
#
#=head2  is_valid_type
#
# Title   : is_valid_type
# Usage   : $self->is_valid_type($feature_hash)
# Function: Validates a ...
# Returns : 1 if valid undef otherwise.
# Args    : A string of text.
#
#=cut
#
#sub is_valid_type {
#
#
#}
#
##-----------------------------------------------------------------------------
#
#=head2  is_valid_start
#
# Title   : is_valid_start
# Usage   : $self->is_valid_start($feature_hash)
# Function: Validates a ...
# Returns : 1 if valid undef otherwise.
# Args    : A string of text.
#
#=cut
#
#sub is_valid_start {
#
#
#}
#
##-----------------------------------------------------------------------------
#
#=head2  is_valid_end
#
# Title   : is_valid_end
# Usage   : $self->is_valid_end($feature_hash)
# Function: Validates a ...
# Returns : 1 if valid undef otherwise.
# Args    : A string of text.
#
#=cut
#
#sub is_valid_end {
#
#
#}
#
##-----------------------------------------------------------------------------
#
#=head2  is_valid_strand
#
# Title   : is_valid_strand
# Usage   : $self->is_valid_strand($feature_hash)
# Function: Validates a ...
# Returns : 1 if valid undef otherwise.
# Args    : A string of text.
#
#=cut
#
#sub is_valid_strand {
#
#
#}
#
##-----------------------------------------------------------------------------
#
#=head2  is_valid_phase
#
# Title   : is_valid_phase
# Usage   : $self->is_valid_phase($feature_hash)
# Function: Validates a ...
# Returns : 1 if valid undef otherwise.
# Args    : A string of text.
#
#=cut
#
#sub is_valid_phase {
#
#
#}
#
##-----------------------------------------------------------------------------
#
#=head2  is_valid_id_attribute
#
# Title   : is_valid_id_attribute
# Usage   : $self->is_valid_id_attribute($feature_hash)
# Function: Validates a id_attribute.
# Returns : 1 if valid undef otherwise.
# Args    : A string of text.
#
#=cut
#
#sub is_valid_id_attribute {
#
#  my $value = shift;
##-----------------------------------------------------------------------------
#
#=head2  is_valid
#
# Title   : is_valid
# Usage   : $self->is_valid($feature_hash)
# Function: Validates a ...
# Returns : 1 if valid undef otherwise.
# Args    : A string of text.
#
#=cut
#
#sub is_valid {
#
#  my $values = shift;
#  if (scalar @{$values} > 1) {
#    my $error_code = 'invalid_ID_attribute_multiple_values';
#    my $all_ids = join ' ', @values;
#    $self->throw($error_code, "($all_ids) $feature_text");
#  }
#  return 1;
#}
#
##-----------------------------------------------------------------------------
#
#=head2  is_valid_name_attribute
#
# Title   : is_valid_name_attribute
# Usage   : $self->is_valid_name_attribute($feature_hash)
# Function: Validates a name_attribute.
# Returns : 1 if valid undef otherwise.
# Args    : A string of text.
#
#=cut
#
#sub is_valid_name_attribute {
#
#  my $values = shift;
#  # No validation of Name attribute
#  return 1;
#
#}
#
##-----------------------------------------------------------------------------
#
#=head2  is_valid_alias_attribute
#
# Title   : is_valid_alias_attribute
# Usage   : $self->is_valid_alias_attribute($feature_hash)
# Function: Validates a alias_attribute.
# Returns : 1 if valid undef otherwise.
# Args    : A string of text.
#
#=cut
#
#sub is_valid_alias_attribute {
#
#  my $values = shift;
#  # No validation of Alias attribute
#  return 1;
#
#}
#
##-----------------------------------------------------------------------------
#
#=head2  is_valid_parent_attribute
#
# Title   : is_valid_parent_attribute
# Usage   : $self->is_valid_parent_attribute($feature_hash)
# Function: Validates a parent_attribute.
# Returns : 1 if valid undef otherwise.
# Args    : A string of text.
#
#=cut
#
#sub is_valid_parent_attribute {
#
#  my $values = shift;
#  # No validation of Parent attribute
#  return 1;
#
#}
#
##-----------------------------------------------------------------------------
#
#=head2  is_valid_target_attribute
#
# Title   : is_valid_target_attribute
# Usage   : $self->is_valid_target_attribute($feature_hash)
# Function: Validates a target_attribute.
# Returns : 1 if valid undef otherwise.
# Args    : A string of text.
#
#=cut
#
#sub is_valid_target_attribute {
#
#  my $values = shift;
#
#  if (scalar @{$values} > 1) {
#    my $error_code = 'invalid_Target_attribute_multiple_values';
#    my $all_values = join ' ', @{$values};
#    $self->throw($error_code, "($all_values) $feature_text");
#  }
#  my $value = $values->[0];
#  if ($value !~ /\S+\s+\d+\s+\d+\s*[\-+]*/) {
#    my $error_code = 'invalid_Target_attribute_value';
#    $self->throw($error_code, "($value) $feature_text");
#  }
#  return 1;
#}
#
##-----------------------------------------------------------------------------
#
#=head2  is_valid_gap_attribute
#
# Title   : is_valid_gap_attribute
# Usage   : $self->is_valid_gap_attribute($feature_hash)
# Function: Validates a gap_attribute.
# Returns : 1 if valid undef otherwise.
# Args    : A string of text.
#
#=cut
#
#sub is_valid_gap_attribute {
#
#  my $values = shift;
#  if (scalar @values > 1) {
#    my $error_code = 'invalid_Gap_attribute_multiple_values';
#    my $all_values = join ' ', @values;
#    $self->throw($error_code, "($all_values) $feature_text");
#  }
#  my $value = $values[0];
#  # TODO: Validate the CIGAR format
#  return 1;
#}
#
##-----------------------------------------------------------------------------
#
#=head2  is_valid_derives_from_attribute
#
# Title   : is_valid_derives_from_attribute
# Usage   : $self->is_valid_derives_from_attribute($feature_hash)
# Function: Validates a derives_from_attribute.
# Returns : 1 if valid undef otherwise.
# Args    : A string of text.
#
#=cut
#
#sub is_valid_derives_from_attribute {
#
#  my $values = shift;
#
#  if (scalar @values > 1) {
#    my $error_code = 'invalid_Derives_from_attribute_multiple_values';
#    my $all_values = join ' ', @values;
#    $self->throw($error_code, "($all_values) $feature_text");
#  }
#  # TODO: Validate the Derives_from relationships
#  return 1;
#}
#
##-----------------------------------------------------------------------------
#
#=head2  is_valid_note_attribute
#
# Title   : is_valid_note_attribute
# Usage   : $self->is_valid_note_attribute($feature_hash)
# Function: Validates a note_attribute.
# Returns : 1 if valid undef otherwise.
# Args    : A string of text.
#
#=cut
#
#sub is_valid_note_attribute {
#
#  my $values = shift;
#  # No validation of Note attribute
#  return 1;
#}
#
##-----------------------------------------------------------------------------
#
#=head2  is_valid_dbxref_attribute
#
# Title   : is_valid_dbxref_attribute
# Usage   : $self->is_valid_dbxref_attribute($feature_hash)
# Function: Validates a dbxref_attribute.
# Returns : 1 if valid undef otherwise.
# Args    : A string of text.
#
#=cut
#
#sub is_valid_dbxref_attribute {
#
#  my $values = shift;
#  return 1;
#
#}
#
##-----------------------------------------------------------------------------
#
#=head2  is_valid_ontology_term_attribute
#
# Title   : is_valid_ontology_term_attribute
# Usage   : $self->is_valid_ontology_term_attribute($feature_hash)
# Function: Validates a ontology_term_attribute.
# Returns : 1 if valid undef otherwise.
# Args    : A string of text.
#
#=cut
#
#sub is_valid_ontology_term_attribute {
#
#  my $values = shift;
#  if (scalar @values > 1) {
#    my $error_code = 'invalid_Is_circular_attribute_multiple_values';
#    my $all_values = join ' ', @values;
#    $self->throw($error_code, "($all_values) $feature_text");
#  }
#  my $value = $values[0];
#  return 1;
#
#}
#
##-----------------------------------------------------------------------------
#
#=head2  is_valid_is_circular_attribute
#
# Title   : is_valid_is_circular_attribute
# Usage   : $self->is_valid_is_circular_attribute($feature_hash)
# Function: Validates a is_circular_attribute.
# Returns : 1 if valid undef otherwise.
# Args    : A string of text.
#
#=cut
#
#sub is_valid_is_circular_attribute {
#
#  my $values = shift;
#  # TODO: Validate Is_circular
#  return 1;
#
#}
#
##-----------------------------------------------------------------------------
#
#=head2  is_valid_variant_seq_attribute
#
# Title   : is_valid_variant_seq_attribute
# Usage   : $self->is_valid_variant_seq_attribute($feature_hash)
# Function: Validates a variant_seq_attribute.
# Returns : 1 if valid undef otherwise.
# Args    : A string of text.
#
#=cut
#
#sub is_valid_variant_seq_attribute {
#
#  my $values = shift;
#  for my $value (@values) {
#    # TODO: Validate the length of the Variant_seq
#    # Compile regex once above
#    # my $Variant_seq_regex = qr/[A-DGHKMNR-WY]+| # Any valid IUPAC Nucleotide
#    #                          ~\d*|            # A ~ optionally followed by an integer
#    #                          [.\-!^\*]        # Any of [.-!^*]
#    #                         /x;
#    if (! $self->is_valid_variant_seq_term($value)) {
#      if ($type eq 'SNV') {
#	# A single nt or .!^*
#	if ($value !~ /^([A-DGHKMNR-WY]|[.!^\*])$/i) {
#	  my $error_code = 'invalid_Variant_seq_attribute_value_for_SNV';
#	  $self->warn($error_code, "($value) $feature_text");
#	}
#      }
#      # A nt string or .!^*
#      elsif ($type eq 'insertion') {
#	if ($value !~ /^([A-DGHKMNR-WY]+|~\d*|[.!^\*])$/i) {
#	  my $error_code = 'invalid_Variant_seq_attribute_value_for_insertion';
#	  $self->warn($error_code, "($value) $feature_text");
#	}
#      }
#      # Any of .!^-
#      elsif ($type eq 'deletion') {
#	if ($value !~ /^[\.\^!-]$/i) {
#	  my $error_code = 'invalid_Variant_seq_attribute_value_for_deletion';
#	  $self->warn($error_code, "($value) $feature_text");
#	}
#      }
#      # A nt string or .!^*-
#      elsif ($type eq 'indel') {
#	if ($value !~ /^([A-DGHKMNR-WY]+|~\d*|[.\-!^])$/i) {
#	  my $error_code = 'invalid_Variant_seq_attribute_value_for_indel';
#	  $self->warn($error_code, "($value) $feature_text");
#	}
#      }
#    }
#    else {
#      my $error_code = 'invalid_Variant_seq_attribute_value';
#      $self->warn($error_code, "($value) $feature_text");
#    }
#  }
#  return 1;
#}
#
##-----------------------------------------------------------------------------
#
#=head2  is_valid_reference_seq_attribute
#
# Title   : is_valid_reference_seq_attribute
# Usage   : $self->is_valid_reference_seq_attribute($feature_hash)
# Function: Validates a reference_seq_attribute.
# Returns : 1 if valid undef otherwise.
# Args    : A string of text.
#
#=cut
#
#sub is_valid_reference_seq_attribute {
#
#  my $values = shift;
#  if (scalar @values > 1) {
#    my $error_code = 'invalid_Reference_seq_attribute_multiple_values';
#    my $all_values = join ', ', @values;
#    $self->warn($error_code, "($all_values) $feature_text");
#  }
#  my $value = $values[0];
#  if (! $self->is_valid_reference_seq($value)) {
#    if ($type eq 'SNV') {
#      if ($value !~ /^([A-DGHKMNR-WY]+|\.)$/) {
#	my $error_code = 'invalid_Reference_seq_attribute_value_for_SNV';
#	$self->warn($error_code, "($value) $feature_text");
#      }
#    }
#    elsif ($type eq 'insertion') {
#      if ($value ne '-') {
#	my $error_code = 'invalid_Reference_seq_attribute_value_insertion';
#	$self->warn($error_code, "($value) $feature_text");
#      }
#    }
#    elsif ($type eq 'deletion') {
#      if ($value !~ /^([A-DGHKMNR-WY]+|[\.~])$/) {
#	my $error_code = 'invalid_Reference_seq_attribute_value_deletion';
#	$self->warn($error_code, "($value) $feature_text");
#      }
#    }
#    elsif ($type eq 'indel') {
#      if ($value !~ /^([A-DGHKMNR-WY]+|[\.\-~])$/) {
#	my $error_code = 'invalid_Reference_seq_attribute_value_indel';
#	$self->warn($error_code, "($value) $feature_text");
#      }
#    }
#  }
#  else {
#    my $error_code = 'invalid_Reference_seq_attribute_value';
#    $self->throw($error_code, "($value) $feature_text");
#  }
#  if ($self->fasta) {
#    my $fasta_seq = $self->fasta->get_seq($seqid, $start, $end);
#    if ($fasta_seq ne $value) {
#      my $error_code = 'Reference_seq_attribute_value_ne_fasta_seq';
#      $self->throw($error_code, "($value ne $fasta_seq) $feature_text");
#    }
#  }
#  return 1;
#}
#
##-----------------------------------------------------------------------------
#
#=head2  is_valid_variant_reads_attribute
#
# Title   : is_valid_variant_reads_attribute
# Usage   : $self->is_valid_variant_reads_attribute($feature_hash)
# Function: Validates a variant_reads_attribute.
# Returns : 1 if valid undef otherwise.
# Args    : A string of text.
#
#=cut
#
#sub is_valid_variant_reads_attribute {
#
#  my $values = shift;
#  for my $value (@values) {
#    if (! $self->is_integer($value)) {
#      my $error_code = 'invalid_Variant_reads_attribute_value';
#      $self->throw($error_code, "($value) $feature_text");
#    }
#  }
#  return 1;
#}
#
##-----------------------------------------------------------------------------
#
#=head2  is_valid_total_reads_attribute
#
# Title   : is_valid_total_reads_attribute
# Usage   : $self->is_valid_total_reads_attribute($feature_hash)
# Function: Validates a total_reads_attribute.
# Returns : 1 if valid undef otherwise.
# Args    : A string of text.
#
#=cut
#
#sub is_valid_total_reads_attribute {
#
#  my $values = shift;
#  if (scalar @values > 1) {
#    my $error_code = 'invalid_Total_reads_attribute_multiple_values';
#    my $all_values = join ' ', @values;
#    $self->throw($error_code, "($all_values) $feature_text");
#  }
#  my $value = $values[0];
#  if (! $self->is_integer($value)) {
#    my $error_code = 'invalid_Total_reads_attribute_value';
#    $self->throw($error_code, "($value) $feature_text");
#  }
#  return 1;
#}
#
##-----------------------------------------------------------------------------
#
#=head2  is_valid_zygosity_attribute
#
# Title   : is_valid_zygosity_attribute
# Usage   : $self->is_valid_zygosity_attribute($feature_hash)
# Function: Validates a zygosity_attribute.
# Returns : 1 if valid undef otherwise.
# Args    : A string of text.
#
#=cut
#
#sub is_valid_zygosity_attribute {
#
#  my $values = shift;
#  for my $value (@values) {
#    if ($value !~ /^(hetero|homo|hemi)zygous$/) {
#      my $error_code = 'invalid_Zygosity_attribute_value';
#      $self->throw($error_code, "($value) $feature_text");
#    }
#  }
#  return 1;
#}
#
##-----------------------------------------------------------------------------
#
#=head2  is_valid_variant_freq_attribute
#
# Title   : is_valid_variant_freq_attribute
# Usage   : $self->is_valid_variant_freq_attribute($feature_hash)
# Function: Validates a variant_freq_attribute.
# Returns : 1 if valid undef otherwise.
# Args    : A string of text.
#
#=cut
#
#sub is_valid_variant_freq_attribute {
#
#  my $values = shift;
#  for my $value (@values) {
#    if (! $self->is_real_number($value)) {
#      my $error_code = 'invalid_Variant_freq_attribute_value';
#      $self->throw($error_code, "($value) $feature_text");
#    }
#  }
#  return 1;
#}
#
##-----------------------------------------------------------------------------
#
#=head2  is_valid_variant_effect_attribute
#
# Title   : is_valid_variant_effect_attribute
# Usage   : $self->is_valid_variant_effect_attribute($feature_hash)
# Function: Validates a variant_effect_attribute.
# Returns : 1 if valid undef otherwise.
# Args    : A string of text.
#
#=cut
#
#sub is_valid_variant_effect_attribute {
#
#  my $values = shift;
#  for my $value (@values) {
#    my ($sequence_variant, $index, $sequence_feature, @featureIDs) = split /\s+/, $value;
#    # TODO: Validate $sequence_variant
#    if (! $self->is_integer($index) || $index > (scalar @{$attributes->{lVariant_seq}} - 1)) {
#      my $error_code = 'invalid_Variant_effect_attribute_index_value';
#      $self->throw($error_code, "($index) $feature_text");
#    }
#    if (! $self->is_valid_variant_effect($sequence_variant)) {
#      my $error_code = 'invalid_Variant_effect_attribute_seq_var_id_value';
#      $self->throw($error_code, "($sequence_variant) $feature_text");
#    }
#    if (! $self->is_valid_sequence_feature($sequence_feature)) {
#      my $error_code = 'invalid_Variant_effect_attribute_seq_feat_id_value';
#      $self->throw($error_code, "($sequence_feature) $feature_text");
#    }
#    # TODO: Validate featureIDs
#  }
#  return 1;
#}
#
##-----------------------------------------------------------------------------
#
#=head2  is_valid_start_range_attribute
#
# Title   : is_valid_start_range_attribute
# Usage   : $self->is_valid_start_range_attribute($feature_hash)
# Function: Validates a start_range_attribute.
# Returns : 1 if valid undef otherwise.
# Args    : A string of text.
#
#=cut
#
#sub is_valid_start_range_attribute {
#
#  my $values = shift;
#  if (scalar @values != 2) {
#    my $error_code = 'invalid_Start_range_attribute_must_have_two_values';
#    my $all_values = join ',', @values;
#    $self->throw($error_code, "($all_values) $feature_text");
#  }
#  if (is_integer($values[0]) && $values[0] > $start) {
#    my $error_code = 'invalid_Start_range_attribute_values_must_contain_start';
#    my $all_values = join ',', @values;
#    $self->throw($error_code, "($all_values, $start) $feature_text");
#  }
#  if (is_integer($values[1]) && $values[1] < $start) {
#    my $error_code = 'invalid_Start_range_attribute_values_must_contain_start';
#    my $all_values = join ',', @values;
#    $self->throw($error_code, "($all_values, $start) $feature_text");
#  }
#  for my $value (@values) {
#    if ($value ne '.' && ! $self->is_integer($value)) {
#      my $error_code = 'invalid_Start_range_attribute_value';
#      $self->throw($error_code, "($value) $feature_text");
#    }
#  }
#  return 1;
#}
#
##-----------------------------------------------------------------------------
#
#=head2  is_valid_end_range_attribute
#
# Title   : is_valid_end_range_attribute
# Usage   : $self->is_valid_end_range_attribute($feature_hash)
# Function: Validates a end_range_attribute.
# Returns : 1 if valid undef otherwise.
# Args    : A string of text.
#
#=cut
#
#sub is_valid_end_range_attribute {
#
#  my $values = shift;
#  if (scalar @values != 2) {
#    my $error_code = 'invalid_End_range_attribute_must_have_2_values';
#    my $all_values = join ',', @values;
#    $self->throw($error_code, "($all_values) $feature_text");
#  }
#  if (is_integer($values[0]) && $values[0] > $end) {
#    my $error_code = 'invalid_End_range_attribute_values_must_contain_end';
#    my $all_values = join ',', @values;
#    $self->throw($error_code, "($all_values, $end) $feature_text");
#  }
#  for my $value (@values) {
#    if ($value ne '.' && ! $self->is_integer($value)) {
#      my $error_code = 'invalid_End_range_attribute_value';
#      $self->throw($error_code, "($value) $feature_text");
#    }
#  }
#  return 1;
#}
#
##-----------------------------------------------------------------------------
#
#=head2  is_valid_phased_attribute
#
# Title   : is_valid_phased_attribute
# Usage   : $self->is_valid_phased_attribute($feature_hash)
# Function: Validates a phased_attribute.
# Returns : 1 if valid undef otherwise.
# Args    : A string of text.
#
#=cut
#
#sub is_valid_phased_attribute {
#
#  my $values = shift;
#  if (scalar @values > 1) {
#    my $error_code = 'invalid_Phased_attribute_multiple_values';
#    my $all_values = join ' ', @values;
#    $self->throw($error_code, "($all_values) $feature_text");
#  }
#  my $value = $values[0];
#  # TODO: Validate the Phased attribute
#  return 1;
#}
#
##-----------------------------------------------------------------------------
#
#=head2  is_valid_genotype_attribute
#
# Title   : is_valid_genotype_attribute
# Usage   : $self->is_valid_genotype_attribute($feature_hash)
# Function: Validates a genotype_attribute.
# Returns : 1 if valid undef otherwise.
# Args    : A string of text.
#
#=cut
#
#sub is_valid_genotype_attribute {
#
#  my $values = shift;
#  for my $value (@values) {
#    if ($value !~ /^\d+:\d+$/) {
#      my $error_code = 'invalid_Genotype_attribute_value';
#      $self->throw($error_code, "($value) $feature_text");
#    }
#  }
#  return 1;
#
#}
#
##-----------------------------------------------------------------------------
#
#=head2  is_valid_individual_attribute
#
# Title   : is_valid_individual_attribute
# Usage   : $self->is_valid_individual_attribute($feature_hash)
# Function: Validates a individual_attribute.
# Returns : 1 if valid undef otherwise.
# Args    : A string of text.
#
#=cut
#
#sub is_valid_individual_attribute {
#
#  my $values = shift;
#  for my $value (@values) {
#    if (! $self->is_integer($value)) {
#      my $error_code = 'invalid_Individual_attribute_value';
#      $self->throw($error_code, "($value) $feature_text");
#    }
#  }
#  return 1;
#}
#
##-----------------------------------------------------------------------------
#
#=head2  is_valid_variant_copy_number_attribute
#
# Title   : is_valid_variant_copy_number_attribute
# Usage   : $self->is_valid_variant_copy_number_attribute($feature_hash)
# Function: Validates a variant_copy_number_attribute.
# Returns : 1 if valid undef otherwise.
# Args    : A string of text.
#
#=cut
#
#sub is_valid_variant_copy_number_attribute {
#
#  my $values = shift;
#  for my $value (@values) {
#    if (! $self->is_integer($value)) {
#      my $error_code = 'invalid_Variant_copy_number_attribute_value';
#      $self->throw($error_code, "($value) $feature_text");
#    }
#  }
#  return 1;
#}
#
##-----------------------------------------------------------------------------
#
#=head2  is_valid_reference_copy_number_attribute
#
# Title   : is_valid_reference_copy_number_attribute
# Usage   : $self->is_valid_reference_copy_number_attribute($feature_hash)
# Function: Validates a reference_copy_number_attribute.
# Returns : 1 if valid undef otherwise.
# Args    : A string of text.
#
#=cut
#
#sub is_valid_reference_copy_number_attribute {
#
#  my $values = shift;
#  if (scalar @values > 1) {
#    my $error_code = 'invalid_Reference_copy_number_attribute_multiple_values';
#    my $all_values = join ' ', @values;
#    $self->throw($error_code, "($all_values) $feature_text");
#  }
#  my $value = $values[0];
#  for my $value (@values) {
#    if (! $self->is_integer($value)) {
#      my $error_code = 'invalid_Reference_copy_number_attribute_value';
#      $self->throw($error_code, "($value) $feature_text");
#    }
#  }
#  return 1;
#}
#
##-----------------------------------------------------------------------------
#
#=head2  is_valid_variant_codon_attribute
#
# Title   : is_valid_variant_codon_attribute
# Usage   : $self->is_valid_variant_codon_attribute($feature_hash)
# Function: Validates a variant_codon_attribute.
# Returns : 1 if valid undef otherwise.
# Args    : A string of text.
#
#=cut
#
#sub is_valid_variant_codon_attribute {
#
#  my $values = shift;
#  for my $value (@values) {
#    if (! $self->is_iupac_nt($value)) {
#      my $error_code = 'invalid_Variant_codon_attribute_value';
#      $self->throw($error_code, "($value) $feature_text");
#    }
#  }
#  return 1;
#}
#
##-----------------------------------------------------------------------------
#
#=head2  is_valid_reference_codon_attribute
#
# Title   : is_valid_reference_codon_attribute
# Usage   : $self->is_valid_reference_codon_attribute($feature_hash)
# Function: Validates a reference_codon_attribute.
# Returns : 1 if valid undef otherwise.
# Args    : A string of text.
#
#=cut
#
#sub is_valid_reference_codon_attribute {
#
#  my $values = shift;
#  for my $value (@values) {
#    if (! $self->is_iupac_nt($value)) {
#      my $error_code = 'invalid_Reference_codon_attribute_value';
#      $self->throw($error_code, "($value) $feature_text");
#    }
#  }
#  return 1;
#}
#
##-----------------------------------------------------------------------------
#
#=head2  is_valid_variant_aa_attribute
#
# Title   : is_valid_variant_aa_attribute
# Usage   : $self->is_valid_variant_aa_attribute($feature_hash)
# Function: Validates a variant_aa_attribute.
# Returns : 1 if valid undef otherwise.
# Args    : A string of text.
#
#=cut
#
#sub is_valid_variant_aa_attribute {
#
#  my $values = shift;
#  for my $value (@values) {
#    if (! $self->is_aa_sequence($value)) {
#      my $error_code = 'invalid_Variant_aa_attribute_value';
#      $self->throw($error_code, "($value) $feature_text");
#    }
#  }
#  return 1;
#}
#
##-----------------------------------------------------------------------------
#
#=head2  is_valid_reference_aa_attribute
#
# Title   : is_valid_reference_aa_attribute
# Usage   : $self->is_valid_reference_aa_attribute($feature_hash)
# Function: Validates a reference_aa_attribute.
# Returns : 1 if valid undef otherwise.
# Args    : A string of text.
#
#=cut
#
#sub is_valid_reference_aa_attribute {
#
#  my $values = shift;
#  for my $value (@values) {
#    if (! $self->is_aa_sequence($value)) {
#      my $error_code = 'invalid_Reference_aa_attribute_value';
#      $self->throw($error_code, "($value) $feature_text");
#    }
#  }
#  return 1;
#}
#
##-----------------------------------------------------------------------------
#
#=head2  is_valid_breakpoint_detail_attribute
#
# Title   : is_valid_breakpoint_detail_attribute
# Usage   : $self->is_valid_breakpoint_detail_attribute($feature_hash)
# Function: Validates a breakpoint_detail_attribute.
# Returns : 1 if valid undef otherwise.
# Args    : A string of text.
#
#=cut
#
#sub is_valid_breakpoint_detail_attribute {
#
#  my $values = shift;
#  my $value = $values[0];
#  my ($seqid, $start_end, $strand) = split /:/, $value;
#  my ($start, $end) = split /\-/, $start_end;
#  if (! $self->is_integer($start)) {
#    my $error_code = 'invalid_start_in_Breakpoint_detail_attribute';
#    $self->warn($error_code, $feature_text);
#  }
#  if ($end && ! $self->is_integer($end)) {
#    my $error_code = 'invalid_end_in_Breakpoint_detail_attribute';
#    $self->warn($error_code, $feature_text);
#  }
#  if ($strand && $strand !~ /^(=|\-)$/) {
#    my $error_code = 'invalid_strand_in_Breakpoint_detail_attribute';
#    $self->warn($error_code, $feature_text);
#  }
#  return 1;
#}
#
##-----------------------------------------------------------------------------
#
#=head2  is_valid_sequence_context_attribute
#
# Title   : is_valid_sequence_context_attribute
# Usage   : $self->is_valid_sequence_context_attribute($feature_hash)
# Function: Validates a sequence_context_attribute.
# Returns : 1 if valid undef otherwise.
# Args    : A string of text.
#
#=cut
#
#sub is_valid_sequence_context_attribute {
#
#  my $values = shift;
#  for my $value (@values) {
#    if (! $self->is_iupac_nt($value)) {
#      my $error_code = 'invalid_Sequence_context_attribute_value';
#      $self->throw($error_code, "($value) $feature_text");
#    }
#  }
#  return 1;
#}
#
##-----------------------------------------------------------------------------
#
#=head2  is_valid_all
#
#
# Title   : is_valid_all
# Usage   : $self->is_valid_all($feature_hash)
# Function: Validates all attributes
# Returns : 1 if valid undef otherwise.
# Args    : A string of text.
#
#=cut
#
#sub is_valid_all {
#
#  my $values = shift;
#  return 1;
#
#}
#
##-----------------------------------------------------------------------------
#
#=head2  is_valid
#
# Title   : is_valid
# Usage   : $self->is_valid($feature_hash)
# Function: Validates a ...
# Returns : 1 if valid undef otherwise.
# Args    : A string of text.
#
#=cut
#
#sub is_valid {
#
#
#}
#
##-----------------------------------------------------------------------------
#
#=head2  is_valid_reference_seq
#
# Title   : is_valid_reference_seq
# Usage   : $self->is_valid_reference_seq($feature_hash)
# Function: Validates a ...
# Returns : 1 if valid undef otherwise.
# Args    : A string of text.
#
#=cut
#
#sub is_valid_reference_seq {
#
#
#}
#
#
##-----------------------------------------------------------------------------
#
#=head2  is_valid_sequence_feature
#
# Title   : is_valid_sequence_feature
# Usage   : $self->is_valid_sequence_feature($feature_hash)
# Function: Validates a ...
# Returns : 1 if valid undef otherwise.
# Args    : A string of text.
#
#=cut
#
#sub is_valid_sequence_feature {
#
#
#}
#
##-----------------------------------------------------------------------------
#
#=head2  is_valid_sequence_alteration
#
# Title   : is_valid_sequence_alteration
# Usage   : $self->is_valid_sequence_alteration($feature_hash)
# Function: Validates a ...
# Returns : 1 if valid undef otherwise.
# Args    : A string of text.
#
#=cut
#
#sub is_valid_sequence_alteration {
#
#
#}
#
##-----------------------------------------------------------------------------
#
#=head2  is_valid_variant_effect
#
# Title   : is_valid_variant_effect
# Usage   : $self->is_valid_variant_effect($feature_hash)
# Function: Validates a ...
# Returns : 1 if valid undef otherwise.
# Args    : A string of text.
#
#=cut
#
#sub is_valid_variant_effect {
#
#
#}
#
##-----------------------------------------------------------------------------
#
#=head2  is_valid_variant_seq
#
# Title   : is_valid_variant_seq
# Usage   : $self->is_valid_variant_seq($feature_hash)
# Function: Validates a ...
# Returns : 1 if valid undef otherwise.
# Args    : A string of text.
#
#=cut
#
#sub is_valid_variant_seq {
#
#
#}

#-----------------------------------------------------------------------------

=head1 DIAGNOSTICS

=over

=item C<< invalid_arguments_to_prepare_args >>

C<GAL::Base::prepare_args> accepts an array, a hash or a reference to
either an array or hash, but it was passed something different.

=item C<< invalid_ipuac_nucleotide_code >>

C<GAL::Base::expand_iupac_nt_codes> was passed a charachter that is
not a valid IUPAC nucleotide code
(http://en.wikipedia.org/wiki/Nucleic_acid_notation).

=item C<< failed_to_load_module >>

C<GAL::Base::load_module> was unable to load (require) the specified
module.  The module may not be installed or it may have a compile time
error.

=item C<< invalid_arguments_to_get_feature_bins >>

C<GAL::Base::get_feature_bins> was called with invalid arguments.  It
must have either a hash with the keys qw(seqid start end) or an object
with those same keys as methods.  The error message will try to give
you some idea of what arguments were passed.

=item C<< invalid_aa_datum_code >>

An invalid amino acid code was passed to
C<GAL::Base::amino_acid_data>.  Single-letter or three-letter amino
acid codes are required.

=back

=head1 CONFIGURATION AND ENVIRONMENT

L<GAL::Base> requires no configuration files or environment variables.

=head1 DEPENDENCIES

=over

=item L<Carp> qw(croak cluck)

=item L<Bio::DB::Fasta>

=back

=head1 INCOMPATIBILITIES

None reported.

=head1 BUGS AND LIMITATIONS

No bugs have been reported.

Please report any bugs or feature requests to:
barry.moore@genetics.utah.edu

=head1 AUTHOR

Barry Moore <barry.moore@genetics.utah.edu>

=head1 LICENCE AND COPYRIGHT

Copyright (c) 2010-2014, Barry Moore <barry.moore@genetics.utah.edu>.
All rights reserved.

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
