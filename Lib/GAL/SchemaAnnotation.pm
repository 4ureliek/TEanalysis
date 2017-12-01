package GAL::SchemaAnnotation;
use base qw(DBIx::Class);
use Scalar::Util qw(weaken);

sub annotation {
  my ($self, $annotation) = @_;
  if (! $self->{annotation} && $annotation) {
    weaken($annotation);
    $self->{annotation} = $annotation;
  }
  return $self->{annotation};
}
1;
