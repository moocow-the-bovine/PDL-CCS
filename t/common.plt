# -*- Mode: CPerl -*-
# File: t/common.plt
# Description: re-usable test subs
use Test::More;
BEGIN { $| = 1; }

# isok($label,@_) -- prints helpful label
sub isok {
  my $label = shift;
  if (@_==1) {
    ok($_[0],$label);
  } elsif (@_==2) {
    is($_[0],$_[1], $label);
  } else {
    die("isok(): expected 1 or 2 non-label arguments, but got ", scalar(@_));
  }
}

# skipok($label,$skip_if_true,@_) -- prints helpful label
sub skipok {
  my ($label,$skip_if_true) = splice(@_,0,2);
  if ($skip_if_true) {
    isok("skip:$label",1);
  } else {
    isok($label,@_);
  }
}

# ulistok($label,\@got,\@expect)
# --> ok() for unsorted lists
sub ulistok {
  my ($label,$l1,$l2) = @_;
  is_deeply([sort @$l1],[sort @$l2],$label);
}

# pdlok($label, $got, $want)
sub pdlok {
  my ($label,$got,$want) = @_;
  isok($label,
       defined($got) && defined($want)
       && $got->ndims==$want->ndims
       && all(pdl([$got->dims])==pdl([$want->dims]))
       && all($want==$got));
}

# pdlapprox($label, $got, $want, $eps=1e-5)
sub pdlapprox {
  my ($label,$got,$want,$eps) = @_;
  $eps = 1e-5 if (!defined($eps));
  isok($label,
       defined($got) && defined($want)
       && $got->ndims==$want->ndims
       && all(pdl([$got->dims])==pdl([$want->dims]))
       && all($want->approx($got,$eps)));
}


#print "common.plt loaded.\n";

1;

