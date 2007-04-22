# -*- Mode: CPerl -*-
# File: t/common.plt
# Description: re-usable test subs for Math::PartialOrder
use Test;
$| = 1;

# isok($label,@_) -- prints helpful label
sub isok {
  my $label = shift;
  print "$label:\n";
  ok(@_);
}

# skipok($label,$skip_if_true,@_) -- prints helpful label
sub skipok {
  my ($label,$skip_if_true) = splice(@_,0,2);
  print "$label:\n";
  skip($skip_if_true,@_);
}

# ulistok($label,\@got,\@expect)
# --> ok() for unsorted lists
sub ulistok {
  my ($label,$l1,$l2) = @_;
  isok($label,join(',',sort(@$l1)),join(',',sort(@$l2)));
}

print "common.plt loaded.\n";

1;

