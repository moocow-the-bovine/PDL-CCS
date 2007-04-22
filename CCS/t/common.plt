# -*- Mode: CPerl -*-
# File: t/common.plt
# Description: re-usable test subs for Math::PartialOrder
use Test;
use PDL;
$| = 1;

#-- basic data
BEGIN {
  our $a = pdl(double, [
			[10,0,0,0,-2],
			[3,9,0,0,0],
			[0,7,8,7,0],
			[3,0,8,7,5],
			[0,8,0,9,9],
			[0,4,0,0,2],
		       ]);
  our $abad   = ($a==0);
  our $agood  = !$abad;
  our $awhich = $a->whichND;
  our $avals  = $a->indexND($awhich);

  our $BAD = pdl(0)->setvaltobad(0);
}

sub matchpdl {
  my ($a,$b) = map {$_->setnantobad} @_[0,1];
  return ($a==$b)->setbadtoval(0) | ($a->isbad & $b->isbad);
}


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

