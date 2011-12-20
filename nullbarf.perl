#!/usr/bin/perl -w

use PDL;

##----------------------------------------------------------------------
sub nullbarf {
  my $a     = sequence(4,3);
  my $vals  = null;

  ## everything works as expected with a null-dimensioned empty pidde:
  my $i_0x0 = pdl([]);
  $a->indexND($i_0x0) = $vals; ##-- ok, expensive null-op

  ## ... but a 2x0 empty piddle pukes
  my $i_2x0 = sequence(2,1)->dice_axis(1,null); ##-- create an empty 2x0 piddle

  ## the next line barf()s with:
  #> Problem with assignment: PDL::Ops::assgn(a,b): Parameter 'b':
  #> Mismatched implicit thread dimension 0: should be 0, is 2
  #> 	 at Basic/Core/Core.pm.PL (i.e. PDL::Core.pm) line 256
  #> 	PDL::Core::barf('PDL=SCALAR(0xa026630)', 'PDL=SCALAR(0xa026690)') called at Basic/Core/Core.pm.PL (i.e. PDL::Core.pm) line 701
  #> 	eval {...} called at Basic/Core/Core.pm.PL (i.e. PDL::Core.pm) line 700
  #> 	__ANON__[Basic/Core/Core.pm.PL (i.e. PDL::Core.pm):712]('PDL=SCALAR(0xa026690)', 'PDL=SCALAR(0xa026630)', undef) called at (eval 160)[/usr/share/perl/5.10/perl5db.pl:638] line 2
  $a->indexND($i_2x0) .= $vals;
}
nullbarf();
