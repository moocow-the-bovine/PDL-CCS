#!/usr/bin/perl -wd

use lib qw(./blib/lib ./blib/arch);
use PDL;
use PDL::CCS::Utils;
use PDL::VectorValued;

BEGIN {
  $, = ' ';
  our $eps=1e-6;
}

##---------------------------------------------------------------------
## isok(): for test script hacking
sub isok {
  my ($label,$bool) = @_;
  print "test($label): ", ($bool ? "ok" : "NOT ok"), "\n";
}

##---------------------------------------------------------------------
## test: data

sub test_ccs_data {
  our $a = pdl([[1,0,0,2],
		[0,3,4,0],
		[5,0,0,6]]);
  our $p = $a;
}


##---------------------------------------------------------------------
## test: ccs: pointer

sub test_ccs_pointer {
  test_ccs_data();
  our $awhich = $a->whichND();
  our $avals  = $a->indexND($awhich);
  our ($aptr0,$awi0) = ccs_pointer($awhich->slice("(0),"));
  our ($aptr1,$awi1) = ccs_pointer($awhich->slice("(1),"));

  ##-- ok, we've got a raw N+1 pointer: build compatible CCS encoded matrix
  our $ccs_rowids = $awhich->slice("(1),")->index($awi0);
  our $ccs_vals   = $avals->index($awi0);
  our (@ccs);
  our ($ptr,$rowids,$nzvals) = @ccs = ($ptr0->slice("0:-2"), $ccs_rowids, $ccs_vals);
}
test_ccs_pointer;


##---------------------------------------------------------------------
## DUMMY
##---------------------------------------------------------------------
foreach $i (0..3) {
  print "--dummy($i)--\n";
}

