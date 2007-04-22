# -*- Mode: CPerl -*-
# t/02_encode.t: test ccs encoding

$TEST_DIR = './t';
#use lib qw(../blib/lib ../blib/arch); $TEST_DIR = '.'; # for debugging

# load common subs
use Test;
do "$TEST_DIR/common.plt";
use PDL;
use PDL::CCS::Utils;
use PDL::VectorValued;

BEGIN { plan tests=>6, todo=>[]; }

##-- setup
our $a = pdl(double, [
		      [10,0,0,0,-2],
		      [3,9,0,0,0],
		      [0,7,8,7,0],
		      [3,0,8,7,5],
		      [0,8,0,9,9],
		      [0,4,0,0,2],
		     ]);

##-- test: encode_pointers
our $awhich = $a->whichND()->vv_qsortvec;
our $avals  = $a->indexND($awhich);
our ($aptr0,$awi0) = ccs_encode_pointers($awhich->slice("(0),"));
our ($aptr1,$awi1) = ccs_encode_pointers($awhich->slice("(1),"));

##-- 1..2
isok("whichND", all($awhich==pdl([[0,0],[0,1],[0,3],[1,1],[1,2],[1,4],[1,5],[2,2],[2,3],[3,2],[3,3],[3,4],[4,0],[4,3],[4,4],[4,5]])));
isok("vals", all($avals==pdl([10,3,3,9,7,8,4,8,8,7,7,9,-2,5,9,2])));

##-- 3..4: ptr0
isok("ccs_encode_pointers:ptr0",    all($aptr0==pdl([0,3,7,9,12,16])));
isok("ccs_encode_pointers:awi0",    all($awi0 ==$awi0->sequence));

##-- 5..6: ptr1
isok("ccs_encode_pointers:ptr1",    all($aptr1==pdl([0,2,4,7,11,14,16])));
our $awi1x = $awhich->slice("(1),")->index($awi1);
isok("ccs_encode_pointers:awi1",    all($awi1x==$awi1x->qsort));


print "\n";
# end of t/02_encode.t

