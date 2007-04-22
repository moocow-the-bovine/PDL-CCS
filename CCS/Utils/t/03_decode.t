# -*- Mode: CPerl -*-
# t/03_encode.t: test ccs pointer-decoding

$TEST_DIR = './t';
#use lib qw(../blib/lib ../blib/arch); $TEST_DIR = '.'; # for debugging

# load common subs
use Test;
do "$TEST_DIR/common.plt";
use PDL;
use PDL::CCS::Utils;
use PDL::VectorValued;

BEGIN { plan tests=>8, todo=>[]; }

##-- setup
our $a = pdl(double, [
		      [10,0,0,0,-2],
		      [3,9,0,0,0],
		      [0,7,8,7,0],
		      [3,0,8,7,5],
		      [0,8,0,9,9],
		      [0,4,0,0,2],
		     ]);

##-- test: decode_pointer
our $awhich = $a->whichND;
our $awhich0 = $awhich->slice("(0)");
our $awhich1 = $awhich->slice("(1)");
our $avals  = $a->indexND($awhich);

##-- 1..2: decode_pointer: dim=0: full
our ($aptr0,$anzi0)     = ccs_encode_pointers($awhich0);
our $aproj0             = sequence(long,$a->dim(0));
our ($aproj0d,$apnzi0d) = ccs_decode_pointer($aptr0,$aproj0);
isok("ccs_decode_pointer:full:dim=0:proj",  all($aproj0d==$awhich0->qsort));
isok("ccs_decode_pointer:full:dim=0:nzi",   all($apnzi0d==$apnzi0d->sequence));

##-- 3..4: decode_pointer: dim=1: full
our ($aptr1,$anzi1)     = ccs_encode_pointers($awhich1);
our $aproj1             = sequence(long,$a->dim(1));
our ($aproj1d,$apnzi1d) = ccs_decode_pointer($aptr1,$aproj1);
isok("ccs_decode_pointer:full:dim=1:proj", all($aproj1d==$awhich1->qsort));
isok("ccs_decode_pointer:full:dim=1:nzi",  all($apnzi1d==$apnzi1d->sequence));

##-- 5..6: decode_pointer: dim=0: partial
$aproj0 = pdl(long,[1,2,4]);
our $aslice0 = $a->dice_axis(0,$aproj0);
($aproj0d,$apnzi0d) = ccs_decode_pointer($aptr0,$aproj0);

our $apnzi      = $anzi0->index($apnzi0d);
our $which_proj = $aproj0d->slice("*1,")->append($awhich->slice("1")->dice_axis(1,$apnzi));
our $vals_proj  = $avals->index($apnzi);

isok("ccs_decode_pointer:partial:dim=0:which", all($which_proj->vv_qsortvec==$aslice0->whichND->vv_qsortvec));
isok("ccs_decode_pointer:partial:dim=0:vals",  all($vals_proj==$aslice0->indexND($which_proj)));

##-- 7..8: decode_pointer: dim=1: partial
$aproj1 = pdl(long,[2,3,5]);
our $aslice1 = $a->dice_axis(1,$aproj1);
($aproj1d,$apnzi1d) = ccs_decode_pointer($aptr1,$aproj1);

$apnzi      = $anzi1->index($apnzi1d);
$which_proj = $aproj1d->slice("*1,")->append($awhich->slice("0")->dice_axis(1,$apnzi))->slice("-1:0");
$vals_proj  = $avals->index($apnzi);

isok("ccs_decode_pointer:partial:dim=1:which", all($which_proj->vv_qsortvec==$aslice1->whichND->vv_qsortvec));
isok("ccs_decode_pointer:partial:dim=1:vals",  all($vals_proj==$aslice1->indexND($which_proj)));


print "\n";
# end of t/03_decode.t

