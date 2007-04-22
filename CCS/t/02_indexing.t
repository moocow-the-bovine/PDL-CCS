# -*- Mode: CPerl -*-
# t/02_indexing.t

$TEST_DIR = './t';
#use lib qw(../blib/lib ../blib/arch); $TEST_DIR = '.'; # for debugging

# load common subs
use Test;
do "$TEST_DIR/common.plt";
use PDL;
use PDL::CCS::Nd;
use PDL::VectorValued;

BEGIN { plan tests=>18, todo=>[]; }

##--------------------------------------------------------------
## missing==0

$ccs = $a->toccs;

##-- 1: which
isok("which:flat", all($ccs->which->qsort == $a->which->qsort));

##-- 2: index (flat)     ------------------> NO PSEUDO-THREADING!
$find = pdl(long(2,4,6,8));
isok("index:flat", all($ccs->index($find)==$a->flat->index($find)));

##-- 3: indexND
$find = pdl(long,[[0,0],[1,0],[1,1]]);
isok("indexND", all($a->indexND($find)==$ccs->indexND($find)));

##-- 4..5: dice_axis
$axisi = pdl(long,[2,4]);
isok("dice_axis(0)", all( matchpdl($a->dice_axis(0,$axisi), $ccs->dice_axis(0,$axisi)->decode) ));
isok("dice_axis(1)", all( matchpdl($a->dice_axis(1,$axisi), $ccs->dice_axis(1,$axisi)->decode) ));

##-- 6..8: at,set
@nzindex = (4,3);
@zindex  = (3,1);
isok("at():nz", $ccs->at(@nzindex)==$a->at(@nzindex));
isok("at:z",    $ccs->at(@zindex)==$a->at(@zindex));
isok("set():nz", all( matchpdl($ccs->set(@nzindex,42)->decode, $a->set(@nzindex,42)) ));

##-- 9..10: reorder
isok("reorder(1,0)",             all($ccs->reorder(1,0)->decode==$a->reorder(1,0)));
isok("post-reorder(1,0):decode", all($ccs->decode==$a));

##-- 11..12: xchg(0,1)
isok("xchg(0,1)",                all($ccs->xchg(0,1)->decode==$a->xchg(0,1)));
isok("post-xchg(0,1):decode",    all($ccs->decode==$a));

##-- 13..14: xchg(0,-1)
isok("xchg(0,-1)",               all($ccs->xchg(0,-1)->decode==$a->xchg(0,-1)));
isok("post-xchg(0,-1):decode",   all($ccs->decode==$a));

##-- 15..16: mv(0,1)
isok("mv(0,1)",                  all($ccs->mv(0,1)->decode==$a->mv(0,1)));
isok("post-mv(0,1):decode",      all($ccs->decode==$a));

##-- 17..18: mv(1,0)
isok("mv(1,0)",                  all($ccs->mv(1,0)->decode==$a->mv(1,0)));
isok("post-mv(1,0):decode",      all($ccs->decode==$a));

print "\n";
# end of t/*.t

