# -*- Mode: Perl -*-
# t/03_ops.t: test ccs native operations

$TEST_DIR = './t';
#use lib qw(../blib/lib ../blib/arch); $TEST_DIR = '.'; # for debugging

# load common subs
use Test;
do "$TEST_DIR/common.plt";
use PDL;
use PDL::CCS;

BEGIN { plan tests=>6, todo=>[]; }

##-- setup
$a = pdl(double, [
		  [10,0,0,0,-2,0],
		  [3,9,0,0,0,3],
		  [0,7,8,7,0,0],
		  [3,0,8,7,5,0],
		  [0,8,0,9,9,13],
		  [0,4,0,0,2,-1],
		 ]);

($ptr,$rowids,$nzvals) = ccsencode($a);

##-- 1: transpose()
($ptrT,$rowidsT,$nzvalsT) = ccstranspose($ptr,$rowids,$nzvals);
$a2 = ccsdecode($ptrT,$rowidsT,$nzvalsT)->xchg(0,1);
isok("transpose()", all($a==$a2));

##-- 2-3: whichND()
($ccols,$crows) = ccswhichND($ptr,$rowids,$nzvals);
($acols,$arows) = whichND($a);
$acoli = $acols->qsorti;
$ccoli = $ccols->qsorti;
isok("whichND() / cols", all($acols->index($acoli) == $ccols->index($ccoli)));
isok("whichND() / rows", all($arows->index($acoli) == $crows->index($ccoli)));

##-- 4: which()
$awhich = which($a)->qsort;
$cwhich = ccswhich($ptr,$rowids,$nzvals)->qsort;
isok("which() / flat", all($awhich==$cwhich));

##-- 5: ccsmult_rv (row vector)
$rv=10**(sequence($a->dim(0))+1);
$nzvals_rv = ccsmult_rv($ptr,$rowids,$nzvals, $rv);
isok("ccsmult_rv()", all(($a * $rv)==ccsdecode($ptr,$rowids,$nzvals_rv)));

##-- 6: ccsmult_cv (col vector)
$cv=10**(sequence($a->dim(1))+1);
$nzvals_cv = ccsmult_cv($ptr,$rowids,$nzvals, $cv);
isok("ccsmult_cv()", all(($a * $cv->slice("*1,"))==ccsdecode($ptr,$rowids,$nzvals_cv)));

print "\n";
# end of t/03_ops.t

