# -*- Mode: CPerl -*-
# t/01_encode.t

$TEST_DIR = './t';
#use lib qw(../blib/lib ../blib/arch); $TEST_DIR = '.'; # for debugging

# load common subs
use Test;
do "$TEST_DIR/common.plt";
use PDL;
use PDL::CCS::Nd;
use PDL::VectorValued;

BEGIN { plan tests=>83, todo=>[]; }

## (i+1)..(i+9): basic properites (missing==0)
sub test_basic {
  my ($label,$a,$ccs,$missing) = @_;

  isok("${label}:defined", defined($ccs));
  isok("${label}:dims",    all(pdl($ccs->dims)==pdl($a->dims)));
  isok("${label}:nelem",   $ccs->nelem==$a->nelem);

  ##-- check missing
  $missing = 0 if (!defined($missing));
  $missing = PDL->topdl($missing);
  if ($missing->isbad) {
    $awhichND = whichND(!isbad($a));
  } else {
    $awhichND = whichND($a!=$missing);
  }

  isok("${label}:_nnz",    $ccs->_nnz==$awhichND->dim(1));
  isok("${label}:whichND", all($ccs->whichND->vv_qsortvec==$awhichND->vv_qsortvec));
  isok("${label}:nzvals",  all(matchpdl($ccs->whichVals, $a->indexND(scalar($ccs->whichND)))));
  isok("${label}:missing:value", matchpdl($ccs->missing, $missing));

  ##-- testdecode
  isok("${label}:decode",  all(matchpdl($ccs->decode,$a)));
  isok("${label}:todense", all(matchpdl($ccs->todense,$a)));
}


##--------------------------------------------------------------
## missing==0

##-- 1*nbasic: newFromDense(): basic properties
$ccs = PDL::CCS::Nd->newFromDense($a);
test_basic("newFromDense:missing=0", $a, $ccs, 0);

##-- 2*nbasic: toccs(): basic properties
$ccs = $a->toccs;
test_basic("toccs:missing=0", $a, $ccs, 0);

##-- 3*nbasic: newFromWhich()
$ccs = PDL::CCS::Nd->newFromWhich($awhich,$avals,missing=>0);
test_basic("newFromWhich:missing=0", $a, $ccs, 0);

##--------------------------------------------------------------
## missing==BAD

##-- 5*nbasic: newFromDense(...BAD): basic properties
our ($abad);
$a = $a->setbadif($abad);
test_basic("newFromDense:missing=BAD:explicit", $a, PDL::CCS::Nd->newFromDense($a,$BAD), $BAD);
test_basic("newFromDense:missing=BAD:implicit", $a, PDL::CCS::Nd->newFromDense($a),      $BAD);

##-- 7*nbasic: toccs(...BAD): basic properties
test_basic("toccs:missing=BAD:explicit", $a, $a->toccs($BAD), $BAD);
test_basic("toccs:missing=BAD:implicit", $a, $a->toccs(),     $BAD);

##-- 9*nbasic: newFromWhich(...BAD)
test_basic("newFromWhich:missing=BAD:explicit", $a, PDL::CCS::Nd->newFromWhich($awhich,$avals,missing=>$BAD), $BAD);
test_basic("newFromWhich:missing=BAD:implicit", $a, PDL::CCS::Nd->newFromWhich($awhich,$avals),               $BAD);

##--------------------------------------------------------------
## global tests
##  (9*nbasic)..((9*nbasic)+2)

## 1..2: PDL->todense, PDL::CCS::Nd->toccs
isok("PDL::todense():no-copy", overload::StrVal($a)   eq overload::StrVal($a->todense));
isok("CCS::toccs():no-copy",   overload::StrVal($ccs) eq overload::StrVal($ccs->toccs));

print "\n";
# end of t/*.t

