# -*- Mode: CPerl -*-
# t/02_encode.t: test ccs encoding
use Test::More;
use strict;
use warnings;

##-- common subs
my $TEST_DIR;
BEGIN {
  use File::Basename;
  use Cwd;
  $TEST_DIR = Cwd::abs_path dirname( __FILE__ );
  eval qq{use lib ("$TEST_DIR/$_/blib/lib","$TEST_DIR/$_/blib/arch");} foreach (qw(../../.. ../.. ..));
  do "$TEST_DIR/common.plt" or  die("$0: failed to load $TEST_DIR/common.plt: $@");
}

##-- common modules
use PDL;
use PDL::CCS::MatrixOps;

##-- common code
sub mccs {
  my $m = shift;
  my $mwhich = $m->whichND->qsortvec;
  my $mvals  = $m->indexND($mwhich);
  return ($mwhich, $mvals);
}

my ($m,$u,$v,$zc,$c);
my $BAD = pdl(0)->setvaltobad(0);

##--------------------------------------------------------------
## test output dims for matmult2d_sdd
sub test_matmult2d_sdd_dims {
  my $a = pdl([[1,0],[0,1],[1,1]]);   # (M=2, N=3)
  my $b = pdl([[0,1,0,1],[1,0,1,0]]); # (O=4, M=2)
  my $c_want = ($a x $b);             # (O=4, M=2)

  my $missing = pdl(0);                # ()           # cannot be (1)!
  my $zc = zeroes($b->dim(0));         # (O=4)
  my ($ixa, $nza) = mccs($a);          # ~(M=2, N=3)

  my $c_got = null;
  ccs_matmult2d_sdd($ixa,$nza,$missing, $b,$zc, $c_got, $a->dim(1));
  pdlok('matmult2d_sdd dims:(c=null):dims', pdl(indx, $c_got->dims), pdl(indx, $c_want->dims));
  pdlok('matmult2d_sdd dims:(c=null):vals', $c_want, $c_got);
}
test_matmult2d_sdd_dims();

##--------------------------------------------------------------
## test matmult2_sdd
sub test_matmult2d_sdd {
  my ($label, $a,$b,$missing) = @_;
  $missing = PDL->topdl($missing)->squeeze;

  my $atmp = $a->pdl;
  $atmp->where($a==0) .= $missing;
  my $c_want = ($atmp x $b);

  my ($ixa, $nza) = mccs($a);
  my $zc = ((zeroes($a->dim(0),1) + $missing) x $b)->flat;
  ccs_matmult2d_sdd($ixa,$nza,$missing, $b,$zc, (my $c_got=null), $a->dim(1));

 SKIP: {
    skip("matmult2d_sdd can't handle missing=BAD", 1)
      if ($missing->isbad);

    pdlok("matmult2d_sdd:$label:missing=$missing", $c_got, $c_want);
  }
}

$m = identity(3)->set(2,2,0);  # [[1,0,0],[0,1,0],[0,0,0]]
$v = zeroes(1,3)->set(0,2,1);  # [[0,0,1]]

foreach my $missing (0,42,$BAD) {
  test_matmult2d_sdd("simple", $m,$v,$missing);
  test_matmult2d_sdd("xvals x yvals", $m->xvals, $v->yvals, $missing);
  test_matmult2d_sdd("yvals x xvals", $m->yvals, $v->xvals, $missing);
  test_matmult2d_sdd("seq x seq", $m->sequence, $v->sequence, $missing);
}


##--------------------------------------------------------------
## test matmult2_zdd
sub test_matmult2d_zdd {
  my ($label, $a,$b) = @_;

  my $c_want = ($a x $b);
  my ($ixa, $nza) = mccs($a);
  ccs_matmult2d_zdd($ixa,$nza, $b, (my $c_got=null), $a->dim(1));

  pdlok("matmult2d_zdd:$label", $c_got, $c_want);
}

$m = identity(3)->set(2,2,0);  # [[1,0,0],[0,1,0],[0,0,0]]
$v = zeroes(1,3)->set(0,2,1);  # [[0,0,1]]

test_matmult2d_zdd("simple", $m,$v);
test_matmult2d_zdd("xvals x yvals", $m->xvals, $v->yvals);
test_matmult2d_zdd("yvals x xvals", $m->yvals, $v->xvals);
test_matmult2d_zdd("seq x seq", $m->sequence, $v->sequence);

##--------------------------------------------------------------
## test vnorm

sub test_vnorm {
  my ($label, $a, $dim) = @_;

  my $vnorm_want = ($a**2);
  $vnorm_want = $vnorm_want->xchg(0,1) if ($dim==0);
  $vnorm_want = $vnorm_want->sumover->sqrt;

  my ($awhich, $avals) = mccs($a);
  my $acols = $awhich->slice("($dim),");
  ccs_vnorm($acols, $avals, (my $vnorm_got=null), $a->dim($dim));

  pdlapprox("vnorm:$label:dim=$dim", $vnorm_got, $vnorm_want, 1e-5);
}
$m = pdl([[1,1,-1],[2,2,-2],[3,2,-3],[4,1,-4]]);
test_vnorm('simple', $m, 0);
test_vnorm('simple', $m, 1);


##--------------------------------------------------------------
done_testing;
