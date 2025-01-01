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
use PDL::CCS::Utils;
use PDL::VectorValued;

##-- setup
my $m = pdl(double, [
                     [10,0,0,0,-2],
                     [3,9,0,0,0],
                     [0,7,8,7,0],
                     [3,0,8,7,5],
                     [0,8,0,9,9],
                     [0,4,0,0,2],
                    ]);

##--------------------------------------------------------------
## common
my $mwhich = $m->whichND()->vv_qsortvec;
my $mvals  = $m->indexND($mwhich);

##--------------------------------------------------------------
## test: xindex1d()
my $a0 = pdl(indx, [1,2,4]);
my $want = $m->dice_axis(0, $a0);
sub test_xindex1d {
  my ($label, @args) = @_;
  print "# xindex1d:$label\n";
  my $nzia = ccs_xindex1d($mwhich, $a0, @args);
  pdlok("ccs_index1d:$label:which[1]", $mwhich->dice_axis(1, $nzia)->slice("1:"), $want->whichND->vv_qsortvec->slice("1:"));
  pdlok("ccs_index1d:$label:vals", $mvals->index($nzia), $m->indexND($mwhich->dice_axis(1, $nzia)));

  my $got = $m->zeroes;
  $got->indexND($mwhich->dice_axis(1, $nzia)) .= $mvals->index($nzia);
  pdlok("ccs_index1d:$label:decoded", $got->dice_axis(0, $a0), $want);
}

test_xindex1d('no-outputs');
test_xindex1d('null-outputs', null, null);

test_xindex1d('prealloc-output-nzia', zeroes(indx, 16), null);
test_xindex1d('prealloc-output-nnza', null, pdl(indx, 16));
test_xindex1d('prealloc-outputs', zeroes(indx, 16), zeroes(indx, 1));

test_xindex1d('minalloc-output-nzia', zeroes(indx, 10), null);
test_xindex1d('minalloc-output-nnza', null, pdl(indx, 10));
test_xindex1d('minalloc-outputs', zeroes(indx, 10), zeroes(indx, 1));

##--------------------------------------------------------------
## test: xindex2d
my $ai = $mwhich->slice("(0),")->uniq;
my $bi = $mwhich->slice("(1),")->uniq;
my $abi_want = sequence(indx, $mwhich->dim(1));

sub test_xindex2d {
  my ($label, @args) = @_;
  print "# xindex2d:$label\n";
  my $abi = ccs_xindex2d($mwhich, $ai, $bi, @args);
  pdlok("ccs_xindex2d:$label:indices", $abi, $abi_want);
}

test_xindex2d('no-outputs');
test_xindex2d('null-outputs', null, null);

test_xindex2d('prealloc-output-ab', $abi_want->zeroes, null);
test_xindex2d('prealloc-output-nab', null, pdl(indx, $abi_want->nelem));
test_xindex2d('prealloc-outputs', $abi_want->zeroes, pdl(indx, $abi_want->nelem));
test_xindex2d('overalloc-outputs', zeroes($abi_want->nelem + 1), pdl(indx, $abi_want->nelem + 1));


##--------------------------------------------------------------
done_testing;
