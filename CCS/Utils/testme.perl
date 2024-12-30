#!/usr/bin/perl -w

use lib qw(./blib/lib ./blib/arch);
use PDL;
use PDL::CCS::Utils;
use PDL::VectorValued;
use Benchmark qw(timethese cmpthese);

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
  our ($ptr,$rowids,$nzvals) = @ccs = ($aptr0->slice("0:-2"), $ccs_rowids, $ccs_vals);
}
#test_ccs_pointer;

##---------------------------------------------------------------------
## test: ccs: xindex2d
sub test_ccs_xindex2d {
  test_ccs_data();
  our $awhich = $a->whichND()->vv_qsortvec;
  our $avals  = $a->indexND($awhich);

  ##-- ok, we've got a raw N+1 pointer: build compatible CCS encoded matrix
  our $ai   = $awhich->slice("(0),")->uniq;
  our $bi   = $awhich->slice("(1),")->uniq;
  our $abi0 = sequence(long,$awhich->dim(1));

  print xvals($awhich->dim(1),1)->glue(1,$awhich->xchg(0,1));
  our $abi = ccs_xindex2d($awhich, $ai,$bi);

  isok("ccs_xindex2d:nelem", $abi->nelem==$abi0->nelem);
  isok("ccs_xindex2d:vals",  all($abi==$abi0));

  ##-- now test with a large random sparse matrix
  srand(0);
  $a = rint(random(100,100)*100)->long;
  (my $tmp=$a->where($a >= 10)) .= 0;
  $awhich = $a->whichND()->vv_qsortvec;
  ##
  $ai = $awhich->slice("(0),")->uniq;
  $bi = $awhich->slice("(1),")->uniq;
  $abi0 = sequence(long,$awhich->dim(1));
  $abi  = ccs_xindex2d($awhich, $ai,$bi);
  #
  isok("ccs_xindex2d:rand:nelem", $abi->nelem==$abi0->nelem);
  isok("ccs_xindex2d:rand:vals",  all($abi==$abi0));
}
test_ccs_xindex2d();

##---------------------------------------------------------------------
## test: ccs: xindex2d (test2)
use PDL::IO::FastRaw;
sub test_ccs_index2d_2 {
  my %mopts = (ReadOnly=>1, Creat=>0);
  my $which = mapfraw("x2d-which.pdl", \%mopts);
  my $ai    = mapfraw("x2d-a.pdl", \%mopts);
  my $bi    = mapfraw("x2d-b.pdl", \%mopts);

  ##-- get "proper" values via vsearchvec
  my $wnd = $ai->slice("*1,")->cat($bi)->clump(2)->xchg(0,1);
  my $abi0      = $wnd->vsearchvec($which);
  my $abi0_mask = ($wnd==$which->dice_axis(1,$abi0))->andover;
  $abi0         = $abi0->where($abi0_mask);

  ##-- xindex2d
  print "xindex2d(which(".join(',',$which->dims)."); ai(".join(',',$ai->dims)."); bi(".join(',',$bi->dims).")\n";
  #print STDERR "Press ENTER to continue: "; $_=<STDIN>;
  my $abi  = ccs_xindex2d($which,$ai,$bi);

  ##-- check
  isok("ccs_xindex2d:2:nelem", $abi->nelem == $abi0->nelem);
  isok("ccs_xindex2d:2:vals",  all($abi==$abi0));

  exit 0;
}
#test_ccs_index2d_2();

##---------------------------------------------------------------------
## bench: ccs_xindex2d
sub bench_ccs_index2d {
  my ($m,$n,$density) = @_;
  $m ||= 100;
  $n ||= $m;
  $density ||= 0.01;

  ##-- create dummy matrix
  my $data = random($m,$n);
  (my $tmp=$data->where($data >= $density)) .= 0;

  ##-- get data
  my $which = $data->whichND->vv_qsortvec;
  my $ai    = $which->slice("(0),")->uniq;
  my $bi    = $which->slice("(1),")->uniq;

  ##-- test subs
  my $cmp_vsearchvec = sub {
    my $abw = $ai->slice("*1,")->cat($bi)->clump(2)->xchg(0,1);
    my $abi = $abw->vsearchvec($which);
    my $abi_mask = ($abw==$which->dice_axis(1,$abi))->andover;
    $abi         = $abi->where($abi_mask);
  };
  my $cmp_xindex2d = sub {
    my $abix = $which->ccs_xindex2d($ai,$bi);
  };
  print "BENCHMARK(m=$m, n=$n, density=$density):\n";
  cmpthese(-1, {'vsearchvec'=>$cmp_vsearchvec, 'xindex2d'=>$cmp_xindex2d});
  ##
  # BENCHMARK(m=100, n=100, density=.01):
  #              Rate vsearchvec   xindex2d
  # vsearchvec 2113/s         --       -73%
  # xindex2d   7859/s       272%         --
  ##
  # BENCHMARK(m=100, n=1000, density=.01):
  #             Rate vsearchvec   xindex2d
  # vsearchvec 131/s         --       -81%
  # xindex2d   704/s       438%         --
  ##
  # BENCHMARK(m=100, n=1000, density=.001):
  #              Rate vsearchvec   xindex2d
  # vsearchvec 1674/s         --       -75%
  # xindex2d   6698/s       300%         --
}
#bench_ccs_index2d(@ARGV); exit 0;

##---------------------------------------------------------------------
## DUMMY
##---------------------------------------------------------------------
foreach $i (0..3) {
  print "--dummy($i)--\n";
}

