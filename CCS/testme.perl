#!/usr/bin/perl -w

use lib qw(./blib/lib ./blib/arch ../blib/lib ../blib/arch);
use PDL;
use PDL::IO::Storable;
use PDL::VectorValued;
#use PDL::CCS::Old;
use PDL::CCS::Utils;
use PDL::CCS::Ufunc;
use PDL::CCS::Ops;
use PDL::CCS::Functions;
use PDL::CCS::Compat;
use PDL::CCS::Nd;
use PDL::IO::Misc;
use PDL::Graphics::PGPLOT;
use Storable qw(store retrieve);
use PDL::IO::Storable;
use Benchmark qw(cmpthese timethese);
use Data::Dumper;

BEGIN {
  $, = ' ';
  our $eps=1e-6;
  select(STDOUT); $|=1;

  our $PDIMS = our $DIMS = $PDL::CCS::Nd::PDIMS;
  our $VDIMS = our $XDIMS = $PDL::CCS::Nd::VDIMS;
  our $WHICH = $PDL::CCS::Nd::WHICH;
  our $VALS  = $PDL::CCS::Nd::VALS;
  our $PTRS  = $PDL::CCS::Nd::PTRS;

  our $P_BYTE = PDL::byte();
  our $P_LONG = PDL::long();

  ##-- flags
  our $CCSND_BAD_IS_MISSING = $PDL::CCS::Nd::CCSND_BAD_IS_MISSING;
  our $CCSND_NAN_IS_MISSING = $PDL::CCS::Nd::CCSND_NAN_IS_MISSING;
  our $CCSND_INPLACE        = $PDL::CCS::Nd::CCSND_INPLACE;
  our $CCSND_FLAGS_DEFAULT  = $PDL::CCS::Nd::CCSND_FLAGS_DEFAULT;

  our $BAD    = PDL->zeroes(1)->setvaltobad(0)->squeeze;

  #$PDL::CCS::Nd::BINOP_BLOCKSIZE_MAX=8; ##-- DEBUG

}

##---------------------------------------------------------------------
## test subs
BEGIN {
  our $N_OK=0;
  our $N_NOTOK=0;
  our $NOTOK_DIE = 1;
  our $TEST_VERBOSE = 1;
}
sub isok {
  my ($label,$bool) = @_;
  if ($TEST_VERBOSE) { print "test($label): ", ($bool ? "ok" : "NOT ok"), "\n"; }
  if ($bool) { ++$N_OK; }
  else       { ++$N_NOTOK; }
  die("isok(): test failed for '$label'!") if ($NOTOK_DIE && !$bool);
}

sub matchpdl {
  my ($a,$b) = map {$_->setnantobad} @_[0,1];
  return ($a==$b)->setbadtoval(0) | ($a->isbad & $b->isbad);
}

sub pdlstr ($) {
  return "$_[0]" if (!ref($_[0]) || $_[0]->ndims==0);
  return '['.join(",", map {pdlstr($_)} $_[0]->dog).']';
}

sub make { system('make'); }

sub test_all {
  my $verbose=shift;
  $TEST_VERBOSE=$verbose if (defined($verbose));

  test_matmult_pre1();
  test_vdims_1();
  test_recode();
  test_nd_binop_sclr_all();
  test_nd_binop_cvrv_all();
  test_nd_binop_mia_all();
  test_nd_unop_all();
  test_nd_ufunc_all();

  our $ntested = $N_OK+$N_NOTOK;
  printf("\ntest_all: $N_OK of $ntested tests passed (%6.2f%%)\n", $ntested>0 ? 100*($N_OK/$ntested) : 0);
}
#test_all(0);
#exit(0);

##---------------------------------------------------------------------
## test data
sub test_data_0 {
  our $p = pdl(double, [
			[10,0,0,0,-2,0],
			[3,9,0,0,0,3],
			[0,7,8,7,0,0],
			[3,0,8,7,5,0],
			[0,8,0,9,9,13],
			[0,4,0,0,2,-1],
		       ]);
  our $a = $p;
}

sub test_data_1 {
  our $p = pdl(double, [
			[10,0,0,0,-2],
			[3,9,0,0,0],
			[0,7,8,7,0],
			[3,0,8,7,5],
			[0,8,0,9,9],
			[0,4,0,0,2],
		       ]);
  our $a = $p;
}

sub test_data_2 {
  our $a = pdl([[1,0,0,2],
		[0,3,4,0],
		[5,0,0,6]]);
}
sub test_data_bg {
  our $bgdata      = Storable::retrieve("bgdata.bin");
  our ($bgw,$bgnz) = @$bgdata{qw(which nzvals)};
  our $bgccs       = PDL::CCS::Nd->newFromWhich($bgw,$bgnz);
  our $bgsyms      = retrieve("bgsyms.bin");
}
#test_data_bg();

##---------------------------------------------------------------------
## bg data: blocksize twiddling
sub test_bg_blocksize {
  test_data_bg;

  $bgd  = $bgccs->double;
  $bgd  = $bgccs;
  $ugd  = $bgd->sumover;
  $bgd /= $ugd; ##-- kira: "Debugger floating point exception" if $bgccs still has an integer type

  $bg0 = $bgccs->copy;
  $bg0->[$VALS] = $bg0->[$VALS]->double;

  cmpthese(4,
	   {
	    'bs=0'    => sub { $PDL::CCS::Nd::BINOP_BLOCKSIZE_MAX=0;     $bg1 = $bg0*$bg0; },
	    'bs=65535'=> sub { $PDL::CCS::Nd::BINOP_BLOCKSIZE_MAX=65535; $bg1 = $bg0*$bg0; },
	    'bs=8192' => sub { $PDL::CCS::Nd::BINOP_BLOCKSIZE_MAX=8192;  $bg1 = $bg0*$bg0; },
	    'bs=4096' => sub { $PDL::CCS::Nd::BINOP_BLOCKSIZE_MAX=1024;  $bg1 = $bg0*$bg0; },
	   });
}
#test_bg_blocksize;

##---------------------------------------------------------------------
## Bigrams: test number of neighbor-links
sub align_neighbors { ##-- like threaded matmult with 'or' instead of '+' for inner() : BROKEN!
  my ($a,$align_dim,$infix) = @_;

  $infix      = '' if (!defined($infix));
  my $outfile = "align_neighbors${infix}.out";
  open(OUT,">$outfile") or die("$0: open failed for '$outfile': $!");
  #print OUT "$,";
  my $b = $a;

  my $nrdims     = 1;
  my $rdpdl      = pdl(long, [ [$align_dim,$align_dim] ]);
  my $report_dim = $align_dim==0 ? 1 : 0;

  ##-- Get & sort relevant indices, vals
  my $ixa    = $a->[$WHICH];
  my $avals  = $a->[$VALS];
  my $nixa   = $ixa->dim(1);
  my ($ixar,$avalsr);
  if ($align_dim==0) {
    ##-- a: relevant dims are a prefix of physical dims, e.g. pre-sorted
    $ixar   = $ixa->slice("0,");
    $avalsr = $avals;
  } else {
    $ixar          = $ixa->slice("$align_dim,");
    my $ixar_sorti = $ixar->qsortveci;
    $ixa           = $ixa->dice_axis(1,$ixar_sorti);
    $ixar          = $ixar->dice_axis(1,$ixar_sorti);
    $avalsr        = $avals->index($ixar_sorti);
  }
  ##
  my $ixb   = $b->[$WHICH];
  my $bvals = $b->[$VALS];
  my $nixb  = $ixb->dim(1);
  my $rb    = $rdpdl->slice("(1)");
  my ($ixbr,$bvalsr);
  if ($align_dim==0) {
    ##-- b: relevant dims are a prefix of physical dims, e.g. pre-sorted
    $ixbr   = $ixb->slice("0,");
    $bvalsr = $bvals;
  } else {
    $ixbr          = $ixb->slice("$align_dim,");
    my $ixbr_sorti = $ixbr->qsortveci;
    $ixb           = $ixb->dice_axis(1,$ixbr_sorti);
    $ixbr          = $ixbr->dice_axis(1,$ixbr_sorti);
    $bvalsr        = $bvals->index($ixbr_sorti);
  }

  ##-- initialize: state vars
  my $blksz  = $nixa > $nixb ? $nixa : $nixb;
  $blksz     = $BINOP_BLOCKSIZE_MIN if ($BINOP_BLOCKSIZE_MIN && $blksz < $BINOP_BLOCKSIZE_MIN);
  $blksz     = $BINOP_BLOCKSIZE_MAX if ($BINOP_BLOCKSIZE_MAX && $blksz > $BINOP_BLOCKSIZE_MAX);
  my $istate = PDL->zeroes($P_LONG,7); ##-- [ nnzai,nnzai_nxt, nnzbi,nnzbi_nxt, nnzci,nnzci_nxt, cmpval ]
  my $ostate = $istate->pdl;

  ##-- initialize: output vectors
  my $nzai    = PDL->zeroes($P_LONG,   $blksz);
  my $nzbi    = PDL->zeroes($P_LONG,   $blksz);
  my $ixc_cur = undef; #PDL->zeroes($P_LONG, 2,$blksz);
  my $ixc_out = PDL->null->long;
  my $nnzc   = 0;
  my $zc     = 0;

  ##-- block-wise variables
  ##   + there are way too many of these...
  my ($nzai_prv,$nzai_pnx, $nzbi_prv,$nzbi_pnx, $nzci_prv,$nzci_pnx,$cmpval_prv);
  my ($nzai_cur,$nzai_nxt, $nzbi_cur,$nzbi_nxt, $nzci_cur,$nzci_nxt,$cmpval);
  my ($nzci_max, $blk_slice, $nnzc_blk,$nnzc_slice_blk);
  my ($nzai_blk,$nzbi_blk,$ixa_blk,$ixb_blk,$ixc_blk,$nzc_blk,$cimask_blk,$ciwhich_blk);

  my ($ixc_e,$ixc_f);
  my $nnzc_prev=0;
  do {
    ##-- align a block of data
    ccs_binop_align_block_mia($ixar,$ixbr,$istate, $nzai,$nzbi,$ostate);

    ##-- parse current alignment algorithm state
    ($nzai_prv,$nzai_pnx, $nzbi_prv,$nzbi_pnx, $nzci_prv,$nzci_pnx,$cmpval_prv) = $istate->list;
    ($nzai_cur,$nzai_nxt, $nzbi_cur,$nzbi_nxt, $nzci_cur,$nzci_nxt,$cmpval)     = $ostate->list;
    $nzci_max = $nzci_cur-1;

    if ($nzci_max >= 0) {
      ##-- construct block output pdls: nzvals
      $blk_slice = "${nzci_prv}:${nzci_max}";
      $nzai_blk  = $nzai->slice($blk_slice);
      $nzbi_blk  = $nzbi->slice($blk_slice);

      ##-- construct block output pdls: ixc_cur
      $ixa_blk = $ixa->slice("($report_dim),")->index($nzai_blk);
      $ixb_blk = $ixb->slice("($report_dim),")->index($nzbi_blk);

      ##-- get minimum indices in $ixa_blk, $ixb_blk
      ($a_b,$b_a) = which_both($ixa_blk < $ixb_blk);
      $ixc_cur    = zeroes(long, 2,$ixa_blk->dim(0));
      $ixc_cur_ab = $ixc_cur->dice_axis(1,$a_b);
      $ixc_cur_ab->slice("(0),") .= $ixa_blk->index($a_b);
      $ixc_cur_ab->slice("(1),") .= $ixb_blk->index($a_b);
      undef($ixc_cur_ab);
      $ixc_cur_ba = $ixc_cur->dice_axis(1,$b_a);
      $ixc_cur_ba->slice("(0),") .= $ixb_blk->index($b_a);
      $ixc_cur_ba->slice("(1),") .= $ixa_blk->index($b_a);
      undef($ixc_cur_ab);

      ##-- make unique: ixc_cur
      ##$ixc_cur = $ixa_blk->slice("$align_dim2")->append($ixb_blk->slice("$align_dim1"));
      ($ixc_f,$ixc_e) = $ixc_cur->vv_qsortvec->rlevec;
      $ixc_e          = $ixc_e->dice_axis(1,$ixc_f->which);
      $ixc_cur        = $ixc_e;
      #@ixc = $ixc_e->flat->list;
      #print OUT map { ($ixc[2*$_],"\t",$ixc[2*$_+1],"\n") } (0..($ixc_e->dim(1)-1));
      ##--
      if (0) { ##-- debug: dump to text file
	print "<DEBUG> dumping block of size=".(1+$nzci_max-$nzci_prv). " to '$outfile' ... ";
	ccs_dump_which($ixc_cur, \*OUT, "%d", " ", "\n");
	print "dumped.\n";
      } else {
	##-- show what we've aligned
	print "<DEBUG> aligned block of size ".sprintf("%8d\n", 1+$nzci_max-$nzci_prv);
      }
      undef($ixc_f);
      undef($ixc_e);
    }

    ##-- add to final output pdl
    if (1) {
      if (!defined($ixc_out) || $ixc_out->isempty) {
	$ixc_out = $ixc_cur->ushort if (defined($ixc_cur));
      } else {
	my $ixc_dim1_old = $ixc_out->dim(1);
	$ixc_out->reshape($ixc_out->dim(0), $ixc_out->dim(1)+$ixc_cur->dim(1));
	$ixc_out->slice(",${ixc_dim1_old}:-1") .= $ixc_cur;
      }
    }
    undef($ixc_cur);

    ##-- possibly allocate for another block
    if ($nzai_cur < $nixa || $nzbi_cur < $nixb) {
      $nzci_nxt -= $nzci_cur;
      $nzci_cur  = 0;

      if ($nzci_nxt > $nzai->dim(0)) {
	$nzai    = $nzai->reshape($nzci_nxt+$blksz);
	$nzbi    = $nzbi->reshape($nzci_nxt+$blksz);
      }

      $istate .= $ostate;
      $istate->set(4, $nzci_cur);
      $istate->set(5, $nzci_nxt);
    }
    $nnzc_prev = $nnzc;
  } while ($nzai_cur < $nixa || $nzbi_cur < $nixb);

  ##-- all done
  close(OUT);

  return $ixc_out;
}

sub test_bg_npairs {
  test_data_bg;

  ##-- ARGH: bug in CCS::Nd->whichND() for all-physical but shuffled dims
  #  our $x = sequence(4,2) * (sequence(4,2)%3);
  #  our $xs = $x->toccs;
  #  isok("sumover", all($x->sumover==$xs->sumover->decode));
  #  $xst  = $xs->xchg(0,1);
  #  $xsts = $xst->sumover;
  #  isok("xchg(0,1)->sumover", all($x->xchg(0,1)->sumover==$xsts->decode));

  our $bg  = $bgccs->clone;

  ##-- get total number of words (& unigram freqs)
  our $uga     = $bgccs->sumover->todense;
  our $nwordsa = $uga->dim(0);

  ##-- hack: elimate hapax BIGRAMS
  $bg->[$VALS]->where($bg->[$VALS]<=1) .= 0;
  $bg = $bg->recode;
  #$bg->[$VALS] .= 1;
  #$bg->missing(0);

  ##-- get pruned unigram freqs (hack: (L+R)/2)
  our $ugp     = ($bg->sumover->todense + $bg->xchg(0,1)->sumover->todense)->double / 2;
  our $nwordsp = $ugp->nnz->sclr;

  ##-- get actual number of words remaining in the vocabulary
  our $words     = $bg->[$WHICH]->slice("(0),")->append($bg->[$WHICH]->slice("(1),"))->flat->qsort;
  our ($wf,$wid) = rle($words);
  $wid = $wid->where($wf);
  $wf  = $wf->where($wf);
  $wid->sever;
  $wf->sever;
  our $nwordsp1 = $wid->nelem;

  if (0) {
    ##-- align on l(r<->r)l , r(l<->l)r: hack: ???
    our $rr = align_neighbors($bg, 0, ".rr");
    our $ll = align_neighbors($bg, 1, ".ll");

    ##-- sort 'em & take union (find all pairs)
    $ll = $ll->vv_qsortvec;
    $rr = $rr->vv_qsortvec;
    $pairs = vv_union($ll,$rr);

    ##-- ... and make pairs unique
    ($pair_f,$pair_e) = $pairs->rlevec;
    $pairs = $pair_e->dice_axis(1,$pair_f->which);
    undef($pair_f);
    undef($pair_e);
    $pairs->sever;

    ##-- save some stuff
    store($ll,'align.ll.bin');
    store($rr,'align.rr.bin');
    store($pairs,'align.pairs.bin');
  } else {
    ##-- just load
    $ll = retrieve 'align.ll.bin';
    $rr = retrieve 'align.rr.bin';
    $pairs = retrieve 'align.pairs.bin';
  }


  ##-- get some stats
  our $npairs  = $pairs->dim(1);
  our $nppairsp = $nwordsp**2;
  our $nppairsa = $nwordsa**2;
  printf("Number of words (pruned/all)    : %10d / %10d (%6.2f%%)\n", $nwordsp, $nwordsa, 100*$nwordsp/$nwordsa);
  printf("Number of pairs (found/pruned^2): %10d / %10d (%6.2f%%)\n", $npairs, $nppairsp, 100.0*$npairs/$nppairsp);
  printf("Number of pairs (found/all^2)   : %10d / %10d (%6.2f%%)\n", $npairs, $nppairsa, 100.0*$npairs/$nppairsa);

  ##--prints
  #Number of words (pruned/all)    :       9681 /      45822 ( 21.13%)
  #Number of pairs (found/pruned^2):    6992524 /   93721761 (  7.46%)
  #Number of pairs (found/all^2)   :    6992524 / 2099655684 (  0.33%)


  ##-- get more word stats
  $ugpp = $ugp->index($wid);
  $ugap = $uga->index($wid);
  printf("Avg word freq (pruned/all)      : %6.2f / %6.2f\n", $ugpp->davg, $uga->davg);
  printf(" + pruned words only            : %6.2f / %6.2f\n", $ugpp->davg, $ugap->davg);
  printf("Log-avg word freq (pruned/all)  : %6.2f / %6.2f\n", $ugpp->log->daverage->exp, $uga->log->average->exp);
  printf(" + pruned words only            : %6.2f / %6.2f\n", $ugpp->log->daverage->exp, $ugap->log->average->exp);

  ##-- prints
  #Avg word freq (pruned/all)      :  18.72 /   7.38
  # + pruned words only            :  18.72 /  29.93
  #Log-avg word freq (pruned/all)  :   2.73 /   1.73
  # + pruned words only            :   2.73 /   6.68

  ##-- look at it
  #points($ugpp->log,$ugap->log); ##-- ok: looks pretty (log-) linear

  print "test_bg_npairs: done.\n";
}
#test_bg_npairs;

##---------------------------------------------------------------------
## CCS::Nd: revised alignment (signature-sensitive)

BEGIN {
  *_ccsnd_parse_signature = \&PDL::CCS::Nd::_ccsnd_parse_signature;
  *_ccsnd_align_dims      = \&PDL::CCS::Nd::_ccsnd_align_dims;
}
sub test_align_v2 {
  our $parsed1 = _ccsnd_parse_signature('int which(Ndims,Nnz); vals(Nnz); [o]out(Nnz);');

  our $a  = zeroes(4,3)->toccs;
  our $br  = ones(4)->toccs;
  our $brv = ones(4)->toccs->dummy(1);
  our $bc  = ones(1,3)->toccs;
  our $bcv = ones(3)->toccs->dummy(0);

  our $sig = _ccsnd_parse_signature('a(); b(); [o]c();');
  our ($dims);
  $dims = _ccsnd_align_dims($sig, [$a,$bcv]);

  ##-- try it with inner()-like

  print "test_align_v2: done.\n";
}
test_align_v2();



##---------------------------------------------------------------------
## CCS: Nd: matrix multiplication, on bigrams
sub test_matmult_bg {
  test_data_bg;

  $bgd = $bgccs->double;
  $bgd->bad_is_missing(1);
  $bgd->nan_is_missing(1);

  our ($NC,$NW) = (32, $bgd->dim(1));
  $phat = random($NC,$NW)->abs;
  $phat->maximum_n_ind($phatmaxi = zeroes(long,4,$phat->dim(1)));
  $phatmask = zeroes(byte,$phat->dims);
  $phatmask->index($phatmaxi->xchg(0,1)) .= 1;
  $phat *= $phatmask;
  $phat /= $phat->sumover->slice("*1,");

  our $phats = $phat->toccs;

  our $bgdr = $bgd->matmult($phats);

  print "bg[NW=$NW,NW] x phat[NC=$NC,NW] --> pdl(", join(',', $bgdr->dims), ")\n";
}
#test_matmult_bg();


##---------------------------------------------------------------------
## CCS: Nd: matrix multiplication (low-level twiddling stuff)
sub test_matmult_pre1 {
  our ($a,$b,$c,$in, $as,$bs,$cs,$ins);

  $a=sequence(4);
  $b=($a+1);
  $as=$a->toccs;
  $bs=$b->toccs;

  ##-- inner(row_vector,row_vector)
  $in  = $a->inner($b);
  $ins = $as->inner($bs);
  $c  = ($a*$b)->sumover;
  $cs = ($as * $bs);
  $cs = $cs->sumover;
  isok("matmult_pre1:(row,row):inner:sparse~dense",  all($ins->todense==$in));
  isok("matmult_pre1:(row,row):inner~sumover(*)", all($c==$in));
  isok("matmult_pre1:(row,row):sparse~dense",     all($cs->todense==$c));

  ##-- inner(matrix,row_vector)
  our ($N,$M,$P) = (4,3,2);
  $a = sequence($N,$M);
  $b = sequence($N)+1;
  $as=$a->toccs;
  $bs=$b->toccs;

  $in = $a->inner($b);
  $ins = $as->inner($bs);
  $c  = ($a*$b)->sumover;
  $cs = ($as*$bs);
  $cs = $cs->sumover;
  isok("matmult_pre1:(mat,row):inner:sparse~dense",  all($ins->todense==$in));
  isok("matmult_pre1:(mat,row):inner~sumover(*)", all($c==$in));
  isok("matmult_pre1:(mat,row):sparse~dense",     all($cs->todense==$c));

  ##-- inner(matrix,col_vector)
  $b = sequence(1,$M)+1;
  $bs=$b->toccs;

  $in = $a->inner($b);
  $ins = $as->inner($bs);
  $c  = ($a*$b)->sumover;
  $cs = ($as*$bs);
  $cs = $cs->sumover;
  isok("matmult_pre1:(mat,col):inner:sparse~dense",  all($ins->todense==$in));
  isok("matmult_pre1:(mat,col):inner~sumover(*)", all($c==$in));
  isok("matmult_pre1:(mat,col):sparse~dense",     all($cs->todense==$c));

  ##-- (debug)
  sub argdims { "a(".dimstr($a).") x b(".dimstr($b).") -> c(".dimstr($c).")\n"; };
  sub dimstr { my $p=shift; return join(',', map {$_==$N?'N':($_==$M?'M':($_==$P?'P':$_))} $p->dims); }

  ##-- matmult: (mat,col) --> col : a(N,M) x b(1,N) --> c(1,M)
  $b = sequence(1,$N)+1;
  $bs = $b->toccs;
  $c  = $a->matmult($b,null);
  $cs = $as->matmult($bs);
  isok("matmult_pre1:(mat,col)->col", all($cs->todense==$c));

  ##-- matmult: (row,mat) --> row : b(M) x a(N,M) -> c(N,1)
  $b = sequence($M)+1;
  $bs = $b->toccs;
  $c  = $b->matmult($a,null);
  $cs = $bs->matmult($as);
  isok("matmult_pre1:(row,mat)->row", all($cs->todense==$c));

  ##-- matmult: (mat,sclr) --> mat : a(N,M) x b() -> c(N,M)
  $b  = pdl(42);
  $c  = $a->matmult($b,null);
  $cs = $as->matmult($b);
  isok("matmult_pre1:(mat,sclr)->mat", all($cs->todense==$c));

  ##-- matmult: (row,col) --> sclr : a(N) x b(1,N) -> c(1,1)
  $a  = sequence($N);
  $b  = sequence(1,$N)+1;
  $as = $a->toccs;
  $bs = $b->toccs;
  $c  = $a->matmult($b,null);
  $cs = $as->matmult($bs);
  isok("matmult_pre1:(row,col)->mat", all($cs->todense==$c));

  ##-- matmult: (col,row) --> mat : b(1,N) x a(N) -> c(N,N)
  $c  = $b->matmult($a,null);
  $cs = $bs->matmult($as);
  isok("matmult_pre1:(col,row)->mat", all($cs->todense==$c));

  ##-- matmult: (mat,mat) --> mat : a(N,M) x b(P,N) -> c(P,M)
  $a = sequence($N,$M);
  $b = sequence($P,$N)+1;
  $as = $a->toccs;
  $bs = $b->toccs;
  $c  = $a->matmult($b,null);
  $cs = $as->matmult($bs);
  isok("matmult_pre1:(mat,mat)->mat", all($cs->todense==$c));

  print "test_matmult_pre1: done.\n";
}
test_matmult_pre1();


##---------------------------------------------------------------------
## CCS: Nd: virtual dims

sub test_vdims_1 {
  test_data_1;
  $as = $a->toccs;

  $asd = $as->decode;
  isok("toccs/decode", all($asd==$a));

  ##-- test: whichND
  $aw  = $a->whichND->vv_qsortvec;
  $asw = $as->whichND;
  isok("ccs:nd:whichND:literal",all($aw==$asw));

  ##-- test: whichND: dummy(0,1)
  $a1  = $a->dummy(0,1);
  $as1 = $as->dummy(0,1);
  $a1w = $a1->whichND;
  $as1w= $as1->whichND;
  $as1v= $as1->whichVals;
  isok("ccs:nd:whichND:+dummy(0,1)",   all($a1w->vv_qsortvec==$as1w->vv_qsortvec));
  isok("ccs:nd:whichVals:+dummy(0,1)", all($as1v==$a1->indexND($as1w)));

  ##-- test: whichND: dummy(0,2)
  $a1  = $a->dummy(0,2);
  $as1 = $as->dummy(0,2);
  $a1w = $a1->whichND;
  $as1w= $as1->whichND;
  $as1v= $as1->whichVals;
  isok("ccs:nd:whichND:+dummy(0,2)", all($a1w->vv_qsortvec==$as1w->vv_qsortvec));
  isok("ccs:nd:whichVals:+dummy(0,2)", all($as1v==$a1->indexND($as1w)));

  ##-- test: whichND: dummy(0,2)+dummy(0,3)
  $a1  = $a->dummy(0,2)->dummy(0,3);
  $as1 = $as->dummy(0,2)->dummy(0,3);
  $a1w = $a1->whichND;
  $as1w= $as1->whichND;
  $as1v= $as1->whichVals;
  isok("ccs:nd:whichND:+dummy(0,2)+dummy(0,3)", all($a1w->vv_qsortvec==$as1w->vv_qsortvec));
  isok("ccs:nd:whichND:+dummy(0,2)+dummy(0,3)", all($as1v==$a1->indexND($as1w)));


  ##-- test: dummy(0,1)
  $a1  = $a->dummy(0,1);
  $as1 = $as->dummy(0,1);
  $as1d= $as1->decode;
  isok("ccs:nd:dummy(0,1)", all($as1d==$a1));

  ##-- test: dummy(0,3)
  $a1 = $a->dummy(0,3);
  $as1 = $as->dummy(0,3);
  $as1d = $as1->decode;
  isok("ccs:nd:dummy(0,3)", all($as1d==$a1));

  ##-- test: dummy(1,3)
  $a1 = $a->dummy(1,3);
  $as1 = $as->dummy(1,3);
  $as1d = $as1->decode;
  isok("ccs:nd:dummy(1,3)", all($as1d==$a1));

  ##-- test: dummy(2,3)
  $a1 = $a->dummy(2,3);
  $as1 = $as->dummy(2,3);
  $as1d = $as1->decode;
  isok("ccs:nd:dummy(2,3)", all($as1d==$a1));

  ##-- test: dummy(-1,3)
  $a1 = $a->dummy(-1,3);
  $as1 = $as->dummy(-1,3);
  $as1d = $as1->decode;
  isok("ccs:nd:dummy(-1,3)", all($as1d==$a1));
}
test_vdims_1;

##---------------------------------------------------------------------
## CCS: Nd: recode

sub test_recode {
  test_data_1;

  ##-- test: assignment
  $as  = $a->toccs;
  $as .= 0;
  isok("assign(as.=0)", all matchpdl($as->decode,$a->zeroes));
  $as .= $a;
  isok("assign(as.=a)", all matchpdl($as->decode,$a));

  ##-- test: ++
  $as++; $a++;
  isok("++", all matchpdl($as->decode,$a));
  $as--; $a--;

  ##-- test: !
  $nas = !$as;
  $na  = !$a;
  isok("!", all matchpdl($nas->decode,$na));

  ##-- test: ==
  $aseq = ($as==$a);
  isok("==", $aseq->all);

  ##-- test: ccs-encode singleton
  our $z0  = pdl(0);
  our $z1  = pdl(1);
  our $z0s = $z0->toccs;
  our $z1s = $z1->toccs;
  isok("ccs:nd:encode:singleton:0",      $z0==$z0s->decode);
  isok("ccs:nd:encode:singleton:1",      $z1==$z1s->decode);
  isok("ccs:nd:encode:singleton:0:bool", (($z0 && $z0s) || (!$z0 && !$z0s)));
  isok("ccs:nd:encode:singleton:1:bool", (($z1 && $z1s) || (!$z1 && !$z1s)));

  ##-- test recode BAD->0
  $as = $a->toccs($BAD);
  $as->missing(0);
  $as->recode;
  isok("ccs:nd:recode(BAD->0)", all matchpdl($as->decode,$a));

  ##-- test recode: empty pdl /0
  $as->[$VALS] .= 0;
  $as->recode();

  $asd = $as->decode;
  isok("ccs:nd:recode(empty:z=0)", all matchpdl($asd,$a->zeroes));

  ##-- test recode: empty/bad
  our $abad = $a->zeroes->setvaltobad(0);
  $as = $abad->toccs($BAD);
  isok("ccs:nd:recode(empty:z=bad)", all matchpdl($as->decode,$abad));

  ##-- test type conversion
  isok("ccs:nd:long()", PDL::CCS::Nd::long()==PDL::long());
  $as     = $a->toccs;
  $aslong = $as->long;
  isok("ccs:nd->long()", $aslong->type == PDL::long());
}
test_recode();

##---------------------------------------------------------------------
## CCS: Nd: binary ops: scalar (correct)

sub test_nd_binop_sclr {
  our ($op_name,$b,$missing,$swap,$flags) = @_;

  ##-- get params
  $op_name = 'add' if (!defined($op_name));
  our $pdl_op = PDL->can($op_name);
  our $ccs_op = PDL::CCS::Nd->can($op_name);
  die("test_nd_binop_sclr(): no PDL op named '$op_name'!") if (!defined($pdl_op));
  die("test_nd_binop_sclr(): no CCS op named '$op_name'!") if (!defined($ccs_op) || $ccs_op eq $pdl_op);

  $missing = 0 if (!defined($missing));
  $missing = pdl($missing);

  $swap = 0 if (!defined($swap));

  ##-- get data
  test_data_1;
  $a = $a->setvaltobad(0);
  $a = $a->setbadtoval($missing) if ($missing->isgood);

  our $as = $a->toccs($missing,$flags);

  ##-- set blocksize?
  #$PDL::CCS::Nd::BINOP_BLOCKSIZE_MAX = 2;

  ##-- guts
  our $cs = $ccs_op->($as,          $b,  $swap);
  our $c  = $pdl_op->($a,  todense($b),  $swap);

  isok("ccs:nd:binop=${op_name},z=${missing},b=sclr($b),swap=$swap:type", $cs->type==$c->type);
  isok("ccs:nd:binop=${op_name},z=${missing},b=sclr($b),swap=$swap:vals", all(matchpdl($cs->decode,$c)));
}
##-- test_nd_binop_sclr(op,b,missing,swap,flags)
#test_nd_binop_sclr('mult',$BAD,0);
test_nd_binop_sclr('mult',42, 0,0,0);
#test_nd_binop_sclr('mult',pdl([[42]])->toccs, 0,0,0);
#test_nd_binop_sclr('divide',0, $BAD,0,0);
#test_nd_binop_sclr('power',0, $BAD,1,0);

sub test_nd_binop_sclr_all {
  my ($sclr,$binop,$missing,$swap);
  foreach $missing ($BAD) {
    foreach $swap (0,1) {
      foreach $binop (
		      qw(plus minus mult divide modulo power),
		      qw(gt ge lt le eq ne spaceship),
		      qw(and2 or2 xor),
		      qw(shiftleft shiftright), #)
		     ) {
	foreach $sclr (0,42) {
	  test_nd_binop_sclr($binop,$sclr,$missing,$swap);
	}
      }
    }
  }
}
#test_nd_binop_sclr_all();


##---------------------------------------------------------------------
## CCS: Nd: binary ops: missing-annihilator: col- & row-vectors

sub test_nd_binop_cvrv_mia {
  my ($op_name,$missing,$swap, $flags) = @_;
  ##-- get params
  $op_name = 'add' if (!defined($op_name));
  our $pdl_op = PDL->can($op_name);
  our $ccs_op = PDL::CCS::Nd->can($op_name);
  die("test_nd_binop(): no PDL op named '$op_name'!") if (!defined($pdl_op));
  die("test_nd_binop(): no CCS op named '$op_name'!") if (!defined($ccs_op) || $ccs_op eq $pdl_op);

  $missing = 0 if (!defined($missing));
  $missing = pdl($missing);

  $swap = 0 if (!defined($swap));

  ##-- get data
  test_data_1;
  $a = $a->setvaltobad(0);
  $a = $a->setbadtoval($missing) if ($missing->isgood);

  our $b0 = sequence(  $a->dim(0))+1;
  our $b1 = sequence(1,$a->dim(1))+1;

  our $as   = $a->toccs($missing,$flags);
  our $bs0  = $b0->toccs($missing,$flags);
  our $bs1  = $b1->toccs($missing,$flags);
  our $bs1v = $b1->flat->toccs($missing,$flags)->dummy(0,1);

  ##-- set blocksize?
  #$PDL::CCS::Nd::BINOP_BLOCKSIZE_MAX = 2;

  ##-- guts
  our $cs0 = $ccs_op->($as,$bs0, $swap);
  our $c0  = $pdl_op->($a, $b0,  $swap);
  isok("ccs:nd:arg2=rv,binop=${op_name},z=${missing},swap=$swap:type", $cs0->type==$c0->type);
  isok("ccs:nd:arg2=rv,binop=${op_name},z=${missing},swap=$swap:vals", all(matchpdl($cs0->decode,$c0)));

  our $cs1  = $ccs_op->($as,$bs1, $swap);
  our $c1   = $pdl_op->($a, $b1,  $swap);
  isok("ccs:nd:arg2=cv,binop=${op_name},z=${missing},swap=$swap:type", $cs1->type==$c1->type);
  isok("ccs:nd:arg2=cv,binop=${op_name},z=${missing},swap=$swap:vals", all(matchpdl($cs1->decode,$c1)));

  our $cs1v = $ccs_op->($as,$bs1v, $swap);
  isok("ccs:nd:arg2=cv(ccs-virtual),binop=${op_name},z=${missing},swap=$swap:type", $cs1v->type==$c1->type);
  isok("ccs:nd:arg2=cv(ccs-virtual),binop=${op_name},z=${missing},swap=$swap:vals", all(matchpdl($cs1v->decode,$c1)));

  our $csrc = $ccs_op->($bs0,$bs1v, $swap);
  our $crc  = $pdl_op->($b0, $b1,   $swap);
  isok("ccs:nd:arg1=rv,arg2=cv(ccs-virtual),binop=${op_name},z=${missing},swap=$swap:type",
       $csrc->type==$crc->type);
  isok("ccs:nd:arg1=rv,arg2=cv(ccs-virtual),binop=${op_name},z=${missing},swap=$swap:vals",
       all(matchpdl($csrc->decode,$crc)));
}
#test_nd_binop_cvrv_mia('plus',$BAD,0);
#test_nd_binop_cvrv_mia('mult',$BAD,0);
test_nd_binop_cvrv_mia('mult',0,0);
#test_nd_binop_cvrv_mia('divide',0,0, 3);

sub test_nd_binop_cvrv_all {
  my ($binop,$missing,$swap);
  foreach $missing ($BAD) {
    foreach $swap (0,1) {
      foreach $binop (
		      qw(plus minus mult divide modulo power),
		      qw(gt ge lt le eq ne spaceship),
		      qw(and2 or2 xor),
		      qw(shiftleft shiftright),
		     )
	{
	  test_nd_binop_cvrv_mia($binop,$missing,$swap);
	}
    }
  }
}
#test_nd_binop_cvrv_all();


##---------------------------------------------------------------------
## CCS: Nd: binary ops: missing-annihilator: full dim-match

sub test_nd_binop_mia {
  my ($op_name,$missing,$swap) = @_;
  ##-- get params
  $op_name = 'add' if (!defined($op_name));
  our $pdl_op = PDL->can($op_name);
  our $ccs_op = PDL::CCS::Nd->can($op_name);
  die("test_nd_binop(): no PDL op named '$op_name'!") if (!defined($pdl_op));
  die("test_nd_binop(): no CCS op named '$op_name'!") if (!defined($ccs_op) || $ccs_op eq $pdl_op);

  $missing = 0 if (!defined($missing));
  $missing = pdl($missing);

  $swap = 0 if (!defined($swap));

  ##-- get data
  test_data_1;
  $a = $a->setvaltobad(0);
  $a = $a->setbadtoval($missing) if ($missing->isgood);
  $a = $a->abs; ##-- use absval so bitwise ops don't overflow

  our $b = $a->rotate(1);
  $b->sever;

  our $as = $a->toccs($missing);
  our $bs = $b->toccs($missing);

  ##-- guts
  our $cs = $ccs_op->($as,$bs, $swap);
  our $c  = $pdl_op->($a, $b,  $swap);

  if (ref($c)) {
    isok("ccs:nd:binop=${op_name},missing=${missing},swap=$swap:type", $cs->type==$c->type);
    isok("ccs:nd:binop=${op_name},missing=${missing},swap=$swap:vals", all(matchpdl($cs->decode,$c)));
  } else {
    isok("ccs:nd:binop=${op_name},missing=${missing},swap=$swap", $cs eq $c);
  }
}
#test_nd_binop_mia('mult',$BAD,0);
#test_nd_binop_mia('and2',$BAD,0);
#test_nd_binop_mia('or2',$BAD,0);
#test_nd_binop_mia('xor',$BAD,0);
test_nd_binop_mia('shiftleft',$BAD,0);

sub test_nd_binop_mia_all {
  my ($binop,$missing,$swap);
  foreach $missing ($BAD) {
    foreach $swap (0,1) {
      foreach $binop (
		      qw(plus minus mult divide modulo power),
		      qw(gt ge lt le eq ne spaceship),
		      qw(and2 or2 xor),
		      qw(shiftleft shiftright),
		     )
	{
	  test_nd_binop_mia($binop,$missing,$swap);
	}
    }
  }
}
#test_nd_binop_mia_all();


##---------------------------------------------------------------------
## CCS: Nd: unary ops

sub test_nd_unop {
  our ($op_name,$missing) = @_;

  ##-- get params
  $op_name = 'sumover' if (!defined($op_name));
  our $ccs_op = PDL::CCS::Nd->can($op_name);
  our $pdl_op = PDL->can($op_name);
  die("test_nd_unop(): no PDL op named '$op_name'!") if (!defined($pdl_op));
  die("test_nd_unop(): no CCS op named '$op_name'!") if (!defined($ccs_op));

  $missing = 0 if (!defined($missing));
  $missing = pdl($missing);

  ##-- get data
  test_data_1;
  $a = $a->setvaltobad(0);
  $a = $a->setbadtoval($missing) if ($missing->isgood);
  our $as = $a->toccs($missing);

  ##-- guts
  our $bs = $ccs_op->($as);
  our $b  = $pdl_op->($a);

  if (ref($b)) {
    isok("ccs:nd:unop=${op_name},missing=${missing}:type", $bs->type==$b->type);
    isok("ccs:nd:unop=${op_name},missing=${missing}:vals", all(matchpdl($bs->decode,$b)));
  } else {
    isok("ccs:nd:unop=${op_name},missing=${missing}", $bs eq $b);
  }
}
#test_nd_unop('abs',0);
#test_nd_unop('bitnot',-1);
test_nd_unop('bitnot',$BAD);

sub test_nd_unop_all {
  my ($unop,$missing);
  foreach $missing (0,-1,$BAD) {
    foreach $unop (qw(bitnot sqrt abs sin cos not exp log log10)) {
      test_nd_unop($unop,$missing);
    }
  }
}
#test_nd_unop_all();

##---------------------------------------------------------------------
## CCS: ufunc: Nd
sub test_nd_ufunc_1 {
  our ($ufunc_name,$missing) = @_;

  ##-- get params
  $ufunc_name = 'sumover' if (!defined($ufunc_name));
  our $ccs_ufunc = PDL::CCS::Nd->can($ufunc_name);
  our $pdl_ufunc = PDL->can($ufunc_name);
  die("test_nd_ufunc(): no PDL ufunc named '$ufunc_name'!") if (!defined($pdl_ufunc));
  die("test_nd_ufunc(): no CCS ufunc named '$ufunc_name'!") if (!defined($ccs_ufunc));

  $missing = 0 if (!defined($missing));
  $missing = pdl($missing);

  ##-- get data
  test_data_1;
  $a = $a->setvaltobad(0);
  $a = $a->setbadtoval($missing) if ($missing->isgood);
  our $ccs = $a->toccs($missing);

  our $ccs2 = $ccs_ufunc->($ccs);
  our $a2   = $pdl_ufunc->($a);

  if (ref($a2)) {
    isok("ccs:nd:ufunc=${ufunc_name},missing=${missing}:type", $ccs2->type==$a2->type);
    isok("ccs:nd:ufunc=${ufunc_name},missing=${missing}:vals", all(matchpdl($ccs2->decode,$a2)));
  } else {
    isok("ccs:nd:ufunc=${ufunc_name},missing=${missing}", $ccs2 eq $a2);
  }
}
#test_nd_ufunc_1('max',0);
#test_nd_ufunc_1('dprod',1);
#test_nd_ufunc_1('sumover',$BAD);
#test_nd_ufunc_1('ngoodover',$BAD);
test_nd_ufunc_1('nbadover',$BAD);

sub test_nd_ufunc_all {
  my ($ufunc,$missing);
  foreach $missing (0,1,$BAD) {
    foreach $ufunc (
		    qw(sumover dsumover prodover dprodover),
		    qw(andover orover bandover borover),
		    qw(maximum minimum),
		    qw(nnz nbadover ngoodover),
		    ##
		    ##-- scalars
		    qw(sum dsum prod dprod max min nbad ngood),
		   )
      {
	test_nd_ufunc_1($ufunc,$missing);
      }
  }
}
test_nd_ufunc_all();

##---------------------------------------------------------------------
## CCS: binops (block-wise alignment / missing-is-annihiliator): with column-vector
sub test_ccs_binops_mia_cv {
  our ($OPNAME,$swap) = @_;

  ##-- set operation sub
  our $OP = PDL->can($OPNAME);
  if (!defined($OP)) {
    warn("$0: no operation '$OPNAME' for PDL!");
    return undef;
  }
  $swap = 0 if (!defined($swap));

  ##-- get data
  test_data_1;
  $a = $a->setvaltobad(0); ##-- use BAD as missing value: ensure annihilator

  our $adims = [$a->dims];
  our $b = sequence($a->type,$a->dim(1))->rotate(1)->setvaltobad(0)->slice("*1,");

  our $ccsa = $a->toccs;
  our @a = our ($ixa,$nza,$za) = ($ccsa->whichND,$ccsa->nzvals,$ccsa->missing);

  our $ccsb = $b->flat->toccs;
  our @b = our ($ixb,$nzb,$zb) = ($ccsb->whichND,$ccsb->nzvals,$ccsb->missing);
  our (@align,$nzai,$nzbi);


  ##-- test: missing-annihilator
  our $nnzc   = ($nza->dim(0) > $nzb->dim(0) ? $nza->dim(0) : $nzb->dim(0));
  our $istate = zeroes(long,7); ##-- [ nnzai,nnzai_nxt, nnzbi,nnzbi_nxt, nnzci,nnzci_nxt, cmpval ]
  our $ostate = $istate->pdl;

  ##-- guts: alignment
  our $ixa1      = $ixa->slice("1,");
  our $ixa1sorti = qsortveci($ixa1);
  $ixa1          = $ixa1->dice_axis(1,$ixa1sorti);
  ccs_binop_align_block_mia($ixa1,$ixb,$istate, $nzai=zeroes(long,$nnzc),$nzbi=zeroes(long,$nnzc),$ostate);

  ##-- parse output state
  our ($nzai_cur,$nzai_nxt, $nzbi_cur,$nzbi_nxt, $nzci_cur,$nzci_nxt,$cmpval) = $ostate->list;

  ##-- trim input pdls
  our $nzci_max = $nzci_cur-1;
  $nzai = $nzai->slice("0:$nzci_max");
  $nzbi = $nzbi->slice("0:$nzci_max");

  ##-- construct output pdl: nzvals
  our $avals = $ccsa->vals;
  our $bvals = $ccsb->vals;
  our $nzc   = zeroes(($avals->type > $bvals->type ? $avals->type : $bvals->type), $nzci_cur);
  our $zc    = $OP->($avals->slice("-1"), $bvals->slice("-1"), $swap);
  $nzc      .= $OP->($avals->index($ixa1sorti)->index($nzai), $bvals->index($nzbi), $swap);
  ##
  ##-- get indices of "good" c() values
  our $cimask  = ($zc->isbad ? ($nzc->isgood) : ($nzc!=$zc));
  our $ciwhich = $cimask->which;
  $nzc = $nzc->index($ciwhich); #->append($zc);

  ##-- construct output pdl: which
  #our $ndims   = $ixa->dim(0);
  #our $cwhich  = zeroes(long, $ndims,$ciwhich->nelem);
  #$cwhich     .= $ixa->dice_axis(1,$nzai->index($ciwhich));
  our $cwhich   = $ixa->dice_axis(1,$ixa1sorti)->dice_axis(1,$nzai->index($ciwhich));

  ##-- re-sort output pdls
  our $cwhichsorti = $cwhich->qsortveci;
  $cwhich          = $cwhich->dice_axis(1,$cwhichsorti);
  $nzc             = $nzc->index($cwhichsorti);

  isok("ccs_align_block_mia:OP=$OPNAME,SWAP=$swap", all matchpdl(ccs_decode($cwhich,$nzc,$zc,[$a->dims]), $OP->($a,$b,$swap)));

  my $foo='bar';
}
sub test_ccs_binops_mia_cv_all {
  test_ccs_binops_mia_cv('plus');
  test_ccs_binops_mia_cv('minus');
  test_ccs_binops_mia_cv('mult');
  test_ccs_binops_mia_cv('divide');
  test_ccs_binops_mia_cv('gt');
  test_ccs_binops_mia_cv('ge');
  test_ccs_binops_mia_cv('lt');
  test_ccs_binops_mia_cv('le');
  test_ccs_binops_mia_cv('eq');
  test_ccs_binops_mia_cv('ne');
  test_ccs_binops_mia_cv('spaceship');
  test_ccs_binops_mia_cv('and2');
  test_ccs_binops_mia_cv('or2');
  test_ccs_binops_mia_cv('xor');
  test_ccs_binops_mia_cv('shiftleft');
  test_ccs_binops_mia_cv('shiftright');
}
#test_ccs_binops_mia_cv_all();

##---------------------------------------------------------------------
## CCS: binops (block-wise alignment / missing-is-annihiliator): #1
sub test_ccs_binops_mia_1 {
  our ($OPNAME,$swap) = @_;

  ##-- set operation sub
  our $OP = PDL->can($OPNAME);
  if (!defined($OP)) {
    warn("$0: no operation '$OPNAME' for PDL!");
    return undef;
  }
  $swap = 0 if (!defined($swap));

  ##-- get data
  test_data_1;
  $a = $a->setvaltobad(0); ##-- use BAD as missing value: ensure annihilator

  our $adims = [$a->dims];
  our $b = $a->rotate(1);
  $b->sever;
  $b->slice("0,0") .= -($a->at(0,0));

  our $ccsa = $a->toccs;
  our @a = our ($ixa,$nza,$za) = ($ccsa->whichND,$ccsa->nzvals,$ccsa->missing);

  our $ccsb = $b->toccs;
  our @b = our ($ixb,$nzb,$zb) = ($ccsb->whichND,$ccsb->nzvals,$ccsb->missing);
  our (@align,$nzai,$nzbi);


  ##-- test: missing-annihilator
  our $nnzc   = ($nza->dim(0) > $nzb->dim(0) ? $nza->dim(0) : $nzb->dim(0));
  our $istate = zeroes(long,7); ##-- [ nnzai,nnzai_nxt, nnzbi,nnzbi_nxt, nnzci,nnzci_nxt, cmpval ]
  our $ostate = $istate->pdl;

  ##-- guts: alignment
  ccs_binop_align_block_mia($ixa,$ixb,$istate, $nzai=zeroes(long,$nnzc),$nzbi=zeroes(long,$nnzc),$ostate);

  ##-- parse output state
  our ($nzai_cur,$nzai_nxt, $nzbi_cur,$nzbi_nxt, $nzci_cur,$nzci_nxt,$cmpval) = $ostate->list;

  ##-- trim input pdls
  our $nzci_max = $nzci_cur-1;
  $nzai = $nzai->slice("0:$nzci_max");
  $nzbi = $nzbi->slice("0:$nzci_max");

  ##-- construct output pdl: nzvals
  our $avals = $ccsa->vals;
  our $bvals = $ccsb->vals;
  our $nzc   = zeroes(($avals->type > $bvals->type ? $avals->type : $bvals->type), $nzci_cur);
  our $zc    = $OP->($avals->slice("-1"), $bvals->slice("-1"), $swap);
  $nzc      .= $OP->($avals->index($nzai), $bvals->index($nzbi), $swap);
  ##
  ##-- get indices of "good" c() values
  our $cimask  = ($zc->isbad ? ($nzc->isgood) : ($nzc!=$zc));
  our $ciwhich = $cimask->which;
  $nzc = $nzc->index($ciwhich); #->append($zc);

  ##-- construct output pdl: which
  #our $ndims   = $ixa->dim(0);
  #our $cwhich  = zeroes(long, $ndims,$ciwhich->nelem);
  #$cwhich     .= $ixa->dice_axis(1,$nzai->index($ciwhich));
  our $cwhich   = $ixa->dice_axis(1,$nzai->index($ciwhich));
  $cwhich->sever;
  isok("ccs_align_block_mia:OP=$OPNAME,SWAP=$swap", all matchpdl(ccs_decode($cwhich,$nzc,$zc,[$a->dims]), $OP->($a,$b,$swap)));

  my $foo='bar';
}
sub test_ccs_binops_mia_1_all {
  test_ccs_binops_mia_1('plus');
  test_ccs_binops_mia_1('minus');
  test_ccs_binops_mia_1('mult');
  test_ccs_binops_mia_1('divide');
  test_ccs_binops_mia_1('gt');
  test_ccs_binops_mia_1('ge');
  test_ccs_binops_mia_1('lt');
  test_ccs_binops_mia_1('le');
  test_ccs_binops_mia_1('eq');
  test_ccs_binops_mia_1('ne');
  test_ccs_binops_mia_1('spaceship');
  test_ccs_binops_mia_1('and2');
  test_ccs_binops_mia_1('or2');
  test_ccs_binops_mia_1('xor');
  test_ccs_binops_mia_1('shiftleft');
  test_ccs_binops_mia_1('shiftright');
}
#test_ccs_binops_mia_1_all();

##---------------------------------------------------------------------
## CCS: binops (block-wise alignment): #2
sub test_ccs_binops_2 {
  test_data_1;

  our $adims = [$a->dims];
  our $b = $a->rotate(1);
  $b->sever;
  $b->slice("0,0") .= -($a->at(0,0));

  our $ccsa = $a->toccs;
  our @a = our ($ixa,$nza,$za) = ($ccsa->whichND,$ccsa->nzvals,$ccsa->missing);

  our $ccsb = $b->toccs;
  our @b = our ($ixb,$nzb,$zb) = ($ccsb->whichND,$ccsb->nzvals,$ccsb->missing);
  our (@align,$nzai,$nzbi);

  ##-- set operation sub
  our $OPNAME = shift || 'plus';
  our $OP = PDL->can($OPNAME);
  if (!defined($OP)) {
    warn("$0: no operation '$OPNAME' for PDL!");
    return undef;
  }

  ##-- test: plus: batch
  our $nnzc   = $a->flat->nnz + $b->flat->nnz;
  our $istate = zeroes(long,7); ##-- [ nnzai,nnzai_nxt, nnzbi,nnzbi_nxt, nnzci,nnzci_nxt, cmpval ]
  our $ostate = $istate->pdl;

  ##-- guts: alignment
  ccs_binop_align_block_mia($ixa,$ixb,$istate, $nzai=zeroes(long,$nnzc),$nzbi=zeroes(long,$nnzc),$ostate);

  ##-- parse output state
  our ($nzai_cur,$nzai_nxt, $nzbi_cur,$nzbi_nxt, $nzci_cur,$nzci_nxt,$cmpval) = $ostate->list;

  ##-- trim input pdls
  our $nzci_max = $nzci_cur-1;
  $nzai = $nzai->slice("0:$nzci_max");
  $nzbi = $nzbi->slice("0:$nzci_max");

  ##-- find mis-alignments
  our $aimask   = ($nzai < 0);
  our $bimask   = ($nzbi < 0);
  our $aiwhich  = $aimask->which;
  our $biwhich  = $bimask->which;

  ##-- construct output pdl: nzvals
  our $avals     = $ccsa->vals;
  our $bvals     = $ccsb->vals;
  our $amissingi = $avals->nelem-1;
  our $bmissingi = $bvals->nelem-1;
  $nzai->index($aiwhich) .= $amissingi;
  $nzbi->index($biwhich) .= $bmissingi;
  ##
  our $nzc = zeroes(($avals->type > $bvals->type ? $avals->type : $bvals->type), $nzci_cur);
  our $zc  = $OP->($avals->slice("-1"), $bvals->slice("-1"), 0);
  $nzc    .= $OP->($avals->index($nzai), $bvals->index($nzbi), 0);
  ##
  ##-- get indices of "good" c() values
  our $cimask  = ($nzc != $zc);
  our $ciwhich = $cimask->which;
  $nzc = $nzc->index($ciwhich); #->append($zc);

  ##-- construct output pdl: which
  $bimask->inplace->not;
  $bimask &= $cimask;
  $bimask &= $aimask;    ##-- bimask is now = ($nzai==$amissingi & $nzbi!=$bmissingi & $nzc!=$zc)
  $aimask->inplace->not;
  $aimask &= $cimask;    ##-- aimask is now = ($nzai!=$amissingi)
  ##-- remap (a|b)i to "ciwhich" indices
  $nzai    = $nzai->index($ciwhich);
  $nzbi    = $nzbi->index($ciwhich);
  $aimask  = $aimask->index($ciwhich);
  $bimask  = $bimask->index($ciwhich);
  $aiwhich = $aimask->which;
  $biwhich = $bimask->which;

  our $ndims   = $ixa->dim(0);
  our $cwhich  = zeroes(long, $ndims,$ciwhich->nelem);
  $cwhich->dice_axis(1,$aiwhich) .= $ixa->dice_axis(1,$nzai->index($aiwhich));
  $cwhich->dice_axis(1,$biwhich) .= $ixb->dice_axis(1,$nzbi->index($biwhich));

  isok("ccs_align_block_mia:OP=$OPNAME", all matchpdl(ccs_decode($cwhich,$nzc,$zc,[$a->dims]), $OP->($a,$b,0)));

  print "test_ccs_binops_2: done.\n";
}
#test_ccs_binops_2('plus');
#test_ccs_binops_2('minus');
#test_ccs_binops_2('mult');
#test_ccs_binops_2('divide');
#test_ccs_binops_2('gt');
#test_ccs_binops_2('ge');
#test_ccs_binops_2('lt');
#test_ccs_binops_2('le');
#test_ccs_binops_2('eq');
#test_ccs_binops_2('ne');
#test_ccs_binops_2('spaceship');
#test_ccs_binops_2('and2');
#test_ccs_binops_2('or2');
#test_ccs_binops_2('xor');
#test_ccs_binops_2('shiftleft');
#test_ccs_binops_2('shiftright');


##---------------------------------------------------------------------
## CCS: binops (hard-coded)

sub test_ccs_binop_plus_byblock_1 {
  return ccs_binop_plus(@_[0..5]) if ($_[6] >= $_[1]->dim(0)+$_[4]->dim(0));

  my ($ixa,$nza,$za, $ixb,$nzb,$zb, $blocksize,$nzctype) = @_;
  $nzctype                   = $nza->type if (!defined($nzctype));
  my $ndims                  = $ixa->dim(0);
  my ($nnzai_max,$nnzbi_max) = ($nza->dim(0), $nzb->dim(0));
  my ($nnzai_from,$nnzbi_from,$nnzci_from) = (0,0,0);
  my ($nnzai_to,$nnzbi_to,$nnzci_to);
  ##
  ##-- initialize collectors
  my ($nnzai_nxt,$nnzbi_nxt,$nnzci_nxt) = map {PDL->pdl($_)->long} ($nnzai_from,$nnzbi_from,$nnzci_from);
  my $ixc = PDL->zeroes(PDL::long(), $ndims,$blocksize);
  my $nzc = PDL->zeroes($nzctype,    $blocksize);
  my $zc  = PDL->zeroes($nzctype,    1);
  do {
    $nnzai_to = $nnzai_from+$blocksize-1;
    $nnzbi_to = $nnzbi_from+$blocksize-1;

    $nnzai_to = $nnzai_max-1 if ($nnzai_to >= $nnzai_max);
    $nnzbi_to = $nnzbi_max-1 if ($nnzbi_to >= $nnzbi_max);

    last if ($nnzai_to < 0 || $nnzbi_to < 0);

    $nnzci_to = $nnzci_from+$blocksize-1;
    if ($ixc->dim(1) < $nnzci_to) {
      $ixc->reshape($ixc->dim(0), $nnzci_to+1);
      $nzc->reshape($nnzci_to+1);
    }

    &PDL::_ccs_binop_plus_int($ixa->slice(",${nnzai_from}:${nnzai_to}"), $nza->slice("${nnzai_from}:${nnzai_to}"), $za,
			      $ixb->slice(",${nnzbi_from}:${nnzbi_to}"), $nzb->slice("${nnzbi_from}:${nnzbi_to}"), $zb,
			      $ixc->slice(",${nnzci_from}:${nnzci_to}"), $nzc->slice("${nnzci_from}:${nnzci_to}"), $zc,
			      $nnzai_nxt, $nnzbi_nxt, $nnzci_nxt);

    $nnzai_from += $nnzai_nxt->sclr;
    $nnzbi_from += $nnzbi_nxt->sclr;
    $nnzci_from += $nnzci_nxt->sclr;
  } while ($nnzai_from < $nnzai_max && $nnzbi_from < $nnzbi_max);

  if ($nnzai_from < $nnzai_max) {
    ##-- gobble leftover data from $a()
    my $ix_gobble = $ixa->slice(",${nnzai_from}:-1");
    my $ngobbled  = $ix_gobble->dim(1);
    $nnzci_to = $nnzci_from+$ngobbled;
    if ($ixc->dim(1) < $nnzci_to) {
      $ixc->reshape($ixc->dim(0), $nnzci_to+1);
      $nzc->reshape($nnzci_to+1);
    }
    &PDL::_ccs_binop_plus_int($ix_gobble, $nza->slice("${nnzai_from}:-1"), $za,
			      $ix_gobble, $zb->slice("*${ngobbled}")->flat, $zb,
			      $ixc->slice(",${nnzci_from}:-1"), $nzc->slice("${nnzci_from}:-1"), $zc,
			      $nnzai_nxt, $nnzbi_nxt, $nnzci_nxt);
  }
  elsif ($nnzbi_from < $nnzbi_max) {
    ##-- gobble leftover data from $b()
    my $ix_gobble = $ixb->slice(",${nnzbi_from}:-1");
    my $ngobbled  = $ix_gobble->dim(1);
    $nnzci_to = $nnzci_from+$ngobbled;
    if ($ixc->dim(1) < $nnzci_to) {
      $ixc->reshape($ixc->dim(0), $nnzci_to+1);
      $nzc->reshape($nnzci_to+1);
    }
    &PDL::_ccs_binop_plus_int($ix_gobble, $za->slice("*${ngobbled}")->flat, $za,
			      $ix_gobble, $nzb->slice(",${nnzbi_from}:-1"), $zb,
			      $ixc->slice(",${nnzci_from}:-1"), $nzc->slice("${nnzci_from}:-1"), $zc,
			      $nnzai_nxt, $nnzbi_nxt, $nnzci_nxt);
  }

  ##-- finally, trim the output pdls & return
  my $trimto = $nnzci_from-1;
  $ixc = $ixc->slice(",0:$trimto");
  $nzc = $nzc->slice("0:$trimto");
  return wantarray ? ($ixc,$nzc,$zc) : $nzc;
}

sub test_binops_1 {
  test_data_1;

  our $adims = [$a->dims];
  our $b = $a->rotate(1);
  $b->sever;

  our $ccsa = $a->toccs;
  our @a = our ($ixa,$nza,$za) = ($ccsa->whichND,$ccsa->nzvals,$ccsa->missing);

  our $ccsb = $b->toccs;
  our @b = our ($ixb,$nzb,$zb) = ($ccsb->whichND,$ccsb->nzvals,$ccsb->missing);

  our @c = our ($ixc,$nzc,$zc);

  ##--------------------------------------
  ## Binary Operations: arithmetic

  ##-- test: binop: plus (batch)
  @c=($ixc,$nzc,$zc)= ccs_binop_plus(@a, @b);
  isok("ccs_binop_plus", all(ccs_decode(@c,$adims) == $a+$b));
  ##
  ##-- test: binop: plus (block-wise: see prototype perl method, above)
  our $blocksize = 7;
  @c=($ixc,$nzc,$zc)= test_ccs_binop_plus_byblock_1(@a, @b, $blocksize);
  isok("ccs_binop_plus_byblock_1", all(ccs_decode(@c,$adims) == $a+$b));
  ##--/block-wise computation

  ##-- test: binop: minus
  @c=($ixc,$nzc,$zc)= ccs_binop_minus(@a,@b);
  isok("ccs_binop_minus", all(ccs_decode(@c,$adims) == $a-$b));

  ##-- test: binop: mult
  @c=($ixc,$nzc,$zc)= ccs_binop_mult(@a,@b);
  isok("ccs_binop_mult", all(ccs_decode(@c,$adims) == $a*$b));

  ##-- test: binop: divide
  @c=($ixc,$nzc,$zc)= ccs_binop_divide(@a,@b);
  isok("ccs_binop_divide", all matchpdl(ccs_decode(@c,$adims), $a/$b));

  ##-- test: binop: modulo
  @c=($ixc,$nzc,$zc)= ccs_binop_modulo(@a,@b);
  isok("ccs_binop_modulo", all(ccs_decode(@c,$adims) == $a%$b));

  ##-- test: binop: power
  @c=($ixc,$nzc,$zc)= ccs_binop_power(@a,@b);
  isok("ccs_binop_power", all(ccs_decode(@c,$adims) == $a**$b));

  ##--------------------------------------
  ## Binary Operations: comparisons
  @c=($ixc,$nzc,$zc)= ccs_binop_gt(@a,@b);
  isok("ccs_binop_gt", all(ccs_decode(@c,$adims) == ($a>$b)));

  @c=($ixc,$nzc,$zc)= ccs_binop_ge(@a,@b);
  isok("ccs_binop_ge", all(ccs_decode(@c,$adims) == ($a>=$b)));

  @c=($ixc,$nzc,$zc)= ccs_binop_lt(@a,@b);
  isok("ccs_binop_lt", all(ccs_decode(@c,$adims) == ($a<$b)));

  @c=($ixc,$nzc,$zc)= ccs_binop_le(@a,@b);
  isok("ccs_binop_le", all(ccs_decode(@c,$adims) == ($a<=$b)));

  @c=($ixc,$nzc,$zc)= ccs_binop_eq(@a,@b);
  isok("ccs_binop_eq", all(ccs_decode(@c,$adims) == ($a==$b)));

  @c=($ixc,$nzc,$zc)= ccs_binop_ne(@a,@b);
  isok("ccs_binop_ne", all(ccs_decode(@c,$adims) == ($a!=$b)));

  @c=($ixc,$nzc,$zc)= ccs_binop_spaceship(@a,@b);
  isok("ccs_binop_spaceship", all(ccs_decode(@c,$adims) == ($a<=>$b)));

  ##--------------------------------------
  ## Binary Operations: bitwise & logical
  @c=($ixc,$nzc,$zc)= ccs_binop_and2(@a,@b);
  isok("ccs_binop_and2", all(ccs_decode(@c,$adims) == ($a&$b)));

  @c=($ixc,$nzc,$zc)= ccs_binop_or2(@a,@b);
  isok("ccs_binop_or2", all(ccs_decode(@c,$adims) == ($a|$b)));

  @c=($ixc,$nzc,$zc)= ccs_binop_xor(@a,@b);
  isok("ccs_binop_xor", all(ccs_decode(@c,$adims) == ($a^$b)));

  @c=($ixc,$nzc,$zc)= ccs_binop_shiftleft(@a,@b);
  isok("ccs_binop_shiftleft", all(ccs_decode(@c,$adims) == ($a<<$b)));

  @c=($ixc,$nzc,$zc)= ccs_binop_shiftright(@a,@b);
  isok("ccs_binop_shiftright", all(ccs_decode(@c,$adims) == ($a>>$b)));

  print "test_binops_1(): done.\n";
}
#test_binops_1();

##---------------------------------------------------------------------
## CCS: Ufuncs

sub test_ufuncs_1 {
  test_data_1;
  our $ccs = $a->toccs;

  our ($ix_out,$nzvals_out,$nnz_out);

  ##-- test: ccs_accum_*: sum
  our $which   = $ccs->whichND;
  our $which0  = $which->slice("(0),");
  our $which1  = $which->slice("1:-1");
  our $whichi1 = $which1->qsortveci;
  $which1      = $which1->dice_axis(1,$whichi1);
  our $nzvals1 = $ccs->_nzvals->index($whichi1);

  ($ix_out,$nzvals_out,$nnz_out) = ccs_accum_sum($which1,$nzvals1, 0,$a->dim(0));
  isok("ccs_accum_sum:missing=0", all(ccs_decode($ix_out,$nzvals_out)==$a->sumover));

  ##-- test: ccs_accum_*: sum (+missing)
  our $missing2 = 42;
  $a2 = $a->pdl;
  $a2->where($a==0) .= $missing2;
  ($ix_out2,$nzvals_out2) = ccs_accum_sum($which1,$nzvals1, $missing2,$a->dim(0));
  isok("ccs_accum_sum:missing=$missing2", all(ccs_decode($ix_out2,$nzvals_out2)==$a2->sumover));

  ##-- test: ccs_accum: prod
  ($ix_out,$nzvals_out,$nnz_out) = ccs_accum_prod($which1,$nzvals1, 0,$a->dim(0));
  isok("ccs_accum_prod:missing=0", all(ccs_decode($ix_out,$nzvals_out)==$a->prodover));
  ($ix_out2,$nzvals_out2) = ccs_accum_prod($which1,$nzvals1, $missing2,$a->dim(0));
  isok("ccs_accum_prod:missing=$missing2", all(ccs_decode($ix_out2,$nzvals_out2)==$a2->prodover));

  ##-- test: ccs_accum: dprod
  our ($la,$la2);
  $la  = $a->long;
  $la2 = $a2->long;
  $lnzvals1 = $nzvals1->long;
  ($ix_out2,$nzvals_out2) = ccs_accum_prod($which1,$lnzvals1, $missing2,$a->dim(0));
  isok("ccs_accum_prod:type", $nzvals_out2->type==$lnzvals1->type);
  ($ix_out2,$nzvals_out2) = ccs_accum_dprod($which1,$lnzvals1, $missing2,$a->dim(0));
  isok("ccs_accum_dprod:type", $lnzvals1->type==long && $nzvals_out2->type==double);
  isok("ccs_accum_dprod:data", all(ccs_decode($ix_out2,$nzvals_out2)==$a2->dprodover));

  ##-- test: ccs_accum: andover
  our ($ba,$ba2,$bwhich,$bwhich0,$bwhich1,$bnzvals1);
  $ba = $a->long;
  $ba->where(($a==0) & ($a->xvals%2==0)) .= 255;
  $ba2 = $ba->pdl;
  $ba2->where($ba==0) .= $missing2;
  $bwhich  = $ba->whichND;
  our $bwhichi1 = $bwhich->slice("-1:0,")->qsortveci;
  $bwhich   = $bwhich->dice_axis(1,$bwhichi1);
  $bwhich0  = $bwhich->slice("0,");
  $bwhich1  = $bwhich->slice("1:-1,");
  $bnzvals1 = $ba->indexND($bwhich);

  ($ix_out,$nzvals_out) = ccs_accum_and($bwhich1,$bnzvals1, 0,$ba->dim(0));
  isok("ccs_accum_and:missing=0", all(ccs_decode($ix_out,$nzvals_out)==$ba->andover));
  ($ix_out,$nzvals_out) = ccs_accum_and($bwhich1,$bnzvals1, $missing2,$ba->dim(0));
  isok("ccs_accum_and:missing=$missing2", all(ccs_decode($ix_out,$nzvals_out)==$ba2->andover));

  ($ix_out,$nzvals_out) = ccs_accum_or($bwhich1,$bnzvals1, 0,$ba->dim(0));
  isok("ccs_accum_or:missing=0", all(ccs_decode($ix_out,$nzvals_out)==$ba->orover));
  ($ix_out,$nzvals_out) = ccs_accum_or($bwhich1,$bnzvals1, $missing2,$ba->dim(0));
  isok("ccs_accum_or:missing=$missing2", all(ccs_decode($ix_out,$nzvals_out)==$ba2->orover));

  ($ix_out,$nzvals_out) = ccs_accum_band($bwhich1,$bnzvals1, 0,$ba->dim(0));
  isok("ccs_accum_band:missing=0", all(ccs_decode($ix_out,$nzvals_out)==$ba->bandover));
  ($ix_out,$nzvals_out) = ccs_accum_band($bwhich1,$bnzvals1, $missing2,$ba->dim(0));
  isok("ccs_accum_band:missing=$missing2", all(ccs_decode($ix_out,$nzvals_out)==$ba2->bandover));

  ($ix_out,$nzvals_out) = ccs_accum_bor($bwhich1,$bnzvals1, 0,$ba->dim(0));
  isok("ccs_accum_bor:missing=0", all(ccs_decode($ix_out,$nzvals_out)==$ba->borover));
  ($ix_out,$nzvals_out) = ccs_accum_bor($bwhich1,$bnzvals1, $missing2,$ba->dim(0));
  isok("ccs_accum_bor:missing=$missing2", all(ccs_decode($ix_out,$nzvals_out)==$ba2->borover));

  ##-- test: ccs_accum: maximum
  ($ix_out,$nzvals_out,$nnz_out) = ccs_accum_maximum($which1,$nzvals1, 0,$a->dim(0));
  isok("ccs_accum_maximum:missing=0", all(ccs_decode($ix_out,$nzvals_out)==$a->maximum));
  ($ix_out2,$nzvals_out2) = ccs_accum_maximum($which1,$nzvals1, $missing2,$a->dim(0));
  isok("ccs_accum_maximum:missing=$missing2", all(ccs_decode($ix_out2,$nzvals_out2)==$a2->maximum));

  ##-- test: ccs_accum: minimum
  ($ix_out,$nzvals_out,$nnz_out) = ccs_accum_minimum($which1,$nzvals1, 0,$a->dim(0));
  isok("ccs_accum_minimum:missing=0", all(ccs_decode($ix_out,$nzvals_out)==$a->minimum));
  ($ix_out2,$nzvals_out2) = ccs_accum_minimum($which1,$nzvals1, $missing2,$a->dim(0));
  isok("ccs_accum_minimum:missing=$missing2", all(ccs_decode($ix_out2,$nzvals_out2)==$a2->minimum));


  ##-- TODO:
  ##  + stats      : average, daverage, medover, oddmedover, pctover, median
  ##  + cumulative": ccs_accum_cumusum, ccs_accum_cumuprod, ..
  print "ufuncs: all done.\n";
}
test_ufuncs_1();

##---------------------------------------------------------------------
## CCS: operations

sub test_ops_1 {
  test_data_1;
  our $ccs = $a->toccs;
  our $xv = 1+sequence($a->type,  $a->dim(0));
  our $yv = 1+sequence($a->type,1,$a->dim(1));
  our @ccs = ccsencode($a);

  ##-- test: add: column vector (missing values are ignored : MISSING~annihilator)
  our ($ccs_xi,$ccs_yi,$newvals,$ccs2);

  ##-- plus: compat
  isok("compat:add_rv",  all( ccsdecode(@ccs[0,1],ccsadd_rv(@ccs,$xv)) == ($a+$xv)*($a!=0) ));
  isok("compat:add_cv",  all( ccsdecode(@ccs[0,1],ccsadd_cv(@ccs,$yv->flat))==($a+$yv)*($a!=0) ));
  isok("compat:mult_rv", all( ccsdecode(@ccs[0,1],ccsmult_rv(@ccs,$xv)) == ($a*$xv)*($a!=0) ));
  isok("compat:mult_cv", all( ccsdecode(@ccs[0,1],ccsmult_cv(@ccs,$yv->flat))==($a*$yv)*($a!=0) ));

  ##-- test: using wrappers in PDL::CCS::Functions
  ##
  ##-- plus
  $newvals = ccs_plus_vector_ma($ccs->[$WHICH]->slice("(0),"), $ccs->nzvals, $xv);
  $ccs2    = $ccs->shadow(which=>$ccs->[$WHICH], vals=>$newvals->append($ccs->missing));
  isok("ccs:plus::xv:missing~annihilator", all(($ccs2->decode==($a+$xv))->or2($a==0,0)));
  $newvals = ccs_plus_vector_ma($ccs->[$WHICH]->slice("(1)"), $ccs->nzvals, $yv->flat);
  $ccs2    = $ccs->shadow(which=>$ccs->[$WHICH], vals=>$newvals->append($ccs->missing));
  isok("ccs:plus:yv:missing~annihilator", all(($ccs2->decode==($a+$yv))->or2($a==0,0)));

  ##-- minus
  $newvals = ccs_minus_vector_ma($ccs->[$WHICH]->slice("(0),"), $ccs->nzvals, $xv);
  $ccs2    = $ccs->shadow(which=>$ccs->[$WHICH], vals=>$newvals->append($ccs->missing));
  isok("ccs:minus::xv:missing~annihilator", all(($ccs2->decode==($a-$xv))->or2($a==0,0)));

  ##-- mult
  $newvals = ccs_mult_vector_ma($ccs->[$WHICH]->slice("(0),"), $ccs->nzvals, $xv);
  $ccs2    = $ccs->shadow(which=>$ccs->[$WHICH], vals=>$newvals->append($ccs->missing));
  isok("ccs:mult::xv:missing~annihilator", all($ccs2->decode==($a*$xv)));
  $newvals = ccs_mult_vector_ma($ccs->[$WHICH]->slice("(1),"), $ccs->nzvals, $yv->flat);
  $ccs2    = $ccs->shadow(which=>$ccs->[$WHICH], vals=>$newvals->append($ccs->missing));
  isok("ccs:mult::xv:missing~annihilator", all($ccs2->decode==($a*$yv)));

  ##-- divide
  $newvals = ccs_divide_vector_ma($ccs->[$WHICH]->slice("(0),"), $ccs->nzvals, $xv);
  $ccs2    = $ccs->shadow(which=>$ccs->[$WHICH], vals=>$newvals->append($ccs->missing));
  isok("ccs:divide::xv:missing~annihilator", all( $ccs2->decode==($a/$xv) ));

  ##-- modulo
  $newvals = ccs_modulo_vector_ma($ccs->[$WHICH]->slice("(0),"), $ccs->nzvals, $xv);
  $ccs2    = $ccs->shadow(which=>$ccs->[$WHICH], vals=>$newvals->append($ccs->missing));
  isok("ccs:modulo::xv:missing~annihilator", all( $ccs2->decode==($a%$xv) ));

  ##-- power
  $newvals = ccs_power_vector_ma($ccs->[$WHICH]->slice("(0),"), $ccs->nzvals, $xv);
  $ccs2    = $ccs->shadow(which=>$ccs->[$WHICH], vals=>$newvals->append($ccs->missing));
  isok("ccs:power::xv:missing~annihilator", all( $ccs2->decode==($a**$xv) ));

  ##---------------------------
  ## Comparisons
  $newvals = ccs_le_vector_ma($ccs->[$WHICH]->slice("(0),"), $ccs->nzvals, $xv);
  $ccs2    = $ccs->shadow(which=>$ccs->[$WHICH], vals=>$newvals->append($ccs->missing));
  isok("ccs:le::xv:missing~annihilator", all( ($ccs2->decode==($a <= $xv))->or2($a==0,0) ));

  $newvals = ccs_lt_vector_ma($ccs->[$WHICH]->slice("(0),"), $ccs->nzvals, $xv);
  $ccs2    = $ccs->shadow(which=>$ccs->[$WHICH], vals=>$newvals->append($ccs->missing));
  isok("ccs:lt::xv:missing~annihilator", all( ($ccs2->decode==($a < $xv))->or2($a==0,0) ));

  $newvals = ccs_ge_vector_ma($ccs->[$WHICH]->slice("(0),"), $ccs->nzvals, $xv);
  $ccs2    = $ccs->shadow(which=>$ccs->[$WHICH], vals=>$newvals->append($ccs->missing));
  isok("ccs:ge::xv:missing~annihilator", all( ($ccs2->decode==($a >= $xv))->or2($a==0,0) ));

  $newvals = ccs_gt_vector_ma($ccs->[$WHICH]->slice("(0),"), $ccs->nzvals, $xv);
  $ccs2    = $ccs->shadow(which=>$ccs->[$WHICH], vals=>$newvals->append($ccs->missing));
  isok("ccs:gt::xv:missing~annihilator", all( ($ccs2->decode==($a > $xv))->or2($a==0,0) ));

  $newvals = ccs_eq_vector_ma($ccs->[$WHICH]->slice("(0),"), $ccs->nzvals, $xv);
  $ccs2    = $ccs->shadow(which=>$ccs->[$WHICH], vals=>$newvals->append($ccs->missing));
  isok("ccs:eq::xv:missing~annihilator", all( ($ccs2->decode==($a == $xv))->or2($a==0,0) ));

  $newvals = ccs_ne_vector_ma($ccs->[$WHICH]->slice("(0),"), $ccs->nzvals, $xv);
  $ccs2    = $ccs->shadow(which=>$ccs->[$WHICH], vals=>$newvals->append($ccs->missing));
  isok("ccs:ne::xv:missing~annihilator", all( ($ccs2->decode==($a != $xv))->or2($a==0,0) ));

  $newvals = ccs_spaceship_vector_ma($ccs->[$WHICH]->slice("(0),"), $ccs->nzvals, $xv);
  $ccs2    = $ccs->shadow(which=>$ccs->[$WHICH], vals=>$newvals->append($ccs->missing));
  isok("ccs:spaceship::xv:missing~annihilator", all( ($ccs2->decode==($a <=> $xv))->or2($a==0,0) ));

  ##---------------------------
  ## Logic / Bitwise
  our $xbv = $xv % 2;

  $newvals = ccs_and2_vector_ma($ccs->[$WHICH]->slice("(0),"), $ccs->nzvals, $xbv);
  $ccs2    = $ccs->shadow(which=>$ccs->[$WHICH], vals=>$newvals->append($ccs->missing));
  isok("ccs:and2::xv:missing~annihilator", all( $ccs2->decode==$a->and2($xbv,0) ));

  $newvals = ccs_or2_vector_ma($ccs->[$WHICH]->slice("(0),"), $ccs->nzvals, $xbv);
  $ccs2    = $ccs->shadow(which=>$ccs->[$WHICH], vals=>$newvals->append($ccs->missing));
  isok("ccs:or2::xv:missing~annihilator", all( ($ccs2->decode==($a->or2($xbv,0)))->or2($a==0,0) ));

  $newvals = ccs_xor_vector_ma($ccs->[$WHICH]->slice("(0),"), $ccs->nzvals, $xbv);
  $ccs2    = $ccs->shadow(which=>$ccs->[$WHICH], vals=>$newvals->append($ccs->missing));
  isok("ccs:xor:xv:missing~annihilator", all( ($ccs2->decode==($a->xor($xbv,0)))->or2($a==0,0) ));

  $newvals = ccs_shiftleft_vector_ma($ccs->[$WHICH]->slice("(0),"), $ccs->nzvals, $xv);
  $ccs2    = $ccs->shadow(which=>$ccs->[$WHICH], vals=>$newvals->append($ccs->missing));
  isok("ccs:shiftleft:xv:missing~annihilator", all(  $ccs2->decode==($a << $xv) ));

  $newvals = ccs_shiftright_vector_ma($ccs->[$WHICH]->slice("(0),"), $ccs->nzvals, $xv);
  $ccs2    = $ccs->shadow(which=>$ccs->[$WHICH], vals=>$newvals->append($ccs->missing));
  isok("ccs:shiftright:xv:missing~annihilator", all( $ccs2->decode==($a >> $xv) ));
}
test_ops_1();


##---------------------------------------------------------------------
## test: non-zero

sub test_nnz {
  our $p = pdl(double, [ [0,1,2], [0,0,1e-7], [0,1,0], [1,1,1] ]);

  our $nnz0 = $p->slice(",(0)")->nnz;
  isok("nnz:slice", $nnz0==2);

  our $nnzf = $p->nnz;
  isok("nnz:flat", all($nnzf==pdl([2,1,1,3])));

  our $nnza = $p->nnza($eps);
  isok("nnza:flat", all($nnza==pdl([2,0,1,3])));
}
#test_nnz; ##-- ok

##---------------------------------------------------------------------
## CCS: encode

sub test_ccsencode_1 {
  test_data_1;
  our $nnz   = $p->flat->nnz;
  our $nrows = $p->dim(1);
  ccsencodefull($p,
		$ptr=zeroes(long,$p->dim(0)),
		$rowids=zeroes(long,$nnz),
		$nzvals=zeroes($p->type,$nnz));

  isok("ccsencodefull:ptr",    all($ptr==pdl([0,3,7,9,12,16])));
  isok("ccsencodefull:rowids", all($rowids==pdl([0,1,3,1,2,4,5,2,3,2,3,4,0,3,4,5,1,4,5])));
  isok("ccsencodefull:nzvals", all($nzvals==pdl([10,3,3,9,7,8,4,8,8,7,7,9,-2,5,9,2,3,13,-1])));
}
#test_ccsencode_1;

sub test_ccsencode_2 {
  test_data_1();

  ##-- get pointer: native perl
  our $pw         = $p->whichND->vv_qsortvec;
  our ($pxf,$pxi) = rle($pw->slice("(0),"));
  $pxf->reshape($pxf->nnz);
  $pxi->reshape($pxf->nelem);

  our $ptr1       = zeroes(long,$p->dim(0)+1);
  $ptr1->index($pxi+1) .= $pxf;
  $ptr1 .= $ptr1->cumusumover;
}
#test_ccsencode_2;

sub test_ccs_encode_pointers {
  test_data_1();

  our $awhich = $a->whichND();
  our $avals  = $a->indexND($awhich);
  our ($aptr0,$awi0) = ccs_encode_pointers($awhich->slice("(0),"));
  our ($aptr1,$awi1) = ccs_encode_pointers($awhich->slice("(1),"));

  ##-- ok, we've got a raw N+1 pointer: build compatible CCS encoded matrix
  our $ccs_rowids = $awhich->slice("(1),")->index($awi0);
  our $ccs_vals   = $avals->index($awi0);
  our (@ccs);
  our ($ptr,$rowids,$nzvals) = @ccs = ($aptr0->slice("0:-2"), $ccs_rowids, $ccs_vals);

  isok("ccs_encode/pointer:ptr",    all($ptr==pdl([0,3,7,9,12,16])));
  isok("ccs_encode/pointer:rowids", all($rowids==pdl([0,1,3,1,2,4,5,2,3,2,3,4,0,3,4,5,1,4,5])));
  isok("ccs_encode/pointer:nzvals", all($nzvals==pdl([10,3,3,9,7,8,4,8,8,7,7,9,-2,5,9,2,3,13,-1])));
}
#test_ccs_encode_pointers;

##---------------------------------------------------------------------
## CCS: Compat (& lots of other stuff)

sub test_ccs_compat {
  test_data_1();

  our $awhich = $a->whichND();
  our $avals  = $a->indexND($awhich);

  ##-- object test
  our $ccs = $a->toccs;

  ##-- object: dimensional shuffling
  isok("obj:reorder(1,0)",             all($ccs->reorder(1,0)->decode==$a->reorder(1,0)));
  isok("obj:post-reorder(1,0):decode", all($ccs->decode==$a));
  isok("obj:xchg(0,1)",                all($ccs->xchg(0,1)->decode==$a->xchg(0,1)));
  isok("obj:post-xchg(0,1):decode",    all($ccs->decode==$a));
  isok("obj:xchg(0,-1)",               all($ccs->xchg(0,-1)->decode==$a->xchg(0,-1)));
  isok("obj:post-xchg(0,-1):decode",   all($ccs->decode==$a));
  isok("obj:mv(0,1)",                  all($ccs->mv(0,1)->decode==$a->mv(0,1)));
  isok("obj:post-mv(0,1):decode",      all($ccs->decode==$a));
  isok("obj:mv(1,0)",                  all($ccs->mv(1,0)->decode==$a->mv(1,0)));
  isok("obj:post-mv(1,0):decode",      all($ccs->decode==$a));

  ##-- object: re-sort : BROKEN
  if (0) {
    ($pix,$pnzix)   = $ccs->ptr(0);
    our ($pix0,$pnzix0) = map {$_->pdl} ($pix,$pnzix);
    $cwhich  = $ccs->whichND;
    our $cwhich0 = $ccs->whichND->pdl;
    isok("resort:pre:pointer", all($cwhich->dice_axis(1,$pnzix)==$cwhich->vv_qsortvec));
    $cvals = $ccs->vals;
    $cnewi  = random($cwhich->dim(1))->qsorti; ##-- randomize: get permutation
    $cwhich .= $cwhich->dice_axis(1,$cnewi); ##-- permute: which
    $cvals->slice("0:-2") .= $cvals->index($cnewi); ##-- permute: nzvals
    $pnzix->index($cnewi) .= $pnzix0;
    our $pnzix1  = $pnzix->pdl;
    our $cwhich1 = $cwhich->pdl;
    isok("resort:randomized:decode",  all($ccs->decode==$a));
    isok("resort:randomized:pointer", all($cwhich->dice_axis(1,$pnzix)==$cwhich->vv_qsortvec));
    isok("resort:randomized:pointer:decode",
	 all($a==ccsdecode($pix->slice("0:-2"), $cwhich->slice("(1),")->index($pnzix), $cvals->index($pnzix))));
    $ccs->sortwhich();
    $cwhich = $ccs->whichND;
    $cvals  = $ccs->vals;
    ($pix,$pnzix) = $ccs->ptr(0);
    isok("resort:randomized+resorted:decode", all($ccs->decode==$a));
    isok("resort:randomized+resorted:pointer", all($cwhich->dice_axis(1,$pnzix)==$cwhich->vv_qsortvec));
    isok("resort:randomized+resorted:pointer:decode",
	 all($a==ccsdecode($pix->slice("0:-2"), $cwhich->slice("(1),")->index($pnzix), $cvals->index($pnzix))));
  }

  ##-- object: at,set
  isok("obj:at:nz", $ccs->at(4,3)==$a->at(4,3));
  isok("obj:at:z",  $ccs->at(3,1)==$a->at(3,1) && $ccs->at(3,1)==$ccs->missing);

  $oldval = $ccs->at(4,3);
  $a->set(4,3, 42);
  $ccs->set(4,3, 42);
  isok("obj:set:nz", $ccs->at(4,3)==$a->at(4,3));
  $ccs->set(4,3, $oldval);
  $a->set(4,3, $oldval);

  ##-- works, but doesn't test well
  #eval { $ccs->set(3,1, 42); };
  #isok("obj:set:z:error",  $@);

  ##-- object: dice_axis
  our $axisi = pdl(long,[2,4]);
  $afound = $a->dice_axis(0,$axisi);
  $cfound = $ccs->dice_axis(0,$axisi);
  isok("obj:dice_axis:0", all($afound==$cfound->decode));
  $afound = $a->dice_axis(1,$axisi);
  $cfound = $ccs->dice_axis(1,$axisi);
  isok("obj:dice_axis:1", all($afound==$cfound->decode));

  ##-- object: indexND
  our $find = pdl(long,[[0,0],[1,0],[1,1]]);
  our $afound = $a->indexND($find);
  our $cfound = $ccs->indexND($find);
  isok("obj:indexND", all($cfound==$afound));

  ##-- object: index (flat): NO (PSEUDO)-THREADING!
  $find = pdl(long,[2,4,6,8]);
  $afound = $a->flat->index($find);
  $cfound = $ccs->index($find);
  isok("obj:index:flat", all($cfound==$afound));

  ##-- object: which (flat)
  isok("obj:which:flat", all($ccs->which==$a->which));

  ##-- object: string
  our $ccs_str = $ccs->string; ##-- ok
  isok("obj:nd:string", $ccs_str);

  ##-- object: todense
  isok("obj:decode", all($ccs->decode==$a));
  isok("obj:toccs",  $ccs->lstring eq $ccs->toccs->lstring);
  isok("a:todense",  ${$a->todense}==$$a);

  ##-- object: basic props
  isok("obj:type", $ccs->type==$a->type);
  isok("obj:dims", all(pdl($ccs->dims)==pdl($a->dims)));
  isok("obj:ndims", $ccs->ndims==$a->ndims);
  isok("obj:nelem", $ccs->nelem==$a->nelem);
  isok("obj:ngood", $ccs->ngood==$a->ngood);
  isok("obj:nnz",   $ccs->nnz==$a->flat->nnz); ##-- needs dimension fix!
  isok("obj:nstored", $ccs->nstored==$a->flat->nnz);
  isok("obj:whichND", all($ccs->whichND->vv_qsortvec==$a->whichND->vv_qsortvec));
  isok("obj:missing", defined($ccs->missing));
  isok("obj:nzvals",  all($ccs->nzvals==$a->indexND($ccs->whichND)));

  ##-- ok, we've got a raw N+1 pointer: build compatible CCS encoded matrix
  our ($aptr0,$awi0) = ccs_encode_pointers($awhich->slice("(0),"));
  our ($aptr1,$awi1) = ccs_encode_pointers($awhich->slice("(1),"));
  our $ccs_rowids = $awhich->slice("(1),")->index($awi0);
  our $ccs_vals   = $avals->index($awi0);
  our (@ccs0);
  @ccs0 = our ($ptr0,$rowids0,$nzvals0) = ($aptr0->slice("0:-2"), $ccs_rowids, $ccs_vals);
  our $ptr_want    =pdl(long,[0,3,7,9,12]);
  our $rowids_want =pdl(long,[0,1,3,1,2,4,5,2,3,2,3,4,0,3,4,5]);
  our $nzvals_want =pdl(long,[10,3,3,9,7,8,4,8,8,7,7,9,-2,5,9,2]);
  isok("ccs_encode/pointer:ptr",    all($ptr0==$ptr_want));
  isok("ccs_encode/pointer:rowids", all($rowids0==$rowids_want));
  isok("ccs_encode/pointer:nzvals", all($nzvals0==$nzvals_want));

  ##-- decode: compat
  our $adc = ccsdecode(@ccs0);
  isok("ccs_decode/compat:0", all($a==$adc));
  $adc = ccsdecodeg(@ccs0);
  isok("ccs_decode/compat:BAD", all($a->where($a)==$adc->where($a)) && all(($a!=0)==$adc->isgood));

  ##-- decode: pointer:
  our ($aptr0d,$aptr0di) = ccs_decode_pointer($aptr0);
  our ($aptr1d,$aptr1di) = ccs_decode_pointer($aptr1);
  isok("ccs_decode/pointer:dim=0", all($aptr0d == $awhich->slice("(0),")->index($awi0)));
  isok("ccs_decode/pointer:dim=1", all($aptr1d == $awhich->slice("(1),")->index($awi1)));

  ##-- decode: whichND:  ----- CONTINUE HERE -----
  isok("ccs_decode:whichND",       all(ccs_decode($awhich,$avals,0) == $a));
  isok("ccs_decode:whichND:dim=0", all(ccs_decode($awhich->dice_axis(1,$awi0),$avals->index($awi0)) == $a));
  isok("ccs_decode:whichND:dim=1", all(ccs_decode($awhich->dice_axis(1,$awi1),$avals->index($awi1)) == $a));

  ##-- project+decode (~ dice_axis along pointer)
  our $aproj0  = pdl(long,[1,2,4]);
  our $aslice0 = $a->dice_axis(0,$aproj0);
  our ($aproj0d,$apnzi0d) = ccs_decode_pointer($aptr0,$aproj0);
  our $apnzi = $awi0->index($apnzi0d);
  our $wproj = $aproj0d->slice("*1,")->append($awhich->slice("1:-1,")->dice_axis(1,$apnzi));
  our $vproj = $avals->index($apnzi);

  our $aprojd = ccs_decode($wproj,$vproj, 0, [$aproj0->nelem, map {$a->dim($_)} (1..($a->ndims-1))]);
  isok("project+decode~dice_axis", all($aprojd==$aslice0));

  our $aprojcd = ccsdecodecols($ptr0,$rowids0,$nzvals0, $aproj0,0,$a->dim(1));
  isok("ccsdecodecols~dice_axis", all($aprojcd==$aslice0));


  ##-- encode: compat
  ($ptr,$rowids,$nzvals) = ccsencode($a);
  isok("ccs_encode/compat:0:ptr",    all($ptr==$ptr_want));
  isok("ccs_encode/compat:0:rowids", all($rowids==$rowids_want));
  isok("ccs_encode/compat:0:nzvals", all($nzvals==$nzvals_want));

  ($ptr,$rowids,$nzvals) = ccsencode_naz($a,0.5);
  isok("ccs_encode/compat:~0:ptr",    all($ptr==$ptr_want));
  isok("ccs_encode/compat:~0:rowids", all($rowids==$rowids_want));
  isok("ccs_encode/compat:~0:nzvals", all($nzvals==$nzvals_want));

  ($ptr,$rowids,$nzvals) = ccsencode_g($a->setvaltobad(0));
  isok("ccs_encode/compat:bad:ptr",    all($ptr==$ptr_want));
  isok("ccs_encode/compat:bad:rowids", all($rowids==$rowids_want));
  isok("ccs_encode/compat:bad:nzvals", all($nzvals==$nzvals_want));

  ($ptr,$rowids,$nzvals) = ccsencode_i($a->which, $a->where($a), $a->dim(0));
  isok("ccs_encode/compat:i:ptr",    all($ptr==$ptr_want));
  isok("ccs_encode/compat:i:rowids", all($rowids==$rowids_want));
  isok("ccs_encode/compat:i:nzvals", all($nzvals==$nzvals_want));

  ($ptr,$rowids,$nzvals) = ccsencode_i2d($a->whichND->xchg(0,1)->dog, $a->where($a!=0));
  isok("ccs_encode/compat:i2d:ptr",    all($ptr==$ptr_want));
  isok("ccs_encode/compat:i2d:rowids", all($rowids==$rowids_want));
  isok("ccs_encode/compat:i2d:nzvals", all($nzvals==$nzvals_want));

  ##-- lookup: compat
  @ccs = ($ptr,$rowids,$nzvals) = ccsencode($a);
  isok("compat:which",    all(ccswhich(@ccs)->qsort == $a->which->qsort));
  isok("compat:whichND",  all(ccswhichND(@ccs)->vv_qsortvec == $a->whichND->vv_qsortvec));
  isok("compat:get:good", all(ccsget(@ccs,$a->which,0) == $a->flat->index($a->which)));
  isok("compat:get:bad",  all(ccsget(@ccs,(($a->which+1)%$a->nelem),0) == $a->flat->index(($a->which+1)%$a->nelem)));
  our ($xi,$yi) = $a->whichND->xchg(0,1)->dog;
  isok("compat:get2d:good", all(ccsget2d(@ccs,$xi,$yi) == $a->index2d($xi,$yi)));
  isok("compat:get2d:bad",  all(ccsget2d(@ccs,$xi,$yi->rotate(1)) == $a->index2d($xi,$yi->rotate(1))));
  our @ccsT = our ($ptrT,$rowidsT,$nzvalsT) = ccstranspose(@ccs);
  isok("compat:transpose", all(ccsdecode(@ccsT)==$a->xchg(0,1)));
}
test_ccs_compat();

##---------------------------------------------------------------------
## DUMMY
##---------------------------------------------------------------------
foreach $i (0..3) {
  print "--dummy($i)--\n";
}

