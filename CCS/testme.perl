#!/usr/bin/perl -wd

use lib qw(./blib/lib ./blib/arch);
use PDL;
use PDL::VectorValued;
#use PDL::CCS::Old;
use PDL::CCS::Utils;
use PDL::CCS::Functions;
use PDL::CCS::Compat;
use PDL::CCS::Nd;
use Storable qw(store retrieve);
use PDL::IO::Storable;
use Data::Dumper;

BEGIN {
  $, = ' ';
  our $eps=1e-6;

  our $DIMS  = $PDL::CCS::Nd::DIMS;
  our $XDIMS = $PDL::CCS::Nd::XDIMS;
  our $WHICH = $PDL::CCS::Nd::WHICH;
  our $VALS  = $PDL::CCS::Nd::VALS;
  our $PTRS  = $PDL::CCS::Nd::PTRS;

  our $P_BYTE = PDL::byte();
  our $P_LONG = PDL::long();
}

##---------------------------------------------------------------------
## test subs
sub isok {
  my ($label,$bool) = @_;
  print "test($label): ", ($bool ? "ok" : "NOT ok"), "\n";
}

sub make { system('make'); }

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
}
#test_data_bg();

##---------------------------------------------------------------------
## CCS: operatopns

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

