#!/usr/bin/perl -w

use lib qw(. ./blib/lib ./blib/arch);
use PDL;
use PDL::CCS;

BEGIN{ $, = ' '; $eps=1e-6; }

##---------------------------------------------------------------------
## N nonzero
use vars qw($p1 $nnz0 $nnzf $nnza);
sub tnnz {
  $p = pdl(double, [ [0,1,2], [0,0,1e-7], [0,1,0], [1,1,1] ]);
  $nnz0 = $p->slice(",(0)")->nnz;    ##-- should be 2
  $nnzf = $p->flat->nnz;             ##-- should be 7 (==sum([2,1,1,3]))
  $nnza   = $p->flat->nnza($eps); ##-- should be 6 (==sum([2,0,1,3]))
}

##---------------------------------------------------------------------
## CCS: encode
use vars qw($ptr $rowids $nzvals $nrows);
sub tccs1 {
  $p = pdl(double, [
		    [10,0,0,0,-2,0],
		    [3,9,0,0,0,3],
		    [0,7,8,7,0,0],
		    [3,0,8,7,5,0],
		    [0,8,0,9,9,13],
		    [0,4,0,0,2,-1],
		   ]);
}

##---------------------------------------------------------------------
## CCS: encode: big
sub tccs1b {
  $a=$p  = random(100,100);
  $a    *= ($a>.9);
}

sub tccs2 {
  $nnz = $p->flat->nnz; ##-- 19
  $nrows = $p->dim(1);
  ccsencodefull($p,
		$ptr=zeroes(long,$p->dim(0)),    ##-- [0,3,7,9,12,16]
		$rowids=zeroes(long,$nnz),       ##-- [0,1,3,1,2,4,5,2,3,2,3,4,0,3,4,5,1,4,5]
		$nzvals=zeroes($p->type, $nnz)); ##-- [10,3,3,9,7,8,4,8,8,7,7,9,-2,5,9,2,3,13,-1]
}
sub tccs2a {
  $nnz = $p->flat->nnza($eps);  ##-- 19
  $nrows = $p->dim(1);
  ccsencodefulla($p,$eps,
		 $ptr=zeroes(long,$p->dim(0)),    ##-- [0,3,7,9,12,16]
		 $rowids=zeroes(long,$nnz),       ##-- [0,1,3,1,2,4,5,2,3,2,3,4,0,3,4,5,1,4,5]
		 $nzvals=zeroes($p->type, $nnz)); ##-- [10,3,3,9,7,8,4,8,8,7,7,9,-2,5,9,2,3,13,-1]
}


sub tccs2c  { ($ptr,$rowids,$nzvals)=ccsencode($p);  $nrows=$p->dim(1); }
sub tccs2ca { ($ptr,$rowids,$nzvals)=ccsencodea($p); $nrows=$p->dim(1); }

##---------------------------------------------------------------------
## CCS: encode: which

sub tccs_enc_w2d {
  if (!defined($a)) { tccs1(); $a=$p; }
  ($awcols,$awrows) = $a->whichND;
  $avals            = $a->index2d($awcols,$awrows);

  ($ptr,$rowids,$nzvals) = ccsencode($a);

  ccsencodefull_i2d($awcols,$awrows,$avals,
		    $wptr   =zeroes(long,     $a->dim(0)),
		    $wrowids=zeroes(long,     $avals->dim(0)),
		    $wnzvals=zeroes($a->type, $avals->dim(0)));

  print "encoding: full: indexND: ", (all(ccsdecode($wptr,$wrowids,$wnzvals)==$a) ? "ok" : "NOT ok"), "\n";

  ($wptr,$wrowids,$wnzvals) = ccsencode_i2d($awcols,$awrows,$avals);
  print "encoding: wrapped: indexND: ", (all(ccsdecode($wptr,$wrowids,$wnzvals)==$a) ? "ok" : "NOT ok"), "\n";
}
#tccs_enc_w2d();

sub tccs_enc_wflat {
  if (!defined($a)) { tccs1(); $a=$p; }
  $awflat = $a->which;
  $avals  = $a->flat->index($awflat);

  ($ptr,$rowids,$nzvals) = ccsencode($a);

  ccsencodefull_i($awflat, $avals,
		  $wptr   =zeroes(long,     $a->dim(0)),
		  $wrowids=zeroes(long,     $avals->dim(0)),
		  $wnzvals=zeroes($a->type, $avals->dim(0)),
		 );
  print "encoding: full: flat: ", (all(ccsdecode($wptr,$wrowids,$wnzvals)==$a) ? "ok" : "NOT ok"), "\n";

  ($wptr,$wrowids,$wnzvals) = ccsencode_i($awflat, $avals, $a->dim(0));
  print "encoding: wrapped: flat: ", (all(ccsdecode($wptr,$wrowids,$wnzvals)==$a) ? "ok" : "NOT ok"), "\n";
}
#tccs_enc_wflat();


##---------------------------------------------------------------------
## CCS: decode
sub tccsd {
  $p2 = zeroes(double, $nrows, $ptr->dim(0));
  ccsdecodefull($ptr,$rowids,$nzvals, $p2);
  print "decoding: ", (all($p==$p2) ? "ok" : "NOT ok"), "\n";
}

sub tccsdc {
  $p2 = ccsdecode($ptr,$rowids,$nzvals);
  print "decoding: ", (all($p==$p2) ? "ok" : "NOT ok"), "\n";
}


##---------------------------------------------------------------------
## CCS: Ops: transpose

sub tp_data {
  $a = pdl([[1,0,0,2],
	    [0,3,4,0],
	    [5,0,0,6]]);
  ($ptr,$rowids,$nzvals) = ccsencode($a);
}

use vars qw($ptrT $rowidsT $nzvalsT);
sub tp_test1 {
  tp_data();
  ccstransposefull($ptr,$rowids,$nzvals, $ptrT=zeroes($a->dim(1)),$rowidsT=zeroes($rowids),$nzvalsT=zeroes($nzvals));

  print
    ("transposefull(): ",
     (all(ccsdecode($ptr,$rowids,$nzvals)==ccsdecode($ptrT,$rowidsT,$nzvalsT)->xchg(0,1)) ? "ok" : "NOT ok"), "\n",
     );
}
sub tp_test2 {
  tp_data();
  ($ptrT,$rowidsT,$nzvalsT) = ccstranspose($ptr,$rowids,$nzvals);

  print
    ("transpose(): ",
     (all(ccsdecode($ptr,$rowids,$nzvals)==ccsdecode($ptrT,$rowidsT,$nzvalsT)->xchg(0,1)) ? "ok" : "NOT ok"), "\n",
     );
}


##---------------------------------------------------------------------
## CCS: access

sub ta_data {
  $a = pdl([[1,0,0,2],
	    [0,3,4,0],
	    [5,0,0,6]]);
  ($ptr, $rowids, $nzvals)  = ccsencode($a);
  ($ptrT,$rowidsT,$nzvalsT) = ccstranspose($ptr,$rowids,$nzvals);
}

sub ta_index_flat {
  ta_data();

  ##-- all present
  $awhich = $a->which;
  $avals  = $a->flat->index($awhich);
  $cvals  = ccsget($ptr,$rowids,$nzvals, $awhich,0);
  print "index (all present): ", (all($cvals==$avals) ? "ok" : "NOT ok"), "\n";

  ##-- some present (zero)
  $allai    = sequence(long,$a->nelem);
  $allavals = $a->flat->index($allai);
  $allcvals = ccsget($ptr,$rowids,$nzvals, $allai,0);
  print "index (some present / zero): ", (all($allavals==$allcvals) ? "ok" : "NOT ok"), "\n";

  ##-- some present (bad)
  $badval    = pdl(0)->setvaltobad(0);
  $allbcvals = ccsget($ptr,$rowids,$nzvals, $allai,$badval);
  print
    "index (some present / bad): ", (all($allbcvals->where($allbcvals->isgood) == $allavals->where($allbcvals->isgood))
					 &&
					 all($allavals->where($allbcvals->isbad) == 0)
					 ? "ok" : "NOT ok"), "\n";
}
#ta_index_flat();


sub ta_index_2d {
  ta_data();

  ##-- all present
  ($acoli,$arowi) = $a->whichND;
  $avals          = $a->index2d($acoli,$arowi);
  $cvals          = ccsget2d($ptr,$rowids,$nzvals, $acoli,$arowi,0);
  print "index2d (all present): ", (all($cvals==$avals) ? "ok" : "NOT ok"), "\n";

  ##-- some present (zero)
  ($allacoli,$allarowi) = ($a->xvals->flat,$a->yvals->flat);
  $allavals = $a->index2d($allacoli,$allarowi);
  $allcvals = ccsget2d($ptr,$rowids,$nzvals, $allacoli,$allarowi,0);
  print "index2d (some present / zero): ", (all($allavals==$allcvals) ? "ok" : "NOT ok"), "\n";

  ##-- some present (bad)
  $badval    = pdl(0)->setvaltobad(0);
  $allbcvals = ccsget2d($ptr,$rowids,$nzvals, $allacoli,$allarowi,$badval);
  print
    "index2d (some present / bad): ", (all($allbcvals->where($allbcvals->isgood) == $allavals->where($allbcvals->isgood))
					 &&
					 all($allavals->where($allbcvals->isbad) == 0)
					 ? "ok" : "NOT ok"), "\n";
}
#ta_index_2d();



sub ta_whichfull {
  ta_data;
  ccswhichfull($ptr,$rowids,$nzvals,    $wcols=zeroes(long,$nzvals->dim(0)),  $wrows=zeroes(long,$nzvals->dim(0)));
  ccswhichfull($ptrT,$rowidsT,$nzvalsT, $wcolsT=zeroes(long,$nzvals->dim(0)), $wrowsT=zeroes(long,$nzvals->dim(0)));

  print
    ("test_whichfull(row-primary): ",
     (all($wrowsT->cat($wcolsT)->xchg(0,1) == scalar($a->whichND)) ? "ok" : "NOT ok"), "\n",
     ##--
     "test_whichfull(col-primary): ",
     (all($wrows ->cat($wcols) ->xchg(0,1) == scalar($a->xchg(0,1)->whichND)) ? "ok" : "NOT ok"), "\n",
    );
}
#ta_whichfull();

##---------------------------------------------------------------------
## CCS: Ops

sub top_data {
  $a = pdl([[1,0,0,2],
	    [0,3,4,0],
	    [5,0,0,6]]);
  our @ccs  = ($ptr, $rowids, $nzvals)  = ccsencode($a);
  our @ccst = ($ptrT,$rowidsT,$nzvalsT) = ccstranspose($ptr,$rowids,$nzvals);
}

use vars qw($cv $nzvals_cv $rv $nzvals_rv);
sub top_mult {
  top_data();
  $nzvals_rv = ccsmult_rv($ptr,$rowids,$nzvals, $rv=10**(sequence($a->dim(0))+1));
  $nzvals_cv = ccsmult_cv($ptr,$rowids,$nzvals, $cv=10**(sequence($a->dim(1))+1));
  print
    ("mult-by-row-vector]: ",
     (all(ccsdecode($ptr,$rowids,$nzvals_rv) == ($a * $rv)) ? "ok" : "NOT ok"), "\n",
     ##--
     "mult-by-col-vector]: ",
     (all(ccsdecode($ptr,$rowids,$nzvals_cv) == ($a * $cv->slice("*1,"))) ? "ok" : "NOT ok"), "\n",
    );

  our $matv           = sequence($a->dims);
  our $nzvals_matv_rv = ccsmult_rv($ptr,$rowids,$nzvals, $matv)->xchg(0,1)->sumover;
  our $nzvals_matv_cv = ccsmult_cv($ptr,$rowids,$nzvals, $matv->xchg(0,1))->xchg(0,1)->sumover;
}
#top_mult();

sub top_sumover {
  top_data;

  $a_sum  = $a->sumover;
  $c_sum  = ccssumover($ptr,$rowids,$nzvals);
  print "sumover: ", (all($a_sum==$c_sum) ? "ok" : "NOT ok"), "\n";

  $at_sum = $a->xchg(0,1)->sumover;
  $ct_sum = ccssumovert($ptr,$rowids,$nzvals);
  print "sumoverT: ", (all($at_sum==$ct_sum) ? "ok" : "NOT ok"), "\n";
}
#top_sumover;

##---------------------------------------------------------------------
## Test: CCS package
##---------------------------------------------------------------------
package PDL::CCS::Obj;
use PDL::Lite;

BEGIN {
  our @ISA = qw(PDL);
  *isa = \&UNIVERSAL::isa;
  our $DIMS   = 0;
  our $MISSING = 1;
  our $PTR    = 2;
  our $ROWIDS = 3;
  our $NZVALS = 4;
  our @CCS    = ($PTR,$ROWIDS,$NZVALS);
}

##--------------------------------------------------------------
## Constructors etc.

## $obj = $class_or_obj->newFromDense($dense2d);
## $obj = $class_or_obj->newFromDense($dense2d,$missing);
sub newFromDense {
  my ($that,$p,$missing) = @_;
  $missing     = (defined($missing)
		  ? PDL->pdl($p->type,$missing)
		  : ($p->badflag
		     ? PDL->pdl($p->type,0)->setvaltobad(0)
		     : PDL->pdl($p->type,0)));
  my $pwhichND = ($missing->isbad->any
		  ? $p->isgood()
		  : ($p != $missing)
		 )->whichND;
  my $pnz      = PDL->pdl($p->indexND($pwhichND));
  return bless([
		[$p->dims],
		$missing,
		PDL::CCS::ccsencode_i2d($pwhichND->xchg(0,1)->dog, $pnz, $p->getdim(0)),
	       ], ref($that)||$that);
}

## DESTROY : avoid PDL inheritance?
sub DESTROY { ; }

## $ccs = $pdl->toccs()
## $ccs = $pdl->toccs($missing)
#BEGIN { *PDL::toccs = \&toccs; }
sub toccs {
  return $_[0] if (isa($_[0],'PDL::CCS::Obj'));
  return PDL::CCS::Obj->newFromDense(@_)
}

## $obj = $obj->clone()
BEGIN { *clone = \&copy; }
sub copy { bless [ [@{$_[0][$DIMS]}], map {PDL->pdl($_)} @{$_[0]}[$MISSING,@CCS] ], ref($_[0]); }

## $obj = $obj->clonei($nzvals)
##  + like clone(), but only clones dimensions & indices
sub clonei {
  bless [
	 [@{$_[0][$DIMS]}],
	 PDL->pdl($_[0][$MISSING]),
	 PDL->pdl($_[0][$PTR]),
	 PDL->pdl($_[0][$ROWIDS]),
	 PDL->topdl($_[1]),
	], ref($_[0]);
}

##--------------------------------------------------------------
## Encoding & Decoding

## $dense = $obj->decode()
sub decode {
  PDL::CCS::ccsdecodecols(@{$_[0]}[@CCS],
			  PDL->sequence(PDL::long(), $_[0][$DIMS][0]),
			  $_[0][$MISSING],
			  $_[0][$DIMS][1]);
}

##--------------------------------------------------------------
## Basic Properties

## $type = $obj->type()
sub type { $_[0][$NZVALS]->type; }

## @dims = $obj->dims()
## ##@dims = $obj->dims(@dims)
sub dims { @{$_[0][$DIMS]}; }

## $ndims = $obj->ndims()
sub ndims { scalar(@{$_[0][$DIMS]}); }

## $dim = $obj->getdim($dimnum)
sub getdim { $_[0][$DIMS][$_[1]]; }

## $nelem = $obj->nelem
sub nelem { PDL->pdl(PDL::long(),@{$_[0][$DIMS]})->dprod; }

##--------------------------------------------------------------
## Low-level CCS access

## $missing = $obj->missing()
## $missing = $obj->missing($missing)
sub missing { $_[0][$MISSING]=$_[1] if (@_>1);  $_[0][$MISSING]; }

## $ptr = $obj->ptr()
## $ptr = $obj->ptr($ptr)
sub ptr { $_[0][$PTR]=$_[1] if (@_ > 1); $_[0][$PTR]; }

## $rowids = $obj->rowids()
## $rowids = $obj->rowids($rowids)
sub rowids { $_[0][$ROWIDS]=$_[1] if (@_ > 1); $_[0][$ROWIDS]; }

## $nzvals = $obj->nzvals();
## $nzvals = $obj->nzvals($nzvals)
sub nzvals { $_[0][$NZVALS]=$_[1] if (@_ > 1); $_[0][$NZVALS]; }

## ($ptr,$rowids,$nzvals) = $obj->ccs()
## ($ptr,$rowids,$nzvals) = $obj->ccs($ptr,$rowids,$nzvals)
sub ccs { @{$_[0]}[@CCS]=@_[1,2,3] if (@_ > 3); @{$_[0]}[@CCS]; }

##--------------------------------------------------------------
## Index Access: which(), whichND(), & friends

## $obj2 = $obj->convert($TYPE)
##  + always copies everything, like PDL method
sub convert { $_[0]->clonei($_[0][$NZVALS]->convert($_[1])); }

## $obj = $obj->_convert($TYPE)
##  + in-place conversion
sub _convert {
  $_[0][$NZVALS] = $_[0][$NZVALS]->convert($_[1])->inplace;
  return $_[0];
}


## $which = $obj->which()
sub which { PDL::CCS::ccswhich(@{$_[0]}[@CCS]); }

## $which2d = $obj->which2d()
sub which2d { PDL::CCS::ccswhich2d(@{$_[0]}[@CCS]); }

## $whichND = $obj->whichND()
sub whichND { PDL::CCS::ccswhichND(@{$_[0]}[@CCS]); }

## $colids = $obj->nzxvals()
sub nzxvals { $_[0]->which2d->slice("(0),"); }

## $rowids = $obj->nzyvals()
sub nzyvals { $_[0][$ROWIDS]; }

## $nzix = $obj->itonzi($ix, [o]$nzix)
##  + BROKEN: missing values aren't recognized:
##        print $ccs->itonzi(sequence($ccs->nelem));
##     - prints:           [  0   0   0   4   4   2   3   3   1   1   1   5  ]
##     - but SHOULD print: [  0  BAD BAD  4  BAD  2   3  BAD  1  BAD BAD  1  ]
sub itonzi { PDL::CCS::ccsitonzi(@{$_[0]}[$PTR,$ROWIDS], $_[1], $_[0][$MISSING], @_[2..$#_]); }

## $nzix = $obj->i2dtonzi($xvals, $yvals, [o]$nzix)
sub i2dtonzi { PDL::CCS::ccsi2dtonzi(@{$_[0]}[$PTR,$ROWIDS], @_[1,2], $_[0][$MISSING], @_[3..$#_]); }

## $nzix = $obj->iNDtonzi($xyvals, [o]$nzix)
##  + currently 2d-only
sub iNDtonzi { $_[0]->i2dtonzi($_[1]->xchg(0,1)->dog, @_[2..$#_]); }

##--------------------------------------------------------------
## Perl Access

## $val = $ccs->at($x,$y)
sub at { PDL::CCS::ccsget2d(@{$_[0]}[@CCS], @_[1,2], $_[0][$MISSING])->sclr; }

## $ccs = $ccs->set($x,$y,$val)
sub set {
  $_[0][$NZVALS]->index($_[0]->i2dtonzi(@_[1,2])) .= $_[3];
  return $_[0];
}

##--------------------------------------------------------------
## Slicing & Dicing

## $objT = $obj->transpose()
sub transpose {
  bless [
	 [@{$_[0][$DIMS]}[1,0,2..$#{$_[0][$DIMS]}]],
	 @{$_[0]}[$MISSING],
	 PDL::CCS::ccstranspose(@{$_[0]}[@CCS]),
	], ref($_[0]);
}

##--------------------------------------------------------------
## inplace stuff (pass to nzvals)

## $inplace = $ccs->inplace()
sub inplace { $_[0][$NZVALS]->inplace(@_[1..$#_]); }

## $bool = $ccs->is_inplace()
sub is_inplace { $_[0][$NZVALS]->is_inplace(@_[1..$#_]); }


##--------------------------------------------------------------
## Operations: Accumulators

## $sumover = $ccs->sumover([o]$sumover)
sub sumover { PDL::CCS::ccssumover(@{$_[0]}[@CCS], $_[0][$DIMS][1], @_[1..$#_]); }

## $sumoverT = $ccs->sumovert([o]$sumovert)
sub sumovert { PDL::CCS::ccssumovert(@{$_[0]}[@CCS], @_[1..$#_]); }

##--------------------------------------------------------------
## Operations: Binary

## $ccs3 = $ccs1->add($ccs2)
## $ccs3 = $ccs1->add($pdl)
##  + TODO: write a "real" ccs_add() routine: @ccs3 = ccs_add(@ccs1, @ccs2)
##    - may require FAST INDEX ACCESS (e.g. assume sorted & do a binary search on indices)
##    - if it works, maybe we can use the same trick for matmult() & friends
sub ccs_add {
}


##--------------------------------------------------------------
## Stringification & Viewing

## $str = $obj->string()
sub string {
  return
    (
     ''
     ."PDL::CCS("     .join(',', @{$_[0][$DIMS]}).")\n"
     .$,." ptr:("    .join(',', $_[0][$PTR]->type,     $_[0][$PTR]->dims)    .")=$_[0][$PTR]\n"
     .$,." rowids:(" .join(',', $_[0][$ROWIDS]->type,  $_[0][$ROWIDS]->dims) .")=$_[0][$ROWIDS]\n"
     .$,." nzvals:(" .join(',', $_[0][$NZVALS]->type,  $_[0][$NZVALS]->dims) .")=$_[0][$NZVALS]\n"
     .$,." missing:(".join(',', $_[0][$MISSING]->type, $_[0][$MISSING]->dims).")=$_[0][$MISSING]\n"
    );
}

## $xpdl = $obj->xpdl()
##  + expanded pdl:
##    [
##     [xi yi BAD nzval BAD missing]
##     ...
##    ]
sub xpdl {
  my $w   = $_[0]->whichND;
  my $z   = PDL->topdl($_[0][$MISSING]);
  my $nz  = $_[0][$NZVALS];
  my $bad = PDL->pdl($nz->type,0)->setvaltobad(0);
  $w->append($bad)->append($nz->slice("*1,"))->append($bad)->append($z->slice("*1,"));
}

## $xstr = $obj->xstring()
##  + string for expanded pdl
sub xstring { $_[0]->xpdl->string(); }

##--------------------------------------------------------------
## Operator overloading
BEGIN {
  package PDL::CCS::Obj;
  use overload (
		"\"\"" => \&string,
	       );
}



package main;

sub ccsobj_data {
  our $al = pdl(long,
		[[1,0,0,2],
		 [0,3,4,0],
		 [5,0,0,6]]);
  our $ad = $al->convert(double);
  our $b  = $al->setvaltobad(0);
  our $a  = $al;
}

# isok($label,@_) -- prints helpful label
use Test::More;
sub isok {
  my $label = shift;
  if (@_==1) {
    ok($_[0],$label);
  } elsif (@_==2) {
    is($_[0],$_[1], $label);
  } else {
    die("isok(): expected 1 or 2 non-label arguments, but got ", scalar(@_));
  }
}

sub test_ccsobj {
  ccsobj_data();
  $ccs = PDL::CCS::Nd->newFromDense($a);

  ##-- pdl() [--> PDL::new()]
  $ccs2 = pdl($ccs);
  isok("pdl(\$ccs)", "$ccs2" eq "$ccs");

  ##-- toccs(): ok
  print "a->toccs/my=", (my $tmp=$a->toccs()); ##--ok
  ##--error:
  ## + (in cleanup) Error - argument is not a recognised data structure at <<HERE>>
  ## + solution: override DESTROY()
  print "a->toccs/lit=", $a->toccs();
  ##--/error
  #print "newFromDense(a)=", @{[PDL::CCS::Obj->newFromDense($a)]}; ##-- also error
  isok("toccs", $a->toccs->string eq $ccs->string);

  ##-- decode()
  our $ccsd = $ccs->decode();
  isok("decode", all($ccsd==$a));

  ##-- encode/decode: bad
  our $bad   = pdl(0)->setvaltobad(0);
  our $ccsb  = $b->toccs();
  our $ccsbd = $ccsb->decode();
  isok("enc/dec:bad:isgood", all($ccsbd->isgood==$b->isgood));
  isok("enc/dec:bad:values", all($ccsbd->where($b->isgood) == $b->where($b->isgood)));

  ##-- transpose
  isok("transpose", all($a->transpose == $ccs->transpose->decode));

  ##-- basic properties: ndims, dims, nelem, type
  isok("ndims", $a->ndims==$ccs->ndims);
  isok("dims", all(pdl($a->dims)==pdl($ccs->dims)));
  isok("nelem", $a->nelem==$ccs->nelem);
  isok("type", $a->type==$ccs->type);

  ##-- which (requires working qsortvec)
  isok("which", all($a->which->qsort == $ccs->which->qsort));
  isok("whichND", all($a->whichND->qsortvec == $ccs->whichND->qsortvec));

  ##-- sumover, sumovert
  isok("sumover",  all($a->sumover == $ccs->sumover));
  isok("sumovert", all($a->xchg(0,1)->sumover == $ccs->sumovert));

  print "done.\n";
}
#test_ccsobj();


##--------------------------------------------------------------
## test: matrix stuff
use PDL::CCS::Nd;
sub test_matstuff {
  my ($M,$N,$O) = (2,3,4);
  my $a = sequence($M,$N);
  my $b = (sequence($O,$M)+1)*10;
  my $c = $a x $b;

  my $az = $a->toccs;
  my $cz = zeroes($O,$N);
  ccs_matmult2d_zdd($az->_whichND,$az->_nzvals, $b, $cz);
  isok("matmult2d_zdd/missing=0", all($cz==$c));
  ##
  my $czs = $cz->zeroes;
  my $abnil = (($az->missing+zeroes($M,1)) x $b)->flat;
  ccs_matmult2d_sdd($az->_whichND,$az->_nzvals,$az->missing, $b, $abnil, $czs);
  isok("matmult2d_sdd/missing=0/func", all($czs==$c));
  ##
  my $czo = $az->matmult2d_sdd($b);
  isok("matmult2d_sdd/missing=0/obj", all($czo==$c));
  my $czoz = $az->matmult2d_zdd($b);
  isok("matmult2d_zdd/missing=0/obj", all($czoz==$c));

  ##-- try to get a handle on non-zero "missing" values
  my $a1 = $a->pdl;
  $a1->where(($a%2)==0) .= 1;
  my $c1 = $a1 x $b;
  my $az1 = $a1->toccs(1);
  my $cz1 = zeroes($O,$N);
  my $abnil1 = (($az1->missing+zeroes($a1->dim(0),1)) x $b)->flat;
  ccs_matmult2d_sdd($az1->_whichND,$az1->_whichVals,$az1->missing, $b, $abnil1, $cz1);
  isok("matmult2d_sdd/missing=1/func", all($cz1==$c1));
  ##
  my $cz1o = $az1->matmult2d_sdd($b);
  isok("matmult2d_sdd/missing=1/obj", all($cz1o==$c1));


  print "$0: test_matstuff() done -- what now?\n";
}
#test_matstuff;

##---------------------------------------------------------------------
## debug: empty index

##----------------------------------------------------------------------
sub testnull {
  my $a = null;
  my $as = $a->toccs;
  my $asd = $as->decode;
  print "testnull: ", (all($a==$asd) ? "ok" : "NOT ok"), "\n";
}
#testnull();

##----------------------------------------------------------------------
sub testgt {
  my $a   = sequence(4,3);
  my $b   = pdl(42);
  my $c   = ($a > $b);
  ##
  my $as  = $a->toccs;
  my $bs  = $b->toccs;
  my $cs  = ($as > $bs);
  ##
  my $asd = $as->decode;
  my $bsd = $bs->decode;
  my $csd = $cs->decode;
  ##
  print "testnull_gt:csd ", (all($csd==$c) ? "ok" : "NOT ok"), "\n";
}
#testnull_gt();

##----------------------------------------------------------------------
sub testnull_min {
  my $a     = sequence(4,3);
  my $vals  = null;

  ## everything works as expected with a null-dimensioned empty pidde:
  my $i_0x0 = pdl([]);
  $a->indexND($i_0x0) = $vals; ##-- ok, expensive null-op

  ## ... but a 2x0 empty piddle pukes
  my $i_2x0 = sequence(2,1)->dice_axis(1,null); ##-- create an empty 2x0 piddle
  $a->indexND($i_2x0) .= $vals;

  ## barf()s with:
  # Problem with assignment: PDL::Ops::assgn(a,b): Parameter 'b':
  # Mismatched implicit thread dimension 0: should be 0, is 2
  # 	 at Basic/Core/Core.pm.PL (i.e. PDL::Core.pm) line 256
  # 	PDL::Core::barf('PDL=SCALAR(0xa026630)', 'PDL=SCALAR(0xa026690)') called at Basic/Core/Core.pm.PL (i.e. PDL::Core.pm) line 701
  # 	eval {...} called at Basic/Core/Core.pm.PL (i.e. PDL::Core.pm) line 700
  # 	__ANON__[Basic/Core/Core.pm.PL (i.e. PDL::Core.pm):712]('PDL=SCALAR(0xa026690)', 'PDL=SCALAR(0xa026630)', undef) called at (eval 160)[/usr/share/perl/5.10/perl5db.pl:638] line 2

}
#testnull_min();

##---------------------------------------------------------------------
## test: qsort

sub test_ccs_qsort {
  my @dims = (5,4);
  my $x  = xvals(@dims) * (2*((xvals(@dims)+yvals(@dims)))-1);
  my $xx = $x->toccs;
  ##
  my $xi = $x->qsorti;
}

##---------------------------------------------------------------------
## test: io

use PDL::CCS::IO::FastRaw;
sub test_io_fraw {
  my $a = sequence(4,3);
  my ($tmp);
  ($tmp=$a->where($a%2)) .= 0;

  my $ccs = $a->toccs;
  $ccs->writefraw("a.raw")
    or die("writefraw a.raw failed: $!");

  my $ccs2 = PDL::CCS::Nd->readfraw("a.raw");
  all($ccs2->decode==$ccs->decode)
    or die("readfraw a.raw inconsistent");

  my $ccs3 = PDL::CCS::Nd->mapfraw("a.raw",{ReadOnly=>1,Creat=>0,Trunc=>0});
  all($ccs3->decode==$ccs->decode)
    or die("mapfraw a.raw inconsistent");

  exit 0;
}
#test_io_fraw();

##---------------------------------------------------------------------
## test vcos

sub test_ccs_vcos {
  my $a = pdl([[1,2,3,4],[1,2,2,1],[-1,-2,-3,-4]])->xchg(0,1);
  my $b = pdl([1,2,3,4]);
  my $ccs = $a->toccs;
  my $vcos = $ccs->vcos_zdd($b);
  my $vcos_want = pdl([1,0.8660254,-1]);
  all($vcos->approx($vcos_want,1e-4))
    or die("vcos test failed : got=$vcos != want=$vcos_want");
  print "vcos ok\n";

  my $b3 = $b->slice(",*3");
  my $vcos3 = $ccs->vcos_zdd($b3);
  my $vcos3_want = $vcos_want->slice(",*3");
  all($vcos3->approx($vcos3_want,1e-4))
    or die("vcos test failed : got=$vcos3 != want=$vcos3_want");
  print "vcos:threaded ok\n";
}
#test_ccs_vcos();

##---------------------------------------------------------------------
## test_ufunc_2_039
## + see http://www.cpantesters.org/cpan/report/9ea0e2bc-a3e5-11eb-8917-5f9a624ac675
sub test_ufunc_2_039 {
  print STDERR "$0: using PDL v$PDL::VERSION\n";
  
  ##-- common data (CCS/t/common.plt)
  our $a = pdl(double, [
			[10,0,0,0,-2],
			[3,9,0,0,0],
			[0,7,8,6,0],
			[3,0,8,7,5],
			[0,8,0,9,7],
			[0,4,0,0,2],
		       ]);
  our $abad   = ($a==0);
  our $agood  = !$abad;
  our $awhich = $a->whichND;
  our $avals  = $a->indexND($awhich);
  our $BAD = pdl(0)->setvaltobad(0);

  ##-- CCS/t/03_ufuncs.t
  my ($ufunc_name, $missing_val) =
    ('bandover',0)
    #('bandover','BAD')
    ;
  print "test_ufunc($ufunc_name, $missing_val)\n";

  my $pdl_ufunc = PDL->can("${ufunc_name}")
    or die("no PDL Ufunc ${ufunc_name} defined!");
  my $ccs_ufunc = PDL::CCS::Nd->can("${ufunc_name}")
    or die("no CCS Ufunc PDL::CCS::Nd::${ufunc_name} defined!");

  $missing_val = 0 if (!defined($missing_val));
  $missing_val = PDL->topdl($missing_val);
  if ($missing_val->isbad) { $a = $a->setbadif($abad); }
  else                     { $a->where($abad) .= $missing_val; $a->badflag(0); }

  ##-- sorting with bad values doesn't work right in PDL-2.015 ; ccs/vv sorts BAD as minimal, PDL sort BAD as maximal: wtf?
  if ($ufunc_name =~ /qsort/ && $missing_val->isbad) {
    my $inf = $^O =~ /MSWin32/i ? (99**99)**99 : 'inf';
    $missing_val = $inf;
    $a->inplace->setbadtoval($inf);
  }

  my $ccs      = $a->toccs($missing_val);
  $ccs->_whichND($ccs->_whichND->ccs_indx()) if ($ccs->_whichND->type != PDL::ccs_indx());
  my $dense_rc = $pdl_ufunc->($a);
  my $ccs_rc   = $ccs_ufunc->($ccs);

  if ($ufunc_name =~ /_ind$/) {
    ##-- hack: adjust $dense_rc for maximum_ind, minimum_ind
    $dense_rc->where( $a->index2d($dense_rc,sequence($a->dim(1))) == $missing_val ) .= -1;
  } elsif ($ufunc_name =~ /qsorti$/) {
    ##-- hack: adjust $dense_rc for qsorti()
    my $ccs_mask = $dense_rc->zeroes;
    $ccs_mask->indexND( scalar($ccs_rc->whichND) ) .= 1;
    $dense_rc->where( $ccs_mask->not ) .= $ccs_rc->missing;
  }
  my $label = "${ufunc_name}:missing=$missing_val";

  ##-- check output type
  my $HAVE_PDL_2_039 = version->parse($PDL::VERSION) >= version->parse("2.039");
  SKIP: {
      skip("${label}:type - bandover, borover for PDL >= v2.039 ",1) if ($HAVE_PDL_2_039);
      isok("${label}:type", $ccs_rc->type, $dense_rc->type);
  }

  ##-- check output values
  SKIP: {
    ##-- RT bug #126294 (ses also analogous tests in CCS/Ufunc/t/01_ufunc.t)
    skip("RT #126294 - PDL::borover() appears to be broken", 1)
      if ($label eq 'borover:missing=BAD' && pdl([10,0,-2])->setvaltobad(0)->borover->sclr != -2);

    pdlok("${label}:vals", $ccs_rc->decode, $dense_rc);
  }
}
test_ufunc_2_039();

# pdlok($label, $got, $want)
sub pdlok {
  my ($label,$got,$want) = @_;
  $got  = PDL->topdl($got) if (defined($got));
  $want = PDL->topdl($want) if (defined($want));
  my $ok = (defined($got) && defined($want)
	    && cmp_dims($got,$want)
	    && all(matchpdl($want,$got))
	   );
  isok(labstr($label,$ok,$got,$want), $ok);
}

# cmp_dims($got_pdl,$expect_pdl)
sub cmp_dims {
  my ($p1,$p2) = @_;
  return $p1->ndims==$p2->ndims && all(pdl(PDL::long(),[$p1->dims])==pdl(PDL::long(),[$p2->dims]));
}

# matchpdl($a,$b) : returns pdl identity check, including BAD
sub matchpdl {
  my ($a,$b) = map {PDL->topdl($_)->setnantobad} @_[0,1];
  return ($a==$b)->setbadtoval(0) | ($a->isbad & $b->isbad) | ($a->isfinite->not & $b->isfinite->not);
}
# matchpdl($a,$b,$eps) : returns pdl approximation check, including BAD
sub matchpdla {
  my ($a,$b) = map {$_->setnantobad} @_[0,1];
  my $eps = $_[2];
  $eps    = 1e-5 if (!defined($eps));
  return $a->approx($b,$eps)->setbadtoval(0) | ($a->isbad & $b->isbad) | ($a->isfinite->not & $b->isfinite->not);
}

sub pdlstr {
  my $a = shift;
  return '(undef)' if (!defined($a));
  my $typ = UNIVERSAL::can($a,'type') ? $a->type : 'NOTYPE';
  my $str = "($typ) $a";
  #$str =~ s/\n/ /g;
  return $str;
}
sub labstr {
  my ($label,$ok,$got,$want) = @_;
  $label .= "\n  :     got=".pdlstr($got)."\n  :  wanted=".pdlstr($want) if (!$ok);
  return $label;
}

##---------------------------------------------------------------------
## DUMMY
##---------------------------------------------------------------------
foreach $i (0..3) {
  print "--dummy($i)--\n";
}

