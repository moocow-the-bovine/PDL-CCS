# -*- Mode: CPerl -*-
# t/02_encode.t: test ccs encoding (compat)

$TEST_DIR = './t';
#use lib qw(../blib/lib ../blib/arch); $TEST_DIR = '.'; # for debugging

# load common subs
use Test;
do "$TEST_DIR/common.plt";
use PDL;
use PDL::CCS;

BEGIN { plan tests=>28, todo=>[]; }

##-- setup
$p = pdl(double, [
		  [10,0,0,0,-2,0],
		  [3,9,0,0,0,3],
		  [0,7,8,7,0,0],
		  [3,0,8,7,5,0],
		  [0,8,0,9,9,13],
		  [0,4,0,0,2,-1],
		 ]);

$nnz = $p->flat->nnz;

$want_ptr=pdl(long,[0,3,7,9,12,16]);
$want_rowids=pdl(long,[0,1,3,1,2,4,5,2,3,2,3,4,0,3,4,5,1,4,5]);
$want_nzvals=pdl(long,[0,1,3,1,2,4,5,2,3,2,3,4,0,3,4,5,1,4,5]);

##-- 1--3: test ccsencodefull()
ccsencodefull($p,
	      $ptr=zeroes(long,$p->dim(0)),
	      $rowids=zeroes(long,$nnz),
	      $nzvals=zeroes($p->type, $nnz));

isok("encodefull():ptr", all($ptr==$want_ptr));
isok("encodefull():rowids", all($rowids==$want_rowids));
isok("encodefull():nzvals", all($nzvals==$want_nzvals));

##-- 4--6: test ccsencode()
($ptr,$rowids,$nzvals) = ccsencode($p);
isok("encodefull():ptr", all($ptr==$want_ptr));
isok("encodefull():rowids", all($rowids==$want_rowids));
isok("encodefull():nzvals", all($nzvals==$want_nzvals));


##-- 7--9: test ccsencodefulla()
$eps=2.5;
$want_ptr_a=pdl(long,[0,3,7,9,12,14]);
$want_rowids_a=pdl(long,[0,1,3,1,2,4,5,2,3,2,3,4,3,4,1,4]);
$want_nzvals_a=pdl(long,[10,3,3,9,7,8,4,8,8,7,7,9,5,9,3,13]);
$nnz = $p->flat->nnza($eps);
ccsencodefulla($p, $eps,
	       $ptra=zeroes(long,$p->dim(0)),
	       $rowidsa=zeroes(long,$nnz),
	       $nzvalsa=zeroes($p->type, $nnz));
isok("encodefulla():ptr", all($ptra==$want_ptr_a));
isok("encodefulla():rowids", all($rowidsa==$want_rowids_a));
isok("encodefulla():nzvals", all($nzvalsa==$want_nzvals_a));

##-- 10--12: : test ccsencodea()
($ptra,$rowidsa,$nzvalsa) = ccsencodea($p,$eps);
isok("encodefulla():ptr",    all($ptra==$want_ptr_a));
isok("encodefulla():rowids", all($rowidsa==$want_rowids_a));
isok("encodefulla():nzvals", all($nzvalsa==$want_nzvals_a));

##-- 13..15 : test ccsencodefull_i2d()
#($pwcols,$pwrows) = $p->whichND; ##-- in pdl-2.4.9_014: WARNING - deprecated list context for whichND (may switch to scalar case soon)
($pwcols,$pwrows) = $p->whichND->xchg(0,1)->dog;
$pwvals           = $p->index2d($pwcols,$pwrows);
$nnz              = $pwvals->nelem;
ccsencodefull_i2d($pwcols,$pwrows,$pwvals,
		  $ptr=zeroes(long,$p->dim(0)),
		  $rowids=zeroes(long,$nnz),
		  $nzvals=zeroes($p->type, $nnz));

isok("encodefull_i2d():ptr",    all($ptr==$want_ptr));
isok("encodefull_i2d():rowids", all($rowids==$want_rowids));
isok("encodefull_i2d():nzvals", all($nzvals==$want_nzvals));

##-- 16..18 : test ccsencode_i2d()
($ptr,$rowids,$nzvals) =  ccsencode_i2d($pwcols,$pwrows,$pwvals);
isok("encode_i2d():ptr",    all($ptr==$want_ptr));
isok("encode_i2d():rowids", all($rowids==$want_rowids));
isok("encode_i2d():nzvals", all($nzvals==$want_nzvals));

##-- 19..21 : test ccsencodefull_i()
$pwhich = $p->which;
$pwvals = $p->flat->index($pwhich);
$nnz    = $pwvals->nelem;
ccsencodefull_i($pwhich, $pwvals,
		$ptr   =zeroes(long,$p->dim(0)),
		$rowids=zeroes(long,$nnz),
		$nzvals=zeroes($p->type, $nnz));

isok("encodefull_i():ptr",    all($ptr==$want_ptr));
isok("encodefull_i():rowids", all($rowids==$want_rowids));
isok("encodefull_i():nzvals", all($nzvals==$want_nzvals));

##-- 22..24 : test ccsencode_i()
our $N = $p->dim(0);
($ptr,$rowids,$nzvals) = ccsencode_i($pwhich, $pwvals, $N);

isok("encode_i():ptr",    all($ptr==$want_ptr));
isok("encode_i():rowids", all($rowids==$want_rowids));
isok("encode_i():nzvals", all($nzvals==$want_nzvals));


##-- 25 : test ccsdecodecols (single col)
our $M = $p->dim(1);
($ptr,$rowids,$nzvals) = ccsencode($p);

$col0 = ccsdecodecols($ptr,$rowids,$nzvals, 0,0);
isok("decodecols(0)", all($col0==$p->slice("0,")));

##-- 26 : test ccsdecodecols (full)
$dense = ccsdecodecols($ptr,$rowids,$nzvals, sequence($p->dim(0)),0);
isok("decodecols(all)", all($dense==$p));


##-- 27 : test decodefull()
$p2 = zeroes($p->type,$p->dims);
ccsdecodefull($ptr,$rowids,$nzvals, $p2);
isok("decodefull()", all($p==$p2));

##-- 28 : test decode()
$p2 = ccsdecode($ptr,$rowids,$nzvals);
isok("decode()", all($p==$p2));

print "\n";
# end of t/*.t

