#!/usr/bin/perl -wd

use lib qw(./blib/lib ./blib/arch);
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
## DUMMY
##---------------------------------------------------------------------
foreach $i (0..100) {
  print "--dummy($i)--\n";
}

