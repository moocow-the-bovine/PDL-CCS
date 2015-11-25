##-*- Mode: CPerl -*-
$TEST_DIR = './t';
#use lib qw(../blib/lib ../blib/arch); $TEST_DIR = '.'; # for debugging

use Test::More tests=>(4+20);
use PDL;
use PDL::CCS;

BEGIN {
  use_ok('PDL::CCS::IO::FastRaw');
  use_ok('PDL::CCS::IO::FITS');
  use_ok('PDL::CCS::IO::MatrixMarket');
  use_ok('PDL::CCS::IO::LDAC');
  $| = 1;
}

##-- basic data
our $a = pdl(double, [
		      [10,0,0,0,-2],
		      [3,9,0,0,0],
		      [0,7,8,7,0],
		      [3,0,8,7,5],
		      [0,8,0,9,9],
		      [0,4,0,0,2],
		     ]);
our $ccs = $a->toccs();

##-- pdl equality
sub pdleq {
  my ($a,$b) = @_;
  return 0 if (!$a->ndims == $b->ndims || !all(pdl(long,[$a->dims])==pdl(long,[$b->dims])));
  if (UNIVERSAL::isa($a,'PDL::CCS::Nd')) {
    return 0 if ($a->_nnz_p != $b->_nnz_p);
    return all($a->_whichND==$b->_whichND) && all($a->_vals==$b->_vals);
  } else {
    return all($a==$b);
  }
}

##-- n..(n+3): i/o testing
sub iotest {
  my ($p, $file, $reader,$writer) = @_;
  my ($q);
  $reader = $p->can($reader) if (!ref($reader));
  $writer = $p->can($writer) if (!ref($writer));
  ok($writer->($p,"$TEST_DIR/$file"), "$file - write");
  ok(defined($q = $reader->("$TEST_DIR/$file")), "$file - read");
  is(ref($q), ref($p), "$file - ref");
  ok(pdleq($p,$q), "$file - data");

  ##-- unlink test data
  unlink($_) foreach (glob("$TEST_DIR/$file*"));
}

##-- 1..4 : raw
iotest($ccs, 'ccs.raw', qw(readfraw writefraw));

##-- 5..8 : fits
iotest($ccs, 'ccs.fits', qw(rfits wfits));

##-- 9..12 : mm/sparse
iotest($ccs, 'ccs.mm', qw(readmm writemm));

##-- 13..16 : mm/dense
iotest($a, 'dense.mm', qw(readmm writemm));

##-- 17..20 : ldac
iotest($a, 'ccs.ldac', qw(readldac writeldac));
