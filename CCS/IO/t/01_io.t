##-*- Mode: CPerl -*-
$TEST_DIR = './t';
#use lib qw(.. ../../.. ../blib/lib ../blib/arch); $TEST_DIR = '.'; # for debugging

use Test::More tests=>(5+(9*6));
use PDL;
use PDL::CCS;

BEGIN {
  use_ok('PDL::CCS::IO::Common');
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

##-- *6: i/o testing
sub iotest {
  my ($p, $file, $reader,$writer, $opts) = @_;
  my ($q);
  $reader = $p->can($reader) if (!ref($reader));
  $writer = $p->can($writer) if (!ref($writer));
  ok(defined($writer), "$file - writer sub");
  ok(defined($reader), "$file - reader sub");
  
  ok($writer->($p,"$TEST_DIR/$file",$opts), "$file - write");
  ok(defined($q = $reader->("$TEST_DIR/$file",$opts)), "$file - read");
  is(ref($q), ref($p), "$file - ref");
  ok(pdleq($p,$q), "$file - data");

  ##-- unlink test data
  #unlink($_) foreach (glob("$TEST_DIR/$file*"));
}

##-- x1 : raw
iotest($ccs, 'ccs.raw', qw(readfraw writefraw));

##-- x2 : fits
iotest($ccs, 'ccs.fits', qw(rfits wfits));

##-- x3-x5 : mm
do {
  iotest($ccs, 'ccs.mm', qw(readmm writemm));			##-- mm: sparse
  iotest($ccs, 'ccs.mm0', qw(readmm writemm), {header=>0});	##-- mm: sparse, no header
  iotest($a, 'dense.mm', qw(readmm writemm));			##-- mm: dense
};

##-- x6-x9 : ldac
do {
  iotest($ccs, 'ccs.ldac', qw(readldac writeldac));				##-- ldac: natural
  iotest($ccs, 'ccs.ldac0', qw(readldac writeldac), {header=>0});		##-- ldac: natural, no-header
  iotest($ccs, 'ccs.ldact', qw(readldac writeldac), {transpose=>1});		##-- ldac: transposed
  iotest($ccs, 'ccs.ldact0', qw(readldac writeldac), {header=>0,transpose=>1});	##-- ldac: transposed, no-header
}
