# -*- Mode: CPerl -*-
# t/01_nnz.t: test n nonzeros

$TEST_DIR = './t';
#use lib qw(../blib/lib ../blib/arch); $TEST_DIR = '.'; # for debugging

# load common subs
use Test;
do "$TEST_DIR/common.plt";
use PDL;
use PDL::CCS::Utils;

BEGIN { plan tests=>4, todo=>[]; }

## 1--4: test nnz
$p = pdl(double, [ [0,1,2], [0,0,1e-7], [0,1,0], [1,1,1] ]);
isok("nnz(0)",     $p->slice(",(0)")->nnz==2);
isok("nnz(flat)",  $p->flat->nnz==7);
isok("nnza(flat)", $p->flat->nnza(1e-6)==6);
isok("nnza(flat,.5)", $p->flat->nnza(.5)==7);


print "\n";
# end of t/01_nnz.t

