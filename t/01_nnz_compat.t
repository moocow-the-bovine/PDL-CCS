# -*- Mode: CPerl -*-
# t/01_nnz_compat.t: test n nonzeros (compat)

use Test::More tests => 5;

$TEST_DIR = './t';
#use lib qw(../blib/lib ../blib/arch); $TEST_DIR = '.'; # for debugging

# load common subs
do "$TEST_DIR/common.plt";
use PDL;
use PDL::CCS;

## 1--4: test nnz
my $p = pdl(double, [ [0,1,2], [0,0,1e-7], [0,1,0], [1,1,1] ]);
isok("nnz(0)",     $p->slice(",(0)")->nnz, 2);
isok("nnz(flat)",  $p->flat->nnz, 7);
isok("nnza(flat,1e-8)", $p->flat->nnza(1e-8), 7);
isok("nnza(flat,1e-5)", $p->flat->nnza(1e-5), 6);
isok("nnza(flat:int,1)",     $p->flat->long->nnza(1), 1);

print "\n";
# end of t/?.t

