# -*- Mode: CPerl -*-
# t/05_binops.t

$TEST_DIR = './t';
#use lib qw(../blib/lib ../blib/arch); $TEST_DIR = '.'; # for debugging

# load common subs
use Test;
do "$TEST_DIR/common.plt";
use PDL;
use PDL::CCS::Nd;

BEGIN {
  my $N_MATOPS = 2;
  my $N_TESTS_PER_MATOP  = 8;
  my $N_RUNS_PER_BLOCK = 6;
  my $N_BLOCKS = 5;
  my $N_HACKS = 3;
  plan(tests=>(
	       $N_BLOCKS*$N_RUNS_PER_BLOCK*$N_TESTS_PER_MATOP*$N_MATOPS
	       +
	       $N_HACKS,
	      ),
       todo=>[]);
}

our ($BAD);

##--------------------------------------------------------------
## hacks

sub test_matmult2d_sdd {
  my ($lab,$a,$b,$az) = @_;  ##-- dense args
  $az = $a->toccs if (!defined($az));
  my $c = $a x $b;       ##-- dense output (desired)
  my $cz = $az->matmult2d_sdd($b);
  isok("${lab}:matmult2d_sdd:obj:missing=".($az->missing->sclr), all($c==$cz));
}
sub test_matmult2d_zdd {
  my ($lab,$a,$b,$az) = @_;  ##-- dense args
  $az = $a->toccs if (!defined($az));
  my $c = $a x $b;       ##-- dense output (desired)
  my $cz = $az->matmult2d_zdd($b);
  isok("${lab}:matmult2d_zdd:obj:missing=".($az->missing->sclr), all($c==$cz));
}
sub test_matmult2d_all {
  my ($M,$N,$O) = (2,3,4);
  my $a = sequence($M,$N);
  my $b = (sequence($O,$M)+1)*10;
  test_matmult2d_sdd('m0',$a,$b, $a->toccs);
  test_matmult2d_zdd('m0',$a,$b, $a->toccs);

  my $a1 = $a->pdl;
  $a1->where(($a%2)==0) .= 1;
  test_matmult2d_sdd('m1',$a,$b, $a->toccs(1));
}
test_matmult2d_all();

##--------------------------------------------------------------
## matrix operation test (manual swap)

##-- i..(i+8): test_matop($label, $op_name, $op_op_or_undef, $swap, $missing_val, $b,$bs)
##   + globals "$a" and "$abad" must always be defined
##   + "$as" is $a->toccs($missing_val);
##   + always tests
##   + for $swap==0
##     $PDL_FUNC->($a,$b) ~ $CCS_FUNC->($as,($b|$bs))
##     ($a OP $b)         ~ ($as OP ($bs|$b))
##   + for $swap==1
##     $PDL_FUNC->($b,$a) ~ $CCS_FUNC->($bs,($a|$as))
##     ($b OP $a)         ~ ($bs OP ($a|$as))
sub test_matop {
  my ($lab, $op_name, $op_op, $swap, $missing_val, $b,$bs) = @_;
  print "test_matop(lab=$lab, name=$op_name, op=", ($op_op||'NONE'), ", swap=$swap, missing=$missing_val)\n";

  my $pdl_func = PDL->can("${op_name}")
    or die("no PDL Ufunc ${op_name} defined!");
  my $ccs_func = PDL::CCS::Nd->can("${op_name}")
    or die("no CCS Ufunc PDL::CCS::Nd::${op_name} defined!");

  $missing_val = 0 if (!defined($missing_val));
  $missing_val = PDL->topdl($missing_val);
  if ($missing_val->isbad) { $a = $a->setbadif($abad); }
  else                     { $a->where($abad) .= $missing_val; $a->badflag(0); }

  my $a = $::a;
  $as = $a->toccs($missing_val);

  $b  = PDL->topdl($b);
  $bs = $b->toccs($missing_val) if (!defined($bs));
  if ($op_name eq 'matmult') {
    if ($lab eq 'mat.mat' && $b->ndims > 1 && $b->dim(1) != 1) {
      ##-- hack: mat.mat
      $b  = $b->xchg(0,1);
      $bs = $bs->xchg(0,1);
    }
    elsif ($lab eq 'mat.rv' && $b->ndims >= 1 && $b->dim(0)==$a->dim(0)) {
      ##-- hack: mat.rv --> rv.mat
      ($a,$as, $b,$bs) = ($b,$bs, $a,$as);
      $b  = $b->xchg(0,1);
      $bs = $bs->xchg(0,1);
      $swap = 0;
    }
    elsif ($lab eq 'mat.cv' && $b->ndims > 1 && $b->dim(0) == 1) {
      ##-- hack: mat.cv
      $a  = $a->xchg(0,1);
      $as = $as->xchg(0,1);
      $swap = 0;
    }
    elsif ($lab eq 'rv.cv') {
      $a  = $a->xchg(0,1);
      $as = $as->xchg(0,1);
      $b  = $b->xchg(0,1);
      $bs = $bs->xchg(0,1);
      $swap = 0;
    }
  }


  ##-- test: function syntax
  my ($dense_rc,$ccs_bs,$ccs_b);
  if (!$swap) {
    $pdl_func->($a,  $b, $dense_rc=null);
    $ccs_bs   = $ccs_func->($as, $bs);
    $ccs_b    = $ccs_func->($as, $b);
  } else {
    $pdl_func->($b,  $a, $dense_rc=null);
    $ccs_bs   = $ccs_func->($bs, $as);
    $ccs_b    = $ccs_func->($bs, $a);
  }

  isok("$lab:${op_name}:func:b=sparse:missing=$missing_val:swap=$swap:type",
       $dense_rc->type==$ccs_bs->type);
  isok("$lab:${op_name}:func:b=sparse:missing=$missing_val:swap=$swap:vals",
       all( matchpdl($dense_rc, $ccs_bs->decode) ));
  isok("$lab:${op_name}:func:b=dense:missing=$missing_val:swap=$swap:type",
       $dense_rc->type==$ccs_b->type);
  isok("$lab:${op_name}:func:b=dense:missing=$missing_val:swap=$swap:vals",
       all( matchpdl($dense_rc, $ccs_b->decode) ));

  if (defined($op_op)) {
    if (!$swap) {
      eval "\$dense_rc = (\$a  $op_op \$b);";
      eval "\$ccs_bs   = (\$as $op_op \$bs);";
      eval "\$ccs_b    = (\$as $op_op \$b);";
    } else {
      eval "\$dense_rc = (\$b  $op_op \$a);";
      eval "\$ccs_bs   = (\$bs $op_op \$as);";
      eval "\$ccs_b    = (\$bs $op_op \$a);";
    }
    isok("$lab:${op_name}:op=$op_op:b=sparse:missing=$missing_val:swap=$swap:type",
	 $dense_rc->type==$ccs_bs->type);
    isok("$lab:${op_name}:op=$op_op:b=sparse:missing=$missing_val:swap=$swap:vals",
	 all( matchpdl($dense_rc,$ccs_bs->decode) ));
    isok("$lab:${op_name}:op=$op_op:b=dense:missing=$missing_val:swap=$swap:type",
	 $dense_rc->type==$ccs_b->type);
    isok("$lab:${op_name}:op=$op_op:b=dense:missing=$missing_val:swap=$swap:vals",
	 all( matchpdl($dense_rc,$ccs_b->decode) ));
  } else {
    isok("$lab:${op_name}:op=NONE:b=sparse:missing=$missing_val:swap=$swap:type (dummy)", 1);
    isok("$lab:${op_name}:op=NONE:b=sparse:missing=$missing_val:swap=$swap:vals (dummy)", 1);
    isok("$lab:${op_name}:op=NONE:b=dense:missing=$missing_val:swap=$swap:type  (dummy)", 1);
    isok("$lab:${op_name}:op=NONE:b=dense:missing=$missing_val:swap=$swap:vals  (dummy)", 1);
  }
}

my @matops = (
	      ##-- Matrix operations
	      'inner',
	      [qw(matmult x)],
	     );

##-- Block 1 : mat * mat (rotated)
$b = $a->flat->rotate(1)->reshape($a->dims);
foreach $missing (0,127,$BAD) {   ##-- *3
  foreach $swap (0,1) {           ##-- *2
    foreach $op (@matops) {       ##-- *NMATOPS
      if (ref($op)) { test_matop('mat.mat', @$op,        $swap, $missing, $b); }
      else          { test_matop('mat.mat', $op, undef,  $swap, $missing, $b); }
    }
  }
}

##-- Block 2 : mat * scalar
$b = PDL->topdl(42);
foreach $missing (0,127,$BAD) {   ##-- *3
  foreach $swap (0,1) {           ##-- *2
    foreach $op (@matops) {       ##-- *NMATOPS
      if (ref($op)) { test_matop('mat.sclr', $op->[0], $op->[1], $swap, $missing, $b); }
      else          { test_matop('mat.sclr', $op,      undef,    $swap, $missing, $b); }
    }
  }
}

##-- Block 3 : mat * row
$b  = sequence($a->dim(0),1)+1;
foreach $missing (0,127,$BAD) {   ##-- *3
  foreach $swap (0,1) {           ##-- *2
    foreach $op (@matops) {         ##-- *NMATOPS
      if (ref($op)) { test_matop('mat.rv', $op->[0], $op->[1], 1,     $missing, $b); } ##-- hack
      else          { test_matop('mat.rv', $op,      undef,    $swap, $missing, $b); }
    }
  }
}

##-- Block 4 : mat * col
$b  = sequence(1,$a->dim(1))+1;
$bs = $b->flat->toccs->dummy(0,1);
foreach $missing (0,127,$BAD) {   ##-- *3
  foreach $swap (0,1) {           ##-- *2
    foreach $op (@matops) {       ##-- *NMATOPS
      if (ref($op)) { test_matop('mat.cv', $op->[0], $op->[1], $swap, $missing, $b,$bs); }
      else          { test_matop('mat.cv', $op,      undef,    $swap, $missing, $b,$bs); }
    }
  }
}

##-- Block 5 : col * row
my @save = ($a,$abad);
$b  = sequence(1,$a->dim(1))+1;
$bs = $b->flat->toccs->dummy(0,1);
$a  = sequence($a->dim(0),1);
$abad = ($a==0);
foreach $missing (0,127,$BAD) {   ##-- *3
  foreach $swap (0,1) {           ##-- *2
    foreach $op (@matops) {       ##-- *NMATOPS
      if (ref($op)) { test_matop('rv.cv', $op->[0], $op->[1], $swap, $missing, $b,$bs); }
      else          { test_matop('rv.cv', $op,      undef,    $swap, $missing, $b,$bs); }
    }
  }
}

($a,$abad) = @save;

print "\n";
# end of t/*.t

