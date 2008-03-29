# -*- Mode: CPerl -*-
# t/01_ufunc.t

$TEST_DIR = './t';
#use lib qw(../blib/lib ../blib/arch); $TEST_DIR = '.'; # for debugging

# load common subs
use Test;
do "$TEST_DIR/common.plt";
use PDL;
use PDL::CCS::Ufunc;
use PDL::VectorValued;

BEGIN { plan tests=>104, todo=>[]; }

##-- basic data
our $a = pdl(double, [
		      [10,0,0,0,-2],
		      [3,9,0,0,0],
		      [0,7,8,7,0],
		      [3,0,8,7,5],
		      [0,8,0,9,9],
		      [0,4,0,0,2],
		     ]);

our $agood    = ($a!=0);
our $abad     = !$agood;
our $awhich   = $a->whichND;
our $awhich1  = $awhich->slice("(1)")->qsort->slice("*1,");
our $awhich1i = $awhich->slice("(1)")->qsorti;
our $avals    = $a->indexND($awhich)->index($awhich1i);

sub matchpdl {
  my ($a,$b) = map {$_->setnantobad} @_[0,1];
  return ($a==$b)->setbadtoval(0) | ($a->isbad & $b->isbad);
}


##-- i..(i+2): test_ufunc($pdl_ufunc_name, $ccs_ufunc_name, $missing_val)
sub test_ufunc {
  my ($pdl_ufunc_name, $ccs_ufunc_name, $missing_val) = @_;
  print "test_ufunc($pdl_ufunc_name, $ccs_ufunc_name, $missing_val)\n";

  my $ccs_ufunc = PDL->can("ccs_accum_${ccs_ufunc_name}")
    or die("no CCS Ufunc ccs_accum_${ccs_ufunc_name} defined!");
  my $pdl_ufunc = PDL->can("${pdl_ufunc_name}")
    or die("no PDL Ufunc ${pdl_ufunc_name} defined!");

  $missing_val = 0 if (!defined($missing_val));
  $missing_val = PDL->topdl($missing_val);
  if ($missing_val->isbad) { $a = $a->setbadif($abad); }
  else                     { $a->where($abad) .= $missing_val; $a->badflag(0); }

  my $dense_rc = $pdl_ufunc->($a);
  my ($which_rc,$nzvals_rc) = $ccs_ufunc->($awhich1, $avals, $missing_val, $a->dim(0));
  my $decoded_rc = $dense_rc->zeroes;
  $decoded_rc   .= $missing_val;
  $decoded_rc->indexND($which_rc) .= $nzvals_rc;

  isok("${pdl_ufunc_name}:missing=$missing_val:type", $nzvals_rc->type==$dense_rc->type);
  isok("${pdl_ufunc_name}:missing=$missing_val:vals", all($decoded_rc==$dense_rc));
}

our $BAD = pdl(0)->setvaltobad(0);

foreach $missing (0,1,31,$BAD) { ## *4
  foreach $pdl_ufunc_name (
			   qw(sumover prodover dsumover dprodover),  ## *13
			   qw(andover orover bandover borover),
			   qw(maximum minimum),
			   qw(nbadover ngoodover), #nnz
			   qw(average),
			  )
    {
      my $ccs_ufunc_name = $pdl_ufunc_name;
      $ccs_ufunc_name =~ s/over$//;
      test_ufunc($pdl_ufunc_name, $ccs_ufunc_name, $missing);
    }
}

print "\n";
# end of t/*.t

