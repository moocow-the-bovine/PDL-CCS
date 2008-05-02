# -*- Mode: CPerl -*-
# t/03_ufuncs.t

$TEST_DIR = './t';
#use lib qw(../blib/lib ../blib/arch); $TEST_DIR = '.'; # for debugging

# load common subs
use Test;
do "$TEST_DIR/common.plt";
use PDL;
use PDL::CCS::Nd;

BEGIN { plan tests=>2*4*15, todo=>[]; }

##--------------------------------------------------------------
## basic test

##-- i..(i+2): test_ufunc($ufunc_name, $missing_val)
sub test_ufunc {
  my ($ufunc_name, $missing_val) = @_;
  print "test_ufunc($ufunc_name, $missing_val)\n";

  my $pdl_ufunc = PDL->can("${ufunc_name}")
    or die("no PDL Ufunc ${ufunc_name} defined!");
  my $ccs_ufunc = PDL::CCS::Nd->can("${ufunc_name}")
    or die("no CCS Ufunc PDL::CCS::Nd::${ufunc_name} defined!");

  $missing_val = 0 if (!defined($missing_val));
  $missing_val = PDL->topdl($missing_val);
  if ($missing_val->isbad) { $a = $a->setbadif($abad); }
  else                     { $a->where($abad) .= $missing_val; $a->badflag(0); }

  my $ccs      = $a->toccs($missing_val);
  my $dense_rc = $pdl_ufunc->($a);
  my $ccs_rc   = $ccs_ufunc->($ccs);

  if ($ufunc_name =~ /_ind/) {
    ##-- hack: adjust $dense_rc for maximum_ind, minimum_ind
    $dense_rc->where( $a->index2d($dense_rc,sequence($a->dim(1))) == $missing ) .= -1;
  }

  isok("${ufunc_name}:missing=$missing_val:type", $dense_rc->type==$ccs_rc->type);
  isok("${ufunc_name}:missing=$missing_val:vals", all( matchpdl($ccs_rc->decode, $dense_rc) ));
}

our ($BAD);
foreach $missing (0,1,255,$BAD) { ##-- *4
  foreach $ufunc (
		  qw(sumover prodover dsumover dprodover),  ## *15
		  qw(andover orover bandover borover),
		  qw(maximum minimum),
		  qw(maximum_ind minimum_ind),
		  qw(nbadover ngoodover), #nnz
		  qw(average),
		 )
    {
      test_ufunc($ufunc,$missing);
    }
}

print "\n";
# end of t/*.t

