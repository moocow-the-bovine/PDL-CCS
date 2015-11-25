##-*- Mode: CPerl -*-
use Test::More tests=>3;

######################### We start with some black magic to print on failure.
#use lib '../blib/lib','../blib/arch';

BEGIN {
  use_ok('PDL::CCS::IO::FastRaw');
  use_ok('PDL::CCS::IO::FITS');
  use_ok('PDL::CCS::IO::MatrixMarket');
  $| = 1;
}


######################### End of black magic.
