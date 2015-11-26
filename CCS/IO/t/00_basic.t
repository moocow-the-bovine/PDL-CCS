##-*- Mode: CPerl -*-
use Test::More tests=>5;

######################### We start with some black magic to print on failure.
#use lib '../blib/lib','../blib/arch';

BEGIN {
  use_ok('PDL::CCS::IO::Common');
  use_ok('PDL::CCS::IO::FastRaw');
  use_ok('PDL::CCS::IO::FITS');
  use_ok('PDL::CCS::IO::MatrixMarket');
  use_ok('PDL::CCS::IO::LDAC');
  $| = 1;
}


######################### End of black magic.
