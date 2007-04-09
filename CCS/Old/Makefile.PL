use PDL::Core::Dev;
use ExtUtils::MakeMaker;
PDL::Core::Dev->import();

$package = ["Old.pd", 'Old', 'PDL::CCS::Old'];
%hash = pdlpp_stdargs($package);

$hash{AUTHOR} = 'Bryan Jurish';
$hash{ABSTRACT} = 'Compressed Column Storage for 2d PDLs';
$hash{PREREQ_PM}{PDL} = '2.4.2'; ##-- ensure PDL-ability
push(@{$hash{LIBS}}, '-lm');
$hash{DIR} = [];
#$hash{INC} .= '';
#$hash{OBJECT} .= '';
$hash{VERSION_FROM} = '../Version.pm';
$hash{realclean}{FILES} .= '*~ *.tmp README.txt';

WriteMakefile(%hash);

sub MY::postamble { pdlpp_postamble($package); }
package main;
