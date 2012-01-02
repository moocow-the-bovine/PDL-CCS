use PDL::Core::Dev;
use ExtUtils::MakeMaker;
PDL::Core::Dev->import();
require "../../pdlmaker.plm";

$package = ["ccsops.pd", 'Ops', 'PDL::CCS::Ops'];
%hash = pdlmaker_init($package);

$hash{AUTHOR} = 'Bryan Jurish';
$hash{ABSTRACT} = 'Low-level binary operations for compressed storage sparse PDLs';
$hash{VERSION_FROM} = '../Version.pm';
$hash{PREREQ_PM}{PDL} = '2.4.2'; ##-- ensure PDL-ability
push(@{$hash{LIBS}}, '-lm');
$hash{DIR} = [];
#$hash{INC} .= '';
#$hash{OBJECT} .= '';
$hash{realclean}{FILES} .= '*~ *.tmp README.txt';

#my $pmfile = $package[0];
#$pmfile =~ s/\.pd$/\.pm/;
#$hash{PM}{$pmfile} = "\$(INST_LIBDIR)/CCS/$pmfile";

WriteMakefile(%hash);
