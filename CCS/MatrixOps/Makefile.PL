use PDL::Core::Dev;
use ExtUtils::MakeMaker;
PDL::Core::Dev->import();

$package = ["ccsmatops.pd", 'MatrixOps', 'PDL::CCS::MatrixOps'];
%hash = pdlpp_stdargs($package);

$hash{AUTHOR} = 'Bryan Jurish';
$hash{ABSTRACT} = 'Low-level matrix operations for compressed storage sparse PDLs';
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

sub MY::postamble { pdlpp_postamble($package); }