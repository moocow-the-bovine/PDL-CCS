use PDL::Core::Dev;
use ExtUtils::MakeMaker;
PDL::Core::Dev->import();

$package = ["CCS.pd",CCS,PDL::CCS];
%hash = pdlpp_stdargs($package);

$hash{AUTHOR} = 'Bryan Jurish';
$hash{ABSTRACT} = 'Compressed Column Storage for 2d PDLs';
$hash{PREREQ_PM}{PDL} = '2.4.2'; ##-- ensure PDL-ability
push(@{$hash{LIBS}}, '-lm');
#$hash{INC} .= '';
$hash{DIR} = [];
#$hash{OBJECT} .= '';
$hash{VERSION_FROM} = 'CCS.pd';
$hash{realclean}{FILES} .= '*~ *.tmp README.txt';

WriteMakefile(%hash);

sub MY::depend {
  package MY;
  my $inherited = shift->SUPER::depend(@_);
  $inherited .= (
		 "\n\n"
		 ."dist: README.txt\n"
		 ."\n"
		 ."create_distdir: README.txt\n"
		 ."\n"
		 ."distcheck: README.txt\n"
		 ."\n"
		 ."manicheck: README.txt\n"
		 ."\n"
		);
}

sub MY::special_targets {
  package MY;
  my $inherited = shift->SUPER::special_targets(@_);
  #$inherited =~ s/^(\.SUFFIXES\s*:.*)$/$1 .pod .pm .man .txt .html;/
  $inherited .=
"
.SUFFIXES: .pm .pod .man .txt .html

.pm.html:
\tpod2html --outfile \$@ \$<

.pm.txt:
\tpod2text \$< \$@

.pod.html:
\tpod2html --outfile \$@ \$<

.pod.txt:
\tpod2text \$< \$@

#README.html: README.pod
#\tpod2html --outfile \$@ \$<
#
#README.txt: README.pod
#\tpod2text \$< \$@
"
}


sub MY::top_targets {
  package MY;
  my $inherited = shift->SUPER::top_targets(@_);
  $inherited;
}

sub MY::postamble {
  pdlpp_postamble($package);
}
