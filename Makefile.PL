use ExtUtils::MakeMaker;
require "pdlmaker.plm";
pdlmaker_init();

##-- put 'realclean_files' in a variable: avoid MakeMaker puking with:
##     ERROR from evaluation of .../ccsutils/Makefile.PL:
##     Modification of a read-only value attempted at /usr/share/perl/5.8/Carp/Heavy.pm line 45.
my $realclean_files = '*~ *.tmp README.txt README.html';
WriteMakefile(
	      NAME   =>'PDL::CCS',
	      AUTHOR => 'Bryan Jurish',
	      ABSTRACT =>'Sparse N-dimensional PDLs with compressed column storage',
	      ##
	      VERSION_FROM => 'CCS/Version.pm',
	      #VERSION => '1.17_001',
	      ##
	      #PM => { (map {$_=>"\$(INST_LIBDIR)/CCS/$_"} <*.pm>), },
	      DIR =>[
		     'CCS',
		    ],
	      realclean=>{ FILES=>$realclean_files, },
	      PREREQ_PM => {
			    'PDL' => '2.4.2',
			    'PDL::VectorValued' => '0.06', #0.06
			   },
	     );
