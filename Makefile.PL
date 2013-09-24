use ExtUtils::MakeMaker;
require "pdlmaker.plm";
pdlmaker_init();

##-- read in user variables
do "Config.PL";
die "$0: reading 'Config.PL' failed: $@" if ($@);


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
			    'PDL' => 0, #'2.4.2',
			    'PDL::VectorValued' => 0, #0.06
			   },
	      CONFIGURE_REQUIRES => {
				     'PDL' => 0,
				     'ExtUtils::MakeMaker'=>0,
				    },
	     );
