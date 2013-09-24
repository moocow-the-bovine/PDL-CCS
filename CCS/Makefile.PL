use ExtUtils::MakeMaker;
require "../pdlmaker.plm";
pdlmaker_init();

WriteMakefile(
	      NAME=>'PDL::CCS::Nd',
	      VERSION_FROM => 'Version.pm',
	      #PM => { (map {$_=>"\$(INST_LIBDIR)/CCS/$_"} <*.pm>), },
	      DIR =>[
		     #'Old',
		     'Utils',
		     'Ufunc',
		     'Ops',
		     'MatrixOps',
		    ], ##-- debug#2
	      PREREQ_PM => {
			    'PDL' => 0,
			    'PDL::VectorValued' => 0, #'0.04',
			   },
	      CONFIGURE_REQUIRES => {
				     'PDL'=>0,
				     'ExtUtils::MakeMaker'=>0,
				    },
	     );
