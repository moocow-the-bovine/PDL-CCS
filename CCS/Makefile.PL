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
			    PDL => '2.4.2',
			    'PDL::VectorValued' => undef, #'0.04',
			   },
	     );
