use ExtUtils::MakeMaker;
require "pdlmaker.plm";
pdlmaker_init();

##-- preqrequistes (for PREREQ_PM)
my %prereq = (
	      'PDL' => 0, #'2.4.2',
	      'PDL::VectorValued::Version' => '1.0.3',
	     );

##-- workaround for bogus cpantesters error messages
## + see http://sourceforge.net/p/pdl/mailman/message/34623263/ ("headaches with indirect PDL-related dependencies on cpantesters", 2015-11-16 13:16:44)
## + workaround: ensure working PDL::VectorValued at configure-time; return UNKNOWN otherwise
use version;
require 'PDL/VectorValued.pm'
  or die("$0: ERROR: could not require PDL/VectorValued.pm");
die("$0: ERROR: unsupported \$PDL::VectorValued::VERSION $PDL::VectorValued::VERSION (need at least $prereq{PDL::VectorValued::Version})")
  if (version->parse($PDL::VectorValued::VERSION) < version->parse($prereq{'PDL::VectorValued::Version'}));

##-- read in user variables
do "Config.PL";
die "$0: reading 'Config.PL' failed: $@" if ($@);

##-- put 'realclean_files' in a variable: avoid MakeMaker puking with:
##     ERROR from evaluation of .../ccsutils/Makefile.PL:
##     Modification of a read-only value attempted at /usr/share/perl/5.8/Carp/Heavy.pm line 45.
my $realclean_files = '*~ *.tmp README.txt README.html CCS/Config.pm';
WriteMakefile(
	      NAME   =>'PDL::CCS',
	      AUTHOR => 'Bryan Jurish',
	      ABSTRACT =>'Sparse N-dimensional PDLs with compressed column storage',
	      ##
	      VERSION_FROM => 'CCS/Version.pm',
	      ##
	      #PM => { (map {$_=>"\$(INST_LIBDIR)/CCS/$_"} <*.pm>), },
	      DIR =>[
		     'CCS',
		    ],
	      realclean=>{ FILES=>$realclean_files, },
	      PREREQ_PM => {%prereq},
	      CONFIGURE_REQUIRES => {
				     %prereq,
				     'ExtUtils::MakeMaker'=>0,
				     'Data::Dumper' => 0,
				    },
	     );

##-- avoid applying 'processPL' rules to 'Config.PL'
sub MY::processPL { return ''; }
