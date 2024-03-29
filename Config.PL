## File: PDL-CCS/Config.PL
## Description: user variables for PDL::CCS package

##-- load cached values?
if (0 && -e "./CCS/Config.pm") {
  require "./CCS/Config.pm";
  if ($@) {
    warn("$0: could not load cache data from './CCS/Config.pm': $@");
  }
  %cconfig = %PDL::CCS::Config::ccsConfig;
}

##--
## $val = cprompt($key, $message)
## $val = cprompt($key, $message, $default)
##   + sets $cconfig{$key}
sub cprompt {
  my ($key, $msg, $default)=@_;
  return $cconfig{$key} if (defined($cconfig{$key}));
  $default = '' if (!defined($default));
  my $answer = ExtUtils::MakeMaker::prompt("<Config.PL> $msg [$default] ? ");
  chomp($answer);
  return $cconfig{$key} = ($answer eq '' ? $default : $answer);
}

require PDL::Core;

$cconfig{INDX_CTYPE} = "PDL_Indx";
$cconfig{INDX_SIG}   = "indx";
$cconfig{INDX_FUNC}  = "indx";
$cconfig{INDX_TYPEDEF} = "typedef $cconfig{INDX_CTYPE} CCS_Indx;  /**< typedef for CCS indices */\n";
$cconfig{INDX_FUNCDEF} = "*ccs_indx = \\&PDL::$cconfig{INDX_FUNC}; ##-- typecasting for CCS indices\n";

##-- figure out what integer types we have available
require PDL::Types;
if (version->parse($PDL::VERSION) >= version->parse("2.037")) {
  ##-- integer types for b*over &c, PDL >= v2.037: elegant (and more correct)
  local $, = ' ';
  $cconfig{INT_TYPE_KEYS} = [map {$_->sym} grep {$_->integer} PDL::Types::types()];
  $cconfig{INT_TYPE_CHRS} = [map {$_->ppsym} grep {$_->integer} PDL::Types::types()];
  $cconfig{INT_TYPE_MAX_IONAME} = (grep {$_->integer} PDL::Types::types())[-1]->ioname;
} else {
  ##-- integer types for b*over &c, PDL < v2.037: functional (and mostly equivalent)
  $cconfig{INT_TYPE_KEYS} = [map {$_->{sym}}
			     sort {$a->{numval} <=> $b->{numval}}
			     grep {$_->{ppsym} =~ /^(?:[BSULQN]|LL|US)$/}
			     values %PDL::Types::typehash
			    ];
  $cconfig{INT_TYPE_CHRS} = [map {$_->{ppsym}} @PDL::Types::typehash{ @{$cconfig{INT_TYPE_KEYS}} }];

  ##-- PDL < v2.037 downcasts to 'indx' if available (but probably should use 'longlong' if it could)
  ## + behavior changed (for the better) apparently due to PDL commit #f892aeb4ae on Basic/Ufunc/ufunc.pd
  #$cconfig{INT_TYPE_MAX_IONAME} = $PDL::Types::typehash{$cconfig{INT_TYPE_KEYS}[-1]}{ioname}; ##-- -> longlong
  $cconfig{INT_TYPE_MAX_IONAME} = 'ccs_indx';
}

##-- save cache file
open(CONFIGPM,">./CCS/Config.pm")
  or die("$0: failed to open ./CCS/Config.pm for writing: $!");
print CONFIGPM <<'EOF';
## Automatically generated, remove to re-configure!

package PDL::CCS::Config;
use PDL qw();
our @ISA = qw(Exporter);
our (%ccsConfig);
our @EXPORT      = qw(ccs_indx);
our @EXPORT_OK   = ('%ccsConfig', 'ccs_indx');
our %EXPORT_TAGS = (config=>['%ccsConfig'], Func=>\@Export, default=>\@EXPORT, all=>\@EXPORT_OK);

EOF

##-- config hash
use Data::Dumper;
$Data::Dumper::Sortkeys=1; # reproducible order of hash keys
print CONFIGPM Data::Dumper->Dump([\%cconfig],['*ccsConfig']), "\n";

##-- type conversion sub
print CONFIGPM << "EOF";
\*PDL::ccs_indx = $cconfig{INDX_FUNCDEF}

1; ##-- be happy
EOF
close CONFIGPM;

1; ##-- return nicely
