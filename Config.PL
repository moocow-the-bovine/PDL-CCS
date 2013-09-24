## File: PDL-CCS/Config.PL
## Description: user variables for PDL::CCS package

##-- load cache
if (-e "./CCS/Config.pm") {
  require "./CCS/Config.pm";
  if ($@) {
    warn("$0: could not load cache data from './CCS/Config.pm': $@");
  }
}

##--
## $val = cprompt($key, $message)
## $val = cprompt($key, $message, $default)
##   + sets $cconfig{$key}
sub cprompt {
  my ($key, $msg, $default)=@_;
  return $cconfig{$key} if (defined($cconfig{$key}));
  $default = '' if (!defined($default));
  print STDERR "$msg [$default] ? ";
  my $answer = <STDIN>;
  chomp($answer);
  return $cconfig{$key} = ($answer eq '' ? $default : $answer);
}

## $USE_PDL_INDX : whether to use built-in PDL_Indx type (PDL >= v2.007)
require PDL::Core;
if (defined(&PDL::indx)) {
  cprompt('USE_PDL_INDX', "use 64-bit indices via PDL_Indx if available (PDL >= v2.007)", 1);
} else {
  print STDERR "Info: 64-bit indices via PDL_Indx disabled for PDL < v2.007\n";
  $cconfig{USE_PDL_INDX} = 0;
}

if ($cconfig{USE_PDL_INDX}) {
  $cconfig{CCS_INDX_CTYPE} = "PDL_Indx";
  $cconfig{CCS_INDX_SIG}   = "indx";
  $cconfig{CCS_INDX_FUNC}  = "indx";
} else {
  $cconfig{CCS_INDX_CTYPE} = "PDL_Long";
  $cconfig{CCS_INDX_SIG}   = "int";
  $cconfig{CCS_INDX_FUNC}  = "long";
}

##-- save cache file
open(CACHE,">./CCS/Config.pm")
  or die("$0: failed to open ./CCS/Config.pm for writing: $!");
print CACHE <<EOF;
## Automatically generated, remove to re-configure!

package PDL::CCS::Config;
our (%Config);

EOF

while (($key,$val)=each(%cconfig)) {
  print CACHE "\$Config{$key} = q($val);\n";
}
print CACHE "\n1; ##-- be happy\n";
close CACHE;

1; ##-- return nicely
