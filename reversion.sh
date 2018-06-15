#!/bin/bash

## + requires perl-reversion from Perl::Version (debian package libperl-version-perl)
## + example call:
##    ./reversion.sh -bump -dryrun

pmfiles=(./CCS.pm `find CCS \( \( -name '*.pd' -o -name '*.pm' \) -a \! -path '*blib*' \) -print | grep -v Old`)

exec perl-reversion "$@" "${pmfiles[@]}"


