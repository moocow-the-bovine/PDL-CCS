#!/bin/bash

## + requires perl-reversion from Perl::Version (debian package libperl-version-perl)
## + example call:
##    ./reversion.sh -bump -dryrun

exec perl-reversion "$@" ./CCS.pm ./CCS/*.pm ./CCS/IO/FastRaw.pm ./CCS/*/*.pd


