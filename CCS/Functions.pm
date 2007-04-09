## File: PDL::CCS::Functions.pm
## Author: Bryan Jurish <moocow@ling.uni-potsdam.de>
## Description: useful perl-level functions for PDL::CCS

package PDL::CCS::Functions;
use PDL::CCS::Version;
use PDL::CCS::Utils;
use PDL;
use strict;

our $VERSION = $PDL::CCS::VERSION;
our @ISA = ('PDL::Exporter');
our @EXPORT_OK =
  (
   ##-- Decoding
   qw(ccs_decode ccs_pointerlen),
  );
our %EXPORT_TAGS =
  (
   Func => [@EXPORT_OK],               ##-- respect PDL conventions (hopefully)
  );


##======================================================================
## pod: headers
=pod

=head1 NAME

PDL::CCS::Functions - Useful perl-level functions for PDL::CCS

=head1 SYNOPSIS

 use PDL;
 use PDL::CCS::Functions;

 ##---------------------------------------------------------------------
 ## ... stuff happens

=cut


##======================================================================
## Decoding
=pod

=head1 Decoding

=cut

##---------------------------------------------------------------
## Decoding: utils
=pod

=head2 ccs_pointerlen

=for sig

  Signature: (int ptr(N+1); int len(N))

Get number of non-missing values for each axis value from a CCS-encoded
offset pointer vector $ptr().

=cut

;#-- emacs

*PDL::ccs_pointerlen = \&ccs_pointerlen;
sub ccs_pointerlen {
  my ($ptr,$len) = @_;
  if (!defined($len)) {
    $len = $ptr->slice("1:-1") - $ptr->slice("0:-2");
  } else {
    $len .= $ptr->slice("1:-1");
    $len -= $ptr->slice("0:-2");
  }
  return $len;
}

##---------------------------------------------------------------
## Decoding: generic
=pod

=head2 ccs_decode

=for sig

  Signature: (int whichnd(Ndims,Nnz); nzvals(Nnz); missing(); \@Dims; [o]a(@Dims))


Decode a CCS-encoded matrix (no dataflow).

=cut

;#-- emacs

*PDL::ccs_decode = \&ccs_decode;
sub ccs_decode {
  my ($aw,$nzvals,$missing,$dims,$a) = @_;
  $missing = $PDL::undefval if (!defined($missing));
  $dims = [ map {$aw->slice("($_),")->max+1} (0..($aw->dim(0)-1))] if (!defined($dims));
  $a    = zeroes($nzvals->type, @$dims) if (!defined($a));
  $a   .= $missing;
  $a->indexND($aw) .= $nzvals;
  return $a;
}
