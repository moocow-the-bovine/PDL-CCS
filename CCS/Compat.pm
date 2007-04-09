## File: PDL::CCS::Compat.pm
## Author: Bryan Jurish <moocow@ling.uni-potsdam.de>
## Description: backwards-compatibility hacks for PDL::CCS

package PDL::CCS::Compat;
use PDL::CCS::Version;
use PDL::CCS::Functions;
use PDL::CCS::Utils;
use PDL;
use strict;

our $VERSION = $PDL::CCS::VERSION;
our @ISA = ('PDL::Exporter');
our @EXPORT_OK =
  (
   ##
   ##-- Encoding
   qw(ccs_encode_compat),
   qw(ccsencode ccsencode_nz ccsencodefull ccsencodefull_nz),
   qw(ccsencodea ccsencode_naz ccsencodefulla ccsencodefull_naz),
   qw(ccsencodeg ccsencode_g ccsencodefullg ccsencodefull_g),
   qw(ccsencodei ccsencode_i ccsencodefulli ccsencodefull_i),
   qw(ccsencodei2d ccsencode_i2d ccsencodefulli2d ccsencodefull_i2d),
   ##
   ##-- Decoding
   qw(_ccsdecodecols ccsdecodecols),
   qw(ccsdecode ccsdecodefull),
   qw(ccsdecode_g ccsdecodeg ccsdecodefull_g ccsdecodefullg),
  );
our %EXPORT_TAGS =
  (
   Func => [@EXPORT_OK],               ##-- respect PDL conventions (hopefully)
  );

##======================================================================
## pod: headers
=pod

=head1 NAME

PDL::CCS::Compat - Backwards-compatibility module for PDL::CCS

=head1 SYNOPSIS

 use PDL;
 use PDL::CCS::Compat;
 ##---------------------------------------------------------------------
 ## ... stuff happens

=cut

##======================================================================
## Encoding
=pod

=head1 Encoding

=cut

##---------------------------------------------------------------
## Encoding: generic

=pod

=head2 ccs_encode_compat

=for sig

  Signature: (int awhich(Nnz,2); avals(Nnz);
              int $N; int $M;
              int [o]ptr(N); int [o]rowids(Nnz); [o]nzvals(Nnz))

Generic wrapper for backwards-compatible ccsencode() variants.

=cut

*ccs_encode_compat = \&PDL::ccs_encode_compat;
sub PDL::ccs_encode_compat {
  my ($aw,$avals,$N,$M,$ptr,$rowids,$nzvals) = @_;

  $N = $aw->slice("(0),")->max+1 if (!defined($N));
  $M = $aw->slice("(1),")->max+1 if (!defined($M));

  my ($ptr1,$awi) = ccs_encode_pointers($aw->slice("(0),"), $N);

  if (defined($ptr))    { $ptr .= $ptr1->slice("0:-2"); }
  else                  { $ptr  = $ptr1->slice("0:-2"); $ptr->sever; }

  if (defined($rowids)) { $rowids .= $aw->slice("(1),")->index($awi); }
  else                  { $rowids  = $aw->slice("(1),")->index($awi); $rowids->sever; }

  if (defined($nzvals)) { $nzvals .= $avals->index($awi); }
  else                  { $nzvals  = $avals->index($awi); $nzvals->sever; }

  return ($ptr,$rowids,$nzvals);
}

##---------------------------------------------------------------
## Encoding: MISSING=zero
=pod

=head2 ccsencode

=head2 ccsencode_nz

=for sig

  Signature: (a(N,M); int [o]ptr(N); int [o]rowids(Nnz); [o]nzvals(Nnz))

Encodes matrix $a() in compressed column format, interpreting zeroes
as "missing" values.

Allocates output vectors if required.

=cut

*ccsencode
  = *ccsencodefull      = *ccsencodefull_nz
  = *PDL::ccsencode     = *PDL::ccsencode_nz
  = *PDL::ccsencodefull = *PDL::ccsencodefull_nz
  = \&ccsencode_nz;
sub ccsencode_nz {
  #my ($a,$ptr,$rowids,$nzvals) = @_;
  my $a  = shift->clump(-2);
  my $aw = $a->whichND;
  return ccs_encode_compat($aw, $a->indexND($aw), $a->dims, @_);
}

##---------------------------------------------------------------
## Encoding: MISSING=ZERO (approx)
=pod

=head2 ccsencodea

=head2 ccsencode_naz

=for sig

  Signature: (a(N,M); eps(); int [o]ptr(N); int [o]rowids(Nnz); [o]nzvals(Nnz))

Encodes matrix $a() in CCS format interpreting approximate zeroes as "missing" values.
This function is just like ccsencode_nz(), but uses the tolerance parameter
$eps() to determine which elements are to be treated as zeroes.

Allocates output vectors if required.

=cut

*ccsencodea
  = *ccsencodefulla      = *ccsencodefull_naz
  = *PDL::ccsencodea     = *PDL::ccsencode_naz
  = *PDL::ccsencodefulla = *PDL::ccsencodefull_naz
  = \&ccsencode_naz;
sub ccsencode_naz {
  #my ($a,$eps,$ptr,$rowids,$nzvals) = @_;
  my $a   = shift->clump(-2);
  my $eps = shift;
  my $aw  = $a->approx(0,$eps)->inplace->not->whichND;          ##-- FIXME: optimize
  return ccs_encode_compat($aw, $a->indexND($aw), $a->dims, @_);
}

##---------------------------------------------------------------
## Encoding: MISSING=BAD
=pod

=head2 ccsencodeg

=head2 ccsencode_g

=for sig

  Signature: (a(N,M); int [o]ptr(N); int [o]rowids(Nnz); [o]nzvals(Nnz))

Encodes matrix $a() in CCS format interpreting BAD values
as "missing".  Requires bad-value support built into PDL.

Allocates output vectors if required.

=cut

*ccsencodeg
  = *ccsencodefullg      = *ccsencodefull_g
  = *PDL::ccsencodeg     = *PDL::ccsencode_g
  = *PDL::ccsencodefullg = *PDL::ccsencodefull_g
  = \&ccsencode_g;
sub ccsencode_g {
  #my ($a,$ptr,$rowids,$nzvals) = @_;
  my $a     = shift->clump(-2);
  my $amask = zeroes(byte,$a->dims);
  $a->isgood($amask);
  my $aw    = $amask->whichND;
  return ccs_encode_compat($aw, $a->indexND($aw), $a->dims, @_);
}

##---------------------------------------------------------------
## Encoding: from flat index
=pod

=head2 ccsencode_i

=for sig

  Signature: (int ix(Nnz); nzvals(Nnz); int $N; int [o]ptr(N); int [o]rowids(Nnz); [o]nzvals_enc(Nnz))

General-purpose CCS encoding method for flat indices.
Encodes values $nzvals() from flat-index locations $ix() into a CCS matrix ($ptr(), $rowids(), $nzvals_enc()).

Allocates output vectors if required.

$N (~ $a-E<gt>dim(0)) must be specified.

=cut

*ccsencodei
  = *ccsencodefulli      = *ccsencodefull_i
  = *PDL::ccsencodei     = *PDL::ccsencode_i
  = *PDL::ccsencodefulli = *PDL::ccsencodefull_i
  = \&ccsencode_i;
sub ccsencode_i {
  #my ($iflat,$avals,$N,$ptr,$rowids,$nzvals) = @_;
  my ($iflat,$avals,$N) = splice(@_,0,3);
  my $aw = ($iflat % $N)->cat($iflat/$N)->xchg(0,1);
  return ccs_encode_compat($aw, $avals, $N, undef, @_);
}


##---------------------------------------------------------------
## Encoding: from 2d index
=pod

=head2 ccsencode_i2d

=for sig

  Signature: (
              int  xvals(Nnz)       ;
              int  yvals(Nnz)       ;
                  nzvals(Nnz)       ;
              int       $N          ;
              int [o]ptr(N)         ;
              int [o]rowids(Nnz)    ;
                  [o]nzvals_enc(Nnz);
             )

General-purpose encoding method.
Encodes values $nzvals() from 2d-index locations ($xvals(), $yvals()) in an $N-by-(whatever) PDL
into a CCS matrix $ptr(), $rowids(), $nzvals_enc().

Allocates output vectors if required.
If $N is omitted, it defaults to the maximum column index given in $xvals().

=cut

*ccsencodei2d
  = *ccsencodefulli2d      = *ccsencodefull_i2d
  = *PDL::ccsencodei2d     = *PDL::ccsencode_i2d
  = *PDL::ccsencodefulli2d = *PDL::ccsencodefull_i2d
  = \&ccsencode_i2d;
sub ccsencode_i2d {
  #my ($whichx,$whichy,$avals,$N,$ptr,$rowids,$nzvals) = @_;
  my ($whichx,$whichy,$avals,$N) = splice(@_, 0, 4);
  my $aw = $whichx->cat($whichy)->xchg(0,1);
  return ccs_encode_compat($aw, $avals, $N, undef, @_);
}

##======================================================================
## Decoding
=pod

=head1 Decoding

=cut

##---------------------------------------------------------------
## Decoding: column-wise

=pod

=head2 ccsdecodecols

=for sig

  Signature: (
              int    ptr    (N)  ;
              int    rowids (Nnz);
                     nzvals (Nnz);
              int    xvals  (I)  ; # default=sequence($N)
                     missing()   ; # default=0
                     M      ()   ; # default=rowids->max+1
                  [o]cols   (I,M); # default=new
              )

Extract dense columns from a CCS-encoded matrix (no dataflow).
Allocates output matrix if required.
If $a(N,M) was the dense source matrix for the CCS-encoding, and
if missing values are zeros, then the
following two calls are equivalent (modulo data flow):

  $cols = $a->dice_axis(1,$col_ix);
  $cols = ccsdecodecols($ptr,$rowids,$nzvals, $col_ix,0);

=cut

*_ccsdecodecols =
  = *PDL::_ccsencodecols = *PDL::_ccsdecodecols
  = \&ccsdecodecols;
sub ccsdecodecols {
  my ($ptr,$rowids,$nzvals, $coli,$missing,$M, $cols) = @_;
  $coli    = sequence(long,$ptr->dim(0)) if (!defined($coli));
  $coli    = pdl(long,$coli)             if (!ref($coli));
  my $ptr1 = zeroes(long,$ptr->nelem+1);
  $ptr1->slice("0:-2") .= $ptr;
  $ptr1->set(-1 => $nzvals->nelem);
  $M = $rowids->max+1 if (!defined($M));
  my ($ptrix,$nzix) = ccs_decode_pointer($ptr1,$coli);
  my $which = $ptrix->cat($rowids->index($nzix))->xchg(0,1);
  if (!defined($cols)) {
    $cols = ccs_decode($which, $nzvals->index($nzix), $missing, [$coli->nelem,$M]);
    $cols->sever; ##-- compat
  } else {
    ccs_decode($which, $nzvals->index($nzix), $missing, [$coli->nelem,$M], $cols);
  }
  return $cols;
}


##---------------------------------------------------------------
## Decoding: MISSING=0
=pod

=head2 ccsdecode

=for sig

  Signature: (int ptr(N); int rowids(Nnz); nzvals(Nnz); $M; [o]dense(N,M))

Decodes compressed column format vectors $ptr(), $rowids(), and $nzvals()
into dense output matrix $a().
Allocates the output matrix if required.

Note that if the original
matrix (pre-encoding) contained trailing rows with no nonzero elements,
such rows will not be allocated by this method (unless you specify either $M or $dense).
In such cases, you might prefer to call ccsdecodecols() directly.

=cut

*PDL::ccsdecodefull = \&ccsdecodefull; ##-- (int ptr(N); int rowids(Nnz); nzvals(Nnz); [o]dense(N,M))
sub ccsdecodefull { ccsdecodecols(@_[0,1,2], undef,0, @_[3..$#_]); }

*PDL::ccsdecode = \&ccsdecode;
sub ccsdecode {
  my ($ptr,$rowids,$nzvals, $M, $dense)=@_;
  if (!defined($dense)) {
    ##-- check for old calling convention (is $M a multi-dim PDL?)
    if (ref($M) && UNIVERSAL::isa($M, 'PDL') && $M->dim(0)==$ptr->dim(0)) {
      $dense = $M;
    } else {
      $M     = $rowids->max+1 if (!defined($M));
      $dense = zeroes($nzvals->type,$ptr->dim(0),$M);
    }
  }
  ccsdecodecols($ptr,$rowids,$nzvals, undef,0,$M, $dense);
  return $dense;
}


##---------------------------------------------------------------
## Decoding: MISSING=BAD
=pod

=head2 ccsdecode_g

=for sig

  Signature: (int ptr(N); int rowids(Nnz); nzvals(Nnz); $M; [o]dense(N,M))

Convenience method.
Like ccsdecode() but sets "missing" values to BAD.

=cut

*ccsdecodefullg = *PDL::ccsdecodefullg = *PDL::ccsdecodefull_g = \&ccsdecodefull_g;
sub ccsdecodefull_g {
  my $badval = pdl($_[2]->type,0)->setvaltobad(0);
  ccsdecodecols(@_[0,1,2], undef,$badval,undef,undef, @_[3..$#_]);
}

*ccsdecodeg = *PDL::ccsdecodeg = *PDL::ccsdecode_g = \&ccsdecode_g;
sub ccsdecode_g {
  my ($ptr,$rowids,$nzvals, $M, $dense)=@_;
  if (!defined($dense)) {
    ##-- check for old calling convention (is $M a multi-dim PDL?)
    if (ref($M) && UNIVERSAL::isa($M, 'PDL') && $M->dim(0)==$ptr->dim(0)) {
      $dense = $M;
    } else {
      $M     = $rowids->max+1 if (!defined($M));
      $dense = zeroes($nzvals->type,$ptr->dim(0),$M);
    }
  }
  my $badval = pdl($nzvals->type,0)->setvaltobad(0);
  ccsdecodecols($ptr,$rowids,$nzvals, undef,$badval,$M, $dense);
  return $dense;
}


##======================================================================
## Index Conversion: TODO
##======================================================================
=pod

=head1 Index Conversion

=cut

##------------------------------------------------------
## ccsitonzi() : index conversion: flat

##------------------------------------------------------
## ccsi2dtonzi() : index conversion: 2d

##------------------------------------------------------
## ccswhich*() : get indices


##------------------------------------------------------
## ccswhichfull() : which(): low-level


##======================================================================
## Lookup: TODO
##======================================================================

=pod

=head1 Lookup

=cut

##------------------------------------------------------
## ccsget() : lookup: flat

##------------------------------------------------------
## ccsget2d() : lookup: 2d

##======================================================================
## Operations: TODO
##======================================================================
=pod

=head1 Native Operations

=cut

##------------------------------------------------------
## ccstranspose() : transposition (convenience)

##------------------------------------------------------
## vector operations, column-vectors

#ccs${opname}_cv()


##------------------------------------------------------
## vector operations: row-vectors

#ccs${opname}_rv()

##------------------------------------------------------
## Ufuncs (accumulators)

#ccs${name}over()

