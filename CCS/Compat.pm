## File: PDL::CCS::Compat.pm
## Author: Bryan Jurish <moocow@ling.uni-potsdam.de>
## Description: backwards-compatibility hacks for PDL::CCS

package PDL::CCS::Compat;
use PDL;
use PDL::VectorValued;
use PDL::CCS::Version;
use PDL::CCS::Functions;
use PDL::CCS::Utils;
use PDL::CCS::Ufunc;
use PDL::CCS::Ops;
use strict;

our $VERSION = $PDL::CCS::VERSION;
our @ISA = ('PDL::Exporter');
our @ccs_binops = (qw(plus minus mult divide modulo power),
		   qw(gt ge lt le eq ne spaceship),
		   qw(and2 or2 xor shiftleft shiftright),
		  );
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
   ##
   ##-- Indexing
   qw(ccsiNDtonzi ccsi2dtonzi ccsitonzi),
   qw(ccswhichND ccswhich2d ccswhichfull ccswhich),
   qw(ccstranspose ccstransposefull),
   ##
   ##-- Lookup
   qw(ccsget ccsget2d),
   ##
   ##-- Operations
   (map {("ccs${_}_cv","ccs${_}_rv")} (@ccs_binops,qw(add diff))),
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

  Signature: (int awhich(2,Nnz); avals(Nnz);
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

*_ccsdecodecols
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
## Index Conversion
##======================================================================
=pod

=head1 Index Conversion

=cut

##------------------------------------------------------
## ccsiNDtonzi() : index conversion: N-dimensional

=pod

=for sig

  Signature: (int ptr(N); int rowids(Nnz); int ind(2,I); int missing(); int [o]nzix(I))

=head2 ccsiNDtonzi

Convert N-dimensional index values $ind() appropriate for a dense matrix (N,M)
into indices $nzix() appropriate for the $rowids() and/or $nzvals() components
of the CCS-encoded matrix ($ptr(),$rowids(),$nzvals()).
Missing values are returned in $nzix() as $missing().

=cut

*PDL::ccsiNDtonzi = \&ccsiNDtonzi;
sub ccsiNDtonzi {
  my ($ptr,$rowids,$ind, $missing, $nzix) = @_;
  my ($ptri,$ptrnzi) = ccs_decode_pointer($ptr->append($rowids->nelem));
  my $ccswnd         = $ptri->cat($rowids->index($ptrnzi))->xchg(0,1)->vv_qsortvec;
  $nzix              = $ind->vsearchvec($ccswnd);
  my $nzix_mask      = ($ind==$ccswnd->dice_axis(1,$nzix))->andover;
  $nzix_mask->inplace->not;
  $nzix->where($nzix_mask) .= $missing;
  return $nzix;
}

##------------------------------------------------------
## ccsi2dtonzi() : index conversion: 2d

=pod

=head2 ccsi2dtonzi

=for sig

  Signaure: (int ptr(N); int rowids(Nnz); int col_ix(I); int row_ix(I); int missing(); int [o]nzix(I))

Convert 2d index values $col_ix() and $row_ix() appropriate for a dense matrix (N,M)
into indices $nzix() appropriate for the $rowids() and/or $nzvals() components
of the CCS-encoded matrix ($ptr(),$rowids(),$nzvals()).
Missing values are returned in $nzix() as $missing().

=cut

*PDL::ccsi2dtonzi = \&ccsi2dtonzi;
sub ccsi2dtonzi {
  my ($ptr,$rowids,$xi,$yi, $missing, $nzix) = @_;
  return ccsiNDtonzi($ptr,$rowids, $xi->cat($yi)->xchg(0,1), $missing,$nzix);
}


##------------------------------------------------------
## ccsitonzi() : index conversion: flat

=pod

=for sig

  Signature: (int ptr(N); int rowids(Nnz); int ix(I); int missing(); int [o]nzix(I))

=head2 ccsitonzi

Convert flat index values $ix() appropriate for a dense matrix (N,M)
into indices $nzix() appropriate for the $rowids() and/or $nzvals() components
of the CCS-encoded matrix ($ptr(),$rowids(),$nzvals()).
Missing values are returned in $nzix() as $missing().

=cut

*PDL::ccsitonzi = \&ccsitonzi;
sub ccsitonzi {
  my ($ptr,$rowids,$ix, $missing, $nzix) = @_;
  my $dummy    = pdl(byte,0)->slice("*".($ptr->dim(0)).",*".($rowids->max+1));
  my ($xi,$yi) = $dummy->one2nd($ix);
  return ccsiNDtonzi($ptr,$rowids, $xi->cat($yi)->xchg(0,1), $missing,$nzix);
}



##------------------------------------------------------
## ccswhichND: get indices (N-dimensional)

=pod

=head2 ccswhichND

=head2 ccswhich2d

=head2 ccswhichfull

=for sig

  Signature: (int ptr(N); int rowids(Nnz); nzvals(Nnz); int [o]which_cols(Nnz); int [o]which_rows(Nnz)',

In scalar context, returns concatenation of $which_cols() and $which_rows(),
similar to the builtin whichND().  Note however that ccswhichND() may return
its index PDLs sorted in a different order than the builtin whichND() method
for dense matrices.  Use the qsort() or qsorti() methods if you need sorted index PDLs.

=cut

*ccswhich2d = *PDL::which2d = *PDL::ccswhichND
  = *ccswhichfull = *PDL::ccswhichfull
  = \&ccswhichND;
sub ccswhichND {
  my ($ptr,$rowids,$nzvals, $which_cols,$which_rows) = @_;
  my ($ptrnzi);
  ($which_cols,$ptrnzi) = ccs_decode_pointer($ptr->append($rowids->nelem),
					     sequence(long, $ptr->nelem),
					     $which_cols
					    );
  $which_rows  = zeroes(long, $rowids->nelem) if (!defined($which_rows));
  $which_rows .= $rowids->index($ptrnzi);
  return wantarray ? ($which_cols,$which_rows) : $which_cols->cat($which_rows)->xchg(0,1);
}


##------------------------------------------------------
## ccswhich(): get indices (flat)

=head2 ccswhich

=for sig

  Signature: (int ptr(N); int rowids(Nnz); nzvals(Nnz); int [o]which(Nnz); int [t]wcols(Nnz)',

Convenience method.
Calls ccswhichfull(), and scales the output PDLs to correspond to a flat enumeration.
The output PDL $which() is B<not> guaranteed to be sorted in any meaningful order.
Use the qsort() method if you need sorted output.

=cut

*PDL::ccswhich = \&ccswhich;
sub ccswhich {
  my ($ptr,$rowids,$nzvals, $which, $wcols) = @_;
  my $nnz = $rowids->dim(0);
  $which = zeroes(long,$nnz) if (!defined($which));
  $wcols = zeroes(long,$nnz) if (!defined($wcols));
  ccswhichfull($ptr,$rowids,$nzvals, $wcols,$which);
  $which *= $ptr->dim(0);
  $which += $wcols;
  return $which;
}

##------------------------------------------------------
## ccstranspose() : transposition (convenience)

=pod

=head2 ccstranspose

=for sig

  Signature: (int ptr(N); int rowids(Nnz); nzvals(Nnz); int [o]ptrT(M); int [o]rowidsT(Nnz); [o]nzvalsT(Nnz)',

Transpose a compressed matrix.

=cut

*ccstransposefull = *PDL::ccstransposefull = *PDL::ccstranspose = \&ccstranspose;
sub ccstranspose {
  my ($ptr,$rowids,$nzvals, $ptrT,$rowidsT,$nzvalsT)=@_;
  my $N   = $ptr->dim(0);
  my $M   = defined($ptrT) ? $ptrT->dim(0) : $rowids->max+1;
  my $wnd = ccswhichND($ptr,$rowids,$nzvals)->slice("1:0,");
  return ccs_encode_compat($wnd,$nzvals,$M,$N, $ptrT,$rowidsT,$nzvalsT);
}


##======================================================================
## Lookup
##======================================================================

=pod

=head1 Lookup

=cut

##------------------------------------------------------
## ccsget2d() : lookup: 2d

=pod

=head2 ccsget2d

=for sig

  Signature: (int ptr(N); int rowids(Nnz); nzvals(Nnz); int xvals(I); int yvals(I); missing(); [o]ixvals(I))

Lookup values in a CCS-encoded PDL by 2d source index (no dataflow).
Pretty much like ccsi2dtonzi(), but returns values instead of indices.
If you know that your index PDLs $xvals() and $yvals() do not refer to any missing
values in the CCS-encoded matrix,
then the following two calls are equivalent (modulo dataflow):

  $ixvals =                ccsget2d   ($ptr,$rowids,$nzvals, $xvals,$yvals,0);
  $ixvals = index($nzvals, ccsi2dtonzi($ptr,$rowids,         $xvals,$yvals,0));

The difference is that only the second incantation will cause subsequent changes to $ixvals
to be propagated back into $nzvals.

=cut

*PDL::ccsget2d = \&ccsget2d;
sub ccsget2d {
  my ($ptr,$rowids,$nzvals, $xi,$yi, $missing, $ixnzvals) = @_;
  my $nzi        = ccsi2dtonzi($ptr,$rowids, $xi,$yi, -1);
  my $nzi_isgood = ($nzi != -1);
  $ixnzvals      = zeroes($nzvals->type, $xi->nelem) if (!defined($ixnzvals));
  $ixnzvals->where( $nzi_isgood) .= $nzvals->index($nzi->where($nzi_isgood));
  $ixnzvals->where(!$nzi_isgood) .= $missing;
  return $ixnzvals;
}


##------------------------------------------------------
## ccsget() : lookup: flat

=pod

=head2 ccsget

=for sig

  Signature: (int ptr(N); int rowids(Nnz); nzvals(Nnz); int ix(I); missing(); [o]ixvals(I))

Lookup values in a CCS-encoded PDL by flat source index (no dataflow).
Pretty much like ccsitonzi(), but returns values instead of indices.
If you know that your index PDL $ix() does not refer to any missing
values in the CCS-encoded matrix,
then the following two calls are equivalent (modulo dataflow):

  $ixvals =                ccsget   ($ptr,$rowids,$nzvals, $ix,0);
  $ixvals = index($nzvals, ccsitonzi($ptr,$rowids,         $ix,0))

The difference is that only the second incantation will cause subsequent changes to $ixvals
to be propagated back into $nzvals.

=cut

*PDL::ccsget = \&ccsget;
sub ccsget {
  my ($ptr,$rowids,$nzvals, $ix, $missing, $ixnzvals) = @_;
  my $nzi        = ccsitonzi($ptr,$rowids, $ix,-1);
  my $nzi_isgood = ($nzi != -1);
  $ixnzvals      = zeroes($nzvals->type, $ix->nelem) if (!defined($ixnzvals));
  $ixnzvals->where( $nzi_isgood) .= $nzvals->index($nzi->where($nzi_isgood));
  $ixnzvals->where(!$nzi_isgood) .= $missing;
  return $ixnzvals;
}


##======================================================================
## Vector Operations
##======================================================================
=pod

=head1 Vector Operations

=cut

##======================================================================
## Vector Operations: Columns
##======================================================================
=pod

=head2 ccs${OP}_cv

=for sig

  Signature: (int ptr(N); int rowids(Nnz); nzvals_in(Nnz);  colvec(M);  [o]nzvals_out(Nnz))

Column-vector operation ${OP} on CCS-encoded PDL.

Should do something like the following
(without decoding the CCS matrix):

 ($colvec ${OP} ccsdecode(\$ptr,\$rowids,\$nzvals))->ccsencode;

Missing values in the CCS-encoded PDL are not affected by this operation.

${OP} is one of the following:

 plus       ##-- addition    (alias: 'add')
 minus      ##-- subtraction (alias: 'diff')
 mult       ##-- multiplication (NOT matrix-multiplication)
 divide     ##-- division    (alias: 'div')
 modulo     ##-- modulo
 power      ##-- potentiation

 gt         ##-- greater-than
 ge         ##-- greater-than-or-equal
 lt         ##-- less-than
 le         ##-- less-than-or-equal
 eq         ##-- equality
 ne         ##-- inequality
 spaceship  ##-- 3-way comparison

 and2       ##-- binary AND
 or2        ##-- binary OR
 xor        ##-- binary XOR
 shiftleft  ##-- left-shift
 shiftright ##-- right-shift

=cut

sub ccs_binop_compat_cv {
  my $ccsop = shift;
  return sub { $ccsop->(@_[1,2,3,4]) };
}
foreach my $op (@ccs_binops) {
  eval "*ccs${op}_cv = *PDL::ccs${op}_cv = ccs_binop_compat_cv(\\&PDL::ccs_${op}_vector_ma);";
}
*ccsadd_cv = *PDL::ccsadd_cv = \&ccsplus_cv;
*ccsdiff_cv = *PDL::ccsdiff_cv = \&ccsminus_cv;
*ccsdiv_cv = *PDL::ccsdiv_cv = \&ccsdivide_cv;

##======================================================================
## Vector Operations: Rows
##======================================================================
=pod

=head2 ccs${OP}_rv

=for sig

  Signature: (int ptr(N); int rowids(Nnz); nzvals_in(Nnz);  rowvec(N);  [o]nzvals_out(Nnz))

Row-vector operation ${OP} on CCS-encoded PDL.
Should do something like the following (without decoding the CCS matrix):

 ($column->slice("*1,") ${OP} ccsdecode($ptr,$rowids,$nzvals))->ccsencode;

Missing values in the CCS-encoded PDL are not effected by this operation.

See ccs${OP}_cv() above for supported operations.

=cut

sub ccs_binop_compat_rv {
  my $ccsop = shift;
  return sub {
    my $ptr = shift;
    my ($ptri,$ptrnzi) = ccs_decode_pointer($ptr->append($_[1]->nelem));
    $ccsop->($ptri, $_[1]->index($ptrnzi), @_[2,3]);
  };
}
foreach my $op (@ccs_binops) {
  eval "*ccs${op}_rv = *PDL::ccs${op}_rv = ccs_binop_compat_rv(\\&PDL::ccs_${op}_vector_ma);";
}
*ccsadd_rv = *PDL::ccsadd_rv = \&ccsplus_rv;
*ccsdiff_rv = *PDL::ccsdiff_rv = \&ccsminus_rv;
*ccsdiv_rv = *PDL::ccsdiv_rv = \&ccsdivide_rv;


##------------------------------------------------------
## Ufuncs (accumulators)

#ccs${name}over()


1; ##-- make perl happy
