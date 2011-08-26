## File: PDL::CCS::Functions.pm
## Author: Bryan Jurish <moocow@cpan.org>
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
   ##
   ##-- Decoding
   qw(ccs_decode ccs_pointerlen),
   ##
   ##-- Vector Operations (compat)
   qw(ccs_binop_vector_mia),
   (map { "ccs_${_}_vector_mia" } (
				   qw(plus minus mult divide modulo power),
				   qw(gt ge lt le eq ne spaceship),
				   qw(and2 or2 xor shiftleft shiftright),
				  )),
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

##======================================================================
## Scalar Operations
=pod

=head1 Scalar Operations

Scalar operations can be performed in parallel directly on C<$nzvals>
(and if applicable on C<$missing> as well):

 $c = 42;

 $nzvals2 = $nzvals  + $c;        $missing2 = $missing  + $c;
 $nzvals2 = $nzvals  - $c;        $missing2 = $missing  - $c;
 $nzvals2 = $nzvals  * $c;        $missing2 = $missing  * $c;
 $nzvals2 = $nzvals  / $c;        $missing2 = $missing  / $c;

 $nzvals2 = $nzvals ** $c;        $missing2 = $missing ** $c;
 $nzvals2 = log($nzvals);         $missing2 = log($missing);
 $nzvals2 = exp($nzvals);         $missing2 = exp(missing);

 $nzvals2 = $nzvals->and2($c,0);  $missing2 = $missing->and($c,0);
 $nzvals2 = $nzvals->or2($c,0);   $missing2 = $missing->or2($c,0);
 $nzvals2 = $nzvals->not();       $missing2 = $missing->not();

Nothing prevents scalar operations from producing new "missing" values (e.g. $nzvals*0),
so you might want to re-encode your compressed data after applying the operation.

=cut


##======================================================================
## Vector Operations
=pod

=head1 Vector Operations

=head2 ccs_OP_vector_mia

=for sig

  Signature: (int whichDimV(Nnz); nzvals(Nnz); vec(V); [o]nzvals_out(Nnz))

A number of row- and column-vector operations may be performed directly
on encoded Nd-PDLs, without the need for decoding to a (potentially huge)
dense temporary.  These operations assume that "missing" values are
annihilators with respect to the operation in question, i.e.
that it holds for all C<$x> in C<$vec> that:

 ($missing __OP__ $x) == $missing

This is in line with the usual PDL semantics if your C<$missing> value is C<BAD>,
but may produce unexpected results when e.g. adding a vector to a sparse PDL with C<$missing>==0.
If you really need to do something like the latter, then you're probably better off decoding to
a dense PDL anyway.

Predefined function names for encoded-PDL vector operations are all of the form:
C<ccs_${OPNAME}_ma>, where ${OPNAME} is the base name of the operation:

 plus       ##-- addition
 minus      ##-- subtraction
 mult       ##-- multiplication (NOT matrix-multiplication)
 divide     ##-- division
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

=head2 \&CODE = ccs_binop_vector_mia($opName, \&PDLCODE);

Returns a generic vector-operation subroutine which reports errors as C<$opName>
and uses \&PDLCODE to perform underlying computation.

=cut

##======================================================================
## Vector Operations: Generic

*PDL::ccs_binop_vector_mia = \&ccs_binop_vector_mia;
sub ccs_binop_vector_mia {
  my ($opName,$pdlCode) = @_;
  return sub {
    my ($wi,$nzvals_in, $vec,$nzvals_out) = @_;
    $nzvals_out = zeroes(($nzvals_in->type > $vec->type ? $nzvals_in->type : $vec->type), $nzvals_in->nelem)
      if (!defined($nzvals_out));
    $pdlCode->($nzvals_in, $vec->index($wi), $nzvals_out, 0);
    return $nzvals_out;
  };
}

##-- Arithmetic
*PDL::ccs_plus_vector_mia   = *ccs_plus_vector_mia   = ccs_binop_vector_mia('plus',\&PDL::plus);      ##-- addition
*PDL::ccs_minus_vector_mia  = *ccs_minus_vector_mia  = ccs_binop_vector_mia('minus',\&PDL::minus);    ##-- subtraction
*PDL::ccs_mult_vector_mia   = *ccs_mult_vector_mia   = ccs_binop_vector_mia('mult',\&PDL::mult);      ##-- multiplication
*PDL::ccs_divide_vector_mia = *ccs_divide_vector_mia = ccs_binop_vector_mia('divide',\&PDL::divide);  ##-- division
*PDL::ccs_modulo_vector_mia = *ccs_modulo_vector_mia = ccs_binop_vector_mia('modulo',\&PDL::modulo);  ##-- modulo
*PDL::ccs_power_vector_mia  = *ccs_power_vector_mia  = ccs_binop_vector_mia('power',\&PDL::power);    ##-- potentiation

##-- Comparison
*PDL::ccs_gt_vector_mia = *ccs_gt_vector_mia = ccs_binop_vector_mia('gt',\&PDL::gt);        ##-- greater-than
*PDL::ccs_ge_vector_mia = *ccs_ge_vector_mia = ccs_binop_vector_mia('ge',\&PDL::ge);        ##-- greater-than-or-equal
*PDL::ccs_lt_vector_mia = *ccs_lt_vector_mia = ccs_binop_vector_mia('lt',\&PDL::lt);        ##-- less-than
*PDL::ccs_le_vector_mia = *ccs_le_vector_mia = ccs_binop_vector_mia('le',\&PDL::le);        ##-- less-than-or-equal
*PDL::ccs_eq_vector_mia = *ccs_eq_vector_mia = ccs_binop_vector_mia('eq',\&PDL::eq);        ##-- equality
*PDL::ccs_ne_vector_mia = *ccs_ne_vector_mia = ccs_binop_vector_mia('ne',\&PDL::ne);        ##-- inequality
*PDL::ccs_spaceship_vector_mia = *ccs_spaceship_vector_mia = ccs_binop_vector_mia('spaceship',\&PDL::spaceship); ##-- <=>

##-- Logic & Bitwise
*PDL::ccs_and2_vector_mia = *ccs_and2_vector_mia = ccs_binop_vector_mia('and',\&PDL::and2);     ##-- logical AND (and2)
*PDL::ccs_or2_vector_mia  = *ccs_or2_vector_mia  = ccs_binop_vector_mia('or',\&PDL::or2);       ##-- logical OR   (or2)
*PDL::ccs_xor_vector_mia  = *ccs_xor_vector_mia  = ccs_binop_vector_mia('xor',\&PDL::xor);      ##-- binary XOR  (xor)
*PDL::ccs_shiftleft_vector_mia = *ccs_shiftleft_vector_mia = ccs_binop_vector_mia('shiftleft',\&PDL::shiftleft);     ##-- <<
*PDL::ccs_shiftright_vector_mia = *ccs_shiftright_vector_mia = ccs_binop_vector_mia('shiftright',\&PDL::shiftright); ##-- >>

1; ##-- make perl happy

