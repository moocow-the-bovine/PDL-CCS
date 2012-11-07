## File: PDL::CCS::Functions.pm
## Author: Bryan Jurish <moocow@cpan.org>
## Description: useful perl-level functions for PDL::CCS

package PDL::CCS::Functions;
use PDL::CCS::Version;
use PDL::CCS::Utils;
use PDL::VectorValued;
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
   ##
   ##-- qsort
   qw(ccs_qsort),
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

  Signature: (int ptr(N+1); int [o]len(N))

Get number of non-missing values for each axis value from a CCS-encoded
offset pointer vector $ptr().

=cut

;#-- emacs

*PDL::ccs_pointerlen = \&ccs_pointerlen;
sub ccs_pointerlen :lvalue {
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
sub ccs_decode  :lvalue {
  my ($aw,$nzvals,$missing,$dims,$a) = @_;
  $missing = $PDL::undefval if (!defined($missing));
  if (!defined($dims)) {
    barf("PDL::CCS::ccs_decode(): whichnd() is empty; you must specify \@Dims!") if ($aw->isempty);
    $dims = [ map {$aw->slice("($_),")->max+1} (0..($aw->dim(0)-1))];
  }
  $a    = zeroes($nzvals->type, @$dims) if (!defined($a));
  $a   .= $missing;

  (my $tmp=$a->indexND($aw)) .= $nzvals; ##-- CPAN tests puke here with "Can't modify non-lvalue subroutine call" in 5.15.x (perl bug #107366)

  ##-- workaround for missing empty pdl support in PDL 2.4.10 release candidates (pdl bug #3462924), fixed in 2.4.9_993
  #$a->indexND($aw) .= $nzvals if (!$nzvals->isempty);
  #if (!$nzvals->isempty) {
  #  my $tmp = $a->indexND($aw);
  #  $tmp .= $nzvals;
  #}

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
  return sub :lvalue {
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

##======================================================================
## Sorting
=pod

=head1 Sorting

=head2 ccs_qsort

=for sig

  Signature: (int whichIn(Ndims,Nnz); nzValsIn(Nnz); missing(); Dim0(); int [o]nzIxOut; int [o]whichOut; int [oca]nzValsOut)

Underlying guts for PDL::CCS::Nd::qsort() and PDL::CCS::Nd::qsorti().
Given an index pdl C<$whichIn>, a corresponding value-pdl C<$nzValsIn>, and a missing value $missing(),
determines indices C<$whichOut> and non-missing values C<$nzValsOut> suitable
for constructing compressed representation of the input data sorted along the 0th dimension, together
with their primary Nnz-indices C<$nzIxOut>.
Note that C<$nzValsOut> will be ALWAYS be returned as an index-reordering
of C<$nzValsIn> with dataflow open; copy or sever() it yourself if required.

qsort() on compressed piddles can then be implemented as:

 $whichnd = $whichOut;
 $nzvals  = $nzValsOut;

and qsorti() as:

 $whichnd = $whichOut;
 $nzvals  = $whichIn->slice("(0),")->index($nzIxOut);

=cut

## $bool = _checkdims(\@dims1,\@dims2,$label);  ##-- match      @dims1 ~ @dims2
## $bool = _checkdims( $pdl1,   $pdl2,$label);  ##-- match $pdl1->dims ~ $pdl2->dims
sub _checkdims {
  #my ($dims1,$dims2,$label) = @_;
  #my ($pdl1,$pdl2,$label) = @_;
  my $d0 = UNIVERSAL::isa($_[0],'PDL') ? pdl(long,$_[0]->dims) : pdl(long,$_[0]);
  my $d1 = UNIVERSAL::isa($_[1],'PDL') ? pdl(long,$_[1]->dims) : pdl(long,$_[0]);
  barf(__PACKAGE__ . "::_checkdims(): dimension mismatch for ".($_[2]||'pdl').": $d0!=$d1")
    if (($d0->nelem!=$d1->nelem) || !all($d0==$d1));
  return 1;
}

sub ccs_qsort {
  my ($whichIn,$nzValsIn,$missing,$dim0, $nzIxOut,$whichOut) = @_;

  ##-- prepare
  $whichOut = $whichIn->zeroes if (!defined($whichOut));
  $whichOut->reshape($whichIn->dims()) if (!$whichOut->isempty);
  _checkdims($whichIn,$whichOut,'ccs_qsort: whichIn~whichOut');
  ##
  $nzIxOut = zeroes(long,$whichIn->dim(1)) if (!defined($nzIxOut));
  $nzIxOut->reshape($whichIn->dim(1)) if ($nzIxOut->isempty);
  _checkdims($nzValsIn,$nzIxOut,'ccs_qsort: nzValsIn~nzIxOut');
  ##
  $dim0 = $whichIn->slice("(0),")->max+1 if (!defined($dim0));

  ##-- collect and sort base data (unsorted indices + values)
  my $wv  = $whichIn->glue(0,$nzValsIn->slice("*1,"));
  $wv->slice("1:-1,")->vv_qsortveci($nzIxOut);

  ##-- enumerate output values, splitting enum-values around $missing()
  $whichOut     .= $whichIn->dice_axis(1,$nzIxOut);
  my $whichOut0  = $whichOut->slice("(0),");
  my $whichOut1  = $whichOut->dim(0)>1 ? $whichOut->slice("1:-1,") : $whichOut0->zeroes; ##-- handle flat sorts without key indices
  my $nzValsOut  = $nzValsIn->index($nzIxOut);
  my ($nzii_l,$nzii_r) = (defined($missing) ? which_both($nzValsOut <= $missing) : ($nzValsOut->xvals,null));
  #
  #$whichOut0 .= -1; ##-- debug
  $whichOut0->index($nzii_l) .= $whichOut1->dice_axis(1,$nzii_l)->enumvec if (!$nzii_l->isempty);
  $whichOut0->index($nzii_r)->slice("-1:0") .= ($dim0 - 1) - $whichOut1->dice_axis(1,$nzii_r)->slice(",-1:0")->enumvec if (!$nzii_r->isempty);

  ##-- all done
  return wantarray ? ($nzIxOut,$whichOut,$nzValsOut) : $nzIxOut;
}

##======================================================================
## Vector Operations: Generic


##======================================================================
## POD: footer
=pod

=head1 ACKNOWLEDGEMENTS

Perl by Larry Wall.

PDL by Karl Glazebrook, Tuomas J. Lukka, Christian Soeller, and others.

=cut


##---------------------------------------------------------------------
=pod

=head1 AUTHOR

Bryan Jurish E<lt>moocow@cpan.orgE<gt>

=head2 Copyright Policy

Copyright (C) 2007-2012, Bryan Jurish. All rights reserved.

This package is free software, and entirely without warranty.
You may redistribute it and/or modify it under the same terms
as Perl itself.

=head1 SEE ALSO

perl(1),
PDL(3perl),
PDL::CCS::Nd(3perl),


=cut


1; ##-- make perl happy
