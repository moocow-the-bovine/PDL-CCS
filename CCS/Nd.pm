## File: PDL::CCS::Nd.pm
## Author: Bryan Jurish <moocow@ling.uni-potsdam.de>
## Description: N-dimensional CCS-encoded pseudo-PDL

package PDL::CCS::Nd;
use PDL::Lite;
use PDL::VectorValued;
use PDL::CCS::Version;
use PDL::CCS::Functions qw(ccs_decode ccs_pointerlen);
use PDL::CCS::Utils     qw(ccs_encode_pointers ccs_decode_pointer);
use PDL::CCS::Ufunc;
use PDL::CCS::Ops;
use UNIVERSAL 'isa';
use Carp;
use strict;

our $VERSION = $PDL::CCS::VERSION;
our @ISA = ('PDL');
our @EXPORT_OK =
  (
   ##-- Encoding/Decoding
   qw(toccs todense),
  );
our %EXPORT_TAGS =
  (
   Func => [@EXPORT_OK],               ##-- respect PDL conventions (hopefully)
  );

##--------------------------------------------------------------
## Operator overloading
use overload (
	      ##-- Binary ops: arithmetic
	      "+" => \&plus_mia,
	      "-" => \&minus_mia,
	      "*" => \&mult_mia,
	      "/" => \&divide_mia,
	      "%" => \&modulo_mia,
	      "**" => \&power_mia,

	      ##-- Binary ops: comparisons
	      ">"  => \&gt_mia,
	      "<"  => \&lt_mia,
	      ">=" => \&ge_mia,
	      "<=" => \&le_mia,
	      "<=>" => \&spaceship_mia,
	      "==" => \&eq_mia,
	      "!=" => \&ne_mia,
	      #"eq" => \&eq_mia

	      ##-- Binary ops: bitwise & logic
	      "|"  => \&or2_mia,
	      "&"  => \&and2_mia,
	      "!"  => \&not_mia,
	      "~"  => \&bitnot_mia,
	      "^"  => \&xor_mia,
	      "<<" => \&shiftleft_mia,
	      ">>" => \&shiftright_mia,

	      ##-- Unary operations
	      "sqrt" => \&sqrt_mia,
	      "abs"  => \&abs_mia,
	      "sin"  => \&sin_mia,
	      "cos"  => \&cos_mia,
	      "log"  => \&log_mia,
	      "exp"  => \&exp_mia,

	      ##-- TODO: assignment (.=), assigning-variants (+=, *=, etc.), matmult (x)

	      ##-- Stringification
	      "\"\"" => \&string,
	     );

##======================================================================
## Globals

our $DIMS    = 0;
our $XDIMS   = 1; ##-- dimension translation pdl: $whichND_logical = $WHICH->dice_axis(0,$XDIMS)
our $WHICH   = 2;
our $VALS    = 3;
our $PTRS    = 4;

our $P_BYTE = PDL::byte();
our $P_LONG = PDL::long();

##======================================================================
## Constructors etc.

## $obj = $class_or_obj->newFromDense($denseND);
## $obj = $class_or_obj->newFromDense($denseND,$missing);
##  + object structure: ARRAY
##     $DIMS    => \@dims,     ##-- $dimi => $dim,
##     $XDIMS   => $xdims,     ##-- pdl(long,$Ndims) : $LOGICAL_DIM => $PHYSICAL_DIM
##                             ##   + s.t. $whichND_logical = $whichND->dice_axis(0,$xdims)
##     $WHICH   => $whichND,   ##-- pdl(long,$Ndims,$Nnz) ~ $dense_orig->whichND
##                             ##   #+ NOT guaranteed to be sorted in any meaningful way
##                             ##   + guaranteed to be sorted as for qsortvec() specs
##                             ##   + NOT changed by dimension-shuffling transformations
##     $VALS    => $vals,      ##-- pdl( ?  ,$Nnz+1)      ~ $dense->where($dense)->append($missing)
##     $PTRS    => \@PTRS,     ##-- array of ccsutils-pointers by original dimension number
##
##  + each element of @PTRS is itself an array:
##     $PTRS[$i] => [ $PTR, $NZI ]
##
sub newFromDense {
  my ($that,$p,$missing) = @_;
  $missing     = (defined($missing)
		  ? PDL->pdl($p->type,$missing)
		  : ($p->badflag
		     ? PDL->pdl($p->type,0)->setvaltobad(0)
		     : PDL->pdl($p->type,0)));
  my $pwhichND = ($missing->isbad->any
		  ? $p->isgood()
		  : ($p != $missing)
		 )->whichND->vv_qsortvec;
  my $pnz  = $p->indexND($pwhichND)->append($missing);
  $pnz->sever;                       ##-- always sever nzvals ?
  my $self = bless [], ref($that)||$that;
  $self->[$DIMS]    = [$p->dims];
  $self->[$XDIMS]   = PDL->sequence($P_LONG,scalar(@{$self->[$DIMS]}));
  $self->[$WHICH]   = $pwhichND;
  $self->[$VALS]    = $pnz;
  $self->[$PTRS]    = [];            ##-- do we really need this ... yes
  return $self;
}

## $obj = $class_or_obj->newFromWhich($whichND,$nzvals,%options);
## $obj = $class_or_obj->newFromWhich($whichND,$nzvals);
##  + %options:
##     dims    => \@dims,   ##-- dimension list; default guessed from $whichND
##     missing => $missing, ##-- default: BAD if $nzvals->badflag, 0 otherwise
##     xdims   => $xdims,   ##-- dimension translation PDL (default: sequence($ndims)) : $LOGICAL_DIM => $PHYSICAL_DIM
##  + no dataflow!
sub newFromWhich {
  my ($that,$wnd,$nzvals,%opts) = @_;
  my $missing = (defined($opts{missing})
		 ? PDL->pdl($nzvals->type,$opts{missing})
		 : ($nzvals->badflag
		    ? PDL->pdl($nzvals->type,0)->setvaltobad(0)
		    : PDL->pdl($nzvals->type,0)));
  my $dims  = ( defined($opts{dims})  ? $opts{dims}  : [($wnd->xchg(0,1)->maximum+1)->list]  );
  my $xdims = ( defined($opts{xdims}) ? $opts{xdims} : PDL->sequence($P_LONG,scalar(@$dims)) );
  $wnd     = $wnd->dice_axis(0,$xdims);
  my $wi   = $wnd->qsortveci;
  $wnd     = $wnd->dice_axis(1,$wi)->pdl();         ##-- copy!
  $nzvals  = $nzvals->index($wi)->append($missing); ##-- copy (b/c append)
  my $self = bless [], ref($that)||$that;
  $self->[$DIMS]    = $dims;
  $self->[$XDIMS]   = PDL->sequence($P_LONG,scalar(@$dims));
  $self->[$WHICH]   = $wnd;
  $self->[$VALS]    = $nzvals;
  $self->[$PTRS]    = [];
  return $self;
}

## DESTROY : avoid PDL inheritanc
sub DESTROY { ; }


## $ccs = $pdl->toccs()
## $ccs = $pdl->toccs($missing)
*PDL::toccs = \&toccs;
sub toccs {
  return $_[0] if (isa($_[0],__PACKAGE__));
  return __PACKAGE__->newFromDense(@_)
}

## $ccs = $ccs->copy()
BEGIN { *clone = \&copy; }
sub copy {
  my $ccs1 = shift;
  my $ccs2 = bless [], ref($ccs1);
  $ccs2->[$DIMS]  = [ @{$_[0]->[$DIMS]} ];
  $ccs2->[$XDIMS] = $ccs1->[$XDIMS]->pdl;
  $ccs2->[$WHICH] = $ccs1->[$WHICH]->pdl;
  $ccs2->[$VALS]  = $ccs1->[$VALS]->pdl;
  $ccs2->[$PTRS]  = [ map {defined($_) ? [map {$_->pdl} @$_] : undef} @{$ccs1->[$PTRS]} ]; ##-- copy pointers?
  #$ccs2->[$MISSING] = $ccs1->[$MISSING]->pdl();
  return $ccs2;
}

## $ccs2 = $ccs->copyShallow()
##  + a very shallow version of copy()
##  + Copied    : @$DIMS, @$PTRS, @{$PTRS->[*]}
##  + Referenced: $XDIMS, $WHICH, $VALS,  $PTRS->[*][*]
sub copyShallow {
  my $ccs = bless [@{$_[0]}], ref($_[0]);
  ##
  ##-- do copy some of it
  $ccs->[$DIMS]  = [ @{$ccs->[$DIMS]} ];
  #$ccs->[$XDIMS] = $ccs->[$XDIMS]->pdl;
  $ccs->[$PTRS]  = [ map {defined($_) ? [@$_] : undef} @{$ccs->[$PTRS]} ];
  $ccs;
}

## $ccs2 = $ccs->shadow(%args)
##  + args:
##     dims  => \@dims,   ##-- default: [@$dims1]
##     xdims => $xdims2,  ##-- default: $xdims1->pdl
##     ptrs  => \@ptrs2,  ##-- default: []
##     which => $which2,  ##-- default: undef
##     vals  => $vals2,   ##-- default: undef
sub shadow {
  my ($ccs,%args) = @_;
  my $ccs2        = bless [], ref($ccs)||$ccs;
  $ccs2->[$DIMS]  = defined($args{dims})  ? $args{dims}  : [@{$ccs->[$DIMS]}];
  $ccs2->[$XDIMS] = defined($args{xdims}) ? $args{xdims} : $ccs->[$XDIMS]->pdl;
  $ccs2->[$PTRS]  = $args{ptrs}  ? $args{ptrs} : [];
  $ccs2->[$WHICH] = $args{which};
  $ccs2->[$VALS]  = $args{vals};
  return $ccs2;
}

##--------------------------------------------------------------
## Maintainence

## $ccs = $ccs->sortwhich()
##  + sorts on $ccs->[$WHICH]
##  + may be DANGEROUS to indexing methods, b/c it alters $VALS
##  + clears pointers
sub sortwhich {
  my $sorti     = $_[0][$WHICH]->qsortveci;
  $_[0][$WHICH] = $_[0][$WHICH]->dice_axis(1,$sorti);
  $_[0][$VALS]  = $_[0][$VALS]->index($sorti->append($_[0][$WHICH]->dim(1)));
#  foreach (grep {defined($_)} @{$_[0][$PTRS]}) {
#    $_->[1]->index($sorti) .= $_->[1]; ##-- dangerous
#  }
  @{$_[0][$PTRS]} = qw() if (! ($sorti==PDL->sequence($P_LONG,$sorti->dims))->all );
  return $_[0];
}


##--------------------------------------------------------------
## Decoding

## $dense = $ccs->decode()
## $dense = $ccs->decode($dense)
sub decode {
  ccs_decode($_[0][$WHICH]->dice_axis(0,$_[0][$XDIMS]),
	     $_[0][$VALS]->slice("0:-2"),
	     $_[0][$VALS]->slice("(-1)"),
	     [@{$_[0][$DIMS]}[$_[0][$XDIMS]->list]],
	     @_[1..$#_]);
}

## $dense = $ccs_or_dense->todense()
*PDL::todense = \&todense;
sub todense { isa($_[0],__PACKAGE__) ? $_[0]->decode() : $_[0]; }


##--------------------------------------------------------------
## Basic Properties

## $type = $obj->type()
sub type { $_[0][$VALS]->type; }

## @dims = $obj->dims()
sub dims { @{$_[0][$DIMS]}[$_[0][$XDIMS]->list]; }

## $dim = $obj->dim($dimi)
*getdim = \&dim;
sub dim { $_[0][$DIMS][$_[0][$XDIMS]->at($_[1])]; }

## $ndims = $obj->ndims()
sub ndims { scalar(@{$_[0][$DIMS]}); }

## $nelem = $obj->nelem
sub nelem { PDL->pdl($P_LONG,$_[0][$DIMS])->dprod; }

## $ngood = $obj->ngood
sub ngood {
  ($_[0][$VALS]->slice("0:-2")->ngood(@_[1..$#_])
   + ($_[0][$VALS]->slice("-1")->isbad()
      ? 0
      : ($_[0]->nelem - $_[0][$WHICH]->dim(1))));
}

## $nnz = $obj->nnz
##   + FIXME: respect dimensions / Ufunc behavor!
sub nnz {
  ($_[0][$VALS]->slice("0:-2")->nnz(@_[1..$#_])
   + ($_[0][$VALS]->at(-1)==0
      ? 0
      : ($_[0]->nelem - $_[0][$WHICH]->dim(1))));
}

##--------------------------------------------------------------
## Low-level CCS access

## $xdims = $obj->xdims()
sub xdims { $_[0][$XDIMS]; }

## $nnz = $obj->nstored()
sub nstored { $_[0][$WHICH]->dim(1); }

## $nmissing = $obj->nmissing()
sub nmissing { 1 + $_[0]->nelem - $_[0]->[$VALS]->nelem; }

## $missing = $obj->missing()
## $missing = $obj->missing($missing)
sub missing {
  $_[0][$VALS]->set(-1,$_[1]) if (@_>1);
  $_[0][$VALS]->slice("-1");
}

## $whichND_stored = $obj->_whichND()
## $whichND_stored = $obj->_whichND($whichND)
sub _whichND {
  $_[0][$WHICH] = $_[1] if (@_>1);
  $_[0][$WHICH];
}

## $nzvals = $obj->nzvals();
## $nzvals = $obj->nzvals($nzvals)
sub nzvals {
  $_[0][$VALS]=$_[1]->append($_[0][$VALS]->slice("-1")) if (@_ > 1);
  $_[0][$VALS]->slice("0:-2");
}

## $vals = $obj->vals()
## $vals = $obj->vals($storedvals)
sub vals {
  $_[0][$VALS]=$_[1] if (@_ > 1);
  $_[0][$VALS];
}


## $ptr           = $obj->ptr($dim); ##-- scalar context
## ($ptr,$pi2nzi) = $obj->ptr($dim); ##-- list context
##   + returns cached value in $ccs->[$PTRS][$dim] if present
##   + caches value in $ccs->[$PTRS][$dim] otherwise
##   + $dim defaults to zero, for compatibility
##   + if $dim is zero, all($pi2nzi==sequence($obj->nstored))
sub ptr {
  my ($ccs,$dim) = @_;
  $dim = 0 if (!defined($dim));
  $ccs->[$PTRS][$dim] = [$ccs->getptr($dim)] if (!defined($ccs->[$PTRS][$dim]));
  return wantarray ? @{$ccs->[$PTRS][$dim]} : $ccs->[$PTRS][$dim][0];
}

## ($ptr,$pi2nzi) = $obj->getptr($dim);
##  + as for ptr(), but does NOT cache anything, and does NOT check the cache
sub getptr { ccs_encode_pointers($_[0][$WHICH]->slice("($_[1]),"), $_[0][$DIMS][$_[1]]); }

##--------------------------------------------------------------
## CONTINUE HERE

## TODO:
##  + type conversion
##  + inplace & dataflow stuff : inplace, sever, ...
##  + index translation (flati|ndi) <-> (nzi)
##  + OPERATIONS:
##    - accumulators: sumover, prodover, ...
##    - unary: bnot, ++, --
##    - binary: logical, arithmetic
##      * with rowVec __OR__ columnVec __OR__ PDL::CCS::Nd
##    - matrix: matmult

##--------------------------------------------------------------
## Dimension Shuffling

## $ccs2 = $ccs->reorder_pdl($dimpdl)
sub reorder_pdl {
  my $ccs2 = $_[0]->copyShallow;
  $ccs2->[$XDIMS] = $ccs2->[$XDIMS]->index($_[1]);
  $ccs2->[$XDIMS]->sever;
  $ccs2;
}

## $ccs2 = $ccs->reorder(@dimlist)
sub reorder { $_[0]->reorder_pdl(PDL->pdl($P_LONG,@_[1..$#_])); }

## $ccs2 = $ccs->xchg($dim1,$dim2)
sub xchg {
  my $dimpdl = PDL->sequence($P_LONG,scalar(@{$_[0][$DIMS]}));
  my $tmp    = $dimpdl->at($_[1]);
  $dimpdl->set($_[1], $dimpdl->at($_[2]));
  $dimpdl->set($_[2], $tmp);
  return $_[0]->reorder_pdl($dimpdl);
}

## $ccs2 = $ccs->mv($dimFrom,$dimTo)
sub mv  {
  my ($d1,$d2) = @_[1,2];
  $d1 = @{$_[0][$DIMS]}+$d1 if ($d1 < 0);
  $d2 = @{$_[0][$DIMS]}+$d2 if ($d2 < 0);
  return $_[0]->reorder($d1 < $d2
			? ((0..($d1-1)), (($d1+1)..$d2), $d1, (($d2+1)..$#{$_[0][$DIMS]}))
			: ((0..($d2-1)), $d1, ($d2..($d1-1)), (($d1+1)..$#{$_[0][$DIMS]}))
		       );
}

## $ccs2 = $ccs->transpose()
sub transpose { return $_[0]->xchg(0,1)->copy; }



##--------------------------------------------------------------
## Indexing

## $nzi = $ccs->indexNDi($ndi)
sub indexNDi {
  my $foundi       = $_[1]->vsearchvec($_[0][$WHICH]);
  my $foundi_mask  = ($_[1]==$_[0][$WHICH]->dice_axis(1,$foundi))->andover;
  $foundi_mask->inplace->not;
  $foundi->where($foundi_mask) .= $_[0][$WHICH]->dim(1);
  return $foundi;
}

## $vals = $ccs->indexND($ndi)
sub indexND { $_[0][$VALS]->index($_[0]->indexNDi($_[1])); }

## $vals = $ccs->index2d($xi,$yi)
sub index2d { $_[0]->indexND($_[1]->cat($_[2])->xchg(0,1)); }

## $vals = $ccs->index($flati)
sub index {
  my ($ccs,$i) = @_;
  my $dummy  = PDL->pdl(0)->slice(join(',', map {"*$_"} @{$ccs->[$DIMS]}));
  my @coords = $dummy->one2nd($i);
  my $ind = PDL->zeroes($P_LONG,scalar(@{$ccs->[$DIMS]}),$i->dims);
  $ind->slice("($_),") .= $coords[$_] foreach (0..$#coords);
  return $ccs->indexND($ind);
}

## $ccs2 = $ccs->dice_axis($axis, $axisi)
##  + returns a new ccs object, should participate in dataflow
sub dice_axis {
  my ($ccs,$axis,$axisi) = @_;
  my ($ptr,$pi2nzi)    = $ccs->ptr($axis);
  my ($ptrix,$pi2nzix) = $ptr->ccs_decode_pointer($axisi);
  my $nzix   = $pi2nzi->index($pi2nzix);
  my $which  = $ccs->[$WHICH]->dice_axis(1,$nzix);
  $which->sever;
  $which->slice("($axis),") .= $ptrix;
  my $nzvals = $ccs->[$VALS]->index($nzix->append($ccs->[$WHICH]->dim(1)));
  ##
  ##-- construct output object
  my $ccs2              = $ccs->shadow();
  $ccs2->[$DIMS][$axis] = $axisi->nelem;
  $ccs2->[$WHICH]       = $which;
  $ccs2->[$VALS]        = $nzvals;
  ##
  ##-- sort output object (if not dicing on 0th dimension)
  return $axis==0 ? $ccs2 : $ccs2->sortwhich();
}

## $onedi = $ccs->n2oned($ndi)
##  + returns a pseudo-index
sub n2oned {
  my $dimsizes = PDL->pdl($P_LONG,[1,@{$_[0][$DIMS]}])->slice("0:-2")->cumuprodover;
  return ($_[1] * $dimsizes)->sumover;
}

## $whichND = $obj->whichND
##  + just returns the literal index PDL: beware of dataflow!
sub whichND { $_[0][$WHICH]->dice_axis(0,$_[0][$XDIMS]); }

## $which = $obj->which()
sub which { $_[0]->n2oned($_[0]->whichND)->qsort; }

## $val = $ccs->at(@index)
sub at { $_[0]->indexND(PDL->pdl($P_LONG,@_[1..$#_]))->sclr; }

## $val = $ccs->set(@index,$value)
sub set {
  my $foundi = $_[0]->indexNDi(PDL->pdl($P_LONG,@_[1..($#_-1)]));
  if ( ($foundi==$_[0][$WHICH]->dim(1))->any ) {
    carp(ref($_[0]).": cannot set() a missing value!")
  } else {
    $_[0][$VALS]->index($foundi) .= $_[$#_];
  }
  return $_[0];
}

##--------------------------------------------------------------
## Ufuncs

## $ufunc_sub = _ufuncsub($subname, \&ccs_accum_sub)
sub _ufuncsub {
  my ($subname,$accumsub) = @_;
  return sub {
    my $ccs = shift;
    ##
    ##-- preparation
    my $which   = $ccs->whichND;
    my $which0  = $which->slice("(0),");
    my $which1  = $which->slice("1:-1,");
    my $sorti   = $which1->qsortveci;
    my $vals    = $ccs->[$VALS];
    my $vals1   = $vals->index($sorti);
    my $missing = $vals->slice("(-1)");
    my @dims    = $ccs->dims;
    ##
    ##-- guts
    my ($which2,$nzvals2) = $accumsub->($which1->dice_axis(1,$sorti),$vals1,
					($missing->isgood ? ($missing,$dims[0]) : (0,0))
				       );
    ##
    ##-- get output pd
    shift(@dims);
    return $ccs->shadow(
			dims  =>[@dims],
			xdims =>PDL->sequence($P_LONG,$#{$ccs->[$DIMS]}),
			which =>$which2,
			vals  =>$nzvals2->append($missing->convert($nzvals2->type)),
		       );
  };
}

foreach my $ufunc (
		   qw(prod dprod sum dsum),
		   qw(and or band bor),
		  )
  {
    eval "*${ufunc}over = _ufuncsub('${ufunc}over', PDL::CCS::Ufunc->can('ccs_accum_${ufunc}'))";
  }
foreach my $ufunc (qw(maximum minimum))
  {
    eval "*${ufunc} = _ufuncsub('${ufunc}', PDL::CCS::Ufunc->can('ccs_accum_${ufunc}'))";
  }

sub sum   { my $z=$_[0]->missing; $_[0][$VALS]->slice("0:-2")->sum   + ($z->isgood ? ($z->sclr *  $_[0]->nmissing) : 0); }
sub dsum  { my $z=$_[0]->missing; $_[0][$VALS]->slice("0:-2")->dsum  + ($z->isgood ? ($z->sclr *  $_[0]->nmissing) : 0); }
sub prod  { my $z=$_[0]->missing; $_[0][$VALS]->slice("0:-2")->prod  * ($z->isgood ? ($z->sclr ** $_[0]->nmissing) : 1); }
sub dprod { my $z=$_[0]->missing; $_[0][$VALS]->slice("0:-2")->dprod * ($z->isgood ? ($z->sclr ** $_[0]->nmissing) : 1); }
sub min   { $_[0][$VALS]->min; }
sub max   { $_[0][$VALS]->max; }

sub any { $_[0][$VALS]->any; }
sub all { $_[0][$VALS]->all; }


##--------------------------------------------------------------
## Unary Operations

## $sub = _unary_op($opname,$pdlsub)
sub _unary_op {
  my ($opname,$pdlsub) = @_;
  return sub { $_[0]->shadow(which=>$_[0][$WHICH]->pdl, vals=>$pdlsub->($_[0][$VALS])); };
}

foreach my $unop (qw(bitnot sqrt abs sin cos not exp log log10))
  {
    eval "*${unop} = _unary_op('${unop}',PDL->can('${unop}'));";
  }

##--------------------------------------------------------------
## Binary Operations: missing-is-annihilator

## \@matchdims_or_undef = _ccs_binop_relevant_dims(\@dims1,\@dims2)
##  + barf()s on mismath
sub _ccs_binop_relevant_dims {
  my ($opname,$dims1,$dims2) = @_;
  my $rdims = [];
  foreach (0..($#$dims1 < $#$dims2 ? $#$dims1 : $#$dims2)) {
    next if ($dims1->[$_] <= 1 || $dims2->[$_] <= 1);
    if ($dims1->[$_] != $dims2->[$_]) {
      barf("PDL::CCS::Nd::",
	   ($opname||'_ccs_binop_relevant_dims'),
	   "(): mismatch on non-trivial dim($_): $dims1->[$_] != $dims2->[$_]");
    }
    push(@$rdims,$_);
  }
  return $rdims;
}

#our $BLOCKSIZE_MIN     =    64;
#our $BLOCKSIZE_MAX     = undef;
##--
our $BLOCKSIZE_MIN     =      1;
our $BLOCKSIZE_MAX     =  65535;
sub _binary_op_mia {
  my ($opname,$pdlsub,$deftype) = @_;

  return sub {
    my ($a,$b,$swap) = @_;
    $swap=0 if (!defined($swap));
    $b = $b->toccs;

    ##-- get relevant dimensions
    my @adims = $a->dims;
    my @bdims = $b->dims;
    my @cdims = (
		 map {
		   $_ <= $#adims && $adims[$_] > 1 ? $adims[$_] : $bdims[$_]
		 } (0..($#adims > $#bdims ? $#adims : $#bdims))
		);
    my $rdims = _ccs_binop_relevant_dims($opname,\@adims,\@bdims);
    my $rdpdl = PDL->pdl($P_LONG,$rdims);

    ##-- get & sort relevant indices
    my $ixa        = $a->[$WHICH]->dice_axis(0,$a->[$XDIMS]);
    my $nixa       = $ixa->dim(1);
    my $ixar       = $ixa->dice_axis(0,$rdpdl);
    my $ixar_sorti = $ixar->qsortveci;
    $ixa           = $ixa->dice_axis(1,$ixar_sorti);
    $ixar          = $ixar->dice_axis(1,$ixar_sorti);
    ##
    my $ixb        = $b->[$WHICH]->dice_axis(0,$b->[$XDIMS]);
    my $nixb       = $ixb->dim(1);
    my $ixbr       = $ixb->dice_axis(0,$rdpdl);
    my $ixbr_sorti = $ixbr->qsortveci;
    $ixb           = $ixb->dice_axis(1,$ixbr_sorti);
    $ixbr          = $ixbr->dice_axis(1,$ixbr_sorti);

    ##-- initialize: values
    my $avals  = $a->[$VALS];
    my $avalsr = $avals->index($ixar_sorti);
    my $bvals  = $b->[$VALS];
    my $bvalsr = $bvals->index($ixbr_sorti);

    ##-- initialize: state vars
    my $blksz  = $nixa > $nixb ? $nixa : $nixb;
    $blksz     = $BLOCKSIZE_MIN if (defined($BLOCKSIZE_MIN) && $blksz < $BLOCKSIZE_MIN);
    $blksz     = $BLOCKSIZE_MAX if (defined($BLOCKSIZE_MAX) && $blksz > $BLOCKSIZE_MAX);
    my $istate = PDL->zeroes($P_LONG,7); ##-- [ nnzai,nnzai_nxt, nnzbi,nnzbi_nxt, nnzci,nnzci_nxt, cmpval ]
    my $ostate = $istate->pdl;

    ##-- initialize: output vectors
    my $nzai   = PDL->zeroes($P_LONG,$blksz);
    my $nzbi   = PDL->zeroes($P_LONG,$blksz);
    my $nzc    = PDL->zeroes((defined($deftype)
			      ? $deftype
			      : ($avals->type > $bvals->type
				 ? $avals->type
				 : $bvals->type)),
			     $blksz);
    my $ixc    = PDL->zeroes($P_LONG, scalar(@cdims), $blksz);
    my $nnzc   = 0;
    my $zc     = $pdlsub->($avals->slice("-1"), $bvals->slice("-1"), $swap)->convert($nzc->type);
    my $zc_isbad = $zc->isbad ? 1 : 0;

    ##-- block-wise variables
    my ($nzai_prv,$nzai_pnx, $nzbi_prv,$nzbi_pnx, $nzci_prv,$nzci_pnx,$cmpval_prv);
    my ($nzai_cur,$nzai_nxt, $nzbi_cur,$nzbi_nxt, $nzci_cur,$nzci_nxt,$cmpval);
    my ($nzci_max, $blk_slice, $nnzc_blk,$nnzc_slice_blk);
    my ($nzai_blk,$nzbi_blk,$ixa_blk,$ixb_blk,$nzc_blk,$cimask_blk,$ciwhich_blk);
    my $nnzc_prev=0;
    do {
      ##-- align a block
      ccs_binop_align_block_mia($ixar,$ixbr,$istate, $nzai,$nzbi,$ostate);

      ##-- parse current I/O state
      ($nzai_prv,$nzai_pnx, $nzbi_prv,$nzbi_pnx, $nzci_prv,$nzci_pnx,$cmpval_prv) = $istate->list;
      ($nzai_cur,$nzai_nxt, $nzbi_cur,$nzbi_nxt, $nzci_cur,$nzci_nxt,$cmpval)     = $ostate->list;
      $nzci_max = $nzci_cur-1;

      if ($nzai_cur >= $nixa && $nzbi_cur >= $nixb) {
	##-- final block: trim pdls
	$nzai = $nzai->slice("0:$nzci_max");
	$nzbi = $nzbi->slice("0:$nzci_max");

	$nzc  = $nzc->slice("0:".($nnzc+$nzci_max));
	$ixc  = $ixc->slice(",0:".($nnzc+$nzci_max));

	last if ($nzci_max < 0);
      }

      ##-- construct block output pdls: nzvals
      $blk_slice = "${nzci_prv}:${nzci_max}";
      $nzai_blk  = $nzai->slice($blk_slice);
      $nzbi_blk  = $nzbi->slice($blk_slice);
      $nzc_blk   = $pdlsub->($avalsr->index($nzai_blk), $bvalsr->index($nzbi_blk), $swap);

      ##-- get indices of "good" c() values
      $cimask_blk   = $zc_isbad ? $nzc_blk->isgood : ($nzc_blk!=$zc);
      $ciwhich_blk  = $cimask_blk->which;
      $nzc_blk      = $nzc_blk->index($ciwhich_blk);

      $nnzc_blk        = $nzc_blk->nelem;
      $nnzc           += $nnzc_blk;
      $nnzc_slice_blk  = "${nnzc_prev}:".($nnzc-1);

      ##-- construct block output pdls: ixc
      $ixa_blk = $ixa->dice_axis(1,$nzai_blk->index($ciwhich_blk));
      $ixb_blk = $ixb->dice_axis(1,$nzbi_blk->index($ciwhich_blk));
      foreach (0..$#cdims) {
	if ($_ <= $#adims && $cdims[$_]==$adims[$_]) {
	  $ixc->slice("($_),$nnzc_slice_blk") .= $ixa_blk->slice("($_),");
	} else {
	  $ixc->slice("($_),$nnzc_slice_blk") .= $ixb_blk->slice("($_),");
	}
      }

      ##-- construct block output pdls: nzc
      $nzc->slice($nnzc_slice_blk) .= $nzc_blk;

      ##-- maybe allocate for another block
      if ($nzai_cur < $nixa || $nzbi_cur < $nixb) {
	$nzci_nxt -= $nzci_cur;
	$nzci_cur  = 0;

	$ixc = $ixc->reshape($ixc->dim(0), $ixc->dim(1)+$nzci_nxt+$blksz);
	$nzc = $nzc->reshape($nzc->dim(0)+$nzci_nxt+$blksz);
	$nzai = $nzai->reshape($nzci_nxt+$blksz) if ($nzci_nxt+$blksz > $nzai->dim(0));
	$nzbi = $nzbi->reshape($nzci_nxt+$blksz) if ($nzci_nxt+$blksz > $nzbi->dim(0));

	$istate .= $ostate;
	$istate->set(4, $nzci_cur);
	$istate->set(5, $nzci_nxt);
      }
      $nnzc_prev = $nnzc;

    } while ($nzai_cur < $nixa || $nzbi_cur < $nixb);

    ##-- set up final output pdl
    my $ixc_sorti = $ixc->qsortveci;
    $nzc          = $nzc->index($ixc_sorti)->append($zc);
    $nzc->sever;

    $ixc          = $ixc->dice_axis(1,$ixc_sorti);
    $ixc->sever;

    return $a->shadow(
		      dims  => [@cdims],
		      xdims => PDL->sequence($P_LONG,scalar(@cdims)),
		      which => $ixc,
		      vals  => $nzc,
		     );
  };
}

foreach my $binop (
		   qw(plus minus mult divide modulo power),
		   qw(gt ge lt le eq ne spaceship),
		  )
  {
    eval "*${binop} = *${binop}_mia = _binary_op_mia('${binop}',PDL->can('${binop}'));";
  }

foreach my $intop (
		   qw(and2 or2 xor shiftleft shiftright),
		  )
  {
    eval "*${intop} = *${intop}_mia = _binary_op_mia('${intop}',PDL->can('${intop}'),\$P_LONG);";
  }


##--------------------------------------------------------------
## General Information

## $density = $ccs->density()
##  + returns PDL density as a scalar (lower is more sparse)
sub density { $_[0][$WHICH]->dim(1) / $_[0]->nelem; }

## $compressionRate = $ccs->compressionRate()
##  + higher is better
##  + negative value indicates that dense storage would be more memory-efficient
##  + pointers aren't included in the statistics: just which,nzvals,missing
sub compressionRate {
  my $ccs     = shift;
  my $dsize   = PDL->pdl($ccs->nelem) * PDL::howbig($ccs->type);
  my $ccssize = (0
		 + PDL->pdl($ccs->[$WHICH]->nelem) * PDL::howbig($ccs->[$WHICH]->type)
		 + PDL->pdl($ccs->[$VALS]->nelem)  * PDL::howbig($ccs->[$VALS]->type)
		 + PDL->pdl($ccs->[$XDIMS]->nelem) * PDL::howbig($ccs->[$XDIMS]->type)
		);
  return (($dsize - $ccssize) / $dsize)->sclr;
}

##--------------------------------------------------------------
## Stringification & Viewing

## $str = $obj->string()
sub string {
  my $whichstr = ''.$_[0][$WHICH]->dice_axis(0,$_[0][$XDIMS])->xchg(0,1);
  $whichstr =~ s/^([^A-Z])/$,  $1/mg;
  #$whichstr =~ s/^\s*?\n\s*//mg;
  chomp($whichstr);
  return
    (
     ''
     .ref($_[0])."(".join(',', @{$_[0][$DIMS]}[$_[0][$XDIMS]->list]).")\n"
     .$,." which:(".join(',', $_[0][$WHICH]->type, $_[0][$WHICH]->dims).")^T=$whichstr\n"
     .$,." vals:(" .join(',', $_[0][$VALS]->type,  $_[0][$VALS]->dims).")=$_[0][$VALS]\n"
    );
}

## $pstr = $obj->lstring()
##  + literal perl-type string (re-blesses the object)
sub lstring {
  my $ccs = shift;
  my $ref = ref($ccs);
  my $str = ''.(bless $ccs, 'ARRAY');
  bless($ccs,$ref);
  return $str;
}


##======================================================================
## AUTOLOAD: pass to nonzero-PDL
##  + doesn't seem to work well
##======================================================================

#our $AUTOLOAD;
#sub AUTOLOAD {
#  my $d = shift;
#  return undef if (!defined($d) || !defined($d->[$VALS]));
#  (my $name = $AUTOLOAD) =~ s/.*:://; ##-- strip qualification
#  my ($sub);
#  if (!($sub=UNIVERSAL::can($d->[$VALS],$name))) {
#    croak( ref($d) , "::$name() not defined for nzvals in ", __PACKAGE__ , "::AUTOLOAD.\n");
#  }
#  return $sub->($d->[$VALS],@_);
#}

1; ##-- make perl happy


##======================================================================
## pod: headers
=pod

=head1 NAME

PDL::CCS::Nd - N-dimensional CCS-encoded PDLs

=head1 SYNOPSIS

 use PDL;
 use PDL::CCS::Nd;

 ##---------------------------------------------------------------------
 ## ... stuff happens

=cut

##======================================================================
## TODO: write pods

