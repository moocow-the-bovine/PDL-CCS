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
our @ISA = ('PDL::Exporter');
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
## Global variables for block-wise computation of binary operations

##-- some (hopefully sensible) defaults
#our $BINOP_BLOCKSIZE_MIN =    64;
#our $BINOP_BLOCKSIZE_MAX = undef; ##-- undef or zero: no maximum

##-- debug/devel defaults
our $BINOP_BLOCKSIZE_MIN  =      1;
our $BINOP_BLOCKSIZE_MAX  =      0;

##======================================================================
## Globals

our $DIMS    = 0;
our $XDIMS   = 1; ##-- dimension translation pdl: $whichND_logical = $WHICH->dice_axis(0,$XDIMS)
our $WHICH   = 2;
our $VALS    = 3;
our $PTRS    = 4;
our $FLAGS   = 5;

##-- flags
our $CCSND_BAD_IS_MISSING = 1;
our $CCSND_NAN_IS_MISSING = 2;
our $CCSND_INPLACE        = 4;
our $FLAGS_DEFAULT        = 0; ##-- default flags

##-- pdl constants
our $P_BYTE = PDL::byte();
our $P_LONG = PDL::long();


##======================================================================
## Constructors etc.

## $obj = $class_or_obj->newFromDense($denseND);
## $obj = $class_or_obj->newFromDense($denseND,$missing);
## $obj = $class_or_obj->newFromDense($denseND,$missing,$flags);
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
##     $FLAGS   => $flags,     ##-- integer holding some flags
##
##  + each element of @PTRS is itself an array:
##     $PTRS[$i] => [ $PTR, $NZI ]
##
sub newFromDense {
  my $that = shift;
  return (bless [], ref($that)||$that)->fromDense(@_);
}

## $obj = $obj->fromDense($denseND,$missing,$flags)
sub fromDense {
  my ($obj,$p,$missing,$flags) = @_;
  $p = PDL->topdl($p);
  $p = $p->slice("*1") if (!$p->dims);
  $missing     = (defined($missing)
		  ? PDL->pdl($p->type,$missing)
		  : ($p->badflag
		     ? PDL->pdl($p->type,0)->setvaltobad(0)
		     : PDL->pdl($p->type,0)));
  $flags = $FLAGS_DEFAULT if (!defined($flags));
  my $pwhichND = ($missing->isbad ? $p->isgood() : ($p != $missing))->whichND->vv_qsortvec;
  my $pnz  = $p->indexND($pwhichND)->append($missing);
  $pnz->sever;                       ##-- always sever nzvals ?
  $obj->[$DIMS]    = [$p->dims];
  $obj->[$XDIMS]   = @{$obj->[$DIMS]} ? PDL->sequence($P_LONG,scalar(@{$obj->[$DIMS]})) : PDL->null->long;
  $obj->[$WHICH]   = $pwhichND;
  $obj->[$VALS]    = $pnz;
  $obj->[$PTRS]    = [];            ##-- do we really need this ... yes
  $obj->[$FLAGS]   = $flags;
  return $obj;
}

## $obj = $class_or_obj->newFromWhich($whichND,$nzvals,%options);
## $obj = $class_or_obj->newFromWhich($whichND,$nzvals);
##  + %options:
##     dims    => \@dims,   ##-- dimension list; default guessed from $whichND
##     missing => $missing, ##-- default: BAD if $nzvals->badflag, 0 otherwise
##     xdims   => $xdims,   ##-- dimension translation PDL (default: sequence($ndims)) : $LOGICAL_DIM => $PHYSICAL_DIM
##     flags   => $flags,   ##-- flags
##  + no dataflow!
sub newFromWhich {
  my $that = shift;
  return bless([],ref($that)||$that)->fromWhich(@_);
}

sub fromWhich {
  my ($obj,$wnd,$nzvals,%opts) = @_;
  my $missing = (defined($opts{missing})
		 ? PDL->pdl($nzvals->type,$opts{missing})
		 : ($nzvals->badflag
		    ? PDL->pdl($nzvals->type,0)->setvaltobad(0)
		    : PDL->pdl($nzvals->type,0)));
  my $dims  = ( defined($opts{dims})  ? $opts{dims}  : [($wnd->xchg(0,1)->maximum+1)->list]  );
  my $xdims = ( defined($opts{xdims}) ? $opts{xdims} : PDL->sequence($P_LONG,scalar(@$dims)) );
  $wnd     = $wnd->dice_axis(0,$xdims);
  my $wi   = $wnd->isempty ? PDL->null->long : $wnd->qsortveci;
  $wnd     = $wnd->dice_axis(1,$wi)->pdl();         ##-- copy!
  $nzvals  = $nzvals->index($wi)->append($missing); ##-- copy (b/c append)
  $obj->[$DIMS]    = $dims;
  $obj->[$XDIMS]   = @{$obj->[$DIMS]} ? PDL->sequence($P_LONG,scalar(@{$obj->[$DIMS]})) : PDL->null->long;
  $obj->[$WHICH]   = $wnd;
  $obj->[$VALS]    = $nzvals;
  $obj->[$PTRS]    = [];
  $obj->[$FLAGS]   = defined($opts{flags}) ? $opts{flags} : $FLAGS_DEFAULT;
  return $obj;
}


## DESTROY : avoid PDL inheritance
sub DESTROY { ; }


## $ccs = $pdl->toccs()
## $ccs = $pdl->toccs($missing)
## $ccs = $pdl->toccs($missing,$flags)
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
  $ccs2->[$DIMS]  = [ @{$ccs1->[$DIMS]} ];
  $ccs2->[$XDIMS] = $ccs1->[$XDIMS]->pdl;
  $ccs2->[$WHICH] = $ccs1->[$WHICH]->pdl;
  $ccs2->[$VALS]  = $ccs1->[$VALS]->pdl;
  $ccs2->[$PTRS]  = [ map {defined($_) ? [map {$_->pdl} @$_] : undef} @{$ccs1->[$PTRS]} ]; ##-- copy pointers?
  $ccs2->[$FLAGS] = $ccs1->[$FLAGS];
  return $ccs2;
}

## $ccs2 = $ccs->copyShallow()
##  + a very shallow version of copy()
##  + Copied    : @$DIMS, @$PTRS, @{$PTRS->[*]}, $FLAGS
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
##     to    => $ccs2,    ##-- default: new
##     dims  => \@dims,   ##-- default: [@$dims1]
##     xdims => $xdims2,  ##-- default: $xdims1->pdl
##     ptrs  => \@ptrs2,  ##-- default: []
##     which => $which2,  ##-- default: undef
##     vals  => $vals2,   ##-- default: undef
##     flags => $flags,   ##-- default: $flags1
sub shadow {
  my ($ccs,%args) = @_;
  my $ccs2        = defined($args{to}) ? $args{to} : bless([], ref($ccs)||$ccs);
  $ccs2->[$DIMS]  = defined($args{dims})  ? $args{dims}  : [@{$ccs->[$DIMS]}];
  $ccs2->[$XDIMS] = defined($args{xdims}) ? $args{xdims} : $ccs->[$XDIMS]->pdl;
  $ccs2->[$PTRS]  = $args{ptrs}  ? $args{ptrs} : [];
  $ccs2->[$WHICH] = $args{which};
  $ccs2->[$VALS]  = $args{vals};
  $ccs2->[$FLAGS] = defined($args{flags}) ? $args{flags} : $ccs->[$FLAGS];
  return $ccs2;
}


##--------------------------------------------------------------
## Maintainence

## $ccs = $ccs->recode()
##  + recodes object, removing any missing values from $nzvals
sub recode {
  my $ccs = shift;
  my $nz = $ccs->[$VALS]->slice("0:-2");
  my $z  = $ccs->[$VALS]->slice("-1");

  ##-- get mask of "real" non-zero values
  my $nzmask = PDL->zeroes($P_BYTE,$nz->nelem);
  my ($nzmask1);
  if ($z->isbad) {
    $nz->isgood($nzmask);
  }
  else {
    $nz->ne($z, $nzmask, 0);
    if ($ccs->[$FLAGS] & $CCSND_BAD_IS_MISSING) {
      $nzmask1 = $nzmask->pdl;
      $nz->isgood($nzmask1);
      $nzmask &= $nzmask1;
    }
  }
  if ($ccs->[$FLAGS] & $CCSND_NAN_IS_MISSING) {
    $nzmask1 = $nzmask->pdl if (!defined($nzmask1));
    $nz->isfinite($nzmask1);
    $nzmask &= $nzmask1;
  }

  ##-- maybe recode
  if (!$nzmask->all) {
    my $nzi = $nzmask->which;
    $ccs->[$WHICH]   = $ccs->[$WHICH]->dice_axis(1,$nzi);
    $ccs->[$VALS]    = $ccs->[$VALS]->index($nzi)->append($z);
    @{$ccs->[$PTRS]} = qw(); ##-- clear pointers
  }

  return $ccs;
}

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
	     ($_[0][$VALS]->dim(0) > 1 ? $_[0][$VALS]->slice("0:-2") : $_[0][$VALS]),
	     $_[0][$VALS]->slice("(-1)"),
	     [@{$_[0][$DIMS]}[$_[0][$XDIMS]->list]],
	     @_[1..$#_]);
}

## $dense = $ccs_or_dense->todense()
*PDL::todense = \&todense;
sub todense { isa($_[0],__PACKAGE__) ? $_[0]->decode() : $_[0]; }


##--------------------------------------------------------------
## PDL API: Basic Properties

## $type = $obj->type()
sub type { $_[0][$VALS]->type; }

## $obj2 = $obj->convert($type)
##  + unlike PDL function, respects 'inplace' flag
sub convert {
  if ($_[0][$FLAGS] & $CCSND_INPLACE) {
    $_[0][$VALS]   = $_[0][$VALS]->convert($_[1]);
    $_[0][$FLAGS] &= ~$CCSND_INPLACE;
    return $_[0];
  }
  $_[0]->shadow(which=>$_[0][$WHICH]->pdl, vals=>$_[0][$VALS]->convert($_[1]));
}

## byte,short,ushort,long,double,...
sub _pdltype_sub {
  my $pdltype = shift;
  return sub { return $pdltype if (!@_); convert(@_,$pdltype); };
}
foreach my $pdltype (qw(byte short ushort long longlong float double)) {
  eval "*${pdltype} = _pdltype_sub(PDL::${pdltype}());";
}

## @dims = $obj->dims()
sub dims { @{$_[0][$DIMS]}[$_[0][$XDIMS]->list]; }

## $dim = $obj->dim($dimi)
*getdim = \&dim;
sub dim { $_[0][$DIMS][$_[0][$XDIMS]->at($_[1])]; }

## $ndims = $obj->ndims()
sub ndims { scalar(@{$_[0][$DIMS]}); }

## $nelem = $obj->nelem
sub nelem { PDL->pdl($P_LONG,$_[0][$DIMS])->dprod; }

## $bool = $obj->isnull
sub isnull { $_[0][$VALS]->isnull; }

## $bool = $obj->isempty
sub isempty { $_[0][$VALS]->isempty; }

##--------------------------------------------------------------
## Low-level CCS access

## $xdims = $obj->xdims()
sub xdims { $_[0][$XDIMS]; }

## $nnz = $obj->nstored()
sub nstored { $_[0][$WHICH]->dim(1); }

## $nnz = $obj->_nnz
##   + returns actual $obj->[$VALS]->dim(0)-1
sub _nnz { $_[0][$VALS]->dim(0)-1; }

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
  return $_[0][$VALS]->index(null) if ($_[0][$VALS]->dim(0)<=1);
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

## $flags = $obj->flags()
## $flags = $obj->flags($flags)
##  + get local flags
sub flags { $_[0][$FLAGS] = $_[1] if (@_ > 1); $_[0][$FLAGS]; }

## $bool = $obj->bad_is_missing()
## $bool = $obj->bad_is_missing($bool)
sub bad_is_missing {
  if (@_ > 1) {
    if ($_[1]) { $_[0][$FLAGS] |=  $CCSND_BAD_IS_MISSING; }
    else       { $_[0][$FLAGS] &= ~$CCSND_BAD_IS_MISSING; }
  }
  $_[0][$FLAGS] & $CCSND_BAD_IS_MISSING;
}

## $bool = $obj->nan_is_missing()
## $bool = $obj->nan_is_missing($bool)
sub nan_is_missing {
  if (@_ > 1) {
    if ($_[1]) { $_[0][$FLAGS] |=  $CCSND_NAN_IS_MISSING; }
    else       { $_[0][$FLAGS] &= ~$CCSND_NAN_IS_MISSING; }
  }
  $_[0][$FLAGS] & $CCSND_NAN_IS_MISSING;
}

## undef = $obj->set_inplace($bool)
##   + sets local inplace flag
sub set_inplace ($$) {
  if ($_[1]) { $_[0][$FLAGS] |=  $CCSND_INPLACE; }
  else       { $_[0][$FLAGS] &= ~$CCSND_INPLACE; }
}

## $bool = $obj->is_inplace()
sub is_inplace ($) { ($_[0][$FLAGS] & $CCSND_INPLACE) ? 1 : 0; }

## $obj = $obj->inplace()
##   + sets local inplace flag
sub inplace ($) { $_[0][$FLAGS] |= $CCSND_INPLACE; $_[0]; }

## $obj = $obj->sever()
##  + severs all sub-pdls
sub sever {
  $_[0][$XDIMS]->sever;
  $_[0][$WHICH]->sever;
  $_[0][$VALS]->sever;
  foreach (grep {defined($_)} (@{$_[0][$PTRS]})) {
    $_->[0]->sever;
    $_->[1]->sever;
  }
  $_[0];
}

## \&code = _twiddle_bad_sub($pdlcode)
sub _twiddle_bad_sub {
  my $pdlsub = shift;
  return sub {
    if ($_[0]->is_inplace) {
      $pdlsub->($_[0][$VALS]->inplace, @_[1..$#_]);
      $_[0]->set_inplace(0);
      return $_[0];
    }
    $_[0]->shadow(
		  which=>$_[0][$WHICH]->pdl,
		  vals=>$pdlsub->($_[0][$VALS],@_[1..$#_]),
		 );
  };
}

## $obj = $obj->setnantobad()
foreach my $badsub (qw(setnantobad setbadtonan setbadtoval setvaltobad)) {
  eval "*${badsub} = _twiddle_bad_sub(PDL->can('$badsub'));";
}

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
##  + always copies
sub transpose { $_[0]->xchg(0,1)->copy; }



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

## $ufunc_sub = _ufuncsub($subname, \&ccs_accum_sub, $allow_bad_missing)
sub _ufuncsub {
  my ($subname,$accumsub,$allow_bad_missing) = @_;
  barf(__PACKAGE__, "::_ufuncsub($subname): no underlying CCS accumulator func!") if (!defined($accumsub));
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
					($allow_bad_missing || $missing->isgood ? ($missing,$dims[0]) : (0,0))
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

*nbadover  = _ufuncsub('nbadover',  PDL::CCS::Ufunc->can('ccs_accum_nbad'), 1);
*ngoodover = _ufuncsub('ngoodover', PDL::CCS::Ufunc->can('ccs_accum_ngood'), 1);
*nnz       = _ufuncsub('nnz', PDL::CCS::Ufunc->can('ccs_accum_nnz'), 1);

sub sum   { my $z=$_[0]->missing; $_[0][$VALS]->slice("0:-2")->sum   + ($z->isgood ? ($z->sclr *  $_[0]->nmissing) : 0); }
sub dsum  { my $z=$_[0]->missing; $_[0][$VALS]->slice("0:-2")->dsum  + ($z->isgood ? ($z->sclr *  $_[0]->nmissing) : 0); }
sub prod  { my $z=$_[0]->missing; $_[0][$VALS]->slice("0:-2")->prod  * ($z->isgood ? ($z->sclr ** $_[0]->nmissing) : 1); }
sub dprod { my $z=$_[0]->missing; $_[0][$VALS]->slice("0:-2")->dprod * ($z->isgood ? ($z->sclr ** $_[0]->nmissing) : 1); }
sub min   { $_[0][$VALS]->min; }
sub max   { $_[0][$VALS]->max; }

sub nbad  { my $z=$_[0]->missing; $_[0][$VALS]->slice("0:-2")->nbad   + ($z->isbad  ? $_[0]->nmissing : 0); }
sub ngood { my $z=$_[0]->missing; $_[0][$VALS]->slice("0:-2")->ngood  + ($z->isgood ? $_[0]->nmissing : 0); }

sub any { $_[0][$VALS]->any; }
sub all { $_[0][$VALS]->all; }


##--------------------------------------------------------------
## Unary Operations

## $sub = _unary_op($opname,$pdlsub)
sub _unary_op {
  my ($opname,$pdlsub) = @_;
  return sub {
    if ($_[0]->is_inplace) {
      $pdlsub->($_[0][$VALS]->inplace);
      $_[0]->set_inplace(0);
      return $_[0];
    }
    $_[0]->shadow(which=>$_[0][$WHICH]->pdl, vals=>$pdlsub->($_[0][$VALS]));
  };
}

foreach my $unop (qw(bitnot sqrt abs sin cos not exp log log10))
  {
    eval "*${unop} = _unary_op('${unop}',PDL->can('${unop}'));";
  }

##--------------------------------------------------------------
## Binary Operations: missing-is-annihilator

## \@matchdims_or_undef = _ccsnd_binop_relevant_dims(\@dims1,\@dims2)
##  + barf()s on mismath
sub _ccsnd_binop_relevant_dims {
  my ($opname,$dims1,$dims2) = @_;
  my $rdims = [];
  foreach (0..($#$dims1 < $#$dims2 ? $#$dims1 : $#$dims2)) {
    next if ($dims1->[$_] <= 1 || $dims2->[$_] <= 1);
    if ($dims1->[$_] != $dims2->[$_]) {
      barf("PDL::CCS::Nd::",
	   ($opname||'_ccsnd_binop_relevant_dims'),
	   "(): mismatch on non-trivial dim($_): $dims1->[$_] != $dims2->[$_]");
    }
    push(@$rdims,$_);
  }
  return $rdims;
}

sub _ccsnd_binary_op_mia {
  my ($opname,$pdlsub,$deftype) = @_;

  return sub {
    my ($a,$b,$swap) = @_;
    $swap=0 if (!defined($swap));

    ##-- check for scalar operations
    if (!ref($b) || $b->nelem==1) {
      if ($a->is_inplace) {
	$pdlsub->($a->[$VALS]->inplace, todense($b), $swap);
	$a->set_inplace(0);
	return $a->recode;
      }
      return $a->shadow(
			which => $a->[$WHICH]->pdl,
			vals  => $pdlsub->($a->[$VALS], todense($b), $swap),
		       )->recode;
    }

    ##-- convert b to CCS
    $b = toccs($b);

    ##-- get relevant dimensions
    my @adims = $a->dims;
    my @bdims = $b->dims;
    my @cdims = (
		 map {
		   $_ <= $#adims && $adims[$_] > 1 ? $adims[$_] : $bdims[$_]
		 } (0..($#adims > $#bdims ? $#adims : $#bdims))
		);
    my $rdims = _ccsnd_binop_relevant_dims($opname,\@adims,\@bdims);
    my $rdpdl = PDL->pdl($P_LONG,$rdims);

    ##-- get & sort relevant indices
    my $ixa        = $a->[$WHICH]->dice_axis(0,$a->[$XDIMS]);
    my $nixa       = $ixa->dim(1);
    my $ixar       = $ixa->dice_axis(0,$rdpdl);
    my $ixar_sorti = $ixar->isempty ? PDL->null->long : $ixar->qsortveci;
    $ixa           = $ixa->dice_axis(1,$ixar_sorti);
    $ixar          = $ixar->dice_axis(1,$ixar_sorti);
    ##
    my $ixb        = $b->[$WHICH]->dice_axis(0,$b->[$XDIMS]);
    my $nixb       = $ixb->dim(1);
    my $ixbr       = $ixb->dice_axis(0,$rdpdl);
    my $ixbr_sorti = $ixbr->isempty ? PDL->null->long : $ixbr->qsortveci;
    $ixb           = $ixb->dice_axis(1,$ixbr_sorti);
    $ixbr          = $ixbr->dice_axis(1,$ixbr_sorti);

    ##-- initialize: values
    my $avals  = $a->[$VALS];
    my $avalsr = $avals->index($ixar_sorti);
    my $bvals  = $b->[$VALS];
    my $bvalsr = $bvals->index($ixbr_sorti);

    ##-- initialize: state vars
    my $blksz  = $nixa > $nixb ? $nixa : $nixb;
    $blksz     = $BINOP_BLOCKSIZE_MIN if ($BINOP_BLOCKSIZE_MIN && $blksz < $BINOP_BLOCKSIZE_MIN);
    $blksz     = $BINOP_BLOCKSIZE_MAX if ($BINOP_BLOCKSIZE_MAX && $blksz > $BINOP_BLOCKSIZE_MAX);
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
    my $nanismissing = ($a->[$FLAGS]&$CCSND_NAN_IS_MISSING);
    my $badismissing = ($a->[$FLAGS]&$CCSND_BAD_IS_MISSING);
    $zc              = $zc->setnantobad() if ($nanismissing && $badismissing);
    my $zc_isbad     = $zc->isbad ? 1 : 0;

    ##-- block-wise variables
    my ($nzai_prv,$nzai_pnx, $nzbi_prv,$nzbi_pnx, $nzci_prv,$nzci_pnx,$cmpval_prv);
    my ($nzai_cur,$nzai_nxt, $nzbi_cur,$nzbi_nxt, $nzci_cur,$nzci_nxt,$cmpval);
    my ($nzci_max, $blk_slice, $nnzc_blk,$nnzc_slice_blk);
    my ($nzai_blk,$nzbi_blk,$ixa_blk,$ixb_blk,$nzc_blk,$cimask_blk,$ciwhich_blk);
    my $nnzc_prev=0;
    do {
      ##-- align a block of data
      ccs_binop_align_block_mia($ixar,$ixbr,$istate, $nzai,$nzbi,$ostate);

      ##-- parse current alignment algorithm state
      ($nzai_prv,$nzai_pnx, $nzbi_prv,$nzbi_pnx, $nzci_prv,$nzci_pnx,$cmpval_prv) = $istate->list;
      ($nzai_cur,$nzai_nxt, $nzbi_cur,$nzbi_nxt, $nzci_cur,$nzci_nxt,$cmpval)     = $ostate->list;
      $nzci_max = $nzci_cur-1;

      if ($nzci_max >= 0) {
	##-- construct block output pdls: nzvals
	$blk_slice = "${nzci_prv}:${nzci_max}";
	$nzai_blk  = $nzai->slice($blk_slice);
	$nzbi_blk  = $nzbi->slice($blk_slice);
	$nzc_blk   = $pdlsub->($avalsr->index($nzai_blk), $bvalsr->index($nzbi_blk), $swap);

	##-- get indices of non-$missing c() values
	$cimask_blk   = $zc_isbad ? $nzc_blk->isgood : ($nzc_blk!=$zc);
	$cimask_blk  &= $nzc_blk->isgood   if (!$zc_isbad && $badismissing);
	$cimask_blk  &= $nzc_blk->isfinite if ($nanismissing);
	if ($cimask_blk->any) {
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
	}
      }

      ##-- possibly allocate for another block
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

    ##-- trim output pdls
    if ($nnzc > 0) {
      ##-- usual case: some values are non-missing
      $ixc = $ixc->slice(",0:".($nnzc-1));
      my $ixc_sorti = $ixc->qsortveci;
      $nzc          = $nzc->index($ixc_sorti)->append($zc);
      $nzc->sever;
      $ixc          = $ixc->dice_axis(1,$ixc_sorti);
      $ixc->sever;
    } else {
      ##-- pathological case: all values are "missing"
      $ixc = $ixc->dice_axis(1,PDL->pdl([]));
      $ixc->sever;
      $nzc = $zc;
    }

    ##-- set up final output object
    my $c = $a->shadow(
		       dims  => [@cdims],
		       xdims => PDL->sequence($P_LONG,scalar(@cdims)),
		       which => $ixc,
		       vals  => $nzc,
		      );
    if ($a->is_inplace) {
      @$a = @$c;
      $a->set_inplace(0);
      return $a;
    }
    return $c;
  };
}

foreach my $binop (
		   qw(plus minus mult divide modulo power),
		   qw(gt ge lt le eq ne spaceship),
		  )
  {
    eval "*${binop} = *${binop}_mia = _ccsnd_binary_op_mia('${binop}',PDL->can('${binop}'));";
  }

foreach my $intop (
		   qw(and2 or2 xor shiftleft shiftright),
		  )
  {
    eval "*${intop} = *${intop}_mia = _ccsnd_binary_op_mia('${intop}',PDL->can('${intop}'),\$P_LONG);";
  }

*rassgn_mia = _ccsnd_binary_op_mia('rassgn', sub { PDL::assgn($_[1],$_[0]); $_[1]; });

## $to = $to->rassgn($from)
##  + calls newFromDense() with $to flags if $from is dense
##  + otherwise, copies $from to $to
##  + calling conventions are the OPPOSITE of PDL::assgn
sub rassgn {
  my ($to,$from) = @_;
  if (!ref($from) || $from->nelem==1) {
    ##-- assignment from a scalar: treat the Nd object as a mask of available values
    $to->[$VALS] .= todense($from);
    return $to;
  }
  if (isa($from,__PACKAGE__)) {
    ##-- assignment from a CCS object: check for full dim match or an empty "$to"
    my $fromdimp = $from->dimpdl;
    my $todimp   = $to->dimpdl;
    if ( $to->[$VALS]->dim(0)<=1 || $todimp->isempty || ($fromdimp==$todimp)->all ) {
      @$to = @{$from->copy};
      return $to;
    }
  }
  ##-- $from is something else: pass it on to 'rassgn_mia': effectively treat $to->[$WHICH] as a mask for $from
  $to->[$FLAGS] |= $CCSND_INPLACE;
  return $to->rassgn_mia($from);
}

## $to = $from->assgn($to)
sub assgn { $_[1]->rassgn($_[0]); }


##--------------------------------------------------------------
## CONTINUE HERE

## TODO:
##  + virtual dimensions: dummy & clump
##  + OPERATIONS:
##    - matrix: matmult
##    - accumulators: (some still missing: statistical, extrema-indices, atan2, ...)


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

##--------------------------------------------------------------
## Operator overloading

use overload (
	      ##-- Binary ops: arithmetic
	      "+" => \&plus_mia,
	      "-" => \&minus_mia,
	      "*" => \&mult_mia,
	      "/" => \&divide_mia,
	      "%" => \&modulo_mia,
	      "**"  => \&power_mia,
	      '+='  => sub { $_[0]->inplace->plus_mia(@_[1..$#_]); },
	      '-='  => sub { $_[0]->inplace->minus_mia(@_[1..$#_]); },
	      '*='  => sub { $_[0]->inplace->mult_mia(@_[1..$#_]); },
	      '%='  => sub { $_[0]->inplace->divide_mia(@_[1..$#_]); },
	      '**=' => sub { $_[0]->inplace->modulo_mia(@_[1..$#_]); },

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
	      "^"  => \&xor_mia,
	      "<<" => \&shiftleft_mia,
	      ">>" => \&shiftright_mia,
	      '|='  => sub { $_[0]->inplace->or2_mia(@_[1..$#_]); },
	      '&='  => sub { $_[0]->inplace->and2_mia(@_[1..$#_]); },
	      '^='  => sub { $_[0]->inplace->xor_mia(@_[1..$#_]); },
	      '<<=' => sub { $_[0]->inplace->shiftleft_mia(@_[1..$#_]); },
	      '>>=' => sub { $_[0]->inplace->shiftright_mia(@_[1..$#_]); },

	      ##-- Unary operations
	      "!"  => \&not,
	      "~"  => \&bitnot,
	      "sqrt" => \&sqrt,
	      "abs"  => \&abs,
	      "sin"  => \&sin,
	      "cos"  => \&cos,
	      "log"  => \&log,
	      "exp"  => \&exp,

	      ##-- assignment & assigning variants
	      ".=" => \&rassgn,


	      ##-- Stringification & casts
	      'bool' => sub {
		my $nelem = $_[0]->nelem;
		return 0 if ($nelem==0);
		croak("multielement ", __PACKAGE__, " pseudo-piddle in conditional expression") if ($nelem!=1);
		$_[0][$VALS]->at(0);
	      },
	      "\"\"" => \&string,
	     );


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

