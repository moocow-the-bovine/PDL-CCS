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
##     dims  => \@dims,   ##-- defaul: [@$dims1]
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
     ."PDL::CCS(".join(',', @{$_[0][$DIMS]}[$_[0][$XDIMS]->list]).")\n"
     .$,." which:(".join(',', $_[0][$WHICH]->type, $_[0][$WHICH]->dims).")^T=$whichstr\n"
     .$,." nzvals:(" .join(',', $_[0][$VALS]->type,  $_[0][$VALS]->dims).")=$_[0][$VALS]\n"
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

