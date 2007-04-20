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
our @ISA = qw();
our @EXPORT_OK =
  (
   ##-- Encoding/Decoding
   qw(toccs todense),
  );
our %EXPORT_TAGS =
  (
   Func => [@EXPORT_OK],               ##-- respect PDL conventions (hopefully)
  );
our @EXPORT = @{$EXPORT_TAGS{Func}};

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

our $PDIMS   = 0;
our $VDIMS   = 1;
our $WHICH   = 2;
our $VALS    = 3;
our $PTRS    = 4;
our $FLAGS   = 5;
our $USER    = 6;

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
##     $PDIMS   => $pdims,     ##-- pdl(long,$NPdims)     : physical dimension sizes : $pdim_i => $dimSize_i
##     $VDIMS   => $vdims,     ##-- pdl(long,$NVdims)     : virtual dimension sizes
##                             ##     + $vdim_i => / -$vdimSize_i   if $vdim_i is dummy
##                             ##                  \  $pdim_i       otherwise
##                             ##     + s.t. $whichND_logical_physical = $whichND->dice_axis(0,$vdims->where($vdims>=0));
##     $WHICH   => $whichND,   ##-- pdl(long,$NPdims,$Nnz) ~ $dense_orig->whichND
##                             ##   + guaranteed to be sorted as for qsortvec() specs
##                             ##   + NOT changed by dimension-shuffling transformations
##     $VALS    => $vals,      ##-- pdl( ?  ,$Nnz+1)      ~ $dense->where($dense)->append($missing)
##     $PTRS    => \@PTRS,     ##-- array of ccsutils-pointers by physical dimension number
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
  my $pdims = PDL->pdl($P_LONG,[$p->dims]);
  $obj->[$PDIMS]   = $pdims;
  $obj->[$VDIMS]   = $pdims->isempty ? $pdims->pdl : $pdims->sequence;
  $obj->[$WHICH]   = $pwhichND;
  $obj->[$VALS]    = $pnz;
  $obj->[$PTRS]    = [];            ##-- do we really need this ... yes
  $obj->[$FLAGS]   = $flags;
  return $obj;
}

## $obj = $class_or_obj->newFromWhich($whichND,$nzvals,%options);
## $obj = $class_or_obj->newFromWhich($whichND,$nzvals);
##  + %options: see $obj->fromWhich()
sub newFromWhich {
  my $that = shift;
  return bless([],ref($that)||$that)->fromWhich(@_);
}

## $obj = $obj->fromWhich($whichND,$nzvals,%options);
## $obj = $obj->fromWhich($whichND,$nzvals);
##  + %options:
##     sorted  => $bool,    ##-- if true, $whichND is assumed to be pre-sorted
##     steal   => $bool,    ##-- if true, $whichND and $nzvals are used literally (implies 'sorted')
##                          ##    + in this case, $nzvals should really be: $nzvals->append($missing)
##     pdims   => $pdims,   ##-- physical dimension list; default guessed from $whichND (alias: 'dims')
##     missing => $missing, ##-- default: BAD if $nzvals->badflag, 0 otherwise
##     vdims   => $vdims,   ##-- virtual dims (default: sequence($nPhysDims)); alias: 'xdims'
##     flags   => $flags,   ##-- flags
sub fromWhich {
  my ($obj,$wnd,$nzvals,%opts) = @_;
  my $missing = (defined($opts{missing})
		 ? PDL->pdl($nzvals->type,$opts{missing})
		 : ($nzvals->badflag
		    ? PDL->pdl($nzvals->type,0)->setvaltobad(0)
		    : PDL->pdl($nzvals->type,0)));
  $opts{pdims} = $opts{dims}  if (!defined($opts{pdims}) && defined($opts{dims}));
  $opts{vdims} = $opts{xdims} if (!defined($opts{vdims}) && defined($opts{xdims}));
  my $pdims = (defined($opts{pdims})
	       ? $opts{pdims}
	       : (defined($opts{dims})
		  ? $opts{dims}
		  : PDL->pdl($P_LONG, [($wnd->xchg(0,1)->maximum+1)->list]) ));
  my $vdims = (defined($opts{vdims})
	       ? $opts{vdims}
	       : (defined($opts{xdims})
		  ? $opts{xdims}
		  : $pdims->sequence));
  ##-- maybe sort & copy
  if (!$opts{steal}) {
    ##-- not stolen: copy or sever
    if (!$opts{sorted}) {
      my $wi   = $wnd->isempty ? PDL->null->long : $wnd->qsortveci;
      $wnd     = $wnd->dice_axis(1,$wi);
      $nzvals  = $nzvals->index($wi);
    }
    $wnd->sever;                         ##-- sever (~ copy)
    $nzvals = $nzvals->append($missing); ##-- copy (b/c append)
  }
  $obj->[$PDIMS]   = $pdims;
  $obj->[$VDIMS]   = $vdims;
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
  $ccs2->[$PDIMS] = $ccs1->[$PDIMS]->pdl;
  $ccs2->[$VDIMS] = $ccs1->[$VDIMS]->pdl;
  $ccs2->[$WHICH] = $ccs1->[$WHICH]->pdl;
  $ccs2->[$VALS]  = $ccs1->[$VALS]->pdl;
  $ccs2->[$PTRS]  = [ map {defined($_) ? [map {$_->pdl} @$_] : undef} @{$ccs1->[$PTRS]} ]; ##-- copy pointers?
  $ccs2->[$FLAGS] = $ccs1->[$FLAGS];
  return $ccs2;
}

## $ccs2 = $ccs->copyShallow()
##  + a very shallow version of copy()
##  + Copied    : $PDIMS, @$PTRS, @{$PTRS->[*]}, $FLAGS
##  + Referenced: $VDIMS, $WHICH, $VALS,  $PTRS->[*][*]
sub copyShallow {
  my $ccs = bless [@{$_[0]}], ref($_[0]);
  ##
  ##-- do copy some of it
  $ccs->[$PDIMS]  = $ccs->[$PDIMS]->pdl;
  #$ccs->[$VDIMS] = $ccs->[$VDIMS]->pdl;
  $ccs->[$PTRS]  = [ map {defined($_) ? [@$_] : undef} @{$ccs->[$PTRS]} ];
  $ccs;
}

## $ccs2 = $ccs->shadow(%args)
##  + args:
##     to    => $ccs2,    ##-- default: new
##     pdims => $pdims2,  ##-- default: $pdims1->pdl  (alias: 'dims')
##     vdims => $vdims2,  ##-- default: $vdims1->pdl  (alias: 'xdims')
##     ptrs  => \@ptrs2,  ##-- default: []
##     which => $which2,  ##-- default: undef
##     vals  => $vals2,   ##-- default: undef
##     flags => $flags,   ##-- default: $flags1
sub shadow {
  my ($ccs,%args) = @_;
  my $ccs2        = defined($args{to}) ? $args{to} : bless([], ref($ccs)||$ccs);
  $ccs2->[$PDIMS] = (defined($args{pdims}) ? $args{pdims} : (defined($args{dims})  ? $args{dims}  : $ccs->[$PDIMS]->pdl));
  $ccs2->[$VDIMS] = (defined($args{vdims}) ? $args{vdims} : (defined($args{xdims}) ? $args{xdims} : $ccs->[$VDIMS]->pdl));
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
  return $_[0] if ($_[0][$WHICH]->isempty);
  my $sorti     = $_[0][$WHICH]->qsortveci;
  $_[0][$WHICH] = $_[0][$WHICH]->dice_axis(1,$sorti);
  $_[0][$VALS]  = $_[0][$VALS]->index($sorti->append($_[0][$WHICH]->dim(1)));
#
#-- DANGEROUS: pointer copy
#  foreach (grep {defined($_)} @{$_[0][$PTRS]}) {
#    $_->[1]->index($sorti) .= $_->[1];
#  }
#--/DANGEROUS: pointer copy
#
  @{$_[0][$PTRS]} = qw() if (! ($sorti==PDL->sequence($P_LONG,$sorti->dims))->all );
  return $_[0];
}


##--------------------------------------------------------------
## Decoding

## $dense = $ccs->decode()
## $dense = $ccs->decode($dense)
sub decode {
  ##-- decode physically stored index+value pairs
  my $dense = ccs_decode($_[0][$WHICH],
			 $_[0]->_nzvals,
			 $_[0]->missing,
			 [ $_[0][$PDIMS] ],
			);

  ##-- map physical dims with reorder()
  my $porder = $_[0][$VDIMS]->where($_[0][$VDIMS]>=0);
  $dense = $dense->reorder($porder->list); #if (($porder!=$_[0][$PDIMS]->sequence)->any);

  ##-- map virtual dims with dummy()
  my @vdims = $_[0][$VDIMS]->list;
  foreach (grep {$vdims[$_]<0} (0..$#vdims)) {
    $dense = $dense->dummy($_, -$vdims[$_]);
  }

  ##-- assign if $dense was specified by the user
  if (defined($_[1])) {
    $_[1] .= $dense;
    return $_[1];
  }

  return $dense;
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

## $dimpdl = $obj->dimpdl()
##  + values in $dimpdl are negative for virtual dimensions
sub dimpdl {
  my $dims  = $_[0][$VDIMS]->pdl;
  my $physi = ($_[0][$VDIMS]>=0)->which;
  $dims->index($physi) .= $_[0][$PDIMS]->index($_[0][$VDIMS]->index($physi));
  return $dims;
}

## @dims = $obj->dims()
sub dims { $_[0]->dimpdl->abs->list; }

## $dim = $obj->dim($dimi)
sub dim { $_[0]->dimpdl->abs->at($_[1]); }
*getdim = \&dim;

## $ndims = $obj->ndims()
sub ndims { $_[0][$VDIMS]->nelem; }
*getndims = \&ndims;

## $nelem = $obj->nelem
sub nelem { $_[0]->dimpdl->abs->dprod; }

## $bool = $obj->isnull
sub isnull { $_[0][$VALS]->isnull; }

## $bool = $obj->isempty
sub isempty { $_[0][$VALS]->isempty; }

##--------------------------------------------------------------
## Low-level CCS access

## $bool = $obj->allmissing
##  + true if no non-missing values are stored
sub allmissing { $_[0][$VALS]->nelem <= 1; }

## $pdims = $obj->pdims()
## $vdims = $obj->vdims()
sub pdims { $_[0][$PDIMS]; }
sub vdims { $_[0][$VDIMS]; }


## $nelem_p = $obj->nelem_p : maximum number of physically addressable elements
## $nelem_v = $obj->nelem_v : maximum number of virtually addressable elements
sub nelem_p { $_[0][$PDIMS]->dprod; }
*nelem_v = \&nelem;

## $v_per_p = $obj->_ccs_nvperp() : number of virtual elements per physical element
sub _ccs_nvperp { $_[0][$VDIMS]->where($_[0][$VDIMS]<0)->abs->dprod; }

## $nstored_p = $obj->nstored_p : actual number of physically stored elements
## $nstored_v = $obj->nstored_v : actual number of phyiscally+virtually stored elements
sub nstored_p { $_[0][$WHICH]->dim(1); }
sub nstored_v { $_[0][$WHICH]->dim(1) * $_[0]->_ccs_nvperp; }
*nstored = \&nstored_v;


## $nnz = $obj->_nnz_p : returns actual  $obj->[$VALS]->dim(0)-1
## $nnz = $obj->_nnz_v : returns virtual $obj->[$VALS]->dim(0)-1
sub _nnz_p { $_[0][$VALS]->dim(0)-1; }
sub _nnz_v { ($_[0][$VALS]->dim(0)-1) * $_[0]->_ccs_nvperp; }
*_nnz = \&_nnz_v;

## $nmissing_p = $obj->nmissing_p()
## $nmissing_v = $obj->nmissing_v()
sub nmissing_p { $_[0]->nelem_p - $_[0]->nstored_p; }
sub nmissing_v { $_[0]->nelem_v - $_[0]->nstored_v; }
*nmissing = \&nmissing_v;


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

## $_nzvals = $obj->_nzvals()
## $_nzvals = $obj->_nzvals($nzvals)
##  + phyiscal storage
sub _nzvals {
  $_[0][$VALS]=$_[1]->append($_[0][$VALS]->slice("-1")) if (@_ > 1);
  return $_[0][$VALS]->index(PDL->null) if ($_[0][$VALS]->dim(0)<=1);
  $_[0][$VALS]->slice("0:-2");
}

## $vals = $obj->_vals()
## $vals = $obj->_vals($storedvals)
##  + physical storage only
sub _vals {
  $_[0][$VALS]=$_[1] if (@_ > 1);
  $_[0][$VALS];
}


## $ptr           = $obj->ptr($dim_p); ##-- scalar context
## ($ptr,$pi2nzi) = $obj->ptr($dim_p); ##-- list context
##   + returns cached value in $ccs->[$PTRS][$dim_p] if present
##   + caches value in $ccs->[$PTRS][$dim_p] otherwise
##   + $dim defaults to zero, for compatibility
##   + if $dim is zero, all($pi2nzi==sequence($obj->nstored))
##   + physical dimensions ONLY
sub ptr {
  my ($ccs,$dim) = @_;
  $dim = 0 if (!defined($dim));
  $ccs->[$PTRS][$dim] = [$ccs->getptr($dim)] if (!defined($ccs->[$PTRS][$dim]));
  return wantarray ? @{$ccs->[$PTRS][$dim]} : $ccs->[$PTRS][$dim][0];
}

## ($ptr,$pi2nzi) = $obj->getptr($dim_p);
##  + as for ptr(), but does NOT cache anything, and does NOT check the cache
##  + physical dimensions ONLY
sub getptr { ccs_encode_pointers($_[0][$WHICH]->slice("($_[1]),"), $_[0][$PDIMS][$_[1]]); }

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
  $_[0][$PDIMS]->sever;
  $_[0][$VDIMS]->sever;
  $_[0][$WHICH]->sever;
  $_[0][$VALS]->sever;
  foreach (grep {defined($_)} (@{$_[0][$PTRS]})) {
    $_->[0]->sever;
    $_->[1]->sever;
  }
  $_[0];
}

## \&code = _setbad_sub($pdlcode)
##  + returns a sub implementing setbadtoval(), setvaltobad(), etc.
sub _setbad_sub {
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
  eval "*${badsub} = _setbad_sub(PDL->can('$badsub'));";
}

##--------------------------------------------------------------
## Dimension Shuffling

## $ccs2 = $ccs->dummy($vdim_index)
## $ccs2 = $ccs->dummy($vdim_index, $vdim_size)
sub dummy {
  my ($ccs,$vdimi,$vdimsize) = @_;
  my @vdims = $ccs->[$VDIMS]->list;
  $vdimsize = 1 if (!defined($vdimsize));
  $vdimi    = 0 if (!defined($vdimi));
  $vdimi    = @vdims + $vdimi + 1 if ($vdimi < 0);
  if ($vdimi < 0) {
    barf(ref($ccs). "::dummy(): negative dimension number ", ($vdimi+@vdims), " exceeds number of dims ", scalar(@vdims));
  }
  splice(@vdims,$vdimi,0,-$vdimsize);
  my $ccs2 = $ccs->copyShallow;
  $ccs2->[$VDIMS] = PDL->pdl($P_LONG,\@vdims);
  return $ccs2;
}

## $ccs2 = $ccs->reorder_pdl($vdim_index_pdl)
sub reorder_pdl {
  my $ccs2 = $_[0]->copyShallow;
  $ccs2->[$VDIMS] = $ccs2->[$VDIMS]->index($_[1]);
  $ccs2->[$VDIMS]->sever;
  $ccs2;
}

## $ccs2 = $ccs->reorder(@vdim_list)
sub reorder { $_[0]->reorder_pdl(PDL->pdl($P_LONG,@_[1..$#_])); }

## $ccs2 = $ccs->xchg($vdim1,$vdim2)
sub xchg {
  my $dimpdl = PDL->sequence($P_LONG,$_[0]->ndims);
  my $tmp    = $dimpdl->at($_[1]);
  $dimpdl->set($_[1], $dimpdl->at($_[2]));
  $dimpdl->set($_[2], $tmp);
  return $_[0]->reorder_pdl($dimpdl);
}

## $ccs2 = $ccs->mv($vDimFrom,$vDimTo)
sub mv  {
  my ($d1,$d2) = @_[1,2];
  my $ndims = $_[0]->ndims;
  $d1 = $ndims+$d1 if ($d1 < 0);
  $d2 = $ndims+$d2 if ($d2 < 0);
  return $_[0]->reorder($d1 < $d2
			? ((0..($d1-1)), (($d1+1)..$d2), $d1,            (($d2+1)..($ndims-1)))
			: ((0..($d2-1)), $d1,            ($d2..($d1-1)), (($d1+1)..($ndims-1)))
		       );
}

## $ccs2 = $ccs->transpose()
##  + always copies
sub transpose {
  if ($_[0]->ndims==1) {
    return $_[0]->dummy(0,1)->copy;
  } else {
    $_[0]->xchg(0,1)->copy;
  }
}



##--------------------------------------------------------------
## PDL API: Indexing

sub slice {
  barf(ref($_[0])."::slice() is not implemented yet (try dummy, dice_axis, indexND, etc.)");
}

## $nzi = $ccs->indexNDi($ndi)
##  + returns Nnz indices for virtual ND-index PDL $ndi
sub indexNDi {
  my ($ccs,$ndi)   = @_;
  ##
  ##-- get physical dims
  my $dims      = $ccs->dimpdl;
  my $whichdimp = ($dims>=0)->which;
  my $pdimi     = $dims->index($whichdimp);
  ##
  $ndi = $ndi->dice_axis(0,$whichdimp)
    if ( $ndi->dim(0)!=$ccs->[$WHICH]->dim(0) || ($pdimi!=PDL->sequence($ccs->[$WHICH]->dim(0)))->any );
  ##
  my $foundi       = $ndi->vsearchvec($ccs->[$WHICH]);
  my $foundi_mask  = ($ndi==$ccs->[$WHICH]->dice_axis(1,$foundi))->andover;
  $foundi_mask->inplace->not;
  $foundi->where($foundi_mask) .= $ccs->[$WHICH]->dim(1);
  return $foundi;
}

## $vals = $ccs->indexND($ndi)
sub indexND { $_[0][$VALS]->index($_[0]->indexNDi($_[1])); }

## $vals = $ccs->index2d($xi,$yi)
sub index2d { $_[0]->indexND($_[1]->cat($_[2])->xchg(0,1)); }

## $vals = $ccs->index($flati)
sub index {
  my ($ccs,$i) = @_;
  my $dummy  = PDL->pdl(0)->slice(join(',', map {"*$_"} $ccs->dims));
  my @coords = $dummy->one2nd($i);
  my $ind = PDL->zeroes($P_LONG,$ccs->ndims,$i->dims);
  $ind->slice("($_),") .= $coords[$_] foreach (0..$#coords);
  return $ccs->indexND($ind);
}

## $ccs2 = $ccs->dice_axis($axis_v, $axisi)
##  + returns a new ccs object, should participate in dataflow
sub dice_axis {
  my ($ccs,$axis_v,$axisi) = @_;
  ##
  ##-- get
  my $ndims = $ccs->ndims;
  $axis_v = $ndims + $axis_v if ($axis_v < 0);
  barf(ref($ccs)."::dice_axis(): axis ".($axis_v<0 ? ($axis_v+$ndims) : $axis_v)." out of range: should be 0<=dim<$ndims")
    if ($axis_v < 0 || $axis_v >= $ndims);
  my $axis  = $ccs->[$VDIMS]->at($axis_v);
  my $asize = $axis < 0 ? -$axis : $ccs->[$PDIMS]->at($axis);
  $axisi    = PDL->topdl($axisi);
  my ($aimin,$aimax) = $axisi->minmax;
  barf(ref($ccs)."::dice_axis(): invalid index $aimin (valid range 0..".($asize-1).")") if ($aimin < 0);
  barf(ref($ccs)."::dice_axis(): invalid index $aimax (valid range 0..".($asize-1).")") if ($aimax >= $asize);
  ##
  ##-- check for virtual
  if ($axis < 0) {
    ##-- we're dicing a virtual axis: ok, but why?
    my $naxisi = $axisi->nelem;
    my $ccs2   = $ccs->copyShallow();
    $ccs2->[$VDIMS] = $ccs->[$VDIMS]->pdl;
    $ccs2->[$VDIMS]->set($axis_v, -$naxisi);
    return $ccs2;
  }
  ##-- ok, we're dicing on a real axis
  my ($ptr,$pi2nzi)    = $ccs->ptr($axis);
  my ($ptrix,$pi2nzix) = $ptr->ccs_decode_pointer($axisi);
  my $nzix   = $pi2nzi->index($pi2nzix);
  my $which  = $ccs->[$WHICH]->dice_axis(1,$nzix);
  $which->sever;
  $which->slice("($axis),") .= $ptrix;
  my $nzvals = $ccs->[$VALS]->index($nzix->append($ccs->[$WHICH]->dim(1)));
  ##
  ##-- construct output object
  my $ccs2 = $ccs->shadow();
  $ccs2->[$PDIMS]->set($axis, $axisi->nelem);
  $ccs2->[$WHICH] = $which;
  $ccs2->[$VALS]  = $nzvals;
  ##
  ##-- sort output object (if not dicing on 0th dimension)
  return $axis==0 ? $ccs2 : $ccs2->sortwhich();
}

## $onedi = $ccs->n2oned($ndi)
##  + returns a pseudo-index
sub n2oned {
  my $dimsizes = PDL->pdl($P_LONG,1)->append($_[0][$VDIMS]->abs)->slice("0:-2")->cumuprodover;
  return ($_[1] * $dimsizes)->sumover;
}

## $whichND = $obj->whichND
##  + just returns the literal index PDL if possible: beware of dataflow!
##  + indices are NOT guaranteed to be returned in any surface-logical order,
##    although physically indexed dimensions should be sorted in physical-lexicographic order
sub whichND {
  my $vpi = ($_[0][$VDIMS]>=0)->which;
  my ($wnd);
  if ( $_[0][$VDIMS]->nelem==$_[0][$PDIMS]->nelem ) {
    if (($vpi==$_[0][$PDIMS]->sequence)->all) { $wnd=$_[0][$WHICH]; }     ##-- all literal & physical
    else { $wnd=$_[0][$WHICH]->dice_axis(0,$_[0][$VDIMS]->index($vpi)); } ##-- all physical, but shuffled
    return wantarray ? $wnd->xchg(0,1)->dog : $wnd;
  }
  ##-- virtual dims are in the game: construct output pdl
  my $ccs = shift;
  my $nvperp = $ccs->_ccs_nvperp;
  my $nv     = $ccs->nstored_v;
  $wnd = PDL->zeroes($P_LONG, $ccs->ndims, $nv);
  $wnd->dice_axis(0,$vpi)->flat .= $ccs->[$WHICH]->slice(",*$nvperp,")->flat;
  my $nzi    = PDL->sequence($P_LONG,$nv);
  my @vdims    = $ccs->[$VDIMS]->list;
  my ($vdimi,);
  foreach (grep {$vdims[$#vdims-$_]<0} (0..$#vdims)) {
    $vdimi = $#vdims-$_;
    $nzi->modulo(-$vdims[$vdimi], $wnd->slice("($vdimi),"), 0);
  }
  return wantarray ? $wnd->xchg(0,1)->dog : $wnd;
}

## $whichVals = $ccs->whichVals()
##  + returns $VALS corresponding to whichND() indices
##  + beware of dataflow!
sub whichVals {
  my $vpi = ($_[0][$VDIMS]>=0)->which;
  return $_[0]->_nzvals() if ( $_[0][$VDIMS]->nelem==$_[0][$PDIMS]->nelem ); ##-- all physical
  ##
  ##-- virtual dims are in the game: construct output pdl
  return $_[0]->_nzvals->slice("*".($_[0]->_ccs_nvperp))->flat;
}

## $which = $obj->which()
##  + not guaranteed to be returned in any meaningful order
sub which { $_[0]->n2oned($_[0]->whichND); }

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
    my $vals    = $ccs->whichVals;
    my $missing = $ccs->missing;
    my @dims    = $ccs->dims;
    my ($which1,$vals1);
    if ($which->dim(0) <= 1) {
      ##-- flat sum
      $which1 = PDL->zeroes($P_LONG,1,$which->dim(1)); ##-- dummy
      $vals1  = $vals;
    } else {
      $which1   = $which->slice("1:-1,");
      my $sorti = $which1->qsortveci;
      $which1   = $which1->dice_axis(1,$sorti);
      $vals1    = $vals->index($sorti);
    }
    ##
    ##-- guts
    my ($which2,$nzvals2) = $accumsub->($which1,$vals1,
					($allow_bad_missing || $missing->isgood ? ($missing,$dims[0]) : (0,0))
				       );
    ##
    ##-- get output pdl
    shift(@dims);
    return $nzvals2->squeeze if (!@dims); ##-- just a scalar: return a plain PDL
    ##
    my $newdims = PDL->pdl($P_LONG,\@dims);
    return $ccs->shadow(
			pdims =>$newdims,
			vdims =>$newdims->sequence,
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

## ($rdpdl,$pdimsc,$vdimsc,$apcp,$bpcp) = _ccsnd_binop_align_dims($pdimsa,$vdimsa, $pdimsb,$vdimsb, $opname)
##  + returns:
##     $rdpdl  : (long,2,$nrdims) : [ [$vdimai,$vdimbi], ...] s.t. $vdimai should align with $vdimbi
##     $pdimsc : (long,$ndimsc)   : physical dim-size pdl for CCS output $c()
##     $vdimsc : (long,$ndimsc)   : virtual dim-size pdl for CCS output $c()
##     $apcp   : (long,2,$nac)    : [ [$apdimi,$cpdimi], ... ] s.t. $cpdimi aligns 1-1 with $apdimi
##     $bpcp   : (long,2,$nbc)    : [ [$bpdimi,$cpdimi], ... ] s.t. $cpdimi aligns 1-1 with $bpdimi
sub _ccsnd_binop_align_dims {
  my ($pdimsa,$vdimsa,$pdimsb,$vdimsb, $opname) = @_;
  $opname = '_ccsnd_binop_relevant_dims' if (!defined($opname));

  ##-- init
  my @pdimsa = $pdimsa->list;
  my @pdimsb = $pdimsb->list;
  my @vdimsa = $vdimsa->list;
  my @vdimsb = $vdimsb->list;

  ##-- get alignment-relevant dims
  my @rdims  = qw();
  my ($vdima,$vdimb, $dimsza,$dimszb);
  foreach (0..($#vdimsa < $#vdimsb ? $#vdimsa : $#vdimsb)) {
    $vdima = $vdimsa[$_];
    $vdimb = $vdimsb[$_];

    ##-- get (virtual) dimension sizes
    $dimsza = $vdima>=0 ? $pdimsa[$vdima] : -$vdima;
    $dimszb = $vdimb>=0 ? $pdimsb[$vdimb] : -$vdimb;

    ##-- check for (virtual) size mismatch
    next if ($dimsza==1 || $dimszb==1);   ##... ignoring (virtual) dims of size 1
    barf( __PACKAGE__ , "::$opname(): dimension size mismatch on dim($_): $dimsza != $dimszb")
      if ($dimsza != $dimszb);

    ##-- dims match: only align if both are physical
    push(@rdims, [$vdima,$vdimb]) if ($vdima>=0 && $vdimb>=0);
  }
  my $rdpdl = PDL->pdl($P_LONG,\@rdims);

  ##-- get output dimension sources
  my @_cdsrc = qw(); ##-- ( $a_or_b_for_dim0, ... )
  foreach (0..($#vdimsa > $#vdimsb ? $#vdimsa : $#vdimsb)) {
    push(@vdimsa, -1) if ($_ >= @vdimsa);
    push(@vdimsb, -1) if ($_ >= @vdimsb);
    $vdima  = $vdimsa[$_];
    $vdimb  = $vdimsb[$_];
    $dimsza = $vdima>=0 ? $pdimsa[$vdima] : -$vdima;
    $dimszb = $vdimb>=0 ? $pdimsb[$vdimb] : -$vdimb;
    if ($vdima>=0) {
      if ($vdimb>=0)  { push(@_cdsrc, $dimsza>=$dimszb ? 0 : 1); } ##-- a:p, b:p --> c:p[max2(sz(a),sz(b)]
      else            { push(@_cdsrc, 0); }                        ##-- a:p, b:v --> c:p[a]
    }
    elsif ($vdimb>=0) { push(@_cdsrc, 1); }                        ##-- a:v, b:p --> c:p[b]
    else              { push(@_cdsrc, $dimsza>=$dimszb ? 0 : 1); } ##-- a:v, b:v --> c:v[max2(sz(a),sz(b))]
  }
  my $_cdsrcp = PDL->pdl($P_LONG,@_cdsrc);

  ##-- get c() dimension pdls
  my @pdimsc = qw();
  my @vdimsc = qw();
  my @apcp  = qw(); ##-- ([$apdimi,$cpdimi], ...)
  my @bpcp  = qw(); ##-- ([$bpdimi,$bpdimi], ...)
  foreach (0..$#_cdsrc) {
    if ($_cdsrc[$_]==0) {
      if ($vdimsa[$_]<0) { $vdimsc[$_]=$vdimsa[$_]; }
      else {
	$vdimsc[$_] = @pdimsc;
	push(@apcp, [$vdimsa[$_],scalar(@pdimsc)]);
	push(@pdimsc, $pdimsa[$vdimsa[$_]]);
      }
    } else {
      if ($vdimsb[$_]<0) { $vdimsc[$_]=$vdimsb[$_]; }
      else {
	$vdimsc[$_] = @pdimsc;
	push(@bpcp, [$vdimsb[$_],scalar(@pdimsc)]);
	push(@pdimsc, $pdimsb [$vdimsb[$_]]);
      }
    }
  }
  my $pdimsc = PDL->pdl($P_LONG,\@pdimsc);
  my $vdimsc = PDL->pdl($P_LONG,\@vdimsc);
  my $apcp   = PDL->pdl($P_LONG,\@apcp);
  my $bpcp   = PDL->pdl($P_LONG,\@bpcp);

  return ($rdpdl,$pdimsc,$vdimsc,$apcp,$bpcp);
}


## \&code = _ccsnd_binary_op_mia($opName, \&pdlSub, $defType)
##  + returns code for wrapping a builtin PDL binary operation \&pdlSub under the name "$opName"
##  + $opName is just used for error reporting
##  + $defType (if specified) is the default output type of the operation (e.g. PDL::long())
sub _ccsnd_binary_op_mia {
  my ($opname,$pdlsub,$deftype) = @_;

  return sub {
    my ($a,$b,$swap) = @_;
    $swap=0 if (!defined($swap));

    ##-- check for & dispatch scalar operations
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

    ##-- align dimensions & determine output sources
    my ($rdpdl,$pdimsc,$vdimsc,$apcp,$bpcp) = _ccsnd_binop_align_dims(@$a[$PDIMS,$VDIMS],
								      @$b[$PDIMS,$VDIMS],
								      $opname);
    my $nrdims = $rdpdl->dim(1);

    ##-- get & sort relevant indices, vals
    my $ixa    = $a->[$WHICH];
    my $avals  = $a->[$VALS];
    my $nixa   = $ixa->dim(1);
    my $ra     = $rdpdl->slice("(0)");
    my ($ixar,$avalsr);
    if ( $rdpdl->isempty ) {
      ##-- a: no relevant dims: align all pairs using a pseudo-dimension
      $ixar   = PDL->zeroes($P_LONG, 1,$nixa);
      $avalsr = $avals;
    } elsif ( ($ra==PDL->sequence($P_LONG,$nrdims))->all ) {
      ##-- a: relevant dims are a prefix of physical dims, e.g. pre-sorted
      $ixar   = $nrdims==$ixa->dim(0) ? $ixa : $ixa->slice("0:".($nrdims-1));
      $avalsr = $avals;
    } else {
      $ixar          = $ixa->dice_axis(0,$ra);
      my $ixar_sorti = $ixar->isempty ? PDL->null->long : $ixar->qsortveci;
      $ixa           = $ixa->dice_axis(1,$ixar_sorti);
      $ixar          = $ixar->dice_axis(1,$ixar_sorti);
      $avalsr        = $avals->index($ixar_sorti);
    }
    ##
    my $ixb   = $b->[$WHICH];
    my $bvals = $b->[$VALS];
    my $nixb  = $ixb->dim(1);
    my $rb    = $rdpdl->slice("(1)");
    my ($ixbr,$bvalsr);
    if ( $rdpdl->isempty ) {
      ##-- b: no relevant dims: align all pairs using a pseudo-dimension
      $ixbr   = PDL->zeroes($P_LONG, 1,$nixb);
      $bvalsr = $bvals;
    } elsif ( ($rb==PDL->sequence($P_LONG,$nrdims))->all ) {
      ##-- b: relevant dims are a prefix of physical dims, e.g. pre-sorted
      $ixbr   = $nrdims==$ixb->dim(0) ? $ixb : $ixb->slice("0:".($nrdims-1));
      $bvalsr = $bvals;
    } else {
      $ixbr          = $ixb->dice_axis(0,$rb);
      my $ixbr_sorti = $ixbr->isempty ? PDL->null->long : $ixbr->qsortveci;
      $ixb           = $ixb->dice_axis(1,$ixbr_sorti);
      $ixbr          = $ixbr->dice_axis(1,$ixbr_sorti);
      $bvalsr        = $bvals->index($ixbr_sorti);
    }


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
    my $ixc    = PDL->zeroes($P_LONG, $pdimsc->nelem, $blksz);
    my $nnzc   = 0;
    my $zc     = $pdlsub->($avals->slice("-1"), $bvals->slice("-1"), $swap)->convert($nzc->type);
    my $nanismissing = ($a->[$FLAGS]&$CCSND_NAN_IS_MISSING);
    my $badismissing = ($a->[$FLAGS]&$CCSND_BAD_IS_MISSING);
    $zc              = $zc->setnantobad() if ($nanismissing && $badismissing);
    my $zc_isbad     = $zc->isbad ? 1 : 0;

    ##-- block-wise variables
    ##   + there are way too many of these...
    my ($nzai_prv,$nzai_pnx, $nzbi_prv,$nzbi_pnx, $nzci_prv,$nzci_pnx,$cmpval_prv);
    my ($nzai_cur,$nzai_nxt, $nzbi_cur,$nzbi_nxt, $nzci_cur,$nzci_nxt,$cmpval);
    my ($nzci_max, $blk_slice, $nnzc_blk,$nnzc_slice_blk);
    my ($nzai_blk,$nzbi_blk,$ixa_blk,$ixb_blk,$ixc_blk,$nzc_blk,$cimask_blk,$ciwhich_blk);
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
	  $ixc_blk = $ixc->slice(",$nnzc_slice_blk");
	  if (!$apcp->isempty) {
	    $ixa_blk = $ixa->dice_axis(1,$nzai_blk->index($ciwhich_blk));
	    $ixc_blk->dice_axis(0,$apcp->slice("(1),")) .= $ixa_blk->dice_axis(0,$apcp->slice("(0),"));
	  }
	  if (!$bpcp->isempty) {
	    $ixb_blk = $ixb->dice_axis(1,$nzbi_blk->index($ciwhich_blk));
	    $ixc_blk->dice_axis(0,$bpcp->slice("(1),")) .= $ixb_blk->dice_axis(0,$bpcp->slice("(0),"));
	  }

	  ##-- construct block output pdls: nzc
	  $nzc->slice($nnzc_slice_blk) .= $nzc_blk;
	}
      }

      ##-- possibly allocate for another block
      if ($nzai_cur < $nixa || $nzbi_cur < $nixb) {
	$nzci_nxt -= $nzci_cur;
	$nzci_cur  = 0;

	if ($nzci_nxt+$blksz > $nzai->dim(0)) {
	  $nzai = $nzai->reshape($nzci_nxt+$blksz);
	  $nzbi = $nzbi->reshape($nzci_nxt+$blksz);
	}
	$ixc = $ixc->reshape($ixc->dim(0), $ixc->dim(1)+$nzai->dim(0));
	$nzc = $nzc->reshape($nzc->dim(0)+$nzai->dim(0));

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
		       pdims => $pdimsc,
		       vdims => $vdimsc,
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

##-- arithmetical & comparison operations
foreach my $binop (
		   qw(plus minus mult divide modulo power),
		   qw(gt ge lt le eq ne spaceship),
		  )
  {
    eval "*${binop} = *${binop}_mia = _ccsnd_binary_op_mia('${binop}',PDL->can('${binop}'));";
  }

##-- integer-only operations
foreach my $intop (
		   qw(and2 or2 xor shiftleft shiftright),
		  )
  {
    eval "*${intop} = *${intop}_mia = _ccsnd_binary_op_mia('${intop}',PDL->can('${intop}'),\$P_LONG);";
  }

## rassgn_mia($to,$from): binary assignment operation with missing-annihilator assumption
##  + argument order is REVERSE of PDL 'assgn()' argument order
*rassgn_mia = _ccsnd_binary_op_mia('rassgn', sub { PDL::assgn($_[1],$_[0]); $_[1]; });

## $to = $to->rassgn($from)
##  + calls newFromDense() with $to flags if $from is dense
##  + otherwise, copies $from to $to
##  + argument order is REVERSED wrt PDL::assgn()
sub rassgn {
  my ($to,$from) = @_;
  if (!ref($from) || $from->nelem==1) {
    ##-- assignment from a scalar: treat the Nd object as a mask of available values
    $to->[$VALS] .= todense($from);
    return $to;
  }
  if (isa($from,__PACKAGE__)) {
    ##-- assignment from a CCS object: copy on a full dim match or an empty "$to"
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
##  + obeys PDL conventions
sub assgn { $_[1]->rassgn($_[0]); }


##--------------------------------------------------------------
## CONTINUE HERE

## TODO:
##  + virtual dimensions: clump
##  + OPERATIONS:
##    - matrix: matmult
##    - accumulators: (some still missing: statistical, extrema-indices, atan2, ...)

##--------------------------------------------------------------
## Matrix operations

## $c = $a->inner($b)
##  + inner product (may produce a large temporary)
sub inner { $_[0]->mult_mia($_[1],0)->sumover; }

## $c = $a->matmult($b)
##  + mostly ganked from PDL::Primitive::matmult
sub matmult {
  barf("Invalid number of arguments for ", __PACKAGE__, "::matmult") if ($#_ < 1);
  my ($a,$b,$c) = @_; ##-- no $c!

  $b=toccs($b); ##-- ensure 2nd arg is a CCS object

  while ($a->getndims < 2) {$a = $a->dummy(-1)} # promote if necessary
  while ($b->getndims < 2) {$b = $b->dummy(-1)}

  if ( ($a->dim(0)==1 && $a->dim(1)==1) || ($b->dim(0)==1 && $b->dim(1)==1) ) {
    if (defined($c)) { @$c = @{$a*$b}; return $c; }
    return $a*$b;
  }

  if ($b->dim(1) != $a->dim(0)) {
    barf(sprintf("Dim mismatch in ", __PACKAGE__ , "::matmult of [%dx%d] x [%dx%d]: %d != %d",
		 $a->dim(0),$a->dim(1),$b->dim(0),$b->dim(1),$a->dim(0),$b->dim(1)));
  }

  my $_c = $a->dummy(1)->inner($b->xchg(0,1)->dummy(2)); ##-- ye olde guttes
  if (defined($c)) { @$c = @$_c; return $c; }

  return $_c;
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
		 + PDL->pdl($ccs->[$VDIMS]->nelem) * PDL::howbig($ccs->[$VDIMS]->type)
		);
  return (($dsize - $ccssize) / $dsize)->sclr;
}

##--------------------------------------------------------------
## Stringification & Viewing

## $dimstr = _dimstr($pdl)
sub _dimstr { return '('.$_[0]->type.', '.join(',',$_[0]->dims).')'; }
sub _pdlstr { return _dimstr($_[0]).'='.$_[0]; }

## $str = $obj->string()
sub string {
  my ($pdims,$vdims,$which,$vals) = @{$_[0]}[$PDIMS,$VDIMS,$WHICH,$VALS];
  my $whichstr  = ''.$which->xchg(0,1);
  $whichstr =~ s/^([^A-Z])/$,  $1/mg;
  chomp($whichstr);
  return
    (
     ''
     .ref($_[0]) . _dimstr($_[0]) ."\n"
     .$,." pdims:" . _pdlstr($pdims) ."\n"
     .$,." vdims:" . _pdlstr($vdims) ."\n"
     .$,." which:" . _dimstr($which)."^T=" . $whichstr . "\n"
     .$,." vals:" . _pdlstr($vals)  ."\n"
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

	      ##-- matrix operations
	      'x' => \&matmult,

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

