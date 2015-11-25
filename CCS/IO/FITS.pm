## File: PDL::CCS::IO::FITS.pm
## Author: Bryan Jurish <moocow@cpan.org>
## Description: PDL::IO::FITS wrappers for PDL::CCS::Nd

package PDL::CCS::IO::FITS;
use PDL::CCS::Version;
use PDL::CCS::Nd;
use PDL;
use PDL::CCS::IO::FastRaw qw(:intern); ##-- for _ccsio_write_header, _ccsio_read_header
use Carp qw(confess);
use strict;

our $VERSION = '1.22.6';
our @ISA = ('PDL::Exporter');
our @EXPORT_OK =
  (
   qw(ccs_wfits ccs_rfits),
  );
our %EXPORT_TAGS =
  (
   Func => [@EXPORT_OK],               ##-- respect PDL conventions (hopefully)
  );



##======================================================================
## pod: headers
=pod

=head1 NAME

PDL::CCS::IO::FITS - PDL::IO::FITS wrappers for PDL::CCS::Nd

=head1 SYNOPSIS

 use PDL;
 use PDL::CCS::Nd;
 use PDL::CCS::IO::FITS;


 $ccs = PDL::CCS::Nd->newFromWhich($which,$nnz);

 ccs_wfits($ccs,$fname);	 # write a pair of FITS files
 $ccs2 = ccs_readfits($fname);   # read a pair of FITS files

=cut


##======================================================================
## I/O utilities
=pod

=head1 I/O Utilities

=cut

##---------------------------------------------------------------
## ccs_wfits
=pod

=head2 ccs_wfits

Write a pair of FITS files using L<PDL::IO::FITS::wfits()|PDL::IO::FITS/wfits>.
Piddles of type L<indx|PDL::Core/indx> will be implicitly converted
to L<long|PDL::Core/long>, since they are not currently supported by L<PDL::IO::FITS|PDL::IO::FITS> in PDL v2.014.

 ccs_wfits($ccs,$fname)
 ccs_wfits($ccs,$fname,\%opts)

Options %opts:

 Header      => $Header,       ##-- default="$fname.hdr"
 ixFile      => $ixFile,       ##-- default="$fname.ix.fits"
 nzFile      => $nzFile,       ##-- default="$fname.nz.fits"

=cut

*PDL::ccs_wfits = *PDL::CCS::Nd::wfits = \&ccs_wfits;
sub ccs_wfits {
  my ($ccs,$fname,$opts) = @_;

  ##-- get filenames
  my $hFile  = $opts->{Header} // "$fname.hdr";
  my $ixFile = $opts->{ixFile} // "$fname.ix.fits";
  my $nzFile = $opts->{nzFile} // "$fname.nz.fits";

  ##-- write header
  _ccsio_write_header({magic=>(ref($ccs)." $VERSION"), pdims=>$ccs->pdims, vdims=>$ccs->vdims, flags=>$ccs->flags}, $hFile)
    or confess("ccs_wfits(): failed to write header-file $hFile: $!");

  ##-- write pdls
  my $ix = $ccs->_whichND;
  $ix    = $ix->long if ($ix->type->ioname eq 'indx'); ##-- hack: treat 'indx' as 'long' until PDL::IO::FITS supports it (PDL v2.014)
  PDL::wfits($ccs->_whichND, $ixFile)
      or confess("ccs_wfits(): failed to write index-file $ixFile: $!");
  PDL::wfits($ccs->_vals,  $nzFile)
      or confess("ccs_wfits(): failed to write values-file $nzFile: $!");

  return 1;
}


##---------------------------------------------------------------
## ccs_rfits
=pod

=head2 ccs_rfits

Read a pair of FITS files using L<PDL::IO::FITS::rfits()|PDL::IO::FITS/rfits()>.

 $ccs = ccs_rfits($fname)
 $ccs = ccs_rfits($fname,\%opts)

Options %opts:

 Header   => $Header,      ##-- default="$fname.hdr"
 ixFile   => $ixFile,      ##-- default="$fname.ix.fits"
 nzFile   => $nzFile,      ##-- default="$fname.nz.fits"
 sorted   => $bool,        ##-- is data on disk already sorted? (default=1)

=cut

*PDL::ccs_rfits = *PDL::CCS::Nd::rfits = \&ccs_rfits;
sub ccs_rfits {
  shift if (UNIVERSAL::isa($_[0],'PDL') || UNIVERSAL::isa($_[0],'PDL::CCS::Nd'));
  my ($fname,$opts) = @_;

  ##-- get filenames
  my $hFile  = $opts->{Header} // "$fname.hdr";
  my $ixFile = $opts->{ixFile} // "$fname.ix.fits";
  my $nzFile = $opts->{nzFile} // "$fname.nz.fits";

  ##-- read header
  my $header = _ccsio_read_header($hFile)
    or confess("ccs_rfits(): failed to read header-file $hFile: $!");

  ##-- read pdls
  defined(my $ix = PDL->rfits($ixFile))
    or confess("ccs_rfits(): failed to read index-file $ixFile: $!");
  defined(my $nz = PDL->rfits($nzFile))
    or confess("ccs_rfits(): failed to read values-file $nzFile: $!");

  ##-- construct and return
  return PDL::CCS::Nd->newFromWhich($ix,$nz,
				    pdims=>$header->{pdims},
				    vdims=>$header->{vdims},
				    flags=>$header->{flags},
				    sorted=>($opts->{sorted}//1),
				    steal=>1);
}



1; ##-- be happy

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

Copyright (C) 2015, Bryan Jurish. All rights reserved.

This package is free software, and entirely without warranty.
You may redistribute it and/or modify it under the same terms
as Perl itself.

=head1 SEE ALSO

perl(1),
PDL(3perl),
PDL::IO::FITS(3perl),
PDL::CCS::Nd(3perl),
...

=cut


1; ##-- make perl happy
