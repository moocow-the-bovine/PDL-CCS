## File: PDL::CCS::IO::LDAC.pm
## Author: Bryan Jurish <moocow@cpan.org>
## Description: LDA-C wrappers for PDL::CCS::Nd

package PDL::CCS::IO::LDAC;
use PDL::CCS::Version;
use PDL::CCS::Config qw(ccs_indx);
use PDL::CCS::Nd;
use PDL::CCS::IO::FastRaw qw(:intern); ##-- for _ccsio_generate_header(), _ccsio_parse_header()
use PDL;
use PDL::IO::Misc;
use Carp qw(confess);
use strict;

our $VERSION = '1.22.6';
our @ISA = ('PDL::Exporter');
our @EXPORT_OK =
  (
   qw(ccs_writeldac ccs_readldac),
  );
our %EXPORT_TAGS =
  (
   Func => [@EXPORT_OK],               ##-- respect PDL conventions (hopefully)
  );

##======================================================================
## pod: headers
=pod

=head1 NAME

PDL::CCS::IO::LDAC - LDA-C format text I/O for PDL::CCS::Nd

=head1 SYNOPSIS

 use PDL;
 use PDL::CCS::Nd;
 use PDL::CCS::IO::LDAC;

 $tdm = PDL::CCS::Nd->newFromWhich($which,$nnz);

 ccs_writeldac($tdm,"tdm.ldac");   # write a sparse matrix market text file
 $tdm2 = ccs_readldac("tdm.ldac"); # read a sparse matrix market text file

=cut


##======================================================================
## I/O utilities
=pod

=head1 I/O Utilities

=cut

##---------------------------------------------------------------
## ccs_writeldac
=pod

=head2 ccs_writeldac

Write a 2d L<PDL::CCS::Nd|PDL::CCS::Nd> (Term x Document) matrix as an LDA-C text file.

 ccs_writeldac($ccs,$filename_or_fh)
 ccs_writeldac($ccs,$filename_or_fh,\%opts)

Options %opts: (none)

=cut

*PDL::ccs_writeldac = *PDL::CCS::Nd::writeldac = \&ccs_writeldac;
sub ccs_writeldac {
  my ($ccs,$file,$opts) = @_;
  $opts = {} if (!defined($opts));
  #$opts->{start} = 0 if (!defined($opts->{start}));

  ##-- sanity check(s)
  confess("ccs_writeldac(): input matrix must be 2d!") if ($ccs->ndims != 2);

  ##-- open output file
  my $fh = _ccsio_open($file,'>')
    or confess("ccs_writeldac(): open failed for output file '$file': $!");
  #binmode($fh,':raw');
  local $,='';

  ##-- convert to lda-c format: use ptr()
  my ($ptr,$pi2nzi) = $ccs->ptr(1);
  my $nd = $ptr->nelem-1;
  my $ix = $ccs->_whichND;
  my $nz = $ccs->_nzvals;
  my ($di,$xi,$xj,$nzi);
  for ($di=0; $di < $ptr->nelem; ++$di) {
    ($xi,$xj) = ($ptr->at($di),$ptr->at($di+1));
    print $fh $xj-$xi;
    for ( ; $xi < $xj; ++$xi) {
      $nzi = $pi2nzi->index($xi)->sclr;
      print $fh ' ', $ix->at(0,$nzi), ":", $nz->at($nzi);
    }
    print $fh "\n";
  }

  ##-- cleanup
  _ccsio_close($file,$fh)
    or confess("ccs_writeldac(): close failed for output file '$file': $!");

  return 1;
}


##---------------------------------------------------------------
## ccs_readldac
=pod

=head2 ccs_readldac

Read a 2d L<PDL::CCS::Nd|PDL::CCS::Nd> (Term x Document) matrix from an LDA-C text file.

 $ccs = ccs_readldac($filename_or_fh)
 $ccs = ccs_readldac($filename_or_fh,\%opts)

Options %opts:

 type => $type,      ##-- value datatype; default = $PDL::IO::Misc::deftype

=cut

*PDL::ccs_readldac = *PDL::CCS::Nd::readldac = \&ccs_readldac;
sub ccs_readldac {
  shift if (UNIVERSAL::isa($_[0],'PDL') || UNIVERSAL::isa($_[0],'PDL::CCS::Nd'));
  my ($file,$opts) = @_;
  $opts = {} if (!defined($opts));
  #$opts->{start} = 1 if (!defined($opts->{start}));
  $opts->{type} = $PDL::IO::Misc::deftype if (!defined($opts->{type}));
  $opts->{type} = PDL->can($opts->{type})->() if (!ref($opts->{type}));

  ##-- open input file
  my $fh = _ccsio_open($file,'<')
    or confess("ccs_readldac(): open failed for input file '$file': $!");

  ##-- get nnz (per doc)
  my $d_nnz = PDL->rcols($fh, [0], { TYPES=>[ccs_indx] });
  my $nnz   = $d_nnz->sum;
  undef($d_nnz);

  ##-- allocate output pdls
  my $nd    = $d_nnz->nelem;
  my $ix    = zeroes(ccs_indx, 2,$nnz);
  my $nz    = zeroes($opts->{type}, $nnz+1);

  ##-- process input
  my ($nzi,$di,$ti,$f);
  for ($nzi=$di=0; $di < $nd && $nzi < $nnz && defined($_=<$fh>); ++$di) {
    chomp;
    next if (/^\s*$/ || /^\s*\D/);
    while (/\b([0-9]+):(\S+)\b/g) {
      ($ti,$f) = ($1,$2);
      $ix->set(0,$nzi => $ti);
      $ix->set(1,$nzi => $di);
      $nz->set($nzi => $f);
      ++$nzi;
    }
    ++$di;
  }
  my $nt = $ix->slice("(0),")->max+1;

  ##-- cleanup
  _ccsio_close($file,$fh)
    or confess("ccs_readldac(): close failed for input file '$file': $!");

  ##-- construct and return
  return PDL::CCS::Nd->newFromWhich($ix,$nz,
				    pdims=>[$nt,$nd],
				    sorted=>0,
				    steal=>1,
				   );
}


1; ##-- be happy

##======================================================================
## POD: footer
=pod

=head1 ACKNOWLEDGEMENTS

Perl by Larry Wall.

PDL by Karl Glazebrook, Tuomas J. Lukka, Christian Soeller, and others.

LDA-C package by by David M. Blei.

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
PDL::IO::Misc(3perl),
PDL::CCS::Nd(3perl),
the LDA-C package documentation at L<http://www.cs.princeton.edu/~blei/lda-c/>
...

=cut


1; ##-- make perl happy
