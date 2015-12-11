## Automatically generated, remove to re-configure!

package PDL::CCS::Config;
use PDL qw();
our @ISA = qw(Exporter);
our (%ccsConfig);
our @EXPORT      = qw(ccs_indx);
our @EXPORT_OK   = ('%ccsConfig', 'ccs_indx');
our %EXPORT_TAGS = (config=>['%ccsConfig'], Func=>\@Export, default=>\@EXPORT, all=>\@EXPORT_OK);

%ccsConfig = (
               'INDX_FUNCDEF' => '*ccs_indx = \\&PDL::indx; ##-- typecasting for CCS indices
',
               'INDX_SIG' => 'indx',
               'INDX_TYPEDEF' => 'typedef PDL_Indx CCS_Indx;  /**< typedef for CCS indices */
',
               'INT_TYPE_KEYS' => [
                                    'PDL_B',
                                    'PDL_IND',
                                    'PDL_L',
                                    'PDL_LL',
                                    'PDL_S',
                                    'PDL_US'
                                  ],
               'INDX_CTYPE' => 'PDL_Indx',
               'INDX_FUNC' => 'indx',
               'INT_TYPE_CHRS' => [
                                    'B',
                                    'L',
                                    'N',
                                    'Q',
                                    'S',
                                    'U'
                                  ],
               'USE_PDL_INDX' => 1
             );

*PDL::ccs_indx = *ccs_indx = \&PDL::indx; ##-- typecasting for CCS indices


1; ##-- be happy
