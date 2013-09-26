## Automatically generated, remove to re-configure!

package PDL::CCS::Config;
use PDL qw();
our @ISA = qw(Exporter);
our (%ccsConfig);
our @EXPORT      = qw(ccs_indx);
our @EXPORT_OK   = ('%ccsConfig', 'ccs_indx');
our %EXPORT_TAGS = (config=>['%ccsConfig'], Func=>\@Export, default=>\@EXPORT, all=>\@EXPORT_OK);

%ccsConfig = (
               'INDX_SIG' => 'int',
               'INT_TYPE_KEYS' => [
                                    'PDL_B',
                                    'PDL_L',
                                    'PDL_LL',
                                    'PDL_S',
                                    'PDL_US'
                                  ],
               'INDX_FUNCDEF' => '*ccs_indx = \\&PDL::long; ##-- typecasting for CCS indices
',
               'INT_TYPE_CHRS' => [
                                    'B',
                                    'L',
                                    'Q',
                                    'S',
                                    'U'
                                  ],
               'INDX_FUNC' => 'long',
               'INDX_TYPEDEF' => 'typedef PDL_Long CCS_Indx;  /**< typedef for CCS indices */
',
               'USE_PDL_INDX' => '',
               'INDX_CTYPE' => 'PDL_Long'
             );

*PDL::ccs_indx = *ccs_indx = \&PDL::long; ##-- typecasting for CCS indices


1; ##-- be happy
