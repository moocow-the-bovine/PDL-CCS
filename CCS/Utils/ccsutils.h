#include "pdl.h"

/* null-detection adapted from PDL_MAYBE_SIZE macro:
 *   https://github.com/PDLPorters/pdl-linearalgebra/blob/f789c4100d04ba9d1b50f8c18249bdef29338496/Real/real.pd#L63-L75
 */
#define CCS_PDL_IS_NULL(pdl) \
  ((pdl)->nvals==0 && ((pdl)->state & PDL_MYDIMS_TRANS))
