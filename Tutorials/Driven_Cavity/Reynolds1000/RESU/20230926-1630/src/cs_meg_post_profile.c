/*----------------------------------------------------------------------------*/
/*
  This file is generated by code_saturne, a general-purpose CFD tool.
*/
/*----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <math.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/

void
cs_meg_post_profiles(const char       *name,
                     int               n_coords,
                     cs_real_t         coords[][3])
{
  /* User-defined coordinates for profile Xvel_YaxisCenterLine
    ---------------------------------------------------------- */

  if (strcmp(name, "Xvel_YaxisCenterLine") == 0) {
    cs_real_t x, y, z;

    for (int p_id = 0; p_id < n_coords; p_id++) {
      cs_real_t s = (cs_real_t)p_id / (cs_real_t)(n_coords-1);

    x = 0.5;
    y = s;
    z = 0;

      coords[p_id][0] = x;
      coords[p_id][1] = y;
      coords[p_id][2] = z;
    }
  }

  /* User-defined coordinates for profile Yvel_XaxisCenterLine
    ---------------------------------------------------------- */

  if (strcmp(name, "Yvel_XaxisCenterLine") == 0) {
    cs_real_t x, y, z;

    for (int p_id = 0; p_id < n_coords; p_id++) {
      cs_real_t s = (cs_real_t)p_id / (cs_real_t)(n_coords-1);

    x = s;
    y = 0.5;
    z = 0;

      coords[p_id][0] = x;
      coords[p_id][1] = y;
      coords[p_id][2] = z;
    }
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS

