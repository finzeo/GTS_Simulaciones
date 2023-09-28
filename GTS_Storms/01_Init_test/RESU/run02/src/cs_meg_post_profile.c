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
  /* User-defined coordinates for profile checkVxProfile
    ---------------------------------------------------- */

  if (strcmp(name, "checkVxProfile") == 0) {
    cs_real_t w = -1.;
    cs_real_t xo = -1.;
    cs_real_t x1 = -1.;
    cs_real_t y0 = -1.;
    cs_real_t y1 = -1.;
    cs_real_t z0 = -1.;
    cs_real_t z1 = -1.;
    cs_real_t x, y, z;

    for (int p_id = 0; p_id < n_coords; p_id++) {
      cs_real_t s = (cs_real_t)p_id / (cs_real_t)(n_coords-1);

    //example: a line segment between (x0, y0, z0) and (x1, y1, z1)
    //(s, the parameter is always between 0 and 1)
    w = 0.32385;
    xo = 0.41 * w;
    x1 = 0.41 * w;
    y0 = 0;
    y1 = 0;
    z0 = - 0.0762;
    z1 = 0.0762;
    x = x1 * s + x0 * (1. - s);
    y = y1 * s + y0 * (1. - s);
    z = z1 * s + z0 * (1. - s);

      coords[p_id][0] = x;
      coords[p_id][1] = y;
      coords[p_id][2] = z;
    }
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS

