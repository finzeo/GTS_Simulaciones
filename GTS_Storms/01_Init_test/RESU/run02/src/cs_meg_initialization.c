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

cs_real_t *
cs_meg_initialization(const char      *zone_name,
                      const cs_lnum_t  n_elts,
                      const cs_lnum_t *elt_ids,
                      const cs_real_t  xyz[][3],
                      const char      *field_name)
{
  cs_real_t *new_vals = NULL;

  /* User defined initialization for variable velocity over zone all_cells
    ---------------------------------------------------------------------- */

  if (strcmp(zone_name, "all_cells") == 0 &&
      strcmp(field_name, "velocity") == 0) {

    const cs_lnum_t vals_size = n_elts * 3;
    BFT_MALLOC(new_vals, vals_size, cs_real_t);

    const cs_real_t vel_inf = cs_notebook_parameter_value_by_name("vel_inf");

    for (cs_lnum_t e_id = 0; e_id < n_elts; e_id++) {
      cs_lnum_t c_id = elt_ids[e_id];
      new_vals[3 * e_id + 0] = vel_inf;
      new_vals[3 * e_id + 1] = 0.;
      new_vals[3 * e_id + 2] = 0.;
    }
  }

  /* User defined initialization for variable turbulence over zone all_cells
    ------------------------------------------------------------------------ */

  if (strcmp(zone_name, "all_cells") == 0 &&
      strcmp(field_name, "turbulence") == 0) {

    const cs_lnum_t vals_size = n_elts * 2;
    BFT_MALLOC(new_vals, vals_size, cs_real_t);

    const cs_real_t vel_inf = cs_notebook_parameter_value_by_name("vel_inf");

    cs_real_t I = -1.;
    cs_real_t dh = -1.;
    cs_real_t l = -1.;
    for (cs_lnum_t e_id = 0; e_id < n_elts; e_id++) {
      cs_lnum_t c_id = elt_ids[e_id];
      I = 0.05;
      // Diametro hidraulico
      dh = 2.510117647;
      // Turbulent length scale (l)
      l = 0.07 * dh;
      new_vals[2 * e_id + 0] = 1.5 * cs_math_pow2((I * vel_inf));
      new_vals[2 * e_id + 1] = (pow(new_vals[2 * e_id + 0], 0.5)) / l;
    }
  }
  return new_vals;
}

/*----------------------------------------------------------------------------*/

END_C_DECLS

