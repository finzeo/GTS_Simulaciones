/*============================================================================
 * This function is called each time step to define physical properties
 *============================================================================*/

/* code_saturne version 8.0 */

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2023 EDF S.A.

  This program is free software; you can redistribute it and/or modify it under
  the terms of the GNU General Public License as published by the Free Software
  Foundation; either version 2 of the License, or (at your option) any later
  version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
  details.

  You should have received a copy of the GNU General Public License along with
  this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
  Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/*----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <math.h>

/*----------------------------------------------------------------------------
 * PLE library headers
 *----------------------------------------------------------------------------*/

#include <ple_coupling.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_source_terms-momentum.c
 *
 * \brief Base examples for additional right-hand side source terms for
 *   momentum equations.
 *
 * See the reference \ref cs_user_source_terms.c for documentation.
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Function called at each time step to define source terms.
 *
 * \param[in, out]  domain   pointer to a cs_domain_t structure
 * \param[in]       f_id     field id of the variable
 * \param[out]      st_exp   explicit source term
 * \param[out]      st_imp   implicit part of the source term
 */
/*----------------------------------------------------------------------------*/

void
cs_user_source_terms(cs_domain_t  *domain,
                     int           f_id,
                     cs_real_t    *st_exp,
                     cs_real_t    *st_imp)
{
  /* field structure */
  const cs_field_t  *f = cs_field_by_id(f_id);

  /* mesh quantities */
  const cs_lnum_t  n_cells = domain->mesh->n_cells;
  const cs_real_t  *cell_f_vol = domain->mesh_quantities->cell_vol;

  /* density */
  const cs_real_t  *cpro_rom = CS_F_(rho)->val;

  /* Example of arbitrary source term for component u:
   *
   *                       S = A * u + B
   *
   *        appearing in the equation under the form
   *
   *                  rho*du/dt = S (+ standard Navier-Stokes terms)
   *
   * In the following example:
   *   A = - rho * CKP
   *   B =   XMMT
   *
   * with:
   *  CKP = 1.0   [1/s]      (return term on velocity)
   *  MMT = 100.0 [kg/m2/s2] (momentum production by volume and time unit)
   *
   * which yields:
   *  st_imp[i][0][0] = cell_f_vol[i] * A = - cell_f_vol[i]*(rho*CKP)
   *  st_exp[i][0]    = cell_f_vol[i] * B =   cell_f_vol[i]*(XMMT)
   */

  if (f == CS_F_(vel)) { /* velocity */

    /* cast to 3D vectors for readability */
    cs_real_3_t    *_st_exp = (cs_real_3_t *)st_exp;
    cs_real_33_t   *_st_imp = (cs_real_33_t *)st_imp;

    /* Density */

    //const cs_real_t ckp = 0.0;
    //const cs_real_t xmmt = 100.0;

    for (cs_lnum_t i = 0; i < n_cells; i++) {
      //_st_imp[i][0][0] = - cell_f_vol[i] * cpro_rom[i] * ckp;
      //_st_exp[i][0]    =   cell_f_vol[i] * cpro_rom[i] * xmmt;
		_st_imp[i][0][0] = 0;
		_st_exp[i][0]    = cell_f_vol[i] * 16.09116261;
    }
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
