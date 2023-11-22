/*============================================================================
 * User definition of physical properties.
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
#include <string.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_physical_properties.c
 *
 * \brief User definition of physical properties.
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Function called at each time step to define physical properties.
 *
 * \param[in, out]  domain   pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_physical_properties(cs_domain_t   *domain)
{
  CS_NO_WARN_IF_UNUSED(domain);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief User definition of enthalpy to temperature conversion.
 *
 * This allows overwriting the solver defaults if necessary.
 *
 * This function may be called on a per-zone basis, so as to allow different
 * conversion relations in zones representing solids or different fluids.
 *
 * \param[in, out]  domain   pointer to a cs_domain_t structure
 * \param[in]       z        zone (volume or boundary) applying to current call
 * \param[in]       z_local  if true, h and t arrays are defined in a compact
 *                           (contiguous) manner for this zone only;
 *                           if false, h and t are defined on the zone's parent
 *                           location (usually all cells or boundary faces)
 * \param[in]       h        enthalpy values
 * \param[in, out]  t        temperature values
 */
/*----------------------------------------------------------------------------*/

void
cs_user_physical_properties_h_to_t(cs_domain_t      *domain,
                                   const cs_zone_t  *z,
                                   bool              z_local,
                                   const cs_real_t   h[],
                                   cs_real_t         t[])
{
  CS_NO_WARN_IF_UNUSED(domain);
  CS_NO_WARN_IF_UNUSED(z);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief User definition of temperature to enthalpy conversion.
 *
 * This allows overwriting the solver defaults if necessary.
 *
 * This function may be called on a per-zone basis, so as to allow different
 * conversion relations in zones representing solids or different fluids.
 *
 * \param[in, out]  domain   pointer to a cs_domain_t structure
 * \param[in]       z        zone (volume or boundary) applying to current call
 * \param[in]       z_local  if true, h and t arrays are defined in a compact
 *                           (contiguous) manner for this zone only;
 *                           if false, h and t are defined on the zone's parent
 *                           location (usually all cells or boundary faces)
 * \param[in]       h        temperature values
 * \param[in, out]  t        enthalpy values
 */
/*----------------------------------------------------------------------------*/

void
cs_user_physical_properties_t_to_h(cs_domain_t      *domain,
                                   const cs_zone_t  *z,
                                   bool              z_local,
                                   const cs_real_t   t[],
                                   cs_real_t         h[])
{
  CS_NO_WARN_IF_UNUSED(domain);
  CS_NO_WARN_IF_UNUSED(z);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief User modification of the Smagorinsky constant for the
 *        dynamic Smagorinsky model.
 *
 * CS = Mij.Lij / Mij.Mij
 *
 * The local averages of the numerator and denominator are done before calling
 * this function, so
 *
 * CS = < Mij.Lij > / < Mij.Mij >
 *
 * In this subroutine, Mij.Lij and Mij.Mij are passed as arguments before
 * the local average.
 *
 * \param[in, out]   domain      pointer to a cs_domain_t structure
 * \param[in]        mijlij      mij.lij before the local averaging
 * \param[in]        mijmij      mij.mij before the local averaging
 */
/*----------------------------------------------------------------------------*/

void
cs_user_physical_properties_smagorinsky_c(cs_domain_t      *domain,
                                          const cs_real_t   mijlij[],
                                          const cs_real_t   mijmij[])
{
  CS_NO_WARN_IF_UNUSED(domain);
}


/*----------------------------------------------------------------------------*/
/*!
 * \brief User function to define a custom law for the thermodynamic pressure.
 *
 * Allows to define a custom law for the constant uniform thermodynamic
 * pressure (whenn \ref cs_velocity_pressure_model_t::idilat = 3 or
 * \ref cs_fluid_properties_t::ipthrm = 1).
 *
 * The density is then updated (in \ref pthrbm.f90) as:
 * \f[\rho^{n+1} =\rho^{n} \cdot \frac{P_{th}^{n+1}}{P_{th}^{n}}\f].
 *
 * \section Usage
 *
 * Here is an extremely basic example where the thermodynamic pressure
 * increases by 10.0 at every iteration.
 * \code{.c}
 * const cs_real_t delta_p = 10.0;
 *
 * td_p = cs_glob_physical_properties->pther + delta_p;
 * \endcode
 *
 * \param[in, out]  td_p  Updated value of the thermodynamic pressure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_physical_properties_td_pressure(cs_real_t  *td_p)
{
  CS_NO_WARN_IF_UNUSED(td_p);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
