/*============================================================================
 * General-purpose user-defined functions called before time stepping, at
 * the end of each time step, and after time-stepping.
 *
 * These can be used for operations which do not fit naturally in any other
 * dedicated user function.
 *============================================================================*/

/* code_saturne version 8.0-beta */

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
#include <stdio.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * PLE library headers
 *----------------------------------------------------------------------------*/

#include <ple_coupling.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_prototypes.h"
#include <limits.h>
#include <cs_physical_constants.h> // For fluid properties
#include <cs_turbulence_model.h> // For turbulence properties

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_extra_operations-boundary_forces.c
 *
 * \brief This function is called at the end of each time step, and has a very
 * general purpose (i.e. anything that does not have another dedicated
 * user function).
 *
 * This is an example of cs_user_extra_operations.c which computes the total
 * force on a boundary zone.
 *
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief This function is called at the end of each time step.
 *
 * It has a very general purpose, although it is recommended to handle
 * mainly postprocessing or data-extraction type operations.
 *
 * \param[in, out]  domain   pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_extra_operations(cs_domain_t     *domain)
{
  
    // Net force calculation
    // Number of boundary faces
    const cs_lnum_t n_b_faces = domain->mesh->n_b_faces;

    // Return a pointer to a field based on its name if present.
    // If no field of the given name is defined, NULL is returned.
    // All forces of all boundary faces, in an array
    cs_field_t *b_forces = cs_field_by_name_try("boundary_forces");

    //if (b_forces != NULL) {
      // Initilization of total force vector (total_b_forces)
      cs_real_3_t total_b_forces = {0., 0., 0.};
      
      // Number and array of boundary faces (to use only certain parts,
      // containing boundary faces that meet a certain criteria)
      cs_lnum_t n_elts, *lst_elts;
      // Allocate memory for n_b_faces elements of type cs_lnum_t,
      // in lst_elts pointer (to fill a list of boundary faces)
      BFT_MALLOC(lst_elts, n_b_faces, cs_lnum_t);
      // Fill a list of boundary faces verifying a given selection criteria
      // For this case, using a group name of mesh entities, only faces on GTS
      // (name set at snappyHexMeshDict)
      // Creo que se pueden usar operadores lógicos
      cs_selector_get_b_face_list("gtsModel", &n_elts, lst_elts);
      // Another example:
      // cs_selector_get_b_face_list 
      // ("normal[0, -1, 0, 0.1] and box[-1000, -1000, -1000, 1000, 0.01, 1000]",
      // &nlelt, lstelt);

      // Loop over GTS boundary faces
      for (cs_lnum_t i_elt = 0; i_elt < n_elts; i_elt++) {
        cs_lnum_t face_id = lst_elts[i_elt];
	// Loop over dimensions
        for (int ii = 0; ii < 3; ii++)
	  // Sum of forces
          total_b_forces[ii] += b_forces->val[face_id*3 + ii];
      }
      // Free resources
      BFT_FREE(lst_elts);

      /* parallel sum */
      cs_parall_sum(3, CS_DOUBLE, total_b_forces);
    //}
  
    // Force coefficients (drag, lift, side force)
    
    // Get user parameters defined in the GUI notebook (dynamic pressure and yaw angle)
    cs_real_t q_inf = cs_notebook_parameter_value_by_name("pressure_dynamic");
    cs_real_t yaw_angle = cs_notebook_parameter_value_by_name("yaw_angle");
    // Projected area of the GTS in the movement direction
    // (to be modified later according to the yaw angle)
    cs_real_t proj_area = -1.0;
    // Tolerance used to match floating-point values written in notebook environment
    double tol_read=1e-6;
    // WARNING: in case of changing this array, modify accordingly the 
    // structure that takes the value of the projected area (proj_area)
    cs_real_t yaw_angle_values[]={0.0,2.5,5.0,7.5,10.0,12.5,14.0};
    int n_yaw_angle_values = sizeof(yaw_angle_values)/sizeof(yaw_angle_values[0]);
    for (int i=0;i<n_yaw_angle_values;i++){
      if( (yaw_angle-yaw_angle_values[i])<tol_read && (yaw_angle-yaw_angle_values[i])>(-tol_read) ){
	switch(i){
	  case 0: proj_area = cs_notebook_parameter_value_by_name("proj_area_0_0");break;
	  case 1: proj_area = cs_notebook_parameter_value_by_name("proj_area_2_5");break;
	  case 2: proj_area = cs_notebook_parameter_value_by_name("proj_area_5_0");break;
	  case 3: proj_area = cs_notebook_parameter_value_by_name("proj_area_7_5");break;
	  case 4: proj_area = cs_notebook_parameter_value_by_name("proj_area_10_0");break;
	  case 5: proj_area = cs_notebook_parameter_value_by_name("proj_area_12_5");break;
	  case 6: proj_area = cs_notebook_parameter_value_by_name("proj_area_14_0");break;
	  default: {fprintf(stderr, "Projected area not found. Exiting...\n"); exit(-1);}
	}
	break;
      }
    }
    if( (proj_area-(-1.0))<tol_read && (proj_area-(-1.0))>(-tol_read) )
      {fprintf(stderr, "Projected area not found. Exiting...\n"); exit(-1);}
      
    // Force coeficients (WARNING: referred to wind tunnel axis, NOT body axis) 
    cs_real_t c_d = total_b_forces[0]/(q_inf*proj_area);
    cs_real_t c_s = total_b_forces[1]/(q_inf*proj_area);
    cs_real_t c_l = total_b_forces[2]/(q_inf*proj_area);
  
    // Net static pressure force calculation
    // Number of boundary faces
    //const cs_lnum_t n_b_faces = domain->mesh->n_b_faces; // Already defined
    // Pointer to array containing surface normals of all boundary faces
    // (L2-norm equals area of the face)
    const cs_real_t *b_f_face_normal =
      domain->mesh_quantities->b_f_face_normal;

    // Initilization of net pressure force vector
    cs_real_3_t total_b_p_forces = {0., 0., 0.};
    // Number and array of boundary faces (to use only certain parts,
    // containing boundary faces that meet a certain criteria)
    cs_lnum_t n_elts_p, *lst_elts_p;
    // Allocate memory for n_b_faces elements of type cs_lnum_t,
    // in lst_elts pointer (to fill a list of boundary faces)
    BFT_MALLOC(lst_elts_p, n_b_faces, cs_lnum_t);
    // Fill a list of boundary faces verifying a given selection criteria
    // For this case, using a group name of mesh entities, only faces on GTS
    // (name set at snappyHexMeshDict)
    // Creo que se pueden usar operadores logicos
    cs_selector_get_b_face_list("gtsModel", &n_elts_p, lst_elts_p);

    /* compute static pressure on selected boundary faces */
    // Array to store static pressure values on selected boundary faces
    cs_real_t *p_b_val;
    BFT_MALLOC(p_b_val, n_elts_p, cs_real_t);
    // Compute pressure on a specific boundary region
    // (A value for each face)
    cs_post_b_pressure(n_elts_p, lst_elts_p, p_b_val);

    // Loop over GTS boundary faces
    for (cs_lnum_t i_elt = 0; i_elt < n_elts_p; i_elt++) {
      cs_lnum_t face_id = lst_elts_p[i_elt];
      for (int ii = 0; ii < 3; ii++)
	// Sum of the values of forces due to static pressure on each face of
	// the selected boundary
        total_b_p_forces[ii] += p_b_val[i_elt]*b_f_face_normal[face_id*3+ii];
    }
    // Free resources
    BFT_FREE(lst_elts_p);
    BFT_FREE(p_b_val);

    /* parallel sum */
    cs_parall_sum(3, CS_DOUBLE, total_b_p_forces);
  
  
  // Get current physical time
  cs_real_t time = cs_glob_time_step->t_cur;
  // Get current iteration number
  int nt = cs_glob_time_step->nt_cur;
  
  // Manage file to write in
  FILE *file = NULL;
  if(cs_glob_rank_id < 1) {
      file = fopen("user_data_forces.csv", "a");
      if(nt==1){
	fprintf(file,"iteration, time, drag coeff, lift coeff, side force coeff, net force [x], net force [y], net force [z], net static pressure force [x], net static pressure force [y], net static pressure force [z] \n");
      }
      fprintf(file,"%d, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g \n",nt,time,c_d,c_l,c_s,
	      total_b_forces[0],total_b_forces[1],total_b_forces[2],
	      total_b_p_forces[0],total_b_p_forces[1],total_b_p_forces[2]);
      if (file != NULL) fclose(file);
    }
    
  //
  // yplus control over GTS surface
  //
  // yplus is a face scalar field
  cs_field_t *yplus_b = cs_field_by_name_try("yplus");
  // Number of boundary faces
  //const cs_lnum_t n_b_faces = domain->mesh->n_b_faces; // Already defined
  // Field of boundary faces center of gravity (cog)
  const cs_real_3_t *cdgfbo = (const cs_real_3_t *)domain->mesh_quantities->b_face_cog;
  // Field of boundary faces surface (surf)
  const cs_real_t *surfbn = cs_glob_mesh_quantities->b_face_surf;
  // Number and array of boundary faces (to use only certain parts,
  // containing boundary faces that meet a certain criteria)
  cs_lnum_t n_elts_y, *lst_elts_y;
  // Allocate memory for n_b_faces elements of type cs_lnum_t,
  // in lst_elts pointer (to fill a list of boundary faces)
  BFT_MALLOC(lst_elts_y, n_b_faces, cs_lnum_t);
  // Fill a list of boundary faces verifying a given selection criteria
  // For this case, using a group name of mesh entities, only faces on GTS
  // (name set at snappyHexMeshDict)
  // Creo que se pueden usar operadores logicos
  cs_selector_get_b_face_list("gtsModel", &n_elts_y, lst_elts_y);
  // yplus variables definition
  cs_real_t yplus_min=INFINITY, yplus_max=-INFINITY, 
    yplus_prom=0, yplus_max_cog[3]={-INFINITY,-INFINITY,-INFINITY},
    gts_area=0;
  // Loop over GTS boundary faces
    for (cs_lnum_t i_elt = 0; i_elt < n_elts_y; i_elt++) {
      cs_lnum_t face_id = lst_elts_y[i_elt];
      // Compute yplus_min
      yplus_min = cs_math_fmin(yplus_min,yplus_b->val[face_id]);
      // Compute yplus_max and save position on yplus_max_cog
      if(yplus_b->val[face_id]>yplus_max)
      {
	  yplus_max = yplus_b->val[face_id];
	  for(int j=0;j<3;j++) yplus_max_cog[j]=cdgfbo[face_id][j];
      }
      // Compute yplus_prom
      yplus_prom += yplus_b->val[face_id]*surfbn[face_id];
      gts_area += surfbn[face_id];
    }
    // Reduce variables over procceses
    cs_parall_min(1,CS_DOUBLE,&yplus_min);
    cs_parall_max_loc_vals(3,&yplus_max,yplus_max_cog);
    cs_parall_sum(1,CS_DOUBLE,&yplus_prom);
    cs_parall_sum(1,CS_DOUBLE,&gts_area);
    yplus_prom = yplus_prom/gts_area;
    // Manage file to write in
    FILE *file_yplus = NULL;
    if(cs_glob_rank_id < 1) {
	file_yplus = fopen("user_data_yplus.csv", "a");
	if(nt==1){
	  fprintf(file_yplus,"iteration, time, yplus_min, yplus_max, yplus_max_pos_x, yplus_max_pos_y, yplus_max_pos_z, yplus_prom \n");
	}
	fprintf(file_yplus,"%d, %g, %g, %g, %g, %g, %g, %g \n",nt,time,yplus_min,yplus_max,
		yplus_max_cog[0],yplus_max_cog[1],yplus_max_cog[2],yplus_prom);
	if (file_yplus != NULL) fclose(file_yplus);
      }
    // Free resources
    BFT_FREE(lst_elts_y);
}

END_C_DECLS
