/*----------------------------------------------------------------------
  PuReMD - Purdue ReaxFF Molecular Dynamics Program

  Copyright (2010) Purdue University
  Hasan Metin Aktulga, hmaktulga@lbl.gov
  Joseph Fogarty, jcfogart@mail.usf.edu
  Sagar Pandit, pandit@usf.edu
  Ananth Y Grama, ayg@cs.purdue.edu

  Please cite the related publication:
  H. M. Aktulga, J. C. Fogarty, S. A. Pandit, A. Y. Grama,
  "Parallel Reactive Molecular Dynamics: Numerical Methods and
  Algorithmic Techniques", Parallel Computing, in press.

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License as
  published by the Free Software Foundation; either version 2 of
  the License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  See the GNU General Public License for more details:
  <http://www.gnu.org/licenses/>.
  ----------------------------------------------------------------------*/

#include "pair_reaxc_sunway.h"
#include "reaxc_types_sunway.h"
#include "reaxc_nonbonded_sunway.h"
#include "reaxc_bond_orders_sunway.h"
#include "reaxc_list_sunway.h"
#include "reaxc_vector_sunway.h"
#include "reaxc_nonbonded_sw64.h"
#include "sunway.h"
#include "gptl.h"
namespace REAXC_SUNWAY_NS{

  static inline int llf(double *dvec){
    double SMALL = 0.0001;
    if (dvec[2] < -SMALL) return 1;
    if (dvec[2] < SMALL){
      if (dvec[1] < -SMALL) return 1;
      if (dvec[1] < SMALL && dvec[0] < -SMALL) return 1;
    }
    return 0;
  }
  static inline int urb(double *dvec){
    double SMALL = 0.0001;
    if (dvec[2] > SMALL) return 1;
    if (dvec[2] > -SMALL){
      if (dvec[1] > SMALL) return 1;
      if (dvec[1] > -SMALL && dvec[0] > SMALL) return 1;
    }
    return 0;
  }
  void vdW_Coulomb_Energy_Full( reax_system *system, control_params *control,
                                simulation_data *data, storage *workspace,
                                reax_list **lists, output_controls *out_control )
  {
    reax_system_c csys;
    system->to_c_sys(&csys);
    GPTLstart("reaxc vdw coul full c");
    vdW_Coulomb_Energy_Full_C(&csys, control, data, workspace, lists, out_control);
    GPTLstop("reaxc vdw coul full c");
    system->from_c_sys(&csys);
    Compute_Polarization_Energy( system, data );
    return;
  }

  void Compute_Polarization_Energy( reax_system *system, simulation_data *data )
  {
    int  i, type_i;
    double q, en_tmp;

    data->my_en.e_pol = 0.0;
    for( i = 0; i < system->n; i++ ) {
      type_i = system->my_atoms[i].type;
      if (type_i < 0) continue;
      q = system->my_atoms[i].q;

      en_tmp = KCALpMOL_to_EV * (system->reax_param.sbp[type_i].chi * q +
                                 (system->reax_param.sbp[type_i].eta / 2.) * SQR(q));
      data->my_en.e_pol += en_tmp;

      /* tally into per-atom energy */
      if( system->pair_ptr->evflag)
        system->pair_ptr->ev_tally(i,i,system->n,1,0.0,en_tmp,0.0,0.0,0.0,0.0);
    }
  }


}
