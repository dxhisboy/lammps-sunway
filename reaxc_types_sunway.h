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

#ifndef __REAX_TYPES_SUNWAY_H_
#define __REAX_TYPES_SUNWAY_H_

#include "lmptype.h"

#include <ctype.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sys/time.h"
#include <time.h>
#include "pair.h"
namespace REAXC_SUNWAY_NS{
#include "reaxc_ctypes_sunway.h"
  using LAMMPS_NS::Pair;

  struct _reax_system
  {
    reax_interaction reax_param;

    rc_bigint        bigN;
    int              n, N, numH;
    int              local_cap, total_cap, gcell_cap, Hcap;
    int              est_recv, est_trans, max_recved;
    int              wsize, my_rank, num_nbrs;
    ivec             my_coords;
    neighbor_proc    my_nbrs[REAX_MAX_NBRS];
    int             *global_offset;
    simulation_box   big_box, my_box, my_ext_box;
    grid             my_grid;
    boundary_cutoff  bndry_cuts;
    reax_atom       *my_atoms;
    double **x;
    class Pair *pair_ptr;
    int my_bonds;
    int mincap;
    double safezone, saferzone;

    int omp_active;
  public:
    void to_c_sys(reax_system_c *csys){
      csys->reax_param    = reax_param;
      csys->bigN          = bigN;
      csys->n             = n;
      csys->N             = N;
      csys->numH          = numH;
      csys->local_cap     = local_cap;
      csys->total_cap     = total_cap;
      csys->gcell_cap     = gcell_cap;
      csys->Hcap          = Hcap;
      csys->est_recv      = est_recv;
      csys->est_trans     = est_trans;
      csys->max_recved    = max_recved;
      csys->wsize         = wsize;
      csys->my_rank       = my_rank;
      csys->num_nbrs      = num_nbrs;
      csys->my_coords[0]  = my_coords[0];
      csys->my_coords[1]  = my_coords[1];
      csys->my_coords[2]  = my_coords[2];
      csys->global_offset = global_offset;
      csys->big_box       = big_box;
      csys->my_box        = my_box;
      csys->my_ext_box    = my_ext_box;
      csys->my_grid       = my_grid;
      csys->bndry_cuts    = bndry_cuts;
      csys->my_atoms      = my_atoms;
      csys->my_bonds      = my_bonds;
      csys->mincap        = mincap;
      csys->safezone      = safezone;
      csys->saferzone     = saferzone;
      csys->x             = (double(*)[3])(void*)x[0];
      if (pair_ptr->vflag_atom)
        csys->vatom         = (double(*)[6])(void*)(pair_ptr->vatom[0]);
      if (pair_ptr->eflag_atom)
      csys->eatom         = pair_ptr->eatom;
      csys->eng_vdwl      = pair_ptr->eng_vdwl;
      csys->eng_coul      = pair_ptr->eng_coul;
      csys->eflag_atom    = pair_ptr->eflag_atom;
      csys->vflag_atom    = pair_ptr->vflag_atom;
      csys->eflag_global  = pair_ptr->eflag_global;
      csys->vflag_global  = pair_ptr->vflag_global;
      csys->evflag        = pair_ptr->evflag;
      memcpy(csys->virial , pair_ptr->virial , sizeof(double) * 6);
      memcpy(csys->my_nbrs, my_nbrs, sizeof(neighbor_proc) * REAX_MAX_NBRS);
    }
    void from_c_sys(reax_system_c *csys){
      reax_param    = csys->reax_param;
      bigN          = csys->bigN;
      n             = csys->n;
      N             = csys->N;
      numH          = csys->numH;
      local_cap     = csys->local_cap;
      total_cap     = csys->total_cap;
      gcell_cap     = csys->gcell_cap;
      Hcap          = csys->Hcap;
      est_recv      = csys->est_recv;
      est_trans     = csys->est_trans;
      max_recved    = csys->max_recved;
      wsize         = csys->wsize;
      my_rank       = csys->my_rank;
      num_nbrs      = csys->num_nbrs;
      my_coords[0]  = csys->my_coords[0];
      my_coords[1]  = csys->my_coords[1];
      my_coords[2]  = csys->my_coords[2];
      global_offset = csys->global_offset;
      big_box       = csys->big_box;
      my_box        = csys->my_box;
      my_ext_box    = csys->my_ext_box;
      my_grid       = csys->my_grid;
      bndry_cuts    = csys->bndry_cuts;
      my_atoms      = csys->my_atoms;
      my_bonds      = csys->my_bonds;
      mincap        = csys->mincap;
      safezone      = csys->safezone;
      saferzone     = csys->saferzone;
      pair_ptr->eng_vdwl = csys->eng_vdwl;
      pair_ptr->eng_coul = csys->eng_coul;

      memcpy(pair_ptr->virial , csys->virial , sizeof(double) * 6);
      memcpy(my_nbrs, csys->my_nbrs, sizeof(neighbor_proc) * REAX_MAX_NBRS);
    }

  };
  typedef _reax_system reax_system;

  /* function pointer defs */
  typedef void (*evolve_function)(reax_system*, control_params*,
                                  simulation_data*, storage*, reax_list**,
                                  output_controls*, mpi_datatypes* );

  typedef void (*interaction_function) (reax_system*, control_params*,
                                        simulation_data*, storage*,
                                        reax_list**, output_controls*);

  typedef void (*print_interaction)(reax_system*, control_params*,
                                    simulation_data*, storage*,
                                    reax_list**, output_controls*);

  typedef double (*lookup_function)(double);

  typedef void (*message_sorter) (reax_system*, int, int, int, mpi_out_data*);
  typedef void (*unpacker) ( reax_system*, int, void*, int, neighbor_proc*, int );

  typedef void (*dist_packer) (void*, mpi_out_data*);
  typedef void (*coll_unpacker) (void*, void*, mpi_out_data*);

}

#endif
