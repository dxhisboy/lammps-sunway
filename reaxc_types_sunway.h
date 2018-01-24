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
