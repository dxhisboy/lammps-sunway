/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Hasan Metin Aktulga, Purdue University
   (now at Lawrence Berkeley National Laboratory, hmaktulga@lbl.gov)

   Please cite the related publication:
   H. M. Aktulga, J. C. Fogarty, S. A. Pandit, A. Y. Grama,
   "Parallel Reactive Molecular Dynamics: Numerical Methods and
   Algorithmic Techniques", Parallel Computing, in press.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(reax/c/sunway,PairReaxCSunway)

#else

#ifndef LMP_PAIR_REAXC_SUNWAY_H
#define LMP_PAIR_REAXC_SUNWAY_H

#include "pair.h"
#include "neighbor.h"
#include "reaxc_types_sunway.h"
namespace LAMMPS_NS {
class PairReaxCSunway : public Pair {
 public:
  PairReaxCSunway(class LAMMPS *);
  ~PairReaxCSunway();
  void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  virtual void init_style();
  virtual void init_list(int, NeighList*);
  double init_one(int, int);
  void *extract(const char *, int &);
  int fixbond_flag, fixspecies_flag;
  int **tmpid;
  double **tmpbo,**tmpr;

  REAXC_SUNWAY_NS::control_params *control;
  REAXC_SUNWAY_NS::reax_system *system;
  REAXC_SUNWAY_NS::output_controls *out_control;
  REAXC_SUNWAY_NS::simulation_data *data;
  REAXC_SUNWAY_NS::storage *workspace;
  REAXC_SUNWAY_NS::reax_list *lists;
  REAXC_SUNWAY_NS::mpi_datatypes *mpi_data;

  bigint ngroup;

 protected:
  double cutmax;
  int nelements;                // # of unique elements
  char **elements;              // names of unique elements
  int *map;
  class FixReaxCSunway *fix_reax;

  double *chi,*eta,*gamma;
  int qeqflag;
  int setup_flag;
  int firstwarn;

  void allocate();
  void setup();
  void create_compute();
  void create_fix();
  void write_reax_atoms();
  void get_distance(REAXC_SUNWAY_NS::rvec, REAXC_SUNWAY_NS::rvec, double *, REAXC_SUNWAY_NS::rvec *);
  //void set_far_nbr(REAXC_SUNWAY_NS::far_neighbor_data *, int, double, REAXC_SUNWAY_NS::rvec);
  int estimate_reax_lists();
  int write_reax_lists();
  int write_reax_full_lists();
  void read_reax_forces(int);

  int nmax;
  void FindBond();
  double memory_usage();

};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Too many ghost atoms

Number of ghost atoms has increased too much during simulation and has exceeded
the size of reax/c arrays.  Increase safe_zone and min_cap in pair_style reax/c
command

*/
