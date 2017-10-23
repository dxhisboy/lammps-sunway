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

#ifndef LMP_NEIGH_LIST_H
#define LMP_NEIGH_LIST_H

#include "pointers.h"
#include "my_page.h"

namespace LAMMPS_NS {

class NeighList : protected Pointers {
 public:
  int index;                   // index of which neigh list this is
                               // also indexes the request it came from
                               // and the npair list of NPair classes

  int bin_method;        // 0 if no binning, else 1-N index into binnames
  int stencil_method;    // 0 if no stencil, else 1-N index into stencilnames
  int pair_method;       // 0 if no pair, else 1-N index into pairnames

  // settings from NeighRequest

  int occasional;                  // 0 if build every reneighbor, 1 if not
  int ghost;                       // 1 if list stores neighbors of ghosts
  int ssa;                         // 1 if list stores Shardlow data
  int copy;                        // 1 if this list is (host) copied from another list
  int copymode;                    // 1 if this is a Kokkos on-device copy
  int dnum;                        // # of doubles per neighbor, 0 if none

  // data structs to store neighbor pairs I,J and associated values

  int inum;                        // # of I atoms neighbors are stored for
  int gnum;                        // # of ghost atoms neighbors are stored for
  int *ilist;                      // local indices of I atoms
  int *numneigh;                   // # of J neighbors for each I atom
  int **firstneigh;                // ptr to 1st J int value of each I atom
  double **firstdouble;            // ptr to 1st J double value of each I atom
  int maxatom;                     // size of allocated per-atom arrays

  int pgsize;                      // size of each page
  int oneatom;                     // max size for one atom
  MyPage<int> *ipage;              // pages of neighbor indices
  MyPage<double> *dpage;           // pages of neighbor doubles, if dnum > 0

  // atom types to skip when building list
  // copied info from corresponding request into realloced vec/array

  int *iskip;         // iskip[i] = 1 if atoms of type I are not in list
  int **ijskip;       // ijskip[i][j] = 1 if pairs of type I,J are not in list

  // settings and pointers for related neighbor lists and fixes

  NeighList *listcopy;          // me = copy list, point to list I copy from
  NeighList *listskip;          // me = skip list, point to list I skip from
  NeighList *listfull;          // me = half list, point to full I derive from

  NeighList *listhistory;       // list storing neigh history
  class Fix *fix_history;       // fix that stores history info

  int respamiddle;              // 1 if this respaouter has middle list
  NeighList *listinner;         // me = respaouter, point to respainner
  NeighList *listmiddle;        // me = respaouter, point to respamiddle

  class Fix *fix_bond;          // fix that stores bond info

  // Kokkos package

  int kokkos;                   // 1 if list stores Kokkos data
  ExecutionSpace execution_space;

  // USER-DPD package and Shardlow Splitting Algorithm (SSA) support

  uint16_t (*ndxAIR_ssa)[8]; // for each atom, last neighbor index of each AIR

  // methods

  NeighList(class LAMMPS *);
  virtual ~NeighList();
  void post_constructor(class NeighRequest *);
  void setup_pages(int, int);           // setup page data structures
  void grow(int,int);                   // grow all data structs
  void print_attributes();              // debug routine
  int get_maxlocal() {return maxatom;}
  bigint memory_usage();
};

}

#endif
