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

#ifdef FIX_CLASS

FixStyle(reax/c/species/sunway,FixReaxCSpeciesSunway)

#else

#ifndef LMP_FIX_REAXC_SPECIES_SUNWAY_H
#define LMP_FIX_REAXC_SPECIES_SUNWAY_H

#include "fix.h"
#include "pointers.h"

#define BUFLEN 1000
namespace REAXC_SUNWAY_NS{
  typedef struct {
    double x, y, z;
  } AtomCoord;
}

namespace LAMMPS_NS {
  //using namespace REAXC_SUNWAY_NS;
class FixReaxCSpeciesSunway : public Fix {
 public:
  FixReaxCSpeciesSunway(class LAMMPS *, int, char **);
  virtual ~FixReaxCSpeciesSunway();
  int setmask();
  virtual void init();
  void init_list(int, class NeighList *);
  void setup(int);
  void post_integrate();
  double compute_vector(int);

 protected:
  int me, nprocs, nmax, nlocal, ntypes, ntotal;
  int nrepeat, nfreq, posfreq;
  int Nmoltype, vector_nmole, vector_nspec;
  int *Name, *MolName, *NMol, *nd, *MolType, *molmap;
  double *clusterID;
  REAXC_SUNWAY_NS::AtomCoord *x0;

  double bg_cut;
  double **BOCut;
  char **tmparg;

  FILE *fp, *pos;
  int eleflag, posflag, multipos, padflag, setupflag;
  int singlepos_opened, multipos_opened;
  char *ele, **eletype, *filepos;

  void Output_ReaxC_Bonds(bigint, FILE *);
  void create_compute();
  void create_fix();
  REAXC_SUNWAY_NS::AtomCoord chAnchor(REAXC_SUNWAY_NS::AtomCoord, REAXC_SUNWAY_NS::AtomCoord);
  virtual void FindMolecule();
  void SortMolecule(int &);
  void FindSpecies(int, int &);
  void WriteFormulas(int, int);
  int CheckExistence(int, int);

  int nint(const double &);
  int pack_forward_comm(int, int *, double *, int, int *);
  void unpack_forward_comm(int, int, double *);
  void OpenPos();
  void WritePos(int, int);
  double memory_usage();

  bigint nvalid;

  class NeighList *list;
  class FixAveAtom *f_SPECBOND;
  class PairReaxCSunway *reaxc;

};
}

#endif
#endif
