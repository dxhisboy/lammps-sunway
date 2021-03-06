/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include <mpi.h>
#include "lammps.h"
#include "input.h"
#include "error.h"
#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include "gptl.h"
#if defined(LAMMPS_TRAP_FPE) && defined(_GNU_SOURCE)
#include <fenv.h>
#endif

using namespace LAMMPS_NS;

/* ----------------------------------------------------------------------
   main program to drive LAMMPS
------------------------------------------------------------------------- */
extern "C"{
  extern void pc_on_sig(int sig);
  extern void cinfo_init();
}
int __rank;
int main(int argc, char **argv)
{
  MPI_Init(&argc,&argv);
  GPTLinitialize();
  cinfo_init();
  MPI_Comm_rank(MPI_COMM_WORLD, &__rank);
  if (__rank == 0){
    signal(SIGUSR1, pc_on_sig);
    printf("%d\n", getpid());
  } else {
    signal(SIGUSR1, SIG_IGN);
  }
// enable trapping selected floating point exceptions.
// this uses GNU extensions and is only tested on Linux
// therefore we make it depend on -D_GNU_SOURCE, too.
  // long fpcr;
  // asm volatile("rfpcr %0\n\t" : "=r"(fpcr));
  // fpcr &= (~3L);
  // asm volatile("wfpcr %0\n\t" : :, "r"(fpcr));
#if defined(LAMMPS_TRAP_FPE) && defined(_GNU_SOURCE)
  fesetenv(FE_NOMASK_ENV);
  fedisableexcept(FE_ALL_EXCEPT);
  feenableexcept(FE_DIVBYZERO);
  feenableexcept(FE_INVALID);
  feenableexcept(FE_OVERFLOW);
#endif
#ifdef LAMMPS_EXCEPTIONS
  try {
    LAMMPS *lammps = new LAMMPS(argc,argv,MPI_COMM_WORLD);
    lammps->input->file();
    delete lammps;
  } catch(LAMMPSAbortException & ae) {
    MPI_Abort(ae.universe, 1);
  } catch(LAMMPSException & e) {
    MPI_Finalize();
    exit(1);
  }
#else
  LAMMPS *lammps = new LAMMPS(argc,argv,MPI_COMM_WORLD);
  lammps->input->file();
  delete lammps;
#endif
  GPTLpr_file("LammpsTime");
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
}
