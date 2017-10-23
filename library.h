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

/*
   C or Fortran style library interface to LAMMPS
   new LAMMPS-specific functions can be added
*/

#include <mpi.h>

/* ifdefs allow this file to be included in a C program */

#ifdef __cplusplus
extern "C" {
#endif

void lammps_open(int, char **, MPI_Comm, void **);
void lammps_open_no_mpi(int, char **, void **);
void lammps_close(void *);
int  lammps_version(void *);
void lammps_file(void *, char *);
char *lammps_command(void *, char *);
void lammps_commands_list(void *, int, char **);
void lammps_commands_string(void *, char *);
void lammps_free(void *);

int lammps_extract_setting(void *, char *);
void *lammps_extract_global(void *, char *);
void lammps_extract_box(void *, double *, double *, 
                        double *, double *, double *, int *, int *);
void *lammps_extract_atom(void *, char *);
void *lammps_extract_compute(void *, char *, int, int);
void *lammps_extract_fix(void *, char *, int, int, int, int);
void *lammps_extract_variable(void *, char *, char *);

void lammps_reset_box(void *, double *, double *, double, double, double);
int lammps_set_variable(void *, char *, char *);
double lammps_get_thermo(void *, char *);

int lammps_get_natoms(void *);
void lammps_gather_atoms(void *, char *, int, int, void *);
void lammps_scatter_atoms(void *, char *, int, int, void *);

// lammps_create_atoms() takes tagint and imageint as args
// ifdef insures they are compatible with rest of LAMMPS
// caller must match to how LAMMPS library is built

#ifdef LAMMPS_BIGBIG
void lammps_create_atoms(void *, int, int64_t *, int *, 
                         double *, double *, int64_t *, int);
#else
void lammps_create_atoms(void *, int, int *, int *, 
                         double *, double *, int *, int);
#endif

#ifdef LAMMPS_EXCEPTIONS
int lammps_has_error(void *);
int lammps_get_last_error_message(void *, char *, int);
#endif

#undef LAMMPS
#ifdef __cplusplus
}
#endif

/* ERROR/WARNING messages:

W: Library error in lammps_gather_atoms

This library function cannot be used if atom IDs are not defined
or are not consecutively numbered.

W: Library error in lammps_scatter_atoms

This library function cannot be used if atom IDs are not defined or
are not consecutively numbered, or if no atom map is defined.  See the
atom_modify command for details about atom maps.

*/
