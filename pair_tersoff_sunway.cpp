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

/* ----------------------------------------------------------------------
   Contributing author: Aidan Thompson (SNL)
------------------------------------------------------------------------- */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pair_tersoff_sunway.h"
#include "atom.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "force.h"
#include "comm.h"
#include "memory.h"
#include "error.h"
#include "sunway.h"
#include "pair_tersoff_sw64.h"
#include "gptl.h"

#include "math_const.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define MAXLINE 1024
#define DELTA 4

/* ---------------------------------------------------------------------- */

PairTersoffSunway::PairTersoffSunway(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 0;
  restartinfo = 0;
  one_coeff = 1;
  manybody_flag = 1;

  nelements = 0;
  elements = NULL;
  nparams = maxparam = 0;
  params = NULL;
  elem2param = NULL;
  map = NULL;
  cutotal = 0;
  maxshort = 10;
  neighshort = NULL;
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

PairTersoffSunway::~PairTersoffSunway()
{
  if (copymode) return;

  if (elements)
    for (int i = 0; i < nelements; i++) delete [] elements[i];
  delete [] elements;
  memory->destroy(params);
  memory->destroy(elem2param);
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(neighshort);
    delete [] map;
  }
  if (cutotal){
    memory->destroy(shortidx);
    memory->destroy(fend);
    memory->destroy(ftmp);
    memory->destroy(fdone);
    memory->destroy(numshort);
    memory->destroy(atom_in);
  }
}

#define THIRD 0.3333333333333333333
void PairTersoffSunway::v_tally3rd(int i, int vflag_global, int vflag_atom,
				double *fi, double *fj, double *drik, double *drjk)
{
  double v[6];

  v[0] = THIRD * (drik[0]*fi[0] + drjk[0]*fj[0]);
  v[1] = THIRD * (drik[1]*fi[1] + drjk[1]*fj[1]);
  v[2] = THIRD * (drik[2]*fi[2] + drjk[2]*fj[2]);
  v[3] = THIRD * (drik[0]*fi[1] + drjk[0]*fj[1]);
  v[4] = THIRD * (drik[0]*fi[2] + drjk[0]*fj[2]);
  v[5] = THIRD * (drik[1]*fi[2] + drjk[1]*fj[2]);

  if (vflag_global){
    virial[0] += v[0];
    virial[1] += v[1];
    virial[2] += v[2];
    virial[3] += v[3];
    virial[4] += v[4];
    virial[5] += v[5];
  }
  if (vflag_atom){
    vatom[i][0] += v[0]; vatom[i][1] += v[1]; vatom[i][2] += v[2];
    vatom[i][3] += v[3]; vatom[i][4] += v[4]; vatom[i][5] += v[5];
  }
}

void PairTersoffSunway::allocate_tmp(int ntotal){
  if (cutotal){
    memory->destroy(fend);
    memory->destroy(ftmp);
    memory->destroy(fdone);
    memory->destroy(shortidx);
    memory->destroy(numshort);
    memory->destroy(atom_in);
  } else
    cutotal = 1;
  while (cutotal < ntotal)
    cutotal *= 2;
  memory->create(shortidx, (cutotal + 1) * maxshort, "pair:shortidx");
  memory->create(fend, cutotal * maxshort + 1, "pair:fend");
  memory->create(ftmp, cutotal + 1, "pair:ftmp");
  memory->create(fdone, cutotal + 8, "pair:fdone");
  memory->create(numshort, cutotal + 1, "pair:numshort");
  memory->create(atom_in, cutotal + 8, "pair:atom_in");
  fend_base = (long)fend;
  ftmp_base = (long)ftmp;
  fdone_base = (long)fdone;
  atom_in_base = (long)atom_in;
}
/* ---------------------------------------------------------------------- */
void PairTersoffSunway::compute(int eflag, int vflag){
  int i,j,k,ii,jj,kk,inum,jnum,knum, gnum, allnum;
  int itype,jtype,ktype,iparam_ij,iparam_ijk, iparam_jik, iparam_jki;
  tagint itag,jtag;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,rsq1,rsq2;
  double delr1[3],delr2[3],fi[3],fj[3],fk[3];
  double zeta_ij,zeta_ji,prefactor_ij, prefactor_ji;
  int *ilist,*jlist,*klist,*numneigh,**firstneigh;
  //short_neigh_t *jlist_short, *klist_short, *jshort, *kshort, *neighshort;;
  double dij[3], dji[3], djk[3], dik[3], r2ij, r2jk, r2ik;
  evdwl = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = vflag_atom = 0;

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int ntotal = atom->nlocal + atom->nghost;
  int newton_pair = force->newton_pair;
  const double cutshortsq = cutmax*cutmax;
  if (ntotal > cutotal)
    allocate_tmp(ntotal);
  inum = list->inum;
  gnum = list->gnum;
  allnum = inum + gnum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  int neightotal;
  int *firstshort;
  
  neightotal = 0;
  for (ii = 0; ii < allnum; ii ++){
    neightotal += numneigh[ii];
  }
  //short_neigh_t *shortlist;

  tersoff_param_t *cparams;
  memory->create(cparams, nelements * nelements * nelements, "pair:cparams");
  int t0, t1, t2;
  for (t0 = 0; t0 < nelements; t0 ++){
    for (t1 = 0; t1 < nelements; t1 ++){
      for (t2 = 0; t2 < nelements; t2 ++){
        int idx_new = (t0 * nelements + t1) * nelements + t2;
        int idx_old = elem2param[t0][t1][t2];
        cparams[idx_new].lam1        = params[idx_old].lam1       ;
        cparams[idx_new].lam2        = params[idx_old].lam2       ;
        cparams[idx_new].lam3        = params[idx_old].lam3       ;
        cparams[idx_new].c           = params[idx_old].c          ;
        cparams[idx_new].d           = params[idx_old].d          ;
        cparams[idx_new].h           = params[idx_old].h          ;
        cparams[idx_new].gamma       = params[idx_old].gamma      ;
        //cparams[idx_new].powerm      = params[idx_old].powerm     ;
        cparams[idx_new].powern      = params[idx_old].powern     ;
        cparams[idx_new].half_powern_inv = 1 / (2.0 * params[idx_old].powern);
        cparams[idx_new].beta        = params[idx_old].beta       ;
        cparams[idx_new].biga        = params[idx_old].biga       ;
        cparams[idx_new].bigb        = params[idx_old].bigb       ;
        cparams[idx_new].bigd        = params[idx_old].bigd       ;
        cparams[idx_new].bigdinv     = 1 / params[idx_old].bigd   ;
        cparams[idx_new].c2divd2     = params[idx_old].c * params[idx_old].c /
          (params[idx_old].d * params[idx_old].d);
        cparams[idx_new].bigr        = params[idx_old].bigr       ;
        //cparams[idx_new].cut         = params[idx_old].cut        ;
        cparams[idx_new].cutsq       = params[idx_old].cutsq      ;
        cparams[idx_new].c1          = params[idx_old].c1         ;
        cparams[idx_new].c4          = params[idx_old].c4         ;
        cparams[idx_new].powermint   = params[idx_old].powermint  ;

      }
    }
  }
  // for (i = 0; i < nparams; i ++){
  //   cparams[i].lam1        = params[i].lam1       ;
  //   cparams[i].lam2        = params[i].lam2       ;
  //   cparams[i].lam3        = params[i].lam3       ;
  //   cparams[i].c           = params[i].c          ;
  //   cparams[i].d           = params[i].d          ;
  //   cparams[i].h           = params[i].h          ;
  //   cparams[i].gamma       = params[i].gamma      ;
  //   cparams[i].powerm      = params[i].powerm     ;
  //   cparams[i].powern      = params[i].powern     ;
  //   cparams[i].half_powern_inv = 1 / (2.0 * params[i].powern);
  //   cparams[i].beta        = params[i].beta       ;
  //   cparams[i].biga        = params[i].biga       ;
  //   cparams[i].bigb        = params[i].bigb       ;
  //   cparams[i].bigd        = params[i].bigd       ;
  //   cparams[i].bigdinv     = 1 / params[i].bigd   ;
  //   cparams[i].c2divd2     = params[i].c * params[i].c / (params[i].d * params[i].d);
  //   cparams[i].bigr        = params[i].bigr       ;
  //   cparams[i].cut         = params[i].cut        ;
  //   cparams[i].cutsq       = params[i].cutsq      ;
  //   cparams[i].c1          = params[i].c1         ;
  //   cparams[i].c2          = params[i].c2         ;
  //   cparams[i].c3          = params[i].c3         ;
  //   cparams[i].c4          = params[i].c4         ;
  //   cparams[i].ielement    = params[i].ielement   ;
  //   cparams[i].jelement    = params[i].jelement   ;
  //   cparams[i].kelement    = params[i].kelement   ;
  //   cparams[i].powermint   = params[i].powermint  ;
  //   cparams[i].Z_i         = params[i].Z_i        ;
  //   cparams[i].Z_j         = params[i].Z_j        ;
  //   cparams[i].ZBLcut      = params[i].ZBLcut     ;
  //   cparams[i].ZBLexpscale = params[i].ZBLexpscale;
  //   cparams[i].c5          = params[i].c5         ;
  //   cparams[i].ca1         = params[i].ca1        ;
  //   cparams[i].ca4         = params[i].ca4        ;
  //   cparams[i].powern_del  = params[i].powern_del ;
  //   cparams[i].c0          = params[i].c0         ;
  // }

  pair_tersoff_compute_param_t pm;
  pm.cutshortsq = cutmax * cutmax;
  pm.fdone_base = fdone_base;
  pm.ftmp_base  = ftmp_base;
  pm.fend_base  = fend_base;
  pm.atom_in_base = atom_in_base;
  pm.rank       = comm->me;
  pm.ilist      = ilist;
  pm.numneigh   = numneigh;
  pm.firstshort = firstshort;
  pm.firstneigh = firstneigh;
  pm.shortidx   = shortidx;
  //pm.shortlist  = shortlist;
  pm.elem2param = elem2param[0][0];
  pm.map        = map;
  pm.params     = cparams;
  pm.x          = (double(*)[3])(void*)x[0];
  pm.f          = (double(*)[3])(void*)f[0];
  if (vflag_atom)
    pm.vatom      = (double(*)[6])(void*)vatom[0];
  if (eflag_atom)
    pm.eatom      = eatom;
  pm.type       = type;
  pm.nlocal     = nlocal;
  pm.nghost     = atom->nghost;
  pm.ntotal     = nlocal + atom->nghost;
  pm.inum       = inum;
  pm.gnum       = gnum;
  pm.ntypes     = atom->ntypes;
  pm.nelements  = nelements;
  pm.nparams    = nparams;
  pm.eng_vdwl   = 0;
  pm.eng_coul   = 0;
  pm.virial[0]  = 0;
  pm.virial[1]  = 0;
  pm.virial[2]  = 0;
  pm.virial[3]  = 0;
  pm.virial[4]  = 0;
  pm.virial[5]  = 0;
  pm.eflag      = eflag;
  pm.vflag      = vflag;
  pm.evflag     = evflag;

  pm.eflag_global = eflag_global;
  pm.vflag_global = vflag_global;
  pm.eflag_atom = eflag_atom;
  pm.vflag_atom = vflag_atom;
  pm.eflag_either = eflag_either;
  pm.vflag_either = vflag_either;
  pm.maxshort     = maxshort;
  pm.numshort     = numshort;
  double fxtmp,fytmp,fztmp;
  int numshorti, numshortj;
  GPTLstart("attractive");
  pair_tersoff_compute_attractive(&pm);
  GPTLstop("attractive");
  if (eflag_global){
    eng_vdwl = pm.eng_vdwl;
    eng_coul = pm.eng_coul;
  }

  if (vflag_global){
    for (i = 0; i < 6; i ++){
      virial[i] = pm.virial[i];
    }
  }
  memory->destroy(cparams);
}

/* ---------------------------------------------------------------------- */

void PairTersoffSunway::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  memory->create(cutsq,n+1,n+1,"pair:cutsq");
  memory->create(neighshort,maxshort,"pair:neighshort");
  map = new int[n+1];
  allocate_tmp(atom->nmax);
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairTersoffSunway::settings(int narg, char **arg)
{
  if (narg != 0) error->all(FLERR,"Illegal pair_style command");
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairTersoffSunway::coeff(int narg, char **arg)
{
  int i,j,n;

  if (!allocated) allocate();

  if (narg != 3 + atom->ntypes)
    error->all(FLERR,"Incorrect args for pair coefficients");

  // insure I,J args are * *

  if (strcmp(arg[0],"*") != 0 || strcmp(arg[1],"*") != 0)
    error->all(FLERR,"Incorrect args for pair coefficients");

  // read args that map atom types to elements in potential file
  // map[i] = which element the Ith atom type is, -1 if NULL
  // nelements = # of unique elements
  // elements = list of element names

  if (elements) {
    for (i = 0; i < nelements; i++) delete [] elements[i];
    delete [] elements;
  }
  elements = new char*[atom->ntypes];
  for (i = 0; i < atom->ntypes; i++) elements[i] = NULL;

  nelements = 0;
  for (i = 3; i < narg; i++) {
    if (strcmp(arg[i],"NULL") == 0) {
      map[i-2] = -1;
      continue;
    }
    for (j = 0; j < nelements; j++)
      if (strcmp(arg[i],elements[j]) == 0) break;
    map[i-2] = j;
    if (j == nelements) {
      n = strlen(arg[i]) + 1;
      elements[j] = new char[n];
      strcpy(elements[j],arg[i]);
      nelements++;
    }
  }

  // read potential file and initialize potential parameters

  read_file(arg[2]);
  setup_params();

  // clear setflag since coeff() called once with I,J = * *

  n = atom->ntypes;
  for (i = 1; i <= n; i++)
    for (j = i; j <= n; j++)
      setflag[i][j] = 0;

  // set setflag i,j for type pairs where both are mapped to elements

  int count = 0;
  for (i = 1; i <= n; i++)
    for (j = i; j <= n; j++)
      if (map[i] >= 0 && map[j] >= 0) {
        setflag[i][j] = 1;
        count++;
      }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairTersoffSunway::init_style()
{
  if (atom->tag_enable == 0)
    error->all(FLERR,"Pair style Tersoff requires atom IDs");
  // if (force->newton_pair == 0)
  //   error->all(FLERR,"Pair style Tersoff requires newton pair on");

  // need a full neighbor list

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
  neighbor->requests[irequest]->ghost= 1;

  if (comm->cutghostuser < (2.0*cutmax + neighbor->skin) )
    comm->cutghostuser = 2.0*cutmax + neighbor->skin;
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairTersoffSunway::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  return cutmax;
}

/* ---------------------------------------------------------------------- */

void PairTersoffSunway::read_file(char *file)
{
  int params_per_line = 17;
  char **words = new char*[params_per_line+1];

  memory->sfree(params);
  params = NULL;
  nparams = maxparam = 0;

  // open file on proc 0

  FILE *fp;
  if (comm->me == 0) {
    fp = force->open_potential(file);
    if (fp == NULL) {
      char str[128];
      sprintf(str,"Cannot open Tersoff potential file %s",file);
      error->one(FLERR,str);
    }
  }

  // read each line out of file, skipping blank lines or leading '#'
  // store line of params if all 3 element tags are in element list

  int n,nwords,ielement,jelement,kelement;
  char line[MAXLINE],*ptr;
  int eof = 0;

  while (1) {
    if (comm->me == 0) {
      ptr = fgets(line,MAXLINE,fp);
      if (ptr == NULL) {
        eof = 1;
        fclose(fp);
      } else n = strlen(line) + 1;
    }
    MPI_Bcast(&eof,1,MPI_INT,0,world);
    if (eof) break;
    MPI_Bcast(&n,1,MPI_INT,0,world);
    MPI_Bcast(line,n,MPI_CHAR,0,world);

    // strip comment, skip line if blank

    if ((ptr = strchr(line,'#'))) *ptr = '\0';
    nwords = atom->count_words(line);
    if (nwords == 0) continue;

    // concatenate additional lines until have params_per_line words

    while (nwords < params_per_line) {
      n = strlen(line);
      if (comm->me == 0) {
        ptr = fgets(&line[n],MAXLINE-n,fp);
        if (ptr == NULL) {
          eof = 1;
          fclose(fp);
        } else n = strlen(line) + 1;
      }
      MPI_Bcast(&eof,1,MPI_INT,0,world);
      if (eof) break;
      MPI_Bcast(&n,1,MPI_INT,0,world);
      MPI_Bcast(line,n,MPI_CHAR,0,world);
      if ((ptr = strchr(line,'#'))) *ptr = '\0';
      nwords = atom->count_words(line);
    }

    if (nwords != params_per_line)
      error->all(FLERR,"Incorrect format in Tersoff potential file");

    // words = ptrs to all words in line

    nwords = 0;
    words[nwords++] = strtok(line," \t\n\r\f");
    while ((words[nwords++] = strtok(NULL," \t\n\r\f"))) continue;

    // ielement,jelement,kelement = 1st args
    // if all 3 args are in element list, then parse this line
    // else skip to next line

    for (ielement = 0; ielement < nelements; ielement++)
      if (strcmp(words[0],elements[ielement]) == 0) break;
    if (ielement == nelements) continue;
    for (jelement = 0; jelement < nelements; jelement++)
      if (strcmp(words[1],elements[jelement]) == 0) break;
    if (jelement == nelements) continue;
    for (kelement = 0; kelement < nelements; kelement++)
      if (strcmp(words[2],elements[kelement]) == 0) break;
    if (kelement == nelements) continue;

    // load up parameter settings and error check their values

    if (nparams == maxparam) {
      maxparam += DELTA;
      params = (Param *) memory->srealloc(params,maxparam*sizeof(Param),
                                          "pair:params");
    }

    params[nparams].ielement = ielement;
    params[nparams].jelement = jelement;
    params[nparams].kelement = kelement;
    params[nparams].powerm = atof(words[3]);
    params[nparams].gamma = atof(words[4]);
    params[nparams].lam3 = atof(words[5]);
    params[nparams].c = atof(words[6]);
    params[nparams].d = atof(words[7]);
    params[nparams].h = atof(words[8]);
    params[nparams].powern = atof(words[9]);
    params[nparams].beta = atof(words[10]);
    params[nparams].lam2 = atof(words[11]);
    params[nparams].bigb = atof(words[12]);
    params[nparams].bigr = atof(words[13]);
    params[nparams].bigd = atof(words[14]);
    params[nparams].lam1 = atof(words[15]);
    params[nparams].biga = atof(words[16]);

    // currently only allow m exponent of 1 or 3

    params[nparams].powermint = int(params[nparams].powerm);

    if (params[nparams].c < 0.0 || params[nparams].d < 0.0 ||
        params[nparams].powern < 0.0 || params[nparams].beta < 0.0 ||
        params[nparams].lam2 < 0.0 || params[nparams].bigb < 0.0 ||
        params[nparams].bigr < 0.0 ||params[nparams].bigd < 0.0 ||
        params[nparams].bigd > params[nparams].bigr ||
        params[nparams].lam1 < 0.0 || params[nparams].biga < 0.0 ||
        params[nparams].powerm - params[nparams].powermint != 0.0 ||
        (params[nparams].powermint != 3 && params[nparams].powermint != 1) ||
        params[nparams].gamma < 0.0)
      error->all(FLERR,"Illegal Tersoff parameter");

    nparams++;
  }

  delete [] words;
}

/* ---------------------------------------------------------------------- */

void PairTersoffSunway::setup_params()
{
  int i,j,k,m,n;

  // set elem2param for all element triplet combinations
  // must be a single exact match to lines read from file
  // do not allow for ACB in place of ABC

  memory->destroy(elem2param);
  memory->create(elem2param,nelements,nelements,nelements,"pair:elem2param");

  for (i = 0; i < nelements; i++)
    for (j = 0; j < nelements; j++)
      for (k = 0; k < nelements; k++) {
        n = -1;
        for (m = 0; m < nparams; m++) {
          if (i == params[m].ielement && j == params[m].jelement &&
              k == params[m].kelement) {
            if (n >= 0) error->all(FLERR,"Potential file has duplicate entry");
            n = m;
          }
        }
        if (n < 0) error->all(FLERR,"Potential file is missing an entry");
        elem2param[i][j][k] = n;
      }


  // compute parameter values derived from inputs

  for (m = 0; m < nparams; m++) {
    params[m].cut = params[m].bigr + params[m].bigd;
    params[m].cutsq = params[m].cut*params[m].cut;

    params[m].c1 = pow(2.0*params[m].powern*1.0e-16,-1.0/params[m].powern);
    params[m].c2 = pow(2.0*params[m].powern*1.0e-8,-1.0/params[m].powern);
    params[m].c3 = 1.0/params[m].c2;
    params[m].c4 = 1.0/params[m].c1;
  }

  // set cutmax to max of all params

  cutmax = 0.0;
  for (m = 0; m < nparams; m++)
    if (params[m].cut > cutmax) cutmax = params[m].cut;
}

/* ---------------------------------------------------------------------- */

void PairTersoffSunway::repulsive(Param *param, double rsq, double &fforce,
                            int eflag, double &eng)
{
  double r,tmp_fc,tmp_fc_d,tmp_exp;

  r = sqrt(rsq);
  tmp_fc = ters_fc(r,param);
  tmp_fc_d = ters_fc_d(r,param);
  tmp_exp = exp(-param->lam1 * r);
  fforce = -param->biga * tmp_exp * (tmp_fc_d - tmp_fc*param->lam1) / r;
  if (eflag) eng = tmp_fc * param->biga * tmp_exp;
}

/* ---------------------------------------------------------------------- */

double PairTersoffSunway::zeta(Param *param, double rsqij, double rsqik,
                         double *delrij, double *delrik)
{
  double rij,rik,costheta,arg,ex_delr;

  rij = sqrt(rsqij);
  rik = sqrt(rsqik);
  costheta = (delrij[0]*delrik[0] + delrij[1]*delrik[1] +
              delrij[2]*delrik[2]) / (rij*rik);

  if (param->powermint == 3) arg = pow(param->lam3 * (rij-rik),3.0);
  else arg = param->lam3 * (rij-rik);

  if (arg > 69.0776) ex_delr = 1.e30;
  else if (arg < -69.0776) ex_delr = 0.0;
  else ex_delr = exp(arg);

  return ters_fc(rik,param) * ters_gijk(costheta,param) * ex_delr;
}


/* ---------------------------------------------------------------------- */

void PairTersoffSunway::force_zeta(Param *param, double rsq, double zeta_ij,
                             double &fforce, double &prefactor,
                             int eflag, double &eng)
{
  double r,fa,fa_d,bij;

  r = sqrt(rsq);
  fa = ters_fa(r,param);
  fa_d = ters_fa_d(r,param);
  bij = ters_bij(zeta_ij,param);
  fforce = 0.5*bij*fa_d / r;
  prefactor = -0.5*fa * ters_bij_d(zeta_ij,param);
  if (eflag) eng = 0.5*bij*fa;
}

/* ----------------------------------------------------------------------
   attractive term
   use param_ij cutoff for rij test
   use param_ijk cutoff for rik test
------------------------------------------------------------------------- */

void PairTersoffSunway::attractive(Param *param, double prefactor,
                             double rsqij, double rsqik,
                             double *delrij, double *delrik,
                             double *fi, double *fj, double *fk)
{
  double rij_hat[3],rik_hat[3];
  double rij,rijinv,rik,rikinv;

  rij = sqrt(rsqij);
  rijinv = 1.0/rij;
  vec3_scale(rijinv,delrij,rij_hat);

  rik = sqrt(rsqik);
  rikinv = 1.0/rik;
  vec3_scale(rikinv,delrik,rik_hat);

  ters_zetaterm_d(prefactor,rij_hat,rij,rik_hat,rik,fi,fj,fk,param);
}

/* ---------------------------------------------------------------------- */

double PairTersoffSunway::ters_fc(double r, Param *param)
{
  double ters_R = param->bigr;
  double ters_D = param->bigd;

  if (r < ters_R-ters_D) return 1.0;
  if (r > ters_R+ters_D) return 0.0;
  return 0.5*(1.0 - sin(MY_PI2*(r - ters_R)/ters_D));
}

/* ---------------------------------------------------------------------- */

double PairTersoffSunway::ters_fc_d(double r, Param *param)
{
  double ters_R = param->bigr;
  double ters_D = param->bigd;

  if (r < ters_R-ters_D) return 0.0;
  if (r > ters_R+ters_D) return 0.0;
  return -(MY_PI4/ters_D) * cos(MY_PI2*(r - ters_R)/ters_D);
}

/* ---------------------------------------------------------------------- */

double PairTersoffSunway::ters_fa(double r, Param *param)
{
  if (r > param->bigr + param->bigd) return 0.0;
  return -param->bigb * exp(-param->lam2 * r) * ters_fc(r,param);
}

/* ---------------------------------------------------------------------- */

double PairTersoffSunway::ters_fa_d(double r, Param *param)
{
  if (r > param->bigr + param->bigd) return 0.0;
  return param->bigb * exp(-param->lam2 * r) *
    (param->lam2 * ters_fc(r,param) - ters_fc_d(r,param));
}

/* ---------------------------------------------------------------------- */

double PairTersoffSunway::ters_bij(double zeta, Param *param)
{
  double tmp = param->beta * zeta;
  if (tmp > param->c1) return 1.0/sqrt(tmp);
  if (tmp > param->c2)
    return (1.0 - pow(tmp,-param->powern) / (2.0*param->powern))/sqrt(tmp);
  if (tmp < param->c4) return 1.0;
  if (tmp < param->c3)
    return 1.0 - pow(tmp,param->powern)/(2.0*param->powern);
  return pow(1.0 + pow(tmp,param->powern), -1.0/(2.0*param->powern));
}

/* ---------------------------------------------------------------------- */

double PairTersoffSunway::ters_bij_d(double zeta, Param *param)
{
  double tmp = param->beta * zeta;
  if (tmp > param->c1) return param->beta * -0.5*pow(tmp,-1.5);
  if (tmp > param->c2)
    return param->beta * (-0.5*pow(tmp,-1.5) *
			  // error in negligible 2nd term fixed 9/30/2015
			  // (1.0 - 0.5*(1.0 +  1.0/(2.0*param->powern)) *
                          (1.0 - (1.0 +  1.0/(2.0*param->powern)) *
                           pow(tmp,-param->powern)));
  if (tmp < param->c4) return 0.0;
  if (tmp < param->c3)
    return -0.5*param->beta * pow(tmp,param->powern-1.0);

  double tmp_n = pow(tmp,param->powern);
  return -0.5 * pow(1.0+tmp_n, -1.0-(1.0/(2.0*param->powern)))*tmp_n / zeta;
}

/* ---------------------------------------------------------------------- */

void PairTersoffSunway::ters_zetaterm_d(double prefactor,
                                  double *rij_hat, double rij,
                                  double *rik_hat, double rik,
                                  double *dri, double *drj, double *drk,
                                  Param *param)
{
  double gijk,gijk_d,ex_delr,ex_delr_d,fc,dfc,cos_theta,tmp;
  double dcosdri[3],dcosdrj[3],dcosdrk[3];

  fc = ters_fc(rik,param);
  dfc = ters_fc_d(rik,param);
  if (param->powermint == 3) tmp = pow(param->lam3 * (rij-rik),3.0);
  else tmp = param->lam3 * (rij-rik);

  if (tmp > 69.0776) ex_delr = 1.e30;
  else if (tmp < -69.0776) ex_delr = 0.0;
  else ex_delr = exp(tmp);

  if (param->powermint == 3)
    ex_delr_d = 3.0*pow(param->lam3,3.0) * pow(rij-rik,2.0)*ex_delr;
  else ex_delr_d = param->lam3 * ex_delr;

  cos_theta = vec3_dot(rij_hat,rik_hat);
  gijk = ters_gijk(cos_theta,param);
  gijk_d = ters_gijk_d(cos_theta,param);
  costheta_d(rij_hat,rij,rik_hat,rik,dcosdri,dcosdrj,dcosdrk);

  // compute the derivative wrt Ri
  // dri = -dfc*gijk*ex_delr*rik_hat;
  // dri += fc*gijk_d*ex_delr*dcosdri;
  // dri += fc*gijk*ex_delr_d*(rik_hat - rij_hat);

  vec3_scale(-dfc*gijk*ex_delr,rik_hat,dri);
  vec3_scaleadd(fc*gijk_d*ex_delr,dcosdri,dri,dri);
  vec3_scaleadd(fc*gijk*ex_delr_d,rik_hat,dri,dri);
  vec3_scaleadd(-fc*gijk*ex_delr_d,rij_hat,dri,dri);
  vec3_scale(prefactor,dri,dri);

  // compute the derivative wrt Rj
  // drj = fc*gijk_d*ex_delr*dcosdrj;
  // drj += fc*gijk*ex_delr_d*rij_hat;

  vec3_scale(fc*gijk_d*ex_delr,dcosdrj,drj);
  vec3_scaleadd(fc*gijk*ex_delr_d,rij_hat,drj,drj);
  vec3_scale(prefactor,drj,drj);

  // compute the derivative wrt Rk
  // drk = dfc*gijk*ex_delr*rik_hat;
  // drk += fc*gijk_d*ex_delr*dcosdrk;
  // drk += -fc*gijk*ex_delr_d*rik_hat;

  vec3_scale(dfc*gijk*ex_delr,rik_hat,drk);
  vec3_scaleadd(fc*gijk_d*ex_delr,dcosdrk,drk,drk);
  vec3_scaleadd(-fc*gijk*ex_delr_d,rik_hat,drk,drk);
  vec3_scale(prefactor,drk,drk);
}

/* ---------------------------------------------------------------------- */

void PairTersoffSunway::costheta_d(double *rij_hat, double rij,
                             double *rik_hat, double rik,
                             double *dri, double *drj, double *drk)
{
  // first element is devative wrt Ri, second wrt Rj, third wrt Rk

  double cos_theta = vec3_dot(rij_hat,rik_hat);

  vec3_scaleadd(-cos_theta,rij_hat,rik_hat,drj);
  vec3_scale(1.0/rij,drj,drj);
  vec3_scaleadd(-cos_theta,rik_hat,rij_hat,drk);
  vec3_scale(1.0/rik,drk,drk);
  vec3_add(drj,drk,dri);
  vec3_scale(-1.0,dri,dri);
}
