#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pair_lj_cut_full.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "update.h"
#include "integrate.h"
#include "respa.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace MathConst;

PairLJCutFull::PairLJCutFull(LAMMPS *lmp) : PairLJCut(lmp)
{
  respa_enable = 0;
}

PairLJCutFull::~PairLJCutFull()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
    memory->destroy(epsilon);
    memory->destroy(sigma);
    memory->destroy(lj1);
    memory->destroy(lj2);
    memory->destroy(lj3);
    memory->destroy(lj4);
    memory->destroy(offset);
  }
}
void PairLJCutFull::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,r2inv,r6inv,forcelj,factor_lj;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;

  if (newton_pair)
    inum = list->inum + list->gnum;
  else
    inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;
      if (newton_pair){
	if (i >= nlocal){
	  if (j >= nlocal) continue;
	  if (x[j][2] > ztmp) continue;
	  if (x[j][2] == ztmp){
	    if (x[j][1] > ytmp) continue;
	    if (x[j][1] == ytmp && x[j][0] > xtmp) continue;
	  }
	}
	if (j >= nlocal){
	  if (i >= nlocal) continue;
	  if (x[j][2] < ztmp) continue;
	  if (x[j][2] == ztmp){
	    if (x[j][1] < ytmp) continue;
	    if (x[j][1] == ytmp && x[j][0] < xtmp) continue;
	  }
	}
      }
      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
	r2inv = 1.0/rsq;
	r6inv = r2inv*r2inv*r2inv;
	forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
	fpair = factor_lj*forcelj*r2inv;

	f[i][0] += delx*fpair;
	f[i][1] += dely*fpair;
	f[i][2] += delz*fpair;

	if (eflag) {
	  evdwl = r6inv*(lj3[itype][jtype]*r6inv-lj4[itype][jtype]) -
	    offset[itype][jtype];
	  evdwl *= factor_lj;
	}

	if (evflag) ev_tally_full(i,evdwl,0.0,fpair,delx,dely,delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

void PairLJCutFull::init_style()
{
  // request regular or rRESPA neighbor lists

  int irequest;

  irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->id = 0;
  neighbor->requests[irequest]->newton = 2;
  neighbor->requests[irequest]->ghost = 1;
  neighbor->requests[irequest]->full = 1;
  neighbor->requests[irequest]->half = 0;
}

  /* ----------------------------------------------------------------------
     neighbor callback to inform pair style of neighbor list to use
     regular or rRESPA
     ------------------------------------------------------------------------- */

void PairLJCutFull::init_list(int id, NeighList *ptr)
{
  if (id == 0) list = ptr;
}

