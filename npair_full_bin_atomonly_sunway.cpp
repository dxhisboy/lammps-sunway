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

#include "npair_full_bin_atomonly_sunway.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "atom.h"
#include "atom_vec.h"
#include "domain.h"
#include "my_page.h"
#include "error.h"
#include "memory.h"
#include "gptl.h"
#include "sunway.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */
extern "C"{
  void npair_full_bin_atomonly_sunway_build_packed(neigh_param_t *);
}
NPairFullBinAtomonlySunway::NPairFullBinAtomonlySunway(LAMMPS *lmp) : NPair(lmp) {}

/* ----------------------------------------------------------------------
   binned neighbor list construction for all neighbors
   every neighbor pair appears in list of both atoms i and j
------------------------------------------------------------------------- */

void NPairFullBinAtomonlySunway::build(NeighList *list)
{
  int i,j,k,n,itype,jtype,ibin, ii, jj;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  int *neighptr;
  GPTLstart("full_build_sunway");
  double **x = atom->x;
  int *type = atom->type;
  int *mask = atom->mask;
  tagint *molecule = atom->molecule;
  int nlocal = atom->nlocal;
  int nall = atom->nlocal + atom->nghost;
  // printf("%d %d\n", atom->nlocal, atom->nghost);
  // printf("%d %d\n", mbins, nall);

  if (includegroup) nlocal = atom->nfirst;

  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  MyPage<int> *ipage = list->ipage;

  int inum = 0;
  ipage->reset();

  //int *binpack, *binpackhead, *binpacktype, **binpackfn, *binpacknn;
  //double **binpackx;
  int *binpacknn;
  int *binpackhead;
  bin_pack_atom_t *binpack;
  memory->create(binpack, nall, "binpack");
  // memory->create(binpack, nall, "binpack");
  memory->create(binpackhead, mbins + 1, "binpackhead");
  // memory->create(binpackx, nall, 3, "binpackx");
  // memory->create(binpacktype, nall, "binpacktype");
  // binpackfn = (int**)memory->smalloc(sizeof(int*) * nall, "binpackfn");
  memory->create(binpacknn, nall, "binpacknn");

  for (i = 0; i < nall; i ++){
    binpack[i].firstneigh = ipage->vget();
    ipage->vgot(ipage->maxchunk);
    binpacknn[i] = 0;
  }
  
  int newhead = 0;
  GPTLstart("binpack");
  for (i = 0; i < mbins; i ++){
    binpackhead[i] = newhead;
    for (j = binhead[i]; j >= 0; j = bins[j]){
      binpack[newhead].id = j;
      binpack[newhead].x[0] = x[j][0];
      binpack[newhead].x[1] = x[j][1];
      binpack[newhead].x[2] = x[j][2];
      binpack[newhead].type = type[j]; 
      newhead ++;
    }
  }
  //printf("%d %d\n", newhead, nall);
  binpackhead[mbins] = newhead;
  // puts("build");
  // fflush(stdout);
  GPTLstop("binpack");

  neigh_param_t pm;
  // pm.binpack      = binpack;
  pm.binpackhead  = binpackhead;
  // pm.binpacktype  = binpacktype;
  // pm.binpackfn    = binpackfn;
  pm.binpacknn    = binpacknn;
  // pm.binpackx     = (double(*)[3])(void*)binpackx[0];
  
  pm.binpack      = binpack;
  pm.cutneighsq   = (double*)(void*)cutneighsq[0];
  pm.nstencil     = nstencil;
  pm.stencil      = stencil;
  pm.nlocal       = atom->nlocal;
  pm.nghost       = atom->nghost;
  pm.ntypes       = atom->ntypes;
  pm.mbins        = mbins;
  pm.maxchunk     = ipage->maxchunk;

  npair_full_bin_atomonly_sunway_build_packed(&pm);
  // for (ibin = 0; ibin < mbins; ibin ++){
  //   for (k = 0; k < nstencil; k ++) {
  //     if (0 <= ibin + stencil[k] && ibin + stencil[k] < mbins){
  //       int js = binpackhead[ibin + stencil[k]];
  //       int je = binpackhead[ibin + stencil[k] + 1];
  //       for (j = js; j < je; j ++){
  //         jtype = binpacktype[j];
  //         for (i = binpackhead[ibin]; i < binpackhead[ibin + 1]; i ++){
  //           if (binpack[i] < nlocal && binpack[i] != binpack[j]){
  //             itype = binpacktype[i];

  //             delx = binpackx[i][0] - binpackx[j][0];
  //             dely = binpackx[i][1] - binpackx[j][1];
  //             delz = binpackx[i][2] - binpackx[j][2];

  //             rsq = delx * delx + dely * dely + delz * delz;

  //             if (rsq <= cutneighsq[itype][jtype]){
  //               binpackfn[i][binpacknn[i] ++] = binpack[j];
  //             }
  //           }
  //         }
  //       }
  //     }
  //   }
  // }
  // puts("packfn");
  // fflush(stdout);
  GPTLstart("binunpack");
  for (i = 0; i < nall; i ++){
    if (binpack[i].id < nlocal){
      firstneigh[binpack[i].id] = binpack[i].firstneigh;
      numneigh[binpack[i].id] = binpacknn[i];
    }
  }
  GPTLstop("binunpack");
  inum = nlocal;
  for (i = 0; i < nlocal; i ++)
    ilist[i] = i;

  list->inum = inum;
  list->gnum = 0;
  memory->destroy(binpack);
  memory->destroy(binpackhead);
  // memory->destroy(binpackx);
  // memory->destroy(binpacktype);
  // memory->sfree(binpackfn);
  memory->destroy(binpacknn);
  GPTLstop("full_build_sunway");
}

// void NPairFullBinAtomonlySunway::build(NeighList *list)
// {
//   int i,j,jj,k,n,itype,jtype,ibin, ii, ioff;
//   double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
//   int *neighptr;
//   GPTLstart("full_build_sunway");
//   double **x = atom->x;
//   int *type = atom->type;
//   int *mask = atom->mask;
//   tagint *molecule = atom->molecule;
//   int nlocal = atom->nlocal;
//   if (includegroup) nlocal = atom->nfirst;

//   int *ilist = list->ilist;
//   int *numneigh = list->numneigh;
//   int **firstneigh = list->firstneigh;
//   MyPage<int> *ipage = list->ipage;
//   int nall = atom->nlocal + atom->nghost;
//   int *binpack, *binpackhead, *binpacktype;
//   double **binpackx;
//   int *neighptr_para[64], neighoffset[64][NEIGH_ISTEP + 32], iilast[64], iistart[64];
//   memory->create(binpack, nall, "binpack");
//   memory->create(binpackhead, mbins + 1, "binpackhead");
//   memory->create(binpackx, nall, 3, "binpackx");
//   memory->create(binpacktype, nall, "binpacktype");
  
//   int newhead = 0;
//   GPTLstart("binpack");
//   for (i = 0; i < mbins; i ++){
//     binpackhead[i] = newhead;
//     for (j = binhead[i]; j >= 0; j = bins[j]){
//       binpack[newhead] = j;
//       binpackx[newhead][0] = x[j][0];
//       binpackx[newhead][1] = x[j][1];
//       binpackx[newhead][2] = x[j][2];
//       binpacktype[newhead] = type[j]; 
//       newhead ++;
//     }
//   }
//   binpackhead[mbins] = newhead;
//   puts("build");
//   GPTLstop("binpack");
//   GPTLstart("build");

//   neigh_param_t pm;
//   pm.binhead      = binhead;
//   pm.binpack      = binpack;
//   pm.binpackhead  = binpackhead;
//   pm.binpacktype  = binpacktype;
//   pm.binpackx     = (double(*)[3])(void*)binpackx[0];
//   pm.cutneighsq   = (double*)(void*)cutneighsq[0];
//   pm.nstencil     = nstencil;
//   pm.stencil      = stencil;
//   pm.neighptr     = neighptr_para;
//   pm.neighoffset  = neighoffset;
//   pm.iilast       = iilast;
//   pm.nlocal       = atom->nlocal;
//   pm.nghost       = atom->nghost;
//   //pm.inum         = 0;
//   pm.x            = (double(*)[3])(void*)x[0];
//   pm.type         = type;
//   pm.ntypes       = atom->ntypes;
//   pm.maxchunk     = ipage->maxchunk;
//   pm.atom2bin     = atom2bin;
//   for (i = 0; i < nlocal; i ++){
//     neighptr = ipage->vget();
//     firstneigh[i] = neighptr;
//     ipage->vgot(ipage->maxchunk);
//   }
//   for (ibin = 0; ibin < mbins; ibin ++){
//     for (ii = 0; ii < binpackhead; ii ++){
      
//     }
//   }
//   //inum = nlocal;
//   // int inum = 0;
//   // ipage->reset();
//   // for (i = 0; i < nlocal; i++) {
//   //   n = 0;
//   //   neighptr = ipage->vget();

//   //   itype = type[i];
//   //   xtmp = x[i][0];
//   //   ytmp = x[i][1];
//   //   ztmp = x[i][2];

//   //   // loop over all atoms in surrounding bins in stencil including self
//   //   // skip i = j

//   //   ibin = atom2bin[i];
//   //   //GPTLstart("stencil scan");    
//   //   for (k = 0; k < nstencil; k++) {

//   //     //for (j = binhead[ibin+stencil[k]]; j >= 0; j = bins[j]) {
//   //     int kbin = ibin + stencil[k];
//   //     //printf("%d %d\n", binpackhead[kbin], binpackhead[kbin + 1]);
//   //     for (jj = binpackhead[kbin]; jj < binpackhead[kbin + 1]; jj ++) {
//   //       j = binpack[jj];
//   //       if (i == j) continue;

//   //       //jtype = type[j];
//   //       jtype = binpacktype[jj];
//   //       if (exclude && exclusion(i,j,itype,jtype,mask,molecule)) continue;

//   //       // delx = xtmp - x[j][0];
//   //       // dely = ytmp - x[j][1];
//   //       // delz = ztmp - x[j][2];

//   //       delx = xtmp - binpackx[jj][0];
//   //       dely = ytmp - binpackx[jj][1];
//   //       delz = ztmp - binpackx[jj][2];

//   //       rsq = delx*delx + dely*dely + delz*delz;

//   //       if (rsq <= cutneighsq[itype][jtype]) neighptr[n++] = j;
//   //     }
//   //   }
//   //   //GPTLstop("stencil scan");
//   //   ilist[inum++] = i;
//   //   firstneigh[i] = neighptr;
//   //   numneigh[i] = n;
//   //   ipage->vgot(n);
//   //   if (ipage->status())
//   //     error->one(FLERR,"Neighbor list overflow, boost neigh_modify one");
//   // }
//   GPTLstop("build");
//   list->inum = atom->nlocal;
//   list->gnum = 0;
//   memory->destroy(binpack);
//   memory->destroy(binpackhead);
//   memory->destroy(binpackx);
//   memory->destroy(binpacktype);
//   puts("build done");
//   fflush(stdout);
//   GPTLstop("full_build_sunway");
// }
