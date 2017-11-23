#include "sunway.h"
#define BNDL_SIZE 16
#ifdef MPE
extern SLAVE_FUN(npair_full_bin_atomonly_sunway_build_packed_para)(neigh_param_t *pm);

void npair_full_bin_atomonly_sunway_build_packed(neigh_param_t *pm){
  if (athread_idle() == 0)
    athread_init();
  athread_spawn(npair_full_bin_atomonly_sunway_build_packed_para, pm);
  athread_join();
  /* int ibins; */
  /* for (ibins = 0; ibins < pm->mbins; ibins += BNDL_SIZE){ */
  /*   int ibine = ibins + BNDL_SIZE; */
  /*   if (ibine > pm->mbins) */
  /*     ibine = pm->mbins; */
  /*   int ibin; */
  /*   for (ibin = ibins; ibin < ibine; ibin ++){ */
  /*     int k; */
  /*     for (k = 0; k < pm->nstencil; k ++){ */
  /*       if (ibin + pm->stencil[k] > 0 && ibin + pm->stencil[k] < pm->mbins){ */
  /*         int js = pm->binpackhead[ibin + pm->stencil[k]]; */
  /*         int je = pm->binpackhead[ibin + pm->stencil[k] + 1]; */
  /*         int j; */
  /*         for (j = js; j < je; j ++){ */
  /*           int jtype = pm->binpack[j].type; */
  /*           int i; */
  /*           for (i = pm->binpackhead[ibin]; i < pm->binpackhead[ibin + 1]; i ++){ */
  /*             if (pm->binpack[i].id < pm->nlocal && */
  /*                 pm->binpack[i].id != pm->binpack[j].id){ */
  /*               int itype = pm->binpack[i].type; */
  /*               double delx = pm->binpack[i].x[0] - pm->binpack[j].x[0]; */
  /*               double dely = pm->binpack[i].x[1] - pm->binpack[j].x[1]; */
  /*               double delz = pm->binpack[i].x[2] - pm->binpack[j].x[2]; */
            
  /*               double rsq = delx * delx + dely * dely + delz * delz; */
  /*               if (rsq <= pm->cutneighsq[itype * (pm->ntypes + 1) + jtype]){ */
  /*                 pm->binpack[i].firstneigh[pm->binpack[i].numneigh ++] = pm->binpack[j].id; */
  /*               } */
  /*             } */
  /*           } */
  /*         } */
  /*       } */
  /*     } */
  /*   } */
  /* } */
}
#endif
#ifdef CPE
#define BIN_PAGESIZE 64
#define NEIGH_PAGESIZE 128
void npair_full_bin_atomonly_sunway_build_packed_para(neigh_param_t *pm){
  pe_init();
  neigh_param_t l_pm;
  pe_get(pm, &l_pm, sizeof(neigh_param_t));
  pe_syn();

  int nstencil = l_pm.nstencil;
  int mbins = l_pm.mbins;
  int stencil[l_pm.nstencil];
  double cutneighsq[l_pm.ntypes + 1][l_pm.ntypes + 1];
  bin_pack_atom_t binpacki[BIN_PAGESIZE], binpackk[BIN_PAGESIZE];
  int binpackheadi[BNDL_SIZE + 1], binpackheadk[BNDL_SIZE + 1];
  int firstneigh[BIN_PAGESIZE][NEIGH_PAGESIZE], numneigh[BIN_PAGESIZE], totalneigh[BIN_PAGESIZE];
  pe_get(l_pm.stencil, stencil, nstencil * sizeof(int));
  pe_get(l_pm.cutneighsq, cutneighsq, sizeof(double) * (l_pm.ntypes + 1) * (l_pm.ntypes + 1));
  pe_syn();
  int ibins;
  for (ibins = BNDL_SIZE * _MYID; ibins < mbins; ibins += BNDL_SIZE * 64){
    int ibine = ibins + BNDL_SIZE;
    if (ibine > mbins)
      ibine = mbins;
    int ibin, ibinpacks, ibinpacke, ibinpackn;
    if (ibine - ibins + 1 > 0) {
      pe_get(l_pm.binpackhead + ibins, binpackheadi, (ibine - ibins + 1) * sizeof(int));
      pe_syn();
    
      ibinpacks = binpackheadi[0];
      ibinpacke = binpackheadi[ibine - ibins];
      ibinpackn = ibinpacke - ibinpacks;
      if (ibinpackn > 0){
        pe_get(l_pm.binpack + ibinpacks, binpacki, sizeof(bin_pack_atom_t) * ibinpackn);
        int i;
        for (i = 0; i < ibinpackn; i ++)
          numneigh[i] = totalneigh[i] = 0;
        pe_syn();
      }
    }
    int k;
    for (k = 0; k < l_pm.nstencil; k ++){
      int kbins = ibins + stencil[k];
      int kbine = ibine + stencil[k];
      if (kbins < 0) kbins = 0;
      if (kbine > mbins) kbine = mbins;
      if (kbins >= kbine) continue;
      int kbinpacks, kbinpacke, kbinpackn;
      if (kbine - kbins + 1 > 0) {
        pe_get(l_pm.binpackhead + kbins, binpackheadk, (kbine - kbins + 1) * sizeof(int));
        pe_syn();

        kbinpacks = binpackheadk[0];
        kbinpacke = binpackheadk[kbine - kbins];
        kbinpackn = kbinpacke - kbinpacks;
        if (kbinpackn > 0){
          pe_get(l_pm.binpack + kbinpacks, binpackk, sizeof(bin_pack_atom_t) * kbinpackn);
          pe_syn();
        }
      }

      for (ibin = ibins; ibin < ibine; ibin ++){
        int iboff = ibin - ibins;
        if (ibin + stencil[k] > 0 && ibin + stencil[k] < mbins){
          int js = binpackheadk[ibin - kbins + stencil[k]];
          int je = binpackheadk[ibin - kbins + stencil[k] + 1];
          int j;
          for (j = js; j < je; j ++){
            bin_pack_atom_t *jatom = binpackk + j - kbinpacks;
            int jtype = jatom->type;
            int i;
            for (i = binpackheadi[iboff]; i < binpackheadi[iboff + 1]; i ++){
              bin_pack_atom_t *iatom = binpacki + i - ibinpacks;
              if (iatom->id < l_pm.nlocal &&
                  iatom->id != jatom->id){
                int itype = iatom->type;
                double delx = iatom->x[0] - jatom->x[0];
                double dely = iatom->x[1] - jatom->x[1];
                double delz = iatom->x[2] - jatom->x[2];
            
                double rsq = delx * delx + dely * dely + delz * delz;
                if (rsq <= cutneighsq[itype][jtype]){
                  int ioff = i - ibinpacks;
                  firstneigh[ioff][numneigh[ioff] ++] = jatom->id;
                  if (numneigh[ioff] >= NEIGH_PAGESIZE){
                    pe_put(iatom->firstneigh + totalneigh[ioff], firstneigh[ioff], sizeof(int) * NEIGH_PAGESIZE);
                    totalneigh[ioff] += numneigh[ioff];
                    numneigh[ioff] = 0;
                    pe_syn();
                  }
                  //l_pm.binpacknn[i] ++;
                  //iatom->firstneigh[l_pm.binpacknn[i] ++] = jatom->id;
                }
              }
            }
          }
        }
      }
    }
    if (ibinpackn > 0){
      int i;
      for (i = 0; i < ibinpackn; i ++){
        if (numneigh[i] > 0){
          pe_put(binpacki[i].firstneigh + totalneigh[i], firstneigh[i], sizeof(int) * numneigh[i]);
          totalneigh[i] += numneigh[i];
          pe_syn();

        }
      }
      pe_put(l_pm.binpacknn + ibinpacks, totalneigh, sizeof(int) * ibinpackn);
      pe_syn();
    }
    
  }
}
#endif
