#include "sunway.h"
#ifdef MPE
extern SLAVE_FUN(npair_full_bin_atomonly_sunway_scan_para)(neigh_param_t *pm);
void npair_full_bin_atomonly_sunway_scan(neigh_param_t *pm){
  if (athread_idle() == 0)
    athread_init();

  athread_spawn(npair_full_bin_atomonly_sunway_scan_para, pm);
  athread_join();
  /* int _MYID; */
  /* for (_MYID = 0; _MYID < 64; _MYID ++){ */
  /*   int i, ii, j, k, jj, itype, jtype; */
  /*   int ibin, kbin; */
  /*   double xi[3]; */
  /*   double delx, dely, delz, rsq; */
  /*   int chunkfull = 0; */
  /*   int neighcnt = 0; */
  /*   int listend = 0; */
  /*   for (ii = pm->iilast[_MYID]; ii < NEIGH_ISTEP; ii ++){ */
  /*     pm->neighoffset[_MYID][ii] = neighcnt; */
  /*     i = _MYID * NEIGH_ISTEP + pm->istart + ii; */
  /*     if (i >= pm->nlocal){ */
  /*       pm->iilast[_MYID] = NEIGH_ISTEP; */
  /*       listend = 1; */
  /*       break; */
  /*     } */
  /*     xi[0] = pm->x[i][0]; */
  /*     xi[1] = pm->x[i][1]; */
  /*     xi[2] = pm->x[i][2]; */
  /*     itype = pm->type[i]; */
  /*     ibin = pm->atom2bin[i]; */
  /*     for (k = 0; k < pm->nstencil; k ++){ */
  /*       kbin = ibin + pm->stencil[k]; */
  /*       //printf("%d %d %d\n", kbin, pm->binpackhead[kbin], pm->binpackhead[kbin + 1]); */
  /*       for (jj = pm->binpackhead[kbin]; jj < pm->binpackhead[kbin + 1]; jj ++){ */
  /*         j = pm->binpack[jj]; */
  /*         if (i == j) continue; */
  /*         jtype = pm->binpacktype[jj]; */
  /*         delx = xi[0] - pm->binpackx[jj][0]; */
  /*         dely = xi[1] - pm->binpackx[jj][1]; */
  /*         delz = xi[2] - pm->binpackx[jj][2]; */
  /*         rsq = delx * delx + dely * dely + delz * delz; */
  /*         if (rsq <= pm->cutneighsq[itype * (pm->ntypes + 1) + jtype]){ */
  /*           //printf("%d %d %d\n", ii, pm->maxchunk, neighcnt); */
  /*           if (neighcnt < pm->maxchunk){ */
  /*             pm->neighptr[_MYID][neighcnt++] = j; */
  /*           } else */
  /*             chunkfull = 1; */
  /*         } */
  /*       } */
  /*       if (chunkfull) */
  /*         break; */
  /*     } */
  /*     if (chunkfull) */
  /*       break; */
  /*   } */
  /*   if (ii == NEIGH_ISTEP) */
  /*     pm->neighoffset[_MYID][ii] = neighcnt; */

  /*   if (listend) */
  /*     continue; */
  /*   pm->iilast[_MYID] = ii; */
  /* } */
}
#endif
#ifdef CPE
#define JPAGE_SIZE 64
#define MAXT2 16
#define NEIGH_BUFSIZE 2048
#define MAX_STENCIL 512
void npair_full_bin_atomonly_sunway_scan_para(neigh_param_t *pm){
  pe_init();
  int i, ii, j, k, jj, itype, jtype;
  int ibin, kbin;
  //double xi[3];
  double delx, dely, delz, rsq;
  int chunkfull = 0;
  int neighcnt = 0;
  int listend = 0;
  int iilast = pm->iilast[_MYID];

  neigh_param_t l_pm;
  pe_get(pm, &l_pm, sizeof(neigh_param_t));
  pe_syn();

  int neighoffset[NEIGH_ISTEP + 32];

  double xi[NEIGH_ISTEP][3];
  double cutneighsq[MAXT2];
  int ti[NEIGH_ISTEP], a2b[NEIGH_ISTEP];
  int ibase = _MYID * NEIGH_ISTEP + l_pm.istart;
  int stencil[MAX_STENCIL];
  pe_get(l_pm.x[ibase], xi[0], sizeof(double) * 3 * NEIGH_ISTEP);
  pe_get(l_pm.type + ibase, ti, sizeof(int) * NEIGH_ISTEP);
  pe_get(l_pm.atom2bin + ibase, a2b, sizeof(int) * NEIGH_ISTEP);
  pe_get(l_pm.cutneighsq, cutneighsq, sizeof(double) * (l_pm.ntypes + 1) * (l_pm.ntypes + 1));
  pe_get(l_pm.stencil, stencil, sizeof(int) * l_pm.nstencil);
  pe_syn();


  double xj[JPAGE_SIZE][3];
  int tj[JPAGE_SIZE];
  int jlist[JPAGE_SIZE];
  int jpage_start, jpage_size;
  int neighbuf[NEIGH_BUFSIZE], neighbufoff = 0, neighptroff = 0;
  int *neighptr = l_pm.neighptr[_MYID];
  for (ii = iilast; ii < NEIGH_ISTEP; ii ++){
    neighoffset[ii] = neighcnt;
    i = _MYID * NEIGH_ISTEP + l_pm.istart + ii;
    if (i >= l_pm.nlocal){
      l_pm.iilast[_MYID] = NEIGH_ISTEP;
      listend = 1;
      break;
    }
    itype = ti[ii];
    ibin = a2b[ii];
    for (k = 0; k < l_pm.nstencil; k ++){
      kbin = ibin + stencil[k];
      int jr[2];
      pe_get(l_pm.binpackhead + kbin, jr, sizeof(int) * 2);
      pe_syn();
      /* int js = l_pm.binpackhead[kbin]; */
      /* int je = l_pm.binpackhead[kbin + 1]; */
      for (jpage_start = jr[0]; jpage_start < jr[1]; jpage_start += JPAGE_SIZE){
        if (jpage_start + JPAGE_SIZE <= jr[1])
          jpage_size = JPAGE_SIZE;
        else
          jpage_size = jr[1] - jpage_start;
        //for (jj = l_pm.binpackhead[kbin]; jj < l_pm.binpackhead[kbin + 1]; jj ++){
        if (jpage_size > 0){
          pe_get(l_pm.binpack + jpage_start, jlist, jpage_size * sizeof(int));
          pe_get(l_pm.binpackx + jpage_start, xj[0], jpage_size * sizeof(double) * 3);
          pe_get(l_pm.binpacktype + jpage_start, tj, jpage_size * sizeof(int));
          pe_syn();
        }
        for (jj = 0; jj < jpage_size; jj ++){
          //j = l_pm.binpack[jj];
          j = jlist[jj];
          if (i == j) continue;
          //jtype = l_pm.binpacktype[jj];
          jtype = tj[jj];
          /* delx = xi[ii][0] - l_pm.binpackx[jj][0]; */
          /* dely = xi[ii][1] - l_pm.binpackx[jj][1]; */
          /* delz = xi[ii][2] - l_pm.binpackx[jj][2]; */
          delx = xi[ii][0] - xj[jj][0];
          dely = xi[ii][1] - xj[jj][1];
          delz = xi[ii][2] - xj[jj][2];

          rsq = delx * delx + dely * dely + delz * delz;
          if (rsq <= cutneighsq[itype * (l_pm.ntypes + 1) + jtype]){
            //printf("%d %d %d\n", ii, l_pm.maxchunk, neighcnt);
            if (neighcnt < l_pm.maxchunk){
              neighbuf[neighbufoff ++] = j;
              neighcnt ++;
              if (neighbufoff == NEIGH_BUFSIZE){
                pe_put(neighptr + neighptroff, neighbuf, sizeof(int) * neighbufoff);
                pe_syn();
                neighptroff += NEIGH_BUFSIZE;
                neighbufoff = 0;
              }
              //l_pm.neighptr[_MYID][neighcnt++] = j;
            } else
              chunkfull = 1;
          }
        }
        if (chunkfull)
          break;
      }
      if (chunkfull)
        break;
    }
    if (chunkfull)
      break;
  }
  if (neighbufoff > 0 && !chunkfull){
    pe_put(neighptr + neighptroff, neighbuf, sizeof(int) * neighbufoff);
    pe_syn();
    //neighptroff += NEIGH_BUFSIZE;
    //neighbufoff = 0;
  }
  if (ii == NEIGH_ISTEP)/* { */
    /*   pe_put(l_pm.neighoffset[_MYID] + ii, &neighcnt, sizeof(int)); */
    /*   pe_syn(); */
    /* } */
    //l_pm.neighoffset[_MYID][ii] = neighcnt;
    neighoffset[ii] = neighcnt;
  pe_put(l_pm.neighoffset[_MYID], neighoffset, sizeof(int) * (NEIGH_ISTEP + 32));
  pe_syn();
  if (!listend)
    l_pm.iilast[_MYID] = ii;
}

#endif
