#include "sunway.h"
#include "pair_reaxc_sw64.h"
#include <math.h>
void write_reax_atoms_and_pack(write_atoms_param_t *pm){
  int i;
  for (i = 0; i < pm->csys->N; i ++){
    pm->csys->my_atoms[i].orig_id    = pm->tag[i];
    pm->csys->my_atoms[i].type       = pm->map[pm->type[i]];
    pm->csys->my_atoms[i].x[0]       = pm->x[i][0];
    pm->csys->my_atoms[i].x[1]       = pm->x[i][1];
    pm->csys->my_atoms[i].x[2]       = pm->x[i][2];
    pm->csys->my_atoms[i].q          = pm->q[i];
    pm->csys->my_atoms[i].num_bonds  = pm->num_bonds[i];
    pm->csys->my_atoms[i].num_hbonds = pm->num_hbonds[i];

    pm->csys->packed_atoms[i].orig_id    = pm->tag[i];
    pm->csys->packed_atoms[i].type       = pm->map[pm->type[i]];
    pm->csys->packed_atoms[i].x[0]       = pm->x[i][0];
    pm->csys->packed_atoms[i].x[1]       = pm->x[i][1];
    pm->csys->packed_atoms[i].x[2]       = pm->x[i][2];
    pm->csys->packed_atoms[i].q          = pm->q[i];
  }
}

void write_reax_lists_c(write_lists_param_t *pm){
  int numall = pm->numall;
  int *ilist = pm->ilist;
  int *numneigh = pm->numneigh;
  int **firstneigh = pm->firstneigh;
  reax_list *far_nbrs = pm->far_nbrs;
  far_neighbor_data_full *far_list = far_nbrs->select.far_nbr_list_full;
  int maxfar = pm->maxfar;
  atom_pack_t *patoms = pm->patoms;
  double nonb_cutsq = pm->nonb_cutsq;
  int ii;
  for (ii = 0; ii < numall; ii ++){
    int i = ilist[ii];
    int *jlist = firstneigh[i];
    int jnum = numneigh[i];
    int jj;
    far_nbrs->index[i] = i * maxfar;
    int knum = 0;
    far_neighbor_data_full *klist = far_list + i * maxfar;
    for (jj = 0; jj < jnum; jj ++){
      int j = jlist[jj];
      j &= NEIGHMASK;
      rvec dvec;
      dvec[0] = patoms[j].x[0] - patoms[i].x[0];
      dvec[1] = patoms[j].x[1] - patoms[i].x[1];
      dvec[2] = patoms[j].x[2] - patoms[i].x[2];
      double d_sqr = dvec[0] * dvec[0] + dvec[1] * dvec[1] + dvec[2] * dvec[2];
      if (d_sqr <= nonb_cutsq){
        far_neighbor_data_full *nbr = klist + knum;
        nbr->nbr = j;
        nbr->d = sqrt(d_sqr);
        nbr->dvec[0] = dvec[0];
        nbr->dvec[1] = dvec[1];
        nbr->dvec[2] = dvec[2];
        nbr->rel_box[0] = 0;
        nbr->rel_box[1] = 0;
        nbr->rel_box[2] = 0;
        nbr->type = patoms[j].type;
        nbr->orig_id = patoms[j].orig_id;
        nbr->q = patoms[j].q;
        knum ++;
      }
    }
    far_nbrs->end_index[i] = i * maxfar + knum;
  }
}
