#include "reaxc_types_sunway.h"
using namespace REAXC_SUNWAY_NS;
void reax_system::to_c_sys(reax_system_c *csys){
  //puts("begin");
  csys->reax_param    = reax_param;
  csys->bigN          = bigN;
  csys->n             = n;
  csys->N             = N;
  csys->numH          = numH;
  csys->local_cap     = local_cap;
  csys->total_cap     = total_cap;
  csys->gcell_cap     = gcell_cap;
  csys->Hcap          = Hcap;
  csys->est_recv      = est_recv;
  csys->est_trans     = est_trans;
  csys->max_recved    = max_recved;
  csys->wsize         = wsize;
  csys->my_rank       = my_rank;
  csys->num_nbrs      = num_nbrs;
  //puts("coord");
  csys->my_coords[0]  = my_coords[0];
  csys->my_coords[1]  = my_coords[1];
  csys->my_coords[2]  = my_coords[2];
  csys->global_offset = global_offset;
  csys->big_box       = big_box;
  csys->my_box        = my_box;
  csys->my_ext_box    = my_ext_box;
  csys->my_grid       = my_grid;
  csys->bndry_cuts    = bndry_cuts;
  csys->my_atoms      = my_atoms;
  csys->packed_atoms  = packed_atoms;
  csys->my_bonds      = my_bonds;
  csys->mincap        = mincap;
  csys->maxfar        = maxfar;
  csys->safezone      = safezone;
  csys->saferzone     = saferzone;
  //puts("xve");
  if (pair_ptr->evflag) {
    csys->x             = (double(*)[3])x[0];//(void*)pair_ptr->atom->x[0];
    //puts("v");
    if (pair_ptr->vflag_atom)
      csys->vatom         = (double(*)[6])(void*)(pair_ptr->vatom[0]);
    //puts("e");
    if (pair_ptr->eflag_atom)
    csys->eatom         = pair_ptr->eatom;
  }
  //puts("flgs");
  csys->eng_vdwl      = pair_ptr->eng_vdwl;
  csys->eng_coul      = pair_ptr->eng_coul;
  csys->eflag_atom    = pair_ptr->eflag_atom;
  csys->vflag_atom    = pair_ptr->vflag_atom;
  csys->eflag_global  = pair_ptr->eflag_global;
  csys->vflag_global  = pair_ptr->vflag_global;
  csys->evflag        = pair_ptr->evflag;
  //puts("cpy");
  memcpy(csys->virial , pair_ptr->virial , sizeof(double) * 6);
  memcpy(csys->my_nbrs, my_nbrs, sizeof(neighbor_proc) * REAX_MAX_NBRS);
}

void reax_system::from_c_sys(reax_system_c *csys){
  reax_param    = csys->reax_param;
  bigN          = csys->bigN;
  n             = csys->n;
  N             = csys->N;
  numH          = csys->numH;
  local_cap     = csys->local_cap;
  total_cap     = csys->total_cap;
  gcell_cap     = csys->gcell_cap;
  Hcap          = csys->Hcap;
  est_recv      = csys->est_recv;
  est_trans     = csys->est_trans;
  max_recved    = csys->max_recved;
  wsize         = csys->wsize;
  my_rank       = csys->my_rank;
  num_nbrs      = csys->num_nbrs;
  my_coords[0]  = csys->my_coords[0];
  my_coords[1]  = csys->my_coords[1];
  my_coords[2]  = csys->my_coords[2];
  global_offset = csys->global_offset;
  big_box       = csys->big_box;
  my_box        = csys->my_box;
  my_ext_box    = csys->my_ext_box;
  my_grid       = csys->my_grid;
  bndry_cuts    = csys->bndry_cuts;
  my_atoms      = csys->my_atoms;
  packed_atoms  = csys->packed_atoms;
  my_bonds      = csys->my_bonds;
  mincap        = csys->mincap;
  maxfar        = csys->maxfar;
  safezone      = csys->safezone;
  saferzone     = csys->saferzone;
  pair_ptr->eng_vdwl = csys->eng_vdwl;
  pair_ptr->eng_coul = csys->eng_coul;

  memcpy(pair_ptr->virial , csys->virial , sizeof(double) * 6);
  memcpy(my_nbrs, csys->my_nbrs, sizeof(neighbor_proc) * REAX_MAX_NBRS);
}

