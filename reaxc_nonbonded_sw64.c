#ifdef MPE
#include <mpi.h>
#endif
#ifdef CPE
#include "STUBS/mpi.h"
#endif
#include <stdio.h>
#include <math.h>
#include "reaxc_defs_sunway.h"
#include "reaxc_ctypes_sunway.h"
#include "reaxc_inlines_sw64.h"
static inline int llf(double *dvec){
  double SMALL = 0.0001;
  if (dvec[2] < -SMALL) return 1;
  if (dvec[2] < SMALL){
    if (dvec[1] < -SMALL) return 1;
    if (dvec[1] < SMALL && dvec[0] < -SMALL) return 1;
  }
  return 0;
}
static inline int urb(double *dvec){
  double SMALL = 0.0001;
  if (dvec[2] > SMALL) return 1;
  if (dvec[2] > -SMALL){
    if (dvec[1] > SMALL) return 1;
    if (dvec[1] > -SMALL && dvec[0] > SMALL) return 1;
  }
  return 0;
}

typedef struct tbp_boost_t{
  double powgi_vdW1;
  double r_vdWinv;
  double rcoreinv;
  double alpha_div_rvdW;
  double a_div_rcore;
  double re6;
  
} tbp_boost_t;
void vdW_Coulomb_Energy_Full_C( reax_system_c *system, control_params *control,
                                simulation_data *data, storage *workspace,
                                reax_list **lists, output_controls *out_control )
{
  int i, j, pj, natoms, nlocal;
  int start_i, end_i, flag;
  rc_tagint orig_i, orig_j;
  double p_vdW1, p_vdW1i;
  double powr_vdW1, powgi_vdW1;
  double tmp, r_ij, fn13, exp1, exp2;
  double Tap, dTap, dfn13, CEvd, CEclmb, de_core;
  double dr3gamij_1, dr3gamij_3;
  double e_ele, e_vdW, e_core, SMALL = 0.0001;
  double e_lg, de_lg, r_ij5, r_ij6, re6;
  rvec temp, ext_press;
  two_body_parameters *twbp;
  far_neighbor_data_full *nbr_pj;
  reax_list *far_nbrs;

  // Tallying variables:
  double pe_vdw, f_tmp, delij[3];

  natoms = system->N;
  nlocal = system->n;
  far_nbrs = (*lists) + FAR_NBRS_FULL;
  p_vdW1 = system->reax_param.gp.l[28];
  p_vdW1i = 1.0 / p_vdW1;
  e_core = 0;
  e_vdW = 0;
  e_lg = de_lg = 0.0;

  int ntypes = system->reax_param.num_atom_types;
  tbp_boost_t tbpb[ntypes][ntypes];
  tbp_boost_t *tbpt;
  two_body_parameters **tbpr = system->reax_param.tbp;
  /* printf("%d\n", ntypes); */
  for (i = 0; i < ntypes; i ++)
    for (j = 0; j < ntypes; j ++){
      tbpb[i][j].powgi_vdW1 = pow(tbpr[i][j].gamma_w, -p_vdW1);
      tbpb[i][j].r_vdWinv = 1.0 / tbpr[i][j].r_vdW;
      tbpb[i][j].alpha_div_rvdW = tbpr[i][j].alpha / tbpr[i][j].r_vdW;
      if(system->reax_param.gp.vdw_type==2 || system->reax_param.gp.vdw_type==3){
        tbpb[i][j].a_div_rcore = tbpr[i][j].acore / tbpr[i][j].rcore;
        tbpb[i][j].rcoreinv = 1.0 / tbpr[i][j].rcore;
      }
      double lgresq = tbpr[i][j].lgre * tbpr[i][j].lgre;
      tbpb[i][j].re6 = lgresq * lgresq * lgresq;
    }
  //system->reax_param.gp.vdw_type = 1;
  for( i = 0; i < natoms; ++i ) {
    if (system->my_atoms[i].type < 0) continue;
    start_i = Start_Index(i, far_nbrs);
    end_i   = End_Index(i, far_nbrs);
    orig_i  = system->my_atoms[i].orig_id;

    for( pj = start_i; pj < end_i; ++pj ) {
      nbr_pj = &(far_nbrs->select.far_nbr_list_full[pj]);
      j = nbr_pj->nbr;
      if (nbr_pj->type < 0) continue;
      orig_j  = nbr_pj->orig_id;
      flag = 1;

      if (i >= nlocal && j >= nlocal) continue;
      if (i >= nlocal){
        if (orig_j > orig_i) continue;
        if (orig_j == orig_i && urb(nbr_pj->dvec)) continue;
      }
      if (j >= nlocal){
        if (orig_i > orig_j) continue;
        if (orig_i == orig_j && llf(nbr_pj->dvec)) continue;
      }
        

      r_ij = nbr_pj->d;
      double rijinv = 1 / r_ij;
      double r2ij = r_ij * r_ij;
      twbp = &(system->reax_param.tbp[ system->my_atoms[i].type ]
               [ nbr_pj->type ]);
      tbpt = tbpb[system->my_atoms[i].type] + nbr_pj->type;
      Tap = workspace->Tap[7] * r_ij + workspace->Tap[6];
      Tap = Tap * r_ij + workspace->Tap[5];
      Tap = Tap * r_ij + workspace->Tap[4];
      Tap = Tap * r_ij + workspace->Tap[3];
      Tap = Tap * r_ij + workspace->Tap[2];
      Tap = Tap * r_ij + workspace->Tap[1];
      Tap = Tap * r_ij + workspace->Tap[0];

      dTap = 7*workspace->Tap[7] * r_ij + 6*workspace->Tap[6];
      dTap = dTap * r_ij + 5*workspace->Tap[5];
      dTap = dTap * r_ij + 4*workspace->Tap[4];
      dTap = dTap * r_ij + 3*workspace->Tap[3];
      dTap = dTap * r_ij + 2*workspace->Tap[2];
      dTap += workspace->Tap[1] * rijinv;

      /*vdWaals Calculations*/
      if(system->reax_param.gp.vdw_type==1 || system->reax_param.gp.vdw_type==3)
        { // shielding
          powr_vdW1 = pow(r_ij, p_vdW1);
          powgi_vdW1 = tbpt->powgi_vdW1;
          //pow( twbp->gamma_w, -p_vdW1);

          fn13 = pow( powr_vdW1 + powgi_vdW1, p_vdW1i );
          //exp1 = exp( twbp->alpha * (1.0 - fn13 * tbpt->r_vdWinv) );
          exp2 = exp( 0.5 * twbp->alpha * (1.0 - fn13 * tbpt->r_vdWinv) );
          exp1 = exp2 * exp2;
          e_vdW = twbp->D * (exp1 - 2.0 * exp2);
          data->my_en.e_vdW += Tap * e_vdW;

          dfn13 = pow( powr_vdW1 + powgi_vdW1, p_vdW1i - 1.0) *
            pow(r_ij, p_vdW1 - 2.0);

          CEvd = dTap * e_vdW -
            Tap * twbp->D * (tbpt->alpha_div_rvdW) * (exp1 - exp2) * dfn13;
        }
      else{ // no shielding
        exp2 = exp( 0.5 * twbp->alpha * (1.0 - r_ij * tbpt->r_vdWinv) );
        exp1 = exp2 * exp2;
        e_vdW = twbp->D * (exp1 - 2.0 * exp2);
        data->my_en.e_vdW += Tap * e_vdW;

        CEvd = dTap * e_vdW -
          Tap * twbp->D * (tbpt->alpha_div_rvdW) * (exp1 - exp2) * rijinv;
      }

      if(system->reax_param.gp.vdw_type==2 || system->reax_param.gp.vdw_type==3)
        { // inner wall
          e_core = twbp->ecore * exp(twbp->acore * (1.0-(r_ij*tbpt->rcoreinv)));
          data->my_en.e_vdW += Tap * e_core;

          de_core = -(tbpt->a_div_rcore) * e_core;
          CEvd += dTap * e_core + Tap * de_core  * rijinv;

          //  lg correction, only if lgvdw is yes
          if (control->lgflag) {
            r_ij5 = r2ij * r2ij * r_ij;
            r_ij6 = r2ij * r2ij * r2ij;
            re6 = tbpt->re6;
            double rrinv = 1 / (r_ij6 + re6);
            e_lg = -(twbp->lgcij * rrinv);
            data->my_en.e_vdW += Tap * e_lg;

            de_lg = -6.0 * e_lg *  r_ij5 * rrinv;
            CEvd += dTap * e_lg + Tap * de_lg  * rijinv;
          }

        }

      /*Coulomb Calculations*/
      dr3gamij_1 = ( r_ij * r_ij * r_ij + twbp->gamma );
      dr3gamij_3 = pow( dr3gamij_1 , 0.33333333333333 );
      double dr3gamij_1inv = 1 / dr3gamij_1;
      double dr3gamij_3inv = 1 / dr3gamij_3;
      tmp = Tap * dr3gamij_3inv;
      data->my_en.e_ele += e_ele =
        C_ele * system->my_atoms[i].q * nbr_pj->q * tmp;

      CEclmb = C_ele * system->my_atoms[i].q * nbr_pj->q *
        ( dTap -  Tap * r_ij * dr3gamij_1inv ) * dr3gamij_3inv;

      /* tally into per-atom energy */
      if( system->evflag || system->vflag_atom) {
        pe_vdw = Tap * (e_vdW + e_core + e_lg);
        /* rvec_ScaledSum( delij, 1., system->my_atoms[i].x, */
        /*                 -1., system->my_atoms[j].x ); */
        f_tmp = -(CEvd + CEclmb);
        /* system->pair_ptr->ev_tally_full(i,pe_vdw,e_ele,f_tmp,delij[0], */
        /*                                 delij[1],delij[2]); */
        //ev_tally_full_sys(i, pe_vdw, e_ele, f_tmp, delij, system);
        ev_tally_full_sys(i, pe_vdw, e_ele, f_tmp, nbr_pj->dvec, system);
      }

      if( control->virial == 0 ) {
        rvec_ScaledAdd( workspace->f[i], -(CEvd + CEclmb), nbr_pj->dvec );
      }
      else { /* NPT, iNPT or sNPT */
        rvec_Scale( temp, CEvd + CEclmb, nbr_pj->dvec );

        rvec_ScaledAdd( workspace->f[i], -1., temp );
        rvec_iMultiply( ext_press, nbr_pj->rel_box, temp );
        rvec_ScaledAdd( data->my_ext_press, 0.5, ext_press );
      }
    }
  }    
}
