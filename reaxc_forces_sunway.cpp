/*----------------------------------------------------------------------
  PuReMD - Purdue ReaxFF Molecular Dynamics Program

  Copyright (2010) Purdue University
  Hasan Metin Aktulga, hmaktulga@lbl.gov
  Joseph Fogarty, jcfogart@mail.usf.edu
  Sagar Pandit, pandit@usf.edu
  Ananth Y Grama, ayg@cs.purdue.edu

  Please cite the related publication:
  H. M. Aktulga, J. C. Fogarty, S. A. Pandit, A. Y. Grama,
  "Parallel Reactive Molecular Dynamics: Numerical Methods and
  Algorithmic Techniques", Parallel Computing, in press.

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License as
  published by the Free Software Foundation; either version 2 of
  the License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  See the GNU General Public License for more details:
  <http://www.gnu.org/licenses/>.
  ----------------------------------------------------------------------*/

#include "pair_reaxc_sunway.h"
#include "reaxc_forces_sunway.h"
#include "reaxc_bond_orders_sunway.h"
#include "reaxc_bonds_sunway.h"
#include "reaxc_hydrogen_bonds_sunway.h"
#include "reaxc_io_tools_sunway.h"
#include "reaxc_list_sunway.h"
#include "reaxc_lookup_sunway.h"
#include "reaxc_multi_body_sunway.h"
#include "reaxc_nonbonded_sunway.h"
#include "reaxc_tool_box_sunway.h"
#include "reaxc_torsion_angles_sunway.h"
#include "reaxc_valence_angles_sunway.h"
#include "reaxc_vector_sunway.h"

#include "gptl.h"
namespace REAXC_SUNWAY_NS{


  interaction_function Interaction_Functions[NUM_INTRS];

  void Dummy_Interaction( reax_system *system, control_params *control,
                          simulation_data *data, storage *workspace,
                          reax_list **lists, output_controls *out_control )
  {
  }


  void Init_Force_Functions( control_params *control )
  {
    Interaction_Functions[0] = BO;
    Interaction_Functions[1] = Bonds; //Dummy_Interaction;
    Interaction_Functions[2] = Atom_Energy; //Dummy_Interaction;
    Interaction_Functions[3] = Valence_Angles; //Dummy_Interaction;
    Interaction_Functions[4] = Torsion_Angles; //Dummy_Interaction;
    if( control->hbond_cut > 0 )
      Interaction_Functions[5] = Hydrogen_Bonds;
    else Interaction_Functions[5] = Dummy_Interaction;
    Interaction_Functions[6] = Dummy_Interaction; //empty
    Interaction_Functions[7] = Dummy_Interaction; //empty
    Interaction_Functions[8] = Dummy_Interaction; //empty
    Interaction_Functions[9] = Dummy_Interaction; //empty
  }


  void Compute_Bonded_Forces( reax_system *system, control_params *control,
                              simulation_data *data, storage *workspace,
                              reax_list **lists, output_controls *out_control,
                              MPI_Comm comm )
  {
    int i;

    /* Implement all force calls as function pointers */
    for( i = 0; i < NUM_INTRS; i++ ) {
      char s[200];
      sprintf(s, "reaxc bonded interaction %d", i);
      GPTLstart(s);
      (Interaction_Functions[i])( system, control, data, workspace,
                                  lists, out_control );
      GPTLstop(s);
    }
  }


  void Compute_NonBonded_Forces( reax_system *system, control_params *control,
                                 simulation_data *data, storage *workspace,
                                 reax_list **lists, output_controls *out_control,
                                 MPI_Comm comm )
  {
    /* van der Waals and Coulomb interactions */
    vdW_Coulomb_Energy_Full( system, control, data, workspace,
                             lists, out_control );
    // if( control->tabulate == 0 )
    // vdW_Coulomb_Energy( system, control, data, workspace,
    //                     lists, out_control );
    // else
    //   Tabulated_vdW_Coulomb_Energy( system, control, data, workspace,
    //                                 lists, out_control );
  }


  void Compute_Total_Force( reax_system *system, control_params *control,
                            simulation_data *data, storage *workspace,
                            reax_list **lists, mpi_datatypes *mpi_data )
  {
    int i, pj;
    reax_list *bonds = (*lists) + BONDS;

    for( i = 0; i < system->N; ++i )
      for( pj = Start_Index(i, bonds); pj < End_Index(i, bonds); ++pj )
        if( i < bonds->select.bond_list[pj].nbr ) {
          if( control->virial == 0 )
            Add_dBond_to_Forces( system, i, pj, workspace, lists );
          else
            Add_dBond_to_Forces_NPT( i, pj, data, workspace, lists );
        }

  }

  void Validate_Lists( reax_system *system, storage *workspace, reax_list **lists,
                       int step, int n, int N, int numH, MPI_Comm comm )
  {
    int i, comp, Hindex;
    reax_list *bonds, *hbonds;

    double saferzone = system->saferzone;

    /* bond list */
    if( N > 0 ) {
      bonds = *lists + BONDS;

      for( i = 0; i < N; ++i ) {
        system->my_atoms[i].num_bonds = MAX(Num_Entries(i,bonds)*2, MIN_BONDS);

        if( i < N-1 )
          comp = Start_Index(i+1, bonds);
        else comp = bonds->num_intrs;

        if( End_Index(i, bonds) > comp ) {
          fprintf( stderr, "step%d-bondchk failed: i=%d end(i)=%d str(i+1)=%d\n",
                   step, i, End_Index(i,bonds), comp );
          MPI_Abort( comm, INSUFFICIENT_MEMORY );
        }
      }
    }


    /* hbonds list */
    if( numH > 0 ) {
      hbonds = *lists + HBONDS;

      for( i = 0; i < N; ++i ) {
        Hindex = system->my_atoms[i].Hindex;
        if( Hindex > -1 ) {
          system->my_atoms[i].num_hbonds =
            (int)(MAX( Num_Entries(Hindex, hbonds)*saferzone, MIN_HBONDS ));

          //if( Num_Entries(i, hbonds) >=
          //(Start_Index(i+1,hbonds)-Start_Index(i,hbonds))*0.90/*DANGER_ZONE*/){
          //  workspace->realloc.hbonds = 1;

          if( Hindex < numH-1 )
            comp = Start_Index(Hindex+1, hbonds);
          else comp = hbonds->num_intrs;

          if( End_Index(Hindex, hbonds) > comp ) {
            fprintf(stderr,"step%d-hbondchk failed: H=%d end(H)=%d str(H+1)=%d\n",
                    step, Hindex, End_Index(Hindex,hbonds), comp );
            MPI_Abort( comm, INSUFFICIENT_MEMORY );
          }
        }
      }
    }
  }
  static void print_bond(int i, bond_data *data, reax_system *system){
    printf("%d %d, d: %.10f, dvec: [%.10f %.10f %.10f], ", 
           i, data->nbr, data->d, data->dvec[0],
           data->dvec[1], data->dvec[2]);
    bond_order_data *bo = &(data->bo_data);
    printf("BO: %.10f, BO_s: %.10f, BO_pi: %.10f, BO_pi2: %.10f, ", 
           bo->BO, bo->BO_s, bo->BO_pi, bo->BO_pi2);
    printf("dBOp: %.10f, dln_BOp_s: %.10f, dln_BOp_pi: %.10f, dln_BOp_pi2: %.10f",
           bo->dBOp, bo->dln_BOp_s, bo->dln_BOp_pi, bo->dln_BOp_pi2);
    int j = data->nbr;
    double *xi = system->my_atoms[i].x;
    double *xj = system->my_atoms[j].x;
    printf("i: [%.10f %.10f %.10f], ", xi[0], xi[1], xi[2]);
    printf("j: [%.10f %.10f %.10f]\n", xj[0], xj[1], xj[2]);
  }

  int cmp_bonds(const void *a, const void *b){
    bond_data *abond = (bond_data*)a;
    bond_data *bbond = (bond_data*)b;
    if (abond->nbr < bbond->nbr)
      return -1;
    if (abond->nbr > bbond->nbr)
      return 1;
    return 0;
  }

  static void sort_bonds(reax_list *bonds, reax_system *system){
    int i, j;
    for (i = 0; i < system->N; i ++){
      int i_start = Start_Index(i, bonds);
      int i_stop = End_Index(i, bonds);
      qsort(bonds->select.bond_list + i_start, i_stop - i_start, sizeof(bond_data), cmp_bonds);
    }
    for (i = 0; i < system->N; i ++){
      int start_i = Start_Index( i, bonds );
      int end_i = End_Index( i, bonds);
      int pj;
      for( pj = start_i; pj < end_i; ++pj ) {
        bond_data *jbond = bonds->select.bond_list + pj;
        j = jbond->nbr;
        int start_j = Start_Index( j, bonds);
        int end_j = End_Index( j, bonds );
        int pk;
        int flag = 0;
        for ( pk = start_j; pk < end_j; ++pk ) {
          bond_data *kbond = bonds->select.bond_list + pk;
          int k = kbond->nbr;
          if (k == i){
            jbond->sym_index = pk;
            flag = 1;
          }
        }
        if (!flag)
          printf("sym not found %d %d\n", i, j);
      }
    }
  }
  static void print_bonds(reax_list *bonds, reax_system *system){
    int i, j;
    for (i = 0; i < system->N; i ++){
      int i_start = Start_Index(i, bonds);
      int i_stop = End_Index(i, bonds);
      int pj;
      for (pj = i_start; pj < i_stop; ++pj){
        j = bonds->select.bond_list[pj].nbr;
        print_bond(i, bonds->select.bond_list + pj, system);
      }
    }
  }
  // void Init_Forces_noQEq( reax_system *system, control_params *control,
  //                         simulation_data *data, storage *workspace,
  //                         reax_list **lists, output_controls *out_control,
  //                         MPI_Comm comm ) {
  //   int i, j, pj;
  //   int start_i, end_i;
  //   int type_i, type_j;
  //   int btop_i, btop_j, num_bonds, num_hbonds;
  //   int ihb, jhb, ihb_top, jhb_top;
  //   int local, flag, renbr;
  //   double cutoff;
  //   reax_list *far_nbrs, *bonds, *hbonds;
  //   single_body_parameters *sbp_i, *sbp_j;
  //   two_body_parameters *twbp;
  //   far_neighbor_data *nbr_pj;
  //   reax_atom *atom_i, *atom_j;

  //   far_nbrs = *lists + FAR_NBRS;
  //   bonds = *lists + BONDS;
  //   hbonds = *lists + HBONDS;

  //   for( i = 0; i < system->n; ++i )
  //     workspace->bond_mark[i] = 0;
  //   for( i = system->n; i < system->N; ++i ) {
  //     workspace->bond_mark[i] = 1000; // put ghost atoms to an infinite distance
  //   }

  //   num_bonds = 0;
  //   num_hbonds = 0;
  //   btop_i = btop_j = 0;
  //   renbr = (data->step-data->prev_steps) % control->reneighbor == 0;
  //   //printf("%d\n", renbr);
  //   for( i = 0; i < system->N; ++i ) {
  //     atom_i = &(system->my_atoms[i]);
  //     type_i  = atom_i->type;
  //     if (type_i < 0) continue;
  //     start_i = Start_Index(i, far_nbrs);
  //     end_i   = End_Index(i, far_nbrs);
  //     btop_i = End_Index( i, bonds );
  //     sbp_i = &(system->reax_param.sbp[type_i]);

  //     if( i < system->n ) {
  //       local = 1;
  //       cutoff = MAX( control->hbond_cut, control->bond_cut );
  //     }
  //     else {
  //       local = 0;
  //       cutoff = control->bond_cut;
  //     }

  //     ihb = -1;
  //     ihb_top = -1;
  //     if( local && control->hbond_cut > 0 ) {
  //       ihb = sbp_i->p_hbond;
  //       if( ihb == 1 )
  //         ihb_top = End_Index( atom_i->Hindex, hbonds );
  //       else ihb_top = -1;
  //     }

  //     /* update i-j distance - check if j is within cutoff */
  //     for( pj = start_i; pj < end_i; ++pj ) {
  //       nbr_pj = &( far_nbrs->select.far_nbr_list[pj] );
  //       j = nbr_pj->nbr;
  //       atom_j = &(system->my_atoms[j]);

  //       if( renbr ) {
  //         if( nbr_pj->d <= cutoff )
  //           flag = 1;
  //         else flag = 0;
  //       }
  //       else{
  //         nbr_pj->dvec[0] = atom_j->x[0] - atom_i->x[0];
  //         nbr_pj->dvec[1] = atom_j->x[1] - atom_i->x[1];
  //         nbr_pj->dvec[2] = atom_j->x[2] - atom_i->x[2];
  //         nbr_pj->d = rvec_Norm_Sqr( nbr_pj->dvec );
  //         if( nbr_pj->d <= SQR(cutoff) ) {
  //           nbr_pj->d = sqrt(nbr_pj->d);
  //           flag = 1;
  //         }
  //         else {
  //           flag = 0;
  //         }
  //       }

  //       if( flag ) {
  //         type_j = atom_j->type;
  //         if (type_j < 0) continue;
  //         sbp_j = &(system->reax_param.sbp[type_j]);
  //         twbp = &(system->reax_param.tbp[type_i][type_j]);

  //         if( local ) {
  //           /* hydrogen bond lists */
  //           if( control->hbond_cut > 0 && (ihb==1 || ihb==2) &&
  //               nbr_pj->d <= control->hbond_cut ) {
  //             // fprintf( stderr, "%d %d\n", atom1, atom2 );
  //             jhb = sbp_j->p_hbond;
  //             if( ihb == 1 && jhb == 2 ) {
  //               hbonds->select.hbond_list[ihb_top].nbr = j;
  //               hbonds->select.hbond_list[ihb_top].scl = 1;
  //               hbonds->select.hbond_list[ihb_top].ptr = nbr_pj;
  //               ++ihb_top;
  //               ++num_hbonds;
  //             }
  //             else if( j < system->n && ihb == 2 && jhb == 1 ) {
  //               jhb_top = End_Index( atom_j->Hindex, hbonds );
  //               hbonds->select.hbond_list[jhb_top].nbr = i;
  //               hbonds->select.hbond_list[jhb_top].scl = -1;
  //               hbonds->select.hbond_list[jhb_top].ptr = nbr_pj;
  //               Set_End_Index( atom_j->Hindex, jhb_top+1, hbonds );
  //               ++num_hbonds;
  //             }
  //           }
  //         }

  //         if( //(workspace->bond_mark[i] < 3 || workspace->bond_mark[j] < 3) &&
  //            nbr_pj->d <= control->bond_cut &&
  //            BOp( workspace, bonds, control->bo_cut,
  //                 i , btop_i, nbr_pj, sbp_i, sbp_j, twbp ) ) {
  //           num_bonds += 2;
  //           ++btop_i;

  //           if( workspace->bond_mark[j] > workspace->bond_mark[i] + 1 )
  //             workspace->bond_mark[j] = workspace->bond_mark[i] + 1;
  //           else if( workspace->bond_mark[i] > workspace->bond_mark[j] + 1 ) {
  //             workspace->bond_mark[i] = workspace->bond_mark[j] + 1;
  //           }
  //         }
  //       }
  //     }

  //     Set_End_Index( i, btop_i, bonds );
  //     if( local && ihb == 1 )
  //       Set_End_Index( atom_i->Hindex, ihb_top, hbonds );
  //   }


  //   workspace->realloc.num_bonds = num_bonds;
  //   workspace->realloc.num_hbonds = num_hbonds;

  //   Validate_Lists( system, workspace, lists, data->step,
  //                   system->n, system->N, system->numH, comm );

  //   // for (i = 0; i < system->N; i ++){
  //   //   int i_start = Start_Index(i, bonds);
  //   //   int i_stop = End_Index(i, bonds);
  //   //   int pj;
  //   //   for (pj = i_start; pj < i_stop; ++pj){
  //   //     j = bonds->select.bond_list[pj].nbr;
  //   //     //if (i == 190 && j == 203){
  //   //     print_bond(i, bonds->select.bond_list + pj, system);
  //   //     // }
  //   //   }
  //   // }

  // }


  int BOp_single( storage *workspace, reax_list *bonds, double bo_cut,
                  int i, int btop_i, far_neighbor_data_full *nbr_pj,
                  single_body_parameters *sbp_i, single_body_parameters *sbp_j,
                  two_body_parameters *twbp ) {
    double r2, C12, C34, C56;
    double Cln_BOp_s, Cln_BOp_pi, Cln_BOp_pi2;
    double BO, BO_s, BO_pi, BO_pi2;
    bond_data *ibond;
    bond_order_data *bo_ij;

    r2 = SQR(nbr_pj->d);

    if( sbp_i->r_s > 0.0 && sbp_j->r_s > 0.0 ) {
      C12 = twbp->p_bo1 * pow( nbr_pj->d / twbp->r_s, twbp->p_bo2 );
      BO_s = (1.0 + bo_cut) * exp( C12 );
    }
    else BO_s = C12 = 0.0;

    if( sbp_i->r_pi > 0.0 && sbp_j->r_pi > 0.0 ) {
      C34 = twbp->p_bo3 * pow( nbr_pj->d / twbp->r_p, twbp->p_bo4 );
      BO_pi = exp( C34 );
    }
    else BO_pi = C34 = 0.0;

    if( sbp_i->r_pi_pi > 0.0 && sbp_j->r_pi_pi > 0.0 ) {
      C56 = twbp->p_bo5 * pow( nbr_pj->d / twbp->r_pp, twbp->p_bo6 );
      BO_pi2= exp( C56 );
    }
    else BO_pi2 = C56 = 0.0;

    /* Initially BO values are the uncorrected ones, page 1 */
    BO = BO_s + BO_pi + BO_pi2;

    if( BO >= bo_cut ) {
      /****** bonds i-j and j-i ******/
      ibond = &( bonds->select.bond_list[btop_i] );

      ibond->nbr = nbr_pj->nbr;
      ibond->d = nbr_pj->d;
      //jbond->d = nbr_pj->d;
      rvec_Copy( ibond->dvec, nbr_pj->dvec );
      //rvec_Scale( jbond->dvec, -1, nbr_pj->dvec );
      ivec_Copy( ibond->rel_box, nbr_pj->rel_box );
      //ivec_Scale( jbond->rel_box, -1, nbr_pj->rel_box );
      ibond->dbond_index = btop_i;
      //jbond->dbond_index = btop_i;
      //ibond->sym_index = btop_j;
      //jbond->sym_index = btop_i;
      //Set_End_Index( j, btop_j+1, bonds );

      bo_ij = &( ibond->bo_data );
      // bo_ji = &( jbond->bo_data );
      // bo_ji->BO = 
      bo_ij->BO = BO;
      // bo_ji->BO_s = 
      bo_ij->BO_s = BO_s;
      // bo_ji->BO_pi = 
      bo_ij->BO_pi = BO_pi;
      // bo_ji->BO_pi2 = 
      bo_ij->BO_pi2 = BO_pi2;

      /* Bond Order page2-3, derivative of total bond order prime */
      Cln_BOp_s = twbp->p_bo2 * C12 / r2;
      Cln_BOp_pi = twbp->p_bo4 * C34 / r2;
      Cln_BOp_pi2 = twbp->p_bo6 * C56 / r2;

      /* Only dln_BOp_xx wrt. dr_i is stored here, note that
         dln_BOp_xx/dr_i = -dln_BOp_xx/dr_j and all others are 0 */
      //printf("%f %f %f\n", ibond->dvec[0], ibond->dvec[1], ibond->dvec[2]);
      //, bo_ij->BO_s);
      rvec_Scale(bo_ij->dln_BOp_s,-bo_ij->BO_s*Cln_BOp_s,ibond->dvec);
      rvec_Scale(bo_ij->dln_BOp_pi,-bo_ij->BO_pi*Cln_BOp_pi,ibond->dvec);
      rvec_Scale(bo_ij->dln_BOp_pi2,
                 -bo_ij->BO_pi2*Cln_BOp_pi2,ibond->dvec);
      // rvec_Scale(bo_ji->dln_BOp_s, -1., bo_ij->dln_BOp_s);
      // rvec_Scale(bo_ji->dln_BOp_pi, -1., bo_ij->dln_BOp_pi );
      // rvec_Scale(bo_ji->dln_BOp_pi2, -1., bo_ij->dln_BOp_pi2 );

      rvec_Scale( bo_ij->dBOp,
                  -(bo_ij->BO_s * Cln_BOp_s +
                    bo_ij->BO_pi * Cln_BOp_pi +
                    bo_ij->BO_pi2 * Cln_BOp_pi2), ibond->dvec );
      //rvec_Scale( bo_ji->dBOp, -1., bo_ij->dBOp );

      rvec_Add( workspace->dDeltap_self[i], bo_ij->dBOp );
      //rvec_Add( workspace->dDeltap_self[j], bo_ji->dBOp );

      bo_ij->BO_s -= bo_cut;
      bo_ij->BO -= bo_cut;
      // bo_ji->BO_s -= bo_cut;
      // bo_ji->BO -= bo_cut;
      workspace->total_bond_order[i] += bo_ij->BO; //currently total_BOp
      //workspace->total_bond_order[j] += bo_ji->BO; //currently total_BOp
      bo_ij->Cdbo = bo_ij->Cdbopi = bo_ij->Cdbopi2 = 0.0;
      //bo_ji->Cdbo = bo_ji->Cdbopi = bo_ji->Cdbopi2 = 0.0;
      // if (i == 150)
      //print_bond(i, ibond);
      return 1;
    }

    return 0;
  }


  void Init_Forces_noQEq_Full( reax_system *system, control_params *control,
                               simulation_data *data, storage *workspace,
                               reax_list **lists, output_controls *out_control,
                               MPI_Comm comm ) {
    int i, j, pj;
    int start_i, end_i;
    int type_i, type_j;
    int btop_i, btop_j, num_bonds, num_hbonds;
    int ihb, jhb, ihb_top, jhb_top;
    int local, flag, renbr;
    double cutoff;
    reax_list *far_nbrs, *bonds, *hbonds;
    single_body_parameters *sbp_i, *sbp_j;
    two_body_parameters *twbp;
    far_neighbor_data_full *nbr_pj;
    reax_atom *atom_i, *atom_j;

    far_nbrs = *lists + FAR_NBRS_FULL;
    bonds = *lists + BONDS;
    hbonds = *lists + HBONDS;

    for( i = 0; i < system->n; ++i )
      workspace->bond_mark[i] = 0;
    for( i = system->n; i < system->N; ++i ) {
      workspace->bond_mark[i] = 1000; // put ghost atoms to an infinite distance
    }

    num_bonds = 0;
    num_hbonds = 0;
    btop_i = btop_j = 0;
    renbr = (data->step-data->prev_steps) % control->reneighbor == 0;

    for( i = 0; i < system->N; ++i ) {
      atom_i = &(system->my_atoms[i]);
      type_i  = atom_i->type;
      if (type_i < 0) continue;
      start_i = Start_Index(i, far_nbrs);
      end_i   = End_Index(i, far_nbrs);
      btop_i = End_Index( i, bonds );
      sbp_i = &(system->reax_param.sbp[type_i]);

      cutoff = control->bond_cut;

      /* update i-j distance - check if j is within cutoff */
      for( pj = start_i; pj < end_i; ++pj ) {
        nbr_pj = &( far_nbrs->select.far_nbr_list_full[pj] );
        j = nbr_pj->nbr;
        atom_j = &(system->my_atoms[j]);

        if( nbr_pj->d <= cutoff )
          flag = 1;
        else flag = 0;


        if( flag ) {
          type_j = atom_j->type;
          if (type_j < 0) continue;
          sbp_j = &(system->reax_param.sbp[type_j]);
          twbp = &(system->reax_param.tbp[type_i][type_j]);

          if(nbr_pj->d <= control->bond_cut &&
             BOp_single( workspace, bonds, control->bo_cut,
                         i , btop_i, nbr_pj, sbp_i, sbp_j, twbp ) ) {
            num_bonds += 1;
            ++btop_i;
            if (j > i){
              if( workspace->bond_mark[j] > workspace->bond_mark[i] + 1 )
                workspace->bond_mark[j] = workspace->bond_mark[i] + 1;
              else 
                if( workspace->bond_mark[i] > workspace->bond_mark[j] + 1 ) {
                  workspace->bond_mark[i] = workspace->bond_mark[j] + 1;
                }
            }
          }
        }
      }

      Set_End_Index( i, btop_i, bonds );
    }

    workspace->realloc.num_bonds = num_bonds;
  }
  // void Init_Forces_noQEq_HB( reax_system *system, control_params *control,
  //                            simulation_data *data, storage *workspace,
  //                            reax_list **lists, output_controls *out_control,
  //                            MPI_Comm comm ) {
  //   int i, j, pj;
  //   int start_i, end_i;
  //   int type_i, type_j;
  //   //int btop_i, btop_j, num_bonds
  //   int num_hbonds;
  //   int ihb, jhb, ihb_top, jhb_top;
  //   int local, flag, renbr;
  //   double cutoff;
  //   reax_list *far_nbrs, *hbonds;
  //   single_body_parameters *sbp_i, *sbp_j;
  //   two_body_parameters *twbp;
  //   far_neighbor_data *nbr_pj;
  //   reax_atom *atom_i, *atom_j;

  //   far_nbrs = *lists + FAR_NBRS;
  //   hbonds = *lists + HBONDS;

  //   num_hbonds = 0;
  //   renbr = (data->step-data->prev_steps) % control->reneighbor == 0;

  //   for( i = 0; i < system->N; ++i ) {
  //     atom_i = &(system->my_atoms[i]);
  //     type_i  = atom_i->type;
  //     if (type_i < 0) continue;
  //     start_i = Start_Index(i, far_nbrs);
  //     end_i   = End_Index(i, far_nbrs);
  //     //btop_i = End_Index( i, bonds );
  //     sbp_i = &(system->reax_param.sbp[type_i]);

  //     if( i < system->n ) {
  //       local = 1;
  //       cutoff = MAX( control->hbond_cut, control->bond_cut );
  //     }
  //     else {
  //       local = 0;
  //       cutoff = control->bond_cut;
  //     }

  //     ihb = -1;
  //     ihb_top = -1;
  //     if( local && control->hbond_cut > 0 ) {
  //       ihb = sbp_i->p_hbond;
  //       if( ihb == 1 )
  //         ihb_top = End_Index( atom_i->Hindex, hbonds );
  //       else ihb_top = -1;
  //     }

  //     /* update i-j distance - check if j is within cutoff */
  //     for( pj = start_i; pj < end_i; ++pj ) {
  //       nbr_pj = &( far_nbrs->select.far_nbr_list[pj] );
  //       j = nbr_pj->nbr;
  //       atom_j = &(system->my_atoms[j]);

  //       if( nbr_pj->d <= cutoff )
  //         flag = 1;
  //       else flag = 0;

  //       if( flag ) {
  //         type_j = atom_j->type;
  //         if (type_j < 0) continue;
  //         sbp_j = &(system->reax_param.sbp[type_j]);
  //         twbp = &(system->reax_param.tbp[type_i][type_j]);

  //         if( local ) {
  //           /* hydrogen bond lists */
  //           if( control->hbond_cut > 0 && (ihb==1 || ihb==2) &&
  //               nbr_pj->d <= control->hbond_cut ) {
  //             // fprintf( stderr, "%d %d\n", atom1, atom2 );
  //             jhb = sbp_j->p_hbond;
  //             if( ihb == 1 && jhb == 2 ) {
  //               hbonds->select.hbond_list[ihb_top].nbr = j;
  //               hbonds->select.hbond_list[ihb_top].scl = 1;
  //               hbonds->select.hbond_list[ihb_top].ptr = nbr_pj;
  //               ++ihb_top;
  //               ++num_hbonds;
  //             }
  //             else if( j < system->n && ihb == 2 && jhb == 1 ) {
  //               jhb_top = End_Index( atom_j->Hindex, hbonds );
  //               hbonds->select.hbond_list[jhb_top].nbr = i;
  //               hbonds->select.hbond_list[jhb_top].scl = -1;
  //               hbonds->select.hbond_list[jhb_top].ptr = nbr_pj;
  //               Set_End_Index( atom_j->Hindex, jhb_top+1, hbonds );
  //               ++num_hbonds;
  //             }
  //           }
  //         }
  //       }
  //     }

  //     //Set_End_Index( i, btop_i, bonds );
  //     if( local && ihb == 1 )
  //       Set_End_Index( atom_i->Hindex, ihb_top, hbonds );
  //   }


  //   workspace->realloc.num_hbonds = num_hbonds;
  //   printf("%d\n", num_hbonds);
  //   Validate_Lists( system, workspace, lists, data->step,
  //                   system->n, system->N, system->numH, comm );
  // }

  void Init_Forces_noQEq_HB_Full( reax_system *system, control_params *control,
                                  simulation_data *data, storage *workspace,
                                  reax_list **lists, output_controls *out_control,
                                  MPI_Comm comm ) {
    int i, j, pj;
    int start_i, end_i;
    int type_i, type_j;
    //int btop_i, btop_j, num_bonds
    int num_hbonds;
    int ihb, jhb, ihb_top;
    int local, flag, renbr;
    double cutoff;
    reax_list *far_nbrs, *hbonds;
    single_body_parameters *sbp_i, *sbp_j;
    two_body_parameters *twbp;
    far_neighbor_data_full *nbr_pj;
    reax_atom *atom_i, *atom_j;

    far_nbrs = *lists + FAR_NBRS_FULL;
    hbonds = *lists + HBONDS;

    num_hbonds = 0;
    renbr = (data->step-data->prev_steps) % control->reneighbor == 0;
    if (control->hbond_cut <= 0)
      return;
    for( i = 0; i < system->n; ++i ) {
      atom_i = &(system->my_atoms[i]);
      type_i  = atom_i->type;
      if (type_i < 0) continue;
      start_i = Start_Index(i, far_nbrs);
      end_i   = End_Index(i, far_nbrs);
      sbp_i = &(system->reax_param.sbp[type_i]);

      local = 1;
      cutoff = control->hbond_cut;

      ihb = -1;
      ihb_top = -1;
      ihb = sbp_i->p_hbond;
      if( ihb == 1 )
        ihb_top = End_Index( atom_i->Hindex, hbonds );
      else continue;

      for( pj = start_i; pj < end_i; ++pj ) {
        nbr_pj = &( far_nbrs->select.far_nbr_list_full[pj] );
        j = nbr_pj->nbr;
        atom_j = &(system->my_atoms[j]);

        if( nbr_pj->d <= cutoff )
          flag = 1;
        else flag = 0;

        if( flag ) {
          type_j = atom_j->type;
          if (type_j < 0) continue;
          sbp_j = &(system->reax_param.sbp[type_j]);
          twbp = &(system->reax_param.tbp[type_i][type_j]);

          if( local ) {
            if( (ihb==1 || ihb==2)) {
              jhb = sbp_j->p_hbond;
              if( ihb == 1 && jhb == 2 ) {
                hbonds->select.hbond_list[ihb_top].nbr = j;
                hbonds->select.hbond_list[ihb_top].scl = 1;
                hbonds->select.hbond_list[ihb_top].ptr = nbr_pj;
                ++ihb_top;
                ++num_hbonds;
              }
            }
          }
        }
      }

      if(ihb == 1 )
        Set_End_Index( atom_i->Hindex, ihb_top, hbonds );
    }

    workspace->realloc.num_hbonds = num_hbonds;
    Validate_Lists( system, workspace, lists, data->step,
                    system->n, system->N, system->numH, comm );
  }


  void Estimate_Storages( reax_system *system, control_params *control,
                          reax_list **lists, int *Htop, int *hb_top,
                          int *bond_top, int *num_3body, MPI_Comm comm )
  {
    int i, j, pj;
    int start_i, end_i;
    int type_i, type_j;
    int ihb, jhb;
    int local;
    double cutoff;
    double r_ij;
    double C12, C34, C56;
    double BO, BO_s, BO_pi, BO_pi2;
    reax_list *far_nbrs;
    single_body_parameters *sbp_i, *sbp_j;
    two_body_parameters *twbp;
    far_neighbor_data *nbr_pj;
    reax_atom *atom_i, *atom_j;

    int mincap = system->mincap;
    double safezone = system->safezone;
    double saferzone = system->saferzone;

    far_nbrs = *lists + FAR_NBRS;
    *Htop = 0;
    memset( hb_top, 0, sizeof(int) * system->local_cap );
    memset( bond_top, 0, sizeof(int) * system->total_cap );
    *num_3body = 0;

    for( i = 0; i < system->N; ++i ) {
      atom_i = &(system->my_atoms[i]);
      type_i  = atom_i->type;
      if (type_i < 0) continue;
      start_i = Start_Index(i, far_nbrs);
      end_i   = End_Index(i, far_nbrs);
      sbp_i = &(system->reax_param.sbp[type_i]);

      if( i < system->n ) {
        local = 1;
        cutoff = control->nonb_cut;
        ++(*Htop);
        ihb = sbp_i->p_hbond;
      }
      else {
        local = 0;
        cutoff = control->bond_cut;
        ihb = -1;
      }

      for( pj = start_i; pj < end_i; ++pj ) {
        nbr_pj = &( far_nbrs->select.far_nbr_list[pj] );
        j = nbr_pj->nbr;
        atom_j = &(system->my_atoms[j]);

        if(nbr_pj->d <= cutoff) {
          type_j = system->my_atoms[j].type;
          if (type_j < 0) continue;
          r_ij = nbr_pj->d;
          sbp_j = &(system->reax_param.sbp[type_j]);
          twbp = &(system->reax_param.tbp[type_i][type_j]);

          if( local ) {
            if( j < system->n || atom_i->orig_id < atom_j->orig_id ) //tryQEq ||1
              ++(*Htop);

            /* hydrogen bond lists */
            if( control->hbond_cut > 0.1 && (ihb==1 || ihb==2) &&
                nbr_pj->d <= control->hbond_cut ) {
              jhb = sbp_j->p_hbond;
              if( ihb == 1 && jhb == 2 )
                ++hb_top[i];
              else if( j < system->n && ihb == 2 && jhb == 1 )
                ++hb_top[j];
            }
          }

          /* uncorrected bond orders */
          if( nbr_pj->d <= control->bond_cut ) {
            if( sbp_i->r_s > 0.0 && sbp_j->r_s > 0.0) {
              C12 = twbp->p_bo1 * pow( r_ij / twbp->r_s, twbp->p_bo2 );
              BO_s = (1.0 + control->bo_cut) * exp( C12 );
            }
            else BO_s = C12 = 0.0;

            if( sbp_i->r_pi > 0.0 && sbp_j->r_pi > 0.0) {
              C34 = twbp->p_bo3 * pow( r_ij / twbp->r_p, twbp->p_bo4 );
              BO_pi = exp( C34 );
            }
            else BO_pi = C34 = 0.0;

            if( sbp_i->r_pi_pi > 0.0 && sbp_j->r_pi_pi > 0.0) {
              C56 = twbp->p_bo5 * pow( r_ij / twbp->r_pp, twbp->p_bo6 );
              BO_pi2= exp( C56 );
            }
            else BO_pi2 = C56 = 0.0;

            /* Initially BO values are the uncorrected ones, page 1 */
            BO = BO_s + BO_pi + BO_pi2;

            if( BO >= control->bo_cut ) {
              ++bond_top[i];
              ++bond_top[j];
            }
          }
        }
      }
    }

    *Htop = (int)(MAX( *Htop * safezone, mincap * MIN_HENTRIES ));
    for( i = 0; i < system->n; ++i )
      hb_top[i] = (int)(MAX( hb_top[i] * saferzone, MIN_HBONDS ));

    for( i = 0; i < system->N; ++i ) {
      *num_3body += SQR(bond_top[i]);
      bond_top[i] = MAX( bond_top[i] * 2, MIN_BONDS );
    }

  }


  void Compute_Forces( reax_system *system, control_params *control,
                       simulation_data *data, storage *workspace,
                       reax_list **lists, output_controls *out_control,
                       mpi_datatypes *mpi_data )
  {
    MPI_Comm comm = mpi_data->world;
    GPTLstart("reaxc init forces noqeq");
    // if (getenv("NOMOD"))
    //   Init_Forces_noQEq( system, control, data, workspace,
    //                      lists, out_control, comm );
    // else{
    Init_Forces_noQEq_Full( system, control, data, workspace,
                            lists, out_control, comm );
    Init_Forces_noQEq_HB_Full( system, control, data, workspace,
                               lists, out_control, comm );
    //}
    //int i;
    sort_bonds(*lists + BONDS, system);
    // print_bonds(*lists + BONDS, system);
    // for (i = 0; i < system->N; i ++)
    //   printf("bmk: %d %d\n", i, workspace->bond_mark[i]);

    GPTLstop("reaxc init forces noqeq");
    /********* bonded interactions ************/
    GPTLstart("reaxc bonded forces");
    Compute_Bonded_Forces( system, control, data, workspace,
                           lists, out_control, mpi_data->world );
    GPTLstop("reaxc bonded forces");
    /********* nonbonded interactions ************/
    GPTLstart("reaxc nonbonded forces");
    Compute_NonBonded_Forces( system, control, data, workspace,
                              lists, out_control, mpi_data->world );
    GPTLstop("reaxc nonbonded forces");
    /*********** total force ***************/
    GPTLstart("reaxc total force");
    Compute_Total_Force( system, control, data, workspace, lists, mpi_data );
    GPTLstop("reaxc total force");

  }

}
