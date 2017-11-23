#include "sunway.h"
#include <stdlib.h>
#include <simd.h>
#include <stdio.h>

#ifdef MPE
extern SLAVE_FUN(fix_nve_initial_integrate_sunway_compute_para)(fix_nve_param_t *pm);
extern SLAVE_FUN(fix_nve_final_integrate_sunway_compute_para)(fix_nve_param_t *pm);


void fix_nve_initial_integrate_sunway_compute(fix_nve_param_t *pm){
  
  
  if (athread_idle() == 0)
    athread_init();


  athread_spawn(fix_nve_initial_integrate_sunway_compute_para, pm);
  //puts("spawned");
  athread_join();
  //puts("joined");


}


void fix_nve_final_integrate_sunway_compute(fix_nve_param_t *pm){
  
  
  if (athread_idle() == 0)
    athread_init();


  athread_spawn(fix_nve_final_integrate_sunway_compute_para, pm);
  //puts("spawned");
  athread_join();
  //puts("joined");


}
#endif



#ifdef CPE

#define C_SIZE_ 512
#define Cach_Size 16
#define Cach_Num 16

void fix_nve_final_integrate_sunway_compute_para(fix_nve_param_t *pm){
  
  pe_init();
  int i,j,k,g,h,ii;
  fix_nve_param_t local_pm;
  pe_get(pm, &local_pm, sizeof(fix_nve_param_t));
  pe_syn();

//printf("fix_nve\n");
  double x[C_SIZE_ * 3], v[C_SIZE_ * 3], f[C_SIZE_ * 3];
  double rmass[C_SIZE_], mass[Cach_Num][Cach_Size];
  int mask[C_SIZE_], type[C_SIZE_], mass_id[Cach_Num];

  int mass_num = -1;
  for(i = 0;i < Cach_Num;i++)
    mass_id[i] = -1;
  int my_id = athread_get_id(-1);

  int nlocal = local_pm.nlocal;
  int groupbit = local_pm.groupbit;
  double dtv = local_pm.dtv; 
  double dtf = local_pm.dtf;


  int C_SIZE = nlocal / 64;
  if(C_SIZE > C_SIZE_)
    C_SIZE = C_SIZE_;



  if(local_pm.rmass)
  {
    for(ii = my_id * C_SIZE; ii < nlocal; ii += 64 * C_SIZE)
    {
      if(C_SIZE + ii > nlocal)
        C_SIZE = nlocal - ii;

      pe_get( local_pm.v[ii], v, sizeof(double) * C_SIZE * 3);
      pe_syn();
      pe_get( local_pm.f[ii], f, sizeof(double) * C_SIZE * 3);
      pe_syn();
      pe_get( &local_pm.rmass[ii], rmass, sizeof(double) * C_SIZE);
      pe_syn();
      pe_get( &local_pm.mask[ii], mask, sizeof(int) * C_SIZE);
      pe_syn();


      for(i = 0;i < C_SIZE; i++)
      {
        if(mask[i] & groupbit)
        {
          int pi_3 = 3 * i;
          double dtfm = dtf / rmass[i];
          v[pi_3 + 0] += dtfm * f[pi_3 + 0];
          v[pi_3 + 1] += dtfm * f[pi_3 + 1];
          v[pi_3 + 2] += dtfm * f[pi_3 + 2];
        }
      }
    }
    
    pe_put( local_pm.v[ii], v, sizeof(double) * C_SIZE * 3);
    pe_syn();
  }
  else
  {
    for(ii = C_SIZE * my_id; ii < nlocal; ii += C_SIZE * 64)
    {
      if(C_SIZE + ii > nlocal)
        C_SIZE = nlocal - ii;


      pe_get( local_pm.v[ii], v, sizeof(double) * C_SIZE * 3);
      pe_syn();
      pe_get( local_pm.f[ii], f, sizeof(double) * C_SIZE * 3);
      pe_syn();
      pe_get( &local_pm.type[ii], type, sizeof(int) * C_SIZE);
      pe_syn();
      pe_get( &local_pm.mask[ii], mask, sizeof(int) * C_SIZE);
      pe_syn();

//
//
      for(i = 0;i < C_SIZE; i++)
      {
        if(mask[i] & groupbit)
        {
          
          int pi_3 = 3 * i;
          int TP = type[i];

          int mass_num_id = (15) & (TP >> 4);

          if(mass_id[mass_num_id] != TP >> 4)
          {

            mass_id[mass_num_id] = TP >> 4;
            pe_get( &local_pm.mass[mass_id[mass_num_id] << 4], mass[mass_num_id], sizeof(double) * Cach_Size);
            pe_syn();

          }

          

          double dtfm = dtf / mass[mass_num_id][type[i]&(15)];
          v[pi_3 + 0] += dtfm * f[pi_3 + 0];
          v[pi_3 + 1] += dtfm * f[pi_3 + 1];
          v[pi_3 + 2] += dtfm * f[pi_3 + 2];
        }
      }
      pe_put( local_pm.v[ii], v, sizeof(double) * C_SIZE * 3);
      pe_syn();

    }

  }

}


void fix_nve_initial_integrate_sunway_compute_para(fix_nve_param_t *pm){
  
  pe_init();
  int i,j,k,g,h,ii;
  fix_nve_param_t local_pm;
  pe_get(pm, &local_pm, sizeof(fix_nve_param_t));
  pe_syn();


  double x[C_SIZE_ * 3], v[C_SIZE_ * 3], f[C_SIZE_ * 3];
  double rmass[C_SIZE_], mass[Cach_Num][Cach_Size];
  int mask[C_SIZE_], type[C_SIZE_], mass_id[Cach_Num];

  int mass_num = -1;
  for(i = 0;i < Cach_Num;i++)
    mass_id[i] = -1;
  int my_id = athread_get_id(-1);

  int nlocal = local_pm.nlocal;
  int groupbit = local_pm.groupbit;
  double dtv = local_pm.dtv; 
  double dtf = local_pm.dtf;


  int C_SIZE = nlocal / 64;
  if(C_SIZE > C_SIZE_)
    C_SIZE = C_SIZE_;



  if(local_pm.rmass)
  {
    for(ii = my_id * C_SIZE; ii < nlocal; ii += 64 * C_SIZE)
    {
      if(C_SIZE + ii > nlocal)
        C_SIZE = nlocal - ii;


      pe_get( local_pm.v[ii], v, sizeof(double) * C_SIZE * 3);
      pe_syn();
      pe_get( local_pm.f[ii], f, sizeof(double) * C_SIZE * 3);
      pe_syn();
      pe_get( local_pm.x[ii], x, sizeof(double) * C_SIZE * 3);
      pe_syn();
      pe_get( &local_pm.rmass[ii], rmass, sizeof(double) * C_SIZE);
      pe_syn();
      pe_get( &local_pm.mask[ii], mask, sizeof(int) * C_SIZE);
      pe_syn();


      for(i = 0;i < C_SIZE; i++)
      {
        if(mask[i] & groupbit)
        {
          int pi_3 = 3 * i;
          double dtfm = dtf / rmass[i];
          v[pi_3 + 0] += dtfm * f[pi_3 + 0];
          v[pi_3 + 1] += dtfm * f[pi_3 + 1];
          v[pi_3 + 2] += dtfm * f[pi_3 + 2];

          x[pi_3 + 0] += dtv * v[pi_3 + 0]; 
          x[pi_3 + 1] += dtv * v[pi_3 + 1]; 
          x[pi_3 + 2] += dtv * v[pi_3 + 2]; 

        }
      }
    }
    
    pe_put( local_pm.v[ii], v, sizeof(double) * C_SIZE * 3);
    pe_syn();
    pe_put( local_pm.x[ii], x, sizeof(double) * C_SIZE * 3);
    pe_syn();


  }
  else
  {
    for(ii = C_SIZE * my_id; ii < nlocal; ii += C_SIZE * 64)
    {
      if(C_SIZE + ii > nlocal)
        C_SIZE = nlocal - ii;


      pe_get( local_pm.v[ii], v, sizeof(double) * C_SIZE * 3);
      pe_syn();
      pe_get( local_pm.f[ii], f, sizeof(double) * C_SIZE * 3);
      pe_syn();
      pe_get( local_pm.x[ii], x, sizeof(double) * C_SIZE * 3);
      pe_syn();
      pe_get( &local_pm.type[ii], type, sizeof(int) * C_SIZE);
      pe_syn();
      pe_get( &local_pm.mask[ii], mask, sizeof(int) * C_SIZE);
      pe_syn();

//
//
      for(i = 0;i < C_SIZE; i++)
      {
        if(mask[i] & groupbit)
        {
          
          int pi_3 = 3 * i;
          int TP = type[i];

          int mass_num_id = (15) & (TP >> 4);

          if(mass_id[mass_num_id] != TP >> 4)
          {

            mass_id[mass_num_id] = TP >> 4;
            pe_get( &local_pm.mass[mass_id[mass_num_id] << 4], mass[mass_num_id], sizeof(double) * Cach_Size);
            pe_syn();

          }

          

          double dtfm = dtf / mass[mass_num_id][type[i]&(15)];
          v[pi_3 + 0] += dtfm * f[pi_3 + 0];
          v[pi_3 + 1] += dtfm * f[pi_3 + 1];
          v[pi_3 + 2] += dtfm * f[pi_3 + 2];

          x[pi_3 + 0] += dtv  * v[pi_3 + 0]; 
          x[pi_3 + 1] += dtv  * v[pi_3 + 1]; 
          x[pi_3 + 2] += dtv  * v[pi_3 + 2]; 

        }
      }
      pe_put( local_pm.v[ii], v, sizeof(double) * C_SIZE * 3);
      pe_syn();
      pe_put( local_pm.x[ii], x, sizeof(double) * C_SIZE * 3);
      pe_syn();


    }

  }

}
#endif
