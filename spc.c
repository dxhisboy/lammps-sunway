#ifdef MPE
#include <athread.h>
#include <signal.h>
#include "sw5_pcr.h"
extern SLAVE_FUN(gid)(long*);
long cgid;
void cinfo_init(){
  if (!athread_idle())
    athread_init();
  athread_spawn(gid, &cgid);
  athread_join();
}
#define CPE_PC_ADDR(cpe_id)   (0x8003002000 | (cgid << 36) | ((cpe_id) << 16))
#define CPE_PCRC_ADDR(cpe_id) (0x8003000480 | (cgid << 36) | ((cpe_id) << 16))
#define CPE_PCR0_ADDR(cpe_id) (0x8003000380 | (cgid << 36) | ((cpe_id) << 16))
#define CPE_PCR1_ADDR(cpe_id) (0x8003000400 | (cgid << 36) | ((cpe_id) << 16))
#define CPE_TA_ADDR(cpe_id)   (0x800300a180 | (cgid << 36) | ((cpe_id) << 16))
inline long set_cpe_pcrc(long event){
  int i;
  for (i = 0; i < 64; i ++){
    *(long*)CPE_PCRC_ADDR(i) = event;
  }
}
inline long set_cpe_pcr0(long event){
  int i;
  for (i = 0; i < 64; i ++){
    *(long*)CPE_PCR0_ADDR(i) = event << 59;
  }
}
inline long set_cpe_pcr1(long event){
  int i;
  for (i = 0; i < 64; i ++){
    *(long*)CPE_PCR1_ADDR(i) = event << 59;
  }
}

inline long get_cpe_pcrc(long id){
  long pcr = *(long*)CPE_PCRC_ADDR(id);
  return pcr & ((1L << 59) - 1);
}
inline long get_cpe_pcr0(long id){
  long pcr = *(long*)CPE_PCR0_ADDR(id);
  return pcr & ((1L << 59) - 1);
}
inline long get_cpe_pcr1(long id){
  long pcr = *(long*)CPE_PCR1_ADDR(id);
  return pcr & ((1L << 59) - 1);
}

inline long get_cpe_pc(int cpe_id){
  //athread_task_info();
  return *(long*)(549806153728L + cgid * 68719476736L + cpe_id * 65536L);
}

void pc_on_sig(int sig){
  printf("writing cpe PCs on signal %d\n", sig);
  int i, j;
  for (i = 0; i < 8; i ++) {
    for (j = 0; j < 8; j ++)
      printf("%3d: %12llx ", i * 8 + j, get_cpe_pc(i * 8 + j));
    puts("");
  }
}
#endif
#ifdef CPE
#include <slave.h>
void gid(long *id){
  int cid;
  asm ("rcsr %0, 0\n\t": "=r"(cid));
  if (_MYID == 0)
    *id = cid >> 6;
}
#endif
