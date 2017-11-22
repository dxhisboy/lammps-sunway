#include <athread.h>
#include <signal.h>
inline long get_cpe_pc(int cpe_id){
  athread_task_info();
  return *(long*)(549806153728L + __cgid * 68719476736L + cpe_id * 65536L);
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


