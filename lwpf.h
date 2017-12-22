#ifndef __LWPF_SW_PCR
#define __LWPF_SW_PCR
//从核性能计数器0的计数事件选择。无复位值。具体定义如下：
#define PC0_CYCLE              0x00 //5'h00:   周期计数
#define PC0_CNT_COND_JMP       0x01 //5'h01:   执行站台执行的条件转移指令计数
#define PC0_CNT_JMP            0x02 //5'h02:   执行站台执行的转移指令计数
#define PC0_CNT_INST           0x03 //5'h03:   执行指令计数
#define PC0_CNT_PIPE0_INST     0x04 //5'h04:   流水线0发射指令计数
#define PC0_CYC_LD_ST_BLOCK    0x05 //5'h05:   因LD/ST缓冲满而导致流水线停顿的周期计数
#define PC0_CYC_GLQ_BLOCK      0x08 //5'h08:   因GLQ队列满而导致流水线停顿的周期计数
#define PC0_CYC_GSQ_BLOCK      0x09 //5'h09:   因GSQ队列满而导致流水线停顿的周期计数
#define PC0_CYC_SBMD_BLOCK     0x0a //5'h0a:   因SBMD接收缓冲满而导致流水线停顿的周期计数
#define PC0_CYC_CHNL_BLOCK     0x0b //5'h0b:   因通道缓冲满而导致流水线停顿的周期计数
#define PC0_CNT_LDM_READ       0x0c //5'h0c:   执行站台发起的LDM读访问次数计数
#define PC0_CNT_LDM_WRITE      0x0d //5'h0d:   执行站台发起的LDM写访问次数计数
#define PC0_CNT_ICACHE_ACCESS  0x0f //5'h0f:   IBOX对L1 ICache的访问次数计数
#define PC0_CYC_MBOX_LDM_BLOCK 0x10 //5'h10:   MBOX访问LDM的请求被阻塞的周期计数
#define PC0_CNT_JMP_PRED_FAIL  0x11 //5'h11:   执行站台判断出JMP指令预测失败的次数
#define PC0_CNT_2_LDM_SCHED    0x12 //5'h12:   LDM同时仲裁上两个（流水线和核外）访问LDM请求的次数计数
#define PC0_CNT_2DATA_1ACCESS  0x13 //5'h13:   LDM两个数据体一次只有一个被访问的次数计数
#define PC0_CNT_ADD_SUB_MUL    0x14 //5'h14:   浮点加减乘、乘加类指令等效操作计数（标量加减乘+1，标量乘加+2，向量加减乘+4，向量乘加+8）

//从核性能计数器1的计数事件选择。无复位值。具体定义如下：
#define PC1_CYCLE              0x00 //5'h00:   周期计数                                                    
#define PC1_CNT_JMP_FAIL       0x01 //5'h01:   执行站台判断出条件转移指令失败的次数                        
#define PC1_CNT_COND_JMP       0x02 //5'h02:   执行站台执行的条件转移指令计数                              
#define PC1_CNT_NOCOND_JMP     0x03 //5'h03:   执行站台执行的无条件转移指令计数                            
#define PC1_CNT_JMP            0x04 //5'h04:   执行站台执行的JMP指令计数                                   
#define PC1_CNT_PIPE1_INST     0x05 //5'h05:   流水线1发射指令计数                                         
#define PC1_CYC_SYN            0x06 //5'h06:   同步指令导致流水线停顿的周期计数                            
#define PC1_CYC_ACCESS_BLOCK   0x07 //5'h07:   因访存指令未处理完而导致流水线停顿的周期计数                
#define PC1_CYC_PUT_BLOCK      0x08 //5'h08:   因PUT缓冲满导致流水线停顿的周期计数                         
#define PC1_CYC_GETR_BLOCK     0x0a //5'h0a:   因GETR缓冲空导致流水线停顿的周期计数                        
#define PC1_CYC_GETC_BLOCK     0x0b //5'h0b:   因GETC缓冲空导致流水线停顿的周期计数                        
#define PC1_CYC_1INST          0x0c //5'h0c:   发射站不能发射两条指令的周期计数                            
#define PC1_CYC_DATA_REL_BLOCK 0x0d //5'h0d:   数据相关导致流水线停顿的周期计数                            
#define PC1_CYC_SYN_JMP_BLOCK  0x0e //5'h0e:   因同步指令和条件转移指令导致流水线停顿的周期计数            
#define PC1_CNT_ICACHE_MISS    0x0f //5'h0f:   L1 ICache脱靶次数计数                                       
#define PC1_CNT_LDM_ACCESS     0x10 //5'h10:   执行站台发起的LDM访问次数                                   
#define PC1_CNT_TBOX_BLK_MBOX  0x11 //5'h11:   MBOX访问LDM的请求被TBOX普通请求阻塞的周期计数               
#define PC1_CNT_ATOM_BLK_MBOX  0x12 //5'h12:   MBOX访问LDM的请求被原子修改写请求阻塞的周期计数             
#define PC1_CNT_ADDR_BLK_MBOX  0x13 //5'h13:   MBOX访问LDM的请求因为地址相关阻塞的周期计数                 
#define PC1_CNT_2DATA_2ACCESS  0x14 //5'h14:   LDM两个数据体一次被同时访问的次数计数                       
#define PC1_CNT_DIV_SQRT       0x15 //5'h15:   浮点除法、平方根类指令等效操作计数（标量指令+1，向量指令+4）

//TA性能计数器的性能计数事件选择
#define PC2_CNT_CTRL_REQ       0x00 //5’h00:   阵列控制网络请求总数           
#define PC2_CNT_GLD            0x01 //5’h01:   gld请求计数                    
#define PC2_CNT_GST            0x02 //5’h02:   gst请求计数                    
#define PC2_CNT_GF_AND_A       0x03 //5’h03:   gf&a请求计数                   
#define PC2_CNT_GUPDATE        0x04 //5’h04:   gupdate请求计数                
#define PC2_CNT_ICACHE_MISS    0x05 //5’h05:   L1 ICache自主脱靶次数（由于有SBMD模式，该统计值与从核核心的统计值可能有所不同）
#define PC2_CNT_ROW_SYN        0x06 //5’h06:   行同步请求计数                 
#define PC2_CNT_COL_SYN        0x07 //5’h07:   列同步请求计数                 
#define PC2_CNT_USER_INT       0x08 //5’h08:   用户中断请求计数               
#define PC2_CNT_DMA_REQ        0x09 //5’h09:   发起DMA请求计数                
#define PC2_CNT_SBMD_START     0x0a //5’h0a:   sbmd_start命令计数             
#define PC2_CNT_SBMD_END       0x0b //5’h0b:   sbmd_end命令计数               
#define PC2_CNT_SBMD_BREAK     0x0c //5’h0c:   sbmd_break命令计数             
#define PC2_CYC_PIPE_MISS      0x0d //5’h0d:   从核流水线脱靶等待的周期计数   
#define PC2_CYC_SBMD           0x0e //5’h0e:   SBMD状态周期计数               
#define PC2_CYC_SELF           0x0f //5’h0f:   自主运行状态周期计数           
#define PC2_CNT_REG_PUT        0x10 //5’h10:   寄存器通信“PUT”请求包总数      
#define PC2_CNT_REG_PUTR       0x11 //5’h11:   寄存器通信行向“PUT”请求包数    
#define PC2_CNT_REG_PUTC       0x12 //5’h12:   寄存器通信列向“PUT”请求包数    
#define PC2_CNT_FORE_LDM_READ  0x13 //5’h13:   外部对LDM读总数                
#define PC2_CNT_FORE_LDM_WRITE 0x14 //5’h14:   外部对LDM写总数                
#define PC2_CNT_REPLY_INC      0x15 //5’h15:   回答字自增1总数                
#define PC2_CNT_INST_LOAD      0x16 //5’h16:   指令装填总数（16B为单位）      
#define PC2_CNT_INST_LOAD_MISS 0x17 //5’h17:   指令自主脱靶装填总数（16B为单位）
#define PC2_CNT_INST_LOAD_SBMD 0x18 //5’h18:   SBMD指令装填总数（16B为单位）  
#define PC2_CNT_MEM_READ_RESP  0x19 //5’h19:   主存读响应总数                 
#define PC2_CNT_MEM_WRITE_RESP 0x1a //5’h1a:   主存写响应总数                 
#define PC2_CNT_OUTER_READ     0x1b //5’h1b:   外部IO读总数（不含LDM访问）    
#define PC2_CNT_OUTER_WRITE    0x1c //5’h1c:   外部IO写总数（不含LDM访问）    
#define PC2_CNT_REG_LDRC       0x1d //5’h1d:   寄存器通信“ld&发送”请求包总数  
#define PC2_CNT_REG_LDR        0x1e //5’h1e:   寄存器通信行向“ld&发送”请求包数
#define PC2_CNT_REG_LDC        0x1f //5’h1f:   寄存器通信列向“ld&发送”请求包数

#define PCRC_ALL_FLOP          (0x3f << 16)
#define PCRC_ALL_PC            (0x1f)
#define PCRC_ALL               (PCRC_ALL_FLOP | PCRC_ALL_PC)
#define PCRC_VDIV_VSQRT        (1 << 21) //浮点向量除法平方根指令计数使能。本位为1，从核性能计数器1的计数使能为1，且计数事件为5’h15（浮点除法、平方根类指令等效操作计数）时，浮点向量除法平方指令参与计数（每次+4）。
#define PCRC_FDIV_FSQRT        (1 << 20) //浮点标量除法平方根指令计数使能。本位为1，从核性能计数器1的计数使能为1，且计数事件为5’h15（浮点除法、平方根类指令等效操作计数）时，浮点标量除法平方指令参与计数（每次+1）。
#define PCRC_VMA               (1 << 19) //浮点向量乘加类指令计数使能。本位为1，从核性能计数器0的计数使能为1，且计数事件为5’h14（浮点加减乘、乘加类指令等效操作计数）时，浮点向量乘加类指令参与计数（每次+8）。
#define PCRC_VADD_VSUB_VMUL    (1 << 18) //浮点向量加减乘指令计数使能。本位为1，从核性能计数器0的计数使能为1，且计数事件为5’h14（浮点加减乘、乘加类指令等效操作计数）时，浮点向量加减乘指令参与计数（每次+4）。
#define PCRC_FMA               (1 << 17) //浮点标量乘加类指令计数使能。本位为1，从核性能计数器0的计数使能为1，且计数事件为5’h14（浮点加减乘、乘加类指令等效操作计数）时，浮点标量乘加类指令参与计数（每次+2）。
#define PCRC_FADD_FSUB_FMUL    (1 << 16) //浮点标量加减乘指令计数使能。本位为1，从核性能计数器0的计数使能为1，且计数事件为5’h14（浮点加减乘、乘加类指令等效操作计数）时，浮点标量加减乘指令参与计数（每次+1）。
#define PCRC_PC2_OVERFLOW      (1 <<  5) //TA性能计数器的溢出中断使能。为1时使能。   
#define PCRC_PC2               (1 <<  4) //TA性能计数器的计数使能。为1时使能。        
#define PCRC_PC1_OVERFLOW      (1 <<  3) //从核性能计数器1的溢出中断使能。为1时使能。
#define PCRC_PC1               (1 <<  2) //从核性能计数器1的计数使能。为1时使能。    
#define PCRC_PC0_OVERFLOW      (1 <<  1) //从核性能计数器0的溢出中断使能。为1时使能。
#define PCRC_PC0               (1 <<  0) //从核性能计数器0的计数使能。为1时使能。

typedef struct perf_config_t{
  long pcr0, pcr1, pcr2, pcrc;
} perf_config_t;

static inline void set_perf_mode(long pcr0, long pcr1, long pcr2, long pcrc){
  pcr0 <<= 59;
  pcr1 <<= 59;
  pcr2 <<= 59;

  asm volatile("wcsr %0, 5\n"
               "wcsr %1, 6\n"
               "wcsr %2, 7\n"
               "wcsr %3, 8\n": :
               "r"(pcr0), "r"(pcr1), "r"(pcr2), "r"(pcrc));
}
static inline long pcr0(){
  long pcr;
  long msk = (1LL << 59) - 1;
  asm volatile("rcsr %0, 5" : "=r"(pcr));
  return pcr & msk;
}
static inline long pcr1(){
  long pcr;
  long msk = (1LL << 59) - 1;
  asm volatile("rcsr %0, 6" : "=r"(pcr));
  return pcr & msk;
}
static inline long pcr2(){
  long pcr;
  long msk = (1LL << 59) - 1;
  asm volatile("rcsr %0, 7" : "=r"(pcr));
  return pcr & msk;
}
static inline long rpcc(){
  long rpc;
  asm volatile("rcsr %0, 4" : "=r"(rpc));
  return rpc;
}
#endif
#ifdef CPE
#ifndef __LWPF_CPE_H
#define __LWPF_CPE_H

#include <slave.h>

//define a enum of kernels
#define _K(x) x
#define K(x) , x
typedef enum {
  LWPF_KERNELS
} lwpf_kernel;
#undef _K
#undef K

//define an array of kernel names
#define _K(x) #x
#define K(x) , #x
#define U(x) lwpf_kernel_names_ ## x
char *LWPF_UNIT[] = {LWPF_KERNELS};
#undef _K
#undef K
#undef U
//define arrays of kernels
#define _K(x) 1
#define K(x) +1
#define U(x) lwpf_kernel_count_ ## x
const long LWPF_UNIT = LWPF_KERNELS;
__thread_local static long KCNT = LWPF_KERNELS;
#undef U
__thread_local static long __rpcc[LWPF_KERNELS];
__thread_local static long __pcr0[LWPF_KERNELS];
__thread_local static long __pcr1[LWPF_KERNELS];
__thread_local static long __pcr2[LWPF_KERNELS];
#undef _K
#undef K

/* static inline void lwpf_start(lwpf_kernel kernel){ */
/*   __rpcc[kernel] -= rpcc(); */
/*   __pcr0[kernel] -= pcr0(); */
/*   __pcr1[kernel] -= pcr1(); */
/*   __pcr2[kernel] -= pcr2(); */
/* } */

/* static inline void lwpf_stop(lwpf_kernel kernel){ */
/*   __rpcc[kernel] += rpcc(); */
/*   __pcr0[kernel] += pcr0(); */
/*   __pcr1[kernel] += pcr1(); */
/*   __pcr2[kernel] += pcr2(); */
/* } */
#define lwpf_start(kernel) {                                    \
    long rpc, pc0, pc1, pc2;                                    \
    asm volatile("rcsr %0, 4\n"                                 \
                 "rcsr %1, 5\n"                                 \
                 "rcsr %2, 6\n"                                 \
                 "rcsr %3, 7\n"                                 \
                 : "=r"(rpc), "=r"(pc0), "=r"(pc1), "=r"(pc2)); \
    __rpcc[kernel] -= rpc;                                      \
    __pcr0[kernel] -= pc0;                                      \
    __pcr1[kernel] -= pc1;                                      \
    __pcr2[kernel] -= pc2;                                      \
  }                                                             \

#define lwpf_stop(kernel) {                                     \
    long rpc, pc0, pc1, pc2;                                    \
    asm volatile("rcsr %0, 4\n"                                 \
                 "rcsr %1, 5\n"                                 \
                 "rcsr %2, 6\n"                                 \
                 "rcsr %3, 7\n"                                 \
                 : "=r"(rpc), "=r"(pc0), "=r"(pc1), "=r"(pc2)); \
    __rpcc[kernel] += rpc;                                      \
    __pcr0[kernel] += pc0;                                      \
    __pcr1[kernel] += pc1;                                      \
    __pcr2[kernel] += pc2;                                      \
  }

#define U(x) lwpf_init_ ## x
void LWPF_UNIT(perf_config_t *conf){
  int i;
  for (i = 0; i < KCNT; i ++){
    __rpcc[i] = 0;
    __pcr0[i] = 0;
    __pcr1[i] = 0;
    __pcr2[i] = 0;
  }  
  set_perf_mode(conf->pcr0, conf->pcr1, conf->pcr2, conf->pcrc);
}
#undef U
#define U(x) lwpf_finish_ ## x
void LWPF_UNIT(long *mem){
  volatile int reply = 0;
  long *mrpcc = mem;
  long *mpcr0 = mem + KCNT * 64;
  long *mpcr1 = mem + KCNT * 2 * 64;
  long *mpcr2 = mem + KCNT * 3 * 64;
  //if (_MYID == 0)
  //  printf("%d\n", KCNT);
  reply = 0;
  athread_put(PE_MODE, __rpcc, mrpcc + KCNT * _MYID, sizeof(long) * KCNT, (void *)&reply, 0, 0);
  athread_put(PE_MODE, __pcr0, mpcr0 + KCNT * _MYID, sizeof(long) * KCNT, (void *)&reply, 0, 0);
  athread_put(PE_MODE, __pcr1, mpcr1 + KCNT * _MYID, sizeof(long) * KCNT, (void *)&reply, 0, 0);
  athread_put(PE_MODE, __pcr2, mpcr2 + KCNT * _MYID, sizeof(long) * KCNT, (void *)&reply, 0, 0);
  while (reply != 4);
}
#undef U
#endif
#endif
#ifdef MPE
#ifndef __LWPF_MPE_H
#define __LWPF_MPE_H
#include <athread.h>

//extern kernel names
#define U(x) extern char *lwpf_kernel_names_ ## x[];
LWPF_UNITS
#undef U

//extern kernel counts
#define U(x) extern long lwpf_kernel_count_ ## x;
LWPF_UNITS
#undef U

//extern init functions
#define U(x) extern SLAVE_FUN(lwpf_init_ ## x)(void *);
LWPF_UNITS
#undef U

//extern finish functions
#define U(x) extern SLAVE_FUN(lwpf_finish_ ## x)();
LWPF_UNITS
#undef U

//create init function
#define U(x) athread_spawn(lwpf_init_ ## x, conf); athread_join();
static inline void lwpf_init(perf_config_t *conf){
  LWPF_UNITS
}
#undef U
inline int longlen(long in){
  if (in == 0)
    return 1;
  int p = in, ret = 0;
  while (p > 0){
    p /= 10;
    ret ++;
  }
  return ret;
}
inline void calc_unit(long max, char *unit, long *div){
  *unit = '1';
  *div = 1;
  if (max > 1000000) {
    *unit = 'K';
    *div = 1000;
  }
  if (max > 1000000000) {
    *unit = 'M';
    *div = 10000000;
  }
  if (max > 1000000000000) {
    *unit = 'G';
    *div = 10000000000;
  }
  if (max > 1000000000000000) {
    *unit = 'T';
    *div = 10000000000000;
  }
}
inline void lwpf_print(FILE *f, long *ctrs, int kcnt, const char *unit, char *kname[]){
  long *crpcc = ctrs;
  long *cpcr0 = ctrs + kcnt * 64;
  long *cpcr1 = ctrs + kcnt * 2 * 64;
  long *cpcr2 = ctrs + kcnt * 3 * 64;

  fprintf(f, "counter staticists for %s\n", unit);
  long rpcc_max = 0;
  long pcr0_max = 0;
  long pcr1_max = 0;
  long pcr2_max = 0;

  long *rpcc_maxk = malloc(sizeof(long) * kcnt);
  long *pcr0_maxk = malloc(sizeof(long) * kcnt);
  long *pcr1_maxk = malloc(sizeof(long) * kcnt);
  long *pcr2_maxk = malloc(sizeof(long) * kcnt);

  long *rpcc_mink = malloc(sizeof(long) * kcnt);
  long *pcr0_mink = malloc(sizeof(long) * kcnt);
  long *pcr1_mink = malloc(sizeof(long) * kcnt);
  long *pcr2_mink = malloc(sizeof(long) * kcnt);

  long *rpcc_avgk = malloc(sizeof(long) * kcnt);
  long *pcr0_avgk = malloc(sizeof(long) * kcnt);
  long *pcr1_avgk = malloc(sizeof(long) * kcnt);
  long *pcr2_avgk = malloc(sizeof(long) * kcnt);

  //printf("%d\n", kcnt);
  int i, j;
  for (i = 0; i < kcnt; i ++){
    rpcc_maxk[i] = 0; rpcc_mink[i] = 1L << 62; rpcc_avgk[i] = 0;
    pcr0_maxk[i] = 0; pcr0_mink[i] = 1L << 62; pcr0_avgk[i] = 0;
    pcr1_maxk[i] = 0; pcr1_mink[i] = 1L << 62; pcr1_avgk[i] = 0;
    pcr2_maxk[i] = 0; pcr2_mink[i] = 1L << 62; pcr2_avgk[i] = 0;

    for (j = 0; j < 64; j ++){
      int ic = j * kcnt + i;
      if (crpcc[ic] > rpcc_maxk[i]) rpcc_maxk[i] = crpcc[ic];
      if (cpcr0[ic] > pcr0_maxk[i]) pcr0_maxk[i] = cpcr0[ic];
      if (cpcr1[ic] > pcr1_maxk[i]) pcr1_maxk[i] = cpcr1[ic];
      if (cpcr2[ic] > pcr2_maxk[i]) pcr2_maxk[i] = cpcr2[ic];
      if (crpcc[ic] < rpcc_mink[i]) rpcc_mink[i] = crpcc[ic];
      if (cpcr0[ic] < pcr0_mink[i]) pcr0_mink[i] = cpcr0[ic];
      if (cpcr1[ic] < pcr1_mink[i]) pcr1_mink[i] = cpcr1[ic];
      if (cpcr2[ic] < pcr2_mink[i]) pcr2_mink[i] = cpcr2[ic];
      rpcc_avgk[i] += crpcc[ic];
      pcr0_avgk[i] += cpcr0[ic];
      pcr1_avgk[i] += cpcr1[ic];
      pcr2_avgk[i] += cpcr2[ic];
    }
    rpcc_avgk[i] = (rpcc_avgk[i] + 32) / 64;
    pcr0_avgk[i] = (pcr0_avgk[i] + 32) / 64;
    pcr1_avgk[i] = (pcr1_avgk[i] + 32) / 64;
    pcr2_avgk[i] = (pcr2_avgk[i] + 32) / 64;

    if (rpcc_maxk[i] > rpcc_max) rpcc_max = rpcc_maxk[i];
    if (pcr0_maxk[i] > pcr0_max) pcr0_max = pcr0_maxk[i];
    if (pcr1_maxk[i] > pcr1_max) pcr1_max = pcr1_maxk[i];
    if (pcr2_maxk[i] > pcr2_max) pcr2_max = pcr2_maxk[i];
  }

  char rpcc_unit, pcr0_unit, pcr1_unit, pcr2_unit;
  long rpcc_div, pcr0_div, pcr1_div, pcr2_div;
  calc_unit(rpcc_max, &rpcc_unit, &rpcc_div);
  calc_unit(pcr0_max, &pcr0_unit, &pcr0_div);
  calc_unit(pcr1_max, &pcr1_unit, &pcr1_div);
  calc_unit(pcr2_max, &pcr2_unit, &pcr2_div);

  char *fmth = "|%16s|%18s(%c)|%18s(%c)|%18s(%c)|%18s(%c)|\n";
  char *fmtl = "|%16s|%7d%7d%7d|%7d%7d%7d|%7d%7d%7d|%7d%7d%7d|\n";
  //puts(fmth);
  //puts(fmtl);
  fprintf(f, fmth, "kernel", "rpcc(avg/max/min)", rpcc_unit, "pcr0(avg/max/min)", pcr0_unit, "pcr1(avg/max/min)", pcr1_unit, "pcr2(avg/max/min)", pcr2_unit);
  char hz[128];
  memset(hz, '-', sizeof(hz));
  hz[21] = 0;
  fprintf(f, "+----------------+%s+%s+%s+%s+\n", hz, hz, hz, hz);
  for (i = 0; i < kcnt; i ++){
    fprintf(f, fmtl, kname[i]
            , rpcc_avgk[i] / rpcc_div, rpcc_maxk[i] / rpcc_div, rpcc_mink[i] / rpcc_div
            , pcr0_avgk[i] / pcr0_div, pcr0_maxk[i] / pcr0_div, pcr0_mink[i] / pcr0_div
            , pcr1_avgk[i] / pcr1_div, pcr1_maxk[i] / pcr1_div, pcr1_mink[i] / pcr1_div
            , pcr2_avgk[i] / pcr2_div, pcr2_maxk[i] / pcr2_div, pcr2_mink[i] / pcr2_div);
  }

  free(rpcc_maxk);free(rpcc_mink);free(rpcc_avgk);
  free(pcr0_maxk);free(pcr0_mink);free(pcr0_avgk);
  free(pcr1_maxk);free(pcr1_mink);free(pcr1_avgk);
  free(pcr2_maxk);free(pcr2_mink);free(pcr2_avgk);
}
#define U(x)                                                            \
  ctrs = malloc(sizeof(long) * lwpf_kernel_count_ ## x * 64 * 4);       \
  athread_spawn(lwpf_finish_ ##x, ctrs);                                \
  athread_join();                                                       \
  lwpf_print(f, ctrs, lwpf_kernel_count_ ## x, #x, lwpf_kernel_names_ ##x); \
  free(ctrs);                                                           \

static inline void lwpf_finish(FILE *f){
  long *ctrs;
  LWPF_UNITS;
  puts("finished");
}
#undef U
#endif
#endif
