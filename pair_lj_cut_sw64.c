#include "sunway.h"
#include "gptl.h"
#include <stdlib.h>
#include <simd.h>

#ifdef MPE
#define LWPF_UNITS U(LJCUT)
#include "lwpf.h"
int r = 0;

extern SLAVE_FUN(pair_lj_cut_sunway_compute_para)(compute_param_t *);
extern SLAVE_FUN(pair_lj_cut_sunway_compute_para_vec)(compute_param_t *);
extern SLAVE_FUN(pair_lj_cut_sunway_compute_a2s)(compute_param_t *);

void pair_lj_cut_sunway_compute(compute_param_t *pm){
  if (athread_idle() == 0)
    athread_init();
  if (r == 0){
    perf_config_t conf;
    conf.pcr0 = PC0_CYC_CHNL_BLOCK;
    conf.pcr1 = PC1_CYCLE;
    conf.pcr2 = PC2_CNT_DMA_REQ;
    conf.pcrc = PCRC_ALL_PC;
    lwpf_init(&conf);
  }
  atom_in_t *atom_in = malloc(sizeof(atom_in_t) * (pm->nlocal + pm->nghost + 32));
  //printf("e:%d ee:%d ve:%d %d %d %d %d\n", pm->eflag_either, pm->vflag_either, pm->eflag_atom, pm->vflag_atom, pm->eflag_global, pm->vflag_global);
  pm->atom_in = (void*)((long)atom_in | 255) + 1;
  GPTLstart("a2s");
  athread_spawn(pair_lj_cut_sunway_compute_a2s, pm);
  athread_join();
  GPTLstop("a2s");
  GPTLstart("comp");
  athread_spawn(pair_lj_cut_sunway_compute_para_vec, pm);
  athread_join();
  /* printf("%e %e %e\n", pm->eng_coul, pm->eng_vdwl, pm->virial[0]); */
  /* printf("%e %e\n", pm->x[0][0], pm->f[0][0]); */
  GPTLstop("comp");
  if (r == 10 && pm->rank == 0){
    lwpf_finish(stdout);
  }
  r ++;
  free(atom_in);
  //exit(1);
}
#endif
#ifdef CPE
#define LWPF_KERNELS _K(COMP) K(CACHE) K(IREAD) K(IWRITE) K(ALL) K(GST) K(JLOOP)
#define LWPF_UNIT U(LJCUT)
#include "lwpf.h"
#include <dma.h>
#include <math.h>
/* #define JPAGE_SIZE 64 */
/* #define JLIST_SIZE 256 */
/* #define IPAGE_SIZE 16 */
inline void ev_tally_full(compute_param_t *pm, int i, double evdwl, double ecoul, double fpair,
                   double delx, double dely, double delz,
                   double *eng_coul, double *eng_vdwl, double *virial,
                   double *eatom, double (*vatom)[6])
{
  double v[6];

  if (pm->eflag_either) {
    if (pm->eflag_global) {
      eng_vdwl[0] += 0.5*evdwl;
      eng_coul[0] += 0.5*ecoul;
    }
    if (pm->eflag_atom) eatom[i] += 0.5 * (evdwl + ecoul);
  }

  if (pm->vflag_either) {
    v[0] = 0.5*delx*delx*fpair;
    v[1] = 0.5*dely*dely*fpair;
    v[2] = 0.5*delz*delz*fpair;
    v[3] = 0.5*delx*dely*fpair;
    v[4] = 0.5*delx*delz*fpair;
    v[5] = 0.5*dely*delz*fpair;

    if (pm->vflag_global) {
      virial[0] += v[0];
      virial[1] += v[1];
      virial[2] += v[2];
      virial[3] += v[3];
      virial[4] += v[4];
      virial[5] += v[5];
    }

    if (pm->vflag_atom) {
      vatom[i][0] += v[0];
      vatom[i][1] += v[1];
      vatom[i][2] += v[2];
      vatom[i][3] += v[3];
      vatom[i][4] += v[4];
      vatom[i][5] += v[5];
    }
  }
}

static inline void reg_reduce_inplace_doublev4(doublev4 *arr, int len){
  int i, j;
  doublev4 tmp;
  for (i = 1; i < 8; i += i){
    if ((_ROW & i) == i){
      for (j = 0; j < len; j ++)
        asm("putc %0, %1": : "r"(arr[j]), "r"(_ROW ^ i));
    }
    if ((_ROW & i) == 0){
      for (j = 0; j < len; j ++){
        asm("getc %0\n" : "=r"(tmp));
        arr[j] += tmp;
      }
    }
    athread_syn(COL_SCOPE, 0xff);
  }
  athread_syn(ARRAY_SCOPE, 0xffff);
  if (_ROW == 0){
    for (i = 1; i < 8; i += i){
      if ((_COL & i) == i){
        for (j = 0; j < len; j ++)
          asm("putr %0, %1": : "r"(arr[j]), "r"(_COL ^ i));
      }
      if ((_COL & i) == 0){
        for (j = 0; j < len; j ++){
          asm("getr %0\n" : "=r"(tmp));
          arr[j] += tmp;
        }
      }
    }
    athread_syn(ROW_SCOPE, 0xff);
  }
}
#define transpose4x4(in0, in1, in2, in3, ot0, ot1, ot2, ot3) { \
    doublev4 o0 = simd_vshff(in1,in0,0x44);                     \
    doublev4 o1 = simd_vshff(in1,in0,0xEE);                     \
    doublev4 o2 = simd_vshff(in3,in2,0x44);                     \
    doublev4 o3 = simd_vshff(in3,in2,0xEE);                     \
    ot0 = simd_vshff(o2,o0,0x88);                               \
    ot1 = simd_vshff(o2,o0,0xDD);                               \
    ot2 = simd_vshff(o3,o1,0x88);                               \
    ot3 = simd_vshff(o3,o1,0xDD);                               \
  }
#define vshuffd_rc(a, b, c, d) (d | (c << 2) | (b << 4) | (a << 6))
#define simd_vsumd(x) {                                 \
    x += simd_vshff(x, x, vshuffd_rc(2, 3, 0, 1));      \
    x += simd_vshff(x, x, vshuffd_rc(1, 0, 3, 2));      \
  } 
#define JLIST_PAGESIZE 64
#define ILIST_PAGESIZE 128

#define JCACHE_HBIT 10
#define JCACHE_SBIT 3
#define JCACHE_LINESIZE (1 << JCACHE_SBIT)
#define JCACHE_LINECNT (1 << (JCACHE_HBIT - JCACHE_SBIT))
#define JCACHE_MMASK (JCACHE_LINESIZE - 1)
#define JCACHE_LMASK (JCACHE_LINECNT - 1)

void pair_lj_cut_sunway_compute_a2s(compute_param_t *pm){
  pe_init();
  compute_param_t l_pm;
  pe_get(pm, &l_pm, sizeof(compute_param_t));
  pe_syn();
  double x[ILIST_PAGESIZE][3];
  int type[ILIST_PAGESIZE];
  atom_in_t atom_in[ILIST_PAGESIZE];
  int ipage_start, i;
  int inum = l_pm.nlocal + l_pm.nghost;
  for (ipage_start = ILIST_PAGESIZE * _MYID; ipage_start < inum; ipage_start += ILIST_PAGESIZE * 64){
    int ipage_end = ipage_start + ILIST_PAGESIZE;
    if (ipage_end > inum)
      ipage_end = inum;
    int ipage_size = ipage_end - ipage_start;
    pe_get(l_pm.x[ipage_start], x, ipage_size * sizeof(double) * 3);
    pe_get(l_pm.type + ipage_start, type, ipage_size * sizeof(int));
    pe_syn();
    for (i = 0; i < ipage_size; i ++){
      atom_in[i].x[0] = x[i][0];
      atom_in[i].x[1] = x[i][1];
      atom_in[i].x[2] = x[i][2];
      atom_in[i].type = type[i];
    }
    pe_put(l_pm.atom_in + ipage_start, atom_in, sizeof(atom_in_t) * ipage_size);
    pe_syn();
  }
}

void pair_lj_cut_sunway_compute_para(compute_param_t *pm){
  //lwpf_start(ALL);
  pe_init();
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,r2inv,r6inv,forcelj,factor_lj;
  int *ilist,*jlist,*numneigh,**firstneigh;

  double (*x)[3] = pm->x;
  double (*f)[3] = pm->f;
  //double xi[IPAGE_SIZE][3], fi[IPAGE_SIZE][3];
  //double ti[IPAGE_SIZE];
  int *type = pm->type;
  int nlocal = pm->nlocal;
  double *special_lj;
  doublev4 eng_virial[2];
  double *eng_vdwl = (double*)(void*)eng_virial;
  double *eng_coul = eng_vdwl + 1;
  double *virial = eng_coul + 1;

  compute_param_t l_pm;
  pe_get(pm, &l_pm, sizeof(compute_param_t));
  double INF = 1.0 / 0.0;
  pe_syn();
  int ntp1 = l_pm.ntypes + 1;
  double cutsq[ntp1][ntp1], lj1[ntp1][ntp1], lj2[ntp1][ntp1], lj3[ntp1][ntp1], lj4[ntp1][ntp1], offset[ntp1][ntp1];
  pe_get(l_pm.cutsq, cutsq, (l_pm.ntypes + 1) * (l_pm.ntypes + 1) * sizeof(double));
  pe_get(l_pm.lj1, lj1, (l_pm.ntypes + 1) * (l_pm.ntypes + 1) * sizeof(double));
  pe_get(l_pm.lj2, lj2, (l_pm.ntypes + 1) * (l_pm.ntypes + 1) * sizeof(double));
  pe_get(l_pm.lj3, lj3, (l_pm.ntypes + 1) * (l_pm.ntypes + 1) * sizeof(double));
  pe_get(l_pm.lj4, lj4, (l_pm.ntypes + 1) * (l_pm.ntypes + 1) * sizeof(double));
  pe_get(l_pm.offset, offset, (l_pm.ntypes + 1) * (l_pm.ntypes + 1) * sizeof(double));
  special_lj = l_pm.special_lj;
  pe_syn();

  eng_virial[0] = 0;
  eng_virial[1] = eng_virial[0];
  inum = pm->inum;
  numneigh = pm->numneigh;
  firstneigh = pm->firstneigh;
  int jlist_buf[JLIST_PAGESIZE];
  atom_in_t jlist_atom[JLIST_PAGESIZE];
  atom_in_t j_cache[JCACHE_LINECNT][JCACHE_LINESIZE];
  int jcache_tag[JCACHE_LINECNT];
  double *xjc, xj[3];
  int tj;
  double xi[ILIST_PAGESIZE][3], fi[ILIST_PAGESIZE][3];
  double ei[ILIST_PAGESIZE], vi[ILIST_PAGESIZE][6], v[6];
  int ti[ILIST_PAGESIZE], *fn[ILIST_PAGESIZE], nn[ILIST_PAGESIZE];
  int ipage_start;
  volatile int cache_reply;
  dma_desc cache_get_desc = 0;
  //memset(&cache_get_desc, 0, sizeof(dma_desc));
  dma_set_mode(&cache_get_desc, PE_MODE);
  dma_set_size(&cache_get_desc, sizeof(atom_in_t) * JCACHE_LINESIZE);
  dma_set_op(&cache_get_desc, DMA_GET);
  dma_set_reply(&cache_get_desc, &cache_reply);
  lwpf_start(COMP);
  for (i = 0; i < JCACHE_LINECNT; i ++)
    jcache_tag[i] = -1;
  for (ipage_start = ILIST_PAGESIZE * _MYID; ipage_start < inum; ipage_start += ILIST_PAGESIZE * 64){
    int ipage_end = ipage_start + ILIST_PAGESIZE;
    if (ipage_end > inum)
      ipage_end = inum;
    int ipage_size = ipage_end - ipage_start;
    lwpf_start(IREAD);
    pe_get(x[ipage_start], xi, ipage_size * sizeof(double) * 3);
    pe_get(type + ipage_start, ti, ipage_size * sizeof(int));
    pe_get(firstneigh + ipage_start, fn, ipage_size * sizeof(int*));
    pe_get(numneigh + ipage_start, nn, ipage_size * sizeof(int));

    pe_syn();
    lwpf_stop(IREAD);
    for (ii = ipage_start; ii < ipage_end; ii ++) {
      i = ii;
      int ioff = i - ipage_start;
      fi[ioff][0] = fi[ioff][1] = fi[ioff][2] = 0;
      ei[ioff] = 0;
      vi[ioff][0] = vi[ioff][1] = vi[ioff][2] = 0;
      vi[ioff][3] = vi[ioff][4] = vi[ioff][5] = 0;
      itype = ti[ioff];
      jlist = fn[ioff];
      jnum = nn[ioff];
      int jpage_start;
      for (jpage_start = 0; jpage_start < jnum; jpage_start += JLIST_PAGESIZE){
        int jpage_size = JLIST_PAGESIZE;
        if (JLIST_PAGESIZE + jpage_start > jnum)
          jpage_size = jnum - jpage_start;
        pe_get(jlist + jpage_start, jlist_buf, jpage_size * sizeof(int));
        pe_syn();
        
        lwpf_start(CACHE);
        for (jj = 0; jj < jpage_size; jj++) {
          j = jlist_buf[jj];
          factor_lj = special_lj[sbmask(j)];
          j &= NEIGHMASK;
          int line = j >> JCACHE_SBIT & JCACHE_LMASK;
          int tag = j >> JCACHE_HBIT;

          if (jcache_tag[line] != tag){
            int mem = j & ~JCACHE_MMASK;
            cache_reply = 0;
            dma(cache_get_desc, l_pm.atom_in + mem, j_cache[line]);
            while (cache_reply != 1);
            /* pe_get(l_pm.atom_in + mem, j_cache[line], sizeof(atom_in_t) * JCACHE_LINESIZE); */
            /* pe_syn(); */
            jcache_tag[line] = tag;
          }
          jlist_atom[jj] = j_cache[line][j & JCACHE_MMASK];
        }
        lwpf_stop(CACHE);
        for (jj = jpage_size; jj < ((jpage_size + 3) & (~3)); jj ++){
          jlist_atom[jj].x[0] = INF;
        }
        lwpf_start(JLOOP);
        for (jj = 0; jj < jpage_size; jj ++){
          jtype = jlist_atom[jj].type;

          delx = xi[ioff][0] - jlist_atom[jj].x[0];
          dely = xi[ioff][1] - jlist_atom[jj].x[1];
          delz = xi[ioff][2] - jlist_atom[jj].x[2];
          rsq = delx*delx + dely*dely + delz*delz;

          if (rsq < cutsq[itype][jtype]) {
            r2inv = 1.0/rsq;
            r6inv = r2inv*r2inv*r2inv;
            forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
            fpair = factor_lj*forcelj*r2inv;

            fi[ioff][0] += delx*fpair;
            fi[ioff][1] += dely*fpair;
            fi[ioff][2] += delz*fpair;
            if (l_pm.evflag) {
              if (l_pm.eflag_either) {
                evdwl = r6inv*(lj3[itype][jtype]*r6inv-lj4[itype][jtype]) -
                  offset[itype][jtype];
                evdwl *= factor_lj;

                if (l_pm.eflag_global) {
                  eng_vdwl[0] += 0.5*evdwl;
                }
                if (l_pm.eflag_atom) ei[i] += 0.5 * evdwl;
              }

              if (l_pm.vflag_either) {
                v[0] = 0.5*delx*delx*fpair;
                v[1] = 0.5*dely*dely*fpair;
                v[2] = 0.5*delz*delz*fpair;
                v[3] = 0.5*delx*dely*fpair;
                v[4] = 0.5*delx*delz*fpair;
                v[5] = 0.5*dely*delz*fpair;

                if (l_pm.vflag_global) {
                  virial[0] += v[0];
                  virial[1] += v[1];
                  virial[2] += v[2];
                  virial[3] += v[3];
                  virial[4] += v[4];
                  virial[5] += v[5];
                }

                if (l_pm.vflag_atom) {
                  vi[i][0] += v[0];
                  vi[i][1] += v[1];
                  vi[i][2] += v[2];
                  vi[i][3] += v[3];
                  vi[i][4] += v[4];
                  vi[i][5] += v[5];
                }
              }

            }
          }
        }
        lwpf_stop(JLOOP);
      }
    }
    lwpf_start(IWRITE);
    pe_put(f[ipage_start], fi[0], sizeof(double) * 3 * ipage_size);
    if (l_pm.eflag_either && l_pm.eflag_atom)
      pe_put(l_pm.eatom + ipage_start, ei, sizeof(double) * ipage_size);
    if (l_pm.vflag_either && l_pm.vflag_atom)
      pe_put(l_pm.vatom[ipage_start], vi[0], sizeof(double) * 6 * ipage_size);
    pe_syn();
    lwpf_stop(IWRITE);
  }
  lwpf_stop(COMP);
  reg_reduce_inplace_doublev4(eng_virial, 2);
  lwpf_start(GST);
  if (_MYID == 0){
    pe_put(&(pm->eng_vdwl), eng_vdwl, sizeof(double) * 8);
    pe_syn();
  }
  lwpf_stop(GST);
  lwpf_stop(ALL);
}

void pair_lj_cut_sunway_compute_para_vec(compute_param_t *pm){
  lwpf_start(ALL);
  pe_init();
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,r2inv,r6inv,forcelj,factor_lj;
  int *ilist,*jlist,*numneigh,**firstneigh;

  double (*x)[3] = pm->x;
  double (*f)[3] = pm->f;
  //double xi[IPAGE_SIZE][3], fi[IPAGE_SIZE][3];
  //double ti[IPAGE_SIZE];
  int *type = pm->type;
  int nlocal = pm->nlocal;
  double *special_lj;
  doublev4 eng_virial[2];
  double *eng_vdwl = (double*)(void*)eng_virial;
  double *eng_coul = eng_vdwl + 1;
  double *virial = eng_coul + 1;

  compute_param_t l_pm;
  pe_get(pm, &l_pm, sizeof(compute_param_t));
  double INF = 1.0 / 0.0;
  pe_syn();
  int ntp1 = l_pm.ntypes + 1;
  double cutsq[ntp1][ntp1], lj1[ntp1][ntp1], lj2[ntp1][ntp1], lj3[ntp1][ntp1], lj4[ntp1][ntp1], offset[ntp1][ntp1];
  pe_get(l_pm.cutsq, cutsq, (l_pm.ntypes + 1) * (l_pm.ntypes + 1) * sizeof(double));
  pe_get(l_pm.lj1, lj1, (l_pm.ntypes + 1) * (l_pm.ntypes + 1) * sizeof(double));
  pe_get(l_pm.lj2, lj2, (l_pm.ntypes + 1) * (l_pm.ntypes + 1) * sizeof(double));
  pe_get(l_pm.lj3, lj3, (l_pm.ntypes + 1) * (l_pm.ntypes + 1) * sizeof(double));
  pe_get(l_pm.lj4, lj4, (l_pm.ntypes + 1) * (l_pm.ntypes + 1) * sizeof(double));
  pe_get(l_pm.offset, offset, (l_pm.ntypes + 1) * (l_pm.ntypes + 1) * sizeof(double));
  special_lj = l_pm.special_lj;
  pe_syn();

  doublev4 lj_v4[4][4], csqoff_v4[4][4], slj_v4;
  slj_v4 = simd_set_doublev4(special_lj[0], special_lj[1], special_lj[2], special_lj[3]);
  for (i = 0; i < ntp1; i ++)
    for (j = 0; j < ntp1; j ++){
      lj_v4[i][j] = simd_set_doublev4(lj1[i][j], lj2[i][j], lj3[i][j], lj4[i][j]);
      csqoff_v4[i][j] = simd_set_doublev4(cutsq[i][j], offset[i][j], 0, 0);
    }
  doublev4 eng_vdwl_v4 = 0;
  doublev4 eng_coul_v4 = 0;
  doublev4 virial_v4[6];
  for (i = 0; i < 6; i ++)
    virial_v4[i] = 0;
  eng_virial[0] = 0;
  eng_virial[1] = eng_virial[0];
  inum = pm->inum;
  numneigh = pm->numneigh;
  firstneigh = pm->firstneigh;
  int jlist_buf[JLIST_PAGESIZE];
  atom_in_t jlist_atom[JLIST_PAGESIZE];

  doublev4 padding_v4_0, padding_v4_1;
  doublev4 jlist_x_v4[JLIST_PAGESIZE >> 2][3];
  doublev4 jlist_lj_v4[JLIST_PAGESIZE >> 2][4];
  doublev4 jlist_cutsq_v4[JLIST_PAGESIZE >> 2];
  doublev4 jlist_offset_v4[JLIST_PAGESIZE >> 2];
  doublev4 jlist_flj_v4[JLIST_PAGESIZE >> 2];
  //int tsb[8];

  atom_in_t j_cache[JCACHE_LINECNT][JCACHE_LINESIZE];
  int jcache_tag[JCACHE_LINECNT];
  double *xjc, xj[3];
  int tj;
  double xi[ILIST_PAGESIZE][3], fi[ILIST_PAGESIZE][3];
  double ei[ILIST_PAGESIZE], vi[ILIST_PAGESIZE][6], v[6];
  int ti[ILIST_PAGESIZE], *fn[ILIST_PAGESIZE], nn[ILIST_PAGESIZE];
  int ipage_start;
  volatile int cache_reply;
  dma_desc cache_get_desc = 0;
  //memset(&cache_get_desc, 0, sizeof(dma_desc));
  dma_set_mode(&cache_get_desc, PE_MODE);
  dma_set_size(&cache_get_desc, sizeof(atom_in_t) * JCACHE_LINESIZE);
  dma_set_op(&cache_get_desc, DMA_GET);
  dma_set_reply(&cache_get_desc, &cache_reply);
  lwpf_start(COMP);
  for (i = 0; i < JCACHE_LINECNT; i ++)
    jcache_tag[i] = -1;
  for (ipage_start = ILIST_PAGESIZE * _MYID; ipage_start < inum; ipage_start += ILIST_PAGESIZE * 64){
    int ipage_end = ipage_start + ILIST_PAGESIZE;
    if (ipage_end > inum)
      ipage_end = inum;
    int ipage_size = ipage_end - ipage_start;
    lwpf_start(IREAD);
    pe_get(x[ipage_start], xi, ipage_size * sizeof(double) * 3);
    pe_get(type + ipage_start, ti, ipage_size * sizeof(int));
    pe_get(firstneigh + ipage_start, fn, ipage_size * sizeof(int*));
    pe_get(numneigh + ipage_start, nn, ipage_size * sizeof(int));

    pe_syn();
    lwpf_stop(IREAD);
    for (ii = ipage_start; ii < ipage_end; ii ++) {
      i = ii;
      int ioff = i - ipage_start;
      doublev4 fi_v4[3];
      fi_v4[0] = fi_v4[1] = fi_v4[2] = 0;
      doublev4 ei_v4 = 0;
      doublev4 vi_v4[6];
      vi_v4[0] = vi_v4[1] = vi_v4[2] = 0;
      vi_v4[3] = vi_v4[4] = vi_v4[5] = 0;
      //fi[ioff][0] = fi[ioff][1] = fi[ioff][2] = 0;
      //ei[ioff] = 0;
      //vi[ioff][0] = vi[ioff][1] = vi[ioff][2] = 0;
      //vi[ioff][3] = vi[ioff][4] = vi[ioff][5] = 0;
      itype = ti[ioff];
      jlist = fn[ioff];
      jnum = nn[ioff];
      int jpage_start;
      for (jpage_start = 0; jpage_start < jnum; jpage_start += JLIST_PAGESIZE){
        int jpage_size = JLIST_PAGESIZE;
        if (JLIST_PAGESIZE + jpage_start > jnum)
          jpage_size = jnum - jpage_start;
        pe_get(jlist + jpage_start, jlist_buf, jpage_size * sizeof(int));
        pe_syn();
        
        lwpf_start(CACHE);
        for (jj = 0; jj < jpage_size; jj++) {
          j = jlist_buf[jj];
          int sbj = sbmask(j);
          j &= NEIGHMASK;
          int line = j >> JCACHE_SBIT & JCACHE_LMASK;
          int tag = j >> JCACHE_HBIT;

          if (jcache_tag[line] != tag){
            int mem = j & ~JCACHE_MMASK;
            //cache_reply = 0;
            //dma(cache_get_desc, l_pm.atom_in + mem, j_cache[line]);
            //while (cache_reply != 1);
            pe_get(l_pm.atom_in + mem, j_cache[line], sizeof(atom_in_t) * JCACHE_LINESIZE);
            pe_syn();

            jcache_tag[line] = tag;
          }
          jlist_atom[jj] = j_cache[line][j & JCACHE_MMASK];
          jlist_atom[jj].sbj = sbj;
        }
        lwpf_stop(CACHE);
        for (jj = jpage_size; jj < ((jpage_size + 3) & (~3)); jj ++){
          jlist_atom[jj].x[0] = 1e8;
          jlist_atom[jj].sbj = 0;
          jlist_atom[jj].type = 0;
        }

        for (jj = 0; jj < jpage_size; jj += 4){
          int jv = jj >> 2;
          doublev4 xj0, xj1, xj2, xj3;
          simd_load(xj0, jlist_atom + jj + 0);
          simd_load(xj1, jlist_atom + jj + 1);
          simd_load(xj2, jlist_atom + jj + 2);
          simd_load(xj3, jlist_atom + jj + 3);
          transpose4x4(xj0, xj1, xj2, xj3,
                       jlist_x_v4[jv][0], jlist_x_v4[jv][1], jlist_x_v4[jv][2],
                       padding_v4_0);

          int *tsb = &padding_v4_0;
          int t0 = tsb[0]; //jlist_atom[jj + 0].type; //simd_vextw0(padding_v4_0);
          int t1 = tsb[2]; //jlist_atom[jj + 1].type; //simd_vextw2(padding_v4_0);
          int t2 = tsb[4]; //jlist_atom[jj + 2].type; //simd_vextw4(padding_v4_0);
          int t3 = tsb[6]; //jlist_atom[jj + 3].type; //simd_vextw6(padding_v4_0);

          int it = itype;
          transpose4x4(lj_v4[it][t0], lj_v4[it][t1], lj_v4[it][t2], lj_v4[it][t3],
                       jlist_lj_v4[jv][0], jlist_lj_v4[jv][1],
                       jlist_lj_v4[jv][2], jlist_lj_v4[jv][3]);

          int sbj0 = tsb[1]; //jlist_atom[jj + 0].sbj; //simd_vextw1(padding_v4_0);
          int sbj1 = tsb[3]; //jlist_atom[jj + 1].sbj; //simd_vextw3(padding_v4_0);
          int sbj2 = tsb[5]; //jlist_atom[jj + 2].sbj; //simd_vextw5(padding_v4_0);
          int sbj3 = tsb[7]; //jlist_atom[jj + 3].sbj; //simd_vextw7(padding_v4_0);
          jlist_flj_v4[jv] = simd_vshff(slj_v4, slj_v4,
                                        sbj0 | (sbj1 << 2) | (sbj2 << 4) | (sbj3 << 6));
          transpose4x4(csqoff_v4[it][t0], csqoff_v4[it][t1],
                       csqoff_v4[it][t2], csqoff_v4[it][t3],
                       jlist_cutsq_v4[jv], jlist_offset_v4[jv], padding_v4_0, padding_v4_1);

        }
        /* for (jj = 0; jj < jpage_size; jj += 4){ */
        /*   int jv = jj >> 2; */
        /*   doublev4 xj0, xj1, xj2, xj3; */
        /*   simd_load(xj0, jlist_atom + jj + 0); */
        /*   simd_load(xj1, jlist_atom + jj + 1); */
        /*   simd_load(xj2, jlist_atom + jj + 2); */
        /*   simd_load(xj3, jlist_atom + jj + 3); */
        /*   transpose4x4(xj0, xj1, xj2, xj3, */
        /*                jlist_x_v4[jv][0], jlist_x_v4[jv][1], jlist_x_v4[jv][2], */
        /*                padding_v4_0); */
        /*   //simd_store(padding_v4_0, tsb); */
        /*   int *tsb = padding_v4_0; */
        /*   int t0 = 1; //jlist_atom[jj + 0].type; //simd_vextw0(padding_v4_0); */
        /*   int t1 = 1; //jlist_atom[jj + 1].type; //simd_vextw2(padding_v4_0); */
        /*   int t2 = 1; //jlist_atom[jj + 2].type; //simd_vextw4(padding_v4_0); */
        /*   int t3 = 1; //jlist_atom[jj + 3].type; //simd_vextw6(padding_v4_0); */

        /*   int it = itype; */
        /*   transpose4x4(lj_v4[it][t0], lj_v4[it][t1], lj_v4[it][t2], lj_v4[it][t3], */
        /*                jlist_lj_v4[jv][0], jlist_lj_v4[jv][1], */
        /*                jlist_lj_v4[jv][2], jlist_lj_v4[jv][3]); */
        /*   int sbj0 = 0; //jlist_atom[jj + 0].sbj; //simd_vextw1(padding_v4_0); */
        /*   int sbj1 = 0; //jlist_atom[jj + 1].sbj; //simd_vextw3(padding_v4_0); */
        /*   int sbj2 = 0; //jlist_atom[jj + 2].sbj; //simd_vextw5(padding_v4_0); */
        /*   int sbj3 = 0; //jlist_atom[jj + 3].sbj; //simd_vextw7(padding_v4_0); */
        /*   jlist_flj_v4[jv] = simd_vshff(slj_v4, slj_v4, */
        /*                                 sbj0 | (sbj1 << 2) | (sbj2 << 4) | (sbj3 << 6)); */
        /*   transpose4x4(csqoff_v4[it][t0], csqoff_v4[it][t1], */
        /*                csqoff_v4[it][t2], csqoff_v4[it][t3], */
        /*                jlist_cutsq_v4[jv], jlist_offset_v4[jv], padding_v4_0, padding_v4_1); */
        /* } */
        lwpf_start(JLOOP);
        int jpage_size_v4 = (jpage_size + 3) >> 2;
        doublev4 xi_v4[3];
        xi_v4[0] = xi[ioff][0];
        xi_v4[1] = xi[ioff][1];
        xi_v4[2] = xi[ioff][2];
        doublev4 v4_1 = 1.0;
        doublev4 v4_0 = 0;
        doublev4 v4_half = 0.5;
        for (jj = 0; jj < jpage_size_v4; jj ++){
          doublev4 delx_v4 = xi_v4[0] - jlist_x_v4[jj][0];
          doublev4 dely_v4 = xi_v4[1] - jlist_x_v4[jj][1];
          doublev4 delz_v4 = xi_v4[2] - jlist_x_v4[jj][2];
          doublev4 rsq_v4 = delx_v4 * delx_v4 + dely_v4 * dely_v4 + delz_v4 * delz_v4;
          
          doublev4 msk_v4 = -simd_vfcmplt(rsq_v4, jlist_cutsq_v4[jj]);
          doublev4 r2inv_v4 = v4_1 / rsq_v4;
          r2inv_v4 = simd_vsellt(msk_v4, r2inv_v4, v4_0);
          doublev4 r6inv_v4 = r2inv_v4 * r2inv_v4 * r2inv_v4;
          doublev4 forcelj_v4 = r6inv_v4 * (jlist_lj_v4[jj][0] * r6inv_v4 - jlist_lj_v4[jj][1]);
          doublev4 fpair_v4 = jlist_flj_v4[jj] * forcelj_v4 * r2inv_v4;
          fpair_v4 = simd_vsellt(msk_v4, fpair_v4, v4_0);
          delx_v4 = simd_vsellt(msk_v4, delx_v4, v4_0);
          
          fi_v4[0] += delx_v4 * fpair_v4;
          fi_v4[1] += dely_v4 * fpair_v4;
          fi_v4[2] += delz_v4 * fpair_v4;
          doublev4 v_fac = simd_vsellt(msk_v4, v4_half, v4_0);
          if (l_pm.evflag) {
            if (l_pm.eflag_either){
              doublev4 evdwl_v4 = r6inv_v4 * (jlist_lj_v4[jj][2] * r6inv_v4 - jlist_lj_v4[jj][3]) - jlist_offset_v4[jj];
              evdwl_v4 *= jlist_flj_v4[jj];
              if (_MYID == 0){
                /* simd_print_doublev4(evdwl_v4); */
                /* simd_print_doublev4(eng_vdwl_v4); */
                /* simd_print_doublev4(jlist_flj_v4[jj]); */
                /* simd_print_doublev4(jlist_lj_v4[jj][2]); */
                /* simd_print_doublev4(jlist_lj_v4[jj][3]); */
                /* simd_print_doublev4(jlist_offset_v4[jj]); */
                /* simd_print_doublev4(r6inv_v4); */
                /* simd_print_doublev4(r2inv_v4); */
                /* simd_print_doublev4(rsq_v4); */
                /* puts(""); */
              }
              if (l_pm.eflag_global)
                eng_vdwl_v4 += v_fac * evdwl_v4;
              if (l_pm.eflag_atom) ei_v4 += v_fac * evdwl_v4;
            }
            if (l_pm.vflag_either) {
              doublev4 v_v4[6];
              v_v4[0] = v_fac * delx_v4 * delx_v4 * fpair_v4;
              v_v4[1] = v_fac * dely_v4 * dely_v4 * fpair_v4;
              v_v4[2] = v_fac * delz_v4 * delz_v4 * fpair_v4;
              v_v4[3] = v_fac * delx_v4 * dely_v4 * fpair_v4;
              v_v4[4] = v_fac * delx_v4 * delz_v4 * fpair_v4;
              v_v4[5] = v_fac * dely_v4 * delz_v4 * fpair_v4;
              if (l_pm.vflag_global) {
                virial_v4[0] += v_v4[0];
                virial_v4[1] += v_v4[1];
                virial_v4[2] += v_v4[2];
                virial_v4[3] += v_v4[3];
                virial_v4[4] += v_v4[4];
                virial_v4[5] += v_v4[5];
              }
              if (l_pm.vflag_atom) {
                vi_v4[0] += v_v4[0];
                vi_v4[1] += v_v4[1];
                vi_v4[2] += v_v4[2];
                vi_v4[3] += v_v4[3];
                vi_v4[4] += v_v4[4];
                vi_v4[5] += v_v4[5];
              }
            }
          }
        }
        lwpf_stop(JLOOP);
      }

      simd_vsumd(fi_v4[0]);
      simd_vsumd(fi_v4[1]);
      simd_vsumd(fi_v4[2]);

      simd_vsumd(vi_v4[0]);
      simd_vsumd(vi_v4[1]);
      simd_vsumd(vi_v4[2]);
      simd_vsumd(vi_v4[3]);
      simd_vsumd(vi_v4[4]);
      simd_vsumd(vi_v4[5]);

      simd_vsumd(ei_v4);

      fi[ioff][0] = fi_v4[0];
      fi[ioff][1] = fi_v4[1];
      fi[ioff][2] = fi_v4[2];

      vi[ioff][0] = vi_v4[0];
      vi[ioff][1] = vi_v4[1];
      vi[ioff][2] = vi_v4[2];
      vi[ioff][3] = vi_v4[3];
      vi[ioff][4] = vi_v4[4];
      vi[ioff][5] = vi_v4[5];

      ei[ioff] = ei_v4;

    }
    lwpf_start(IWRITE);
    pe_put(f[ipage_start], fi[0], sizeof(double) * 3 * ipage_size);
    if (l_pm.eflag_either && l_pm.eflag_atom)
      pe_put(l_pm.eatom + ipage_start, ei, sizeof(double) * ipage_size);
    if (l_pm.vflag_either && l_pm.vflag_atom)
      pe_put(l_pm.vatom[ipage_start], vi[0], sizeof(double) * 6 * ipage_size);
    pe_syn();
    lwpf_stop(IWRITE);
  }
  lwpf_stop(COMP);
  simd_vsumd(eng_vdwl_v4);
  simd_vsumd(eng_coul_v4);
  simd_vsumd(virial_v4[0]);
  simd_vsumd(virial_v4[1]);
  simd_vsumd(virial_v4[2]);
  simd_vsumd(virial_v4[3]);
  simd_vsumd(virial_v4[4]);
  simd_vsumd(virial_v4[5]);

  *eng_vdwl = eng_vdwl_v4;
  *eng_coul = eng_coul_v4;
  virial[0] = virial_v4[0];
  virial[1] = virial_v4[1];
  virial[2] = virial_v4[2];
  virial[3] = virial_v4[3];
  virial[4] = virial_v4[4];
  virial[5] = virial_v4[5];
  /* for (i = 0; i < 64; i ++){ */
  /*   if (_MYID == i){ */
  /*     simd_print_doublev4(eng_vdwl_v4); */
  /*     printf("%d %f\n", i, *eng_vdwl); */
  /*   } */
  /*   athread_syn(ARRAY_SCOPE, 0xffff); */
  /* } */
  reg_reduce_inplace_doublev4(eng_virial, 2);
  lwpf_start(GST);
  if (_MYID == 0){
    pe_put(&(pm->eng_vdwl), eng_vdwl, sizeof(double) * 8);
    pe_syn();
  }
  lwpf_stop(GST);
  lwpf_stop(ALL);
}

#endif
