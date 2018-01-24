#ifndef PAIR_LJ_CUT_SW64_H_
#define PAIR_LJ_CUT_SW64_H_
#ifdef __cplusplus
extern "C"{
#endif
  typedef struct atom_in_t{
    double x[3];
    int type, sbj;
  } atom_in_t;
  typedef struct pair_lj_cut_compute_param_t{
    int *ilist, *numneigh, **firstneigh;
    //vars to be put back
    double (*x)[3], (*f)[3], (*vatom)[6], *eatom;
    atom_in_t *atom_in;
    double *cutsq, *lj1, *lj2, *lj3, *lj4, *offset;
    int *type;
    int nlocal, nghost, ntotal, inum, ntypes, rank;
    double special_lj[4];
    //vars to be copy back
    double eng_vdwl, eng_coul, virial[6];
    int eflag, vflag, evflag;
    int eflag_global, vflag_global;
    int eflag_atom, vflag_atom;
    int eflag_either, vflag_either;
  } pair_lj_cut_compute_param_t;
#ifdef __cplusplus
}
#endif
#endif
#ifdef COMPUTE_FUNC
{
  //lwpf_start(ALL);
  pe_init();
  volatile int cache_reply;            
  doublev4 v4_1 = 1.0;
  doublev4 v4_0 = 0;
  doublev4 v4_half = 0.5;

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

  pair_lj_cut_compute_param_t l_pm;
  pe_get(pm, &l_pm, sizeof(pair_lj_cut_compute_param_t));
  //double INF = 1.0 / 0.0;
  doublev4 INF_v4 = simd_set_doublev4(1e8, 1e8, 1e8, 0);
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
  doublev4 jlist_x_v4[3];
  doublev4 jlist_lj_v4[4];
  doublev4 jlist_cutsq_v4;
  doublev4 jlist_offset_v4;
  doublev4 jlist_flj_v4;

  atom_in_t j_cache[JCACHE_LINECNT][JCACHE_LINESIZE];
  int jcache_tag[JCACHE_LINECNT];
  double *xjc, xj[3];
  int tj;
  double xi[ILIST_PAGESIZE][3], fi[ILIST_PAGESIZE][3];
  double ei[ILIST_PAGESIZE], vi[ILIST_PAGESIZE][6], v[6];
  int ti[ILIST_PAGESIZE], *fn[ILIST_PAGESIZE], nn[ILIST_PAGESIZE];
  int ipage_start;

  volatile dma_desc cache_get_desc = 0;

  dma_set_mode(&cache_get_desc, PE_MODE);
  dma_set_size(&cache_get_desc, sizeof(atom_in_t) * JCACHE_LINESIZE);
  dma_set_op(&cache_get_desc, DMA_GET);
  dma_set_reply(&cache_get_desc, &cache_reply);

  //lwpf_start(COMP);
  for (i = 0; i < JCACHE_LINECNT; i ++)
    jcache_tag[i] = -1;
  for (ipage_start = ILIST_PAGESIZE * _MYID; ipage_start < inum; ipage_start += ILIST_PAGESIZE * 64){
    int ipage_end = ipage_start + ILIST_PAGESIZE;
    if (ipage_end > inum)
      ipage_end = inum;
    int ipage_size = ipage_end - ipage_start;
    //lwpf_start(IREAD);
    pe_get(x[ipage_start], xi, ipage_size * sizeof(double) * 3);
    pe_get(type + ipage_start, ti, ipage_size * sizeof(int));
    pe_get(firstneigh + ipage_start, fn, ipage_size * sizeof(int*));
    pe_get(numneigh + ipage_start, nn, ipage_size * sizeof(int));

    pe_syn();
    //lwpf_stop(IREAD);
    for (ii = ipage_start; ii < ipage_end; ii ++) {
      i = ii;
      int ioff = i - ipage_start;
      doublev4 fi_v4[3];
      fi_v4[0] = fi_v4[1] = fi_v4[2] = 0;
      doublev4 ei_v4 = 0;
      doublev4 vi_v4[6];
      vi_v4[0] = vi_v4[1] = vi_v4[2] = 0;
      vi_v4[3] = vi_v4[4] = vi_v4[5] = 0;

.      doublev4 xi_v4[3];
      xi_v4[0] = xi[ioff][0];
      xi_v4[1] = xi[ioff][1];
      xi_v4[2] = xi[ioff][2];

      itype = ti[ioff];
      jlist = fn[ioff];
      jnum = nn[ioff];
      int jpage_start, jpage_size;
      for (jpage_start = 0; jpage_start < jnum; jpage_start += JLIST_PAGESIZE){
        jpage_size = JLIST_PAGESIZE;
        if (JLIST_PAGESIZE + jpage_start > jnum)
          jpage_size = jnum - jpage_start;
        pe_get(jlist + jpage_start, jlist_buf, jpage_size * sizeof(int));
        pe_syn();
        
        int jpage_size_to = jpage_size & ~3;

        for (jj = 0; jj < jpage_size_to; jj+=4) {
          int sbj[4];
          int jjj;
          int flj_msk = 0;
          doublev4 xjv[4];
          //lwpf_start(CACHE);
          for (jjj = 0; jjj < 4; jjj ++){
            j = jlist_buf[jj + jjj];
            sbj[jjj] = sbmask(j);
            j &= NEIGHMASK;
            int line = j >> JCACHE_SBIT & JCACHE_LMASK;
            int tag = j >> JCACHE_HBIT;
            flj_msk |= sbj[jjj] << (jjj * 2);
            if (jcache_tag[line] != tag){
              int mem = j & ~JCACHE_MMASK;
              cache_reply = 0;
              dma(cache_get_desc, l_pm.atom_in + mem, j_cache[line]);
              while (cache_reply != 1);
              jcache_tag[line] = tag;
            }
            simd_load(xjv[jjj], j_cache[line] + (j & JCACHE_MMASK));
          }
          
          //lwpf_stop(CACHE);
          //lwpf_start(VECPACK);
          {
            transpose4x4(xjv[0], xjv[1], xjv[2], xjv[3],
                         jlist_x_v4[0], jlist_x_v4[1], jlist_x_v4[2],
                         padding_v4_0);

            /* int t0, t1, t2, t3; */
            /* asm("vextw %4, 0, %0\n" */
            /*     "vextw %4, 2, %1\n" */
            /*     "vextw %4, 4, %2\n" */
            /*     "vextw %4, 6, %3\n" */
            /*     : "=r"(t0), "=r"(t1), "=r"(t2), "=r"(t3) : "r"(padding_v4_0)); */
            int *tsb = &padding_v4_0;
            int t0 = tsb[0];
            int t1 = tsb[2];
            int t2 = tsb[4];
            int t3 = tsb[6];

            int it = itype;
            transpose4x4(lj_v4[it][t0], lj_v4[it][t1], lj_v4[it][t2], lj_v4[it][t3],
                         jlist_lj_v4[0], jlist_lj_v4[1],
                         jlist_lj_v4[2], jlist_lj_v4[3]);
            transpose4x4_2x4(csqoff_v4[it][t0], csqoff_v4[it][t1],
                             csqoff_v4[it][t2], csqoff_v4[it][t3],
                             jlist_cutsq_v4, jlist_offset_v4);
            jlist_flj_v4 = simd_vshff(slj_v4, slj_v4, flj_msk);
          }
          //lwpf_stop(VECPACK);
          //lwpf_start(VECCOMP);
          doublev4 delx_v4 = xi_v4[0] - jlist_x_v4[0];
          doublev4 dely_v4 = xi_v4[1] - jlist_x_v4[1];
          doublev4 delz_v4 = xi_v4[2] - jlist_x_v4[2];
          doublev4 rsq_v4 = delx_v4 * delx_v4 + dely_v4 * dely_v4 + delz_v4 * delz_v4;
          
          doublev4 msk_v4 = simd_vfcmple(jlist_cutsq_v4, rsq_v4);
          doublev4 r2inv_v4 = v4_1 / rsq_v4;
          doublev4 r6inv_v4 = r2inv_v4 * r2inv_v4 * r2inv_v4;
          doublev4 forcelj_v4 = r6inv_v4 * (jlist_lj_v4[0] * r6inv_v4 - jlist_lj_v4[1]);
          doublev4 fpair_v4 = jlist_flj_v4 * forcelj_v4 * r2inv_v4;
          fpair_v4 = simd_vselle(msk_v4, fpair_v4, v4_0);
          
          fi_v4[0] += delx_v4 * fpair_v4;
          fi_v4[1] += dely_v4 * fpair_v4;
          fi_v4[2] += delz_v4 * fpair_v4;
          if (l_pm.evflag) {
            doublev4 v_fac = simd_vselle(msk_v4, v4_half, v4_0);
            if (l_pm.eflag_either){
              doublev4 evdwl_v4 = r6inv_v4 * (jlist_lj_v4[2] * r6inv_v4 - jlist_lj_v4[3]) - jlist_offset_v4;
              evdwl_v4 *= jlist_flj_v4;
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
          //lwpf_stop(VECCOMP);
        }
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

      for (jj = jpage_size & ~3; jj < jpage_size; jj++) {
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
          jcache_tag[line] = tag;
        }
        
        delx = xi[ioff][0] - j_cache[line][j & JCACHE_MMASK].x[0];
        dely = xi[ioff][1] - j_cache[line][j & JCACHE_MMASK].x[1];
        delz = xi[ioff][2] - j_cache[line][j & JCACHE_MMASK].x[2];
        jtype = j_cache[line][j & JCACHE_MMASK].type;
        rsq = delx * delx + dely * dely + delz * delz;
        
        if (rsq < cutsq[itype][jtype]){
          r2inv = 1 / rsq;
          r6inv = r2inv * r2inv * r2inv;
          forcelj = r6inv * (lj1[itype][jtype] * r6inv - lj2[itype][jtype]);
          fpair = factor_lj * forcelj * r2inv;
          
          fi[ioff][0] += delx * fpair;
          fi[ioff][1] += dely * fpair;
          fi[ioff][2] += delz * fpair;
          if (l_pm.evflag) {
            if (l_pm.eflag_either){
              double evdwl = r6inv * (lj3[itype][jtype] * r6inv - lj4[itype][jtype]) - offset[itype][jtype];
              evdwl *= factor_lj;
              if (l_pm.eflag_global)
                *eng_vdwl += 0.5 * evdwl;
              if (l_pm.eflag_atom) ei[ioff] += 0.5 * evdwl;
            }
            if (l_pm.vflag_either) {
              double v[6];
              v[0] = 0.5 * delx * delx * fpair;
              v[1] = 0.5 * dely * dely * fpair;
              v[2] = 0.5 * delz * delz * fpair;
              v[3] = 0.5 * delx * dely * fpair;
              v[4] = 0.5 * delx * delz * fpair;
              v[5] = 0.5 * dely * delz * fpair;
              if (l_pm.vflag_global) {
                virial[0] += v[0];
                virial[1] += v[1];
                virial[2] += v[2];
                virial[3] += v[3];
                virial[4] += v[4];
                virial[5] += v[5];
              }
              if (l_pm.vflag_atom) {
                vi[ioff][0] += v[0];
                vi[ioff][1] += v[1];
                vi[ioff][2] += v[2];
                vi[ioff][3] += v[3];
                vi[ioff][4] += v[4];
                vi[ioff][5] += v[5];
              }
            }
          }
        }
      }
    }
    //lwpf_start(IWRITE);
    pe_put(f[ipage_start], fi[0], sizeof(double) * 3 * ipage_size);
    if (l_pm.eflag_either && l_pm.eflag_atom)
      pe_put(l_pm.eatom + ipage_start, ei, sizeof(double) * ipage_size);
    if (l_pm.vflag_either && l_pm.vflag_atom)
      pe_put(l_pm.vatom[ipage_start], vi[0], sizeof(double) * 6 * ipage_size);
    pe_syn();
    //lwpf_stop(IWRITE);
  }
  //lwpf_stop(COMP);
  simd_vsumd(eng_vdwl_v4);
  simd_vsumd(eng_coul_v4);
  simd_vsumd(virial_v4[0]);
  simd_vsumd(virial_v4[1]);
  simd_vsumd(virial_v4[2]);
  simd_vsumd(virial_v4[3]);
  simd_vsumd(virial_v4[4]);
  simd_vsumd(virial_v4[5]);

  *eng_vdwl += eng_vdwl_v4;
  *eng_coul += eng_coul_v4;
  virial[0] += virial_v4[0];
  virial[1] += virial_v4[1];
  virial[2] += virial_v4[2];
  virial[3] += virial_v4[3];
  virial[4] += virial_v4[4];
  virial[5] += virial_v4[5];

  reg_reduce_inplace_doublev4(eng_virial, 2);
  if (_MYID == 0){
    pe_put(&(pm->eng_vdwl), eng_vdwl, sizeof(double) * 8);
    pe_syn();
  }

  //lwpf_stop(ALL);
}
#endif
