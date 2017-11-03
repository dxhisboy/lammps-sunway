void FNF(compute_param_t *pm){
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

  atom_in_t j_cache[JCACHE_LINECNT][JCACHE_LINESIZE];
  int jcache_tag[JCACHE_LINECNT];
  double *xjc, xj[3];
  int tj;
  double xi[ILIST_PAGESIZE][3], fi[ILIST_PAGESIZE][3];
  double ei[ILIST_PAGESIZE], vi[ILIST_PAGESIZE][6], v[6];
  int ti[ILIST_PAGESIZE], *fn[ILIST_PAGESIZE], nn[ILIST_PAGESIZE];
  int ipage_start;
  volatile int cache_reply;
  dma_desc cache_get_desc;
  memset(&cache_get_desc, 0, sizeof(dma_desc));
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
        //lwpf_start(JLOOP);
        for (jj = 0; jj < jpage_size; jj++) {
          j = jlist_buf[jj];
          factor_lj = special_lj[sbmask(j)];
          j &= NEIGHMASK;
          int line = j >> JCACHE_SBIT & JCACHE_LMASK;
          int tag = j >> JCACHE_HBIT;

          //lwpf_start(CACHE);
          if (jcache_tag[line] != tag){
            int mem = j & ~JCACHE_MMASK;
            cache_reply = 0;
            dma(cache_get_desc, l_pm.atom_in + mem, j_cache[line]);
            while (cache_reply != 1);
            jcache_tag[line] = tag;
          }
          //lwpf_stop(CACHE);
          atom_in_t *jc = j_cache[line] + (j & JCACHE_MMASK);
          xj[0] = jc->x[0];
          xj[1] = jc->x[1];
          xj[2] = jc->x[2];
          jtype = jc->type;

          delx = xi[ioff][0] - xj[0];
          dely = xi[ioff][1] - xj[1];
          delz = xi[ioff][2] - xj[2];
          rsq = delx*delx + dely*dely + delz*delz;

          if (rsq < cutsq[itype][jtype]) {
            r2inv = 1.0/rsq;
            r6inv = r2inv*r2inv*r2inv;
            forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
            fpair = factor_lj*forcelj*r2inv;

            fi[ioff][0] += delx*fpair;
            fi[ioff][1] += dely*fpair;
            fi[ioff][2] += delz*fpair;

#if defined(EATOM) || defined(EGLOBAL)
            evdwl = r6inv*(lj3[itype][jtype]*r6inv-lj4[itype][jtype]) -
              offset[itype][jtype];
            evdwl *= factor_lj;
#endif
            //            if (l_pm.evflag) {
#if defined(EATOM) || defined(EGLOBAL) || defined(VATOM) || defined(VGLOBAL)
            //if (l_pm.eflag_either) {
#if defined(EATOM) || defined(EGLOBAL)
            //if (l_pm.eflag_global) {
#if defined(EGLOBAL)
            eng_vdwl[0] += 0.5*evdwl;
            eng_coul[0] += 0.5*0;
#endif
            //  }
            //if (l_pm.eflag_atom)
#if defined(EATOM)
            ei[i] += 0.5 * (evdwl + 0);
#endif
#endif
            //}
            //if (l_pm.vflag_either) {
#if defined(VATOM) || defined(VGLOBAL)
            v[0] = 0.5*delx*delx*fpair;
            v[1] = 0.5*dely*dely*fpair;
            v[2] = 0.5*delz*delz*fpair;
            v[3] = 0.5*delx*dely*fpair;
            v[4] = 0.5*delx*delz*fpair;
            v[5] = 0.5*dely*delz*fpair;

            //if (l_pm.vflag_global) {
#if defined(VGLOBAL)
            virial[0] += v[0];
            virial[1] += v[1];
            virial[2] += v[2];
            virial[3] += v[3];
            virial[4] += v[4];
            virial[5] += v[5];
            //}
#endif
            //if (l_pm.vflag_atom) {
#if defined(VATOM)
            vi[i][0] += v[0];
            vi[i][1] += v[1];
            vi[i][2] += v[2];
            vi[i][3] += v[3];
            vi[i][4] += v[4];
            vi[i][5] += v[5];
            //}
#endif
#endif
            //}
#endif
          }
        }
        //lwpf_stop(JLOOP);
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
  reg_reduce_inplace_doublev4(eng_virial, 2);
  //lwpf_start(GST);
  if (_MYID == 0){
    pe_put(&(pm->eng_vdwl), eng_vdwl, sizeof(double) * 8);
    pe_syn();
  }
}
#undef FN
#undef FNEG
#undef FNEA
#undef FNVG
#undef FNF
