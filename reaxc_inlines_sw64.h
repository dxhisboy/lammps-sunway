static inline void rvec_Copy( rvec dest, rvec src )
{
  dest[0] = src[0], dest[1] = src[1], dest[2] = src[2];
}


static inline void rvec_Scale( rvec ret, double c, rvec v )
{
  ret[0] = c * v[0], ret[1] = c * v[1], ret[2] = c * v[2];
}


static inline void rvec_Add( rvec ret, rvec v )
{
  ret[0] += v[0], ret[1] += v[1], ret[2] += v[2];
}


static inline void rvec_ScaledAdd( rvec ret, double c, rvec v )
{
  ret[0] += c * v[0], ret[1] += c * v[1], ret[2] += c * v[2];
}


static inline void rvec_ScaledSum( rvec ret, double c1, rvec v1 ,double c2, rvec v2 )
{
  ret[0] = c1 * v1[0] + c2 * v2[0];
  ret[1] = c1 * v1[1] + c2 * v2[1];
  ret[2] = c1 * v1[2] + c2 * v2[2];
}


double rvec_Dot( rvec v1, rvec v2 )
{
  return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
}


static inline void rvec_iMultiply( rvec r, ivec v1, rvec v2 )
{
  r[0] = v1[0] * v2[0];
  r[1] = v1[1] * v2[1];
  r[2] = v1[2] * v2[2];
}


static inline void rvec_Cross( rvec ret, rvec v1, rvec v2 )
{
  ret[0] = v1[1] * v2[2] - v1[2] * v2[1];
  ret[1] = v1[2] * v2[0] - v1[0] * v2[2];
  ret[2] = v1[0] * v2[1] - v1[1] * v2[0];
}


static inline double rvec_Norm_Sqr( rvec v )
{
  return SQR(v[0]) + SQR(v[1]) + SQR(v[2]);
}


static inline double rvec_Norm( rvec v )
{
  return sqrt( SQR(v[0]) + SQR(v[1]) + SQR(v[2]) );
}


static inline void rvec_MakeZero( rvec v )
{
  v[0] = v[1] = v[2] = 0.000000000000000e+00;
}


static inline void rtensor_MatVec( rvec ret, rtensor m, rvec v )
{
  int i;
  rvec temp;

  if( ret == v )
    {
      for( i = 0; i < 3; ++i )
        temp[i] = m[i][0] * v[0] + m[i][1] * v[1] + m[i][2] * v[2];

      for( i = 0; i < 3; ++i )
        ret[i] = temp[i];
    }
  else
    {
      for( i = 0; i < 3; ++i )
        ret[i] = m[i][0] * v[0] + m[i][1] * v[1] + m[i][2] * v[2];
    }
}


static inline void rtensor_MakeZero( rtensor t )
{
  t[0][0] = t[0][1] = t[0][2] = 0;
  t[1][0] = t[1][1] = t[1][2] = 0;
  t[2][0] = t[2][1] = t[2][2] = 0;
}


static inline void ivec_MakeZero( ivec v )
{
  v[0] = v[1] = v[2] = 0;
}


static inline void ivec_Copy( ivec dest, ivec src )
{
  dest[0] = src[0], dest[1] = src[1], dest[2] = src[2];
}


static inline void ivec_Scale( ivec dest, double C, ivec src )
{
  dest[0] = (int)(C * src[0]);
  dest[1] = (int)(C * src[1]);
  dest[2] = (int)(C * src[2]);
}


static inline void ivec_Sum( ivec dest, ivec v1, ivec v2 )
{
  dest[0] = v1[0] + v2[0];
  dest[1] = v1[1] + v2[1];
  dest[2] = v1[2] + v2[2];
}

static inline int Start_Index( int i, reax_list *l )
{
  return l->index[i];
}

static inline int End_Index( int i, reax_list *l )
{
  return l->end_index[i];
}

static inline void ev_tally_full(int i, double evdwl, double ecoul, double fpair,
                                 double delx, double dely, double delz,
                                 int eflag_global, int vflag_global,
                                 int eflag_atom, int vflag_atom,
                                 double *eng_vdwl, double *eng_coul,
                                 double *virial,
                                 double *eatom, double (*vatom)[6]){
  double v[6];

  if (eflag_global || eflag_atom) {
    if (eflag_global) {
      eng_vdwl[0] += 0.5*evdwl;
      eng_coul[0] += 0.5*ecoul;
    }
    if (eflag_atom) eatom[i] += 0.5 * (evdwl + ecoul);
  }

  if (vflag_global || vflag_atom) {
    v[0] = 0.5*delx*delx*fpair;
    v[1] = 0.5*dely*dely*fpair;
    v[2] = 0.5*delz*delz*fpair;
    v[3] = 0.5*delx*dely*fpair;
    v[4] = 0.5*delx*delz*fpair;
    v[5] = 0.5*dely*delz*fpair;

    if (vflag_global) {
      virial[0] += v[0];
      virial[1] += v[1];
      virial[2] += v[2];
      virial[3] += v[3];
      virial[4] += v[4];
      virial[5] += v[5];
    }

    if (vflag_atom) {
      vatom[i][0] += v[0];
      vatom[i][1] += v[1];
      vatom[i][2] += v[2];
      vatom[i][3] += v[3];
      vatom[i][4] += v[4];
      vatom[i][5] += v[5];
    }
  }

}

static inline void ev_tally_full_sys(int i, double evdwl, double ecoul, double fpair,
                                     double *del, reax_system_c *sys){
  double v[6];

  if (sys->eflag_global || sys->eflag_atom) {
    if (sys->eflag_global) {
      sys->eng_vdwl += 0.5*evdwl;
      sys->eng_coul += 0.5*ecoul;
    }
    if (sys->eflag_atom) sys->eatom[i] += 0.5 * (evdwl + ecoul);
  }

  if (sys->vflag_global || sys->vflag_atom) {
    v[0] = 0.5*del[0]*del[0]*fpair;
    v[1] = 0.5*del[1]*del[1]*fpair;
    v[2] = 0.5*del[2]*del[2]*fpair;
    v[3] = 0.5*del[0]*del[1]*fpair;
    v[4] = 0.5*del[0]*del[2]*fpair;
    v[5] = 0.5*del[1]*del[2]*fpair;

    if (sys->vflag_global) {
      sys->virial[0] += v[0];
      sys->virial[1] += v[1];
      sys->virial[2] += v[2];
      sys->virial[3] += v[3];
      sys->virial[4] += v[4];
      sys->virial[5] += v[5];
    }

    if (sys->vflag_atom) {
      sys->vatom[i][0] += v[0];
      sys->vatom[i][1] += v[1];
      sys->vatom[i][2] += v[2];
      sys->vatom[i][3] += v[3];
      sys->vatom[i][4] += v[4];
      sys->vatom[i][5] += v[5];
    }
  }

}
