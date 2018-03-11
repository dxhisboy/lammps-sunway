extern inline double p_lnd(double in);
extern inline double p_expd(double in);
extern inline double p_powd(double base, double power);
extern inline double p_sind(double in);
extern inline double p_cosd(double in);
extern inline double p_sinnpi_pid(double in);
extern inline double p_cosnpi_pid(double in);

extern inline doublev4 simd_vfloord(doublev4 in);
extern inline doublev4 simd_vceild(doublev4 in);
extern inline doublev4 simd_vroundd(doublev4 in);
extern inline doublev4 simd_vltod(doublev4 in);
extern inline doublev4 simd_vlnd(doublev4 in);
extern inline doublev4 simd_vln1_2d(doublev4 in);
extern inline doublev4 simd_vexpd(doublev4 in);
extern inline doublev4 simd_vpowd(doublev4 base, doublev4 power);
extern inline doublev4 simd_vsind(doublev4 in);
extern inline doublev4 simd_vcosd(doublev4 in);
extern inline doublev4 simd_vsinnpi_pid(doublev4 in);
extern inline doublev4 simd_vcosnpi_pid(doublev4 in);
extern inline doublev4 simd_vinv_sqrtd(doublev4 in);

extern __thread_local long exp_ln_coef[24];
extern __thread_local long sin_cos_coef[24];
extern __thread_local long inv_sqrt_coef[4];

