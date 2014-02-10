#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#define NUM_PCLS (8192 * 8)  // 2048 * 30 * 30
#define ALIGNMENT 64

#define ALIGNED(X) __assume_aligned(X, ALIGNMENT)

inline double time_sec()
{
  static struct timeval tv;

  gettimeofday(&tv, NULL);

  return (tv.tv_sec + tv.tv_usec * (double) 1e-6);
}

int main(void)
{
  /* Static variables */                                  
  int num_pcls = NUM_PCLS;
  int pidx, c, i, j;
  double time1, time2;
  double *ptr;
  double qdto2mc;
  double abs_pos[3];
  double rel_pos[3];
  double cm1_pos[3];
  double w0[3], w1[3], weight[4];
  double Om[3];
  double t[3], orig[3];
  double omsq, denom;
  double udotOm;
  double avg[3];
  double *weights[8];
  double *Bxl[8], *Byl[8], *Bzl[8];
  double *Exl[8], *Eyl[8], *Ezl[8];
  double *field_components[8][6];
  
  // Member variables (random values)
  double xstart = 0, ystart = 0, zstart = 0;
  double inv_dx = .25, inv_dy = .25, inv_dz = .25;
  double dto2 = .4;
  double cx = 1, cy = 2, cz = 3;

  
  /* Alloc */
  double *x = ptr = _mm_malloc(sizeof(double) * NUM_PCLS, ALIGNMENT);
  if (!ptr) { printf("%d: Malloc failed.\n", __LINE__); exit(1); }
  double *y = ptr = _mm_malloc(sizeof(double) * NUM_PCLS, ALIGNMENT);
  if (!ptr) { printf("%d: Malloc failed.\n", __LINE__); exit(1); }
  double *z = ptr = _mm_malloc(sizeof(double) * NUM_PCLS, ALIGNMENT);
  if (!ptr) { printf("%d: Malloc failed.\n", __LINE__); exit(1); }
  double *u = ptr = _mm_malloc(sizeof(double) * NUM_PCLS, ALIGNMENT);
  if (!ptr) { printf("%d: Malloc failed.\n", __LINE__); exit(1); }
  double *v = ptr = _mm_malloc(sizeof(double) * NUM_PCLS, ALIGNMENT);
  if (!ptr) { printf("%d: Malloc failed.\n", __LINE__); exit(1); }
  double *w = ptr = _mm_malloc(sizeof(double) * NUM_PCLS, ALIGNMENT);
  if (!ptr) { printf("%d: Malloc failed.\n", __LINE__); exit(1); }
  double *_xavg = ptr = _mm_malloc(sizeof(double) * NUM_PCLS, ALIGNMENT);
  if (!ptr) { printf("%d: Malloc failed.\n", __LINE__); exit(1); }
  double *_yavg = ptr = _mm_malloc(sizeof(double) * NUM_PCLS, ALIGNMENT);
  if (!ptr) { printf("%d: Malloc failed.\n", __LINE__); exit(1); }
  double *_zavg = ptr = _mm_malloc(sizeof(double) * NUM_PCLS, ALIGNMENT);
  if (!ptr) { printf("%d: Malloc failed.\n", __LINE__); exit(1); }
  ALIGNED(x);
  ALIGNED(y);
  ALIGNED(z);
  ALIGNED(u);
  ALIGNED(v);
  ALIGNED(w);
  ALIGNED(_xavg);
  ALIGNED(_yavg);
  ALIGNED(_zavg);
  
  for (i = 0; i < 8; i++) {
    weights[i] = ptr = _mm_malloc(sizeof(double) * NUM_PCLS, ALIGNMENT);
    if (!ptr) { printf("%d: Malloc failed.\n", __LINE__); exit(1); }
    Bxl[i] = ptr = _mm_malloc(sizeof(double) * NUM_PCLS, ALIGNMENT);
    if (!ptr) { printf("%d: Malloc failed.\n", __LINE__); exit(1); }
    Byl[i] = ptr = _mm_malloc(sizeof(double) * NUM_PCLS, ALIGNMENT);
    if (!ptr) { printf("%d: Malloc failed.\n", __LINE__); exit(1); }
    Bzl[i] = ptr = _mm_malloc(sizeof(double) * NUM_PCLS, ALIGNMENT);
    if (!ptr) { printf("%d: Malloc failed.\n", __LINE__); exit(1); }
    Exl[i] = ptr = _mm_malloc(sizeof(double) * NUM_PCLS, ALIGNMENT);
    if (!ptr) { printf("%d: Malloc failed.\n", __LINE__); exit(1); }
    Eyl[i] = ptr = _mm_malloc(sizeof(double) * NUM_PCLS, ALIGNMENT);
    if (!ptr) { printf("%d: Malloc failed.\n", __LINE__); exit(1); }
    Ezl[i] = ptr = _mm_malloc(sizeof(double) * NUM_PCLS, ALIGNMENT);
    if (!ptr) { printf("%d: Malloc failed.\n", __LINE__); exit(1); }
  }
  ALIGNED(weights[0]); ALIGNED(weights[1]); ALIGNED(weights[2]); ALIGNED(weights[3]); ALIGNED(weights[4]); ALIGNED(weights[5]); ALIGNED(weights[6]); ALIGNED(weights[7]);
  ALIGNED(Bxl[0]); ALIGNED(Bxl[1]); ALIGNED(Bxl[2]); ALIGNED(Bxl[3]); ALIGNED(Bxl[4]); ALIGNED(Bxl[5]); ALIGNED(Bxl[6]); ALIGNED(Bxl[7]);
  ALIGNED(Byl[0]); ALIGNED(Byl[1]); ALIGNED(Byl[2]); ALIGNED(Byl[3]); ALIGNED(Byl[4]); ALIGNED(Byl[5]); ALIGNED(Byl[6]); ALIGNED(Byl[7]);
  ALIGNED(Bzl[0]); ALIGNED(Bzl[1]); ALIGNED(Bzl[2]); ALIGNED(Bzl[3]); ALIGNED(Bzl[4]); ALIGNED(Bzl[5]); ALIGNED(Bzl[6]); ALIGNED(Bzl[7]);
  ALIGNED(Exl[0]); ALIGNED(Exl[1]); ALIGNED(Exl[2]); ALIGNED(Exl[3]); ALIGNED(Exl[4]); ALIGNED(Exl[5]); ALIGNED(Exl[6]); ALIGNED(Exl[7]);
  ALIGNED(Eyl[0]); ALIGNED(Eyl[1]); ALIGNED(Eyl[2]); ALIGNED(Eyl[3]); ALIGNED(Eyl[4]); ALIGNED(Eyl[5]); ALIGNED(Eyl[6]); ALIGNED(Eyl[7]);
  ALIGNED(Ezl[0]); ALIGNED(Ezl[1]); ALIGNED(Ezl[2]); ALIGNED(Ezl[3]); ALIGNED(Ezl[4]); ALIGNED(Ezl[5]); ALIGNED(Ezl[6]); ALIGNED(Ezl[7]);
  
  for (c = 0; c < 8; c++) {
    for (j = 0; j < 6; j++) {
      field_components[c][j] = ptr = _mm_malloc(sizeof(double) * NUM_PCLS, ALIGNMENT);
      ALIGNED(field_components[c][j]);
      if (!ptr) { printf("%d: Malloc failed.\n", __LINE__); exit(1); }
    }
  }

  /* Init arrays */
  for (pidx = 0; pidx < num_pcls; pidx++)
  {
    x[pidx] = (double) random();
    y[pidx] = (double) random();
    z[pidx] = (double) random();
    u[pidx] = (double) random();
    v[pidx] = (double) random();
    w[pidx] = (double) random();
    _xavg[pidx] = (double) random();
    _yavg[pidx] = (double) random();
    _zavg[pidx] = (double) random();
  }
  for (c = 0; c < 8; c++) {
    for (j = 0; j < 6; j++) {
      for (pidx = 0; pidx < num_pcls; pidx++)
      field_components[c][j][pidx] = (double) random();
    }
  }

                           
  /* Prepare weights */
  for (pidx = 0; pidx < num_pcls; pidx++)
  {
    abs_pos[0] = _xavg[pidx];
    abs_pos[1] = _yavg[pidx];
    abs_pos[2] = _zavg[pidx];
    // xstart marks start of domain excluding ghosts
    rel_pos[0] = abs_pos[0] - xstart;
    rel_pos[1] = abs_pos[1] - ystart;
    rel_pos[2] = abs_pos[2] - zstart;
    // cell position minus 1 (due to ghost cells)
    cm1_pos[0] = rel_pos[0] * inv_dx;
    cm1_pos[1] = rel_pos[1] * inv_dy;
    cm1_pos[2] = rel_pos[2] * inv_dz;
    // index of interface to right of cell 
    // NOT NEEDED HERE
    //const int ix = cx + 1;
    //const int iy = cy + 1;
    //const int iz = cz + 1;
    // fraction of the distance from the right of the cell
    w1[0] = cx - cm1_pos[0];
    w1[1] = cy - cm1_pos[1];
    w1[2] = cz - cm1_pos[2];
    // fraction of distance from the left
    w0[0] = 1-w1[0];
    w0[1] = 1-w1[1];
    w0[2] = 1-w1[2];  
    weight[0] = w0[0]*w0[1];
    weight[1] = w0[0]*w1[1];
    weight[2] = w1[0]*w0[1];
    weight[3] = w1[0]*w1[1];    
    weights[0][pidx] = weight[0]*w0[2]; // weight000
    weights[1][pidx] = weight[0]*w1[2]; // weight001
    weights[2][pidx] = weight[1]*w0[2]; // weight010
    weights[3][pidx] = weight[1]*w1[2]; // weight011
    weights[4][pidx] = weight[2]*w0[2]; // weight100
    weights[5][pidx] = weight[2]*w1[2]; // weight101
    weights[6][pidx] = weight[3]*w0[2]; // weight110
    weights[7][pidx] = weight[3]*w1[2]; // weight111    
  }
  
  //printf("Weights prepared\n");
  //fflush(stdout);

  /* Init B{x,y,z}l, E{x,y,z}l */
  for (c = 0; c < 8; c++) {
    for (pidx = 0; pidx < num_pcls; pidx++) {
      Bxl[c][pidx] = 0.0;
      Byl[c][pidx] = 0.0;
      Bzl[c][pidx] = 0.0;
      Exl[c][pidx] = 0.0;
      Eyl[c][pidx] = 0.0;
      Ezl[c][pidx] = 0.0;
    }
  }


{
  double *p1 = Bxl[0];
  double *p2 = weights[0];
  double *p3 = field_components[0][0];

#if 0
  // Warm up cache 
#pragma vector aligned
  for (pidx = 0; pidx < num_pcls; pidx++) {
    //Bxl[0][pidx] = weights[0][pidx] * field_components[0][0][pidx];
    p1[pidx] = p2[pidx] * p3[pidx];
  }
#endif

#pragma omp parallel 
{}

   time1 = time_sec();

#pragma omp parallel for  
#pragma vector aligned
  for (pidx = 0; pidx < num_pcls; pidx++) {
    //Bxl[0][pidx] = weights[0][pidx] * field_components[0][0][pidx];
    p1[pidx] = p2[pidx] * p1[pidx];  // Important!!!
  }

  time2 = time_sec();
  printf("Time   : %f\n", time2 - time1);
  printf("GFlops : %f\n", num_pcls / (time2 - time1) / 1e9);
}

  
#if 0
  time1 = time_sec();

  /* Calc B{x,y,z}l[0..7], E{x,y,z}l[0..7] for all pcls */
#pragma ivdep
  for (pidx = 0; pidx < num_pcls; pidx++) {
    Bxl[0][pidx] = weights[0][pidx] * field_components[0][0][pidx];
    Byl[0][pidx] = weights[0][pidx] * field_components[0][1][pidx];
    Bzl[0][pidx] = weights[0][pidx] * field_components[0][2][pidx];
    Exl[0][pidx] = weights[0][pidx] * field_components[0][3][pidx];
    Eyl[0][pidx] = weights[0][pidx] * field_components[0][4][pidx];
    Ezl[0][pidx] = weights[0][pidx] * field_components[0][5][pidx];

    Bxl[1][pidx] = weights[1][pidx] * field_components[1][0][pidx];
    Byl[1][pidx] = weights[1][pidx] * field_components[1][1][pidx];
    Bzl[1][pidx] = weights[1][pidx] * field_components[1][2][pidx];
    Exl[1][pidx] = weights[1][pidx] * field_components[1][3][pidx];
    Eyl[1][pidx] = weights[1][pidx] * field_components[1][4][pidx];
    Ezl[1][pidx] = weights[1][pidx] * field_components[1][5][pidx];

    Bxl[2][pidx] = weights[2][pidx] * field_components[2][0][pidx];
    Byl[2][pidx] = weights[2][pidx] * field_components[2][1][pidx];
    Bzl[2][pidx] = weights[2][pidx] * field_components[2][2][pidx];
    Exl[2][pidx] = weights[2][pidx] * field_components[2][3][pidx];
    Eyl[2][pidx] = weights[2][pidx] * field_components[2][4][pidx];
    Ezl[2][pidx] = weights[2][pidx] * field_components[2][5][pidx];

    Bxl[3][pidx] = weights[3][pidx] * field_components[3][0][pidx];
    Byl[3][pidx] = weights[3][pidx] * field_components[3][1][pidx];
    Bzl[3][pidx] = weights[3][pidx] * field_components[3][2][pidx];
    Exl[3][pidx] = weights[3][pidx] * field_components[3][3][pidx];
    Eyl[3][pidx] = weights[3][pidx] * field_components[3][4][pidx];
    Ezl[3][pidx] = weights[3][pidx] * field_components[3][5][pidx];

    Bxl[4][pidx] = weights[4][pidx] * field_components[4][0][pidx];
    Byl[4][pidx] = weights[4][pidx] * field_components[4][1][pidx];
    Bzl[4][pidx] = weights[4][pidx] * field_components[4][2][pidx];
    Exl[4][pidx] = weights[4][pidx] * field_components[4][3][pidx];
    Eyl[4][pidx] = weights[4][pidx] * field_components[4][4][pidx];
    Ezl[4][pidx] = weights[4][pidx] * field_components[4][5][pidx];

    Bxl[5][pidx] = weights[5][pidx] * field_components[5][0][pidx];
    Byl[5][pidx] = weights[5][pidx] * field_components[5][1][pidx];
    Bzl[5][pidx] = weights[5][pidx] * field_components[5][2][pidx];
    Exl[5][pidx] = weights[5][pidx] * field_components[5][3][pidx];
    Eyl[5][pidx] = weights[5][pidx] * field_components[5][4][pidx];
    Ezl[5][pidx] = weights[5][pidx] * field_components[5][5][pidx];

    Bxl[6][pidx] = weights[6][pidx] * field_components[6][0][pidx];
    Byl[6][pidx] = weights[6][pidx] * field_components[6][1][pidx];
    Bzl[6][pidx] = weights[6][pidx] * field_components[6][2][pidx];
    Exl[6][pidx] = weights[6][pidx] * field_components[6][3][pidx];
    Eyl[6][pidx] = weights[6][pidx] * field_components[6][4][pidx];
    Ezl[6][pidx] = weights[6][pidx] * field_components[6][5][pidx];

    Bxl[7][pidx] = weights[7][pidx] * field_components[7][0][pidx];
    Byl[7][pidx] = weights[7][pidx] * field_components[7][1][pidx];
    Bzl[7][pidx] = weights[7][pidx] * field_components[7][2][pidx];
    Exl[7][pidx] = weights[7][pidx] * field_components[7][3][pidx];
    Eyl[7][pidx] = weights[7][pidx] * field_components[7][4][pidx];
    Ezl[7][pidx] = weights[7][pidx] * field_components[7][5][pidx];
#if 1
    Bxl[0][pidx] = Bxl[0][pidx] + Bxl[1][pidx] + Bxl[2][pidx] + Bxl[3][pidx] + Bxl[4][pidx] + Bxl[5][pidx] + Bxl[6][pidx] + Bxl[7][pidx];
    Byl[0][pidx] = Byl[0][pidx] + Byl[1][pidx] + Byl[2][pidx] + Byl[3][pidx] + Byl[4][pidx] + Byl[5][pidx] + Byl[6][pidx] + Byl[7][pidx];
    Bzl[0][pidx] = Bzl[0][pidx] + Bzl[1][pidx] + Bzl[2][pidx] + Bzl[3][pidx] + Bzl[4][pidx] + Bzl[5][pidx] + Bzl[6][pidx] + Bzl[7][pidx];

    Exl[0][pidx] = Exl[0][pidx] + Exl[1][pidx] + Exl[2][pidx] + Exl[3][pidx] + Exl[4][pidx] + Exl[5][pidx] + Exl[6][pidx] + Exl[7][pidx];
    Eyl[0][pidx] = Eyl[0][pidx] + Eyl[1][pidx] + Eyl[2][pidx] + Eyl[3][pidx] + Eyl[4][pidx] + Eyl[5][pidx] + Eyl[6][pidx] + Eyl[7][pidx];
    Ezl[0][pidx] = Ezl[0][pidx] + Ezl[1][pidx] + Ezl[2][pidx] + Ezl[3][pidx] + Ezl[4][pidx] + Ezl[5][pidx] + Ezl[6][pidx] + Ezl[7][pidx];
#endif
  }
  
#if 0
  /* Sum B{x,y,z}l[0..7], E{x,y,z}l[0..7] for all pcls */
#pragma nofusion
  for (pidx = 0; pidx < num_pcls; pidx++) {
    Bxl[0][pidx] = Bxl[0][pidx] + Bxl[1][pidx] + Bxl[2][pidx] + Bxl[3][pidx] + Bxl[4][pidx] + Bxl[5][pidx] + Bxl[6][pidx] + Bxl[7][pidx];
    Byl[0][pidx] = Byl[0][pidx] + Byl[1][pidx] + Byl[2][pidx] + Byl[3][pidx] + Byl[4][pidx] + Byl[5][pidx] + Byl[6][pidx] + Byl[7][pidx];
    Bzl[0][pidx] = Bzl[0][pidx] + Bzl[1][pidx] + Bzl[2][pidx] + Bzl[3][pidx] + Bzl[4][pidx] + Bzl[5][pidx] + Bzl[6][pidx] + Bzl[7][pidx];
    
    Exl[0][pidx] = Exl[0][pidx] + Exl[1][pidx] + Exl[2][pidx] + Exl[3][pidx] + Exl[4][pidx] + Exl[5][pidx] + Exl[6][pidx] + Exl[7][pidx];
    Eyl[0][pidx] = Eyl[0][pidx] + Eyl[1][pidx] + Eyl[2][pidx] + Eyl[3][pidx] + Eyl[4][pidx] + Eyl[5][pidx] + Eyl[6][pidx] + Eyl[7][pidx];
    Ezl[0][pidx] = Ezl[0][pidx] + Ezl[1][pidx] + Ezl[2][pidx] + Ezl[3][pidx] + Ezl[4][pidx] + Ezl[5][pidx] + Ezl[6][pidx] + Ezl[7][pidx];
  }
#endif

  time2 = time_sec();
  printf("time: %f\n", time2 - time1);
#endif

  /* Do what is left */
  for (pidx = 0; pidx < num_pcls; pidx++) {
    Om[0] = qdto2mc * Bxl[0][pidx];
    Om[1] = qdto2mc * Byl[0][pidx];
    Om[2] = qdto2mc * Bzl[0][pidx];
    
    // end interpolation
    omsq = (Om[0] * Om[0] + Om[1] * Om[1] + Om[2] * Om[2]);
    denom = 1.0 / (1.0 + omsq);    
    
    // solve the position equation
    t[0] = u[pidx] + qdto2mc * Exl[0][pidx];
    t[1] = v[pidx] + qdto2mc * Eyl[0][pidx];
    t[2] = w[pidx] + qdto2mc * Ezl[0][pidx];
    //const pfloat udotb = ut * Bxl + vt * Byl + wt * Bzl;
    udotOm = t[0] * Om[0] + t[1] * Om[1] + t[2] * Om[2];
    
    // solve the velocity equation
    avg[0] = (t[0] + (t[1] * Om[2] - t[2] * Om[1] + udotOm * Om[0])) * denom;
    avg[1] = (t[1] + (t[2] * Om[0] - t[0] * Om[2] + udotOm * Om[1])) * denom;
    avg[2] = (t[2] + (t[0] * Om[1] - t[1] * Om[0] + udotOm * Om[2])) * denom;
    
    // update average position
    _xavg[pidx] = x[pidx] + avg[0] * dto2;
    _yavg[pidx] = y[pidx] + avg[1] * dto2;
    _zavg[pidx] = z[pidx] + avg[2] * dto2;

    /*
     * Need outer loop to do this here.
     * So we skip it.
     */
#if 0
    // if it is the last iteration, update the position and velocity
    // (hopefully this will not compromise vectorization...)
    if(niter==NiterMover)
    {
      x[pidx] = xorig + uavg * dt;
      y[pidx] = yorig + vavg * dt;
      z[pidx] = zorig + wavg * dt;
      u[pidx] = 2.0 * uavg - uorig;
      v[pidx] = 2.0 * vavg - vorig;
      w[pidx] = 2.0 * wavg - worig;
    }
#endif
  }
  
  /* Free */
  _mm_free(x);
  _mm_free(y);
  _mm_free(z);
  _mm_free(u);
  _mm_free(v);
  _mm_free(w);
  _mm_free(_xavg);
  _mm_free(_yavg);
  _mm_free(_zavg);

  for (i = 0; i < 8; i++) {
    _mm_free(weights[i]);
    _mm_free(Bxl[i]);
    _mm_free(Byl[i]);
    _mm_free(Bzl[i]);
    _mm_free(Exl[i]);
    _mm_free(Eyl[i]);
    _mm_free(Ezl[i]);
  }

  for (c = 0; c < 8; c++) {
    for (j = 0; j < 6; j++) {
      _mm_free(field_components[c][j]);
    }
  }

  printf("Finished successfully.\n");
  fflush(stdout);

  return 0;                         
}

#if 0
for(int pidx=bucket_offset; pidx<bucket_end; pidx++)
{
  // copy the particle
  const pfloat xorig = x[pidx];
  const pfloat yorig = y[pidx];
  const pfloat zorig = z[pidx];
  const pfloat uorig = u[pidx];
  const pfloat vorig = v[pidx];
  const pfloat worig = w[pidx];

  // compute weights for field components
  //
  double weights[8];
  const double abs_xpos = _xavg[pidx];
  const double abs_ypos = _yavg[pidx];
  const double abs_zpos = _zavg[pidx];
  // xstart marks start of domain excluding ghosts
  const double rel_xpos = abs_xpos - xstart;
  const double rel_ypos = abs_ypos - ystart;
  const double rel_zpos = abs_zpos - zstart;
  // cell position minus 1 (due to ghost cells)
  const double cxm1_pos = rel_xpos * inv_dx;
  const double cym1_pos = rel_ypos * inv_dy;
  const double czm1_pos = rel_zpos * inv_dz;
  // index of interface to right of cell
  const int ix = cx + 1;
  const int iy = cy + 1;
  const int iz = cz + 1;
  // fraction of the distance from the right of the cell
  const double w1x = cx - cxm1_pos;
  const double w1y = cy - cym1_pos;
  const double w1z = cz - czm1_pos;
  // fraction of distance from the left
  const double w0x = 1-w1x;
  const double w0y = 1-w1y;
  const double w0z = 1-w1z;
  const double weight00 = w0x*w0y;
  const double weight01 = w0x*w1y;
  const double weight10 = w1x*w0y;
  const double weight11 = w1x*w1y;
  weights[0] = weight00*w0z; // weight000
  weights[1] = weight00*w1z; // weight001
  weights[2] = weight01*w0z; // weight010
  weights[3] = weight01*w1z; // weight011
  weights[4] = weight10*w0z; // weight100
  weights[5] = weight10*w1z; // weight101
  weights[6] = weight11*w0z; // weight110
  weights[7] = weight11*w1z; // weight111

  pfloat Exl = 0.0;
  pfloat Eyl = 0.0;
  pfloat Ezl = 0.0;
  pfloat Bxl = 0.0;
  pfloat Byl = 0.0;
  pfloat Bzl = 0.0;

  // would expanding this out help to vectorize?
  for(int c=0; c<8; c++)
  {
    Bxl += weights[c] * field_components[c][0];
    Byl += weights[c] * field_components[c][1];
    Bzl += weights[c] * field_components[c][2];
    Exl += weights[c] * field_components[c][3];
    Eyl += weights[c] * field_components[c][4];
    Ezl += weights[c] * field_components[c][5];
  }
  const double Omx = qdto2mc*Bxl;
  const double Omy = qdto2mc*Byl;
  const double Omz = qdto2mc*Bzl;

  // end interpolation
  const pfloat omsq = (Omx * Omx + Omy * Omy + Omz * Omz);
  const pfloat denom = 1.0 / (1.0 + omsq);
  // solve the position equation
  const pfloat ut = uorig + qdto2mc * Exl;
  const pfloat vt = vorig + qdto2mc * Eyl;
  const pfloat wt = worig + qdto2mc * Ezl;
  //const pfloat udotb = ut * Bxl + vt * Byl + wt * Bzl;
  const pfloat udotOm = ut * Omx + vt * Omy + wt * Omz;
  // solve the velocity equation
  const pfloat uavg = (ut + (vt * Omz - wt * Omy + udotOm * Omx)) * denom;
  const pfloat vavg = (vt + (wt * Omx - ut * Omz + udotOm * Omy)) * denom;
  const pfloat wavg = (wt + (ut * Omy - vt * Omx + udotOm * Omz)) * denom;
  // update average position
  _xavg[pidx] = xorig + uavg * dto2;
  _yavg[pidx] = yorig + vavg * dto2;
  _zavg[pidx] = zorig + wavg * dto2;

  // if it is the last iteration, update the position and velocity
  // (hopefully this will not compromise vectorization...)
  if(niter==NiterMover)
  {
    x[pidx] = xorig + uavg * dt;
    y[pidx] = yorig + vavg * dt;
    z[pidx] = zorig + wavg * dt;
    u[pidx] = 2.0 * uavg - uorig;
    v[pidx] = 2.0 * vavg - vorig;
    w[pidx] = 2.0 * wavg - worig;
  }
}
#endif
