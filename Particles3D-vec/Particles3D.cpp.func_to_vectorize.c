#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <assert.h>

#define CACHE_SIZE (10*1024*1024) // 1 MB
#define ITERS      100000
#define NUM_PCLS   2048 
#define ALIGNMENT  64
#define BLOCKSIZE  32

#define pfloat double
#define ALIGNED(X) __assume_aligned(X, ALIGNMENT)

inline double time_sec();
void flush_cache();
void move_bucket_old();
void move_bucket_new();
void move_bucket_new_blocked();

/* 
 * Global variables
 */

// Array to flush caches
char volatile *cache;

// Member variables (random values)
double xstart = 0, ystart = 0, zstart = 0;
double inv_dx = .25, inv_dy = .25, inv_dz = .25;
double dto2 = .4;
double cx = 1, cy = 2, cz = 3;
double qdto2mc = .123;

// Particles
int num_pcls = NUM_PCLS;
double *x, *y, *z, *u, *v, *w;
double *restrict _xavg, *restrict _yavg, *restrict _zavg;

double *field_components[8];


int main(void)
{
  int pidx, c, i, j;
  double *ptr;

  /* Array to flush cache */
  cache = _mm_malloc(sizeof(char) * CACHE_SIZE, ALIGNMENT);
  if (!cache) { printf("%d: Malloc failed.\n", __LINE__); exit(1); }
  
  /* Particles in bucket */
  x = ptr = _mm_malloc(sizeof(double) * NUM_PCLS, ALIGNMENT);
  if (!ptr) { printf("%d: Malloc failed.\n", __LINE__); exit(1); }
  y = ptr = _mm_malloc(sizeof(double) * NUM_PCLS, ALIGNMENT);
  if (!ptr) { printf("%d: Malloc failed.\n", __LINE__); exit(1); }
  z = ptr = _mm_malloc(sizeof(double) * NUM_PCLS, ALIGNMENT);
  if (!ptr) { printf("%d: Malloc failed.\n", __LINE__); exit(1); }
  u = ptr = _mm_malloc(sizeof(double) * NUM_PCLS, ALIGNMENT);
  if (!ptr) { printf("%d: Malloc failed.\n", __LINE__); exit(1); }
  v = ptr = _mm_malloc(sizeof(double) * NUM_PCLS, ALIGNMENT);
  if (!ptr) { printf("%d: Malloc failed.\n", __LINE__); exit(1); }
  w = ptr = _mm_malloc(sizeof(double) * NUM_PCLS, ALIGNMENT);
  if (!ptr) { printf("%d: Malloc failed.\n", __LINE__); exit(1); }
  _xavg = ptr = _mm_malloc(sizeof(double) * NUM_PCLS, ALIGNMENT);
  if (!ptr) { printf("%d: Malloc failed.\n", __LINE__); exit(1); }
  _yavg = ptr = _mm_malloc(sizeof(double) * NUM_PCLS, ALIGNMENT);
  if (!ptr) { printf("%d: Malloc failed.\n", __LINE__); exit(1); }
  _zavg = ptr = _mm_malloc(sizeof(double) * NUM_PCLS, ALIGNMENT);
  if (!ptr) { printf("%d: Malloc failed.\n", __LINE__); exit(1); }
  ALIGNED(x); ALIGNED(y); ALIGNED(z);
  ALIGNED(u); ALIGNED(v); ALIGNED(w);
  ALIGNED(_xavg); ALIGNED(_yavg); ALIGNED(_zavg);
  
  /* Fields for bucket */
  for (c = 0; c < 8; c++) {
      field_components[c] = ptr = _mm_malloc(sizeof(double) * 6, ALIGNMENT);
      ALIGNED(field_components[c]);
      if (!ptr) { printf("%d: Malloc failed.\n", __LINE__); exit(1); }
    }

  /* Init particles */
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
  /* Init fields */
  for (c = 0; c < 8; c++) {
    for (j = 0; j < 6; j++) {
      field_components[c][j] = (double) random();
    }
  }

  //move_bucket_old();
  move_bucket_new_blocked();
  /*
  flush_cache();
  flush_cache();
  move_bucket_old();
  flush_cache();
  //move_bucket_new();
  flush_cache();
  move_bucket_new_blocked();
  flush_cache();
  flush_cache();
  move_bucket_old();
  flush_cache();
  flush_cache();
  move_bucket_new_blocked();
  flush_cache();
  flush_cache();
  move_bucket_old();
*/
  
  printf("Finished successfully.\n");

  return 0;                         
}

void move_bucket_new_blocked()
{
    int iters = ITERS;
    int blocksize = BLOCKSIZE; // number particles per block
    int vec_length = 8;
    int num_threads = omp_get_max_threads();
    int num_vecs, rest_num_pcls;
    int tid_num_pcls[num_threads];
    int tid_start_pidx[num_threads];
    int tid_end_pidx[num_threads];
    double time[2];

    // Num particles must be multiple of 128
    assert(num_pcls%128 == 0);

    // Make sure that we run at most 2 threads 
    // as the calculation is only valid for < 3 threads
    assert(num_threads < 3);

    // Create threads in advance for timing
    #pragma omp parallel 
    {}

    time[0] = time_sec();

    num_vecs = num_pcls / vec_length;
    rest_num_pcls = num_pcls % vec_length;

    tid_num_pcls[0] = (num_vecs/2 + num_vecs%2) * vec_length;
    tid_num_pcls[1] = num_vecs/2 * vec_length + rest_num_pcls;

    tid_start_pidx[0] = 0;
    tid_start_pidx[1] = tid_num_pcls[0];

    tid_end_pidx[0] = tid_start_pidx[0] + tid_num_pcls[0]-1;
    tid_end_pidx[1] = tid_start_pidx[1] + tid_num_pcls[1]-1;

    // If only one thread, then it gets all particles
    if (1 == num_threads) {
        tid_num_pcls[0] = num_pcls;
        tid_end_pidx[0] = num_pcls - 1;
    }
  
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        int start_pidx = tid_start_pidx[tid];
        int end_pidx = tid_end_pidx[tid];
        int num_pcls = tid_num_pcls[tid];
        int num_blocks, block, block_pidx;
        int i, j, pidx;
        // For weights calculation
        double abs_pos[3];
        double rel_pos[3];
        double cm1_pos[3];
        double w0[3], w1[3], weight[4];
        double *weights[8];
        __attribute__ ((align(64))) double weights_0[BLOCKSIZE];
        __attribute__ ((align(64))) double weights_1[BLOCKSIZE];
        __attribute__ ((align(64))) double weights_2[BLOCKSIZE];
        __attribute__ ((align(64))) double weights_3[BLOCKSIZE];
        __attribute__ ((align(64))) double weights_4[BLOCKSIZE];
        __attribute__ ((align(64))) double weights_5[BLOCKSIZE];
        __attribute__ ((align(64))) double weights_6[BLOCKSIZE];
        __attribute__ ((align(64))) double weights_7[BLOCKSIZE];
        // For calculating BxyzExyz_l
        double *field_components_0, *field_components_1, *field_components_2;
        double *field_components_3, *field_components_4, *field_components_5;
        double *field_components_6, *field_components_7;
        double *ptr;
        double *BxyzExyz_l[6];
        __attribute__ ((align(64))) double Bxl[BLOCKSIZE];
        __attribute__ ((align(64))) double Byl[BLOCKSIZE];
        __attribute__ ((align(64))) double Bzl[BLOCKSIZE];
        __attribute__ ((align(64))) double Exl[BLOCKSIZE];
        __attribute__ ((align(64))) double Eyl[BLOCKSIZE];
        __attribute__ ((align(64))) double Ezl[BLOCKSIZE];
        // For what is left
        double qdto2mc;
        double omsq, denom;
        double t[3], Om[3], avg[3];
        double udotOm;
        // Timing
        double time[6];


        //printf("tid %d: pcls: [ %d , %d ] total: %d\n", tid, start_pidx, end_pidx, num_pcls);

        //time[0] = time_sec();


        BxyzExyz_l[0] = Bxl; BxyzExyz_l[1] = Byl; BxyzExyz_l[2] = Bzl;
        BxyzExyz_l[3] = Exl; BxyzExyz_l[4] = Eyl; BxyzExyz_l[5] = Ezl;

        field_components_0 = field_components[0]; ALIGNED(field_components_0);
        field_components_1 = field_components[1]; ALIGNED(field_components_1);
        field_components_2 = field_components[2]; ALIGNED(field_components_2);
        field_components_3 = field_components[3]; ALIGNED(field_components_3);
        field_components_4 = field_components[4]; ALIGNED(field_components_4);
        field_components_5 = field_components[5]; ALIGNED(field_components_5);
        field_components_6 = field_components[6]; ALIGNED(field_components_6);
        field_components_7 = field_components[7]; ALIGNED(field_components_7);

        // Number particle blocks for this thread
        num_blocks = (num_pcls-1) / BLOCKSIZE + 1;
        //printf("tid %d: block_pidx: %d\n", tid, block_pidx);

        //time[1] = time_sec();

        // For all blocks
        for (block = 0; block < num_blocks; block++)
        {
            block_pidx = start_pidx + block * BLOCKSIZE;
            //printf("tid %d: block_pidx: %d\n", tid, block_pidx);
            
            /* Compute weights for field components */
            for (i = 0, pidx = block_pidx; i < BLOCKSIZE; i++, pidx++)
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
                weights_0[i] = weight[0]*w0[2]; // weight000
                weights_1[i] = weight[0]*w1[2]; // weight001
                weights_2[i] = weight[1]*w0[2]; // weight010
                weights_3[i] = weight[1]*w1[2]; // weight011
                weights_4[i] = weight[2]*w0[2]; // weight100
                weights_5[i] = weight[2]*w1[2]; // weight101
                weights_6[i] = weight[3]*w0[2]; // weight110
                weights_7[i] = weight[3]*w1[2]; // weight111    
            }

            //printf("%d: weights done\n", block); fflush(stdout);

            //time[2] = time_sec();

            // For Bxyz Exyz
            for (i = 0; i < 6; i++) {
                ptr = BxyzExyz_l[i]; // iterates over B{x,y,z}l E{x,y,z}l
                ALIGNED(ptr);

                // For all particles in block
                #pragma ivdep
                #pragma noprefetch
                for (j = 0; j < BLOCKSIZE; j++)
                {
                    // For all 8 components of one particle
                    ptr[j] = 0.0;
                    ptr[j] += weights_0[j] * field_components_0[i];
                    ptr[j] += weights_1[j] * field_components_1[i];
                    ptr[j] += weights_2[j] * field_components_2[i];
                    ptr[j] += weights_3[j] * field_components_3[i];
                    ptr[j] += weights_4[j] * field_components_4[i];
                    ptr[j] += weights_5[j] * field_components_5[i];
                    ptr[j] += weights_6[j] * field_components_6[i];
                    ptr[j] += weights_7[j] * field_components_7[i];
                }
            }
            //printf("%d: BE done\n", block); fflush(stdout);

            //time[3] = time_sec();

            /* Do what is left */
            //#pragma vector nontemporal (_xavg,_yavg,_zavg)
            for (i = 0, pidx = block_pidx; i < BLOCKSIZE; i++, pidx++) {
                Om[0] = qdto2mc * Bxl[pidx];
                Om[1] = qdto2mc * Byl[pidx];
                Om[2] = qdto2mc * Bzl[pidx];

                // end interpolation
                omsq = (Om[0] * Om[0] + Om[1] * Om[1] + Om[2] * Om[2]);
                denom = 1.0 / (1.0 + omsq);    

                // solve the position equation
                t[0] = u[pidx] + qdto2mc * Exl[pidx];
                t[1] = v[pidx] + qdto2mc * Eyl[pidx];
                t[2] = w[pidx] + qdto2mc * Ezl[pidx];
                //const pfloat udotb = ut * Bxl + vt * Byl + wt * Bzl;
                udotOm = t[0] * Om[0] + t[1] * Om[1] + t[2] * Om[2];

                // solve the velocity equation
                avg[0] = (t[0] + (t[1] * Om[2] - t[2] * Om[1] + udotOm * Om[0])) * denom;
                avg[1] = (t[1] + (t[2] * Om[0] - t[0] * Om[2] + udotOm * Om[1])) * denom;
                avg[2] = (t[2] + (t[0] * Om[1] - t[1] * Om[0] + udotOm * Om[2])) * denom;

                //printf("%d: rest: pidx: %d\n", block, pidx); fflush(stdout);
                // update average position
                _xavg[pidx] = x[pidx] + avg[0] * dto2;
                _yavg[pidx] = y[pidx] + avg[1] * dto2;
                _zavg[pidx] = z[pidx] + avg[2] * dto2;
                

#if 0
                /*
                 * Need outer loop to do this here.
                 * So we skip it.
                 */

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
            //printf("block %d finished\n", block); fflush(stdout);
        }
        //time[4] = time_sec();
        //printf("tid: %d  malloc: %f  weights: %f BE: %f rest: %f\n", tid, time[1]-time[0], time[2]-time[1], time[3]-time[2], time[4]-time[3]);
    }
    time[1] = time_sec();
    printf("move_bucket_new_blocked: total: %f\n", time[1] - time[0]); 
}

void move_bucket_new()
{
    int iters = ITERS;
    int vec_length = 8;
    int num_threads = omp_get_max_threads();
    int num_vecs, rest_num_pcls;
    int tid_num_pcls[num_threads];
    int tid_start_pidx[num_threads];
    int tid_end_pidx[num_threads];
    double time[2];

    // Make sure that we run at most 2 threads 
    // as the calculation is only valid for < 3 threads
    assert(num_threads < 3);

    num_vecs = num_pcls / vec_length;
    rest_num_pcls = num_pcls % vec_length;

    tid_num_pcls[0] = (num_vecs/2 + num_vecs%2) * vec_length;
    tid_num_pcls[1] = num_vecs/2 * vec_length + rest_num_pcls;

    tid_start_pidx[0] = 0;
    tid_start_pidx[1] = tid_num_pcls[0];

    tid_end_pidx[0] = tid_start_pidx[0] + tid_num_pcls[0]-1;
    tid_end_pidx[1] = tid_start_pidx[1] + tid_num_pcls[1]-1;

    // If only one thread, then it gets all particles
    if (1 == num_threads) {
        tid_num_pcls[0] = num_pcls;
        tid_end_pidx[0] = num_pcls - 1;
    }
  
    // Create threads in advance for timing
    #pragma omp parallel 
    {}

    time[0] = time_sec();

    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        int start_pidx = tid_start_pidx[tid];
        int end_pidx = tid_end_pidx[tid];
        int num_pcls = tid_num_pcls[tid];
        int i, j, pidx;
        // For weights calculation
        double abs_pos[3];
        double rel_pos[3];
        double cm1_pos[3];
        double w0[3], w1[3], weight[4];
        double *weights[8];
        double *restrict weights_0;
        double *restrict weights_1; 
        double *restrict weights_2;
        double *restrict weights_3;
        double *restrict weights_4;
        double *restrict weights_5;
        double *restrict weights_6;
        double *restrict weights_7;
        // For calculating BxyzExyz_l
        double *field_components_0, *field_components_1, *field_components_2;
        double *field_components_3, *field_components_4, *field_components_5;
        double *field_components_6, *field_components_7;
        double *ptr;
        double *BxyzExyz_l[6];
        double *Bxl, *Byl, *Bzl;
        double *Exl, *Eyl, *Ezl;
        // For what is left
        double qdto2mc;
        double omsq, denom;
        double t[3], Om[3], avg[3];
        double udotOm;
        // Timing
        double time[6];


        //printf("tid %d: pcls: [ %d , %d ] total: %d\n", tid, start_pidx, end_pidx, num_pcls);

        time[0] = time_sec();

        // Weights for my particles for their 8 components
        for (i = 0; i < 8; i++)
        {
            weights[i] = ptr =_mm_malloc(sizeof(double) * num_pcls, ALIGNMENT);
            if (!ptr) { printf("%d: Malloc failed.\n", __LINE__); exit(1); }
        }
        weights_0 = weights[0]; ALIGNED(weights_0);
        weights_1 = weights[1]; ALIGNED(weights_1);
        weights_2 = weights[2]; ALIGNED(weights_2);
        weights_3 = weights[3]; ALIGNED(weights_3);
        weights_4 = weights[4]; ALIGNED(weights_4);
        weights_5 = weights[5]; ALIGNED(weights_5);
        weights_6 = weights[6]; ALIGNED(weights_6);
        weights_7 = weights[7]; ALIGNED(weights_7);

        // Arrays for BxyzExyz_l
        Bxl = ptr = _mm_malloc(sizeof(double) * num_pcls, ALIGNMENT);
        if (!ptr) { printf("%d: Malloc failed.\n", __LINE__); exit(1); }
        Byl = ptr = _mm_malloc(sizeof(double) * num_pcls, ALIGNMENT);
        if (!ptr) { printf("%d: Malloc failed.\n", __LINE__); exit(1); }
        Bzl = ptr = _mm_malloc(sizeof(double) * num_pcls, ALIGNMENT);
        if (!ptr) { printf("%d: Malloc failed.\n", __LINE__); exit(1); }
        Exl = ptr = _mm_malloc(sizeof(double) * num_pcls, ALIGNMENT);
        if (!ptr) { printf("%d: Malloc failed.\n", __LINE__); exit(1); }
        Eyl = ptr = _mm_malloc(sizeof(double) * num_pcls, ALIGNMENT);
        if (!ptr) { printf("%d: Malloc failed.\n", __LINE__); exit(1); }
        Ezl = ptr = _mm_malloc(sizeof(double) * num_pcls, ALIGNMENT);
        if (!ptr) { printf("%d: Malloc failed.\n", __LINE__); exit(1); }
        ALIGNED(Bxl); ALIGNED(Byl); ALIGNED(Bzl); 
        ALIGNED(Exl); ALIGNED(Eyl); ALIGNED(Ezl); 
        BxyzExyz_l[0] = Bxl; BxyzExyz_l[1] = Byl; BxyzExyz_l[2] = Bzl;
        BxyzExyz_l[3] = Exl; BxyzExyz_l[4] = Eyl; BxyzExyz_l[5] = Ezl;


        field_components_0 = field_components[0]; ALIGNED(field_components_0);
        field_components_1 = field_components[1]; ALIGNED(field_components_1);
        field_components_2 = field_components[2]; ALIGNED(field_components_2);
        field_components_3 = field_components[3]; ALIGNED(field_components_3);
        field_components_4 = field_components[4]; ALIGNED(field_components_4);
        field_components_5 = field_components[5]; ALIGNED(field_components_5);
        field_components_6 = field_components[6]; ALIGNED(field_components_6);
        field_components_7 = field_components[7]; ALIGNED(field_components_7);

        time[1] = time_sec();

        /* Compute weights for field components */
        for (i = 0, pidx = start_pidx; i < num_pcls; i++, pidx++)
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
            weights_0[i] = weight[0]*w0[2]; // weight000
            weights_1[i] = weight[0]*w1[2]; // weight001
            weights_2[i] = weight[1]*w0[2]; // weight010
            weights_3[i] = weight[1]*w1[2]; // weight011
            weights_4[i] = weight[2]*w0[2]; // weight100
            weights_5[i] = weight[2]*w1[2]; // weight101
            weights_6[i] = weight[3]*w0[2]; // weight110
            weights_7[i] = weight[3]*w1[2]; // weight111    
        }

        time[2] = time_sec();

        // For Bxyz Exyz
        for (i = 0; i < 6; i++) {
            ptr = BxyzExyz_l[i]; // iterates over B{x,y,z}l E{x,y,z}l
            ALIGNED(ptr);

            // For all particles
            // Note: Might be better to do this in batches of particles for better cache re-use
            #pragma ivdep
            for (j = 0; j < num_pcls; j++)
            {
                // For all 8 components of one particle
                ptr[j] = 0.0;
                ptr[j] += weights_0[j] * field_components_0[i];
                ptr[j] += weights_1[j] * field_components_1[i];
                ptr[j] += weights_2[j] * field_components_2[i];
                ptr[j] += weights_3[j] * field_components_3[i];
                ptr[j] += weights_4[j] * field_components_4[i];
                ptr[j] += weights_5[j] * field_components_5[i];
                ptr[j] += weights_6[j] * field_components_6[i];
                ptr[j] += weights_7[j] * field_components_7[i];
            }
        }

        time[3] = time_sec();

        /* Do what is left */
        for (i = 0, pidx = start_pidx; i < num_pcls; i++, pidx++) {
            Om[0] = qdto2mc * Bxl[pidx];
            Om[1] = qdto2mc * Byl[pidx];
            Om[2] = qdto2mc * Bzl[pidx];

            // end interpolation
            omsq = (Om[0] * Om[0] + Om[1] * Om[1] + Om[2] * Om[2]);
            denom = 1.0 / (1.0 + omsq);    

            // solve the position equation
            t[0] = u[pidx] + qdto2mc * Exl[pidx];
            t[1] = v[pidx] + qdto2mc * Eyl[pidx];
            t[2] = w[pidx] + qdto2mc * Ezl[pidx];
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

#if 0
            /*
             * Need outer loop to do this here.
             * So we skip it.
             */

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
        time[4] = time_sec();
        printf("tid: %d  malloc: %f  weights: %f BE: %f rest: %f\n", tid, time[1]-time[0], time[2]-time[1], time[3]-time[2], time[4]-time[3]);
    }
    time[1] = time_sec();
    printf("move_bucket_new: total: %f\n", time[1] - time[0]); 
}

void move_bucket_old()
{
    int i, c, pidx;
    double time[2];

    // Create threads in advance for timing
    #pragma omp parallel 
    {}

    time[0] = time_sec();

    #pragma omp parallel for private(i, c)
    for(pidx = 0; pidx < num_pcls; pidx++)
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
        for(c=0; c<8; c++)
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
    time[1] = time_sec();
    printf("move_bucket_old:         total: %f\n", time[1]-time[0]);
}

inline double time_sec()
{
  static struct timeval tv;

  gettimeofday(&tv, NULL);

  return (tv.tv_sec + tv.tv_usec * (double) 1e-6);
}

void flush_cache()
{
    int i;

    for (i = 0; i < CACHE_SIZE; i++)
        cache[i] = (cache[i] + i) % 256;
}

