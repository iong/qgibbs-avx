#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <x86intrin.h>
#ifdef _OPENMP
#include <omp.h>
#endif

static int NGAUSS, Natom, *nnb, *nbidx, nnbmax;
float bl, bl2;

static float  LJA[16] __attribute__((aligned(32)));
static float  LJC[16] __attribute__((aligned(32)));



#ifdef __AVX__
static __m256 _mm256_flip_sign_ps(__m256 x)
{
	return _mm256_xor_ps(x, _mm256_set1_ps(-0.0f));
}
#else
static __m128 _mm_flip_sign_ps(__m128 x)
{
	return _mm_xor_ps(x, _mm_set1_ps(-0.0f));
}
#endif

#ifdef __FMA4__
#include "gaussian_average_fma4.h"
#elif __AVX__
#include "gaussian_average_acc.h"
#elif (__SSE4_2__) || (__SSE4_1__)
#include "gaussian_average_sse4_1.h"
#else
#error "Not supported!"
#endif



static void min_image(int n, int stride, float *dq)
{
	int	i, k;

	for (k=0; k<3; k++) {
		for (i=k*stride; i<k*stride+n; i++) {
			if (dq[i] > bl2) {
				dq[i] =  dq[i] - bl;
			}
			else if (-dq[i] > bl2) {
				dq[i] = bl + dq[i];
			}
		}
	}
}


/* vgw.f90 already allocates UPV and UPM to accomodate the maximum number of threads.
 */
void gaussian_average_acc(double *y, double *Uout, double *UPV, double *UPM)
{
	float E3[6] = {1.0, 0.0, 0.0, 1.0, 0.0, 1.0};
	double U=0.0;

	float *dq, *GC, *U0, *UX0, *UXX0;
	int	tid=0;
	int	i;
#ifdef _OPENMP
	int	nthreads = omp_get_num_threads();
	tid = omp_get_thread_num();
#endif

	memset(UPV+3*Natom*tid, 0, Natom*3*sizeof(double));
	memset(UPM+6*Natom*tid, 0, Natom*6*sizeof(double));
#pragma omp master
	{*Uout = 0.0;}

	posix_memalign((void **)&dq  , 32, nnbmax*3*sizeof(float));
	posix_memalign((void **)&GC  , 32, nnbmax*6*sizeof(float));
	posix_memalign((void **)&U0 ,  32, nnbmax*sizeof(float));
	posix_memalign((void **)&UX0 , 32, nnbmax*3*sizeof(float));
	posix_memalign((void **)&UXX0, 32, nnbmax*6*sizeof(float));

	//printf("%d: %x, %x\n", tid, expav, dq);

#pragma omp for schedule(dynamic, 12)
	for (i=0; i<Natom-1; i++) {
		int j, k;
		int NN1 = nnb[i];	
		int NN18 = 8*((NN1 + 7) / 8);

		for (j=0; j<NN1; j++) {
			int ioffset = 3*i;
			int joffset = 3*(nbidx[Natom*i+j] - 1);

			for (k=0; k<3; k++)
				dq[k*nnbmax + j] = y[k + ioffset] - y[k + joffset];

			ioffset = 3*Natom + 6*i;
			joffset = 3*Natom + 6*(nbidx[Natom*i + j]-1);
			for (k=0; k<6; k++)
				GC[k*nnbmax + j] = y[k + ioffset] + y[k + joffset];
		}
		min_image(NN1, nnbmax, dq);

		for (k=0; k<3; k++) 
			for (j=NN1; j<NN18; j++)
				dq[k*nnbmax + j] = 1.0;

		for (k=0; k<6; k++) 
			for (j=NN1; j<NN18; j++)
				GC[k*nnbmax + j] = E3[k];

#if defined (__AVX__) || defined (__FMA4__)
		for (j=0; j < NN1; j += 8) {
#else
		for (j=0; j < NN1; j += 4) {
#endif
			//printf ("%d\n", j);
			vgw_kernel(j, nnbmax, dq, GC, U0, UX0, UXX0);
		}

		for (j=0; j<NN1; j++) {
			int idest = Natom*tid + i;
			int jdest = Natom*tid + nbidx[Natom*i+j] - 1;

			U += U0[j];

			for (k=0; k<3; k++) UPV[k + 3*idest] += UX0[j+k*nnbmax];
			for (k=0; k<3; k++) UPV[k + 3*jdest] -= UX0[j+k*nnbmax];

			for (k=0; k<6; k++) UPM[k + 6*idest] += UXX0[j+k*nnbmax];
			for (k=0; k<6; k++) UPM[k + 6*jdest] += UXX0[j+k*nnbmax];
		}
	}

#pragma omp critical
	{*Uout += U;}

#ifdef _OPENMP
#pragma omp for schedule(static)
	for (i=0; i<Natom; i++) {
		int j, k;
		for (j=1; j<nthreads; j++) {
			for (k=0; k<3; k++) UPV[k+3*i] += UPV[k+3*(i + Natom*j)];
			for (k=0; k<6; k++) UPM[k+6*i] += UPM[k+6*(i + Natom*j)];
		}
	}
#endif

	free(dq);
	free(GC);
	free(U0);
	free(UX0);
	free(UXX0);

// #pragma omp master
	//printf("%lg\n", *Uout);

}


void gaussian_average_acc_init(int iNatom, int *innb, int *inbidx, int innbmax, float *iLJA, float *iLJC, int iNGAUSS, double ibl)
{
	int i;

	Natom = iNatom;
	nnb = innb;
	nbidx = inbidx;
	nnbmax = innbmax;


	NGAUSS = iNGAUSS;
	bl = ibl;
	bl2 = bl / 2.0;

	for (i=0; i<NGAUSS; i++) {
		LJA[i] = iLJA[i];
		LJC[i] = iLJC[i];
	}

}

void gaussian_average_acc_cleanup()
{
}


