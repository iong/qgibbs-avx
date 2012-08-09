#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <immintrin.h>

#ifdef _OPENMP
#include <omp.h>
#endif

static float  expav[8] __attribute__((aligned(32)));
#pragma omp threadprivate(expav)

static int NGAUSS, Natom, *nnb, *nbidx, nnbmax;
float bl, bl2;

static float  LJA[16] __attribute__((aligned(32)));
static float  LJC[16] __attribute__((aligned(32)));


static inline __m256 vmuladd3e(__m256 rx1, __m256 ry1, __m256 rx2, __m256 ry2, __m256 rx3, __m256 ry3)
{
	rx1 = _mm256_mul_ps(rx1, ry1);
	rx2 = _mm256_mul_ps(rx2, ry2);
	rx3 = _mm256_mul_ps(rx3, ry3);

	rx1 = _mm256_add_ps(rx1, rx2);
	rx1 = _mm256_add_ps(rx1, rx3);

	return rx1;
}

static inline __m256 invdetinv_avx(__m256 *a, __m256 *ia)
{
	ia[0] = _mm256_sub_ps(_mm256_mul_ps(a[3], a[5]), _mm256_mul_ps(a[4], a[4]));
	ia[1] = _mm256_sub_ps(_mm256_mul_ps(a[2], a[4]), _mm256_mul_ps(a[1], a[5]));
	ia[2] = _mm256_sub_ps(_mm256_mul_ps(a[1], a[4]), _mm256_mul_ps(a[2], a[3]));
	ia[3] = _mm256_sub_ps(_mm256_mul_ps(a[0], a[5]), _mm256_mul_ps(a[2], a[2]));
	ia[4] = _mm256_sub_ps(_mm256_mul_ps(a[2], a[1]), _mm256_mul_ps(a[0], a[4]));
	ia[5] = _mm256_sub_ps(_mm256_mul_ps(a[0], a[3]), _mm256_mul_ps(a[1], a[1]));

	__m256 m1 =_mm256_mul_ps(ia[0], a[0]);
	__m256 detia = _mm256_mul_ps(ia[1], a[1]);
	detia = _mm256_add_ps(detia, m1);
	detia = _mm256_add_ps(detia, _mm256_mul_ps(ia[2], a[2]));
	detia = _mm256_rcp_ps(detia);

	ia[0] = _mm256_mul_ps(ia[0],detia);
	ia[1] = _mm256_mul_ps(ia[1],detia);
	ia[2] = _mm256_mul_ps(ia[2],detia);
	ia[3] = _mm256_mul_ps(ia[3],detia);
	ia[4] = _mm256_mul_ps(ia[4],detia);
	ia[5] = _mm256_mul_ps(ia[5],detia);

	return detia;
}

static void vgw_kernel_avx(int first, int stride, float *q, float *GC, float *U, float *UX0, float *UXX0)
{
    int IG, i, j;
    __m256 dq[3], gc[6], a[6], deta;

    for (i=0; i<3; i++) dq[i] = _mm256_load_ps( q + first + i*stride);
    for (i=0; i<6; i++) gc[i] = _mm256_load_ps(GC + first + i*stride);

    deta = invdetinv_avx(gc, a);

    __m256 ux0[3], uxx0[6];
    for (i=0; i<3; i++) ux0[i] = _mm256_setzero_ps();
    for (i=0; i<6; i++) uxx0[i] = _mm256_setzero_ps();
    __m256 Ulocal =  _mm256_setzero_ps();


    for(IG=0; IG<NGAUSS; IG++) {
	__m256 lja_ig = _mm256_broadcast_ss(LJA+IG);
	__m256 ag[6];
	ag[0] = _mm256_add_ps(a[0], lja_ig);
	ag[1] = a[1];
	ag[2] = a[2];
	ag[3] = _mm256_add_ps(a[3], lja_ig);
	ag[4] = a[4];
	ag[5] = _mm256_add_ps(a[5], lja_ig);

	__m256 idetag, z[6], Zq[3], qZq;
	idetag = invdetinv_avx(ag, z);

	__m256 lja_igsq = _mm256_sub_ps(_mm256_setzero_ps(), _mm256_mul_ps(lja_ig, lja_ig));

	for (j=0; j<6; j++) z[j] = _mm256_mul_ps(lja_igsq, z[j]);
	z[0] = _mm256_add_ps(z[0], lja_ig);
	z[3] = _mm256_add_ps(z[3], lja_ig);
	z[5] = _mm256_add_ps(z[5], lja_ig);

	//Zq = matmul(Z, Q12) ! R = -2.0*Zq
	Zq[0]=vmuladd3e(z[0], dq[0], z[1], dq[1], z[2], dq[2]);
	Zq[1]=vmuladd3e(z[1], dq[0], z[3], dq[1], z[4], dq[2]);
	Zq[2]=vmuladd3e(z[2], dq[0], z[4], dq[1], z[5], dq[2]);

	//qZq = dot_product(Q12, Zq) 
	qZq = vmuladd3e(Zq[0], dq[0], Zq[1], dq[1], Zq[2], dq[2]);

	_mm256_store_ps(expav, qZq);
	for (i=0; i<8; i++) expav[i] = expf(-expav[i]);

	__m256 two = _mm256_set_ps(2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f, 2.0f);
	__m256 v0 = _mm256_mul_ps(_mm256_mul_ps(two, _mm256_broadcast_ss(LJC+IG)), _mm256_mul_ps(_mm256_load_ps(expav), _mm256_sqrt_ps(_mm256_mul_ps(deta, idetag))));

	Ulocal = _mm256_add_ps(Ulocal, v0);
/*
	__m256 thresh = _mm256_set_ps(1e3, 1e3, 1e3, 1e3, 1e3, 1e3, 1e3, 1e3);
*/
	for (j=0; j<3; j++) {
	    ux0[j] = _mm256_sub_ps(ux0[j], _mm256_mul_ps(v0, Zq[j]));
	    /*
	    if (_mm256_movemask_ps(_mm256_cmp_ps(_mm256_andnot_ps(mask, ux0[j]), thresh, _CMP_GT_OS)) & 0xff) {
		    printf("%d, %d, %d\n", first, IG, j);
	    }
	    */
	}

	uxx0[0] = _mm256_add_ps(uxx0[0], _mm256_mul_ps(v0, _mm256_sub_ps(_mm256_mul_ps(two, _mm256_mul_ps(Zq[0], Zq[0])), z[0])));
	uxx0[1] = _mm256_add_ps(uxx0[1], _mm256_mul_ps(v0, _mm256_sub_ps(_mm256_mul_ps(two, _mm256_mul_ps(Zq[1], Zq[0])), z[1])));
	uxx0[2] = _mm256_add_ps(uxx0[2], _mm256_mul_ps(v0, _mm256_sub_ps(_mm256_mul_ps(two, _mm256_mul_ps(Zq[2], Zq[0])), z[2])));
	uxx0[3] = _mm256_add_ps(uxx0[3], _mm256_mul_ps(v0, _mm256_sub_ps(_mm256_mul_ps(two, _mm256_mul_ps(Zq[1], Zq[1])), z[3])));
	uxx0[4] = _mm256_add_ps(uxx0[4], _mm256_mul_ps(v0, _mm256_sub_ps(_mm256_mul_ps(two, _mm256_mul_ps(Zq[2], Zq[1])), z[4])));
	uxx0[5] = _mm256_add_ps(uxx0[5], _mm256_mul_ps(v0, _mm256_sub_ps(_mm256_mul_ps(two, _mm256_mul_ps(Zq[2], Zq[2])), z[5])));

    }
    __m256 half = _mm256_set_ps(0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f);
    _mm256_store_ps(U+first, _mm256_mul_ps(half, Ulocal));

    for (i=0; i<3; i++) _mm256_store_ps( UX0+first + i*stride,  ux0[i]);
    for (i=0; i<6; i++) _mm256_store_ps(UXX0+first + i*stride, uxx0[i]);
}

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
void gaussian_average_avx(double *y, double *Uout, double *UPV, double *UPM)
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

#pragma omp for schedule(dynamic, 16)
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

		for (j=0; j < NN1; j += 8) {
			//printf ("%d\n", j);
			vgw_kernel_avx(j, nnbmax, dq, GC, U0, UX0, UXX0);
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


void gaussian_average_avx_init(int iNatom, int *innb, int *inbidx, int innbmax, float *iLJA, float *iLJC, int iNGAUSS, double ibl)
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

void gaussian_average_avx_cleanup()
{
}


