#include "gmx_math_x86_sse4_1_single.h"


static inline __m128 vmuladd3e(__m128 rx1, __m128 ry1, __m128 rx2, __m128 ry2, __m128 rx3, __m128 ry3)
{
	rx1 = _mm_mul_ps(rx1, ry1);
	rx2 = _mm_mul_ps(rx2, ry2);
	rx3 = _mm_mul_ps(rx3, ry3);

	rx1 = _mm_add_ps(rx1, rx2);
	rx1 = _mm_add_ps(rx1, rx3);

	return rx1;
}

static inline __m128 invdetinv(__m128 *a, __m128 *ia)
{
	ia[0] = _mm_sub_ps(_mm_mul_ps(a[3], a[5]), _mm_mul_ps(a[4], a[4]));
	ia[1] = _mm_sub_ps(_mm_mul_ps(a[2], a[4]), _mm_mul_ps(a[1], a[5]));
	ia[2] = _mm_sub_ps(_mm_mul_ps(a[1], a[4]), _mm_mul_ps(a[2], a[3]));
	ia[3] = _mm_sub_ps(_mm_mul_ps(a[0], a[5]), _mm_mul_ps(a[2], a[2]));
	ia[4] = _mm_sub_ps(_mm_mul_ps(a[2], a[1]), _mm_mul_ps(a[0], a[4]));
	ia[5] = _mm_sub_ps(_mm_mul_ps(a[0], a[3]), _mm_mul_ps(a[1], a[1]));

	__m128 m1 =_mm_mul_ps(ia[0], a[0]);
	__m128 detia = _mm_mul_ps(ia[1], a[1]);
	detia = _mm_add_ps(detia, m1);
	detia = _mm_add_ps(detia, _mm_mul_ps(ia[2], a[2]));
	detia = _mm_rcp_ps(detia);

	ia[0] = _mm_mul_ps(ia[0],detia);
	ia[1] = _mm_mul_ps(ia[1],detia);
	ia[2] = _mm_mul_ps(ia[2],detia);
	ia[3] = _mm_mul_ps(ia[3],detia);
	ia[4] = _mm_mul_ps(ia[4],detia);
	ia[5] = _mm_mul_ps(ia[5],detia);

	return detia;
}

static void vgw_kernel(int first, int stride, float *q, float *GC, float *U, float *UX0, float *UXX0)
{
    int IG, i, j;
    __m128 dq[3], gc[6], a[6], deta;

    for (i=0; i<3; i++) dq[i] = _mm_load_ps( q + first + i*stride);
    for (i=0; i<6; i++) gc[i] = _mm_load_ps(GC + first + i*stride);

    deta = invdetinv(gc, a);

    __m128 ux0[3], uxx0[6];
    for (i=0; i<3; i++) ux0[i] = _mm_setzero_ps();
    for (i=0; i<6; i++) uxx0[i] = _mm_setzero_ps();
    __m128 Ulocal =  _mm_setzero_ps();


    for(IG=0; IG<NGAUSS; IG++) {
	__m128 lja_ig = _mm_set1_ps(LJA[IG]);
	__m128 ag[6];
	ag[0] = _mm_add_ps(a[0], lja_ig);
	ag[1] = a[1];
	ag[2] = a[2];
	ag[3] = _mm_add_ps(a[3], lja_ig);
	ag[4] = a[4];
	ag[5] = _mm_add_ps(a[5], lja_ig);

	__m128 idetag, z[6], Zq[3], qZq;
	idetag = invdetinv(ag, z);

	__m128 lja_igsq = _mm_sub_ps(_mm_setzero_ps(), _mm_mul_ps(lja_ig, lja_ig));

	for (j=0; j<6; j++) z[j] = _mm_mul_ps(lja_igsq, z[j]);
	z[0] = _mm_add_ps(z[0], lja_ig);
	z[3] = _mm_add_ps(z[3], lja_ig);
	z[5] = _mm_add_ps(z[5], lja_ig);

	//Zq = matmul(Z, Q12) ! R = -2.0*Zq
	Zq[0]=vmuladd3e(z[0], dq[0], z[1], dq[1], z[2], dq[2]);
	Zq[1]=vmuladd3e(z[1], dq[0], z[3], dq[1], z[4], dq[2]);
	Zq[2]=vmuladd3e(z[2], dq[0], z[4], dq[1], z[5], dq[2]);

	//qZq = dot_product(Q12, Zq) 
	qZq = vmuladd3e(Zq[0], dq[0], Zq[1], dq[1], Zq[2], dq[2]);

	__m128 two = _mm_set1_ps(2.0f);
	__m128 v0 = _mm_mul_ps(two, gmx_mm_exp_ps(_mm_flip_sign_ps(qZq)));
        v0 = _mm_mul_ps(_mm_mul_ps(v0, _mm_set1_ps(LJC[IG])), _mm_sqrt_ps(_mm_mul_ps(deta, idetag)));

	Ulocal = _mm_add_ps(Ulocal, v0);
	for (j=0; j<3; j++) {
	    ux0[j] = _mm_sub_ps(ux0[j], _mm_mul_ps(v0, Zq[j]));
	}

	uxx0[0] = _mm_add_ps(uxx0[0], _mm_mul_ps(v0, _mm_sub_ps(_mm_mul_ps(two, _mm_mul_ps(Zq[0], Zq[0])), z[0])));
	uxx0[1] = _mm_add_ps(uxx0[1], _mm_mul_ps(v0, _mm_sub_ps(_mm_mul_ps(two, _mm_mul_ps(Zq[1], Zq[0])), z[1])));
	uxx0[2] = _mm_add_ps(uxx0[2], _mm_mul_ps(v0, _mm_sub_ps(_mm_mul_ps(two, _mm_mul_ps(Zq[2], Zq[0])), z[2])));
	uxx0[3] = _mm_add_ps(uxx0[3], _mm_mul_ps(v0, _mm_sub_ps(_mm_mul_ps(two, _mm_mul_ps(Zq[1], Zq[1])), z[3])));
	uxx0[4] = _mm_add_ps(uxx0[4], _mm_mul_ps(v0, _mm_sub_ps(_mm_mul_ps(two, _mm_mul_ps(Zq[2], Zq[1])), z[4])));
	uxx0[5] = _mm_add_ps(uxx0[5], _mm_mul_ps(v0, _mm_sub_ps(_mm_mul_ps(two, _mm_mul_ps(Zq[2], Zq[2])), z[5])));

    }
    _mm_store_ps(U+first, _mm_mul_ps(_mm_set1_ps(0.5f), Ulocal));

    for (i=0; i<3; i++) _mm_store_ps( UX0+first + i*stride,  ux0[i]);
    for (i=0; i<6; i++) _mm_store_ps(UXX0+first + i*stride, uxx0[i]);
}


