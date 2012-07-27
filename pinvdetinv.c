#include <stdio.h>
#include <x86intrin.h>

void pinvdetinv_avx(int N, int stride, float * a, float * inva, float * detinva)
{
	int	i, N8;
	//printf("%d, %d, %x, %x, %x\n", N, stride, (int)a, (int)inva, (int)detinva);
	N8 = ((N+7)/8) * 8;
	for (i=N; i<N8; i++) a[i] = 1.0;
	for (i=N; i<N8; i++) a[i +   stride] = 0.0;
	for (i=N; i<N8; i++) a[i + 2*stride] = 0.0;
	for (i=N; i<N8; i++) a[i + 3*stride] = 1.0;
	for (i=N; i<N8; i++) a[i + 4*stride] = 0.0;
	for (i=N; i<N8; i++) a[i + 5*stride] = 1.0;

	for (i = 0; i<N8; i += 8) {
		//printf("%d\n", i);
		__m256 a0 = _mm256_load_ps(&(a[i]));
		__m256 a1 = _mm256_load_ps(&(a[i +   stride]));
		__m256 a2 = _mm256_load_ps(&(a[i + 2*stride])); 
		__m256 a3 = _mm256_load_ps(&(a[i + 3*stride])); 
		__m256 a4 = _mm256_load_ps(&(a[i + 4*stride])); 
		__m256 a5 = _mm256_load_ps(&(a[i + 5*stride]));

		__m256 ia0 = _mm256_sub_ps(_mm256_mul_ps(a3, a5), _mm256_mul_ps(a4, a4));
		__m256 ia1 = _mm256_sub_ps(_mm256_mul_ps(a2, a4), _mm256_mul_ps(a1, a5));
		__m256 ia2 = _mm256_sub_ps(_mm256_mul_ps(a1, a4), _mm256_mul_ps(a2, a3));
		__m256 ia3 = _mm256_sub_ps(_mm256_mul_ps(a0, a5), _mm256_mul_ps(a2, a2));
		__m256 ia4 = _mm256_sub_ps(_mm256_mul_ps(a2, a1), _mm256_mul_ps(a0, a4));
		__m256 ia5 = _mm256_sub_ps(_mm256_mul_ps(a0, a3), _mm256_mul_ps(a1, a1));

		__m256 detia = _mm256_add_ps(_mm256_mul_ps(ia0, a0), _mm256_mul_ps(ia1, a1));
		detia = _mm256_add_ps(detia, _mm256_mul_ps(ia2, a2));
		detia = _mm256_rcp_ps(detia);

		ia0 = _mm256_mul_ps(ia0, detia);
		ia1 = _mm256_mul_ps(ia1, detia);
		ia2 = _mm256_mul_ps(ia2, detia);
		ia3 = _mm256_mul_ps(ia3, detia);
		ia4 = _mm256_mul_ps(ia4, detia);
		ia5 = _mm256_mul_ps(ia5, detia);


		_mm256_store_ps(&(detinva[i]), detia);
		_mm256_store_ps(&(inva[i]), ia0); 
		_mm256_store_ps(&(inva[i +   stride]), ia1); 
		_mm256_store_ps(&(inva[i + 2*stride]), ia2); 
		_mm256_store_ps(&(inva[i + 3*stride]), ia3); 
		_mm256_store_ps(&(inva[i + 4*stride]), ia4); 
		_mm256_store_ps(&(inva[i + 5*stride]), ia5);
	}
}
