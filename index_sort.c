#include <stdlib.h>

static double *x;

static int xcomp(const void *ap, const void *bp)
{
	int a = * (int *) ap - 1;
	int b = * (int *) bp - 1;
	
	return (x[a] > x[b]) - (x[b] > x[a]);
}

void index_sort(size_t N, int *idx, double *data)
{
	x = data;
	qsort(idx, N, sizeof(int), &xcomp);
}
