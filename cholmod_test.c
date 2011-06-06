#include <stdio.h>
#include <cholmod.h>

int main (void)
{
	FILE *f;
    cholmod_sparse *A, *Ls ;
    cholmod_dense *x, *b, *r ;
    cholmod_factor *L ;
    double one [2] = {1,0}, m1 [2] = {-1,0} ;	    /* basic scalars */
    cholmod_common c ;
	double *Ax, *Lx, detA;
	int	*Ai, *Ap, *Lp;
	
	size_t	i, N = 8, nnz=15;
		//int ia[ 9] = { 1, 5, 8, 10, 12, 15, 17, 18, 19 };
	int ia[9]={0, 2, 4, 6, 8, 10, 12, 14, 15	};
	/*int ja[18] = { 1, 3, 6, 7,
		2, 3, 5,
		3, 8,
		4, 7,
		5, 6, 7,
		6, 8,
			7,
		8 };*/
	int ja[15] = {0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7};
	/*
	double a[18] = { 7.0, 1.0, 2.0, 7.0,
		-4.0, 8.0, 2.0,
		1.0, 5.0,
		7.0, 9.0,
		5.0, 1.0, 5.0,
		-1.0, 5.0,
		11.0,
		5.0 };
	 */
	
	double a[15] = {
		4.0/3.0,
        -2.0/3.0,
		2.13333333333333,
		-1.2,
		3.08571428571429,
		-1.71428571428571,
		4.06349206349206,
		-20.0/9.0,
		5.05050505050505,
		-2.72727272727273,
		6.04195804195804,
		-3.23076923076922,
		7.03589743589742,
		-3.73333333333333,
		4.26666666666666
	};
	
    cholmod_start (&c) ;
	/* convert to packed LL' when done */
    c.final_asis = 0 ;
    c.final_super = 0 ;
    c.final_ll = 1 ;
    c.final_pack = 1 ;
    c.final_monotonic = 1 ;
	
	A = cholmod_allocate_sparse(N, N, nnz, 1, 1, -1, CHOLMOD_REAL, &c);
	
	Ai = A->i;
	Ap = A->p;
	Ax = A->x;
	
	for (i=0; i<=N; i++) {
		Ap[i] = ia[i];
	}
	for (i=0; i<nnz; i++) {
		Ai[i] = ja[i];
		Ax[i] = a[i];
	}
	
	f=fopen("L.mtx", "w+");
		//cholmod_write_sparse(f, A, NULL, NULL, &c);
	fclose(f);
	
    cholmod_print_sparse (A, "A", &c) ;		    /* print the matrix */
    if (A == NULL || A->stype == 0)		    /* A must be symmetric */
    {
		cholmod_free_sparse (&A, &c) ;
		cholmod_finish (&c) ;
		return (0) ;
    }
    b = cholmod_ones (A->nrow, 1, A->xtype, &c) ;   /* b = ones(n,1) */
    L = cholmod_analyze (A, &c) ;		    /* analyze */
    cholmod_factorize (A, L, &c) ;		    /* factorize */
	Ls = cholmod_factor_to_sparse(L, &c);
	cholmod_print_sparse (Ls, "L", &c);
	cholmod_write_sparse(stdout, Ls, NULL, NULL, &c);
	
	Lx = Ls->x;
	Lp = Ls->p;
	detA = 1.0;
	for (i=0; i<N; i++) {
		detA = detA*Lx[Lp[i]];
	}
	detA = detA * detA;
	printf("det(A) = %lg\n", detA);
    x = cholmod_solve (CHOLMOD_A, L, b, &c) ;	    /* solve Ax=b */
    r = cholmod_copy_dense (b, &c) ;		    /* r = b */
    cholmod_sdmult (A, 0, m1, one, x, r, &c) ;	    /* r = r-Ax */
    printf ("norm(b-Ax) %8.1e\n",
			cholmod_norm_dense (r, 0, &c)) ;	    /* print norm(r) */
	
    cholmod_free_factor (&L, &c) ;
    cholmod_free_sparse (&A, &c) ;
    cholmod_free_dense (&r, &c) ;
    cholmod_free_dense (&x, &c) ;
    cholmod_free_dense (&b, &c) ;
	
    cholmod_finish (&c) ;
    return (0) ;
}
