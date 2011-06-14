#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <mkl_spblas.h>

#include <cholmod.h>

static cholmod_common c1;
static cholmod_sparse *mA, *mB;


double cholmod_logdet(double *Gb, int *Gbia, int *Gbja, int Natom)
{
    static cholmod_common c;
    static cholmod_sparse *A, *Ls;
    static cholmod_factor *L ;
    
    int     i, N, nnz, nnzb, s;
    int     *Ai, *Ap, *Lp;
    double    *Ax, *Lx;

	int	ldabsr, mblk, info;
	int job[6] = { 1, 0, 1, 0, 0, 1 };
	char matdescra[6]="GxxFxx";

    double logDetA;

	FILE	*f;

    cholmod_start (&c);
    c.final_asis = 0 ;
    c.final_super = 0 ;
    c.final_ll = 1 ;
    c.final_pack = 1 ;
    c.final_monotonic = 1 ;
	c.nmethods = 1 ;
	c.method [0].ordering = CHOLMOD_NATURAL ;
	c.postorder = 0 ;

	N = 3*Natom;
	nnzb = Gbia[Natom];
	nnz = 9*nnzb;
    A = cholmod_allocate_sparse(N, N, nnz, 1, 1, -1, CHOLMOD_REAL, &c);

	mblk = 3;
	ldabsr = 9;
	mkl_dcsrbsr(job, &Natom, &mblk, &ldabsr, A->x, A->i, A->p, Gb, Gbja, Gbia, &info);

    L = cholmod_analyze (A, &c) ;
    cholmod_factorize (A, L, &c) ;

	if (c.status == CHOLMOD_NOT_POSDEF) {
		f = fopen("G_cholmod.mtx", "w+");
		cholmod_write_sparse(f, A, NULL, NULL, &c);
		fclose(f);
		exit(EXIT_FAILURE);
		return 0.0;
	}

    Ls = cholmod_factor_to_sparse(L, &c);

    cholmod_band_inplace(0,0,1, Ls, &c);

    
    Lx = Ls->x;
    Lp = Ls->p;
    logDetA = 0.0;
    s = 1;
    for (i=0; i<3*Natom; i++) {
        if (Lx[Lp[i]] < 0.0) {
            s = -s;
        }
        logDetA = logDetA + 2.0*log(fabs(Lx[Lp[i]]));
    }


    if (s<0) {
        fprintf(stderr, "det(G) < 0!\n");
        exit(EXIT_FAILURE);
    }

    cholmod_free_factor (&L, &c) ;
    cholmod_free_sparse (&A, &c) ;
    cholmod_free_sparse (&Ls, &c) ;

    cholmod_finish (&c) ;

    return logDetA;
}

/*
void init_cholmod(int N, int nnz, void **Ax, void **Ai, void **Ap, void **Bx, void **Bi, void **Bp)
{
    cholmod_start (&mc);

    mA = cholmod_allocate_sparse(N, N, nnz, 1, 1, -1, CHOLMOD_REAL, &mc);
    mB = cholmod_allocate_sparse(N, N, nnz, 1, 1, -1, CHOLMOD_REAL, &mc);

    *Ax = A->x;
    *Ai = A->i;
    *Ap = A->p;
    *Bx = B->x;
    *Bi = B->i;
    *Bp = B->p;
}

void destroy_cholmod()
{
	cholmod_free_sparse(&mA, &mc);
	cholmod_free_sparse(&mB, &mc);

	cholmod_finish(&mc);
}
*/
