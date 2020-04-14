#include<stdint.h>
#include<memory>
#include<omp.h>
#include<graphblas.h>

const int32_t N = 32;

int main
(
    int		argc,
    char      **argv
)
{ 
    GrB_Matrix A[N], B[N], C[N], D[N], E[N];  // A, B, C, D, E are "arrays" of matrices

    // Create all input matrices.
    for (int i=0; i<N; i++) 
    {
	GrB_Matrix_new(A[i]);
	GrB_Matrix_new(B[i]);
	GrB_Matrix_new(D[i]);
    }

#pragma omp parallel
    {
	if (0 == omp_get_thread_num())
	{
	    int p = omp_get_num_threads();
	    fprintf(stdout, "Running on %3d thread(s)\n", p);
	}

#pragma omp for schedule(dynamic)
	for (int i=0; i<N; i++)
	    GrB_mxm(C[i], A[i], B[i]);

#pragma omp for schedule(dynamic)
	for (int i=0; i<N; i++)
	    GrB_mxm(E[i], C[i], D[i]);

    }

    for (int I=0; I<N; I++)
    {
	for (int i=0; i<DIM; i++) for (int j=0; j<DIM; j++)
	{
	    double S = 0;
	    for (int k=0; k<DIM; k++)
	    {
		double R = 0;
		for (int l=0; l<DIM; l++)
		    R += (*A[I])(i,l)*(*B[I])(l,k);
		S += R*(*D[I])(k,j);
	    }
	    if (S != (*E[I])(i,j)) 
	    {
		printf("FAIL: S = %f, E[%d](%d,%d) = %f\n", S, I, i, j, (*E[I])(i,j));
		exit(-1);
	    }
	}
    }
    printf("PASS\n");

    return 0;
}
