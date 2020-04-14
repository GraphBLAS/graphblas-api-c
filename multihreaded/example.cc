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
    GrB_Matrix 	A[N], B[N], C[N], D[N], E[N], F[N], G[N];	// A, B, C, D, E are "arrays" of matrices
    double	S[N][DIM][DIM];					// Results for checking

#pragma omp parallel
    {
	if (0 == omp_get_thread_num())
	{
	    int p = omp_get_num_threads();
	    fprintf(stdout, "Running on %3d thread(s)\n", p);
	}

	// Create all input matrices.
#pragma omp for schedule(dynamic)
	for (int i=0; i<N; i++) 
	{
	    GrB_Matrix_new(A[i]);
	    GrB_Matrix_new(B[i]);
	    GrB_Matrix_new(D[i]);
	    GrB_Matrix_new(F[i]);
	    GrB_Matrix_new(G[i]);
	}

	// Compute the reference results for testing E = A*B*D;
#pragma omp for schedule(dynamic)
	for (int I=0; I<N; I++)
	{
	    for (int i=0; i<DIM; i++) for (int j=0; j<DIM; j++)
	    {
		S[I][i][j] = 0;
		for (int k=0; k<DIM; k++)
		{
		    double R = 0;
		    for (int l=0; l<DIM; l++)
			R += (*A[I])(i,l)*(*B[I])(l,k);
		    S[I][i][j] += R*(*D[I])(k,j);
		}
	    }
	}

#pragma omp for schedule(dynamic)
	for (int i=0; i<N; i++)
	    GrB_mxm(C[i], A[i], B[i]);

#pragma omp for schedule(dynamic)
	for (int i=0; i<N; i++)
	{
	    GrB_mxm(A[i], F[i], G[i]);
	    GrB_mxm(E[i], C[i], D[i]);
	}

#pragma omp for schedule(dynamic)
	for (int I=0; I<N; I++)
	{
	    for (int i=0; i<DIM; i++) for (int j=0; j<DIM; j++)
	    {
		if (S[I][i][j] != (*E[I])(i,j)) 
		{
		    printf("FAIL: S = %f, E[%d](%d,%d) = %f\n", S[I][i][j], I, i, j, (*E[I])(i,j));
		    exit(-1);
		}
	    }
	}

#pragma omp for schedule(dynamic)
	for (int I=0; I<N; I++)
	{
	    for (int i=0; i<DIM; i++) for (int j=0; j<DIM; j++)
	    {
		double S = 0;
		for (int k=0; k<DIM; k++)
		    S += (*F[I])(i,k) * (*G[I])(k,j);
		if (S != (*A[I])(i,j))
		{
		    printf("FAIL: S = %f, A[%d](%d,%d) = %f\n", S, I, i, j, (*A[I])(i,j));
		    exit(-1);
		}
	    }
	}
    }

    printf("PASS\n");

    return 0;
}
