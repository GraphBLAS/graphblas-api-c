#include<stdlib.h>
#include<stdint.h>
#include<assert.h>

typedef enum
{
    GrB_SUCCESS
} GrB_Info;

typedef uint64_t	GrB_Index;

const int32_t	DIM = 32;

class Matrix
{
    private:

    public:

	virtual double operator()(GrB_Index i, GrB_Index j) const = 0;
};

class ConcreteMatrix : public Matrix
{
    private:
	double	val[DIM][DIM];

    public:
	ConcreteMatrix()
	{
	    // The critical section is necessary because "random()" is not thread safe!
	    // There is a random_r version, but the interface seems to be a pain :-)
#pragma omp critical
	    {
		for (int i=0; i<DIM; i++) for (int j=0; j<DIM; j++)
		    val[i][j] = random() % DIM;
	    }
	}
	virtual double operator()(GrB_Index i, GrB_Index j) const { assert(i<DIM); assert(j<DIM); return val[i][j]; }
};

class MatrixTimesMatrix : public Matrix
{
    private:
	std::shared_ptr<Matrix>	_A;
	std::shared_ptr<Matrix>	_B;

    public:
	MatrixTimesMatrix(std::shared_ptr<Matrix> A, std::shared_ptr<Matrix> B) : _A(A), _B(B) { }
	virtual double operator()(GrB_Index i, GrB_Index j) const;
};

double MatrixTimesMatrix::operator()
(
    GrB_Index i, 
    GrB_Index j
) const
{
    double S = 0;
    for (int k=0; k<DIM; k++) S += (*_A)(i,k) * (*_B)(k,j);
    return S;
}

typedef std::shared_ptr<Matrix>	GrB_Matrix;

GrB_Info GrB_Matrix_new(
    GrB_Matrix& A
)
{
    A = std::shared_ptr<Matrix>((new ConcreteMatrix()));
    return GrB_SUCCESS;
}

GrB_Info GrB_mxm(
    GrB_Matrix&	C,
    GrB_Matrix	A,
    GrB_Matrix 	B
)
{	
    std::shared_ptr<Matrix> CC(new MatrixTimesMatrix(A,B));
    C = CC;
    return GrB_SUCCESS;
}
