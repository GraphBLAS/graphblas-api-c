#include<stdlib.h>
#include<stdint.h>
#include<assert.h>

bool ForceError = false;

static const char *IndexOutOfBoundsMsg = "Index out of bounds";

class Error
{
};

class APIError : public Error
{
};

class ExecutionError : public Error
{
};

class IndexOutOfBounds : public ExecutionError
{
};

static uint32_t max
(
    uint32_t 	a,
    uint32_t	b
)
{
    return (a > b) ? a : b;
}

typedef enum
{
    GrB_SUCCESS,
    GrB_INDEX_OUT_OF_BOUNDS
} GrB_Info;

typedef uint64_t	GrB_Index;

const int32_t	DIM = 32;

class Matrix
{
    private:

	const char 	*_error = 0;

    protected:
	uint32_t 	 _depth;

    public:

	virtual double operator()(GrB_Index i, GrB_Index j) const = 0;
	const char *error() const { return _error; }
	const uint32_t depth() const { return _depth; }
	void setError(const char* error) { _error = error; }
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
	    _depth = 0;
	}

	ConcreteMatrix(const Matrix& A)
	{
	    for (int i=0; i<DIM; i++) for (int j=0; j<DIM; j++)
		val[i][j] = A(i,j);
	}

	virtual double operator()(GrB_Index i, GrB_Index j) const 
	{ 
	    assert(i<DIM); assert(j<DIM); 
	    if (ForceError) throw IndexOutOfBounds();
	    return val[i][j]; 
	}
};

class MatrixTimesMatrix : public Matrix
{
    private:
	std::shared_ptr<Matrix>	_A;
	std::shared_ptr<Matrix>	_B;

    public:
	MatrixTimesMatrix(std::shared_ptr<Matrix> A, std::shared_ptr<Matrix> B) : _A(A), _B(B) { _depth = max(A->depth(), B->depth()) + 1; }
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

GrB_Info GrB_Matrix_error
(
    const char		**error,
    const GrB_Matrix	  A
)
{
    *error = A->error();
    return GrB_SUCCESS;
}

GrB_Info GrB_Matrix_wait
(
    GrB_Matrix&	A
)
{
    try
    {
	std::shared_ptr<Matrix> AA(new ConcreteMatrix(*A));
	A = AA;
	return GrB_SUCCESS;
    }
    catch (IndexOutOfBounds e)
    {
	A->setError(IndexOutOfBoundsMsg);
	return GrB_INDEX_OUT_OF_BOUNDS;
    }
}
