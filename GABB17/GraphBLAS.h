#ifndef GRAPHBLAS_H
#define GRAPHBLAS_H

#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>

typedef uint32_t GrB_Info;
typedef uint32_t GrB_Type;
typedef uint32_t GrB_Vector;
typedef uint32_t GrB_Matrix;
typedef uint32_t GrB_Descriptor;
typedef uint32_t GrB_Index;
typedef uint32_t GrB_BinaryOp;
typedef uint32_t GrB_UnaryOp;
typedef uint32_t GrB_Monoid;
typedef uint32_t GrB_Semiring;
typedef uint32_t GrB_Field;
typedef uint32_t GrB_Value;

const GrB_BinaryOp GrB_PLUS_FP32 = 0;
const GrB_BinaryOp GrB_TIMES_FP32 = 0;
const GrB_BinaryOp GrB_MINV_FP32 = 0;
const GrB_BinaryOp GrB_IDENTITY_BOOL = 0;
const GrB_BinaryOp GrB_TIMES_INT32 = 0;
const GrB_BinaryOp GrB_PLUS_INT32 = 0;

const GrB_Index *GrB_ALL = 0;

const GrB_Type GrB_FP32 = 0;
const GrB_Type GrB_INT32 = 0;
const GrB_Type GrB_BOOL = 0;
const uint32_t GrB_NULL = 0;

enum
{
  GrB_SUCCESS,
  GrB_REPLACE,
  GrB_OUTP,
  GrB_MASK,
  GrB_INP0,
  GrB_TRAN,
  GrB_SCMP
};

GrB_Info GrB_Descriptor_new(GrB_Descriptor*);
GrB_Info GrB_Descriptor_set(GrB_Descriptor, const GrB_Field, const GrB_Value);
GrB_Info GrB_Monoid_INT32_new(GrB_Monoid*, const GrB_Type, const GrB_BinaryOp, int32_t);
GrB_Info GrB_Monoid_FP32_new (GrB_Monoid*, const GrB_Type, const GrB_BinaryOp, float);
GrB_Info GrB_Semiring_new(GrB_Semiring*, const GrB_Monoid, const GrB_BinaryOp);

GrB_Info GrB_Matrix_build_INT32(GrB_Matrix*, const GrB_Matrix, const GrB_BinaryOp, const GrB_Index*, const GrB_Index*, int32_t*, const GrB_Index, const GrB_BinaryOp, const GrB_Descriptor);
GrB_Info GrB_Vector_new(GrB_Vector*, const GrB_Type, const GrB_Index n);
GrB_Info GrB_Matrix_new(GrB_Matrix*, const GrB_Type, const GrB_Index m, const GrB_Index n); 
GrB_Info GrB_Matrix_nrows(GrB_Index*, const GrB_Matrix);
GrB_Info GrB_Matrix_nvals(GrB_Index*, const GrB_Matrix);
GrB_Info GrB_Matrix_free(GrB_Matrix);
GrB_Info GrB_Matrix_Monoid_eWiseAdd (GrB_Matrix*, const GrB_Matrix, const GrB_BinaryOp, const GrB_Monoid  , const GrB_Matrix, const GrB_Matrix, const GrB_Descriptor);
GrB_Info GrB_Matrix_Monoid_eWiseMult(GrB_Matrix*, const GrB_Matrix, const GrB_BinaryOp, const GrB_Monoid  , const GrB_Matrix, const GrB_Matrix, const GrB_Descriptor);
GrB_Info GrB_Matrix_BinaryOp_reduce (GrB_Vector*, const GrB_Vector, const GrB_BinaryOp, const GrB_BinaryOp, const GrB_Matrix                  , const GrB_Descriptor);
GrB_Info GrB_mxm(GrB_Matrix*, const GrB_Matrix, const GrB_BinaryOp, const GrB_Semiring, const GrB_Matrix, const GrB_Matrix, const GrB_Descriptor);
GrB_Info GrB_Matrix_apply(GrB_Matrix*, const GrB_Matrix, const GrB_BinaryOp, const GrB_UnaryOp, const GrB_Matrix, const GrB_Descriptor);
GrB_Info GrB_Matrix_extract(GrB_Matrix*, const GrB_Matrix, const GrB_BinaryOp, const GrB_Matrix, const GrB_Index*, const GrB_Index, const GrB_Index*, const GrB_Index, const GrB_Descriptor);
GrB_Info GrB_Matrix_assign_constant_FP32(GrB_Matrix*, const GrB_Matrix, const GrB_BinaryOp, float, const GrB_Index*, const GrB_Index, const GrB_Index*, const GrB_Index, const GrB_Descriptor);

#define GrB_Matrix_build(C,Mask,accum,rowIDs,colIDs,values,n,dup,desc) _Generic((values), int32_t* : GrB_Matrix_build_INT32)(C,Mask,accum,rowIDs,colIDs,values,n,dup,desc)
#define GrB_free(obj) _Generic((obj), GrB_Matrix : GrB_Matrix_free)(obj)
#define GrB_free_all(...)
#define GrB_reduce(u,mask,accum,op,A,desc) _Generic((A), GrB_Matrix : _Generic((op), GrB_BinaryOp : GrB_Matrix_BinaryOp_reduce))(u,mask,accum,op,A,desc)
#define GrB_eWiseMult(C,Mask,accum,op,A,B,desc) _Generic((A), GrB_Matrix : _Generic((op), GrB_Monoid : GrB_Matrix_Monoid_eWiseMult))(C,Mask,accum,op,A,B,desc)
#define GrB_eWiseAdd(C,Mask,accum,op,A,B,desc) _Generic((A), GrB_Matrix : _Generic((op), GrB_Monoid : GrB_Matrix_Monoid_eWiseAdd ))(C,Mask,accum,op,A,B,desc)
#define GrB_Monoid_new(monoid,d1,op,identity) _Generic((identity), int32_t : GrB_Monoid_INT32_new, float : GrB_Monoid_FP32_new)(monoid,d1,op,identity)
#define GrB_apply(C,Mask,accum,op,A,desc) _Generic((A), GrB_Matrix : GrB_Matrix_apply)(C,Mask,accum,op,A,desc)
#define GrB_extract(C,Mask,accum,A,rows,m,cols,n,desc) _Generic((A), GrB_Matrix : GrB_Matrix_extract)(C,Mask,accum,A,rows,m,cols,n,desc)
#define GrB_assign(C,Mask,accum,val,rows,m,cols,n,desc) _Generic((val), float : GrB_Matrix_assign_constant_FP32)(C,Mask,accum,val,rows,m,cols,n,desc)

#endif // GRAPHBLAS_H
