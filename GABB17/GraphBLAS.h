#ifndef GRAPHBLAS_H
#define GRAPHBLAS_H

#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>

typedef uint32_t GrB_Info;
typedef uint32_t GrB_Type;
typedef struct { int x[1]; } *GrB_Vector;
typedef struct { int x[2]; } *GrB_Matrix;
typedef struct { int x[3]; } *GrB_Descriptor;
typedef uint32_t GrB_Index;
typedef struct { int x[4]; } *GrB_BinaryOp;
typedef struct { int x[5]; } *GrB_UnaryOp;
typedef struct { int x[6]; } *GrB_Monoid;
typedef struct { int x[7]; } *GrB_Semiring;
typedef uint32_t GrB_Field;
typedef uint32_t GrB_Value;

const GrB_UnaryOp GrB_MINV_FP32 = 0;
const GrB_UnaryOp GrB_IDENTITY_BOOL = 0;
const GrB_BinaryOp GrB_PLUS_FP32 = 0;
const GrB_BinaryOp GrB_TIMES_FP32 = 0;
const GrB_BinaryOp GrB_MAX_FP32 = 0;
const GrB_BinaryOp GrB_TIMES_INT32 = 0;
const GrB_BinaryOp GrB_PLUS_INT32 = 0;

const GrB_Index *GrB_ALL = 0;

const GrB_Type GrB_FP32 = 0;
const GrB_Type GrB_INT32 = 0;
const GrB_Type GrB_BOOL = 0;
// const uint32_t GrB_NULL = 0;
#define GrB_NULL 0

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
GrB_Info GrB_Vector_size(GrB_Index*, const GrB_Vector);
GrB_Info GrB_Vector_new(GrB_Vector*, const GrB_Type, const GrB_Index n);
GrB_Info GrB_Vector_dup(GrB_Vector*, const GrB_Vector);
GrB_Info GrB_Matrix_new(GrB_Matrix*, const GrB_Type, const GrB_Index m, const GrB_Index n); 
GrB_Info GrB_Matrix_nrows(GrB_Index*, const GrB_Matrix);
GrB_Info GrB_Matrix_nvals(GrB_Index*, const GrB_Matrix);
GrB_Info GrB_Matrix_free(GrB_Matrix);
GrB_Info GrB_Vector_free(GrB_Vector);
GrB_Info GrB_Monoid_free(GrB_Monoid);
GrB_Info GrB_Semiring_free(GrB_Semiring);
GrB_Info GrB_Matrix_Monoid_eWiseAdd (GrB_Matrix*, const GrB_Matrix, const GrB_BinaryOp, const GrB_Monoid  , const GrB_Matrix, const GrB_Matrix, const GrB_Descriptor);
GrB_Info GrB_Vector_Monoid_eWiseAdd (GrB_Vector*, const GrB_Vector, const GrB_BinaryOp, const GrB_Monoid  , const GrB_Vector, const GrB_Vector, const GrB_Descriptor);
GrB_Info GrB_Matrix_BinaryOp_eWiseAdd (GrB_Matrix*, const GrB_Matrix, const GrB_BinaryOp, const GrB_BinaryOp  , const GrB_Matrix, const GrB_Matrix, const GrB_Descriptor);
GrB_Info GrB_Vector_BinaryOp_eWiseAdd (GrB_Vector*, const GrB_Vector, const GrB_BinaryOp, const GrB_BinaryOp  , const GrB_Vector, const GrB_Vector, const GrB_Descriptor);
GrB_Info GrB_Matrix_Monoid_eWiseMult(GrB_Matrix*, const GrB_Matrix, const GrB_BinaryOp, const GrB_Monoid  , const GrB_Matrix, const GrB_Matrix, const GrB_Descriptor);
GrB_Info GrB_Matrix_BinaryOp_reduce (GrB_Vector*, const GrB_Vector, const GrB_BinaryOp, const GrB_BinaryOp, const GrB_Matrix                  , const GrB_Descriptor);
GrB_Info GrB_mxm(GrB_Matrix*, const GrB_Matrix, const GrB_BinaryOp, const GrB_Semiring, const GrB_Matrix, const GrB_Matrix, const GrB_Descriptor);
GrB_Info GrB_mxv(GrB_Vector*, const GrB_Vector, const GrB_BinaryOp, const GrB_Semiring, const GrB_Matrix, const GrB_Vector, const GrB_Descriptor);
GrB_Info GrB_Matrix_apply(GrB_Matrix*, const GrB_Matrix, const GrB_BinaryOp, const GrB_UnaryOp, const GrB_Matrix, const GrB_Descriptor);
GrB_Info GrB_Matrix_extract(GrB_Matrix*, const GrB_Matrix, const GrB_BinaryOp, const GrB_Matrix, const GrB_Index*, const GrB_Index, const GrB_Index*, const GrB_Index, const GrB_Descriptor);
GrB_Info GrB_Matrix_assign_constant_FP32(GrB_Matrix*, const GrB_Matrix, const GrB_BinaryOp, float, const GrB_Index*, const GrB_Index, const GrB_Index*, const GrB_Index, const GrB_Descriptor);
GrB_Info GrB_Vector_assign_constant_FP32(GrB_Vector*, const GrB_Vector, const GrB_BinaryOp, float, const GrB_Index*, const GrB_Index, const GrB_Descriptor);
GrB_Info GrB_Vector_assign_constant_FP64(GrB_Vector*, const GrB_Vector, const GrB_BinaryOp, double, const GrB_Index*, const GrB_Index, const GrB_Descriptor);

#define GrB_GET_7_OR_8_OR_9(_1,_2,_3,_4,_5,_6,_7,_8,_9,   NAME,...) NAME
#define GrB_Matrix_build(C,Mask,accum,rowIDs,colIDs,values,n,dup,desc) _Generic((values), int32_t* : GrB_Matrix_build_INT32)(C,Mask,accum,rowIDs,colIDs,values,n,dup,desc)
#define GrB_free(obj) _Generic((obj), GrB_Matrix : GrB_Matrix_free, GrB_Semiring : GrB_Semiring_free, GrB_Vector : GrB_Vector_free, GrB_Monoid : GrB_Monoid_free)(obj)
#define GrB_free_all(...)
#define GrB_reduce(u,mask,accum,op,A,desc) _Generic((A), GrB_Matrix : _Generic((op), GrB_BinaryOp : GrB_Matrix_BinaryOp_reduce))(u,mask,accum,op,A,desc)
#define GrB_eWiseMult(C,Mask,accum,op,A,B,desc) _Generic((A), GrB_Matrix : _Generic((op), GrB_Monoid : GrB_Matrix_Monoid_eWiseMult))(C,Mask,accum,op,A,B,desc)
#define GrB_eWiseAdd(C,Mask,accum,op,A,B,desc) _Generic((A), GrB_Matrix : _Generic((op), GrB_Monoid : GrB_Matrix_Monoid_eWiseAdd, GrB_BinaryOp : GrB_Matrix_BinaryOp_eWiseAdd ), GrB_Vector : _Generic((op), GrB_Monoid : GrB_Vector_Monoid_eWiseAdd, GrB_BinaryOp : GrB_Vector_BinaryOp_eWiseAdd))(C,Mask,accum,op,A,B,desc)
#define GrB_Monoid_new(monoid,d1,op,identity) _Generic((identity), int32_t : GrB_Monoid_INT32_new, float : GrB_Monoid_FP32_new)(monoid,d1,op,identity)
#define GrB_apply(C,Mask,accum,op,A,desc) _Generic((A), GrB_Matrix : GrB_Matrix_apply)(C,Mask,accum,op,A,desc)
#define GrB_extract(C,Mask,accum,A,rows,m,cols,n,desc) _Generic((A), GrB_Matrix : GrB_Matrix_extract)(C,Mask,accum,A,rows,m,cols,n,desc)
#define GrB_assign_7(u,mask,accum,val,indices,n,desc) _Generic((val), float : GrB_Vector_assign_constant_FP32, double : GrB_Vector_assign_constant_FP64)(u,mask,accum,val,indices,n,desc)
#define GrB_assign_9(C,Mask,accum,val,rows,m,cols,n,desc) _Generic((val), float : GrB_Matrix_assign_constant_FP32)(C,Mask,accum,val,rows,m,cols,n,desc)
#define GrB_assign(...) GrB_GET_7_OR_8_OR_9(__VA_ARGS__,GrB_assign_9,GrB_assign_8,GrB_assign_7)(__VA_ARGS__)

#endif // GRAPHBLAS_H
