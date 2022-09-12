# A more concrete query proposal

## New Get and Set Functions

### **`GrB_<TYPE>_Option_Set(field::GrB_Option_Field, ...)`**

- `GrB_Global_Option_Set`
  
  - Non-Optional Reserved Fields:
  
  - Optional Reserved Fields:
    
    - `Print_1Based`::Bool
    
    - `DefaultStorageOrder_Hint`::`int32_t` $\in$ {`GrB_BYROW`, `GrB_BYCOL`}
    
    - `nthreads`::Integer

- `GrB_<Matrix | Vector | Scalar>_Option_Set`
  
  - Optional Reserved Fields:
    
    - `Name`::Null Terminated 128-Byte `char*`
    
    - `StorageOrder_Hint`::`int32_t` $\in$ {`GrB_BYROW`, `GrB_BYCOL`}

- `GrB_<OPERATOR>_Option_Set`
  
  - Non-Optional Reserved Fields:
    
    - `Name`::Null Terminated 128-Byte `char*`

- `GrB_Type_Option_Set`
  
  - Non-Optional Reserved Fields:
    
    - `Name`::Null Terminated 128-Byte `char*`

### **`GrB_<TYPE>_Option_Get(field::GrB_Option_Field, ...)`**

- `GrB_Global_Option_get`
  
  - - `Library_Name`::Null Terminated 128-Byte `char*`
    
    - `Library_Version`::{3 Integer Pointers | length 3 integer array}
    
    - `API_Version`::{3 Integers Pointers | length 3 integer array}
    
    - `Blocking_Mode`::`int32_t*`$\in$ {`GrB_BLOCKING`, `GrB_NONBLOCKING`}

- `GrB_<Matrix | Vector | Scalar>_Option_Get`
  
  - Non-Optional Reserved Fields:
    
    - `StorageOrder`::`int32_t*` $\in$ {`GrB_BYROW`, `GrB_BYCOL`}
    
    - `EltypeString`::Null-terminated 128-Byte `char*`
      
      - GraphBLAS Type (`GrB_Int32`) or `Name` for UDTs
      
      - This could be replaced by querying `Name` of `GrB_Type` returned by  `GrB_Matrix_Eltype(A)`. However Tim Davis would prefer this *as well*.
    
    - `Size`::`int64_t*`
      
      - Total size including all direct internal structures and header
  
  - Optional Reserved Fields:
    
    - `Name`::Null Terminated 128-Byte `char*`
    
    - `StorageOrder_Hint`::`int32_t*` $\in$ {`GrB_BYROW`, `GrB_BYCOL`}
    
    - `Format`::`int32_t*` $\in$ format enum
      
      - Optional since it *is* opacity breaking. However I believe this will be common enough, for performance debugging and more, to include as an *optional* field.

- `GrB_<OPERATOR>_Option_Get`
  
  - Non-Optional Reserved Fields:
    
    - `IsBuiltin`::`Bool*`
    
    - `DomainString`::Null-terminated `char*` of TBD length and format
      
      - Or 1/2/3 `char*`s depending on the operator (1 for monoid, 3 for binop/semiring, 2 for unaryop)

- `GrB_Type_Option_Get`
  
  - Non-Optional Reserved Fields:
    
    - `IsBuiltin`::`Bool*`
    
    - `Name`::Null Terminated 128-Byte `char*`
    
    - `Size`::`int64_t*`

**TODO: Descriptor**

## New Query Functions (not covered by get / set)

- `GrB_<Matrix | Vector | Scalar>_ElementType_(A, type::GrB_Type*)`

- `GrB_BinaryOp_Domain(op::BinaryOp, input1::GrB_Type*, input2::GrB_Type*, input3::GrB_Type*)`
  
  - And siblings

- Query:
  
  - Monoid and Multiplicative Op of Semiring
  
  - BinaryOp of a Monoid
  
  - Identity of a Monoid::`GrB_Scalar`

## String Handling

- Input
  
  - User passes a NULL-terminated string to the library. This string is then copied into internal buffers. Where a max length is specified the library must check that the max length is not exceeded.

- Output
  
  - For Fields with a known-size string the user must provide a `char*` of that size to be filled by the library
  
  - For Fields without a statically-size string (All of these will be non-reserved) the library must provide a `<Field>_StringSize` to be queried before querying the Field itself. 
