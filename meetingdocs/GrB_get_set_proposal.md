# Proposal for GrB_get and GrB_set

For the following objects:

- GrB_Scalar, GrB_Vector, GrB_Matrix
- GrB_UnaryOp, GrB_IndexUnaryOp, GrB_BinaryOp, GrB_Monoid, GrB_Semiring
- GrB_Descriptor, GrB_Type
- One unique Global “object”
  - The Global object will hold options like 

propose the creation of two methods:

```c
GrB_Info GrB_<TYPE>_set(GrB_<TYPE>	x, GrB_Field field, ...)
```

and

```c
GrB_Info GrB_<TYPE>_get(GrB_<TYPE> x, GrB_Field field, ...)
```

(The Global version will not take the argument `x`).

These functions will be variadic, the type of the variadic arguments determined by the `GrB_Field field` argument (and optionally `<TYPE>`).

The `GrB_Descriptor_set` function already exists, and will be subsumed by this spec if possible. If not, then `GrB_Desc_set` may be used.

### Return Values

The return values held by `GrB_Info` for this operation are:

- `GrB_SUCCESS`: An option has been successfully `set` or `get`.
- `GrB_INVALID_VALUE`: One or more parameters were not supported or were malformed. Particular objects will typically only support a subset of `GrB_Field` options.
- `GrB_OUT_OF_MEMORY`
- `GrB_UNINITIALIZED_OBJECT`
- `GrB_PANIC`

### GrB_Field

Proposed values for the GrB_Field enum, **REQUIRED**:

- [ ] `GrB_NAME::(NULL terminated 256 length char*)`: All objects except Global
- [ ] `GrB_IMPLEMENTATION_NAME::(NULL terminated 256 length char*)`: Gloabl only
- [ ] `GrB_SPEC_VERSION::(3 UInt64s)`: Global only
  - Not clear if this should be a length 3 array, or 3 arguments to variadic
- [ ] `GrB_IMPLEMENTATION_VERSION::(3 UInt64s)`: Global only
  - Same issue as directly above
- [ ] Existing `GrB_Desc_Field` values

`GrB_IMPLEMENTATION_NAME` and `GrB_IMPLEMENTATION_VERSION` ***must*** exist to enable disambiguating `GrB_Field` extensions

Proposed values for the GrB_Field enum, questionable or optional:

- [ ] `GrB_ORIENTATION::(GrB_Orientation enum, GrB_BYROW, GrB_BYCOL, GrB_DEFAULT)`: Global and arrays
  - Optional, allows users to specify row or column orientation if available. 
  - Important since some languages use row major, and others column major arrays. The overhead of conversion during import and export may be high if there is a mismatch.
  - Breaks opacity
- [ ] `GrB_NTHREADS::(UInt64)`: Global
  - Optional since a CUDA implementation would not want
- [ ] `GrB_NAME_SIZE::(UInt64)`: All objects except global
  - Potentially required, see below

### Undecided Issues

- What should be a get/set and what should have its own function?

  - I would argue that any computation inducing *or* `GrB` object returning function must be its own function. 
  - Anything returning a primitive or `char*` in an *immediate* fashion should be a get/set
  - [ ] Decision

- String handling

  - I am *strongly* in favor of zero ownership transfer (implicit and explicit), even if it makes a mismatch with the error handling
  - Reason: If `T` is a GrB_Type (or GrB_Matrix) with a name, and I `GrB_get(T, GrB_NAME, ...)` then:
    - If I have a reference to an internal GraphBLAS string which outlives `T` I may segfault. 
    - If instead I provide a buffer to be filled then there is *no* chance of segfaulting
  - I am also generally in favor of a size *and* a null terminator. This size may be runtime (`GrB_get(T, GrB_NAME_SIZE, ...)`) or static, but the null termination is valuable for a lot of subsequent string processing tools. It also enables the library to return `GrB_NAME_SIZE` as static, overestimating the necessary buffer.
  - Decision on ownership:
    - [ ] Expose internal string
    - [ ] Fill buffer
  - Decision on Null termination:
    - [ ] Null termination required
    - [ ] Non-terminated
  - Decision on `GrB_NAME_SIZE`:
    - [ ] Required
    - [ ] Static size of []

- When and where can type name be set:

  - Decision on number of times type name can be set:
    - [ ] Once (return `GrB_INVALID_VALUE`)
    - [ ] Multiple times
  - Decision on where it can be set:
    - [ ] Immediately after type creation
    - [ ] As part of `GrB_Type` ctor

  

  

