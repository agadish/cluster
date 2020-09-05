/*
 * @file spmat_array.h
 * @purpose Sparse matrix implemented using arrays
 */
#ifndef __SPMAT_ARRAY_H__
#define __SPMAT_ARRAY_H__

/* Includes **************************************************************************************/
#include "matrix.h"

/* Functions Declarations ************************************************************************/
/* Allocates a new arrays sparse matrix of size n with nnz non-zero elements */
result_t
SPMAT_ARRAY_allocate(int n, int nnz, matrix_t **mat_out);


#endif /* __SPMAT_ARRAY_H__ */

