/*
 * @file spmat_array.h
 * @purpose Sparse matrix implemented using arrays
 */
#ifndef __SPMAT_ARRAY_H__
#define __SPMAT_ARRAY_H__

/* Includes ******************************************************************/
#include "matrix.h"
#include "submatrix.h"


/* Functions Declarations ****************************************************/
/* Allocates a new arrays sparse matrix of size n with nnz non-zero elements */
result_t
SPMAT_ARRAY_allocate(int n, int nnz, matrix_t **mat_out);

double
SPMAT_ARRAY_matrix_vector_sandwich(const matrix_t *mat, const double *v);

double
SUBMAT_SPMAT_ARRAY_get_1norm(const submatrix_t *submatrix);

double
SUBMAT_SPMAT_ARRAY_mult_vmv(const submatrix_t *submatrix, const double *vector);



#endif /* __SPMAT_ARRAY_H__ */

