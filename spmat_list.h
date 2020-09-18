/*
 * @file spmat_list.h
 * @purpose Sparse matrix implemented using linked lists
 */
#ifndef __SPMAT_LIST_H__
#define __SPMAT_LIST_H__

/* Includes **************************************************************************************/
#include <stddef.h>
#include <stdio.h>

#include "matrix.h"


/* Functions Declarations ************************************************************************/
/* Allocates a new linked-lists sparse matrix of size n */
result_t
SPMAT_LIST_allocate(int n, matrix_t **mat);

result_t
SPMAT_LIST_initialise_rows_numbers(matrix_t *mat);

result_t
SPMAT_LIST_divide_matrix(const matrix_t *matrix,
                         const double * vector_s,
                         matrix_t **matrix1_out,
                         matrix_t **matrix2_out);

void
SPMAT_LIST_print(const char *matrix_name, matrix_t *mat_in);

result_t
SPMAT_LIST_get_1norm(const matrix_t *matrix, double *norm_out);

result_t
SPMAT_LIST_decrease_rows_sums_from_diag(matrix_t *matrix);

result_t
SPMAT_LIST_write_neighbors(const matrix_t *matrix, FILE *file);

result_t
SPMAT_LIST_initialise_rows_numbers(matrix_t *mat);

/**
 * @purpose multiplying row vector, matrix, and same vector as col vector
 * @param mat input matrix
 * @param v input vector

 * @return result of multiplication
 *
 */
double
SPMAT_LIST_matrix_vector_sandwich(const matrix_t *mat, const double *v);


#endif /* __SPMAT_LIST_H__ */
