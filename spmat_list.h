/*
 * @file spmat_list.h
 * @purpose Sparse matrix implemented using linked lists
 */
#ifndef __SPMAT_LIST_H__
#define __SPMAT_LIST_H__

/* Includes ******************************************************************/
#include <stddef.h>
#include <stdio.h>

#include "matrix.h"
#include "submatrix.h"
#include "common.h"


/* Functions Declarations ****************************************************/
/* Allocates a new linked-lists sparse matrix of size n */
result_t
SPMAT_LIST_allocate(int n, matrix_t **mat);

void
SPMAT_LIST_print(const char *matrix_name, matrix_t *mat_in);

result_t
SPMAT_LIST_write_neighbors(const matrix_t *matrix, FILE *file);

result_t
SPMAT_LIST_transpose(const matrix_t *matrix, matrix_t **transposed_out);

double
SUBMAT_SPMAT_LIST_get_1norm(const submatrix_t *submatrix,
                            double *temp_n_sized_vector);

void
SUBMAT_SPMAT_LIST_mult(const submatrix_t *submatrix,
                       const double *vector,
                       double *result);

double
SUBMAT_SPMAT_LIST_calculate_q(const submatrix_t *submatrix,
                              const double *s_vector);

result_t
SUBMAT_SPMAT_LIST_split(submatrix_t *smat,
                        const double *s_vector,
                        submatrix_t *split1,
                        submatrix_t *split2);

#endif /* __SPMAT_LIST_H__ */
