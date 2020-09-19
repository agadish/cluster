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
#include "common.h"


/* Functions Declarations ************************************************************************/
/* Allocates a new linked-lists sparse matrix of size n */
result_t
SPMAT_LIST_allocate(int n, bool_t should_init_rows_numbers, matrix_t **mat);

void
SPMAT_LIST_print(const char *matrix_name, matrix_t *mat_in);

result_t
SPMAT_LIST_write_neighbors(const matrix_t *matrix, FILE *file);


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
