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

/*
 * Calculate the 1-norm of a given submatrix
 *
 * @param submatrix The submatrix
 * @param tmp_rows_sums Temp buffer with size n
 */
double
SUBMAT_SPMAT_LIST_get_1norm(const submatrix_t *smat,
                            double *tmp_row_sums);

/*
 * Multiply the submatrix with a given vector, to a pre-allocated buffer
 *
 * @param submatrix The submatrix
 * @param vector Buffer to multiply with
 * @param result pre-allocated buffer
 */
void
SUBMAT_SPMAT_LIST_mult(const submatrix_t *submatrix,
                       const double *vector,
                       double *result);

/**
 * Calculate the Q of the submatrix with a given vector using the formula
 * learned in class
 */
double
SUBMAT_SPMAT_LIST_calculate_q(const submatrix_t *submatrix,
                              const double *s_vector);

/**
 * Split a submatrix into two submatrices accordingly to a given s-vector
 */
result_t
SUBMAT_SPMAT_LIST_split(submatrix_t *smat,
        const double * vector_s,
        int *temp_s_indexes,
        submatrix_t **matrix1_out,
        submatrix_t **matrix2_out);

/**
 * Calculate the the improved formula Q score within algorithm 4
 */
double
SUBMAT_SPMAT_LIST_calc_q_score(const submatrix_t *smat,
                       const double *vector,
                       int row);


#endif /* __SPMAT_LIST_H__ */
