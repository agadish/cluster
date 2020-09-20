/**
 * @file matrix_raw.c
 * @purpose 
 */

/* Includes **************************************************************************************/
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

#include "matrix.h"
#include "common.h"
#include "vector.h"
#include "debug.h"

/* Macros ****************************************************************************************/
#define RAW_ARRAY(matrix) ((double *)((matrix)->private))
#define MATRIX_AT(m, i, j) (RAW_ARRAY(m)[((i) * (m)->n) + (j)])

/* Functions Declarations ************************************************************************/
/**
 * @purpose Free a matrix which was previously created by matrix_raw_open or MATRIX_RAW_allocate
 *
 * @param matrix The matrix to free
 *
 * @remark Safe to call with NULL
 */
static
void
matrix_raw_free(matrix_t * matrix);

static
result_t
matrix_raw_add_row(matrix_t *mat, const double *row, int i);

/**
 * @purpose Given a matrix A and a column-vector b, Computes Ab
 * @param matrix input Matirx A [nxn]
 * @param vector input Vector b [nx1]
 * @param multiplication output Result vector [nx1]
 *
 * @return One of result_t values
 */
static
void
matrix_raw_mat_vector_multiply(const matrix_t *matrix,
                               const double *vector,
                               double * multiplication);

static
double
matrix_raw_get_1norm(const matrix_t *matrix);

/**
 * @purpose Given a matrix A and a line-vector b, Computes bA
 * @param matrix input Matirx A [nxn]
 * @param vector input Vector b [nx1]
 * @param multiplication output Result vector [nx1]
 *
 * @return One of result_t values
 */
#if 0
static
void
matrix_raw_vector_mat_multiply(const matrix_t *matrix,
                               const double *vector,
                               double * multiplication);
#endif /* 0 */

static
result_t
matrix_raw_split(matrix_t *matrix,
                  const double *vector_s,
                  int *temp_s_indexes,
                  matrix_t **matrix1_out,
                  matrix_t **matrix2_out);


/* Virtual Table *************************************************************/
const matrix_vtable_t MATRIX_RAW_VTABLE = {
    .add_row = matrix_raw_add_row,
    .free = matrix_raw_free,
    .mult = matrix_raw_mat_vector_multiply,
    .mult_vmv = NULL,
    .get_1norm = matrix_raw_get_1norm,
    .split = matrix_raw_split,
};


/* Functions *****************************************************************/
result_t
MATRIX_RAW_allocate(int n, matrix_t ** matrix_out)
{
    result_t result = E__UNKNOWN;
    matrix_t * matrix = NULL;
    size_t array_length = 0;

    /* 0. Input validation */
    if (NULL == matrix_out) {
        result = E__NULL_ARGUMENT;
        goto l_cleanup;
    }

    if (0 > n) {
        result = E__INVALID_SIZE;
        goto l_cleanup;
    }

    /* 1. Allocate matrix_t */
    matrix = (matrix_t *)malloc(sizeof(*matrix));
    if (NULL == matrix) {
        result = E__MALLOC_ERROR;
        goto l_cleanup;
    }
    /* 2. Initialize */
    (void)memset(matrix, 0, sizeof(*matrix));

    /* 3. Allocate data */
    array_length = sizeof(double) * n * n;
    matrix->vtable = &MATRIX_RAW_VTABLE;
    matrix->n = n;
    matrix->type = MATRIX_TYPE_RAW;
    /* DEBUG_PRINT("allocating matrix n=%d, size=%lu: addr %p to %p", */
    /*         n, */
    /*         array_length, */
    /*         matrix->private, */
    /*         (void *)((uint8_t* )matrix->private + array_length)); */
    matrix->private = malloc(array_length);
    if (NULL == matrix->private) {
        result = E__MALLOC_ERROR;
        goto l_cleanup;
    }

    /* Success */
    *matrix_out = matrix;

    result = E__SUCCESS;

l_cleanup:
    if (E__SUCCESS != result) {
        matrix_raw_free(matrix);
    }

    return result;
}

static
void
matrix_raw_free(matrix_t * matrix)
{
    if (NULL != matrix) {
        FREE_SAFE(matrix->private);
        FREE_SAFE(matrix);
    }
}


static
void
matrix_raw_mat_vector_multiply(const matrix_t *matrix,
                               const double *vector,
                               double * multiplication)
{
    double *current_row = NULL;
    int row = 0;

    for (row = 0 ; row < matrix->n ; ++row) {
        current_row = &MATRIX_AT(matrix, row, 0);
        multiplication[row] = VECTOR_scalar_multiply(current_row, vector, matrix->n);
    }
}

#if 0
static
void
matrix_raw_vector_mat_multiply(const matrix_t *matrix,
                               const double *vector,
                               double * multiplication)
{
    double column_sum = 0.0;
    int row = 0;
    int column = 0;

    /* TODO: Optimize */
    for (column = 0 ; column < matrix->n ; ++column) {
        column_sum = 0.0;
        for (row = 0 ; row < matrix->n ; ++row) {
            column_sum += vector[column] * MATRIX_AT(matrix, row, column);
        }
        multiplication[column] = column_sum;
    }
}
#endif /* 0 */

static
result_t
matrix_raw_add_row(matrix_t *mat, const double *row, int i)
{
    result_t result = E__UNKNOWN;

    double *destination_row = NULL;

    if ((NULL == mat) || (NULL == row)) {
        result = E__NULL_ARGUMENT;
        goto l_cleanup;

    }

    if (!MATRIX_IS_VALID_ROW_INDEX(mat, i)) {
        result = E__INVALID_ROW_INDEX;
        goto l_cleanup;
    }

    destination_row = &MATRIX_AT(mat, i, 0);
    (void)memcpy(destination_row, row, sizeof(*row) * mat->n);

    result = E__SUCCESS;
l_cleanup:

    return result;
} 

static
double
matrix_raw_get_1norm(const matrix_t *matrix)
{
    double max_col_sum = 0.0;
    double current_col_sum = 0.0;
    double *row = 0;
    double *last_row = 0;
    double *col = 0;
    double *last_col = 0;
    
    /* Go over each column */
    last_col = &MATRIX_AT(matrix, matrix->n, 0);
    for (col = &MATRIX_AT(matrix, 0, 0) ; last_col != col ; col += matrix->n) {
        /* Sum up all the rows within this column */
        last_row = col + matrix->n;
        for (row = col ; row < last_row ; ++row) {
            current_col_sum += fabs(*row);
        }

        /* Update max col sum */
        max_col_sum = MAX(max_col_sum, current_col_sum);
    }

    return max_col_sum;
}
static
result_t
matrix_raw_split(matrix_t *matrix,
                  const double * vector_s,
                  int * temp_s_indexes,
                  matrix_t **matrix1_out,
                  matrix_t **matrix2_out)
{
    result_t result = E__UNKNOWN;
    int matrix1_n = 0;
    int matrix2_n = 0;
    matrix_t *matrix1 = NULL;
    matrix_t *matrix2 = NULL;

    /* 0. Input validation */
    if ((NULL == matrix) || (NULL == vector_s) || (NULL == matrix1_out) ||
        (NULL == matrix2_out))
    {
        result = E__NULL_ARGUMENT;
        goto l_cleanup;
    }

    /* 1. Calculate s-index */
    matrix1_n = VECTOR_create_s_indexes(vector_s,
                                     matrix->n,
                                     temp_s_indexes);

    /* 2. Create matrixes */
    /* 2.1. matrix 1 */
    result = MATRIX_create_matrix(matrix1_n, MATRIX_TYPE_RAW, &matrix1);
    if (E__SUCCESS != result) {
        goto l_cleanup;
    }

    /* 2.2. matrix 2 */
    matrix2_n = matrix->n - matrix1_n;
    result = MATRIX_create_matrix(matrix2_n, MATRIX_TYPE_RAW, &matrix2);
    if (E__SUCCESS != result) {
        goto l_cleanup;
    }

    result = E__SUCCESS;
l_cleanup:
    if (E__SUCCESS != result) {
    }

    return result;
}
