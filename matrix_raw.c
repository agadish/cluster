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
#include "list_node.h"

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
void
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

/**
 * @purpose Given a matrix A and a line-vector b, Computes bA
 * @param matrix input Matirx A [nxn]
 * @param vector input Vector b [nx1]
 * @param multiplication output Result vector [nx1]
 *
 * @return One of result_t values
 */
static
void
matrix_raw_vector_mat_multiply(const matrix_t *matrix,
                               const double *vector,
                               double * multiplication);



/* Functions ************************************************************************************/
result_t
MATRIX_RAW_allocate(int n, matrix_t ** matrix_out)
{
    result_t result = E__UNKNOWN;
    matrix_t * matrix = NULL;

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
    (void)memset(&matrix, 0, sizeof(matrix));
    matrix->n = n;
    matrix->add_row = matrix_raw_add_row;
    matrix->free = matrix_raw_free;
    matrix->mult = matrix_raw_mat_vector_multiply;
    matrix->rmult = matrix_raw_vector_mat_multiply;

    /* 3. Allocate data */
    matrix->private = malloc(sizeof(double) * n * n);
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


static
void
matrix_raw_add_row(matrix_t *mat, const double *row, int i)
{
    result_t result = E__UNKNOWN;

    double *destination_row = NULL;

    if ((NULL == mat) || (NULL == row)) {
        result = E__NULL_ARGUMENT;
        goto l_cleanup;

    }

    if (!SPMAT_IS_VALID_ROW_INDEX(mat, i)) {
        result = E__INVALID_ROW_INDEX;
        goto l_cleanup;
    }

    destination_row = &MATRIX_AT(mat, i, 0);
    (void)memcpy(destination_row, row, sizeof(*row) * mat->n);

    result = E__SUCCESS;
l_cleanup:
    (void)result;

    return;
} 

