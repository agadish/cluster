/**
 * @file matrix_raw.c
 * @purpose 
 */

/* Includes ******************************************************************/
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


/* Macros ********************************************************************/
#define GET_DATA(matrix) ((matrix_raw_data_t *)((matrix)->private))
#define GET_INDEXES(matrix) (GET_DATA((matrix))->indexes)
#define GET_ARRAY(matrix) (GET_DATA((matrix))->array)
#define MATRIX_AT(m, i, j) (GET_ARRAY(m)[((i) * (m)->n) + (j)])


/* Structs *******************************************************************/
/* spmat_row_t * matrix->private: n-sized array of spmat_row_t */
typedef struct matrix_raw_data_s {
    /* n*n array contains the data */
    double *array;

    /* n array contains the indexes */
    int *indexes;
} matrix_raw_data_t;


/* Functions Declarations ****************************************************/
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

static
result_t
matrix_raw_decrease_row_sum_from_diag(matrix_t *matrix);


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
matrix_raw_divide(matrix_t *matrix,
                  const double *vector_s,
                  int *temp_s_indexes,
                  matrix_t **matrix1_out,
                  matrix_t **matrix2_out);

static
double
matrix_raw_mult_vmv(const matrix_t *mat, const double *v);


static
result_t
matrix_raw_write_neighbors(const matrix_t *matrix, FILE *file);

static
void
matrix_raw_initialise_rows_numbers(matrix_t *mat);

static
void
matrix_raw_get_row(matrix_t *mat, double *row, int i);


/* Virtual Table *************************************************************/
const matrix_vtable_t MATRIX_RAW_VTABLE = {
    .add_row = matrix_raw_add_row,
    .get_row = matrix_raw_get_row,
    .free = matrix_raw_free,
    .mult = matrix_raw_mat_vector_multiply,
    .mult_vmv = matrix_raw_mult_vmv,
    .get_1norm = matrix_raw_get_1norm,
    .decrease_rows_sums_from_diag = matrix_raw_decrease_row_sum_from_diag,
    .divide = matrix_raw_divide,
    .write_neighbors = matrix_raw_write_neighbors,
    .initialise_row_numbers = matrix_raw_initialise_rows_numbers,
};


/* Functions *****************************************************************/
result_t
MATRIX_RAW_allocate(int n, matrix_t ** matrix_out)
{
    result_t result = E__UNKNOWN;
    matrix_t * matrix = NULL;
    size_t array_length = 0;
    matrix_raw_data_t *data = NULL;
    double *array = NULL;
    int *indexes = NULL;

    /* 0. Input validation */
    if (NULL == matrix_out) {
        result = E__NULL_ARGUMENT;
        goto l_cleanup;
    }

    if (0 > n) {
        result = E__INVALID_SIZE;
        goto l_cleanup;
    }

    /* 1. Allocate indexes */
    indexes = (int *)malloc(sizeof(*indexes) * n);
    if (NULL == indexes) {
        result = E__MALLOC_ERROR;
        goto l_cleanup;
    }
    DEBUG_PRINT("matrix->data->indexes is %p", (void *)indexes);

    /* 1. Allocate array */
    array_length = sizeof(double) * n * n;
    array = malloc(array_length);
    if (NULL == array) {
        result = E__MALLOC_ERROR;
        goto l_cleanup;
    }
    (void)memset(array, 0, array_length);
    DEBUG_PRINT("matrix->data->array is %p", (void *)array);

    /* 4. Allocate matrix */
    matrix = (matrix_t *)malloc(sizeof(*matrix));
    if (NULL == matrix) {
        result = E__MALLOC_ERROR;
        goto l_cleanup;
    }
    DEBUG_PRINT("matrix is %p", (void *)matrix);

    /* 3. Allocate data */
    data = (matrix_raw_data_t *)malloc(sizeof(*data));
    if (NULL == data) {
        result = E__MALLOC_ERROR;
        goto l_cleanup;
    }
    DEBUG_PRINT("matrix->data is %p", (void *)data);


    /* 4. Assign values */
    data->array = array;
    data->indexes = indexes;
    matrix->n = n;
    matrix->type = MATRIX_TYPE_RAW;
    matrix->vtable = &MATRIX_RAW_VTABLE;
    matrix->private = (void *)data;

    /* Success */
    *matrix_out = matrix;

    result = E__SUCCESS;

l_cleanup:
    if (E__SUCCESS != result) {
        FREE_SAFE(matrix);
        FREE_SAFE(data);
        FREE_SAFE(array);
        FREE_SAFE(indexes);
    }

    return result;
}

static
void
matrix_raw_free(matrix_t * matrix)
{
    if (NULL != matrix) {
        if (NULL != matrix->private) {
            FREE_SAFE(GET_INDEXES(matrix));
            FREE_SAFE(GET_ARRAY(matrix));
            FREE_SAFE(matrix->private);
        }
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
    int col = 0;
    int row = 0;
    
    /* Go over each column */
    for (col = 0 ; col < matrix->n ; ++col) {
        /* Sum up all the rows within this column */
        current_col_sum = 0.0;
        for (row = 0 ; row < matrix->n ; ++row) {
            current_col_sum += fabs(MATRIX_AT(matrix, row, col));
        }

        /* Update max col sum */
        max_col_sum = MAX(max_col_sum, current_col_sum);
    }

    return max_col_sum;
}

static
result_t
matrix_raw_decrease_row_sum_from_diag(matrix_t *matrix)
{
    result_t result = E__UNKNOWN;
    double current_row_sum = 0.0;
    int row = 0;
    int col = 0;
    
    /* 0. Input validation */
    if (NULL == matrix) {
        result = E__NULL_ARGUMENT;
        goto l_cleanup;
    }

    /* Go over each col */
    for (row = 0 ; row < matrix->n ; ++row) {
        /* Sum up the row */
        for (col = 0 ; col < matrix->n ; ++col) {
            current_row_sum += MATRIX_AT(matrix, row, col);
        }

        /* Decrease row sum from diag */
        MATRIX_AT(matrix, row, row) -= current_row_sum;
    }

    result = E__SUCCESS;
l_cleanup:
    
    return result;
}

static
result_t
matrix_raw_divide(matrix_t *matrix,
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

static
double
matrix_raw_mult_vmv(const matrix_t *mat, const double *v)
{
    int i = 0;
    double current_mult = 0.0;
    double result = 0.0;

    for (i = 0 ; i < mat->n ; ++i) {
        current_mult = VECTOR_scalar_multiply(&MATRIX_AT(mat, i, 0), v, mat->n);
        result += v[i] * current_mult;
    }

    return result;
}

static
result_t
matrix_raw_write_neighbors(const matrix_t *matrix, FILE *file)
{
    result_t result = E__UNKNOWN;
    int *current_index = NULL;
    int *last_index = NULL;
    size_t result_write = 0;

    if ((NULL == matrix) || (NULL == file)) {
        result = E__NULL_ARGUMENT;
        goto l_cleanup;
    }

    current_index = GET_INDEXES(matrix);
    last_index = GET_INDEXES(matrix) + matrix->n;

    for (; last_index > current_index ; ++current_index) {
        result_write = fwrite((void *)current_index,
                              sizeof(*current_index),
                              1,
                              file);
        if (1 != result_write) {
            result = E__FWRITE_ERROR;
            goto l_cleanup;
        }
    }
    
    result = E__SUCCESS;
l_cleanup:

    return result;
}

static
void
matrix_raw_initialise_rows_numbers(matrix_t *mat)
{
    int i = 0;

    for (i = 0 ; i < mat->n ; ++i) {
        GET_INDEXES(mat)[i] = i;
    }
}

static
void
matrix_raw_get_row(matrix_t *mat, double *row, int i)
{
    (void)memcpy(row, &MATRIX_AT(mat, i, 0), sizeof(*row) * mat->n);
}
