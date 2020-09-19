/*
 * @file spmat_array.c
 * @purpose Sparse matrix implemented using arrays
 */
/* Includes ******************************************************************/
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "results.h"
#include "matrix.h"
#include "spmat_array.h"
#include "common.h"
#include "submatrix.h"

/* Macros ********************************************************************/
#define GET_ARRAY(mat) ((spmat_array_t *)((mat)->private))

/* Structs *******************************************************************/
typedef struct spmat_array_s {
    /* 3 arrays representing a sparse matrix: non zero elements array, their
     * indexes array, and an array of the the first non zero element in every
     * row
     */
    double *values;
    int *colind;
    int *rowptr;
    int last_index;
} spmat_array_t;


/* Functions Declarations ****************************************************/
static
result_t
spmat_array_add_row(matrix_t *A, const double *row, int i);

static
void
spmat_array_free(matrix_t *A);

static
void
spmat_array_mult(const matrix_t *A, const double *v, double *result);


/* Global Variables **********************************************************/
const matrix_vtable_t SPMAT_ARRAY_VTABLE = {
    .add_row = spmat_array_add_row,
    .free = spmat_array_free,
    .mult = spmat_array_mult,
    .mult_vmv = NULL,
    .get_1norm = NULL,
    .decrease_rows_sums_from_diag = NULL
};


/* Functions *****************************************************************/
result_t
SPMAT_ARRAY_allocate(int n, int nnz, matrix_t **mat_out)
{
    result_t result = E__UNKNOWN;
    matrix_t * mat = NULL;
    double *values = NULL;
    int *colind = NULL;
    int *rowptr = NULL;
    int i = 0;
    spmat_array_t *spmat_array_data = NULL;

    mat = (matrix_t *)malloc(sizeof(*mat));
    if (NULL == mat) {
        result = E__MALLOC_ERROR;
        goto l_cleanup;
    }
    mat->private = NULL;
    mat->vtable = &SPMAT_ARRAY_VTABLE;
    mat->n = n;
    mat->type = MATRIX_TYPE_SPMAT_ARRAY;

    spmat_array_data = (spmat_array_t *)malloc(sizeof(spmat_array_t));
    if (NULL == spmat_array_data) {
        result = E__MALLOC_ERROR;
        goto l_cleanup;
    }
    mat->private = spmat_array_data;
    (void)memset(spmat_array_data, 0, sizeof(spmat_array_t));
    spmat_array_data->last_index = 0;

    values = (double*)malloc(sizeof(double)*nnz);
    if (NULL == values) {
        result = E__MALLOC_ERROR;
        goto l_cleanup;
    }
    spmat_array_data->values = values;
    (void)memset(values, 0, sizeof(double)*nnz);

    colind = (int*)malloc(sizeof(int)*nnz);
    if (NULL == colind) {
        result = E__MALLOC_ERROR;
        goto l_cleanup;
    }
    spmat_array_data->colind = colind;
    (void)memset(colind, 0, sizeof(int)*nnz);

    rowptr = ((int*)malloc(sizeof(int)*(n+1)));
    if (NULL == rowptr) {
        result = E__MALLOC_ERROR;
        goto l_cleanup;
    }
    spmat_array_data->rowptr = rowptr;
    for (i = 0; i<n+1; i++){
        rowptr[i] = nnz;
    }

    /* Success */
    result = E__SUCCESS;
    *mat_out = mat;

l_cleanup:
    if (E__SUCCESS != result) {
        spmat_array_free(mat);
        mat = NULL;
    }

    return result;
}

double
SPMAT_ARRAY_matrix_vector_sandwich(const matrix_t *mat, const double *v)
{
    int row;
    int n;
    int i;
    int first_index_row;
    int non_zero_in_row;
    double row_sum = 0.0;
    double result = 0.0;
    spmat_array_t *spmat_array_data;
    double *values;
    int *colind, *rowptr;


    spmat_array_data = (spmat_array_t *)mat->private;
    n = mat->n;
    values = spmat_array_data->values;
    colind = spmat_array_data->colind;
    rowptr = spmat_array_data->rowptr;

    for (row = 0 ; row < n; ++row) {
        first_index_row = rowptr[row];
        non_zero_in_row = rowptr[row + 1] - rowptr[row];
        if (0 != non_zero_in_row) {
            for (i = 0 ; i < non_zero_in_row ; ++i){
                row_sum += (values[first_index_row + i] * v[colind[first_index_row + i]]);
            }
        }

        result += row_sum * v[row];
        row_sum = 0.0;
    }

    return result;
}

static
void
spmat_array_free(matrix_t *mat)
{
    spmat_array_t *spmat_array_data = NULL;

    if (NULL != mat) {
        spmat_array_data = (spmat_array_t *)mat->private;
        if (NULL != spmat_array_data){
            FREE_SAFE(spmat_array_data->values);
            FREE_SAFE(spmat_array_data->colind);
            FREE_SAFE(spmat_array_data->rowptr);
        }
        FREE_SAFE(spmat_array_data);
        mat->private = NULL;
        FREE_SAFE(mat);
    }
}

static
result_t
spmat_array_add_row(matrix_t *mat, const double *data, int row_index)
{
    result_t result = E__UNKNOWN;
    bool_t first_found = FALSE;
    int j, n;
    int k;
    double number;
    spmat_array_t *spmat_array_data;
    double *values;
    int *colind, *rowptr;

    /* 0. Input validation */
    if ((NULL == mat) || (NULL == mat->private) || (NULL == data)) {
        result = E__NULL_ARGUMENT;
        goto l_cleanup;
    }

    if (!MATRIX_IS_VALID_ROW_INDEX(mat, row_index)) {
        result = E__INVALID_ROW_INDEX;
        goto l_cleanup;
    }

    spmat_array_data = (spmat_array_t *)mat->private;
    n = mat->n;
    values = spmat_array_data->values;
    colind = spmat_array_data->colind;
    rowptr = spmat_array_data->rowptr;

    for (j = 0; j < n; ++j) {
        number = data[j];
        if (0 != number) {
            if (FALSE == first_found) {
                first_found = TRUE;
                for (k = row_index;
                     ((k >= 0) && (rowptr[k] == rowptr[n])) ;
                     --k)
                {
                    rowptr[k] = spmat_array_data->last_index;
                }
            }
            values[spmat_array_data->last_index] = number;
            colind[spmat_array_data->last_index] = j;
            ++spmat_array_data->last_index;
        }
    }

    result = E__SUCCESS;
l_cleanup:

    return result;
}

static
void
spmat_array_mult(const matrix_t *mat, const double *v, double *result)
{
    int row;
    int n;
    int i;
    int first_index_row;
    int non_zero_in_row;
    double sum = 0.0;
    spmat_array_t *spmat_array_data;
    double *values;
    int *colind, *rowptr;


    spmat_array_data = (spmat_array_t *)mat->private;
    n = mat->n;
    values = spmat_array_data->values;
    colind = spmat_array_data->colind;
    rowptr = spmat_array_data->rowptr;

    for (row = 0 ; row < n; ++row) {
        first_index_row = rowptr[row];
        non_zero_in_row = rowptr[row + 1] - rowptr[row];
        if (0 != non_zero_in_row) {
            for (i = 0 ; i < non_zero_in_row ; ++i){
                sum += (values[first_index_row + i] *
                        v[colind[first_index_row + i]]);
            }
        }

        result[row] = sum;
        sum = 0.0;
    }
}

#if 0
double
SUBMAT_SPMAT_ARRAY_get_1norm(const submatrix_t *submatrix)
{
    const matrix_t *matrix = NULL;
    const spmat_array_t *arr = NULL;
    int i = 0;
    double current_row_sum = 0.0;
    double max_row_sum = 0.0;
    bool_t is_last_in_row = FALSE;

    /* The 1-norm is the max column abs sum. We will go over the transpoed
     * matrix' rows */
    matrix = submatrix->transposed;
    arr = GET_ARRAY(matrix);

    for (i = 0 ; i < matrix->n ; ++i) {
        current_row_sum += fabs(arr->values[i]);

        is_last_in_row = (arr->rowptr[i + 1] > arr->rowptr[i]);
        if (is_last_in_row) {
            max_row_sum = MAX(current_row_sum, max_row_sum);
            current_row_sum = 0;
        }
    }

    return max_row_sum;
}

double
SUBMAT_SPMAT_ARRAY_mult_vmv(const submatrix_t *submatrix, const double *vector)
{
    int i = 0;
    double total_mult = 0.0;
    double row_result = 0.0;
    const matrix_t *matrix = NULL;
    const spmat_array_t *arr = NULL;
    bool_t is_last_in_row = FALSE;

    matrix = submatrix->original;
    arr = GET_ARRAY(matrix);

    for (i = 0 ; i < matrix->n ; ++i) {
        row_result += arr->values[i] * vector[arr->colind[i]];
        is_last_in_row = (row_scanner[1] > row_scanner[0]);
        if (is_last_in_row) {
            total_mult = 
            max_row_sum = MAX(current_row_sum, max_row_sum);
            current_row_sum = 0;
        }
    }

}
#endif
