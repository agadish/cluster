/**
 * @file matrix.c
 * @purpose A generic matrix
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
#include "spmat_list.h"

/* Functions ************************************************************************************/
result_t
MATRIX_create_matrix(int n, matrix_type_t type, matrix_t **matrix_out)
{
    result_t result = E__UNKNOWN;
    matrix_t *mat = NULL;

    /* 2. Create the requested matrix_t according to the input dimensions */
    switch (type)
    {
    case MATRIX_TYPE_SPMAT_LIST:
        result = SPMAT_LIST_allocate(n, &mat);
        if (E__SUCCESS != result) {
            goto l_cleanup;
        }
        break;
    default:
        result = E__UNKNOWN_MATRIX_IMPLEMNTATION;
        goto l_cleanup;
    }

    /* Success */
    *matrix_out = mat;

    result = E__SUCCESS;
l_cleanup:
    if (E__SUCCESS != result) {
        /* On failure free the smpat */
        MATRIX_FREE_SAFE(mat);
    }

    return result;
}

#ifdef NEED_COL_VECTOR_TRANSPOSE
result_t
MATRIX_col_vector_transpose(matrix_t *vector_in, matrix_t **vector_out)
{
    result_t result = E__UNKNOWN;
    matrix_t * res_vector = NULL;
    int length = vector_in->header.rows;
    int col = 0; /* For iterations */

    result = MATRIX_allocate(length, 1, &vector);
    if (E__SUCCESS != result) {
        goto l_cleanup;
    }

    for (col = 0 ; col < length ; ++col) {
        MATRIX_AT(res_vector, 0, col) = MATRIX_AT(vector_in, col, 0);
    }

    *vector_out = res_vector;

    result = E__SUCCESS;

l_cleanup:
    if (E__SUCCESS != result) {
        FREE_SAFE(res_vector);
    }

    return result;
}
#endif

result_t
MATRIX_add_diag(matrix_t *matrix, double onenorm)
{
    result_t result = E__UNKNOWN;
    double * temp_row_vector = NULL;
    int i = 0;

    /* 0. Input validation */
    if (NULL == matrix) {
        result = E__NULL_ARGUMENT;
        goto l_cleanup;
    }

    /* 1. Allocate temp row vector */
    temp_row_vector = (double *)malloc(sizeof(*temp_row_vector) * matrix->n);
    if (NULL == temp_row_vector) {
        result = E__MALLOC_ERROR;
        goto l_cleanup;
    }
       
    /* 1.2. Zero vector */
    (void)memset(temp_row_vector, 0, sizeof(*temp_row_vector) * matrix->n);

    /* 2. Add (i,i) for each row */
    for (i = 0 ; i < matrix->n ; ++i) {
        /* 2.1. Set i-th cell as 1-norm */
        temp_row_vector[i] = onenorm;

        result = MATRIX_ADD_ROW(matrix, temp_row_vector, i);
        if (E__SUCCESS != result) {
            goto l_cleanup;
        }

        /* 2.3. Restore i-th cell to 0 */
        temp_row_vector[i] = 0.0;
    }

    /* Success */
    result = E__SUCCESS;
l_cleanup:

    FREE_SAFE(temp_row_vector);

    return result;
}

