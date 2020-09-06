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
#include "spmat_array.h"
#include "matrix_raw.h"


/* Functions ************************************************************************************/
result_t
MATRIX_create_matrix(int n, matrix_type_t type, matrix_t **matrix_out)
{
    result_t result = E__UNKNOWN;
    matrix_t *mat = NULL;
    int nnz = 0;

    /* 2. Create the requested matrix_t according to the input dimensions */
    switch (type)
    {
    case MATRIX_TYPE_SPMAT_LIST:
        result = SPMAT_LIST_allocate(n, &mat);
        if (E__SUCCESS != result) {
            goto l_cleanup;
        }
        break;
    case MATRIX_TYPE_SPMAT_ARRAY:
/* TODO: get nonzero values */
        /* The allocator needs the count of non-zero values on the matrix */
#if 0
        result = LAZY_MATRIX_count_nonzero_values(input_matrix, &nnz);
        if (E__SUCCESS != result) {
            goto l_cleanup;
        }
#endif
        nnz = -1;

        result = SPMAT_ARRAY_allocate(n, nnz, &mat);
        if (E__SUCCESS != result) {
            goto l_cleanup;
        }
        break;
    case MATRIX_TYPE_RAW:
        result = MATRIX_RAW_allocate(n, n, &mat);
        if (E__SUCCESS != result) {
            goto l_cleanup;
        }
        break;
    default:
        result = E__UNKNOWN_MATRIX_IMPLEMNTATION;
        goto l_cleanup;
    }

#if 0
    /* 3. Read input into the spmat */
    /* 3.1. Rewind to the first line */
    result = LAZY_MATRIX_rewind(input_matrix);
    if (E__SUCCESS != result) {
        goto l_cleanup;
    }

    /* 3.2. Read */
    for (current_row = 0 ; n > current_row ; ++current_row) {
        result = LAZY_MATRIX_read_next_line(input_matrix);
        if (E__SUCCESS != result) {
            goto l_cleanup;
        }

        mat->add_row(mat, input_matrix->current_line, current_row);
    }
#endif

    /* Success */
    *matrix_out = mat;

    result = E__SUCCESS;
l_cleanup:
    if (E__SUCCESS != result) {
        /* On failure free the smpat */
        if (NULL != mat) {
            mat->free(mat);
            mat = NULL;
        }
    }
#if 0
    LAZY_MATRIX_close(input_matrix);
#endif

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

