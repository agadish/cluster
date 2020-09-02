/**
 * @file matrix.c
 * @purpose 
 */

/* Includes **************************************************************************************/
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <math.h>

#include "matrix.h"
#include "common.h"
#include "eigen.h"



/* Functions ************************************************************************************/
result_t
MATRIX_calculate_eigen(spmat *input,
                const matrix_t *b_vector,
                matrix_t **eigen_out,
				matrix_t **prev_vector)
{
    result_t result = E__UNKNOWN;
    matrix_t * vector_res = NULL;
    matrix_t * prev_vector_res = NULL;
    matrix_t * temp = NULL;

    if ((NULL == input) || (NULL == b_vector) || (NULL == eigen_out)) {
        result = E__NULL_ARGUMENT;
        goto l_cleanup;
    }

    /* 1. Allocate vectors */
    /* 1.1. Create two vectors: for result, and previous result */
    result = MATRIX_allocate(input->n, 1, &vector_res);
    if (E__SUCCESS != result) {
        goto l_cleanup;
    }
    /* 1.2. Create eigen vector */
    result = MATRIX_allocate(input->n, 1, &prev_vector_res);
    if (E__SUCCESS != result) {
        goto l_cleanup;
    }

    /* 2. Multiply matrix and b-vector for the first time */
    input->mult(input, b_vector->data, vector_res->data);
    result = MATRIX_normalize(vector_res);
    if (E__SUCCESS != result) {
        goto l_cleanup;
    }

    /* 3. Multiply the matrix and vector until the previous result is close enough */
    do {
        /* Swap */
        temp = vector_res;
        vector_res = prev_vector_res;
        prev_vector_res = temp;

        input->mult(input, prev_vector_res->data, vector_res->data);
        result = MATRIX_normalize(vector_res);
        if (E__SUCCESS != result) {
            goto l_cleanup;
        }
    } while (!MATRIX_is_close(vector_res, prev_vector_res, 0.00001));

    /* Success */
    *eigen_out = vector_res;
    *prev_vector = prev_vector_res;

    result = E__SUCCESS;

l_cleanup:
    if (E__SUCCESS != result) {
        MATRIX_free(vector_res);
        MATRIX_free(prev_vector_res);
    }

    return result;
}

