/**
 * @file eigen.c
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
#include "vector.h"
#include "config.h"


/* Functions ************************************************************************************/
result_t
MATRIX_calculate_eigen(const matrix_t *input,
                       double *b_vector,
                       double **eigen_out)
{
    result_t result = E__UNKNOWN;
    double * vector_res = NULL;
    double * prev_vector_res = NULL;
    double * temp = NULL;

    if ((NULL == input) || (NULL == b_vector) || (NULL == eigen_out)) {
        result = E__NULL_ARGUMENT;
        goto l_cleanup;
    }

    /* 1. Allocate vectors */
    /* 1.1. Create two vectors: for result, and previous result */
    prev_vector_res = (double *)malloc(sizeof(*prev_vector_res) * input->n);
    if (NULL == prev_vector_res) {
        result = E__MALLOC_ERROR;
        goto l_cleanup;
    }

    /* 2. Multiply the matrix and vector until the previous result is close enough */
    /* 2.1. The first initialization works as the given b_vector was the loop's vector_res */
    vector_res = b_vector;
    do {
        /* Swap */
        temp = vector_res;
        vector_res = prev_vector_res;
        prev_vector_res = temp;

        input->mult(input, prev_vector_res, vector_res);
        result = VECTOR_normalize(vector_res, input->n);
        if (E__SUCCESS != result) {
            goto l_cleanup;
        }
    } while (!VECTOR_is_close(vector_res, prev_vector_res, input->n, EPSILON));

    /* Success */
    *eigen_out = vector_res;

    result = E__SUCCESS;

l_cleanup:
    if (E__SUCCESS != result) {
        FREE_SAFE(vector_res);
    }

    return result;
}

