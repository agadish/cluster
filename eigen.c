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
    double * original_vector_res = NULL;
    double * vector_res = NULL;
    double * temp = NULL;

    if ((NULL == input) || (NULL == b_vector) || (NULL == eigen_out)) {
        result = E__NULL_ARGUMENT;
        goto l_cleanup;
    }

    /* 1. Allocate vectors */
    /* 1.1. Create two vectors: for result, and previous result */
    vector_res = (double *)malloc(sizeof(*vector_res) * input->n);
    if (NULL == vector_res) {
        result = E__MALLOC_ERROR;
        goto l_cleanup;
    }
    original_vector_res = vector_res;


    /* 2. Multiply the matrix and vector until the previous result is close enough */
    /* 2.1. Swap before first iteration */
    vector_res = b_vector;
    b_vector = original_vector_res;

    /* 2.2. Do power iterations */
    do {
        /* Swap */
        temp = vector_res;
        vector_res = b_vector;
        b_vector = temp;

        input->mult(input, b_vector, vector_res);
        result = VECTOR_normalize(vector_res, input->n);
        if (E__SUCCESS != result) {
            goto l_cleanup;
        }
    } while (!VECTOR_is_close(vector_res, b_vector, input->n, EPSILON));

    /* 3. Make sure the result is in prev_vector_res */
    if (original_vector_res == b_vector) {
        (void)memcpy(original_vector_res, vector_res, sizeof(*vector_res) * input->n);
    }

    /* Success */
    *eigen_out = original_vector_res;

    result = E__SUCCESS;

l_cleanup:
    if (E__SUCCESS != result) {
        FREE_SAFE(original_vector_res);
    }

    return result;
}

