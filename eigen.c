/**
 * @file eigen.c
 * @purpose 
 */

/* Includes ******************************************************************/
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
#include "submatrix.h"
#include "spmat_list.h"


/* Functions *****************************************************************/
result_t
EIGEN_calculate_eigen(const submatrix_t *smat,
                      double *b_vector,
                      double onenorm,
                      double *eigen)
{
    result_t result = E__UNKNOWN;
    double *original_vector_res = NULL;
    double *vector_res = NULL;
    double *temp = NULL;
    int n = 0;

    if ((NULL == smat) || (NULL == b_vector)) {
        result = E__NULL_ARGUMENT;
        goto l_cleanup;
    }

    /* 1. Allocate vectors */
    /* 1.1. Create two vectors: for result, and previous result */
    n = smat->adj->original->n;
    original_vector_res = eigen;

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

        SUBMAT_SPMAT_LIST_mult(smat, b_vector, onenorm, vector_res);
        result = VECTOR_normalize(vector_res, n);
        if (E__SUCCESS != result) {
            goto l_cleanup;
        }
    } while (!VECTOR_is_close(vector_res, b_vector, n, EPSILON));

    /* 3. Make sure the result is in prev_vector_res */
    if (original_vector_res == b_vector) {
        (void)memcpy(original_vector_res, vector_res, sizeof(*vector_res) * n);
    }

    /* Success */
    result = E__SUCCESS;

l_cleanup:

    return result;
}

