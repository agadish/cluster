/*
 * @file algorithm1.c
 * @purpose Compute algorithm 1 from the assignment
 */
#ifndef __ALGO1_C__
#define __ALGO1_C__


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
#include "results.h"

/* Functions ************************************************************************************/
result_t
calculate_leading_eigenvalue(matrix_t *matrix,
					   matrix_t *eigen_vector,
					   double *eigen_value_out)
{
    result_t result = E__UNKNOWN;
    double denominator;
    double numerator;
    double res;

    if ((NULL == matrix) || (NULL == eigen_vector) || (NULL == eigen_value_out)) {
        result = E__NULL_ARGUMENT;
        goto l_cleanup;
    }

    /* need to add multiplication of matrices */

    res = numerator / denominator;
    *eigen_value_out = res;

    result = E__SUCCESS;
l_cleanup:

    return result;
}

result_t
algorithm1(matrix_t *input,
		matrix_t **group1,
		matrix_t **group2)
{
    result_t result = E__UNKNOWN;
    matrix_t *leading_eigen = NULL;
    matrix_t *prev_eigen = NULL;
    matrix_t *s_vector = NULL;
    matrix_t *transpose= NULL;
    double leading_eigenvalue;
    matrix_t Bs;
    matrix_t sTBs;
    int i = 0; /* For iterations */

    /* 1. Calculate leading eigenvector */
    result = MATRIX_calculate_eigen(input, &leading_eigen, &prev_eigen);
    if (RESULT_SUCCESS != result) {
        goto l_cleanup;
    }

    /* 2. Calculate corresponding eigenvalue */
    result = calculate_leading_eigenvalue(&leading_eigen, &prev_eigen, &leading_eigenvalue);
    if (RESULT_SUCCESS != result) {
        goto l_cleanup;
    }

    /* 3. Check divisibility #1 */
    if (leading_eigenvalue <= 0) {
    	//Network is indivisable
    	*group1 = NULL;
    	*group2 = NULL;
    	return result;
    }

    /* 4. Generate S vector */
    result = MATRIX_allocate(1, input->n, &s_vector);
    if (E__SUCCESS != result) {
        goto l_cleanup;
    }

    for (i = 0 ; i < input->n; ++row) {
    	if(MATRIX_AT(leading_eigen, i, 0) > 0) {
    		MATRIX_AT(s_vector, i, 0) = 1;
    	} else {
    		MATRIX_AT(s_vector, i, 0) = -1;
    	}
    }

    /* 5. Check divisibility #2 */
    /* 5.1 Calculating sTBs */
    MATRIX_row_vector_transpose(s_vector, &s_transpose);
    MATRIX_mat_vector_multiply(input, s_vector, &Bs)
    MATRIX_mat_vector_multiply(s_transpose, Bs, &sTBs)

    //sTBs has only one value at (0,0)

	/* 5.2 Actually Checking divisibility */
    if (MATRIX_AT(sTBs, 0, 0) <= 0) {
    	//Network is indivisable
    	*group1 = NULL;
    	*group2 = NULL;
    	return result;
    }

    /* 6. Generate divised group */
    //If we are here the netwrk is divisible

    /////////// ?? Im not sure vectors is the right representations of the groups ???


    l_cleanup:
		/// ??
        if (E__SUCCESS != result) {
        	////??
        }

        return result;
}


#endif /* __ALGO1_C__ */
