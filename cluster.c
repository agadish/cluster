/*
 * @file cluster.c
 * @purpose Compute algorithm 1 from the assignment
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
#include "results.h"
#include "vector.h"


/* Functions Declarations ***********************************************************************/
static
result_t
cluster_calculate_leading_eigenvalue(const matrix_t *matrix,
        const double *eigen_vector,
        double *eigen_value_out);


/* Functions ************************************************************************************/
static
result_t
cluster_calculate_leading_eigenvalue(const matrix_t *matrix,
        const double *eigen_vector,
        double *eigen_value_out)
{
    result_t result = E__UNKNOWN;
    double total_ratio = 0.0;
    size_t number_of_ratios = 0;
    double eigen_value = 0.0;
    double *av = NULL;
    int i = 0;

    /* 0. Input validation */
    if ((NULL == matrix) || (NULL == eigen_vector) || (NULL == eigen_value_out)) {
        result = E__NULL_ARGUMENT;
        goto l_cleanup;
    }

    /* 1. Multiply the matrice with the vector */
    /* 1.1. Allocate result */
    av = (double *)malloc(sizeof(*eigen_vector) * matrix->n);
    if (NULL == av) {
        result = E__MALLOC_ERROR;
        goto l_cleanup;
    }

    /* 1.2. Multiply */
    matrix->mult(matrix, eigen_vector, av);

    /* 2. Calculate the avarage ratio between elements of Av and the eigen */
    /* 2.1. Sum all ratios */
    for (i = 0 ; i < matrix->n ; ++i) {
        /* Av=gv ---> g = (Av)/v
         * We'll calculate the avarage ratio. v[i] must not be zero */
        if (0 != eigen_vector[i]) {
            ++number_of_ratios;
            total_ratio += av[i] / eigen_vector[i];
        }
    }

    /* 2.2. Calculate the avarage, or set as 0 */
    if (0 < number_of_ratios) {
        eigen_value = total_ratio / number_of_ratios;
    } else {
        eigen_value = 0;
    }

    /* Success */
    *eigen_value_out = eigen_value;

    result = E__SUCCESS;
l_cleanup:

    FREE_SAFE(av);
    return result;
}

result_t
CLUSTER_divide(matrix_t *input,
               matrix_t **group1,
               matrix_t **group2)
{
    result_t result = E__UNKNOWN;
    double *leading_eigen = NULL;
    double *b_vector = NULL;
    double *s_vector = NULL;
    double leading_eigenvalue = 0.0;
    double *bs = NULL;
    double stbs = 0.0;
    size_t s_ones = 0;
    int i = 0;

    /* 1. Calculate leading eigenvector */
    /* 1.1. Create b-vector */
    result = VECTOR_random_vector(input->n, &b_vector);
    if (E__SUCCESS != result) {
        goto l_cleanup;
    }

    /* 1.2. Create b-vector */
    result = MATRIX_calculate_eigen(input, b_vector, &leading_eigen);
    if (E__SUCCESS != result) {
        goto l_cleanup;
    }
    /* 1.3. Free the b-vector */
    FREE_SAFE(b_vector);

    /* 2. Calculate corresponding eigenvalue */
    result = cluster_calculate_leading_eigenvalue(input, leading_eigen, &leading_eigenvalue);
    if (E__SUCCESS != result) {
        goto l_cleanup;
    }

    /* 3. Check divisibility #1 */
    if (0 >= leading_eigenvalue) {
        /* 4.1. Network is indivisable */
        *group1 = NULL;
        *group2 = NULL;
        goto l_cleanup;
    }

    /* 4. Generate S vector */
    /* 4.1. Allocate */
    s_vector = (double *)malloc(sizeof(*s_vector) * input->n);
    if (NULL == s_vector) {
        result = E__MALLOC_ERROR;
        goto l_cleanup;
    }

    /* 5.1. Calculate */
    for (i = 0 ; i < input->n; ++i) {
        if (0 < leading_eigen[i]) {
            s_vector[i] = 1;
            ++s_ones;
        } else {
            s_vector[i] = -1;
        }
    }

    /* 6. Calculating stbs */
    /* 6.1. Allocate stbs */
    bs = (double *)malloc(sizeof(*bs) * input->n);
    if (NULL == bs) {
        result = E__MALLOC_ERROR;
        goto l_cleanup;
    }

    /* 6.2. Bs: Multiply B with s */
    input->mult(input, s_vector, bs);

    /* 6.3. sTBs: Multiply s_transposed with Bs */
    stbs = VECTOR_scalar_multiply(s_vector, bs, input->n);

    /* 6.4. Free bs */
    FREE_SAFE(bs);

    /* 7. Check divisibility #2 */
    if (0 >= stbs) {
        /* 7.1. Network is indivisable */
        *group1 = NULL;
        *group2 = NULL;
        return result;
    }

    /* 8. Generate divised group */
    /* If we are here the netwrk is divisible */

    /* Im not sure vectors is the right representations of the groups ??? */


l_cleanup:
    if (E__SUCCESS != result) {
    }
    FREE_SAFE(b_vector);
    FREE_SAFE(bs);
    FREE_SAFE(s_vector);

    return result;
}

