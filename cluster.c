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
#include "spmat_list.h"
#include "debug.h"


/* Functions Declarations ************************************************************************/
/**
 * @purpose finding leading eigenvalue of Sparse Matrix
 * @param leading_vector The input leading eigenvector
 * @param prev_vector previous candidate for leading eigenvector from Power Iterations
 * @param eigen_value The calculated leading eigenvalue (output)
 *
 * @return One of result_t values
 *
 * @remark The returned eigen must be freed using MATRIX_free
 */
STATIC
result_t
cluster_calculate_leading_eigenvalue(const matrix_t *matrix,
                                     const double *eigen_vector,
                                     double *eigen_value_out);


/* Functions *************************************************************************************/
STATIC
result_t
cluster_calculate_leading_eigenvalue(const matrix_t *matrix,
        const double *eigen_vector,
        double *eigen_value_out)
{
    result_t result = E__UNKNOWN;
    double *av = NULL;
    double eigen_norm = 0.0;
    double eigen_value_non_normalized = 0.0;
    double eigen_value = 0.0;

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

    /* 1.2. Multiply Av */
    matrix->mult(matrix, eigen_vector, av);

    /* 2. Multiply v with Av */
    eigen_value_non_normalized = VECTOR_scalar_multiply(eigen_vector, av, matrix->n);
    eigen_norm = VECTOR_scalar_multiply(eigen_vector, eigen_vector, matrix->n);

    /* 2.2. Calculate the avarage, or set as 0 */
    if (0 < eigen_norm) {
        eigen_value = eigen_value_non_normalized / eigen_norm;
    } else {
        eigen_value = 0;
    }
    DEBUG_PRINT("eigen value: %f/%f=%f", eigen_value_non_normalized, eigen_norm, eigen_value);

    /* Success */
    *eigen_value_out = eigen_value;

    result = E__SUCCESS;
l_cleanup:

    FREE_SAFE(av);
    return result;
}

result_t
CLUSTER_divide(matrix_t *input,
               matrix_t **group1_out,
               matrix_t **group2_out)
{
    result_t result = E__UNKNOWN;
    double *leading_eigen = NULL;
    double *b_vector = NULL;
    double *s_vector = NULL;
    double leading_eigenvalue = 0.0;
    double *bs = NULL;
    double stbs = 0.0;
    size_t s_ones = 0;
    double onenorm = 0.0;
    int i = 0;
    matrix_t *group1 = NULL;
    matrix_t *group2 = NULL;

    /* 1. Matrix shifting */
    /* 1.1. Calculate 1-norm shifting */
    result = SPMAT_LIST_get_1norm(input, &onenorm);
    if (E__SUCCESS != result) {
        goto l_cleanup;
    }
    DEBUG_PRINT("1norm=%f", onenorm);

    /* 1.2. Add 1-norm to the diag */
    result = MATRIX_add_diag(input, onenorm);
    if (E__SUCCESS != result) {
        goto l_cleanup;
    }

    /* 2. Calculate leading eigenvector */
    /* 2.1. Create b-vector */
    result = VECTOR_random_vector(input->n, &b_vector);
    if (E__SUCCESS != result) {
        goto l_cleanup;
    }

    /* 2.2. Calculate eigen */
    result = MATRIX_calculate_eigen(input, b_vector, &leading_eigen);
    if (E__SUCCESS != result) {
        goto l_cleanup;
    }
    /* 2.2. Free the b-vector */
    FREE_SAFE(b_vector);

    /* 3. Calculate eigenvalue */
    /* 3.1. Calculate eigenvalue plus 1-norm */ 
    result = cluster_calculate_leading_eigenvalue(input, leading_eigen, &leading_eigenvalue);
    if (E__SUCCESS != result) {
        goto l_cleanup;
    }

    /* 3.2. Decrease 1-norm to the diag */
    result = MATRIX_add_diag(input, -onenorm);
    if (E__SUCCESS != result) {
        goto l_cleanup;
    }

    /* 3.2. Decrease 1-norm */ 
    leading_eigenvalue -= onenorm;

    /* 3. Check divisibility #1 */
    if (0 >= leading_eigenvalue) {
        /* 4.1. Network is indivisable */
        DEBUG_PRINT("Network is indivisable! leading eigenvalue is %f", leading_eigenvalue);
        *group1_out = NULL;
        *group2_out = NULL;
        result = E__SUCCESS;
        goto l_cleanup;
    }

    /* 4. Calculate S-vector */
    /* 4.1. Allocate s-vector */
    s_vector = (double *)malloc(sizeof(*s_vector) * input->n);
    if (NULL == s_vector) {
        result = E__MALLOC_ERROR;
        goto l_cleanup;
    }

    /* 4.2. Normalize leading eigenvector */
    result = VECTOR_normalize(leading_eigen, input->n);
    if (E__SUCCESS != result) {
        goto l_cleanup;
    }

    /* 4.1. Calculate s-vector */
    for (i = 0 ; i < input->n; ++i) {
        if (0 < leading_eigen[i]) {
            s_vector[i] = 1;
            ++s_ones;
        } else {
            s_vector[i] = -1;
        }
    }

    /* 5. Calculating stbs */
    /* 5.1. Allocate stbs */
    bs = (double *)malloc(sizeof(*bs) * input->n);
    if (NULL == bs) {
        result = E__MALLOC_ERROR;
        goto l_cleanup;
    }

    /* 5.2. Bs: Multiply B with s */
    input->mult(input, s_vector, bs);

    /* 5.3. sTBs: Multiply s_transposed with Bs */
    stbs = VECTOR_scalar_multiply(s_vector, bs, input->n);

    /* 5.4. Free bs */
    FREE_SAFE(bs);

    /* 6. Check divisibility #2 */
    if (0 >= stbs) {
        /* 7.1. Network is indivisable */
        *group1_out = NULL;
        *group2_out = NULL;
        return result;
    }

    /* 7. Generate divised group */
    /* If we are here the netwrk is divisible */
    result = SPMAT_LIST_divide_matrix(input, s_vector, &group1, &group2);
    if (E__SUCCESS != result) {
        goto l_cleanup;
    }

    /* Success */
    *group1_out = group1;
    *group2_out = group2;
    result = E__SUCCESS;

l_cleanup:
    if (E__SUCCESS != result) {
        MATRIX_FREE_SAFE(group1);
        MATRIX_FREE_SAFE(group2);
    }
    FREE_SAFE(b_vector);
    FREE_SAFE(leading_eigen);
    FREE_SAFE(bs);
    FREE_SAFE(s_vector);

    return result;
}

result_t
CLUSTER_divide_repeatedly(matrix_t *matrix)
{
    result_t result = E__UNKNOWN;
    double *line_vector_tmp = NULL;
    int i = 0;

    /* 0. Input validation */
    if (NULL == matrix) {
        result = E__NULL_ARGUMENT;
        goto l_cleanup;
    }

    result = SPMAT_LIST_decrease_rows_sums_from_diag(matrix);
    if (E__SUCCESS != result) {
        goto l_cleanup;
    }

    for (i = 0 ; i < matrix->n ; ++i) {
        line_vector_tmp[i] = mat
    }

    

    /* Success */
    result = E__SUCCESS;
l_cleanup:

    return result;
}
