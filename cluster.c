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
#include "list.h"
#include "division_file.h"


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
STATIC
result_t
cluster_optimize_division_iteration(matrix_t *matrix,
                              double *s_vector,
                              double *improve,
                              int *indices,
                              double *temp,
                              double *delta_q_out);

STATIC
result_t
cluster_optimize_division(matrix_t *matrix,
                          double *s_vector);

/**
 * @purpose divide a network to two groups
 * @param input Matrix to divide
 * @param s_vector A pre-allocated input->n sized s-vector
 *
 * @return One of result_t values, E__UNDIVISIBLE_NETWORK if network is undivisible
 *
 * @remark The returned eigen must be freed using MATRIX_free
 */
STATIC
result_t
cluster_divide(matrix_t *input,
               double *s_vector);

STATIC
result_t
cluster_sub_divide(matrix_t *input,
                   double *s_vector);

STATIC
result_t
cluster_sub_divide_optimized(matrix_t *input,
                           double *s_vector);

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
    MATRIX_MULT(matrix, eigen_vector, av);

    /* 2. Multiply v with Av */
    eigen_value_non_normalized = VECTOR_scalar_multiply(eigen_vector, av, matrix->n);
    eigen_norm = VECTOR_scalar_multiply(eigen_vector, eigen_vector, matrix->n);

    /* 2.2. Calculate the avarage, or set as 0 */
    if (IS_POSITIVE(eigen_norm)) {
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

STATIC
result_t
cluster_divide(matrix_t *input,
               double *s_vector)
{
    result_t result = E__UNKNOWN;
    double *leading_eigen = NULL;
    double *b_vector = NULL;
    double leading_eigenvalue = 0.0;
    double *bs = NULL;
    double stbs = 0.0;
    size_t s_ones = 0;
    double onenorm = 0.0;
    int i = 0;

    /* 1. Matrix shifting */
    /* 1.1. Calculate 1-norm shifting */
    onenorm = MATRIX_GET_1NORM(input);
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
        result = E__UNDIVISIBLE_NETWORK;
        goto l_cleanup;
    }

    /* 4. Calculate S-vector */
    /* 4.1. Allocate s-vector */
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
    MATRIX_MULT(input, s_vector, bs);

    /* 5.3. sTBs: Multiply s_transposed with Bs */
    stbs = VECTOR_scalar_multiply(s_vector, bs, input->n);

    /* 5.4. Free bs */
    FREE_SAFE(bs);

    /* 6. Check divisibility #2 */
    if (0 >= stbs) {
        /* 7.1. Network is indivisable */
        result = E__UNDIVISIBLE_NETWORK;
        goto l_cleanup;
    }

    /* Success */
    result = E__SUCCESS;

l_cleanup:
    FREE_SAFE(b_vector);
    FREE_SAFE(leading_eigen);
    FREE_SAFE(bs);

    return result;
}

/* Algorithm 2 */
STATIC
result_t
cluster_sub_divide(matrix_t *input,
                   double *s_vector)
{
    result_t result = E__UNKNOWN;

    result = MATRIX_DECREASE_ROWS_SUMS_FROM_DIAG(input);
    if (E__SUCCESS != result) {
        goto l_cleanup;
    }

    result = cluster_divide(input, s_vector);
    if ( E__SUCCESS != result) {
        goto l_cleanup;
    }

    result = E__SUCCESS;
l_cleanup:

    return result;
}

STATIC
result_t
cluster_sub_divide_optimized(matrix_t *input,
                           double *s_vector)
{
    result_t result = E__UNKNOWN;

    result = cluster_sub_divide(input, s_vector);
    if (E__SUCCESS != result) {
        goto l_cleanup;
    }

    result = cluster_optimize_division(input, s_vector);
    if (E__SUCCESS != result) {
        goto l_cleanup;
    }

    result = E__SUCCESS;
l_cleanup:

    return result;
}

result_t
CLUSTER_divide_repeatedly(matrix_t *initial_matrix, division_file_t *output_file)
{
    result_t result = E__UNKNOWN;
    result_t division_result = E__UNKNOWN;
    matrix_t **p_group = NULL;
    matrix_t *current_matrix = NULL;
    size_t p_group_length = 0;
    double *s_vector = NULL;
    matrix_t *group1 = NULL;
    matrix_t *group2 = NULL;
    int *temp_s_indexes = NULL;

    /* 0. Input validation */
    if ((NULL == initial_matrix) || (NULL == output_file)) {
        result = E__NULL_ARGUMENT;
        goto l_cleanup;
    }

    /* 1. Initializations */
    p_group = (matrix_t **)malloc(initial_matrix->n * sizeof(*p_group));
    if (NULL == p_group) {
        result = E__MALLOC_ERROR;
        goto l_cleanup;
    }

    s_vector = (double *)malloc(initial_matrix->n * sizeof(*s_vector));
    if (NULL == s_vector) {
        result = E__MALLOC_ERROR;
        goto l_cleanup;
    }

    temp_s_indexes = (int *)malloc(initial_matrix->n * sizeof(*temp_s_indexes));
    if (NULL == temp_s_indexes) {
        result = E__MALLOC_ERROR;
        goto l_cleanup;
    }

    /* 1. Initailize p-group */
    p_group[0] = initial_matrix;
    ++p_group_length;
    /* Prevent input from being double freed */
    initial_matrix = NULL;

    for (; 0 < p_group_length ; --p_group_length)
    {
        /* Take next matrix */
        current_matrix = p_group[p_group_length - 1];
        /* (void)printf("\n\n-------------%s: diving matrix\n", __func__); */
        /* SPMAT_LIST_print("Original", current_matrix); */

        division_result = cluster_sub_divide_optimized(current_matrix,
                                                       s_vector);
        if (E__SUCCESS != division_result) {
            if (E__UNDIVISIBLE_NETWORK == division_result) {
                /* Matrix is undivisibe - write it */
                result = DIVISION_FILE_write_matrix(output_file, current_matrix);
                if (E__SUCCESS != result) {
                    goto l_cleanup;
                }
                /* printf("no division"); */

                MATRIX_FREE_SAFE(current_matrix);
                p_group[p_group_length - 1] = NULL;

                /* Get next matrix from the p-group */
                continue;
            } else {
                result = division_result;
                goto l_cleanup;
            }
        }

        /* If we are here the netwrk is divisible */
        result = MATRIX_DIVIDE(current_matrix,
                               s_vector,
                               temp_s_indexes,
                               &group1,
                               &group2);
        if (E__SUCCESS != result) {
            goto l_cleanup;
        }
        MATRIX_FREE_SAFE(current_matrix);
        p_group[p_group_length - 1] = NULL;
        /* printf("divide result-matrix1=%p matrix2=%p\n", (void *)group1, (void *)group2); */

        if (0 == group1->n) {
            result = DIVISION_FILE_write_matrix(output_file, group2);
            /* SPMAT_LIST_print("group2 trivial", group2); */
            if (E__SUCCESS != result) {
                goto l_cleanup;
            }

            MATRIX_FREE_SAFE(group1);
            MATRIX_FREE_SAFE(group2);
            continue;
        }
        if (0 == group2->n) {
            /* SPMAT_LIST_print("group1 trivial", group1); */
            result = DIVISION_FILE_write_matrix(output_file, group1);
            if (E__SUCCESS != result) {
                goto l_cleanup;
            }
            MATRIX_FREE_SAFE(group1);
            MATRIX_FREE_SAFE(group2);
            continue;
        }

        /* SPMAT_LIST_print("Matrix1", group1); */
        /* SPMAT_LIST_print("Matrix2", group2); */

        if (1 == group1->n) {
            result = DIVISION_FILE_write_matrix(output_file, group1);
            if (E__SUCCESS != result) {
                goto l_cleanup;
            }
            MATRIX_FREE_SAFE(group1);
        } else {
            p_group[p_group_length - 1] = group1;
            ++p_group_length;
        }
        
        if (1 == group2->n) {
            result = DIVISION_FILE_write_matrix(output_file, group2);
            if (E__SUCCESS != result) {
                goto l_cleanup;
            }
            MATRIX_FREE_SAFE(group2);
        } else {
            p_group[p_group_length - 1] = group2;
            ++p_group_length;
        }

    }

    /* Success */
    result = E__SUCCESS;
l_cleanup:

    /* Will be freed in case of failure before being inserted to p_group */
    MATRIX_FREE_SAFE(initial_matrix);
    while (p_group_length > 0) {
        --p_group_length;
        MATRIX_FREE_SAFE(p_group[p_group_length]);
    }

    FREE_SAFE(p_group);
    FREE_SAFE(s_vector);
    FREE_SAFE(temp_s_indexes);

    return result;
}

double
cluster_calculate_q(const matrix_t *matrix,
                    const double *s_vector,
                    double *temp)
{
    double q = 0.0;
    MATRIX_MULT(matrix, s_vector, temp);
    q = VECTOR_scalar_multiply(s_vector, temp, matrix->n);

    return q;
}

STATIC
result_t
cluster_optimize_division(matrix_t *matrix,
                          double *s_vector)
{
    result_t result = E__UNKNOWN;
    double delta_q = 0.0;
    double *temp = NULL;
    int *indices = NULL;
    double *improve = NULL;

    /* 0. Input validation */

    if ((NULL == matrix) || (NULL == s_vector)) {
        result = E__NULL_ARGUMENT;
        goto l_cleanup;
    }

    /* 1. Allocate memory */
    temp = (double *)malloc(sizeof(*temp) * matrix->n);
    if (NULL == temp) {
        result = E__MALLOC_ERROR;
        goto l_cleanup;
    }

    indices = (int *)malloc(sizeof(*indices) * matrix->n);
    if (NULL == indices) {
        result = E__MALLOC_ERROR;
        goto l_cleanup;
    }

    improve = (double *)malloc(sizeof(*improve) * matrix->n);
    if (NULL == improve) {
        result = E__MALLOC_ERROR;
        goto l_cleanup;
    }

    /* 2. Do iterations as long as there's improvement */
    do {
        result = cluster_optimize_division_iteration(matrix,
                s_vector,
                improve,
                indices,
                temp,
                &delta_q);
        if (E__SUCCESS != result) {
            goto l_cleanup;
        }
    } while (IS_POSITIVE(delta_q));

    result = E__SUCCESS;
l_cleanup:

    FREE_SAFE(temp);
    FREE_SAFE(improve);
    FREE_SAFE(indices);

    return result;
}

STATIC
result_t
cluster_optimize_division_iteration(matrix_t *matrix,
                                    double *s_vector,
                                    double *improve,
                                    int *indices,
                                    double *temp,
                                    double *delta_q_out)
{
    result_t result = E__UNKNOWN;
    list_t *unmoved_scores = NULL;
    node_t *scanner = NULL;
    node_t *max_unmoved = NULL;
    double q_0 = 0.0;
    double delta_q = 0.0;
    double q_score = 0.0;
    int i = 0;
    int k = 0;
    /* A negative value will never be the max since the last improvement is always 0 */
    double max_improvement_value = -1.0;
    int max_improvement_index = 0;

    /* 1. Initialize list */
    result = LIST_range(matrix->n, &unmoved_scores);
    if (E__SUCCESS != result) {
        goto l_cleanup;
    }

    for (i = 0 ; i < matrix->n ; ++i) {
        /* 2. Calculate Q_0 */
        q_0 = cluster_calculate_q(matrix, s_vector, temp);

        /* 3. Computing DeltaQ for the move of each unmoved vertex */
        for (scanner = unmoved_scores->first ; NULL != scanner ; scanner = scanner->next) {
            /* Calculate score when moving k */
            k = scanner->index;
            s_vector[k] *= -1;
            q_score = cluster_calculate_q(matrix, s_vector, temp);
            s_vector[k] *= -1;
            scanner->value = q_score - q_0;

            /* Update max score */
            if (NULL == max_unmoved) {
                max_unmoved = scanner;
            } else if (max_unmoved->value < scanner->value) {
                max_unmoved = scanner;
            }
        }

        /* 4. Move vertex max_score_index with a maximal score */
        s_vector[max_unmoved->index] = -s_vector[max_unmoved->index];
        indices[i] = max_unmoved->index;
        if (0 == i) {
            improve[i] = max_unmoved->value;
        } else {
            improve[i] = max_unmoved->value + improve[i - 1];
        }

        result = LIST_remove_node(unmoved_scores, max_unmoved);
        if (E__SUCCESS != result) {
            goto l_cleanup;
        }
        max_unmoved = NULL;

        /* 5. Update max improvement */
        if (max_improvement_value < improve[i]) {
            max_improvement_index = i;
            max_improvement_value = improve[i];
        }
    }

    /* 6. Apply the max improvement to the s-vector */
    for (i = matrix->n - 1 ; i > max_improvement_index ; --i) {
        s_vector[indices[i]] *= -1;
    }

    if (max_improvement_index == matrix->n - 1) {
        delta_q = 0.0;
    } else {
        delta_q = improve[max_improvement_index];
    }

    *delta_q_out = delta_q;

    /* Success */
    result = E__SUCCESS;
l_cleanup:

    LIST_destroy(unmoved_scores);
    unmoved_scores = NULL;

    return result;
}


