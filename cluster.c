/*
 * @file cluster.c
 * @purpose Compute algorithm 1 from the assignment
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
#include "results.h"
#include "vector.h"
#include "spmat_list.h"
#include "spmat_array.h"
#include "debug.h"
#include "list.h"
#include "division_file.h"
#include "submatrix.h"


/* Functions Declarations ****************************************************/
/**
 * @purpose finding leading eigenvalue of Sparse Matrix
 * @param leading_vector The input leading eigenvector
 * @param prev_vector Previous candidate for leading eigenvector from Power
 *                    Iterations
 * @param eigen_value The calculated leading eigenvalue (output)
 *
 * @return One of result_t values
 *
 * @remark The returned eigen must be freed using MATRIX_free
 */
static
result_t
cluster_calculate_leading_eigenvalue(const submatrix_t *matrix,
                                     const double *eigen_vector,
                                     double *eigen_value_out);
#if 0
static
result_t
cluster_optimize_division_iteration(submatrix_t *matrix,
                                    double *s_vector,
                                    double *improve,
                                    int *indices,
                                    double *delta_q_out);
#endif
static
result_t
cluster_optimize_division_iteration2(submatrix_t *smat,
                                    double *s_vector,
                                    double *improve,
                                    int *indices,
                                    double *delta_q_out);

static
result_t
cluster_optimize_division(submatrix_t *smat,
                          double *s_vector);

/**
 * @purpose divide a network to two groups
 * @param input Matrix to divide
 * @param s_vector A pre-allocated input->n sized s-vector
 *
 * @return One of result_t values, E__UNDIVISIBLE_NETWORK if network is
 *         undivisible
 *
 * @remark The returned eigen must be freed using MATRIX_free
 */
static
result_t
cluster_divide(submatrix_t *smat,
               double *temp_b_vector,
               double *temp_eigen_vector,
               double *s_vector);

static
result_t
cluster_sub_divide_optimized(submatrix_t *smat,
                             double *temp_b_vector,
                             double *temp_eigen_vector,
                             double *s_vector);

static
result_t
cluster_create_submatrix(const adjacency_matrix_t *adj, submatrix_t **smat_out);


/* Functions *****************************************************************/
static
result_t
cluster_calculate_leading_eigenvalue(const submatrix_t *matrix,
                                     const double *eigen_vector,
                                     double *eigen_value_out)
{
    result_t result = E__UNKNOWN;
    double eigen_value_numerator = 0.0;
    double eigen_value_denominator = 0.0;
    double eigen_value = 0.0;

    /* 0. Input validation */
    if ((NULL == matrix) || (NULL == eigen_vector) ||
        (NULL == eigen_value_out)) {
        result = E__NULL_ARGUMENT;
        goto l_cleanup;
    }

    /* 1. Multiply the matrice with the vector */
    /* 2.1. Caculate eigen norm */
    eigen_value_numerator = SUBMAT_SPMAT_LIST_calculate_q(matrix,
                                                          eigen_vector);
    eigen_value_denominator = VECTOR_scalar_multiply(eigen_vector,
                                                     eigen_vector,
                                                     matrix->adj->original->n);

    /* 2.2. Calculate the avarage, or set as 0 */
    if (IS_POSITIVE(eigen_value_denominator)) {
        eigen_value = eigen_value_numerator / eigen_value_denominator;
    } else {
        eigen_value_denominator = 0;
    }
    DEBUG_PRINT("eigen value: %f/%f=%f",
                eigen_value_numerator,
                eigen_value_denominator,
                eigen_value);

    /* Success */
    *eigen_value_out = eigen_value;

    result = E__SUCCESS;
l_cleanup:

    return result;
}

static
result_t
cluster_divide(submatrix_t *smat,
               double *temp_b_vector,
               double *temp_eigen_vector,
               double *s_vector)
{
    result_t result = E__UNKNOWN;
    double leading_eigenvalue = 0.0;
    double stbs = 0.0;
    size_t s_ones = 0;
    double onenorm = 0.0;
    int i = 0;
    int n = 0;


    /* 1. Calculate leading eigenvector */
    n = smat->g_length;

    /* 1.1. Calculate the 1-norm of the matrix, using b-vector as temp vector */
    onenorm = SUBMAT_SPMAT_LIST_get_1norm(smat, temp_b_vector);

    /* 1.2. Randomize b-vector */
    VECTOR_random_vector(n, temp_b_vector);

    /* 1.3. Add the 1-norm to the diag for the eigen calculation */
    smat->add_to_diag = onenorm;

    /* 1.4. Calculate eigen vector */
    result = EIGEN_calculate_eigen(smat,
                                   temp_b_vector,
                                   temp_eigen_vector);
    if (E__SUCCESS != result) {
        goto l_cleanup;
    }

    /* 3. Calculate eigenvalue */
    /* 3.1. Calculate eigenvalue plus 1-norm */ 
    result = cluster_calculate_leading_eigenvalue(smat,
                                                  temp_eigen_vector,
                                                  &leading_eigenvalue);
    if (E__SUCCESS != result) {
        goto l_cleanup;
    }

    /* 3.2. Decrease 1-norm from the result, restore the diag variable */ 
    smat->add_to_diag = 0.0;
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
    result = VECTOR_normalize(temp_eigen_vector, n);
    if (E__SUCCESS != result) {
        goto l_cleanup;
    }

    /* 4.1. Calculate s-vector */
    for (i = 0 ; i < n; ++i) {
        if (0 < temp_eigen_vector[i]) {
            s_vector[i] = 1;
            ++s_ones;
        } else {
            s_vector[i] = -1;
        }
    }

    /* 5. Calculating stbs */
    stbs = SUBMAT_SPMAT_LIST_calculate_q(smat, s_vector);

    /* 6. Check divisibility #2 */
    if (0 >= stbs) {
        /* 7.1. Network is indivisable */
        result = E__UNDIVISIBLE_NETWORK;
        goto l_cleanup;
    }

    /* Success */
    result = E__SUCCESS;

l_cleanup:

    return result;
}


static
result_t
cluster_sub_divide_optimized(submatrix_t *smat,
                             double *temp_b_vector,
                             double *temp_eigen_vector,
                             double *s_vector)
{
    result_t result = E__UNKNOWN;

    result = cluster_divide(smat, temp_b_vector, temp_eigen_vector,s_vector);
    if (E__SUCCESS != result) {
        goto l_cleanup;
    }

    result = cluster_optimize_division(smat, s_vector);
    if (E__SUCCESS != result) {
        goto l_cleanup;
    }

    result = E__SUCCESS;
l_cleanup:

    return result;
}

static
result_t
cluster_create_submatrix(const adjacency_matrix_t *adj, submatrix_t **smat_out)
{
    result_t result = E__UNKNOWN;
    submatrix_t *smat = NULL;
    int i = 0;

    /* 1. Create submatrix */
    result = SUBMATRIX_create(adj, &smat);
    if (E__SUCCESS != result) {
        goto l_cleanup;
    }

    /* 2. Initialise g-vector */
    for (i = 0 ; i < smat->g_length ; ++i) {
        smat->g[i] = i;
    }

    /* Success */
    *smat_out = smat;

    result = E__SUCCESS;
l_cleanup:
    if (E__SUCCESS != result) {
        SUBMATRIX_FREE_SAFE(smat);
    }

    return result;
}

result_t
CLUSTER_divide_repeatedly(adjacency_matrix_t *adj, division_file_t *output_file)
{
    result_t result = E__UNKNOWN;
    result_t division_result = E__UNKNOWN;
    submatrix_t *smat = NULL;
    submatrix_t **p_group = NULL;
    submatrix_t *current_matrix = NULL;
    size_t p_group_length = 0;
    double *s_vector = NULL;
    submatrix_t *group1 = NULL;
    submatrix_t *group2 = NULL;
    double *temp_b_vector = NULL;
    double *temp_eigen_vector = NULL;

    /* 0. Input validation */
    if ((NULL == adj) || (NULL == output_file)) {
        result = E__NULL_ARGUMENT;
        goto l_cleanup;
    }

    /* 1. Initializations */
    p_group = (submatrix_t **)malloc(adj->original->n * sizeof(*p_group));
    if (NULL == p_group) {
        result = E__MALLOC_ERROR;
        goto l_cleanup;
    }

    s_vector = (double *)malloc(adj->original->n * sizeof(*s_vector));
    if (NULL == s_vector) {
        result = E__MALLOC_ERROR;
        goto l_cleanup;
    }

    temp_b_vector = (double *)malloc(adj->original->n * sizeof(*temp_b_vector));
    if (NULL == temp_b_vector) {
        result = E__MALLOC_ERROR;
        goto l_cleanup;
    }

    temp_eigen_vector = (double *)malloc(adj->original->n *
                                         sizeof(*temp_eigen_vector));
    if (NULL == temp_eigen_vector) {
        result = E__MALLOC_ERROR;
        goto l_cleanup;
    }

    result = cluster_create_submatrix(adj, &smat);
    if (E__SUCCESS != result) {
        goto l_cleanup;
    }

    /* 1. Initailize p-group */
    p_group[0] = smat;
    ++p_group_length;
    /* Prevent input from being double freed */

    for (; 0 < p_group_length ; --p_group_length)
    {
        /* Take next matrix */
        current_matrix = p_group[p_group_length - 1];
        /* (void)printf("\n\n-------------%s: diving matrix\n", __func__); */
        /* SPMAT_LIST_print("Original", current_matrix); */

        division_result = cluster_sub_divide_optimized(current_matrix,
                                                       temp_b_vector,
                                                       temp_eigen_vector,
                                                       s_vector);
        if (E__SUCCESS != division_result) {
            if (E__UNDIVISIBLE_NETWORK == division_result) {
                /* Matrix is undivisibe - write it */
                result = DIVISION_FILE_write_matrix(output_file,
                                                    current_matrix->g,
                                                    current_matrix->g_length);
                if (E__SUCCESS != result) {
                    goto l_cleanup;
                }
                /* printf("no division"); */

                SUBMATRIX_FREE_SAFE(current_matrix);
                p_group[p_group_length - 1] = NULL;

                /* Get next matrix from the p-group */
                continue;
            } else {
                result = division_result;
                goto l_cleanup;
            }
        }

        /* Network is divisible */
        result = SUBMAT_SPMAT_LIST_split(current_matrix,
                                         s_vector,
                                         &group1,
                                         &group2);
        if (E__SUCCESS != result) {
            goto l_cleanup;
        }
        SUBMATRIX_FREE_SAFE(current_matrix);
        p_group[p_group_length - 1] = NULL;
        /* printf("divide result-matrix1=%p matrix2=%p\n", (void *)group1, (void *)group2); */

        if (0 == group1->g_length) {
            result = DIVISION_FILE_write_matrix(output_file,
                                                group2->g,
                                                group2->g_length);
            /* SPMAT_LIST_print("group2 trivial", group2); */
            if (E__SUCCESS != result) {
                goto l_cleanup;
            }

            SUBMATRIX_FREE_SAFE(group1);
            SUBMATRIX_FREE_SAFE(group2);
            continue;
        }
        if (0 == group2->g_length) {
            /* SPMAT_LIST_print("group1 trivial", group1); */
            result = DIVISION_FILE_write_matrix(output_file,
                                                group1->g,
                                                group1->g_length);
            if (E__SUCCESS != result) {
                goto l_cleanup;
            }
            SUBMATRIX_FREE_SAFE(group1);
            SUBMATRIX_FREE_SAFE(group2);
            continue;
        }

        /* SPMAT_LIST_print("Matrix1", group1); */
        /* SPMAT_LIST_print("Matrix2", group2); */

        if (1 == group1->g_length) {
            result = DIVISION_FILE_write_matrix(output_file,
                                                group1->g,
                                                group1->g_length);
            if (E__SUCCESS != result) {
                goto l_cleanup;
            }
            SUBMATRIX_FREE_SAFE(group1);
        } else {
            p_group[p_group_length - 1] = group1;
            ++p_group_length;
        }
        
        if (1 == group2->g_length) {
            result = DIVISION_FILE_write_matrix(output_file,
                                                group2->g,
                                                group2->g_length);
            if (E__SUCCESS != result) {
                goto l_cleanup;
            }
            SUBMATRIX_FREE_SAFE(group2);
        } else {
            p_group[p_group_length - 1] = group2;
            ++p_group_length;
        }

    }

    /* Success */
    result = E__SUCCESS;
l_cleanup:

    /* Will be freed in case of failure before being inserted to p_group */
    while (p_group_length > 0) {
        --p_group_length;
        SUBMATRIX_FREE_SAFE(p_group[p_group_length]);
    }


    FREE_SAFE(p_group);
    FREE_SAFE(s_vector);
    FREE_SAFE(temp_b_vector);
    FREE_SAFE(temp_eigen_vector);

    return result;
}


static
result_t
cluster_optimize_division(submatrix_t *smat,
                          double *s_vector)
{
    result_t result = E__UNKNOWN;
    double delta_q = 0.0;
    int *indices = NULL;
    double *improve = NULL;

    /* 0. Input validation */

    if ((NULL == smat) || (NULL == s_vector)) {
        result = E__NULL_ARGUMENT;
        goto l_cleanup;
    }

    /* 1. Allocate memory */

    indices = (int *)malloc(sizeof(*indices) * smat->g_length);
    if (NULL == indices) {
        result = E__MALLOC_ERROR;
        goto l_cleanup;
    }

    improve = (double *)malloc(sizeof(*improve) * smat->g_length);
    if (NULL == improve) {
        result = E__MALLOC_ERROR;
        goto l_cleanup;
    }

    /* 2. Do iterations as long as there's improvement */
    do {
        result = cluster_optimize_division_iteration2(smat,
                                                     s_vector,
                                                     improve,
                                                     indices,
                                                     &delta_q);
        if (E__SUCCESS != result) {
            goto l_cleanup;
        }
    } while (IS_POSITIVE(delta_q));

    result = E__SUCCESS;
l_cleanup:

    FREE_SAFE(improve);
    FREE_SAFE(indices);

    return result;
}

#if 0
static
result_t
cluster_optimize_division_iteration(submatrix_t *smat,
                                    double *s_vector,
                                    double *improve,
                                    int *indices,
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
    result = LIST_range(smat->g_length, &unmoved_scores);
    if (E__SUCCESS != result) {
        goto l_cleanup;
    }

    for (i = 0 ; i < smat->g_length ; ++i) {
        /* 2. Calculate Q_0 */
        q_0 = SUBMAT_SPMAT_LIST_calculate_q(smat, s_vector);

        /* 3. Computing DeltaQ for the move of each unmoved vertex */
        for (scanner = unmoved_scores->first ; NULL != scanner ; scanner = scanner->next) {
            /* Calculate score when moving k */
            k = scanner->index;
            s_vector[k] *= -1;
            q_score = SUBMAT_SPMAT_LIST_calculate_q(smat, s_vector);
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
    for (i = smat->g_length - 1 ; i > max_improvement_index ; --i) {
        s_vector[indices[i]] *= -1;
    }

    if (max_improvement_index == smat->g_length - 1) {
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
#endif

static
result_t
cluster_optimize_division_iteration2(submatrix_t *smat,
                                    double *s_vector,
                                    double *improve,
                                    int *indices,
                                    double *delta_q_out)
{
    result_t result = E__UNKNOWN;
    list_t *unmoved_scores = NULL;
    node_t *scanner = NULL;
    node_t *max_unmoved = NULL;
    int i = 0;
    int k = 0;
    double delta_q = 0.0;
    /* A negative value will never be the max since the last improvement is always 0 */
    double max_improvement_value = -1.0;
    int max_improvement_index = 0;

    /* 1. Initialize list */
    result = LIST_range(smat->g_length, &unmoved_scores);
    if (E__SUCCESS != result) {
        goto l_cleanup;
    }

    for (i = 0 ; i < smat->g_length ; ++i) {
        /* 3. Computing DeltaQ for the move of each unmoved vertex */
        for (scanner = unmoved_scores->first ; NULL != scanner ; scanner = scanner->next) {
            /* Calculate score when moving k */
            k = scanner->index;
            s_vector[k] *= -1;
            scanner->value = SUBMAT_SPMAT_LIST_calc_q_score(smat, s_vector, k);
            s_vector[k] *= -1;

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
    for (i = smat->g_length - 1 ; i > max_improvement_index ; --i) {
        s_vector[indices[i]] *= -1;
    }

    if (max_improvement_index == smat->g_length - 1) {
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


