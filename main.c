/**
 *
 * @file main.c
 * @purpose Cluster main logic
 */

#include <stdio.h>

#include "results.h"
#include "matrix.h"
#include "adjacency_matrix.h"
#include "common.h"
#include "spmat_list.h"
#include "cluster.h"
#include "config.h"


/* Enums *****************************************************************************************/
enum clsuter_args_e {
    ARG_PROGRAM_NAME,
    ARG_INPUT_ADJACENCY,
    ARG_OUTPUT_GRAPH,
    ARG_COUNT
};


/* Functions *************************************************************************************/
int main(int argc, const char * argv[])
{
    result_t result = E__UNKNOWN;
    adjacency_matrix_t *adj_matrix = NULL;
    matrix_t *mod_matrix = NULL;
    matrix_t *group1 = NULL;
    matrix_t *group2 = NULL;

    /* 1. Input validation */
    if (ARG_COUNT != argc) {
        (void)fprintf(stderr, "Usage: %s INPUT_ADJACENCY OUTPUT_MATRICES\n", argv[0]);

        result = E__INVALID_CMDLINE_ARGS;
        goto l_cleanup;
    }

    /* 2. Open adjacency matrix */
    result = ADJACENCY_MATRIX_open(argv[ARG_INPUT_ADJACENCY], &adj_matrix);
    if (E__SUCCESS != result) {
        goto l_cleanup;
    }

    /* 3. Calculate modularity matrix */
    result = ADJACENCY_MATRIX_calculate_modularity(adj_matrix,
                                                   MOD_MATRIX_TYPE,
                                                   &mod_matrix);
    if (E__SUCCESS != result) {
        goto l_cleanup;
    }


    SPMAT_LIST_print("mod_matrix", mod_matrix);

    /* 4. Divide */
    result = CLUSTER_divide(mod_matrix, &group1, &group2);
    if (E__SUCCESS != result) {
    }
    SPMAT_LIST_print("group1", group1);
    SPMAT_LIST_print("group2", group2);


    result = E__SUCCESS;
l_cleanup:
    if (NULL != adj_matrix) {
        ADJACENCY_MATRIX_free(adj_matrix);
        adj_matrix = NULL;
    }

    MATRIX_FREE_SAFE(mod_matrix);
    MATRIX_FREE_SAFE(group1);
    MATRIX_FREE_SAFE(group2);

    return (int)result;
}

