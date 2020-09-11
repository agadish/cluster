/**
 * @file main.c
 * @purpose Cluster main logic
 */

#include <stdio.h>

#include "results.h"
#include "matrix.h"
#include "adjacency_matrix.h"
#include "common.h"
#include "cluster.h"


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
                                                   MATRIX_TYPE_RAW,
                                                   &mod_matrix);
    if (E__SUCCESS != result) {
        goto l_cleanup;
    }


    /* 4. Divide */
    result = CLUSTER_divide(mod_matrix, &group1, &group2);
    if (E__SUCCESS != result) {
    }


    result = E__SUCCESS;
l_cleanup:
    if (NULL != adj_matrix) {
        ADJACENCY_MATRIX_free(adj_matrix);
        adj_matrix = NULL;
    }

    MATRIX_FREE_SAFE(mod_matrix);

    return (int)result;
}

