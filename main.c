/**
 *
 * @file main.c
 * @purpose Cluster main logic
 */

#include <stdio.h>
#include <time.h>

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
    matrix_t *group1 = NULL;
    matrix_t *group2 = NULL;
    division_file_t *division_file = NULL;
    clock_t start = 0;
    clock_t end = 0;

    /* 1. Input validation */
    if (ARG_COUNT != argc) {
        (void)fprintf(stderr, "Usage: %s INPUT_ADJACENCY OUTPUT_MATRICES\n", argv[0]);

        result = E__INVALID_CMDLINE_ARGS;
        goto l_cleanup;
    }

    start = clock();

    /* 2. Open adjacency matrix */
    result = ADJACENCY_MATRIX_open(argv[ARG_INPUT_ADJACENCY], &adj_matrix);
    if (E__SUCCESS != result) {
        goto l_cleanup;
    }


    /* 3. Create output file */
    result = DIVISION_FILE_open(argv[ARG_OUTPUT_GRAPH], &division_file);
    if (E__SUCCESS != result) {
        goto l_cleanup;
    }

    /* 5. Divide. Note: mod_matrix is freed by divide */
    result = CLUSTER_divide_repeatedly(adj_matrix, division_file);
    if (E__SUCCESS != result) {
        goto l_cleanup;
    }

    /* 6. Finalize output file */
    result = DIVISION_FILE_finalize(division_file);
    if (E__SUCCESS != result) {
        goto l_cleanup;
    }

    end = clock();
    printf("OUR PROGRAM: took %f sec\n", (double)(end - start) / CLOCKS_PER_SEC);


    result = E__SUCCESS;
l_cleanup:
    if (NULL != adj_matrix) {
        ADJACENCY_MATRIX_free(adj_matrix);
        adj_matrix = NULL;
    }

    DIVISION_FILE_close(division_file);
    division_file = NULL;

    MATRIX_FREE_SAFE(group1);
    MATRIX_FREE_SAFE(group2);

    return (int)result;
}

