/**
 * @file adjacency_matrix.c
 * @purpose 
 */

/* Includes **************************************************************************************/
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>

#include "matrix.h"
#include "common.h"
#include "adjacency_matrix.h"


/* Functions Declarations ***********************************************************************/
/*
 *
 */
static result_t
adjacency_matrix_read_neighbors_line(FILE * file, adjacency_matrix_t *matrix, int line_index);


/* Functions ************************************************************************************/
static result_t
adjacency_matrix_read_neighbors_line(FILE * file,
                                     adjacency_matrix_t *adj_matrix,
                                     int line_index)
{
    result_t result = E__UNKNOWN;
    size_t result_fread = 0;
    int number_of_edges = 0;
    int i = 0;
    int current_edge = 0;

    /* 1. Read number of edges n */
    result_fread = fread((void *)&number_of_edges, 1, sizeof(number_of_edges), file);
    if (1 > result_fread) {
        result = E__FREAD_ERROR;
        goto l_cleanup;
    }

    /* 2. Define the (i,i) as 0 */
    MATRIX_AT(adj_matrix->matrix, line_index, line_index) = 0;

    /* 3. Read all edges */
    for (i = 0 ; i < number_of_edges ; ++i) {
        /* 3.1. Read current edge */
        result_fread = fread((void *)&current_edge, 1, sizeof(current_edge), file);
        if (1 > result_fread) {
            result = E__FREAD_ERROR;
            goto l_cleanup;
        }

        /* 3.2. Save number of neighbors */
        adj_matrix->neighbors[line_index] = number_of_edges;
        adj_matrix->M += number_of_edges;

        /* 3.3. Assign 1 to edge */
        MATRIX_AT(adj_matrix->matrix, line_index, current_edge) = 1.0;
    }

    result = E__SUCCESS;
l_cleanup:

    return result;
}

result_t
ADJACENCY_MATRIX_open(const char * path, adjacency_matrix_t **matrix_out)
{
    result_t result = E__UNKNOWN;
    adjacency_matrix_t * matrix = NULL;
    int matrix_n = 0 ;
    FILE * file = NULL;
    size_t result_fread = 0;
    int i = 0;

    /* 1. Allocate adj matrix */
    matrix = (adjacency_matrix_t *)malloc(sizeof(*matrix));
    if (NULL == matrix) {
        result = E__MALLOC_ERROR;
        goto l_cleanup;
    }
    (void)memset(matrix, 0, sizeof(*matrix));

    /* 2. Open matrix file */
    file = fopen(path, "rb");
    if (NULL == file) {
        result = E__FOPEN_ERROR;
        goto l_cleanup;
    }

    /* 3. Read matrix n */
    result_fread = fread((void *)&matrix_n, 1, sizeof(matrix_n), file);
    if (1 > result_fread) {
        result = E__FREAD_ERROR;
        goto l_cleanup;
    }

    /* 4. Allocations matrix */
    /* 4.1. Allocate matrix */
    result = MATRIX_allocate(matrix_n, matrix_n, &matrix->matrix);
    if (E__SUCCESS != result) {
        goto l_cleanup;
    }

    /* 4.2. Allocate neighbors array */
    matrix->neighbors = (int *)malloc(sizeof(*(matrix->neighbors)) * matrix_n);
    if (NULL == matrix->neighbors) {
        result = E__MALLOC_ERROR;
        goto l_cleanup;
    }

    /* 5. Initialize each line */
    for (i = 0 ; i < matrix_n ; ++i) {
        result = adjacency_matrix_read_neighbors_line(file, matrix, i);
        if (E__SUCCESS != result) {
            goto l_cleanup;
        }
    }

    /* Success */
    *matrix_out = matrix;

    result = E__SUCCESS;
l_cleanup:
    if (E__SUCCESS != result) {
        ADJACENCY_MATRIX_free(matrix);
    }
    FCLOSE_SAFE(file);

    return result;
}

void
ADJACENCY_MATRIX_free(adjacency_matrix_t *adjacency_matrix)
{
    if (NULL != adjacency_matrix) {
        MATRIX_free(adjacency_matrix->matrix);
        adjacency_matrix->matrix = NULL;

        FREE_SAFE(adjacency_matrix->neighbors);
        FREE_SAFE(adjacency_matrix);
    }
}


result_t
ADJACENCY_MATRIX_calculate_modularity(adjacency_matrix_t *adj, matrix_t **mod_matrix_out)
{
    result_t result = E__UNKNOWN;
    matrix_t *mod_matrix = NULL;
    int row = 0;
    int col = 0;
    double expected_edges = 0.0;
    double b_value = 0.0;

    /* 0. Input validation */
    if ((NULL == adj) || (NULL == mod_matrix_out)) {
        result = E__NULL_ARGUMENT;
        goto l_cleanup;
    }

    /* 1. Allocate modulation matrix */
    result = MATRIX_allocate(adj->matrix->header.rows, adj->matrix->header.columns, &mod_matrix);
    if (E__SUCCESS != result) {
        goto l_cleanup;
    }

    /* 2. Calculate each cell */
    for (row = 0 ; row < adj->matrix->header.rows ; ++row) {
        for (col = 0 ; col < adj->matrix->header.columns ; ++col) {
            expected_edges = (adj->neighbors[row] * adj->neighbors[col]) / adj->M;
            b_value = MATRIX_AT(adj->matrix, row, col) - expected_edges;
            MATRIX_AT(mod_matrix, row, col) = b_value;
        }
    }

    /* Success */
    *mod_matrix_out = mod_matrix;

    result = E__SUCCESS;
l_cleanup:
    if (E__SUCCESS != result) {
        MATRIX_free(mod_matrix);
        mod_matrix = NULL;
    }

    return result;
}
