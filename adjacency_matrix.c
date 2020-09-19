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
#include "config.h"
#include "matrix_raw.h"
#include "debug.h"
#include "spmat_list.h"

/* Functions Declarations ***********************************************************************/
/*
 *
 */
static result_t
adjacency_matrix_read_neighbors_line(FILE * file,
                                     adjacency_matrix_t *matrix,
                                     int line_index,
                                     double *tmp_buffer);


/* Functions ************************************************************************************/
static result_t
adjacency_matrix_read_neighbors_line(FILE *file,
                                     adjacency_matrix_t *adj_matrix,
                                     int line_index,
                                     double *tmp_neighbors_buffer)
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


    /* 2. Save number of neighbors */
    adj_matrix->neighbors[line_index] = number_of_edges;
    adj_matrix->M += number_of_edges;

    /* 3. Zero the neighbors buffer */
    memset(tmp_neighbors_buffer, 0, sizeof(*tmp_neighbors_buffer) * adj_matrix->matrix->n);

    /* 4. Read all edges */
    for (i = 0 ; i < number_of_edges ; ++i) {
        /* 4.1. Read current edge */
        result_fread = fread((void *)&current_edge, 1, sizeof(current_edge), file);
        if (1 > result_fread) {
            result = E__FREAD_ERROR;
            goto l_cleanup;
        }

        /* 4.2. Assign 1 to edge */
        tmp_neighbors_buffer[current_edge] = 1.0;
    }

    /* 4. Set row in matrix */
    result = MATRIX_ADD_ROW(adj_matrix->matrix, tmp_neighbors_buffer, line_index);
    if (E__SUCCESS != result) {
        goto l_cleanup;
    }
    

    result = E__SUCCESS;
l_cleanup:

    return result;
}

result_t
ADJACENCY_MATRIX_open(const char *path, adjacency_matrix_t **matrix_out)
{
    result_t result = E__UNKNOWN;
    adjacency_matrix_t *matrix = NULL;
    int matrix_n = 0 ;
    FILE *file = NULL;
    size_t result_fread = 0;
    double *tmp_neighbors_buffer = NULL;
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

    /* 4. Allocations */
    /* 4.1. Allocate matrix */
    result = MATRIX_create_matrix(matrix_n, MATRIX_TYPE_RAW, &matrix->matrix);
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
    /* 5.1. Allocate temporary neighbors buffer */
    tmp_neighbors_buffer = (double *)malloc(sizeof(*tmp_neighbors_buffer) * matrix_n);
    if (NULL == tmp_neighbors_buffer) {
        result = E__MALLOC_ERROR;
        goto l_cleanup;
    }

    for (i = 0 ; i < matrix_n ; ++i) {
        result = adjacency_matrix_read_neighbors_line(file, matrix, i, tmp_neighbors_buffer);
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
    FREE_SAFE(tmp_neighbors_buffer);

    return result;
}

void
ADJACENCY_MATRIX_free(adjacency_matrix_t *adjacency_matrix)
{
    if (NULL != adjacency_matrix) {
        MATRIX_FREE_SAFE(adjacency_matrix->matrix);
        FREE_SAFE(adjacency_matrix->neighbors);

        FREE_SAFE(adjacency_matrix);
    }
}


result_t
ADJACENCY_MATRIX_calculate_modularity(adjacency_matrix_t *adj,
                                      matrix_type_t mod_matrix_type,
                                      matrix_t **mod_matrix_out)
{
    result_t result = E__UNKNOWN;
    matrix_t *mod_matrix = NULL;
    int row = 0;
    int col = 0;
    double expected_edges = 0.0;
    double b_value = 0.0;
    double *current_row = NULL;
    double tmp = 0.0;

    /* 0. Input validation */
    if ((NULL == adj) || (NULL == mod_matrix_out)) {
        result = E__NULL_ARGUMENT;
        goto l_cleanup;
    }

    /* 1. Allocate modulation matrix */
    result = MATRIX_create_matrix(adj->matrix->n,
                                  mod_matrix_type,
                                  &mod_matrix);
    if (E__SUCCESS != result) {
        goto l_cleanup;
    }

    /* 2. Fill the modulation matrix */
    /* 2.1. Allocate buffer for current row */
    current_row = (double *)malloc(sizeof(*current_row) * adj->matrix->n);
    if (NULL == current_row) {
        result = E__MALLOC_ERROR;
        goto l_cleanup;
    }

    /* 2.2. Go over each row */
    for (row = 0 ; row < adj->matrix->n ; ++row) {
        /* 2.2.1. Calculate each column */
        MATRIX_GET_ROW(adj->matrix, current_row, row);

        for (col = 0 ; col < adj->matrix->n ; ++col) {
            expected_edges = ((double)(adj->neighbors[row]) * (double)(adj->neighbors[col])) / (double)(adj->M);
            /* DEBUG_PRINT("getting matrix at %d,%d (addr=%p)", */
            /*         row, */
            /*         col, */
            /*         (void *)&(MATRIX_RAW_AT(adj->matrix, row, col))); */
            b_value = tmp - expected_edges;
            current_row[col] -= b_value;
        }

        /* 2.2.2. ADd to mod matrix */
        result = MATRIX_ADD_ROW(mod_matrix, current_row, row);
        if (E__SUCCESS != result) {
            goto l_cleanup;
        }
    }

    /* 3. Initialise rows numbers */
    MATRIX_INITIALISE_ROW_NUMBERS(mod_matrix);

    /* Success */
    *mod_matrix_out = mod_matrix;

    result = E__SUCCESS;
l_cleanup:
    if (E__SUCCESS != result) {
        MATRIX_FREE_SAFE(mod_matrix);
    }
    FREE_SAFE(current_row);

    return result;
}
