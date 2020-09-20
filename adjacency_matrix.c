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
#include "debug.h"
#include "spmat_list.h"

/* Functions Declarations ***********************************************************************/
/*
 *
 */
static
result_t
adjacency_matrix_read_neighbors_line(FILE * file,
                                     adjacency_t *adj,
                                     matrix_t *matrix,
                                     int line_index,
                                     double *tmp_buffer);

static
void
adjacency_matrix_calculate_neighbors_div_M(adjacency_t *adj);



/* Functions ************************************************************************************/
static
result_t
adjacency_matrix_read_neighbors_line(FILE *file,
                                     adjacency_t *adj,
                                     matrix_t *matrix,
                                     int line_index,
                                     double *tmp_neighbors_buffer)
{
    result_t result = E__UNKNOWN;
    size_t result_fread = 0;
    int number_of_edges = 0;
    int i = 0;
    int current_edge = 0;

    /* 1. Read number of edges n */
    result_fread = fread((void *)&number_of_edges,
                         1,
                         sizeof(number_of_edges),
                         file);
    if (1 > result_fread) {
        result = E__FREAD_ERROR;
        goto l_cleanup;
    }

    /* 2. Save number of neighbors */
    adj->neighbors[line_index] = number_of_edges;
    adj->M += number_of_edges;

    /* 3. Zero the neighbors buffer */
    (void)memset(tmp_neighbors_buffer,
                 0,
                 sizeof(*tmp_neighbors_buffer) * adj->n);

    /* 4. Read all edges */
    for (i = 0 ; i < number_of_edges ; ++i) {
        /* 4.1. Read current edge */
        result_fread = fread((void *)&current_edge,
                             1,
                             sizeof(current_edge),
                             file);
        if (1 > result_fread) {
            result = E__FREAD_ERROR;
            goto l_cleanup;
        }

        /* 4.2. Assign 1 to edge */
        tmp_neighbors_buffer[current_edge] = 1.0;
    }

    /* 4. Set row in matrix */
    result = MATRIX_ADD_ROW(matrix,
                            tmp_neighbors_buffer,
                            line_index);
    if (E__SUCCESS != result) {
        goto l_cleanup;
    }
    

    result = E__SUCCESS;
l_cleanup:

    return result;
}

result_t
ADJACENCY_MATRIX_open(const char *path,
                      adjacency_t **adj_out,
                      matrix_t **matrix_out)
{
    result_t result = E__UNKNOWN;
    adjacency_t *adj = NULL;
    matrix_t *matrix = NULL;
    int matrix_n = 0 ;
    FILE *file = NULL;
    size_t result_fread = 0;
    double *tmp_neighbors_buffer = NULL;
    int i = 0;

    /* 1. Allocate adj adj */
    adj = (adjacency_t *)malloc(sizeof(*adj));
    if (NULL == adj) {
        result = E__MALLOC_ERROR;
        goto l_cleanup;
    }
    (void)memset(adj, 0, sizeof(*adj));

    /* 2. Open adj file */
    file = fopen(path, "rb");
    if (NULL == file) {
        result = E__FOPEN_ERROR;
        goto l_cleanup;
    }

    /* 3. Read adj n */
    result_fread = fread((void *)&matrix_n, 1, sizeof(matrix_n), file);
    if (1 > result_fread) {
        result = E__FREAD_ERROR;
        goto l_cleanup;
    }

    /* 4. Allocations */
    /* 4.1. Allocate adj */
    result = MATRIX_create_matrix(matrix_n,
                                  MATRIX_TYPE_SPMAT_LIST,
                                  &matrix);
    if (E__SUCCESS != result) {
        goto l_cleanup;
    }

    /* 4.2. Allocate neighbors array */
    adj->n = matrix_n;
    adj->neighbors = (int *)malloc(sizeof(*(adj->neighbors)) * matrix_n);
    if (NULL == adj->neighbors) {
        result = E__MALLOC_ERROR;
        goto l_cleanup;
    }

    adj->neighbors_div_M = (double *)malloc(sizeof(*(adj->neighbors_div_M)) *
                                         matrix_n);
    if (NULL == adj->neighbors_div_M) {
        result = E__MALLOC_ERROR;
        goto l_cleanup;
    }

    /* 5. Initialize each line */
    /* 5.1. Allocate temporary neighbors buffer */
    tmp_neighbors_buffer = (double *)malloc(sizeof(*tmp_neighbors_buffer) *
                                            matrix_n);
    if (NULL == tmp_neighbors_buffer) {
        result = E__MALLOC_ERROR;
        goto l_cleanup;
    }

    for (i = 0 ; i < matrix_n ; ++i) {
        result = adjacency_matrix_read_neighbors_line(file,
                                                      adj,
                                                      matrix,
                                                      i,
                                                      tmp_neighbors_buffer);
        if (E__SUCCESS != result) {
            goto l_cleanup;
        }
    }

    /* 6. Calculate neighbosr div M for optimization */
    adjacency_matrix_calculate_neighbors_div_M(adj);

    /* Success */
    *adj_out = adj;
    *matrix_out = matrix;

    result = E__SUCCESS;
l_cleanup:
    if (E__SUCCESS != result) {
        ADJACENCY_MATRIX_free(adj);
    }
    FCLOSE_SAFE(file);
    FREE_SAFE(tmp_neighbors_buffer);

    return result;
}

void
ADJACENCY_MATRIX_free(adjacency_t *adj)
{
    FREE_SAFE(adj->neighbors_div_M);
    FREE_SAFE(adj->neighbors);
    adj->n = 0;
    adj->M = 0;
    FREE_SAFE(adj);
}


static
void
adjacency_matrix_calculate_neighbors_div_M(adjacency_t *adj)
{
    int i = 0;

    for (i = 0 ; i < adj->n ; ++i) {
        adj->neighbors_div_M[i] = (double)adj->neighbors[i] / adj->M;
    }
}
