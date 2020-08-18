/**
 * @file lazy_matrix.c
 * @purpose 
 */

/* Includes **************************************************************************************/
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

#include "lazy_matrix.h"
#include "common.h"


/* Functions ************************************************************************************/
result_t
LAZY_MATRIX_open(const char * path,
                 lazy_matrix_t **matrix_out)
{
    result_t result = E__UNKNOWN;
    FILE * file = NULL;
    lazy_matrix_t * matrix = NULL;
    size_t result_fread = 0;

    if ((NULL == path) || (NULL == matrix_out)) {
        result = E__NULL_ARGUMENT;
        goto l_cleanup;
    }

    matrix = (lazy_matrix_t *)malloc(sizeof(*matrix));
    if (NULL == matrix) {
        result = E__MALLOC_ERROR;
        goto l_cleanup;
    }
    (void)memset(matrix, 0, sizeof(*matrix));

    file = fopen(path, "rb");
    if (NULL == file) {
        result = E__FOPEN_ERROR;
        goto l_cleanup;
    }
    matrix->file = file;

    /* Read matrix dimensions */
    result_fread = fread((void *)&matrix->header, 1, sizeof(matrix->header), file);
    if (sizeof(matrix->header) > result_fread) {
        result = E__FREAD_ERROR;
        goto l_cleanup;
    }

    matrix->current_line = (double *)malloc(sizeof(*matrix->current_line) * matrix->header.columns);
    if (NULL == matrix->current_line) {
        result = E__MALLOC_ERROR;
        goto l_cleanup;
    }

    /* Success */
    *matrix_out = matrix;

    result = E__SUCCESS;
l_cleanup:
    if (E__SUCCESS != result) {
        LAZY_MATRIX_close(matrix);
    }

    return result;
}

result_t
LAZY_MATRIX_rewind(lazy_matrix_t * matrix)
{
    result_t result = E__UNKNOWN;
    int result_fseek = -1;

    if (NULL == matrix) {
        result = E__NULL_ARGUMENT;
        goto l_cleanup;
    }

    result_fseek = fseek(matrix->file, sizeof(matrix->header), SEEK_SET);
    if (-1 == result_fseek) {
        result = E__FSEEK_ERROR;
        goto l_cleanup;
    }

    result = E__SUCCESS;
l_cleanup:

    return result;
}

void
LAZY_MATRIX_close(lazy_matrix_t * matrix)
{
    if (NULL != matrix) {
        FCLOSE_SAFE(matrix->file);
        FREE_SAFE(matrix->current_line);
    }

    FREE_SAFE(matrix);
}

result_t
LAZY_MATRIX_read_next_line(lazy_matrix_t * matrix)
{
    result_t result = E__UNKNOWN;
    size_t result_fread = 0;

    if (NULL == matrix) {
        result = E__NULL_ARGUMENT;
        goto l_cleanup;
    }

    result_fread = fread(matrix->current_line,
                         sizeof(double),
                         matrix->header.columns,
                         matrix->file);
    if ((size_t)matrix->header.columns > result_fread) {
        result = E__FREAD_ERROR;
        goto l_cleanup;
    }

    result = E__SUCCESS;
l_cleanup:

    return result;
}

result_t
LAZY_MATRIX_count_nonzero_values(lazy_matrix_t *matrix, int *count_out)
{
    result_t result = E__UNKNOWN;
    int nonzero_count = 0;
    int row = 0;
    int column = 0;

    for (row = 0 ; matrix->header.rows > row ; ++row) {
        result = LAZY_MATRIX_read_next_line(matrix);
        if (E__SUCCESS != result) {
            goto l_cleanup;
        }

        for (column = 0 ; matrix->header.columns > column ; ++column) {
            if (0.0 != matrix->current_line[column]) {
                ++nonzero_count;
            }
        }
    }

    /* Success */
    *count_out = nonzero_count;

    result = E__SUCCESS;
l_cleanup:

    return result;
}
