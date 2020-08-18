/**
 * @file matrix.c
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

#include "matrix.h"
#include "common.h"



/* Functions ************************************************************************************/
result_t
MATRIX_allocate(int rows, int columns, matrix_t ** matrix_out)
{
    result_t result = E__UNKNOWN;
    matrix_t * matrix = NULL;
    size_t total_cells = 0;

    if (NULL == matrix_out) {
        result = E__NULL_ARGUMENT;
        goto l_cleanup;
    }

    if ((0 > rows) || (0 > columns)) {
        result = E__INVALID_SIZE;
        goto l_cleanup;
    }

    matrix = (matrix_t *)malloc(sizeof(*matrix));
    if (NULL == matrix) {
        result = E__MALLOC_ERROR;
        goto l_cleanup;
    }
    matrix->header.rows = rows;
    matrix->header.columns = columns;

    total_cells = rows * columns;
    matrix->data = (double *)malloc(sizeof(*matrix->data) * total_cells);
    if (NULL == matrix->data) {
        result = E__MALLOC_ERROR;
        goto l_cleanup;
    }

    /* Success */
    *matrix_out = matrix;

    result = E__SUCCESS;

l_cleanup:
    if (E__SUCCESS != result) {
        MATRIX_free(matrix);
    }

    return result;
}


result_t
MATRIX_open(const char * path, matrix_t **matrix_out)
{
    result_t result = E__UNKNOWN;
    matrix_t * matrix = NULL;
    FILE * file = NULL;
    size_t cells_count = 0;
    size_t result_fread = 0;

    /* 1. Allocate matrix */
    matrix = (matrix_t *)malloc(sizeof(*matrix));
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

    /* 3. Read matrix dimensions */
    result_fread = fread((void *)&matrix->header, 1, sizeof(matrix->header), file);
    if (sizeof(matrix->header) > result_fread) {
        result = E__FREAD_ERROR;
        goto l_cleanup;
    }

    /* 4. Read matrix data */
    /* 4.1. Allocate */
    cells_count = MATRIX_COUNT(matrix);
    matrix->data = (double *)malloc(sizeof(*matrix->data) * cells_count);
    if (NULL == matrix->data) {
        result = E__MALLOC_ERROR;
        goto l_cleanup;
    }

    /* 4.2. Read */
    result_fread = fread(matrix->data, sizeof(*matrix->data), cells_count, file);
    if ((size_t)cells_count > result_fread) {
        result = E__FREAD_ERROR;
        goto l_cleanup;
    }

    /* Success */
    *matrix_out = matrix;

    result = E__SUCCESS;
l_cleanup:
    if (E__SUCCESS != result) {
        MATRIX_free(matrix);
    }
    FCLOSE_SAFE(file);

    return result;
}

void
MATRIX_free(matrix_t * matrix)
{
    if (NULL != matrix) {
        FREE_SAFE(matrix->data);
        FREE_SAFE(matrix);
    }
}

result_t
MATRIX_div(matrix_t * vector, double value)
{
    result_t result = E__UNKNOWN;
    double * buffer = vector->data;
    double * buffer_end = NULL;

    if (NULL == vector) {
        result = E__NULL_ARGUMENT;
        goto l_cleanup;
    }
    if (0.0 == value) {
        result = E__ZERO_DIV_ERROR;
        goto l_cleanup;
    }

    buffer_end = vector->data + MATRIX_COUNT(vector);
    for ( ; buffer < buffer_end - 1; buffer += 2) {
        buffer[0] /= value;
        buffer[1] /= value;
    }

    if (buffer < buffer_end) {
        buffer[0] /= value;
    }

    result = E__SUCCESS;
l_cleanup:

    return result;
}

result_t
MATRIX_normalize(matrix_t *matrix)
{
    result_t result = E__UNKNOWN;
    double magnitude = 0.0;

    if (NULL == matrix) {
        result = E__NULL_ARGUMENT;
        goto l_cleanup;
    }

    magnitude = MATRIX_calculate_magnitude(matrix);

    result = MATRIX_div(matrix, magnitude);
    if (E__SUCCESS != result) {
        goto l_cleanup;
    }

    result = E__SUCCESS;
l_cleanup:

    return result;
}

double
MATRIX_calculate_magnitude(matrix_t *matrix)
{
    double * current = matrix->data;
    const double * end = matrix->data + MATRIX_COUNT(matrix);
    double magnitude = 0.0;

    for ( ; end - 1 > current ; current +=2) {
        magnitude += ((current[0] * current[0]) + (current[1] * current[1]));
    }

    for ( ; end > current ; ++current) {
        magnitude += (current[0] * current[0]);
    }

    magnitude = sqrt(magnitude);

    return magnitude;
}

bool_t
MATRIX_is_close(matrix_t * prev_vector, matrix_t * curr_vector, double epsilon){

    bool_t result = TRUE;
    int row = 0;

    for (row = 0 ; row < prev_vector->header.rows; ++row) {
        if (fabs((prev_vector->data[row]) - (curr_vector->data[row])) >= epsilon) {
            result = FALSE;
            break;
        }
    }
    return result;
}

result_t
MATRIX_random_vector(const int length, matrix_t ** vector_out)
{
    result_t result = E__UNKNOWN;
    matrix_t * vector = NULL;
    int row = 0;
    int random_int = 0;

    if ((NULL == vector_out) || (0 > length)) {
        result = E__NULL_ARGUMENT;
        goto l_cleanup;
    }

    srand(time(NULL));

    result = MATRIX_allocate(length, 1, &vector);
    if (E__SUCCESS != result) {
        goto l_cleanup;
    }

    for (row = 0 ; row < length ; ++row) {
        random_int = (rand() % 1000);
        MATRIX_AT(vector, row, 0) = (double)random_int;
    }

    *vector_out = vector;

    result = E__SUCCESS;
l_cleanup:
    if (E__SUCCESS != result) {
        FREE_SAFE(vector);
    }

    return result;
}
