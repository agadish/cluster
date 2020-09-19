/**
 * @file division.c
 * @purpose The output file of the spmat division
 */

/** Includes *************************************************************************************/
#include <stdio.h>

#include "results.h"
#include "matrix.h"
#include "division_file.h"
#include "spmat_list.h"


/** Structs **************************************************************************************/
struct division_file_s {
    FILE *file;
    int number_of_matrices;
};


/** Functions Declarations ***********************************************************************/
result_t
DIVISION_FILE_open(const char *path, division_file_t **division_file_out)
{
    result_t result = E__UNKNOWN;
    division_file_t *division_file = NULL;
    int result_fseek = -1;

    /* 0. Input validation */
    if ((NULL == path) || (NULL == division_file_out)) {
        result = E__NULL_ARGUMENT;
        goto l_cleanup;
    }

    /* 1. Allocate struct and open file */
    division_file = (division_file_t *)malloc(sizeof(*division_file));
    if (NULL == division_file) {
        result = E__MALLOC_ERROR;
        goto l_cleanup;
    }

    division_file->number_of_matrices = 0;
    division_file->file = fopen(path, "wb");
    if (NULL == division_file->file) {
        result = E__FOPEN_ERROR;
        goto l_cleanup;
    }

    /* 2. Initialize writing position */
    result_fseek = fseek(division_file->file,
                         sizeof(division_file->number_of_matrices),
                         SEEK_SET);
    if (-1 == result_fseek) {
        result = E__FSEEK_ERROR;
        goto l_cleanup;
    }

    /* Success */
    *division_file_out = division_file;

    result = E__SUCCESS;
l_cleanup:
    if (E__SUCCESS != result) {
        DIVISION_FILE_close(division_file);
        division_file = NULL;
    }

    return result;
}

result_t
DIVISION_FILE_write_matrix(division_file_t *division_file, const matrix_t *matrix)
{
    result_t result = E__UNKNOWN;
    size_t result_write = 0;

    /* 0. Input validation */
    if ((NULL == division_file) || (NULL == matrix)) {
        result = E__NULL_ARGUMENT;
        goto l_cleanup;
    }

    if (0 == matrix->n) {
        /* Nothing to write */
        result = E__SUCCESS;
        goto l_cleanup;
    }

    /* 1. Write n to file */
    result_write = fwrite(&matrix->n,
                          sizeof(matrix->n),
                          1,
                          division_file->file);
    if (1 != result_write) {
        result = E__FWRITE_ERROR;
        goto l_cleanup;
    }

    /* 2. Write neighbors to file */
    result = MATRIX_WRITE_NEIGHBORS(matrix, division_file->file);
    if (E__SUCCESS != result) {
        goto l_cleanup;
    }

    /* 3. Increase matrix count */
    ++division_file->number_of_matrices;

    /* Success */
    result = E__SUCCESS;
l_cleanup:

    return result;
}

result_t
DIVISION_FILE_finalize(division_file_t *division_file)
{
    result_t result = E__UNKNOWN;
    int result_fseek = -1;
    size_t result_write = 0;

    if (NULL == division_file) {
        result = E__NULL_ARGUMENT;
        goto l_cleanup;
    }

    /* 1. Write final size to file */
    /* 1.1. Seek to beginning */
    result_fseek = fseek(division_file->file, 0, SEEK_SET);
    if (-1 == result_fseek) {
        result = E__FSEEK_ERROR;
        goto l_cleanup;
    }

    /* 1.2. Write matrix count */
    result_write = fwrite((void *)&division_file->number_of_matrices,
                          sizeof(division_file->number_of_matrices),
                          1,
                          division_file->file);
    if (1 != result_write) {
        result = E__FWRITE_ERROR;
        goto l_cleanup;
    }

    /* Success */
    result = E__SUCCESS;
l_cleanup:

    return result;
}


void
DIVISION_FILE_close(division_file_t *division_file)
{
    if (NULL != division_file)
    {
        FCLOSE_SAFE(division_file->file);
        FREE_SAFE(division_file);
    }
}
