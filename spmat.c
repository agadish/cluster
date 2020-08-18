/* Includes **************************************************************************************/
#include <string.h>
#include <stddef.h>

#include "results.h"
#include "lazy_matrix.h"
#include "spmat.h"
#include "spmat_list.h"
#include "spmat_array.h"


/* Functions *************************************************************************************/
result_t
SPMAT_FACTORY_create_spmat(const char * input_matrix_path,
                           const char * spmat_type,
                           spmat **A_out)
{
    result_t result = E__UNKNOWN;
    lazy_matrix_t * input_matrix = NULL;
    spmat *mat = NULL;
    int rows = 0;
    int current_row = 0;
    int nnz = 0;

    /* 1. Open input matrix */
    result = LAZY_MATRIX_open(input_matrix_path, &input_matrix);
    if (E__SUCCESS != result) {
        goto l_cleanup;
    } 
    rows = input_matrix->header.rows;

    /* 2. Create the requested spmat according to the input dimensions */
    if (0 == strcmp(spmat_type, SPMAT_TYPE_LIST)) {
        /* SPMAT: List */
        result = SPMAT_LIST_allocate(rows, &mat);
        if (E__SUCCESS != result) {
            goto l_cleanup;
        }
    } else if (0 == strcmp(spmat_type, SPMAT_TYPE_ARRAY)) {
        /* SPMAT: Array */
        /* The allocator needs the count of non-zero values on the matrix */
        result = LAZY_MATRIX_count_nonzero_values(input_matrix, &nnz);
        if (E__SUCCESS != result) {
            goto l_cleanup;
        }

        result = SPMAT_ARRAY_allocate(rows, nnz, &mat);
        if (E__SUCCESS != result) {
            goto l_cleanup;
        }
    } else {
        result = E__UNKNOWN_SPMAT_IMPLEMNTATION;
        goto l_cleanup;
    }

    /* 3. Read input into the spmat */
    /* 3.1. Rewind to the first line */
    result = LAZY_MATRIX_rewind(input_matrix);
    if (E__SUCCESS != result) {
        goto l_cleanup;
    }

    /* 3.2. Read */
    for (current_row = 0 ; rows > current_row ; ++current_row) {
        result = LAZY_MATRIX_read_next_line(input_matrix);
        if (E__SUCCESS != result) {
            goto l_cleanup;
        }

        mat->add_row(mat, input_matrix->current_line, current_row);
    }

    /* Success */
    *A_out = mat;

    result = E__SUCCESS;
l_cleanup:
    if (E__SUCCESS != result) {
        /* On failure free the smpat */
        if (NULL != mat) {
            mat->free(mat);
            mat = NULL;
        }
    }
    LAZY_MATRIX_close(input_matrix);

    return result;
}

