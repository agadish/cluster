/*
 * @file spmat_list.c
 * @purpose Sparse matrix implemented using linked lists
 */
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include "list_node.h"
#include "results.h"
#include "matrix.h"
#include "spmat_list.h"
#include "common.h"


/* Structs ***************************************************************************************/
typedef struct spmat_list_s {
    /* n-sized array of linked lists */
    node_t **rows;
} spmat_list_t;


/* Macros ****************************************************************************************/
#define SPMAT_LIST_IS_VALID_ROW_INDEX(A, i) (((A)->n > i) && (0 <= i))
#define GET_ROW(list, row_index) (((node_t **)((list)->private))[(row_index)])


/* Functions Declarations ************************************************************************/
static
void
spmat_list_add_row(matrix_t *A, const double *row, int i);

static
void
spmat_list_free(matrix_t *A);

static
void
spmat_list_mult(const matrix_t *A, const double *v, double *result);

/**
 * @remark vector_s must be valid with given length, and values 1 or -1
 */
static
result_t
spmat_list_create_s_indexes(const double * vector_s,
                            int length,
                            int **s_indexes_out,
                            int *matrix1_n_out);
static
result_t
spmat_list_reduce_row(const node_t *original_row,
                       const double * vector_s,
                       double relevant_vector_s_value,
                       const int * s_indexes,
                       node_t **row_out);


/* Functions *************************************************************************************/
result_t
SPMAT_LIST_allocate(int n, matrix_t **mat_out)
{
    result_t result = E__UNKNOWN;

    matrix_t * mat = NULL;
    node_t **rows = NULL;
    int lists_array_size = 0;

    /* 0. Input validation */
    if (NULL == mat) {
        result = E__NULL_ARGUMENT;
        goto l_cleanup;
    }

    if (0 > n) {
        result = E__INVALID_SIZE;
        goto l_cleanup;
    }

    /* 1. Allocate matrix_t */
    mat = (matrix_t *)malloc(sizeof(*mat));
    if (NULL == mat) {
        result = E__MALLOC_ERROR;
        goto l_cleanup;
    }

    /* 2. Initialize */
    (void)memset(&mat, 0, sizeof(mat));
    mat->n = n;
    mat->add_row = spmat_list_add_row;
    mat->free = spmat_list_free;
    mat->mult = spmat_list_mult;
    /* TODO: Implement */
    mat->rmult = NULL;

    mat->private = NULL;

    /* 3. Rows array */
    lists_array_size = n * sizeof(node_t *);
    rows = malloc(lists_array_size);
    if (NULL == rows) {
        result = E__MALLOC_ERROR;
        goto l_cleanup;
    }
    (void)memset(rows, 0, lists_array_size);
    mat->private = (void *)rows;

    /* Success */
    *mat_out = mat;

    result = E__SUCCESS;
l_cleanup:
    if (E__SUCCESS != result) {
        spmat_list_free(mat);
        mat = NULL;
    }

    return result;
}

static
void
spmat_list_free(matrix_t *mat)
{
    node_t **rows = NULL;
    int i = 0;

    if (NULL != mat) {
        rows = (node_t **)mat->private;
        if (NULL != mat->private) {
            for (i = 0 ; i < mat->n ; ++i) {
                LIST_NODE_destroy(rows[i]);
                rows[i] = NULL;
            }
            FREE_SAFE(mat->private);
            rows = NULL;
            mat->private = NULL;
        }
        FREE_SAFE(mat);
    }
}

void
spmat_list_add_row(matrix_t *mat, const double *row, int i)
{
    result_t result = E__UNKNOWN;

    node_t *new_row_end = NULL;
    node_t *new_row_current = NULL;
    int column = 0;
    node_t **rows = NULL;

    /* 1. Input validation */
    if ((NULL == mat) || (NULL == mat->private) || (NULL == row)) {
        result = E__NULL_ARGUMENT;
        goto l_cleanup;
    }

    if (!SPMAT_LIST_IS_VALID_ROW_INDEX(mat, i)) {
        result = E__INVALID_ROW_INDEX;
        goto l_cleanup;
    }

    rows = (node_t **)mat->private;

    /* 2. Destory current row value */
    if (NULL != rows[i]) {
        LIST_NODE_destroy(rows[i]);
        rows[i] = NULL;
    }

    /* 3. Create a new row */
    for (column = 0 ; mat->n > column ; ++column) {
        if (0 != row[column]) {
            /* 3.1. Found a non-zero value, create a node for it */
            result = LIST_NODE_create(row[column], column, &new_row_current);
            if (E__SUCCESS != result) {
                goto l_cleanup;
            }
            
            /* 3.2.1. Is this the first node in the row? */
            if (NULL == new_row_end) {
                rows[i] = new_row_current;
            } else {
                new_row_end->next = new_row_current;
            }
            new_row_end = new_row_current;
        }
    }

    result = E__SUCCESS;
l_cleanup:
    if (E__SUCCESS != result) {
        LIST_NODE_destroy(rows[i]);
        rows[i] = NULL;
    }

    return;
} 

static
void
spmat_list_mult(const matrix_t *mat, const double *v, double *multiplication_result)
{
    result_t result = E__UNKNOWN;
    node_t **rows = NULL;
    int i = 0;

    if ((NULL == mat) || (NULL == v) || (NULL == multiplication_result)) {
        /* Remark: this function signature shouldn't be void */
        result = E__NULL_ARGUMENT;
        goto l_cleanup;
    }

    rows = (node_t **)mat->private;
    for (i = 0 ; i < mat->n ; ++i) {
        multiplication_result[i] = LIST_NODE_scalar_multiply(rows[i], v);
    }

    result = E__SUCCESS;
l_cleanup:
    UNUSED_ARG(result);

    return;
}

static
result_t
spmat_list_create_s_indexes(const double * vector_s,
                            int length,
                            int **s_indexes_out,
                            int *matrix1_n_out)
{
    result_t result = E__UNKNOWN;
    int *s_indexes = NULL;
    int i = 0;
    int index_1 = 0;
    int index_2 = 0;

    s_indexes = (int *)malloc(sizeof(*s_indexes) * length);
    if (NULL == s_indexes) {
        result = E__MALLOC_ERROR;
        goto l_cleanup;
    }

    /* Go over the s-vector */
    for (i = 0 ; i < length ; ++i) {
        if (1 == vector_s[i]) {
            s_indexes[i] = index_1;
            ++index_1;
        } else if (-1 == vector_s[i]) {
            s_indexes[i] = index_2;
            ++index_2;
        } 
    }

    /* Success */
    *s_indexes_out = s_indexes;
    *matrix1_n_out = index_1;
    result = E__SUCCESS;
l_cleanup:
    if (E__SUCCESS != result) {
        FREE_SAFE(s_indexes);
    }

    return result;
}

/**
 *
 * @param relevant_value 1 or -1 accordingly to the values belong to the given reduced row
 */
static
result_t
spmat_list_reduce_row(const node_t *original_row,
                       const double * vector_s,
                       double relevant_vector_s_value,
                       const int * s_indexes,
                       node_t **row_out)
{
    result_t result = E__UNKNOWN;
    const node_t *original_row_scanner = NULL;
    node_t *row = NULL;
    node_t *row_end = NULL;
    double scanned_s_value = 0.0; /* 1 or -1 */
    int scanned_index = 0.0; /* 0 ... n */

    /* 1. Check if original row is zeroes */
    if (NULL == original_row) {
        result = E__SUCCESS;
        goto l_cleanup;
    }

    /* Go over nodes in the given row */
    for (original_row_scanner = original_row ;
            NULL != original_row_scanner ;
            original_row_scanner = original_row_scanner->next) {
        /* We're looking at relevant_value = -1:
         *      V    V  V    V  V
         * 1 1 -1 1 -1 -1 1 -1 -1..... 1
         *
         * The original row contains:
         * * * 5 *  0  *  *  0  2 ... *
         *
         * The linked list is:
         *        [ BLAH  ]    [       ]     [ BLAH  ]     [       ]
         *  HEAD: [ val X ] -->[ val 5 ] --> [ val X ] --> [ val 2 ] --> NULL
         *        [ ind X ]    [ ind 2 ]     [ ind X ]     [ ind 8 ]    
         *
         * We will create:
         *        [       ]     [       ]
         *  HEAD: [ val 5 ] --> [ val 2 ] --> NULL
         *        [ ind 0 ]     [ ind 1 ]
         *
         */
        scanned_s_value = vector_s[original_row_scanner->index];
        if (scanned_s_value != relevant_vector_s_value) {
            continue;
        }

        scanned_index = s_indexes[original_row_scanner->index];
        result = LIST_NODE_append(&row_end, original_row_scanner->value, scanned_index);
        if (E__SUCCESS != result) {
            goto l_cleanup;
        }

        if (NULL == row) {
            row = row_end;
        }
    }


    /* Success */
    *row_out = row;
    result = E__SUCCESS;
l_cleanup:
    if (E__SUCCESS != result) {
        LIST_NODE_destroy(row);
    }

    return result;
}

result_t
SPMAT_LIST_divide_matrix(const matrix_t *matrix,
                         const double * vector_s,
                         matrix_t **matrix1_out,
                         matrix_t **matrix2_out)
{
    result_t result = E__UNKNOWN;
    matrix_t *matrix1 = NULL;
    matrix_t *matrix2 = NULL;
    int i = 0;
    int matrix1_n = 0;
    double scanned_s_value = 0.0;
    int *s_indexes = NULL;
    int scanned_s_index = 0;
    node_t **relevant_row_pointer = NULL;

    /* 0. Input validation */
    /* Null arguments */
    if ((NULL == matrix) ||
            (NULL == vector_s) ||
            (NULL == matrix1_out) ||
            (NULL == matrix2_out)) {
        result = E__NULL_ARGUMENT;
        goto l_cleanup;
    }

    /* Input matrix must be spmat list */
    if (MATRIX_TYPE_SPMAT_LIST != matrix->type) {
        result = E__INVALID_MATRIX_TYPE;
        goto l_cleanup;
    }

    /* 1. Create s-indexes vector, get matrix1's length */
    result = spmat_list_create_s_indexes(vector_s, matrix->n, &s_indexes, &matrix1_n);
    if (E__SUCCESS != result) {
        goto l_cleanup;
    }

    /* 2. Create matrixes as spmat lists */
    /* 2.1. Matrix 1 */
    result = MATRIX_create_matrix(matrix1_n, MATRIX_TYPE_SPMAT_LIST, &matrix1);
    if (E__SUCCESS != result) {
        goto l_cleanup;
    }

    /* 2.2. Matrix 2 */
    result = MATRIX_create_matrix(matrix->n - matrix1_n, MATRIX_TYPE_SPMAT_LIST, &matrix2);
    if (E__SUCCESS != result) {
        goto l_cleanup;
    }
    
    /* 3. Go over each row of the original matrix */
    for (i = 0 ; i < matrix->n ; ++i) {
        /* 3.1. Get the relevant s-value (1=matrix1, -1=matrix2) */
        scanned_s_value = vector_s[i]; /* 1 or -1 */
        /* 3.2. Get the row index within the selected matrix */
        scanned_s_index = s_indexes[i];

        /* 3.3. Get the actual row to be added */
        if (1.0 == scanned_s_value) {
            relevant_row_pointer = &GET_ROW(matrix1, scanned_s_index);
        } else if (-1.0 == scanned_s_value) {
            relevant_row_pointer = &GET_ROW(matrix2, scanned_s_index);
        } else {
            result = E__INVALID_S_VECTOR;
            goto l_cleanup;
        }

        /* 3.4. Add the filtered values in the row */
        result = spmat_list_reduce_row(GET_ROW(matrix, i),
                                        vector_s,
                                        scanned_s_value,
                                        s_indexes,
                                        relevant_row_pointer);
        if (E__SUCCESS != result) {
            goto l_cleanup;
        }
    }

    /* Success */
    *matrix1_out = matrix1;
    *matrix2_out = matrix2;

    result = E__SUCCESS;
l_cleanup:

    FREE_SAFE(s_indexes);

    if (E__SUCCESS != result) {
        if (NULL != matrix1) {
            MATRIX_FREE(matrix1);
            matrix1 = NULL;
        }
        if (NULL != matrix2) {
            MATRIX_FREE(matrix2);
            matrix2 = NULL;
        }
    }

    return result;
}

