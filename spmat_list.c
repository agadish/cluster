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


/* Functions *************************************************************************************/
result_t
SPMAT_LIST_allocate(int n, matrix_t **mat_out)
{
    result_t result = E__UNKNOWN;

    matrix_t * mat = NULL;
    node_t **rows = NULL;
    size_t lists_array_size = 0;

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

/**
 *
 * @param relevant_value 1 or -1 accordingly to the values belong to the given reduced line
 */
static
result_t
spmat_list_reduce_line(const node_t *original_line,
                       const double * vector_s,
                       size_t vector_s_length,
                       double relevant_value,
                       node_t **reduced_list_out)
{
    result_t result = E__UNKNOWN;
    node_t *original_line_scanner = NULL;
    node_t *reduced_list_start = NULL; /* Will be returned */
    node_t *reduced_list_end = NULL; /* Will stick to its end */
    node_t *reduced_last_node = NULL; /* Temporary for newly created last node */
    size_t s_scanner = 0;
    int reduced_index = 0;

    /* 1. Check if original line is zeroes */
    if (NULL == original_line) {
        result = E__SUCCESS;
        goto l_cleanup;
    }

    original_line_scanner = original_line;

    /* Scan the whole s vector:
     * 1 1 -1 1 -1 -1 1 -1 -1 ..... 1 */
    for (s_scanner = 0 ;
        ((s_scanner < vector_s_length) && (NULL != original_line_scanner)) ;
        ++s_scanner)
    {
        /* We're looking at relevant_value = -1:
         *      V    V  V    V  V
         * 1 1 -1 1 -1 -1 1 -1 -1..... 1
         *
         * The original line contains:
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
        if (vector_s[s_scanner] == relevant_value) {
            /* original_line_scanner->index is ind 2 */
            if (original_line_scanner->index == s_scanner) {
                /* We encountered the ind 2 or ind 5 */
                /* Create a new node */
                result = LIST_NODE_create(original_line_scanner->value,
                                          reduced_index,
                                          &reduced_last_node);
                if (E__SUCCESS != result) {
                    goto l_cleanup;
                }

                /* Append the new node */
                if (NULL != reduced_list_end) {
                    /* First node in list */
                    reduced_list_start = reduced_last_node;
                } else {
                    /* Link with the list's end */
                    reduced_list_end->next = reduced_last_node;
                }

                reduced_list_end = reduced_last_node;
                
                /* Increase original_line_scanner */
                original_line_scanner = original_line_scanner->next;
            }

            ++reduced_index;
        }
    }

    /* Success */
    *reduced_list_out = reduced_list_start;
    result = E__SUCCESS;
l_cleanup:

    return result;
}

result_t
SPMAT_LIST_divide_matrix(const matrix_t *matrix,
                         const double * vector_s,
                         matrix_t **matrix1_out,
                         matrix_t **matrix2_out,
                         size_t ones_count)
{
    result_t result = E__UNKNOWN;
    matrix_t *matrix1 = NULL;
    matrix_t *matrix2 = NULL;
    size_t i = 0;
    size_t matrix1_i = 0;
    size_t matrix2_i = 0;
    node_t **rows = NULL;
    /* Original line in matrix */
    const node_t *scanned_row = NULL;
    /* Points to newly initialised row either of matrix1 or matrix2 */
    node_t *new_row = NULL;
    int new_s_value = 0;
    int scanned_s_value = 0;

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

    /* 1. Create matrixes as spmat lists */
    result = MATRIX_create_matrix(ones_count, MATRIX_TYPE_SPMAT_LIST, &matrix1);
    if (E__SUCCESS != result) {
        goto l_cleanup
    }

    result = MATRIX_create_matrix(matrix->n - ones_count, MATRIX_TYPE_SPMAT_LIST, &matrix2);
    if (E__SUCCESS != result) {
        goto l_cleanup
    }
    
    /* 2. Split matrix */
    rows = (node_t **)mat->private;
    for (i = 0 ; i < matrix->n ; ++i) {
        scanned_s_value = vector_s[i];

        /* 2.1. Find matching row */
        switch (scanned_s_value)
        {
        case 1:
            /* Row belongs to matrix1 */
            new_row = matrix1->rows[matrix1_i];
            ++matrix1_i;
            break;
        case -1:
            /* Row belongs to matrix2 */
            new_row = matrix2->rows[matrix2_i];
            ++matrix2_i;
            break;
        default:
            /* Error */
            result = E__INVALID_S_VECTOR;
            goto l_cleanup;
        }

        /* 2.2. Split matching row */
        for (scanned_row = rows[i] ; NULL != scanned_row ; scanned_row = scanned_row->next) {
            new_s_value = s_vector[scanned_row->index];
            if (scanned_s_value == 

        }
    }


    result = E__SUCCESS;
l_cleanup:

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

