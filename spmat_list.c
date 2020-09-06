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

