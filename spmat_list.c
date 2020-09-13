/*
 * @file spmat_list.c
 * @purpose Sparse matrix implemented using linked lists
 */

/* Includes **************************************************************************************/
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "list_node.h"
#include "results.h"
#include "matrix.h"
#include "spmat_list.h"
#include "common.h"
#include "debug.h"


/* Structs ***************************************************************************************/
/* spmat_row_t * matrix->private: n-sized array of spmat_row_t */
typedef struct spmat_row_s {
    /* Begin of row's linked list */
    node_t *begin;
    double sum;
} spmat_row_t;

/* spmat_row_t * matrix->private: n-sized array of spmat_row_t */
typedef struct spmat_data_s {
    /* Begin of row's linked list */
    spmat_row_t *rows;
    double *columns_sum;
} spmat_data_t;

/* Macros ****************************************************************************************/
#define GET_SPMAT_DATA(matrix) ((spmat_data_t *)((matrix)->private))

#define GET_ROWS_ARRAY(matrix) (GET_SPMAT_DATA(matrix)->rows)

#define GET_ROW(matrix, row_index) (GET_ROWS_ARRAY(matrix)[row_index])


/* Functions Declarations ************************************************************************/
static
result_t
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
spmat_list_reduce_row(const spmat_row_t *original_row,
                       const double * vector_s,
                       double relevant_vector_s_value,
                       const int * s_indexes,
                       spmat_row_t *row_out);

/* Functions *************************************************************************************/
result_t
SPMAT_LIST_allocate(int n, matrix_t **mat_out)
{
    result_t result = E__UNKNOWN;
    matrix_t *mat = NULL;
    spmat_data_t *spmat_data = NULL;
    spmat_row_t *rows_array = NULL;
    double *columns_sum = NULL;
    int rows_array_size = 0;

    /* 0. Input validation */
    if (NULL == mat_out) {
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
    (void)memset(mat, 0, sizeof(*mat));
    mat->add_row = spmat_list_add_row;
    mat->free = spmat_list_free;
    mat->mult = spmat_list_mult;
    mat->private = NULL;

    /* 3. spmat struct */
    spmat_data = (spmat_data_t *)malloc(sizeof(*spmat_data));
    if (NULL == spmat_data) {
        result = E__MALLOC_ERROR;
        goto l_cleanup;
    }
    (void)memset(spmat_data, 0, sizeof(*spmat_data));

    /* 4. Rows array */
    rows_array_size = n * sizeof(*rows_array);
    rows_array = (spmat_row_t *)malloc(rows_array_size);
    if (NULL == rows_array) {
        result = E__MALLOC_ERROR;
        goto l_cleanup;
    }
    (void)memset(rows_array, 0, rows_array_size);
    spmat_data->rows = rows_array;

    /* 5. Columns sums */
    columns_sum = (double *)malloc(sizeof(*columns_sum) * n);
    if (NULL == columns_sum) {
        result = E__MALLOC_ERROR;
        goto l_cleanup;
    }
    (void)memset(columns_sum, 0, sizeof(*columns_sum) * n);
    spmat_data->columns_sum = columns_sum;

    mat->private = (void *)spmat_data;

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
    spmat_data_t *spmat_data = NULL;
    spmat_row_t *rows_array = NULL;
    int i = 0;

    if (NULL != mat) {
        spmat_data = GET_SPMAT_DATA(mat);
        if (NULL != spmat_data) {
            rows_array = spmat_data->rows;
            if (NULL != rows_array) {
                for (i = 0 ; i < mat->n ; ++i) {
                    LIST_NODE_destroy(rows_array[i].begin);
                    rows_array[i].begin = NULL;
                }

                FREE_SAFE(spmat_data->rows);
                rows_array = NULL;
            }
            
            FREE_SAFE(spmat_data->columns_sum);
            FREE_SAFE(spmat_data);
            mat->private = NULL;
        }
        FREE_SAFE(mat);
    }
}

result_t
spmat_list_add_row(matrix_t *mat, const double *values, int row_index)
{
    result_t result = E__UNKNOWN;
    spmat_data_t *matrix_data = NULL;
    spmat_row_t *row = NULL;
    node_t *prev_node = NULL;
    node_t *next_node = NULL;
    int col = 0;
    bool_t is_first_in_line = FALSE;

    /* 0. Input validation */
    if ((NULL == mat) || (NULL == mat->private) || (NULL == values)) {
        result = E__NULL_ARGUMENT;
        goto l_cleanup;
    }

    if (!MATRIX_IS_VALID_ROW_INDEX(mat, row_index)) {
        result = E__INVALID_ROW_INDEX;
        goto l_cleanup;
    }

    /* 1. Initialisations */
    matrix_data = GET_SPMAT_DATA(mat);
    row = &GET_ROW(mat, row_index);
    prev_node = NULL;
    next_node = row->begin;
    is_first_in_line = TRUE;

#if 0
    DEBUG_PRINT("called with row %d, values: ", row_index);
    for (col =  0 ; mat->n > col ; ++col) {
        printf("%f ", values[col]);
    }
    printf("\n");
#endif /* 0 */

    /* 2. Add row to matrix */
    for (col = 0 ; mat->n > col ; ++col) {
        if (0 != values[col]) {
            /* 2.1. Update insertion point */
            while ((NULL != next_node) && (col > next_node->index)) {
                is_first_in_line = FALSE;
                prev_node = next_node;
                next_node = next_node->next;
            }

            /* 2.2. Find insertion place */
            if (NULL == next_node) {
                /* DEBUG_PRINT("row %d: Inserting last node at col %d value %f", row_index, col, values[col]); */
                /* 2.2.1. Next node is NULL - we simply add it */
                result = LIST_NODE_append(&prev_node, values[col], col);
                if (E__SUCCESS != result) {
                    goto l_cleanup;
                }

                /* 2.2.2. If prev node was NULL we are the row's beginning */
                if (is_first_in_line) {
                    /* DEBUG_PRINT("row %d was NULL", row_index); */
                    is_first_in_line = FALSE;
                    row->begin = prev_node;
                }
            /* 2.3. Next node can be either equal or larger*/
            } else {
                if (col == next_node->index) {
                    /* 2.3.1. Found existing node - increase its value */
                    /* DEBUG_PRINT("row %d: Adding %f to existing node at col %d with prev value %f", row_index, values[col], col, next_node->value); */
                    next_node->value += values[col];
                    is_first_in_line = FALSE;

                    /* 2.3.2. Increase next node */
                    prev_node = next_node;
                    next_node = next_node->next;
                } else if (col < next_node->index) {
                    /* DEBUG_PRINT("row %d: Inserting node in the middle node at col %d value %f", row_index, col, values[col]); */
                    /* 2.3.3. Add new node
                     *        Note: row->begin already exists */
                    result = LIST_NODE_append(&prev_node, values[col], col);
                    if (E__SUCCESS != result) {
                        goto l_cleanup;
                    }

                    /* 2.3.3. Set the new node's next as next_node */
                    if (is_first_in_line) {
                        is_first_in_line = FALSE;
                        row->begin = prev_node;
                    }
                    prev_node->next = next_node;
                }
            }

            row->sum += values[col];
            matrix_data->columns_sum[col] += fabs(values[col]);
        }
    }

    result = E__SUCCESS;
l_cleanup:

    return result;
} 

static
void
spmat_list_mult(const matrix_t *mat, const double *v, double *multiplication_result)
{
    result_t result = E__UNKNOWN;
    spmat_row_t *rows_array = NULL;
    int i = 0;

    if ((NULL == mat) || (NULL == v) || (NULL == multiplication_result)) {
        /* Remark: this function signature shouldn't be void */
        result = E__NULL_ARGUMENT;
        goto l_cleanup;
    }

    rows_array = GET_ROWS_ARRAY(mat);
    for (i = 0 ; i < mat->n ; ++i) {
        multiplication_result[i] = LIST_NODE_scalar_multiply(rows_array[i].begin, v);
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
spmat_list_reduce_row(const spmat_row_t *original_row,
                       const double * vector_s,
                       double relevant_vector_s_value,
                       const int * s_indexes,
                       spmat_row_t *row_out)
{
    result_t result = E__UNKNOWN;
    const node_t *scanner = NULL;
    node_t *begin = NULL;
    node_t *row_end = NULL;
    double scanned_s_value = 0.0; /* 1 or -1 */
    int scanned_index = 0.0; /* 0 ... n */
    double sum = 0.0;

    /* 1. Check if original row is zeroes */
    if (NULL == original_row) {
        result = E__SUCCESS;
        goto l_cleanup;
    }

    /* 2. Go over nodes in the given row */
    for (scanner = original_row->begin ;
            NULL != scanner ;
            scanner = scanner->next) {
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
        scanned_s_value = vector_s[scanner->index];
        /* 2.1. Skip nodes with different s-value */
        if (scanned_s_value != relevant_vector_s_value) {
            continue;
        }

        /* 2.2. Append nodes with our s-value to the end of our new row */
        scanned_index = s_indexes[scanner->index];
        sum += scanner->value;
        result = LIST_NODE_append(&row_end, scanner->value, scanned_index);
        if (E__SUCCESS != result) {
            goto l_cleanup;
        }

        /* 2.3. First node only: set as the row's begin */
        if (NULL == begin) {
            begin = row_end;
        }
    }

    /* Success */
    row_out->begin = begin;
    row_out->sum = sum;

    result = E__SUCCESS;
l_cleanup:

    if (E__SUCCESS != result) {
        LIST_NODE_destroy(begin);
        begin = NULL;
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
    spmat_row_t *relevant_row_pointer = NULL;

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
        result = spmat_list_reduce_row(&GET_ROW(matrix, i),
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
        MATRIX_FREE_SAFE(matrix1);
        MATRIX_FREE_SAFE(matrix2);
    }

    return result;
}

void
SPMAT_LIST_print(const char *matrix_name, matrix_t *mat_in)
{
	spmat_row_t *relevant_row_pointer = NULL;
    const node_t *scanner = NULL;
	int row = 0;
	int col = 0;
	int last_scanner_index = 0;

    if (NULL == mat_in) {
        DEBUG_PRINT("matrix is null");
        return;
    }

    if (MATRIX_TYPE_SPMAT_LIST != mat_in->type) {
        DEBUG_PRINT("invalid matrix type (got %d)", mat_in->type);
        return;
    }

    if (NULL != mat_in) {
        printf("%s:\n----------------------\n", matrix_name);
        for (row = 0; row < mat_in->n; row++){
            col = 0;
            last_scanner_index = 0;
            relevant_row_pointer = &GET_ROW(mat_in, row);
            printf("(");
            for (scanner = relevant_row_pointer->begin ;
                    NULL != scanner ;
                    scanner = scanner->next) {

                for ( ; scanner->index > col ; ++col){
                    printf("%5.2f ", 0.0);
                }

                printf("%5.2f ", scanner->value);
                col++;

                if (NULL != scanner) {
                    last_scanner_index = col;
                }

            }

            for ( ; last_scanner_index < mat_in->n ; ++last_scanner_index) {
                printf("%5.2f ", 0.0);
            }
            printf("\n");
        }
    }

}

result_t
SPMAT_LIST_get_1norm(const matrix_t *matrix, double *norm_out)
{
    result_t result = E__SUCCESS;
    double norm = 0.0;
    const spmat_data_t *data = NULL;
    const double *columns_sum = NULL;
    int i = 0;

    /* 0. Input validation */
    if ((NULL == matrix) || (NULL == norm_out)) {
        result = E__NULL_ARGUMENT;
        goto l_cleanup;
    }

    /* Go over all the columns' norms */
    data = GET_SPMAT_DATA(matrix);
    columns_sum = data->columns_sum;
    for (i = 0 ; i < matrix->n ; ++i) {
        norm = fabs(MAX(norm, columns_sum[i]));
    }

    /* Success */
    *norm_out = norm;

    result = E__SUCCESS;
l_cleanup:

    return result;
}


result_t
SPMAT_LIST_decrease_rows_sums_from_diag(const matrix_t *matrix)
{
    result_t result = E__UNKNOWN;
    double *line_vector_tmp = NULL;
    int i = 0;
    spmat_row_t *row = NULL;
    double row_sum = 0.0;

    /* 0. Input validation */
    if (NULL == matrix) {
        result = E__NULL_ARGUMENT;
        goto l_cleanup;
    }

    /* 1. Allocate temp buffer */
    line_vector_tmp = (double *)malloc(sizeof(*line_vector_tmp) * matrix->n));
    if (NULL == line_vector_tmp) {
        result = E__MALLOC_ERROR;
        goto l_cleanup;
    }
    (void)memset(line_vector_tmp, 0, sizeof(*line_vector_tmp) * matrix->n);

    for (i = 0 ; i < matrix->n ; ++i) {
        line_vector_tmp[i] = matrix;
        row_sum = GET_ROW(matrix, i).sum;
        result = MATRIX_ADD_ROW(matrix, line_vector_tmp, i);
        if (E__SUCCESS != result) {
            goto l_cleanup;
        }
        line_vector_tmp[i] = 0.0;
    }

    /* Success */
    result = E__SUCCESS;
l_cleanup:

    FREE_SAFE(line_vector_tmp);

    return result;
}
