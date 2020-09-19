/*
 * @file spmat_list.c
 * @purpose Sparse matrix implemented using linked lists
 */

/* Includes ******************************************************************/
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "list.h"
#include "results.h"
#include "matrix.h"
#include "spmat_list.h"
#include "common.h"
#include "debug.h"
#include "vector.h"


/* Structs *******************************************************************/
/* spmat_row_t * matrix->private: n-sized array of spmat_row_t */
typedef struct spmat_row_s {
    /* Begin of row's linked list */
    list_t *list;
    double sum;
    int index;
} spmat_row_t;

/* spmat_row_t * matrix->private: n-sized array of spmat_row_t */
typedef struct spmat_data_s {
    /* Begin of row's linked list */
    spmat_row_t *rows;
    double *columns_sum;
    double *temp_buffer; /* n-sized buffer required for various calculations */
} spmat_data_t;


/* Macros ********************************************************************/
#define GET_SPMAT_DATA(matrix) ((spmat_data_t *)((matrix)->private))

#define GET_ROWS_ARRAY(matrix) (GET_SPMAT_DATA(matrix)->rows)

#define GET_ROW(matrix, row_index) (GET_ROWS_ARRAY(matrix)[row_index])


/* Functions Declarations ************************************************************************/
/**
 * @purpose updating a row in a sparse matrix implemented by linked lists
 * @param A input Matrix
 * @param row input updated values for the row
 * @param i index of the row we want to update inside the matrix
 *
 * @return One of result_t values
 *
 */
static
result_t
spmat_list_add_row(matrix_t *A, const double *row, int i);

/**
 * @purpose free allocated memory for a matrix
 * @param A input Matrix
 *
 */
static
void
spmat_list_free(matrix_t *A);

/**
 * @purpose multiplying a matrix with vector
 * @param A input Matrix
 * @param v input vector
 * @param result output vector
 *
 */
static
void
spmat_list_mult(const matrix_t *A, const double *v, double *result);

/**
 * @purpose creating a list of incrementing indexes for each one of the groups
 *          according to values of s vector
 * @param vector_s input s vector
 * @param length length of s vector
 * @param s_indexes_out output result list
 * @param matrix1_n_out length of first group after division
 *
 * @return One of result_t values
 * @remark vector_s must be valid with given length, and values 1 or -1
 */
static
result_t
spmat_list_create_s_indexes(const double * vector_s,
                            int length,
                            int **s_indexes_out,
                            int *matrix1_n_out);
/**
 * @purpose reducing a row inside a sparse matrix to values only
 *          relevant to a specific group after division
 * @param original_row input row list from sparse matrix
 * @param vector_s input s vector describing division
 * @param relavant_vector_s_value inoput 1 or -1, describing which group we are building
 * @param s_indexes input list of incrementing indexes for each one of the groups
 * @param row_out output new row
 *
 * @return One of result_t values
 *
 */
static
result_t
spmat_list_reduce_row(const spmat_row_t *original_row,
                       const double * vector_s,
                       double relevant_vector_s_value,
                       const int * s_indexes,
                       spmat_row_t *row_out);

static
double
spmat_list_get_1norm(const matrix_t *matrix);

static 
result_t
spmat_list_decrease_rows_sums_from_diag(matrix_t *matrix);

static
void
spmat_list_initialise_rows_numbers(matrix_t *mat);

/**
 * @purpose multiplying row vector, matrix, and same vector as col vector
 * @param mat input matrix
 * @param v input vector

 * @return result of multiplication
 *
 */
static
double
spmat_list_matrix_vector_sandwich(const matrix_t *mat, const double *v);

/**
 * @see matrix_divide_f on matrix.h
 */
static
result_t
spmat_list_divide_matrix(matrix_t *matrix,
        const double * vector_s,
        int *temp_s_index,
        matrix_t **matrix1_out,
        matrix_t **matrix2_out);


/* Virtual Table *************************************************************/
const matrix_vtable_t SPMAT_LIST_VTABLE = {
    .add_row = spmat_list_add_row,
    .free = spmat_list_free,
    .mult = spmat_list_mult,
    .mult_vmv = spmat_array_matrix_vector_sandwich,
    .get_1norm = spmat_list_get_1norm,
    .decrease_rows_sums_from_diag = spmat_list_decrease_rows_sums_from_diag,
    .divide = spmat_list_divide_matrix
};


/* Functions *****************************************************************/
result_t
SPMAT_LIST_allocate(int n, bool_t should_initialise_row_numbers, matrix_t **mat_out)
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
    mat->private = NULL;
    mat->vtable = &SPMAT_LIST_VTABLE;
    mat->n = n;
    mat->type = MATRIX_TYPE_SPMAT_LIST;

    /* 3. spmat struct */
    spmat_data = (spmat_data_t *)malloc(sizeof(*spmat_data));
    if (NULL == spmat_data) {
        result = E__MALLOC_ERROR;
        goto l_cleanup;
    }
    (void)memset(spmat_data, 0, sizeof(*spmat_data));

    /* 4. Rows array */
    /* 4.1. Allocate */
    rows_array_size = n * sizeof(*rows_array);
    rows_array = (spmat_row_t *)malloc(rows_array_size);
    if (NULL == rows_array) {
        result = E__MALLOC_ERROR;
        goto l_cleanup;
    }
    /* 4.2. Initialise as NULL, sum 0.0, index 0 */
    (void)memset(rows_array, 0, rows_array_size);

    /* 4.3. Assign rows */
    spmat_data->rows = rows_array;

    /* 5. Columns sums */
    columns_sum = (double *)malloc(sizeof(*columns_sum) * n);
    if (NULL == columns_sum) {
        result = E__MALLOC_ERROR;
        goto l_cleanup;
    }
    (void)memset(columns_sum, 0, sizeof(*columns_sum) * n);
    spmat_data->columns_sum = columns_sum;

    /* Temp buffer */
    spmat_data->temp_buffer = (double *)malloc(sizeof(*columns_sum) * n);
    if (NULL == spmat_data->temp_buffer) {
        result = E__MALLOC_ERROR;
        goto l_cleanup;
    }

    /* 7. Validate matrix */
    mat->private = (void *)spmat_data;

    /* 8. Initialise row's indexes in increasing order */
    if (should_initialise_row_numbers) {
        spmat_list_initialise_rows_numbers(mat);
    }
    
    DEBUG_PRINT("%s: addr %p n=%d\n", __func__, (void *)mat, n);
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

    DEBUG_PRINT("%s: addr %p n=%d\n", __func__, (void *)mat, mat->n);
    if (NULL != mat) {
        spmat_data = GET_SPMAT_DATA(mat);
        if (NULL != spmat_data) {
            rows_array = spmat_data->rows;
            if (NULL != rows_array) {
                for (i = 0 ; i < mat->n ; ++i) {
                    LIST_destroy(rows_array[i].list);
                    rows_array[i].list = NULL;
                }

                FREE_SAFE(spmat_data->rows);
                rows_array = NULL;
            }
            
            FREE_SAFE(spmat_data->columns_sum);
            FREE_SAFE(spmat_data->temp_buffer);
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
    if (NULL != row->list) {
        next_node = row->list->first;
    }

    /* 2. Add row to matrix */
    for (col = 0 ; mat->n > col ; ++col) {
        if (0 != values[col]) {
            /* 2.1. Update insertion point */
            while ((NULL != next_node) && (col > next_node->index)) {
                prev_node = next_node;
                next_node = next_node->next;
            }

            /* 2.2. Find insertion place */
            if (NULL == next_node) {
                /* DEBUG_PRINT("row %d: Inserting last node at col %d value %f", row_index, col, values[col]); */
                /* 2.2.1. Next node is NULL - we simply add it */
                if (NULL == row->list) {
                    result = LIST_create(&row->list);
                    if (E__SUCCESS != result) {
                        goto l_cleanup;
                    }
                }
                result = LIST_insert(row->list, prev_node, values[col], col);
                if (E__SUCCESS != result) {
                    goto l_cleanup;
                }
            /* 2.3. Next node can be either equal or larger*/
            } else {
                if (col == next_node->index) {
                    /* 2.3.1. Found existing node - increase its value */
                    next_node->value += values[col];

                    /* 2.3.2. Increase next node */
                    prev_node = next_node;
                    next_node = next_node->next;
                } else if (col < next_node->index) {
                    /* DEBUG_PRINT("row %d: Inserting node in the middle node at col %d value %f", row_index, col, values[col]); */
                    /* 2.3.3. Add new node
                     *        Note: row->list already exists */
                    if (NULL == row->list) {
                        result = LIST_create(&row->list);
                        if (E__SUCCESS != result) {
                            goto l_cleanup;
                        }
                    }
                    result = LIST_insert(row->list, prev_node, values[col], col);
                    if (E__SUCCESS != result) {
                        goto l_cleanup;
                    }

                    /* 2.3.3. Set the new node's next as next_node */
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
        multiplication_result[i] = LIST_scalar_multiply(rows_array[i].list, v);
    }

    result = E__SUCCESS;
l_cleanup:

    UNUSED_ARG(result);
    return;
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
    list_t *reduced_list = NULL;
    double scanned_s_value = 0.0; /* 1 or -1 */
    int scanned_index = 0.0; /* 0 ... n */
    double sum = 0.0;

    /* 1. Check if original row is zeroes */
    if (NULL == original_row) {
        result = E__SUCCESS;
        goto l_cleanup;
    }

    /* 2. Create new list */
    result = LIST_create(&reduced_list);
    if (E__SUCCESS != result) {
        goto l_cleanup;
    }

    /* 3. Go over nodes in the given row */
    for (scanner = original_row->list->first ;
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
        result = LIST_insert(reduced_list, NULL, scanner->value, scanned_index);
        if (E__SUCCESS != result) {
            goto l_cleanup;
        }
    }

    /* Success */
    row_out->list = reduced_list;
    row_out->sum = sum;
    row_out->index = original_row->index;

    result = E__SUCCESS;
l_cleanup:
    if (E__SUCCESS != result) {
        LIST_destroy(reduced_list);
        reduced_list = NULL;
    }

    return result;
}

static
result_t
spmat_list_divide_matrix(matrix_t *matrix,
        const double * vector_s,
        int *temp_s_indexes,
        matrix_t **matrix1_out,
        matrix_t **matrix2_out)
{
    result_t result = E__UNKNOWN;
    matrix_t *matrix1 = NULL;
    matrix_t *matrix2 = NULL;
    int i = 0;
    int matrix1_n = 0;
    double scanned_s_value = 0.0;
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
    matrix1_n = VECTOR_create_s_indexes(vector_s,
                                     matrix->n,
                                     temp_s_indexes);

    /* 2. Create matrixes as spmat lists */
    /* 2.1. Matrix 1 */
    result = SPMAT_LIST_allocate(matrix1_n, FALSE, &matrix1);
    if (E__SUCCESS != result) {
        goto l_cleanup;
    }

    /* 2.2. Matrix 2 */
    result = SPMAT_LIST_allocate(matrix->n - matrix1_n, FALSE, &matrix2);
    if (E__SUCCESS != result) {
        goto l_cleanup;
    }

    /* 3. Go over each row of the original matrix */
    for (i = 0 ; i < matrix->n ; ++i) {
        /* 3.1. Get the relevant s-value (1=matrix1, -1=matrix2) */
        scanned_s_value = vector_s[i]; /* 1 or -1 */
        /* 3.2. Get the row index within the selected matrix */
        scanned_s_index = (int)temp_s_indexes[i];

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
                temp_s_indexes,
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
        printf("%s:\n----------------------\n",
                (NULL == matrix_name) ? "matrix" : matrix_name);
        for (row = 0; row < mat_in->n; row++){
            col = 0;
            last_scanner_index = 0;
            relevant_row_pointer = &GET_ROW(mat_in, row);
            printf("(");
            for (scanner = relevant_row_pointer->list->first ;
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

static
double
spmat_list_get_1norm(const matrix_t *matrix)
{
    double norm = 0.0;
    const spmat_data_t *data = NULL;
    const double *columns_sum = NULL;
    int i = 0;

    /* Go over all the columns' norms */
    data = GET_SPMAT_DATA(matrix);
    columns_sum = data->columns_sum;
    for (i = 0 ; i < matrix->n ; ++i) {
        norm = fabs(MAX(norm, columns_sum[i]));
    }

    return norm;
}


static
result_t
spmat_list_decrease_rows_sums_from_diag(matrix_t *matrix)
{
    result_t result = E__UNKNOWN;
    double *line_vector_tmp = NULL;
    int i = 0;
    double row_sum = 0.0;

    /* 0. Input validation */
    if (NULL == matrix) {
        result = E__NULL_ARGUMENT;
        goto l_cleanup;
    }

    line_vector_tmp = GET_SPMAT_DATA(matrix)->temp_buffer;
    (void)memset(line_vector_tmp, 0, sizeof(*line_vector_tmp) * matrix->n);
    for (i = 0 ; i < matrix->n ; ++i) {
        row_sum = GET_ROW(matrix, i).sum;
        line_vector_tmp[i] = row_sum;
        result = MATRIX_ADD_ROW(matrix, line_vector_tmp, i);
        if (E__SUCCESS != result) {
            goto l_cleanup;
        }
        line_vector_tmp[i] = 0.0;
    }

    /* Success */
    result = E__SUCCESS;
l_cleanup:

    return result;
}


result_t
SPMAT_LIST_write_neighbors(const matrix_t *matrix, FILE *file)
{
    result_t result = E__UNKNOWN;
    size_t result_write = 0;
    int i = 0;
    int index = -1;

    if ((NULL == matrix) || (NULL == file)) {
        result = E__NULL_ARGUMENT;
        goto l_cleanup;
    }

    /* Go over the matrix lines */
    for (i = 0 ; i < matrix->n ; ++i) {
        /* Get the line's index */
        index = GET_ROW(matrix, i).index;

        /* Write the line's index */
        result_write = fwrite(&index, sizeof(index), 1, file);
        if (1 != result_write) {
            result = E__FWRITE_ERROR;
            goto l_cleanup;
        }
    }

    result = E__SUCCESS;
l_cleanup:

    return result;
}

void
spmat_list_initialise_rows_numbers(matrix_t *mat)
{
    int i = 0;

    for (i = 0 ; i < mat->n ; ++i) {
        GET_ROW(mat, i).index = i;
    }
}

static
double
spmat_array_matrix_vector_sandwich(const matrix_t *mat, const double *v)
{
    spmat_row_t *rows_array = NULL;
    double row_sum = 0.0;
    double result = 0.0;
    int row = 0;

    rows_array = GET_ROWS_ARRAY(mat);
    for (row = 0 ; row < mat->n ; ++row) {
        if (rows_array[i].list == NULL){
        	continue;
        } else {
            for (scanner = rows_array[i].list->first ; NULL != scanner ; scanner = scanner->next) {
                row_sum += (scanner->value * v[scanner->index]);
            }
        }

        result += row_sum * v[row]
        row_sum = 0.0;
    }

    return result;
}
