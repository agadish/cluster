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
/* spmat_row_t * matrix->private: n-sized array of spmat_row_t */
typedef struct spmat_row_s {
    /* Begin of row's linked list */
    node_t *begin;
    double sum;
} spmat_row_t;


/* Macros ****************************************************************************************/
#define GET_ROWS_ARRAY(list) ((spmat_row_t *)((list)->private))

#define GET_ROW(list, row_index) (GET_ROWS_ARRAY(list)[row_index])


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

static
result_t
spmat_list_init_row(spmat_row_t *row, const double *data, int data_length);

/* Functions *************************************************************************************/
result_t
SPMAT_LIST_allocate(int n, matrix_t **mat_out)
{
    result_t result = E__UNKNOWN;
    matrix_t *mat = NULL;
    spmat_row_t *rows_array = NULL;
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
    mat->n = n;
    mat->add_row = spmat_list_add_row;
    mat->free = spmat_list_free;
    mat->mult = spmat_list_mult;
    /* TODO: Implement */
    mat->rmult = NULL;
    mat->private = NULL;

    /* 3. Rows array */
    rows_array_size = n * sizeof(*rows_array);
    rows_array = (spmat_row_t *)malloc(rows_array_size);
    if (NULL == rows_array) {
        result = E__MALLOC_ERROR;
        goto l_cleanup;
    }
    (void)memset(rows_array, 0, rows_array_size);
    mat->private = (void *)rows_array;

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
    spmat_row_t *rows_array = NULL;
    int i = 0;

    if (NULL != mat) {
        rows_array = GET_ROWS_ARRAY(mat);
        if (NULL != rows_array) {
            for (i = 0 ; i < mat->n ; ++i) {
                LIST_NODE_destroy(rows_array[i].begin);
                rows_array[i].begin = NULL;
            }

            FREE_SAFE(rows_array);
            mat->private = NULL;
        }
        FREE_SAFE(mat);
    }
}

result_t
spmat_list_add_row(matrix_t *mat, const double *data, int row_index)
{
    result_t result = E__UNKNOWN;
    spmat_row_t *rows_array = NULL;

    /* 1. Input validation */
    if ((NULL == mat) || (NULL == mat->private) || (NULL == data)) {
        result = E__NULL_ARGUMENT;
        goto l_cleanup;
    }

    if (!MATRIX_IS_VALID_ROW_INDEX(mat, row_index)) {
        result = E__INVALID_ROW_INDEX;
        goto l_cleanup;
    }

    rows_array = GET_ROWS_ARRAY(mat);

    /* 2. Destory current row data */
    if (NULL != rows_array[row_index].begin) {
        LIST_NODE_destroy(rows_array[row_index].begin);
        rows_array[row_index].begin = NULL;
        rows_array[row_index].sum = 0;
    }

    /* 3. Initialise row with new data */
    result = spmat_list_init_row(&rows_array[row_index], data, mat->n);
    if (E__SUCCESS != result) {
        goto l_cleanup;
    }

    result = E__SUCCESS;
l_cleanup:

    return result;
} 

static
result_t
spmat_list_init_row(spmat_row_t *row, const double *data, int data_length)
{
    result_t result = E__UNKNOWN;
    int i = 0;
    node_t *last_node = NULL;
    node_t *current_node = NULL;

    /* 3. Create a new row */
    for (i = 0 ; data_length > i ; ++i) {
        if (0 != data[i]) {
            /* 3.1. Found a non-zero value, create a node for it */
            result = LIST_NODE_create(data[i], i, &current_node);
            if (E__SUCCESS != result) {
                goto l_cleanup;
            }
            
            /* 3.2. Link with last node */
            if (NULL == last_node) {
                /* 3.2.1. This node is the beginning of the row */
                row->begin = current_node;
            } else {
                /* 3.2.2. Link to last node */
                last_node->next = current_node;
            }
            /* 3.2.1. Is this the first node in the row? */
            last_node = current_node;
            row->sum += data[i];
        }
    }

    result = E__SUCCESS;
l_cleanup:

    if (E__SUCCESS != result) {
        LIST_NODE_destroy(row->begin);
        row->begin = NULL;
    }
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
SPMAT_LIST_print(matrix_t **mat_in){
	spmat_row_t *relevant_row_pointer = NULL;
    const node_t *scanner = NULL;
	int length = 0;
	int row = 0;
	int col = 0;
	int last_scanner_index = 0;

    if (NULL != mat) {
        length = mat_in->n;
        if (NULL != rows_array) {
        	printf("MATRIX");
        	for (row = 0; row < n; row++){
        		col = 0;
        		last_scanner_index = 0;
        		relevant_row_pointer = &GET_ROW((*mat_in), row);
        		printf("||");
        	    for (scanner = relevant_row_pointer->begin ;
        	            NULL != scanner ;
        	            scanner = scanner->next) {

        	    	while (scanner->index > col){
        	    		printf(" 0 ,");
        	    	}

        	    	printf(" %lf ,", scanner-> value );
        	    	col++;

        	    	if (NULL != scanner){ last_scanner_index = scanner->index; }

        	    }

        	    while (last_scanner_index < length) { printf(" 0 ,"); }
        	    printf("|| \n");
        	}
        }
    }

}

