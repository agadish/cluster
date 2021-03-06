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
#include "submatrix.h"


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
} spmat_data_t;


/* Macros ********************************************************************/
#define GET_SPMAT_DATA(matrix) ((spmat_data_t *)((matrix)->private))

#define GET_ROWS_ARRAY(matrix) (GET_SPMAT_DATA(matrix)->rows)

#define GET_ROW(matrix, row_index) (GET_ROWS_ARRAY(matrix)[row_index])

#define SPMAT_GET_EXPECTED_VALUE(smat, i, j) (                              \
    (smat->adj->neighbors[(i)] * smat->adj->neighbors_div_M[(j)])           \
)


/* Functions Declarations ****************************************************/
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

/**
 * @see matrix_split_f on matrix.h
 */

static
void
spmat_list_get_rows_sums(const submatrix_t *smat,
                         double *vector);

static
double
submat_spmat_list_mult_row_with_s(const submatrix_t *submatrix,
                                  int row_g,
                                  const double *s_vector);

/**
 * Calculate the multiplication result of the row_g'th row of B matrix,
 * a B matrix (withno hat!) with a given s_vector.
 * Algorithm is improved: is calculates sequences of zeroes effeciently
 *
 * @param smat The submatrix
 * @param row_g The row index to multiply
 * @param s_vector s_vector The s vector to multiply with
 *
 * @return The multiplication 
 */
static
double
submat_spmat_list_mult_row_with_s_no_hat_improved(const submatrix_t *smat,
                                                  int row_g,
                                                  const double *s_vector);


/* Virtual Table *************************************************************/
const matrix_vtable_t SPMAT_LIST_VTABLE = {
    .add_row = spmat_list_add_row,
    .free = spmat_list_free,
    .mult = spmat_list_mult,
    .mult_vmv = NULL,
};


/* Functions *****************************************************************/
result_t
SPMAT_LIST_allocate(int n, matrix_t **mat_out)
{
    result_t result = E__UNKNOWN;
    matrix_t *mat = NULL;
    spmat_data_t *spmat_data = NULL;
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

    /* 5. Validate matrix */
    mat->private = (void *)spmat_data;
    
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
                result = LIST_insert(row->list, next_node, values[col], col);
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
                    result = LIST_insert(row->list, next_node, values[col], col);
                    if (E__SUCCESS != result) {
                        goto l_cleanup;
                    }

                    /* 2.3.3. Set the new node's next as next_node */
                    prev_node->next = next_node;
                }
            }

            row->sum += values[col];
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
    if (NULL == original_row->list) {
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
        result = LIST_insert(reduced_list,
                             NULL,
                             scanner->value,
                             scanned_index);
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

/* TODO: Remove temp_s_indexes */
result_t
SUBMAT_SPMAT_LIST_split(submatrix_t *smat,
        const double * vector_s,
        int *temp_s_indexes,
        submatrix_t **matrix1_out,
        submatrix_t **matrix2_out)
{
    result_t result = E__UNKNOWN;
    matrix_t *orig = NULL;
    matrix_t *matrix1 = NULL;
    matrix_t *matrix2 = NULL;
    submatrix_t *smat1 = NULL;
    submatrix_t *smat2 = NULL;
    int i = 0;
    int matrix1_n = 0;
    double scanned_s_value = 0.0;
    spmat_row_t *relevant_row_pointer = NULL;

    /* 0. Input validation */
    /* Null arguments */
    if ((NULL == smat) ||
            (NULL == vector_s) ||
            (NULL == matrix1_out) ||
            (NULL == matrix2_out)) {
        result = E__NULL_ARGUMENT;
        goto l_cleanup;
    }

    /* 1. Create s-indexes vector, get matrix1's length */
    matrix1_n = VECTOR_create_s_indexes(vector_s,
                                        smat->g_length,
                                        temp_s_indexes);

    result = SPMAT_LIST_allocate(matrix1_n, &matrix1);
    if (E__SUCCESS != result) {
        goto l_cleanup;
    }

    /* 2. Create matrixes as spmat lists */
    /* 2.1. smat 1 */
    result = SUBMATRIX_create(smat->adj, matrix1, &smat1);
    if (E__SUCCESS != result) {
        goto l_cleanup;
    }
    matrix1 = NULL;

    result = SPMAT_LIST_allocate(smat->g_length - matrix1_n, &matrix2);
    if (E__SUCCESS != result) {
        goto l_cleanup;
    }

    /* 2.2. smat 2 */
    result = SUBMATRIX_create(smat->adj, matrix2, &smat2);
    if (E__SUCCESS != result) {
        goto l_cleanup;
    }
    matrix2 = NULL;

    orig = smat->orig;
    /* 3. Go over each row of the original smat */
    for (i = 0 ; i < orig->n ; ++i) {
        /* 3.1. Get the relevant s-value (1=matrix1, -1=matrix2) */
        scanned_s_value = vector_s[i]; /* 1 or -1 */

        /* 3.3. Split g vector */
        if (1.0 == scanned_s_value) {
            relevant_row_pointer = &GET_ROW(smat1->orig, smat1->g_length);
            smat1->g[smat1->g_length] = smat->g[i];
            ++smat1->g_length;
        } else if (-1.0 == scanned_s_value) {
            relevant_row_pointer = &GET_ROW(smat2->orig, smat2->g_length);
            smat2->g[smat2->g_length] = smat->g[i];
            ++smat2->g_length;
        } else {
            result = E__INVALID_S_VECTOR;
            goto l_cleanup;
        }

        /* 3.4. Add the filtered values in the row */
        result = spmat_list_reduce_row(&GET_ROW(orig, i),
                vector_s,
                scanned_s_value,
                temp_s_indexes,
                relevant_row_pointer);
        if (E__SUCCESS != result) {
            goto l_cleanup;
        }
    }

    /* Success */
    *matrix1_out = smat1;
    *matrix2_out = smat2;

    result = E__SUCCESS;
l_cleanup:

    if (E__SUCCESS != result) {
        MATRIX_FREE_SAFE(matrix1);
        MATRIX_FREE_SAFE(matrix2);
        SUBMATRIX_FREE_SAFE(smat1);
        SUBMATRIX_FREE_SAFE(smat2);
    }

    return result;
}

void
SPMAT_LIST_print(const char *matrix_name, matrix_t *mat_in)
{
    spmat_row_t *relevant_row_pointer = NULL;
    const node_t *scanner = NULL;
    const list_t *list = NULL;
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
            list = relevant_row_pointer->list;
            if (NULL != list) {
                for (scanner = list->first ;
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
            }

            for ( ; last_scanner_index < mat_in->n ; ++last_scanner_index) {
                printf("%5.2f ", 0.0);
            }
            printf("\n");
        }
    }

}

static
void
spmat_list_get_rows_sums(const submatrix_t *smat,
                         double *vector)
{
    int row_g = 0;
    int row_i = 0;
    int col_g = 0;
    int col_i = 0;
    double current_sum = 0.0;
    double current_cell = 0.0;
    list_t *current_list = NULL;
    node_t *s = NULL;

    for (row_g = 0 ; row_g < smat->g_length ; ++row_g) {
        current_sum = 0.0;
        row_i = smat->g[row_g];
        current_list = GET_ROW(smat->orig, row_g).list;
        if (NULL != current_list) {
            s = current_list->first;
        }
        for (col_g = 0 ; col_g < smat->g_length ; ++col_g) {
            col_i = smat->g[col_g];
            /* Skip irrelevant values from original matrix */
            while ((NULL != s) && (s->index < col_g)) {
                s = s->next;
            }
            if ((NULL != s) && (s->index == col_g)) {
                current_cell = s->value;
            } else {
                current_cell = 0.0;
            }

            current_cell -= SPMAT_GET_EXPECTED_VALUE(smat, row_i, col_i);
            /* DEBUG_PRINT("mod[%d,%d]=%f", row_g, col_g, current_cell); */
            current_sum += current_cell;
        }
        vector[row_g] = current_sum;
    }
}

double
SUBMAT_SPMAT_LIST_get_1norm(const submatrix_t *smat,
                            double *tmp_row_sums)
{
    double norm = 0.0;
    double current_row_norm = 0.0;
    int trow_g = 0;
    int trow_i = 0;
    int tcol_g = 0;
    int tcol_i = 0;
    list_t *l = NULL;
    node_t *s = NULL;
    double cell_value = 0.0;
    double a = 0.0;
    double diag_add = 0.0;
    double expected_value;


    /* TODO: Handle zero sized */
    if (smat->g_length == 0){ 
        printf("WARNING WARNING WARNING DAMN SUBMAT_SPMAT_LIST_get_1norm got zero sized mat\n");
    }

    /* Note: The adjacency matrix is symmetric, therefore 1-norm can be done on
     *       either max row sum or max column sum */

    spmat_list_get_rows_sums(smat, tmp_row_sums);

    for (trow_g = 0 ; trow_g < smat->g_length ; ++trow_g) {
        /* Go over the sub rows */
        current_row_norm = 0.0;
        trow_i = smat->g[trow_g];
        l = GET_ROW(smat->orig, trow_g).list;
        if (NULL == l) {
            /* Zeroes row cannot increase norm */
            continue;
        }

        s = l->first;
        for (tcol_g = 0 ; tcol_g < smat->g_length ; ++tcol_g) {
            /* Go over the sub columns */
            tcol_i = smat->g[tcol_g];

            /* 1. Set s as the next member with greater/equal index */
            for (; (NULL != s) && (s->index) < tcol_g ; s = s->next);
                /* DEBUG_PRINT("increasing s_index from %d (tcol_g is %d)", s->index, tcol_g); */

            /* 2. Check if row is over */
            if ((NULL == s) || (s->index > tcol_g)) {
                /* if (NULL != s) { */
                /*     DEBUG_PRINT("s=%p, tcol_g=%d s->index=%d", (void *)s, tcol_g, s->index); */
                /* } else { */
                /*     DEBUG_PRINT("s=%p, tcol_g=%d", (void *)s, tcol_g); */
                /* } */
                a = 0.0;
            } else {
                a = s->value;
            }

            /* 3. Check if the found index is what we're searching for */
            expected_value = SPMAT_GET_EXPECTED_VALUE(smat, trow_i, tcol_i);
            if (tcol_g == trow_g) {
                /* DEBUG_PRINT("lol found diag"); */
                diag_add = smat->add_to_diag - tmp_row_sums[tcol_g];
            } else {
                diag_add = 0.0;
            }
            cell_value = a - expected_value + diag_add;
            current_row_norm += fabs(cell_value);
            /* DEBUG_PRINT("1norm[%d,%d]:a=%f,expected=%f,add=%f", tcol_g, trow_g, a, expected_value, diag_add); */
            /* DEBUG_PRINT("cell[%d,%d]=%f", tcol_g, trow_g, cell_value); */
        }

        norm = MAX(norm, current_row_norm);
    }

    /* Success */
    return norm;
}

static
double
submat_spmat_list_mult_row_with_s_no_hat_improved(const submatrix_t *smat,
                                                  int row_g,
                                                  const double *s_vector)
{
    double result = 0.0;
    list_t *l = NULL;
    node_t *s = NULL;
    int col_i = 0;
    int row_i = 0;
    double expected_value = 0.0;
    int prev_g = 0;
    int current_g = 0;
    double kj_zeroes_sums = 0.0;
    double values_sum = 0.0;

    row_i = smat->g[row_g];

    l = GET_ROW(smat->orig, row_g).list;
    if (NULL != l) {
        /* If line is NULL, continue with calculation */
        for (s = l->first ; NULL != s ; s = s->next) {
            current_g = s->index;

            /* Add previous 0.0 cells */
            kj_zeroes_sums += VECTOR_scalar_multiply_with_s(
                &smat->adj->neighbors_div_M[prev_g],
                &s_vector[prev_g],
                current_g - prev_g
            );
            /* If current_g is last, the loop will break and the fix will be
             * applied */
            prev_g = current_g + 1;

            col_i = smat->g[current_g];
            expected_value = SPMAT_GET_EXPECTED_VALUE(smat, row_i, col_i);
            values_sum += (s->value - expected_value) * (s_vector[current_g] > 0 ? 1 : -1);
        }
    }
    current_g = smat->g_length;

    kj_zeroes_sums += VECTOR_scalar_multiply_with_s(
        &smat->adj->neighbors_div_M[prev_g],
        &s_vector[prev_g],
        MAX(current_g - prev_g, 0)
    );

    kj_zeroes_sums *= (double)smat->adj->neighbors[row_i];

    result = values_sum - kj_zeroes_sums + (smat->add_to_diag * s_vector[row_g]);
    result += (smat->add_to_diag) * s_vector[row_g];

    return result;
}

static
double
submat_spmat_list_mult_row_with_s(const submatrix_t *smat,
                                  int row_g,
                                  const double *s_vector)
{
    double result = 0.0;
    list_t *l = NULL;
    node_t *s = NULL;
    int col_g = 0;
    int col_i = 0;
    int row_i = 0;
    double a = 0.0;
    double expected_value = 0.0;
    bool_t is_zero = TRUE;
    double row_sum = 0.0;

    row_i = smat->g[row_g];

    l = GET_ROW(smat->orig, row_g).list;
    if (NULL != l) {
        /* If line is NULL, continue with calculation */
        s = l->first;
    }

    /* Go over the columns */
    for (col_g = 0 ;  col_g < smat->g_length ; ++col_g) {
        col_i = smat->g[col_g];

        expected_value = SPMAT_GET_EXPECTED_VALUE(smat, row_i, col_i);

        /* Skip irrelevant values from original matrix */
        while ((NULL != s) && (s->index < col_g)) {
            s = s->next;
        }

        /* If column is non-zero take s-value */
        is_zero = (NULL == s) || (s->index > col_g);
        if (is_zero) {
            a = 0.0;
        } else {
            a = s->value;
            s = s->next;
        }
        row_sum += (a - expected_value);

        result += (a - expected_value) * s_vector[col_g];
    }

    result += (smat->add_to_diag - row_sum) * s_vector[row_g];

    return result;
}


double
SUBMAT_SPMAT_LIST_calculate_q(const submatrix_t *submatrix,
                              const double *s_vector)
{
    double current_row_mul = 0.0;
    double mult_vmv = 0.0;
    int row_g = 0;

    /* Multiply each row with s-vector */
    for (row_g = 0 ; row_g < submatrix->g_length ; ++row_g) {
        current_row_mul = 0.0;

        /* Add to result v[row] times M[row, :]*v */
        current_row_mul = submat_spmat_list_mult_row_with_s(submatrix,
                                                            row_g,
                                                            s_vector);
        mult_vmv += (s_vector[row_g] * current_row_mul);
    }

    return mult_vmv;
}

#if 0
result_t
SUBMAT_SPMAT_LIST_split(submatrix_t *smat,
                        const double *s_vector,
                        submatrix_t **split1_out,
                        submatrix_t **split2_out)
{
    result_t result = E__UNKNOWN;
    int i = 0;
    int split1_length = 0;
    int split2_length = 0;
    submatrix_t *split1 = NULL;
    submatrix_t *split2 = NULL;

    /* 0. Input validation */
    /* 0.1. Null arguments */
    if ((NULL == smat) || (NULL == s_vector) ||
            (NULL == split1_out) || (NULL == split2_out)) {
        result = E__NULL_ARGUMENT;
        goto l_cleanup;
    }

    /* 0.2. Input matrix must be spmat list */
    if (MATRIX_TYPE_SPMAT_LIST != smat->orig->type) {
        result = E__INVALID_MATRIX_TYPE;
        goto l_cleanup;
    }

    result = SUBMATRIX_create(smat->adj, &split1);
    if (E__SUCCESS != result) {
        goto l_cleanup;
    }

    result = SUBMATRIX_create(smat->adj, &split2);
    if (E__SUCCESS != result) {
        goto l_cleanup;
    }

    /* XXX: The g-vector is allocated for size n, the g-vector is
     *      non-ascending therefore we will never overflow this vector */
    /* 2. Update g vector */
    for (i = 0 ; i < smat->g_length ; ++i) {
        if (1.0 == s_vector[i]) {
            split1->g[split1_length] = smat->g[i];
            ++split1_length;
        } else { /* Equals -1.0 */
            split2->g[split2_length] = smat->g[i];
            ++split2_length;
        }
    }
    /* 3. Update g lengths */
    split1->g_length = split1_length;
    split2->g_length = split2_length;

    /* Success */
    *split1_out = split1;
    *split2_out = split2;
    result = E__SUCCESS;
l_cleanup:
    if (E__SUCCESS != result) {
        SUBMATRIX_FREE_SAFE(split1);
        SUBMATRIX_FREE_SAFE(split2);
    }

    return result;
}
#endif

void
SUBMAT_SPMAT_LIST_mult(const submatrix_t *submatrix,
                       const double *vector,
                       double *result)
{
    int row_g = 0;
    double current_row_mul = 0.0;

    for (row_g = 0 ; row_g < submatrix->g_length ; ++row_g) {
        current_row_mul = submat_spmat_list_mult_row_with_s(submatrix,
                                                            row_g,
                                                            vector);
        result[row_g] = current_row_mul;
    }
}

double
SUBMAT_SPMAT_LIST_calc_q_score(const submatrix_t *smat,
                               const double *vector,
                               int row_g)
{
	double q_part1 = 0.0;
	double expected_value = 0.0;
	double q_score = 0.0;
    int row_i = 0;

    q_part1 = submat_spmat_list_mult_row_with_s_no_hat_improved(smat,
                                                                row_g,
                                                                vector);
    row_i = smat->g[row_g];
    expected_value = SPMAT_GET_EXPECTED_VALUE(smat, row_i, row_i);
    q_score = 4 * (vector[row_g] * q_part1 + expected_value);

    return q_score;
}
