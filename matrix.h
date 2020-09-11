/**
 * @file matrix.h
 * @purpose Reprent a matrix data structure, implement mathematical operations and de/serialization
 */
#ifndef __MATRIX_H__
#define __MATRIX_H__

/* Includes **************************************************************************************/
#include <stddef.h>

#include "results.h"
#include "common.h"


/* Macros ****************************************************************************************/
/* Number of elements on the matrix */
#define MATRIX_COUNT(m) ((size_t)((m)->n * (m)->n))

/* Get the element at a given cow and column */
#define MATRIX_IS_VALID_ROW_INDEX(A, i) (((A)->n > (i)) && (0 <= (i)))

#define MATRIX_FREE(m) ((m)->free((m)))
#define MATRIX_ADD_ROW(m, row, i) ((m)->add_row((m), (row), (i)))


/* Enums *****************************************************************************************/
typedef enum matrix_type_e {
    MATRIX_TYPE_RAW = 0,
    MATRIX_TYPE_SPMAT_LIST,
    MATRIX_TYPE_SPMAT_ARRAY,
    MATRIX_TYPE_MAX
} matrix_type_t;


/* Typedefs ***************************************************************************************/
typedef struct matrix_s matrix_t;

/* Adds row i the matrix. Called before any other call, exactly n times in order (i = 0 to n-1) */
typedef result_t (*matrix_add_row_f)(matrix_t *matrix, const double *row, int i);

/*
 * Gets the i-th row if the matrix.
 *
 * @param matrix The matrix to get whose row
 * @param row_index The index of the row
 * @param row A pre-allocated n-sizezd buffer which the row will be written to
 */
typedef void (*matrix_get_row_f)(matrix_t *matrix, int row_index, double *row);

/* Frees all resources used by A */
typedef void (*matrix_free_f)(matrix_t *matrix);

/* Multiplies matrix A by vector v, into result (result is pre-allocated) */
typedef void (*matrix_mult_f)(const matrix_t *matrix, const double *vector, double *result);


/* Structs ***************************************************************************************/
/* A square matrix implementation */
struct matrix_s {
	int	n; /* Matrix size (n*n) */
    matrix_type_t type;
    matrix_add_row_f add_row;
#if 0
    matrix_get_row_f get_row;
#endif
    matrix_free_f free;
    matrix_mult_f mult; /* Multiple the matrix with a column vector */
    matrix_mult_f rmult; /* Multiple a line vector with the matrix */
	void *private;
};


/* Functions Declarations *************************************************************/
/* 
 * @purpose Create a matrix accordingly to the given type
 * 
 * @param n The matrix's size nxn
 * @param type The matrix's implementation
 * @param matrix_out The matrix creadet
 *
 * @return One of result_t values
 */
result_t
MATRIX_create_matrix(int n, matrix_type_t type, matrix_t **matrix_out);

/**
 * @remark vector_out must be freed later by MATRIX_free
 */
result_t
MATRIX_random_vector(const int length, double ** vector_out);

#ifdef NEED_COL_VECTOR_TRANSPOSE
/*
 * @purpose Create a "transpose" vector to row vector
 *
 * @param vector_in column vector.
 * @param vector_out row vector corresponding to input vector
 *
 * @return One of result_t values
 */
result_t
MATRIX_col_vector_transpose(matrix_t *vector_in, matrix_t **vector_out);
#endif


#endif /* __MATRIX_H__ */
