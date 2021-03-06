/**
 * @file matrix.h
 * @purpose Abstract matrix class: fields and virtual table
 */
#ifndef __MATRIX_H__
#define __MATRIX_H__

/* Includes ******************************************************************/
#include <stddef.h>

#include "results.h"
#include "common.h"


/* Macros ********************************************************************/
/* Number of elements on the matrix */
#define MATRIX_COUNT(m) ((size_t)((m)->n * (m)->n))

/* Get the element at a given cow and column */
#define MATRIX_IS_VALID_ROW_INDEX(A, i) (((A)->n > (i)) && (0 <= (i)))

/* Vtable macros */
#define MATRIX_VTABLE(m) ((m)->vtable)
#define MATRIX_FREE(m) (MATRIX_VTABLE(m)->free((m)))

#define MATRIX_ADD_ROW(m, row, i) \
    MATRIX_VTABLE((m))->add_row((m), (row), (i))

#define MATRIX_FREE_SAFE(m) do {                            \
    if (NULL != (m)) {                                      \
        MATRIX_FREE(m);                                     \
        (m) = NULL;                                         \
    }                                                       \
} while (0)


/* Enums *********************************************************************/
typedef enum matrix_type_e {
    MATRIX_TYPE_SPMAT_LIST,
    MATRIX_TYPE_MAX
} matrix_type_t;


/* Typedefs ******************************************************************/
typedef struct matrix_s matrix_t;

/*
 * Increase the values of the row in the matrix with a given row
 **/
typedef result_t (*matrix_add_row_f)(matrix_t *matrix,
                                     const double *row,
                                     int i);

/* Frees all resources used by A */
typedef void (*matrix_free_f)(matrix_t *matrix);

/* Multiplies matrix A by vector v, into result (result is pre-allocated) */
typedef void (*matrix_mult_f)(const matrix_t *matrix,
                              const double *vector,
                              double *result);

/* Given a vector v and matrix A,
 * calculate the v-transposed multiply M multiply v */
typedef double (*matrix_mult_vmv_f)(const matrix_t *matrix,
                                          const double *vector);

/**
 * Calcualte the 1-norm of the matrix
 *
 * @param matrix The matrix to calculate whose norm. Must be valid!
 *
 * @return The 1-norm of the matrix
 *
 * @remark If matrix is invalid the function's behavior will be unexpected
 */
typedef double (*matrix_get_1norm_f)(const matrix_t *matrix);


/* Structs *******************************************************************/
/**
 * Virtual table. The functions that every matrix module must implement
 **/
typedef struct matrix_vtable_s {
    matrix_add_row_f add_row;
    matrix_free_f free;
    matrix_mult_f mult; /* Calculate M*v */
    matrix_mult_vmv_f mult_vmv; /* Calculate v^T*M*v */
} matrix_vtable_t;

/* A square matrix implementation */
struct matrix_s {
    int n; /* Matrix size (n*n) */
    matrix_type_t type;
    const matrix_vtable_t *vtable;
    void *private;
};


/* Functions Declarations ****************************************************/
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

result_t
MATRIX_add_diag(matrix_t *matrix, double onenorm);


#endif /* __MATRIX_H__ */
