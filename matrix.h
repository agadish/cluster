/**
 * @file matrix.h
 * @purpose Reprent a matrix data structure, implement mathematical operations and de/serialization
 */
#ifndef __MATRIX_H__
#define __MATRIX_H__

/* Includes **************************************************************************************/
#include "results.h"
#include "common.h"


/* Macros ****************************************************************************************/
/* Number of elements on the matrix */
#define MATRIX_COUNT(m) ((size_t)((m)->header.rows * (m)->header.columns))

/* Get the element at a given cow and column */
#define MATRIX_AT(m, row, col) ((m)->data[((row) * (m)->header.columns) + (col)])


/* Structs ***************************************************************************************/
/* The matrix dimensions, from its binary format */
typedef struct matrix_header_s {
    int columns;
    int rows;
} matrix_header_t;

/* The matrix dimensions and data */
typedef struct matrix_s {
    matrix_header_t header;
    double * data;
} matrix_t;


/* Functions Declarations *************************************************************/
/*
 * @purpose Create a new matrix with the given dimensions
 *
 * @param rows The number of rows in the matrix
 * @param columns The number of columns in the matrix
 * @param matrix_out An external pointer which will point the new matrix
 *
 * @return One of result_t values
 *
 * @remark matrix must be freed using MATRIX_free
 */
result_t
MATRIX_allocate(int rows, int columns, matrix_t ** matrix_out);

/*
 * @purpose Create a new matrix by a given encoded matrix file path
 *
 * @param matrix_path The path to the matrix file. must be valid matrix!
 * @param matrix_out An external pointer which will point the new matrix
 *
 * @return One of result_t values
 *
 * @remark matrix must be freed using MATRIX_free
 */
result_t
MATRIX_open(const char * path, matrix_t **matrix_out);

/**
 * @purpose Free a matrix which was previously created by MATRIX_open or MATRIX_allocate
 *
 * @param matrix The matrix to free
 *
 * @remark Safe to call with NULL
 */
void
MATRIX_free(matrix_t * matrix);

/**
 * @purpose Write a matrix to a given file
 * @param matrix The matrix to write
 * @param path The output path
 *
 * @return One of result_t values
 */
result_t
MATRIX_write_to_file(const matrix_t * matrix, const char * path);

/**
 * @remark vector_out must be freed later by MATRIX_free
 */
result_t
MATRIX_random_vector(const int length, matrix_t ** vector_out);

/**
 * @purpose divides vector by scalar umber
 * @param vector The vector to divide
 * @param value The value to divide with
 *
 * @return One of result_t values
 */
result_t
MATRIX_div(matrix_t * vector, double value);


/**
 * @purpose Normalize each line of the matrix using the regular scalar multiplication
 * @param matrix Matrix to normalize
 *
 * @return One of result_t values
 */
result_t
MATRIX_normalize(matrix_t *matrix);

/**
 * @remark matrix must be valid
 */
double
MATRIX_calculate_magnitude(matrix_t *matrix);

/**
 * @purpose Checks if all the differences between the matrix' values are smaller than epsilon
 * @param matrix_a The first matrix
 * @param matrix_b The second matrix
 * @param epsilon The largest difference allowed
 *
 * @return TRUE if close enough, otherwise FALSE
 *
 * @remark Matrixes must have the same size
 */
bool_t
MATRIX_is_close(matrix_t * matrix_a, matrix_t * matrix_b, double epsilon);

/**
 * @remark vector_out must be freed later by MATRIX_free
 */
result_t
MATRIX_random_vector(const int length, matrix_t ** vector_out);

/**
 * Open parse adjacency matrix
 */
result_t
MATRIX_open_adjacency(const char * path, matrix_t **matrix_out);


#endif /* __MATRIX_H__ */
