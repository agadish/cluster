/**
 * @file lazy_matrix.h
 * @purpose A lazy matrix is a representation of a matrix file, on which
 *          only 1 line is held in the memory simultaneously.
 *          You can read lines using LAZY_MATRIX_read_next_line
 */
#ifndef __LAZY_MATRIX_H__
#define __LAZY_MATRIX_H__

/* Includes **************************************************************************************/
#include <stdio.h>

#include "results.h"
#include "matrix.h"


/* Structs ***************************************************************************************/
/* The matrix dimensions, from its binary format */
typedef struct lazy_matrix_s {
    matrix_header_t header;
    double * current_line;
    FILE * file;
} lazy_matrix_t;


/* Functions Declarations *************************************************************/
/**
 * @purpose Opens a lazy matrix. A lazy matrix is a representation of a matrix file, on which
 *          only 1 line is held in the memory simultaneously.
 *          You can read lines using LAZY_MATRIX_read_next_line
 * @param path A path to a matrix file
 * @param matrix_out The created lazy matrix
 *
 * @return One of result_t values
 *
 * @remark Must be closed using LAZY_MATRIX_close
 */
result_t
LAZY_MATRIX_open(const char * path,
                 lazy_matrix_t **matrix_out);

/**
 * @purpose Closes a lazy matrix
 * @param matrix The matrix to closed
 */
void
LAZY_MATRIX_close(lazy_matrix_t * matrix);


/**
 * @purpose Seeks the matrix's file thus the next call to LAZY_MATRIX_read_next_line will read the
 *          first line.
 *
 * @param matrix Matrix to read
 *
 * @return One of result_t values
 *
 * @remark This function will change the matrix's current_line buffer to store the next line buffer
 */
result_t
LAZY_MATRIX_rewind(lazy_matrix_t * matrix);

/**
 * @purpose Read the next line from the matrix to its internal buffer
 *
 * @param matrix Matrix to read
 *
 * @return One of result_t values
 *
 * @remark This function will change the matrix's current_line buffer to store the next line buffer
 */
result_t
LAZY_MATRIX_read_next_line(lazy_matrix_t * matrix);

/**
 * @purpose Count how many values the matrix has which are non-zero
 *
 * @param matrix Matrix to normalize
 * @param count_out The returned number of nonzero values in the matrix
 *
 * @return One of result_t values
 *
 * @remark matrix must be valid
 */
result_t
LAZY_MATRIX_count_nonzero_values(lazy_matrix_t *matrix, int *count_out);

#endif /* __LAZY_MATRIX_H__ */
