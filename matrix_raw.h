/**
 * @file matrix_raw.h
 * @purpose A matrix implemented using a single array
 */
#ifndef __MATRIX_RAW_H__
#define __MATRIX_RAW_H__

/* Includes **************************************************************************************/
#include "results.h"
#include "common.h"
#include "matrix.h"

/* Macros ****************************************************************************************/
/**
 * @remark Can only be used with MATRIX_TYPE_RAW matrices
 */
#define MATRIX_RAW_GET_ARRAY(matrix) ((double *)((matrix)->private))

/**
 * @remark Can only be used with MATRIX_TYPE_RAW matrices
 */
#define MATRIX_RAW_AT(m, i, j) (MATRIX_RAW_GET_ARRAY(m)[((i) * (m)->n) + (j)])


/* Functions Declarations *************************************************************/
/**
 * @purpose Create a new matrix with the given dimensions
 *
 * @param rows The number of rows in the matrix
 * @param columns The number of columns in the matrix
 * @param matrix_out An external pointer which will point the new matrix
 *
 * @return One of result_t values
 *
 * @remark matrix must be freed using MATRIX_RAW_free
 */
result_t
MATRIX_RAW_allocate(int rows, int columns, matrix_t ** matrix_out);


#endif /* __MATRIX_RAW_H__ */
