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

/* Functions Declarations *************************************************************/
/**
 * @purpose Create a new matrix with the given dimensions
 *
 * @param n The number of rows/columns of the matrix
 * @param matrix_out An external pointer which will point the new matrix
 *
 * @return One of result_t values
 *
 * @remark matrix must be freed using MATRIX_RAW_free
 */
result_t
MATRIX_RAW_allocate(int n, matrix_t ** matrix_out);

/** XXX: The following line DOESNT cause error although the implementaiton is different 
MATRIX_RAW_allocate(int rows, int columns, matrix_t ** matrix_out); */


#endif /* __MATRIX_RAW_H__ */
