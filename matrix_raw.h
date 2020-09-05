/**
 * @file matrix.h
 * @purpose Reprent a matrix data structure, implement mathematical operations and de/serialization
 */
#ifndef __MATRIX_RAW_H__
#define __MATRIX_RAW_H__

/* Includes **************************************************************************************/
#include "results.h"
#include "common.h"

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
 * @remark matrix must be freed using MATRIX_RAW_free
 */
result_t
MATRIX_RAW_allocate(int rows, int columns, matrix_t ** matrix_out);


#endif /* __MATRIX_RAW_H__ */
