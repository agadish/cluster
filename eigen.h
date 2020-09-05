/**
 * @file eigen.h
 * @purpose Implement eigen logic of a mtrix
 */
#ifndef __EIGEN_H__
#define __EIGEN_H__

/* Includes **************************************************************************************/
#include "results.h"
#include "matrix.h"
#include "matrix.h"


/* Macros ****************************************************************************************/

/**
 * @purpose fivides eigenvector as close as epsilon
 * @param spmat The input matrix
 * @param b_vecor A pre allocated b_vector
 * @param eigen_out The calculated eigen value
 *
 * @return One of result_t values
 *
 * @remark The returned eigen must be freed using MATRIX_free
 */
result_t
MATRIX_calculate_eigen(matrix_t *input,
                       const matrix_t *b_vector,
                       matrix_t **eigen_out,
					   matrix_t **prev_vector);

#endif /* __EIGEN_H__ */
