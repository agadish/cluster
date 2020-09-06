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


/* Functions *************************************************************************************/
/**
 * @purpose fivides eigenvector as close as epsilon
 * @param spmat The input matrix
 * @param b_vecor A pre initialized b_vector. Will be modified!
 * @param eigen_out The calculated eigen value
 *
 * @return One of result_t values
 *
 * @remark The returned eigen must be freed using MATRIX_free
 * @remark v_vector will be modified
 */
result_t
MATRIX_calculate_eigen(const matrix_t *input,
                       double *b_vector,
                       double **eigen_out);

#endif /* __EIGEN_H__ */
