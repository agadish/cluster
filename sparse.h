/**
 * @file eigen.h
 * @purpose Implement eigen logic of a mtrix
 */
#ifndef __EIGEN_H__
#define __EIGEN_H__

/* Includes **************************************************************************************/
#include "results.h"
#include "matrix.h"


/* Macros ****************************************************************************************/
/**
 * @purpose Normalize each line of the eigen using the regular scalar multiplication
 * @param eigen Eigen to normalize
 *
 * @return One of result_t values
 */
result_t
EIGEN_normalize(matrix_t *eigen);

/**
 * @remark eigen must be valid
 */
double
EIGEN_calculate_magnitude(matrix_t *eigen);

/**
 * @purpose Checks if all the differences between the eigen' values are smaller than epsilon
 * @param eigen_a The first eigen
 * @param eigen_b The second eigen
 * @param epsilon The largest difference allowed
 *
 * @return TRUE if close enough, otherwise FALSE
 *
 * @remark Eigenes must have the same size
 */
bool_t
EIGEN_is_close(matrix_t * eigen_a, matrix_t * eigen_b, double epsilon);

/**
 * @remark vector_out must be freed later by EIGEN_free
 */
result_t
EIGEN_random_vector(const int length, eigen_t ** vector_out);




#endif /* __EIGEN_H__ */
