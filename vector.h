/** 
 * @file vector.h
 * @purpose Common vector functions
 */
#ifndef __VECTOR_H__
#define _VECTOR_H__

/* Includes **************************************************************************************/
#include <stddef.h>

#include "results.h"
#include "common.h"


/* Functions Declarations ************************************************************************/
/**
 * @purpose Calculate the scalar multiplication between two vectors
 *
 * @param l1 First vector - must be valid n-sized double array!
 * @param l2 Second vector - must be valid n-sized double array!
 * @param n The length of the vectors
 *
 * @return The scalar multiplication result
 * @remark The vectors must be valid
 */
double
VECTOR_scalar_multiply(const double * l1, const double * l2, size_t n);

/**
 * @remark vector_out must be freed later by MATRIX_free
 */
result_t
VECTOR_random_vector(size_t length, double ** vector_out);


/**
 * @purpose Normalize a vector using scalar multiplication
 * @param vector The vector to normalize
 * @param length The vector's length
 *
 * @return One of result_t values
 */
result_t
VECTOR_normalize(double *vector, size_t length);

/**
 * @purpose Checks if all the differences between the matrix' values are smaller than epsilon
 * @param vector_a The first vector
 * @param vector_b The second vector
 * @param length The vectors' length
 * @param epsilon The largest difference allowed
 *
 * @return TRUE if close enough, otherwise FALSE
 *
 * @remark Matrixes must have the same size
 */
bool_t
VECTOR_is_close(const double * vector_a,
                const double * vector_b,
                size_t length,
                double epsilon);

#endif /* __VECTOR_H__ */
