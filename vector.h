/** 
 * @file vector.h
 * @purpose Common vector functions
 */
#ifndef __VECTOR_H__
#define _VECTOR_H__

#include "results.h"

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
VECTOR_scalar_multiply(const double * l1, const double * l2, int n);

/**
 * @remark vector_out must be freed later by MATRIX_free
 */
result_t
VECTOR_random_vector(const int length, double ** vector_out);

#endif /* __VECTOR_H__ */
