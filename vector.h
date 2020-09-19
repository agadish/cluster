/** 
 * @file vector.h
 * @purpose Common vector functions
 */
#ifndef __VECTOR_H__
#define __VECTOR_H__

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
 * Create a random vector
 *
 * @param length The vector's length
 * @param vector A pre-allocated lenth-sized vector
 */
void
VECTOR_random_vector(size_t length, double *vector);


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

/**
 * @brief Convert an s-vector to s-index vector
 *
 * Given an s-vector that maps node n to a group (-1.0 or 1.0), this
 * function will calculate a s-index vectro that maps a node n to its index
 * within its new group.
 * Both groups' new indexes will be written to the one array.
 *
 * @param matrix The matrix to divide
 * @param s_vector The given s-vector
 * @param length The length of the s-vector
 * @param s_indexes A pre-allocated buffer which the indexes of both groups
 *                  will be written to
 * @param matrix1_n_out Will return the 1st group's lenght (=number of 1.0's)
 *
 * @return One of result_t values
 *
 * @remark vector_s and s_indexes must be valid vectors with the given length
 */
int
VECTOR_create_s_indexes(const double * vector_s,
                        int length,
                        int *s_indexes);

#endif /* __VECTOR_H__ */
