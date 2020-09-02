/*
 * @file algorithm1.h
 * @purpose Compute algorithm 1 from the assignment
 */
#ifndef __ALGO1_H__
#define __ALGO1_H__


 /* Includes **************************************************************************************/
#include "matrix.h"
#include "common.h"
#include "eigen.h"
#include "results.h"
#include "spmat.h"

/* Macros ****************************************************************************************/

/**
 * @purpose finding leading eigenvalue of Sparse Matrix
 * @param leading_vector The input leading eigenvector
 * @param prev_vector previous candidate for leading eigenvector from Power Iterations
 * @param eigen_value The calculated leading eigenvalue (output)
 *
 * @return One of result_t values
 *
 * @remark The returned eigen must be freed using MATRIX_free
 */
result_t
calculate_leading_eigenvalue(matrix_t **leading_vector,
					   matrix_t **prev_vector
					   double *eigen_value);

/**
 * @purpose divide a network to two groups
 * @param input modularity Sparse Matrix
 * @param group1 vector representing the first group from division algorithm 1
 * @param group2 vector representing the second group from division algorithm 1
 *
 * @return One of result_t values
 *
 * @remark The returned eigen must be freed using MATRIX_free
 */
result_t
algorithm1(spmat *input,
		matrix_t **group1,
		matrix_t **group2);







#endif /* __ALGO1_H__ */
