/*
 * @file algorithm1.h
 * @purpose Compute algorithm 1 from the assignment
 */
#ifndef __CLUSTER_H__
#define __CLUSTER_H__


 /* Includes **************************************************************************************/
#include "matrix.h"
#include "common.h"
#include "results.h"


/* Functions Declarations ************************************************************************/
/**
 * @purpose divide a network to two groups
 * @param input Matrix to divide
 * @param group1 vector representing the first group from division algorithm 1
 * @param group2 vector representing the second group from division algorithm 1
 *
 * @return One of result_t values
 *
 * @remark The returned eigen must be freed using MATRIX_free
 */
result_t
CLUSTER_divide(matrix_t *input,
               matrix_t **group1_out,
               matrix_t **group2_out);

#if 0
result_t
CLUSTER_divide_repeatedly(matrix_t *matrix);
#endif

result_t
CLUSTER_divide_alg4(matrix_t *matrix,
                    matrix_t **group1_out,
                    matrix_t **group2_out);


#endif /* __CLUSTER_H__ */
