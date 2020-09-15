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
result_t
CLUSTER_divide_repeatedly(matrix_t *matrix);

#if 0
result_t
CLUSTER_divide_repeatedly(matrix_t *matrix);
#endif

result_t
CLUSTER_divide_alg4(matrix_t *matrix,
                    matrix_t **group1_out,
                    matrix_t **group2_out);


#endif /* __CLUSTER_H__ */
