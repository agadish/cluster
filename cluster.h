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
#include "division_file.h"


/* Functions Declarations ************************************************************************/
result_t
CLUSTER_divide_repeatedly(matrix_t *matrix, division_file_t *output_file);


#endif /* __CLUSTER_H__ */
