/**
 * @file adjacency_matrix.h
 * @purpose Reprent a adjacency_matrix data structure, implement mathematical operations and de/serialization
 */
#ifndef __ADJACENCY_MATRIX_H__
#define __ADJACENCY_MATRIX_H__

/* Includes **************************************************************************************/
#include "results.h"
#include "common.h"


/* Macros ****************************************************************************************/
/* Number of elements on the adjacency_matrix */


/* Structs ***************************************************************************************/
/* The adjacency_matrix dimensions and data */
typedef struct adjacency_matrix_s {
    matrix_t *matrix;
    int * neighbors;
    int M;
} adjacency_matrix_t;


/* Functions Declarations *************************************************************/
/*
 * @purpose Create a new adjacency matrix by a given encoded adj-matrix file path
 *
 * @param path The path to the adjacency matrix file. must be valid!
 * @param adj_matrix_out An external pointer which will point the new matrix
 *
 * @return One of result_t values
 *
 * @remark adj_matrix must be freed using ADJACENCY_MATRIX_free
 */
result_t
ADJACENCY_MATRIX_open(const char * path, adjacency_matrix_t **adj_matrix_out);

/**
 * @purpose Free an adjacency matrix which was previously created by
 *          ADJACENCY_MATRIX_open
 *
 * @param matrix The matrix to free
 *
 * @remark Safe to call with NULL
 */
void
ADJACENCY_MATRIX_free(adjacency_matrix_t *adjacency_matrix);

/**
 * @purpose Calculate the modularity matrix of a given adjacency matrix
 */
result_t
ADJACENCY_MATRIX_calculate_modularity(adjacency_matrix_t *adj, matrix_t **mod_matrix_out);


#endif /* __ADJACENCY_MATRIX_H__ */
