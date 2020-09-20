/**
 * @file adjacency_matrix.h
 * @purpose Reprent a adjacency_matrix data structure.
 *          Implement mathematical operations and de/serialization
 */
#ifndef __ADJACENCY_MATRIX_H__
#define __ADJACENCY_MATRIX_H__

/* Includes **************************************************************************************/
#include "results.h"
#include "common.h"
#include "matrix.h"


/* Macros ****************************************************************************************/
/* Number of elements on the adjacency_matrix */


/* Structs ***************************************************************************************/
/**
 * @brief The adjacency data
 * @param neighbors neighbors buffer length
 * @param neighbors Mapping array from vertice index to its neighbors count
 * @param M The total neighbors count (equals edges count times 2)
 */
typedef struct adjacency_s {
    int n;
    int *neighbors;
    int M;
} adjacency_t;


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
ADJACENCY_MATRIX_open(const char *path,
                      adjacency_t **adj_out,
                      matrix_t **matrix_out);

/**
 * @purpose Free an adjacency matrix which was previously created by
 *          ADJACENCY_MATRIX_open
 *
 * @param matrix The matrix to free
 *
 * @remark Safe to call with NULL
 */
void
ADJACENCY_MATRIX_free(adjacency_t *adjacency);


#endif /* __ADJACENCY_MATRIX_H__ */
