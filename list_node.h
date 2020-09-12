/*
 * @file list_node.h
 * @purpose Linked list node
 */
#ifndef __LIST_NODE_H__
#define __LIST_NODE_H__

#include "results.h"

/* Structs ***************************************************************************************/
typedef struct node_s {
    double value;
    int index;
    struct node_s *next;
} node_t;


/* Functions Declarations ************************************************************************/
/**
 * @purpose Create a node_t with the given index and value
 * @param value The value
 * @param index The index
 * @param node_out The created node
 *
 * @return One of result_t values
 */
result_t
LIST_NODE_create(double value, int index, node_t **node_out);

/*
 * @purpose Extract a row with a given index from the Spmat List to a pre-allocated buffer
 * @param A The Spmat List
 * @param row_index The index of the row
 * @param row A pre allocate buffer to store the row. must contain (sizeof(double) * A->n) bytes
 *
 * @return One of result_t values
 */

void
LIST_NODE_destroy(node_t *row);

/**
 * @purpose Perform scalar multiplication with the node and a given vector
 * @param row An SPMat List's line row. Must be valid
 * @param v A buffer (with SPMat List's n values) to multiply with
 * 
 * @return The scalar multiplication result
 */
double
LIST_NODE_scalar_multiply(node_t *row, const double *v);

result_t
LIST_NODE_append(node_t **last_node, double value, int index);


#endif /* __LIST_NODE_H__ */

