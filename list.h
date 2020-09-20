/*
 * @file list.h
 * @purpose Linked list node
 */
#ifndef __LIST_H__
#define __LIST_H__

#include <stddef.h>

#include "results.h"

/* Structs ***************************************************************************************/
typedef struct node_s node_t;

typedef struct list_s {
    struct node_s *first;
    struct node_s *last;
} list_t;

struct node_s {
    double value;
    int index;
    struct node_s *prev;
    struct node_s *next;
};


/* Functions Declarations ************************************************************************/
/**
 * @purpose Create an empty list_t
 * @param list_out The new list
 *
 * @return One of result_t values
 */
result_t
LIST_create(list_t **list_out);

/*
 * @purpose Destroy a list_t
 * @param list The list to destroy
 *
 * @return One of result_t values
 */
void
LIST_destroy(list_t *list);

/**
 * Create a new node with the given value and index.
 * The node will be inserted before next_node, or if NULL is given it will be
 * appended to the end of the list
 *
 * @param list The list
 * @param node_next An existing node which will the the new one's next
 * @param value The value
 * @param index The index
 *
 * @return One of result_t values
 */
result_t
LIST_insert(list_t *list, node_t *next_node, double value, int index);

/*
 * Create a list from 0 (including) to count - 1 (including)
 */
result_t
LIST_range(size_t count, list_t **list_out);

/**
 * @purpose Perform scalar multiplication with the node and a given vector
 * @param row A linked list. Nodes' indexes must not overflow v length
 * @param v A vector
 * 
 * @return The scalar multiplication result
 */
double
LIST_scalar_multiply(list_t *row, const double *v);

/**
 * @purpose remove a node from a linked list
 * @param list - the linked list
 * @param node - the node we want to delete
 *
 * @return The scalar multiplication result
 */
result_t
LIST_remove_node(list_t *list, node_t *node);


#endif /* __LIST_H__ */

