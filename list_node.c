/*
 * @file list_node.c
 * @purpose Linked list node
 */

/* Includes **************************************************************************************/
#include <stdlib.h>
#include <stddef.h>

#include "common.h"
#include "results.h"
#include "list_node.h"


/* Functions *************************************************************************************/
result_t
LIST_NODE_create(double value, int index, node_t **node_out)
{
    result_t result = E__UNKNOWN;
    node_t *node = NULL;

    node = (node_t *)malloc(sizeof(*node));
    if (NULL == node) {
        result = E__MALLOC_ERROR;
        goto l_cleanup;
    }
    node->value = value;
    node->index = index;
    node->next = NULL;

    *node_out = node;

    result = E__SUCCESS;
l_cleanup:
    if (E__SUCCESS != result) {
        FREE_SAFE(node);
    }

    return result;
}

void
LIST_NODE_destroy(node_t *list)
{
    node_t *prev_node = NULL;
    node_t *node = list;

    while (NULL != node) {
        prev_node = node;
        node = node->next;
        FREE_SAFE(prev_node);
    }
}

double
LIST_NODE_scalar_multiply(node_t *row, const double *v)
{
    double multiplication = 0.0;
    node_t * scanner = NULL;
    for (scanner = row ; NULL != scanner ; scanner = scanner->next) {
        multiplication += (scanner->value * v[scanner->index]);
    }

    return multiplication;
}


result_t
LIST_NODE_append(node_t **last_node, double value, int index)
{
    result_t result = E__UNKNOWN;
    node_t *new_last_node = NULL;

    if (NULL == last_node)
    {
        result = E__NULL_ARGUMENT;
        goto l_cleanup;
    }

    result = LIST_NODE_create(value, index, &new_last_node);
    if (E__SUCCESS != result) {
        goto l_cleanup;
    }

    /* Set as the last node's next, if exists */
    if (NULL != (*last_node)) {
        /* First node in list */
        (*last_node)->next = new_last_node;
    }

    *last_node = new_last_node;

    result = E__SUCCESS;
l_cleanup:

    if (E__SUCCESS != result) {
        LIST_NODE_destroy(new_last_node);
        new_last_node = NULL;
    }

    return result;
}

result_t
LIST_NODE_range(size_t count, node_t **node_out)
{
    result_t result = E__UNKNOWN;
    node_t *head = NULL;
    node_t *tail = NULL;
    size_t i = 0;

    if (NULL == node_out) {
        result = E__NULL_ARGUMENT;
    }

    for (i = 0 ; i < count ; ++i) {
        result = LIST_NODE_append(&tail, 0.0, i);
        if (E__SUCCESS != result) {
            goto l_cleanup;
        }
        
        if (NULL == head) {
            head = tail;
        }
    }

    *node_out = head;
    /* Success */
    result = E__SUCCESS;
l_cleanup:
    if (E__SUCCESS != result) {
        LIST_NODE_destroy(head);
        head = NULL;
    }

    return result;
}
