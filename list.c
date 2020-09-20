/*
 * @file list.c
 * @purpose Linked list node
 */

/* Includes **************************************************************************************/
#include <stdlib.h>
#include <stddef.h>

#include "common.h"
#include "results.h"
#include "list.h"


/* Functions Declarations ****************************************************/
static
void
node_link(node_t *first, node_t *second);


/* Functions *****************************************************************/
static
void
node_link(node_t *first, node_t *second)
{
    if (NULL != first) {
        first->next = second;
    }

    if (NULL != second) {
        second->prev = first;
    }
}

result_t
LIST_create(list_t **list_out)
{
    result_t result = E__UNKNOWN;
    list_t *list = NULL;

    list = (list_t *)malloc(sizeof(*list));
    if (NULL == list) {
        result = E__MALLOC_ERROR;
        goto l_cleanup;
    }

    list->first = NULL;
    list->last = NULL;

    *list_out = list;

    result = E__SUCCESS;
l_cleanup:
    if (E__SUCCESS != result) {
        FREE_SAFE(list);
    }

    return result;
}

void
LIST_destroy(list_t *list)
{
    node_t *prev_node = NULL;
    node_t *node = NULL;

    if (NULL != list) {
        node = list->first;
        while (NULL != node) {
            prev_node = node;
            node = node->next;
            FREE_SAFE(prev_node);
        }
        FREE_SAFE(list);
    }
}

result_t
LIST_insert(list_t *list, node_t *next_node, double value, int index)
{
    result_t result = E__UNKNOWN;
    node_t *new_node = NULL;
    node_t *prev_node = NULL;

    if (NULL == list)
    {
        result = E__NULL_ARGUMENT;
        goto l_cleanup;
    }

    new_node = (node_t *)malloc(sizeof(*new_node));
    if (NULL == new_node) {
        result = E__MALLOC_ERROR;
        goto l_cleanup;
    }
    new_node->value = value;
    new_node->index = index;

    if (NULL == next_node) {
        prev_node = list->last;
        list->last = new_node;
    } else {
        prev_node = next_node->prev;
    }

    if (NULL == prev_node) {
        list->first = new_node;
    }

    node_link(prev_node, new_node);
    node_link(new_node, next_node);


    result = E__SUCCESS;
l_cleanup:

    if (E__SUCCESS != result) {
        LIST_destroy(list);
        new_node = NULL;
    }

    return result;
}

double
LIST_scalar_multiply(list_t *list, const double *v)
{
    double multiplication = 0.0;
    node_t * scanner = NULL;
    for (scanner = list->first ; NULL != scanner ; scanner = scanner->next) {
        multiplication += (scanner->value * v[scanner->index]);
    }

    return multiplication;
}


result_t
LIST_range(size_t count, list_t **list_out)
{
    result_t result = E__UNKNOWN;
    list_t *list = NULL;
    size_t i = 0;

    if (NULL == list_out) {
        result = E__NULL_ARGUMENT;
    }

    result = LIST_create(&list);
    if (E__SUCCESS != result) {
        goto l_cleanup;
    }

    for (i = 0 ; i < count ; ++i) {
        result = LIST_insert(list, NULL, 0.0, i);
        if (E__SUCCESS != result) {
            goto l_cleanup;
        }
    }

    /* Success */
    *list_out = list;

    result = E__SUCCESS;
l_cleanup:

    if (E__SUCCESS != result) {
        LIST_destroy(list);
        list = NULL;
    }

    return result;
}

result_t
LIST_remove_node(list_t *list, node_t *node)
{
    result_t result = E__UNKNOWN;
    node_t *prev = NULL;
    node_t *next = NULL;

    if ((NULL == list) || (NULL == node)) {
        result = E__NULL_ARGUMENT;
        goto l_cleanup;
    }

    prev = node->prev;
    next = node->next;

    if (list->first == node) {
        list->first = next;
    }
    if (list->last == node){ 
        list->last = prev;
    }

    node_link(prev, next);

    FREE_SAFE(node);

    result = E__SUCCESS;
l_cleanup:

    return result;
}
