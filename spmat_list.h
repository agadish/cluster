/*
 * @file spmat_list.h
 * @purpose Sparse matrix implemented using linked lists
 */
#ifndef __SPMAT_LIST_H__
#define __SPMAT_LIST_H__

/* Includes **************************************************************************************/
#include <stddef.h>
#include <stdio.h>

#include "matrix.h"
#include "common.h"


/* Functions Declarations ************************************************************************/
/* Allocates a new linked-lists sparse matrix of size n */
result_t
SPMAT_LIST_allocate(int n, matrix_t **mat);

void
SPMAT_LIST_print(const char *matrix_name, matrix_t *mat_in);

#endif /* __SPMAT_LIST_H__ */
