/*
 * @file spmat_list.h
 * @purpose Sparse matrix implemented using linked lists
 */
#ifndef __SPMAT_LIST_H__
#define __SPMAT_LIST_H__

/* Includes **************************************************************************************/
#include "matrix.h"

/* Functions Declarations ************************************************************************/
/* Allocates a new linked-lists sparse matrix of size n */
result_t
SPMAT_LIST_allocate(int n, matrix_t **mat);


#endif /* __SPMAT_LIST_H__ */
