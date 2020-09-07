/*
 * @file spmat_list.h
 * @purpose Sparse matrix implemented using linked lists
 */
#ifndef __SPMAT_LIST_H__
#define __SPMAT_LIST_H__

/* Includes **************************************************************************************/
#include <stddef.h>

#include "matrix.h"


/* Functions Declarations ************************************************************************/
/* Allocates a new linked-lists sparse matrix of size n */
result_t
SPMAT_LIST_allocate(int n, matrix_t **mat);


result_t
SPMAT_LIST_divide_matrix(const matrix_t *matrix,
                         const double * vector_s,
                         matrix_t **matrix1_out,
                         matrix_t **matrix2_out,
                         size_t ones_count);



#endif /* __SPMAT_LIST_H__ */
