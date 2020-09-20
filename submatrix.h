/**
 * @file submatrix.h
 * @purpose Abstract submatrix class: fields and virtual table
 */
#ifndef __SUBMATRIX_H__
#define __SUBMATRIX_H__

/* Includes ******************************************************************/
#include <stddef.h>

#include "results.h"
#include "common.h"
#include "adjacency_matrix.h"


/* Macros ********************************************************************/
/* Number of elements on the submatrix */
#define SUBMATRIX_COUNT(m) ((size_t)((m)->n * (m)->n))

/* Vtable macros */
#define SUBMATRIX_VTABLE(m) ((m)->vtable)
#define SUBMATRIX_FREE(m) (SUBMATRIX_free((m)))

#define SUBMATRIX_FREE_SAFE(m) do { \
    if (NULL != (m)) {              \
        SUBMATRIX_FREE(m);          \
        (m) = NULL;                 \
    }                               \
} while (0)


/* Typedefs ******************************************************************/
typedef struct submatrix_s submatrix_t;


/* Structs *******************************************************************/

/*
 * A representation of a submatrix given the whole matrix and subindexes.
 * The submatrix has virtual values, given by the following furmula:
 *
 * [ ^  ]     [    ]     ki*kj            (  g )             
 * | B  |  =  | A  |  -  -----  - delta * ( f  )  +  delta *  add_to_diag 
 * [ ij ]     [ ij ]       M         ij   (  i )        ij
 */
struct submatrix_s {
    const adjacency_t *adj;
    matrix_t *orig;
    /* Count of subindexes */
    int g_length;
    /* Array of subindexes of this matrix */
    int *g;
    /* A constant value to add to the diag */
    double add_to_diag;
};


/* Functions Declarations ****************************************************/
/* 
 * @purpose Create a submatrix that wraps a given matrxi
 * 
 * @param full_matrix 
 * @param type The submatrix's implementation
 * @param submatrix_out The submatrix creadet
 *
 * @return One of result_t values
 */
result_t
SUBMATRIX_create(const adjacency_t *adj,
                 matrix_t *matrix,
                 submatrix_t **smat_out);

result_t
SUBMATRIX_update(submatrix_t *submatrix, int *g, int g_length);

/*
 * @remark The original, transpoed and g-vector are not freed!
 */
void
SUBMATRIX_free(submatrix_t *submatrix_out);


#endif /* __SUBMATRIX_H__ */
