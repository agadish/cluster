/**
 * @file submatrix.c
 * @purpose A generic submatrix
 */

/* Includes ******************************************************************/
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

#include "submatrix.h"
#include "common.h"
#include "spmat_list.h"
#include "spmat_array.h"


/* Functions *****************************************************************/
result_t
SUBMATRIX_create(const adjacency_t *adj,
                 matrix_t *matrix,
                 submatrix_t **smat_out)
{
    result_t result = E__UNKNOWN;
    int *g = NULL;
    submatrix_t *smat = NULL;

    if ((NULL == adj) || (NULL == smat_out))
    {
        result = E__NULL_ARGUMENT;
        goto l_cleanup;
    }

    /* 1. Allocate g-vector with length n */
    g = (int *)malloc(adj->n * sizeof(*g));
    if (NULL == g) {
        result = E__MALLOC_ERROR;
        goto l_cleanup;
    }   

    smat = (submatrix_t *)malloc(sizeof(*smat));
    if (NULL == smat) {
        result = E__MALLOC_ERROR;
        goto l_cleanup;
    }

    smat->adj = adj;
    smat->g = g;
    smat->g_length = 0;
    smat->add_to_diag = 0.0;
    smat->orig = matrix;

    /* result = SPMAT_LIST_transpose(orig, &smat->transposed); */
    /* if (E__SUCCESS != result) { */
    /*     goto l_cleanup; */
    /* } */

    *smat_out = smat;

    result = E__SUCCESS;
l_cleanup:
    if (E__SUCCESS != result) {
        /* On failure free the smpat */
        SUBMATRIX_FREE_SAFE(smat);
    }

    return result;
}

void
SUBMATRIX_free(submatrix_t *smat)
{
    if (NULL != smat) {
        FREE_SAFE(smat->g);
        MATRIX_FREE_SAFE(smat->orig);
        smat->g_length = 0;
        FREE_SAFE(smat->orig);
    }
    FREE_SAFE(smat);
}
