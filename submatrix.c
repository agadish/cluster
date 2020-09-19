/**
 * @file submatrix.c
 * @purpose A generic submatrix
 */

/* Includes **************************************************************************************/
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


/* Functions ************************************************************************************/
result_t
SUBMATRIX_create(const adjacency_matrix_t *adj,
                 submatrix_t **smat_out)
{
    result_t result = E__UNKNOWN;
    int *g = NULL;
    int g_length = 0;
    submatrix_t *smat = NULL;

    if ((NULL == adj) || (NULL == g) || (NULL == smat_out))
    {
        result = E__NULL_ARGUMENT;
        goto l_cleanup;
    }

    if (0 > g_length) {
        result = E__INVALID_SIZE;
        goto l_cleanup;
    }

    /* 1. Allocate g-vector with length n */
    g_length = adj->original->n;
    g = (int *)malloc(g_length * sizeof(*g));
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
    smat->g_length = g_length;

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
        smat->g_length = 0;
        FREE_SAFE(smat->g);
    }
    FREE_SAFE(smat);
}
