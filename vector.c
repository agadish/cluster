/** 
 * @file vector.c
 * @purpose Common vector functions
 */

/* Includes **************************************************************************************/
#include <time.h>
#include <stdlib.h>

#include "vector.h"
#include "results.h"
#include "common.h"

#define OPTIMIZE_SCALAR_MULTIPLICATION


/* Functions *************************************************************************************/
double
VECTOR_scalar_multiply(const double * l1, const double * l2, int n)
{
    double result = 0.0;
    const double * l1_end = l1 + n;

#ifdef OPTIMIZE_SCALAR_MULTIPLICATION
    for ( ; l1 < l1_end - 1 ; l1 += 2, l2 += 2) {
        result += ((l1[0] * l2[0]) + (l1[1] * l2[1]));
    }
#endif /* OPTIMIZE_SCALAR_MULTIPLICATION */

    /* If OPTIMIZE_SCALAR_MULTIPLICATION is defined, this will run up to 3 iterations */
    for ( ; l1 < l1_end ; ++l1, ++l2) {
        result += (l1[0] * l2[0]);
    }

    return result;
}

result_t
VECTOR_random_vector(const int length, double ** vector_out)
{
    result_t result = E__UNKNOWN;
    double * vector = NULL;
    int i = 0;
    int random_int = 0;

    if ((NULL == vector_out) || (0 > length)) {
        result = E__NULL_ARGUMENT;
        goto l_cleanup;
    }

    srand(time(NULL));

    vector = (double *)malloc(sizeof(*vector) * length);
    if (NULL == vector) {
        result = E__MALLOC_ERROR;
        goto l_cleanup;
    }

    for (i = 0 ; i < length ; ++i) {
        random_int = (rand() % 1000);
        vector[i] = (double)random_int;
    }

    *vector_out = vector;

    result = E__SUCCESS;
l_cleanup:
    if (E__SUCCESS != result) {
        FREE_SAFE(vector);
    }

    return result;
}

