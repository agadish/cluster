/** 
 * @file vector.c
 * @purpose Common vector functions
 */

/* Includes **************************************************************************************/
#include <time.h>
#include <stdlib.h>
#include <math.h>

#include "vector.h"
#include "results.h"
#include "common.h"

#define OPTIMIZE_VECTOR_OPERATIONS

/* Functions Declarations ************************************************************************/
static
double
vector_calculate_magnitude(double *vector, size_t length);

/**
 * @purpose divides vector by scalar umber
 * @param vector The vector to divide
 * @param length The vector length
 * @param value The value to divide with
 *
 * @remark vector must be valid buffer, value must be non-zero
 */
static
void
vector_div(double *vector, size_t length, double value);


/* Functions *************************************************************************************/
double
VECTOR_scalar_multiply(const double * l1, const double * l2, size_t n)
{
    double result = 0.0;
    const double * l1_end = l1 + n;

#ifdef OPTIMIZE_VECTOR_OPERATIONS
    for ( ; l1 < l1_end - 1 ; l1 += 2, l2 += 2) {
        result += ((l1[0] * l2[0]) + (l1[1] * l2[1]));
    }
#endif /* OPTIMIZE_VECTOR_OPERATIONS */

    /* If OPTIMIZE_VECTOR_OPERATIONS is defined, this will run up to 1 iteration */
    for ( ; l1 < l1_end ; ++l1, ++l2) {
        result += (l1[0] * l2[0]);
    }

    return result;
}

result_t
VECTOR_random_vector(size_t length, double ** vector_out)
{
    result_t result = E__UNKNOWN;
    double * vector = NULL;
    size_t i = 0;
    int random_int = 0;

    if (NULL == vector_out) {
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

static
double
vector_calculate_magnitude(double *vector, size_t length)
{
    double norm_square = 0.0;
    double magnitude = 0.0;

    norm_square = VECTOR_scalar_multiply(vector, vector, length);
    magnitude = sqrt(norm_square);

    return magnitude;
}

static
void
vector_div(double *vector, size_t length, double value)
{
    double * i = NULL;
    double * vector_end = NULL;

    /* 1. Initialize vector start and end */
    i = vector;
    vector_end = vector + length;

    /* 2. Divide using optimization */
#ifdef OPTIMIZE_VECTOR_OPERATIONS
    for ( ; i < vector_end - 1; i += 2) {
        i[0] /= value;
        i[1] /= value;
    }
#endif /* OPTIMIZE_VECTOR_OPERATIONS */

    for (; i < vector_end ; ++i) {
        i[0] /= value;
    }
}

result_t
VECTOR_normalize(double *vector, size_t length)
{
    result_t result = E__UNKNOWN;
    double magnitude = 0.0;

    if (NULL == vector) {
        result = E__NULL_ARGUMENT;
        goto l_cleanup;
    }

    magnitude = vector_calculate_magnitude(vector, length);
    vector_div(vector, length, magnitude);

    result = E__SUCCESS;
l_cleanup:

    return result;
}

bool_t
VECTOR_is_close(const double * vector_a,
                const double * vector_b,
                size_t length,
                double epsilon)
{
    bool_t result = TRUE;
    size_t row = 0;

    for (row = 0 ; row < length ; ++row) {
        if (fabs(vector_a[row] - vector_b[row]) >= epsilon) {
            result = FALSE;
            break;
        }
    }
    return result;
}

