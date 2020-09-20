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

/**
 * @purpose calculate magnitude of vector
 * @param vector input
 * @param length The vector length
 *
 * @remark vector must be valid buffer, value must be non-zero
 */
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

double
VECTOR_scalar_multiply_with_s(const double * l1, const double * s, size_t n)
{
    double result = 0.0;
    const double * l1_end = l1 + n;

#ifdef OPTIMIZE_VECTOR_OPERATIONS
    for ( ; l1 < l1_end - 1 ; l1 += 2, s += 2) {
        result += ((l1[0] * ((s[0] > 0) ? 1 : -1)) +
                (l1[1] * ((s[1] > 0) ? 1 : -1)));
    }
#endif /* OPTIMIZE_VECTOR_OPERATIONS */

    /* If OPTIMIZE_VECTOR_OPERATIONS is defined, this will run up to 1 iteration */
    for ( ; l1 < l1_end ; ++l1, ++s) {
        result += (l1[0] * (s[0] > 0 ? 1 : -1));
    }

    return result;
}

int
VECTOR_scalar_multiply_int_with_s(const int * l1, const double * s, size_t n)
{
    int result = 0.0;
    const int * l1_end = l1 + n;

#ifdef OPTIMIZE_VECTOR_OPERATIONS
    for ( ; l1 < l1_end - 1 ; l1 += 2, s += 2) {
        result += ((l1[0] * (s[0] > 0 ? 1 : -1)) + (l1[1] * (s[1]) > 0 ? 1 : -1));
    }
#endif /* OPTIMIZE_VECTOR_OPERATIONS */

    /* If OPTIMIZE_VECTOR_OPERATIONS is defined, this will run up to 1 iteration */
    for ( ; l1 < l1_end ; ++l1, ++s) {
        result += (l1[0] * (s[0] > 0 ? 1 : -1));
    }

    return result;
}

void
VECTOR_random_vector(size_t length, double *vector)
{
    size_t i = 0;
    int random_int = 0;
    srand(time(NULL));

    for (i = 0 ; i < length ; ++i) {
        random_int = (rand() % 1000);
        vector[i] = (double)random_int;
    }
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


int
VECTOR_create_s_indexes(const double * vector_s,
                        int length,
                        int *s_indexes)
{
    int i = 0;
    int index_1 = 0;
    int index_2 = 0;

    /* Go over the s-vector */
    for (i = 0 ; i < length ; ++i) {
        if (1 == vector_s[i]) {
            s_indexes[i] = index_1;
            ++index_1;
        } else if (-1 == vector_s[i]) {
            s_indexes[i] = index_2;
            ++index_2;
        } 
    }

    return index_1;
}
