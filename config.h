/**
 * @file config.h
 * @purpose Configuration constants
 */

#ifndef __CONFIG_H__
#define __CONFIG_H__

/* Constants *************************************************************************************/
#include "matrix.h"


/* Constants *************************************************************************************/
#define EPSILON (0.00001)

#ifndef ADJ_MATRIX_TYPE
#define ADJ_MATRIX_TYPE (MATRIX_TYPE_SPMAT_LIST)
#endif /* ADJ_MATRIX_TYPE */

#ifndef MOD_MATRIX_TYPE
#define MOD_MATRIX_TYPE (MATRIX_TYPE_SPMAT_LIST)
#endif /* MOD_MATRIX_TYPE */


#endif /* __CONFIG_H__ */

