/*
 * @file common.h
 * @purpose Commonly used constants and macros
 */
#ifndef __COMMON_H__
#define __COMMON_H__
/* Includes **************************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <stddef.h>


/* Constants *************************************************************************************/
#define FALSE (0)
#define TRUE (!FALSE)


/* Typedefs **************************************************************************************/
typedef uint8_t bool_t;


/* Macros ****************************************************************************************/
/** Close a FILE* and set it as NULL */
#define FCLOSE_SAFE(f) do {     \
    if (NULL != (f)) {          \
        (void)fclose((f));      \
        (f) = NULL;             \
    }                           \
} while (0)

/** Free a pointer and set it as NULL */
#define FREE_SAFE(p) do {       \
    free((p));                  \
    (p) = NULL;                 \
} while (0)

#define UNUSED_ARG(a) ((void)(a))

#define MAX(a, b) (((a) > (b)) ? (a) : (b))


#endif /* __COMMON_H__ */
