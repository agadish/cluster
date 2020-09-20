/**
 * @file debug
 * @purpose Debug macros
 */

#ifndef __DEBUG_H__
#define __DEBUG_H__

#include <stdio.h>

#ifdef __DEBUG__
#define DEBUG_PRINT(...) do {                                               \
    (void)fprintf(stderr, "%s:%d:%s: ", __FILE__, __LINE__, __func__);      \
    (void)fprintf(stderr, __VA_ARGS__);                                     \
    (void)fprintf(stderr, "\n");                                            \
} while (0)
#else
#define DEBUG_PRINT(...) 
#endif

#endif /* __DEBUG_H__ */

