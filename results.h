/**
 * @file results.h
 * @purpose Define the return values of the project
 */

#ifndef __RESULTS_H__
#define __RESULTS_H__

typedef enum result_e { 
    E__UNKNOWN = -1, 
    E__SUCCESS = 0, 
    E__NULL_ARGUMENT, 
    E__FOPEN_ERROR, 
    E__FREAD_ERROR, 
    E__FWRITE_ERROR, 
    E__FSEEK_ERROR, 
    E__MALLOC_ERROR, 
    E__INVALID_SIZE, 
    E__INVALID_CMDLINE_ARGS, 
    E__ZERO_DIV_ERROR, 
    E__INVALID_ROW_INDEX, 
    E__UNKNOWN_MATRIX_IMPLEMNTATION,
    E__INVALID_MATRIX_TYPE,
    E__INVALID_S_VECTOR,
    E__ROW_ALREADY_IN_USE,
} result_t; 

#endif /* __RESULTS_H__ */
