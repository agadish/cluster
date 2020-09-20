/**
 * @file division_file.h
 * @purpose The output file of the spmat division
 */

#ifndef __DIVISION_FILE_H__
#define __DIVISION_FILE_H__

/** Includes *************************************************************************************/
#include "results.h"
#include "matrix.h"


/** Structs **************************************************************************************/
typedef struct division_file_s division_file_t;


/** Functions Declarations ***********************************************************************/
/**
 * @purpose open file
 * @param path- path of input file
 * @param division_file_out path of output file
 *
 * @return one of retrun_t values
 */
result_t
DIVISION_FILE_open(const char *path, division_file_t **division_file_out);

/**
 * @purpose write a matrix to the output file
 * @param division_file- path of output file
 * @param indexes - indexes describing written matrix
 * @param length - matrix length
 *
 * @return one of retrun_t values
 */
result_t
DIVISION_FILE_write_matrix(division_file_t *division_file,
                           const int *indexes,
                           int length);

/**
 * @purpose finalize the output file to include correct values
 * @param division_file- path of output file
 * @return one of retrun_t values
 */
result_t
DIVISION_FILE_finalize(division_file_t *division_file);

/**
 * @purpose close the output file after writting it
 * @param division_file- path of output file
 * @return one of retrun_t values
 */
void
DIVISION_FILE_close(division_file_t *division_file);




#endif /* __DIVISION_FILE_H__ */

