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
 * @remark On error, the file won't be created. If already exists it may be overriden
 */
result_t
DIVISION_FILE_open(const char *path, division_file_t **division_file_out);

result_t
DIVISION_FILE_write_matrix(division_file_t *division_file,
                           const int *indexes,
                           int length);

result_t
DIVISION_FILE_finalize(division_file_t *division_file);

void
DIVISION_FILE_close(division_file_t *division_file);




#endif /* __DIVISION_FILE_H__ */

