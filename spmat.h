/*
 * @file spmat.h
 * @purpose Sparse matrix definition and factory method
 */
#ifndef __SPMAT_H__
#define __SPMAT_H__

/* Constants *************************************************************************************/
#define SPMAT_TYPE_LIST ("-list")
#define SPMAT_TYPE_ARRAY ("-array")


/* Structs ***************************************************************************************/
typedef struct _spmat {
	/* Matrix size (n*n) */
	int	 n;

	/* Adds row i the matrix. Called before any other call,
	 * exactly n times in order (i = 0 to n-1) */
	void (*add_row)(struct _spmat *A, const double *row, int i);

	/* Frees all resources used by A */
	void (*free)(struct _spmat *A);

	/* Multiplies matrix A by vector v, into result (result is pre-allocated) */
	void (*mult)(const struct _spmat *A, const double *v, double *result);

	/* Private field for inner implementation.
	 * Should not be read or modified externally */
	void *private;
} spmat;


/* Functions *************************************************************************************/
result_t
SPMAT_FACTORY_create_spmat(const char * input_matrix_path,
                           const char * spmat_type,
                           spmat **A_out);

#endif /* __SPMAT_H__ */
