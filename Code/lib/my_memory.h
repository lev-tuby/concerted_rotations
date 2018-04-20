#ifndef HDR_MY_MEM
#define HDR_MY_MEM
/**
 * @file
 * @brief Header file contain functions for memory allocation/deallocation
 */

#include "messages.h"
#include <stdio.h>

/** @brief Allocation of 1D FILE array */
FILE **F2t(int n1);
/** @brief Allocation of 2D FILE array */
FILE ***F3t(int n1, int n2);


/** @brief Allocation of 1D char array */
char *c1t(int n1);
/** @brief Allocation of 2D char array */
char **c2t(int n1, int n2);
/** @brief Allocation of 3D char array */
char ***c3t(int n1, int n2, int n3);
/** @brief Allocation of 4D char array */
char ****c4t(int n1, int n2, int n3, int n4);
/** @brief Deallocation of #c4t() */
void free_c4t(char ****p);
/** @brief Deallocation of #c3t() */
void free_c3t(char ***p);
/** @brief Deallocation of #c2t() */
void free_c2t(char **p);
/** @brief Deallocation of #c1t() */
void free_c1t(char *p);


/** @brief Allocation of 1D short array */
short *s1t(int n1);
/** @brief Allocation of 2D short array */
short **s2t(int n1, int n2);
/** @brief Allocation of 3D short array */
short ***s3t(int n1, int n2, int n3);
/** @brief Allocation of 4D short array */
short ****s4t (int n1, int n2, int n3, int n4);
/** @brief Deallocation of #s4t() */
void free_s4t(short ****p);
/** @brief Deallocation of #s3t() */
void free_s3t(short ***p);
/** @brief Deallocation of #s2t() */
void free_s2t(short **p);
/** @brief Deallocation of #s1t() */
void free_s1t(short *p);


/** @brief Allocation of 1D int array */
int *i1t(int n1);
/** @brief Allocation of 2D int array */
int **i2t(int n1, int n2);
/** @brief Allocation of 3D int array */
int ***i3t(int n1, int n2, int n3);
/** @brief Allocation of 4D int array */
int ****i4t (int n1, int n2, int n3, int n4);
/** @brief Deallocation of #i4t() */
void free_i4t(int ****p);
/** @brief Deallocation of #i3t() */
void free_i3t(int ***p);
/** @brief Deallocation of #i2t() */
void free_i2t(int **p);
/** @brief Deallocation of #i1t() */
void free_i1t(int *p);


/** @brief Allocation of 1D float array */
float *f1t(int n1);
/** @brief Allocation of 2D float array */
float **f2t(int n1, int n2);
/** @brief Allocation of 3D float array */
float ***f3t(int n1, int n2, int n3);
/** @brief Allocation of 4D float array */
float ****f4t(int n1, int n2, int n3, int n4);
/** @brief Allocation of 5D float array */
float *****f5t(int n1, int n2, int n3, int n4, int n5);
/** @brief Deallocation of #f5t() */
void free_f5t(float *****p);
/** @brief Deallocation of #f4t() */
void free_f4t(float ****p);
/** @brief Deallocation of #f3t() */
void free_f3t(float ***p);
/** @brief Deallocation of #f2t() */
void free_f2t(float **p);
/** @brief Deallocation of #f1t() */
void free_f1t(float *p);

/** @brief Allocation of 1D double array */
double *d1t(int n1);
/** @brief Allocation of 2D double array */
double **d2t(int n1, int n2);
/** @brief Allocation of 3D double array */
double ***d3t(int n1, int n2, int n3);
/** @brief Allocation of 4D double array */
double ****d4t(int n1, int n2, int n3, int n4);
/** @brief Deallocation of #d4t() */
void free_d4t(double ****p);
/** @brief Deallocation of #d3t() */
void free_d3t(double ***p);
/** @brief Deallocation of #d2t() */
void free_d2t(double **p);
/** @brief Deallocation of #d1t() */
void free_d1t(double *p);

/**
 * @brief Function skip tu end of line in FILE
 */
void readeol(FILE *fp);

/**
 * @brief Function calculate \f$a^n\f$ for int base
 */
int myipow(int a, int n);
/**
 * @brief Function calculate \f$a^n\f$ for float base
 */
float myfpow(float a, int n);
/**
 * @brief Function calculate \f$a^n\f$ for double base
 */
double mydpow(double a, int n);



/**
 * @brief Print out double array with dimension (n, m)
 */
void pdarray(int n, int m, double **a);
/**
 * @brief Print out double array with dimension (n)
 */
void pdvector(int n, double *a);
/**
 * @brief Print out float array with dimension (n, m)
 */
void pfarray(int n, int m, float **a);
/**
 * @brief Print out float array with dimension (n)
 */
void pfvector(int n, float *a);
/**
 * @brief Print out int array with dimension (n, m)
 */
void piarray(int n, int m, int **a);
/**
 * @brief Print out int array with dimension (n)
 */
void pivector(int n, int *a);


/**
 * @brief Set all values to zero in double array with dimension(n, m)
 */
void zdarray(int n, int m, double **a);
/**
 * @brief Set all values to zero in double array with dimension(n)
 */
void zdvector(int n, double *a);
/**
 * @brief Set all values to zero in float array with dimension(n, m)
 */
void zfarray(int n, int m, float **a);
/**
 * @brief Set all values to zero in float array with dimension(n)
 */
void zfvector(int n, float *a);
/**
 * @brief Set all values to zero in int array with dimension(n, m)
 */
void ziarray(int n, int m, int **a);
/**
 * @brief Set all values to zero in int array with dimension(n)
 */
void zivector(int n, int *a);

// --ADDED LT ---
/**
 * @brief Function read 2D array of doubles from file
 *
 * If function fail to read all numbers end up with #failed().
 *
 * Function assume that numbers are each on new line.
 * a[0][0] = 1. line
 * a[0][1] = 2. line
 *  ...    = ...
 *  ...    = ...
 * a[1][0] = #m-th line
 *  ...    = ...
 *  ...    = ...
 * a[n][m] = (n*m)-th line
 */
void rdarray			(int n, int m, double **a, FILE *F);

/**
 * @brief Function read 2D diagonaly symmetric array of doubles from file
 *
 * If function fail to read all numbers end up with #failed().
 *
 * Function assume that numbers are each on new line.
 * a[0][0] = 1. line
 * a[0][1] = 2. line
 * a[0][2] = 3. line
 *  ...    = ...
 *  ...    = ...
 * a[0][n] = n-th line
 * a[1][1] = (n+1)-th line
 *  ...    = ...
 *  ...    = ...
 * a[1][n] = (n+n-1)-th line
 *  ...    = ...
 *  ...    = ...
 * well it reads just one half of matrix + diagonal
 */
void rdarray_symm (int n, double **a, FILE *F);

/**
 * @brief Function read 1D array of doubles from file
 *
 * If function fail to read all numbers end up with #failed().
 *
 * Function assume that numbers are each on new line.
 * a[0] = 1. line
 * a[1] = 2. line
 * a[2] = 3. line
 *  ... = ...
 *  ... = ...
 */
void rdvector 		(int n, double  *a, FILE *F);
#endif //HDR_MY_MEM
