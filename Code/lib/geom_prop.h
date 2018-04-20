#ifndef __GEOM_PROPS__
#define __GEOM_PROPS__
/**
 * @file
 * @brief Header file for functions to calculate geometrical properties
 * @todo We might merge my_geom and this ... since they are sort of related anyway ...
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>

#include "my_geom.h"
#include "my_memory.h"

/**
 * @brief Symetry tensor data
 */
typedef struct {
	int size;           /**< Dimension size. */
	gsl_matrix *matrix; /**< Actual symetry tensor. */
	gsl_vector *evals;  /**< Eigenvalues. */
	gsl_matrix *evecs;  /**< Eigenvectors. */
} SYMM_TENS;

/**
 * @brief Function generate inertia tensor from point coordinates and their masses
 */
void inertia_tensor 			( double **coord, double *mass,int N, gsl_matrix* inrt );

/**
 * @brief Function generate #SYMM_TENS tensor structure from point coordinates and their masses
 *
 * Generation of inertia tensor inside function is done by inertia_tensor().
 * Eigenvectors and values are then calculated from inertia tensor.
 */
SYMM_TENS * inertia_axes 	(	double **coord, double *mass, int N);

/**
 * @brief Function calculate asphericity
 *
 * @todo Add more info.
 */
double asphericity				(	gsl_vector *semiaxis );

/**
 * @brief Function calculate prolateness
 *
 * @todo Add more info.
 */
double prolateness				( gsl_vector *semiaxis );

/**
 * @brief Function calculate radius of gyration
 */
double gyration_radius(double **coord, int N);

/**
 * @brief Function rotate set of N points along given center and rotation axis by given angle
 */
void rotate ( double **coord, int N, double center[3],double axis[3], double theta);

/**
 * @brief Function create orthonormal basis with respect to 3D vector
 *
 * Metod use Gram-Schmidt orthogonalization.
 */
gsl_matrix * graham_schmidt ( gsl_vector *v);

#endif //__GEOM_PROPS__
