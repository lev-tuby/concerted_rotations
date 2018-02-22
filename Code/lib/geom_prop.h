#ifndef _GEOM_PROPS_
#define _GEOM_PROPS_

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
//symm tensor
typedef struct {
	int size;
	gsl_matrix *matrix;
	gsl_vector *evals;
	gsl_matrix *evecs;
} SYMM_TENS;

void inertia_tensor 			( double **coord, double *mass,int N, gsl_matrix* inrt );
SYMM_TENS * inertia_axes 	(	double **coord, double *mass, int N);
double asphericity				(	gsl_vector *semiaxis );
double prolateness				( gsl_vector *semiaxis );
//
double gyration_radius(double **coord, int N);
//
void rotate ( double **coord, int N, double center[3],double axis[3], double theta);
//
gsl_matrix * graham_schmidt ( gsl_vector *v);

#endif
