#ifndef __HDR_GEOM__
#define __HDR_GEOM__
/**
 * @file
 * @brief Source file for all functions for work vectors and geometric object. Adapted from the <a href="https://github.com/luca-tubiana/KymoKnot"> KymoKnot github repository </a> by Luca Tubiana, Guido Polles, Enzo Orlandini and Cristian Micheletti.
 * Vectors are implemented as arrays without boundary checks.
 */

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

/** @brief Return values of function that should return TRUE/FALSE.*/
#ifndef TRUE
#define TRUE 1
#endif
/** @brief Return values of function that should return TRUE/FALSE.*/
#ifndef FALSE
#define FALSE 0
#endif

/**
 * @brief Computes the cross product of 2 vectors in 3D
 */
void vecprod_d (const double *a, const double *b, double *c);

/**
 * @brief Computes the scalar product of 2 vectors in dim D
 *
 */
double scal_d (const double *a, const double *b, int dim);

/**
 * @brief Computes the norm of a dim D dimensional vector
 *
 */
double norm_d (const double *a, int dim);

/**
 * @brief Normalizes a D dimensional vector
 *
 */
void normalize_d (double *a, int dim);

/*******************************/
/**
 * @brief Calculates the Euclidean distance between 2 points in dim D space
 *
 */
double dist_d (const double *a, const double *b, int dim);

/**
 * @brief Calculates the minimal distance between two line segments in 3D
 *
 */
double dist_segments(const double p0[3], const double p1[3], const double p2[3], const double p3[3]);

/**
 * @brief Calculate intersection of line segment and a triangle
 */
int LineFacet(const double *p1, const double *p2, const double *pa, const double *pb, const double *pc, double *p);


/**
 * @brief Computes the center of mass of given number of 3D points
 */
void center_of_mass( int N, double ** coord, double *cm );

/**
 * @brief Function find square of farthest point from CM of points
 *
 */
double find_radius(int N,double ** coord,double center[3]);

/**
 * @brief Computes the angle \f$\angle ABC\f$ between points *A, *B, and *C
 */
double angle_ABC(double *A,double *B, double *C);

/**
 * @brief Rotates a set of coordinates
 *
 */
void rotate ( double **coord, int N, double center[3],double axis[3], double theta);

/**
 * @brief Creates an orthonormal basis with respect to 3D vector using the Graham-Schmidt approach.
 */
gsl_matrix * gram_schmidt ( gsl_vector *v);



#endif //__HDR_GEOM__
