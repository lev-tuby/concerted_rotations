#ifndef __HDR_GEOM__
#define __HDR_GEOM__
/**
 * @file
 * @brief Header file for all functions for work vectors and geometric object
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/** @brief Return values of function that should return TRUE/FALSE.*/
#ifndef TRUE
#define TRUE 1
#endif
/** @brief Return values of function that should return TRUE/FALSE.*/
#ifndef FALSE
#define FALSE 0
#endif

/**
 * @brief Function calculate cross product of 2 vectors in 3D
 */
void vecprod_d (const double *a, const double *b, double *c);

/**
 * @brief Function calculate scalar product of 2 vectors in dim D
 *
 * @note Function assume both vectors are of dimensin dim. Otherwise it might cause SEGFAULT.
 */
double scal_d (const double *a, const double *b, int dim);

/**
 * @brief Function calculate norm of dim D dimensional vector
 *
 * @note Function assume both vectors are of dimensin dim. Otherwise it might cause SEGFAULT.
 */
double norm_d (const double *a, int dim);

/**
 * @brief Function normalize dim D dimensional vector
 *
 * @note Function assume both vectors are of dimensin dim. Otherwise it might cause SEGFAULT.
 */
void normalize_d (double *a, int dim);

/*******************************/
/**
 * @brief Function calculate Euclidean distance between 2 points in dim D space
 *
 * @note Function assume both vectors are of dimensin dim. Otherwise it might cause SEGFAULT.
 */
double dist_d (const double *a, const double *b, int dim);

/**
 * @brief Function calculate minimal distance between two line segments in 3D
 *
 * First segment is defined as vector p1 - p0 and second as p3 - p2.
 * @note Vectors have to be in 3D.
 */
double dist_segments(const double p0[3], const double p1[3], const double p2[3], const double p3[3]);

/**
 * @brief Function calculate intersection of line segment and vertex
 *
 * Determine whether or not the line segment p1,p2
 * Intersects the 3 vertex facet bounded by pa,pb,pc
 * Return true/false and the intersection point p
 * The equation of the line is p = p1 + mu (p2 - p1)
 * The equation of the plane is a x + b y + c z + d = 0
 * n.x x + n.y y + n.z z + d = 0
 */
int LineFacet(const double *p1, const double *p2, const double *pa, const double *pb, const double *pc, double *p);

/**
 * @brief Function calculate intersection of line segment and vertex
 *
 * Determine whether or not the line segment p1,p2
 * Intersects the 3 vertex facet bounded by pa,pb,pc
 * Return true/false and the intersection point p
 * The equation of the line is p = p1 + mu (p2 - p1)
 * The equation of the plane is a x + b y + c z + d = 0
 * n.x x + n.y y + n.z z + d = 0
 */
int LineFacetOLD( double *p1, double *p2, double *pa, double *pb, double *pc, double *p);

/**
 * @brief Function calculate center of mass of given number of 3D points
 */
void center_of_mass( int N, double ** coord, double *cm );

/**
 * @brief Function find square of farthest point from CM of points
 *
 * The following is to find the farthest vertex from the center of mass (cm),
 * by first computing the cm coordinates, then computing the distances of
 * ALL the vertices to the cm, and finally finding the maximum amongst the
 * distances. To avoid sqrt, which is repeated Len times in each run, it is 
 * better to compute the squared distances then find the maximum amongst them...
 *
 * @note However, remember to return the sqrt of the found squared radius!
 */
double find_radius(int N,double ** coord,double center[3]);

/**
 * @brief Function calculate angle \f$\angle ABC\f$ between points *A, *B, and *C
 */
double angle_ABC(double *A,double *B, double *C);

#endif //__HDR_GEOM__
