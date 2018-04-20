#ifndef __MY_QUATERNIONS_HDR
#define __MY_QUATERNIONS_HDR
/**
 * @file
 * @brief Header file for all functions for work with histogram structure
 * Well not much else to say ...
 */

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_blas.h>

/**
 * @brief Function calculate rotation matrix in 3D from quaternion
 */
int rotation3D_build 		( gsl_matrix * rot, const gsl_vector *quat);

/**
 * @brief Function calculate rotation-translation matrix in 3D from quaternion and translation vector
 *
 * builds  a matrix
 *   /  Rot tr \
 *   \  0   1  /
 */
int rototransl3D_build 	( gsl_matrix * rt, const gsl_vector *quat, const gsl_vector * transl);

/**
 * @brief Function calculate quaternion for rotation in 3D
 *
 * @param[in,out]   *q                Calculated quaternion.
 * @param[in]       *v                Vector around which rotation is performed.
 * @param[in]        theta            Angle by which rotation being performed.
 *
 * @return GSL err code
 */
int quaternion_build 		(gsl_vector * q, const gsl_vector *v, const double theta );

/**
 * @brief Function calculate multiplication of rotation-translation matrix
 *
 * Multiply two rototranslations rt1 and rt2 giving rt2(rt1(.))
 */
int rototrans3D_mult 		( gsl_matrix *rt21, const gsl_matrix *rt2, const gsl_matrix *rt1);

/**
 * @brief Function calculate multiplication of two quaternions
 *
 * Returns the multiplied quaternion q21=q2*q1.
 */
int quaternion_mult 		( gsl_vector *q21, const gsl_vector *q2, const gsl_vector *q1);

/**
 * @brief Function rotate vector via rotation matrix
 */
int rotate3D_vector_m			( gsl_vector *rV, const gsl_matrix *r, const gsl_vector *v);

/**
 * @brief Function rotate and translate vector via rotation-translation matrix
 */
int rototranslate3D_vector_m (gsl_vector *rtV, const gsl_matrix *rt, const gsl_vector *V);
#endif //__MY_QUATERNIONS_HDR
