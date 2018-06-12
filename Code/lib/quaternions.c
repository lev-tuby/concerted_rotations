/**
 * @file
 * @brief Source file contain function for manipulation with quaternions and rotations
 * Added 20180201 -- rototranslations
 */

#include "quaternions.h"
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <math.h>

/**
 * @brief Function calculate rotation matrix in 3D from quaternion
 *
 * @param[in,out]   *rot              Minimal value of histogram in 1D histogram.
 * @param[in]       *quat             Quaternion vector.
 *
 * @return GSL err code
 */
int rotation3D_build 		( gsl_matrix * rot, const gsl_vector *quat)
{
	double r[3][3];
	double q0,q1,q2,q3;
	if(quat->size!=4 || rot->size1 !=3 || rot->size2 !=3  ) {
		return GSL_EBADLEN;
	}
	q0=gsl_vector_get(quat,0);
	q1=gsl_vector_get(quat,1);
	q2=gsl_vector_get(quat,2);
	q3=gsl_vector_get(quat,3);
	//rotation.
	r[0][0]=1-2*(q2*q2+q3*q3);
	r[1][1]=1-2*(q3*q3+q1*q1);
	r[2][2]=1-2*(q1*q1+q2*q2);

	r[0][1]=2*(q1*q2-q0*q3);
	r[0][2]=2*(q1*q3+q0*q2);

	r[1][0]=2*(q1*q2+q0*q3);
	r[1][2]=2*(q2*q3-q0*q1);

	r[2][0]=2*(q1*q3-q0*q2);
	r[2][1]=2*(q2*q3+q0*q1);
	//set the rotation matrix
	for(int i=0;i<3;i++) {
		gsl_matrix_set(rot,i,0,r[i][0]);
		gsl_matrix_set(rot,i,1,r[i][1]);
		gsl_matrix_set(rot,i,2,r[i][2]);
	}
	return GSL_SUCCESS;
}

/**
 * @brief Function calculate rotation-translation matrix in 3D from quaternion and translation vector
 *
 * builds  a matrix
 *   /  Rot tr \
 *   \  0   1  /
 *
 * @param[in,out]   *rt               Minimal value of histogram in 1D histogram.
 * @param[in]       *quat             Quaternion vector.
 * @param[in]       *transl           Translation vector.
 *
 * @return GSL err code
 */
int rototransl3D_build ( gsl_matrix * rt, const gsl_vector *quat, const gsl_vector * transl)
{
	double t[4];
	if(quat->size!=4 || rt->size1 !=4 || rt->size2 !=4 ) {
		return GSL_EBADLEN;
	}
	else if (transl->size != 3 ) {
		return GSL_EBADLEN;
	}

	gsl_matrix_set_zero(rt);

	t[0]=gsl_vector_get(transl,0);
	t[1]=gsl_vector_get(transl,1);
	t[2]=gsl_vector_get(transl,2);
	t[3]=1;

	gsl_vector_const_view t_v=gsl_vector_const_view_array(t,4);
	gsl_matrix_view rot_v = gsl_matrix_submatrix ( rt,0,0,3,3);

	rotation3D_build 		( &rot_v.matrix,quat);
	gsl_matrix_set_col	(rt,3,&t_v.vector);
	return GSL_SUCCESS;
}

/**
 * @brief Function calculate quaternion for rotation in 3D
 *
 * @param[in,out]   *q                Calculated quaternion.
 * @param[in]       *v                Vector around which rotation is performed.
 * @param[in]        theta            Angle (in radians) by which rotation being performed.
 *
 * @return GSL err code
 */
int quaternion_build (gsl_vector * q, const gsl_vector *v, const double theta )
{
	// q = (cos(theta/2), \vec{v}*sin(theta/2)
	if(q->size!=4 ) {
		return GSL_EBADLEN;
	}
	if (v->size != 3 ) {
		return GSL_EBADLEN;
	}
	double norm_v=gsl_blas_dnrm2(v);
	double th_2=0.5*theta;
	gsl_vector_set(q,0,cos(th_2));
	gsl_vector_set(q,1,sin(th_2)*gsl_vector_get(v,0)/norm_v);
	gsl_vector_set(q,2,sin(th_2)*gsl_vector_get(v,1)/norm_v);
	gsl_vector_set(q,3,sin(th_2)*gsl_vector_get(v,2)/norm_v);
	return GSL_SUCCESS;
}

/**
 * @brief Function calculate multiplication of rotation-translation matrix
 *
 * Multiply two rototranslations rt1 and rt2 giving rt2(rt1(.))
 *
 * @param[in,out]   *rt21             Result of matrix multiplication.
 * @param[in]       *rt2              First rotationaly-translational matrix.
 * @param[in]       *rt1              Second rotationaly-translational matrix.
 *
 * @return GSL err code
 */
int rototrans3D_mult ( gsl_matrix *rt21, const gsl_matrix *rt2, const gsl_matrix *rt1)
{
	if(rt1->size1 !=4 || rt1->size2 !=4 ) {
		return GSL_EBADLEN;
	}
	if(rt2->size1 !=4 || rt2->size2 !=4 ) {
		return GSL_EBADLEN;
	}
	if(rt21->size1 !=4 || rt21->size2 !=4 ) {
		return GSL_EBADLEN;
	}
	int error;
	error=gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.,rt1,rt2,0.,rt21);
	return error;
}

/**
 * @brief Function calculate multiplication of two quaternions
 *
 * Returns the multiplied quaternion q21=q2*q1.
 *
 * @param[in,out]   *q21              Result of matrix multiplication.
 * @param[in]       *q2               First quaternion in 3D.
 * @param[in]       *q1               Second quaternion in 3D.
 *
 * @return GSL err code
 */
int quaternion_mult ( gsl_vector *q21, const gsl_vector *q2, const gsl_vector *q1)
{
	if(q1->size!=4 || q2->size !=4 || q21->size !=4) {
		return GSL_EBADLEN;
	}
	int error;
	double q1cmp[8];
	double q2cmp[8];
	double q21cmp[8];
	gsl_complex c_1=gsl_complex_rect(1,0);
	gsl_complex c_0=gsl_complex_rect(0,0);

	q1cmp[0]= gsl_vector_get(q1,0); //re
	q1cmp[1]= gsl_vector_get(q1,1); //im
	q1cmp[2]= gsl_vector_get(q1,2); //..etc
	q1cmp[3]= gsl_vector_get(q1,3);
	q1cmp[4]= -gsl_vector_get(q1,2);
	q1cmp[5]= gsl_vector_get(q1,3);
	q1cmp[6]= gsl_vector_get(q1,0);
	q1cmp[7]=-gsl_vector_get(q1,1);

	q2cmp[0]= gsl_vector_get(q2,0);
	q2cmp[1]= gsl_vector_get(q2,1);
	q2cmp[2]= gsl_vector_get(q2,2);
	q2cmp[3]= gsl_vector_get(q2,3);
	q2cmp[4]= -gsl_vector_get(q2,2);
	q2cmp[5]= gsl_vector_get(q2,3);
	q2cmp[6]= gsl_vector_get(q2,0);
	q2cmp[7]=-gsl_vector_get(q2,1);

	gsl_matrix_complex_const_view q1_v=gsl_matrix_complex_const_view_array(q1cmp,2,2);
	gsl_matrix_complex_const_view q2_v=gsl_matrix_complex_const_view_array(q2cmp,2,2);
	gsl_matrix_complex_view q21_v=gsl_matrix_complex_view_array(q21cmp,2,2);

	error=gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,c_1,&q2_v.matrix,&q1_v.matrix,c_0,&q21_v.matrix);

	if(error!=GSL_SUCCESS) { return error; }
	gsl_vector_set(q21,0,GSL_REAL(gsl_matrix_complex_get(&q21_v.matrix,0,0)));
	gsl_vector_set(q21,1,GSL_IMAG(gsl_matrix_complex_get(&q21_v.matrix,0,0)));
	gsl_vector_set(q21,2,GSL_REAL(gsl_matrix_complex_get(&q21_v.matrix,0,1)));
	gsl_vector_set(q21,3,GSL_IMAG(gsl_matrix_complex_get(&q21_v.matrix,0,1)));
	return GSL_SUCCESS;
}
/*
int quaternion_mult ( gsl_vector *q21, const gsl_vector *q2, const gsl_vector *q1)
{
	//returns the multiplied quaternion q21=q2*q1
	if(q1->size!=4 || q2->size !=4 || q21->size !=4) {
		return GSL_EBADLEN;
	}
	int error;
	double q1cmp[8];
	double q2cmp[8];
	double q21cmp[8];
	gsl_complex c_1=gsl_complex_rect(1,0);
	gsl_complex c_0=gsl_complex_rect(0,0);

	q1cmp[0]= gsl_vector_get(q1,0); //re
	q1cmp[1]= gsl_vector_get(q1,1); //im
	q1cmp[2]= gsl_vector_get(q1,2); //..etc
	q1cmp[3]= gsl_vector_get(q1,3);
	q1cmp[4]= -gsl_vector_get(q1,2);
	q1cmp[5]= gsl_vector_get(q1,3);
	q1cmp[6]= gsl_vector_get(q1,0);
	q1cmp[7]=-gsl_vector_get(q1,1);

	q2cmp[0]= gsl_vector_get(q2,0);
	q2cmp[1]= gsl_vector_get(q2,1);
	q2cmp[2]= gsl_vector_get(q2,2);
	q2cmp[3]= gsl_vector_get(q2,3);
	q2cmp[4]=-gsl_vector_get(q2,2);
	q2cmp[5]= gsl_vector_get(q2,3);
	q2cmp[6]= gsl_vector_get(q2,0);
	q2cmp[7]=-gsl_vector_get(q2,1);

	gsl_matrix_complex_const_view q1_v=gsl_matrix_complex_const_view_array(q1cmp,2,2);
	gsl_matrix_complex_const_view q2_v=gsl_matrix_complex_const_view_array(q2cmp,2,2);
	gsl_matrix_complex_view q21_v=gsl_matrix_complex_view_array(q21cmp,2,2);

	error=gsl_blas_zgemm(CblasNoTrans,CblasNoTrans,c_1,&q2_v.matrix,&q1_v.matrix,c_0,&q21_v.matrix);

	if(error!=GSL_SUCCESS) { return error; }
	gsl_vector_set(q21,0,GSL_REAL(gsl_matrix_complex_get(&q21_v.matrix,0,0)));
	gsl_vector_set(q21,1,GSL_IMAG(gsl_matrix_complex_get(&q21_v.matrix,0,0)));
	gsl_vector_set(q21,2,GSL_REAL(gsl_matrix_complex_get(&q21_v.matrix,0,1)));
	gsl_vector_set(q21,3,GSL_IMAG(gsl_matrix_complex_get(&q21_v.matrix,0,1)));
	return GSL_SUCCESS;
}
*/

/**
 * @brief Function rotate vector via rotation matrix
 *
 * @param[in,out]   *rV               3D vector rotated by rotation matrix.
 * @param[in]       *r                Rotation matrix in 3D.
 * @param[in]       *V                3D vector to be rotated.
 *
 * @return GSL err code
 */
int rotate3D_vector_m			( gsl_vector *rV, const gsl_matrix *r, const gsl_vector *V)
{
	if (r->size1!=3 || r->size2!=3 || V->size!=3 || rV->size !=3 ) {
		return GSL_EBADLEN;
	}
	int error=gsl_blas_dgemv(CblasNoTrans,1.,r,V,0.,rV);
	return error;
}

/**
 * @brief Function rotate and translate vector via rotation-translation matrix
 *
 * @param[in,out]   *rtV              3D vector rotated and translated by rotation-translation matrix.
 * @param[in]       *rt               Rotation-translation matrix in 3D.
 * @param[in]       *V                3D vector to be rotated and translated.
 *
 * @return GSL err code
 */
int rototranslate3D_vector_m (gsl_vector *rtV, const gsl_matrix *rt, const gsl_vector *V)
{
	if (rt->size1!=4 || rt->size2!=4 || V->size!=3 || rtV->size !=3 ) {
		return GSL_EBADLEN;
	}
	double tmp_V[4]={0.,0.,0.,0.};
	double tmp_rtV[4];
	tmp_V[0]=gsl_vector_get(V,0);
	tmp_V[1]=gsl_vector_get(V,1);
	tmp_V[2]=gsl_vector_get(V,2);
	tmp_V[3]=1;
	gsl_vector_const_view tmp_V_v = gsl_vector_const_view_array(tmp_V,4);
	gsl_vector_view tmp_rtV_v = gsl_vector_view_array(tmp_rtV,4);
	int error=gsl_blas_dgemv(CblasNoTrans,1.,rt,&tmp_V_v.vector,0.,&tmp_rtV_v.vector);
	gsl_vector_set(rtV,0,tmp_rtV[0]);
	gsl_vector_set(rtV,1,tmp_rtV[1]);
	gsl_vector_set(rtV,2,tmp_rtV[2]);
	return error;
}

