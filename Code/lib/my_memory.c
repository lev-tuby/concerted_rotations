/**
 * @file
 * @brief Source file contain functions for memory allocation/deallocation
 */

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include "my_memory.h"
#include "messages.h"

/**
 * @brief Allocation of 1D FILE array
 *
 * @param[in]        n1              Number of files in array.
 *
 * @return pointer to array of FILEs
 */
FILE **F2t (int n1)
{
  FILE **p;
  if ((p = (FILE **) malloc ((size_t) n1 * sizeof (FILE *))) == NULL)
      failed ("F2t: failed n1");
  return p;
}

/**
 * @brief Allocation of 2D FILE array
 *
 * @param[in]        n1              Size of first dimension.
 * @param[in]        n2              Size of first dimension.
 *
 * @return pointer to matrix (n1, n2) of FILEs
 */
FILE ***F3t (int n1, int n2)
{
  FILE ***p;
  int i;
  if ((p = (FILE ***) malloc ((size_t) n1 * sizeof (FILE **))) == NULL)
      failed ("F3t: failed n1");
  if ((p[0] = (FILE **) malloc ((size_t) n1 * n2 * sizeof (FILE *))) == NULL)
      failed ("F3t: failed n2");
  for (i = 0; i < n1 - 1; i++)
    p[i + 1] = p[i] + n2;
  return p;
}

/**
 * @brief Allocation of 1D char array
 *
 * @param[in]        n1              Size of first dimension.
 *
 * @return pointer to array of char
 */
char *c1t (int n1)
{
  char *p;
  if ((p = (char *) malloc ((size_t) n1 * sizeof (char))) == NULL)
      failed ("c1t: failed");
  return p;
}

/**
 * @brief Allocation of 2D char array
 *
 * @param[in]        n1              Size of first dimension.
 * @param[in]        n2              Size of second dimension.
 *
 * @return pointer to matrix (n1, n2) of FILEs
 */
char **c2t (int n1, int n2)
{
  char **p;
  int i;
  if ((p = (char **) malloc ((size_t) n1 * sizeof (char *))) == NULL)
      failed ("c2t: failed n1");
  if ((p[0] = (char *) malloc ((size_t) n1 * n2 * sizeof (char))) == NULL)
      failed ("c2t: failed n2");
  for (i = 0; i < n1 - 1; i++)
    p[i + 1] = p[i] + n2;
  return p;
}

/**
 * @brief Allocation of 3D char array
 *
 * @param[in]        n1              Size of first dimension.
 * @param[in]        n2              Size of second dimension.
 * @param[in]        n3              Size of third dimension.
 *
 * @return pointer to char matrix
 */
char ***c3t (int n1, int n2, int n3)
{
  char ***p;
  int i, j;
  if ((p = (char ***) malloc ((size_t) n1 * sizeof (char **))) == NULL)
      failed ("c3t: failed n1");
  if ((p[0] = (char **) malloc ((size_t) n1 * n2 * sizeof (char *))) == NULL)
      failed ("c3t: failed n2");
  if ((p[0][0] = (char *) malloc ((size_t) n1 * n2 * n3 * sizeof (char))) == NULL)
      failed ("c3t: failed n3");
  for (i = 0; i < n2 - 1; i++)
    p[0][i + 1] = p[0][i] + n3;
  for (i = 0; i < n1 - 1; i++)
    {
      p[i + 1] = p[i] + n2;
      p[i + 1][0] = p[i][0] + n2 * n3;
      for (j = 0; j < n2 - 1; j++)
	p[i + 1][j + 1] = p[i + 1][j] + n3;
    }
  return p;
}

/**
 * @brief Allocation of 4D char array
 *
 * @param[in]        n1              Size of first dimension.
 * @param[in]        n2              Size of second dimension.
 * @param[in]        n3              Size of third dimension.
 * @param[in]        n4              Size of fourth dimension.
 *
 * @return pointer to char matrix
 */
char ****c4t (int n1, int n2, int n3, int n4)
{
  char ****p, *a;
  int i, j, k;
  if ((p = (char ****) malloc ((size_t) n1 * sizeof (char ***))) == NULL)
      failed ("c4t: failed n1");
  if ((p[0] = (char ***) malloc ((size_t) n1 * n2 * sizeof (char **))) == NULL)
      failed ("c4t: failed n2");
  if ((p[0][0] = (char **) malloc ((size_t) n1 * n2 * n3 * sizeof (char *))) == NULL)
      failed ("c4t: failed n3");
  if ((p[0][0][0] = (char *) malloc ((size_t) n1 * n2 * n3 * n4 * sizeof (char ))) == NULL)
      failed ("c4t: failed n4");

  for (i = 0; i < n3 - 1; i++)
    p[0][0][i + 1] = p[0][0][i] + n4;

  for (i = 0; i < n2 - 1; i++)
    {
      p[0][i + 1] = p[0][i] + n3;
      p[0][i + 1][0] = p[0][i][0] + n3 * n4;
      for (j = 0; j < n3 - 1; j++)
	p[0][i + 1][j + 1] = p[0][i + 1][j] + n4;
    }

  for (i = 0; i < n1 - 1; i++)
    {
      p[i + 1] = p[i] + n2;
      p[i + 1][0] = p[i][0] + n2 * n3;
      p[i + 1][0][0] = p[i][0][0] + n2 * n3 * n4;

      for (j = 0; j < n3 - 1; j++)
	p[i + 1][0][j + 1] = p[i + 1][0][j] + n4;

      for (j = 0; j < n2 - 1; j++)
	{
	  p[i + 1][j + 1] = p[i + 1][j] + n3;
	  p[i + 1][j + 1][0] = p[i + 1][j][0] + n3 * n4;
	  for (k = 0; k < n3 - 1; k++)
	    p[i + 1][j + 1][k + 1] = p[i + 1][j + 1][k] + n4;
	}
    }

  for (i = 0, a = p[0][0][0]; i < n1 * n2 * n3 * n4; i++)
    *a++ = 0.0;

  return p;
}

/**
 * @brief Deallocation of #c4t() 
 *
 * @param[in,out]  ****p              Array.
 *
 * @return \c void
 */
void free_c4t(char ****p){

  free(p[0][0][0]);
  free(p[0][0]);
  free(p[0]);
  free(p);
}

/**
 * @brief Deallocation of #c3t() 
 *
 * @param[in,out]   ***p              Array.
 *
 * @return \c void
 */
void free_c3t(char ***p){

  free(p[0][0]);
  free(p[0]);
  free(p);
}

/**
 * @brief Deallocation of #c2t() 
 *
 * @param[in,out]    **p              Array.
 *
 * @return \c void
 */
void free_c2t(char **p){

  free(p[0]);
  free(p);
}

/**
 * @brief Deallocation of #c1t() 
 *
 * @param[in,out]     *p              Array.
 *
 * @return \c void
 */
void free_c1t(char *p){

  free(p);
}


/*********************************/

/**
 * @brief Allocation of 1D short array
 *
 * @param[in]        n1              Size of first dimension.
 *
 * @return pointer to short array
 */
short *s1t (int n1)
{
  short *p, *a;
  int i;
  if ((p = (short *) malloc ((size_t) n1 * sizeof (short))) == NULL)
      failed ("s1t: failed");
  for (i = 0, a = p; i < n1; i++)
    *a++ = 0;
  return p;
}

/**
 * @brief Allocation of 2D short array
 *
 * @param[in]        n1              Size of first dimension.
 * @param[in]        n2              Size of second dimension.
 *
 * @return pointer to short array
 */
short **s2t (int n1, int n2)
{
  short **p, *a;
  int i;
  if ((p = (short **) malloc ((size_t) n1 * sizeof (short *))) == NULL)
      failed ("s2t: failed n1");
  if ((p[0] = (short *) malloc ((size_t) n1 * n2 * sizeof (short))) == NULL)
      failed ("s2t: failed n2");
  for (i = 0; i < n1 - 1; i++)
    p[i + 1] = p[i] + n2;
  for (i = 0, a = p[0]; i < n1 * n2; i++)
    *a++ = 0;
  return p;
}

/**
 * @brief Allocation of 3D short array
 *
 * @param[in]        n1              Size of first dimension.
 * @param[in]        n2              Size of second dimension.
 * @param[in]        n3              Size of third dimension.
 *
 * @return pointer to short array
 */
short ***s3t (int n1, int n2, int n3)
{
  short ***p, *a;
  int i, j;
  if ((p = (short ***) malloc ((size_t) n1 * sizeof (short **))) == NULL)
      failed ("s3t: failed n1");
  if ((p[0] = (short **) malloc ((size_t) n1 * n2 * sizeof (short *))) == NULL)
      failed ("s3t: failed n2");
  if ((p[0][0] = (short *) malloc ((size_t) n1 * n2 * n3 * sizeof (short))) == NULL)
      failed ("s3t: failed n3");
  for (i = 0; i < n2 - 1; i++)
    p[0][i + 1] = p[0][i] + n3;
  for (i = 0; i < n1 - 1; i++)
    {
      p[i + 1] = p[i] + n2;
      p[i + 1][0] = p[i][0] + n2 * n3;
      for (j = 0; j < n2 - 1; j++)
	p[i + 1][j + 1] = p[i + 1][j] + n3;
    }
  for (i = 0, a = p[0][0]; i < n1 * n2 * n3; i++)
    *a++ = 0;
  return p;
}

/**
 * @brief Allocation of 4D short array
 *
 * @param[in]        n1              Size of first dimension.
 * @param[in]        n2              Size of second dimension.
 * @param[in]        n3              Size of third dimension.
 * @param[in]        n4              Size of fourth dimension.
 *
 * @return pointer to short array
 */
short ****s4t (int n1, int n2, int n3, int n4)
{
  short ****p, *a;
  int i, j, k;
  if ((p = (short ****) malloc ((size_t) n1 * sizeof (short ***))) == NULL)
      failed ("s4t: failed n1");
  if ((p[0] = (short ***) malloc ((size_t) n1 * n2 * sizeof (short **))) == NULL)
      failed ("s4t: failed n2");
  if ((p[0][0] = (short **) malloc ((size_t) n1 * n2 * n3 * sizeof (short *))) == NULL)
      failed ("s4t: failed n3");
  if ((p[0][0][0] = (short *) malloc ((size_t) n1 * n2 * n3 * n4 * sizeof (short ))) == NULL)
      failed ("s4t: failed n4");

  for (i = 0; i < n3 - 1; i++)
    p[0][0][i + 1] = p[0][0][i] + n4;

  for (i = 0; i < n2 - 1; i++)
    {
      p[0][i + 1] = p[0][i] + n3;
      p[0][i + 1][0] = p[0][i][0] + n3 * n4;
      for (j = 0; j < n3 - 1; j++)
	p[0][i + 1][j + 1] = p[0][i + 1][j] + n4;
    }

  for (i = 0; i < n1 - 1; i++)
    {
      p[i + 1] = p[i] + n2;
      p[i + 1][0] = p[i][0] + n2 * n3;
      p[i + 1][0][0] = p[i][0][0] + n2 * n3 * n4;

      for (j = 0; j < n3 - 1; j++)
	p[i + 1][0][j + 1] = p[i + 1][0][j] + n4;

      for (j = 0; j < n2 - 1; j++)
	{
	  p[i + 1][j + 1] = p[i + 1][j] + n3;
	  p[i + 1][j + 1][0] = p[i + 1][j][0] + n3 * n4;
	  for (k = 0; k < n3 - 1; k++)
	    p[i + 1][j + 1][k + 1] = p[i + 1][j + 1][k] + n4;
	}
    }

  for (i = 0, a = p[0][0][0]; i < n1 * n2 * n3 * n4; i++)
    *a++ = 0.0;

  return p;
}

/**
 * @brief Deallocation of #s4t() 
 *
 * @param[in,out]  ****p              Array.
 *
 * @return \c void
 */
void free_s4t(short ****p){

  free(p[0][0][0]);
  free(p[0][0]);
  free(p[0]);
  free(p);
}

/**
 * @brief Deallocation of #s3t() 
 *
 * @param[in,out]   ***p              Array.
 *
 * @return \c void
 */
void free_s3t(short ***p){

  free(p[0][0]);
  free(p[0]);
  free(p);
}

/**
 * @brief Deallocation of #s2t() 
 *
 * @param[in,out]    **p              Array.
 *
 * @return \c void
 */
void free_s2t(short **p){

  free(p[0]);
  free(p);
}

/**
 * @brief Deallocation of #s1t() 
 *
 * @param[in,out]     *p              Array.
 *
 * @return \c void
 */
void free_s1t(short *p){

  free(p);
}
/*********************************/


/**
 * @brief Allocation of 1D int array
 *
 * @param[in]        n1              Size of first dimension.
 *
 * @return pointer to int array
 */
int *i1t (int n1)
{
  int *p, *a;
  int i;
  if ((p = (int *) malloc ((size_t) n1 * sizeof (int))) == NULL)
      failed ("i1t: failed");
  for (i = 0, a = p; i < n1; i++)
    *a++ = 0;
  return p;
}

/**
 * @brief Allocation of 2D int array
 *
 * @param[in]        n1              Size of first dimension.
 * @param[in]        n2              Size of second dimension.
 *
 * @return pointer to int array
 */
int **i2t (int n1, int n2)
{
  int **p, *a;
  int i;
  if ((p = (int **) malloc ((size_t) n1 * sizeof (int *))) == NULL)
      failed ("i2t: failed n1");     
  if ((p[0] = (int *) malloc ((size_t) n1 * n2 * sizeof (int))) == NULL)
      failed ("i2t: failed n2");
  for (i = 0; i < n1 - 1; i++)
    p[i + 1] = p[i] + n2;
  for (i = 0, a = p[0]; i < n1 * n2; i++)
    *a++ = 0;
  return p;
}

/**
 * @brief Allocation of 3D int array
 *
 * @param[in]        n1              Size of first dimension.
 * @param[in]        n2              Size of second dimension.
 * @param[in]        n3              Size of third dimension.
 *
 * @return pointer to int array
 */
int ***i3t (int n1, int n2, int n3)
{
  int ***p, *a;
  int i, j;
  if ((p = (int ***) malloc ((size_t) n1 * sizeof (int **))) == NULL)
      failed ("i3t: failed n1");
  if ((p[0] = (int **) malloc ((size_t) n1 * n2 * sizeof (int *))) == NULL)
      failed ("i3t: failed n2");
  if ((p[0][0] = (int *) malloc ((size_t) n1 * n2 * n3 * sizeof (int))) == NULL)
      failed ("i3t: failed n3");
  for (i = 0; i < n2 - 1; i++)
    p[0][i + 1] = p[0][i] + n3;
  for (i = 0; i < n1 - 1; i++)
    {
      p[i + 1] = p[i] + n2;
      p[i + 1][0] = p[i][0] + n2 * n3;
      for (j = 0; j < n2 - 1; j++)
	p[i + 1][j + 1] = p[i + 1][j] + n3;
    }
  for (i = 0, a = p[0][0]; i < n1 * n2 * n3; i++)
    *a++ = 0;
  return p;
}

/**
 * @brief Allocation of 4D int array
 *
 * @param[in]        n1              Size of first dimension.
 * @param[in]        n2              Size of second dimension.
 * @param[in]        n3              Size of third dimension.
 * @param[in]        n4              Size of fourth dimension.
 *
 * @return pointer to int array
 */
int ****i4t (int n1, int n2, int n3, int n4)
{
  int ****p, *a;
  int i, j, k;
  if ((p = (int ****) malloc ((size_t) n1 * sizeof (int ***))) == NULL)
      failed ("i4t: failed n1");
  if ((p[0] = (int ***) malloc ((size_t) n1 * n2 * sizeof (int **))) == NULL)
      failed ("i4t: failed n2");
  if ((p[0][0] = (int **) malloc ((size_t) n1 * n2 * n3 * sizeof (int *))) == NULL)
      failed ("i4t: failed n3");
  if ((p[0][0][0] = (int *) malloc ((size_t) n1 * n2 * n3 * n4 * sizeof (int ))) == NULL)
      failed ("i4t: failed n4");

  for (i = 0; i < n3 - 1; i++)
    p[0][0][i + 1] = p[0][0][i] + n4;

  for (i = 0; i < n2 - 1; i++)
    {
      p[0][i + 1] = p[0][i] + n3;
      p[0][i + 1][0] = p[0][i][0] + n3 * n4;
      for (j = 0; j < n3 - 1; j++)
	p[0][i + 1][j + 1] = p[0][i + 1][j] + n4;
    }

  for (i = 0; i < n1 - 1; i++)
    {
      p[i + 1] = p[i] + n2;
      p[i + 1][0] = p[i][0] + n2 * n3;
      p[i + 1][0][0] = p[i][0][0] + n2 * n3 * n4;

      for (j = 0; j < n3 - 1; j++)
	p[i + 1][0][j + 1] = p[i + 1][0][j] + n4;

      for (j = 0; j < n2 - 1; j++)
	{
	  p[i + 1][j + 1] = p[i + 1][j] + n3;
	  p[i + 1][j + 1][0] = p[i + 1][j][0] + n3 * n4;
	  for (k = 0; k < n3 - 1; k++)
	    p[i + 1][j + 1][k + 1] = p[i + 1][j + 1][k] + n4;
	}
    }

  for (i = 0, a = p[0][0][0]; i < n1 * n2 * n3 * n4; i++)
    *a++ = 0.0;

  return p;
}



/**
 * @brief Deallocation of #i4t() 
 *
 * @param[in,out]  ****p              Array.
 *
 * @return \c void
 */
void free_i4t(int ****p){

  free(p[0][0][0]);
  free(p[0][0]);
  free(p[0]);
  free(p);
}

/**
 * @brief Deallocation of #i4t() 
 *
 * @param[in,out]   ***p              Array.
 *
 * @return \c void
 */
void free_i3t(int ***p){

  free(p[0][0]);
  free(p[0]);
  free(p);
}

/**
 * @brief Deallocation of #i4t() 
 *
 * @param[in,out]    **p              Array.
 *
 * @return \c void
 */
void free_i2t(int **p){

  free(p[0]);
  free(p);
}

/**
 * @brief Deallocation of #i4t() 
 *
 * @param[in,out]     *p              Array.
 *
 * @return \c void
 */
void free_i1t(int *p){

  free(p);
}

/*********************************/

/**
 * @brief Allocation of 1D float array
 *
 * @param[in]        n1              Size of first dimension.
 *
 * @return pointer to float array
 */
float *f1t (int n1)
{
  float *p, *a;
  int i;
  if ((p = (float *) malloc ((size_t) n1 * sizeof (float))) == NULL)
      failed ("f1t: failed");
  for (i = 0, a = p; i < n1; i++)
    *a++ = 0;
  return p;
}

/**
 * @brief Allocation of 2D float array
 *
 * @param[in]        n1              Size of first dimension.
 * @param[in]        n2              Size of second dimension.
 *
 * @return pointer to float array
 */
float **f2t (int n1, int n2)
{
  float **p, *a;
  int i;
  if ((p = (float **) malloc ((size_t) n1 * sizeof (float *))) == NULL)
      failed ("f2t: failed n1");
  if ((p[0] = (float *) malloc ((size_t) n1 * n2 * sizeof (float))) == NULL)
      failed ("f2t: failed n2");
  for (i = 0; i < n1 - 1; i++)
    p[i + 1] = p[i] + n2;
  for (i = 0, a = p[0]; i < n1 * n2; i++)
    *a++ = 0;
  return p;
}

/**
 * @brief Allocation of 3D float array
 *
 * @param[in]        n1              Size of first dimension.
 * @param[in]        n2              Size of second dimension.
 * @param[in]        n3              Size of third dimension.
 *
 * @return pointer to float array
 */
float ***f3t (int n1, int n2, int n3)
{
  float ***p, *a;
  int i, j;
  if ((p = (float ***) malloc ((size_t) n1 * sizeof (float **))) == NULL)
      failed ("f3t: failed n1");
  if ((p[0] = (float **) malloc ((size_t) n1 * n2 * sizeof (float *))) == NULL)
      failed ("f3t: failed n2");
  if ((p[0][0] = (float *) malloc ((size_t) n1 * n2 * n3 * sizeof (float))) == NULL)
      failed ("f3t: failed n3");
  for (i = 0; i < n2 - 1; i++)
    p[0][i + 1] = p[0][i] + n3;
  for (i = 0; i < n1 - 1; i++)
    {
      p[i + 1] = p[i] + n2;
      p[i + 1][0] = p[i][0] + n2 * n3;
      for (j = 0; j < n2 - 1; j++)
	p[i + 1][j + 1] = p[i + 1][j] + n3;
    }
  for (i = 0, a = p[0][0]; i < n1 * n2 * n3; i++)
    *a++ = 0;
  return p;
}

/**
 * @brief Allocation of 4D float array
 *
 * @param[in]        n1              Size of first dimension.
 * @param[in]        n2              Size of second dimension.
 * @param[in]        n3              Size of third dimension.
 * @param[in]        n4              Size of fourth dimension.
 *
 * @return pointer to float array
 */
float ****f4t (int n1, int n2, int n3, int n4)
{
  float ****p, *a;
  int i, j, k;
  if ((p = (float ****) malloc ((size_t) n1 * sizeof (float ***))) == NULL)
      failed ("f4t: failed n1");
  if ((p[0] = (float ***) malloc ((size_t) n1 * n2 * sizeof (float **))) == NULL)
      failed ("f4t: failed n2");
  if ((p[0][0] = (float **) malloc ((size_t) n1 * n2 * n3 * sizeof (float *))) == NULL)
      failed ("f4t: failed n3");
  if ((p[0][0][0] = (float *) malloc ((size_t) n1 * n2 * n3 * n4 * sizeof (float))) == NULL)
      failed ("f4t: failed n4");

  for (i = 0; i < n3 - 1; i++)
    p[0][0][i + 1] = p[0][0][i] + n4;

  for (i = 0; i < n2 - 1; i++)
    {
      p[0][i + 1] = p[0][i] + n3;
      p[0][i + 1][0] = p[0][i][0] + n3 * n4;
      for (j = 0; j < n3 - 1; j++)
	p[0][i + 1][j + 1] = p[0][i + 1][j] + n4;
    }

  for (i = 0; i < n1 - 1; i++)
    {
      p[i + 1] = p[i] + n2;
      p[i + 1][0] = p[i][0] + n2 * n3;
      p[i + 1][0][0] = p[i][0][0] + n2 * n3 * n4;

      for (j = 0; j < n3 - 1; j++)
	p[i + 1][0][j + 1] = p[i + 1][0][j] + n4;

      for (j = 0; j < n2 - 1; j++)
	{
	  p[i + 1][j + 1] = p[i + 1][j] + n3;
	  p[i + 1][j + 1][0] = p[i + 1][j][0] + n3 * n4;
	  for (k = 0; k < n3 - 1; k++)
	    p[i + 1][j + 1][k + 1] = p[i + 1][j + 1][k] + n4;
	}
    }

  for (i = 0, a = p[0][0][0]; i < n1 * n2 * n3 * n4; i++)
    *a++ = 0.0;

  return p;
}

/**
 * @brief Allocation of 5D float array
 *
 * @param[in]        n1              Size of first dimension.
 * @param[in]        n2              Size of second dimension.
 * @param[in]        n3              Size of third dimension.
 * @param[in]        n4              Size of fourth dimension.
 * @param[in]        n5              Size of fifth dimension.
 *
 * @return pointer to float array
 */
float *****f5t (int n1, int n2, int n3, int n4, int n5)
{
  float *****p, *a;
  int i, j, k, l;
  if ((p = (float *****) malloc ((size_t) n1 * sizeof (float ****))) == NULL)
      failed ("f5t: failed n1");
  if ((p[0] = (float ****) malloc ((size_t) n1 * n2 * sizeof (float ***))) == NULL)
      failed ("f5t: failed n2");
  if ((p[0][0] =
    (float ***) malloc ((size_t) n1 * n2 * n3 * sizeof (float **))) == NULL)
      failed ("f5t: failed n3");
  if ((p[0][0][0] =
  (float **) malloc ((size_t) n1 * n2 * n3 * n4 * sizeof (float *))) == NULL)
      failed ("f5t: failed n4");
  if ((p[0][0][0][0] =
       (float *) malloc ((size_t) n1 * n2 * n3 * n4 * n5 * sizeof (float))) == NULL)
      failed ("f5t: failed n5");

  for (i = 0; i < n4 - 1; i++)
    p[0][0][0][i + 1] = p[0][0][0][i] + n5;

  for (i = 0; i < n3 - 1; i++)
    {
      p[0][0][i + 1] = p[0][0][i] + n4;
      p[0][0][i + 1][0] = p[0][0][i][0] + n4 * n5;
      for (j = 0; j < n4 - 1; j++)
	p[0][0][i + 1][j + 1] = p[0][0][i + 1][j] + n5;
    }

  for (i = 0; i < n2 - 1; i++)
    {
      p[0][i + 1] = p[0][i] + n3;
      p[0][i + 1][0] = p[0][i][0] + n3 * n4;
      p[0][i + 1][0][0] = p[0][i][0][0] + n3 * n4 * n5;

      for (j = 0; j < n4 - 1; j++)
	p[0][i + 1][0][j + 1] = p[0][i + 1][0][j] + n5;

      for (j = 0; j < n3 - 1; j++)
	{
	  p[0][i + 1][j + 1] = p[0][i + 1][j] + n4;
	  p[0][i + 1][j + 1][0] = p[0][i + 1][j][0] + n4 * n5;
	  for (k = 0; k < n4 - 1; k++)
	    p[0][i + 1][j + 1][k + 1] = p[0][i + 1][j + 1][k] + n5;
	}
    }

  for (i = 0; i < n1 - 1; i++)
    {
      p[i + 1] = p[i] + n2;
      p[i + 1][0] = p[i][0] + n2 * n3;
      p[i + 1][0][0] = p[i][0][0] + n2 * n3 * n4;
      p[i + 1][0][0][0] = p[i][0][0][0] + n2 * n3 * n4 * n5;

      for (j = 0; j < n4 - 1; j++)
	p[i + 1][0][0][j + 1] = p[i + 1][0][0][j] + n5;

      for (j = 0; j < n3 - 1; j++)
	{
	  p[i + 1][0][j + 1] = p[i + 1][0][j] + n4;
	  p[i + 1][0][j + 1][0] = p[i + 1][0][j][0] + n4 * n5;
	  for (k = 0; k < n4 - 1; k++)
	    p[i + 1][0][j + 1][k + 1] = p[i + 1][0][j + 1][k] + n5;
	}

      for (j = 0; j < n2 - 1; j++)
	{
	  p[i + 1][j + 1] = p[i + 1][j] + n3;
	  p[i + 1][j + 1][0] = p[i + 1][j][0] + n3 * n4;
	  p[i + 1][j + 1][0][0] = p[i + 1][j][0][0] + n3 * n4 * n5;

	  for (k = 0; k < n3 - 1; k++)
	    p[i + 1][j + 1][0][k + 1] = p[i + 1][j + 1][0][k] + n5;

	  for (k = 0; k < n3 - 1; k++)
	    {
	      p[i + 1][j + 1][k + 1] = p[i + 1][j + 1][k] + n4;
	      p[i + 1][j + 1][k + 1][0] = p[i + 1][j + 1][k][0] + n4 * n5;
	      for (l = 0; l < n4 - 1; l++)
		p[i + 1][j + 1][k + 1][l + 1] = p[i + 1][j + 1][k + 1][l] + n5;
	    }
	}
    }

  for (i = 0, a = p[0][0][0][0]; i < n1 * n2 * n3 * n4 * n5; i++)
    *a++ = 0.0;

  return p;
}

/**
 * @brief Deallocation of #f5t() 
 *
 * @param[in,out] *****p              Array.
 *
 * @return \c void
 */
void free_f5t(float *****p){

  free(p[0][0][0][0]);
  free(p[0][0][0]);
  free(p[0][0]);
  free(p[0]);
  free(p);
}

/**
 * @brief Deallocation of #f4t() 
 *
 * @param[in,out]  ****p              Array.
 *
 * @return \c void
 */
void free_f4t(float ****p){

  free(p[0][0][0]);
  free(p[0][0]);
  free(p[0]);
  free(p);
}

/**
 * @brief Deallocation of #f3t() 
 *
 * @param[in,out]   ***p              Array.
 *
 * @return \c void
 */
void free_f3t(float ***p){

  free(p[0][0]);
  free(p[0]);
  free(p);
}

/**
 * @brief Deallocation of #f2t() 
 *
 * @param[in,out]    **p              Array.
 *
 * @return \c void
 */
void free_f2t(float **p){

  free(p[0]);
  free(p);
}

/**
 * @brief Deallocation of #f1t() 
 *
 * @param[in,out]     *p              Array.
 *
 * @return \c void
 */
void free_f1t(float *p){

  free(p);
}

/*********************************/

/**
 * @brief Allocation of 1D double array
 *
 * @param[in]        n1              Size of first dimension.
 *
 * @return pointer to double array
 */
double *d1t (int n1)
{
  double *p, *a;
  int i;
  if ((p = (double *) malloc ((size_t) n1 * sizeof (double))) == NULL)
      failed ("d1t: failed");
  for (i = 0, a = p; i < n1; i++)
    *a++ = 0;
  return p;
}

/**
 * @brief Allocation of 2D double array
 *
 * @param[in]        n1              Size of first dimension.
 * @param[in]        n2              Size of second dimension.
 *
 * @return pointer to double array
 */
double **d2t (int n1, int n2)
{
  double **p, *a;
  int i;
  if ((p = (double **) malloc ((size_t) n1 * sizeof (double *))) == NULL)
      failed ("d2t: failed n1");
  if ((p[0] = (double *) malloc ((size_t) n1 * n2 * sizeof (double))) == NULL)
      failed ("d2t: failed n2");
  for (i = 0; i < n1 - 1; i++)
    p[i + 1] = p[i] + n2;
  for (i = 0, a = p[0]; i < n1 * n2; i++)
    *a++ = 0;
  return p;
}

/**
 * @brief Allocation of 3D double array
 *
 * @param[in]        n1              Size of first dimension.
 * @param[in]        n2              Size of second dimension.
 * @param[in]        n3              Size of third dimension.
 *
 * @return pointer to double array
 */
double ***d3t (int n1, int n2, int n3)
{
  double ***p, *a;
  int i, j;
  if ((p = (double ***) malloc ((size_t) n1 * sizeof (double **))) == NULL)
      failed ("d3t: failed n1");
  if ((p[0] = (double **) malloc ((size_t) n1 * n2 * sizeof (double *))) == NULL)
      failed ("d3t: failed n2");
  if ((p[0][0] = (double *) malloc ((size_t) n1 * n2 * n3 * sizeof (double))) == NULL)
      failed ("d3t: failed n3");
  for (i = 0; i < n2 - 1; i++)
    p[0][i + 1] = p[0][i] + n3;
  for (i = 0; i < n1 - 1; i++)
    {
      p[i + 1] = p[i] + n2;
      p[i + 1][0] = p[i][0] + n2 * n3;
      for (j = 0; j < n2 - 1; j++)
	p[i + 1][j + 1] = p[i + 1][j] + n3;
    }
  for (i = 0, a = p[0][0]; i < n1 * n2 * n3; i++)
    *a++ = 0;
  return p;
}

/**
 * @brief Allocation of 4D double array
 *
 * @param[in]        n1              Size of first dimension.
 * @param[in]        n2              Size of second dimension.
 * @param[in]        n3              Size of third dimension.
 * @param[in]        n4              Size of fourth dimension.
 *
 * @return pointer to double array
 */
double ****d4t (int n1, int n2, int n3, int n4)
{
  double ****p, *a;
  int i, j, k;
  if ((p = (double ****) malloc ((size_t) n1 * sizeof (double ***))) == NULL)
      failed ("d4t: failed n1");
  if ((p[0] = (double ***) malloc ((size_t) n1 * n2 * sizeof (double **))) == NULL)
      failed ("d4t: failed n2");
  if ((p[0][0] = (double **) malloc ((size_t) n1 * n2 * n3 * sizeof (double *))) == NULL)
      failed ("d4t: failed n3");
  if ((p[0][0][0] = (double *) malloc ((size_t) n1 * n2 * n3 * n4 * sizeof (double ))) == NULL)
      failed ("d4t: failed n4");

  for (i = 0; i < n3 - 1; i++)
    p[0][0][i + 1] = p[0][0][i] + n4;

  for (i = 0; i < n2 - 1; i++)
    {
      p[0][i + 1] = p[0][i] + n3;
      p[0][i + 1][0] = p[0][i][0] + n3 * n4;
      for (j = 0; j < n3 - 1; j++)
	p[0][i + 1][j + 1] = p[0][i + 1][j] + n4;
    }

  for (i = 0; i < n1 - 1; i++)
    {
      p[i + 1] = p[i] + n2;
      p[i + 1][0] = p[i][0] + n2 * n3;
      p[i + 1][0][0] = p[i][0][0] + n2 * n3 * n4;

      for (j = 0; j < n3 - 1; j++)
	p[i + 1][0][j + 1] = p[i + 1][0][j] + n4;

      for (j = 0; j < n2 - 1; j++)
	{
	  p[i + 1][j + 1] = p[i + 1][j] + n3;
	  p[i + 1][j + 1][0] = p[i + 1][j][0] + n3 * n4;
	  for (k = 0; k < n3 - 1; k++)
	    p[i + 1][j + 1][k + 1] = p[i + 1][j + 1][k] + n4;
	}
    }

  for (i = 0, a = p[0][0][0]; i < n1 * n2 * n3 * n4; i++)
    *a++ = 0.0;

  return p;
}



/**
 * @brief Deallocation of #d4t() 
 *
 * @param[in,out]  ****p              Array.
 *
 * @return \c void
 */
void free_d4t(double ****p){

  free(p[0][0][0]);
  free(p[0][0]);
  free(p[0]);
  free(p);
}

/**
 * @brief Deallocation of #d3t() 
 *
 * @param[in,out]   ***p              Array.
 *
 * @return \c void
 */
void free_d3t(double ***p){

  free(p[0][0]);
  free(p[0]);
  free(p);
}

/**
 * @brief Deallocation of #d2t() 
 *
 * @param[in,out]    **p              Array.
 *
 * @return \c void
 */
void free_d2t(double **p){

  free(p[0]);
  free(p);
}

/**
 * @brief Deallocation of #d1t() 
 *
 * @param[in,out]     *p              Array.
 *
 * @return \c void
 */
void free_d1t(double *p){

  free(p);
}

/**********************************/


/**
 * @brief Function skip tu end of line in FILE
 *
 * @param[in,out]   *fp              Pointer to file which will got to end of line.
 *
 * @return \c void
 */
void readeol (FILE * fp)
{
  char s;
  while ((s = getc (fp)) != EOF)
    if (s == '\n')
      return;
}

/**
 * @brief Function calculate \f$a^n\f$ for int base
 *
 * @param[in]        a              Base.
 * @param[in]        n              Power.
 *
 * @return \f$a^n\f$
 */
int myipow (int a, int n)
{
  int b = a;
  while (--n > 0)
    a *= b;
  return a;
}

/**
 * @brief Function calculate \f$a^n\f$ for float base
 *
 * @param[in]        a              Base.
 * @param[in]        n              Power.
 *
 * @return \f$a^n\f$
 */
float myfpow (float a, int n)
{
  float b = a;
  while (--n > 0)
    a *= b;
  return a;
}

/**
 * @brief Function calculate \f$a^n\f$ for double base
 *
 * @param[in]        a              Base.
 * @param[in]        n              Power.
 *
 * @return \f$a^n\f$
 */
double mydpow (double a, int n)
{
  double b = a;
  while (--n > 0)
    a *= b;
  return a;
}

/**
 * @brief Print out double array with dimension (n, m)
 *
 * @param[in]        n              First dimension of array.
 * @param[in]        m              Second dimension of array.
 * @param[in]      **a              Array.
 *
 * @return \c void
 */
void pdarray (int n, int m, double **a)
{
  int i, j;
  printf ("\n");
  for (i = 0; i < n; i++)
    {
      for (j = 0; j < m; j++)
	printf ("%8.4f ", a[i][j]);
      printf ("\n");
    }
}

/**
 * @brief Print out double array with dimension (n)
 *
 * @param[in]        n              First dimension of array.
 * @param[in]      **a              Array.
 *
 * @return \c void
 */
void pdvector (int n, double *a)
{
  int i;
  printf ("\n");
  for (i = 0; i < n; i++)
    printf ("%8.4f ", a[i]);
  printf ("\n");
}

/**
 * @brief Print out float array with dimension (n, m)
 *
 * @param[in]        n              First dimension of array.
 * @param[in]        m              Second dimension of array.
 * @param[in]      **a              Array.
 *
 * @return \c void
 */
void pfarray (int n, int m, float **a)
{
  int i, j;
  printf ("\n");
  for (i = 0; i < n; i++)
    {
      for (j = 0; j < m; j++)
	printf ("%8.4f ", a[i][j]);
      printf ("\n");
    }
}

/**
 * @brief Print out float array with dimension (n)
 *
 * @param[in]        n              First dimension of array.
 * @param[in]      **a              Array.
 *
 * @return \c void
 */
void pfvector (int n, float *a)
{
  int i;
  printf ("\n");
  for (i = 0; i < n; i++)
    printf ("%8.4f ", a[i]);
  printf ("\n");
}

/**
 * @brief Print out int array with dimension (n, m)
 *
 * @param[in]        n              First dimension of array.
 * @param[in]        m              Second dimension of array.
 * @param[in]      **a              Array.
 *
 * @return \c void
 */
void piarray (int n, int m, int **a)
{
  int i, j;
  printf ("\n");
  for (i = 0; i < n; i++)
    {
      for (j = 0; j < m; j++)
	printf ("%8d ", a[i][j]);
      printf ("\n");
    }
}

/**
 * @brief Print out int array with dimension (n)
 *
 * @param[in]        n              First dimension of array.
 * @param[in]      **a              Array.
 *
 * @return \c void
 */
void pivector (int n, int *a)
{
  int i;
  printf ("\n");
  for (i = 0; i < n; i++)
    printf ("%8d ", a[i]);
  printf ("\n");
}

/**
 * @brief Set all values to zero in double array with dimension(n, m)
 *
 * @param[in]        n              First dimension of array.
 * @param[in]        m              Second dimension of array.
 * @param[in,out]  **a              Array.
 *
 * @return \c void
 */
void zdarray (int n, int m, double **a)
{
  int i, j;
  for (i = 0; i < n; i++)
    for (j = 0; j < m; j++)
      a[i][j] = 0.0;
}

/**
 * @brief Set all values to zero in double array with dimension(n)
 *
 * @param[in]        n              First dimension of array.
 * @param[in,out]  **a              Array.
 *
 * @return \c void
 */
void zdvector (int n, double *a)
{
  int i;
  for (i = 0; i < n; i++)
    a[i] = 0.0;
}

/**
 * @brief Set all values to zero in float array with dimension(n, m)
 *
 * @param[in]        n              First dimension of array.
 * @param[in]        m              Second dimension of array.
 * @param[in,out]  **a              Array.
 *
 * @return \c void
 */
void zfarray (int n, int m, float **a)
{
  int i, j;
  for (i = 0; i < n; i++)
    for (j = 0; j < m; j++)
      a[i][j] = 0.0;
}

/**
 * @brief Set all values to zero in float array with dimension(n)
 *
 * @param[in]        n              First dimension of array.
 * @param[in,out]  **a              Array.
 *
 * @return \c void
 */
void zfvector (int n, float *a)
{
  int i;
  for (i = 0; i < n; i++)
    a[i] = 0.0;
}

/**
 * @brief Set all values to zero in int array with dimension(n, m)
 *
 * @param[in]        n              First dimension of array.
 * @param[in]        m              Second dimension of array.
 * @param[in,out]  **a              Array.
 *
 * @return \c void
 */
void ziarray (int n, int m, int **a)
{
  int i, j;
  for (i = 0; i < n; i++)
    for (j = 0; j < m; j++)
      a[i][j] = 0;
}

/**
 * @brief Set all values to zero in int array with dimension(n)
 *
 * @param[in]        n              First dimension of array.
 * @param[in,out]  **a              Array.
 *
 * @return \c void
 */
void zivector (int n, int *a)
{
  int i;
  for (i = 0; i < n; i++)
    a[i] = 0;
}



/*********************************/





/**
 * @brief Function read 2D array of doubles from file
 *
 * --ADDED LT ---
 * If function fail to read all numbers end up with #failed().
 *
 * Function assume that numbers are each on new line.
 * a[0][0] = 1. line
 * a[0][1] = 2. line
 *  ...    = ...
 *  ...    = ...
 * a[1][0] = m-th line
 *  ...    = ...
 *  ...    = ...
 * a[n][m] = (n*m)-th line
 *
 * @param[in]          n              First dimension.
 * @param[in]          m              Second dimension.
 * @param[in,out]    **a              2D double array.
 * @param[in,out]     *F              File from where values are readed.
 *
 * @return \c void
 */
void rdarray(int n, int m, double **a, FILE *F)
{
	int i,j;
	int ret;
	char msg[1024];

	for(i = 0; i < n ; i++ )
	{
		for( j = 0; j < m ; j++)
		{
			ret=fscanf (F,"%lf",&a[i][j]);
			if(ret!=1)
			{
				sprintf(msg,
						"Failed reading double in _CATENR_rdarray_symm. ret= %d, i=%d, j=%d,n=%d\n",
						ret,i,j,n);
				failed(msg);
			}
		}
	}
}

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
 *
 * @param[in]          n              First and second dimension.
 * @param[in,out]    **a              2D double array.
 * @param[in,out]     *F              File from where values are readed.
 *
 * @return \c void
 */
void rdarray_symm (int n, double **a, FILE *F)
{
	int i,j;
	int ret;
	char msg[1024];

	for(i = 0; i < n ; i++ )
	{
		for( j = i; j < n ; j++)
		{
			ret=fscanf (F,"%lf",&a[i][j]);
			if(ret!=1)
			{
				sprintf(msg,
						"Failed reading double in _CATENR_rdarray_symm. ret= %d, i=%d, j=%d,n=%d\n",
						ret,i,j,n);
				failed(msg);
			}
			a[j][i]=a[i][j];
		}
	}
}

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
 *
 * @param[in]          n              Length of vector.
 * @param[in,out]    **a              1D double vector.
 * @param[in,out]     *F              File from where values are readed.
 *
 * @return \c void
 */
void rdvector (int n, double *a, FILE *F)
{
	int i;
	int ret;
	char msg[1024];

	for(i = 0; i < n ; i++ )
	{
		ret=fscanf (F,"%lf",&a[i]);
		if(ret!=1)
		{
				sprintf(msg,
						"Failed reading double in _CATENR_rdvector. ret= %d, i=%d,n=%d\n",
						ret,i,n);
				failed(msg);
		}
	}
}

