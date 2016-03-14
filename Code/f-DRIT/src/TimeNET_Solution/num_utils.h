/*
 * num_utils.h
 *
 *  Created on: Mar 14, 2016
 *      Author: sushma_n
 */

#ifndef TIMENET_SOLUTION_NUM_UTILS_H_
#define TIMENET_SOLUTION_NUM_UTILS_H_

/************************************************************

  Filename  :   num_utils.h

  Purpose   :   header file for num_utils.c

  Author    :   Christian Kelling

  Date      :   May  1993

************************************************************/



#define  SUCCEED    0
#define  FAIL       1

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define SQR(a) (a == 0.0 ? 0.0 : a*a)

#define TINY  1.0e-16
#define TOL   1.0e-5

#ifndef min
  #define min(a,b) (((a)>(b))?(b):(a))
#endif
#ifndef max
  #define max(a,b) (((a)>(b))?(a):(b))
#endif


void nrerror();          /* default function for errors      */
int *ivector();          /* int and double                   */
void free_ivector();
void free_dvector();     /* functions for freeing memory     */
double *dvector();       /* functions for reserving memory   */
int **imatrix();
void free_imatrix();
double **dmatrix();      /* for vectors and matrices of type */
void free_dmatrix();
double pythag();
void fpoly();            /* function for polynomial functions*/
double **mult_matr();    /* calculates the product of two matrices */
double **transpose();
void ludcmp();           /* LU (lower/upper decomposition of a matrix */
void lubksb();           /* solve A*x=B for an output matix of ludcmp */
void solve_LES();        /* solving the set of linear equations using
			    LU decomposition and backward substiution */
double **inverse();      /* computes the inverse of a matrix */
void svbksb();           /* prepares matrices for next       */
                         /* function                         */
void svdcmp();           /* singular value decomposition     */
void svdfit();           /* least square method              */
double *fit_poly_norm ();  /* least squares polynomial fit  for normal
			      distributed data */
double *fit_multi_lin_norm ();  /* least squares fit  for normal
				   distributed multidimensional data */










#endif /* TIMENET_SOLUTION_NUM_UTILS_H_ */
