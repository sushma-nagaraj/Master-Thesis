/*
 * stat_utils.cpp
 *
 *  Created on: Mar 14, 2016
 *      Author: sushma_n
 */

/*************************************************************

   Filename :  stat_utils.c

   Purpose:    performs several statistical basic functions used in the
               DSPN Simulation modul
	       >
	       >

   Author :    Christian Kelling

   Date :      May 1993

*************************************************************/

#include <stdio.h>
#include <math.h>
#include "stat_utils.h"

/* windows */
#ifndef M_PI
#define M_PI 3.14159265
#endif

/* constant definitions for normal_distr */
/* see [Abramowitz, Stegun]: "Handbook of mathematical functions" */
/* page 933, Formula 26.2.23 */
static double C_const[3] = {	2.515517,
	                      	0.802853,
        	              	0.010328
                	      	};

static double D_const[4] = {	1.0,
				1.432788,
                      		0.189269,
                      		0.001308
                      		};

/* constant definition for student_distr */
/* see [Abramowitz, Stegun]: "Handbook of mathematical functions" */
/* page 990, Table 26.10 */
/* STUDENT [degree of freedom: 1..30][1-alpha: 0.90; 0.95; 0.98; 0.99; 0.995; 0.999] */
static double STUDENT_const[30][6] = {{6.314, 12.706, 31.821, 63.657, 127.321, 636.619},
				      {2.920,  4.303,  6.965,  9.925,  14.089,  31.598},
				      {2.353,  3.182,  4.541,  5.841,   7.453,  12.924},
				      {2.132,  2.776,  3.747,  4.604,   5.598,   8.610},
				      {2.015,  2.571,  3.365,  4.032,   4.773,   6.869},

				      {1.943,  2.447,  3.143,  3.707,   4.317,   5.959},
				      {1.895,  2.365,  2.998,  3.499,   4.029,   5.408},
				      {1.860,  2.306,  2.896,  3.355,   3.833,   5.041},
				      {1.833,  2.262,  2.821,  3.250,   3.690,   4.781},
				      {1.812,  2.228,  2.764,  3.169,   3.581,   4.587},

				      {1.796,  2.201,  2.718,  3.106,   3.497,   4.437},
				      {1.782,  2.179,  2.681,  3.055,   3.428,   4.318},
				      {1.771,  2.160,  2.650,  3.012,   3.372,   4.221},
				      {1.761,  2.145,  2.624,  2.977,   3.326,   4.140},
				      {1.753,  2.131,  2.602,  2.947,   3.286,   4.073},

				      {1.746,  2.120,  2.583,  2.921,   3.252,   4.015},
				      {1.740,  2.110,  2.567,  2.898,   3.223,   3.965},
				      {1.734,  2.101,  2.552,  2.878,   3.197,   3.922},
				      {1.729,  2.093,  2.539,  2.861,   3.174,   3.883},
				      {1.725,  2.086,  2.528,  2.845,   3.153,   3.850},

				      {1.721,  2.080,  2.518,  2.831,   3.135,   3.819},
				      {1.717,  2.074,  2.508,  2.819,   3.119,   3.792},
				      {1.714,  2.069,  2.500,  2.807,   3.104,   3.768},
				      {1.711,  2.064,  2.492,  2.797,   3.090,   3.745},
				      {1.708,  2.060,  2.485,  2.787,   3.078,   3.725},

				      {1.706,  2.056,  2.479,  2.779,   3.067,   3.707},
				      {1.703,  2.052,  2.473,  2.771,   3.057,   3.690},
				      {1.701,  2.048,  2.467,  2.763,   3.047,   3.674},
				      {1.699,  2.045,  2.462,  2.756,   3.038,   3.659},
				      {1.697,  2.042,  2.457,  2.750,   3.030,   3.646}
				      };


/*************************************************************************

  function:	means

  parameters:	sequence -> pointer to linked list of data
		count	 -> number of data to calculate mean value

  returns:	mean value of data

  purpose:	calculates mean value for a number of data

  author:	Christian Kelling
		Technical University of Berlin

  date:		May 1993

*************************************************************************/

double means(int count, SampSeq *sequence)
{
  int n;
  double mean;

  mean = 0.0;
  n = count;
  while (n--)
    {
      mean += sequence->Value;
      if (sequence->Next == NULL)
	{
	  printf ("Error in function means, counter = %d",n);
	  printf ("insufficient length of sequence");
	}
      sequence = sequence->Next;
    }
  return (mean/count);
} /* means */

/*************************************************************************

  function:	normal_distr

  parameters:	p -> probability value

  returns:	p-quantile

  purpose:	calculates the p-quantile (percentile) of the normal-
		distribution applying a numerical approximation.
		See [Abramowitz, Stegun]: "Handbook of mathematical functions",
		page 933, formula 26.2.23
		this function assumes the one-sided questionary !!!

  author:	Rainer Hagenau
		Technical University of Berlin

  date:		May 1991

*************************************************************************/

double normal_distr(double p)
{
  double t, zaehler, nenner;
  int i;

  t = sqrt(log(1.0/pow(p,2.0)));
  for (zaehler = 0.0, i=0; i<3; i++)
    zaehler += C_const[i] * pow(t,(double)i);
  for (nenner = 0.0, i=0; i<4; i++)
    nenner += D_const[i] * pow(t,(double)i);
  return(t - zaehler/nenner);
} /* normal_distr */


/*************************************************************************

  function:	student_distr

  parameters:	free  -> degree of freedom
		level -> probability value

  returns:	p-quantile

  purpose:	calculates the p-quantile (percentile) of the t-distribution
		with n := free degrees of freedom.
		For values of free <= 30 the function uses values of a table
		and calculates by lineare interpolation for values outside the
		table. If the degree of freedom is greater than 30, it will
		approximate by using the normal-distribution.

  author:	Rainer Hagenau
		Technical University of Berlin

  date:		May 1991

*************************************************************************/

double student_distr(int free, double level)
{
  /* values of table applied by handbook	*/
  /* level := 1-alpha			*/
  int i;
  double x, delta;
  void nrerror(char []);

  if  (free <= 0)
    nrerror ("t_student: invalid degree of freedom.");
  if (free > 30)
    return(normal_distr((1-level)/2.0));  /* this is the transformation from
	     the two-sided questionary (assumed by the student_distr) into the
	     two-sided questionary assumed by normal_distr */
  else
    {
      switch((int)fabs(level*1000))
	{
	case 900:
	  return(STUDENT_const[free-1][0]);
	case 950:
	  return(STUDENT_const[free-1][1]);
	case 980:
	  return(STUDENT_const[free-1][2]);
	case 990:
	  return(STUDENT_const[free-1][3]);
	case 995:
	  return(STUDENT_const[free-1][4]);
	case 999:
	  return(STUDENT_const[free-1][5]);

	default:
	  if (level < 0.9 || level > 0.999)
	    printf("student_distr: unknown level [%f] by degree of freedom %d",
		   level, free);
	  else
	    {
	      /* evaluating by linear interpolation */
	      if (level > 0.995)
		{
		  i = 4;
		  delta = (level-0.995)/0.004;
		}
	      else if (level > 0.990)
		{
		  i = 3;
		  delta = (level-0.990)/0.005;
		}
	      else if (level > 0.980)
		{
		  i = 2;
		  delta = (level-0.980)/0.01;
		}
	      else if (level > 0.950)
		{
		  i = 1;
		  delta = (level-0.950)/0.03;
		}
	      else
		{
		  i = 0;
		  delta = (level-0.900)/0.05;
		}
	      x = (STUDENT_const[free-1][i+1] -
		   STUDENT_const[free-1][i]) * delta;
	      x += STUDENT_const[free-1][i];
	      return(x);
	    }
	  break;
	}
    }
	return(0.0);
} /* student_distr */


/*************************************************************************

  function:	periodogram

  parameters:	xarray -> Array of observations x1, x2, .. xn
		n      -> number of observations
		j      -> index of the discrete Fourier transforms Ax(j)

  returns:	value of periodogram

  purpose:	calculates the periodogram for a number of observations.
		See [Pawlikowsky]: Equations 45-46, p.138 .

  author:	Rainer Hagenau
		Technical University of Berlin

  date:		May 1991

*************************************************************************/

double periodogram (double *xarray, int j, int n)
{
  double angle, real, imag;
  int s;

  real = 0.0;
  imag = 0.0;
   for (s=1; s<=n; s++)
    {
      angle = -2.0 * M_PI * (s-1) * j/n;
      real += xarray[s] * cos(angle);
      imag += xarray[s] * sin(angle);
    }
  return((real*real + imag*imag)/n);
} /* periodogram */


/*************************************************************************

  function:	SpectralVarAnalysis

  parameters:	testsamplesarray -> array of data to be tested
		arraylen	 -> number of data

  returns:	estimated variance of a sequence of observations

  purpose:	This function determines the variance of the
		that is needed by Schrubens stationarity test
		and by the Soectral Variance Analysis
		The least squares extrapolition procedure is
		from: Press, W.H. et al."Numerical Recipes in C,
		  Cambridge, University Press,1992

		For further information to used functions (fit_poly_norm etc.)
		see num_utils.c

  author:	Christian Kelling
		Technical University of Berlin

  date:		May 1993

*************************************************************************/

double SpectralVarAnalysis (double *testsamplesarray, int arraylen)
{
  int i, nap, npol;
  double normcon, **x, *y, *a, *period;
  double periodogram(), *dvector(), **dmatrix();
  double *fit_poly_norm();

  nap = 25;
  npol = 3;
  normcon = 0.882;
  x=dmatrix(1,nap,1,1);
  y=dvector(1,nap);
  a=dvector(1,npol);
  period = dvector(1,2*nap);
  for (i=1;i<=(2*nap);i++)
    period[i] = periodogram(testsamplesarray, i, arraylen);
  for (i=1;i<=nap;i++)
    {
      x[i][1] = (4.0*i-1)/(2.0*arraylen);
      y[i] = log((period[2*i-1]+period[2*i])/2.0)+0.27;
    }
  a = fit_poly_norm (x, y, nap, npol);
  return (normcon*exp(a[1])/arraylen);
} /* SpectralVarAnalysis */


