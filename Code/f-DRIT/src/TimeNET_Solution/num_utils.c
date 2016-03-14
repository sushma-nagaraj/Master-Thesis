/*
 * num_utils.c
 *
 *  Created on: Mar 14, 2016
 *      Author: sushma_n
 */
/************************************************************
   Filename :  num_utils.c

   Purpose:    performs several computations used in the DSPN Simulation modul
                 > computes LU decomposition of a matrix (ludcmp)
		 > computes least squares interpolation (svdfit) using
		   sigular value decomposition
	       see Numerical Recipes in C; Press, W.H.,
	       cambridge university press, 1992

   Author :    Christian Kelling

   Date :      May 1993

*************************************************************/

#include <stdio.h>
#include <math.h>
#include "num_utils.h"

#include <stdlib.h>

//extern char *er_malloc();

/*****************************************************************************

utility-functions for error-handling and

                  memory managment

*****************************************************************************/


/*****************************************************************************

  standard error handler

*****************************************************************************/

void nrerror(char *error_text)
{
	fprintf(stderr,"%s\n",error_text);
	exit(FAIL);
}


/*****************************************************************************

  allocates an int vector with subscript range v[nl..nh]

*****************************************************************************/

int *ivector(int nl, int nh)
{
	int *v;

	v=(int *)er_malloc((unsigned) (nh-nl+1)*sizeof(int));
	if (!v) nrerror("allocation failure in ivector()");
	return v-nl;
}


/*****************************************************************************

  free an int vector allocated with ivector

*****************************************************************************/

void free_ivector(int *v, int nl, int nh)
{
	free((char*) (v+nl));
}


/*****************************************************************************

  allocates a double vector with subscript range v[nl..nh]

*****************************************************************************/

double *dvector(int nl, int nh)
{
	double *v;

	v=(double *)er_malloc((unsigned) (nh-nl+1)*sizeof(double));
	if (!v) nrerror("allocation failure in dvector()");
	return v-nl;
}

/*****************************************************************************

  free a double vector allocated with ivector

*****************************************************************************/

void free_dvector(double *v, int nl, int nh)
{
	free((char*) (v+nl));
}

/*****************************************************************************

  allocates an int matrix with subscript range m[nrl..nrh][ncl..nch]

*****************************************************************************/

int **imatrix(int nrl, int nrh, int ncl, int nch)
{
	int i,**m;

	m=(int **)er_malloc((unsigned) (nrh-nrl+1)*sizeof(int*));
	if (!m) nrerror("allocation failure 1 in imatrix()");
	m -= nrl;

	for(i=nrl;i<=nrh;i++) {
		m[i]=(int *)er_malloc((unsigned) (nch-ncl+1)*sizeof(int));
		if (!m[i]) nrerror("allocation failure 2 in imatrix()");
		m[i] -= ncl;
	}
	return m;
}


/*****************************************************************************

  free an int matrix allocated with imatrix

*****************************************************************************/

void free_imatrix(int **m, int nrl, int nrh, int ncl, int nch)
{
	int i;

	for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
	free((char*) (m+nrl));
}


/*****************************************************************************

  allocates an double matrix with subscript range m[nrl..nrh][ncl..nch]

*****************************************************************************/

double **dmatrix(int nrl, int nrh, int ncl, int nch)
{
	int i;
	double **m;

	m=(double **) er_malloc((unsigned) (nrh-nrl+1)*sizeof(double*));
	if (!m) nrerror("allocation failure 1 in dmatrix()");
	m -= nrl;

	for(i=nrl;i<=nrh;i++) {
		m[i]=(double *) er_malloc((unsigned) (nch-ncl+1)*sizeof(double));
		if (!m[i]) nrerror("allocation failure 2 in dmatrix()");
		m[i] -= ncl;
	}
	return m;
}


/*****************************************************************************

  free a double matrix allocated with dmatrix

*****************************************************************************/

void free_dmatrix(double **m, int nrl, int nrh, int ncl, int nch)
{
	int i;

	for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
	free((char*) (m+nrl));
}


/*****************************************************************************

  computes sqrt(a^2+b^2) with overflow

*****************************************************************************/

double pythag(double a, double b)
{
  double absa,absb;
  absa = fabs(a);
  absb = fabs(b);
  if (absa > absb)
    return absa*sqrt(1.0+SQR(absb/absa));
  else
    return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb)));
}


/*****************************************************************************

  end of utility-functions

*****************************************************************************/


/*****************************************************************************

  service functions used in the function svdfit to determine the so called
  basis funtions X(x) which can be arbitrary functions of x

  we implemented the functions:

  fpoly (fitting for a polynomial of degree n-1 with coefficients p[1..np])

  fmultdim (fitting to the simple multidimensional case X[k](x[i])=x[i][k]
            k= 1...number_of_variates; i=1...number_of_samples)
	    see: Press, Numerical Recipes in C (1992), p. 680
*****************************************************************************/

void fpoly(double *x, double *p, int np)
{
	int j;

	p[1]=1.0;
	for (j=2;j<=np;j++) p[j]=x[1]*p[j-1];
}

void fmultdim(double *x, double *p, int np)
{
  int i;

  p[0] = 1.0;
  for (i=1;i<=np;i++)
      p[i] = x[i];

}

/*****************************************************************************

  end of service functions for svdfit

*****************************************************************************/


/*************************************************************************

  function:	mult_matr

  parameters:	a -> pointer to pointer to an m*n-matrix a
                b -> pointer to pointer to an  n*p-matrix b
                m,n,p -> dimensions of the system

  returns:      pointer to pointer to the matrix a*b

  purpose:	Computes the product Y=AxB

  author:	Christian Kelling
		Technical University of Berlin

  date:		May 1993

*************************************************************************/

double **mult_matr(double **a, double **b, int m, int n, int p)
{
  int i, j, k;
  double **y;

  y = dmatrix(1,m,1,p);
  for (i=1;i<=m;i++)
    for (j=1;j<=p;j++) y[i][j] = 0.0;
  for (i=1;i<=m;i++)
    for (j=1;j<=p;j++)
      for (k=1;k<=n;k++) y[i][j] += a[i][k]*b[k][j];
  return y;
}


/*************************************************************************

  function:	transpose

  parameters:	a -> pointer to pointer to an m*n-matrix a
                m,n -> dimensions of the system

  returns:      pointer to pointer to the transpose of matrix a

  purpose:      transpose a matrix A

  author:	Christian Kelling
		Technical University of Berlin

  date:		May 1993

*************************************************************************/

double **transpose(double **a, int m, int n)
{
  int i, j;
  double **t, **dmatrix();

  t = dmatrix(1,n,1,m);
  for (i=1;i<=m;i++)
    for (j=1;j<=n;j++)
      t[j][i] = a[i][j];
  return t;
}


/*************************************************************************

  function:	ludcmp

  parameters:	a -> pointer to pointer to an matrix
                n -> dimension of a

  returns:      indx[1..n] -> pointer to the row permutation effected by
                              partial pivoting
                d -> pointer to double d, which denotes wether:
		     +1 number of row interchanges is even
		     -1 number of row interchanges is odd

  purpose:	replaces the matrix a by its LU (lower/upper) decomposition
                This routine is used in combination with lubksb to solve
		linear equations or invert a matrix

  author:	see: Press, Numerical Recipes in C (1992) pp. 46--47
                Christian Kelling
		Technical University of Berlin

  date:		May 1993

*************************************************************************/

void ludcmp (double **a, int n, int *indx, double *d)
{
  int i, j, k, imax = 0;
  double big, dum, sum, temp;
  double *vv;
  void free_dvector();

  vv = dvector(1,n);
  *d = 1.0;
  for (i=1; i<=n; i++) /* loop over rows to get the implicit scaling info*/
    {
      big = 0.0;
      for (j=1; j<=n; j++)
	if ((temp = fabs(a[i][j])) > big) big = temp;
      if (big == 0.0) nrerror("Singular matrix in ludcmp");
      vv[i] = 1.0/big;  /* save the scaling*/
    }
  for (j=1; j<=n; j++)  /* loop over columns of Crouts method(see Press) */
    {
      for (i=1; i<j; i++)
	{
	  sum = a[i][j];
	  for (k=1; k<i; k++) sum -= a[i][k]*a[k][j];
	  a[i][j] = sum;
	}
      big = 0.0;  /* init search for the largest pivot element */
      for (i=j; i<=n; i++)
	{
	  sum = a[i][j];
	  for (k=1; k<j; k++) sum -= a[i][k]*a[k][j];
	  a[i][j] = sum;
	  if ((dum = vv[i]*fabs(sum)) >= big)  /* is the figure of merit */
	    {  				       /* for the pivot better */
	      big = dum;		       /* than the best so far */
	      imax = i;
	    }
	}
      if (j != imax)  /* do we need interchange rows ? */
	{
	  for (k=1; k<=n; k++)  /* interchange rows */
	    {
	      dum = a[imax][k];
	      a[imax][k] = a[j][k];
	      a[j][k] = dum;
	    }
	  *d = -(*d);  /* even/odd interchanges */
	  vv[imax] = vv[j];  /* interchange the scale factor */
	}
      indx[j] = imax;
      if (a[j][j] == 0.0) a[j][j] = TINY;
      /* if the pivot element is zero, the matrix is singular */
      if ( j != n)  /* now finally divide by the pivot element */
	{
	  dum = 1.0/a[j][j];
	  for (i=j+1;i<=n; i++) a[i][j] *= dum;
	}
    }  /* go back for the next column in the reduction */
  free_dvector (vv,1,n);
}


/*************************************************************************

  function:	lubksb

  parameters:	a -> pointer to pointer to the LU decomposition of a
                     matrix a, determined by ludcmp
	        indx[1..n] -> pointer to the row permutation vector,
		              returned by ludcmp
		b[1..n] -> pointer to the right hand side vector B
                n -> dimension of the system

  returns:      b[1..n] -> solution vector X

  purpose:	Solves the set of n linear equations A*X = B, using the
                backward substitution algorithm

  author:	see: Press, Numerical Recipes in C (1992), pp. 47
                Christian Kelling
		Technical University of Berlin

  date:		May 1993

*************************************************************************/

void lubksb (double **a, int n, int *indx, double *b)
{
  int  i, ii=0, ip, j;
  double sum;

  for (i=1; i<=n; i++)
    {
      ip = indx[i];
      sum = b[ip];
      b[ip] = b[i];
      if (ii)  /* when ii is positive, it becomes the index of the first
		  non-vanishing element, then forward substitution */
	for (j=ii; j<=i-1; j++) sum -= a[i][j]*b[j];
      else  /* a non zero element was encountered, do the sums in loop above */
	if (sum) ii = i;
      b[i] = sum;
    }
  for (i=n; i>=1; i--)  /* backward substitution */
    {
      sum = b[i];
      for (j=i+1; j<=n; j++) sum -= a[i][j]*b[j];
      b[i] = sum/a[i][i];  /* store an element of the solution vector */
    }
}

/*************************************************************************

  function:	solve_LES

  parameters:	a -> pointer to pointer to an matrix a
		b[1..n] -> pointer to the right hand side vector B
                n -> dimension of the system

  returns:      b[1..n] -> pointer to the solution vector X

  purpose:	Solves the set of n linear equations A*X = B, using the
                procedures ludcmp and lubksb

  author:	see: Press, Numerical Recipes in C (1992), pp. 48
                Christian Kelling
		Technical University of Berlin

  date:		May 1993

*************************************************************************/

void solve_LES (double **a, double *b, int n)
{
  double d;
  int *indx, *ivector();
  void ludcmp(), lubksb();

  indx = ivector(1,n);

  ludcmp(a,n,indx,&d);
  lubksb(a,n,indx,b);

}

/*************************************************************************

  function:	inverse

  parameters:	a -> pointer to pointer to an matrix a
                n -> dimension of the system

  returns:      y -> pointer to pointer to the inverse of a (matrix y)

  purpose:	Computes the inverse of an matrix a
                procedures ludcmp and lubksb

  author:	see: Press, Numerical Recipes in C (1992) pp. 48
                Christian Kelling
		Technical University of Berlin

  date:		May 1993

*************************************************************************/

double **inverse (double **a, int n)
{
  double  **y, d, *col, *dvector();
  int  *indx, i, j, *ivector();
  void  ludcmp(), lubksb();

  col = dvector (1,n);
  indx = ivector (1,n);
  y = dmatrix (1,n,1,n);

  ludcmp (a,n,indx,&d);
  for (j=1;j<=n;j++)
    {
      for (i=1;i<=n;i++) col[i] = 0.0;
      col[j] = 1.0;
      lubksb(a,n,indx,col);
      for (i=1;i<=n;i++) y[i][j] = col[i];
    }
  free_dvector (col,1,n);
  free_ivector (indx,1,n);
  return y;
}

/*************************************************************************

  function:	svdcmp

  parameters:	 a -> pointer to pointer to an m*n matrix a
                 m,n -> dimensions of the matrix

  returns:      a -> pointer to pointer to the result matrix U of the singular
                     value decomposition
		v -> pointer to pointer to the n*n result matrix V (not the
		     transpose!) of the singular value decomposition
		w -> vector of the diagonal matrix of singular values W as a
		     vector [n], all other elements of W are zero

  purpose:	computes the singular value decomposition of a given matrix A
                A=U*W*(transpose)V

  author:	see: Press, Numerical Recipies in C (1992), pp.59--70
                Christian Kelling
		Technical University of Berlin

  date:		May 1993

*************************************************************************/

void svdcmp(double **a, int m, int n, double *w, double **v)
{
  int flag,i,its,j,jj,k,l,nm;
  double c,f,h,s,x,y,z;
  double anorm=0.0,g=0.0,scale=0.0;
  double *rv1,*dvector(),pythag();
  void nrerror(),free_dvector();

  l = 0;
  nm = 0;

  rv1=dvector(1,n);

  for (i=1;i<=n;i++) {
    l=i+1;
    rv1[i]=scale*g;
    g=s=scale=0.0;
    if (i <= m)
      {
	for (k=i;k<=m;k++) scale += fabs(a[k][i]);
	if (scale)
	  {
	    for (k=i;k<=m;k++)
	      {
		a[k][i] /= scale;
		s += a[k][i]*a[k][i];
	      }
	    f=a[i][i];
	    g = -SIGN(sqrt(s),f);
	    h=f*g-s;
	    a[i][i]=f-g;
	    for (j=l;j<=n;j++)
	      {
		for (s=0.0,k=i;k<=m;k++) s += a[k][i]*a[k][j];
		f=s/h;
		for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
	      }
	    for (k=i;k<=m;k++) a[k][i] *= scale;
	  }
      }
    w[i]=scale*g;
    g=s=scale=0.0;
    if (i <= m && i != n)
      {
	for (k=l;k<=n;k++) scale += fabs(a[i][k]);
	if (scale)
	  {
	    for (k=l;k<=n;k++)
	      {
		a[i][k] /= scale;
		s += a[i][k]*a[i][k];
	      }
	    f=a[i][l];
	    g = -SIGN(sqrt(s),f);
	    h=f*g-s;
	    a[i][l]=f-g;
	    for (k=l;k<=n;k++) rv1[k]=a[i][k]/h;
	    for (j=l;j<=m;j++)
	      {
		for (s=0.0,k=l;k<=n;k++) s += a[j][k]*a[i][k];
		for (k=l;k<=n;k++) a[j][k] += s*rv1[k];
	      }
	    for (k=l;k<=n;k++) a[i][k] *= scale;
	  }
      }
    anorm=max(anorm,(fabs(w[i])+fabs(rv1[i])));
  }
  for (i=n;i>=1;i--)
    {
      if (i < n)
	{
	  if (g)
	    {
	      for (j=l;j<=n;j++)
		  v[j][i]=(a[i][j]/a[i][l])/g;
	      for (j=l;j<=n;j++)
		{
		  for (s=0.0,k=l;k<=n;k++) s += a[i][k]*v[k][j];
		  for (k=l;k<=n;k++) v[k][j] += s*v[k][i];
		}
	    }
	  for (j=l;j<=n;j++) v[i][j]=v[j][i]=0.0;
	}
      v[i][i]=1.0;
      g=rv1[i];
      l=i;
    }
  for (i=min(m,n);i>=1;i--)
    {
      l=i+1;
      g=w[i];
      for (j=l;j<=n;j++) a[i][j]=0.0;
      if (g)
	{
	  g=1.0/g;
	  for (j=l;j<=n;j++)
	    {
	      for (s=0.0,k=l;k<=m;k++) s += a[k][i]*a[k][j];
	      f=(s/a[i][i])*g;
	      for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
	    }
	  for (j=i;j<=m;j++) a[j][i] *= g;
	}
      else
	for (j=i;j<=m;j++) a[j][i]=0.0;
      ++a[i][i];
    }
  for (k=n;k>=1;k--)
    {
      for (its=1;its<=30;its++)
       {
	  flag=1;
	  for (l=k;l>=1;l--)
	    {
	      nm=l-1;
	      if (fabs(rv1[l])+anorm == anorm)
		{
		  flag=0;
		  break;
		}
	      if (fabs(w[nm])+anorm == anorm) break;
	    }
	  if (flag)
	    {
	      c=0.0;
	      s=1.0;
	      for (i=l;i<=k;i++)
		{
		  f=s*rv1[i];
		  rv1[i]=c*rv1[i];
		  if (fabs(f)+anorm == anorm) break;
		  g=w[i];
		  h=pythag(f,g);
		  w[i]=h;
		  h=1.0/h;
		  c=g*h;
		  s=(-f*h);
		  for (j=1;j<=m;j++)
		    {
		      y=a[j][nm];
		      z=a[j][i];
		      a[j][nm]=y*c+z*s;
		      a[j][i]=z*c-y*s;
		    }
		}
	    }
	  z=w[k];
	  if (l == k)
	    {
	      if (z < 0.0)
		{
		  w[k] = -z;
		  for (j=1;j<=n;j++)
		    {
		      v[j][k]=(-v[j][k]);
		    }
		}
	      break;
	    }
	  if (its == 30)
	    nrerror("No convergence in 30 SVDCMP iterations");
	  x=w[l];
	  nm=k-1;
	  y=w[nm];
	  g=rv1[nm];
	  h=rv1[k];
	  f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
	  g=pythag(f,1.0);
	  f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
	  c=s=1.0;
	  for (j=l;j<=nm;j++)
	    {
	      i=j+1;
	      g=rv1[i];
	      y=w[i];
	      h=s*g;
	      g=c*g;
	      z=pythag(f,h);
	      rv1[j]=z;
	      c=f/z;
	      s=h/z;
	      f=x*c+g*s;
	      g=g*c-x*s;
	      h=y*s;
	      y *= c;
	      for (jj=1;jj<=n;jj++)
		{
		  x=v[jj][j];
		  z=v[jj][i];
		  v[jj][j]=x*c+z*s;
		  v[jj][i]=z*c-x*s;
		}
	      z=pythag(f,h);
	      w[j]=z;
	      if (z)
		{
		  z=1.0/z;
		  c=f*z;
		  s=h*z;
		}
	      f=c*g+s*y;
	      x=c*y-s*g;
	      for (jj=1;jj<=m;jj++)
		{
		  y=a[jj][j];
		  z=a[jj][i];
		  a[jj][j]=y*c+z*s;
		  a[jj][i]=z*c-y*s;
		}
	    }
	  rv1[l]=0.0;
	  rv1[k]=f;
	  w[k]=x;
	}  /* end for its */
    } /* end for k  */
  free_dvector(rv1,1,n);
}

/*************************************************************************

  function:	svbksb

  parameters:	 u -> pointer to pointer to an m*n matrix U
                 v -> pointer to pointer to an n*n matrix V
		 w -> array [n] for diagonal matrix W
		 b -> array [m] for the right hand side input B
                 m,n -> dimensions of the matrix

  returns:      x -> array [n] for the solution vector x

  purpose:	solves A*X=B for a given vector x using the singular value
                decomposition of svdcmp

  author:	see: Press, Numerical Recipies in C (1992), pp.64
                Christian Kelling
		Technical University of Berlin

  date:		May 1993

*************************************************************************/

void svbksb(double **u, double *w, double **v, int m, int n, double *b, double *x)
{
  int jj,j,i;
  double s,*tmp,*dvector();
  void free_dvector();

  tmp=dvector(1,n);
  for (j=1;j<=n;j++)
    {  /* calculate U(transpose)B */
      s=0.0;
      if (w[j])
	{
	  for (i=1;i<=m;i++) s += u[i][j]*b[i];
	  s /= w[j];
	}
      tmp[j]=s;
    }
  for (j=1;j<=n;j++)
    {
      s=0.0;
      for (jj=1;jj<=n;jj++) s += v[j][jj]*tmp[jj];
      x[j]=s;
    }
  free_dvector(tmp,1,n);
}


/*************************************************************************

  function:	svdfit

  parameters:	x -> array[ndata] of x-data points
                y -> array[ndata] of y-data points
		sig -> array[ndata] of individual standard deviations
		ndata -> number of data points
                ma -> number of fit parameters

		u -> workspace for the singular value decomposition
		     (ndata*ma)-matrix U of the (ndata*ma) input matrix
		v -> workspace for the singular value decomposition
		     (ma*ma)-matrix V of the (ndata*ma) input matrix
		w -> workspace for the singular value decomposition
                     (ma*ma)-matrix W of the (ndata*ma) input matrix, stored
		     as vector[ma] (all off-diagonal elements are zero)

  returns:      a -> array[ma] for the fit parameters a
                chisq -> minimized chi sqare values

  purpose:	determine the coefficient of the fitting function
                y = sum(i)(a[i]*afunc(i)(x))
		the user supplies a routine funcs(x,afunc,ma) that returns
		the basis functions in the array afunc

  author:	see: Press, Numerical Recipies in C (1992), pp.676--679
        	Christian Kelling
		Technical University of Berlin

  remarks:      he definition of x has been slightly modified in order to
                make possible multidimensional fits

  date:		May 1993

**********************************************************************v***/


void svdfit(double **x, double *y, double *sig, int ndata, double *a, int ma,
			double **u, double **v, double *w, double *chisq,
			void (*funcs)(double *,double *,int))
{
  int j,i;
  double wmax,tmp,thresh,sum,*b,*afunc,*dvector();
  void svdcmp(),svbksb(),free_dvector();

  b=dvector(1,ndata);
  afunc=dvector(1,ma);

  for (i=1;i<=ndata;i++)
    {  /* accumulate coefficients of the fitting matrix */
      (*funcs)(x[i],afunc,ma);
      tmp=1.0/sig[i];
      for (j=1;j<=ma;j++) u[i][j]=afunc[j]*tmp;
      b[i]=y[i]*tmp;
    }
  svdcmp(u,ndata,ma,w,v);
  wmax=0.0;
  for (j=1;j<=ma;j++)
    if (w[j] > wmax) wmax=w[j];
  thresh=TOL*wmax;
  for (j=1;j<=ma;j++)
    if (w[j] < thresh) w[j]=0.0;
  svbksb(u,w,v,ndata,ma,b,a);
  *chisq=0.0;
  for (i=1;i<=ndata;i++)
    {
      (*funcs)(x[i],afunc,ma);
      for (sum=0.0,j=1;j<=ma;j++) sum += a[j]*afunc[j];
      *chisq += (tmp=(y[i]-sum)/sig[i],tmp*tmp);
    }
  free_dvector(afunc,1,ma);
  free_dvector(b,1,ndata);
}


/*************************************************************************

  function:	fit_poly_norm

  parameters:	x -> pointer to pointer to an array[ndata][nvariates]
                     of x-data points (nvariates=1 for polynomial fit)
                y -> array[ndata] of y-data points
		ndata -> number of data points
		npol -> degree of the fitting polynom +1

  returns:      pointer to the  array[npol] for the fit parameters a

  purpose:	fits a polynom of degree npol-1 to the data points
                the function assumes the individual standard deviations of
		the data points as 1.0 (this is normal distribution)

  author:      	Christian Kelling
		Technical University of Berlin

  date:		May 1993

**************************************************************************/

double  *fit_poly_norm (double **x, double *y, int ndata, int npol)
{
  double chisq, *sig, *w, **u, **v, *a;
  double *dvector(), **dmatrix();
  int  i;
  void svdfit();

  a = dvector(1,npol);
  w = dvector(1,npol);
  sig = dvector(1,ndata);
  u = dmatrix(1,ndata,1,npol);
  v = dmatrix(1,npol,1,npol);
  for (i=1;i<=ndata;i++) sig[i] = 1.0;
  svdfit(x, y, sig, ndata, a, npol, u, v, w, &chisq, fpoly);
  return a;
}


/*************************************************************************

  function:	fit_multi_lin_norm

  parameters:	x -> pointer to pointer to an array[ndata][nvariates]
                     of x-data vectors
                y -> array[ndata] of y-data points
		ndata -> number of data points
		nap -> dimension of the multidimensional fitting function

  returns:      pointer to the array[nap] for the fit parameters a

  purpose:	fits a multidimensional function to the data vectors
                the function assumes the individual standard deviations of
		the data vectors as 1.0 (this is normal distribution)

  author:      	Christian Kelling
		Technical University of Berlin

  date:		May 1993

**************************************************************************/

double *fit_multi_lin_norm (double **x, double *y, int ndata, int nap)
{
  double chisq, *sig, *w, **u, **v, *a;
  double *dvector(), **dmatrix();
  int  i;
  void svdfit();

  a = dvector(1,nap);
  w = dvector(1,nap);
  sig = dvector(1,ndata);
  u = dmatrix(1,ndata,1,nap);
  v = dmatrix(1,ndata,1,nap);
  for (i=1;i<=ndata;i++) sig[i] = 1.0;
  svdfit(x, y, sig, ndata, a, nap, u, v, w, &chisq, fmultdim);
  printf("chisquare: %f\n",chisq);
  return a;
}




