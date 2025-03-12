#include <stdlib.h>
#include <stdio.h>
#include <math.h>

int main(){
  double A[10];
  double B[10];
  for (int i=0;i<10;i++){
    A[i]=i*1.0;
    B[i]=i*2.0;
  }
  double yval,ypval;
  spline_linear_val ( 10, A, B,11.4, &yval, &ypval );

  printf("%g\n",yval);
}









/******************************************************************************/

void spline_linear_val ( int ndata, double tdata[], double ydata[], 
  double tval, double *yval, double *ypval )

/******************************************************************************/
/*
  Purpose:

    SPLINE_LINEAR_VAL evaluates a piecewise linear spline at a point.

  Discussion:

    Because of the extremely simple form of the linear spline,
    the raw data points ( TDATA(I), YDATA(I)) can be used directly to
    evaluate the spline at any point.  No processing of the data
    is required.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    24 February 2004

  Author:

    John Burkardt

  Parameters:

    Input, int NDATA, the number of data points defining the spline.

    Input, double TDATA[NDATA], YDATA[NDATA], the values of the independent
    and dependent variables at the data points.  The values of TDATA should
    be distinct and increasing.

    Input, double TVAL, the point at which the spline is to be evaluated.

    Output, double *YVAL, *YPVAL, the value of the spline and its first
    derivative dYdT at TVAL.  YPVAL is not reliable if TVAL is exactly
    equal to TDATA(I) for some I.
*/
{
  int left;
  int right;
/*
  Find the interval [ TDATA(LEFT), TDATA(RIGHT) ] that contains, or is
  nearest to, TVAL.
*/
  r8vec_bracket ( ndata, tdata, tval, &left, &right );
/*
  Now evaluate the piecewise linear function.
*/
  *ypval = ( ydata[right-1] - ydata[left-1] ) 
         / ( tdata[right-1] - tdata[left-1] );

  *yval = ydata[left-1] +  ( tval - tdata[left-1] ) * (*ypval);

  return;
}

/******************************************************************************/

void r8vec_bracket ( int n, double x[], double xval, int *left,
  int *right )

/******************************************************************************/
/*
  Purpose:

    R8VEC_BRACKET searches a sorted array for successive brackets of a value.

  Discussion:

    An R8VEC is a vector of R8's.

    If the values in the vector are thought of as defining intervals
    on the real line, then this routine searches for the interval
    nearest to or containing the given value.

    It is always true that RIGHT = LEFT+1.

    If XVAL < X[0], then LEFT = 1, RIGHT = 2, and
      XVAL   < X[0] < X[1];
    If X(1) <= XVAL < X[N-1], then
      X[LEFT-1] <= XVAL < X[RIGHT-1];
    If X[N-1] <= XVAL, then LEFT = N-1, RIGHT = N, and
      X[LEFT-1] <= X[RIGHT-1] <= XVAL.

    For consistency, this routine computes indices RIGHT and LEFT
    that are 1-based, although it would be more natural in C and
    C++ to use 0-based values.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    31 May 2009

  Author:

    John Burkardt

  Parameters:

    Input, int N, length of input array.

    Input, double X[N], an array that has been sorted into ascending order.

    Input, double XVAL, a value to be bracketed.

    Output, int *LEFT, *RIGHT, the results of the search.
*/
{
  int i;

  for ( i = 2; i <= n - 1; i++ )
  {
    if ( xval < x[i-1] )
    {
      *left = i - 1;
      *right = i;
      return;
    }

   }

  *left = n - 1;
  *right = n;

  return;
}