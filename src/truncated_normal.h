/* https://people.sc.fsu.edu/~jburkardt/cpp_src/truncated_normal/truncated_normal.html
Distributed under the MIT license.
Norman Johnson, Samuel Kotz, Narayanaswamy Balakrishnan,
Continuous Univariate Distributions,
Second edition,
Wiley, 1994,
ISBN: 0471584940,
LC: QA273.6.J6. */
# include <cmath>
# include <cstdlib>
# include <cstring>
# include <iomanip>
# include <iostream>
# include <RcppArmadillo.h>
using namespace std;

double normal_01_cdf ( double x );
double normal_01_cdf_inv ( double cdf );
double normal_01_pdf ( double x );
double truncated_normal_a_cdf ( double x, double mu, double sigma, double a );
double truncated_normal_a_cdf_inv ( double cdf, double mu, double sigma, double a );
double truncated_normal_a_pdf ( double x, double mu, double sigma, double a );

double normal_01_cdf ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    NORMAL_01_CDF evaluates the Normal 01 CDF.
//
//  Licensing:
//
//    This code is distributed under the MIT license.
//
//  Modified:
//
//    10 February 1999
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    A G Adams,
//    Areas Under the Normal Curve,
//    Algorithm 39,
//    Computer j.,
//    Volume 12, pages 197-198, 1969.
//
//  Parameters:
//
//    Input, double X, the argument of the CDF.
//
//    Output, double CDF, the value of the CDF.
//
{
  double a1 = 0.398942280444;
  double a2 = 0.399903438504;
  double a3 = 5.75885480458;
  double a4 = 29.8213557808;
  double a5 = 2.62433121679;
  double a6 = 48.6959930692;
  double a7 = 5.92885724438;
  double b0 = 0.398942280385;
  double b1 = 3.8052E-08;
  double b2 = 1.00000615302;
  double b3 = 3.98064794E-04;
  double b4 = 1.98615381364;
  double b5 = 0.151679116635;
  double b6 = 5.29330324926;
  double b7 = 4.8385912808;
  double b8 = 15.1508972451;
  double b9 = 0.742380924027;
  double b10 = 30.789933034;
  double b11 = 3.99019417011;
  double cdf;
  double q;
  double y;
//
//  |X| <= 1.28.
//
  if ( fabs ( x ) <= 1.28 )
  {
    y = 0.5 * x * x;

    q = 0.5 - fabs ( x ) * ( a1 - a2 * y / ( y + a3 - a4 / ( y + a5
      + a6 / ( y + a7 ) ) ) );
//
//  1.28 < |X| <= 12.7
//
  }
  else if ( fabs ( x ) <= 12.7 )
  {
    y = 0.5 * x * x;

    q = exp ( - y ) * b0 / ( fabs ( x ) - b1
      + b2  / ( fabs ( x ) + b3
      + b4  / ( fabs ( x ) - b5
      + b6  / ( fabs ( x ) + b7
      - b8  / ( fabs ( x ) + b9
      + b10 / ( fabs ( x ) + b11 ) ) ) ) ) );
//
//  12.7 < |X|
//
  }
  else
  {
    q = 0.0;
  }
//
//  Take account of negative X.
//
  if ( x < 0.0 )
  {
    cdf = q;
  }
  else
  {
    cdf = 1.0 - q;
  }

  return cdf;
}
//****************************************************************************80


double normal_01_pdf ( double x )

//****************************************************************************80
//
//  Purpose:
//
//    NORMAL_01_PDF evaluates the Normal 01 PDF.
//
//  Discussion:
//
//    The Normal 01 PDF is also called the "Standard Normal" PDF, or
//    the Normal PDF with 0 mean and standard deviation 1.
//
//    PDF(X) = exp ( - 0.5 * X^2 ) / sqrt ( 2 * PI )
//
//  Licensing:
//
//    This code is distributed under the MIT license.
//
//  Modified:
//
//    18 September 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the PDF.
//
//    Output, double PDF, the value of the PDF.
//
{
  double pdf;
  const double r8_pi = 3.14159265358979323;

  pdf = exp ( -0.5 * x * x ) / sqrt ( 2.0 * r8_pi );

  return pdf;
}
//****************************************************************************80

// [[Rcpp::export]]
double truncated_normal_a_cdf ( double x, double mu, double sigma, double a )

//****************************************************************************80
//
//  Purpose:
//
//    TRUNCATED_NORMAL_A_CDF evaluates the lower truncated Normal CDF.
//
//  Licensing:
//
//    This code is distributed under the MIT license.
//
//  Modified:
//
//    24 January 2017
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the CDF.
//
//    Input, double MU, SIGMA, the mean and standard deviation of the
//    parent Normal distribution.
//
//    Input, double A, the lower truncation limit.
//
//    Output, double TRUNCATED_NORMAL_A_CDF, the value of the CDF.
//
{
  double alpha;
  double alpha_cdf;
  double cdf;
  double xi;
  double xi_cdf;

  if ( x < a )
  {
    cdf = 0.0;
  }
  else
  {
    alpha = ( a - mu ) / sigma;
    xi = ( x - mu ) / sigma;

    alpha_cdf = normal_01_cdf ( alpha );
    xi_cdf = normal_01_cdf ( xi );

    cdf = ( xi_cdf - alpha_cdf ) / ( 1.0 - alpha_cdf );
  }
  
  return cdf;
}

// [[Rcpp::export]]
double truncated_normal_a_pdf ( double x, double mu, double sigma, double a )

//****************************************************************************80
//
//  Purpose:
//
//    TRUNCATED_NORMAL_A_PDF evaluates the lower truncated Normal PDF.
//
//  Licensing:
//
//    This code is distributed under the MIT license.
//
//  Modified:
//
//    24 January 2017
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the PDF.
//
//    Input, double MU, SIGMA, the mean and standard deviation of the
//    parent Normal distribution.
//
//    Input, double A, the lower truncation limit.
//
//    Output, double TRUNCATED_NORMAL_A_PDF, the value of the PDF.
//
{
  double alpha;
  double alpha_cdf;
  double pdf;
  double xi;
  double xi_pdf;

  if ( x < a )
  {
    pdf = 0.0;
  }
  else
  {
    alpha = ( a - mu ) / sigma;
    xi = ( x - mu ) / sigma;

    alpha_cdf = normal_01_cdf ( alpha );
    xi_pdf = normal_01_pdf ( xi );

    pdf = xi_pdf / ( 1.0 - alpha_cdf ) / sigma;
  }
  
  return pdf;
}
//****************************************************************************80



