/* #include <stdlib.h> */
/* #include <math.h> */
/* #include "linpack_d.h" */
/* #include <stdio.h> */

/*
// double _BernoulliSlow(double x, int degre) {
//   double y;
//   if      (degre == 0) y = 1.;
//   else if (degre == 1) y = x - 0.5;
//   else if (degre == 2) y = x*x - x + 1./6.;
//   else if (degre == 3) y = x*x*x - 1.5 * x*x + 0.5 * x;
//   else if (degre == 4) y = x*x*x*x - 2. * x*x*x + x*x - 1./30.;
//   else if (degre == 5) y = x*x*x*x*x - 2.5 * x*x*x*x + 5./3. * x*x*x - 1./6. * x;
//   else if (degre == 6) y = x*x*x*x*x*x - 3 * x*x*x*x*x + 2.5 * x*x*x*x - 0.5 * x*x + 1./42.;
//   else y = 0.;
//   return(y);
// }
*/

double _Bernoulli(double x, int degre) {
  double y;
  if      (degre == 0) y = 1.;
  else if (degre == 1) y = x - 0.5;
  else if (degre == 2) y = x * x - x + 1./6.;
  else if (degre == 3) y = x * (x * (x - 1.5) + 0.5);
  else if (degre == 4) y = x * (x * (x * (x - 2.0) + 1.0)) - 1./30.;
  else if (degre == 5) y = x * (x * (x * (x * (x - 2.5) + 5./3.)) - 1./6.);
  else if (degre == 6) y = x * (x * (x * (x * (x * (x - 3.0) + 2.5)) - 0.5)) + 1./42.;
  else y = 0.;
  return(y);
}

double _KernelBernoulli(double x, double y, int degre) {
  double v, s, tgamma(double), pow(double, double), floor(double);
  int d;
  for(s=1.,d=1;d<=degre;d++) s = s + _Bernoulli(x, d) * _Bernoulli(y, d) / pow(tgamma(d+1.),2.);
  v = x - y; v = v - floor(v);
  if (degre%2) s = s + _Bernoulli(v, 2 * degre) / tgamma(2. * degre + 1.);
  else         s = s - _Bernoulli(v, 2 * degre) / tgamma(2. * degre + 1.);
  return(s);
}

void KernelBernoulli(double *X1, double *X2, int n1, int n2, int d, int degre, double *K) {
  int i,j,k;
  double s, *px1, *px2;
  for(i=0;i<n1;i++) {
    px2 = X2;
    for(j=0;j<n2;j++) {
      px1 = X1 + i * d;
      for(s=1.,k=0; k<d; k++) s = s * _KernelBernoulli(*px1++,*px2++,degre);
      K[i+j*n1] = s;
    }
  }
}

void KernelBernoulliAni(double *X1, double *X2, int n1, int n2, int d, int *degre, double *K) {
  int i,j,k;
  double s, *px1, *px2;
  for(i=0;i<n1;i++) {
    px2 = X2;
    for(j=0;j<n2;j++) {
      px1 = X1 + i * d;
      for(s=1.,k=0; k<d; k++) s = s * _KernelBernoulli(*px1++,*px2++,degre[k]);
      K[i+j*n1] = s;
    }
  }
}


/* Noyau passe-band */
double _KernelPB(double x, double y, double omega) {
  double s, d, sin(double), fabs(double);
  d = x - y;
  double pi = 3.14159265358979323846;
  if (fabs(d) > 1.e-16)
    s = sin(omega * d) / (d * pi);
  else
    s = omega / pi;
  return(s);
}

void KernelPB(double *X1, double *X2, int n1, int n2, int d, double omega, double *K) {
  int i,j,k;
  double s, *px1, *px2, dx, sin(double), fabs(double), sqrt(double);
  double pi = 3.14159265358979323846;
  for(i=0;i<n1;i++) {
    px2 = X2;
    for(j=0;j<n2;j++) {
      px1 = X1 + i * d;
      for(s=0.,k=0; k<d; k++) {
	dx = *px1++ - *px2++;
	s = s + dx * dx;
      }
      s = sqrt(s) / d;
      if (fabs(s) > 1.e-16)
	K[i+j*n1] = sin(omega * s) / (s * pi);
      else
	K[i+j*n1] = omega / pi;
    }
  }
}

void KernelPolynomial(double *X1, double *X2, int n1, int n2, int d, int degre, double *K) {
  int i,j,k;
  double s, *px1, *px2, pow(double, double);
  s=0;
  for(i=0;i<n1;i++) {
    px2 = X2;
    for(j=0;j<n2;j++) {
      px1 = X1 + i * d;
      for(s=1.,k=0; k<d; k++) s = s + *px1++ * *px2++;
      K[i+j*n1] = pow(s, (double) degre);
    }
  }
}

/*
// double _KernelSpline(double x, double y, int degre) {
//   double z, u, v;
//   if (x<y)		{u=x; v=y;}
//   else			{u=y; v=x;}
//   if      (degre == 0) 	z = 1.;
//   else if (degre == 1) 	z = u;
//   else if (degre == 2) 	z = u * u * (v / 2. - u / 6.);
//   else 			z = 0.;
//   return(z);
// }
*/

double _KernelSpline(double x, double y, int degre) {
  double z, u, v;
  if (x<y)		{u=x; v=y;}
  else			{u=y; v=x;}
  if      (degre == 0) 	z = 1.;
  else if (degre == 1) 	z = 1 + u;
  else if (degre == 2) 	z = 1 + u * v + u * u * (v - u / 3. ) / 2.;
  else 			z = 0.;
  return(z);
}


void KernelSpline(double *X1, double *X2, int n1, int n2, int d, int degre, double *K) {
  int i,j,k;
  double s, *px1, *px2;
  for(i=0;i<n1;i++) {
    px2 = X2;
    for(j=0;j<n2;j++) {
      px1 = X1 + i * d;
      for(s=1.,k=0; k<d; k++) s = s * _KernelSpline(*px1++,*px2++,degre);
      K[i+j*n1] = s;
    }
  }
}

