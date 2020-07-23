#include "kernels.h"

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





double Bernoulli_Kprod(double *x, double *y, double *k,int n,double p) {
  int i;
  for(i=0;i<n;i++){
   p = p*_KernelBernoulli(x[i],y[i],k[i]);
    }
  return(p);
}

void add(PyObject* a, int dim_a)
{ 

    int *array = NULL;
    if ( PyList_Check( int_list ) )
    {
        int nInts = PyList_Size( int_list );
        array = malloc( nInts * sizeof( *array ) );
        for ( int ii = 0; ii < nInts; ii++ )
        {
            PyObject *oo = PyList_GetItem( int_list, ii );
            if ( PyInt_Check( oo ) )
            {
                array[ ii ] = ( int ) PyInt_AsLong( oo );
            }
        }
    }

    a[0] = a[0] + a[0]; 
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
