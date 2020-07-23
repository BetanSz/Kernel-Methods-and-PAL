%module kernels

%{
#define SWIG_FILE_WITH_INIT
#include "kernels.h"
%}

double _Bernoulli(double x, int degre);
double _KernelBernoulli(double x, double y, int degre);
double _KernelSpline(double x, double y, int degre);
double _KernelPB(double x, double y, double omega);
double Bernoulli_Kprod(double *x, double *y, double *k,int n,double p);
void KernelBernoulli(double *X1, double *X2, int n1, int n2, int d, int degre, double *K);
void add(PyObject* a, int dim_a) ;
