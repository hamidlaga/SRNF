#ifndef MYVECTOR_H
#define MYVECTOR_H

#include <iostream>
#include <cstring>
#include <cmath>
#include <math.h>
#include "myTemplate.h"
// #include "surfaceChristoffel.h"

using namespace std;

#define SQRT5	2.23606797749978969641
#define MYLARGE 10000.0
// #define EPS	0.00000000001
// # define PI 3.14159265359

// Prototypes: -----------------------------------------------------------------
void pointProd(double *w, const double *u, const double *v, int n);
double dotProd(const double *u, const double *v, int d);
void innerS2(double *s, const double *u, const double *v, const double *sinPhi, int N);
void innerR3(double *w, const double *u, const double *v, int N);

void gradInverseQ(double *gradC, const double *Bu, const double *Bv, const double *fu, const double *fv, const double *nf, const double *w, const double *r1, const double *sinPhi, double dA, int N, int nB);
void surf2qnew(double *q, double *nf, double *r, double *r1, const double *fu, const double *fv, int N);

void gradInverseQnew(double *gradC, const double *Bu, const double *Bv, const double *fu, const double *fv, const double *nf, const double *w, const double *r1, const double *sinPhi, double dA, int N, int nB);

#endif