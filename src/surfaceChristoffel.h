#ifndef SURFACECHRISTOFFEL_H
#define SURFACECHRISTOFFEL_H

#include <cstring>
#include <cmath>
#include <stdio.h>
#include<stdlib.h>
#include <math.h>
#include <stddef.h>
#include <time.h>
using namespace std;

// #include "f2c.h"
// #include "clapack.h"

#include "myTemplate.h"

#define PI 3.14159265
#define EPS 0.00000000000000000000000000000000001

// prototypes
void LUdecomp(double *L, double *U, const double *A, int n);
void invMatSqrt(double *Ainv, const double *A, int n);
void Determinant(double det, const double *a,int n);
void Transpose(double *a,int n);
void CoFactor(const double *a, int n, double *b);

void invertMat(double *Ainv, const double* A, int N);
void invMat(double *Ainv, const double* A, int N);
void intfS2(double *fint, const double *f, const double *phi, int n1, int n2);
void crossBasis(double *Cij, double *B, int d, int n, int nB);
void normalSurf(double *nf, const double *Cij, const double *a, int d, int n, int nB);
void dnormalSurf(double *nf, double *nb, double *g, double *S, double *S1, const double *Cij, const double *a, const double *phi, int n1, int n2, int d, int n, int nB);
void paramSurfaceTensor(double *Gl, const double *Cij, const double *phi, double *nf, double *nb, const double *S, const double *S1, int n1, int n2, int d, int n, int nB);
void getGu(double *Gu, const double *Gl, const double *ginv, int nB);

void findgrad(double *dfdu, double *dfdv, const double *f, int n, int t);
void findgrad2D(double *dfdu, double *dfdv, const double *f, int n, int t, int d);
void multfact_image(double *multfact, const double *dfdu, const double *dfdv, int n, int t, int d);
void surface_to_q(double *q, const double *f, const double *multfact, int n, int t, int d);
void Calculate_Distance(double *H, const double *q1, const double *q2, int n, int t, int d);
void findphistar(double *w, const double *q, const double *b, int n, int t, int d, int K);
void findupdategam(double *gamupdate, const double *v, const double *w, const double *b, int n, int t, int d, int K);
void updategam(double *gamnew, const double *gamupdate, const double *gamid, double eps, int n, int t, int D);
void Apply_gam_gamid(double *gamcum, const double *gamid, const double *gaminc, int m, int n);
void Apply_Gamma_Surf(double *Fnew, const double *F, const double *gam, int n, int d);
bool check_crossing(double *C, const double *f, int n, int t, int D);
int ReparamSurf(double *Fnew, double *gampnew, double *H, 
                const double *Ft, const double *Fm, const double *gamp, 
                const double *b, const double *gamid,
                const int n, const int t, const int d, const int D, const int K, 
				const double eps, const int itermax, const double tol);

void InvtGamma(double *gaminv, const double *gam, const double *gamid, int n);
void jacob_image(double *A, const double *F, int n, int t);
// #endif // !UNITSQUARE_H

// subroutine =============================================================
// void invertMat(double *Ainv, const double* A, integer N)
// {
//     integer *IPIV;
//     integer INFO;
//     integer LWORK = N*N;
//     double *WORK = new double[LWORK];
//     
//     IPIV = new integer[N];
//     
//     for (int i=0; i<N*N; i++) {
//       Ainv[i] = A[i];
//     }
// 
//     dgetrf_(&N,&N,Ainv,&N,IPIV,&INFO);
//     dgetri_(&N,Ainv,&N,IPIV,WORK,&LWORK,&INFO);
// 
//     delete [] IPIV;
//     delete [] WORK;
// }

void invMatSqrt(double *Ainv, const double *A, int n) {
    int i, j, k;
    double sum;
    
    for (i=0; i<n; i++) {
        Ainv[i*n+i] = 1/A[i*n+i];
    }
    
    for (i=1; i<n; i++) {
        for (j=1; j<n; j++) {
            sum = 0.0;
            for (k=0; k<i+1; k++) {
                sum += A[j*n+i]*Ainv[j*n+i];
                Ainv[i*n+j] = -sum;
            }
        }
    }

}

void LUdecomp(double *L, double *U, const double *A, int n) {
    double *a;
    int k,i,j;
    
    a = new double[n*n];
    
    for (i=0; i<n*n; i++) {
        a[i] = A[i];
    }
    
    for (k=0; k<n; k++) {
        L[k*n+k] = 1;
        for (i=k+1; i<n; i++) {
            L[k*n+i] = a[k*n+i]/a[k*n+k];
            for (j=k+1; j<n; j++) {
                a[j*n+i] -= L[k*n+i]*a[j*n+k];
            }
        }
        for (j=k; j<n; j++) {
            U[j*n+k] = a[j*n+k];
        }
    }
    delete [] a;
}

/*
   Recursive definition of determinate using expansion by minors.
*/
void Determinant(double *det, const double *a,int n) {
    int i,j,j1,j2;
    double *m, tmp;

    det[0] = 0;
    
    if (n < 1) { /* Error */

    } 
    else if (n == 1) { /* Shouldn't get used */
    	det[0] = a[0];
    } 
    else if (n == 2) {
    	det[0] = a[0]*a[3] - a[1]*a[2];
    } 
    else {
        det[0] = 0;
        m = new double[(n-1)*(n-1)];
        for (j1=0 ; j1<n ; j1++) {
            for (i=1;i<n;i++) {
                j2 = 0;
                for (j=0;j<n;j++) {
                   if (j == j1) continue;
                   m[j2*(n-1)+i-1] = a[j*n+i];
                   j2++;
                }
            }
        Determinant(&tmp,m,n-1);
        det[0] += pow(-1.0,j1+2.0) * a[j1*n] * tmp;
        }
        delete [] m;
    }
}

/*
   Transpose of a square matrix, do it in place
*/
void Transpose(double *a,int n)
{
   int i,j;
   double tmp;

   for (i=1;i<n;i++) {
      for (j=0;j<i;j++) {
         tmp = a[j*n+i];
         a[j*n+i] = a[i*n+j];
         a[i*n+j] = tmp;
      }
   }
}


/*
   Find the cofactor matrix of a square matrix
*/
void CoFactor(const double *a,int n,double *b) //b is adjoing matrix of a
{
   int i,j,ii,jj,i1,j1;
   double det;
   double *c;

   c = new double[n*n];

   for (j=0;j<n;j++) {
      for (i=0;i<n;i++) {

         /* Form the adjoint a_ij */
         i1 = 0;
         for (ii=0;ii<n;ii++) {
            if (ii == i)
               continue;
            j1 = 0;
            for (jj=0;jj<n;jj++) {
               if (jj == j)
                  continue;
               c[j1*n+i1] = a[jj*n+ii];
               j1++;
            }
            i1++;
         }

         /* Calculate the determinate */
         Determinant(&det, c,n-1);

         /* Fill in the elements of the cofactor */
         b[j*n+i] = pow(-1.0,i+j+2.0) * det;
      }
   }
   delete [] c;
}


// void invMat(double *Ainv, const double *A, int N) {
//     double Adet, *Ajoint;
//     int i, j, ii;
//     
//     dgesv(a,ainv,n,idmat)
// }


void intfS2(double *fint, const double *f, const double *phi, int n1, int n2) {
    // res = [n2,n1];
    double dtheta, dphi, ds;
    int n = n1*n2;
    
    dtheta = 2*PI/(n1-1);
    dphi = PI/(n2-1);
    ds = dtheta*dphi;
    
    fint[0] = 0.0;
    for (int i=0; i<n; ++i) {
        fint[0] += f[i]*sin(phi[i])*ds;
    }
}

void normalSurf(double *nf, const double *Cij, const double *a, int d, int n, int nB) {
    int i, j, k, l, N=d*n;
    
    for (i=0; i<N; ++i) { 
        nf[i] = 0.0;
        for (l=0; l<nB*nB; ++l) {
            j = l/nB;
            k = l%nB;
            nf[i] += a[j]*a[k]*Cij[N*l+i];
        }
    }
}

void dnormalSurf(double *nf, double *nb, double *g, double *S, double *S1, const double *Cij, const double *a, const double *phi, int n1, int n2, int d, int n, int nB) {
    double ds, T1, T2;
    int i, j, k, l, N=d*n;
      
    ds = (double) PI*2*PI/((n1-1)*(n2-1));
    
    for (i=0; i<(N*nB); ++i) {
        nb[i] = 0.0;
    }
    
    for (i=0; i<N; ++i) { 
        nf[i] = 0.0;
        for (l=0; l<nB*nB; ++l) {
            j = l/nB;
            k = l%nB;
            nf[i] += a[j]*a[k]*Cij[N*l+i];
            nb[N*k+i] += a[j]*(Cij[N*(nB*j+k)+i]+Cij[N*(nB*k+j)+i]);
        }
    }
    
    for (i=0; i<n; ++i) {
        S[i] = (nf+d*i,nf+d*i,d);
        if (S[i] < EPS)
            S[i] = EPS;
        S1[i] = pow(S[i],0.5);
    }
    
    for (l=0; l<nB*nB; ++l) {
        k = l%nB; // row index
        j = l/nB; // column index
        if (k<=j) {
//             continue;
//         g[l] = 0.0;
            T1 = 0.0;
            T2 = 0.0;
            for (i=0; i<n; ++i) {
                T1 += dot3(nb+N*k+d*i, nb+N*j+d*i)/S1[i]*sin(phi[i]);
                T2 += dot3(nf+d*i, nb+N*k+d*i)*dot3(nf+d*i, nb+N*j+d*i)/(S[i]*S1[i])*sin(phi[i]);
            }
            g[l] = (T1 - 3.0/4*T2)*ds;
        }
        else
           g[k*nB+j] = g[l];
    }
}

void paramSurfaceTensor(double *Gl, const double *Cij, const double *phi, double *nf, double *nb, const double *S, const double *S1, int n1, int n2, int d, int n, int nB) {
	double *gm, *Ckm, *Chm, T1, T2, T3, T4, ds;
    int i, j, k, h, m, l, N=d*n;
    
    ds = (double) PI*2*PI/((n1-1)*(n2-1));
    
    gm = new double[nB*nB*nB];   
    Ckm = new double[d];
    Chm = new double[d];
        
   // compute gm
    for (l=0; l<nB*nB*nB; ++l) {
        m = l/(nB*nB); //3rd dim
        h = (l-m*nB*nB)/nB; // column
        k = (l-m*nB*nB)%nB; // row

        if (h>=k) {
//             continue;
//         gm[l] = 0.0;
            T1 = 0.0;
            T2 = 0.0;
            T3 = 0.0;
            T4 = 0.0;
            for (i=0; i<n; ++i) {    
                add3(Ckm, Cij+N*(nB*m+k)+d*i,Cij+N*(nB*k+m)+d*i);
                add3(Chm, Cij+N*(nB*m+h)+d*i,Cij+N*(nB*h+m)+d*i);

                T1 += (dot3(nb+N*h+d*i,Ckm) + dot3(nb+N*k+d*i,Chm))/S1[i]*sin(phi[i]);
                T2 += dot3(nb+N*k+d*i,nb+N*h+d*i)*dot3(nf+d*i, nb+N*m+d*i)/(S[i]*S1[i])*sin(phi[i]);
                T3 += (dot3(nb+N*k+d*i,nb+N*m+d*i) + dot3(nf+d*i,Ckm))*dot3(nf+d*i,nb+N*h+d*i)/(S[i]*S1[i])*sin(phi[i]);
                T3 += (dot3(nb+N*h+d*i,nb+N*m+d*i) + dot3(nf+d*i,Chm))*dot3(nf+d*i,nb+N*k+d*i)/(S[i]*S1[i])*sin(phi[i]);
                T4 += dot3(nf+d*i,nb+N*k+d*i)*dot3(nf+d*i,nb+N*h+d*i)*dot3(nf+d*i,nb+N*m+d*i)/(S[i]*S[i]*S1[i])*sin(phi[i]);
            }
            gm[l] = (T1 - T2 - 3.0/4.0*T3 + 9.0/4.0*T4)*ds;
        }
        else
            gm[m*nB*nB+k*nB+h] = gm[l];
    }
    delete [] Ckm;
    delete [] Chm;
    
    // compute Gl
    for (l=0; l<nB*nB*nB; ++l) {
        k = l/(nB*nB); //3rd dim
        j = (l-k*nB*nB)/nB; // column
        i = (l-k*nB*nB)%nB; // row
        Gl[l] = 1.0/2.0*(gm[j*nB*nB+k*nB+i] + gm[i*nB*nB+k*nB+j] - gm[l]); //[k*nB*nB+j*nB+i]
    }  
    delete [] gm;
}

void getGu(double *Gu, const double *Gl, const double *ginv, int nB) {
    double T;
    int i,j,k,l,m;
    
    for (l=0; l<nB*nB*nB; l++) {
        k = l/(nB*nB); //3rd dim
        j = (l-k*nB*nB)/nB; // column
        i = l-k*nB*nB-j*nB; // row
        T = 0.0;
        for (m=0; m<nB; m++) {
            T += ginv[m*nB+k]*Gl[m*nB*nB+j*nB+i]; // use dot(,) template?
        }
        Gu[l] = T;
    }
}

//

void findgrad(double *dfdu, double *dfdv, const double *f, int n, int t) {
	double du, dv;
    int i, j, k, N = n*t;

	du = 1.0/(t-1);
	dv = 1.0/(n-1);

    for (i = 0; i < 1; ++i) {
		// j = 0, k = 0
		dfdu[n*(t*i + 0) + 0] = fdiff(f + n*(t*i + 0) + 0, du, n);
		dfdv[n*(t*i + 0) + 0] = fdiff(f + n*(t*i + 0) + 0, dv, 1);

		// k = 0
		for (j = 1; j < n-1; ++j) {
			dfdu[n*(t*i + 0) + j] = fdiff(f + n*(t*i + 0) + j, du, n);
			dfdv[n*(t*i + 0) + j] = cdiff(f + n*(t*i + 0) + j, dv, 1);
		}

		// j = n-1, k = 0
		dfdu[n*(t*i + 0) + n-1] = fdiff(f + n*(t*i + 0) + n-1, du, n);
		dfdv[n*(t*i + 0) + n-1] = bdiff(f + n*(t*i + 0) + n-1, dv, 1);

		for (k = 1; k < t-1; ++k) {
			// j = 0
			dfdu[n*(t*i + k) + 0] = cdiff(f + n*(t*i + k) + 0, du, n);
			dfdv[n*(t*i + k) + 0] = fdiff(f + n*(t*i + k) + 0, dv, 1);

			for (j = 1; j < n-1; ++j) {
				dfdu[n*(t*i + k) + j] = cdiff(f + n*(t*i + k) + j, du, n);
				dfdv[n*(t*i + k) + j] = cdiff(f + n*(t*i + k) + j, dv, 1);
			}

			// j = n-1
			dfdu[n*(t*i + k) + n-1] = cdiff(f + n*(t*i + k) + n-1, du, n);
			dfdv[n*(t*i + k) + n-1] = bdiff(f + n*(t*i + k) + n-1, dv, 1);
		}

		// j = 0, k = t-1
		dfdu[n*(t*i + t-1) + 0] = bdiff(f + n*(t*i + t-1) + 0, du, n);
		dfdv[n*(t*i + t-1) + 0] = fdiff(f + n*(t*i + t-1) + 0, dv, 1);

		// k = t-1
		for (j = 1; j < n-1; ++j) {
			dfdu[n*(t*i + t-1) + j] = bdiff(f + n*(t*i + t-1) + j, du, n);
			dfdv[n*(t*i + t-1) + j] = cdiff(f + n*(t*i + t-1) + j, dv, 1);
		}

		// j = n-1, k = t-1
		dfdu[n*(t*i + t-1) + n-1] = bdiff(f + n*(t*i + t-1) + n-1, du, n);
		dfdv[n*(t*i + t-1) + n-1] = bdiff(f + n*(t*i + t-1) + n-1, dv, 1);
	}

}

// ------------------------------------------------------------------------
void findgrad2D(double *dfdu, double *dfdv, const double *f, int n, int t, int d) {
	double du, dv;
    int i, j, k, N = n*t;

	du = 1.0/(t-1);
	dv = 1.0/(n-1);
    
    for (i = 0; i < d; ++i) {
		// j = 0, k = 0
		dfdu[n*(t*i + 0) + 0] = fdiff(f + n*(t*i + 0) + 0, du, n);
		dfdv[n*(t*i + 0) + 0] = fdiff(f + n*(t*i + 0) + 0, dv, 1);

		// k = 0
		for (j = 1; j < n-1; ++j) {
			dfdu[n*(t*i + 0) + j] = fdiff(f + n*(t*i + 0) + j, du, n);
			dfdv[n*(t*i + 0) + j] = cdiff(f + n*(t*i + 0) + j, dv, 1);
		}

		// j = n-1, k = 0
		dfdu[n*(t*i + 0) + n-1] = fdiff(f + n*(t*i + 0) + n-1, du, n);
		dfdv[n*(t*i + 0) + n-1] = bdiff(f + n*(t*i + 0) + n-1, dv, 1);

		for (k = 1; k < t-1; ++k) {
			// j = 0
			dfdu[n*(t*i + k) + 0] = cdiff(f + n*(t*i + k) + 0, du, n);
			dfdv[n*(t*i + k) + 0] = fdiff(f + n*(t*i + k) + 0, dv, 1);

			for (j = 1; j < n-1; ++j) {
				dfdu[n*(t*i + k) + j] = cdiff(f + n*(t*i + k) + j, du, n);
				dfdv[n*(t*i + k) + j] = cdiff(f + n*(t*i + k) + j, dv, 1);
			}

			// j = n-1
			dfdu[n*(t*i + k) + n-1] = cdiff(f + n*(t*i + k) + n-1, du, n);
			dfdv[n*(t*i + k) + n-1] = bdiff(f + n*(t*i + k) + n-1, dv, 1);
		}

		// j = 0, k = t-1
		dfdu[n*(t*i + t-1) + 0] = bdiff(f + n*(t*i + t-1) + 0, du, n);
		dfdv[n*(t*i + t-1) + 0] = fdiff(f + n*(t*i + t-1) + 0, dv, 1);

		// k = t-1
		for (j = 1; j < n-1; ++j) {
			dfdu[n*(t*i + t-1) + j] = bdiff(f + n*(t*i + t-1) + j, du, n);
			dfdv[n*(t*i + t-1) + j] = cdiff(f + n*(t*i + t-1) + j, dv, 1);
		}

		// j = n-1, k = t-1
		dfdu[n*(t*i + t-1) + n-1] = bdiff(f + n*(t*i + t-1) + n-1, du, n);
		dfdv[n*(t*i + t-1) + n-1] = bdiff(f + n*(t*i + t-1) + n-1, dv, 1);
	}

}

// ------------------------------------------------------------------------
void multfact_image(double *multfact, const double *dfdu, const double *dfdv, int n, int t,int d) {
	int N = n*t;
    
//     multfact = new double[N];
    
    if (d < 3) {
        for (int i = 0; i < N; ++i) {
            multfact[i] = abs(dfdu[N*0+i]*dfdv[N*1+i] - dfdu[N*1+i]*dfdv[N*0+i]);
        }
    }
    else if (d < 4) {
        for (int i = 0; i < N; ++i) {
            multfact[i] = pow(pow(dfdu[N*0+i]*dfdv[N*1+i] - dfdu[N*1+i]*dfdv[N*0+i],2) + 
                    pow(dfdu[N*0+i]*dfdv[N*2+i] - dfdu[N*2+i]*dfdv[N*0+i],2) + 
                    pow(dfdu[N*1+i]*dfdv[N*2+i] - dfdu[N*2+i]*dfdv[N*1+i],2), 1/2);
        }
    }
    else {
        for (int i = 0; i < N; ++i) {
            multfact[i] = pow(pow(dfdu[N*0+i]*dfdv[N*1+i] - dfdu[N*1+i]*dfdv[N*0+i],2) + 
                    pow(dfdu[N*0+i]*dfdv[N*2+i] - dfdu[N*2+i]*dfdv[N*0+i],2) + 
                    pow(dfdu[N*0+i]*dfdv[N*3+i] - dfdu[N*3+i]*dfdv[N*0+i],2) + 
                    pow(dfdu[N*1+i]*dfdv[N*2+i] - dfdu[N*2+i]*dfdv[N*1+i],2) + 
                    pow(dfdu[N*1+i]*dfdv[N*3+i] - dfdu[N*3+i]*dfdv[N*1+i],2) + 
                    pow(dfdu[N*2+i]*dfdv[N*3+i] - dfdu[N*3+i]*dfdv[N*2+i],2), 1/2);
        }
    }
}

// ------------------------------------------------------------------------
void surface_to_q(double *q, const double *f, const double *multfact, int n, int t, int d) {
	int N = n*t;

	for (int k = 0; k < d; ++k) {
		for (int i = 0; i < N; ++i) {
			q[N*k + i] = sqrt(multfact[i])*f[N*k + i];
		}
	}
}

// ------------------------------------------------------------------------
void Calculate_Distance(double *H, const double *q1, const double *q2, int n, int t, int d) {
	int N = n*t*d;
	double tmp, du, dv;

	du = 1.0/(n-1);
	dv = 1.0/(t-1);

	*H = 0;

	for (int i = 0; i < N; ++i) {
		tmp = q1[i] - q2[i];
		*H += tmp*tmp;
	}

	*H = sqrt((*H)*du*dv);
}

// ------------------------------------------------------------------------
void findphistar(double *w, const double *q, const double *b, int n, int t, int d, int K) {
	int D = 2;
	double du, dv, dbxdu, dbydv, divb, *dqdu, *dqdv, *expr1, *expr2;

    du = 1.0/(t-1);
	dv = 1.0/(n-1);
    
    expr1 = new double[d];
    expr2 = new double[d];
    
	dqdu = new double[n*t*d];
	dqdv = new double[n*t*d];

	findgrad2D(dqdu, dqdv, q, n, t, d);
    
    memset(w,0,n*t*d*K*sizeof(double));

	for (int j = 0; j < K; ++j) {
		// k = 0, i = 0
		dbxdu = fdiff(b + n*(t*(D*j + 0) + 0) + 0, du, n);
		dbydv = fdiff(b + n*(t*(D*j + 1) + 0) + 0, dv, 1);
		divb = 0.5*(dbxdu + dbydv);
        
        for (int ii = 0; ii < d; ++ii) {
            expr1[ii] = divb*q[n*(t*ii + 0) + 0];
            expr2[ii] = dqdu[n*(t*ii + 0) + 0]*b[n*(t*(D*j + 0) + 0) + 0] + dqdv[n*(t*ii + 0) + 0]*b[n*(t*(D*j + 1) + 0) + 0];
            w[n*(t*(d*j + ii) + 0) + 0] = expr1[ii] + expr2[ii];
        }

		// k = 0
		for (int i = 1; i < n-1; ++i) {
			dbxdu = fdiff(b + n*(t*(D*j + 0) + 0) + i, du, n);
			dbydv = cdiff(b + n*(t*(D*j + 1) + 0) + i, dv, 1);
			divb = 0.5*(dbxdu + dbydv);
            
            for (int ii = 0; ii < d; ++ii) {
                expr1[ii] = divb*q[n*(t*ii + 0) + i];
                expr2[ii] = dqdu[n*(t*ii + 0) + i]*b[n*(t*(D*j + 0) + 0) + i] + dqdv[n*(t*ii + 0) + i]*b[n*(t*(D*j + 1) + 0) + i];
                w[n*(t*(d*j + ii) + 0) + i] = expr1[ii] + expr2[ii];
            }
		}

		// i = n-1, k = 0
		dbxdu = fdiff(b + n*(t*(D*j + 0) + 0) + n-1, du, n);
		dbydv = bdiff(b + n*(t*(D*j + 1) + 0) + n-1, dv, 1);
		divb = 0.5*(dbxdu + dbydv);
        
        for (int ii = 0; ii < d; ++ii) {
            expr1[ii] = divb*q[n*(t*ii + 0) + n-1];
            expr2[ii] = dqdu[n*(t*ii + 0) + n-1]*b[n*(t*(D*j + 0) + 0) + n-1] + dqdv[n*(t*ii + 0) + n-1]*b[n*(t*(D*j + 1) + 0) + n-1];
            w[n*(t*(d*j + ii) + 0) + n-1] = expr1[ii] + expr2[ii];
        }

		for (int k = 1; k < t-1; ++k) {
			// i = 0
			dbxdu = cdiff(b + n*(t*(D*j + 0) + k) + 0, du, n);
			dbydv = fdiff(b + n*(t*(D*j + 1) + k) + 0, dv, 1);
			divb = 0.5*(dbxdu + dbydv);
            
            for (int ii = 0; ii < d; ++ii) {
                expr1[ii] = divb*q[n*(t*ii + k) + 0];
                expr2[ii] = dqdu[n*(t*ii + k) + 0]*b[n*(t*(D*j + 0) + k) + 0] + dqdv[n*(t*ii + k) + 0]*b[n*(t*(D*j + 1) + k) + 0];
                w[n*(t*(d*j + ii) + k) + 0] = expr1[ii] + expr2[ii];
            }

			for (int i = 1; i < n-1; ++i) {
				dbxdu = cdiff(b + n*(t*(D*j + 0) + k) + i, du, n);
				dbydv = cdiff(b + n*(t*(D*j + 1) + k) + i, dv, 1);
				divb = 0.5*(dbxdu + dbydv);
                
                for (int ii = 0; ii < d; ++ii) {
                    expr1[ii] = divb*q[n*(t*ii + k) + i];
                    expr2[ii] = dqdu[n*(t*ii + k) + i]*b[n*(t*(D*j + 0) + k) + i] + dqdv[n*(t*ii + k) + i]*b[n*(t*(D*j + 1) + k) + i];
                    w[n*(t*(d*j + ii) + k) + i] = expr1[ii] + expr2[ii];
                }
			}

			// i = n-1
			dbxdu = cdiff(b + n*(t*(D*j + 0) + k) + n-1, du, n);
			dbydv = bdiff(b + n*(t*(D*j + 1) + k) + n-1, dv, 1);
			divb = 0.5*(dbxdu + dbydv);
            
            for (int ii = 0; ii < d; ++ii) {
                expr1[ii] = divb*q[n*(t*ii + k) + n-1];
                expr2[ii] = dqdu[n*(t*ii + k) + n-1]*b[n*(t*(D*j + 0) + k) + n-1] + dqdv[n*(t*ii + k) + n-1]*b[n*(t*(D*j + 1) + k) + n-1];
                w[n*(t*(d*j + ii) + k) + n-1] = expr1[ii] + expr2[ii];
            }
		}

		// i = 0, k = t-1
		dbxdu = bdiff(b + n*(t*(2*j + 0) + t-1) + 0, du, n);
		dbydv = fdiff(b + n*(t*(2*j + 1) + t-1) + 0, dv, 1);
		divb = 0.5*(dbxdu + dbydv);
        
        for (int ii = 0; ii < d; ++ii) {
            expr1[ii] = divb*q[n*(t*ii + t-1) + 0];
            expr2[ii] = dqdu[n*(t*ii + t-1) + 0]*b[n*(t*(D*j + 0) + t-1) + 0] + dqdv[n*(t*ii + t-1) + 0]*b[n*(t*(D*j + 1) + t-1) + 0];
            w[n*(t*(d*j + ii) + t-1) + 0] = expr1[ii] + expr2[ii];
        }

		// k = t-1
		for (int i = 1; i < n-1; ++i) {
			dbxdu = bdiff(b + n*(t*(D*j + 0) + t-1) + i, du, n);
			dbydv = cdiff(b + n*(t*(D*j + 1) + t-1) + i, dv, 1);
			divb = 0.5*(dbxdu + dbydv);
            
            for (int ii = 0; ii < d; ++ii) {
                expr1[ii] = divb*q[n*(t*ii + t-1) + i];
                expr2[ii] = dqdu[n*(t*ii + t-1) + i]*b[n*(t*(D*j + 0) + t-1) + i] + dqdv[n*(t*ii + t-1) + i]*b[n*(t*(D*j + 1) + t-1) + i];
                w[n*(t*(d*j + ii) + t-1) + i] = expr1[ii] + expr2[ii];
            }
		}

		// i = n-1, k = t-1
		dbxdu = bdiff(b + n*(t*(D*j + 0) + t-1) + n-1, du, n);
		dbydv = bdiff(b + n*(t*(D*j + 1) + t-1) + n-1, dv, 1);
		divb = 0.5*(dbxdu + dbydv);
        
        for (int ii = 0; ii < d; ++ii) {
            expr1[ii] = divb*q[n*(t*ii + t-1) + n-1];
            expr2[ii] = dqdu[n*(t*ii + t-1) + n-1]*b[n*(t*(D*j + 0) + t-1) + n-1] + dqdv[n*(t*ii + t-1) + n-1]*b[n*(t*(D*j + 1) + t-1) + n-1];
            w[n*(t*(d*j + ii) + t-1) + n-1] = expr1[ii] + expr2[ii];
        }
	}

	delete [] dqdu;
	delete [] dqdv;
    delete [] expr1;
    delete [] expr2;
}


// ------------------------------------------------------------------------
void findupdategam(double *gamupdate, const double *v, const double *w, const double *b, int n, int t, int d, int K) {
	int N = n*t, D = 2;
	double innp, du, dv;

	du = 1.0/(n-1);
	dv = 1.0/(t-1);

	memset(gamupdate,0,n*t*D*sizeof(double));

    for (int k = 0; k < K; ++k) {
        innp = 0;
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < d; ++j) {
                innp += v[N*j + i]*w[N*(d*k + j) + i];
            }
        }
        
        innp *= du*dv;

        for (int i = 0; i < N; ++i) {
            gamupdate[N*0 + i] += innp*b[D*N*k + N*0 + i];
            gamupdate[N*1 + i] += innp*b[D*N*k + N*1 + i];
        }
    }

}

// ------------------------------------------------------------------------
void updategam(double *gamnew, const double *gamupdate, const double *gamid, double eps, int n, int t, int D) {
	int N = n*t;

	for (int k = 0; k < D; ++k) {
		for (int i = 0; i < N; ++i) {
			gamnew[N*k + i] = gamid[N*k + i] + eps*gamupdate[N*k + i];
		}
	}
}

// ------------------------------------------------------------------------
void Apply_gam_gamid(double *gamcum, const double *gamid, const double *gaminc, int m, int n) {
	int k;
	double u0, v0, ndu, ndv, t, *y, *D1, *D2;

	u0 = 0;
	v0 = 0;
	ndu = 1;
	ndv = 1;

	y = new double[n];
	D1 = new double[n];
    D2 = new double[m];

	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < n; ++j) {
			y[j] = gamid[m*(n*0+j)+i]; // for each row (i-th)
		}

		spline(D1, y, n);

		for (int j = 0; j < n; ++j) {
			lookupspline(t, k, gaminc[m*(n*0+j)+i]-u0, ndu, n);
			gamcum[m*(n*0+j)+i] = evalspline(t, D1+k, y+k);
		}	
	}
    
	for (int j = 0; j < n; ++j) {
		spline(D2, gamid + m*(n*1+j)+0, m); // for each column (j-th)

		for (int i = 0; i < m; ++i) {
			lookupspline(t, k, gaminc[m*(n*1+j)+i]-v0, ndv, m);
			gamcum[m*(n*1+j)+i] = evalspline(t, D2+k, gamid + m*(n*1+j)+k);
		}
	}

	delete [] y;
	delete [] D1;
    delete [] D2;
}

//-------------------------------------------------------------------------
void Apply_Gamma_Surf(double *Fnew, const double *F, const double *gam, int m, int n, int d) {
	int j, N = m*n;
	double *Du, *Dv, *zu, u0, v0, ndu, ndv, u, v;

	Dv = new double[N];
	Du = new double[n];
	zu = new double[n];

	ndu = 1;
	ndv = 1;
	u0 = 0;
	v0 = 0;

	for (int k = 0; k < d; ++k) {
		interp2(Dv, F+(N*k+0), m, n); //col

		for (int i = 0; i < N; ++i) {
			lookupspline(v, j, gam[N*1+i]-v0, ndv, m); // col
			evalinterp2(v, Du, zu, Dv+j, F+(N*k+j), m, n); // row

			lookupspline(u, j, gam[N*0+i]-u0, ndu, n); // row
			Fnew[N*k + i] = evalspline(u, Du+j, zu+j); // row
		}
	}

	delete [] Dv;
	delete [] Du;
	delete [] zu;
}

// ------------------------------------------------------------------------
bool check_crossing(double *C, const double *f, int n, int t, int D) {
	bool is_diffeo = true;
    int N = n*t;
    double c, *dfdu, *dfdv;
    
    dfdu = new double[n*t*D];
    dfdv = new double[n*t*D];
    
    findgrad2D(dfdu,dfdv,f,n,t,D);
    
    *C = 0;
	for (int i=0; i<N; ++i) {
        c = dfdu[N*0 + i]*dfdv[N*1 + i] - dfdu[N*1 + i]*dfdv[N*0 + i];
        
        if (c < 0) {
            *C += 1;
        }
	}
    
    if (C > 0)
        is_diffeo = false;
    
    delete [] dfdu;
    delete [] dfdv;
    
    return is_diffeo;
}

//-------------------------------------------------------------------------
int ReparamSurf(double *Fnew, double *gampnew, double *H, 
                const double *Ft, const double *Fm, const double *gamp, 
                const double *b, const double *gamid,
                const int n, const int t, const int d, const int D, const int K, 
				double eps, const int itermax, const double tol)
{
// 	const double tol = 0.0004; // tolerance level for Hdiff
    bool is_diffeo;
	int iter = 0, N = n*t;
	double *qt, *qm, *v, *w, *gamupdate, *gamnew, *gampold, *Fold, Hdiff, epsTmp;
    double *Ft1, *Ft2, *Fm1, *Fm2, *multfactt, *multfactm, *C;
    
	qt = new double[n*t*d];
	qm = new double[n*t*d];
	v = new double[n*t*d];
	w = new double[n*t*d*K];
	gamupdate = new double[n*t*D];
	gamnew = new double[n*t*D];
    gampold = new double[n*t*D];
//     gampnew = new double[n*t*D];
	Fold = new double[n*t*d];
    Ft1 = new double[n*t*d];
    Ft2 = new double[n*t*d];
    Fm1 = new double[n*t*d];
    Fm2 = new double[n*t*d];
    multfactt = new double[n*t];
    multfactm = new double[n*t];
    C = new double[1];
    
    Apply_Gamma_Surf(Fold,Fm,gamp,n,t,d);
    
    // compute q functions for Ft and Fm
    findgrad2D(Ft1,Ft2,Ft,n,t,d);
    findgrad2D(Fm1,Fm2,Fnew,n,t,d);
    multfact_image(multfactt,Ft1,Ft2,n,t,d);
    multfact_image(multfactm,Fm1,Fm2,n,t,d);
	surface_to_q(qt,Ft,multfactt,n,t,d);
	surface_to_q(qm,Fold,multfactm,n,t,d);

    // compute initial energy
    Calculate_Distance(H+iter,qt,qm,n,t,d);

    Hdiff = 100;
    
    for (int i=0; i<(n*t*D); i++) {
        gampold[i] = gamp[i];
    }
//     for (int i=0; i<(n*t*d); i++) {
//         Fnew[i] = Fold[i];
//     }
    
    // main iteration
	for (iter = 1; iter < itermax && Hdiff > tol; ++iter) {
        // basis to tangent space
        findphistar(w,qm,b,n,t,d,K);

        // find v = q1-q2
		for (int i = 0; i < N; ++i) {
            for (int jj = 0; jj < d; ++jj) {
                v[N*jj+i] = qt[N*jj+i]-qm[N*jj+i];
            }
		}

        findupdategam(gamupdate,v,w,b,n,t,d,K);

        epsTmp = eps;
        do {  
            updategam(gamnew,gamupdate,gamid,epsTmp,n,t,D);

            Apply_gam_gamid(gampnew,gampold,gamnew,n,t);

            //Apply_Gamma_Surf(Fnew,Fold,gamnew,n,d); // apply incremental deformation to last time object
            Apply_Gamma_Surf(Fnew,Fm,gampnew,n,t,d); // apply cumulative deformation to original object

            findgrad2D(Fm1,Fm2,Fnew,n,t,d);
            multfact_image(multfactm,Fm1,Fm2,n,t,d);
            surface_to_q(qm,Fnew,multfactm,n,t,d);

            Calculate_Distance(H+iter,qt,qm,n,t,d);
            
            is_diffeo = check_crossing(C,gampnew,n,t,D);
            
            if (is_diffeo)
                epsTmp *= 1.2;
            else
                epsTmp *= 0.67;
            
        } while (H[iter] >= H[iter-1] && epsTmp > eps*0.001);

        // update iteration or break out
        if (H[iter] < H[iter-1]) { // && is_diffeo
            for (int i=0; i<(n*t*D); i++) {
                gampold[i] = gampnew[i];
            }
            for (int i=0; i<(n*t*d); i++) {
                Fold[i] = Fnew[i];
            }

            Hdiff = (H[iter-1]-H[iter])/H[iter-1];
        }
        else {
            for (int i=0; i<(n*t*D); i++) {
                gampnew[i] = gampold[i];
            }
            for (int i=0; i<(n*t*d); i++) {
                Fnew[i] = Fold[i];
            }
            break;
        }
    }

	delete [] qt;
	delete [] qm;
	delete [] v;
	delete [] w;
	delete [] gamupdate;
	delete [] gamnew;
	delete [] gampold;
	delete [] Fold;
    delete [] Ft1;
    delete [] Ft2;
    delete [] Fm1;
    delete [] Fm2;
    delete [] multfactt;
    delete [] multfactm;
    delete [] C;

	return iter;
}

// inverse transform ------------------------------------------------------
void InvtGamma(double *gaminv, const double *gam, const double *gamid, int n) {
	int i, j, k;
	double ndu, ndv, t, *y, *D;

	ndu = 1;
	ndv = 1;

	y = new double[n];
	D = new double[n];
    
	for (i = 0; i < n; ++i) {
		for (j = 0; j < n; ++j) {
			y[j] = gam[n*(n*0+j)+i];
		}

		spline(D, y, n);

		for (j = 0; j < n; ++j) {
			lookupspline(t, k, gamid[n*(n*0+j)+i], ndu, n);
			gaminv[n*(n*0+j)+i] = evalspline(t, D+k, y+k);
		}	
	}

	for (j = 0; j < n; ++j) {
		spline(D, gam + n*(n*1+j)+0, n);

		for (i = 0; i < n; ++i) {
			lookupspline(t, k, gamid[n*(n*1+j)+i], ndv, n);
			gaminv[n*(n*1+j)+i] = evalspline(t, D+k, gam + n*(n*1+j)+k);
		}
	}

	delete [] y;
	delete [] D; 
}

// ------------------------------------------------------------------------
void jacob_image(double *A, const double *F, int n, int t) {
	int j, k, N = n*t;
	double dfdu[2], dfdv[2], c, du, dv;

	du = 1.0/(n-1);
	dv = 1.0/(t-1);

	// j = 0, k = 0
	fdiff2(dfdu, F + n*(t*0 + 0) + 0, du, n, N);
	fdiff2(dfdv, F + n*(t*0 + 0) + 0, dv, 1, N);
	jacob(c,dfdu,dfdv);
	A[n*0 + 0] = sqrt(abs(c));

	// k = 0
	for (j = 1; j < n-1; ++j) {
		fdiff2(dfdu, F + n*(t*0 + 0) + j, du, n, N);
		cdiff2(dfdv, F + n*(t*0 + 0) + j, dv, 1, N);
		jacob(c,dfdu,dfdv);
		A[n*0 + j] = sqrt(abs(c));
	}

	// j = n-1, k = 0
	fdiff2(dfdu, F + n*(t*0 + 0) + n-1, du, n, N);
	bdiff2(dfdv, F + n*(t*0 + 0) + n-1, dv, 1, N);
	jacob(c,dfdu,dfdv);
	A[n*0 + n-1] = sqrt(abs(c));

	for (k = 1; k < t-1; ++k) {
		// j = 0
		cdiff2(dfdu, F + n*(t*0 + k) + 0, du, n, N);
		fdiff2(dfdv, F + n*(t*0 + k) + 0, dv, 1, N);
		jacob(c,dfdu,dfdv);
		A[n*k + 0] = sqrt(abs(c));

		for (j = 1; j < n-1; ++j) {
			cdiff2(dfdu, F + n*(t*0 + k) + j, du, n, N);
			cdiff2(dfdv, F + n*(t*0 + k) + j, dv, 1, N);
			jacob(c,dfdu,dfdv);
			A[n*k + j] = sqrt(abs(c));
		}

		// j = n-1
		cdiff2(dfdu, F + n*(t*0 + k) + n-1, du, n, N);
		bdiff2(dfdv, F + n*(t*0 + k) + n-1, dv, 1, N);
		jacob(c,dfdu,dfdv);
		A[n*k + n-1] = sqrt(abs(c));
	}

	// j = 0, k = t-1
	bdiff2(dfdu, F + n*(t*0 + t-1) + 0, du, n, N);
	fdiff2(dfdv, F + n*(t*0 + t-1) + 0, dv, 1, N);
	jacob(c,dfdu,dfdv);
	A[n*(t-1) + 0] = sqrt(abs(c));

	// k = t-1
	for (j = 1; j < n-1; ++j) {
		bdiff2(dfdu, F + n*(t*0 + t-1) + j, du, n, N);
		cdiff2(dfdv, F + n*(t*0 + t-1) + j, dv, 1, N);
		jacob(c,dfdu,dfdv);
		A[n*(t-1) + j] = sqrt(abs(c));
	}

	// j = n-1, k = t-1
	bdiff2(dfdu, F + n*(t*0 + t-1) + n-1, du, n, N);
	bdiff2(dfdv, F + n*(t*0 + t-1) + n-1, dv, 1, N);
	jacob(c,dfdu,dfdv);
	A[n*(t-1) + n-1] = sqrt(abs(c));

}

#endif