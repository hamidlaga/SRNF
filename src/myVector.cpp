#include "myVector.h"
#include "mex.h"

// Subroutines ============================================================
void gradInverseQ(double *gradC, const double *Bu, const double *Bv, 
                  const double *fu, const double *fv, 
                  const double *nf, const double *w, 
                  const double *r1, const double *sinPhi, double dA, int N, int nB) {
    
    double *wdnf, *r5, nb[3], intgrand = 0.0, wdnb=0.0, ndnb=0.0;
    
    wdnf = new double[N];
    r5 = new double[N];
    
    #pragma omp parallel for
    for (int j=0; j<N; j++) {
        wdnf[j] = w[0*N+j]*nf[0*N+j] + w[1*N+j]*nf[1*N+j] + w[2*N+j]*nf[2*N+j];
        r5[j] = pow(r1[j],5);
        if (r5[j] > MYLARGE)
            r5[j] = 0.0;
    }
    
    #pragma omp parallel for
    for (int i=0; i<nB; i++) {   
        intgrand = 0.0;
        for (int j=0; j<N; j++) {
			//             nb[0] = fu[1*N+j]*Bv[i*d*N+2*N+j]-fu[2*N+j]*Bv[i*d*N+1*N+j]+Bu[i*d*N+1*N+j]*fv[2*N+j]-Bu[i*d*N+2*N+j]*fv[1*N+j];
			//             nb[1] = -fu[0*N+j]*Bv[i*d*N+2*N+j]+fu[2*N+j]*Bv[i*d*N+0*N+j]-Bu[i*d*N+0*N+j]*fv[2*N+j]+Bu[i*d*N+2*N+j]*fv[0*N+j];
			//             nb[2] = fu[0*N+j]*Bv[i*d*N+1*N+j]-fu[1*N+j]*Bv[i*d*N+0*N+j]+Bu[i*d*N+0*N+j]*fv[1*N+j]-Bu[i*d*N+1*N+j]*fv[0*N+j]; 
			
            crossLong(nb,fu+0*N+j,fv+0*N+j,Bu+i*3*N+0*N+j,Bv+i*3*N+0*N+j,N);
            
            wdnb = dot3Long(w+j,nb,N);
            ndnb = dot3Long(nf+j,nb,N);
			//             wdnb = w[0*N+j]*nb[0] + w[1*N+j]*nb[1] + w[2*N+j]*nb[2];
			//             ndnb = nf[0*N+j]*nb[0] + nf[1*N+j]*nb[1] + nf[2*N+j]*nb[2];
			
            intgrand += (2*wdnb*r1[j] - ndnb*wdnf[j]*r5[j])*sinPhi[j];
        }
        gradC[i] = intgrand*dA;
    }
    delete [] wdnf;
    delete [] r5;
}

void gradInverseQnew(double *gradC, const double *Bu, const double *Bv, 
                     const double *fu, const double *fv, 
                     const double *nf, const double *w, const double *r1, 
                     const double *sinPhi, double dA, int N, int nB) {
    double *wdnf, *r5, nb[3], intgrand=0.0, wdnb=0.0, ndnb=0.0;
    
    wdnf = new double[N];
    r5   = new double[N];
    
    //  mexPrintf("Here .....\n");
    #pragma omp parallel for
    for (int j=0; j<N; j++) {
        wdnf[j] = dot3(w+j*3,nf+j*3);
		//         wdnf[j] = w[0*N+j]*nf[0*N+j] + w[1*N+j]*nf[1*N+j] + w[2*N+j]*nf[2*N+j];
        r5[j] = pow(r1[j], 5); //5);
        // mexPrintf("%f\n", r5[j]);
        if (r5[j] > MYLARGE)
            r5[j] = 0.0;
    }
    
    
    #pragma omp parallel for
    for (int i=0; i<nB; i++) {   
        intgrand = 0.0;
        for (int j=0; j<N; j++) {
			
            crossDiff(nb,fu+j*3, fv+j*3, Bu+i*3*N+j*3, Bv+i*3*N+j*3);
            
            wdnb = dot3(w+j*3,nb);
            ndnb = dot3(nf+j*3,nb);
			
            intgrand += (2*wdnb*r1[j] - ndnb*wdnf[j]*r5[j])*sinPhi[j];
        }
        gradC[i] = intgrand*dA;
    }
    
    delete [] wdnf;
    delete [] r5;
}

void surf2q(double *q, double *nf, double *r, double *r1, const double *fu, const double *fv, int N) {
    
    #pragma omp parallel for
    for (int j=0; j<N; j++) {
        nf[0*N+j] = fu[1*N+j]*fv[2*N+j] - fu[2*N+j]*fv[1*N+j];
        nf[1*N+j] = - fu[0*N+j]*fv[2*N+j] + fu[2*N+j]*fv[0*N+j];
        nf[2*N+j] = fu[0*N+j]*fv[1*N+j] - fu[1*N+j]*fv[0*N+j];
        
        r[j] = pow(nf[0*N+j]*nf[0*N+j] + nf[1*N+j]*nf[1*N+j] + nf[2*N+j]*nf[2*N+j],.25);
		
        r1[j] = 1/r[j];
        if (r1[j] > MYLARGE)
            r1[j] = 0.0;
		
        
        q[0*N+j] = nf[0*N+j]*r1[j];
        q[1*N+j] = nf[1*N+j]*r1[j];
        q[2*N+j] = nf[2*N+j]*r1[j];
        
    }
}

void surf2qnew(double *q, double *nf, double *r, double *r1, const double *fu, const double *fv, int N) {
    
    #pragma omp parallel for
    for (int j=0; j<N; j++) {
        cross(nf+j*3+0,fu+j*3+0,fv+j*3+0);
//         printf("%f %f %f\n", fu[3*j], fu[3*j+1],fu[3*j+2]);
//         printf("%f %f %f\n", fv[3*j], fv[3*j+1],fv[3*j+2]);
//         getchar();
        r[j] = dot3(nf+j*3+0,nf+j*3+0);
       // printf("%f\n", r[j]);
        r[j] = pow(*(r+j),.25);
		
        r1[j] = 1/r[j];
		if (r1[j] > MYLARGE)
            r1[j] = 0.0;
		 
        q[j*3+0] = nf[j*3+0]*r1[j];
        q[j*3+1] = nf[j*3+1]*r1[j];
        q[j*3+2] = nf[j*3+2]*r1[j];
    }
}

void innerR3(double *w, const double *u, const double *v, int N)
{
    #pragma omp parallel for
    for (int j=0; j<N; j++) {
        w[j] = dot3(u+j*3+0,v+j*3+0);
    }
}


void innerS2(double *s, const double *u, const double *v, const double *sinPhi, int N)
{
    s[0] = 0.0;
    #pragma omp parallel for
    for (int j=0; j<N; j++) {
        s[0] += dot3(u+j*3, v+j*3)*sinPhi[j];
    }
}

// void pointProd(double *w, const double *u, const double *v, int n)
// {
//     for (int i = 0; i < n; i++) {
// 		w[i] = dot(u+i,v+i);
// 	}
//     
// }

double dotProd(const double *u, const double *v, int d)
{
	double innp = 0.0;
	
	for (int i = 0; i < d; i++) {
		innp += u[i]*v[i];
	}
	
	return innp;
}



// void GramSchmitd(double *x, int &n, int d)
// {
// 	// n=# of basis; d= dimension for each basis element
// 	double innp;
//     int cnt = 0, k;
// 	    
//     cnt = n; //initial # of vectors
//     
//     // 1st vector
//     k = 0;
//     innp = sqrt(InProd(x,x,d));
//     for (int i=0; i<d; i++) {
//         x[d*k + i] /= innp;
//     }
//     
//     // 2nd to the last
//     k = 1;
// 	do{
//         // kth vector
//         for (int j=0; j<k; j++) {
//             innp = InProd(x+d*k, x+d*j,d);
// 
//             for (int i=0; i<d; i++) {
//                 x[d*k+i] -= innp*x[d*j + i];
//             }
//         }
// 
// 		innp = sqrt(InProd(x+d*k,x+d*k,d));
// 		if (innp>0.0000000001) {
// 			for (int i=0; i<d; i++) {
// 				x[d*k + i] /= innp;
// 			}
//             k += 1;
// 		}
//         else {
//             for (int i=0; i<d; i++) {
// 				x[d*k + i] = x[d*cnt + i];
// 			}
//             cnt -= 1;
//         }
// 	} while (k < cnt);
//     
//     n = cnt;
// }
// 
// double innerSquare(const double *u, const double *v, int n1, int n2, int d)
// {
//     int N = n1*n2;
// 	double innp = 0.0, du, dv;
// 
//     du = 1.0/(n1-1);
// 	dv = 1.0/(n2-1);
//     
// 	for (int i = 0; i < N*d; i++) {
// 		innp += u[i]*v[i];
// 	}
// 
//     innp *= du*dv;
// 	return innp;
// }
// 


// void GramSchmitdSquare(double *x, int &n, int n1, int n2, int d)
// {
// 	// n=# of basis; d= dimension for each basis element
// 	double innp;
//     int cnt = 0, k, N = n1*n2*d, i, j;
//     
//     cnt = n-1;
//     
//     // 1st vector
//     k = 0;
//     innp = sqrt(innerSquare(x,x,n1,n2,d));
//     for (i=0; i<N; i++) {
//         x[N*k + i] /= innp;
//     }
//     
//     // 2nd to the last
//     k = 1;
// 	do{
//         // kth vector
//         for (j=0; j<k; j++) {
//             innp = innerSquare(x+N*k, x+N*j,n1,n2,d);
// 
//             for (i=0; i<N; i++) {
//                 x[N*k+i] -= innp*x[N*j + i];
//             }
//         }
// 
// 		innp = sqrt(innerSquare(x+N*k,x+N*k,n1,n2,d));
// 		if (innp>EPS) {
// 			for (i=0; i<N; i++) {
// 				x[N*k + i] /= innp;
// 			}
//             k += 1;
// 		}
//         else {
//             for (i=0; i<N; i++) {
// 				x[N*k + i] = x[N*cnt + i];
// 			}
//             cnt -= 1;
//         }
// 	} while (k <= cnt);
//     
//     n = k;
// }
// 
// void GramSchmitdSphere(double *x, int &n, double *phi, int n1, int n2, int d)
// {
// 	// n=# of basis; d= dimension for each basis element
// 	double innp;
//     int cnt = 0, k, N = n1*n2*d, i, j;
//     
//     cnt = n-1;
//     
//     // 1st vector
//     k = 0;
//     innp = sqrt(innerSphere(x,x,phi,n1,n2,d));
//     for (i=0; i<N; i++) {
//         x[N*k + i] /= innp;
//     }
//     
//     // 2nd to the last
//     k = 1;
// 	do{
//         // kth vector
//         for (j=0; j<k; j++) {
//             innp = innerSphere(x+N*k, x+N*j,phi,n1,n2,d);
// 
//             for (i=0; i<N; i++) {
//                 x[N*k+i] -= innp*x[N*j + i];
//             }
//         }
// 
// 		innp = sqrt(innerSphere(x+N*k,x+N*k,phi,n1,n2,d));
// 		if (innp>EPS) {
// 			for (i=0; i<N; i++) {
// 				x[N*k + i] /= innp;
// 			}
//             k += 1;
// 		}
//         else {
//             for (i=0; i<N; i++) {
// 				x[N*k + i] = x[N*cnt + i];
// 			}
//             cnt -= 1;
//         }
// 	} while (k <= cnt);
//     
//     n = k;
// }