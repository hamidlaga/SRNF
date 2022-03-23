#include <cmath>
#include <cstring>
#include "Surface.h"
#include "ClosedIan.h"

void ApplyTran(double *qnew, const double *q, const double *T, int n, int t) {
	int i, k, N = n*t;

	for (k = 0; k < 3; ++k) {
		for (i = 0; i < N; ++i) {
			qnew[N*k + i] = q[N*k + i] + T[k];
		}
	}
}

void area_surf_closed(double *Snorm, double *A, double *A_tmp1, double *A_tmp2, const double *F, int n, int t) {
	int i, j, k, N = n*t;
	double dfdu[3], dfdv[3], c[3], dtheta, dphi, sint;

	*A = 0;

	dphi = (2*M_PI)/(n-1);
	dtheta = M_PI*(n-0.02)/((n-1)*(n-1+0.02));

	sint = sin(M_PI*(0.01+0)/(n-1+0.02));

	// j = 0, k = 0
	fdiff3(dfdu, F + n*(t*0 + 0) + 0, dtheta, n, N);
	fdiff3(dfdv, F + n*(t*0 + 0) + 0, dphi, 1, N);
	div3(dfdv,sint);
	cross(c,dfdu,dfdv);
	for (i = 0; i < 3; ++i) {
		Snorm[n*(t*i + 0) + 0] = c[i];
	}
	A_tmp1[n*0 + 0] = L2norm3(c);
	A_tmp2[n*0 + 0] = sqrt(A_tmp1[n*0 + 0]);
	*A += sint*A_tmp1[n*0 + 0];

	// k = 0
	for (j = 1; j < n-1; ++j) {
		fdiff3(dfdu, F + n*(t*0 + 0) + j, dtheta, n, N);
		cdiff3(dfdv, F + n*(t*0 + 0) + j, dphi, 1, N);
		div3(dfdv,sint);
		cross(c,dfdu,dfdv);
		for (i = 0; i < 3; ++i) {
			Snorm[n*(t*i + 0) + j] = c[i];
		}	
		A_tmp1[n*0 + j] = L2norm3(c);
		A_tmp2[n*0 + j] = sqrt(A_tmp1[n*0 + j]);
		*A += sint*A_tmp1[n*0 + j];
	}

	// j = n-1, k = 0
	fdiff3(dfdu, F + n*(t*0 + 0) + n-1, dtheta, n, N);
	bdiff3(dfdv, F + n*(t*0 + 0) + n-1, dphi, 1, N);
	div3(dfdv,sint);
	cross(c,dfdu,dfdv);
	for (i = 0; i < 3; ++i) {
		Snorm[n*(t*i + 0) + n-1] = c[i];
	}
	A_tmp1[n*0 + n-1] = L2norm3(c);
	A_tmp2[n*0 + n-1] = sqrt(A_tmp1[n*0 + n-1]);
	*A += sint*A_tmp1[n*0 + n-1];

	for (k = 1; k < n-1; ++k) {

		sint = sin(M_PI*(0.01+k)/(n-1+0.02));

		// j = 0
		cdiff3(dfdu, F + n*(t*0 + k) + 0, dtheta, n, N);
		fdiff3(dfdv, F + n*(t*0 + k) + 0, dphi, 1, N);
		div3(dfdv,sint);
		cross(c,dfdu,dfdv);
		for (i = 0; i < 3; ++i) {
			Snorm[n*(t*i + k) + 0] = c[i];
		}	
		A_tmp1[n*k + 0] = L2norm3(c);
		A_tmp2[n*k + 0] = sqrt(A_tmp1[n*k + 0]);
		*A += sint*A_tmp1[n*k + 0];

		for (j = 1; j < n-1; ++j) {
			cdiff3(dfdu, F + n*(t*0 + k) + j, dtheta, n, N);
			cdiff3(dfdv, F + n*(t*0 + k) + j, dphi, 1, N);
			div3(dfdv,sint);
			cross(c,dfdu,dfdv);
			for (i = 0; i < 3; ++i) {
				Snorm[n*(t*i + k) + j] = c[i];
			}
			A_tmp1[n*k + j] = L2norm3(c);
			A_tmp2[n*k + j] = sqrt(A_tmp1[n*k + j]);
			*A += sint*A_tmp1[n*k + j];
		}

		// j = n-1
		cdiff3(dfdu, F + n*(t*0 + k) + n-1, dtheta, n, N);
		bdiff3(dfdv, F + n*(t*0 + k) + n-1, dphi, 1, N);
		div3(dfdv,sint);
		cross(c,dfdu,dfdv);
		for (i = 0; i < 3; ++i) {
			Snorm[n*(t*i + k) + n-1] = c[i];
		}
		A_tmp1[n*k + n-1] = L2norm3(c);
		A_tmp2[n*k + n-1] = sqrt(A_tmp1[n*k + n-1]);
		*A += sint*A_tmp1[n*k + n-1];
	}

	sint = sin(M_PI*(0.01+n-1)/(n-1+0.02));

	// j = 0, k = n-1
	bdiff3(dfdu, F + n*(t*0 + n-1) + 0, dtheta, n, N);
	fdiff3(dfdv, F + n*(t*0 + n-1) + 0, dphi, 1, N);
	div3(dfdv,sint);
	cross(c,dfdu,dfdv);
	for (i = 0; i < 3; ++i) {
		Snorm[n*(t*i + n-1) + 0] = c[i];
	}	
	A_tmp1[n*(n-1) + 0] = L2norm3(c);
	A_tmp2[n*(n-1) + 0] = sqrt(A_tmp1[n*(n-1) + 0]);
	*A += sint*A_tmp1[n*(n-1) + 0];

	// k = n-1
	for (j = 1; j < n-1; ++j) {
		bdiff3(dfdu, F + n*(t*0 + n-1) + j, dtheta, n, N);
		cdiff3(dfdv, F + n*(t*0 + n-1) + j, dphi, 1, N);
		div3(dfdv,sint);
		cross(c,dfdu,dfdv);
		for (i = 0; i < 3; ++i) {
			Snorm[n*(t*i + n-1) + j] = c[i];
		}	
		A_tmp1[n*(n-1) + j] = L2norm3(c);
		A_tmp2[n*(n-1) + j] = sqrt(A_tmp1[n*(n-1) + j]);
		*A += sint*A_tmp1[n*(n-1) + j];
	}

	// j = n-1, k = n-1
	bdiff3(dfdu, F + n*(t*0 + n-1) + n-1, dtheta, n, N);
	bdiff3(dfdv, F + n*(t*0 + n-1) + n-1, dphi, 1, N);
	div3(dfdv,sint);
	cross(c,dfdu,dfdv);
	for (i = 0; i < 3; ++i) {
		Snorm[n*(t*i + n-1) + n-1] = c[i];
	}	
	A_tmp1[n*(n-1) + n-1] = L2norm3(c);
	A_tmp2[n*(n-1) + n-1] = sqrt(A_tmp1[n*(n-1) + n-1]);
	*A += sint*A_tmp1[n*(n-1) + n-1];

	*A *= dphi*dtheta;
}

void surface_to_q(double *q, const double *f, const double *multfact, int n, int t) {
	int i, k, N = n*t;

	for (k = 0; k < 3; ++k) {
		for (i = 0; i < N; ++i) {
				q[N*k + i] = f[N*k + i]/multfact[i];
		}
	}
}

void findgrad_closed(double *dfdu, double *dfdv, const double *f, const double *Theta, int n, int t) {
	int i, j, k;
	double dtheta, dphi;

	dtheta = (M_PI*(n-0.02))/((n-1)*(n-1+0.02));
	dphi = (2*M_PI)/(n-1);

	for (i = 0; i < 3; ++i) {
		// j = 0, k = 0
		dfdu[n*(t*i + 0) + 0] = fdiff(f + n*(t*i + 0) + 0, dtheta, n);
		dfdv[n*(t*i + 0) + 0] = fdiff(f + n*(t*i + 0) + 0, dphi, 1);
		dfdv[n*(t*i + 0) + 0] /= sin(Theta[n*0 + 0]);

		// k = 0
		for (j = 1; j < n-1; ++j) {
			dfdu[n*(t*i + 0) + j] = fdiff(f + n*(t*i + 0) + j, dtheta, n);
			dfdv[n*(t*i + 0) + j] = cdiff(f + n*(t*i + 0) + j, dphi, 1);
			dfdv[n*(t*i + 0) + j] /= sin(Theta[n*0 + j]);
		}

		// j = n-1, k = 0
		dfdu[n*(t*i + 0) + n-1] = fdiff(f + n*(t*i + 0) + n-1, dtheta, n);
		dfdv[n*(t*i + 0) + n-1] = bdiff(f + n*(t*i + 0) + n-1, dphi, 1);
		dfdv[n*(t*i + 0) + n-1] /= sin(Theta[n*0 + n-1]);

		for (k = 1; k < n-1; ++k) {
			// j = 0
			dfdu[n*(t*i + k) + 0] = cdiff(f + n*(t*i + k) + 0, dtheta, n);
			dfdv[n*(t*i + k) + 0] = fdiff(f + n*(t*i + k) + 0, dphi, 1);
			dfdv[n*(t*i + k) + 0] /= sin(Theta[n*k + 0]);

			for (j = 1; j < n-1; ++j) {
				dfdu[n*(t*i + k) + j] = cdiff(f + n*(t*i + k) + j, dtheta, n);
				dfdv[n*(t*i + k) + j] = cdiff(f + n*(t*i + k) + j, dphi, 1);
				dfdv[n*(t*i + k) + j] /= sin(Theta[n*k + j]);
			}

			// j = n-1
			dfdu[n*(t*i + k) + n-1] = cdiff(f + n*(t*i + k) + n-1, dtheta, n);
			dfdv[n*(t*i + k) + n-1] = bdiff(f + n*(t*i + k) + n-1, dphi, 1);
			dfdv[n*(t*i + k) + n-1] /= sin(Theta[n*k + n-1]);
		}

		// j = 0, k = n-1
		dfdu[n*(t*i + n-1) + 0] = bdiff(f + n*(t*i + n-1) + 0, dtheta, n);
		dfdv[n*(t*i + n-1) + 0] = fdiff(f + n*(t*i + n-1) + 0, dphi, 1);
		dfdv[n*(t*i + n-1) + 0] /= sin(Theta[n*(n-1) + 0]);

		// k = n-1
		for (j = 1; j < n-1; ++j) {
			dfdu[n*(t*i + n-1) + j] = bdiff(f + n*(t*i + n-1) + j, dtheta, n);
			dfdv[n*(t*i + n-1) + j] = cdiff(f + n*(t*i + n-1) + j, dphi, 1);
			dfdv[n*(t*i + n-1) + j] /= sin(Theta[n*(n-1) + j]);
		}

		// j = n-1, k = n-1
		dfdu[n*(t*i + n-1) + n-1] = bdiff(f + n*(t*i + n-1) + n-1, dtheta, n);
		dfdv[n*(t*i + n-1) + n-1] = bdiff(f + n*(t*i + n-1) + n-1, dphi, 1);
		dfdv[n*(t*i + n-1) + n-1] /= sin(Theta[n*(n-1) + n-1]);
	}
}

void updategam(double *gamnew, const double *gamupdate, const double *gamid, double eps, int n, int t) {
	int i, N = n*t;
	double th1 = M_PI*0.01/(n-1+0.02), thend = M_PI*(0.01+n-1)/(n-1+0.02);

	for (i = 0; i < N; ++i) {
			gamnew[N*0 + i] = gamid[N*0 + i] + eps*gamupdate[N*0 + i];

			if (gamnew[N*0 + i] < th1)
				gamnew[N*0 + i]  = th1 + (th1 - gamnew[N*0 + i]);
			else if (gamnew[N*0 + i] > thend)
				gamnew[N*0 + i] = thend - (gamnew[N*0+ i] - thend);
	}

	for (i = 0; i < N; ++i) {
			gamnew[N*1 + i] = gamid[N*1 + i] + eps*gamupdate[N*1 + i];
	}		

}

void center_surface_closed(double *Fnew, const double *F, const double *multfact, double A, const double *Theta, int n, int t) {
	double dphi, dtheta, sums[3] = { 0, 0, 0 };
	int i, k, N = n*t;

	dphi = 2*M_PI/(n-1);
	// XXX: Off by one?
	dtheta = M_PI*(n+1 - 0.02)/(n*(n+0.02));

	for (k = 0; k < 3; ++k) {
		for (i = 0; i < N; ++i) {
			sums[k] += F[N*k + i]*multfact[i]*sin(Theta[i]);
		}
	}

	mult3(sums,dtheta*dphi/A);

	for (k = 0; k < 3; ++k) {
		for (i = 0; i < N; ++i) {
			Fnew[N*k + i] = F[N*k + i] - sums[k];
		}
	}

}

void scale_surf(double *Fnew, const double *F, double A, int n, int t) {
	int i, N = n*t*3;
	double tmp = sqrt(A);

	for (i = 0; i < N; ++i) {
		Fnew[i] = F[i]/tmp;
	}
}

void findphistarclosed(double *w, const double *q2, const double *Psi, const double *b, const double *Theta, int n, int t, int a) {
	int i, j, k, N = n*t;
	double du, dv, dbxdu, dbydv, expr11[3], expr21[3], extraterm, divb, *dq2du, *dq2dv;

	du = (M_PI*(n-0.02))/((n-1)*(n-1+0.02));
	dv = (2*M_PI)/(n-1);

	dq2du = new double[n*t*3];
	dq2dv = new double[n*t*3];

	findgrad_closed(dq2du,dq2dv,q2,Theta,n,t);

for (j = 0; j < a; ++j) {
		
		// k = 0, i = 0
		dbxdu = fdiff(b + n*(t*(2*j+0)+0)+0, du, n);
		dbydv = fdiff(b + n*(t*(2*j+1)+0)+0, dv, 1);
		dbydv /= sin(Theta[n*0+0]);
		extraterm = cos(Theta[n*0+0])*b[n*(t*(2*j+0)+0)+0]/sin(Theta[n*0+0]);
		divb = 0.5*(dbxdu + extraterm + dbydv);
		expr11[0] = divb*q2[n*(t*0+0)+0];
		expr11[1] = divb*q2[n*(t*1+0)+0];
		expr11[2] = divb*q2[n*(t*2+0)+0];
		expr21[0] = dq2du[n*(t*0+0)+0]*b[n*(t*(2*j+0)+0)+0] + dq2dv[n*(t*0+0)+0]*b[n*(t*(2*j+1)+0)+0];
		expr21[1] = dq2du[n*(t*1+0)+0]*b[n*(t*(2*j+0)+0)+0] + dq2dv[n*(t*1+0)+0]*b[n*(t*(2*j+1)+0)+0];
		expr21[2] = dq2du[n*(t*2+0)+0]*b[n*(t*(2*j+0)+0)+0] + dq2dv[n*(t*2+0)+0]*b[n*(t*(2*j+1)+0)+0];
		w[n*(t*(3*j+0)+0)+0] = expr11[0] + expr21[0];
		w[n*(t*(3*j+1)+0)+0] = expr11[1] + expr21[1];
		w[n*(t*(3*j+2)+0)+0] = expr11[2] + expr21[2];


		// k = 0
		for (i = 1; i < n-1; ++i) {
			dbxdu = fdiff(b + n*(t*(2*j+0)+0)+i, du, n);
			dbydv = cdiff(b + n*(t*(2*j+1)+0)+i, dv, 1);
			dbydv /= sin(Theta[n*0+i]);
			extraterm = cos(Theta[n*0+i])*b[n*(t*(2*j+0)+0)+i]/sin(Theta[n*0+i]);
			divb = 0.5*(dbxdu + extraterm + dbydv);
			expr11[0] = divb*q2[n*(t*0+0)+i];
			expr11[1] = divb*q2[n*(t*1+0)+i];
			expr11[2] = divb*q2[n*(t*2+0)+i];
			expr21[0] = dq2du[n*(t*0+0)+i]*b[n*(t*(2*j+0)+0)+i] + dq2dv[n*(t*0+0)+i]*b[n*(t*(2*j+1)+0)+i];
			expr21[1] = dq2du[n*(t*1+0)+i]*b[n*(t*(2*j+0)+0)+i] + dq2dv[n*(t*1+0)+i]*b[n*(t*(2*j+1)+0)+i];
			expr21[2] = dq2du[n*(t*2+0)+i]*b[n*(t*(2*j+0)+0)+i] + dq2dv[n*(t*2+0)+i]*b[n*(t*(2*j+1)+0)+i];
			w[n*(t*(3*j+0)+0)+i] = expr11[0] + expr21[0];
			w[n*(t*(3*j+1)+0)+i] = expr11[1] + expr21[1];
			w[n*(t*(3*j+2)+0)+i] = expr11[2] + expr21[2];
		}

		// k = 0, i = n-1
		dbxdu = fdiff(b + n*(t*(2*j+0)+0)+n-1, du, n);
		dbydv = bdiff(b + n*(t*(2*j+1)+0)+n-1, dv, 1);
		dbydv /= sin(Theta[n*0+n-1]);
		extraterm = cos(Theta[n*0+n-1])*b[n*(t*(2*j+0)+0)+n-1]/sin(Theta[n*0+n-1]);
		divb = 0.5*(dbxdu + extraterm + dbydv);
		expr11[0] = divb*q2[n*(t*0+0)+n-1];
		expr11[1] = divb*q2[n*(t*1+0)+n-1];
		expr11[2] = divb*q2[n*(t*2+0)+n-1];
		expr21[0] = dq2du[n*(t*0+0)+n-1]*b[n*(t*(2*j+0)+0)+n-1] + dq2dv[n*(t*0+0)+n-1]*b[n*(t*(2*j+1)+0)+n-1];
		expr21[1] = dq2du[n*(t*1+0)+n-1]*b[n*(t*(2*j+0)+0)+n-1] + dq2dv[n*(t*1+0)+n-1]*b[n*(t*(2*j+1)+0)+n-1];
		expr21[2] = dq2du[n*(t*2+0)+n-1]*b[n*(t*(2*j+0)+0)+n-1] + dq2dv[n*(t*2+0)+n-1]*b[n*(t*(2*j+1)+0)+n-1];
		w[n*(t*(3*j+0)+0)+n-1] = expr11[0] + expr21[0];
		w[n*(t*(3*j+1)+0)+n-1] = expr11[1] + expr21[1];
		w[n*(t*(3*j+2)+0)+n-1] = expr11[2] + expr21[2];

		for (k = 1; k < n-1; ++k) {
			// i = 0
			dbxdu = cdiff(b + n*(t*(2*j+0)+k)+0, du, n);
			dbydv = fdiff(b + n*(t*(2*j+1)+k)+0, dv, 1);
			dbydv /= sin(Theta[n*k+0]);
			extraterm = cos(Theta[n*k+0])*b[n*(t*(2*j+0)+k)+0]/sin(Theta[n*k+0]);
			divb = 0.5*(dbxdu + extraterm + dbydv);
			expr11[0] = divb*q2[n*(t*0+k)+0];
			expr11[1] = divb*q2[n*(t*1+k)+0];
			expr11[2] = divb*q2[n*(t*2+k)+0];
			expr21[0] = dq2du[n*(t*0+k)+0]*b[n*(t*(2*j+0)+k)+0] + dq2dv[n*(t*0+k)+0]*b[n*(t*(2*j+1)+k)+0];
			expr21[1] = dq2du[n*(t*1+k)+0]*b[n*(t*(2*j+0)+k)+0] + dq2dv[n*(t*1+k)+0]*b[n*(t*(2*j+1)+k)+0];
			expr21[2] = dq2du[n*(t*2+k)+0]*b[n*(t*(2*j+0)+k)+0] + dq2dv[n*(t*2+k)+0]*b[n*(t*(2*j+1)+k)+0];
			w[n*(t*(3*j+0)+k)+0] = expr11[0] + expr21[0];
			w[n*(t*(3*j+1)+k)+0] = expr11[1] + expr21[1];
			w[n*(t*(3*j+2)+k)+0] = expr11[2] + expr21[2];

			for (i = 1; i < n-1; ++i) {
				dbxdu = cdiff(b + n*(t*(2*j+0)+k)+i, du, n);
				dbydv = cdiff(b + n*(t*(2*j+1)+k)+i, dv, 1);
				dbydv /= sin(Theta[n*k+i]);
				extraterm = cos(Theta[n*k+i])*b[n*(t*(2*j+0)+k)+i]/sin(Theta[n*k+i]);
				divb = 0.5*(dbxdu + extraterm + dbydv);
				expr11[0] = divb*q2[n*(t*0+k)+i];
				expr11[1] = divb*q2[n*(t*1+k)+i];
				expr11[2] = divb*q2[n*(t*2+k)+i];
				expr21[0] = dq2du[n*(t*0+k)+i]*b[n*(t*(2*j+0)+k)+i] + dq2dv[n*(t*0+k)+i]*b[n*(t*(2*j+1)+k)+i];
				expr21[1] = dq2du[n*(t*1+k)+i]*b[n*(t*(2*j+0)+k)+i] + dq2dv[n*(t*1+k)+i]*b[n*(t*(2*j+1)+k)+i];
				expr21[2] = dq2du[n*(t*2+k)+i]*b[n*(t*(2*j+0)+k)+i] + dq2dv[n*(t*2+k)+i]*b[n*(t*(2*j+1)+k)+i];
				w[n*(t*(3*j+0)+k)+i] = expr11[0] + expr21[0];
				w[n*(t*(3*j+1)+k)+i] = expr11[1] + expr21[1];
				w[n*(t*(3*j+2)+k)+i] = expr11[2] + expr21[2];
			}

			// i = n-1
			dbxdu = cdiff(b + n*(t*(2*j+0)+k)+n-1, du, n);
			dbydv = bdiff(b + n*(t*(2*j+1)+k)+n-1, dv, 1);
			dbydv /= sin(Theta[n*k+n-1]);
			extraterm = cos(Theta[n*k+n-1])*b[n*(t*(2*j+0)+k)+n-1]/sin(Theta[n*k+n-1]);
			divb = 0.5*(dbxdu + extraterm + dbydv);
			expr11[0] = divb*q2[n*(t*0+k)+n-1];
			expr11[1] = divb*q2[n*(t*1+k)+n-1];
			expr11[2] = divb*q2[n*(t*2+k)+n-1];
			expr21[0] = dq2du[n*(t*0+k)+n-1]*b[n*(t*(2*j+0)+k)+n-1] + dq2dv[n*(t*0+k)+n-1]*b[n*(t*(2*j+1)+k)+n-1];
			expr21[1] = dq2du[n*(t*1+k)+n-1]*b[n*(t*(2*j+0)+k)+n-1] + dq2dv[n*(t*1+k)+n-1]*b[n*(t*(2*j+1)+k)+n-1];
			expr21[2] = dq2du[n*(t*2+k)+n-1]*b[n*(t*(2*j+0)+k)+n-1] + dq2dv[n*(t*2+k)+n-1]*b[n*(t*(2*j+1)+k)+n-1];
			w[n*(t*(3*j+0)+k)+n-1] = expr11[0] + expr21[0];
			w[n*(t*(3*j+1)+k)+n-1] = expr11[1] + expr21[1];
			w[n*(t*(3*j+2)+k)+n-1] = expr11[2] + expr21[2];
		}

		// k = n-1, i = 0
		dbxdu = bdiff(b + n*(t*(2*j+0)+n-1)+0, du, n);
		dbydv = fdiff(b + n*(t*(2*j+1)+n-1)+0, dv, 1);
		dbydv /= sin(Theta[n*(n-1)+0]);
		extraterm = cos(Theta[n*(n-1)+0])*b[n*(t*(2*j+0)+n-1)+0]/sin(Theta[n*(n-1)+0]);
		divb = 0.5*(dbxdu + extraterm + dbydv);
		expr11[0] = divb*q2[n*(t*0+n-1)+0];
		expr11[1] = divb*q2[n*(t*1+n-1)+0];
		expr11[2] = divb*q2[n*(t*2+n-1)+0];
		expr21[0] = dq2du[n*(t*0+n-1)+0]*b[n*(t*(2*j+0)+n-1)+0] + dq2dv[n*(t*0+n-1)+0]*b[n*(t*(2*j+1)+n-1)+0];
		expr21[1] = dq2du[n*(t*1+n-1)+0]*b[n*(t*(2*j+0)+n-1)+0] + dq2dv[n*(t*1+n-1)+0]*b[n*(t*(2*j+1)+n-1)+0];
		expr21[2] = dq2du[n*(t*2+n-1)+0]*b[n*(t*(2*j+0)+n-1)+0] + dq2dv[n*(t*2+n-1)+0]*b[n*(t*(2*j+1)+n-1)+0];
		w[n*(t*(3*j+0)+n-1)+0] = expr11[0] + expr21[0];
		w[n*(t*(3*j+1)+n-1)+0] = expr11[1] + expr21[1];
		w[n*(t*(3*j+2)+n-1)+0] = expr11[2] + expr21[2];

		// k = n-1
		for (i = 1; i < n-1; ++i) {
			dbxdu = bdiff(b + n*(t*(2*j+0)+n-1)+i, du, n);
			dbydv = cdiff(b + n*(t*(2*j+1)+n-1)+i, dv, 1);
			dbydv /= sin(Theta[n*(n-1)+i]);
			extraterm = cos(Theta[n*(n-1)+i])*b[n*(t*(2*j+0)+n-1)+i]/sin(Theta[n*(n-1)+i]);
			divb = 0.5*(dbxdu + extraterm + dbydv);
			expr11[0] = divb*q2[n*(t*0+n-1)+i];
			expr11[1] = divb*q2[n*(t*1+n-1)+i];
			expr11[2] = divb*q2[n*(t*2+n-1)+i];
			expr21[0] = dq2du[n*(t*0+n-1)+i]*b[n*(t*(2*j+0)+n-1)+i] + dq2dv[n*(t*0+n-1)+i]*b[n*(t*(2*j+1)+n-1)+i];
			expr21[1] = dq2du[n*(t*1+n-1)+i]*b[n*(t*(2*j+0)+n-1)+i] + dq2dv[n*(t*1+n-1)+i]*b[n*(t*(2*j+1)+n-1)+i];
			expr21[2] = dq2du[n*(t*2+n-1)+i]*b[n*(t*(2*j+0)+n-1)+i] + dq2dv[n*(t*2+n-1)+i]*b[n*(t*(2*j+1)+n-1)+i];
			w[n*(t*(3*j+0)+n-1)+i] = expr11[0] + expr21[0];
			w[n*(t*(3*j+1)+n-1)+i] = expr11[1] + expr21[1];
			w[n*(t*(3*j+2)+n-1)+i] = expr11[2] + expr21[2];
		}

		// k = n-1, i = n-1
		dbxdu = bdiff(b + n*(t*(2*j+0)+n-1)+n-1, du, n);
		dbydv = bdiff(b + n*(t*(2*j+1)+n-1)+n-1, dv, 1);
		dbydv /= sin(Theta[n*(n-1)+n-1]);
		extraterm = cos(Theta[n*(n-1)+n-1])*b[n*(t*(2*j+0)+n-1)+n-1]/sin(Theta[n*(n-1)+n-1]);
		divb = 0.5*(dbxdu + extraterm + dbydv);
		expr11[0] = divb*q2[n*(t*0+n-1)+n-1];
		expr11[1] = divb*q2[n*(t*1+n-1)+n-1];
		expr11[2] = divb*q2[n*(t*2+n-1)+n-1];
		expr21[0] = dq2du[n*(t*0+n-1)+n-1]*b[n*(t*(2*j+0)+n-1)+n-1] + dq2dv[n*(t*0+n-1)+n-1]*b[n*(t*(2*j+1)+n-1)+n-1];
		expr21[1] = dq2du[n*(t*1+n-1)+n-1]*b[n*(t*(2*j+0)+n-1)+n-1] + dq2dv[n*(t*1+n-1)+n-1]*b[n*(t*(2*j+1)+n-1)+n-1];
		expr21[2] = dq2du[n*(t*2+n-1)+n-1]*b[n*(t*(2*j+0)+n-1)+n-1] + dq2dv[n*(t*2+n-1)+n-1]*b[n*(t*(2*j+1)+n-1)+n-1];
		w[n*(t*(3*j+0)+n-1)+n-1] = expr11[0] + expr21[0];
		w[n*(t*(3*j+1)+n-1)+n-1] = expr11[1] + expr21[1];
		w[n*(t*(3*j+2)+n-1)+n-1] = expr11[2] + expr21[2];
	}


	delete [] dq2du;
	delete [] dq2dv;
}

void Calculate_Distance_Closed(double *H, const double *q1, const double *q2, const double *Theta, int n, int t) {
	int i, k, N = n*t;
	double dtheta, dphi, tmp;

	--n;
		
	dtheta = M_PI*(n+1-0.02)/(n*(n+0.02));
	dphi = 2*M_PI/n;

	*H = 0;

	for (k = 0; k < 3; ++k) {
		for (i = 0; i < N; ++i) {
			tmp = q1[N*k + i] - q2[N*k + i];
			*H += sin(Theta[i])*tmp*tmp;
		}
	}

	*H = sqrt((*H)*dtheta*dphi);
}

void findupdategamclosed(double *gamupdate, const double *v, const double *w, const double *b, const double *Theta, int n, int t, int a) {
	int i, j, k, N = n*t;
	double dtheta, dphi, innp;
		
	--n;

	dtheta = M_PI*(n+1-0.02)/(n*(n+0.02));
	dphi = 2*M_PI/n;

	memset(gamupdate,0,N*2*sizeof(double));

	for (k = 0; k < a; ++k) {
		innp = 0;
		for (j = 0; j < 3; ++j) {
			for (i = 0; i < N; ++i) {
				innp += w[N*(3*k + j) + i]*v[N*j + i]*sin(Theta[i]);
			}
		}

		innp *= dtheta*dphi;

		for (j = 0; j < 2; ++j) {
			for (i = 0; i < N; ++i) {
				gamupdate[N*j + i] += innp*b[N*(2*k + j) + i];
			}
		}
	}
}

void optimal_rot_surf(double *Ot, const double *q1, const double *q2, const double *Theta, int n, int t) {
	int i, N = n*t;
	double dphi, dtheta, sint, tmp, det, A[3][3], U[3][3], V[3][3], s[3];

	dphi = 2*M_PI/(n-1);
	// XXX: Off by one?
	dtheta = M_PI*(n+1 - 0.02)/(n*(n+0.02));

	memset(A,0,sizeof(A));

	for (i = 0; i < N; ++i) {
		sint = sin(Theta[i]);
		A[0][0] += q1[N*0 + i]*q2[N*0 + i]*sint;
		A[1][0] += q1[N*0 + i]*q2[N*1 + i]*sint;
		A[2][0] += q1[N*0 + i]*q2[N*2 + i]*sint;

		A[0][1] += q1[N*1 + i]*q2[N*0 + i]*sint;
		A[1][1] += q1[N*1 + i]*q2[N*1 + i]*sint;
		A[2][1] += q1[N*1 + i]*q2[N*2 + i]*sint;

		A[0][2] += q1[N*2 + i]*q2[N*0 + i]*sint;
		A[1][2] += q1[N*2 + i]*q2[N*1 + i]*sint;
		A[2][2] += q1[N*2 + i]*q2[N*2 + i]*sint;
	}

	tmp = dphi*dtheta;
	A[0][0] *= tmp;
	A[0][1] *= tmp;
	A[0][2] *= tmp;

	A[1][0] *= tmp;
	A[1][1] *= tmp;
	A[1][2] *= tmp;

	A[2][0] *= tmp;
	A[2][1] *= tmp;
	A[2][2] *= tmp;

	det = det3((double *)A);

	svd3((double *)U,s,(double *)V,(double *)A);

	if (det > 0) {
		Ot[3*0 + 0] = U[0][0]*V[0][0] + U[1][0]*V[1][0] + U[2][0]*V[2][0];
		Ot[3*1 + 0] = U[0][0]*V[0][1] + U[1][0]*V[1][1] + U[2][0]*V[2][1];
		Ot[3*2 + 0] = U[0][0]*V[0][2] + U[1][0]*V[1][2] + U[2][0]*V[2][2];

		Ot[3*0 + 1] = U[0][1]*V[0][0] + U[1][1]*V[1][0] + U[2][1]*V[2][0];
		Ot[3*1 + 1] = U[0][1]*V[0][1] + U[1][1]*V[1][1] + U[2][1]*V[2][1];
		Ot[3*2 + 1] = U[0][1]*V[0][2] + U[1][1]*V[1][2] + U[2][1]*V[2][2];

		Ot[3*0 + 2] = U[0][2]*V[0][0] + U[1][2]*V[1][0] + U[2][2]*V[2][0];
		Ot[3*1 + 2] = U[0][2]*V[0][1] + U[1][2]*V[1][1] + U[2][2]*V[2][1];
		Ot[3*2 + 2] = U[0][2]*V[0][2] + U[1][2]*V[1][2] + U[2][2]*V[2][2];
	}
	else {
		Ot[3*0 + 0] = U[0][0]*V[0][0] + U[1][0]*V[1][0] - U[2][0]*V[2][0];
		Ot[3*1 + 0] = U[0][0]*V[0][1] + U[1][0]*V[1][1] - U[2][0]*V[2][1];
		Ot[3*2 + 0] = U[0][0]*V[0][2] + U[1][0]*V[1][2] - U[2][0]*V[2][2];

		Ot[3*0 + 1] = U[0][1]*V[0][0] + U[1][1]*V[1][0] - U[2][1]*V[2][0];
		Ot[3*1 + 1] = U[0][1]*V[0][1] + U[1][1]*V[1][1] - U[2][1]*V[2][1];
		Ot[3*2 + 1] = U[0][1]*V[0][2] + U[1][1]*V[1][2] - U[2][1]*V[2][2];

		Ot[3*0 + 2] = U[0][2]*V[0][0] + U[1][2]*V[1][0] - U[2][2]*V[2][0];
		Ot[3*1 + 2] = U[0][2]*V[0][1] + U[1][2]*V[1][1] - U[2][2]*V[2][1];
		Ot[3*2 + 2] = U[0][2]*V[0][2] + U[1][2]*V[1][2] - U[2][2]*V[2][2];
	}
}

void Apply_gam_gamid_closed(double *gamcum, const double *gamid, const double *gaminc, int n) {
	int i, j, k;
	double ndth, ndph, t, *D, *y, th0, ph0;

	D = new double[n];
	y = new double[n];

	
	ndth = M_PI*(n-1)/(n-1+0.02);
	ndph = 2*M_PI;
	th0  = M_PI*0.01/(n-1+0.02);
	ph0  = 0;

	// THE matlab version is as follows
	/*
n=n-1;

th = pi*[.01:n+1-.01]/(n+.02);
ph = 2*pi*[0:n]'/(n);

for i=1:n+1
    gamcum(i,:,1) = spline(th,gamid(i,:,1),gaminc(i,:,1));
end

for j=1:n+1
    gamcum(:,j,2) = spline(ph,gamid(:,j,2),gaminc(:,j,2));
end
*/
	
	// (n, t, d) ==> n*(t*k + j) + i
	for (i = 0; i < n; ++i) {
		for (j = 0; j < n; ++j) {
			y[j] = gamid[n*(n*0+j)+i];
		}

		spline(D, y, n);

		for (j = 0; j < n; ++j) {
			lookupspline(t, k, gaminc[n*(n*0+j)+i]-th0, ndth, n);
			gamcum[n*(n*0+j)+i] = evalspline(t, D+k, y+k);
		}
	}

	for (j = 0; j < n; ++j) {
		spline(D, gamid+(n*(n*1+j)+0), n);

		for (i = 0; i < n; ++i) {
			lookupspline(t, k, gaminc[n*(n*1+j)+i]-ph0, ndph, n);
			gamcum[n*(n*1+j)+i] = evalspline(t, D+k, gamid+(n*(n*1+j)+k));
								
		}
	}

	delete [] D;
	delete [] y;
}

// I am using this one
void Apply_gam_gamid_closed(double *gamcum, const double *gamid, const double *gaminc, 
							const double *Theta, const double *Phi, int n) {
	int i, j, k;
	double ndth, ndph, t, *D, *y, th0, ph0;

	D = new double[n];
	y = new double[n];

	/*
	ndth = M_PI*(n-1)/(n-1+0.02);
	ndph = 2*M_PI;
	th0  = M_PI*0.01/(n-1+0.02);
	ph0  = 0;
	*/
	ndth = Theta[n*(n-1)+0] - Theta[0];
	ndph = Phi[n-1] - Phi[0];
	th0  = Theta[0];
	ph0  = Phi[0];
	
	// (n, t, d) ==> n*(t*k + j) + i
	for (i = 0; i < n; ++i) {
		for (j = 0; j < n; ++j) {
			y[j] = gamid[n*(n*0+j)+i];
		}

		spline(D, y, n);

		for (j = 0; j < n; ++j) {
			lookupspline(t, k, gaminc[n*(n*0+j)+i]-th0, ndth, n);
			gamcum[n*(n*0+j)+i] = evalspline(t, D+k, y+k);
		}
	}

	for (j = 0; j < n; ++j) {
		spline(D, gamid+(n*(n*1+j)+0), n);

		for (i = 0; i < n; ++i) {
			// y[i] = gamid[n*(n*1+j)+i];
			lookupspline(t, k, gaminc[n*(n*1+j)+i]-ph0, ndph, n);
			gamcum[n*(n*1+j)+i] = evalspline(t, D+k, gamid+(n*(n*1+j)+k));								
		}
	}

	delete [] D;
	delete [] y;
}

void Apply_Gamma_Surf_Closed(double *Fnew, const double *F, const double *Theta, const double *Phi, const double *gam, int n) {
	double *Dx, *Dy, *zx, ndth, ndph, t, u, th0, ph0;
	int i, j, k, N = n*n;

	Dy = new double[N];
	Dx = new double[n];
	zx = new double[n];

	ndth = Theta[n*(n-1)+0] - Theta[0];
	ndph = Phi[n-1] - Phi[0];
	th0  = Theta[0];
	ph0  = Phi[0];

	for (k = 0; k < 3; ++k) {
		interp2(Dy, F+(N*k+0), n, n);

		for (i = 0; i < N; ++i) {
			lookupspline(u, j, gam[N + i]-ph0, ndph, n);
			evalinterp2(u, Dx, zx, Dy+j, F+(N*k+j), n, n);

			lookupspline(t, j, gam[i]    -th0, ndth, n);
			Fnew[N*k + i] = evalspline(t, Dx+j, zx+j);
		}
	}

	delete [] Dy;
	delete [] Dx;
	delete [] zx;
}

int ReParamclosed(double *Fnew, double *H, double *gampnew, const double *f1, const double *f2new, const double *Snorm1, const double *Snorm2, const double *A_tmp21, 
				  const double *A_tmp22new, const double *b, const double *gamp, const double *gamid, const double *Theta, 
				  const double *Phi, const double *Psi, int itermax, double eps, int n, int t, int a) {
	
	double *Snormnew, *q1, *q2, *v, *w, *gamnew, *gamupdate, *gamptmp, *Ftmp, Anew, *multfactnew, *sqrtmultfactnew, Hdiff;
	int i, j, iter, N = n*t*3, Np=n*t;

	int memsize = n*t*2;  // original was n*t*3
	memcpy(gampnew, gamp, memsize*sizeof(double)); // here it assumes that gamp is of nxtx 3 (problem here !!)

	// return 0;
	double th0  = Theta[0];
	double ph0  = Phi[0];

	// double x, y, z;
	// double* S = new double[N];
	// double* Snew = new double[N];
	

	q1 = new double[n*t*3];
	q2 = new double[n*t*3];
	v  = new double[n*t*3];
	w  = new double[n*t*3*a];
	
	gamupdate = new double[n*t*2];
	gamnew    = new double[n*t*2];   // why here they are of 2 dim
	gamptmp   = new double[memsize]; // why here it is of 3
	Ftmp	  = new double[n*t*3];
	Snormnew  = new double[n*t*3];
	multfactnew = new double[n*t];
	sqrtmultfactnew = new double[n*t];

	surface_to_q(q1, Snorm1, A_tmp21, n, t);
	surface_to_q(q2, Snorm2, A_tmp22new, n, t);

	Calculate_Distance_Closed(H+0,q1,q2,Theta,n,t);

	memcpy(Fnew,f2new, N*sizeof(double));
	Hdiff = 100;

	for (iter=1; iter < itermax && Hdiff > 0.0001; ++iter) {

		findphistarclosed(w,q2,Psi,b,Theta,n,t,a);

		for (i = 0; i < N; ++i) {
			v[i] = q1[i] - q2[i];
		}

		findupdategamclosed(gamupdate,v,w,b,Theta,n,t,a);

		updategam(gamnew,gamupdate,gamid,eps,n,t);		// 2dims
		
		Apply_gam_gamid_closed(gamptmp, gampnew, gamnew, Theta, Phi, n);
		memcpy(gampnew, gamptmp, memsize*sizeof(double));

		for (i = 0; i < Np; ++i) {	
			if (gamnew[Np*1 + i] < 0)
				gamnew[Np*1 + i] = gamnew[Np*1 + i] + 2*M_PI;
			else if (gamnew[Np*1 + i] > 2*M_PI)
				gamnew[Np*1 + i] = gamnew[Np*1+ i] - 2*M_PI;	
		}

		// Apply_Gamma_Surf_Closed(gamptmp,gampnew,Theta,Phi,gamnew,n);  // CHECK THIS CHECK THIS 
		
			
		/*
		for (i = 0; i < Np; ++i) {	
			if (gamptmp[Np*1 + i] < 0)
				gamptmp[Np*1 + i] = gamptmp[Np*1 + i] + 2*M_PI;
			else if (gamptmp[Np*1 + i] > 2*M_PI)
				gamptmp[Np*1 + i] = gamptmp[Np*1+ i] - 2*M_PI;	
		}

		*/
		/*
		for (i = 0; i < Np; ++i) {	
			if (gamptmp[i] < th0)
				gamptmp[i] =  th0 + (th0 - gamptmp[i]);
			else if (gamptmp[i] > Theta[n-1])
				gamptmp[i] = Theta[n-1] - (gamptmp[i] - Theta[n-1]);	
		}

		*/
		
		

		Apply_Gamma_Surf_Closed(Ftmp,Fnew,Theta,Phi,gamnew,n);
		memcpy(Fnew,Ftmp,n*t*3*sizeof(double));

		for (j = 0; j < 3; ++j) {
			for (i = 0; i < t; ++i) {
				Fnew[n*(t*j + i) + 0] = Fnew[n*(t*j + i) + n-1];
			}
		}

		area_surf_closed(Snormnew, &Anew, multfactnew, sqrtmultfactnew, Fnew, n, t);
		surface_to_q(q2, Snormnew, sqrtmultfactnew, n, t);

		Calculate_Distance_Closed(H+iter, q1, q2, Theta, n, t);

		if (H[iter] > H[iter-1]) {
			iter = 0;
			break;
		}

		if (iter > 1)
			Hdiff = (H[iter-2]-H[iter-1])/H[iter-2];
	}
	
	delete [] q1;
	delete [] q2;
	delete [] v;
	delete [] w;
	delete [] gamupdate;
	delete [] gamnew;
	delete [] gamptmp;
	delete [] Ftmp;
	delete [] multfactnew;
	delete [] sqrtmultfactnew;
	delete [] Snormnew;

	// delete [] S;
	
	return iter - 1;
}

void findoptimalparamet(double *f1new, double *f2new, double *A2new, double *A_tmp12new, double *A_tmp22new, double *A_tmp21, 
						double *optdodeca, double *optrot, const double *f1, const double *f2, const double *Theta, const double *Phi,
						int scale, int n, int t) {
	double A1, A2, *Snorm2new, *Snorm1, *Snorm2, *A_tmp11, *A_tmp12, *A_tmp22, *q1, *q2n, *f2tmp, *f2n, *f2n2, x[3], y[3], *gamnew, tmp, Hmin = 0;
	int i, k, N = n*t, s;

	A_tmp11 = new double[n*t];
	A_tmp12 = new double[n*t];
	A_tmp22 = new double[n*t];
	Snorm1 = new double[n*t*3];
	Snorm2 = new double[n*t*3];
	Snorm2new = new double[n*t*3];
	q1 = new double[n*t*3];
	q2n = new double[n*t*3];
	f2tmp = new double[n*t*3];
	f2n = new double[n*t*3];
	f2n2 = new double[n*t*3];
	gamnew = new double[n*t*2];

	area_surf_closed(Snorm1,&A1,A_tmp11,A_tmp21,f1,n,t);

	if (scale == 1) {
		center_surface_closed(f1new,f1,A_tmp11,A1,Theta,n,t);
		scale_surf(f1new,f1new,A1,n,t);
		area_surf_closed(Snorm1,&A1,A_tmp11,A_tmp21,f1new,n,t);
	}
	else {
		memcpy(f1new,f1,n*t*3*sizeof(double));
	}

	surface_to_q(q1,Snorm1,A_tmp21,n,t);
	area_surf_closed(Snorm2,&A2,A_tmp12,A_tmp22,f2,n,t);

	if (scale == 1) {
		center_surface_closed(f2tmp,f2,A_tmp12,A2,Theta,n,t);
		scale_surf(f2tmp,f2tmp,A2,n,t);
	}
	else {
		memcpy(f2tmp,f2,n*t*3*sizeof(double));
	}

	for (k = 0; k < 60; ++k) {
		for (i = 0; i < N; ++i) {
			spherical_to_cart(x+0,x+1,x+2,Theta[i],Phi[i],1);
			matvec3(y,DodecaElements[k][0]+0,x);
			cartesian_to_sph(gamnew+(N*0+i), gamnew+(N*1+i), &tmp, y[0], y[1], y[2]);
			gamnew[N*1+i] += (gamnew[N*1+i] < 0)*(2*M_PI);

			gamnew[N*0+i] = (gamnew[N*0+i] > 0)*gamnew[N*0+i];
			s = (gamnew[N*0+i] < M_PI);
			gamnew[N*0+i] = s*gamnew[N*0+i] + (!s)*M_PI;

			gamnew[N*1+i] = (gamnew[N*1+i] > 0)*gamnew[N*1+i];
			s = (gamnew[N*1+i] < 2*M_PI);
			gamnew[N*1+i] = s*gamnew[N*1+i] + (!s)*(2*M_PI);
		}

		Apply_Gamma_Surf_Closed(f2n,f2tmp,Theta,Phi,gamnew,n);

		for (i = 0; i < N; ++i) {
			f2n2[N*0+i] = DodecaElements[k][0][0]*f2n[N*0+i] + DodecaElements[k][0][1]*f2n[N*1+i] + DodecaElements[k][0][2]*f2n[N*2+i];
			f2n2[N*1+i] = DodecaElements[k][1][0]*f2n[N*0+i] + DodecaElements[k][1][1]*f2n[N*1+i] + DodecaElements[k][1][2]*f2n[N*2+i];
			f2n2[N*2+i] = DodecaElements[k][2][0]*f2n[N*0+i] + DodecaElements[k][2][1]*f2n[N*1+i] + DodecaElements[k][2][2]*f2n[N*2+i];
		}

		area_surf_closed(Snorm2,&A2,A_tmp12,A_tmp22,f2n2,n,t);

		if (scale == 1) {	
			center_surface_closed(f2n,f2n2,A_tmp12,A2,Theta,n,t);
			scale_surf(f2n,f2n,A2,n,t);
		}
		else {
			memcpy(f2n,f2n2,n*t*3*sizeof(double));
		}

		surface_to_q(q2n,Snorm2,A_tmp22,n,t);

		Calculate_Distance_Closed(&tmp,q1,q2n,Theta,n,t);
		if (tmp < Hmin || k == 0) {
			area_surf_closed(Snorm2new,A2new,A_tmp12new,A_tmp22new,f2n,n,t);
			Hmin = tmp;
			memcpy(f2new,f2n,n*t*3*sizeof(double));
			memcpy(optdodeca,gamnew,n*t*2*sizeof(double));
			memcpy(optrot,DodecaElements[k][0]+0,9*sizeof(double));
		}

	}

	delete [] A_tmp11;
	delete [] A_tmp12;
	delete [] A_tmp22;
	delete [] Snorm1;
	delete [] Snorm2;
	delete [] q1;
	delete [] q2n;
	delete [] f2tmp;
	delete [] f2n;
	delete [] f2n2;
	delete [] gamnew;
}

