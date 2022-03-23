#ifndef CLOSEDIAN_H
#define CLOSEDIAN_H

#include <cmath>

extern const double DodecaElements[60][3][3];

void ApplyTran(double *qnew, const double *q, const double *T, int n, int t);
void area_surf_closed(double *Snorm, double *A, double *A_tmp1, double *A_tmp2, const double *F, int n, int t);
void surface_to_q(double *q, const double *f, const double *multfact, int n, int t);
void findgrad_closed(double *dfdu, double *dfdv, const double *f, const double *Theta, int n, int t);
void updategam(double *gamnew, const double *gamupdate, const double *gamid, double eps, int n, int t);
void center_surface_closed(double *Fnew, const double *F, const double *multfact, double A, const double *Theta, int n, int t);
void scale_surf(double *Fnew, const double *F, double A, int n, int t);
void findphistarclosed(double *w, const double *q2, const double *Psi, const double *b, const double *Theta, int n, int t, int a);
void Calculate_Distance_Closed(double *H, const double *q1, const double *q2, const double *Theta, int n, int t);
void findupdategamclosed(double *gamupdate, const double *v, const double *w, const double *b, const double *Theta, int n, int t, int a);
void optimal_rot_surf(double *Ot, const double *q1, const double *q2, const double *Theta, int n, int t);
void Apply_gam_gamid_closed(double *gamcum, const double *gamid, const double *gaminc, int n);
void Apply_gam_gamid_closed(double *gamcum, const double *gamid, const double *gaminc, const double *Theta, const double *Phi, int n);
void Apply_Gamma_Surf_Closed(double *Fnew, const double *F, const double *Theta, const double *Phi, const double *gam, int n);
int ReParamclosed(double *Fnew, double *H, double *gampnew, const double *f1, const double *f2new, const double *Snorm1, const double *Snorm2, const double *A_tmp21, 
				  const double *A_tmp22new, const double *b, const double *gamp, const double *gamid, const double *Theta, 
				  const double *Phi, const double *Psi, int itermax, double eps, int n, int t, int a);
void findoptimalparamet(double *f1new, double *f2new, double *A2new, double *A_tmp12new, double *A_tmp22new, double *A_tmp21, 
						double *optdodeca, double *optrot, const double *f1, const double *f2, const double *Theta, const double *Phi,
						int scale, int n, int t);
void spherical_to_cart(double *x, double *y, double *z, double theta, double phi, double r=1.0);
void cartesian_to_sph(double *theta, double *phi, double *r, double x, double y, double z);

inline void spherical_to_cart(double *x, double *y, double *z, double theta, double phi, double r) {
	double rsint = r*sin(theta);
	*x = rsint*sin(phi);
	*y = rsint*cos(phi);
	*z = r*cos(theta);
}

inline void cartesian_to_sph(double *theta, double *phi, double *r, double x, double y, double z) {
	*r = sqrt(x*x+y*y+z*z);
	*theta = acos(z/(*r));
	*phi = atan2(x,y);
}

#endif // !CLOSED_H
