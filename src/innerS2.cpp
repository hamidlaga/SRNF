#include <cstring>
#include <cmath>
#include "mex.h"
#include "matrix.h"
#include "myTemplate.h"
#include "myVector.h"

// main function ==========================================================
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	double *fint;
    const double *u, *v, *sinPhi;
    const int *Fd;
	int N;
	
	if (nrhs != 4)
		mexErrMsgTxt("usage: fint = innerS2(u,v,sinPhi,N)");

	if (!mxIsDouble(prhs[0]))
		mexErrMsgTxt("Double precision arguments expected.");

//     dims = mxGetNumberOfDimensions(prhs[0]);
    Fd = mxGetDimensions(prhs[0]);

    u = mxGetPr(prhs[0]);
    v = mxGetPr(prhs[1]);
    sinPhi = mxGetPr(prhs[2]);
    N = (int)mxGetScalar(prhs[3]);
    
	if (nlhs != 1)
		mexErrMsgTxt("Expected one returns.");

	plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
	fint = mxGetPr(plhs[0]);
    
	innerS2(fint,u,v,sinPhi,N);

}
