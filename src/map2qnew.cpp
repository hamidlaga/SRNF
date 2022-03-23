#include "mex.h"
#include "matrix.h"
#include "myTemplate.h"
#include "myVector.h"
#include "surfaceChristoffel.h"

// main function ==========================================================
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    double *q, *nf, *r, *r1;
	int N;
	const double *fu, *fv;
    const int *Fd;
	
	if (nrhs != 3)
		mexErrMsgTxt("usage: [q,nf,r] = map2qnew(fu,fv,d,N);");
	if (!mxIsDouble(prhs[0]))
		mexErrMsgTxt("Double precision arguments expected.");

    if (nlhs != 4)
		mexErrMsgTxt("Expected four(4) returns.");
 
    fu = mxGetPr(prhs[0]);
    fv = mxGetPr(prhs[1]);
    N = (int)mxGetScalar(prhs[2]);
    
    
    Fd = mxGetDimensions(prhs[0]);
    if ( (Fd[0]*Fd[1]) != (3*N) ) {
        mexPrintf("1st dim: %d.\n", Fd[0]*Fd[1]); //
        mexPrintf("2nd dim: %d.\n", 3*N); //
        mexErrMsgTxt("Dimension of fu mismatches with N.");
    }
    
	plhs[0] = mxCreateDoubleMatrix(3,N,mxREAL);
	q = mxGetPr(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix(3,N,mxREAL);
	nf = mxGetPr(plhs[1]);
    plhs[2] = mxCreateDoubleMatrix(1,N,mxREAL);
	r = mxGetPr(plhs[2]);
    plhs[3] = mxCreateDoubleMatrix(1,N,mxREAL);
	r1 = mxGetPr(plhs[3]);

    // mexPrintf("1st dim: %d.\n", Fd[0]); //d = 3
    // mexPrintf("2nd dim: %d.\n", Fd[1]); //N = prod(res);
    // mexErrMsgTxt("TEST.");
    
	surf2qnew(q,nf,r,r1,fu,fv,N);

}
