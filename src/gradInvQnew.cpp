#include "mex.h"
#include "matrix.h"
#include "myTemplate.h"
#include "myVector.h"
#include "surfaceChristoffel.h"

// main function ==========================================================
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    double *gradC, dA;
	const double *Bu, *Bv, *fu, *fv, *nf, *w, *r1, *sinPhi;
    int N, nB;
    const int *Fd;

	if (nrhs != 11)
		mexErrMsgTxt("usage: gradC = gradInvQ(Bu,Bv,fu,fv,nf,r1,w,sinPhi,dA,N,nB)");

	if (!mxIsDouble(prhs[0]))
		mexErrMsgTxt("Double precision arguments expected.");
    
//     if (nlhs != 1)
// 		mexErrMsgTxt("Expected one returns.");

    Bu = mxGetPr(prhs[0]);
    Bv = mxGetPr(prhs[1]);
    fu = mxGetPr(prhs[2]);
    fv = mxGetPr(prhs[3]);
    nf = mxGetPr(prhs[4]);
    w = mxGetPr(prhs[5]);
    r1 = mxGetPr(prhs[6]);
    sinPhi = mxGetPr(prhs[7]);
    dA = mxGetScalar(prhs[8]);
    N = (int)mxGetScalar(prhs[9]);
    nB = (int)mxGetScalar(prhs[10]);
    
    Fd = mxGetDimensions(prhs[0]);
    
    if ( Fd[0] != 3 ) {
        mexPrintf("1st dim = 3: %d.\n", Fd[0]); //
        mexPrintf("2nd dim = N: %d.\n", Fd[1]); //
        mexErrMsgTxt("First 2 dimension of Bu mismatches.");
    }
    if (Fd[2] != nB)
        mexErrMsgTxt("3rd dimension of Bu mismatches with nB.");

    //     mexPrintf("d: %d .\n", Fd[0]); 
    //     mexPrintf("N: %d .\n", Fd[1]); 
    //     mexPrintf("N: %d .\n", Fd[2]); 
    //     mexErrMsgTxt("TEST.");
    
	plhs[0] = mxCreateDoubleMatrix(nB,1,mxREAL);
	gradC = mxGetPr(plhs[0]);

	gradInverseQnew(gradC,Bu,Bv,fu,fv,nf,w,r1,sinPhi,dA,N,nB);
}
