/* disloc_mex.c -- MEX interface to disloc.c

   Record of revisions:

   Date          Programmer            Description of Change
   ====          ==========            =====================
   10/28/2000    Peter Cervelli        Original Code

*/

#include <mex.h>

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
     double *pOutput;
     double *pModel, *pCoords, *pOptions;
     double nu;
     int i, m, n, NumStat, NumDisl;

     /*Check argument syntax*/

     if (nrhs != 3 || nlhs > 1)
          {
               mexPrintf("disloc 1.2    10/28/2000\nUsage: u=disloc(model,stations,nu)\n");
               return;
          }

     /*Check model vector*/

          m=mxGetM(prhs[0]);
          NumDisl=mxGetN(prhs[0]);

          if (m != 10)
               mexErrMsgTxt("First argument must be a 10xn matrix containing n dislocations, stored columnwise.");

          pModel = mxGetPr(prhs[0]);

     /*Check station coordinate matrix*/

          m=mxGetM(prhs[1]);
          n=mxGetN(prhs[1]);

          if (m != 2)
               mexErrMsgTxt("Second argument must be a 2xn matrix containing n station coordinates, stored columnwise.");

          NumStat=n;
          pCoords = mxGetPr(prhs[1]);

     /*Check Poisson's ratio*/

          m=mxGetM(prhs[2]);
          n=mxGetN(prhs[2]);

          if (m != 1 || n != 1)
               mexErrMsgTxt("Third argument must be a scalar (Poisson's ratio).");

          nu =  mxGetScalar(prhs[2]);

     /*Create output array and call main function*/

          plhs[0] = mxCreateDoubleMatrix(3, NumStat, mxREAL);
          pOutput = mxGetPr(plhs[0]);

          Disloc(pOutput, pModel, pCoords, nu, NumStat, NumDisl);
}

