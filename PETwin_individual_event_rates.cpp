//
//
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
//#include <string.h>
#include "mex.h"

///////////////////////////////////////////////////////////////////////////////
/// MEX entry point function
///////////////////////////////////////////////////////////////////////////////
void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    char *easic, *echan;
    double *era;
    const mwSize *NA, *NC,*NAC;
    int aa, cc, era_index;
    register int k;
    
    if(nrhs!=3) {
        mexErrMsgTxt("PETwin_individual_event_rates takes 3 arguments:\n\
                * 1  easic[NE] - int8 (char) vector of asic numbers\n\
                * 2  echan[NE] - int8 (char) vector of chanel numbers\n\
                * 3  era[NA*NC] - double vector of event rates\n");
                
                return;
    }
    
    easic  =  (char *)mxGetPr(prhs[0]);
    echan =  (char *)mxGetPr(prhs[1]);

    era  = (double *)mxGetPr(prhs[2]);
    
    NA = mxGetDimensions(prhs[0]);
    NC = mxGetDimensions(prhs[1]);
    NAC = mxGetDimensions(prhs[2]);
    
    if(NA[1]!=1) { mexErrMsgTxt("Asic array must be a vector\n");  return; }
    if(NC[1]!=1) { mexErrMsgTxt("Channel array must be a vector\n");  return; }
    if(NC[0]!=NA[0]) { mexErrMsgTxt("Same number of Asics and Channel must be\n");  return; }
    
    for(k=0; k<(NA[0]); k++) {

        aa = (int) easic[k];
        cc = (int) echan[k];
        
        era_index = aa*32 + cc;
//         printf("## k=%d:     A=%d C=%d I=%d\n",k,aa,cc,era_index);
        
        era[era_index]++;
    }
    
    
} // main
