//You can include any C libraries that you normally use
//#define CHAR16_T
#include "mex.h"   //--This one is required
#include <math.h>

#include <cstdio>
#include <cstdlib>
#include "image.h"
#include "misc.h"
#include "pnmfile.h"
#include "segment-image.h"
#include "2dmem_template_class.h"

/* [IJCV04] Efficient Graph-Based Image Segmentation
 *
 * Usage:
 * [segment nSeg] = EGBIS(ucRGBImg, sigma, k, min_size);
 *
 * Input:
 * ucRGBImg - [uint8] | unsigned char | input image
 *		[height*width*3]
 * sigma -  double | float | to smooth the image.
 * k - double | float | constant for threshold function.
 * min_size - int32 | int | minimum component size (enforced by post-processing stage).
 *
 *
 * Output:
 * segment - [int32] | TC2DMem<int> | segment
 *		[height*width]
 * nSeg - int32 | int | num. of the segment
 *
 */

void CheckInput (int nrhs, const mxArray *prhs[]){
/*
 * Input:
 * ucRGBImg - [uint8] | unsigned char | input image
 *		[height*width*3]
 * sigma -  double | float | to smooth the image.
 * k - double | float | constant for threshold function.
 * min_size - int32 | int | minimum component size (enforced by post-processing stage).
 */
    if(nrhs<4)
        mexErrMsgTxt("EGBIS requires 4 Inputs.");

    if(!mxIsClass(prhs[0], "uint8"))		//check ucRGBImg
        mexErrMsgTxt("in[0] ucRGBImg must be 'uint8'.");

    if(!mxIsClass(prhs[1], "double"))		//check sigma
        mexErrMsgTxt("in[1] sigma must be 'double'.");

    if(!mxIsClass(prhs[2], "double"))		//check k
        mexErrMsgTxt("in[2] k must be 'double'.");

    if(!mxIsClass(prhs[3], "int32"))		//check min_size
        mexErrMsgTxt("in[3] min_size must be 'int32'.");
}

void TransferInputData (int nrhs, const mxArray *prhs[],
						image<rgb>** ppUCRGBImg,
						float& fSigma,
						float& fK,
						int& nMinSize){
/*
 * Input:
 * ucRGBImg - [uint8] | unsigned char | input image
 *		[height*width*3]
 * sigma -  double | float | to smooth the image.
 * k - double | float | constant for threshold function.
 * min_size - int32 | int | minimum component size (enforced by post-processing stage).
 */
	int nHeight, nWidth;

	int nPatch;
	int nYLabel, nHLabel, nFeature;
	int i, j;

	const mwSize* pDims;

	unsigned char* pUC;
	double* pDouble;

	//check input
	CheckInput(nrhs, prhs);

	//[0] ucRGBImg
	pDims = mxGetDimensions(prhs[0]);
	nHeight = (int)(pDims[0]);
	nWidth = (int)(pDims[1]);

	(*ppUCRGBImg) = new image<rgb>(nWidth, nHeight);
	pUC = (unsigned char*)mxGetData(prhs[0]);
	//r
	for(j=0; j<nWidth; j++){	//nWidth
		for(i=0; i<nHeight; i++){	//nHeight
			(*ppUCRGBImg)->access[i][j].r = (*pUC);
			pUC++;
		}
	}
	//g
	for(j=0; j<nWidth; j++){	//nWidth
		for(i=0; i<nHeight; i++){	//nHeight
			(*ppUCRGBImg)->access[i][j].g = (*pUC);
			pUC++;
		}
	}
	//b
	for(j=0; j<nWidth; j++){	//nWidth
		for(i=0; i<nHeight; i++){	//nHeight
			(*ppUCRGBImg)->access[i][j].b = (*pUC);
			pUC++;
		}
	}

	//[1] sigma
	fSigma = (float)(((double*)mxGetData(prhs[1]))[0]);
	
	//[2] k
	fK = (float)(((double*)mxGetData(prhs[2]))[0]);

	//[3] min_size
	nMinSize = ((int*)mxGetData(prhs[3]))[0];
}

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	//Input:
	image<rgb>* pUCRGBImg;
	float fSigma;
	float fK;
	int nMinSize;

	//Output:
	TC2DMem<int> segment;
	int nSeg;

	//temp
	image<rgb>* pSeg;
	int i, j;
	int* pInt;

	//transfer input data
	TransferInputData(nrhs, prhs, &pUCRGBImg, fSigma, fK, nMinSize);

	pSeg = segment_image(pUCRGBImg, fSigma, fK, nMinSize, &nSeg, segment);

	delete pUCRGBImg;
	delete pSeg;

	//output
	//[0] segment
	plhs[0] = mxCreateNumericMatrix(segment.m_row, segment.m_col, mxINT32_CLASS, mxREAL);
	pInt = (int*)mxGetData(plhs[0]);

	for(j=0; j<segment.m_col; j++){
		for(i=0; i<segment.m_row; i++){
			(*pInt) = segment[i][j];
			pInt++;
		}
	}

	//[1] nSeg
	plhs[1] = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL);
	*((int*)mxGetData(plhs[1])) = nSeg;
}
