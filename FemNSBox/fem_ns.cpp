
// matlab header
#include <mex.h>

// 
#include <cmath>
#include <cstring>

////////////////////////////////////////////////////////////////////////////////
// The computational routine
////////////////////////////////////////////////////////////////////////////////

void arrayProduct(double x, double *y, double *z, mwSize n)
{
	mwSize i;
	/* multiply each element y by x */
	for (i = 0; i < n; i++) {
		z[i] = x * y[i];
	}
}


////////////////////////////////////////
// gauss

const int ngp = 7;

constexpr double w1 = 0.225000000000000;
constexpr double w2 = 0.132394152788506;
constexpr double w3 = 0.125939180544827;
constexpr double a1 = 0.333333333333333;
constexpr double a2 = 0.059715871789770;
constexpr double a3 = 0.797426985353087;
constexpr double b2 = (1 - a2) / 2;
constexpr double b3 = (1 - a3) / 2;

// Gauss points
static double gp[ngp][3] = {
	a1,a1,1 - a1 - a1,
	a2,b2,1 - a2 - b2,
	b2,a2,1 - b2 - a2,
	b2,b2,1 - b2 - b2,
	a3,b3,1 - a3 - b3,
	b3,a3,1 - b3 - a3,
	b3,b3,1 - b3 - b3,
};
// Gauss weights
static double w[ngp] = { w1,w2,w2,w2,w3,w3,w3 };

////////////////////////////////////////////////////////////////////////////////
// Global data
static double sh1[ngp][3]; // P1 Shape functions
static double sh1r[ngp][3]; // P1 Shape functions r-derivatives
static double sh1s[ngp][3]; // P1 Shape functions s-derivatives
static double sh2[ngp][6]; // P2 Shape functions
static double sh2r[ngp][6]; // P2 Shape functions r-derivatives
static double sh2s[ngp][6]; // P2 Shape functions s-derivatives


void init_shape() {
	for (int i = 0; i < ngp; i++) {
		// P1 Shape functions
		sh1[i][0] = gp[i][2];
		sh1[i][1] = gp[i][0];
		sh1[i][2] = gp[i][1];
		// P1 Shape functions r-derivatives
		sh1r[i][0] = -1;
		sh1r[i][1] = 1;
		sh1r[i][2] = 0;
		// P1 Shape functions s-derivatives
		sh1s[i][0] = -1;
		sh1s[i][1] = 0;
		sh1s[i][2] = 1;

		// P2 Shape functions
		sh2[i][0] = 1 - 3 * gp[i][0] - 3 * gp[i][1] + 2 * gp[i][0] * gp[i][0] +
			4 * gp[i][0] * gp[i][1] + 2 * gp[i][1] * gp[i][1];
		sh2[i][1] = -gp[i][0] + 2 * gp[i][0] * gp[i][0];
		sh2[i][2] = -gp[i][1] + 2 * gp[i][1] * gp[i][1];
		sh2[i][3] = 4 * gp[i][0] * gp[i][1];
		sh2[i][4] = 4 * gp[i][1] - 4 * gp[i][0] * gp[i][1] - 4 * gp[i][1] * gp[i][1];
		sh2[i][5] = 4 * gp[i][0] - 4 * gp[i][0] * gp[i][1] - 4 * gp[i][0] * gp[i][0];
		// P2 Shape functions r-derivatives
		sh2r[i][0] = -3 + 4 * gp[i][0] + 4 * gp[i][1];
		sh2r[i][1] = -1 + 4 * gp[i][0];
		sh2r[i][2] = 0;
		sh2r[i][3] = 4 * gp[i][1];
		sh2r[i][4] = -4 * gp[i][1];
		sh2r[i][5] = 4 - 8 * gp[i][0] - 4 * gp[i][1];
		// P2 Shape functions s-derivatives
		sh2s[i][0] = -3 + 4 * gp[i][0] + 4 * gp[i][1];
		sh2s[i][1] = 0;
		sh2s[i][2] = -1 + 4 * gp[i][1];
		sh2s[i][3] = 4 * gp[i][0];
		sh2s[i][4] = 4 - 8 * gp[i][1] - 4 * gp[i][0];
		sh2s[i][5] = -4 * gp[i][0];
	}
}

void localKL() {
	
}



////////////////////////////////////////////////////////////////////////////////
// The gateway function
////////////////////////////////////////////////////////////////////////////////

void mexFunction(
	int nlhs, mxArray *plhs[],
	int nrhs, const mxArray *prhs[])
{
	if (nlhs <= 0) return;

	init_shape();

	{
		mwSize dims[] = { 1 };
		const char *fieldnames[] = { "ngp" };
		const size_t hoge = _countof(fieldnames);
		mxArray *s = mxCreateStructArray(1, dims, _countof(fieldnames), fieldnames);

		mxSetField(s, 0, "ngp", mxCreateDoubleScalar(ngp));

		{
			mxArray *a = mxCreateDoubleMatrix(ngp, 3, mxREAL);
			double *p = mxGetPr(a);
			for (int j = 0; j < 3; j++) {
				for (int i = 0; i < ngp; i++) {
					*p = gp[i][j];
					++p;
				}
			}

			mxAddField(s, "gp");
			mxSetField(s, 0, "gp", a);
		}

		{
			mxArray *a = mxCreateDoubleMatrix(ngp, 3, mxREAL);
			double *p = mxGetPr(a);
			for (int j = 0; j < 3; j++) {
				for (int i = 0; i < ngp; i++) {
					*p = sh1[i][j];
					++p;
				}
			}

			int fieldidx = mxAddField(s, "shape");
			mxSetField(s, 0, "shape", a);
		}


		plhs[0] = s;

	}

	return;

	double multiplier;              /* input scalar */
	double *inMatrix;               /* 1xN input matrix */
	size_t ncols;                   /* size of matrix */
	double *outMatrix;              /* output matrix */

									/* check for proper number of arguments */
	if (nrhs != 2) {
		mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs", "Two inputs required.");
	}
	if (nlhs != 1) {
		mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs", "One output required.");
	}
	/* make sure the first input argument is scalar */
	if (!mxIsDouble(prhs[0]) ||
		mxIsComplex(prhs[0]) ||
		mxGetNumberOfElements(prhs[0]) != 1) {
		mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notScalar", "Input multiplier must be a scalar.");
	}

	/* make sure the second input argument is type double */
	if (!mxIsDouble(prhs[1]) ||
		mxIsComplex(prhs[1])) {
		mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble", "Input matrix must be type double.");
	}

	/* check that number of rows in second input argument is 1 */
	if (mxGetM(prhs[1]) != 1) {
		mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notRowVector", "Input must be a row vector.");
	}

	/* get the value of the scalar input  */
	multiplier = mxGetScalar(prhs[0]);

	/* create a pointer to the real data in the input matrix  */
	inMatrix = mxGetPr(prhs[1]);

	/* get dimensions of the input matrix */
	ncols = mxGetN(prhs[1]);

	/* create the output matrix */
	plhs[0] = mxCreateDoubleMatrix(1, (mwSize)ncols, mxREAL);

	/* get a pointer to the real data in the output matrix */
	outMatrix = mxGetPr(plhs[0]);

	/* call the computational routine */
	arrayProduct(multiplier, inMatrix, outMatrix, (mwSize)ncols);
}





void mexFunction1(int nlhs, mxArray *plhs[],
	int nrhs, const mxArray *prhs[])
{
	double multiplier;              /* input scalar */
	double *inMatrix;               /* 1xN input matrix */
	size_t ncols;                   /* size of matrix */
	double *outMatrix;              /* output matrix */

									/* check for proper number of arguments */
	if (nrhs != 2) {
		mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs", "Two inputs required.");
	}
	if (nlhs != 1) {
		mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs", "One output required.");
	}
	/* make sure the first input argument is scalar */
	if (!mxIsDouble(prhs[0]) ||
		mxIsComplex(prhs[0]) ||
		mxGetNumberOfElements(prhs[0]) != 1) {
		mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notScalar", "Input multiplier must be a scalar.");
	}

	/* make sure the second input argument is type double */
	if (!mxIsDouble(prhs[1]) ||
		mxIsComplex(prhs[1])) {
		mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notDouble", "Input matrix must be type double.");
	}

	/* check that number of rows in second input argument is 1 */
	if (mxGetM(prhs[1]) != 1) {
		mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notRowVector", "Input must be a row vector.");
	}

	/* get the value of the scalar input  */
	multiplier = mxGetScalar(prhs[0]);

	/* create a pointer to the real data in the input matrix  */
	inMatrix = mxGetPr(prhs[1]);

	/* get dimensions of the input matrix */
	ncols = mxGetN(prhs[1]);

	/* create the output matrix */
	plhs[0] = mxCreateDoubleMatrix(1, (mwSize)ncols, mxREAL);

	/* get a pointer to the real data in the output matrix */
	outMatrix = mxGetPr(plhs[0]);

	/* call the computational routine */
	arrayProduct(multiplier, inMatrix, outMatrix, (mwSize)ncols);
}





