#include "head.h"

void nrerror(char error_text[])
//Numerical Recipes standard error handler
{
	cout << "Numerical Recipes run-time error..." << endl;
	cout << error_text;
	cout << "...now exiting to system" << endl;
	exit(0);
}

int **imatrix(long nrow, long ncol)
//allocate a int matrix with subscript range m[0..nrow][0..ncol]
{
	char temp[50];
	long i;
	int **m;
	//allocate pointers to rows
	m = new int*[nrow];
	if (!m) {
		strcpy_s(temp, "allocation failure 1 in imatrix(,)");
		nrerror(temp);
	}
	//allocate rows and set pointers to them
	m[0] = new int[nrow*ncol];
	if (!m[0]) {
		strcpy_s(temp, "allocation failure 2 in imatrix(,)");
		nrerror(temp);
	}
	for (i = 1; i < nrow; i++) m[i] = m[i - 1] + ncol;
	//return pointer to array of pointers to rows
	return m;
}

double **dmatrix(long nrow, long ncol)
//allocate a double matrix with subscript range m[0..nrow][0..ncol]
{
	char temp[50];
	long i;
	double **m;
	//allocate pointers to rows
	m = new double*[nrow];
	if (!m) {
		strcpy_s(temp, "allocation failure 1 in dmatrix(,)");
		nrerror(temp);
	}
	//allocate rows and set pointers to them
	m[0] = new double[nrow*ncol];
	if (!m[0]) {
		strcpy_s(temp, "allocation failure 2 in dmatrix(,)");
		nrerror(temp);
	}
	for (i = 1; i < nrow; i++) m[i] = m[i - 1] + ncol;
	//return pointer to array of pointers to rows
	return m;
}

complex<double> **cmatrix(long nrow, long ncol)
//allocate a complex<double> matrix with subscript range m[0..nrow][0..ncol]
{
	char temp[50];
	long i;
	complex<double> **m;
	//allocate pointers to rows
	m = new complex<double>*[nrow];
	if (!m) {
		strcpy_s(temp, "allocation failure 1 in cmatrix(,)");
		nrerror(temp);
	}
	//allocate rows and set pointers to them
	m[0] = new complex<double>[nrow*ncol];
	if (!m[0]) {
		strcpy_s(temp, "allocation failure 2 in cmatrix(,)");
		nrerror(temp);
	}
	for (i = 1; i < nrow; i++) m[i] = m[i - 1] + ncol;
	//return pointer to array of pointers to rows
	return m;
}
