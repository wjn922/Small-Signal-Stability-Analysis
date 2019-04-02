#include "head.h"

void matrix_add(int m, int n, double **a, double **b, double **c)
//c[1..m][1..n]=a[1..m][1..n]+b[1..m][1..n]
{
	int i, j;
	for (i = 1; i <= m; i++)
		for (j = 1; j <= n; j++)
			c[i][j] = a[i][j] + b[i][j];
}

void matrix_sub(int m, int n, double **a, double **b, double **c)
//c[1..m][1..n]=a[1..m][1..n]-b[1..m][1..n]
{
	int i, j;
	for (i = 1; i <= m; i++)
		for (j = 1; j <= n; j++)
			c[i][j] = a[i][j] - b[i][j];
}

void matrix_mul(int m, int l, int n, double **a, double **b, double **c)
//c[1..m][1..n]=a[1..m][1..l]*b[1..l][1..n]
{
	int i, j, k;
	double d;
	for (i = 1; i <= m; i++) {
		for (j = 1; j <= n; j++) {
			d = 0.0;
			for (k = 1; k <= l; k++) d += a[i][k] * b[k][j];
			c[i][j] = d;
		}
	}
}

void matrix_Crout(int n, double **a)
//a[1..n][1..n]=L(单位下三角阵)*U(上三角阵)，存入 a
{
	int i, j, k;
	double c, d;
	for (j = 1; j <= n; j++) {
		for (i = 1; i <= j; i++) {
			d = a[i][j];
			for (k = 1; k <= i - 1; k++) d -= a[i][k] * a[k][j];
			a[i][j] = d;
		}
		c = a[j][j];
		for (i = j + 1; i <= n; i++) {
			d = a[i][j];
			for (k = 1; k <= j - 1; k++) d -= a[i][k] * a[k][j];
			a[i][j] = d / c;
		}
	}
}