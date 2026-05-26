/*
   diag-c.h
   global declarations for the Diag routines in C
   this file is part of Diag
   last modified 20 Aug 15 th
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h> 
#include "matr_sqrt.h"

static inline RealType Sq(cComplexType x) { return Re(x*Conjugate(x)); }
static inline RealType sq(cRealType x) { return x*x; }
static inline RealType min(cRealType a, cRealType b) { return (a < b) ? a : b; }
static inline RealType max(cRealType a, cRealType b) { return (a > b) ? a : b; }
static inline int imin(cint a, cint b) { return (a < b) ? a : b; }
static inline int imax(cint a, cint b) { return (a > b) ? a : b; }

int nsweeps;

/*
   TakagiFactor.c
   computes the Takagi factorization of a complex symmetric matrix
   code adapted from the "Handbook" routines
   (Wilkinson, Reinsch: Handbook for Automatic Computation, p. 202)
   this file is part of the Diag library
   last modified 20 Aug 15 th
*/

/*
   TakagiFactor factorizes a complex symmetric n-by-n matrix
   Input:	n, A = n-by-n matrix, complex symmetric
		(only the upper triangle of A needs to be filled),
   Output:	d = vector of diagonal values,
		U = transformation matrix, unitary (U^-1 = U^+),
   these fulfill
	d = U^* A U^+,  A = U^T d U,  U^* A = d U  (UCOLS=0),
	d = U^+ A U^*,  A = U d U^T,  A U^* = U d  (UCOLS=1).
*/

void TakagiFactor(cint n, ComplexType *A, cint ldA,
  RealType *d, ComplexType *U, cint ldU, cint sort)
{
  int p, q;
  cRealType red = .04/(n*n*n*n);
  ComplexType ev[n][2];

  for( p = 0; p < n; ++p ) {
    ev[p][0] = 0;
    ev[p][1] = A(p,p);
  }

  for( p = 0; p < n; ++p ) {
    memset(&U(p,0), 0, n*sizeof(ComplexType));
    U(p,p) = 1;
  }

  for( nsweeps = 1; nsweeps <= 50; ++nsweeps ) {
    RealType thresh = 0;
    for( q = 1; q < n; ++q )
      for( p = 0; p < q; ++p )
        thresh += Sq(A(p,q));
    if( !(thresh > SYM_EPS) ) goto done;

    thresh = (nsweeps < 4) ? thresh*red : 0;

    for( q = 1; q < n; ++q )
      for( p = 0; p < q; ++p ) {
        cComplexType Apq = A(p,q);
        cRealType off = Sq(Apq);
        cRealType sqp = Sq(ev[p][1]);
        cRealType sqq = Sq(ev[q][1]);
        if( nsweeps > 4 && off < SYM_EPS*(sqp + sqq) )
          A(p,q) = 0;
        else if( off > thresh ) {
          RealType t, invc;
          ComplexType f;
          int j;

          t = .5*absr(sqp - sqq);
          if( t > DBL_EPS )
            f = sign(1, sqp - sqq)*
              (ev[q][1]*Conjugate(Apq) + Conjugate(ev[p][1])*Apq);
          else
            f = (sqp == 0) ? 1 : Sqrt(ev[q][1]/ev[p][1]);
          t += sqrt(t*t + Sq(f));
          f /= t;

          ev[p][1] = A(p,p) + (ev[p][0] += Apq*Conjugate(f));
          ev[q][1] = A(q,q) + (ev[q][0] -= Apq*f);

          t = Sq(f);
          invc = sqrt(t + 1);
          f /= invc;
          t /= invc*(invc + 1);

          for( j = 0; j < p; ++j ) {
            cComplexType x = A(j,p);
            cComplexType y = A(j,q);
            A(j,p) = x + (Conjugate(f)*y - t*x);
            A(j,q) = y - (f*x + t*y);
          }

          for( j = p + 1; j < q; ++j ) {
            cComplexType x = A(p,j);
            cComplexType y = A(j,q);
            A(p,j) = x + (Conjugate(f)*y - t*x);
            A(j,q) = y - (f*x + t*y);
          }

          for( j = q + 1; j < n; ++j ) {
            cComplexType x = A(p,j);
            cComplexType y = A(q,j);
            A(p,j) = x + (Conjugate(f)*y - t*x);
            A(q,j) = y - (f*x + t*y);
          }

          A(p,q) = 0;

          for( j = 0; j < n; ++j ) {
            cComplexType x = UL(p,j);
            cComplexType y = UL(q,j);
            UL(p,j) = x + (f*y - t*x);
            UL(q,j) = y - (Conjugate(f)*x + t*y);
          }
        }
      }

    for( p = 0; p < n; ++p ) {
      ev[p][0] = 0;
      A(p,p) = ev[p][1];
    }
  }

  fputs("Bad convergence in TakagiFactor\n", stderr);

done:

/* make the diagonal elements nonnegative */

  for( p = 0; p < n; ++p ) {
    cComplexType App = A(p,p);
    d[p] = Abs(App);
    if( d[p] > DBL_EPS && d[p] != Re(App) ) {
      cComplexType f = Sqrt(App/d[p]);
      for( q = 0; q < n; ++q ) UL(p,q) *= f;
    }
  }

  if( sort == 0 ) return;

/* sort the eigenvalues */

  for( p = 0; p < n - 1; ++p ) {
    int j = p;
    RealType t = d[p];
    for( q = p + 1; q < n; ++q )
      if( sort*(t - d[q]) > 0 ) t = d[j = q];
    if( j == p ) continue;
    d[j] = d[p];
    d[p] = t;
    for( q = 0; q < n; ++q ) {
      cComplexType x = UL(p,q);
      UL(p,q) = UL(j,q);
      UL(j,q) = x;
    }
  }
}

/*Here follows the use of Schur decomposition for square root of polarizability tensor*/

typedef struct {
	int m, n;
	double complex ** v;
} mat_t, *mat;

mat matrix_new(int m, int n)
{
	mat x = malloc(sizeof(mat_t));
	x->v = malloc(sizeof(double complex*) * m);
	x->v[0] = calloc(m * n, sizeof(double complex));
	for (int i = 0; i < m; i++)
		x->v[i] = x->v[0] + n * i;
	x->m = m;
	x->n = n;
	return x;
}

void matrix_delete(mat m)
{
	free(m->v[0]);
	free(m->v);
	free(m);
}

void matrix_hermitian(mat m)
{
	for (int i = 0; i < m->m; i++) {
		for (int j = 0; j < i; j++) {
			double complex t = conj(m->v[i][j]);
			m->v[i][j] = conj(m->v[j][i]);
			m->v[j][i] = t;
		}
		m->v[i][i] = conj(m->v[i][i]);
	}
}

mat matrix_copy(int n, double complex a[][n], int m)
{
	mat x = matrix_new(m, n);
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			x->v[i][j] = a[i][j];
	return x;
}

mat matrix_mul(mat x, mat y)
{
	if (x->n != y->m) return 0;
	mat r = matrix_new(x->m, y->n);
	for (int i = 0; i < x->m; i++)
		for (int j = 0; j < y->n; j++)
			for (int k = 0; k < x->n; k++)
				r->v[i][j] += x->v[i][k] * y->v[k][j];
	return r;
}

mat matrix_minor(mat x, int d)
{
	mat m = matrix_new(x->m, x->n);
	for (int i = 0; i < d; i++)
		m->v[i][i] = 1.0;
	for (int i = d; i < x->m; i++)
		for (int j = d; j < x->n; j++)
			m->v[i][j] = x->v[i][j];
	return m;
}

double complex *vmadd(double complex a[], double complex b[], double complex s, double complex c[], int n)
{
	for (int i = 0; i < n; i++)
		c[i] = a[i] + s * b[i];
	return c;
}

mat vmul(double complex v[], int n)
{
	mat x = matrix_new(n, n);
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			x->v[i][j] = -2.0 * v[i] * conj(v[j]);
	for (int i = 0; i < n; i++)
		x->v[i][i] += 1.0;
	return x;
}

double vnorm(double complex x[], int n)
{
	double sum = 0;
	for (int i = 0; i < n; i++)
		sum += (conj(x[i]) * x[i]);
	return sqrt(sum);
}

double complex* vdiv(double complex x[], double d, double complex y[], int n)
{
	for (int i = 0; i < n; i++) y[i] = x[i] / d;
	return y;
}

double complex* mcol(mat m, double complex *v, int c)
{
	for (int i = 0; i < m->m; i++)
		v[i] = m->v[i][c];
	return v;
}

void householder(mat m, mat *R, mat *Q)
{
	mat q[m->m];
	mat z = m, z1;

	for (int k = 0; k < m->n && k < m->m - 1; k++) {
		double complex e[m->m], x[m->m], a;

		z1 = matrix_minor(z, k);
		if (z != m) matrix_delete(z);
		z = z1;

		mcol(z, x, k);
		a = vnorm(x, m->m);
		if (creal(z->v[k][k]) < 0) a = -a;

		for (int i = 0; i < m->m; i++)
			e[i] = (i == k) ? 1.0 : 0.0;

		vmadd(x, e, a, e, m->m);
		vdiv(e, vnorm(e, m->m), e, m->m);

		q[k] = vmul(e, m->m);

		z1 = matrix_mul(q[k], z);
		if (z != m) matrix_delete(z);
		z = z1;
	}

	matrix_delete(z);

	*Q = q[0];
	*R = matrix_mul(q[0], m);

	for (int i = 1; i < m->n && i < m->m - 1; i++) {
		z1 = matrix_mul(q[i], *Q);
		if (i > 1) matrix_delete(*Q);
		*Q = z1;
		matrix_delete(q[i]);
	}

	matrix_delete(q[0]);

	z = matrix_mul(*Q, m);
	matrix_delete(*R);
	*R = z;

	matrix_hermitian(*Q);
}

void SchurDecomposition(int n,double complex _A[][n],double complex _U[][n],double complex _D[][n])
{
    mat Q, R;
    mat A = matrix_copy(n, _A, n);
    mat U = matrix_new(n, n);
		int max_iter = 1e3;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            U->v[i][j] = (i == j);
    for (int iter = 0; iter < max_iter; iter++) {
        householder(A, &R, &Q);
        mat tmpA = matrix_mul(R, Q);
        matrix_delete(A);
        A = tmpA;
        mat tmpU = matrix_mul(U, Q);
        matrix_delete(U);
        U = tmpU;
        double off = 0;
        for (int i = 1; i < n; i++)
            for (int j = 0; j < i; j++)
                off += (conj(A->v[i][j]) * A->v[i][j]);
        if (sqrt(off) < __DBL_EPSILON__)
            break;
    }
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            _D[i][j] = A->v[i][j];
            _U[i][j] = U->v[i][j];
        }
    matrix_delete(A);
    matrix_delete(U);
}

void TriSqrt(double complex D[3][3], double complex S[3][3]) {
	/*This function is used for Schur decomposition (see p. 136 of https://doi.org/10.1137/1.9780898717778)*/
    double complex tmp[3][3];
    for(int i = 0; i < 3; i++)
        for(int j = 0; j < 3; j++) S[i][j] = 0.0;
    for(int i = 0; i < 3; i++) {
        S[i][i] = csqrt(D[i][i]);
        for(int j = i + 1; j < 3; j++) {
            tmp[i][j] = 0.0; 
            for(int k = i + 1; k < j; k++)
                tmp[i][j] += S[i][k] * S[k][j];
            S[i][j] = (D[i][j] - tmp[i][j]) / (S[i][i] + S[j][j]);
        }
    }
}