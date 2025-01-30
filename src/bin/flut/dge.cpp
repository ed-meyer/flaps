//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

// test the interval versions of dgefa, dgesl


#include "blas.h"
#include "interval.h"
#include "trace.h"

/* Table of constant values */

static int c__1 = 1;

int
dgefa_(double *a, int *lda, int *n, int *ipvt, int *info) {
    /* System generated locals */
    int a_dim1, a_offset, i__1, i__2, i__3;

    /* Local variables */
    int j, k, l;
    double t;
    int kp1, nm1;
    //!! extern /* Subroutine */ int dscal_(int *, double *, double *, 
	    //!! int *), daxpy_(int *, double *, double *, int 
	    //!! *, double *, int *);
    //!! extern int idamax_(int *, double *, int *);


/*     dgefa factors a double precision matrix by gaussian elimination. */

/*     dgefa is usually called by dgeco, but it can be called */
/*     directly with a saving in time if  rcond  is not needed. */
/*     (time for dgeco) = (1 + 9/n)*(time for dgefa) . */

/*     on entry */

/*        a       double precision(lda, n) */
/*                the matrix to be factored. */

/*        lda     int */
/*                the leading dimension of the array  a . */

/*        n       int */
/*                the order of the matrix  a . */

/*     on return */

/*        a       an upper triangular matrix and the multipliers */
/*                which were used to obtain it. */
/*                the factorization can be written  a = l*u  where */
/*                l  is a product of permutation and unit lower */
/*                triangular matrices and  u  is upper triangular. */

/*        ipvt    int(n) */
/*                an int vector of pivot indices. */

/*        info    int */
/*                = 0  normal value. */
/*                = k  if  u(k,k) .eq. 0.0 .  this is not an error */
/*                     condition for this subroutine, but it does */
/*                     indicate that dgesl or dgedi will divide by zero */
/*                     if called.  use  rcond  in dgeco for a reliable */
/*                     indication of singularity. */

/*     linpack. this version dated 08/14/78 . */
/*     cleve moler, university of new mexico, argonne national lab. */

/*     subroutines and functions */

/*     blas daxpy,dscal,idamax */

/*     internal variables */



/*     gaussian elimination with partial pivoting */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --ipvt;

    /* Function Body */
    *info = 0;
    nm1 = *n - 1;
    if (nm1 < 1) {
	goto L70;
    }
    i__1 = nm1;
    for (k = 1; k <= i__1; ++k) {
	kp1 = k + 1;

/*        find l = pivot index */

	i__2 = *n - k + 1;
	l = idamax_(&i__2, &a[k + k * a_dim1], &c__1) + k - 1;
	ipvt[k] = l;

/*        zero pivot implies this column already triangularized */

	if (a[l + k * a_dim1] == 0.) {
	    goto L40;
	}

/*           interchange if necessary */

	if (l == k) {
	    goto L10;
	}
	t = a[l + k * a_dim1];
	a[l + k * a_dim1] = a[k + k * a_dim1];
	a[k + k * a_dim1] = t;
L10:

/*           compute multipliers */

	t = -1. / a[k + k * a_dim1];
	i__2 = *n - k;
	dscal_(&i__2, &t, &a[k + 1 + k * a_dim1], &c__1);

/*           row elimination with column indexing */

	i__2 = *n;
	for (j = kp1; j <= i__2; ++j) {
	    t = a[l + j * a_dim1];
	    if (l == k) {
		goto L20;
	    }
	    a[l + j * a_dim1] = a[k + j * a_dim1];
	    a[k + j * a_dim1] = t;
L20:
	    i__3 = *n - k;
	    daxpy_(&i__3, &t, &a[k + 1 + k * a_dim1], &c__1, &a[k + 1 + j * 
		    a_dim1], &c__1);
/* L30: */
	}
	goto L50;
L40:
	*info = k;
L50:
/* L60: */
	;
    }
L70:
    ipvt[*n] = *n;
    if (a[*n + *n * a_dim1] == 0.) {
	*info = *n;
    }
    return 0;
} /* dgefa_ */

#ifdef __cplusplus
//!! 	}
#endif

#ifdef __cplusplus
//!! extern "C" {
#endif
//!! #include "f2c.h"

/* Table of constant values */

//!! static int c__1 = 1;

int
dgesl_(double *a, int *lda, int *n, int *ipvt, Interval* b) {
    /* System generated locals */
    int a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    int k, l;
    Interval t;
    int kb, nm1;


/*     dgesl solves the double precision system */
/*     a * x = b  or  trans(a) * x = b */
/*     using the factors computed by dgeco or dgefa. */

/*     on entry */

/*        a       double precision(lda, n) */
/*                the output from dgeco or dgefa. */

/*        lda     int */
/*                the leading dimension of the array  a . */

/*        n       int */
/*                the order of the matrix  a . */

/*        ipvt    int(n) */
/*                the pivot vector from dgeco or dgefa. */

/*        b       double precision(n) */
/*                the right hand side vector. */

/*        job     int */
/*                = 0         to solve  a*x = b , */
/*                = nonzero   to solve  trans(a)*x = b  where */
/*                            trans(a)  is the transpose. */

/*     on return */

/*        b       the solution vector  x . */

/*     error condition */

/*        a division by zero will occur if the input factor contains a */
/*        zero on the diagonal.  technically this indicates singularity */
/*        but it is often caused by improper arguments or improper */
/*        setting of lda .  it will not occur if the subroutines are */
/*        called correctly and if dgeco has set rcond .gt. 0.0 */
/*        or dgefa has set info .eq. 0 . */

/*     to compute  inverse(a) * c  where  c  is a matrix */
/*     with  p  columns */
/*           call dgeco(a,lda,n,ipvt,rcond,z) */
/*           if (rcond is too small) go to ... */
/*           do 10 j = 1, p */
/*              call dgesl(a,lda,n,ipvt,c(1,j),0) */
/*        10 continue */

/*     linpack. this version dated 08/14/78 . */
/*     cleve moler, university of new mexico, argonne national lab. */

/*     subroutines and functions */

/*     blas daxpy,ddot */

/*     internal variables */


    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --ipvt;
    --b;

    /* Function Body */
    nm1 = *n - 1;

/*        job = 0 , solve  a * x = b */
/*        first solve  l*y = b */

    if (nm1 < 1) {
	goto L30;
    }
    i__1 = nm1;
    for (k = 1; k <= i__1; ++k) {
	l = ipvt[k];
	t = b[l];
	if (l == k) {
	    goto L10;
	}
	b[l] = b[k];
	b[k] = t;
L10:
	i__2 = *n - k;
	//!!  daxpy_(&i__2, &t, &a[k + 1 + k * a_dim1], &c__1, &b[k + 1], &c__1);
	for (int i=0; i<i__2; i++)
		b[k+1+i] += t*a[k+1+k*a_dim1+i];
/* L20: */
    }
L30:

/*        now solve  u*x = y */

    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
	b[k] /= a[k + k * a_dim1];
	t = -b[k];
	i__2 = k - 1;
	//!! daxpy_(&i__2, &t, &a[k * a_dim1 + 1], &c__1, &b[1], &c__1);
	for (int i=0; i<i__2; i++)
		b[1+i] += t*a[k*a_dim1+1+i];
/* L40: */
    }
    return 0;
} /* dgesl_ */

#ifdef __cplusplus
	//!! }
#endif

#ifdef MAIN

#include <iostream>
#include <vector>
#include "fptype.h"

using namespace std;

void
makeab(vector<double>& a, vector<Interval>& b) {
	size_t n = b.size();
	assert(a.size() == n*n);
	for (size_t i=0; i<n; i++) {
		for (size_t j=0; j<n; j++) {
			double x = flaps::xrand();
			if (i != j)
				x *= 0.1;
			a[i+j*n] = x;
		}
		b[i] = flaps::xrand();
	}
}

int
main(int argc, char** argv) {
	T_(Trace trc(0,"dge");)
	int n{5};
	vector<double> a(n*n);
	vector<Interval> b(n);
	vector<int> ipvt(n);
	makeab(a, b);
	vector<double> af{a};
	vector<Interval> x{b};
	int info;
	T_(trc.dprintm(n,n,n,a,"a");)
	T_(trc.dprint("b",b);)
	// dgefa_(double *a, int *lda, int *n, int *ipvt, int *info);
	dgefa_(&af[0], &n, &n, &ipvt[0], &info);
	// dgesl_(double *a, int *lda, int *n, int *ipvt, Interval* b);
	dgesl_(&af[0], &n, &n, &ipvt[0], &x[0]);

	// check
	vector<Interval> resid(n, Interval(0.0));
	for (int i=0; i<n; i++) {
		for (int j=0; j<n; j++) {
			resid[i] += a[i+j*n]*x[j];
		}
		resid[i] -= b[i];
	}

	T_(trc.dprint("residual: ",resid);)

}
#endif // MAIN
