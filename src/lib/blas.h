//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

#ifndef blas_h
#define blas_h

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "Ad.h"
#include "trace.h"
// A collection of various ways to access the Basic Linear Algebra Subroutines
// (BLAS) from C++, plus some additional linear-algebraic functions:
// 1) double-valued functions with names like blas_copy
// 2) templated 
//   real_rep:  create the real (2m,2n) representation of a complex (m,n) matrix

// Templatized version of the BLAS
// Specializations in blas.c

template<typename Type>
void
blas_copy (int n, Type const* x, int incx, Type* y, int incy) {
	int ix = 0;
	int iy = 0;
	for (int i=0; i<n; i++) {
		y[iy] = x[ix];
		iy += incy;
		ix += incx;
	}
}

template<typename Type>
void
blas_axpy (int n, Type const& a, Type const* x, int incx, Type* y, int incy) {
// y = a*x + y
	int ix = 0;
	int iy = 0;
	for (int i=0; i<n; i++) {
		y[iy] += a*x[ix];
		iy += incy;
		ix += incx;
	}
}

template<typename Type>
void
blas_scal (int n, Type const& a, Type* x, int incx) {
// x <= a*x
	int ix = 0;
	for (int i=0; i<n; i++) {
		x[ix] *= a;
		ix += incx;
	}
}

template<typename Type>
void
blas_dot (int n, Type const* x, int incx, Type const* y, int incy, Type& z) {
// z = x'*y
	Trace trc(3,"blas_dot");
	int ix = 0;
	int iy = 0;
	z = 0.0;
	for (int i=0; i<n; i++) {
		z += x[ix]*y[iy];
		trc.dprint("i ",i," x ",x[ix],", y ",y[iy],", z ",z);
		ix += incx;
		iy += incy;
	}
}

template<typename Type>
void
blas_cdotc (int n, std::complex<Type> const* x, int incx, std::complex<Type> const* y, int incy, std::complex<Type>& z) {
// z = x(*)*y
	int ix = 0;
	int iy = 0;
	z = 0.0;
	for (int i=0; i<n; i++) {
		z += conj(x[ix])*y[iy];
		ix += incx;
		iy += incy;
	}
}

template<typename Type>
void
blas_cdotu (int n, std::complex<Type> const* x, int incx, std::complex<Type> const* y, int incy, std::complex<Type>& z) {
// z = x(*)*y
	int ix = 0;
	int iy = 0;
	z = 0.0;
	for (int i=0; i<n; i++) {
		z += x[ix]*y[iy];
		ix += incx;
		iy += incy;
	}
}


template<typename Type>
double
blas_snrm2 (int n, Type const* a, int inca) {
// 2-norm for real data types: sqrt of the dot product
	Type rval(0.0);
	blas_dot (n, a, inca, a, inca, rval);
	return sqrt(rval);
}
// specialization for double
template<>
double
blas_snrm2 (int n, double const* a, int inca);


template<typename Type>
double
blas_scnrm2 (int n, std::complex<Type> const* a, int inca) {
// 2-norm for complex data types: sqrt of the conjugated dot product
	double rval(0.0);
	int ia = 0;
	for (int i=0; i<n; i++) {
		Type ar = a[i].real();
		Type ai = a[i].imag();
		rval += ar*ar + ai*ai;
		ia += inca;
	}
	return sqrt(rval);
}
// specialization for double
template<>
double
blas_scnrm2 (int n, std::complex<double> const* a, int inca);

template<typename Type>
int
blas_isamax (int n, Type const* a, int inca) {
// 1-based component number of the largest abs element
	int rval = 1;
	int ia = 0;
	double amax(0.0);
	for (int i=0; i<n; i++) {
		if (std::abs(a[i]) > amax) {
			amax = std::abs(a[i]);
			rval = i+1;
		}
		ia += inca;
	}
	return rval;
}

template<typename Type>
int
blas_icamax (int n, std::complex<Type> const* a, int inca) {
// this function is unnecessary: identical to isamax - replace
// both with blas_iamax?
	int rval = 1;
	int ia = 0;
	double amax(0.0);
	for (int i=0; i<n; i++) {
		if (std::abs(a[i]) > amax) {
			amax = std::abs(a[i]);
			rval = i+1;
		}
		ia += inca;
	}
	return rval;
}

template<typename Type>
double
blas_sasum (int n, Type const* a, int inca) {
	double rval(0.0);
	int ia = 0;
	for (int i=0; i<n; i++) {
		rval += std::abs(a[i]);
		ia += inca;
	}
	return rval;
}

template<typename Type>
double
blas_scasum (int n, std::complex<Type> const* a, int inca) {
	double rval(0.0);
	int ia = 0;
	for (int i=0; i<n; i++) {
		rval += std::abs(a[i]);
		ia += inca;
	}
	return rval;
}



inline
bool
lsame (char const* a, char const* b) { return tolower(*a) == tolower(*b); }

template<typename Type>
void
conjugate (std::complex<Type>* r, std::complex<Type> const* z) {
	Type zi = z->imag();
	*r = std::complex<Type>(z->real(), -zi);
}


template<typename Type>
int
blas_sgemv(char const* trans, int m, int n, Type alpha, Type const* a,
		int lda, Type const* x, int incx, Type beta, Type* y, int incy) {

    int a_offset, i__1, i__2;

    /* Local variables */
    int i__, j, ix, iy, jx, jy, kx, ky, info;
    Type temp;
    int lenx, leny;

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  SGEMV  performs one of the matrix-vector operations */

/*     y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y, */

/*  where alpha and beta are scalars, x and y are vectors and A is an */
/*  m by n matrix. */

/*  Arguments */
/*  ========== */

/*  TRANS  - CHARACTER*1. */
/*           On entry, TRANS specifies the operation to be performed as */
/*           follows: */

/*              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y. */

/*              TRANS = 'T' or 't'   y := alpha*A**T*x + beta*y. */

/*              TRANS = 'C' or 'c'   y := alpha*A**T*x + beta*y. */

/*           Unchanged on exit. */

/*  M      - int. */
/*           On entry, M specifies the number of rows of the matrix A. */
/*           M must be at least zero. */
/*           Unchanged on exit. */

/*  N      - int. */
/*           On entry, N specifies the number of columns of the matrix A. */
/*           N must be at least zero. */
/*           Unchanged on exit. */

/*  ALPHA  - REAL            . */
/*           On entry, ALPHA specifies the scalar alpha. */
/*           Unchanged on exit. */

/*  A      - REAL             array of DIMENSION ( LDA, n ). */
/*           Before entry, the leading m by n part of the array A must */
/*           contain the matrix of coefficients. */
/*           Unchanged on exit. */

/*  LDA    - int. */
/*           On entry, LDA specifies the first dimension of A as declared */
/*           in the calling (sub) program. LDA must be at least */
/*           max( 1, m ). */
/*           Unchanged on exit. */

/*  X      - REAL             array of DIMENSION at least */
/*           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n' */
/*           and at least */
/*           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise. */
/*           Before entry, the incremented array X must contain the */
/*           vector x. */
/*           Unchanged on exit. */

/*  INCX   - int. */
/*           On entry, INCX specifies the increment for the elements of */
/*           X. INCX must not be zero. */
/*           Unchanged on exit. */

/*  BETA   - REAL            . */
/*           On entry, BETA specifies the scalar beta. When BETA is */
/*           supplied as zero then Y need not be set on input. */
/*           Unchanged on exit. */

/*  Y      - REAL             array of DIMENSION at least */
/*           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n' */
/*           and at least */
/*           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise. */
/*           Before entry with BETA non-zero, the incremented array Y */
/*           must contain the vector y. On exit, Y is overwritten by the */
/*           updated vector y. */

/*  INCY   - int. */
/*           On entry, INCY specifies the increment for the elements of */
/*           Y. INCY must not be zero. */
/*           Unchanged on exit. */

/*  Further Details */
/*  =============== */

/*  Level 2 Blas routine. */
/*  The vector and matrix arguments are not referenced when N = 0, or M = 0 */

/*  -- Written on 22-October-1986. */
/*     Jack Dongarra, Argonne National Lab. */
/*     Jeremy Du Croz, Nag Central Office. */
/*     Sven Hammarling, Nag Central Office. */
/*     Richard Hanson, Sandia National Labs. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */

/*     Test the input parameters. */

    /* Parameter adjustments */
    a_offset = 1 + lda;
    a -= a_offset;
    --x;
    --y;

    /* Function Body */
    info = 0;
    if (! lsame(trans, "N") && ! lsame(trans, "T") && ! lsame(trans, "C")) {
	info = 1;
    } else if (m < 0) {
	info = 2;
    } else if (n < 0) {
	info = 3;
    } else if (lda < std::max(1,m)) {
	info = 6;
    } else if (incx == 0) {
	info = 8;
    } else if (incy == 0) {
	info = 11;
    }
    if (info != 0) {
		// errorMessage("blas_sgemv %s", info);
		return info;
    }

/*     Quick return if possible. */

	if ((m == 0) || (n == 0) ||
		((alpha == Type(0.0)) && (beta - Type(1.0)) == Type(0.0))) {
		return 0;
	}

/*     Set  LENX  and  LENY, the lengths of the vectors x and y, and set */
/*     up the start points in  X  and  Y. */

    if (lsame(trans, "N")) {
	lenx = n;
	leny = m;
    } else {
	lenx = m;
	leny = n;
    }
    if (incx > 0) {
	kx = 1;
    } else {
	kx = 1 - (lenx - 1) * incx;
    }
    if (incy > 0) {
	ky = 1;
    } else {
	ky = 1 - (leny - 1) * incy;
    }

/*     Start the operations. In this version the elements of A are */
/*     accessed sequentially with one pass through A. */

/*     First form  y := beta*y. */

	if (beta != Type(1.0)) {
		if (incy == 1) {
			if (beta == Type(0.0)) {
				i__1 = leny;
				for (i__ = 1; i__ <= i__1; ++i__) {
					y[i__] = Type(0.0);
/* L10: */
				}
			} else {
				i__1 = leny;
				for (i__ = 1; i__ <= i__1; ++i__) {
					y[i__] = beta * y[i__];
/* L20: */
				}
			}
		} else {
			iy = ky;
			if (beta == Type(0.0)) {
				i__1 = leny;
				for (i__ = 1; i__ <= i__1; ++i__) {
					y[iy] = Type(0.0);
					iy += incy;
/* L30: */
				}
			} else {
				i__1 = leny;
				for (i__ = 1; i__ <= i__1; ++i__) {
					y[iy] = beta * y[iy];
					iy += incy;
/* L40: */
				}
			}
		}
	}
	if (alpha == Type(0.0)) {
		return 0;
	}
	if (lsame(trans, "N")) {

/*        Form  y := alpha*A*x + y. */

	jx = kx;
	if (incy == 1) {
	    i__1 = n;
	    for (j = 1; j <= i__1; ++j) {
		if (x[jx] != Type(0.0)) {
		    temp = alpha * x[jx];
		    i__2 = m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			y[i__] += temp * a[i__ + j * lda];
/* L50: */
		    }
		}
		jx += incx;
/* L60: */
	    }
	} else {
	    i__1 = n;
	    for (j = 1; j <= i__1; ++j) {
		if (x[jx] != Type(0.0)) {
		    temp = alpha * x[jx];
		    iy = ky;
		    i__2 = m;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			y[iy] += temp * a[i__ + j * lda];
			iy += incy;
/* L70: */
		    }
		}
		jx += incx;
/* L80: */
	    }
	}
    } else {

/*        Form  y := alpha*A**T*x + y. */

	jy = ky;
	if (incx == 1) {
	    i__1 = n;
	    for (j = 1; j <= i__1; ++j) {
		temp = Type(0.0);
		i__2 = m;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    temp += a[i__ + j * lda] * x[i__];
/* L90: */
		}
		y[jy] += alpha * temp;
		jy += incy;
/* L100: */
	    }
	} else {
	    i__1 = n;
	    for (j = 1; j <= i__1; ++j) {
		temp = Type(0.0);
		ix = kx;
		i__2 = m;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    temp += a[i__ + j * lda] * x[ix];
		    ix += incx;
/* L110: */
		}
		y[jy] += alpha * temp;
		jy += incy;
/* L120: */
	    }
	}
    }

    return 0;
} /* blas_sgemv_ */

// specialization for double
template<>
int
blas_sgemv<double> (char const* trans, int m, int n, double alpha, double const* a,
		int lda, double const* x, int incx, double beta, double *y, int incy);



template<typename Type>
int
blas_cgemv(char const* trans, int m, int n, std::complex<Type> alpha,
		std::complex<Type> const* a, int lda, std::complex<Type> const* x,
		int incx, std::complex<Type> beta, std::complex<Type>* y, int incy) {

    int a_offset, i__1, i__2, i__3, i__4, i__5;
	 std::complex<Type> q__1, q__2, q__3;

    /* Local variables */
    int i__, j, ix, iy, jx, jy, kx, ky, info;
	 std::complex<Type> temp;
    int lenx, leny;
    bool noconj = false;

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  CGEMV performs one of the matrix-vector operations */

/*     y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,   or */

/*     y := alpha*A**H*x + beta*y, */

/*  where alpha and beta are scalars, x and y are vectors and A is an */
/*  m by n matrix. */

/*  Arguments */
/*  ========== */

/*  TRANS  - CHARACTER*1. */
/*           On entry, TRANS specifies the operation to be performed as */
/*           follows: */

/*              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y. */

/*              TRANS = 'T' or 't'   y := alpha*A**T*x + beta*y. */

/*              TRANS = 'C' or 'c'   y := alpha*A**H*x + beta*y. */

/*           Unchanged on exit. */

/*  M      - int. */
/*           On entry, M specifies the number of rows of the matrix A. */
/*           M must be at least zero. */
/*           Unchanged on exit. */

/*  N      - int. */
/*           On entry, N specifies the number of columns of the matrix A. */
/*           N must be at least zero. */
/*           Unchanged on exit. */

/*  ALPHA  - COMPLEX         . */
/*           On entry, ALPHA specifies the scalar alpha. */
/*           Unchanged on exit. */

/*  A      - COMPLEX          array of DIMENSION ( LDA, n ). */
/*           Before entry, the leading m by n part of the array A must */
/*           contain the matrix of coefficients. */
/*           Unchanged on exit. */

/*  LDA    - int. */
/*           On entry, LDA specifies the first dimension of A as declared */
/*           in the calling (sub) program. LDA must be at least */
/*           max( 1, m ). */
/*           Unchanged on exit. */

/*  X      - COMPLEX          array of DIMENSION at least */
/*           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n' */
/*           and at least */
/*           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise. */
/*           Before entry, the incremented array X must contain the */
/*           vector x. */
/*           Unchanged on exit. */

/*  INCX   - int. */
/*           On entry, INCX specifies the increment for the elements of */
/*           X. INCX must not be zero. */
/*           Unchanged on exit. */

/*  BETA   - COMPLEX         . */
/*           On entry, BETA specifies the scalar beta. When BETA is */
/*           supplied as zero then Y need not be set on input. */
/*           Unchanged on exit. */

/*  Y      - COMPLEX          array of DIMENSION at least */
/*           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n' */
/*           and at least */
/*           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise. */
/*           Before entry with BETA non-zero, the incremented array Y */
/*           must contain the vector y. On exit, Y is overwritten by the */
/*           updated vector y. */

/*  INCY   - int. */
/*           On entry, INCY specifies the increment for the elements of */
/*           Y. INCY must not be zero. */
/*           Unchanged on exit. */

/*  Further Details */
/*  =============== */

/*  Level 2 Blas routine. */
/*  The vector and matrix arguments are not referenced when N = 0, or M = 0 */

/*  -- Written on 22-October-1986. */
/*     Jack Dongarra, Argonne National Lab. */
/*     Jeremy Du Croz, Nag Central Office. */
/*     Sven Hammarling, Nag Central Office. */
/*     Richard Hanson, Sandia National Labs. */

/*  ===================================================================== */

/*     .. Parameters .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */

/*     Test the input parameters. */

    /* Parameter adjustments */
    a_offset = 1 + lda;
    a -= a_offset;
    --x;
    --y;

    /* Function Body */
   info = 0;
	if (! lsame(trans, "N") && ! lsame(trans, "T")
	    && ! lsame(trans, "C")) {
		info = 1;
	} else if (m < 0) {
	info = 2;
    } else if (n < 0) {
	info = 3;
    } else if (lda < std::max(1,m)) {
	info = 6;
    } else if (incx == 0) {
	info = 8;
    } else if (incy == 0) {
	info = 11;
    }
	if (info != 0) {
		// errorMessage("blas_cgemv %d",info);
		return info;
	}

/*     Quick return if possible. */

	if ((m == 0) || (n == 0) || ((alpha == std::complex<Type>(Type(0.0))) 
	    && (beta == std::complex<Type>(Type(1.0))))) {
		return 0;
    }

    noconj = lsame(trans, "T");

/*     Set  LENX  and  LENY, the lengths of the vectors x and y, and set */
/*     up the start points in  X  and  Y. */

	if (lsame(trans, "N")) {
		lenx = n;
		leny = m;
	} else {
		lenx = m;
		leny = n;
	}
	if (incx > 0) {
	kx = 1;
    } else {
	kx = 1 - (lenx - 1) * incx;
    }
    if (incy > 0) {
	ky = 1;
    } else {
	ky = 1 - (leny - 1) * incy;
    }

/*     Start the operations. In this version the elements of A are */
/*     accessed sequentially with one pass through A. */

/*     First form  y := beta*y. */

	if ((beta.real() - Type(1.0)) != 0.0 ||
			(beta.imag() != Type(0.0))) {
		if (incy == 1) {
			if ((beta.real() == 0.0) && (beta.imag() == 0.0)) {
				i__1 = leny;
				for (i__ = 1; i__ <= i__1; ++i__) {
					i__2 = i__;
					y[i__2] = std::complex<Type>(Type(0.0));
/* L10: */
				}
			} else {
				i__1 = leny;
				for (i__ = 1; i__ <= i__1; ++i__) {
					i__2 = i__;
					i__3 = i__;
					q__1 = beta * y[i__3];
					y[i__2] = q__1;
/* L20: */
				}
			}
	} else {
		iy = ky;
		if (std::abs(beta) == 0.0) {
			i__1 = leny;
			for (i__ = 1; i__ <= i__1; ++i__) {
				i__2 = iy;
				y[i__2] = std::complex<Type>(Type(0.0));
				iy += incy;
/* L30: */
			}
		} else {
			i__1 = leny;
			for (i__ = 1; i__ <= i__1; ++i__) {
				i__2 = iy;
				i__3 = iy;
				q__1 = beta * y[i__3];
				y[i__2] = q__1;
				iy += incy;
/* L40: */
				}
			}
		}
	}

	if (std::abs(alpha) == 0.0) {
		return 0;
	}

	if (lsame(trans, "N")) {

/*        Form  y := alpha*A*x + y. */

		jx = kx;
		if (incy == 1) {
			i__1 = n;
			for (j = 1; j <= i__1; ++j) {
				i__2 = jx;
				if (std::abs(x[i__2]) != 0.0) {
					i__2 = jx;
					q__1 = alpha * x[i__2];
					temp = q__1;
					i__2 = m;
					for (i__ = 1; i__ <= i__2; ++i__) {
						i__3 = i__;
						i__4 = i__;
						i__5 = i__ + j * lda;
						q__2 = temp * a[i__5];
						q__1 = y[i__4] + q__2;
						y[i__3] = q__1;
/* L50: */
					}
				}
				jx += incx;
/* L60: */
			}
		} else {
			i__1 = n;
			for (j = 1; j <= i__1; ++j) {
				i__2 = jx;
				if (std::abs(x[i__2]) != 0.0) {
					i__2 = jx;
					q__1 = alpha * x[i__2];
					temp = q__1;
					iy = ky;
					i__2 = m;
					for (i__ = 1; i__ <= i__2; ++i__) {
						i__3 = iy;
						i__4 = iy;
						i__5 = i__ + j * lda;
						q__2 = temp * a[i__5];
						q__1 = y[i__4] + q__2;
						y[i__3] = q__1;
						iy += incy;
/* L70: */
					}
				}
				jx += incx;
/* L80: */
			}
		}
	} else {

/*        Form  y := alpha*A**T*x + y  or  y := alpha*A**H*x + y. */

		jy = ky;
		if (incx == 1) {
			i__1 = n;
			for (j = 1; j <= i__1; ++j) {
				temp = std::complex<Type>(Type(0.0));
				if (noconj) {
					i__2 = m;
					for (i__ = 1; i__ <= i__2; ++i__) {
						i__3 = i__ + j * lda;
						i__4 = i__;
						q__2 = a[i__3] * x[i__4];
						q__1 = temp + q__2;
						temp = q__1;
/* L90: */
					}
				} else {
					i__2 = m;
					for (i__ = 1; i__ <= i__2; ++i__) {
						conjugate(&q__3, &a[i__ + j * lda]);
						i__3 = i__;
						q__2 = q__3 * x[i__3];
						q__1 = temp + q__2;
						temp = q__1;
/* L100: */
					}
				}
				i__2 = jy;
				i__3 = jy;
				q__2 = alpha * temp;
				q__1 = y[i__3] + q__2;
				y[i__2] = q__1;
				jy += incy;
/* L110: */
			}
		} else {
			i__1 = n;
			for (j = 1; j <= i__1; ++j) {
				temp = std::complex<Type>(Type(0.0));
				ix = kx;
				if (noconj) {
					i__2 = m;
					for (i__ = 1; i__ <= i__2; ++i__) {
						i__3 = i__ + j * lda;
						i__4 = ix;
						q__2 = a[i__3] * x[i__4];
						q__1 = temp + q__2;
						temp = q__1;
						ix += incx;
/* L120: */
					}
				} else {
					i__2 = m;
					for (i__ = 1; i__ <= i__2; ++i__) {
						conjugate(&q__3, &a[i__ + j * lda]);
						i__3 = ix;
						q__2 = q__3 * x[i__3];
						q__1 = temp + q__2;
						temp = q__1;
						ix += incx;
/* L130: */
					}
				}
				i__2 = jy;
				i__3 = jy;
				q__2 = alpha * temp;
				q__1 = y[i__3] + q__2;
				y[i__2] = q__1;
				jy += incy;
/* L140: */
			}
		}
	}

    return 0;

} /* blas_cgemv_ */


// specializations
int
blas_sgemm(char const* transa, char const* transb, int m, int n, int k,
		double alpha, double const* a, int lda, double const* b, int ldb,
		double beta, double* c__, int ldc);

int
blas_cgemm (char const* transa, char const* transb, int m, int n, int k,
		std::complex<double> alpha, std::complex<double> const* a, int lda,
		std::complex<double> const* b, int ldb, std::complex<double> beta,
		std::complex<double> *c, int ldc);

template<>
int
blas_cgemv<double> (char const* trans, int m, int n, std::complex<double> alpha,
		std::complex<double> const* a, int lda, std::complex<double> const* x, int incx,
		std::complex<double> beta, std::complex<double> *y, int incy);


double clobber_constant();

// Declarations for linking with Fortran library libblas.so
extern "C" {
// level 1 BLAS
int idamax_(const int* n, double* x, const int* incx);
void drot_(const int* n, double* x, const int* incx,
		double* y, const int* incy, const double* c, const double* s);
void dcopy_(int const* n, double const* x, int const* incx,
		double *y, int const* incy);
void dswap_(int *n, double *x, int *incx, double *y, int *incy);
void dscal_(int *n, const double *t, double *x, int *incx);
void daxpy_(int *n, const double *t, const double *x, int *incx,
		double *y, int *incy);
double ddot_(int *n, const double *x, int *incx, const double *y, int *incy);
double dnrm2_(int *n, const double *x, int *incx);
double dznrm2_(int *n, const std::complex<double> *x, int *incx);
// level 2 blas
double dtrsv_(const char* uplo, const char* trans, const char* diag,
		int *n, const double* a, int *lda, double* x, int* incx);
void dlartg_(double* f, double* g, double* cs, double* sn, double* r);
void dgemv_(const char* trans, int const* m, int const* n,
		double const* alpha, double const* a, int const* lda,
		double const* x, int const* incx,
		double const* beta, double const* y, int const* incy);
void zgemv_(const char* trans, int const* m, int const* n,
		std::complex<double> const* alpha, std::complex<double> const* a,
		int const* lda, std::complex<double> const* x, int const* incx,
		std::complex<double> const* beta,
		std::complex<double> const* y, int const* incy);
void dger_(const int* m, const int* n, const double* alpha,
		const double* x, const int* incx, const double* y, const int* incy,
		double* a, const int* lda);
// level 3 blas
void dgemm_(const char* transa, const char* transb, int const* m,
		int const* n, int const* k, double const* alpha, double const* a,
		int const* lda, double const* b, int const* ldb,
		double const* beta, double* c, int const* ldc);
void zgemm_(const char* transa, const char* transb, int const* m,
		int const* n, int const* k, std::complex<double> const* alpha,
		std::complex<double> const* a, int const* lda,
		std::complex<double> const* b, int const* ldb,
		std::complex<double> const* beta, std::complex<double>* c, int const* ldc);
} // extern "C"

namespace blas {

// Mixed-mode blas
template<typename Dtx, typename Dty, typename Dtz>
void
dot(int n, Dtx const* x, int incx, Dty const* y, int incy, Dtz& c) noexcept {
	c = Dtz();
	int ix{0};
	int iy{0};
	for (int i=0; i<n; i++) {
		c += x[ix]*y[iy];
		ix += incx;
		iy += incy;
	}
}

#ifdef NEVER // no Ad specializations
// specializations
// Ad-Ad
template<>
void
dot(int n, Ad const* x, int incx, Ad const* y, int incy, Ad& z) noexcept ;
// double-Ad
template<>
void
dot(int n, double const* x, int incx, Ad const* y, int incy, Ad& z) noexcept ;
// Ad-double
template<>
void
dot(int n, Ad const* x, int incx, double const* y, int incy, Ad& z) noexcept ;
// double-complex<Ad>
template<>
void
dot(int n, double const* x, int incx, std::complex<Ad> const* y,
		int incy, std::complex<Ad>& z) noexcept ;
#endif // NEVER // no Ad specializations

// scal
template<typename Dta, typename Dtx>
void
scal (int n, Dta const& a, Dtx* x, int incx)  noexcept {
// x <= a*x
	int ix = 0;
	for (int i=0; i<n; i++) {
		x[ix] *= a;
		ix += incx;
	}
}
// specializations
// complex<Ad>, complex<Ad>
template<>
void
scal(int n, std::complex<Ad> const& a, std::complex<Ad>* x, int incx) noexcept ;

// axpy
template<typename Dta, typename Dtx, typename Dty>
void
axpy(int n, Dta const& a, Dtx const* x, int incx, Dty* y, int incy)  noexcept {
// y = a*x + y
	int ix{0};
	int iy{0};
	for (int i=0; i<n; i++) {
		y[iy] += a*x[ix];
		iy += incy;
		ix += incx;
	}
}
// specializations
template<>
void
axpy(int n, std::complex<Ad> const& a, std::complex<Ad> const* x, int incx,
		std::complex<Ad>* y, int incy) noexcept ;

// gemv: matrix-vector multiplication
template<typename Dtalpha, typename Dta, typename Dtx, typename Dtbeta, typename Dty>
void
gemv(std::string const& trans, int m, int n, Dtalpha const& alpha, Dta const* a, int lda,
		Dtx const* x, int incx, Dtbeta const& beta, Dty* y, int incy) {
// y = alpha*A*x + beta*y
// A:  (m,n)
// x:  n
// y:  m
	if (trans != "n" && trans != "N")
		throw std::runtime_error("transpose not implemented in blas::gemv");
	int iy{0};
	for (int i=0; i<m; i++) {
		y[iy] *= beta;
		int ix{0};
		for (int j=0; j<n; j++) {
			y[iy] += alpha*a[i+j*lda]*x[ix];
			ix += incx;
		}
		iy += incy;
	}
}
// specializations
// gemv (double, complex<Ad>, complex<Ad>, double, complex<Ad>)
template<>
void
gemv(std::string const& trans, int m, int n, double const& alpha, std::complex<Ad> const* a,
		int lda, std::complex<Ad> const* x, int incx, double const& beta,
		std::complex<Ad>* y, int incy);
// double-Ad gemv
template<>
void
gemv(std::string const& trans, int m, int n, double const& alpha, double const* a,
		int lda, Ad const* x, int incx, double const& beta, Ad* y, int incy);

// blas::gemm: 
template<typename Dtalpha, typename Dta, typename Dtb, typename Dtbeta, typename Dtc>
void
gemm(std::string const& transa, std::string const& transb, int m, int n, int p,
	Dtalpha const& alpha, Dta const* a, int lda,
		Dtb const* b, int ldb, Dtbeta const& beta, Dtc* c, int ldc) noexcept {
// c = alpha*A*B + beta*C
// A:  (m,p)
// B:  (p,n)
// y:  m
	for (int i=0; i<m; i++) {
		for (int j=0; j<n; j++) {
			Dtc cij{0.0};
			for (int k=0; k<p; k++)
				cij += a[i+k*lda]*b[k+j*ldb];
			c[i+j*ldc] = cij;
		}
	}
}

// Copy elements from 2D array A to 2D array B based on row and
// column keys: B = beta*B + alpha*A
// kind of like a generalized version of saxpy
// 4 different datatypes are possible:
//    Rkey   row keys for A and B
//    Ckey   column keys for A and B
//    Dta    A datatype
//    Dtb    B datatype with Dtb(Dta) constructor
// There are 2 versions: one takes Dta* and Dtb*, the other
// takes vector<Dta> and vector<Dtb>. The vector version does a
// bit more checking.
// Both A and B arrays must be in column major format with leading
// dimensions the number of row keys in each. The number of columns
// must be the number of column keys.

template<typename Rkey, typename Ckey, typename Dta, typename Dtb>
void
copy (const std::vector<Rkey>& rkeya, const std::vector<Ckey>& ckeya,
	Dta alpha, const Dta* a, const std::vector<Rkey>& rkeyb,
	const std::vector<Ckey>& ckeyb, Dtb beta, Dtb* b) {
// insert rectangular matrix A into rectangular matrix B using row & column keys
// pointer version
	size_t lda = rkeya.size();
	size_t ldb = rkeyb.size();

	auto rbegin{rkeyb.begin()};
	auto rend{rkeyb.end()};
	auto cbegin{ckeyb.begin()};
	auto cend{ckeyb.end()};
	for (size_t ja=0; ja<ckeya.size(); ja++) {
		auto idx = find(cbegin, cend, ckeya[ja]);
		// ok if a ckey of A is not in B
		if (idx == cend) continue;
		int jb = idx - cbegin;
		for (size_t ia=0; ia<rkeya.size(); ia++) {
			auto idx = find(rbegin, rend, rkeya[ia]);
			// ok if a rkey of A is not in B
			if (idx == rend) continue;
			int ib = idx - rbegin;
			b[ib+jb*ldb] = beta*b[ib+jb*ldb] + alpha*a[ia+ja*lda];
		}
	}
}

template<typename Rkey, typename Ckey, typename Dta, typename Dtb>
void
copy (const std::vector<Rkey>& rkeya, const std::vector<Ckey>& ckeya,
	Dta alpha, const std::vector<Dta>& a, const std::vector<Rkey>& rkeyb,
	const std::vector<Ckey>& ckeyb, Dtb beta, std::vector<Dtb>& b) {
// insert rectangular matrix A into rectangular matrix B using row & column keys
// vector version
	size_t lda = rkeya.size();
	size_t nca = ckeya.size();
	size_t ldb = rkeyb.size();
	size_t ncb = ckeyb.size();
	size_t na = a.size();
	size_t nb = b.size();
	// just check sizes...
	if (lda*nca != a.size()) {
		std::cerr << "size mismatch for A: " << lda << " rkeys, " << nca
				<< " ckeys, A.size = " << na << std::endl;
		throw std::runtime_error("blas::copy: size mismatch for A");
	}
	if (ldb*ncb != b.size()) {
		std::cerr << "size mismatch for B: " << ldb << " rkeys, " << ncb
				<< " ckeys, B.size = " << nb << std::endl;
		throw std::runtime_error("blas::copy: size mismatch for B");
	}
	// ... and call the pointer version
	return blas::copy(rkeya, ckeya, alpha, &a[0], rkeyb, ckeyb, beta, &b[0]);
}

//------------------------------------------------------------------
// Various non-BLAS linear algebraic functions

template<typename Dta, typename Dtb>
void
append (const std::vector<Dta>& a, std::vector<Dtb>& b) {
// append na-vector "a" to (na,m) matrix "b" store column-major in
// a vector with na*m elements on entry, na*(m+1) elements on return
	size_t na = a.size();
	size_t m = b.size()/na;
	if (m*na != b.size())
		throw std::runtime_error(vastr("blas::append called with ",na,"-vector a, ",
					b.size(),"-vector b"));

	b.resize(na*(m+1));
	// copy item-wise in case Dta != Dtb
	Dtb* bp = &b[na*m];
	for (size_t i=0; i<na; i++)
		*bp++ = a[i];
}

// append

// create the real (2m,2n) representation of a complex (m,n) matrix
void
real_rep (int m, int n, std::complex<double>* c, int ldc, double* r, int ldr);


// Various normalization schemes XXX some overlap here?
template<typename Type>
double
normalize_inf(size_t n, Type* a, size_t incx, bool normalize=false, double sf=1.0) {
// Compute the infinity norm of a vector and optionally
// scale the vector so that the infinity norm is "sf"
	double rval{0.0};
	Type largest{0.0};
	for (size_t i=0; i<n; i++) {
		double ai = std::abs(a[i]);
		if (ai > rval) {
			rval = ai;
			largest = a[i];
		}
	}
	if (normalize) {
		if (rval > 0.0) {
			Type t = sf/largest;
			for (size_t i=0; i<n; i++)
				a[i] *= t;
		}
	}
	return rval;
}

template<typename Type>
double
normalize2(size_t n, Type* a, size_t incx, bool normalize=false, double sf=1.0) {
// Compute the 2-norm of a vector and optionally
// scale the vector so that the infinity norm is "sf"
	double rval = blas_snrm2(n, a, incx);  // XXX only real datatypes
	if (normalize) {
		if (rval > 0.0) {
			double t = sf/rval;
			// cant use blas_scal: "a" may not be double
			for (size_t i=0; i<n; i++) {
				a[i] *= t;
			}
		}
	}
	return rval;
}

template<typename Type>
void
scatter(size_t n, Type const* x, int const* indx, Type* y) {
/*------------------------------------------------------------------
 * C/C++ version of ssctr_()
 * Scatter (copy) the components of a sparse vector stored in
 * compressed form into specified components of a vector in full-
 * storage form. If y is initially the zero vector, then ssctr()
 * performs
 *		y := x,
 * which is represented by the following FORTRAN:
 *		Y(INDX(I)) = X(I)   for  I = 1, ..., NZ
 *
 * This routine is a C implementation of the sparse BLAS routine
 * described in
 * Dodson, D.S., and Grimes, R.G., "Sparse Extensions to the FORTRAN
 * Basic Linear Algebra Subprograms", ACM TOMS V.17 No. 2, June 1991,
 * pp. 253-263
 *
 * Note: this routine should be replaced with the blas routine
 * of the same name and same interface
 *------------------------------------------------------------------*/

	while (n--)
		y[*indx++ - 1] = *x++;
}

} // namespace blas

#endif // blas_h
