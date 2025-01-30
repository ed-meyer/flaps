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

#include "Ad.h"	// for specializations

// Templated implementations of BLAS functions for mixed-modes with
// specializations in blas.cpp
// Some additional linear-algebraic functions:
//   real_rep:  create the real (2m,2n) representation of a complex (m,n) matrix

double clobber_constant();

// Declarations for linking with Fortran library libblas.so
// used in specializations of BLAS
extern "C" {
// level 1 BLAS
int idamax_(const int* n, const double* x, const int* incx);
int izamax_(const int* n, const std::complex<double>* x, const int* incx);
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

// #ifdef NEVER // try templated
double snrm2(int n, const double* x, int incx);
double scnrm2(int n, const std::complex<double>* x, int incx);
int isamax(int n, const double* x, int incx);
int icamax(int n, const std::complex<double>* x, int incx);
// #else // NEVER // try templated

// Mixed-mode blas
template<typename Dtx, typename Dty>
void
copy(int n, Dtx const* x, int incx, Dty* y, int incy) noexcept {
	int ix{0};
	int iy{0};
	for (int i=0; i<n; i++) {
		y[iy] = x[ix];
		ix += incx;
		iy += incy;
	}
}

template<typename Dtx, typename Dty, typename Dtz>
void
dot(int n, Dtx const* x, int incx, Dty const* y, int incy, Dtz& c) noexcept {
	c = Dtz(0.0);
	int ix{0};
	int iy{0};
	for (int i=0; i<n; i++) {
		c += x[ix]*y[iy];
		ix += incx;
		iy += incy;
	}
}
// specializations
template<>
void
dot(int n, const double* x, int incx, const Ad* y, int incy, Ad& c) noexcept;

template<typename Dtx, typename Dty, typename Dtz>
void
dotc(int n, const std::complex<Dtx>* x, int incx, Dty const* y, int incy, Dtz& c) noexcept {
	c = Dtz(0);	// initialize to zero
	int ix{0};
	int iy{0};
	for (int i=0; i<n; i++) {
		std::complex<Dtx> xc(x[ix].real(), -x[ix].imag());
		c += xc*y[iy];
		ix += incx;
		iy += incy;
	}
}


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

// axpy: y += a*x
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
// specialization: complex<Ad> += complex<Ad>*complex<Ad>
template<>
void
axpy(int n, const std::complex<Ad>& a, const std::complex<Ad>* x, int incx,
			std::complex<Ad>* y, int incy) noexcept;

//------------------------------------------------------------------
// gemv: matrix-vector multiplication
template<typename Dtalpha, typename Dta, typename Dtx, typename Dtbeta, typename Dty>
void
gemvx(std::string const& trans, int m, int n, Dtalpha const& alpha, Dta const* a, int lda,
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

	blas::scal(m, beta, y, incy);
	const Dta* ap = a;
	for (int j=0; j<n; j++) {
		blas::axpy(m, alpha*x[j], ap, 1, y, incy);
		ap += lda;
	}
}

// specializations
// gemv <double, double, double, double, double>
template<>
void
gemv(const std::string& trans, int m, int n, const double& alpha, const double* a,
		int lda, const double* x, int incx, const double& beta, double* y, int incy);

// gemv (complex<double>)
template<>
void
gemv(const std::string& trans, int m, int n, const std::complex<double>& alpha,
		const std::complex<double>* a, int lda,
		const std::complex<double>* x, int incx,
		const std::complex<double>& beta,
		std::complex<double>* y, int incy);

//------------------------------------------------------------------
// gemm XXX transposed matrices not implemented
template<typename Dtalpha, typename Dta, typename Dtb,
	typename Dtbeta, typename Dtc>
void
gemm(std::string const& transa, std::string const& transb, int m, int n, int p,
	Dtalpha const& alpha, Dta const* a, int lda,
		Dtb const* b, int ldb, Dtbeta const& beta, Dtc* c, int ldc) {
// c = alpha*A*B + beta*C
// A:  (m,p)
// B:  (p,n)
// y:  m
	// transposed matrices not implemented
	if (tolower(transa[0]) != 'n' || tolower(transb[0]) != 'n')
			std::cerr << "transposed matrices not implemented in blas::gemm\n";

	for (int i=0; i<m; i++) {
		for (int j=0; j<n; j++) {
			Dtc cij{0.0};
			for (int k=0; k<p; k++)
				cij += a[i+k*lda]*b[k+j*ldb];
			c[i+j*ldc] = cij;
		}
	}
}

// double specialization
template<>
void
gemm(std::string const& transa, std::string const& transb, int m, int n, int p,
	double const& alpha, double const* a, int lda,
		double const* b, int ldb, double const& beta, double* c, int ldc);

// complex<double> specialization
template<>
void
gemm(std::string const& transa, std::string const& transb, int m, int n, int p,
	std::complex<double> const& alpha, std::complex<double> const* a, int lda,
		std::complex<double> const* b, int ldb, std::complex<double> const& beta,
		std::complex<double>* c, int ldc);

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

template<typename Dtx>
double
nrm2(int n, const Dtx* x, int incx) {
	double rval;
	blas::dot(n,x,incx,x,incx,rval);
	return sqrt(rval);
}
template<typename Dtx>
double
nrm2(int n, const std::complex<Dtx>* x, int incx) {
	std::complex<Dtx> rval;
	blas::dotc(n,x,incx,x,incx,rval);
	return sqrt(rval.real());
}

template<typename Type>
double
normalize2(size_t n, Type* a, size_t incx, bool normalize=false, double sf=1.0) {
// Compute the 2-norm of a vector and optionally
// scale the vector so that the infinity norm is "sf"
	double rval = nrm2(n, a, incx);;
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
