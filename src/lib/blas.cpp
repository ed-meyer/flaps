//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

#include "config.h"
#include "blas.h"
#include "message.h"
#include "text.h"
#include "trace.h"
#include "vastr.h"

using namespace std;

namespace blas {

// dot specializations
template<>
void
dot(int n, const double* x, int incx, const Ad* y, int incy, Ad& c) noexcept {
// dot specialization for double*Ad (used by vspline)
	c.zero();
	int ix{0};
	int iy{0};
	for (int i=0; i<n; i++) {
		double xi = x[ix];
		c.data_[0] += xi*y[iy].data_[0];
		for (int j=1; j<Ad::ndata; j++)
			c.data_[j] += xi*y[iy].data_[j];
		ix += incx;
		iy += incy;
	}
}

// axpy specializations
template<>
void
axpy(int n, const complex<Ad>& a, const complex<Ad>* x, int incx,
	complex<Ad>* y, int incy) noexcept {
// for complex<Ad> += complex<Ad>*complex<Ad>
	complex<Ad> ax(0.0);
	int ix{0};
	int iy{0};
   Ad wk(0.0);
   for (int i=0; i<n; i++) {
      Ad::multcad(a, x[ix], ax, wk);
      y[iy] += ax;
		ix += incx;
		iy += incy;
   }
}


// nrm2 specializations
double
snrm2 (int n, double const* a, int inca) { return dnrm2_(&n, a, &inca); }
double
scnrm2 (int n, const complex<double>* a, int inca) { return dznrm2_(&n, a, &inca); }

// amax specializations
int
isamax(int n, const double* x, int incx) { return idamax_(&n, x, &incx); }
int
icamax(int n, const complex<double>* x, int incx) { return izamax_(&n, x, &incx); }

// gemv specializations
template<>
void
gemv(const string& trans, int m, int n, const double& alpha, const double* a,
		int lda, const double* x, int incx, const double& beta, double *y, int incy) {
// specialization: gemv(double) call libblas fcn
	dgemv_ (trans.c_str(), &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);
}

template<>
void
gemv(const string& trans, int m, int n, const complex<double>& alpha,
	const complex<double>* a, int lda,
	const complex<double>* x, int incx,
	const complex<double>& beta, complex<double>* y, int incy) {
// specialization: gemv(complex<double>) call libblas fcn
	zgemv_ (trans.c_str(), &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);
}

template<>
void
gemv(const string& trans, int m, int n, const double& alpha,
	const complex<Ad>* a, int lda, const complex<Ad>* x, int incx,
	const double& beta, complex<Ad>* y, int incy) {
// specialization for flut: alpha=1, beta=0, complex<Ad>
	if (trans != "n" && trans != "N")
		throw std::runtime_error("transpose not implemented in blas::gemv");
	if (alpha != 1.0 || beta != 0.0)
		throw std::runtime_error("special gemv only for alpha=1, beta=0");

	const complex<Ad>* ap = a; // column pointer
	int jx{0};
   for (int j=0; j<n; j++) {
      blas::axpy(m, x[jx], ap, 1, y, incy);
      ap += lda;
		jx += incx;
   }
}

// gemm double specialization: call libblas version (allows transposed matrices)
template<>
void
gemm(std::string const& transa, std::string const& transb, int m, int n, int p,
	double const& alpha, double const* a, int lda,
		double const* b, int ldb, double const& beta, double* c, int ldc) {
	dgemm_(transa.c_str(), transb.c_str(), &m, &n, &p, &alpha, a, &lda, b, &ldb,
			&beta, c, &ldc);
}

// gemm complex<double> specialization: call libblas version (allows transposed matrices)
template<>
void
gemm(std::string const& transa, std::string const& transb, int m, int n, int p,
	std::complex<double> const& alpha, std::complex<double> const* a, int lda,
		std::complex<double> const* b, int ldb, std::complex<double> const& beta,
		std::complex<double>* c, int ldc) {
	zgemm_(transa.c_str(), transb.c_str(), &m, &n, &p, &alpha, a, &lda, b, &ldb,
			&beta, c, &ldc);
}

void
real_rep (int m, int n, complex<double>* c, int ldc, double* r, int ldr) {
/*------------------------------------------------------------------
 * creates the real (2m,2n) representation of a complex (m,n) matrix
 *
 *    Each element of the complex matrix becomes a 2 by 2 block of
 *    the real matrix, and each 2 by 2 block has the form
 *           ! a   -b !
 *           ! b    a !
 *    where the complex element is complex(a,b)
 * input:
 *  m      number of rows in the complex matrix c
 *  n      number of columns in complex matrix c
 *  c      complex n by m matrix contained in a complex
 *         array dimensioned ldc by m in the calling routine
 *  ldc    row dimension of c in the calling routine
 *  ldr    row dimension of real matrix r in the calling routine
 *
 * output:
 *  r      (double) 2*m by 2*n matrix containing the real
 *         representation of complex matrix c.
 *------------------------------------------------------------------*/
	// treat C as a (2m,n) double matrix in a (2*ldc,n) double matrix
	int ldc2 = 2*ldc;
	for (int j=0; j<n; j++) {
		// assuming column-major (fortran) storage, set pointer to the
		// top of column j of C, columns 2*j and 2*j+1 of R
		double* r1 = r + 2*j*ldr;
		double* r2 = r1 + ldr;
		double* cp = (double*)c + j*ldc2;
		for (int i=0; i<m; i++) {
			double cr = *cp++;
			double ci = *cp++;
			*r1++ = cr;
			*r1++ = ci;
			*r2++ = -ci;
			*r2++ = cr;
		}
	}
}

} // namespace blas


#ifdef MAIN
#undef MAIN

#include "Ad.h"
#include "tyme.h"

void
init() {
	int nder = 2;  // number of Ad derivative parameters
	// set nder parameter names as Ad parameters, then
	// Ad constructors will use the number of parameters
	// to allocate space for derivatives
	T_(Trace trc(0,"init");)
	vector<string> mypar;
	for (int i=0; i<nder; i++) {
		mypar.push_back(vastr("par", i+1));
	}
	Ad::initialize(mypar);
}

void
test_dot(int n) {
	T_(Trace trc(0,"test_dot");)

	vector<complex<Ad>> x(n);
	vector<complex<Ad>> y(n);
	Ad xv{2.0};
	xv.der(0,1.0);
	xv.der(1,0.0);
	Ad yv{2.0};
	yv.der(0,1.0);
	yv.der(1,0.0);
	complex<Ad> cxv{xv, xv};
	complex<Ad> cyv{yv, yv};
	for (int i=0; i<n; i++) {
		x[i] = cxv;
		y[i] = cyv;
	}
	T_(trc.dprint("testing blas::dot(complex<Ad>) with x = ", x);)
	T_(trc.dprint("testing blas::dot(complex<Ad>) with y = ", y);)

	complex<Ad> result;
	{
		Tyme t("blas::dot");
		blas::dot(n, &x[0], 1, &y[0], 1, result);
		cerr << "blas::dot result: " << result << endl;
	}
}

void
test_axpy(int n) {
	T_(Trace trc(0,"test_axpy");)

	vector<complex<Ad>> x(n);
	vector<complex<Ad>> y(n);
	complex<Ad> a{Ad(2.0)};
	Ad::real(a).der(0,1.0);	// a is par1
	Ad::real(a).der(1,0.0);
	complex<Ad> p{Ad(3.0), Ad(4.0)};
	for (int i=0; i<n; i++) {
		x[i] = p;
	}
	T_(trc.dprint("testing blas::axpy(complex<Ad>) with a = ", a);)
	T_(trc.dprint(" x = ", x);)
	T_(trc.dprint(" y = ", y);)

	blas::axpy(n, a, &x[0], 1, &y[0], 1);
	T_(trc.dprint("result: ", y);)

	complex<Ad> zero{Ad(0.0)};
	for (int i=0; i<n; i++)
		y[i] = zero;
	blas::axpy(n, a, &x[0], 1, &y[0], 1);
	T_(trc.dprint("blas::axpy result: ", y);)
}

void
test_gemv(int n) {
	T_(Trace trc(0,"test_gemv");)

	vector<complex<Ad>> x(n);
	vector<complex<Ad>> y(n);
	vector<complex<Ad>> A(n*n);
	for (int i=0; i<n; i++)
		A[i+n*i] = 1.0;
	complex<Ad> alpha{Ad(2.0)};
	Ad::real(alpha).der(0,1.0);	// a is par1
	Ad::real(alpha).der(1,0.0);
	complex<Ad> p{Ad(3.0), Ad(4.0)};
	for (int i=0; i<n; i++) {
		x[i] = p;
	}
	double beta{3.0};

	T_(trc.dprint("testing blas::gemv(complex<Ad>) with A = ", A);)
	T_(trc.dprint(" alpha = ",alpha);)
	T_(trc.dprint(" x = ", x);)
	T_(trc.dprint(" beta = ",beta);)
	T_(trc.dprint(" y = ", y);)

	// save y for gemvx
	vector<complex<Ad>> ys = y;

	blas::gemv("n", n, n, alpha, A.data(), n, &x[0], 1, beta, &y[0], 1);
	T_(trc.dprint("gemv: ", y);)

	// re-do with gemvx
	blas::gemv("n", n, n, alpha, A.data(), n, &x[0], 1, beta, &ys[0], 1);
	T_(trc.dprint("gemvx: ", ys);)
}

int
main(int argc, char** argv) {
	char* dbg = getenv("DEBUG");
	if (dbg != nullptr) Trace::debug(std::stoi(dbg));
	T_(Trace trc(0,"blas testing");)
	// what to test:
	bool axpy{false};
	bool dot{false};
	bool gemv{true};

	int n{3};	// default n
	if (argc > 1)
		n = stoi(argv[1]);

	// initialize Ad
	init();

	if (axpy)
		test_axpy(n);

	if (dot)
		test_dot(n);

	if (gemv)
		test_gemv(n);
}
#endif // MAIN
