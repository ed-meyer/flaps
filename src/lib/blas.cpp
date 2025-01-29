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
#include "Ad.h"
#include "blas.h"
#include "message.h"
#include "text.h"
#include "trace.h"
#include "vastr.h"

using namespace std;


// Specializations of templatized versions of the BLAS in blas.h

template<>
double
blas_snrm2 (int n, double const* a, int inca) { return dnrm2_(&n, a, &inca); }

template<>
double
blas_scnrm2 (int n, complex<double> const* a, int inca) { return dznrm2_(&n, a, &inca); }

template<>
int
blas_sgemv<double> (char const* trans, int m, int n, double alpha, double const* a,
		int lda, double const* x, int incx, double beta, double *y, int incy) {
// specialization: sgemv(double)
	dgemv_ (trans, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);
	return 0;
}

template<>
int
blas_cgemv<double> (char const* trans, int m, int n, complex<double> alpha,
		complex<double> const* a, int lda, complex<double> const* x, int incx,
		complex<double> beta, complex<double> *y, int incy) {
// specialization: cgemv(complex<double>)
	zgemv_ (trans, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy);
	return 0;
}

int
blas_sgemm (char const* transa, char const* transb, int m, int n, int k,
		double alpha, double const* a, int lda, double const* b, int ldb,
		double beta, double *c, int ldc) {
	dgemm_(transa, transb, &m, &n, &k, &alpha, a, &lda, b, &ldb,
			&beta, c, &ldc);
	return 0;
}

int
blas_cgemm (char const* transa, char const* transb, int m, int n, int k,
		complex<double> alpha, complex<double> const* a, int lda,
		complex<double> const* b, int ldb, complex<double> beta,
		complex<double> *c, int ldc) {
	zgemm_(transa, transb, &m, &n, &k, &alpha, a, &lda, b, &ldb,
			&beta, c, &ldc);
	return 0;
}

namespace blas {

#ifdef NEVER // no specializations
// specialization: dot (Ad,Ad)
template<>
void
dot(int n, Ad const* x, int incx, Ad const* y, int incy, Ad& z) noexcept {
// specialization: Ad = Ad * Ad
	// initialize: zero z
	int nd = Ad::ndata();
	double* zd = z.data();
	for (int i=0; i<nd; i++)
		zd[i] = 0.0;
		
	int ix{0};
	int iy{0};
	for (int i=0; i<n; i++) {
		double const* xd = x[ix].data();
		double const* yd = y[iy].data();
		zd[0] += xd[0]*yd[0];
		// derivatives
		for (int j=1; j<nd; j++)
			zd[j] += xd[0]*yd[j] + xd[j]*yd[0];
		ix += incx;
		iy += incy;
	}
}

// specialization: dot (double,Ad)
template<>
void
dot(int n, double const* x, int incx, Ad const* y, int incy, Ad& z) noexcept {
// specialization: Ad = double * Ad
	// initialize: zero z
	int nd = Ad::ndata();
	double* zd = z.data();
	for (int i=0; i<nd; i++)
		zd[i] = 0.0;
		
	int ix{0};
	int iy{0};
	for (int i=0; i<n; i++) {
		double const* yd = y[iy].data();
		double xi = x[ix];
		for (int j=0; j<nd; j++)
			zd[j] += xi*yd[j];
		ix += incx;
		iy += incy;
	}
}

// specialization dot (Ad,double)
template<>
void
dot(int n, Ad const* x, int incx, double const* y, int incy, Ad& z) noexcept {
// commutative
	return blas::dot(n,y,incy,x,incx,z);
}

// specialization dot (double,complex<Ad>)
template<>
void
dot(int n, double const* x, int incx, complex<Ad> const* y, int incy, complex<Ad>& z) noexcept {
// specialization: Ad = double * complex<Ad>
	// initialize: zero z
	int nd = Ad::ndata();
	double* zrd = Ad::real(z).data();
	double* zid = Ad::imag(z).data();
	for (int i=0; i<nd; i++) {
		zrd[i] = 0.0;
		zid[i] = 0.0;
	}
		
	int ix{0};
	int iy{0};
	for (int i=0; i<n; i++) {
		double const* yrd = Ad::real(y[iy]).data();
		double const* yid = Ad::imag(y[iy]).data();
		double xi = x[ix];
		for (int j=0; j<nd; j++) {
			zrd[j] += xi*yrd[j];
			zid[j] += xi*yid[j];
		}
		ix += incx;
		iy += incy;
	}
}
#endif // NEVER // no specializations

// scal (complex<Ad>,complex<Ad>)
template<>
void
scal(int n, complex<Ad> const& a, complex<Ad>* x, int incx) noexcept {
// Reduced Ad copy constructor version of mixed-mode scal
	int nd = Ad::ndata();
	int ix{0};
	double const* ard = Ad::real(a).data();
	double const* aid = Ad::imag(a).data();
	double arv = ard[0];
	double aiv = aid[0];
	for (int i=0; i<n; i++) {
		double* xrd = Ad::real(x[ix]).data();
		double* xid = Ad::imag(x[ix]).data();
		double xrv = xrd[0];
		double xiv = xid[0];
		// x[i] := a*x[i] = (ar*xr - ai*xi, ai*xr + ar*xi)
		//    x[i].real = ar*xr - ai*xi
		//    x[i].imag = ai*xr + ar*xi
		// derivatives:
		for (int j=1; j<nd; j++) {
			xrd[j] = (arv*xrd[j] + ard[j]*xrv) - (aiv*xid[j] + aid[j]*xiv);
			xid[j] = (aiv*xrd[j] + aid[j]*xrv) + (arv*xid[j] + ard[j]*xiv);
		}
		xrd[0] = arv*xrv - aiv*xiv;
		xid[0] = aiv*xrv + arv*xiv;
		ix += incx;
	}
}

// axpy (complex<Ad>,complex<Ad>)
template<>
void
axpy(int n, complex<Ad> const& a, complex<Ad> const* x, int incx,
		complex<Ad>* y, int incy) noexcept {
// Reduced Ad copy constructor version of mixed-mode axpy
	int nd = Ad::ndata();
	int ix{0};
	int iy{0};
	double const* ard = Ad::real(a).data();
	double const* aid = Ad::imag(a).data();
	double arv = ard[0];
	double aiv = aid[0];
	for (int i=0; i<n; i++) {
		double const* xrd = Ad::real(x[ix]).data();
		double const* xid = Ad::imag(x[ix]).data();
		double* yrd = Ad::real(y[iy]).data();
		double* yid = Ad::imag(y[iy]).data();
		double xrv = xrd[0];
		double xiv = xid[0];
		// y[i] += a*x[i] = (ar*xr - ai*xi, ai*xr + ar*xi)
		//    y[i].real = ar*xr - ai*xi
		//    y[i].imag = ai*xr + ar*xi
		for (int j=1; j<nd; j++) {
			yrd[j] += (arv*xrd[j] + ard[j]*xrv) - (aiv*xid[j] + aid[j]*xiv);
			yid[j] += (aiv*xrd[j] + aid[j]*xrv) + (arv*xid[j] + ard[j]*xiv);
		}
		yrd[0] += arv*xrv - aiv*xiv;
		yid[0] += aiv*xrv + arv*xiv;
		ix += incx;
		iy += incy;
	}
}

// cgemv (complex<Ad>, complex<Ad>)
template<>
void
gemv(string const& trans, int m, int n, double const& alpha, complex<Ad> const* a,
		int lda, complex<Ad> const* x, int incx, double const& beta,
		complex<Ad>* y, int incy) {
// specialization: gemv(double, complex<Ad>, complex<Ad>, double, complex<Ad>)
// XXX alpha must be 1, beta must be zero
// Reduced Ad copy constructor version of mixed-mode gemv
	// limited version: no transpose, alpha, beta doubles 1 & 0
	if (trans != "n" && trans != "N")
		flaps::error ("blas::gemv: transpose not implemented");
		// throw runtime_error("blas::gemv: transpose not implemented");
	if (alpha != 1.0 || beta != 0.0)
		flaps::error ("blas::gemv: alpha!=1, beta!=0 not implemented");
		// throw runtime_error("blas::gemv: alpha!=1, beta!=0 not implemented");

	int nd = Ad::ndata();
	int iy{0};
	for (int i=0; i<m; i++) {
		double* yrd = Ad::real(y[iy]).data();
		double* yid = Ad::imag(y[iy]).data();
		// initialize y[i]
		for (int k=0; k<nd; k++) {
			yrd[k] = 0.0;
			yid[k] = 0.0;
		}
		int ix{0};
		for (int j=0; j<n; j++) {
			double const* ard = Ad::real(a[i+j*lda]).data();
			double const* aid = Ad::imag(a[i+j*lda]).data();
			double arv = ard[0];
			double aiv = aid[0];
			double const* xrd = Ad::real(x[ix]).data();
			double const* xid = Ad::imag(x[ix]).data();
			double xrv = xrd[0];
			double xiv = xid[0];
			// y[i] += a[i,j]*x[j] = (ar*xr - ai*xi, ai*xr + ar*xi)
			//    y[i].real = ar*xr - ai*xi
			//    y[i].imag = ai*xr + ar*xi
			// derivatives:
			for (int k=1; k<nd; k++) {
				yrd[k] += (arv*xrd[k] + ard[k]*xrv) - (aiv*xid[k] + aid[k]*xiv);
				yid[k] += (aiv*xrd[k] + aid[k]*xrv) + (arv*xid[k] + ard[k]*xiv);
			}
			// values:
			yrd[0] += arv*xrv - aiv*xiv;
			yid[0] += aiv*xrv + arv*xiv;
			ix += incx;
		}
		iy += incy;
	}
}

// sgemv (double,Ad)
template<>
void
gemv(string const& trans, int m, int n, double const& alpha, double const* a,
		int lda, Ad const* x, int incx, double const& beta, Ad* y, int incy) {
// zero Ad copy constructor version of mixed-mode gemv
//   y = alpha*A*x + beta*y
// A: (m,n) double in a (lda,n) array
// x: (n) Ad
// y: (m) Ad
	// limited version: no transpose, alpha, beta doubles 1 & 0
	if (trans != "n" && trans != "N")
		flaps::error ("blas::gemv: transpose not implemented");
		// throw runtime_error("blas::gemv: transpose not implemented");
	if (alpha != 1.0 || beta != 0.0)
		flaps::error("blas::gemv: alpha!=1, beta!=0 not implemented");
		// throw runtime_error("blas::gemv: alpha!=1, beta!=0 not implemented");

	int nd = Ad::ndata();
	int iy{0};
	for (int i=0; i<m; i++) {
		double* yd = y[iy].data();
		// initialize y[i]
		for (int k=0; k<nd; k++)
			yd[k] = 0.0;

		int ix{0};
		for (int j=0; j<n; j++) {
			double aij = a[i+j*lda];
			double const* xd = x[ix].data();
			// y[i] += a[i,j]*x[j]
			for (int k=0; k<nd; k++)
				yd[k] += aij*xd[k];
			ix += incx;
		}
		iy += incy;
	}
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
#include "tyme.h"

void
init() {
	int nder = 2;  // number of Ad derivative parameters
	// set nder parameter names as Ad parameters, then
	// Ad constructors will use the number of parameters
	// to allocate space for derivatives
	vector<string> mypar;
	for (int i=0; i<nder; i++) {
		mypar.push_back(vastr("par", i+1));
	}
	Ad::initialize(mypar);
	cerr << Ad::nder() << " automatic differentiation parameters: "
		<< Ad::toString() << endl;
}

void
test_dot(int n) {
	Trace trc(1,"test_dot");

	vector<complex<Ad>> x(n);
	vector<complex<Ad>> y(n);
	Ad xv{2.0};
	xv.deriv(0,1.0);
	xv.deriv(1,0.0);
	Ad yv{2.0};
	yv.deriv(0,1.0);
	yv.deriv(1,0.0);
	complex<Ad> cxv{xv, xv};
	complex<Ad> cyv{yv, yv};
	for (int i=0; i<n; i++) {
		x[i] = cxv;
		y[i] = cyv;
	}
	trc.dprint("testing blas::dot(complex<Ad>) with x = ", x);
	trc.dprint("testing blas::dot(complex<Ad>) with y = ", y);

	complex<Ad> result;
	{
		Tyme t("blas::dot");
		blas::dot(n, &x[0], 1, &y[0], 1, result);
		cerr << "blas::dot result: " << result << endl;
	}

	{
		Tyme t("blas_dot");
		blas_dot(n, &x[0], 1, &y[0], 1, result);
		cerr << "blas_dot result: " << result << endl;
	}
}

void
test_axpy(int n) {

	vector<complex<Ad>> x(n);
	vector<complex<Ad>> y(n);
	complex<Ad> a{Ad(2.0)};
	Ad::real(a).deriv(0,1.0);	// a is par1
	Ad::real(a).deriv(1,0.0);
	complex<Ad> p{Ad(3.0), Ad(4.0)};
	for (int i=0; i<n; i++) {
		x[i] = p;
	}
	cerr << "testing blas::axpy(complex<Ad>) with a = " << a << endl;
	cerr << " x = " << x << endl;
	cerr << " y = " << y << endl;

	blas::axpy(n, a, &x[0], 1, &y[0], 1);
	cerr << "result: " << y << endl;

	complex<Ad> zero{Ad(0.0)};
	for (int i=0; i<n; i++)
		y[i] = zero;
	blas_axpy(n, a, &x[0], 1, &y[0], 1);
	cerr << "blas_axpy result: " << y << endl;
}

int
main(int argc, char** argv) {
	Trace trc(1,"blas testing");
	bool axpy{false};
	bool dot{true};

	int n{3};	// default n
	if (argc > 1)
		n = stoi(argv[1]);

	// initialize Ad
	init();

	if (axpy)
		test_axpy(n);

	if (dot)
		test_dot(n);
}
#endif // MAIN
