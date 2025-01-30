//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

// Cover routines for LAPACK subroutines to allow them
// to be called from C or C++.
// Fortran names differ with architectures; here we assume
// the fortran name is lower case with a trailing underscore
// There are 2 ways to call fortran subroutines, e.g.
//    lapack::dqr1up(m,n,k,Q,ldq,R,ldr,u,v)
//    dqr1up_(&m,&n,&k,Q,&ldq,R,&ldr,u,v)
// (assuming Q, R, u, and v are double*)
//
// Note: I decided against using lapacke (the C interface that
// now comes with lapack) for a few reasons:
//   - they require e.g. lapack_double instead of double
//   - needs LAPACKE_name
//   - takes an extra arg LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR but
//     I always use column major
//   - lapacke may not be on all machines


#include "config.h"
#include "blas.h"
#include "lapack.h"
#include "matrix.h"
#include "trace.h"

#define IJ(i,j,ld1) ((i) + (ld1)*(j))

using namespace std;

namespace lapack {

int
dposv(const char* uplo, int n, int nrhs, double* a, int lda,
	double* b, int ldb) {
	int info;
	dposv_(uplo, &n, &nrhs, a, &lda, b, &ldb, &info);
	return info;
}

int
dpotrf(const char* uplo, int n, double* a, int lda) {
// Cholesky factorization of a symmetric positive-definite matrix
	int info;
	dpotrf_(uplo, &n, a, &lda, &info);
	return info;
}

int
dpotrs(const char* uplo, int n, int nrhs, double* a,
	int lda, double* b, int ldb) {
// Solution of linear equations A*X = B using the factorization
// from dpotrf
	int info;
	dpotrs_(uplo, &n, &nrhs, a, &lda, b, &ldb, &info);
	return info;
}

int dgetrf (int m, int n, double* a, int lda, int* ipiv) {
	int info;
	dgetrf_(&m, &n, a, &lda, ipiv, &info);
	return info;
}

int
dgesv (int n, int nrhs, double* a, int lda, int* ipiv, double* b, int ldb) {
	int info;
	dgesv_(&n, &nrhs, a, &lda, ipiv, b, &ldb, &info);
	return info;
}

int
dgesvx (char const* fact, char const* trans, int n, int nrhs,
	double* a, int lda, double* af, int ldaf, int* ipiv, 
	char equed[2], double* r, double* c, double* b, int ldb,
	double* x, int ldx, double* rcond, double* ferr, double* berr) {

	int info;

	vector<double> work(4*n, 0.0);
	vector<int> iwork(n, 0);

	dgesvx_ (fact, trans, &n, &nrhs, a, &lda, af, &ldaf, ipiv,
		equed, r, c, b, &ldb, x, &ldx, rcond, ferr, berr,
		&work[0], &iwork[0], &info);

	return info;
}

int zgesvx (char const* fact, char const* trans, int n, int nrhs,
	complex<double>* a, int lda, complex<double>* af, int ldaf, int* ipiv, 
	char equed[2], double* r, double* c, complex<double>* b, int ldb,
	complex<double>* x, int ldx, double* rcond, double* ferr, double* berr) {

	int info;
	vector<double> rwork(2*n, 0.0);
	vector<complex<double>> work(2*n, 0.0);

	zgesvx_ (fact, trans, &n, &nrhs, a, &lda, af, &ldaf, ipiv,
		equed, r, c, b, &ldb, x, &ldx, rcond, ferr, berr,
		&work[0], &rwork[0], &info);
	return info;
}

int
dgeqrf (int m, int n, double* a, int lda, double* tau) {
	int info;
	int lwork = n*m;
	vector<double> work(lwork, 0.0);
	dgeqrf_(&m, &n, a, &lda, tau, &work[0], &lwork, &info);
	return info;
}

int
dgeqp3 (int m, int n, double* a, int lda, int* jpvt, double* tau) {

	int info;
	// workspace query
	int lwork{-1};
	double query;
	dgeqp3_ (&m, &n, a, &lda, jpvt, tau, &query, &lwork, &info);
	lwork = query;
	assert(lwork > 0);
	assert(info == 0);

	vector<double> work(lwork, 0.0);

	dgeqp3_ (&m, &n, a, &lda, jpvt, tau, &work[0], &lwork, &info);

	return info;
}

double
dtrcon(char const* norm, char const* uplo, char const* diag, int n,
		double* a, int lda) {
	int info;
	vector<double> work(3*n);
	vector<int> iwork(n);
	double rcond;

	dtrcon_ (norm, uplo, diag, &n, a, &lda, &rcond,
			&work[0], &iwork[0], &info);

	return rcond;
}

int
dsyev(char const* jobz, char const* uplo, int n, double* a, int lda, double* w) {
// Solution of the standard symmetric eigenvalue problem Ax = \lambda x
	int info;
	// workspace query
	int lwork{-1};
	double query;
	dsyev_(jobz, uplo, &n, a, &lda, w, &query, &lwork, &info);
	lwork = query;
	assert(lwork > 0);
	assert(info == 0);
	vector<double> work(lwork, 0.0);
	dsyev_(jobz, uplo, &n, a, &lda, w, &work[0], &lwork, &info);
	return info;
}

int
dsygv(int itype, char const* jobz, char const* uplo, int n, double* a,
	int lda, double* b,  int ldb, double* w) {
// Solution of the generalized symmetric eigenvalue problem
//   Ax = \lambda Bx  (itype=1),
// with B positive definite
	int info;
	// workspace query
	int lwork{-1};
	double query;
	dsygv_(&itype, jobz, uplo, &n, a, &lda, b, &ldb, w, &query, &lwork, &info);
	lwork = query;
	assert(lwork > 0);
	assert(info == 0);
	vector<double> work(lwork, 0.0);
	dsygv_(&itype, jobz, uplo, &n, a, &lda, b, &ldb, w, &work[0], &lwork, &info);
	// if (info != 0)
	// 	throw runtime_error(vastr("dsygv returned error ",info));
	return info;
}


int
dgeev (char const* jobvl, char const* jobvr, int n, double* a, int lda,
	double* wr, double* wi, double* vl, int ldvl, double* vr, int ldvr) {

	int info;
	// workspace query
	int lwork{-1};
	double query;
	dgeev_ (jobvl, jobvr, &n, a, &lda, wr, wi, vl, &ldvl,
				vr, &ldvr, &query, &lwork, &info);
	lwork = query;
	assert(lwork > 0);
	assert(info == 0);

	vector<double> work(lwork, 0.0);

	if (tolower(jobvl[0]) == 'n')
		vl = &work[0];
	if (tolower(jobvr[0]) == 'n')
		vr = &work[0];

	dgeev_ (jobvl, jobvr, &n, a, &lda, wr, wi, vl, &ldvl,
				vr, &ldvr, &work[0], &lwork, &info);

	return info;
}

int dgeevx (char const* balanc, char const* jobvl, char const* jobvr,
	char const* sense, int n, double* a, int lda, double* wr, double* wi,
	double* vl, int ldvl, double* vr, int ldvr,
	int& ilo, int& ihi, double* scale, double& abnrm, double* rconde, double* rcondv) {

	vector<int> iwork(2*n-2);

	int info;
	// workspace query
	int lwork{-1};
	double query;

	dgeevx_(balanc, jobvl, jobvr, sense, &n, a, &lda, wr, wi, vl, &ldvl,
		vr, &ldvr, &ilo, &ihi, scale, &abnrm, rconde, rcondv, &query,
		&lwork, iwork.data(), &info);

	lwork = query;
	if (info != 0)
		throw runtime_error(vastr("dgeevx query returned info ",info));
	if (lwork <= 0)
		throw runtime_error(vastr("dgeevx returned query ",query));

	vector<double> work(lwork, 0.0);

	dgeevx_(balanc, jobvl, jobvr, sense, &n, a, &lda, wr, wi, vl, &ldvl,
		vr, &ldvr, &ilo, &ihi, scale, &abnrm, rconde, rcondv, work.data(),
		&lwork, iwork.data(), &info);

	return info;
}

int
dggev (char const* jobvl, char const* jobvr, int n, double* a, int lda,
	double* b, int ldb, double* alphar, double* alphai, double* beta,
	double* vl, int ldvl, double* vr, int ldvr) {
	T_(Trace trc(1,"dggev");)

	// workspace query
	int info{0};
	int lwork{-1};
	double query;
	dggev_ (jobvl, jobvr, &n, a, &lda, b, &ldb, alphar, alphai, beta,
				vl, &ldvl, vr, &ldvr, &query, &lwork, &info);
	lwork = query;
	assert(lwork > 0);
	assert(info == 0);
	T_(trc.dprint("using lwork ",lwork);)

	vector<double> work(lwork, 0.0);

	dggev_ (jobvl, jobvr, &n, a, &lda, b, &ldb, alphar, alphai, beta,
				vl, &ldvl, vr, &ldvr, &work[0], &lwork, &info);
	return info;
}

int
zggev (char const* jobvl, char const* jobvr, int n, complex<double>* a, int lda,
	complex<double>* b, int ldb, complex<double>* alpha, complex<double>* beta,
	complex<double>* vl, int ldvl, complex<double>* vr, int ldvr) {
	T_(Trace trc(1,"zggev");)

	// workspace query
	int info = 0;
	int lwork{-1};
	complex<double> query;
	vector<double> rwork(8*n, 0.0);
	zggev_ (jobvl, jobvr, &n, a, &lda, b, &ldb, alpha, beta,
				vl, &ldvl, vr, &ldvr, &query, &lwork, &rwork[0], &info);
	lwork = query.real();
	assert(lwork > 0);
	assert(info == 0);
	T_(trc.dprint("using lwork ",lwork);)

	vector<complex<double> > work(lwork, complex<double>(0.0));

	zggev_ (jobvl, jobvr, &n, a, &lda, b, &ldb, alpha, beta,
				vl, &ldvl, vr, &ldvr, &work[0], &lwork, &rwork[0], &info);
	return info;
}

int
zggevx (char const* balanc, char const* jobvl, char const* jobvr,
	char const* sense, int n, complex<double>* a, int lda,
	complex<double>* b, int ldb, complex<double>* alpha, complex<double>* beta,
	complex<double>* vl, int ldvl, complex<double>* vr,
	int ldvr, double* rconde, double* rcondv) {

	T_(Trace trc(1,"zggevx");)
	int info = 0;
	int ilo, ihi;
	double abnrm, bbnrm;

	vector<double> rwork(6*n, 0.0);
	vector<int> iwork(n+2, 0);
	vector<int> bwork(n, 0);   // XXX Fortran logical
	vector<double> rscale(n, 0.0);
	vector<double> lscale(n, 0.0);

	// workspace query
	int lwork{-1};
	complex<double> query;
	zggevx_ (balanc, jobvl, jobvr, sense, &n, a, &lda,
			b, &ldb, alpha, beta, vl, &ldvl, vr, &ldvr,
			&ilo, &ihi, &lscale[0], &rscale[0], &abnrm, &bbnrm, rconde, rcondv,
			&query, &lwork, &rwork[0], &iwork[0], &bwork[0], &info);
	lwork = query.real();
	assert(lwork > 0);
	assert(info == 0);
	T_(trc.dprint("using lwork ",lwork);)

	vector<complex<double>> work(lwork, complex<double>(0.0));

	T_(trc.dprint("zggevx n ",n,", sense ",sense);)
	
	// Generalized eigenvalue routines have had problems in the
	// past when scaling is done (see mail re deprecated routine
	// CGEGV)

	zggevx_ (balanc, jobvl, jobvr, sense, &n, a, &lda,
			b, &ldb, alpha, beta, vl, &ldvl, vr, &ldvr,
			&ilo, &ihi, &lscale[0], &rscale[0], &abnrm, &bbnrm, rconde, rcondv,
			&work[0], &lwork, &rwork[0], &iwork[0], &bwork[0], &info);
	T_(trc.dprint("ZGGEVX returned info ",info);)
	return info;
}

int
dormqr (char const* side, char const* trans, int m, int n, int k,
		double const* a, int lda, double const* tau, double* c, int ldc) {

	int info;

	// workspace query
	int lwork{-1};
	double query;
	dormqr_ (side, trans, &m, &n, &k, a, &lda,
		tau, c, &ldc, &query, &lwork, &info);
	lwork = query;
	assert(lwork > 0);
	assert(info == 0);

	vector<double> work(lwork, 0.0);

	dormqr_ (side, trans, &m, &n, &k, a, &lda,
		tau, c, &ldc, &work[0], &lwork, &info);
	return info;
}

void
dtrsm (char const* side, char const* uplo, char const* transa,
		char const* diag, int m, int n, double alpha, double *a,
		int lda, double *b, int ldb) {
	dtrsm_ (side, uplo, transa, diag, &m, &n, &alpha, a, &lda, b, &ldb);
}

int
dgesvd (char const* jobu, char const* jobvt, int m, int n, double* a, int lda,
	double* s, double* u, int ldu, double* vt, int ldvt) {

	int info;
	// workspace query
	int lwork{-1};
	double query;
	dgesvd_ (jobu, jobvt, &m, &n, a, &lda,
				s, u, &ldu, vt, &ldvt,
				&query, &lwork, &info);
	lwork = query;
	assert(lwork > 0);
	assert(info == 0);

	vector<double> work(lwork, 0.0);

	dgesvd_ (jobu, jobvt, &m, &n, a, &lda,
				s, u, &ldu, vt, &ldvt,
				&work[0], &lwork, &info);

	return info;
}

string
eigstring (complex<double> w, int n, complex<double>* vr) {
// Create a string with a one-line summary of a complex eigenvalue/eigenvector
	ostringstream os;
	os << std::setw(15) << std::showpoint << std::left
		<< std::setprecision(8) << w.real();
	os << std::setw(15) << std::showpoint << std::left
		<< std::setprecision(8) << w.imag();
	// normalize vr: largest component 1.0
	vector<complex<double> > v(n);
	blas::copy(n, vr, 1, &v[0], 1);
	blas::normalize_inf(n, &v[0], 1, true, 1.0);
	os << flaps::summarize(v, 80);
	return os.str();
}

string
eigstring (double w, int n, double* vr) {
// Create a string with a one-line summary of a real eigenvalue/eigenvector
	ostringstream os;
	os << std::setw(15) << std::showpoint << std::left
		<< std::setprecision(8) << w;
	// normalize vr: largest component 1.0
	vector<double> v(n);
	blas::copy(n, vr, 1, &v[0], 1);
	blas::normalize_inf(n, &v[0], 1, true, 1.0);
	os << flaps::summarize(v, 80);
	return os.str();
}

static void
filter (int n, complex<double>* a, string const& title) {
// delete very small elements from matrix "a"
	double eps{sqrt(std::numeric_limits<double>::epsilon())};
	vector<complex<double> > zeroed(n*n, complex<double>(0.0));
	for (int i=0; i<n; i++) {
		double rownorm = blas::scnrm2(n, &a[IJ(i,0,n)], n);
		for (int j=0; j<n; j++) {
			double colnorm = blas::scnrm2(n, &a[IJ(0,j,n)], 1);
			double aa = abs(a[IJ(i,j,n)]);
			if (aa < eps*(rownorm+colnorm)) {
				zeroed[IJ(i,j,n)] = a[IJ(i,j,n)];
				a[IJ(i,j,n)] = complex<double>(0.0);
			}
		}
	}
}

void
sortEigen (size_t n, size_t m, complex<double> *s, complex<double>* vl, complex<double>* vr) {
/*
 * Sort a set of m eigenvalues and their corresponding right
 * and left eigenvectors by increasing imag part
 * Use a selection sort (see Numerical Recipes in C p 330)
 */
	int i, j;
	vector<complex<double>> zl(n, complex<double>(0.0));
	vector<complex<double>> zr(n, complex<double>(0.0));

	for (j=1; j<(int)m; j++) {
		complex<double> t = s[j];
		if (vl)
			blas::copy (n, &vl[IJ(0,j,n)], 1, &zl[0], 1);
		if (vr)
			blas::copy (n, &vr[IJ(0,j,n)], 1, &zr[0], 1);
		i = j - 1;
		// sort criteria:
		// while(i >= 0 && abs(s[i].imag()) > abs(t.imag())) 
		while(i >= 0 && abs(s[i]) > abs(t)) {
			s[i+1] = s[i];
			if (vl)
				blas::copy (n, &vl[IJ(0,i,n)], 1, &vl[IJ(0,i+1,n)], 1);
			if (vr)
				blas::copy (n, &vr[IJ(0,i,n)], 1, &vr[IJ(0,i+1,n)], 1);
			i--;
		}
		s[i+1] = t;
		if (vl)
			blas::copy (n, &zl[0], 1, &vl[IJ(0,i+1,n)], 1);
		if (vr)
			blas::copy (n, &zr[0], 1, &vr[IJ(0,i+1,n)], 1);
	}
}

static void
insert(int m, int nb, int ib, int jb, const double* a,
	double scale, double* B) {
// insert (m,m) matrix "a" into 0b block (ib,jb) of (nb,nb) matrix "B"
	for (int i=0; i<m; i++)
		for (int j=0; j<m; j++)
			B[IJ(i+m*ib,j+m*jb,nb)] = scale*a[IJ(i,j,m)];
}

static void
tinsert(int m, int nb, int ib, int jb, const double* a,
	double scale, double* B) {
// insert the transpose of (m,m) matrix "a" into 0b block (ib,jb) of (nb,nb) matrix "B"
	for (int i=0; i<m; i++)
		for (int j=0; j<m; j++)
			B[IJ(i+m*ib,j+m*jb,nb)] = scale*a[IJ(j,i,m)];
}

static void
cinsert(size_t m, size_t nb, size_t ib, size_t jb, complex<double> const* a, complex<double>* B, complex<double> scale) {
	for (size_t i=0; i<m; i++) {
		for (size_t j=0; j<m; j++) {
			B[IJ(i+m*ib,j+m*jb,nb)] = scale*a[IJ(i,j,m)];
		}
	}
}

double
cmatrixNorm2(size_t n, complex<double> const* a) {
	double rval = blas::scnrm2(n*n, a, 1);
	return rval;
}

int
perfind (int n, const vector<double>& A, const vector<double>& B,
		double w, double* x) {
// compute a performance index for the generalized symmetric eigenvalue
// problem Ax = \lambda B x
// where 0 is bad and 14 is about as good as possible
	T_(Trace trc(1,"perfind(A,B)");)
	T_(trc.dprintm(n,1,n,x,"eigenvector");)
	int nsq{n*n};
	vector<double> D{A};
	for (int i=0; i<nsq; i++)
		D[i] -= w*B[i];
	T_(trc.dprintm(n,n,n,D,"D matrix");)
	double alpha{1.0};
	double beta{0.0};
	vector<double> res(n, 0.0);
	blas::gemv("n",n,n,alpha,&D[0],n,&x[0],1,beta,&res[0],1);
	double rnorm = blas::snrm2(n, &res[0], 1);
	double dnorm = blas::snrm2(nsq, &D[0], 1);
	double eps{sqrt(std::numeric_limits<double>::epsilon())};
	if (rnorm < eps)
		return 14;
	if (dnorm < eps)
		return 0.0;
	double t = -log10(rnorm/dnorm);
	int rval{0};
	if (t < 1)
		rval = 0;
	else
		rval = t;
	T_(trc.dprint("residual ",rnorm,", D ",dnorm,", t ",t,", eps ",eps,", pi ",rval);)
	return rval;
}

int
perfind (int n, const vector<double>& m, const vector<double>& g,
		const vector<double>& k, double w, complex<double>* x) {
// compute a performance index for the free-vibration gyro problem
// where 0 is bad and 14 is about as good as possible
	T_(Trace trc(1,"perfind");)
	int rval{0};

	T_(trc.dprintm(n,1,n,x,"eigenvector");)
	int nsq{n*n};
	vector<complex<double>> D(nsq, 0.0);
	complex<double> iw(0.0, w);
	double ww{w*w};
	for (int i=0; i<nsq; i++)
		D[i] = k[i];
	for (int i=0; i<nsq; i++)
		D[i] += iw*g[i];
	for (int i=0; i<nsq; i++)
		D[i] += -ww*m[i];
	T_(trc.dprintm(n,n,n,D,"D matrix");)
	complex<double> alpha(1.0);
	complex<double> beta(0.0);
	vector<complex<double>> res(n, 0.0);
	blas::gemv("n",n,n,alpha,&D[0],n,&x[0],1,beta,&res[0],1);
	double rnorm = blas::scnrm2(n, &res[0], 1);
	double dnorm = blas::scnrm2(nsq, &D[0], 1);
	double eps{std::numeric_limits<double>::epsilon()};
	T_(trc.dprint("residual ",rnorm,", D ",dnorm,", eps ",eps);)
	if (dnorm < eps)
		rval = 0.0;
	else {
		double t = rnorm/dnorm;
		if (t < eps)
			rval = 14;
		else
			t = -log10(t);
		if (t > 1)
			rval = t;
	}
	T_(trc.dprint("pi ",rval);)
	return rval;
}

int
gyroeig (int n, const vector<double>& m, const vector<double>& g,
		const vector<double>& k, vector<double>& w,
		vector<complex<double>>& x, vector<int>& pi) {
// Solution of the nth-order free-vibration problem with gyroscopic terms:
//   (\omega^2 m + i\omega g + k)x = 0
// using the method proposed in [1], where the nth order problem is cast
// into a 2nth order problem with real, symmetric matrices:
//   \omega^2 I y = Ky
// where
//    I = | m  0 |      K = | g'm^{-1}g + k    g'm^{-1}k |
//        | 0  k |          |   km^{-1}g        km^{-1}k |
// are real and symmetric.
// Reference:
// @article{meirovitch1974new,
//   title={A new method of solution of the eigenvalue problem
//         for gyroscopic systems},
//   author={Meirovitch, Leonard},
//   journal={AiAA Journal},
//   volume={12},
//   number={10},
//   pages={1337--1342},
//   year={1974}
// }
	T_(Trace trc(1,"gyroeig");)

	int n2{2*n};
	assert(m.size() == n*n);
	assert(m.size() == g.size());
	assert(m.size() == k.size());
	T_(trc.dprintm(n,n,n,m,"gyroeig_m");)
	T_(trc.dprintm(n,n,n,g,"gyroeig_g");)
	T_(trc.dprintm(n,n,n,k,"gyroeig_k");)

	// create I
	vector<double> I(n2*n2, 0.0);
	insert(n, n2, 0, 0, &m[0], 1.0, &I[0]);
	insert(n, n2, 1, 1, &k[0], 1.0, &I[0]);

	// create K
	vector<double> K(n2*n2, 0.0);
	// Cholesky factorization of the mass matrix
	vector<double> mi{m};
	int info = lapack::dpotrf("u", n, &mi[0], n);
	if (info != 0)
		throw runtime_error("singular mass matrix");
	// g'm^{-1}g + k  (gmig)
	vector<double> mig{g};
	info = lapack::dpotrs("u", n, n, &mi[0], n, &mig[0], n);
	vector<double> gmig{k};
	blas::gemm("t","n",n,n,n,1.0,&g[0],n,&mig[0],n, 1.0, &gmig[0], n);
	insert(n, n2, 0, 0, &gmig[0], 1.0, &K[0]);
	// g'm^{-1}k
	vector<double> mik{k};
	info = lapack::dpotrs("u", n, n, &mi[0], n, &mik[0], n);
	vector<double> gmik{k};
	blas::gemm("t","n",n,n,n,1.0,&g[0],n,&mik[0],n, 0.0, &gmik[0], n);
	insert(n, n2, 0, 1, &gmik[0], 1.0, &K[0]);
	tinsert(n, n2, 1, 0, &gmik[0], 1.0, &K[0]);
	// km^{-1}k
	vector<double> kmik{k};
	blas::gemm("n","n",n,n,n,1.0,&k[0],n,&mik[0],n, 0.0, &kmik[0], n);
	insert(n, n2, 1, 1, &kmik[0], 1.0, &K[0]);

	T_(trc.dprintm(n2,n2,n2,K,"K matrix");)
	T_(trc.dprintm(n2,n2,n2,I,"I matrix");)
	// save K and I - destroyed in dsygv
	vector<double> Kp{K};
	vector<double> Ip{I};

	// solve the eigenvalue problem \omega^2 I y = Ky
	vector<double> ev(n2, 0.0);
	double eps{sqrt(std::numeric_limits<double>::epsilon())};
	info = lapack::dsygv(1, "v", "u", n2, &K[0], n2, &I[0], n2, &ev[0]);
	if (info != 0) {
		// failure: try to figure out which matrix is the culprit
		T_(trc.dprint("dsygv returned ",info);)
		// (k,m) eigenvalues
		vector<double> w(n,0.0);
		vector<double> mass{m};
		vector<double> stif{k};
		info = lapack::dsygv(1, "v", "u", n, &stif[0], n, &mass[0], n, &w[0]);
		if (info != 0)
			throw runtime_error(vastr("solution of the free-vibration (K,M) "
				"eigenvalue problem failed: dsygv returned ",info));
		T_(trc.dprint("eigenvalues of (stif,mass):",w);)
		if (w[0] < -eps*w[n-1]) {
			// mass matrix eigenvalues
			mass = m;
			info = dsyev("v", "u", n, &mass[0], n, &w[0]);
			if (w[0] <= 0.0)
				throw runtime_error(vastr("solution of the free-vibration (K,G,M) "
					"problem failed: the mass matrix is not positive definite:",w));
			// stiffness matrix eigenvalues
			stif = k;
			info = dsyev("v", "u", n, &stif[0], n, &w[0]);
			if (w[0] <= 0.0) {
				T_(trc.dprint("eigenvalues of stif:",w);)
				throw runtime_error(vastr("solution of the free-vibration (K,G,M) "
					"problem failed: the stiffness matrix is not positive definite: ",
					w[0]," ",flaps::summarize(n,&stif[0],80)));
			}
		}
		throw runtime_error(vastr("solution of the free-vibration (K,G,M) "
			"eigenvalue problem failed: eigenvalues of the (K,M) problem:",w));
	}
	// check results from dsygv
	vector<int> pis;
	vector<double> zero(n2*n2, 0.0);
	for (int j=0; j<n2; j++) {
		pis.push_back(perfind(n2,Kp,Ip,ev[j],&K[j*n2]));
	}
	T_(trc.dprint("perf ind for dsygv:",pis);)
	T_(trc.dprint("eigenvalues:\n",ev);)
	T_(trc.dprintm(n2,n2,n2,K,"eigenvectors");)

	// there are 2n eigenvalues and vectors; there will be 2 of
	// each eigenvalue - the 2 corresponding eigenvectors are the
	// real and imaginary parts of x
	x = vector<complex<double>>(n*n, 0.0);
	w.clear();
	pi.clear();
	int nsig{3}; // pairs of eigenvalues from gyroeig should match
	             // to within 'nsig' significant figures
	for (int j=0; j<n; j++) {
		if (!is_equal(ev[2*j],ev[2*j+1], nsig))
			throw runtime_error(vastr("gyroeig w[",2*j,"] 0b = ",ev[2*j],
				" does not equal w[",2*j+1,"] = ",ev[2*j+1]));
		for (int i=0; i<n; i++)
			x[i+j*n] = complex<double>(K[i+n2*2*j], K[i+n2*(2*j+1)]);
		blas::normalize_inf(n, &x[j*n], 1, true);
		// the eigenvalues are \omega^2 - return \omega but only one of
		// (+\omega, -\omega) is a solution: determine which by computing
		// a performance index
		w.push_back(sqrt(ev[2*j]));
		int p = perfind(n,m,g,k,w[j],&x[j*n]);
		if (p < 8) {
			// -\omega is the solution, but we can conjugate the eigenvector
			// instead of using negative frequency
			for (int i=0; i<n; i++)
				x[i+j*n] = complex<double>(K[i+n2*2*j], -K[i+n2*(2*j+1)]);
			blas::normalize_inf(n, &x[j*n], 1, true);
			p = perfind(n,m,g,k,w[j],&x[j*n]);
		}
		pi.push_back(p);
	}
	T_(trc.dprint("frequencies (rad/s): ",w);)
	T_(trc.dprint("performance indices\n",pi);)
	return 0;
}

int
poly (int n, vector<double*>& A, vector<complex<double>>& w, vector<complex<double>>& x) {
// Solution of the polynomial eigenvalue problem
//   (A_0 + \lambda A_1 + ... + \lambda^r A_r) x = 0
// by putting it into generalized eigenvalue form:
//   \bar{A} \bar{x} = \lambda \bar{B} \bar{x}
// where
//   \bar{A} = [ A_0  0 ]
//             [  0   I ]
//   \bar{B} = [ -A_1  -A_2 ]
//             [  I    0  ]
//   \bar{x} = {  x }
//             { \lambda x }
//
	T_(Trace trc(1,"poly"," n ",n,", ",A.size()," matrices");)
	int r = A.size() - 1;

	// Abar and Bbar are (r*n, r*n)
	int rn = r*n;
	vector<double> Abar(rn*rn, 0.0);
	vector<double> Bbar(rn*rn, 0.0);

	// identity matrices: Abar diagonals (except 0,0),
	// Bbar subdiagonals: (1,0) (2,1)...
	vector<double> identity(n*n, 0.0);
	for (int j=0; j<n; j++)
		identity[IJ(j,j,n)] = 1.0;

	for (int j=0; j<r-1; j++) {
		insert(n, rn, j+1, j, &identity[0], -1.0, &Bbar[0]);
		insert(n, rn, j+1, j+1, &identity[0], 1.0, &Abar[0]);
	}
	// A_0: Abar(0,0)
	insert(n, rn, 0, 0, A[0], 1.0, &Abar[0]);
	// A_j j=1,r go in the first row of Bbar
	for (int j=1; j<=r; j++)
		insert(n, rn, 0, j-1, A[j], -1.0, &Bbar[0]);

	vector<double> alphar(rn, 0.0);
	vector<double> alphai(rn, 0.0);
	vector<double> beta(rn, 0.0);
	vector<double> vlbar(rn*rn, 0.0);
	vector<double> vrbar(rn*rn, 0.0);

	T_(trc.dprintm(rn,rn,rn,Abar,"poly Abar");)
	T_(trc.dprintm(rn,rn,rn,Bbar,"poly Bbar");)

	int info = dggev ("n", "V", rn, &Abar[0], rn, &Bbar[0], rn, &alphar[0],
		&alphai[0], &beta[0], &vlbar[0], rn, &vrbar[0], rn);

	if (info != 0)
		throw runtime_error(vastr("error ",info," returned from dggev"));

	T_(trc.dprint("alphar",alphar,"alphai",alphai);)
	T_(trc.dprint("beta",beta);)
	T_(trc.dprintm(rn,rn,rn,vrbar,"vr");)

	// dggev returns w = (alphar+ialphai)/beta
	// Save the first n elements of the left and right
	// eigenvectors in vr unless they are zero,
	// in which case ignore the eigenvalue also
	int neig{0};

	double eps = sqrt(std::numeric_limits<double>::epsilon());
	complex<double> big(1.0/eps, 1.0/eps);
	w = vector<complex<double>>(rn, complex<double>(0.0));
	x = vector<complex<double>>(rn*rn, complex<double>(0.0));
	for (int j=0; j<rn; j++) {
		double vnorm = blas::snrm2(n, &vrbar[IJ(0,j,rn)], 1);
		if (vnorm < eps) continue;
		if (beta[j] == 0.0)
			beta[j] = eps;
		// if the jth eigenvalue is real (alphai[j]==0) col j of vr is
		// the (real) eigenvector
		if (alphai[j] == 0.0 && vnorm > eps) {
			if (alphar[j] >= 0.0) {
				for (int i=0; i<n; i++)
					x[IJ(i,neig,n)] = complex<double>(vrbar[IJ(i,j,rn)]);
				w[neig] = complex<double>(alphar[j])/beta[j];
			} else
				continue;
		} else {
			double vnormi = blas::snrm2(n, &vrbar[IJ(0,j+1,rn)], 1);
			if (vnorm > eps || vnormi > eps) {
				for (int i=0; i<n; i++)
					x[IJ(i,neig,n)] = complex<double>(vrbar[IJ(i,j,rn)],vrbar[IJ(i,j+1,rn)]);
				w[neig] = complex<double>(alphar[j],alphai[j])/beta[j];
				j++;
			} else
				continue;
		}
		T_(trc.dprint("eig[",neig,"]: ",w[neig]," vector: ", flaps::summarize(n,&x[IJ(0,neig,n)],80));)
		neig++;
	}

	sortEigen (n, n, &w[0], nullptr, &x[0]);

	T_(trc.dprint("returning ",neig," eigenpair");)
	return neig;
}

int
polyeig (size_t n, int maxeig, vector<complex<double>*>& A,
		complex<double>* w, complex<double>* vl, complex<double>* vr, int* pil, int* pir,
		int nrefine) {
// Solution of the polynomial eigenvalue problem
//   (A_0 + \lambda A_1 + ... + \lambda^r A_r) x = 0
// by putting it into generalized eigenvalue form:
//   \bar{A} \bar{x} = \lambda \bar{B} \bar{x}
// where for r=2
//   \bar{A} = [ A_0  0 ]
//             [  0   I ]
//   \bar{B} = [ -A_1  -A_2 ]
//             [  I    0  ]
//   \bar{x} = {  y }
//             { sy }
//
	T_(Trace trc(1,"polyeig"," n ",n,", ",A.size()," matrices");)
	size_t r = A.size() - 1;
	ostringstream os;
	size_t j;
	int neig{0};
	complex<double> one{1.0};
	double gamma{1.0};

	// first and last point must not be null
	if (A[0] == nullptr) {
		os << "attempt to solve a polynomial eigenvalue problem "
			<< "with A[0] null";
		T_(trc.dprint("throwing exception: ",os.str());)
		throw runtime_error(os.str());
	}
	if (A[r] == nullptr) {
		os << "attempt to solve a polynomial eigenvalue problem with A["
			<< r << "] null";
		T_(trc.dprint("throwing exception: ",os.str());)
		throw runtime_error(os.str());
	}

	// special case: 3 matrices w/middle one null:
	//    A_0 + s^2 A_2 = 0
	// zggevx solves A_0 - \lambda A_2 so we change the sign on alpha
	//
	if (A.size() == 3 && A[1] == nullptr) {
		vector<complex<double> > alpha(n);
		vector<complex<double> > beta(n);
		vector<double> rconde(n);
		vector<double> rcondv(n);
		double eps(std::numeric_limits<double>::epsilon());
		T_(trc.dprintm(n,n,n,A[0],"polyeig A");)
		T_(trc.dprintm(n,n,n,A[2],"polyeig B");)
		// int info = zggevx ("N", "V", "V", "B", n, A[0], n, A[2], n,
	 	// 	&alpha[0], &beta[0], vl, n, vr, n, &rconde[0], &rcondv[0]);
		int info = zggev ("V", "V", n, A[0], n, A[2], n,
				&alpha[0], &beta[0], vl, n, vr, n);
		// throw an exception if zggevx returned info != 0
		if (info != 0) {
			if (info < 0)
				os << "argument " << info << " to zggevx had an illegal value";
			else if (info <= (int)n) {
				os << "zggevx failed (eigenvalue " << info << ")\n";
				os << strArray(n,n,n,&A[0][0]) << endl;
				os << strArray(n,n,n,&A[2][0]) << endl;
			} else {
				os << "error " << info << " returned from zggevx";
			}
			T_(trc.dprint("throwing exception: ",os.str());)
			throw runtime_error(os.str());
		}
		for (size_t i=0; i<n; i++) {
			if (abs(beta[i]) > eps) {
				alpha[i] = -alpha[i];
				w[i] = alpha[i]/beta[i];
				w[i] = sqrt(w[i]);
				if (w[i].imag() < 0.0)
					w[i] = -w[i];
			} else {
				w[i] = 1.0/eps;
			}
			T_(trc.dprint("w[",i,"] = sqrt(",alpha[i],'/',beta[i],") = ", w[i], "  ",flaps::summarize(n,&vr[IJ(0,i,n)]));)
		}
		sortEigen (n, n, w, vl, vr);
		return n;
	}

	/*
	 * scale the matrices if there are three
	 * using the scheme proposed in "Fan, H-Y, Lin, W-W, and Van Dooren, P.,
	 * Normwise scaling of second order polynomial matrices, SIAM Journal on Matrix
	 * Analysis and Applications, 2006; 26:252-256" and used in
	 * "Higham, N.J., Mackey, D.S., Tisseur, F., and Garvey, S.D.,
	 * Scaling, sensitivity and stability in the numerical
	 * solution of quadratic eigenvalue problems, Int. J. Numer. Meth. Engng.,
	 *	73:344-360, 2008"
	 */
	vector<complex<double>> scale(A.size(), complex<double>(1.0));
	char* normscale = getenv("SCALE");
	if (A.size() == 3 && normscale != nullptr) {
		double k = cmatrixNorm2(n, A[0]);
		double d = cmatrixNorm2(n, A[1]);
		double m = cmatrixNorm2(n, A[2]);
		double eps(std::numeric_limits<double>::epsilon());
		T_(trc.dprint("norm(A) = ",k,", norm(B) = ",d,", norm(C) = ",m);)
		if (abs(m) > eps) {
			gamma = sqrt(k/m);
			double delta = 2.0/(k + d*gamma);
			scale[0] = delta;
			scale[1] = gamma*delta;
			scale[2] = gamma*gamma*delta;
			T_(trc.dprint("gamma ",gamma,", delta ",delta);)
		}
	}

	T_(trc.dprintm(n,n,n,A[0],"polyeig A[0]");)
	T_(trc.dprintm(n,n,n,A[1],"polyeig A[1]");)
	T_(trc.dprintm(n,n,n,A[2],"polyeig A[2]");)

	// Abar and Bbar are (r*n, r*n)
	size_t rn = r*n;
	vector<complex<double>> Abar(rn*rn, complex<double>(0.0));
	vector<complex<double>> Bbar(rn*rn, complex<double>(0.0));
	vector<int> pilbar(rn, 0);
	vector<int> pirbar(rn, 0);

	// identity matrices: Abar diagonals (except 0,0),
	// Bbar subdiagonals
	vector<complex<double>> identity(n*n, complex<double>(0.0));
	for (j=0; j<n; j++)
		identity[IJ(j,j,n)] = complex<double>(1.0);

	for (j=0; j<r-1; j++) {
		cinsert(n, rn, j+1, j, &identity[0], &Bbar[0], one);
		cinsert(n, rn, j+1, j+1, &identity[0], &Abar[0], one);
		// cinsert(n, rn, j+1, j+1, A[j+2], &Abar[0], one);
	}
	// A_0: Abar(0,0)
	cinsert(n, rn, 0, 0, A[0], &Abar[0], scale[0]);
	// A_j j=1,r go in the first row of Bbar
	for (j=1; j<=r; j++) {
		if (A[j]) {
			cinsert(n, rn, 0, j-1, A[j], &Bbar[0], -scale[j]);
			// cinsert(n, rn, 0, j-1, A[j], &Bbar[0], scale[j]);
		}
	}

	vector<complex<double>> alpha(rn, complex<double>(0.0));
	vector<complex<double>> beta(rn, complex<double>(0.0));
	vector<complex<double>> vlbar(rn*rn, complex<double>(0.0));
	vector<complex<double>> vrbar(rn*rn, complex<double>(0.0));
	vector<double> rconde(rn, 0.0);
	vector<double> rcondv(rn, 0.0);

	// Filter small off-diagonals - cause cggevx to choke
	filter(rn, &Abar[0], "filtered Abar");
	filter(rn, &Bbar[0], "filtered Bbar");

	T_(trc.dprintm(rn,rn,rn,Abar,"polyeig Abar");)
	T_(trc.dprintm(rn,rn,rn,Bbar,"polyeig Bbar");)

	int info = zggevx ("B", "V", "V", "B", rn, &Abar[0], rn, &Bbar[0], rn,
	 		&alpha[0], &beta[0], &vlbar[0], rn, &vrbar[0], rn, &rconde[0], &rcondv[0]);

	if (info != 0) {
		if (info < 0)
			os << "argument " << info << " to zggevx had an illegal value";
		else if (info <= (int)n) {
			os << "the QZ algorithm failed (eigenvalue " << info << ')';
			os << "A matrix:\n" << strArray(rn,rn,rn,&Abar[0]) << endl;
			os << "B matrix:\n" << strArray(rn,rn,rn,&Bbar[0]) << endl;
		} else {
			os << "error " << info << " returned from zggevx";
		}
		T_(trc.dprint("throwing exception: ",os.str());)
		throw runtime_error(os.str());
	}

	T_(trc.dprint("alpha",alpha);)
	T_(trc.dprint("beta",beta);)
	T_(trc.dprintm(rn,rn,rn,vrbar,"vr");)
	T_(trc.dprintm(rn,rn,rn,vlbar,"vl");)

	// zggevx returns w = alpha/beta
	// Save the first n elements of the left and right
	// eigenvectors in vl, vr unless they are zero,
	// in which case ignore the eigenvalue also
	if (maxeig <= 0)
		maxeig = rn;
	neig = 0;

	double eps = sqrt(std::numeric_limits<double>::epsilon());
	complex<double> big(1.0/eps, 1.0/eps);
	// size_t rm1n = (r-1)*n;
	size_t rm1n = 0;
	for (j=0; j<(size_t)maxeig; j++) {
		double vlnorm = blas::scnrm2(n, &vlbar[IJ(rm1n,j,rn)], 1);
		double vrnorm = blas::scnrm2(n, &vrbar[IJ(rm1n,j,rn)], 1);
		T_(trc.dprint("alpha[",j,"]=",alpha[j],", beta=",beta[j], " vlnorm ",vlnorm," vrnorm ",vrnorm);)
		if (vlnorm > eps && vrnorm > eps) {
			blas::copy(n, &vlbar[IJ(rm1n,j,rn)], 1, &vl[IJ(0,neig,n)], 1);
			blas::copy(n, &vrbar[IJ(rm1n,j,rn)], 1, &vr[IJ(0,neig,n)], 1);
			if (abs(beta[j]) > eps) {
				w[neig] = gamma*alpha[j]/beta[j];
			} else {
				w[neig] = big;
			}
			T_(trc.dprint("eig[",j,"]: ",w[neig]," vector: ", flaps::summarize(n,&vr[IJ(0,neig,n)],80));)
			neig++;
		} else {
			T_(trc.dprint("ignoring eig[",j,"]: zero vector(s)");)
		}
	}

	// sort the eigenvalues prior to refining
	sortEigen (n, neig, w, vl, vr);

	T_(for (int j=0; j<neig; j++))
		T_(trc.dprint("sorted eig[",j,"]: ",w[j]," vector: ", flaps::summarize(n,&vr[IJ(0,j,n)],80));)

	T_(trc.dprint("returning ",neig," eigenpair");)
	return neig;
}

#ifdef NEVER // no Ad specializations
// specializations of triprod
template<>
void
triprod(int n, int m, const double* T, const complex<Ad>* A, complex<Ad>* tat) noexcept {
// triple-product tat = T'AT with T(n,m), A(n,n), tat(m,m)
// assuming Dta*Dtt = Dta e.g. double * complex -> complex
	vector<complex<Ad>> at(n);
	for (int i=0; i<m; i++) {
		for (int j=0; j<m; j++) {
			for (int k=0; k<n; k++)
				blas::dot(n, &T[j*n], 1, &A[k], n, at[k]);
			blas::dot(n, &T[i*n], 1, &at[0], 1, tat[i+j*m]);
		}
	}
}
#endif // NEVER // no Ad specializations

} // namespace lapack

//------------------------------------------------------------------
// dqr1up: rank-one update of a qr factorization
// Translated by f2c, then cleaned up to not require f2c.h or libf2c.a

extern "C"
int
dqr1up_(int *m, int *n, int *k, double* q, int* ldq, double* r__,
		int *ldr, double *u, double *v, double *work);

int
dqr1up(int m, int n, int k, double* q, int ldq, double* r, int ldr,
		double* u, double* v) {
// cover function for dqr1up_
	vector<double> work(2*k, 0.0);
	return dqr1up_(&m, &n, &k, q, &ldq, r, &ldr, u, v, &work[0]);
}


/* Copyright (C) 2008, 2009  VZLU Prague, a.s., Czech Republic */
/* Author: Jaroslav Hajek <highegg@gmail.com> */
/* This file is part of qrupdate. */
/* qrupdate is free software; you can redistribute it and/or modify */
/* it under the terms of the GNU General Public License as published by */
/* the Free Software Foundation; either version 3 of the License, or */
/* (at your option) any later version. */
/* This program is distributed in the hope that it will be useful, */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the */
/* GNU General Public License for more details. */
/* You should have received a copy of the GNU General Public License */
/* along with this software; see the file COPYING.  If not, see */
/* <http://www.gnu.org/licenses/>. */


using integer = int;
using doublereal = double;
using logical = int;

int
dqrder_(integer *m, integer *n, doublereal *q, integer* ldq,
		doublereal* r__, integer* ldr, integer* j, doublereal* w);
int
dqrqh_(integer *, integer *, doublereal*, integer *,
		doublereal *, doublereal *);
int
dqhqr_(integer*, integer*, doublereal*, integer*, doublereal*, doublereal*);
int
dqrtv1_(integer *, doublereal *, doublereal *);
int
dqrot_(const char *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *);
int
dch1up_(integer *, doublereal *, integer *, doublereal *, doublereal *);
int
dgqvec_(integer *, integer *, doublereal *, integer *, doublereal *);


int
dqr1up_(integer *m, integer *n, integer *k, doublereal *
	q, integer *ldq, doublereal *r__, integer *ldr, doublereal *u, 
	doublereal *v, doublereal *w) {
    /* System generated locals */
    integer q_dim1, q_offset, r_dim1, r_offset, i__1, i__2, i__3;
    doublereal d__1;
	// Table of constant values
	int c__1 = 1;

    /* Local variables */
    int i__;
    double ru, ruu;
    int info;
    int full;


/* purpose:      updates a QR factorization after rank-1 modification */
/*               i.e., given a m-by-k orthogonal Q and m-by-n upper */
/*               trapezoidal R, an m-vector u and n-vector v, */
/*               this subroutine updates Q -> Q1 and R -> R1 so that */
/*               Q1*R1 = Q*R + u*v', and Q1 is again orthonormal */
/*               and R1 upper trapezoidal. */
/*               (real version) */
/* arguments: */
/* m (in)        number of rows of the matrix Q. */
/* n (in)        number of columns of the matrix R. */
/* k (in)        number of columns of Q, and rows of R. Must be */
/*               either k = m (full Q) or k = n < m (economical form). */
/* Q (io)        on entry, the orthogonal m-by-k matrix Q. */
/*               on exit, the updated matrix Q1. */
/* ldq (in)      the leading dimension of Q. ldq >= m. */
/* R (io)        on entry, the upper trapezoidal m-by-n matrix R.. */
/*               on exit, the updated matrix R1. */
/* ldr (in)      the leading dimension of R. ldr >= k. */
/* u (io)        the left m-vector. On exit, if k < m, u is destroyed. */
/* v (io)        the right n-vector. On exit, v is destroyed. */
/* w (out)       a workspace vector of size 2*k */

/* quick return if possible. */
    /* Parameter adjustments */
    q_dim1 = *ldq;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    r_dim1 = *ldr;
    r_offset = 1 + r_dim1;
    r__ -= r_offset;
    --u;
    --v;
    --w;

    /* Function Body */
    if (*k == 0 || *n == 0) {
	return 0;
    }
/* check arguments. */
    info = 0;
    if (*m < 0) {
	info = 1;
    } else if (*n < 0) {
	info = 2;
    } else if (*k != *m && (*k != *n || *n > *m)) {
	info = 3;
    } else if (*ldq < *m) {
	info = 5;
    } else if (*ldr < *k) {
	info = 7;
    }
    if (info != 0) {
	// xerbla_("DQR1UP", &info, (ftnlen)6);
	throw runtime_error(vastr("DQR1UP ",info));
	return 0;
    }
	 if (*k == *m)
		 full = 1;
	else
		full = 0;
    // full = *k == *m;
/* in the non-full case, we shall need the norm of u. */
    if (! full) {
	ru = dnrm2_(m, &u[1], &c__1);
    }
/* form Q'*u. In the non-full case, form also u - Q*Q'u. */
    i__1 = *k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	w[i__] = ddot_(m, &q[i__ * q_dim1 + 1], &c__1, &u[1], &c__1);
	if (! full) {
	    d__1 = -w[i__];
	    daxpy_(m, &d__1, &q[i__ * q_dim1 + 1], &c__1, &u[1], &c__1);
	}
    }
/* generate rotations to eliminate Q'*u. */
    dqrtv1_(k, &w[1], &w[*k + 1]);
/* apply rotations to R. */
    dqrqh_(k, n, &r__[r_offset], ldr, &w[*k + 1], &w[2]);
/* apply rotations to Q. */
    dqrot_("B", m, k, &q[q_offset], ldq, &w[*k + 1], &w[2]);
/* update the first row of R. */
    daxpy_(n, &w[1], &v[1], &c__1, &r__[r_dim1 + 1], ldr);
/* retriangularize R. */
    dqhqr_(k, n, &r__[r_offset], ldr, &w[*k + 1], &w[1]);
/* apply rotations to Q. */
/* Computing MIN */
    i__2 = *k, i__3 = *n + 1;
    i__1 = std::min(i__2,i__3);
    dqrot_("F", m, &i__1, &q[q_offset], ldq, &w[*k + 1], &w[1]);
/* in the full case, we're finished */
    if (full) {
	return 0;
    }
/* compute relative residual norm */
    ruu = dnrm2_(m, &u[1], &c__1);
    // ru *= dlamch_("e", (ftnlen)1);
	 ru *= std::numeric_limits<double>::epsilon();
    if (ruu <= ru) {
	return 0;
    }
/* update the orthogonal basis. */
    dscal_(n, &ruu, &v[1], &c__1);
    d__1 = 1. / ruu;
    dscal_(m, &d__1, &u[1], &c__1);
    dch1up_(n, &r__[r_offset], ldr, &v[1], &w[*k + 1]);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	drot_(m, &q[i__ * q_dim1 + 1], &c__1, &u[1], &c__1, &w[*k + i__], &v[
		i__]);
    }
    return 0;
} /* dqr1up_ */

int
dch1up_(integer *n, doublereal *r__, integer *ldr, doublereal *u, doublereal *w)
{
    /* System generated locals */
    integer r_dim1, r_offset, i__1, i__2;

    /* Local variables */
    int i__, j;
    double t, ui, rr;

/* purpose:      given an upper triangular matrix R that is a Cholesky */
/*               factor of a symmetric positive definite matrix A, i.e. */
/*               A = R'*R, this subroutine updates R -> R1 so that */
/*               R1'*R1 = A + u*u' */
/*               (real version) */
/* arguments: */
/* n (in)        the order of matrix R */
/* R (io)        on entry, the upper triangular matrix R */
/*               on exit, the updated matrix R1 */
/* ldr (in)      leading dimension of R. ldr >= n. */
/* u (io)        the vector determining the rank-1 update */
/*               on exit, u contains the rotation sines */
/*               used to transform R to R1. */
/* w (out)       cosine parts of rotations. */

    /* Parameter adjustments */
    r_dim1 = *ldr;
    r_offset = 1 + r_dim1;
    r__ -= r_offset;
    --u;
    --w;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* apply stored rotations, column-wise */
	ui = u[i__];
	i__2 = i__ - 1;
	for (j = 1; j <= i__2; ++j) {
	    t = w[j] * r__[j + i__ * r_dim1] + u[j] * ui;
	    ui = w[j] * ui - u[j] * r__[j + i__ * r_dim1];
	    r__[j + i__ * r_dim1] = t;
	}
/* generate next rotation */
	dlartg_(&r__[i__ + i__ * r_dim1], &ui, &w[i__], &u[i__], &rr);
	r__[i__ + i__ * r_dim1] = rr;
    }
    return 0;
} /* dch1up_ */

int
dgqvec_(integer* m, integer* n, doublereal* q, integer* ldq, doublereal *u) {
    /* System generated locals */
    integer q_dim1, q_offset, i__1;
    doublereal d__1;

    /* Local variables */
    int i__, j;
    double r__;
    int info;
	// Table of constant values
	int c__1 = 1;

/* purpose:      given an orthogonal m-by-n matrix Q, n < m, generates */
/*               a vector u such that Q'*u = 0 and norm(u) = 1. */
/* arguments: */
/* m (in)        number of rows of matrix Q. */
/* n (in)        number of columns of matrix Q. */
/* Q (in)        the orthogonal matrix Q. */
/* ldq (in)      leading dimension of Q. */
/* u (out)       the generated vector. */

/* quick return if possible. */
    /* Parameter adjustments */
    q_dim1 = *ldq;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    --u;

    /* Function Body */
    if (*m == 0) {
	return 0;
    }
    if (*n == 0) {
	u[1] = 1.;
	i__1 = *m;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    u[i__] = 0.;
	}
	return 0;
    }
/* check arguments. */
    info = 0;
    if (*m < 0) {
	info = 1;
    } else if (*n < 0) {
	info = 2;
    } else if (*ldq < *m) {
	info = 4;
    }
    if (info != 0) {
	// xerbla_("DGQVEC", &info, (ftnlen)6);
	throw runtime_error(vastr("DGQVEC ",info));
	return 0;
    }
    j = 1;
L10:
/* probe j-th canonical unit vector. */
    i__1 = *m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	u[i__] = 0.;
    }
    u[j] = 1.;
/* form u - Q*Q'*u */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	r__ = ddot_(m, &q[i__ * q_dim1 + 1], &c__1, &u[1], &c__1);
	d__1 = -r__;
	daxpy_(m, &d__1, &q[i__ * q_dim1 + 1], &c__1, &u[1], &c__1);
    }
    r__ = dnrm2_(m, &u[1], &c__1);
    if (r__ == 0.) {
	++j;
	if (j > *n) {
/* this is fatal, and in theory, it can't happen. */
	    throw runtime_error("fatal: impossible condition in DGQVEC");
	} else {
	    ++j;
	    goto L10;
	}
    }
    d__1 = 1. / r__;
    dscal_(m, &d__1, &u[1], &c__1);
    return 0;
} /* dgqvec_ */


int
dqhqr_(integer *m, integer *n, doublereal *r__, integer *
	ldr, doublereal *c__, doublereal *s) {
    /* System generated locals */
    integer r_dim1, r_offset, i__1, i__2;

    /* Local variables */
    int i__, j;
    double t;
    int ii, info;

/* purpose:      given an m-by-n upper Hessenberg matrix R, this */
/*               subroutine updates R to upper trapezoidal form */
/*               using min(m-1,n) Givens rotations. */
/*               (real version) */
/* arguments: */
/* m (in)        number of rows of the matrix R */
/* n (in)        number of columns of the matrix R */
/* R (io)        on entry, the upper Hessenberg matrix R */
/*               on exit, the updated upper trapezoidal matrix */
/* ldr (in)      leading dimension of R, >= m */
/* c(out)        rotation cosines, size at least min(m-1,n) */
/* s(out)        rotation sines, size at least min(m-1,n) */

/* quick return if possible. */
    /* Parameter adjustments */
    r_dim1 = *ldr;
    r_offset = 1 + r_dim1;
    r__ -= r_offset;
    --c__;
    --s;

    /* Function Body */
    if (*m == 0 || *m == 1 || *n == 0) {
	return 0;
    }
/* check arguments. */
    info = 0;
    if (*m < 0) {
	info = 1;
    } else if (*n < 0) {
	info = 2;
    } else if (*ldr < *m) {
	info = 4;
    }
    if (info != 0) {
	// xerbla_("DQHQR", &info, (ftnlen)5);
	throw runtime_error(vastr("DQHQR ",info));
	return 0;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* apply stored rotations, column-wise */
	t = r__[i__ * r_dim1 + 1];
	ii = std::min(*m,i__);
	i__2 = ii - 1;
	for (j = 1; j <= i__2; ++j) {
	    r__[j + i__ * r_dim1] = c__[j] * t + s[j] * r__[j + 1 + i__ * 
		    r_dim1];
	    t = c__[j] * r__[j + 1 + i__ * r_dim1] - s[j] * t;
	}
	if (ii < *m) {
/* generate next rotation */
	    dlartg_(&t, &r__[ii + 1 + i__ * r_dim1], &c__[i__], &s[i__], &r__[
		    ii + i__ * r_dim1]);
	    r__[ii + 1 + i__ * r_dim1] = 0.;
	} else {
	    r__[ii + i__ * r_dim1] = t;
	}
    }
    return 0;
} /* dqhqr_ */


int
dqrdec_(integer *m, integer *n, integer *k, doublereal *
	q, integer *ldq, doublereal *r__, integer *ldr, integer *j, 
	doublereal *w) {
    /* System generated locals */
    integer q_dim1, q_offset, r_dim1, r_offset, i__1, i__2;

    /* Local variables */
    int i__, info;
	// Table of constant values
	int c__1 = 1;

/* purpose:      updates a QR factorization after deleting */
/*               a column. */
/*               i.e., given an m-by-k orthogonal matrix Q, an k-by-n */
/*               upper trapezoidal matrix R and index j in the range */
/*               1:n+1, this subroutine updates the matrix Q -> Q1 and */
/*               R -> R1 so that Q1 remains orthogonal, R1 is upper */
/*               trapezoidal, and Q1*R1 = [A(:,1:j-1) A(:,j+1:n)], */
/*               where A = Q*R. */
/*               (real version) */
/* arguments: */
/* m (in)        number of rows of the matrix Q. */
/* n (in)        number of columns of the matrix R. */
/* k (in)        number of columns of Q, and rows of R. Must be */
/*               either k = m (full Q) or k = n < m (economical form, */
/*               basis dimension will decrease). */
/* Q (io)        on entry, the unitary m-by-k matrix Q. */
/*               on exit, the updated matrix Q1. */
/* ldq (in)      leading dimension of Q. ldq >= m. */
/* R (io)        on entry, the original matrix R. */
/*               on exit, the updated matrix R1. */
/* ldr (in)      leading dimension of R. ldr >= k. */
/* j (in)        the position of the deleted column in R. */
/*               1 <= j <= n. */
/* w (o)         a workspace vector of size k-j. */

/* quick return if possible. */
    /* Parameter adjustments */
    q_dim1 = *ldq;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    r_dim1 = *ldr;
    r_offset = 1 + r_dim1;
    r__ -= r_offset;
    --w;

    /* Function Body */
    if (*m == 0 || *n == 0 || *j == *n) {
	return 0;
    }
/* check arguments. */
    info = 0;
    if (*m < 0) {
	info = 1;
    } else if (*n < 0) {
	info = 2;
    } else if (*k != *m && (*k != *n || *n >= *m)) {
	info = 3;
    } else if (*ldq < *m) {
	info = 5;
    } else if (*ldr < *k) {
	info = 7;
    } else if (*j < 1 || *j > *n + 1) {
	info = 8;
    }
    if (info != 0) {
	// xerbla_("DQRDEC", &info, (ftnlen)6);
	throw runtime_error(vastr("DQRDEC ",info));
	return 0;
    }
/* delete the j-th column. */
    i__1 = *n - 1;
    for (i__ = *j; i__ <= i__1; ++i__) {
	dcopy_(k, &r__[(i__ + 1) * r_dim1 + 1], &c__1, &r__[i__ * r_dim1 + 1],
		 &c__1);
    }
/* retriangularize. */
    if (*j < *k) {
   	i__1 = *k + 1 - *j;
   	i__2 = *n - *j;
   	dqhqr_(&i__1, &i__2, &r__[*j + *j * r_dim1], ldr,
				&w[1], &r__[*n *r_dim1 + 1]);
/* apply rotations to Q. */
   	i__1 = std::min(*k,*n) + 1 - *j;
	   dqrot_("F", m, &i__1, &q[*j * q_dim1 + 1], ldq,
				&w[1], &r__[*n *r_dim1 + 1]);
    }
    return 0;
} /* dqrdec_ */

int
dqrder_(integer* m, integer* n, doublereal* q, integer* ldq,
		doublereal* r__, integer* ldr, integer* j, doublereal* w) {
    /* System generated locals */
    integer q_dim1, q_offset, r_dim1, r_offset, i__1, i__2;

    /* Local variables */
    int i__, k, info;
	// Table of constant values
	int c__1 = 1;

/* purpose:      updates a QR factorization after deleting a row. */
/*               i.e., given an m-by-m orthogonal matrix Q, an m-by-n */
/*               upper trapezoidal matrix R and index j in the range */
/*               1:m, this subroutine updates Q ->Q1 and an R -> R1 */
/*               so that Q1 is again orthogonal, R1 upper trapezoidal, */
/*               and Q1*R1 = [A(1:j-1,:); A(j+1:m,:)], where A = Q*R. */
/*               (real version) */

/* arguments: */
/* m (in)        number of rows of the matrix Q. */
/* n (in)        number of columns of the matrix R. */
/* Q (io)        on entry, the orthogonal matrix Q. */
/*               on exit, the updated matrix Q1. */
/* ldq (in)      leading dimension of Q. ldq >= m. */
/* R (io)        on entry, the original matrix R. */
/*               on exit, the updated matrix R1. */
/* ldr (in)      leading dimension of R. ldr >= m. */
/* j (in)        the position of the deleted row. */
/* w (out)       a workspace vector of size 2*m. */

/* quick return if possible */
    /* Parameter adjustments */
    q_dim1 = *ldq;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    r_dim1 = *ldr;
    r_offset = 1 + r_dim1;
    r__ -= r_offset;
    --w;

    /* Function Body */
    if (*m == 1) {
	return 0;
    }
/* check arguments */
    info = 0;
    if (*m < 1) {
	info = 1;
    } else if (*j < 1 || *j > *m) {
	info = 7;
    }
    if (info != 0) {
	// xerbla_("DQRDER", &info, (ftnlen)6);
	throw runtime_error(vastr("DQRDER ",info));
	return 0;
    }
/* eliminate Q(j,2:m). */
    dcopy_(m, &q[*j + q_dim1], ldq, &w[1], &c__1);
    dqrtv1_(m, &w[1], &w[*m + 1]);
/* apply rotations to Q. */
    dqrot_("B", m, m, &q[q_offset], ldq, &w[*m + 1], &w[2]);
/* form Q1. */
    i__1 = *m - 1;
    for (k = 1; k <= i__1; ++k) {
	if (*j > 1) {
	    i__2 = *j - 1;
	    dcopy_(&i__2, &q[(k + 1) * q_dim1 + 1], &c__1, &q[k * q_dim1 + 1],
		     &c__1);
	}
	if (*j < *m) {
	    i__2 = *m - *j;
	    dcopy_(&i__2, &q[*j + 1 + (k + 1) * q_dim1], &c__1, &q[*j + k * 
		    q_dim1], &c__1);
	}
    }
/* apply rotations to R. */
    dqrqh_(m, n, &r__[r_offset], ldr, &w[*m + 1], &w[2]);
/* form R1. */
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	i__2 = *m - 1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    r__[i__ + k * r_dim1] = r__[i__ + 1 + k * r_dim1];
	}
    }
    return 0;
} /* dqrder_ */

int
dqrinc_(integer *m, integer *n, integer *k, doublereal *
	q, integer *ldq, doublereal *r__, integer *ldr, integer *j, 
	doublereal *x, doublereal *w) {
    /* System generated locals */
    integer q_dim1, q_offset, r_dim1, r_offset, i__1, i__2;
    doublereal d__1;

    /* Local variables */
    int i__, k1;
    double rx;
    int info;
    int full;
	// Table of constant values
	int c__1 = 1;

/* purpose:      updates a QR factorization after inserting a new */
/*               column. */
/*               i.e., given an m-by-k orthogonal matrix Q, an m-by-n */
/*               upper trapezoidal matrix R and index j in the range */
/*               1:n+1, this subroutine updates the matrix Q -> Q1 and */
/*               R -> R1 so that Q1 is again orthogonal, R1 upper */
/*               trapezoidal, and Q1*R1 = [A(:,1:j-1); x; A(:,j:n)], */
/*               where A = Q*R. */
/*               (real version) */
/* arguments: */
/* m (in)        number of rows of the matrix Q. */
/* n (in)        number of columns of the matrix R. */
/* k (in)        number of columns of Q, and rows of R. Must be */
/*               either k = m (full Q) or k = n <= m (economical form, */
/*               basis dimension will increase). */
/* Q (io)        on entry, the orthogonal m-by-k matrix Q. */
/*               on exit, the updated matrix Q1. */
/* ldq (in)      leading dimension of Q. ldq >= m. */
/* R (io)        on entry, the original matrix R. */
/*               on exit, the updated matrix R1. */
/* ldr (in)      leading dimension of R. ldr >= min(m,n+1). */
/* j (in)        the position of the new column in R1 */
/* x (in)        the column being inserted */
/* w (out)       a workspace vector of size k. */

/* quick return if possible. */
    /* Parameter adjustments */
    q_dim1 = *ldq;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    r_dim1 = *ldr;
    r_offset = 1 + r_dim1;
    r__ -= r_offset;
    --x;
    --w;

    /* Function Body */
    if (*m == 0) {
	return 0;
    }
/* check arguments. */
    info = 0;
    if (*m < 0) {
	info = 1;
    } else if (*n < 0) {
	info = 2;
    } else if (*k != *m && (*k != *n || *n >= *m)) {
	info = 3;
    } else if (*ldq < *m) {
	info = 5;
    } else /* if(complicated condition) */ {
/* Computing MIN */
	i__1 = *m, i__2 = *k + 1;
	if (*ldr < std::min(i__1,i__2)) {
	    info = 7;
	} else if (*j < 1 || *j > *n + 1) {
	    info = 8;
	}
    }
    if (info != 0) {
	// xerbla_("DQRINC", &info, (ftnlen)6);
	throw runtime_error(vastr("DQRINC ",info));
	return 0;
    }
	 if (*k == *m)
		 full = 1;
	 else
		 full = 0;
    // full = *k == *m;
/* insert empty column at j-th position. */
    i__1 = *j;
    for (i__ = *n; i__ >= i__1; --i__) {
	dcopy_(k, &r__[i__ * r_dim1 + 1], &c__1, &r__[(i__ + 1) * r_dim1 + 1],
		 &c__1);
    }
/* insert Q'*u into R. In the nonfull case, form also u-Q*Q'*u. */
    if (full) {
	k1 = *k;
	i__1 = *k;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    r__[i__ + *j * r_dim1] = ddot_(m, &q[i__ * q_dim1 + 1], &c__1, &x[
		    1], &c__1);
	}
    } else {
	k1 = *k + 1;
/* zero last row of R */
	i__1 = *n + 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    r__[k1 + i__ * r_dim1] = 0.;
	}
	dcopy_(m, &x[1], &c__1, &q[k1 * q_dim1 + 1], &c__1);
	i__1 = *k;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    r__[i__ + *j * r_dim1] = ddot_(m, &q[i__ * q_dim1 + 1], &c__1, &q[
		    k1 * q_dim1 + 1], &c__1);
	    d__1 = -r__[i__ + *j * r_dim1];
	    daxpy_(m, &d__1, &q[i__ * q_dim1 + 1], &c__1, &q[k1 * q_dim1 + 1],
		     &c__1);
	}
/* get norm of the inserted column */
	rx = dnrm2_(m, &q[k1 * q_dim1 + 1], &c__1);
	r__[k1 + *j * r_dim1] = rx;
	if (rx == 0.) {
/* in the rare case when rx is exact zero, we still need to provide */
/* a valid orthogonal unit vector. The details are boring, so handle */
/* that elsewhere. */
	    dgqvec_(m, k, &q[q_offset], ldq, &q[k1 * q_dim1 + 1]);
	} else {
/* otherwise, just normalize the added column. */
	    d__1 = 1. / rx;
	    dscal_(m, &d__1, &q[k1 * q_dim1 + 1], &c__1);
	}
    }
/* maybe we're finished. */
    if (*j > *k) {
	return 0;
    }
/* eliminate the spike. */
    i__1 = k1 + 1 - *j;
    dqrtv1_(&i__1, &r__[*j + *j * r_dim1], &w[1]);
/* apply rotations to R(j:k,j:n). */
    if (*j <= *n) {
	i__1 = k1 + 1 - *j;
	i__2 = *n + 1 - *j;
	dqrqh_(&i__1, &i__2, &r__[*j + (*j + 1) * r_dim1], ldr, &w[1], &r__[*
		j + 1 + *j * r_dim1]);
    }
/* apply rotations to Q(:,j:k). */
    i__1 = k1 + 1 - *j;
    dqrot_("B", m, &i__1, &q[*j * q_dim1 + 1], ldq, &w[1], &r__[*j + 1 + *j * 
	    r_dim1]);
/* zero spike. */
    i__1 = k1;
    for (i__ = *j + 1; i__ <= i__1; ++i__) {
	r__[i__ + *j * r_dim1] = 0.;
    }
    return 0;
} /* dqrinc_ */

int
dqrinr_(integer *m, integer *n, doublereal *q, integer *
	ldq, doublereal *r__, integer *ldr, integer *j, doublereal *x, 
	doublereal *w) {
    /* System generated locals */
    integer q_dim1, q_offset, r_dim1, r_offset, i__1, i__2;

    /* Local variables */
    int i__, k, info;
	// Table of constant values
	int c__1 = 1;

/* purpose:      updates a QR factorization after inserting a new */
/*               row. */
/*               i.e., given an m-by-m unitary matrix Q, an m-by-n */
/*               upper trapezoidal matrix R and index j in the range */
/*               1:m+1, this subroutine updates Q -> Q1  and R -> R1 */
/*               so that Q1 is again unitary, R1 upper trapezoidal, */
/*               and Q1*R1 = [A(1:j-1,:); x; A(j:m,:)], where A = Q*R. */
/*               (real version) */
/* arguments: */
/* m (in)        number of rows of the matrix Q. */
/* n (in)        number of columns of the matrix R. */
/* Q (io)        on entry, the unitary matrix Q. */
/*               on exit, the updated matrix Q1. */
/* ldq (in)      leading dimension of Q. ldq >= m+1. */
/* R (io)        on entry, the original matrix R. */
/*               on exit, the updated matrix R1. */
/* ldr (in)      leading dimension of R. ldr >= m+1. */
/* j (in)        the position of the new row in R1 */
/* x (io)        on entry, the row being added */
/*               on exit, x is destroyed. */
/* w (out)       a workspace vector of size min(m,n). */

/* check arguments */
    /* Parameter adjustments */
    q_dim1 = *ldq;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    r_dim1 = *ldr;
    r_offset = 1 + r_dim1;
    r__ -= r_offset;
    --x;
    --w;

    /* Function Body */
    info = 0;
    if (*n < 0) {
	info = 2;
    } else if (*j < 1 || *j > *m + 1) {
	info = 7;
    }
    if (info != 0) {
	// xerbla_("DQRINR", &info, (ftnlen)6);
	throw runtime_error(vastr("DQRINR ",info));
	return 0;
    }
/* permute the columns of Q1 and rows of R1 so that c the new row ends */
/* up being the topmost row of R1. */
    for (i__ = *m; i__ >= 1; --i__) {
	if (*j > 1) {
	    i__1 = *j - 1;
	    dcopy_(&i__1, &q[i__ * q_dim1 + 1], &c__1, &q[(i__ + 1) * q_dim1 
		    + 1], &c__1);
	}
	q[*j + (i__ + 1) * q_dim1] = 0.;
	if (*j <= *m) {
	    i__1 = *m + 1 - *j;
	    dcopy_(&i__1, &q[*j + i__ * q_dim1], &c__1, &q[*j + 1 + (i__ + 1) 
		    * q_dim1], &c__1);
	}
    }
/* set up the 1st column */
    i__1 = *j - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	q[i__ + q_dim1] = 0.;
    }
    q[*j + q_dim1] = 1.;
    i__1 = *m + 1;
    for (i__ = *j + 1; i__ <= i__1; ++i__) {
	q[i__ + q_dim1] = 0.;
    }
/* set up the new matrix R1 */
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	if (k < *m) {
	    r__[*m + 1 + k * r_dim1] = 0.;
	}
	for (i__ = std::min(*m,k); i__ >= 1; --i__) {
	    r__[i__ + 1 + k * r_dim1] = r__[i__ + k * r_dim1];
	}
	r__[k * r_dim1 + 1] = x[k];
    }
/* retriangularize R */
    i__1 = *m + 1;
    dqhqr_(&i__1, n, &r__[r_offset], ldr, &w[1], &x[1]);
/* apply rotations to Q */
    i__1 = *m + 1;
    i__2 = std::min(*m,*n) + 1;
    dqrot_("F", &i__1, &i__2, &q[q_offset], ldq, &w[1], &x[1]);
    return 0;
} /* dqrinr_ */

int
dqrot_(const char *dir, integer *m, integer *n, doublereal *q, 
	integer *ldq, doublereal *c__, doublereal *s) {
    /* System generated locals */
    integer q_dim1, q_offset, i__1;

    /* Local variables */
    int i__;
    int fwd;
    int info;
	// Table of constant values
	int c__1 = 1;

/* purpose:      Apply a sequence of inv. rotations from right */

/* arguments: */
/* dir (in)      if 'B' or 'b', rotations are applied from backwards */
/*               if 'F' or 'f', from forwards. */
/* m (in)        number of rows of matrix Q */
/* n (in)        number of columns of the matrix Q */
/* Q (io)        on entry, the matrix Q */
/*               on exit, the updated matrix Q1 */
/* ldq (in)      the leading dimension of Q */
/* c (in)        n-1 rotation cosines */
/* s (in)        n-1 rotation sines */

/* quick return if possible */
    /* Parameter adjustments */
    q_dim1 = *ldq;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    --c__;
    --s;

    /* Function Body */
    if (*m == 0 || *n == 0 || *n == 1) {
	return 0;
    }
/* check arguments. */
    info = 0;
    // fwd = lsame_(dir, "F", (ftnlen)1, (ftnlen)1);
	 fwd = 0;
	 if (tolower(dir[0]) == 'f')
		fwd = 1;
    // if (! (fwd || lsame_(dir, "B", (ftnlen)1, (ftnlen)1)))
	// info = 1;
	if (!(fwd == 1 || tolower(dir[0]) == 'b')) {
		info = 1;
    } else if (*m < 0) {
	info = 2;
    } else if (*n < 0) {
	info = 3;
    } else if (*ldq < *m) {
	info = 5;
    }
    if (info != 0) {
	// xerbla_("DQROT", &info, (ftnlen)5);
	throw runtime_error(vastr("DQROT ",info));
	return 0;
    }
    if (fwd) {
	i__1 = *n - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    drot_(m, &q[i__ * q_dim1 + 1], &c__1, &q[(i__ + 1) * q_dim1 + 1], 
		    &c__1, &c__[i__], &s[i__]);
	}
    } else {
	for (i__ = *n - 1; i__ >= 1; --i__) {
	    drot_(m, &q[i__ * q_dim1 + 1], &c__1, &q[(i__ + 1) * q_dim1 + 1], 
		    &c__1, &c__[i__], &s[i__]);
	}
    }
    return 0;
} /* dqrot_ */

int
dqrqh_(integer *m, integer *n, doublereal *r__, integer *
	ldr, doublereal *c__, doublereal *s) {
    /* System generated locals */
    integer r_dim1, r_offset, i__1, i__2;

    /* Local variables */
    int i__, j;
    double t;
    int ii, info;

/* purpose:      brings an upper trapezoidal matrix R into upper */
/*               Hessenberg form using min(m-1,n) Givens rotations. */
/*               (real version) */
/* arguments: */
/* m (in)        number of rows of the matrix R */
/* n (in)        number of columns of the matrix R */
/* R (io)        on entry, the upper Hessenberg matrix R */
/*               on exit, the updated upper trapezoidal matrix */
/* ldr (in)      leading dimension of R, >= m */
/* c(in)         rotation cosines, size at least min(m-1,n) */
/* s(in)         rotation sines, size at least min(m-1,n) */

/* quick return if possible. */
    /* Parameter adjustments */
    r_dim1 = *ldr;
    r_offset = 1 + r_dim1;
    r__ -= r_offset;
    --c__;
    --s;

    /* Function Body */
    if (*m == 0 || *m == 1 || *n == 0) {
	return 0;
    }
/* check arguments. */
    info = 0;
    if (*m < 0) {
	info = 1;
    } else if (*n < 0) {
	info = 2;
    } else if (*ldr < *m) {
	info = 4;
    }
    if (info != 0) {
	// xerbla_("DQRQH", &info, (ftnlen)5);
	throw runtime_error(vastr("DQRQH ",info));
	return 0;
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MIN */
	i__2 = *m - 1;
	ii = std::min(i__2,i__);
/* apply stored rotations, column-wise */
	t = r__[ii + 1 + i__ * r_dim1];
	for (j = ii; j >= 1; --j) {
	    r__[j + 1 + i__ * r_dim1] = c__[j] * t - s[j] * r__[j + i__ * 
		    r_dim1];
	    t = c__[j] * r__[j + i__ * r_dim1] + s[j] * t;
	}
	r__[i__ * r_dim1 + 1] = t;
    }
    return 0;
} /* dqrqh_ */

int
dqrshc_(integer *m, integer *n, integer *k, doublereal *
	q, integer *ldq, doublereal *r__, integer *ldr, integer *i__, integer 
	*j, doublereal *w) {
    /* System generated locals */
    integer q_dim1, q_offset, r_dim1, r_offset, i__1, i__2;

    /* Local variables */
    int l, jj, kk, info;
	// Table of constant values
	int c__1 = 1;

/* purpose:      updates a QR factorization after circular shift of */
/*               columns. */
/*               i.e., given an m-by-k orthogonal matrix Q, an k-by-n */
/*               upper trapezoidal matrix R and index j in the range */
/*               1:n+1, this subroutine updates the matrix Q -> Q1 and */
/*               R -> R1 so that Q1 is again orthogonal, R1 upper */
/*               trapezoidal, and */
/*               Q1*R1 = A(:,p), where A = Q*R and p is the permutation */
/*               [1:i-1,shift(i:j,-1),j+1:n] if i < j  or */
/*               [1:j-1,shift(j:i,+1),i+1:n] if j < i. */
/*               (real version) */
/* arguments: */
/* m (in)        number of rows of the matrix Q. */
/* n (in)        number of columns of the matrix R. */
/* k (in)        number of columns of Q1, and rows of R1. Must be */
/*               either k = m (full Q) or k = n <= m (economical form). */
/* Q (io)        on entry, the unitary m-by-k matrix Q. */
/*               on exit, the updated matrix Q1. */
/* ldq (in)      leading dimension of Q. ldq >= m. */
/* R (io)        on entry, the original matrix R. */
/*               on exit, the updated matrix R1. */
/* ldr (in)      leading dimension of R. ldr >= k. */
/* i (in)        the first index determining the range (see above) */
/* j (in)        the second index determining the range (see above) */
/* w (o)         a workspace vector of size 2*k. */

/* quick return if possible. */
    /* Parameter adjustments */
    q_dim1 = *ldq;
    q_offset = 1 + q_dim1;
    q -= q_offset;
    r_dim1 = *ldr;
    r_offset = 1 + r_dim1;
    r__ -= r_offset;
    --w;

    /* Function Body */
    if (*m == 0 || *n == 1) {
	return 0;
    }
    info = 0;
/* check arguments. */
    if (*m < 0) {
	info = 1;
    } else if (*n < 0) {
	info = 2;
    } else if (*k != *m && (*k != *n || *n > *m)) {
	info = 3;
    } else if (*i__ < 1 || *i__ > *n) {
	info = 6;
    } else if (*j < 1 || *j > *n) {
	info = 7;
    }
    if (info != 0) {
	// xerbla_("DQRSHC", &info, (ftnlen)6);
	throw runtime_error(vastr("DQRSHC ",info));
	return 0;
    }
    if (*i__ < *j) {
/* shift columns */
	dcopy_(k, &r__[*i__ * r_dim1 + 1], &c__1, &w[1], &c__1);
	i__1 = *j - 1;
	for (l = *i__; l <= i__1; ++l) {
	    dcopy_(k, &r__[(l + 1) * r_dim1 + 1], &c__1, &r__[l * r_dim1 + 1],
		     &c__1);
	}
	dcopy_(k, &w[1], &c__1, &r__[*j * r_dim1 + 1], &c__1);
/* retriangularize */
	if (*i__ < *k) {
	    kk = std::min(*k,*j);
	    i__1 = kk + 1 - *i__;
	    i__2 = *n + 1 - *i__;
	    dqhqr_(&i__1, &i__2, &r__[*i__ + *i__ * r_dim1], ldr, &w[*k + 1], 
		    &w[1]);
/* apply rotations to Q. */
	    i__1 = kk + 1 - *i__;
	    dqrot_("F", m, &i__1, &q[*i__ * q_dim1 + 1], ldq, &w[*k + 1], &w[1]);
	}
    } else if (*j < *i__) {
/* shift columns */
	dcopy_(k, &r__[*i__ * r_dim1 + 1], &c__1, &w[1], &c__1);
	i__1 = *j + 1;
	for (l = *i__; l >= i__1; --l) {
	    dcopy_(k, &r__[(l - 1) * r_dim1 + 1], &c__1, &r__[l * r_dim1 + 1],
		     &c__1);
	}
	dcopy_(k, &w[1], &c__1, &r__[*j * r_dim1 + 1], &c__1);
/* retriangularize */
	if (*j < *k) {
/* Computing MIN */
	    i__1 = *j + 1;
	    jj = std::min(i__1,*n);
	    kk = std::min(*k,*i__);
/* eliminate the introduced spike. */
	    i__1 = kk + 1 - *j;
	    dqrtv1_(&i__1, &r__[*j + *j * r_dim1], &w[*k + 1]);
/* apply rotations to R */
	    i__1 = kk + 1 - *j;
	    i__2 = *n - *j;
	    dqrqh_(&i__1, &i__2, &r__[*j + jj * r_dim1], ldr, &w[*k + 1], &
		    r__[*j + 1 + *j * r_dim1]);
/* apply rotations to Q */
	    i__1 = kk + 1 - *j;
	    dqrot_("B", m, &i__1, &q[*j * q_dim1 + 1], ldq, &w[*k + 1], &r__[*
		    j + 1 + *j * r_dim1]);
/* zero spike. */
	    i__1 = kk;
	    for (l = *j + 1; l <= i__1; ++l) {
		r__[l + *j * r_dim1] = 0.;
	    }
	}
    }
    return 0;
} /* dqrshc_ */

int
dqrtv1_(integer *n, doublereal *u, doublereal *w) {
    int i__;
    double t, rr;

/* purpose:      generates a sequence of n-1 Givens rotations that */
/*               eliminate all but the first element of a vector u. */
/* arguments: */
/* n (in)        the length of the vector u */
/* u (io)        on entry, the vector u. */
/*               on exit, u(2:n) contains the rotation sines, u(1) */
/*               contains the remaining element. */
/* w (o)         on exit, w contains the rotation cosines. */

/* quick return if possible. */
    /* Parameter adjustments */
    --w;
    --u;

    /* Function Body */
    if (*n <= 0) {
	return 0;
    }
    rr = u[*n];
    for (i__ = *n - 1; i__ >= 1; --i__) {
	dlartg_(&u[i__], &rr, &w[i__], &u[i__ + 1], &t);
	rr = t;
    }
    u[1] = rr;
    return 0;
} /* dqrtv1_ */

#ifdef MAIN
#undef MAIN

#include "matrix.h"

int
main (int argc, char **argv) {
	int rval{0};
	int nr{0}, nc{0};
	int info{0};
	bool is_complex;

	switch (argc) {
		case 1: {
			cerr << "Usage: " << argv[0] << " file\n";
			cerr << "        the file extension determines the format\n";
			rval = 1;
			break;
		}
		case 2: {
			vector<double> a = MM::importer(argv[1], nr, nc, is_complex);
			vector<double> w(nr,0.0);
			vector<double> A{a}; // copy
			// eigenvectors are returned in A
			info = lapack::dsyev("v","u",nr,&A[0],nr,&w[0]);
			if (info != 0) {
				cerr << "dsyev failed: info = " << info << endl;
				rval = 1;
			} else {
				cout << "eigenvalues:\n" << w << endl;
				// test triprod
				vector<double> tat(nr*nr, 0.0);
				triprod(nr, nr, &A[0], &a[0], &tat[0]);
				Dprint::dprintm(nr,nr,nr,&tat[0],"TAT");
			}
			break;
		}
		case 3: {
			cerr << "not implemented\n";
			break;
		}
		case 4: {
			vector<double> k = MM::importer(argv[1], nr, nc, is_complex);
			vector<double> g = MM::importer(argv[2], nr, nc, is_complex);
			// scale g by spin
			double spin{100.0};
			blas::scal(nr*nr, spin, &g[0], 1);
			vector<double> m = MM::importer(argv[3], nr, nc, is_complex);
			vector<double> w(nr,0.0);
			vector<complex<double>> x(nr*nr,0.0);
			vector<int> pi(nr,0);
			info = lapack::gyroeig (nr, m, g, k, w, x, pi);
			if (info != 0) {
				cerr << "gyroeig failed: info = " << info << endl;
				rval = 1;
			} else {
				cout << w.size() << " eigenvalues:\n" << w << endl;
			}
			break;
		}
		default:
			cerr << "usage: " << argv[0] << " stif.mm gyro.mm mass.mm\n";
	}
	return rval;
}
#endif // MAIN
