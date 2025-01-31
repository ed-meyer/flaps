//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
#ifndef lapack_h
#define lapack_h

#include <complex>
#include <limits>
#include <vector>

#include "blas.h"

// C/C++ covers for LAPACK fortran subroutines:
//   - these take int instead of int*
//   - workspace arguments are left off - optimal workspace
//     determined internally

// dqr1up is not part of lapack but modifies QR factorization
// from lapack
int dqr1up(int m, int n, int k, double* q, int ldq,
		double* r, int ldr, double* u, double* v);

namespace lapack {

// Cholesky factor & solve
int dposv(const char* uplo, int n, int nrhs,double* a, int lda, double* b, int ldb);
// Cholesky factorization
int dpotrf(const char* uplo, int n, double* a, int lda);
// solution using Cholesky factorization
int dpotrs(const char* uplo, int n, int nrhs, double* a,
		int lda, double* b, int ldb);

// LU factorization of nonsymmetric matrix
int dgetrf (int m, int n, double* a, int lda, int* ipiv);
// linear equation drivers
int dgesv (int n, int nrhs, double* a, int lda, int* ipiv, double* b, int ldb);
int dgesvx (char const* fact, char const* trans, int n, int nrhs,
	double* a, int lda, double* af, int ldaf, int* ipiv, 
	char equed[2], double* r, double* c, double* b, int ldb,
	double* x, int ldx, double* rcond, double* ferr, double* berr);

int zgesvx (char const* fact, char const* trans, int n, int nrhs,
	std::complex<double>* a, int lda, std::complex<double>* af, int ldaf, int* ipiv, 
	char equed[2], double* r, double* c, std::complex<double>* b, int ldb,
	std::complex<double>* x, int ldx, double* rcond, double* ferr, double* berr);

int dgeqrf (int m, int n, double* a, int lda, double* tau);
int dgeqp3 (int m, int n, double* a, int lda, int* jpvt, double* tau);

double dtrcon(char const* norm, char const* uplo, char const* diag,
		int n, double* a, int lda);

int dsyev(char const* jobz, char const* uplo, int n, double* a, int lda, double* w);

// symmetric generalized eigenvalues
int dsygv(int itype, const char* jobz, const char* uplo, int n, double* a,
		int lda, double* b, int ldb, double* w);

int dgeev (char const* jobvl, char const* jobvr, int n, double* a, int lda,
	double* wr, double* wi, double* vl, int ldvl, double* vr, int ldvr);

int dgeevx (char const* balanc, char const* jobvl, char const* jobvr,
		char const* sense, int n, double* a, int lda,
	double* wr, double* wi, double* vl, int ldvl, double* vr, int ldvr,
	int& ilo, int& ihi, double* scale, double& abnrm, double* rconde, double* rcondv);

int dormqr (char const* side, char const* trans, int m, int n,
		int k, double const* a, int lda, double const* tau, double* c, int ldc);

void dtrsm (char const* side, char const* uplo, char const* transa,
		char const* diag, int m, int n, double alpha, double *a,
		int lda, double *b, int ldb);

int dgesvd (char const* jobu, char const* jobvt, int m, int n, double* a,
		int lda, double* s, double* u, int ldu, double* vt, int ldvt);

int dggev (char const* jobvl, char const* jobvr, int n, double* a, int lda,
		double* b, int ldb, double* alphar, double* alphai, double* beta,
		double* vl, int ldvl, double* vr, int ldvr);

int zggev (char const* jobvl, char const* jobvr,
		int n, std::complex<double>* a, int lda,
		std::complex<double>* b, int ldb, std::complex<double>* alpha,
		std::complex<double>* beta, std::complex<double>* vl, int ldvl,
		std::complex<double>* vr, int ldvr);

int zggevx (char const* balanc, char const* jobvl, char const* jobvr,
		char const* sense, int n, std::complex<double>* a, int lda,
		std::complex<double>* b, int ldb, std::complex<double>* alpha,
		std::complex<double>* beta, std::complex<double>* vl, int ldvl,
		std::complex<double>* vr, int ldvr, double* rconde, double* rcondv);


std::string
eigstring (std::complex<double> w, int n, std::complex<double>* vr);
std::string
eigstring (double w, int n, double* vr);

// eigensolution for free-vibration with gyro
int gyroeig (int n, const std::vector<double>& m, const std::vector<double>& g,
		const std::vector<double>& k, std::vector<double>& w,
		std::vector<std::complex<double>>& x, std::vector<int>& pi);

int poly (int n, std::vector<double*>& A, std::vector<std::complex<double>>& w,
		std::vector<std::complex<double>>& x);

int polyeig (size_t n, int maxeig, std::vector<std::complex<double>*>& A,
		std::complex<double>* w, std::complex<double>* vl,
		std::complex<double>* vr, int* pil, int* pir, int nrefine);

} // namespace lapack


// declarations of the Fortran LAPACK routines assuming the
// convention is lower-case with a trailing underscore
extern "C" {

// Cholesky factorization, solution
void dposv_ (const char* uplo, int* n, int* nrhs, double* a,
		int* lda, double* b, int* ldb, int* info);
void dpotrf_ (const char* uplo, int* n, double* a, int* lda, int* info);
void dpotrs_ (const char* uplo, int* n, int* nrhs, double* a,
		int* lda, double* b, int* ldb, int* info);

// LU factorization
void dgetrf_ (int* m, int* n, double* a, int* lda, int* ipiv, int* info);

int dgesv_ (int* n,int*nrhs, double* a, int* lda, int* ipiv, double* b, int* ldb, int* info);
void dgesvx_ (char const* fact, char const* trans,
			int *n, int *nrhs, double *a, int *lda, double *af, int *ldaf,
			int *ipiv, char const* equed, double *r, double *c, double *b,
			int *ldb, double *x, int *ldx, double *rcond, double *ferr,
			double *berr, double *work, int *iwork, int *info);

int zgesvx_ (char const* fact, char const* trans, int* n, int* nrhs,
	std::complex<double>* a, int* lda, std::complex<double>* af, int* ldaf, int* ipiv, 
	const char* equed, double* r, double* c, std::complex<double>* b, int* ldb,
	std::complex<double>* x, int* ldx, double* rcond, double* ferr, double* berr,
	std::complex<double>* work, double* rwork, int* info);

void dgeqrf_ (int *m, int *n, double *a, int *lda, double* tau,
		double* work, int* lwork, int* info);

void dgeqp3_ (int* m, int* n, double* a, int* lda, int* jpvt, double* tau,
		double* work, int* lwork, int* info);

void dtrcon_ (char const* norm, char const* uplo, char const* diag, int* n,
				double* a, int* lda, double* rcond, double* work, int* iwork, int* info);

void dsyev_(char const*, char const*, int* n, double *a, int *lda,
				double* w, double* work, int* lwork, int* info);
// generalized symmetric eigensolution
void dsygv_(int* itype, const char* jobz, const char* uplo, int* n, double* a,
	int* lda, double* b, int* ldb, double* w, double* work, int* lwork, int* info);

void dgeev_(char const*, char const*, int* n, double *a, int *lda,
				double* wr, double* wi, double* vl, int* ldv, double* vr, int* ldvr,
				double* work, int* lwork, int* info);

void dgeevx_(char const*, char const*, char const*, char const*,
		int* n, double *a, int *lda, double* wr, double* wi, double* vl, int* ldv,
		double* vr, int* ldvr, int* ilo, int* ihi, double* scale, double* abnrm,
		double* rconde, double* rcondv, double* work, int* lwork, int* iwork, int* info);

void dggev_(char const* jobvl, char const* jobvr, int* n,
		double* a, int* lda, double* b, int* ldb, double* alphar, double* alphai, double* beta,
		double* vl, int* ldvl, double* vr, int* ldvr, double* work, int* lwork, int* info);

void zggev_(char const* jobvl, char const* jobvr, int *n,
		std::complex<double> *a, int *lda, std::complex<double> *b,
		int *ldb, std::complex<double> *alpha, std::complex<double> *beta,
		std::complex<double> *vl, int *ldvl, std::complex<double> *vr, int *ldvr,
		std::complex<double> *work, int *lwork, double* rwork, int *info);

void zggevx_(char const* balanc, char const* jobvl, char const* jobvr,
				char const* sense, int *n, std::complex<double> *a, int *lda,
				std::complex<double> *b, int *ldb, std::complex<double> *alpha,
				std::complex<double> *beta, std::complex<double> *vl, int *ldvl,
				std::complex<double> *vr, int *ldvr, int *ilo, int *ihi,
				double* lscale, double* rscale, double* abnrm, double* bbnrm,
				double* rconde, double* rcondv, std::complex<double> *work, int *lwork,
				double* rwork, int* iwork, int* bwork, int *info);

void dormqr_ (char const* side, char const* trans, int *m, int *n, int *k,
		double const* a, int *lda, double const* tau, double const* c,
		int *ldc, double *work, int *lwork, int *info);

void dtrsm_ (char const* side, char const* uplo, char const* transa,
		char const* diag, int *m, int *n, double *alpha, double *a,
		int *lda, double *b, int *ldb);

void dgesvd_ (char const* jobu, char const* jobvt, int *m, int *n,
		double *a, int *lda, double *s, double *u, int *ldu, double *vt,
		int *ldvt, double *work, int *lwork, int *info);

}  // extern "C"


// Some extra functions not part of Lapack but we put in lapack
// namespace so it's obvious where to find them
namespace lapack {

// maxpi determines how much refining is done: 100=none, just get pi's
void refine (size_t n, size_t m, std::vector<std::complex<double> const*>& A,
		std::complex<double>* w, std::complex<double>* vl, std::complex<double>* vr, int* pil,
		int* pir, int maxpi, std::complex<double>* uncert);

void sortEigen (size_t n, size_t m, std::complex<double> *s, std::complex<double>* vl, std::complex<double>* vr);
void printEigen (FILE* stream, char const* title, int n, int m,
		std::complex<double>* s, std::complex<double>* x, int* pi);

void plotCeigen(std::string const& path, std::string const& runid, int m, std::complex<double>* eig);

template<typename Type>
int
PerformanceIndex (size_t nr, size_t nc, Type const* A, Type const* x, Type const* b, Type* residual, bool computeResidual) {
/*
 * compute a performance index for the linear equation Ax = b
 * b is optional - ignored if a null pointer so this function can
 * be used either for the linear equation or polynomial eigenvalue
 * problem
 *   Ax = A_0 + sA_1 + s^2A_2 + ... = 0
 */
	Type trmax(0.0);
	Type tol(sqrt(std::numeric_limits<Type>::epsilon()));
	Type eps(sqrt(std::numeric_limits<Type>::epsilon()));
	size_t i, j;

	/*
	 * compute left, right residuals
	 */
	std::vector<Type> rr(nr, Type(0.0));

	if (!residual || computeResidual) {
		for (i=0; i<nr; i++) {
			Type t = Type(0.0);
			for (j=0; j<nc; j++) {
				t += A[i+(j*nr)]*x[j];
			}
			rr[i] = t;
			if (b)
				rr[i] -= b[i];
			if (residual && computeResidual)
				residual[i] = rr[i];
		}
	} else {
		for (i=0; i<nr; i++)
			rr[i] = residual[i];
	}

	/*
	 * compute the size of the interval an exact
	 * residual can lie in and the ratio of the
	 * actual residual to this interval (for each row)
	 */
	for (i=0; i<nr; i++) {
		double rrmin(0.0);
		for (j=0; j<nc; j++) {
			double xj = abs(x[j]);
			if (xj < tol) xj = tol;
			double aij = abs(A[i+(j*nr)]);
			if (aij < tol) aij = tol;
			rrmin += aij*xj;
		}
		//if (b)
		//	rrmin += abs(b[i]);
		if (rrmin < 1.0)
			rrmin = 1.0;
		rrmin *= eps;
		double tr = abs(rr[i])/rrmin;
		if (tr > trmax)
			trmax = tr;
	}

	int rval = 0;
	if (trmax > 1.0)
		rval = (int)log10((double)trmax);
	return rval;
}

template<typename Dtt, typename Dta>
void
triprod(int n, int m, const Dtt* T, const Dta* A, Dta* tat) noexcept {
// triple-product tat = T'AT with T(n,m), A(n,n), tat(m,m)
// assuming Dta*Dtt = Dta e.g. double * complex -> complex
#ifdef NEVER // use dot
	auto IJ = [](int i, int j, int ld1) { return i + j*ld1; };
	for (int i=0; i<m; i++) {
		for (int j=0; j<m; j++) {
			Dta taij{0.0};
			for (int k=0; k<n; k++) {
				Dta atkj{0.0};
				for (int l=0; l<n; l++)
					atkj += A[IJ(k,l,n)]*T[IJ(l,j,n)];
				taij += T[IJ(k,i,n)]*atkj;
			}
			tat[IJ(i,j,m)] = taij;
		}
	}
#else // NEVER // use dot
	std::vector<Dta> at(n);
	for (int i=0; i<m; i++) {
		for (int j=0; j<m; j++) {
			for (int k=0; k<n; k++)
				blas::dot(n, &T[j*n], 1, &A[k], n, at[k]);
			blas::dot(n, &T[i*n], 1, &at[0], 1, tat[i+j*m]);
		}
	}
#endif // NEVER // use dot
}

#ifdef NEVER // no Ad specializations
// specializations of triprod: double, complex<Ad>
template<>
void
triprod(int n, int m, const double* T, const std::complex<Ad>* A, std::complex<Ad>* tat) noexcept;
#endif // NEVER // no Ad specializations

}	// namespace lapack

#endif // lapack_h
