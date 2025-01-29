//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

#include <cassert>
#include <limits>

#include "config.h"
#include "exim.h"
#include "lapack.h"
#include "matrix.h"
#include "qr.h"
#include "svd.h"
#include "trace.h"
#include "transpose.h"

using namespace std;

double
rcondeps() {
// Returns a number to be used to determine rank deficiencies,
// e.g. if an SVD yields s[i]/s[0] that is a rank deficiency
	//!! double eps = 2048.0*std::numeric_limits<double>::epsilon();
	double eps = 256.0*std::numeric_limits<double>::epsilon();
	// double eps = 128.0*std::numeric_limits<double>::epsilon();
	// double eps = 8.0*std::numeric_limits<double>::epsilon();
	// double eps = std::numeric_limits<double>::epsilon();
	return eps;
}

bool
QR::
usesvd(int act) {
	Trace trc(1,"usesvd");
	static bool svd{false};
	bool rval = svd;

	trc.dprint("act ",act);
	if (act > 0) {
		svd = true;
	} else if (act == 0) {
		svd = false;
	}
	return rval;
}

QR::
QR (size_t nrow, size_t ncol, double const* A) {
/*------------------------------------------------------------------
 * QR constructor
 * Factor the (nrow,ncol) real matrix A using an QR factorization
 * of the transposed matrix; this is equivalent to an
 * LQ factorization of the original matrix, except that the QR
 * factorization allows pivoting and complete orthogonal factorization
 * in the rank-deficient case.
 * Solutions using this factorization are minimum 2-norm solutions,
 * as in Algorithm 5.6.2 in
 *   Golub, G.H.; Van Loan, C.F. Matrix Computations, 4th ed; The
 *   Johns Hopkins University Press: Baltimore, MD, USA, 2013
 *
 * Returns a QR structure with:
 *  af       an (ncol,nrow) double array of QR factors
 *  rcond    reciprocal condition number estimate
 *------------------------------------------------------------------*/
	Trace trc(1,"QR constructor (",nrow,",",ncol,")");
	int info;
	double eps = rcondeps();

	// save the dimensions of the original matrix A
	nr = nrow;
	nc = ncol;

	// the user may have requested SVD instead of QR
	if (usesvd()) {
		svd = new SVD(nr, nc, A);
		rcond = svd->rcond;
		rankdef = svd->rankdef;
		return;
	}
	// Since we are factorizing the transpose of A we
	// set m and n to correspond to the LAPACK usage for A'
	// that is, we factorize an (m,n) matrix
	int m = nc, n = nr, lda = nc;

	// transpose the matrix to the af member for factorization
	pivots.resize(n, 0);  // initialize pivots to zero for dgeqp3
	tau.resize(n);     // min(m,n)
	af.resize(m*n);
	flaps::transpose(nr, nc, &A[0], &af[0]);

//!! #define USESGEQP3
#ifdef USESGEQP3
	info = lapack::dgeqp3(m, n, &af[0], lda, &pivots[0], &tau[0]);
#else
	// no pivoting in dgeqrf so simply set them to default
	for (int i=0; i<n; i++) pivots[i] = i+1;
	info = lapack::dgeqrf (m, n, &af[0], lda, &tau[0]);
#endif
	assert(info == 0);
	
	// estimate the reciprocal condition number of the upper
	// triangular factor (the R in QR) - this is the rcond of
	// the original matrix since the condition of Q is 1
	rcond = lapack::dtrcon("I", "U", "N", n, &af[0], lda);
	trc.dprint("rcond ",rcond);

	// if the matrix is rank deficient do an SVD
	rankdef = 0;
	if (rcond < eps) {
		svd = new SVD(nr, nc, A);
		rankdef = svd->rankdef;
		trc.dprint("rcond by qr, svd: ",rcond,", ",svd->rcond);
		MM::exporter("rank-def_A.mm","rank def",A,nr,nc);
	} else {
		svd = nullptr;
	}
}

QR::
~QR() {
	if (svd)
		delete svd;
}

void
QR::
solve (vector<double> const& b, vector<double>& x) {
/*------------------------------------------------------------------
 * Solve a system of linear equations with one right-hand-side
 *
 * Returns an nc-vector x which is the solution to Ax = b,
 * throws runtime_error if an error occured.
 *------------------------------------------------------------------*/
	Trace trc(3,"QR::solve");
	int n = nr;
	size_t i;

	// if the matrix A is rank-deficient there will be an
	// SVD member - use it
	if (svd != nullptr) {
		svd->solve(b, x);
		return;
	}

	// We have an orthog factorization of
	//    A(t)P   = QR
	//    (m,n)(n,n) = (m,m)(n,n)
	// where P is a permutation matrix (pivots).
	// Then
	//    A(t) = QRP(t)
	//    A = PR(t)Q(t)
	// and the min-norm soln to Ax = b is
	//    x = QR(-t)P(t)b
	//    n = (n,m)m
	int m = nc;
	int lda = nc;
	int rank = nr - rankdef;

	// mult P(t)b => x : n
	assert(!pivots.empty());
	for (i=0; i<nr; i++)
		x[i] = b[pivots[i] - 1];

	// mult R(-t)*x => x
	double alpha{1.0};
	lapack::dtrsm("Left", "Upper", "Trans", "Non-unit",
		rank, 1, alpha, &af[0], lda, &x[0], m);

	// zero last m - rank ?
	for (i = rank; i < nc; i++)
		x[i] = 0.0;

	// mult Q*x => x
	int info = lapack::dormqr("Left", "No transpose", m, 1, n,
				&af[0], lda, &tau[0], &x[0], m);
	assert(info == 0);
}

int
QR::
performanceIndex(const double* A, const double* x, const double* b) {
// compute a performance index for the solution
	int rval = lapack::PerformanceIndex(nr, nc, A, x, b, (double*)nullptr, false);
	return rval;
}


vector<double>
QR::
nullProj (vector<double> const& proj, double scale) {
	Trace trc(2,"QR::nullProj");
	// The orig matrix is (nr,nc), the transposed matrix is (m,n)
	int m = nc;
	size_t i;

	// if the matrix A is rank-deficient there will be an
	// SVD member - use it
	if (svd != nullptr) {
		return svd->nullProj(proj);
	}

	vector<double> rval(nc, 0.0);

	// For square matrices we assume that the jacobian has been
	// augmented by adding the tangent as the last row so we
	// solve Jt = I_n i.e. solve for the null vector who's inner
	// product with the tangent is 1.
	// For rectangular matrices which were factored with QR
	// the transpose of the matrix was factored:
	//    A^t = QR
	// so we want to project "proj" onto the orthogonal complement
	// of the column space of A^t; the projection onto the column
	// space is the solution to A^tx = b and the orthogonal complement
	// is A^tx - b or the residual of the least squares problem
	// where b = proj, p is the projection onto the nullspace:
	//    p = A^tx - b
	//    x = R^{-1}Q^tb
	//    p = QRR^{-1}Q^tb - b = QQ^tb - b
	// p: nc vector
	// Note: we assume that the matrix has not been equilibrated!
	if (nr == nc) {
		vector<double> rhs(nr, 0.0);
		rhs[nr-1] = 1.0;
		solve(rhs, rval);
	} else {
		// initialize rval = proj
		for(i=0; i<nc; i++)
			rval[i] = proj[i];
		assert (nc > nr);
		assert (!af.empty());
		assert (!tau.empty());
		trc.dprint("projecting: ",rval);

		int rank = nr - rankdef;

		// Solving for the residual of the least-squares problem A^tx = b
		// where b = proj (see Lapack section 2.4.2.1)
		// First c = (c_1 c_2)^t = Q^t b  with unblocked Q'rval
		int info = lapack::dormqr("Left", "Transpose", m, 1, rank,
				&af[0], m, &tau[0], &rval[0], m);
		assert(info == 0);

		// zero the first rank elements of rval so that rval = (0 c_2)^t
		for (i=0; i < (size_t)rank; i++)
			rval[i] = 0.0;

		// ... then r = Q (0 c_2)^t = (Q_2 Q_2)^t b i.e. the projection
		// of b onto the nullspace of A
		info = lapack::dormqr("Left", "No Transpose", m, 1, rank,
				&af[0], m, &tau[0], &rval[0], m);
		assert(info == 0);
	}

	// check that the projection is in the same direction as the
	// projected vector
	blas_dot (m, &rval[0], 1, &proj[0], 1, scale);
	if (scale < 0.0) {
		scale = -1.0;
		blas_scal (m, scale, &rval[0], 1);
	}
	trc.dprint("returning ",rval);
	return rval;
}

double
QR::
determinant (double& expon) const {
// Given the QR factorization of the transpose of an (nr,nc) matrix,
// compute the determinant as the product of the diagonals of R
// times the determinant of Q (+/- 1)
// Returns the determinant in the form d*10^expon
//
// det(Q) is (-1)^p where p is the number
// of Householder reflections: the number of non-zero elements in tau
// See http://icl.cs.utk.edu/lapack-forum/viewtopic.php?f=2&t=1741
//
// http://icl.cs.utk.edu/lapack-forum/viewtopic.php?f=2&t=1741
// Re: determinant in QR factorization
// Postby Julien Langou B; Sat Feb 06, 2010 2:01 pm
// 
// For a QR factorization of an m-by-n matrix with m <= n, you have IN GENERAL
// n-1 Householder reflections (all the columns except the last), with m > n,
// you have IN GENERAL n Householder reflections. And that could be as simple as that.
// 
// But it happens sometimes that we do not need to a Householder reflection for a
// given column (if it is already in the correct form). For LAPACK 3.1 and before,
// this was when the updated i-th column of A is already the column of the
// triangular matrix R (i.e. the zeros are already here, no need to apply a
// Householder reflection). For LAPACK 3.2 (and after), this is when the updated
// i-th column of A is already the column of the triangular matrix R (i.e. the
// zeros are already here) AND R(i,i) is real nonnegative (if we have the zeros
// but R(i,i) is negative (or has a nonzero imaginary part), we would do a
// reflection to get a real nonnegative diagonal element for R(i,i)).
// 
// All this to say that you should scan TAU from 1 to n (m > n ) or 1 to n-1
// (m <= n) and count the nonzero elements. Say you get k nonzeros in TAU. Then (-1)^k.
// 
	Trace trc(2,"QR::determinant");
	double rval = 1.0;
	double exp = 0.0;

	// if an svd was used return 1
	if (svd != nullptr) {
		trc.dprint("quick return: svd");
		return rval;
	}

	trc.dprint("tau: ",tau);

	int nref = 0;
	for (size_t i=0; i<nr; i++) {
		if (tau[i] != 0.0)
			nref++;
		rval = af[IJ(i,i,nc)]*rval;
		//!! if (pivots[i] != (int)(i+1))
		//!! 	rval = -rval;
		if (abs(rval) == 0.0) break;
		// adjust exponent so that 1 <= det <= 10
		while (abs(rval) < 1.0) {
			rval *= 10.0;
			exp -= 1.0;
		}
		while (abs(rval) > 10.0) {
			rval /= 10.0;
			exp += 1.0;
		}
	}

	trc.dprint(nref," Householder reflections");

	if (nref%2) {
		rval = -rval;
		trc.dprint("odd-number of nonzero taus: changed to ",rval);
	}

	trc.dprint("unscaled determinant: ",rval," 10^",exp);

	expon = exp;
	return rval;
}

double
QR::
augJacDet(vector<double> const& jac, vector<double> const& tan,
		double& expon) const {
// Special function for monitoring the sign of the determinant of
// the Jacobian matrix augmented with the associated tangent vector.
// See notes in bifurcation().
// Compute the determinant of the QR factorized matrix
// augmented with "tan", iff:
// 1) nc == nr+1
// 2) tan is not empty
	Trace trc(2,"augJacDet");
	double rval = 0.0;

	expon = 0.0;

	// if the matrix A is rank-deficient there will be an
	// SVD member, and the determinant must be computed
	// with eigenvalues
	if (svd != nullptr) {
		return augeigdet(nr, nc, jac, tan, expon);
	}

	// quick return if nc != nr+1...
	if (nc != nr+1) {
		trc.dprint("returning 0: nc ",nc," != nr+1 ",nr+1);
		return 0.0;
	}
	// ... or if a tangent wasn't provided
	if (tan.empty()) {
		trc.dprint("returning 0: no previous tangent");
		return 0.0;
	}
	assert(tan.size() == nc);

	// af is (nc,nr) and contains R in it's upper triangle;
	// calculate its determinant
	rval = this->determinant(expon);

	// nc == nr+1 and we are expected to calculate the determinant
	// of the square matrix [ J^t t ], the Jacobian augmented by
	// the tangent. Because the tan was derived from the null vectors
	// the product of the Q matrix and the tan must be +/- 1
	// If it is -1 change the determinant sign
	// copy tan to a temp array (c), overwritten by dormqr
	vector<double> c(tan);
	int m = nc;
	int n = nr;
	int info = lapack::dormqr("Left", "Transpose", m, 1, n,
			&af[0], m, &tau[0], &c[0], m);
	assert(info == 0);
	trc.dprint("c[nc-1] = ",c[nc-1]);
	if (c[nc-1] < 0.0) {
		rval = -rval;
		trc.dprint("Q^t*tan < 0 so changing sign: ",rval);
	}

	// if there was pivoting in the factorization change
	// the sign once per pivot
	for (size_t i=0; i<nr; i++)
		if (pivots[i] != (int)(i+1))
			rval = -rval;
	
	trc.dprint("returning ",rval," * 10^",expon);
	return rval;
}

double
eigdeterm(int n, vector<double> const& A, double& expon) {
// compute the determinant of a real (n,n) matrix in
// the form d*10^expon
	Trace trc(3,"eigdeterm");

	vector<double> as = A;
	vector<double> wr(n);
	vector<double> wi(n);
	double v;
	int info = lapack::dgeev("n","n",n,&as[0],n,
				&wr[0],&wi[0],&v,n, &v,n);
	if (info != 0)
		throw runtime_error(vastr("eigdeterm: dgeev returned ",info));

	trc.dprintvv(wr,wi,"eigenvalues");

	// the determinant is the product of the eigenvalues
	// adjust det and expon to keep 1 < abs(det) < 10
	complex<double> det(1.0, 0.0);
	expon = 0.0;
	for (int i=0; i<n; i++) {
		complex<double> ei(wr[i],wi[i]);
		trc.dprint("eig ",i," = ",ei);
		det *= ei;
		if (abs(det) == 0.0)
			break;
		while(abs(det) < 1.0) {
			trc.dprint("det ",det," < 1: mult by 10, expon ",expon);
			det *= 10.0;
			expon -= 1.0;
		}
		while(abs(det) > 10.0) {
			trc.dprint("det ",det," > 10: divide by 10, expon ",expon);
			det /= 10.0;
			expon += 1.0;
		}
	}
	trc.dprint("returning ",det,"*10^",expon);
	double rval = det.real();
	return rval;
}

double
augeigdet(size_t nr, size_t nc, vector<double> const& J,
		vector<double> const& tan, double& expon) {
// compute the determinant of an (m,m+1) matrix augmented
// by placing "tan" in the last (m+1st) row and computing
// the eigenvalues of the square matrix. The determinant
// is the product of the eigenvalues.
	Trace trc(2,"augeigdet");
	size_t i, j;
	expon = 0.0;
	if (nc != nr+1) {
		trc.dprint("quick return: nc (",nc,") != nr (",nr,")+1");
		return 0.0;
	}
	// copy the input (nr,nc) matrix into ja and the
	// tan into the last row
	vector<double> ja(nc*nc, 0.0);
	assert(tan.size() == nc);
	for (j=0; j<nc; j++) {
		for (i=0; i<nr; i++) {
			ja[IJ(i,j,nc)] = J[IJ(i,j,nr)];
		}
		ja[IJ(nr,j,nc)] = tan[j];
	}

	trc.dprintm(nc,nc,nc,ja,"augmented Jacobian");

	// compute the determinant at the product of the eigenvalues
	double rval = eigdeterm(nc, ja, expon);
	return rval;
}

#ifdef MAIN

#undef MAIN

#include "exim.h"

// defaults
bool Square = false;		// only test square matrices?
bool Underdetermined = false;
int Nsol = 10;				// number of solutions
bool Scaling = false;	// row, column scaling
int Scalingrange = 6;	// scale factor limit power of 10
double Rcond = 1.0e-5;  // will be zero if Rankdef > 0
bool Usennull = true;  // use Nnull to determine nc = nr + Nnull
int Nnull = 1;				// nc = nr + Nnull
int Rankdef = 0;        // number of rank deficiencies
double Maxelem{1000.0}; // max element size
// int Maxsize = 300;		// largest matrix dimension
int Maxsize = 30;		// largest matrix dimension
int Minsize = 8;			// smallest matrix dimension
#define SEED 47

static void
initRand() {
	static int first = 1;
	if (first)
		srand(SEED);
	first = 0;
}

static double
randReal (double range) {
	Trace trc(3,"randReal");
	double rval;

	trc.dprint("randReal range<",range,"> RAND_MAX<",RAND_MAX,">");

	initRand();

	rval = (double)rand()/(double)RAND_MAX;
	trc.dprint("rand/RAND_MAX <",rval,">");
	rval = range*(2.0*rval - 1.0);
	trc.dprint("returning ",rval);
	return rval;
}

static double
xrand() {
	int k = rand();
	double rval;

	rval = (double)k/(double)RAND_MAX;
	rval = 2.0*rval - 1.0;
	return rval;
}

static int
irand(int range) {
// Random integer in the range 1 to range
	Trace trc(1,"irand");
	int rval;
	double x;

	trc.dprint("range 1-",range);

	initRand();

	x = (double)rand()/(double)RAND_MAX;
	trc.dprint("rand/RAND_MAX: ",x);
	x = ((double)range)*x;
	rval = (int)x;
	if (rval <= 0)
		rval = 1;
	trc.dprint("returning ",rval);
	return rval;
}


static void
norm (int n, double *x) {
	double t;

	t = blas_snrm2 (n, x, 1);
	if (t != 0.0)
		t = 1.0/t;
	blas_scal (n, t, x, 1);
}

static void
multMv (int nr, int nc, double *a, double *x, double *b, double *c) {
/* multiply Ax - b = c accumulating in double */
	int i, j;
	double t;

	for (i=0; i<nr; i++) {
		t = 0.0;
		for (j=0; j<nc; j++)
			t += a[IJ(i,j,nr)]*x[j];
		if (b)
			t -= (b[i]);
		c[i] = t;
	}
}

vector<double>
newA (int nr, int nc, double rcond, double maxelem) {
/*
 * Create a random double matrix with a specified size and
 * reciprocal condition number.
 * Method:
 *   1) create a random (nr,nc) matrix, do a singular-value
 *      decomposition:  A = U \Sigma V^t
 *   2) change the singular values to be evenly spaced
 *      between rcond*maxelem and maxelem
 *   3) form the return matrix as A = U \Sigma V^t where
 *      \Sigma is the new set of singular values.
 */
	Trace trc(2,"newA");
	vector<double> rval(nr*nc, 0.0);
	double smin, smax, sk, dels;
	int lwork, info;
	int i, j, k, mn;

	trc.dprint("nr, nc: ",nr,", ",nc,", rcond ",rcond);

	mn = std::min(nr, nc);

	// form a random (nr,nc) matrix
	for (i=0; i<nr*nc; i++) {
		rval[i] = xrand();
	}
	trc.dprintm(nr,nc,nr,rval,"random A");

	// Do an SVD
	vector<double> s(mn);
	vector<double> u(nr*nr);
	vector<double> vt(nc*nc);
	lwork = std::max(3*mn + std::max(nr,nc), 5*mn-4);
	lwork *=4;
	vector<double> work(lwork);

	// info = ftn_dgesvd ("A", "A", nr, nc, &rval[0], nr,
	// 	&s[0], &u[0], nr, &vt[0], nc);
	dgesvd_("A", "A", &nr, &nc, &rval[0], &nr,
	 	&s[0], &u[0], &nr, &vt[0], &nc, &work[0], &lwork, &info);

	if (info != 0) {
		string exc = vastr("error in dgesvd: info = ",info);
		trc.dprint("throwing exception: ",exc);
		throw runtime_error(exc);
	}

	trc.dprintm(nr,nr,nr,u,"U");
	trc.dprintm(nc,nc,nc,vt,"Vt");

	vector<double> utu(nr*nr);
	blas_sgemm("t","n",nr,nr,nr,1.0,&u[0],nr,&u[0],nr,0.0,&utu[0],nr);
	trc.dprintm(nr,nr,nr,utu,"U(t)*U");

	for (i=0; i<nr; i++)
		for (j=0; j<nc; j++)
			rval[IJ(i,j,nr)] = 0.0;

	// Create the return matrix by forming the product
	//     U \Sigma V^t
	// where the singular values are evenly spaced between rcond*maxelem
	// and maxelem

	smin = rcond*maxelem;
	smax = maxelem;
	if (mn > 1) {
		dels = (smax - smin)/((double)(mn-1));
	} else {
		dels = 0.0;
	}
	for (k=0; k<mn; k++) {
		sk = smin + k*dels;
		if (k < Rankdef)
			sk = 0.0;
		trc.dprint("setting singular value ",k," 0b to ",sk);
		for (i=0; i<nr; i++) {
			for (j=0; j<nc; j++) {
				rval[IJ(i,j,nr)] = rval[IJ(i,j,nr)] + u[IJ(i,k,nr)]*sk*vt[IJ(k,j,nc)];
			}
		}
	}

	// do an SVD on the return matrix to ensure the singular-value
	// distribution is correct
	vector<double> as = rval;
	trc.dprintm(nr,nc,nr,as,"new A");
	// info = ftn_dgesvd ("N", "N", nr, nc, &as[0], nr,
	// 	&s[0], &u[0], nr, &vt[0], nc);
	dgesvd_("N", "N", &nr, &nc, &as[0], &nr,
	 	&s[0], &u[0], &nr, &vt[0], &nc, &work[0], &lwork, &info);
	if (info != 0) {
		string exc = vastr("error in dgesvd on new A: info = ",info);
		trc.dprint("throwing exception: ",exc);
		throw runtime_error(exc);
	}
	trc.dprint("singular values of return A",s);
	double newrcond = s[mn-1]/s[0];
	if (Rankdef == 0 && !is_equal(newrcond,rcond,2)) {
		cout << "rcond(new A) (" << newrcond << " != input rcond (" << rcond << ")\n";
	}

	// Row, column scaling
	if (Scaling) {
		vector<double> rscale(nr);
		vector<double> cscale(nc);
		for (i=0; i<nr; i++)
			rscale[i] = pow(10.0, (double)irand(Scalingrange));
		for (i=0; i<nc; i++)
			cscale[i] = pow(10.0, (double)irand(Scalingrange));
		trc.dprint("row scale factors: ",rscale);
		trc.dprint("col scale factors: ",cscale);
		for (i=0; i<nr; i++)
			blas_scal (nc, rscale[i], &rval[IJ(i,0,nr)], nr);
		for (j=0; j<nc; j++)
			blas_scal (nr, cscale[j], &rval[IJ(0,j,nr)], 1);
	}
	return rval;
}

vector<double>
newb (int nr, int nc, vector<double>& a) {
	vector<double> rval(nr, 0.0);
	int i;

	if (nr == nc) {
		for (i=0; i<nr; i++)
			rval[i] = (double)randReal(1.0);
	} else {
		vector<double> x(nc, 0.0);
		for (i=0; i<nc; i++)
			x[i] = (double)randReal(1.0);
		blas_sgemv("n",nr,nc,1.0,&a[0],nr,&x[0],1,0.0,&rval[0],1);
	}
	return rval;
}

void
testFile (char const* file) {
	Trace trc(2,"testFile");
	int pi;
	Matrix* ap = MM::importer(file);
	double* a = ap->elem();
	size_t nr = ap->rsize();
	size_t nc = ap->csize();

	// only use the real part of complex matrices
	if (ap->is_complex()) {
		double* ra = new double[nr*nc];
		blas_copy(nr*nc, a, 2, ra, 1);
		delete ap;
		a = ra;
	}
	// the QR constructor factorizes a
	QR af(nr, nc, &a[0]);

	// If we have a rectangular matrix compute the nullspace
	if (nc == nr+1) {
		vector<double> dir(nc, 0.0);
		dir[0] = -1.0;
		vector<double> tan = af.nullProj (dir);
		double expon;
		double det = af.augJacDet(af.af, tan, expon);
		trc.dprint("determiant: ",det,"*10^",expon);
	}

	trc.dprint("rcond: ",af.rcond);
	vector<double> b(nr, 0.0);
	for (size_t i=0; i<nr; i++)
		b[i] = xrand();

	vector<double> x(nc, 0.0);
	// af.solve(1, &b[0], &x[0], pi);
	af.solve(b, x);
	pi = af.performanceIndex(&a[0], &x[0], &b[0]);
	trc.dprint("before refine: pi = ",pi);
	// pi = af.refine(&x[0],&b[0]);
	// T("after refine: pi = %d\n", pi);
}

int
main(int argc, char** argv) {
	Trace trc(1,argv[0]);
	int nr, nc;
	int k;
	int pi, nullpi;
	double resNorm, nullResNorm;
	double det{0.0};
	double expon{0.0};

	if (argc > 1) {
		testFile(argv[1]);
		exit(0);
	}


	printf ("      Exact rcond: %7.0e\n", Rcond);
	if (Square)
		printf ("      Square matrices only\n");
	if (Scaling)
		printf ("      Row and Column scaling: [%7.0e : %7.0e]\n",
			1.0, pow(10.0, Scalingrange));
	printf ("      ----------------------------\n");

	printf ("                    residual        null vector"
		"              rank\n");
	printf ("  size              norm (pi)     residual norm    "
		"rcond     def     determinant\n");

	for (k=0; k<Nsol; k++) {
		nr = irand(Maxsize);
		nr = std::max(nr, Minsize);
		trc.dprint("solution ",k+1,": nr = ",nr);
		if (Square) {
			nc = nr;
		} else if (Usennull) {
			nc = nr + Nnull;
		} else {
			nc = irand(Maxsize);
			nc = std::max(nc, Minsize);
			if (Underdetermined && nr > nc) {
				int n = nr;
				nr = nc;
				nc = n;
			}
		}
		vector<double> x(nc);
		vector<double> res(nr);

		// Create a random A matrix
		vector<double> a = newA(nr,nc,Rcond,Maxelem);

		//... and a random rhs "b"
		vector<double> b = newb(nr,nc,a);

		trc.dprintm(nr,nc,nr,a,"a = ");
		trc.dprint("b = ",b);
		string shape("ls");
		if (nr < nc)
			shape = "ud";

		// Factor the matrix...
		QR af(nr, nc, &a[0]);

		trc.dprint("rcond ",af.rcond,", pivots",af.pivots);
		// ... and solve Ax = b
		af.solve(b, x);
		pi = af.performanceIndex(&a[0], &x[0], &b[0]);
		trc.dprint("solution",x);

		// pi = af->refine(x, b);
		// T("pi = %d\n",pi);

		multMv (nr, nc, &a[0], &x[0], &b[0], &res[0]);

		resNorm = blas_snrm2 (nr, &res[0], 1);
		trc.dprint("residual",res);
		// If we have a rectangular matrix compute the nullspace
		if (nc > nr) {
			vector<double> dir(nc, 0.0);
			dir[0] = -1.0;
			vector<double> xv = af.nullProj (dir);
			norm(nc, &xv[0]);						/* normalize the null vector */
			multMv (nr, nc, &a[0], &xv[0], nullptr, &res[0]);
			nullResNorm = blas_snrm2 (nr, &res[0], 1);
			trc.dprint("null vector",xv);
			trc.dprint("residual of null vector",res);
			dir[0] = 0.0;
			nullpi = af.performanceIndex(&a[0], &xv[0], &dir[0]);
		} else {
			nullResNorm = 0.0;
			nullpi = 0;
		}
		// If nc == nr+1 compute the determinant
		if (nc == nr+1) {
			vector<double> dir(nc, 0.0);
			dir[0] = -1.0;
			vector<double> tan = af.nullProj (dir);
			double tn = blas_snrm2(nc, &tan[0], 1);
			blas_scal(nc, 1.0/tn, &tan[0], 1);
			det = af.augJacDet(af.af, tan, expon);
			double eigexp;
			double eigdet = augeigdet(nr, nc, a, tan, eigexp);
			if (!is_equal(eigdet,det,3) || !is_equal(eigexp,expon,3)) {
				cerr << "qr determinant (" << det << "*10^" << expon
					<< " != eig det (" << eigdet << "*10^" << eigexp << endl;
			}
		} else {
			det = 0.0;
			expon = 0.0;
		}

		printf (" (%2d,%2d) %s    %7.0e (%2d)      %7.0e (%d)     "
			"%7.0e    %3d     %7.0e 10^%d\n",
			nr, nc, shape.c_str(), resNorm, pi,
			nullResNorm, nullpi,
			af.rcond, af.rankdef, det, (int)expon);
	}
}

#endif
