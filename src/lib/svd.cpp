//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

// minimum 2-norm solutions of underdetermined linear systems
// as described in Algorithm 5.6.2 in
//    Golub, G.H.; Van Loan, C.F. Matrix Computations, 4th ed; The
//    Johns Hopkins University Press: Baltimore, MD, USA, 2013

#include "config.h"
#include "lapack.h"
#include "matrix.h"
#include "qr.h"
#include "svd.h"
#include "trace.h"
#include "transpose.h"

using namespace std;

SVD::
SVD (int nrow, int ncol, double const* A) {
/*------------------------------------------------------------------
 * SVD constructor
 * Factor the (nrow,ncol) real matrix A using an SVD factorization
 *   A = \Matrix{U} \Matrix{\Sigma} \Matrix{V}^T
 *          (m,m)            (m,n)      (n,n)
 *
 * Returns an SVD structure
 *------------------------------------------------------------------*/
	T_(Trace trc(2,"SVD constructor");)
	nr = nrow;
	nc = ncol;
	double eps = 100.0*std::numeric_limits<double>::epsilon();

	T_(trc.dprint("m ",nr,", n ",nc);)

	int mn = std::min(nr, nc);

	// Copy the matrix to another spot for factorization
	a = A;
	vector<double> af(nr*nc);
	blas::copy (nr*nc, &A[0], 1, &af[0], 1);

	s.resize(mn);
	u.resize(nr*nr);
	vt.resize(nc*nc);

	// Factorize using LAPACK routine DGESVD into U*S*V(t)
	// with s sorted in descending order
	int info = lapack::dgesvd("A", "A", nr, nc, &af[0], nr,
	 	&s[0], &u[0], nr, &vt[0], nc);

	if (info != 0) {
		string exc = vastr("SVD failed: info ",info);
		// T_(trc.dprint("throwing exception: ",exc);)
		T_(trc.dprintm(nr,nc,nr,&af[0],"throwing exception: ",exc);)
		throw runtime_error(exc);
	}

	// The reciprocal condition number is simply the ratio
	// of the smallest to largest singular values (s[mn-1]/s[0])
	if (s[0] > eps)
		rcond = s[mn-1]/s[0];
	else
		rcond = eps;

	rankdef = this->rank_def();

	T_(trc.dprint("rcond = ",rcond,", rankdef = ",rankdef);)
}

ostream&
operator<<(ostream& s, const SVD& t) {
	size_t m = t.nr;
	size_t n = t.nc;
	s << "left singular vectors (U):\n" << strArray(m,m,m,&t.u[0]) << endl;
	s << "singular values (s): {\n" << t.s << "}\n";
	s << "right singular vectors (Vt):\n" << strArray(n,n,n,&t.vt[0]) << endl;
	return s;
}

void
SVD::
solve (vector<double> const& b, vector<double>& x) {
/*------------------------------------------------------------------
 * Solve a system of linear equations with one right-hand-side.
 * If nr < nc the solution is a minimum-norm solution to the
 * underdetermined system. Ref: golub2013matrix pg 300
 *    Ax = b
 *    A: (nr,nc) = U s V^t
 *    x: nc-vector
 *    b: nr-vector
 * throws runtime_error if an error occured.
 *------------------------------------------------------------------*/
	T_(Trace trc(2,"SVD::solve");)
	int j;

	T_(trc.dprintm(nr,1,nr,b,"b");)

	assert((int)b.size() >= nr && nr > 0);
	assert((int)x.size() >= nc && nc > 0);

	for (j=0; j<nc; j++)
		x[j] = 0.0;

	// Ax = (U s V^t)x = b
	//  A: m,n
	//  U: m,m
	//  s: m,n
	//  V^t: n,n
	// x = V s^{-1} U^tb = \sum_{i=1}^{r} \frac{u_i^t b}{\sigma_i} v_i
	// where r = rank(A). To be full rank all nr s[i] must be > tol
	int mn = std::min(nr, nc);
	// double eps = 8.0*std::numeric_limits<double>::epsilon();
	double eps = std::numeric_limits<double>::epsilon();
	double tol = eps*s[0];
	int rank = 0;
	for (j=0; j<mn; j++) {
		if (s[j] < tol)
			break;
		double ub;
		blas::dot(nr, &u[IJ(0,j,nr)], 1, &b[0], 1, ub);
		T_(trc.dprint("u_",j,"^t b = ",ub);)
		ub /= s[j];
		blas::axpy (nc, ub, &vt[j], nc, &x[0], 1);
		rank++;
	}
	T_(trc.dprint("rankdef by rank_def ",rankdef,", and here: ",nr-rank);)
	rankdef = nr - rank;
	// the last nc-nr+rankdef rows of vt are null vectors
	if (rankdef == 1) {
		int nnull = nc - nr + rankdef;
		vector<double> nv(nc);
		for (int i=1; i<=nnull; i++) {
			blas::copy(nc, &vt[nc-i], nc, &nv[0], 1);
			cerr << "null vector " << i << ", s=" << s[nc-1] <<
				": " << flaps::summarize(nv, 80) << endl;
			double alpha{1.0};
			double beta{0.0};
			vector<double> y(nr);
			blas::gemv("n",nr,nc,alpha,&a[0],nr,&nv[0],1,beta,&y[0],1);
			cerr << "norm J*nv = " << blas::snrm2(nr, &y[0], 1) << endl;
		}
	}

	T_(trc.dprint(rankdef," rank deficiencies");)
}

vector<Interval>
SVD::
solve (vector<Interval> const& b) {
/*------------------------------------------------------------------
 * Solve a system of linear equations with one right-hand-side.
 * Special version of SVD::solve to handle interval rhs.
 *    Ax = b
 *    A: (nr,nc) = U s V^t
 *    x: nc-vector
 *    b: nr-vector
 * throws runtime_error if an error occured.
 *------------------------------------------------------------------*/
	T_(Trace trc(2,"SVD::solve(Interval)");)

	T_(trc.dprint("b",b);)
	
	//!! if (nr != nc)
		//!! throw runtime_error(vastr("only square matrices allowed, input is (",nr,',',nc,")"));
	if ((int)b.size() < nr)
		throw runtime_error(vastr("SVD: rhs is ",b.size(),", matrix order is ",nr));

	// invert the matrix to minimize the number of interval
	// ops required: A^{-1} = V \Sigma^{-1} U^T
	// U: (nr,nr)
	// V: (nc,nc)
	// \Sigma: (nr,nc)
	vector<double> tmp(nr*nr,0.0);
	flaps::transpose(nr,nr,&u[0], &tmp[0]);
	for (int i=0; i<nr; i++) {
		double sinv = 1.0/s[i];
		for (int j=0; j<nr; j++)
			tmp[i+nr*j] *= sinv;
	}
	// multiply the first nr columns of V by tmp
	// a inverse is (nc,nr)
	vector<double> ainv(nc*nr,0.0);
	for (int i=0; i<nc; i++) {
		for (int j=0; j<nr; j++) {
			double t{0.0};
			for (int k=0; k<nr; k++) {
				t += vt[k+i*nc]*tmp[k+j*nr];
			}
			ainv[i+j*nc] = t;
		}
	}
	// verify ainv
	vector<double> ident(nr*nr, 0.0);
	for (int i=0; i<nc; i++) {
		for (int j=0; j<nr; j++) {
			double t{0.0};
			for (int k=0; k<nc; k++)
				t += a[i+k*nr]*ainv[k+j*nc];
			ident[i+j*nr] = t;
		}
	}
	T_(trc.dprint("a*ainv = ",ident);)

	// x = A^{-1} b
	vector<Interval> rval(nc);
	for (int i=0; i<nc; i++) {
		Interval t(0.0);
		for (int j=0; j<nr; j++) {
			t += ainv[i+j*nc]*b[j];
		}
		rval[i] = t;
	}
	T_(trc.dprint("solution x:",rval);)
	return rval;
}

std::vector<double>
SVD::
nullProj(const vector<double>& proj) {
	T_(Trace trc(2,"SVD::nullProj");)
// given an svd factorization A = U \Sigma V^t project a
// vector "proj" onto the nullspace:
//   p = V V^t proj
// The last nc-nr+rankdef columns of V (rows of vt) are null vectors
	vector<double> rval(nc, 0.0);
	T_(int nnull = nc - nr + rankdef;)
	T_(trc.dprint(nnull, " null vectors, start at V[",nc-nnull,"] 0b");)
	double vtp;
	int start = (int)(nr - rankdef);
	int end = (int)nc;
	for (int i=start; i<end; i++) {
		blas::dot(nc, &vt[i], nc, &proj[0], 1, vtp);
		T_(trc.dprint("row ",i," vtp ",vtp);)
		if (vtp != 0.0)
			blas::axpy(nc, vtp, &vt[i], nc, &rval[0], 1);
		else {
			blas::copy(nc, &vt[i], nc, &rval[0], 1);
			string msg = vastr("attempt to project ",proj," onto svd null vector:",rval);
			flaps::warning(msg);
		}
	}
	// check that the projection is in the same direction as the
	// projected vector
	double t;
	blas::dot (nc, &rval[0], 1, &proj[0], 1, t);
	if (t < 0.0) {
		t = -1.0;
		blas::scal (nc, t, &rval[0], 1);
	}
	T_(trc.dprint("returning ",flaps::summarize(rval));)
	return rval;
}

int
SVD::
perfIndex (double* x, double* dir) { return 0; }


int
SVD::
rank_def() {
/*------------------------------------------------------------------
 * Estimate the rank deficiency of a matrix given its SVD: U(t)*A*V = S
 * Two ways of doing this:
 * 1) by estimating the change in s with roundoff-sized
 *    perturbations in A; if these perturbations could result in
 *    a zero s we declare that to be a rank deficiency. To do this
 *    we take each element of A to be a parameter, then use an
 *    approximation to the partial of s wrt a parameter, and sum
 *    all partials times the matrix element perturbation.
 *    For small s the derivative of s wrt a parameter is approx
 *      s' = u'(t)*A*v + u(t)*A'v + u(t)*A*v'
 *         = u'(t)*u*s + u(t)*A'*v + s*v(t)*A*v'
 *         ~ u(t)*A'*v
 * 2) the more traditional way is to simply look at ratios
 *    s(i)/s[0] where s[0] is the largest singular value. If
 *    a ratio is less than machine epsilon that is a rank
 *    deficiency.
 *------------------------------------------------------------------*/
	T_(Trace trc(2,"SVD::rank_def");)
	int ns = std::min(nr, nc);
	int rank;
	//!! double eps = sqrt(std::numeric_limits<double>::epsilon());
	double eps = rcondeps();

	T_(trc.dprint("m ",nr,", n ",nc);)
//!! #ifdef NEVER
	// method (1): estimate the interval each s lies in,
	// decrement rank if interval contains zero
	int rank1 = ns;
	for (int k=ns-1; k>=0; k--) {
		double t = 0.0;
		for (int i=0; i<nr; i++) {
			for (int j=0; j<nc; j++) {
				t += std::abs(u[IJ(i,k,nr)]*a[IJ(i,j,nr)]*vt[IJ(k,j,nc)]);
			}
		}
		t *= eps;
		if (s[k] == 0.0 || t > s[k]) {
			rank1--;
		}
	}
//!! #else
	// method (2): the number of s which are < eps*largest s
	eps = s[0]*rcondeps();
	T_(trc.dprint("eps = ",s[0],"*",rcondeps()," = ",eps);)
	rank = 0;
	for (int i=0; i<ns; i++) {
		if (s[i] < eps)
			break;
		rank++;
	}
//!! #endif // NEVER
	T_(trc.dprint("rank by method 1 ",rank1,", and by method 2 ",rank);)
	T_(trc.dprint("returning rank def ",ns-rank);)
	return ns-rank;
}


#ifdef MAIN

#undef MAIN

#include "exim.h"

bool Square = true;			/* only test square matrices? */
int Nsol = 1;					/* number of solutions */
bool Scaling = false;			/* row, column scaling */
int Scalingrange = 6;		/* scale factor limit power of 10 */
double Rcond = 1.0e-5;
// #define USENNULL
#define NNULL 1								/* nc = nr + NNULL */
#define UNDERDETERMINED 1  // nc >= nr
int Rankdef = 0;				// number of rank deficiencies
// int Maxsize = 600;		/* largest matrix dimension */
int Maxsize = 10;		/* largest matrix dimension */
int Minsize = 10;			/* smallest matrix dimension */
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
// returns a random number in the range [-range:range]
	initRand();
	double rval = (double)rand()/(double)RAND_MAX;
	rval = range*(2.0*rval - 1.0);
	return rval;
}

static double xrand() {
// returns a random number in the range [-1:1]
// XXX duplicates randReal()
	int k = rand();
	double rval = (double)k/(double)RAND_MAX;
	rval = 2.0*rval - 1.0;
	return rval;
}

static int
irand(int range) {
// Random integer in the range 1 to range
	T_(Trace trc(1,"irand");)
	int rval;
	double x;

	T_(trc.dprint("range: ",range);)

	initRand();

	x = (double)rand()/(double)RAND_MAX;
	T_(trc.dprint("rand/RAND_MAX ",x);)
	x = ((double)range)*x;
	rval = (int)x;
	if (rval <= 0)
		rval = 1;
	T_(trc.dprint("returning ",rval);)
	return rval;
}


static void
norm (int n, double* x) {
// scale input "x" so it's 2-norm is 1.0
	double t = blas::snrm2 (n, x, 1);
	if (t != 0.0)
		t = 1.0/t;
	blas::scal (n, t, x, 1);
}

static bool
isResOk (int nr, int nc, int exp, vector<double>& a, vector<double>& x, vector<double>& b) {
/*------------------------------------------------------------------
 * Answers the question "is the residual (Ax-b) in the computed
 * solution small enough?"
 *------------------------------------------------------------------*/
	T_(Trace trc(1,"isResOk");)
	long double rmin, rmax;
	long double ax;
	long double tlo, thi;
	int i, j;
	double eps;
	double machineEps(std::numeric_limits<double>::epsilon());

	eps = pow(2.0, (double)exp)*machineEps;
	tlo = 1.0 - eps;
	thi = 1.0 + eps;

	T_(trc.dprint("eps ",eps,", tlo ",tlo,", thi ",thi);)

	for (i=0; i<nr; i++) {
		rmin = rmax = 0.0;
		for (j=0; j<nc; j++) {
			ax = a[IJ(i,j,nr)]*x[j];
			if (ax < 0.0) {
				rmin += thi*ax;
				rmax += tlo*ax;
			} else {
				rmin += tlo*ax;
				rmax += thi*ax;
			}
		}
		if (rmin > b[i] || rmax < b[i]) {
			T_(trc.dprint("returning false: ",rmin," <= ",b[i]," <= ",rmax," failed");)
			return false;
		} else {
			T_(trc.dprint(rmin," <= ",b[i]," <= ",rmax);)
		}
	}
	return true;
}

static void
multMv (int nr, int nc, vector<double>& a, vector<double>& x,
		vector<double>& b, vector<double>& c) {
// multiply Ax - b = c accumulating in long double
	int i, j;
	long double t;

	for (i=0; i<nr; i++) {
		t = 0.0;
		for (j=0; j<nc; j++)
			t += (long double)(a[IJ(i,j,nr)])*(long double)(x[j]);
		t -= (long double)(b[i]);
		c[i] = (double)t;
	}
}

vector<double>
newA (int nr, int nc, double rcond) {
// Create a random Real matrix with a specified size and
// reciprocal condition number.
// Method:
//   1) create a random (nr,nc) matrix, do a singular-value
//      decomposition:  A = U \Sigma V^t
//   2) change the singular values to be evenly spaced
//      between rcond and 1.0
//   3) form the return matrix as A = U \Sigma V^t where
//      \Sigma is the new set of singular values.
	T_(Trace trc(1,"newA");)
	vector<double> rval(nr*nc, 0.0);
	int i, j, k;

	T_(trc.dprint("nr ",nr,", nc ",nc,", rcond ",rcond);)

	int mn = std::min(nr, nc);

	// form a random (nr,nc) matrix
	for (i=0; i<nr; i++) {
		for (j=0; j<nc; j++) {
			rval[IJ(i,j,nr)] = xrand();
		}
	}
	T_(trc.dprintm(nr,nc,nr,&rval[0],"random A");)

	// SVD the random matrix
	vector<double> s(mn, 0.0);
	vector<double> u(nr*nr, 0.0);
	vector<double> vt(nc*nc, 0.0);
	int info = lapack::dgesvd("A", "A", nr, nc, &rval[0], nr,
	 	&s[0], &u[0], nr, &vt[0], nc);

	if (info != 0) {
		throw runtime_error(vastr("svd failed: info ",info));
	}

	T_(trc.dprintm(nr,nr,nr,&u[0],"U");)
	T_(trc.dprintm(nc,nc,nc,&vt[0],"Vt");)

	vector<double> utu(nr*nr, 0.0);
	blas::gemm("t","n",nr,nr,nr,1.0,&u[0],nr,&u[0],nr,0.0,&utu[0],nr);

	for (i=0; i<nr; i++)
		for (j=0; j<nc; j++)
			rval[IJ(i,j,nr)] = 0.0;

	// Create the return matrix by forming the product
	//   U \Sigma V^t
	// where the singular values are evenly spaced between rcond
	// and 1.0
	double smin = rcond;
	double smax = std::max(1.0, 10.0*rcond);
	double dels{0.0};
	if (mn > 1) {
		dels = (smax - smin)/((double)(mn-1));
	}
	for (k=0; k<mn; k++) {
		double sk = smin + k*dels;
		if (k < Rankdef)
			sk = 0.0;
		T_(trc.dprint("setting singular value ",k," to ",sk);)
		for (i=0; i<nr; i++) {
			for (j=0; j<nc; j++) {
				rval[IJ(i,j,nr)] = rval[IJ(i,j,nr)] + u[IJ(i,k,nr)]*sk*vt[IJ(k,j,nc)];
			}
		}
	}

	vector<double> as(nr*nc);
	blas::copy (nr*nc, &rval[0], 1, &as[0], 1);
	T_(trc.dprintm(nr,nc,nr,&as[0],"new A");)
	info = lapack::dgesvd("N", "N", nr, nc, &as[0], nr,
				&s[0], &u[0], nr, &vt[0], nc);
	if (info != 0) cerr << "info= " << info << endl;
	T_(trc.dprint("singular values of return A:",s);)
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
		blas::gemv("n",nr,nc,1.0,&a[0],nr,&x[0],1,0.0,&rval[0],1);
	}
	return rval;
}

void
testFile (char const* file) {
// read a matrix from a Matrix Market file and svd it
	int nr, nc;
	bool iscmplx;
	vector<double> va = MM::importer(file, nr, nc, iscmplx);
	double* a = &va[0];

	// if the matrix is complex copy the real part only
	if (iscmplx) {
		double* ta = new double[nr*nc];
		blas::copy(nr*nc, a, 2, ta, 1);
		a = ta;
	}
	SVD af(nr, nc, a);
	vector<double> b(nr, 0.0);
	for (int i=0; i<nr; i++)
		b[i] = xrand();

	vector<double> x(nc, 0.0);
	af.solve(b, x);
}

int
main(int argc, char** argv) {
	T_(Trace trc(1,argv[0]);)
	int nr, nc, nrnc;
	int i, k;
	double resNorm, nullResNorm;
	double factorTime{0.0};
	SVD* af;
	ostringstream os;

	if (argc > 1) {
		testFile(argv[1]);
		exit(0);
	}

	if (Nsol > 5)
		debug(std::max(debug(), 1));		/* MAIN */


	printf ("      Exact rcond: %7.0e\n", Rcond);
	if (Square)
		printf ("      Square matrices only\n");
	if (Scaling)
		printf ("      Row and Column scaling: [%7.0e : %7.0e]\n",
			1.0, pow(10.0, Scalingrange));
	printf ("      ----------------------------\n");

	printf ("                    residual    null vector"
		"            rank        solution\n");
	printf ("  size              norm       residual norm   "
		"rcond    def          time\n");

	for (k=0; k<Nsol; k++) {
		nr = irand(Maxsize);
		nr = std::max(nr, Minsize);
		T_(trc.dprint("solution ",k+1,": nr = ",nr);)
		if (Square) {
			nc = nr;
		} else {
#ifdef USENNULL
			nc = nr + NNULL;
#else
			nc = irand(Maxsize);
			nc = std::max(nc, Minsize);
#ifdef UNDERDETERMINED
			if (nr > nc) {
				int tmp = nr;
				nr = nc;
				nc = tmp;
			}
#endif // UNDERDETERMINED
#endif // USENNULL
		}
		nrnc = nr*nc;
		vector<double> as(nrnc, 0.0);
		vector<double> x(nc, 0.0);
		vector<double> bs(nr, 0.0);
		vector<double> res(nr, 0.0);

		// Create a random A matrix...
		vector<double> a = newA(nr,nc,Rcond);

		// ... and a random rhs "b"
		vector<double> b = newb(nr,nc,a);

		// save copies of A and b for computing the residual
		blas::copy (nrnc, &a[0], 1, &as[0], 1);
		blas::copy (nr, &b[0], 1, &bs[0], 1);

		string shape("ls");
		if (nr < nc)
			shape = "ud";
		else if (nr == nc)
			shape = "sq";

		// Factor the matrix...
		try {
			af = new SVD(nr, nc, &a[0]);
		} catch(runtime_error& s) {
			cerr << "factorization failed: " << s << endl;
			continue;
		}
		T_(trc.dprint(" rcond = ",af->rcond);)

		// ... and solve Ax = b
		af->solve(b, x);
		T_(trc.dprint("solution: \n",x);)

		multMv (nr, nc, as, x, bs, res);

		resNorm = blas::snrm2 (nr, &res[0], 1);
		T_(trc.dprint("residual: \n",res);)
		bool resok = isResOk (nr, nc, 4, as, x, bs);
		// If we have a rectangular matrix compute the nullspace
		if (nc > nr) {
			vector<double> dir(nc, 0.0);
			for (auto& di : dir)
				di = xrand();
			vector<double> xv = af->nullProj (dir);
			norm(nc, &xv[0]);						/* normalize the null vector */
			vector<double> zero(nr, 0.0);
			multMv (nr, nc, as, xv, zero, res);
			nullResNorm = blas::snrm2 (nr, &res[0], 1);
			T_(trc.dprint("null vector:\n",xv);)
			T_(trc.dprint("residual of null vector:\n",res);)
			for (i=0; i<nc; i++)
				dir[i] = 0.0;
		} else {
			nullResNorm = 0.0;
		}

		printf (" (%3d,%3d) %s    %7.0e (%s)      %7.0e     "
			"%7.0e   %3d           %g\n",
			nr, nc, shape.c_str(), resNorm, resok?"ok":"!!",
			nullResNorm, af->rcond, af->rankdef, factorTime);

	}
}

#endif
