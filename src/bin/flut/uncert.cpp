//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
#include <string>
#include <vector>

#include "blas.h"
#include "interval.h"
#include "lapack.h"
#include "uncert.h"

using namespace std;
//!! int dgefa_(double *a, int *lda, int *n, int * ipvt, int *info);
//!! int dgesl_(double *a, int *lda, int *n, int * ipvt, Interval* b);
#include "dge.cpp"

double
vector_width(vector<Interval> const& a) {
// compute the sqrt of the sum of widths squared
	double rval{0.0};
	for (auto& x : a) {
		double w = x.width();
		rval += w*w;
	}
	return sqrt(rval);
}

Uncert::
Uncert (const Flutcurve& curve, const vector<string>& uncnames,
	const vector<double>& factors, const vector<string>& constpar) {
// Uncert constructor
	Trace trc(2,"Uncert::constructor");
	// create a copy of curve.params to use for constructing hi & lo
	pset altpl(curve.params);
	altpl.desc("altpl");

	// in this alt pset, set constpar Fixed...
	for (auto& pi : constpar) {
		Par* constpp = altpl.findp(pi);
		if (constpp == nullptr)
			throw runtime_error(vastr("uncertain parameter ",pi," has not been defined"));
		constpp->set_fixed();
	}

	// ... and set the uncertain parameters indep, verify they are AD
	for (auto& nm : uncnames) {
		Par* pp = altpl.findp(nm);
		if (pp == nullptr)
			throw runtime_error(vastr("uncertain parameter ",nm," has not been defined"));
		pp->set_indep();
		if (Ad::find(nm) < 0)
			throw runtime_error(vastr("uncertain parameter ",nm," is not AD"));
		uncpar.push_back(nm);
	}
	// create equations for them
	altpl.make_equations();
	// evaluate altpl parameters & set AD derivatives to 1.0
	// XXX cannot set AD par - only one set
	altpl.eval();

	// if factors were given set uncfac
	if (!factors.empty()) {
		if (factors.size() != uncpar.size()) {
			throw runtime_error(vastr(uncpar.size()," uncertainty parameters, but ",
				factors.size()," factors"));
		}
		uncfac = factors;
	}

	// create the 3 curves: uncertainty trace with altpl, hi and lo with gpset
	hi = new Flutcurve(curve.aid(), curve.cid()+"_hi");
	lo = new Flutcurve(curve.aid(), curve.cid()+"_lo");
	unc = new Flutcurve(curve.aid(), curve.cid()+"_unc", &altpl);

	// the interval curves (hi & lo) have the uncertain params
	// as constants at their limits
	// XXX should we just get the values in eval?
	for (size_t i=0; i<uncpar.size(); i++) {
		Par* lopar = lo->params.findp(uncpar[i]);
		double p = lopar->value();
		double delta = uncfac[i]*p;
		lopar->valuef(p-delta);
		Par* hipar = hi->params.findp(uncpar[i]);
		hipar->valuef(p+delta);
	}
}

void
Uncert::
update(vector<Interval>& xu, double coord) {
// update the uncertainty curve (unc) with xu, the update
// curves hi and lo with corresponding parameters
// Note the size of the Ivar in unc is different from those in hi & lo,
// but the params are the same
	Trace trc(2,"Uncert::update");
#ifdef NEVER // interval scheme
	unc->update(xu);
	for (auto& pi : unc->params) {
		double x = pi.second->value();
		Par* phi = hi->params.findp(pi.first);
		double xhi = phi->value();
		if (x > xhi)
			phi->valuef(x);
		Par* plo = lo->params.findp(pi.first);
		double xlo = plo->value();
		if (x < xlo)
			plo->valuef(x);
	}
#else // NEVER : interval scheme
	// lower bounds
	vector<Par*>& lo_par = lo->ivar_par();
	const vector<double>& lo_scale = lo->ivar_scale();
	assert(xu.size() == lo_par.size());
	for (size_t i=0; i<lo_par.size(); i++) {
		double x = xu[i].lower()*lo_scale[i];
		lo_par[i]->value(x);
		trc.dprint("set ",lo_par[i]->name,".lower to ",x);
	}
	// coord
	Par* coordp = lo->params.findp("coord");
	if (coordp != nullptr)
		coordp->value(coord);
	// upper bounds
	vector<Par*>& hi_par = hi->ivar_par();
	const vector<double>& hi_scale = hi->ivar_scale();
	for (size_t i=0; i<hi_par.size(); i++) {
		double x = xu[i].upper()*hi_scale[i];
		hi_par[i]->value(x);
		trc.dprint("set ",lo_par[i]->name,".upper to ",x);
	}
	// coord
	coordp = hi->params.findp("coord");
	if (coordp != nullptr)
		coordp->value(coord);
	// mark all parameters out-of-date...
	hi->params.touch();
	lo->params.touch();
	// ... then set the autodiff derivative parameters derivative spot to 1,
	// and evaluate all derived parameters
	// XXX params should be evaluated in intervals: need eval(lopar,hipar)
	hi->params.eval();
	lo->params.eval();
#endif // NEVER : interval scheme
}

vector<Interval>
Uncert::
compute_f(Point& pt) {
	Trace trc(1,"compute_f");
	// update curve lo and compute the function...
	lo->update(pt);
	size_t nx = pt.x.size();
	vector<double> jac(nx*nx, 0.0);
	size_t nf;
	vector<double> flo(nx, 0.0);
	lo->fcnjac(pt.x, flo, jac, nf);
	// ...update curve lo and compute the function
	hi->update(pt);
	vector<double> fhi(nx, 0.0);
	hi->fcnjac(pt.x, fhi, jac, nf);
	// form the interval return vector; make it nx long with the
	// last element zero to allow a square system
	//!! vector<Interval> rval(nx, 0.0);
	vector<Interval> rval(nf, 0.0);
	for (size_t i=0; i<nf; i++)
		rval[i] = Interval(flo[i], fhi[i]);

	trc.dprint("f_u width ",vector_width(rval),rval);
	return rval;
}

void
Uncert::
eval(Point& pt) {
// estimate the uncertainty at "pt" due to given uncertainties in
// certain parameters (uncpar)
	Trace trc(2,"Uncert::eval coord ",pt.coord);
	Flutcurve* curve{pt.fcv};
	size_t nx = curve->get_nx();
	size_t nf;
	// compute the jacobian at current point on "curve",
	// augmented with the tangent
	vector<double> fx(nx, 0.0);
	vector<double> jac(nx*nx, 0.0);
	vector<int> ipiv(nx,0);
	//!! curve->fcnjac(pt.x, fx, jac, pt.tangent, 0.0, nf);
	curve->fcnjac(pt.x, fx, jac, nf);
	SVD jsvd(nf,nx,&jac[0]);
	trc.dprint("rcond = ",jsvd.rcond);
	// LU factorize the Jacobian
	if (abs(pt.coord) == 1)
		trc.dprintm(nx,nx,nx,&jac[0],"jacobian on curve");
#ifdef NEVER // use svd inverse
	int info;
	int n{nx};
	dgefa_(&jac[0],&n,&n,&ipiv[0],&info);
	trc.dprint("dgefa info ",info);
	if (abs(pt->coord) == 1)
		trc.dprintm(nx,nx,nx,&jac[0],"jacobian factors");
#endif // NEVER : use svd inverse

#ifdef NEVER // use compute_f
	// convert x to xu
	size_t nxu = unc->get_nx();
	vector<Par*>& xpar = curve->ivar_par();
	//!! Tangent tanu(nxu, 0.0);
	Ivar xu(nxu, 0.0);
	for (size_t i=0; i<nx; i++) {
		int idx = unc->ivar_index(xpar[i]->name);
		assert(idx >= 0);
		//!! tanu[idx] = tan[i];
		xu[idx] = pt->x[i];
	}
	//!! trc.dprintm(nxu,1,nxu,tanu,"tanu: curve tangent in unc coords:");

	// set th uncert parameters in xu to the (fixed) values in curve->params
	for (auto nm : uncpar) {
		Par* cp = curve->findp(nm);
		int idx = unc->ivar_index(nm);
		assert(idx >= 0);
		double scale = unc->ivar_scale()[idx];
		xu[idx] = cp->value()/scale;
	}
	trc.dprint("unc ivar xu",xu);

	// compute the jacobian in the unc cs to get df/d(uncpar)
	size_t nfu;
	fx.resize(nxu);
	vector<double> jac_u(nxu*nxu, 0.0);
	unc->fcnjac(xu, fx, jac_u, nfu);
	trc.dprint("nf = ",nf,", nfu = ",nfu);

	// Now extract each uncpar column, create an Interval vector as
	//    f_u = {[-delta*f', delta*f']}
	// ie an approximation to the uncertainty in f assuming f(x) = 0
	// Sum them together to form the overall uncertainty in f(x)
	// get the nominal values and limits for all uncertain parameters
	// in continuation units (CU)
	// NOTE: dimension f_u nf but nf = nfu+1
	vector<Interval> f_u(nf, Interval(0.0,0.0));
	double fwidth{0.0};
	for (size_t k=0; k<uncpar.size(); k++) {
		int j = unc->ivar_index(uncpar[k]);
		double delta = abs(uncfac[k]*xu[j]);
		trc.dprint("delta ",k," = ",delta);
		for (size_t i=0; i<nfu; i++) {
			double df = delta*jac_u[IJ(i,j,nfu)];
			f_u[i] += Interval(-df, df);
			fwidth += 4.0*df*df;
		}
	}
	fwidth = sqrt(fwidth);
	trc.dprint("f_u width ",fwidth,f_u);
#endif // NEVER : use compute_f

	vector<Interval> f_u = compute_f(pt);

	// solve J*h = f_u
	vector<Interval> h = jsvd.solve(f_u);
	trc.dprint("h width ",vector_width(h),h);
	//!! dgesl_(&jac[0], &n, &n, &ipiv[0], &f_u[0]);
	//!! trc.dprint("h width ",vector_width(f_u),f_u);
	// add the uncertainty h to x to form x_u
	vector<Interval> x_u(nx, Interval(0.0));
	for (size_t i=0; i<nx; i++)
		x_u[i] = pt.x[i] + h[i];
		//!! x_u[i] = pt->x[i] + f_u[i];
	trc.dprint("x_u", x_u);
	// update the hi and lo curves with x_u
	this->update(x_u, pt.coord);
	// then add the parameters to the solns
	if (pt.coord > 0.0) {
		lo->back_solns();
		hi->back_solns();
	} else {
		lo->front_solns();
		hi->front_solns();
	}
}


#ifdef NEVER // use dge.c
#ifdef __cplusplus
//!! extern "C" {
#endif
// #include "f2c.h"

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
#endif // NEVER : use dge.c
