//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

#ifndef PAC_H
#define PAC_H 1
// Pseudo-Arclength Continuation
// References - ref number used throughout pac.h and pac.cpp:
//
//   1) PITCON source code version 6.1 available on netlib as dpcon61.f
//
//   2) Rheinboldt, W.C., Numerical Analysis of Parameterized
//      Nonlinear Equations, The University of Arkansas Lecture
//      Notes in the Mathematical Sciences Vol. 7,
//      John Wiley & Sons, 1986.
// 
//   3) Rheinboldt, W.C., and Burkardt, J.V., "A Locally Parameterized
//      Continuation Process", ACM Trans Math Software, Vol 9 No. 2,
//      June 1983, pp 215-235. (rheinboldt1983locally.pdf)
//
//   4) Den Heijer, C. and Rheinboldt, W.C., "On Steplength Algorithms
//      for a Class of Continuation Methods", SIAM J. Numer. Analysis,
//      Vol. 18 No. 5, Oct 1981, pp 925-948
//
//   5) Rheinboldt, W.C. and Burkardt, J.V., "Algorithm 596: a program
//      for a locally parameterized continuation process", ACM Trans
//      Math Software, Vol 9 No. 2, June 1983, pp 236-241
//
//   6) A Fortran90 version of PITCON is available at
//        people.math.sc.edu/Burkardt/f_src/pitcon7/pitcon7.html
//
// Public interfaces:
//    Pac::sigue(Pac& to, const vect& constraint): extend a Pac to "to" using the
//       tangent and stepsize in "this"
//    homotopy(vect& x0, Fjac fjac, vect* results, vect* stepdata):
//       solve a system of n equations in n unknowns, given an initial guess at
//       a solution (x0) and a function that evaluates f(x)
//    nleig(complex<double>& s, vector<complex<double>>& x, Dfcn dmat,
//    		vect* results, vector<pair<string,string>>* params, vect* stepdata):
//    	given an estimate of a solution to the nonlinear eigenvalue problem
//    	D(s)x = 0, refine the estimate using homotopy
//    start.trace (Processfcn process) trace a curve beginning with
//       "start", continuing until "process" returns false.
//
// Aliases for functions used by public interface functions:
//    Fjacfcn: int fjac(const vect& x, vect& f, vect& jac)
//    Processfcn:  int process(Pac& from, Pac& to, Issue* issue)
// if the next step is to be constrained, set to.confcn
#include <complex>
#include <functional>
#include <string>
#include <utility>
#include <vector>

class Pac;

enum class Kind { angle, noconv, divergence, slow, stalled, zerotan, user_request };
class Issue {
public:
	std::string msg;
	Kind kind;
	Issue(std::string m, Kind k) : msg(m), kind(k) {}
	bool converged() {
		if (kind == Kind::angle || kind == Kind::stalled || kind == Kind::zerotan) return true;
		return false;
	}
	bool is_angle() { return kind == Kind::angle; }
	bool is_noconv() { return kind == Kind::noconv; }
	bool is_divergence() { return kind == Kind::divergence; }
	bool is_slow() { return kind == Kind::slow; }
	bool is_stalled() { return kind == Kind::stalled; }
	bool is_zerotan() { return kind == Kind::zerotan; }
};
std::ostream&
operator<<(std::ostream& s, const Issue* t);

// shorthand for vector of doubles
using vect = std::vector<double>;
// alias for the function that computes f(x) and the Jacobian
using Fjacfcn = std::function<int(const vect& x, vect& f, vect& jac)>;
// alias for a function that returns the correction needed for
// a constraint and it's vector
using Constraintfcn =
	std::function<double(const std::vector<double>& x, std::vector<double>& c)>;

// alias for functions which compute d(Dx)/dx, d(Dx)/dsigma, d(Dx)/dfreq for nleig
// solving the nonlinear eigenvalue problem D(s,x)x = 0 for an eigenpair (s,x)
using Dfcn =
	std::function<int(std::complex<double> s, const std::vector<std::complex<double>>& x,
		std::vector<double>& Dx,				// 2n double D*x
		std::vector<double>& Dsigma,		// 2n double d(Dx)/dsigma
		std::vector<double>& Dfreq,			// 2n double d(Dx)/dfreq
		std::vector<double>& Dxx)>;			// 2n,2n double d(Dx)/dx

// alias for processing results at each step of trace()
// returns int: +/0/- means continue/quit/retry
using Processfcn = std::function<int(Pac& from,Pac& to,Issue* issue)>;

// solve f(x) = 0
bool homotopy(vect& x0, Fjacfcn fjac, vect* results=nullptr,
		std::vector<std::pair<std::string,std::string>>* params=nullptr,
		vect* stepdata=nullptr);


class Det {
public:
	double coef{0.0};
	double exp{0.0};
	Det () = default;
	Det (double d, double e) : coef(d), exp(e) {}
	void equalize(Det const& p);	// make my exponent = p.exp
};

std::ostream&
operator<<(std::ostream& s, const Det& t);

//------------------------------------------------------------------
// class Step: computation of stepsize for pseudo-arclength continuation
//
// Public interface:
//   Issue* rheinboldt (Pac const& from, double& stepsize);
//
// a Step contains several items related to the computation
// of "stepsize", the optimal stepsize based on estimated radius
// of convergence as documented in references in reverse chronological order:
//   1) PITCON source code version 6.1 available on netlib as dpcon61.f
//
//   2) Rheinboldt, W.C., Numerical Analysis of Parameterized
//      Nonlinear Equations, The University of Arkansas Lecture
//      Notes in the Mathematical Sciences Vol. 7,
//      John Wiley & Sons, 1986.
// 
//   3) Rheinboldt, W.C., and Burkardt, J.V., "A Locally Parameterized
//      Continuation Process", ACM Trans Math Software, Vol 9 No. 2,
//      June 1983, pp 215-235.
//
//   4) Den Heijer, C. and Rheinboldt, W.C., "On Steplength Algorithms
//      for a Class of Continuation Methods", SIAM J. Numer. Analysis,
//      Vol. 18 No. 5, Oct 1981, pp 925-948
//
// Names used in Step for some variables from ref. 3:
//   deltak      (\delta_k eqn. 4.3) total corrector distance
//   hnorm       (4.10) last corrector distance
//   omegatilde  (\tilde{\omega} eqn. 4.10) hnorm/deltak
//   thetak      (\theta_k table I) convergence quality
//   wk          (2-norm of w^k, eqn. 4.16) curvature between previous
//               and current point
//   alphak      (\alpha_k, eqn. 4.16) angle between prev and current tangents
//   gammak      (\gamma_k, eqn 4.17a,b) predicted curvature between current
//               point and next point
//   epsk        (4.20) truncated ideal starting error
//   hk1         (h_k^1, eqn. 4.21) tentative predicted stepsize
//   hk2         (h_k^2, eqn. 4.22, 4.23) final predicted stepsize
//------------------------------------------------------------------
class Pac_variable {
public:
	std::string name;
	std::string desc;
	Pac_variable(std::string nm, std::string d) : name(nm), desc(d) {}
};

class Fv : public Pac_variable {
public:
	double value;
	Fv(std::string nm, std::string d, double v=0.0) :
		Pac_variable(nm,d), value(v) {}
};

class Iv : public Pac_variable {
public:
	int value;
	Iv(std::string nm, std::string d, int v=0) :
		Pac_variable(nm,d), value(v) {}
};

std::ostream&
operator<<(std::ostream& s, const Fv& t);

std::ostream&
operator<<(std::ostream& s, const Iv& t);

class Step {
public:
	// equation numbers refer to reference 3
	static std::vector<std::pair<std::string, std::string>> parnames;
#ifdef NEVER // new Fv
	Fv deltak{"deltak","corrector dist"};	// dist between predicted & corrected (4.3)
	Fv rcond{"rcond","recip cond of Jacobian"};	// repciprocal cond of Jacobian
	Fv hnorm{"hnorm","last correction norm"};	// norm of last correction
	Fv omegatilde{"omegatilde","hnorm/deltak"};		// ratio hnorm/deltak
	Fv quality{"quality","thetak 4.5"};				// thetak (4.5 & table 1)
	Fv secant{"secant","secant norm"};				// secant norm twix this and prev solution
	Fv alphak{"alphak","tan angle (rad)"};     	// angle (rad) between this and prev tangents
	Fv wk{"wk","curvature"};						// curvature between current and prev
	Fv gammak{"gammak","predicted curvature"};  	// predicted curvature on next step
	Fv epsk{"epsk","starting error"};				// truncated ideal starting error
	Fv fxnorm{"fxnorm","fx norm"};					// 2-norm of f(x)
	Fv xnorm{"xnorm","x norm"};						// 2-norm of x
	Fv secssratio{"secssratio","secant/stepsize"};	// secant/stepsize (should be near 1)
	Fv stepratio{"stepratio","stepsize/prev"};		// stepsize/prev stepsize
	Fv hk1{"hk1","1st est stepsize"};			// first est of stepsize
	Fv hk2{"hk2","2nd est stepsize"};			// second est of stepsize, before truncation
	Fv hk{"hk","predicted stepsize"};			// predicted optimal stepsize
	// integer variables (Iv)
	Iv niter{"niter","# iterations"};				// number of iterations to converge
	Iv conv_crit{"conv_crit","convergence crit"};	// see Pac::convergence()
   Iv red_conv{"red_conv","# reducts: converg."};	// no convergence
   Iv red_angle{"red_angle","# reducts: angle"};	// large angle between tangents
   Iv red_dsc{"red_dsc","# reductions: dsc"};		// Dsc: determinant sign change
	Iv svd{"svd","qr or svd 0/1"};				// qr/svd: 0/1
#else // NEVER // new Fv
	// parameter variables:
	double deltak{0.0};				// dist between predicted & corrected (4.3)
	double rcond{0.0};				// repciprocal cond of Jacobian
	double hnorm{0.0};				// norm of last correction
	double omegatilde{0.0};			// ratio hnorm/deltak
	double quality{0.0};				// thetak (4.5 & table 1)
	double secant{0.0};				// norm of secant between this and prev solution
	double alphak{0.0};     		// angle (rad) between this and prev tangents
	double wk{0.0};					// curvature between current and prev
	double gammak{0.0};	   		// predicted curvature on next step
	double epsk{0.0};					// truncated ideal starting error
	double fxnorm{0.0};				// 2-norm of f(x)
	double xnorm{0.0};				// 2-norm of x
	double secssratio{0.0};			// secant/stepsize (should be near 1)
	double stepratio{0.0};			// stepsize/prev stepsize
	double hk1{0.0};					// first est of stepsize
	double hk2{0.0};					// second est of stepsize, before truncation
	double hk{0.0};					// predicted optimal stepsize
	// these are actually integers:
	double coord{0.0};				// steps from the origin: duplicates Pac::coord
	double niter{0.0};				// number of iterations to converge
	double conv_crit{0.0};			// see Pac::convergence()
   double red_conv{0.0};			// reduced stepsize: no convergence
   double red_angle{0.0};			// reduced stepsize: angle between tangents
   double red_dsc{0.0};				// reduced stepsize: determinant sign change
	double svd{0.0};					// qr/svd: 0/1
	double endofdata;					// end of the doubles
	const double* cbegin() const { return &deltak; }
	const double* cend() const { return &endofdata; }
	double* begin() { return &deltak; }
	double* end() { return &endofdata; }
	int ndata() { return end()-begin(); }
#endif // NEVER // new Fv


	// constructors
	Step() {}
	// default copy constructor, assignment op ok

	double convergence_quality (int maxiter);
	double dtau(Pac const& from, Pac& to, int stepred);

	// 2 ways of computing convergence quality:
	static double coqual (int niter, double omegatilde);
	static double table1(int niter, double omegatilde);

}; // Step

// Specs for Pac, possibly set by the caller
class Pacspecs {
public:
	// factors controlling step size
	double curvaturefac{5.0};	// increase curvature, decrease stepsize
	int ecc{0};						// -1/0/1: never/default/always use minf
	double ecc_eps{1.0e-13};	// used in ecc()
	double epsabs{0.0};				// error tolerance for fxnorm
	double growthfac{3.0};
	double initial_stepsize{0.001};
	double maxangle{0.5235};		// angle between tans: 30 deg
	double maxcurvature{100.0};		// max curvature
	int maxiter{10};				// Newton - use 20 for Modified Newton
	int maxsteps{1000};
	double maxstepsize{0.2};
	std::vector<double> minf;		// smallest possible value for f_j j=1,nf
	double minstepsize{1.0e-12};
	double reductionfac{3.0};
};

std::ostream&
operator<<(std::ostream& s, const Step& t);

class Pac {
public:
	std::vector<double> f;
	std::vector<double> x;
	std::vector<double> jac;
	std::vector<double> tan;
	Fjacfcn fjac{nullptr};
	Constraintfcn confcn{nullptr};
	int coord{0};
	// "stepsize" is initialized (in Pac::sigue) to Step::hk, the predicted stepsize
	// for next point, and may be reduced for various reasons
	double stepsize{0.0};	// default: correct the start point first
	Det det;
	Pacspecs specs;			// factors controlling stepsize
	Step step;					// stepsize data
	double projnorm;			// norm of projection of tan onto nullspace (corrector)

	// Constructors
	Pac() = default;
	Pac(int nf, int nx, Fjacfcn fj);
	// XXX alternative: set starting x and tan
	Pac(const vect& x, const vect& tan, int nf, Fjacfcn fj);

	// extend "this" to "to" with an optional constraint
	Issue* sigue(Pac& to, vect* projectee=nullptr);
	Issue* corrector();
	Issue* corrector(Constraintfcn constraint);
	// compute extended convergence criteria
	double ecc();
	Issue* converged(double prev_fxnorm, double prev_hnorm);
	Issue* convergev(double prev_fxnorm, double prev_hnorm);
	// compute the optimal stepsize based on data collected during continuation,
	// and the previous and current point
	Issue* rheinboldt (Pac const& from);
	// trace a curve starting at "start" until "process" returns false,
	// put results in "results" and "stepdata"
	bool trace(Processfcn process, vect* projectee=nullptr, int direction=1);
	// compute the determinant of the Jacobian augmented with the tangent
	Det determinant();
	// check the Jacobian using finite-differences
	Issue* check_jac(double* xmin, double* xmax, Issue* reason,
			std::vector<std::string>* pnames=nullptr, double tol=0.1);
};	// Pac

std::ostream&
operator<<(std::ostream& s, const Pac& t);


// given an estimate of an eigenvalue/eigenvector for the complex nonlinear
// eigenvalue problem D(s)x = 0, refine the estimate using homotopy
int
nleig(std::complex<double>& s, std::vector<std::complex<double>>& x, Dfcn dmat,
		vect* results=nullptr,
		std::vector<std::pair<std::string,std::string>>* params=nullptr,
		vect* stepdata=nullptr);

// check the Jacobian by finite-differencing (vector version)
Issue*
checkjac(vect& x, vect& f, vect& jac, Fjacfcn fjac, double* xmin, double* xmax,
		Issue* reason, std::vector<std::string>* pnames=nullptr, double tol=0.1);
	
#endif // PAC_H
