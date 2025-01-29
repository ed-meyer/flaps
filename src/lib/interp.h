//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

#ifndef Interp_h
#define Interp_h

#include <string>
#include <vector>

#include "Ad.h"
#include "einspline.h"
#include "fio.h"

// Various types of interpolation
// class interpolant  base class for all interpolation types
// class corner       smooths the transition from 2 interpolants
// class plinear      piecewise linear interpolation
// class spline       cubic spline interpolation (1x,1y)
// class vspline      cubic spline interpolation on a vector of doubles (1x,ny)
// class vspline2     cubic spline interpolation on a vector of doubles (2x,ny)
// class vspline3     cubic spline interpolation on a vector of doubles (3x,ny)
// class Interp       container for any interpolants

// abstract base class for all interpolation schemes
class interpolant : public Fio {

public:
	interpolant() {}
	virtual ~interpolant() = default;
	// things every interpolant should have
	virtual bool contains(double t) const = 0;
	virtual void limits(double& lo, double& hi) const = 0;
	virtual double eval(double t) const = 0;
	virtual Ad eval(Ad const& t) const = 0;
	virtual double deriv(double t) const = 0;

	// Fio required member & functions
	static bool regd;
	static std::string id() { return std::string("interpolant"); }
	virtual std::string vid() const = 0;
	static Fio* get (Receiver& s) { return nullptr; }
	virtual void put (Sender& s) const = 0;
};

class corner : public interpolant {
public:
	double x0, x1;  // 2 breakpoints
	std::vector<double> c;  // 4 cubic coeff

	corner(double x0, double y0, double yp0, double x1, double y1, double yp1);
	corner(Receiver& s);
	~corner() = default;

	bool contains(double t) const;
	void limits(double& lo, double& hi) const;
	double eval(double x) const;
	Ad eval(Ad const& x) const;
	double deriv(double t) const;

	// Fio required member & functions
	static bool regd;
	static std::string id() { return std::string("corner"); }
	std::string vid() const { return id(); }
	static Fio* get (Receiver& s) { return new corner(s); }
	void put (Sender& s) const;
}; // corner

// Piecewise linear
class plinear : public interpolant {
public:
	std::vector<double> x;
	std::vector<double> y;

	// constructors
	plinear(const std::vector<double>& xs, const std::vector<double>& ys) :
		x(xs), y(ys) {}
	plinear(Receiver& s);
	~plinear() = default;

	// things every interpolant should have
	bool contains(double t) const;
	void limits(double& lo, double& hi) const;

	double eval(double x) const;
	Ad eval(Ad const& x) const;
	double deriv(double t) const;


	// Fio required member & functions
	static bool regd;
	static std::string id() { return std::string("plinear"); }
	std::string vid() const { return id(); }
	static Fio* get (Receiver& s) { return new plinear(s); }
	void put (Sender& s) const;

	// printing
	friend std::ostream& operator<<(std::ostream& s, const plinear& t);
}; // plinear

class spline : public interpolant {
	// spline interpolation of a scalar wrt 1 scalar variable
	std::vector<double> xbp;  // nx breakpoints
	std::vector<double> ybp;
	std::vector<double> c0;	// nx coeff
	std::vector<double> c1;	// nx-1 coeff
	std::vector<double> c2;
	std::vector<double> c3;

public:
	spline() {}
	spline(const int n, const double *x, const double *y, double rho=0.0);
	spline(const std::vector<double>& x, const std::vector<double>& y, double rho=0.0);
	spline(Receiver& s);
	~spline() = default;

	// things every interpolant should have
	bool contains(double t) const;
	void limits(double& lo, double& hi) const;

	double eval (double t) const;
	Ad eval (Ad const& t) const;
	double deriv (double t) const;
	int interval(double t) const;  // which interval is t in?

	// Fio required member & functions
	static bool regd;
	static std::string id() { return std::string("spline"); }
	std::string vid() const { return id(); }
	static Fio* get (Receiver& s) { return new spline(s); }
	void put (Sender& s) const;
};

class vspline : public interpolant {
// B-spline interpolation of vectors with respect to a single real
// parameter, with serializable coefficients. If the data is complex<double>
// the real and imag parts are interpolated separately; i.e. construct
// with a y array of double the length.
// Coefficients are such that
//   f(t) = c0[i] + d*(c1[i] + d*(c2[i] + d*c3[i]))
// where d = t-bp[i] and bp[i] <= t <= bp[i+1]
// Note the coefficients are stored reversed from cubgcv
	std::vector<double> bp;  // nx breakpoints
	std::vector<std::vector<double> > yorig; // original (unsmoothed) data points
	// XXX c0 should also be (nx-1,ny)?
	std::vector<std::vector<double> > c0;  // (nx,ny) array of coeff
	std::vector<std::vector<double> > c1;	// (nx-1,ny)
	std::vector<std::vector<double> > c2;	// (nx-1,ny)
	std::vector<std::vector<double> > c3;	// (nx-1,ny)
	int extrap{0};
	double rho{0.0};
public:
	// constructors
	vspline() { extrap = 0; }
	// x: nx doubles
	// y: nx ny-vectors of doubles
	vspline(std::vector<double> const& x,
			std::vector<std::vector<double> > const& y,
			int extrap = 0, double rho = 0.0);
	// fetch constructor
	vspline(Receiver& s);

	~vspline() = default;

	// things every interpolant should have
	bool contains(double t) const;
	void limits(double& lo, double& hi) const;
	// Note: special eval function - vector<complex<Ad>> arg
	bool eval (Ad const& t, std::vector<std::complex<Ad>>& result) const;
	double eval (double t) const {
		throw std::runtime_error("eval(double) not implemented"); }
	Ad eval (Ad const& t) const {
		throw std::runtime_error("eval(Ad) not implemented"); }
	double deriv (double t) const {
		throw std::runtime_error("no deriv(double) for vspline"); }
	
	size_t getnx() const { return bp.size(); }
	size_t getny() const { return c0[0].size(); }


	int segment(double t) const;  // which segment is t in?
	int segment(Ad const& t) const;

	friend std::ostream& operator<<(std::ostream& s, const vspline& t);

	// Fio required member & functions
	static bool regd;
	static std::string id() { return std::string("vspline"); }
	std::string vid() const { return id(); }
	static Fio* get (Receiver& s) { return new vspline(s); }
	void put (Sender& s) const;

};  // class vspline

class vspline2 : public interpolant {
// B-spline interpolation of vectors with respect to 2 real
// parameters, with serializable coefficients. If the data is complex<double>
// the real and imag parts are interpolated separately; i.e. construct
// with a data array of double the length.
public:
	NUBasis x_basis;
	NUBasis y_basis;
	std::vector<NUBspline_2d> splines;
	size_t ndata{0};
	// constructors
	// x: nx doubles
	// y: ny doubles
	// data:  nx*ny vectors of doubles of length nrow*ncol (real matrices)
	//        or 2*nrow*ncol (complex)
	vspline2() = default;
	vspline2(std::vector<double> const& x,
			std::vector<double> const& y,
			std::vector<std::vector<double> > const& data);
	// fetch constructor
	vspline2(Receiver& s);

	~vspline2() = default;

	// things every interpolant should have
	bool contains(double t) const;
	void limits(double& lo, double& hi) const;
	void limitx(std::vector<double> grid, Ad& x, double& lo, double& hi) const;
	// Note: special eval function - vector<complex<Ad>> arg
	bool eval (std::vector<Ad> const& t, std::vector<std::complex<Ad>>& result) ;
	double eval (double t) const {
		throw std::runtime_error("eval(double) not implemented"); }
	Ad eval (Ad const& t) const {
		throw std::runtime_error("eval(Ad) not implemented"); }
	double deriv (double t) const {
		throw std::runtime_error("no deriv(double) for vspline2"); }
	
	size_t getnx() const { return x_basis.grid.size(); }
	size_t getny() const { return y_basis.grid.size(); }
	size_t getndata() const { return ndata; }

	friend std::ostream& operator<<(std::ostream& s, const vspline2& t);

	// Fio required member & functions
	static bool regd;
	static std::string id() { return std::string("vspline2"); }
	std::string vid() const { return id(); }
	static Fio* get (Receiver& s) { return new vspline2(s); }
	void put (Sender& s) const;

};  // class vspline2

class vspline3 : public interpolant {
// B-spline interpolation of vectors with respect to 3 real
// parameters, with serializable coefficients. If the data is complex<double>
// the real and imag parts are interpolated separately; i.e. construct
// with a data array of double the length.
public:
	NUBasis x_basis;
	NUBasis y_basis;
	NUBasis z_basis;
	std::vector<NUBspline_3d> splines;
	size_t ndata{0};
	// constructors
	// x: nx doubles
	// y: ny doubles
	// z: nz doubles
	// data:  nx*ny*nz vectors of doubles of length nrow*ncol (real matrices)
	//        or 2*nrow*ncol (complex)
	vspline3() = default;
	vspline3(std::vector<double> const& x,
			std::vector<double> const& y,
			std::vector<double> const& z,
			std::vector<std::vector<double> > const& data);
	// fetch constructor
	vspline3(Receiver& s);

	//!! ~vspline3();  // impl in interp.c to avoid missing vtable
	~vspline3() = default;

	// things every interpolant should have
	bool contains(double t) const;
	void limits(double& lo, double& hi) const;
	void limitx(std::vector<double> grid, Ad& x, double& lo, double& hi) const;
	// Note: special eval function - vector<complex<Ad>> arg
	bool eval (std::vector<Ad> const& t, std::vector<std::complex<Ad>>& result) ;
	double eval (double t) const {
		throw std::runtime_error("eval(double) not implemented"); }
	Ad eval (Ad const& t) const {
		throw std::runtime_error("eval(Ad) not implemented"); }
	double deriv (double t) const {
		throw std::runtime_error("no deriv(double) for vspline"); }
	
	size_t getnx() const { return x_basis.grid.size(); }
	size_t getny() const { return y_basis.grid.size(); }
	size_t getnz() const { return z_basis.grid.size(); }
	size_t getndata() const { return ndata; }

	friend std::ostream& operator<<(std::ostream& s, const vspline3& t);

	// Fio required member & functions
	static bool regd;
	static std::string id() { return std::string("vspline3"); }
	std::string vid() const { return id(); }
	static Fio* get (Receiver& s) { return new vspline3(s); }
	void put (Sender& s) const;

};  // class vspline3

// an Interp is a collection of interpolants for doing a different type
// of interpolation in each interval 
class Interp : public Fio {
	std::string iid_;
	std::vector<std::string> parnames_;
	std::string dataname_;
	std::vector<interpolant*> segments;
public:
	// constructors
	Interp() = default;
	Interp(const std::string& iid, const std::string& xn, const std::string& yn) :
		iid_(iid), parnames_{xn}, dataname_(yn) {}
	Interp(const std::string& iid, const std::vector<std::string>& pnames,
			const std::string& dname) : iid_(iid), parnames_(pnames), dataname_(dname) {}
	~Interp() = default;

	Interp(Receiver& s);					// de-serialize constructor
	Interp(const std::string& iid);	// fetch constructor

	std::string iid() const { return iid_; }
	// dependent parameter names:
	std::vector<std::string> dependson() const { return parnames_; }
	std::string dataname() const { return dataname_; }

	// XXX always push a *copy* of the interpolant in case it goes out
	// of scope?
	void add(interpolant* t) { segments.push_back(t); }
	void add(std::vector<interpolant*> t);
	void clear() { segments.clear(); }

	void limits(double& lo, double& hi) const;

	double eval (double t) const;
	Ad eval (Ad const& t) const;
	// Note: special eval function for vspline - vector<complex<Ad>> arg
	bool eval (std::vector<Ad> const& t, std::vector<std::complex<Ad>>& result);


	// Fio required member & functions
	void store();
	static Interp* fetch(const std::string& iid);
	static bool regd;
	static std::string id() { return std::string("Interp"); }
	std::string vid() const { return id(); }
	static Fio* get (Receiver& s) { return new Interp(s); }
	void put (Sender& s) const;

	void plot(const std::string& plotfile, int nstep);

	friend std::ostream& operator<<(std::ostream& s, Interp& t);
}; // Interp

// convenience functions
//
// piecewise-linear interpolation with smoothed corners
Interp*
plsmooth(const std::vector<double>& x, const std::vector<double>& y,
		double width, const std::string& iid="",
		const std::string& xname="", const std::string& yname="");

// create corners between each internal junction in a piecewise-linear
// interpolant
std::vector<corner*>
plinear_smooth(const plinear& pl, double width);

#endif // Interp_h
