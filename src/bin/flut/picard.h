//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

#ifndef Picard_h
#define Picard_h

// Public interface: picard()
using Userfcn = int (*)(std::vector<double> const& x, std::vector<double>& f, std::vector<double>& jac);

int
picard (std::vector<double> const& xmin,
		std::vector<double> const& xmax,
		Userfcn fcn, int inint, int& nDtd);

//------------------------------------------------------------------
// classes used to implement topdeg()
// a Point is an n-vector of ints with a marker telling if
// it is an initial arbitrary enclosing point
class Point : public std::vector<double> {
public:
	// std::vector<int> idx;
	// std::vector<double> x;
	std::vector<double> f;
	std::vector<double> jac;
	// int sgndet;
	bool isarb;
	// constructors
	Point (std::vector<double> xs) : std::vector<double>(xs), isarb(false) {}
	Point (std::vector<double> xs, bool arb) : std::vector<double>(xs), isarb(arb) {}
	// Point (std::vector<double> xs, std::vector<int> idxs) :
	// 	std::vector<double>(xs) { idx = idxs; }

	bool operator==(const Point& p) const;
	std::string toString() const;
};


// a Face is a set of n Points (n-vectors of doubles)
// part of an n-dimensional Simplex, describing an
// (n-1)-dimensional hypersurface. Any n Points of
// a Simplex form a Face; a boundary Face is one where
// one coordinate of all n Point's are at a specified
// value, e.g xmin or xmax
class Face : public std::vector<Point> {
public:
	Face() {}
	Face(std::vector<Point> points) : std::vector<Point>(points) {}

	void detX(double* det);
	bool operator==(const Face& f) const;
	void permute();
	void plot(int faceno, std::string const& file, bool append) const;
	// 3D faces (2D surfaces in 3D) can be visualized in gnuplot:
	void plotgnu(int faceno) const;
	static void makegnucmd(int nbf);
	std::string toString() const;
};

// a Simplex comprises n+1 n-dimensional Points forming
// an nD hypersurface, and n+1 (n-1)D hypersurface Faces
class Simplex : public std::vector<Point> {
public:
	// constructor
	Simplex (std::vector<Point> const& points);

	std::vector<Face> faces() const;
	vector<Face> boundaryFaces (std::vector<double> const& lower,
			std::vector<double> const& upper);

	int orient(int dir);  // orient to have direction "dir"
	double orientation(int& expon) const;
	void permute();
	void plot(int simpno, std::string const& file, bool append) const;
#ifdef NEVER // unused
	double td(std::vector<double> const& xmin,
			std::vector<double> const& xmax,
		int (*fcn)(std::vector<double> const& x, std::vector<double>& f, int& sd),
		std::vector<int> nintvl, bool ignsign);
#endif // NEVER : unused
	std::vector<Face> cube_bndry();
};

// an n-dimensional Cube comprises n! Simplex's and has
// 2n boundary faces;
// these faces comprise 2n! Simplex boundary faces
class Cube : public std::vector<Simplex> {
public:
	// std::vector<double> xmin; 
	// std::vector<double> xmax; 
	double topdeg;
	Cube() : topdeg(0.0) {}
	Cube(std::vector<Simplex> s) {
		for (auto sp : s)
			push_back(sp);
		topdeg = 0.0;
	}

	// order() returns n
	size_t order() { return (*this)[0].size()-1; }

	// all faces, internal and boundary
	std::vector<Face> faces();

	// the 2n (n-1)-dimensional boundary Faces for this Cube
	std::vector<Face> boundary();

	// the limits of this cube
	vector<Point> limits();

	// compute the TD for this Cube:
	double td(std::vector<std::vector<double> > pts,
		Userfcn fcn, int& nadded);

#ifdef NEVER // unused
	// compute the sum for the (n+1)-D Picard problem
	double tdpicard(std::vector<double> const& xmin,
			std::vector<double> const& xmax,
			int (*fcn)(std::vector<double> const& x, std::vector<double>& f, int& sd),
		std::vector<int> nintvl);
#endif // NEVER : unused

	void plot(std::string const& file);
};  // Cube

std::ostream&
operator<<(std::ostream& s, const Point& t);
std::ostream&
operator<<(std::ostream& s, const Face& t);
std::ostream&
operator<<(std::ostream& s, const Simplex& t);
std::ostream&
operator<<(std::ostream& s, const Cube& t);

void
plotSimplices(std::vector<Simplex> simplices, std::string const& file);

Cube
proto(int n);

template<typename Type>
Type
determ (size_t n, Type* A, int& expon) {
// returns the determinant of (n,n) matrix A
	Trace trc(2,"determ",DatatypeTraits<Type>::id());
	Type rval{0};
	size_t i, j, kj;
	size_t nm1 = n - 1;

	expon = 0;
	static int depth = 0;
	if (depth == 0)
	depth++;

	if (n == 1) {
		rval = A[0];
	} else if (n == 2) {
		rval = (A[0]*A[3] - A[1]*A[2]);
	} else if (n == 3) {
		rval = A[0]*(A[4]*A[8] - A[5]*A[7]) - A[3]*(A[1]*A[8] - A[2]*A[7])
			+ A[6]*(A[1]*A[5] - A[2]*A[4]);
	} else if (n == 4) {
	 	Type d31 = A[5]*(A[10]*A[15] - A[11]*A[14]) - A[9]*(A[6]*A[15] - A[7]*A[14])
	 		+ A[13]*(A[6]*A[11] - A[7]*A[10]);
	 	Type d32 = A[1]*(A[10]*A[15] - A[11]*A[14]) - A[9]*(A[2]*A[15] - A[3]*A[14])
	 		+ A[13]*(A[2]*A[11] - A[3]*A[10]);
	 	Type d33 = A[1]*(A[6]*A[15] - A[7]*A[14]) - A[5]*(A[2]*A[15] - A[3]*A[14])
	 		+ A[13]*(A[2]*A[7] - A[3]*A[6]);
	 	Type d34 = A[1]*(A[6]*A[11] - A[7]*A[10]) - A[5]*(A[2]*A[11] - A[3]*A[10])
	 		+ A[9]*(A[2]*A[7] - A[3]*A[6]);
		rval = A[0]*d31 - A[4]*d32 + A[8]*d33 - A[12]*d34;
	} else {
		Type scale{1};
		std::vector<Type> M(nm1*nm1);
		for (j=0; j<n; j++) {
			size_t k = 0;
			for (kj=0; kj<n; kj++) {
				if (kj != j) {
					for (i=1; i<n; i++) {
						M[IJ(i-1,k,nm1)] = A[IJ(i,kj,n)];
					}
					k++;
				}
			}
			// rval += scale*A[IJ(0,j,n)]*determ(nm1, &M[0], exi);
			int exi;
			Type di = scale*A[IJ(0,j,n)]*determ(nm1, &M[0], exi);
			if (j == 0) {
				rval = di;
				expon = exi;
			} else {
				if (exi != expon)
					di *= pow(10, exi-expon); 
				rval += di;
			}
			scale = -scale;
		}
	}

	// adjust rval*10^expon so that 1 < rval < 10
	// only if Type is a floating-point type
	std::string tname(typeid(Type).name());
	trc.dprint("Type is ",tname);

	if (rval != 0 && (abs(Type(0.4)) != abs(Type(0)))) {
		while (abs(rval) < 1) {
			rval *= 10;
			expon -= 1;
		}
		while (abs(rval) > 10) {
			rval /= 10;
			expon += 1;
		}
	}

	depth--;
	if (depth == 0) {
		delete time;
	}
	trc.dprint("returning ",rval,"*10^",expon);
	return rval;
}

#endif // Picard_h
