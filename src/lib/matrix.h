//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

#ifndef Matrix_h
#define Matrix_h


#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

class Matrix;

inline size_t
IJ(size_t i, size_t j, size_t ld1) { return i + j*ld1; }
inline size_t
IJK(size_t i, size_t j, size_t k, size_t ld1, size_t ld2) {
	return i + ld1*(j + ld2*k); }
inline size_t
IJKL(size_t i, size_t j, size_t k, size_t l, size_t ld1, size_t ld2, size_t ld3) {
	return (i + ld1*(j + ld2*(k + ld3*l))); }

#include "Ad.h"
#include "blas.h"
#include "fptype.h"
#include "Par.h"
#include "Pz.h"
#include "pset.h"

// a ComplexMatrix is just a Matrix with is_complex()==true
using ComplexMatrix = Matrix;

int checksum(size_t n, char const* buf);

// functions for converting n-dimensional arrays to an index
std::vector<size_t> vec2mdim (size_t vector_index, std::vector<size_t> const& ld);

size_t indices2vi (std::vector<size_t> const& indices,
		std::vector<size_t> const& dim);

std::string
matrix_list_summary(const std::vector<Matrix*>& list);

// a Matrix stores it's values in a vector<double> (_data)
// If the matrix is complex the length of _data is 2*nr*nc,
// and it should be cast to complex<double>, otherwise it is nr*nc
class Matrix : public Fio {
protected:
	size_t nr;
	size_t nc;
	std::string _mid;
	std::string _desc;
	std::string _nonSymmDesc;
	std::vector<double> _data;
public:
	std::vector<Pz*> pz;

	// Constructors
	Matrix() { }
	Matrix (std::string const& id, std::string const& des,
			size_t m, size_t n, bool cmplx);
	Matrix(std::string const& mid);  // fetch constructor (will throw)

	// serializer constructor
	Matrix (Receiver& s);

	Matrix (const Matrix& a);  // copy constructor

	// clone or factory method ("virtual constructor")
	// virtual constructor must have body in Matrix.c for Fio vtable
	Matrix* clone() const { return new Matrix(*this); }
	// destructor
	~Matrix();

	Matrix& operator=(const Matrix& rhs); // assignment op

	std::vector<std::string> dependson(pset& plt) const;
	// this Matrix is constant iff it has no pz's or all dependent
	// parameters in ps are constant
	bool is_constant(pset& ps);

	// convert a real Matrix to complex
	void cast_complex();

	// Fio stuff
	//    get: must be static so that an instance is not needed
	//    id: also static for same reason (implemented in Matrix.c)
	//    vid: virtual function returning id()

	static bool regd;  // has been registered?
	static std::string id();
	std::string vid() const { return id(); }

	static Fio* get (Receiver& s) { return new Matrix(s); }
	void put (Sender& s) const;

	// fetch takes a regular-expression and returns all that match;
	// returns empty if no matches in catalog.
	static std::vector<Matrix*> fetch (const std::string& mid);
	void store();

	// XXX remove
	std::string basedatatype() const { return "n/a"; }
	std::string datatype() const { return "n/a"; }
	size_t eleSize() const { return 0; }

	// evaluate the matrix, putting the result into an (nr,nc) vector<complex<Ad>>
	bool eval(pset& plt, std::vector<std::complex<Ad>>& result) const;

	bool is_complex() const {
		if (nr*nc*2 == _data.size())
			return true;
		return false;
	}

	// accessors for matrix values or derivatives:
	double* elem() { return &_data[0]; }
	const double* elem() const { return &_data[0]; }
	std::complex<double>* celem() { return (std::complex<double>*)&_data[0]; }
	const std::complex<double>* celem() const {
		return (std::complex<double>*)&_data[0]; }
	// wrt: returns deriv wrt parameter "wrt" or matrix values if empty
	std::vector<double> values(std::string wrt, pset& ps);
	// cvalues always returns a complex array
	std::vector<double> cvalues(std::string wrt, pset& ps=gpset::get());

	// a reference to the vector<double> the matrix is stored as:
	std::vector<double>& data() { return _data; }
	std::vector<double> const& const_data() const { return _data; }

	double operator()(size_t i, size_t j=1) const {
	// return the (i,j) 1b element of a real matrix
		assert(nr*nc == _data.size());  // real matrices only
		assert(i>0 && i<=nr);
		assert(j>0 && j<=nc);
		return _data[IJ(i,j,nr)]; }
	// accessors for dimensions, mid, etc
	// the dimensions may be changed but you are responsible for insuring
	// the data array is changed also
	size_t rsize() const { return nr; }
	size_t rsize(int newr) {
		size_t rval = nr;
		nr = newr;
		return rval;
	}
	size_t csize() const { return nc; }
	size_t csize(int newc) {
		size_t rval = nc;
		nc = newc;
		return rval;
	}

	// get/set mid...
	std::string mid() const { return _mid; }
	bool mid(const std::string& m) { _mid = m; return true; }
	// ... and desc
	std::string const desc() const { return _desc; }
	bool desc(std::string const& des) { _desc = des; return true; }

	void plot (pset& pl, std::vector<Elem> const&,
			std::vector<std::string>& params,
			std::string& fixedparam, size_t nstep);
	void plot_apf (pset& pl, std::string const& filename,
			std::vector<Elem> const& ele, size_t nstep, std::string& param);

	void inflate(double eps);	// specializations in Matrix.c

	void transpose();

	// maintain a singlton collection of Matrix*
	//   collection()      return a reference to *the* collection
	//   insert(Matrix*)   add a Matrix
	//   find_desc(string) return a pointer to Matrix with description
	//   find_mid(string)  return a pointer to Matrix with mid
	//   inventory()       return (string) list of all Matrix*
	//   fetch_mids        returns a vector of pairs of (desc,mid)
	//   store_mids        save a list of all matrix mids and their
	//                     descriptions
	static std::vector<Matrix*>& collection();
	static void insert(Matrix* m);
	// find a matrix either by description or mid
	static Matrix* find_desc(std::string const& description);
	static Matrix* find_mid (std::string const& id);
	static std::string inventory() {return matrix_list_summary(collection());}
	static  std::vector<std::pair<std::string,std::string> >
		fetch_mids (std::string const& id);
	static void store_mids(std::string const& id);

#ifdef NEVER // Fma
	// fetch the matrices needed for visualization
	static void vzmatrices(const std::string& vzid, Matrix*& nodes,
			Matrix*& coords, Matrix*& conn, Matrix*& gct);
#endif // NEVER // Fma

	// operator<< prints a summary of the matrix
	friend std::ostream& operator<<(std::ostream& s, const Matrix& a);

};  // Matrix

// General-purpose template functions

template<typename Type>
void
relativeDiff (size_t nr, size_t nc, Type* a, Type* b, double* diff) {
/*
 * Compute the relative difference between two matrices.
 * Here the relative difference is taken as the magnitude
 * of the difference between elements divided by the average
 * of the infinity norms of the element's row and column; thus
 * this number will be smaller than a strict relative difference
 * (magnitude of the difference divided by average magnitude
 * of the elements) but it will emphasize the larger elements
 * in the matrices
 */
	double eps = sqrt(std::numeric_limits<double>::epsilon());
	size_t i, j;

	std::vector<double> rownorm;
	std::vector<double> colnorm;
	for (i=0; i<nr; i++) {
		rownorm.push_back(normalize_inf(nc, &a[IJ(i,0,nr)], nr));
	}
	for (j=0; j<nc; j++) {
		colnorm.push_back(normalize_inf(nr, &a[IJ(0,j,nr)], (size_t)1));
	}

	for (i=0; i<nr; i++) {
		for (j=0; j<nc; j++) {
			Type aij = a[IJ(i,j,nr)];
			Type bij = b[IJ(i,j,nr)];
			double norm = (rownorm[i] + colnorm[j])/2.0;
			double tol = eps*norm;
			double dij = std::abs(aij-bij);
			double x;
			if (dij < tol) {
				x = 0.0;
			} else {
				x = dij/std::max(std::max(std::abs(aij),std::abs(bij)), std::max(norm, 0.01));
			}

			diff[IJ(i,j,nr)] = x;
		}
	}
}

template<typename Type>
void
getElement(std::string const& desc, Elem const& el, Type* val);

// extract the (row,col) 0b element of matrix "desc" in Matrix::collection
std::complex<Ad>
matvij (pset& plt, std::string const& desc, int row, int col);

// compare 2 matrix ids
int
midcmp (std::string const& a, std::string const& b);

// templated functions for printing and summarizing vectors and matrices
namespace flaps {

// IndexReal & LargestMag are for summarizing vectors
class IndexReal {
public:
	int index; // one-based (1b)
	double v;
	IndexReal() : index(0), v(0.0) {}
	IndexReal(int i, double val) : index(i), v(val) {}
	bool operator<(IndexReal const& b) const { return (v < b.v); }
	bool operator>(IndexReal const& b) const { return (v > b.v); }
};

template<typename Type>
class LargestMag : public std::vector<IndexReal> {
public:
	LargestMag(size_t n, Type const* a) {
		for (size_t i=0; i<n; i++) {
			double x{std::abs(a[i])};
			push_back(IndexReal(i+1, x));
		}
		std::sort(begin(), end(), std::greater<IndexReal>());
	}
}; // LargestMag

template<typename DT>
std::string
summarize (size_t n, DT const* x, unsigned int maxchar=80) {
// Summarize a vector in a string with the largest element
// number followed by the second-largest, etc. Element
// numbers are followed by the abs value of the element in
// parentheses except for the first (largest) element
// which only gets the magnitude if it is not 1.0
// For example if the return value is:
//     11  122 (0.46)  131 (0.44)    8 (0.41)   15 (0.30)  130 (0.19)  121 (0.19)   13 (0.16)   19 (0.13)
// this means element 11 is the largest and has an absolute
// value of 1.0, element 122 is the 2nd largest and it's magnitude
// is 0.46 times the magnitude of the largest, and so on.
//
// Specializations are in matrix.c
	if (n == 0 || maxchar < 2) return std::string("");

	std::ostringstream ost;

	flaps::LargestMag<DT> ax(n, &x[0]);

	// Write the largest component ...
	int ndig = n<10000 ? (n<1000 ? (n<100 ? (n<10 ? 1 : 2) : 3) : 4) : 5;
	if (std::abs(ax[0].v-1.0) < 0.001)
		ost << std::setw(ndig) << ax[0].index;
	else
		ost << "[" << ax[0].v << "]" << std::setw(ndig) << ax[0].index;

	// If the largest value in the vector is zero quit
	if (ax[0].v < std::numeric_limits<double>::epsilon())
		return ost.str();

	// ... all other components depend on having enough room
	Form fm;
	fm.precision(2).fixed();
	for (size_t i=1; i<ax.size(); i++) {
		double val = ax[i].v/ax[0].v;
		//if (val < 0.01) break;
		ost << "  " << std::setw(ndig) << ax[i].index << " (" << fm(val) << ")";
		// ost << "  " << std::setw(ndig) << ax[i].index << " (" << val << ")";
		if (ost.str().size() >= maxchar)
			break;
	}
	return ost.str();
}

template<typename DT>
std::string
summarize (std::vector<DT> const& x, unsigned int maxchar=80) {
// summarize a vector: see comments above
	if (x.empty())
		return "empty";
	return flaps::summarize(x.size(), &x[0], maxchar);
}

template<typename DT>
std::string
summarize (std::vector<DT> const& x, unsigned int maxchar, bool isComplex) {
	if (x.empty())
		return "empty";
	if (isComplex) {
		std::complex<DT>* cd = (std::complex<DT>*)&x[0];
		size_t n = x.size()/2;
		return flaps::summarize(n, cd, maxchar);
	}
	return flaps::summarize(x.size(), &x[0], maxchar);
}

// summarize Specializations:
std::string
summarize (std::vector<int> const& x, unsigned int maxchar);
std::string
summarize (size_t n, int const* x, unsigned int maxchar);

} // namespace flaps

std::string
strArray(int nr, int nc, int lda, const std::complex<double>* a);

template<typename Type>
std::string
strArray(int nr, int nc, int lda, const Type* a) {
// write a matrix to a string
// XXX double, complex specializations in matrix.c

	std::ostringstream os;
	Form sci(6);
	sci.width(13);
	for (int i=0; i<nr; i++) {
		for (int j=0; j<nc; j++) {
			os << sci(a[IJ(i,j,lda)]) << " ";
			// os << a[IJ(i,j,lda)] << " ";
		}
		os << std::endl;
	}
	return os.str();
}

#ifdef NEVER // needs work: how to return path(s), input ostream, path
namespace flaps {
void
mprint (int nr,int nc,int lda,const std::vector<Ad>& a,const std::string& title);

void
mprint (int nr,int nc,int lda,
		const std::vector<std::complex<Ad>>& a,const std::string& title);

void
mprint (int nr, int nc, int lda, const int* a, const std::string& title);

void
mprint (int nr, int nc, int lda, const double* a, const std::string& title);

void
mprint (int nr, int nc, int lda,
		const std::complex<double>* a, const std::string& title);

} // namespace flaps
#endif // NEVER // needs work: how to return path(s), input ostream, path

#endif  // Matrix_h
