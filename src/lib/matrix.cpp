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
#include "exim.h"
#include "lexer.h"
#include "matrix.h"
#include "Pz.h"
#include "trace.h"
#include "transpose.h"


using namespace std;

// register classes for Fio::get/put serialization
bool Matrix::regd = Fio::register_class(Matrix::id(), Matrix::get);


string
Matrix::
id() {
// fio serialization identifier
	return "Matrix";
}

int
checksum(size_t n, char const* buf) {
	long int rval = 0;
	for (size_t i=0; i<n; i++)
		rval += buf[i];
	return (int)rval;
}

vector<size_t>
vec2mdim (size_t indx, vector<size_t> const& ld) {
// Convert a 0b index (indx) into a vector into a set of
// m indices into the equivalent m-dimensional matrix,
// given the number of elements in each dimension (ld).
// Useful for taking all possible combinations of a set
// of m vectors, where "ld" contains the length of each
// vector.
	vector<size_t> rval;
	if (ld.empty())
		return rval;
	for (size_t i=0; i<ld.size()-1; i++) {
		rval.push_back(indx%ld[i]);
		indx /= ld[i];
	}
	rval.push_back(indx);
	return rval;
}

size_t
indices2vi (std::vector<size_t> const& indices,
		std::vector<size_t> const& dim) {
// Performs the reverse operation from vec2mdim()
// Given a set of indices into a p-dimensional matrix
// and the size of each dimension (dim), returns
// the equivalent index into a 1-dimensional array.
// This allows a chunk of memory to be treated as
// a p-dimensional array.
// For example if p=2, indices[0]=i and indices[1]=j
// rval = i + j*dim[0]
// Note that only p-1 elements of dim are referenced
	size_t rval = 0;
	// start from the last dimension and work backwards
	for (size_t i=indices.size()-1; i>0; i--) {
		rval += indices[i];
		rval *= dim[i-1];
	}
	rval += indices[0];
	return rval;
}

/*------------------------------------------------------------------
 * Matrix implementation
 *------------------------------------------------------------------*/

Matrix::
Matrix (std::string const& id, std::string const& des,
		size_t m, size_t n, bool cmplx) {
// general constructor
	nr = m;
	nc = n;
	_desc = des;
	_mid = id;

	size_t nn = m*n;
	if (cmplx) {
		nn *= 2;
	}
	_data = vector<double>(nn,0.0);
}

Matrix::
Matrix (Matrix const& a) {
// copy constructor
	nr = a.nr;
	nc = a.nc;
	_desc = a.desc();
	_mid = a.mid();

	_data = a.const_data();
	for (size_t i=0; i<a.pz.size(); i++) {
		pz.push_back(a.pz[i]->clone());
	}
}

Matrix::
Matrix(std::string const& mid) {
// fetch constructor
// Throws runtime_error exception if the matrix is not available -
// use Matrix::fetch() if you want a nullptr instead of exception
	T_(Trace trc(2,"Matrix fetch constructor");)

	if (mid.empty()) {
		string exc{vastr("attempt to fetch a matrix with a blank name")};
		throw runtime_error(exc);
	}

	vector<string> mids = fio::catalog(mid);
	if (mids.empty()) {
		string exc = vastr ("matrix \"",mid,"\" is not available");
		throw runtime_error(exc);
	}

	Receiver file(mid);
	Fio* op = Fio::get(file);
	Matrix* rval = dynamic_cast<Matrix*>(op);
	// check it is the right type
	if (rval == nullptr) {
		string exc = vastr("reading \"",mid,
			"\", expecting a Matrix, got a ",op->vid());
		T_(trc.dprint("throwing exception: ",exc);)
		throw runtime_error(exc);
	}
	// XXX need move assignment
	*this = *rval;
	delete rval;
	T_(trc.dprint("returning ",*this);)
}

vector<Matrix*>
Matrix::
fetch(const string& re) {
// Fetch one or more matrices based on regular-expression "re";
// returns empty if only one matrix matches (e.g. if "re" is not
// a regular expression) and it is not available, otherwise if
// more than one match and any one is not available a runtime_error
// exception is thrown.
	T_(Trace trc(2,"Matrix::fetch ",re);)
	vector<Matrix*> rval;
	// get a list of the matrices to fetch
	vector<string> mids = fio::catalog(re);
	T_(trc.dprint("catalog returned ",mids);)
	for (auto& mid : mids) {
		// mid may not be a Matrix (e.g. Fma) so catch exceptions
		try {
			rval.push_back(new Matrix(mid));
		} catch (runtime_error& s) {
#ifdef NEVER // new behavior
			if (mids.size() == 1)
				return rval;
			else
				throw s;
#else // NEVER // new behavior
			T_(trc.dprint("not adding ",mid," to rval: ",s.what());)
#endif // NEVER // new behavior
		}
	}
	return rval;
}


Matrix::
Matrix (Receiver& s) {
// Matrix deserialize constructor (receiver)
// _data holds the matrix values, which are modified by the
// parameterizations (Pz).
	T_(Trace trc(2,"Matrix::receiver constructor");)
	s.serialize (_mid);
	s.serialize (_desc);
	s.serialize (nr);
	s.serialize (nc);
	_data.clear();
	s.serialize(_data);
	T_(trc.dprintv(_data,_mid+" deserialized");)
	// check for complex
	// compare with size of data (d) read
	size_t nn = nr*nc;
	if (_data.size() != nn && _data.size() != 2*nn) {
		throw runtime_error(vastr("deserializing a (", nr, ", ", nc,
				") Matrix: ", _data.size(), " elements"));
	}
	// Read a vector of Pz*
	std::vector<Fio*> tmp;
	getVFio(s, tmp);
	for (size_t i=0; i<tmp.size(); i++) {
		Pz* p = dynamic_cast<Pz*>(tmp[i]);
		if (p)
			pz.push_back(p);
		else {
			throw runtime_error(vastr("expecting a Pz, got a \"",tmp[i]->vid(),'\"'));
		}
	}
	T_(trc.dprint("returning ",*this);)
}

Matrix::
~Matrix() {
// destructor
	for (size_t i=0; i<pz.size(); i++)
		delete pz[i];
}

void
Matrix::
put (Sender& s) const {
// Matrix serializer (sender)
// _data holds the matrix values, which are modified by the
// parameterizations (Pz).
	T_(Trace trc(2,"Matrix::put ", *this);)

	s.serialize (_mid);
	s.serialize(_desc);
	s.serialize (nr);
	s.serialize (nc);
	// serialize the _data
	s.serialize(_data);
	// serialize each Pz
	std::vector<Fio*> tmp;
	for (size_t i=0; i<pz.size(); i++)
		tmp.push_back(static_cast<Fio*>(pz[i]));
	putVFio(s, tmp);
}

void
Matrix::
store() {
	Sender file(mid());
	Fio::put(file, *this);
}

Matrix&
Matrix::
operator=(Matrix const& rhs) {
	T_(Trace trc(2,"Matrix::operator=(&)");)
	if (this != &rhs) {
		nr = rhs.nr;
		nc = rhs.nc;
		_mid = rhs.mid();
		_desc = rhs.desc();
		_data = rhs.const_data();

		for (size_t i=0; i<rhs.pz.size(); i++)
			pz.push_back(rhs.pz[i]->clone());
	}
	return *this;
}

bool
Matrix::
eval(pset& plt, vector<complex<Ad>>& result) const {
// Evalutate this matrix and put the result into "result" which
// may be larger than "this" and is complex where "this" may be real
	T_(Trace trc(2,"Matrix::eval ",this->mid());)
	bool rval(true);
	//!! size_t nr = result.nr;
	//!! size_t nc = result.nc;
	size_t m = this->rsize();
	size_t n = this->csize();
	//!! bool samesize = ((nr == m) && (nc == n));
	bool samesize = (m*n == result.size());

	T_(trc.dprint(this->pz.size()," pzs, samesize? ",samesize?"y":"n");)

	// zero the result array
	size_t nt = 2*result.size()*Ad::ndata();
	//!! std::fill(result.begin(), result.end(), complex<Ad>(0.0));
	// XXX need contiguous elements
	for (size_t i=0; i<result.size(); i++) {
		Ad::real(result[i]).zero();
		Ad::imag(result[i]).zero();
	}	

	// ok if "result" is larger but not smaller
	if (!samesize) {
		if (result.size() < m*n)
			throw runtime_error(vastr("evaluating ",this->mid(),
					": result (", result.size(), ") is smaller than matrix (",
					m, ", ", n, ")"));
		// smaller matrices only for square
		if (m != n)
			throw runtime_error("only square matrices allowed with extra eqns");
	}

	// copy the matrix to result, then if there are no pz we are done
	// note: _data is real but may represent complex
	// Only need to do this if there are no generative Pz (e.g. IntPz, RFAPz)
	size_t nr = sqrt(result.size());
	bool needcopy{true};
	for (auto pi : this->pz) {
		if (dynamic_cast<IntPz*>(pi) != nullptr || dynamic_cast<RFAPz*>(pi) != nullptr)
			needcopy = false;
	}

	if (needcopy) {
		if (this->is_complex()) {
			//!! result.copy_complex (&_data[0], m, n, nr, nc);
			const complex<double>* cp = this->celem();
			if (!samesize) {
				for (size_t j=0; j<n; j++)
					for (size_t i=0; i<m; i++) {
						Ad::real(result[IJ(i,j,nr)]).value(cp[IJ(i,j,m)].real());
						Ad::imag(result[IJ(i,j,nr)]).value(cp[IJ(i,j,m)].imag());
					}
			} else {
				for (size_t i=0; i<m*n; i++) {
					Ad::real(result[i]).value(cp[i].real());
					Ad::imag(result[i]).value(cp[i].imag());
				}
			}
		} else {
			const double* rp = this->elem();
			if (!samesize) {
				for (size_t j=0; j<n; j++)
					for (size_t i=0; i<m; i++) {
						//!! result[IJ(i,j,nr)] = rp[IJ(i,j,m)];
						Ad::real(result[IJ(i,j,nr)]).value(rp[IJ(i,j,m)]);
					}
			} else {
				for (size_t i=0; i<m*n; i++)
					Ad::real(result[i]).value(rp[i]);
			}
		}
	}

	// evaluate all Pz
	if (this->pz.empty()) {
		T_(trc.dprint("no pzs: copied data to result");)
		return rval;
	}

	// XXX all Pz eval expect complex result - add ones for real matrices
	// that expect real result?
	try {
		if (samesize) {
			for (size_t i=0; i<this->pz.size(); i++) {
				this->pz[i]->eval(plt, result, nr, nc);
			}
		} else {
			// "result" is larger than "this": create a temp array
			// to hold the evaluated matrix...
			vector<complex<Ad>> work(m*n);
			for (size_t i=0; i<this->pz.size(); i++) {
				this->pz[i]->eval(plt, work, m, n);
			}
			// ...then copy work to result a column at a time
			T_(trc.dprint("copying ",n," columns, ",nt," doubles/column");)
			for (size_t j=0; j<n; j++)
				for (size_t i=0; i<m; i++)
					result[IJ(i,j,nr)] = work[IJ(i,j,m)];
		}
	} catch (const std::range_error& s) {
		throw;  // rethrow
	} catch (const std::runtime_error& s) {
		throw runtime_error(vastr("evaluating ",this->mid(),": ",s.what()));
	}

	static int visit{0};
	if (visit++ == 0)
		T_(trc.dprintm(nr,nc,nr,result,this->mid());)
	return rval;
}

std::vector<double>
Matrix::
values (std::string wrt, pset& ps) {
// Returns matrix values (wrt=="0" or empty), or the derivative
// of the matrix with respect to parameter "wrt"
	T_(Trace trc(2,"Matrix::values");)

	size_t m = rsize();
	size_t n = csize();
	size_t mn = m*n;
	if (this->is_complex()) mn *= 2;
	std::vector<double> rval(mn);

	// if the matrix is constant just return the _data array or
	// zero if wrt != "0"
	if (this->pz.empty()) {
		T_(trc.dprint("quick return: constant");)
		if (wrt.empty() || wrt == "0")
			return _data;
		else
			return rval;
	}

	// eval() always expects a complex array
	vector<complex<Ad>> work(m*n);
	this->eval(ps, work);

	// partial wrt "wrt" = sum(\partial A/\partial p_i \partial p_i/\partial w)
	// extract the desired values
	if (this->is_complex()) {
		extract (work, wrt, (complex<double>*)&rval[0]);
	} else {
		std::vector<double> tmp(2*m*n);
		extract (work, wrt, (complex<double>*)&tmp[0]);
		blas::copy (m*n, &tmp[0], 2, &rval[0], 1);
	}
	return rval;
}

std::vector<double>
Matrix::
cvalues (std::string wrt, pset& ps) {
// Returns matrix values or the derivative of the matrix
// with respect to parameter "wrt" cast to complex
	vector<double> vals = this->values(wrt, ps);
	if (!this->is_complex()) {
		vector<double> rval(2*vals.size(), 0.0);
		blas::copy(vals.size(), &vals[0], 1, &rval[0], 2);
		return rval;
	}
	return vals;
}

bool
Matrix::
is_constant(pset& ps) {
// Is this Matrix constant? i.e. are all it's dependent parameters
// constant in ps, or it has no dependent parameters

	// if there are no pz's or only a Pz_const it is constant
	if (pz.empty() || (pz.size() == 1 &&
			dynamic_cast<Pz_const*>(pz[0]) != nullptr)) {
		return true;
	}
	// check all dependent parameters in ps
	std::vector<std::string> dep = this->dependson(ps);
	for (auto& di : dep) {
		Par* pp = ps.findp(di);
		// ok if a parameter has not been defined
		if (pp != nullptr && !pp->is_constant()) {
			return false;
		}
	}
	return true;
}

void
Matrix::
cast_complex() {
// cast a matrix from real to complex
// XXX what about the pz's?
	if (is_complex())
		return;
	// throw exception if attempting to cast parameterized matrix
	// if (!pz.empty()) {
	// 	std::ostringstream os;
	// 	os << "attempt to cast parameterized matrix to complex: "
	// 		<< summary();
	// 	throw Err(os.str());
	// }
	std::vector<double> newdata(2*_data.size(), 0.0);
	for (size_t i=0; i<_data.size(); i++)
		newdata[2*i] = _data[i];
	_data = newdata;
	assert(this->nr*this->nc == _data.size()/2);
}

vector<string>
Matrix::
dependson(pset& plt) const {
// returns a list of the names of all parameters this matrix is
// a function of
	T_(Trace trc(1,"Matrix::dependson");)
	vector<string> rval;
	for (auto& pzi : this->pz) {
		vector<string> di = pzi->dependson(plt);
		for (auto& dj : di) {
			if (std::find(rval.begin(), rval.end(), dj) == rval.end())
				rval.push_back(dj);
		}
	}
	return rval;
}

void
Matrix::
transpose() {
// Transpose a real or complex matrix "in-place".
	size_t m = this->rsize();
	size_t n = this->csize();
	if (this->is_complex()) {
		vector<complex<double> > cdata(m*n);
		flaps::transpose(m, n, this->celem(), &cdata[0]);
		blas::copy(m*n, &cdata[0], 1, this->celem(), 1);
	} else {
		vector<double> rdata(m*n);
		flaps::transpose(m, n, this->elem(), &rdata[0]);
		blas::copy(m*n, &rdata[0], 1, this->elem(), 1);
	}
	this->nr = n;
	this->nc = m;
}

void
Matrix::
plot_apf (pset& pl, string const& filename, vector<Elem> const& ele,
		size_t nstep, string& plotpar) {
// Create an Apf file of specified elements of a matrix that
// is a function of parameter "plotpar"
	T_(Trace trc(1,"Matrix::plot_apf");)
	size_t i, j;
	size_t nr = this->rsize();
	size_t nc = this->csize();
	vector<Ad> work(2*nr*nc);
	std::ostringstream os;

	T_(trc.dprint("plotting ",ele.size()," of ",*this," against ",plotpar);)

	//!! Par* fpsplotpar = gpset::find(plotpar);
	Par* fpsplotpar = pl.findp(plotpar);
	assert(fpsplotpar != nullptr);
	assert(fpsplotpar->has_min());
	assert(fpsplotpar->has_max());
	assert(nstep > 1);

	double xmin = fpsplotpar->min();
	double xmax = fpsplotpar->max();
	double xdel =  (xmax - xmin)/(nstep-1);

	// each ele is a separate curve
	vector<string> names;
	names.push_back(fpsplotpar->name);
	names.push_back("real");
	if (this->is_complex())
		names.push_back("imag");
	vector<complex<Ad>> result(nr*nc);
	vector<complex<double> > val(nr*nc);
	bool append{false};
	for (i=0; i<ele.size(); i++) {
		vector<double> x;
		vector<double> y;
		vector<double> z;
		for (j=0; j<nstep; j++) {
			double xval = xmin + j*xdel;
			x.push_back(xval);
			fpsplotpar->value(xval);
			//!! :this->eval(gpset::get(), result);
			this->eval(pl, result);
			extract(result, "", &val[0]);
			complex<double> cv = val[ele[i].index(nr)];
			y.push_back(cv.real());
			if (this->is_complex())
				z.push_back(cv.imag());
		}
		vector<vector<double>> xs;
		xs.push_back(x);
		xs.push_back(y);
		if (this->is_complex())
			xs.push_back(z);
		// string cid{vastr(this->mid(), ele[i])};
		string cid{vastr(ele[i].row+1,'.',ele[i].col+1)};
		Apf::plot(filename, cid, xs, names, append);
		append = true;
	}
} // plot_apf

void
Matrix::
plot (pset& pl, std::vector<Elem> const& ele, vector<string>& params,
			string& fixedparam, size_t nstep) {
// Plot specified elements of this matrix varying "params" parameters
// and holding "fixedparam" parameters constant
// If there is only one parameter name in "params" the data is
// written in ASCII plot-file format (.apf)
	T_(Trace trc(1,"Matrix::plot");)
	std::ostringstream os;

	T_(trc.dprint("fixed par: ",fixedparam.empty() ? "none" :fixedparam);)

	vector<string> toplot{params};

	// if none given use all parameters in the matrix pz's
	if (toplot.empty()) {
		vector<string> dep = this->dependson(pl);
		for (auto& depi : dep) {
			if (depi != fixedparam) {
				toplot.push_back(depi);
			}
		}
	}
	if (toplot.empty()) {
		T_(trc.dprint("quick return: no dependent parameters and none given");)
		return;
	}

	// the fixed parameter may have altvals: plot at each of them
	size_t nfixed{1};
	Par* stdfixedparam{nullptr};
	if (!fixedparam.empty()) {
		stdfixedparam = pl.findp(fixedparam);
		nfixed = stdfixedparam->altval.size();
		if (nfixed == 0)
			nfixed = 1;
	}

	// Change the params in the pl to Indep so the matrix is not constant
	for (auto& nm : toplot) {
		Par* pp = pl.findp(nm);
		pp->set_indep();
	}

	// if 3D plotting reduce the number of steps
	if (params.size() == 2)
		nstep = (int)sqrt((double)nstep);

	// a different file for each fixed altval?
	for (size_t j=0; j<nfixed; j++) {
		// base file name: mid.apf
		// or mid.j.apf if nfixed > 1
		string filename(this->mid());
		if (nfixed > 1) {
			os.str("");
			os << '.' << fixedparam << '_' << stdfixedparam->altval[j];
			filename += os.str();
		}

		if (!fixedparam.empty() && nfixed > 1) {
			stdfixedparam->valuef(stdfixedparam->altval[j]);
			pl.touch();
			pl.eval();
		}

		if (toplot.size() == 1) {
			plot_apf(pl, filename, ele, nstep, toplot[0]);
		} else if (params.size() == 2) {
			// XXX gnuplot
			// if (getenv("PLOTDX"))
				// plotDx(basefilename, ele[i], nstep, params);
			// else
				// plotGeomview(basefilename, ele[i], nstep, params);
			cerr << "multiple-parameter plotting is disabled\n";
		}
	}
}  // plot

void
Matrix::
inflate(double eps) { /* do nothing */ }

// end of Matrix implementation
// end of Matrix specializations

ostream&
operator<<(ostream& s, const Matrix& a) {
// print a summary of "a": mid (desc) nr by nc (complex) parameterized by ...
	ostringstream ts;
	int nr = a.rsize();
	int nc = a.csize();
	if (a.desc().empty())
		s << a.mid() << " " << nr << " by " << nc;
	else
		s << a.mid() << " (" << a.desc() << ") " << nr << " by " << nc;
	if (a.is_complex())
		s << " (complex)";
	else
		s << " (double)";
	// parameterizations
	if (a.pz.size() == 1) {
		s << " parameterized by " << *a.pz[0];
	} else if (a.pz.size() > 1) {
		s << ' ' << a.pz.size() << " parameterizations:\n";
		ts.str("");
		for (size_t i=0; i<a.pz.size(); i++) {
			ts << i+1 << ") " << *a.pz[i] << endl;
		}
		s << stringIndent(4, ts.str());
	}
#ifdef NEVER // need better way to display matrices, not with operator<<?
	// old operator<< just for small matrices
	const vector<double>& matval = a.const_data();
	if (a.is_complex()) {
		s << strArray(nr,nc,nr, (complex<double>*)&matval[0], a.mid());
	} else {
		s << strArray(nr,nc,nr, &matval[0], a.mid());
	}
#endif // NEVER : need better way to display matrices
	return s;
}


complex<Ad>
matvij (pset& plt, std::string const& desc, int row, int col) {
/*------------------------------------------------------------------
 * returns an element of a matrix given the description
 * ("mass", "stif", "gaf", etc).
 * Intended to be called by user-written C++ functions.
 *    row    1b row number
 *    col    1b column number
 *------------------------------------------------------------------*/
	T_(Trace trc(2,"matvij");)
	complex<Ad> rval;

	T_(trc.dprint("desc ",desc,", row ",row,", col ",col);)

	Matrix* mp = Matrix::find_desc(desc);
	if (mp == nullptr) {
		if (plt.monitoring()) {
			T_(trc.dprint("not throwing exception: monitoring");)
			return rval;
		} else {
			throw runtime_error(vastr("matrix \"", desc,
			"\" has not been defined; available matrices:\n", Matrix::inventory()));
		}
	}

	int nr = (int)mp->rsize();
	int nc = (int)mp->csize();

	if (row < 1 || row > nr) {   // 1b
		throw runtime_error(vastr("row ", row, "1b is out of the range 1 to ",
				nr, " (matvij)"));
	}
	if (col < 1 || col > nc) {
		throw runtime_error(vastr("column ", col, "1b is out of the range 1 to ",
			nc, " (matvij)"));
	}

	// evaluate the matrix...
	vector<complex<Ad>> work(nr*nc);
	mp->eval(plt, work);

	// ...then extract the return value
	rval = work[row-1 + nr*(col-1)];

	T_(trc.dprint("returning ",rval);)
	return rval;
}

//------------------------------------------------------------------
// Matrix::collection(): maintain a set of Matrix*
//------------------------------------------------------------------

vector<Matrix*> _collection;

vector<Matrix*>&
Matrix::
collection() { return _collection; }

Matrix*
Matrix::
find_desc(const string& desc) {
	for (auto mp : _collection)
		if (mp->_desc == desc)
			return mp;
	return nullptr;
}

Matrix*
Matrix::
find_mid(const string& mid) {
// Search for a matrix pointer in Matrix::collection(); if not
// there try to fetch it and add it to collection()
	for (auto mp : _collection)
		if (mp->_mid == mid)
			return mp;

	vector<Matrix*> mp = fetch(mid);
	if (!mp.empty()) {
		Matrix::insert(mp[0]);
		return mp[0];
	}
	return nullptr;
}

void
Matrix::
insert(Matrix* m) {
// Add a matrix pointer to Matrix::collection() if it is
// not already there
	for (auto mp : _collection)
		if (mp->mid() == m->mid())
			return;
	_collection.push_back(m);
}

string
matrix_list_summary(const vector<Matrix*>& list) {
// returns a string with a summary of each matrix in "list" 1/line
	ostringstream os;
	for (size_t i=0; i<list.size(); i++)
		os << i+1 << ") " << *list[i] << endl;
	return os.str();
}

// fetch/store a list of the matrix mids:
vector<pair<string,string> >
Matrix::
fetch_mids (string const& id) {
// Fetch a list of the mids of all matrices in Matrix::collection()
// Return pairs of (desc,mid), e.g. (mass,MHH)
	T_(Trace trc(1,"Matrix::fetch_mids");)
	vector<string> pairs;
	fio::fetch(vastr("mids,",id), pairs);
	vector<pair<string,string> > rval;
	// split these string into (desc,mid)
	for (auto& s : pairs) {
		vector<string> toks = string2tok(s, "=");
		if (toks.size() != 2)
			throw runtime_error(vastr("fetch_mids: bad pair <",s,">"));
		rval.push_back(make_pair(stripwhitespace(toks[0]),stripwhitespace(toks[1])));
		T_(trc.dprint("got ",toks[0]," = ",toks[1]);)
	}
	return rval;
}

void
Matrix::
store_mids (string const& id) {
// Store a list of the mids of all matrices in Matrix::collection()
	if (_collection.empty())
		return;

	vector<string> ml;
	for (auto mp : _collection) {
		ml.push_back(vastr(mp->desc(),"=",mp->mid()));
	}
	fio::store(vastr("mids,",id), ml);
}

#ifdef NEVER // Fma
void
Matrix::
vzmatrices(const string& vzid, Matrix*& nodes,
         Matrix*& coords, Matrix*& conn, Matrix*& gct) {
// Fetch the 4 matrices necessary for amviz, throws runtime_error
// if any of the 4 are missing
	T_(Trace trc(1,"vzmatrices");)

	nodes = coords = conn = gct = nullptr;
	// vzid may be empty - just fetch without extension
	string re;
	if (vzid.empty()) {
		re = "nodes|coords|conn|gct";
	} else {
		// extension must be _vzid: _ is the only char acceptable to
		// matlab, unix
		re = vastr(".*_",vzid);
	}
	vector<Matrix*> mats = Matrix::fetch(re);
	if (mats.empty()) {
		string exc{"no visualization matrices available"};
		if (!vzid.empty())
			exc += vastr(" with vzid \"",vzid,"\"");
		throw runtime_error(exc);
	}

	// the matrix name can be anything if it has a _vzid extension,
	// but the description must be "nodes", "coords", etc
	for (auto mp : mats) {
		if (mp->desc() == "nodes")
			nodes = mp;
		else if (mp->desc() == "coords")
			coords = mp;
		else if (mp->desc() == "conn")
			conn = mp;
		else if (mp->desc() == "gct")
			gct = mp;
		else
			flaps::warning("ignoring ",mp->mid(),", desc = ",mp->desc());
	}
	if (nodes == nullptr)
		throw runtime_error("nodes matrix missing from viz data");
	if (coords == nullptr)
		throw runtime_error("coords matrix missing from viz data");
	if (conn == nullptr)
		throw runtime_error("conn matrix missing from viz data");
	if (gct == nullptr)
		throw runtime_error("gct matrix missing from viz data");
}
#endif // NEVER // Fma

//------------------------------------------------------------------
// midcmp: compare 2 mids, 
static std::string midName(std::string const& mid);
static std::string midAttributes(std::string const& mid);

static bool
isEll (string const& str) {
// is str an ellipsis (...)?
	return (str == "...");
}

int
midcmp (string const& a, string const& b) {
/*------------------------------------------------------------------
 * Compare two matrix id's (mid). Either mid may be
 * a (POSIX extended) regular expression.
 * Returns:
 *   <0   if a < b
 *    0   if a == b
 *   >0   if a > b
 * so this function can be used with the STL sort algorithm
 * to sort mids.
 *
 * Attributes
 * ----------
 * mids may have attributes as in "gaf,rf=4.5,mach=0.4".
 * The attributes do not have to be in the same order.
 * Values associated with attributes may be strings, integers,
 * or floating-point numbers. When comparing floats they
 * are not compared as string, rather they are converted to
 * true floating point numbers and compared within machine
 * precision, so for example "gaf,rf=4.5,mach=0.4" is equal
 * to "gaf,mach=.40,rf=.45e+1"
 *
 * Regular-expressions
 * -------------------
 * gaf,.*        any number of attributes (at least one)
 * gaf,rf=*,... at least attribute "rf" with any value,
 *               possibly followed by any other attributes
 * gaf,...       any number of attributes (including none)
 * pk(1|2|3).*
 * stab,id=(pk1,pk2,pk3),...
 *------------------------------------------------------------------*/
	T_(Trace trc(3,"midcmp ",a," =? ",b);)
	size_t i, j;
	ostringstream os;
	int rval = 0;

	// entire strings match?
	if (a == b) {
		T_(trc.dprint("quick return: match");)
		return 0;
	}

	// strip quotes, split a and b into name and attributes
	string amid(a);
	stripquotes(amid);
	string bmid(b);
	stripquotes(bmid);
	string aname = midName(amid);
	string bname = midName(bmid);
	string aAtt = midAttributes(amid);
	string bAtt = midAttributes(bmid);

	// check the names...
	if (!aname.empty() && !bname.empty()) {
		// allow either name (but not both) to be a regular-expression
		bool ameta = aname.find_first_of("|*[]?") != string::npos;
		bool bmeta = bname.find_first_of("|*[]?") != string::npos;
		if (ameta) {
			if (bmeta) {
				T_(trc.dprint("both names have metachar");)
				rval = 0;
			}
			regex are(aname, regex_constants::icase);
			if (!regex_match(bname, are)) {
				T_(trc.dprint(bname," does not match re ",aname);)
				rval = 1;
			}
		} else if (bmeta) {
			regex bre(bname, regex_constants::icase);
			if (!regex_match(aname, bre)) {
				T_(trc.dprint(aname," does not match re ",bname);)
				rval = 1;
			}
		} else {
			// no regular-expressions - just a straight compare
			rval = aname.compare(bname);
		}
	}
	// return if the names don't match
	if (rval != 0 || (aAtt.empty() && bAtt.empty())) {
		T_(trc.dprint("returning ",rval," names differ or no attributes");)
		return rval;
	}

	// Names match - now we need to check the attributes
	// Parse each mid as a set of options
	// Watch out for options on any of the options,
	// as in "modes,set=1{rows=1}" which is a mistake
	// They do not have to have the same number of
	// attributes if one has ellipses (...)
	rval = aAtt.compare(bAtt);
	if (rval == 0) {
		T_(trc.dprint("returning 0: attributes match exactly");)
		return rval;
	}

	vector<Tok*> alist = flaps::lexer(aAtt);
	vector<Tok*> blist = flaps::lexer(bAtt);
	// check for curly braces in the rhs
	for (auto& ai : alist)
		if (!(ai->lopt.empty() && ai->ropt.empty())) {
			string exc = vastr("illegal matrix id: \"",a,
					"\" should not have curly braces");
			throw runtime_error(exc);
		}

	for (auto bi : blist)
		if (!(bi->lopt.empty() && bi->ropt.empty())) {
			string exc = vastr("illegal matrix id: \"",b,
					"\" should not have curly braces");
			throw runtime_error(exc);
		}

	bool aellipses = false;
	bool bellipses = false;
	for (auto& ai : alist)
		if (isEll(ai->lhs))
			aellipses = true;

	for (auto& bi : blist)
		if (isEll(bi->lhs))
			bellipses = true;

	// If both have ellipses we can return 0 (true) since any
	// combination of attributes will match...
	if (aellipses && bellipses) {
		T_(trc.dprint("returning 0: both have ellipses");)
		return 0;
	}
	// ... on the other hand if neither have ellipses the number
	// of attributes must match
	rval = 0;
	if (!aellipses && !bellipses) {
		if (!aname.empty() && !bname.empty()) {
			if (alist.size() != blist.size()) {
				T_(trc.dprint("returning ",rval,": different number of attributes");)
				return (alist.size() > blist.size() ? 1 : -1);
			}
		}
	}

	// If b has ellipses, switch alist and blist so that alist
	// is always the only one that can have ellipses to make the code
	// simpler
	if (bellipses) {
		vector<Tok*> tmp = alist;
		alist = blist;
		blist = tmp;
		aellipses = true;
		bellipses = false;
	}
	// for each attribute in A check against each in B; if
	// one is found that matches break out of the inner loop.
	// If the inner loops finishes => not found
	// Watch for metacharacters: att=* or att=.* means match any
	rval = 0;
	for (i=0; i<alist.size(); i++) {
		if (isEll(alist[i]->lhs))
			break;
		string arhs;
		if (!alist[i]->svec.empty())
			arhs = alist[i]->svec[0];
		for (j=0; j<blist.size(); j++) {
			string brhs;
			if (!blist[i]->svec.empty())
				brhs = blist[i]->svec[0];
			T_(trc.dprintn(alist[i]," =? ",blist[j]);)
			if (alist[i]->lhs == blist[j]->lhs) {
				if (compare(arhs, brhs)) {
					// found a match - quit this loop
					T_(trc.dprint(" yes");)
					break;
				}
			}
		}
		// did the inner loop finish (ie not found)?
		// but it is ok if stage was in a but not b
		if (j == blist.size()) {
			T_(trc.dprint(alist[i]->lhs," is an attribute in a but not b");)
			rval = 1;
			break;
		}
	}

	T_(trc.dprint("returning ",rval);)
	return rval;
}

string
midName(string const& mid) {
// returns the "name" part of a mid i.e. without attributes
// The mid might have no name, only attributes
	string rval;
	if (mid.empty())
		return mid;
	if (isEll(mid))
		return rval;
	string::size_type idx = mid.find(',');
	if (idx == string::npos) {
		// no commas but watch out for attribute only, e.g. mach=0
		idx = mid.find('=');
		if (idx == string::npos) {
			return mid;
		} else {
			return rval;
		}
	} else {
		rval = mid.substr(0,idx);
		// check that this is not an attribute
		idx = rval.find('=');
		if (idx != string::npos)
			rval.clear();
	}
	return rval;
}

string
midAttributes(string const& mid) {
// returns the "attributes" part of a mid
// Ok if "mid" has no name - only attributes
	string rval;
	if (mid.empty())
		return rval;
	string::size_type idx = mid.find(',');
	if (idx != string::npos) {
		rval = mid.substr(idx+1);
		string name = mid.substr(0,idx);
		idx = name.find('=');
		if (idx != string::npos)
			rval = mid;
	} else {
		// no comma: only attributes if it has = or is ellipsis
		idx = mid.find('=');
		if (idx != string::npos || isEll(mid))
			rval = mid;
	}
	return rval;
}

// printing function for Matrix
namespace flaps {
// Specializations for flaps::summarize()
std::string
summarize (size_t n, int const* x, unsigned int maxchar) {
// specialization for vectors of integers
	std::ostringstream os;
	for (size_t i=0; i<n; i++)
		os << " " << x[i];
	return os.str();
}

std::string
summarize (std::vector<int> const& x, unsigned int maxchar) {
	return flaps::summarize(x.size(), &x[0], maxchar);
}
}  // namespace flaps

// strArray(... double) is a template in matrix.h
string
strArray(int nr, int nc, int lda, const complex<double>* a) {
// write a complex matrix to a string
	ostringstream os;
	Form sci(6);
	sci.width(13);
	for (int i=0; i<nr; i++) {
		for (int j=0; j<nc; j++) {
			os << '(' << sci(a[IJ(i,j,lda)].real()) << ", "
				<< sci(a[IJ(i,j,lda)].imag()) << ") ";
		}
		os << endl;
	}
	return os.str();
}

#ifdef NEVER // needs work: how to return path(s), input ostream, path
// XXX modify matview to groke Ad?
// various specializations of mprint
void
mprint (ostream* s, int nr,int nc,int lda,const std::vector<Ad>& a,const std::string& title) {
// vector<Ad> print: print the values and derivatives
	std::vector<double> av(nr*nc);
	std::vector<std::string> toprint{""};
	std::vector<std::string> const& adparnames = Ad::adnames();
	for (auto& nm : adparnames)
		toprint.push_back(nm);
	// "values" array
	extract(a, "", av.data());
	// check the width of nc elements
	std::string width = strArray(nc,1,1,av.data());
	for (auto& nm : adparnames) {
		extract(a, nm, &av[0]);
		std::string width = strArray(nc,1,1,&av[0]);
		std::string ti{vastr(title," partial wrt ",nm)};
		if (nm.empty())
			ti = title;
		if ((int)width.size() > page_width()) {
			std::string path = vastr(getftmp(),'/',ti,".mm");
			MM::exporter(path, ti, &av[0], nr, nc);
		} else {
			std::ostream* lf = Trace::logfile();
			*lf << ti << std::endl << strArray (nr, nc, lda, &av[0]);
		}
	}
}

void
mprint (int nr,int nc,int lda,const std::vector<std::complex<Ad>>& a,const std::string& title) {
// vector<complex<Ad>> print: print the values and derivatives
	std::vector<std::complex<double> > av(nr*nc);
	std::vector<std::string> toprint{""};
	std::vector<std::string> const& adparnames = Ad::adnames();
	for (auto& nm : adparnames)
		toprint.push_back(nm);
	for (auto& nm : toprint) {
		extract(a, nm, &av[0]);
		std::string width = strArray(nc,1,1,&av[0]);
		std::string ti{vastr(title," partial wrt ",nm)};
		if (nm.empty())
			ti = title;
		if ((int)width.size() > page_width()) {
			std::string path = vastr(getftmp(),'/',ti,".mm");
			MM::exporter(path, ti, &av[0], nr, nc);
		} else {
			std::ostream* lf = Trace::logfile();
			*lf << ti << std::endl << strArray (nr, nc, lda, &av[0]);
		}
	}
}

void
mprint (int nr, int nc, int lda, const int* a, const std::string& title) {
	std::ostream* lf = Trace::logfile();
	*lf << title << std::endl << strArray (nr, nc, lda, a);
}

void
mprint (int nr, int nc, int lda, const double* a, const std::string& title) {
	std::string width = strArray(nc,1,1,a);
	if ((int)width.size() > page_width()) {
		std::string path = vastr(getftmp(),'/',title,".mm");
		MM::exporter(path, title, a, nr, nc);
	} else {
	 	std::ostream* lf = Trace::logfile();
	 	*lf << title << std::endl << strArray (nr, nc, lda, a);
	}
}

void
mprint (int nr, int nc, int lda, const std::complex<double>* a, const std::string& title) {
	std::string width = strArray(nc,1,1,a);
	if ((int)width.size() > page_width()) {
		std::string path = vastr(getftmp(),'/',title,".mm");
		MM::exporter(path, title, a, nr, nc);
	} else {
		std::ostream* lf = Trace::logfile();
		*lf << title << std::endl << strArray (nr, nc, lda, a);
	}
}
#endif // NEVER // needs work: how to return path(s), input ostream, path


#ifdef MAIN
#undef MAIN

#include "fio.h"  // Ftmpdir

static void
usage(char const* prog) {
	T("Usage: %s regex mid:               test midCmp\n",prog);
	T("       %s -o string string:        test stringCmpIC\n",prog);
	exit(1);
}
int
main(int argc, char **argv) {
// this main tests midcmp


   if (argc < 3)
		usage(argv[0]);

	if (argv[1][0] == '-') {
		// midCmpIC
		if (argv[1][1] == 'o') {
			if (argc != 5)
				usage(argv[0]);
			cout << argv[2] << argv[4] << argv[3]
				<< " -> " << midCmpIC(argv[2], argv[3]) << endl;
			exit(0);
		} else {
			usage(argv[0]);
		}
	}

	cout << "does \"" << argv[1] << "\" == \"" << argv[2] << "\"?... ";
	try {
		int cmp = midcmp(argv[1], argv[2]);
		if (cmp == 0) {
			cout << "yes\n";
		} else {
			if (cmp < 0)
				cout << "no (" << cmp << "): \"" << argv[1] << "\" < \"" << argv[2] << "\"\n";
			else
				cout << "no (" << cmp << "): \"" << argv[1] << "\" > \"" << argv[2] << endl;
		}
		cout << "midName(" << argv[1] << ") = " << midName(argv[1]) << endl;
		cout << "midAttributes(" << argv[1] << ") = " << midAttributes(argv[1]) << endl;
		cout << "midName(" << argv[2] << ") = " << midName(argv[2]) << endl;
		cout << "midAttributes(" << argv[2] << ") = " << midAttributes(argv[2]) << endl;
	} catch (runtime_error& s) {
		cerr << "Error: " << s.what() << std::endl;
	}
}


Matrix*
makeMatrix(int n, double val, string const& mid) {

	bool iscmplx{false};
	Matrix* rval = new Matrix (mid, "makeMatrix", n, n, iscmplx);
	double* rp = rval->elem();
	for (int i=0; i<n; i++) {
		for (int j=0; j<n; j++) {
			rp[IJ(i,j,n)] = 0.0;
		}
		rp[IJ(i,i,n)] = val;
	}
	return rval;
}


int
main (int argc, char **argv) {
	Ftmpdir ftmp;
	int n{3};
	string mid{"test"};
	putEnv("KEEPFTMP=1");

	Matrix* rp = makeMatrix(n, 2.0, mid);
	cout << "Created:\n" << strArray(n,n,n,rp->elem()) << endl;
	rp->store();
	Matrix* newmat = Matrix::fetch(mid);
	cout << "Fetched:\n" << strArray(n,n,n,newmat->elem()) << endl;

	try {
		ostringstream os;
		Matrix *ap, *bp;
		vector<Matrix*> matlist;
		Par* vtasp = gpset::find("vtas");
		vector<string> params;
		int i;
		/*
		 * Create some matrices for interpolation wrt vtas...
		 */
		params.push_back("vtas");
		vector<double> parvals(1);
		for (i=0; i<4; i++) {
			os.str("");
			os << "mat" << i;
			string mid{os.str()};
			ap = makeMatrix(n, (double)i,mid);
			// add a NoOpPz
			// parvals[0] = 100.0*i;
			// ap->pz.push_back(new NoOpPz(params, parvals, n, n, ap->data()));
			// cast it to Complex...
			// ComplexMatrix* cp = castReal2Complex(ap);
			matlist.push_back(ap);
		}
		/*
		 * ... then interpolate
		 */
		IntPz* pz = new IntPz(matlist, 0, 0.0);
		/*
		 * ... and associate it with a new matrix
		 * XXX How do we ensure that the parameters in pz
		 * correspond to the gm?
		 */
		mid = "interp";
		bp = makeMatrix(n, -1.0, mid);
		bp->pz.push_back(pz);

		// eval at 150
		vtasp->value(150.0);
		vector<complex<Ad>> result(n,n);
		bp->eval(gpset::get(), result, n, n);
		cout << "evaluated at 150:\n" << strArray(n,n,n,result) << endl;

		// store, fetch it
		bp->store();
		Matrix* newmat = Matrix::fetch(mid);
		cout << "Fetched:\n" << strArray(n,n,n,newmat->elem()) << endl;
		// add it to the collection
		Matrix::insert(newmat);
		Matrix* cp = Matrix::find_mid(mid);
		if (cp == nullptr) {
			string exc = vastr("find_mid(",mid,") failed");
			throw runtime_error(exc);
		}
		// eval at 250
		vtasp->value(250.0);
		cp->eval(gpset::get(), result, n, n);
		cout << "evaluated at 250:\n" << strArray(n,n,n,result) << endl;

		cout << "Inventory:\n" << Matrix::inventory();

	} catch (runtime_error& s) {
		error("Error: ",s);
	}
}
 

#endif // MAIN
