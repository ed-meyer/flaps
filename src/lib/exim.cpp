//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

/*------------------------------------------------------------------
 * Import/export files of various formats
 *   Matlab mat file: namespace Matlab (.mat)
 *   Nastran output4: namespace Output4 (.op4)
 *   Universal file:  namespace UF (.uf)
 *   Matrix Market:   namespace MM (.mm)
 *   xfig files:      namespace Fig (.fig)
 *   ASCII plot:      namespace Apf (.apf)
 *   LTI:             namespace LTI
 *------------------------------------------------------------------*/


#include <cmath>
#include <cstdio>
#include <cstring>
#include <iomanip>
#include <unistd.h> // access
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include "config.h"
#include "exim.h"
#include "fio.h"
#include "matrix.h"
#include "text.h"

using namespace std;

// private utilities for reading text files:
// void read_line(ifstream& file, string& s, string comment=""):
//         read the next line into "s"
// int line_number(): returns the current line number

class EOF_exc : public std::runtime_error {
public:
	EOF_exc(const std::string& msg) : runtime_error(msg) {}
};

string
get_path(const string& name, const string& ext) {
// Check name for read access on the working directory (FWD) and
// on the temporary directory (FTMP), possibly with the extension (.ext)
	Trace trc(2,"get_path");

	// if "name" is already a full path just check it
	if (access(name.c_str(), R_OK) != -1)
		return name;

	string cwd = get_cwd();
	// a lambda for checking access on a specified directory
	auto check_access = [&ext](const string& name, const string& wd) {
		string path = wd + "/" + name;
		if (!ext.empty() && rsubstr(path, ext.size()) != ext) {
			// try with .ext first to avoid taking "name" when there is a "name.ext"
			string expath = path + ext;
			if (access(expath.c_str(), R_OK) != -1)
				return expath;
		}

		if (access(path.c_str(), R_OK) != -1)
			return path;
		return string("");
	};
	string rval = check_access(name, getftmp());
	if (rval.empty())
		rval = check_access(name, getfwd());
	return rval;
}

namespace Output4 {
bool header(Ascii& in, string& mid, int& nr, int& nc,
		int& dtype, int& wordsPerElement, int& formatwidth);
Matrix* data (Ascii& in, string& mid, int nr, int nc,
				int dtype, int wordsPerElement, int formatwidth);
} // namespace Output4

// export functions: XXX should go in namespace Output4
static void
writeComplexFull (ostream& file, size_t nr, size_t nc, complex<double>* data);
static void
writeRealFull (ostream& file, size_t nr, size_t nc, double* data);

static
vector<string>
tokenize(const string& line, size_t width) {
// split a string into tokens, each of length "width"
	Trace trc(1,"tokenize");
	size_t len{line.size()};
	size_t start{0};
	vector<string> rval;
	trc.dprint("line ",line,", len ",len,", width ",width);
	while (start < len) {
		rval.push_back(line.substr(start, width));
		start += width;
	}
	trc.dprint("returning ",rval.size()," tokens");
	return rval;
}

vector<double>
floatTok(string const& line) {
	vector<double> rval;
// figure out how wide each float field is;
// there may be no space between floats, as in
// 2.235581580E-04-3.141353039E+02-4.173338717E-04 6.362884978E+02-1.926432754E-04
	string::size_type start = 0;
	string::size_type end = 0;
	double x;
	while(end != string::npos) {
		if (!str2double(line, x, start, end)) {
			throw runtime_error(vastr("expecting real number, got ", line));
		}
		rval.push_back(x);
		start = end;
	}
	return rval;
}

static vector<Matrix*>
splitQ(Matrix* ap) {
// Nastran stores QHH matrices as (n,nk*n) where each block of n
// columns represents one gaf matrix at a particular k-value
	Trace trc(1,"splitQ ",ap->summary());
	size_t nr = ap->rsize();
	size_t nc = ap->csize();
	vector<Matrix*> rval;
	vector<double>& aval = ap->data();
	complex<double>* matval = (complex<double>*)&aval[0];
	ostringstream os;

	if (nc%nr != 0) {
		flaps::warning("the number of columns in ", ap->mid(), " (",
			nc, ") is not a multiple of the number of rows (", nr, ')');
		return rval;
	}
	size_t nq = nc/nr;
	// quick return if only one matrix
	if (nq == 1) {
		rval.push_back(ap);
		trc.dprint("quick return: only one matrix");
		return rval;
	}

	// treat the input Q as an (nsq2,nq) real matrix
	for (size_t j=0; j<nq; j++) {
		os.str("");
		os << ap->mid() << j+1;
		Matrix* newq = new Matrix(os.str(), "Gaf", nr, nr, true);
		blas_copy (nr*nr, &matval[IJ(0,j,nr*nr)], 1,
				(complex<double>*)&newq->data()[0], 1);
		trc.dprint("new Matrix <", newq->summary());
		rval.push_back(newq);
	}
	return rval;
}

vector<Matrix*>
Output4::
importer(const string& path, const string& output) {
//------------------------------------------------------------------
// Read a set of matrices from a NASTRAN output4 file
// Aero matrices come as one big matrix stacked columnwise -
// it will be split into square matrices
//------------------------------------------------------------------
	Trace trc(1,"Output4:importer");
	vector<Matrix*> rval;
	int wordsPerElement;
	int formatwidth;
	int nr, nc;
	int dtype;
	string mid;
	Ascii in(path);

	try {
		while(true) {
			if (!Output4::header(in, mid, nr, nc, dtype, wordsPerElement, formatwidth))
				break;

			Matrix* ap = Output4::data (in, mid, nr, nc, dtype,
					wordsPerElement, formatwidth);

			 // look for a matrix beginning with 'Q' - break
			 // it into square matrices
			if (ap->mid()[0] == 'Q') {
				vector<Matrix*> aml = splitQ(ap);
				for (auto ap : aml)
					rval.push_back(ap);
			} else {
				rval.push_back(ap);
				trc.dprintm(nr,nc,nr,ap->elem(),vastr("read:",ap->mid()));
			}
		}
	// only catch EOF, let others go upstream
	} catch (EOF_exc& s) {
		trc.dprint("caught EOF exception\n");
	}
	return rval;
}

bool
Output4::
header (Ascii& in, string& mid, int& nr, int& nc,
		int& dtype, int& wordsPerElement, int& formatwidth) {
//------------------------------------------------------------------
// Read the first record for a NASTRAN matrix:
// ncol nr form type dmapname
//------------------------------------------------------------------
	Trace trc(1,"header");
	string formatbuf;
	int form;
	string line;

	// Ascii returns false on EOF
	if (!in(line)) {
		trc.dprint("returning false: EOF");
		return false;
	}

	// split the line into 8-char strings...
	vector<string> toks = tokenize(line, 8);
	if (toks.size() < 6)
		throw runtime_error(vastr("line ",in.line_number(),
				": expecting > 5 items, got \"", line,"\""));
	// ... then convert
	nc = stoi(toks[0]);
	nr = stoi(toks[1]);
	form = stoi(toks[2]);
	dtype = stoi(toks[3]);
	mid = stripwhitespace(toks[4]);
	formatbuf = toks[5];
	if (toks.size() > 6)
		formatbuf += toks[6];
	trc.dprint("got nc ",nc," nr ",nr," form ",form, " dtype ",dtype," mid ",mid," format ",formatbuf);

	// Negative nc means BIGMAT=TRUE; we don't worry about
	// handling this here - just change signs on nc
	if (nr < 0)
		nr = -nr;

	if (nr <= 0 || nc <= 0) {
		string exc{vastr("bad matrix header on line ",in.line_number(),": nr = ",
			nr, ", nc = ", nc)};
		trc.dprint("throwing exception ",exc);
		throw runtime_error(exc);
	}

	// Check header line for a format (e.g. 1P,3E22.15)
	// Only get the width of float fields from this format
	// XXX unused - floatTok() does not need it
	regex re{R"(E([^.]*)\.)", regex_constants::icase};
	string format(formatbuf);
	formatwidth = 0;
	smatch mch;
	if (regex_search(format, mch, re) && mch.size() > 1) {
		str2int(mch[1], formatwidth);
		trc.dprint("got format width = ",formatwidth);
	} else {
		string exc{vastr("no format found on header line ", in.line_number(),
			": \"", format, "\"")};
		trc.dprint("throwing exception ",exc);
		throw runtime_error(exc);
	}


	// Note: the NASTRAN manual says that double-precision types
	// have two words/element; this doesn't make sense for ASCII files
	// but wordsPerElement is used in computing the number of elements
	// in a column as (L-1)/wordsPerElement where L is the
	// length of the string (number of tokens).
	switch (dtype) {
		case 1:						/* single precision */
			wordsPerElement = 1;
			break;
		case 2:						/* double-precision */
			wordsPerElement = 2;
			break;
		case 3:                // single precision complex
			wordsPerElement = 2;
			break;
		case 4:
			wordsPerElement = 4;       // double precision complex
			break;
		default:
			string exc{vastr("unrecognized datatype on line ",
				in.line_number(), " (", dtype, ')')};
			trc.dprint("throwing exception: ",exc);
			throw runtime_error(exc);
	}
	return true;
}

static void
toFull (size_t nr, size_t nc, int* ri, int* cs, double* nz, double* gm) {
// convert a sparse double matrix to a full double matrix
	size_t j;
	for (j=0; j<nr*nc; j++)
			gm[j] = 0.0;
	for (j=0; j<nc; j++) {
		int first = cs[j];
		int last = cs[j+1];
		for (int k=first-1; k<last-1; k++) {
			int i = ri[k] - 1;
			gm[IJ(i,j,nr)] = nz[k];
		}
	}
}

static void
toFull (size_t nr, size_t nc, int* ri, int* cs,
		complex<double>* nz, complex<double>* gm) {
// convert a sparse complex<double> matrix to a full complex<double> matrix
	size_t j;
	for (j=0; j<nr*nc; j++)
			gm[j] = complex<double>(0.0);
	for (j=0; j<nc; j++) {
		int first = cs[j];
		int last = cs[j+1];
		for (int k=first-1; k<last-1; k++) {
			int i = ri[k] - 1;
			gm[IJ(i,j,nr)] = nz[k];
		}
	}
}


Matrix*
Output4::
data (Ascii& in, string& mid, int nr, int nc,
		int dtype, int wordsPerElement, int formatwidth) {
/*------------------------------------------------------------------
 * Reads the data portion of a NASTRAN matrix in an output4 file
 *------------------------------------------------------------------*/
	Trace trc(1,"Output4::data");
	Matrix *rval;
	char *endp;
	int colno, irow, nw, nws, nwr, nel;
	int i;
	int nnz = 0;
	int nextcolstart;
	vector<int> ri;
	int nRealCol;
	double* fp{nullptr};
	vector<double> nz;
	string line;
	string desc{"Nastran"};
	vector<string> quotes;
	bool iscomplex = (dtype > 2);

	trc.dprint("mid ",mid," nr ",nr," nc ",nc," words/element ",wordsPerElement);

	// cs: column-start indices into ri, nz
	vector<int> cs(nc+1, 1);

	// Allocate space to hold one column: data and row indices
	if (iscomplex) {
		nRealCol = 2*nr;
	} else {
		nRealCol = nr;
	}
	vector<double> column(nRealCol, 0.0);
	vector<int> columnri(nr, 0);

	/*
	 * Read one column at a time into "column", append it to nz array,
	 * and append the corresponding row indices in columnri to ri
	 */
	vector<string> toks;
	try {			// XXX Ascii won't throw on EOF
		while (1) {
			fill(column.begin(), column.end(), 0.0);
			if (!in(line))
				throw runtime_error(vastr("premature EOF line ",in.line_number()));
			toks = tokenize(line, 8);
			assert(toks.size() > 2);
			colno = stoi(toks[0]);
			irow = stoi(toks[1]);
			nw = stoi(toks[2]);
			trc.dprint("got colno ",colno," irow ",irow," nw ",nw);
			if (irow-1+nw > nRealCol) {
				string exc{vastr("column ", colno, " of ", mid, " has ",
					nw, " elements, but the matrix was declared to have ",
					nr, " rows")};
				trc.dprint("throwing exception: ",exc);
				throw runtime_error(exc);
			}
			nel = nw;
			if (iscomplex) nel /= 2;
			// colno == nc+1 => last "dummy" column
			if (colno > nc) {
				 // read until we have read nw floats
				while(nw > 0) {
					if (!in(line))
						throw runtime_error(vastr("premature EOF line ",in.line_number()));
					toks = tokenize(line, formatwidth);
					nw -= toks.size();
				}
				break;
			}
			if (colno <= 0) {
				string exc{vastr("bad column number on line ",in.line_number(),
						": ", colno)};
				trc.dprint("throwing exception: ",exc);
				throw runtime_error(exc);
			}
			if (irow < 0 || irow > nr) {
				string exc{vastr("bad row number on line ", in.line_number(),
						": ", irow, " (matrix has ", nr, " rows)")};
				trc.dprint("throwing exception: ",exc);
				throw runtime_error(exc);
			}
			// Read the floats: non-sparse. Question: is it possible
			// to have mixed sparse and non-sparse columns in one file?
			if (irow > 0) {
				if (iscomplex)
					fp = &column[2*(irow-1)];
				else
					fp = &column[irow-1];
				while (nw > 0) {
					if (!in(line))
						throw runtime_error(vastr("premature EOF line ",in.line_number()));
					toks = tokenize(line, formatwidth);
					for (auto& t : toks) {
						*fp++ = stod(t);
						nw--;
					}
				}
				if (iscomplex) {
					for (i=0; i<nel; i++) {
						int cindx = (irow-1)*2 + 2*i;
						nz.push_back(column[cindx]);
						nz.push_back(column[cindx + 1]);
					}
				} else {
					for (i=0; i<nel; i++) {
						nz.push_back(column[irow-1+i]);
					}
				}
				for (i=0;  i<nel; i++)
					ri.push_back(irow+i);
				nnz = ri.size();

				nextcolstart = ri.size() + 1;	// one past current number of nz
				for (i=colno; i<=nc; i++)  // set next and all subsequent cs
					cs[i] = nextcolstart;
				trc.dprint("nnz ",nnz,", cs[",colno,"] = ",cs[colno]);

			// Read the data values for sparse: nw *words* comprising
			// "strings" of words, which if BIGMAT=false are
			// string header (IS)    (one word)
			// (L-1)/wordsPerElement elements
			} else {
				int strhdr{0};
				int ncon{65536};
				for (nwr = 0; nwr < nw; ) {
					if (!in(line))
						throw runtime_error(vastr("premature EOF line ",in.line_number()));
					// Determine if this is BIGMAT or not: BIGMAT format has
					// two integers on this line, non-BIGMAT has the two packed
					// into one
					vector<string> t = string2tok(line, " ", quotes);
					if (t.size() < 2) {
						strhdr = strtol(t[0].c_str(), &endp, 10);
						nws = strhdr/ncon;
						irow = strhdr - nws*ncon;
					} else {
						nwr++;		// for the extra integer
						nws = strtol(t[0].c_str(), &endp, 10);
						irow = strtol(t[1].c_str(), &endp, 10);
					}

					nwr += nws;
					// nel is the actual number of elements in this string - it is
					// the number of *words* (nws-1) divided by the number of
					// words per element
					nel = (nws - 1)/wordsPerElement;
					if (nel < 0 || nel > nr) {
						string exc{vastr("bad string header (",line,") on line ",
							in.line_number(), ": nw=", nel)};
						trc.dprint("throwing exception: ",exc);
						throw runtime_error(exc);
					}
					trc.dprint("strhdr ",strhdr," nws ",nws," nel ",nel," irow ",irow);
					if (irow < 0 || irow > nr) {
						string exc{vastr("bad row number on line ", in.line_number(),
							": ", irow, " (should be < ", nr,')')};
						trc.dprint("throwing exception: ",exc);
						throw runtime_error(exc);
					}

					for (i=0; i<nel; ) {
						if (!in(line))
							throw runtime_error(vastr("premature EOF line ",in.line_number()));
						toks = tokenize(line, formatwidth);
						for (auto& t : toks) {
							column[i] = stod(t);
							columnri[i] = irow+i;
							i++;
						}
					}
					if (nnz+nel > 0) {
						if (iscomplex) {
							for (i=0; i<nel; i++) {
								nz.push_back(column[2*i]);
								nz.push_back(column[2*i+1]);
							}
						} else {
							for (i=0; i<nel; i++) {
								nz.push_back(column[i]);
							}
						}
						for (i=0;  i<nel; i++)
							ri.push_back(columnri[i]);
						nnz += nel;
						for (i=colno; i<=nc; i++)
							cs[i] = nnz+1;	// start of next column
					}
				}
			}
		}
	} catch (runtime_error& s) {
		string exc{vastr("premature end-of-file reading \"", mid, "\"")};
		trc.dprint("throwing exception: ",exc);
		throw runtime_error(exc);
	}

	trc.dprint("row indices",ri);
	trc.dprint("col starts",cs);
	trc.dprint("non-zeros",nz);
	/*
	 * Check that the lower-right element is present
	 */
	nnz = ri.size();
	if (ri[nnz-1] != nr) {
		nnz++;
		cs[nc] = nnz+1;
		ri.push_back(nr);
		if (iscomplex) {
			nz.push_back(0.0);
			nz.push_back(0.0);
		} else {
			nz.push_back(0.0);
		}
	}

	// expand the sparse rep to full
	if (iscomplex) {
		Matrix* mp = new Matrix(mid, desc, nr, nc, iscomplex);
		toFull (nr, nc, &ri[0], &cs[0], (complex<double>*)&nz[0],
				(complex<double>*)&mp->data()[0]);
		rval = mp;
	} else {
		Matrix* mp = new Matrix(mid, desc, nr, nc, iscomplex);
		toFull (nr, nc, &ri[0], &cs[0], &nz[0], &mp->data()[0]);
		rval = mp;
	}

	trc.dprint("returning ",*rval);
	return rval;
}

bool
Output4::
exporter (const string& path, const vector<Matrix*>& ml, bool append) {
//------------------------------------------------------------------
// Write a set of matrices to an output4-formatted file
//------------------------------------------------------------------
	Trace trc(1,"Output4::exporter");
	bool rval = true;
	string format;
	string exc;

	trc.dprint("matrix: ",ml[0]->summary()," to ",path);

	// open the file
	std::ios_base::openmode om{std::ios::trunc};
	if (append)
		om = std::ios::app;
	ofstream file(path, om);

	for (auto mp : ml) {
		size_t nr = mp->rsize();
		size_t nc = mp->csize();
		int form{1};  // square
		if (nr != nc)
			form = 2;       // rectangular

		int type{1};
		if (mp->is_complex())
			type = 3;

		// header
		file << setw(8) << nc << setw(8) << nr << setw(8)
			<< form << setw(8) << type << setw(8) << mp->mid()
			<< " 1P,4E25.16" << endl;
		// fprintf (fp, "%8d%8d%8d%8d%-8s 1P,4E%s\n",
		// 		(int)ml[im]->csize(), (int)ml[im]->rsize(),
		// 		form, type, ml[im]->mid().c_str(), format.c_str());

		if (!mp->is_complex()) {
			writeRealFull (file, nr, nc, mp->elem());
		} else {
			writeComplexFull (file, nr, nc, mp->celem());
		}
	}
	return rval;
}

static void
writeRealFull (ostream& file, size_t nr, size_t nc, double* data) {
	Trace trc(1,"writeRealFull");
	size_t first, last;
	size_t floatsPerLine = 4;
	Form sci(16);
	sci.width(25);

	for (size_t j = 0; j < nc; j++) {
		// find the first non-zero in this column...
		for (first=0; first<nr; first++) {
			if (data[IJ(first,j,nr)] != 0.0)
				break;
		}
		if (first == nr)
			continue;
		// ... and the last (plus one)
		for (last=nr; last>first; last--) {
			if (data[IJ(last-1,j,nr)] != 0.0)
				break;
		}
		if (last == first)
			continue;
		file << setw(8) << j+1 << setw(8) << first+1 << setw(8) << last-first << endl;
		// fprintf (fp, "%8d%8d%8d\n", (int)(j+1), (int)(first+1), (int)(last-first));
		for (size_t i=first; i<last; i++) {
			file << sci(data[IJ(i,j,nr)]);
			// fprintf (fp, fmt, data[IJ(i,j,nr)]);
			if ((i-first+1)%floatsPerLine == 0)
				file << endl;
				// fputc ('\n', fp);
		}
		if ((last-first)%floatsPerLine != 0)
			file << endl;
			// fputc ('\n', fp);
	}
	// sentinal column:
	file << setw(8) << nc+1 << setw(8) << 1 << setw(8) << 1 << endl;
	// fprintf (fp, "%8d%8d%8d\n", (int)(nc+1), 1, 1);
	file << sci(1.0) << endl;
	// fprintf (fp, fmt, 1.0);
	// fputc ('\n', fp);
}

static void
writeComplexFull (ostream& file, size_t nr, size_t nc, complex<double>* data) {
	Trace trc(1,"writeComplexFull");
	size_t first, last;
	Form sci(16);
	size_t elemPerLine = 2;

	sci.width(25);

	for (size_t j = 0; j < nc; j++) {
		// find the first non-zero in this column...
		for (first=0; first<nr; first++) {
			if (std::abs(data[IJ(first,j,nr)]) != 0.0)
				break;
		}
		if (first == nr)
			continue;
		// ... and the last (plus one)
		for (last=nr; last>first; last--) {
			if (std::abs(data[IJ(last-1,j,nr)]) != 0.0)
				break;
		}
		if (last == first)
			continue;
		file << setw(8) << j+1 << setw(8) << first+1 << setw(8) << 2*(last-first) << endl;
		// fprintf (fp, "%8d%8d%8d\n", (int)(j+1), (int)(first+1), (int)(last-first));
		for (size_t i=first; i<last; i++) {
			file << sci(data[IJ(i,j,nr)].real()) << sci(data[IJ(i,j,nr)].imag());
			// fprintf (fp, fmt, data[IJ(i,j,nr)]);
			if ((i-first+1)%elemPerLine == 0)
				file << endl;
				// fputc ('\n', fp);
		}
		if ((last-first)%elemPerLine != 0)
			file << endl;
			// fputc ('\n', fp);
	}
	// sentinal column:
	file << setw(8) << nc+1 << setw(8) << 1 << setw(8) << 1 << endl;
	// fprintf (fp, "%8d%8d%8d\n", (int)(nc+1), 1, 1);
	file << sci(1.0) << sci(0.0) << endl;
	// fprintf (fp, fmt, 1.0);
	// fputc ('\n', fp);
}


namespace UF {
vector<int> connectivity(Ascii& in);
vector<complex<double>> disp(Ascii& in, vector<int>& gct_node_numbers);
Matrix* coord(Ascii& in, vector<int>& node_numbers); // XXX deprecated
vector<double> coordinates(Ascii& in, vector<int>& node_numbers);
void ignor(Ascii& in);

vector<int>
penlift2segments(const vector<int>& penlift) {
// convert connectivities in penlift format to segment format:
// pairs of ints defining 2 nodes where a line is to be drawn
// XXX deprecated: use penlift2segidx()
	Trace trc(2,"penlift2segments");
	vector<int> rval;
	int prev{0};
   for (size_t i=0; i<penlift.size(); i++) {
      if (penlift[i] == 0) {
         prev = 0;
         continue;
      }
      if (prev > 0) {
         rval.push_back(prev);
         rval.push_back(penlift[i]);
         prev = penlift[i];
      } else if (penlift[i+1] > 0) {
         rval.push_back(penlift[i++]);
         rval.push_back(penlift[i]);
			prev = penlift[i];
		}
	}
	trc.dprint("connectivities in segment format",rval);
	return rval;
}

vector<int>
penlift2segidx(const vector<int>& penlift, const vector<int>& nodes) {
// convert connectivity in "penlift" format (tracelines from a Universal
// file) to "segment index" (segidx) format: 0-based indexes into the
// node number array, so the node numbers for (0b) ordinal segment j are
//   nodes[segidx[2*j]] and nodes[segidx[2*j+1]]
// and the coordinates of segment j are
//   coords[3*segidx[2*j]],coords[3*segidx[2*j]+1],coords[3*segidx[2*j]+2]
// and
//   coords[3*segidx[2*j+1]],coords[3*segidx[2*j+1]+1],coords[3*segidx[2*j+1]+2]
// seems complicated but a more efficient way of visualizing with OpenGL
	Trace trc(2,"UF::penlift2segidx");

	// first convert penlift to segments: pairs of node numbers
	vector<int> segments = UF::penlift2segments(penlift);

	// then convert the node numbers to indexes into the nodes array
	vector<int> segidx;
	for (auto ci : segments) {
      auto idx = find(nodes.begin(), nodes.end(), ci);
      if (idx == nodes.end())
         throw runtime_error(vastr("node ",ci," in the connectivity"
            " has not been defined"));
      segidx.push_back((int)(idx - nodes.begin()));	// zero-based index
   }
   trc.dprint("segment indices (segidx)",segidx);
	return segidx;
}

vector<Matrix*>
importer (const string& path, const string& vzid) {
// XXX deprecated since Fma?
// Two importers for UF - this one just calls the vector version
// and creates a vector of 4 Matrix* with mids:
//    "nodes"
//    "coords"
//    "conn"
//    "gct"
// import from a Universal file-formatted file. Recognized data formats:
//   15   node numbers and coordinates:
//               returned as 2 matrices: a (3*nnodes,1) real matrix "coordinates",
//               and an (nnodes,1) real matrix "node_numbers", the node numbers
//               corresponding to each group of (x,y,z) coordinates
//   55   modes: returned as a (3*nnodes,nmodes) complex Matrix ("gct")
//               with rows ordered dx,dy,dz for each node in the same
//               order as "coordinates"
//   82   tracelines:
//            returned as a (2,nsegment) real Matrix "segidx" with each column a
//            pair of 0-based indexes into the "node_numbers" array, that
//            is, the ordinal node number. Because the "gct" and "coordinates"
//            arrays are in the same order, the start of the coordinates
//            for the jth node is 3*j. Each pair of nodes define a line
//            segment used to draw a picture in amvz.
	Trace trc(1,"UF::importer ",path);
	vector<Matrix*> rval;
	vector<int> nodenumbers;
	vector<int> connectivity;
	vector<double> coordinates;
	vector<complex<double>> gctransform;

	// call the vector version
	UF::importer (path, nodenumbers, coordinates, connectivity, gctransform);

	int nnodes = nodenumbers.size();

	// nodenumbers and connectivity are int - convert to double
	string ext;
	if (!vzid.empty())
		ext = vastr(".",vzid);
	string mid = "nodes" + ext;
	string desc{"nodes"};
	Matrix* nodes = new Matrix(mid, desc, nnodes, 1, false);
	double* dp = nodes->elem();
	for (int i=0; i<nnodes; i++)
		dp[i] = nodenumbers[i];
	rval.push_back(nodes);
	// conn
	mid = "conn" + ext;
	desc = "conn";
	int nconn = connectivity.size();
	Matrix* conn = new Matrix(mid, desc, nconn, 1, false);
	dp = conn->elem();
	for (int i=0; i<nconn; i++)
		dp[i] = connectivity[i];
	rval.push_back(conn);
	// coords
	mid = "coords" + ext;
	desc = "coords";
	int n3 = 3*nnodes;
	Matrix* coords = new Matrix(mid, desc, n3, 1, false);
	coords->data() = coordinates;
	rval.push_back(coords);
	// gct
	int nmodes = gctransform.size()/n3;
	mid = "gct" + ext;
	desc = "gct";
	Matrix* gct = new Matrix(mid, desc, n3, nmodes, true);
	complex<double>* cp = gct->celem();
	for (int i=0; i<n3*nmodes; i++)
		cp[i] = gctransform[i];
	rval.push_back(gct);

	trc.dprint("returning ",rval.size()," matrices");
	return rval;
}

bool
importer (const string& path, vector<int>& nodenumbers, vector<double>& coords,
	vector<int>& segments, vector<complex<double>>& gct) {
// import from a Universal file-formatted file. Recognized data formats:
//   15  node numbers and coordinates:
//        nodenumbers: (nnodes ints)
//        coords:      (3*nnodes doubles) x, y, z coordinates for each node
//   55  gct: (3*nnodes*nmodes) vector of complex<double> in "fortran"
//       storage (column major) with rows ordered dx,dy,dz for each node
//       in the same order as "coords"
//   82  tracelines in penlift format - sets of one or more node numbers
//          separated by zeros (penlift)
	Trace trc(1,"UF::importer(vectors) ",path);
	vector<string> toks;
	vector<string> attributes;

	if (path.empty())
		throw runtime_error("UF::importer called with empty path");

	string actualpath = get_path(expenv(path), ".uf");
	if (actualpath.empty())
		throw runtime_error(vastr("UF::importer: ",path," is not available"));
	Ascii in(actualpath);

	string line;
	// collect all displacement vectors in "modes", copy to gct
	vector<vector<complex<double>>> modes;
	vector<int> gct_node_numbers;   // from DS 55

	try {
		while(in(line, "", false)) {
			toks = string2tok(line, " ");
			if (toks.empty())
				continue;
			if (toks[0] == "-1") {
				in(line, "", false);
				// the next line is the data type: 15, 55, 82 are all we recognize
				toks = string2tok(line, " ");
				if (toks[0] == "15") {
					coords = UF::coordinates(in, nodenumbers);
				} else if (toks[0] == "55") {
					// collect all 55 DS in vector "modes", create vector later
					vector<complex<double>> mp = UF::disp(in, gct_node_numbers);
					if (!mp.empty())
						modes.push_back(mp);
				} else if (toks[0] == "82") {
					// connectivity (tracelines)
					segments = UF::connectivity(in);
				} else if (toks[0] == "151") {
					UF::ignor(in);
				} else if (toks[0] == "164") {
					UF::ignor(in);
				}
			} else {
				throw runtime_error(vastr("unrecognized dataset type ", toks[0]));
			}
		}
	} catch (EOF_exc& s) {
		trc.dprint("caught eof exception, line ",in.line_number(),": ",s);
	} catch (std::exception& s) {
		flaps::warning("caught exception: ",s.what());
		throw runtime_error(vastr("caught unknown exception: ",s.what()));
	}

	// assemble all modes into one column-major array (gct)
	if (!modes.empty()) {
		size_t nr = modes[0].size();
		size_t nc = modes.size();

		gct = vector<complex<double>>(nr*nc);
		complex<double>* cp = &gct[0];
		for (size_t j=0; j<nc; j++) {
			if (modes[j].size() != nr)
				throw runtime_error(vastr("dataset 55 ",j+1," has ",modes[j].size(),
					" elements, others have ",nr));
			blas_copy(nr, &modes[j][0], 1, cp, 1);
			cp += nr;
		}
		// re-order the rows to match nodenumbers from DS 15
		bool need_reorder{false};
		for (size_t i=0; i<nodenumbers.size(); i++) {
			if (gct_node_numbers[i] != nodenumbers[i]) {
				need_reorder = true;
				break;
			}
		}
		if (need_reorder) {
			complex<double>* cp = &gct[0];
			vector<complex<double>> reordered(nr*nc);
			// put rows into reordered
			for (size_t i=0; i<nodenumbers.size(); i++) {
				size_t k;
				for (k=0; k<gct_node_numbers.size(); k++) {
					if (gct_node_numbers[k] == nodenumbers[i]) {
						blas_copy(nc, &cp[3*k], nr, &reordered[3*i], nr);
						blas_copy(nc, &cp[3*k+1], nr, &reordered[3*i+1], nr);
						blas_copy(nc, &cp[3*k+2], nr, &reordered[3*i+2], nr);
						break;
					}
				}
				if (k == gct_node_numbers.size()) {
					string exc{vastr("node ",nodenumbers[i]," was not found in DS 55")};
					throw runtime_error(exc);
				}
			}
			gct = reordered;
		}
	}

	return true;
}

Matrix*
coord(Ascii& in, vector<int>& node_numbers) {
// Read Universal Dataset 15 (nodal data). Returns a Matrix named
// "coordinates" that is dimensioned (3,nnode) where nnode is the
// number of nodes and the 3 rows are the x, y, and z coordinates.
// Also returns the node numbers as an nnode-vector of integers.
	Trace trc(1,"UF::coord");
	string s;
	vector<string> toks;
	vector<double> coords;

	// each line has the FORMAT(4I10,1P3E13.5): 7 numbers:
	// 1) (int) node label
	// 2) (int) definition coordinate system number
	// 3) (int) displacement coordinate system number
	// 4) (int) color
	// 5-7) coordinates of node
	while(in(s, "", false)) {
		toks = string2tok(s, " ");
		if (toks.empty())
			continue;
		if (toks[0] == "-1") {
			break;
		}
		if (toks.size() != 7) {
			throw runtime_error(vastr("line ",in.line_number(),
				": expected 7 items, got ", toks.size()));
		}
		int number;
		str2int(toks[0], number);
		node_numbers.push_back(number);
		double x;
		str2double(toks[4], x);
		double y;
		str2double(toks[5], y);
		double z;
		str2double(toks[6], z);
		coords.push_back(x);
		coords.push_back(y);
		coords.push_back(z);
		// nodes.push_back(Node(number, x, y, z));
	}

	int nnodes = node_numbers.size();
	if (nnodes == 0) {
		throw runtime_error(vastr("line ",in.line_number(),": no nodal data found"));
	}
	string mid{"coordinates"};
	Matrix* rval = new Matrix(mid,"(x,y,z)", 3*nnodes, 1, false);
	rval->data() = coords;
	trc.dprint("returning coordinates",*rval);
	return rval;
}

vector<double>
coordinates(Ascii& in, vector<int>& nodenumbers) {
// Read Universal Dataset 15 (nodal data). Returns a (3,nnode) vector
// of doubles where nnode = nodenumbers.size().
	Trace trc(1,"UF::coord");
	string s;
	vector<string> toks;
	vector<double> coords;

	// each line has the FORMAT(4I10,1P3E13.5): 7 numbers:
	// 1) (int) node label
	// 2) (int) definition coordinate system number
	// 3) (int) displacement coordinate system number
	// 4) (int) color
	// 5-7) coordinates of node
	while(in(s, "", false)) {
		toks = string2tok(s, " ");
		if (toks.empty())
			continue;
		if (toks[0] == "-1") {
			break;
		}
		if (toks.size() != 7) {
			throw runtime_error(vastr("line ",in.line_number(),
				": expected 7 items, got ", toks.size()));
		}
		int number;
		str2int(toks[0], number);
		nodenumbers.push_back(number);
		double x;
		str2double(toks[4], x);
		double y;
		str2double(toks[5], y);
		double z;
		str2double(toks[6], z);
		coords.push_back(x);
		coords.push_back(y);
		coords.push_back(z);
		// nodes.push_back(Node(number, x, y, z));
	}

	int nnodes = nodenumbers.size();
	if (nnodes == 0) {
		throw runtime_error(vastr("line ",in.line_number(),": no nodal data found"));
	}
	trc.dprint("returning coordinates",coords);
	return coords;
}

vector<complex<double>>
disp(Ascii& in, vector<int>& gct_node_numbers) {
// import displacements (modes) from a Universal Dataset 55
	Trace trc(1,"UF::disp");
	string line;
	vector<string> toks;
	bool need_node_numbers = gct_node_numbers.empty();

	// lines 1-5: titles
	// Take the first line as the desc member of the return matrix
	string desc;
	in(desc, "", false);
	// ... then skip the next 4 descriptions
	for (int i=0; i<4; i++)
		in(line, "", false);
	/*
	 * These are the only values in record 6 we recognize:
	 * Record 6: 6I10:
	 *   Field 1 - Model Type
	 ---     1: Structural
	 *   Field 2 - Analysis Type
	 *       2: Normal Mode
	 ---     3: Complex Eigenvalue, first order
	 *   Field 3 - Data Characteristices
	 ---     2: 3 DOF Global Translation Vector
	 *       3: 6 DOF Global Translation, Rotation Vector
	 *   Field 4 - Specific Data Type
	 ---     8: Displacement
	 *   Field 5 - Data Type
	 *       2: Real
	 ---     5: Complex
	 *   Field 6 - Number of data values per node (NDV)
	 */
	in(line, "", false);
	toks = string2tok(line, " ");
	if (toks.size() != 6) {
		string exc{vastr("line ", in.line_number(),
			": expecting 6 values, got ", toks.size())};
		trc.dprint("throwing exception: ",exc);
		throw runtime_error(exc);
	}
	// int model_type = string2Int(toks[0]);
	// int analysis_type = string2Int(toks[1]);
	// int data_characteristic = string2Int(toks[2]);
	// int specific_data_type = string2Int(toks[3]);
	int data_type;
	str2int(toks[4], data_type);
	size_t ndv;
	int kdv;
	str2int(toks[5], kdv);
	ndv = kdv;

	// XXX we should also accept Data Type 2 (real modes)
	if (data_type != 5) {
		string exc{vastr("can only treat complex modes (line ",
			in.line_number(), ')')};
		trc.dprint("throwing exception: ",exc);
		throw runtime_error(exc);
	}
	if (ndv != 3 && ndv != 6) {
		string exc{vastr("the number of data values per node must be",
			" either 3 or 6, not ", ndv," (line ",in.line_number(), ')')};
		trc.dprint("throwing exception: ",exc);
		throw runtime_error(exc);
	}

	// line 7: 4 integers (ignoring)
	in(line, "", false);
	toks = string2tok(line, " ");
	if (toks.size() != 4) {
		string exc{vastr("line ", in.line_number(),
			": expecting 4 values, got ", toks.size())};
		trc.dprint("throwing exception: ",exc);
		throw runtime_error(exc);
	}

	// line 8: 6 floats (ignoring)
	in(line, "", false);
	toks = string2tok(line, " ");
	if (toks.size() != 6) {
		string exc{vastr("line ", in.line_number(),
			": expecting 6 values, got ", toks.size())};
		trc.dprint("throwing exception: ",exc);
		throw runtime_error(exc);
	}
	
	int nfloats = ndv;
	bool cmplx{false};
	vector<complex<double>> rval;
	if (data_type == 5) {
		cmplx = true;
		nfloats *= 2;
	}
	// lines 9 on: node number, values
	while(in(line, "", false)) {
		toks = string2tok(line, " ");
		if (toks.empty())
			continue;
		if (toks[0] == "-1") {
			break;
		}
		if (toks.size() != 1) {
			string exc{vastr("line ", in.line_number(),
				": expected 1 items, got ", toks.size())};
			trc.dprint("throwing exception: ",exc);
			throw runtime_error(exc);
		}
		int nodenumber;
		str2int(toks[0], nodenumber);
		// check that node numbers are the same for each gc transform column
		if (need_node_numbers) {
			gct_node_numbers.push_back(nodenumber);
		} else if (nodenumber != gct_node_numbers[rval.size()/3]) {
			string exc{vastr("line ", in.line_number(), ": node ",nodenumber,
					" is out of order - should be ", gct_node_numbers[rval.size()/3])};
			trc.dprint("throwing exception: ",exc);
			throw runtime_error(exc);
		}

		// next line: ndv floats (2*ndv if complex)
		in(line, "", false);
		toks = string2tok(line, " ");
		if (toks.size() != (size_t)nfloats) {
			string exc{vastr("line ", in.line_number(),
				": expected ", nfloats,
				" items, got ", toks.size())};
			trc.dprint("throwing exception: ",exc);
			throw runtime_error(exc);
		}
		// only take the first 3 displacements, ignore rotations
		int k{0};
		double xi{0.0};
		double xr;
		xr = stod(toks[k++]);
		//!! str2double(toks[k++], xr);
		if (cmplx)
			xi = stod(toks[k++]);
			//!! str2double(toks[k++], xi);
		rval.push_back(complex<double>(xr,xi));
		xr = stod(toks[k++]);
		//!! str2double(toks[k++], xr);
		if (cmplx)
			xi = stod(toks[k++]);
			//!! str2double(toks[k++], xi);
		rval.push_back(complex<double>(xr,xi));
		xr = stod(toks[k++]);
		//!! str2double(toks[k++], xr);
		if (cmplx)
			xi = stod(toks[k++]);
			//!! str2double(toks[k++], xi);
		rval.push_back(complex<double>(xr,xi));
	}
	if (rval.empty()) {
		string exc{vastr("line ", in.line_number(), ": no modal data found")};
		trc.dprint("throwing exception: ",exc);
		throw runtime_error(exc);
	}

	trc.dprint("returning ",rval.size()," displacements");
	return rval;
}

vector<int>
connectivity(Ascii& in) {
// Read a Universal file dataset 82: trace lines in pen-lift format
	Trace trc(1,"UF::connectivity");
	string line;
	vector<string> toks;

	// first line: trace_line_number nnodes color
	in(line, "", false);
	toks = string2tok(line, " ");
	int k;
	str2int(toks[1], k);
	size_t nnodes = k;	// number of nodes + penlifts (zeros)

	// second line: identification
	in(line, "", false);
	string id = replace_char(line, ' ', '_');  // replace blanks

	vector<int> penlift;
	size_t nnew{0};
	while(in(line, "", false)) {
		toks = string2tok(line, " ");
		if (toks.empty())
			continue;
		if (toks[0] == "-1")
			break;
		for (auto& tok : toks) {
			str2int(tok, k);
			penlift.push_back(k);
			nnew++;
		}
	}

	// the number of lines (nnodes) may not match the number in
	// the file - ignore nnodes
	if (nnew != nnodes) {
		flaps::warning("reading ", id, " line ", in.line_number(),
			": expected ", nnodes, " values, got ", nnew);
	}
	trc.dprint("returning penlift conn:",penlift);

	return penlift;
}

void
ignor(Ascii& in) {
// skip this Universal dataset - i.e. read lines until the next -1
	string line;
	vector<string> toks;
	while(in(line, "", false)) {
		toks = string2tok(line, " ");
		if (toks.empty())
			continue;
		if (toks[0] == "-1") {
			break;
		}
	}
}

bool
exporter (const string& path, const vector<Matrix*>& ml, bool append) {
// XXX deprecate: use Fma?
// export a set of matrices to a Universal file. Recognized matrix names
// in ml:
//   nodes        (nnodes) integer node numbers
//   coords  (3*nnodes,1) matrix of (x,y,z) at each node
//   gct          (3*nnodes,nmodes) real or complex matrix with dof
//                in the same order as coordinates
//   conn         real matrix of integers: sets of integers separated
//                by zero (penlift)
// These matrices are written to "path" with UF data types 15, 55, and
// 82, respectively.
// This version just calls the vector version
	Trace trc(1,"UF::exporter(Matrix)");

	// get pointers to each of the 4 matrices
	Matrix* nodes_m{nullptr};
	Matrix* coords_m{nullptr};
	Matrix* conn_m{nullptr};
	Matrix* gct_m{nullptr};
	for (auto mp : ml) {
		string mid{mp->mid()};
		if (mid == "coords")
			coords_m = mp;
		else if (mid == "gct")
			gct_m = mp;
		else if (mid == "conn")
			conn_m = mp;
		else if (mid == "nodes")
			nodes_m = mp;
		else
			throw runtime_error(vastr("cannot export matrix ",mid,
				" to a Universal file"));
	}

	// sanity checks
	// assumed to have 3 displacements per node, i.e. row size is 3*nnodes
	vector<pair<string,Matrix*>> matrices {
		{"nodes", nodes_m},
		{"coords", coords_m},
		{"gct", gct_m},
		{"conn", conn_m} };
	for (auto& mat : matrices)
		if (mat.second == nullptr)
			throw runtime_error(vastr("matrix \"", mat.first,
				"\" was not included in the matrices passed to UF:exporter"));
	// coordinates may be dimensioned (3,nnodes) or (3*nnodes,1)
	int nnodes = nodes_m->rsize()*nodes_m->csize();
	if (coords_m->rsize()*coords_m->csize() != (size_t)3*nnodes)
		throw runtime_error(vastr("mismatched dimensions: ",nnodes,
			" nodes numbers, ", coords_m," coords"));
	// gct must be (3*nnodes, nc)
	int nr = gct_m->rsize();
	int nc = gct_m->csize();
	if (nr != 3*nnodes)
		throw runtime_error(vastr("dimension of gct (",nr,',',nc,
				") do not agree with coordinates (", 3*nnodes,")"));

	// 15: nodes are double in nodes_m - cast to ints
	vector<int> nodes;
	for (auto i : nodes_m->data())
		nodes.push_back(i);
	// 82: conn in "penlift" format - doubles in "conn_m"
	//     must be cast to integer
	vector<int> conn;
	for (auto i : conn_m->data())
		conn.push_back(i);

	// displacements (55) if real convert to complex
	int nrgct = gct_m->rsize();
	int ncgct = gct_m->csize();
	vector<complex<double>> gct(nrgct*ncgct);
	if (!gct_m->is_complex())
		blas_copy(nrgct*ncgct, gct_m->elem(), 1,
			reinterpret_cast<double*>(&gct[0]), 2);
	else
		blas_copy(nrgct*ncgct, gct_m->celem(), 1, &gct[0], 1);

	const vector<double>& coords = coords_m->data();

	// then just call the vector version:
	return UF::exporter(path, nodes, coords, conn, gct);
}

bool
exporter (const string& path, const vector<int> nodes,
	const vector<double>& coords, const vector<int>& conn,
	const vector<double>& gct) {
// Real gct version: cast to complex, call complex version
	vector<complex<double>> cgct(gct.size());
	for (size_t i=0; i<gct.size(); i++)
		cgct[i] = gct[i];
	return UF::exporter(path,nodes,coords,conn,cgct);
}

bool
exporter (const string& path, const vector<int> nodes,
	const vector<double>& coords, const vector<int>& conn,
	const vector<complex<double>>& gct) {
// export to a Universal-file formatted file.
// Input
//   path      path of the file to be written to
//   nodes     (nnodes int) node numbers
//   coords    (3*nnodes) matrix of (x,y,z) at each node
//   conn      int vector of node numbers describing line segments
//             in penlift format: a zero indicates a penlift.
//   gct       (3*nnodes,nmodes complex) column-major matrix with rows
//             in the same order as "coords"
// These are written to "path" with UF data types 15, 82, and
// 55, respectively.
	Trace trc(1,"UF::exporter");

	// open the .uf file
	string actualpath{expenv(path)};
	if (rsubstr(actualpath, 3) != ".uf")
		actualpath += ".uf";
	ofstream file(actualpath);
	if (!file)
		throw runtime_error(vastr("could not open ",actualpath));

	// sanity checks
	// assumed to have 3 displacements per node, i.e. row size is 3*nnodes
	int nnodes = nodes.size();
	int n3 = coords.size();
	if (n3 != 3*nnodes)
		throw runtime_error("there in not 3 coordinates per node");
	int ngct = gct.size();
	int nc = ngct/n3;
	if (ngct != n3*nc)
		throw runtime_error(vastr("gct size (",ngct,") wrong: ",
			"not a multiple of coords (",n3,")"));

	// 82: connectivity in "penlift" format
	file << "    -1\n";
	file << "    82\n";
	// record 1 (3I10): trace_line_number number_of_nodes color
	file << "         1" << setw(10) << conn.size() << "         1\n";
	file << "Flaps " << VERSION << " " << get_uname("version") << endl;
	for (size_t i=0; i<conn.size(); i++) {
		file << setw(10) << conn[i];
		if ((i+1)%8 == 0)
			file << endl;
	}
	if (conn.size()%8 != 0)
		file << endl;
	file << "    -1\n";

	// node numbers and coordinates (15)
	Form sci(5);
	sci.width(13);
	file << "    -1\n";
	file << "    15\n";
	for (int i=0; i<nnodes; i++) {
		// node number
		file << setw(10) << nodes[i];
		// ignore freedoms
		file << setw(10) << 0;
		file << setw(10) << 0;
		file << setw(10) << 1;
		int k = 3*i;
		file << sci(coords[k]);
		file << sci(coords[k+1]);
		file << sci(coords[k+2]);
		file << endl;
	}
	file << "    -1\n";

	// displacements (55)
	// record 6 fields:
	int model_type{1};
	int analysis_type{3};
	int data_characteristic{2};
	int specific_data_type{8};
	int data_type{5};  // complex
	int values_per_node{3};
	// record 7 fields:
	int number_of_integer_data_values{2};
	int number_of_real_data_values{6};
	int load_case{1};
	// type 55 for each mode
	for (int k=0; k<nc; k++) {
		file << "    -1\n";
		file << "    55\n";
		// records 1-5: id
		for (int i=0; i<5; i++)
			file << "gct" << k+1 << endl;
		// record 6
		file << setw(10) << model_type;
		file << setw(10) << analysis_type;
		file << setw(10) << data_characteristic;
		file << setw(10) << specific_data_type;
		file << setw(10) << data_type;
		file << setw(10) << values_per_node << endl;
		// record 7
		file << setw(10) << number_of_integer_data_values;
		file << setw(10) << number_of_real_data_values;
		file << setw(10) << load_case;
		file << setw(10) << k+1 << endl;
		// record 8: just put all zeros
		for (int i=0; i<6; i++)
			file << sci(0.0);
		file << endl;
		// records 9 & 10 repeat for each node
		for (int j=0; j<nnodes; j++) {
			file << setw(10) << nodes[j] << endl;
			// gct is (3,nnodes,nc)
			for (int i=0; i<3; i++) {
				file << sci(gct[IJK(i,j,k,3,nnodes)].real());
				file << sci(gct[IJK(i,j,k,3,nnodes)].imag());
			}
			file << endl;
		}
		file << "    -1\n";
	}

	// print user info
	flaps::info ("Export to Universal file ",actualpath,"\n    ",nnodes,
			" nodes\n    ",conn.size()," connectivities\n",
			"    3 diplacements at ",nnodes," nodes");

	return true;
}
} // namespace UF


// namespace Matlab {

// Reference:
// 1) "MATLAB MAT-file Format" available as a pdf file:
//    https://www.mathworks.com/help/pdf_doc/matlab/matfile_format.pdf
//    March 2021 Version 9.10 (Release 2021a)

// MAT-file data types: ref. p 1-6
int miINT8 = 1;
int miUINT8 = 2;
int miINT16 = 3;
int miUINT16 = 4;
int miINT32 = 5;
int miUINT32 = 6;
int miSINGLE = 7;
int miDOUBLE = 9;
int miINT64 = 12;
int miUINT64 = 13;
int miMATRIX = 14;
int miCOMPRESSED = 15;
int miUTF8 = 16;
int miUTF16 = 17;
int miUTF32 = 18;

static bool
bigendian () {
// is this machine big-endian?
	constexpr int n{1};
	if (*(char*)&n == 1)
		return false;
	return true;
}

bool
check_matfile(string const& path) {
	bool rval = false;
// Test to see if "path" is a matlab file by checking the header and
// endianness bytes; return true if they conform with this machines endianness,
// throws an exception otherwise: system_error if the file is not available,
// a runtime_error exception if the file is corrupt or the wrong endianness

	// open the file read-only using fio class Fd: RAII, throws
	// exception if open fails
	Fd pathfd (path, O_RDONLY);
	int fd = pathfd();

	char buf[256];
	// header (128 bytes)
	read(fd, buf, 128);

	// the first bytes should be:
	string first{"MATLAB 5.0"};
	string head(buf, 10);
	if (head != first) {
		throw runtime_error(vastr("bad header: should start with ",first,
					", but it is ",head));
	}
	// bytes 127,128 (126,127 zero-based) indicate endianness -
	// if they are IM the file was created on a little-endian machine;
	// if they are MI it was created on a big-endian machine.
	// Return false if these bytes are not M and I (meaning this is
	// not a matlab file), and print a warning if they are but in the
	// wrong order.
	// We could do byte-swapping to correct wrong endianness but
	// XXX don't bother for now
	if (bigendian()) {
		if (buf[126] == 'M' && buf[127] == 'I') {
			rval = true;
		} else if (buf[127] == 'M' && buf[126] == 'I') {
			throw std::system_error(errno, std::generic_category(),
				vastr(path," was created on a little-endian machine but "
				"this is a big-endian machine"));
		}
	} else {
		if (buf[126] == 'I' && buf[127] == 'M') {
			rval = true;
		} else if (buf[126] == 'M' && buf[127] == 'I') {
			throw std::system_error(errno, std::generic_category(),
				vastr(path," was created on a big-endian machine but "
				"this is a little-endian machine"));
		}
	}
	return rval;
}

static int32_t
align64(int32_t nbytes) {
// Returns nbytes rounded up to a multiple of 8 so it aligns
// to a 64-bit (8 byte) boundary.
// Page 1-6 of ref. 1 states: All data that is uncompressed must be
// aligned on 64-bit boundaries.
	int32_t rval = ceil((nbytes*8.0)/64.0)*8;
	return rval;
}

template<typename MT, typename T>
bool
readdt(int fd, vector<T>& data, int32_t nbytes) {
// template to read nbytes from fd of data type MT data and
// convert each item to data type T, returning them in vector "data"
	int nitems = nbytes/sizeof(MT);
	assert(nitems*sizeof(MT) == (size_t)nbytes);

	// initialize output vector
	data = vector<T>(nitems, 0);

	// if nbytes is not a multiple of 8 we must use
	// a buffer with 8-byte alignment
	int32_t nb = align64(nbytes);

	// if type MT and T are the same, read directly into "data"
	if (nb == nbytes && typeid(MT) == typeid(T)) {
		read(fd, &data[0], nbytes);
	} else {
		vector<MT> buf(nb/sizeof(MT));
		read(fd, &buf[0], nb);
		for (int i=0; i<nitems; i++)
			data[i] = buf[i];
	}
	return true;
}

template<typename T>
bool
readdtype (int fd, int32_t data_type, int32_t nbytes, vector<T>& data) {
// template to read nbytes from fd of Matlab data type "data_type" and
// convert each item to data type T, returning them in vector "data"
	switch (data_type) {
		case 1:
			return readdt<int8_t,T>(fd, data, nbytes);
		case 2:
			return readdt<uint8_t,T>(fd, data, nbytes);
		case 3:
			return readdt<int16_t,T>(fd, data, nbytes);
		case 4:
			return readdt<uint16_t,T>(fd, data, nbytes);
		case 5:
			return readdt<int32_t,T>(fd, data, nbytes);
		case 6:
			return readdt<uint32_t,T>(fd, data, nbytes);
		case 7:
			return readdt<float,T>(fd, data, nbytes);
		case 8:
			throw runtime_error("data type 8 is reserved");
		case 9:
			return readdt<double,T>(fd, data, nbytes);
		case 10:
			throw runtime_error("data type 10 is reserved");
		case 11:
			throw runtime_error("data type 11 is reserved");
		case 12:
			return readdt<int64_t,T>(fd, data, nbytes);
		case 13:
			return readdt<uint64_t,T>(fd, data, nbytes);
		case 14:
			throw runtime_error("cannot import MATLAB arrays");
		case 15:
			throw runtime_error("cannot import compressed MATLAB data yet");
		case 16:
			return readdt<char,T>(fd, data, nbytes);
		case 17:
		case 18:
			throw runtime_error("cannot import UTF-16 or -32 encoded MATLAB data");
	}
	return 0.0;
}

static bool
read_tag(int fd, int32_t& data_type, int32_t& nbytes, int32_t& data) {
// Read an 8-byte tag from file descriptor "fd" and extract the data
// type and number of bytes. This special function is to allow for
// Matlab "small data" where data_type, nbytes, and the data are all
// packed into the 8-byte tag.
// If this is a small-data tag, returns true and "data" is set to
// the 4 bytes of data; otherwise true is returned and "data" is unchanged.
// If eof is encountered an EOF_exc exception is thrown.
	Trace trc(1,"read_tag");
	unsigned char buf[8];

	// watch for EOF
	if (read(fd, buf, 8) == 0) {
		throw EOF_exc("reading a tag");
	}

	// check to see if this is small data: taking the first 4 bytes as
	// int32_t, the 2 high-order bytes are "nbytes" and the 2 low-order
	// bytes are "data_type"; if it is not small the 2 high-order bytes will
	// be zero since then these will be the high-order bytes of an int32_t
	// data type which ranges from 1 to 18.
	bool small_data{false};
	if (bigendian()) {
		if (buf[0] != 0 || buf[1] != 0)
			small_data = true;
	} else {
		if (buf[2] != 0 || buf[3] != 0)
			small_data = true;
	}
	// If this is a small-data tag the 2 high-order bytes are nbytes and
	// the 2 low-order bytes are data_type...
	if (small_data) {
		int16_t nb;
		int16_t dt;
		if (bigendian()) {
			memcpy(&nb, &buf[0], 2);
			memcpy(&dt, &buf[2], 2);
		} else {
			memcpy(&nb, &buf[2], 2);
			memcpy(&dt, &buf[0], 2);
		}
		trc.dprint("small data format: data type ",dt,", nbytes ",nb,", data ", buf[4],buf[5],buf[6],buf[7]);
		data_type = dt;
		nbytes = nb;
		// small data is the 2nd 4 bytes
		memcpy (&data, &buf[4], 4);
	// ... otherwise the first 4 bytes are the data type and the
	// second 4 bytes are the number of bytes.
	} else {
		memcpy(&data_type, buf, 4);
		memcpy(&nbytes, &buf[4], 4);
	}
	if (data_type <= 0 || nbytes <= 0) {
		ostringstream os;
		os << "bad tag: ";
		os << showbase << internal << setfill('0');
		for (size_t i=0; i<8; i++)
			os << std::hex << " " << buf[i];
		throw runtime_error(os.str());
	}
	trc.dprint("returning ",small_data,", data type ",data_type,", nbytes ",nbytes);
	return small_data;
}

static Matrix* read_Matrix(int fd);

bool
period2underscore(string& str) {
// replace ALL periods in str with double underscore
	Trace trc(2,"period2underscore");
	bool rval{false};
	trc.dprint("str: ",str);
	string::size_type idx = str.find(".");
	while(idx != string::npos) {
		str.replace(idx,1,"__");
		idx = str.find(".");
		rval = true;
	}
	trc.dprint("output str: ",str);
	return rval;
}

bool
underscore2period(string& str) {
// replace ALL double-underscore in str with periods
	Trace trc(2,"underscore2period");
	bool rval{false};
	trc.dprint("str: ",str);
	string::size_type idx = str.find("__");
	while(idx != string::npos) {
		str.replace(idx,2,".");
		idx = str.find("__");
		rval = true;
	}
	trc.dprint("output str: ",str);
	return rval;
}

vector<Matrix*>
Matlab::
importer(const string& path, const string& output) {
// read matrices from Matlab file "path"
// If a matrix name contains double underscore (__) change it to .
// output is unused
	Trace trc(1,"Matlab::importer ",path);
	vector<Matrix*> rval;

	// first check that the file exists, is a matlab file, and
	// has the same endianness as this machine
	if (!check_matfile(path)) {
		throw runtime_error(vastr(path," is not a Matlab file"));
	}

	// open the file read-only, using class Fd for RAII
	Fd pathfd (path, O_RDONLY);
	int fd = pathfd();

	// header (128 bytes) - checks already done in check_matfile()
	char buf[128];
	read(fd, buf, 128);

	// read matrices until EOF - read_Matrix will return nullptr on EOF
	Matrix* mp;
	while((mp = read_Matrix(fd)) != nullptr)
		rval.push_back(mp);
	trc.dprint("returning ",rval.size()," matrices");
	return rval;
}	// Matlab::importer

static Matrix*
read_Matrix(int fd) {
// read a Matlab array of a numeric data type (1-13)
	Trace trc(1,"Matlab::read_Matrix ");
	size_t nr, nc;
	int32_t small_data; // for small data tags
	Matrix* rval{nullptr};

	// first tag must be Array Data Element (data_type == miMATRIX)
	int32_t data_type;
	int32_t nbytes;
	try {
		read_tag(fd, data_type, nbytes, small_data);
		trc.dprint("data type ",data_type,", nbytes ",nbytes);
		if (data_type != miMATRIX) {
			if (data_type == miCOMPRESSED) {
				throw runtime_error("cannot treat MATLAB compressed MAT files: "
						"use the -nocompression option when saving in MATLAB");
			} else {
				throw runtime_error(vastr("cannot treat MATLAB data type ",data_type));
			}
		}
	} catch (EOF_exc& s) {
		trc.dprint("caught eof exception");
		return nullptr;
	}

	// Array Flags subelement tag
	read_tag(fd, data_type, nbytes, small_data);
	trc.dprint("Array flags data type ",data_type,", nbytes ",nbytes);
	// ... and the Array Flags data: 2 uint32_t:
	// 1) flags, class in the 2 low-order bytes
	// 2) unused
	// read as 8 unint8_t so we can just pick out the 2 bytes
	vector<uint8_t> flags;
	readdtype (fd, miUINT8, nbytes, flags);
	uint8_t af_class;
	uint8_t af_flags;
	if (bigendian()) {
		af_class = flags[3];
		af_flags = flags[2];
	} else {
		af_class = flags[0];
		af_flags = flags[1];
	}
	trc.dprint("af class ",(int)af_class,", af flags ",(int)af_flags);
	bool is_complex{false};
	if (af_flags >> 3 == 1)
		is_complex = true;

	// Dimensions Array subelement tag: n int32 (n dimensions)
	read_tag(fd, data_type, nbytes, small_data);
	trc.dprint("Dimensions array data type ",data_type,", nbytes ",nbytes);
	// the number of dimensions is nbytes/sizeof(int32)
	int ndim = nbytes/4;
	trc.dprint("ndim ",ndim,", complex? ",is_complex);
	// read each dimension
	vector<int32_t> dim(ndim);
	readdtype (fd, miINT32, nbytes, dim);
	// set nr, nc: ref p 1-11 says "All numeric arrays have at least
	// two dimensions."
	if (dim.size() == 2) {
		nr = dim[0];
		nc = dim[1];
	} else {
		throw runtime_error(vastr("cannot treat ",ndim," dimensional matrices"));
	}

	// Array Name subelement tag - the ref says these are miINT8 so
	// we cast to char
	string name;
	if (read_tag(fd, data_type, nbytes, small_data)) {
		name = string(reinterpret_cast<char*>(&small_data), nbytes);
	} else {
		vector<char> charvec;
		readdtype (fd, data_type, nbytes, charvec);
		name = string(&charvec[0], nbytes);
	}
	// replace double-underscores with periods
	underscore2period(name);
	trc.dprint("name: ",name);

	// check to see if this matrix exists - if so read it and
	// take its parameterizations
	vector<string> cat = fio::catalog(name);
	vector<Pz*> old_pz;
	if (!cat.empty()) {
		Matrix* existing = new Matrix(name);
		old_pz = existing->pz;
	}
	// Real Part subelement tag - possibly small data?
#ifdef NEVER // for amviz: give matrix desc=mid
	Matrix* real_part = new Matrix(name, "matlab", nr, nc, false);
#else // NEVER // for amviz: give matrix desc=mid
	Matrix* real_part = new Matrix(name, name, nr, nc, false);
#endif // NEVER // for amviz: give matrix desc=mid
	if (read_tag(fd, data_type, nbytes, small_data)) {
		*real_part->elem() = small_data;
	} else {
		trc.dprint("real part data type ",data_type,", nbytes ",nbytes);
		readdtype (fd, data_type, nbytes, real_part->data());
	}

	if (is_complex) {
#ifdef NEVER // for amviz: give matrix desc=mid
		Matrix* cmat = new Matrix(name, "matlab", nr, nc, true);
#else // NEVER // for amviz: give matrix desc=mid
		Matrix* cmat = new Matrix(name, name, nr, nc, true);
#endif // NEVER // for amviz: give matrix desc=mid
		complex<double>* cp = cmat->celem();
		// imag Part subelement tag: watch for small data
		if (read_tag(fd, data_type, nbytes, small_data)) {
			double ci = small_data;
			double* rp = real_part->elem();
			cp[0] = complex<double>(rp[0], ci);
		} else {
			trc.dprint("imag part data type ",data_type,", nbytes ",nbytes);
			// copy the real_part into cmat...
			blas_copy (nr*nc, real_part->elem(), 1, cmat->elem(), 2);
			// ... so we can re-use real_part to read the imaginary part...
			readdtype (fd, data_type, nbytes, real_part->data());
			// ... then copy it into cmat
			blas_copy (nr*nc, real_part->elem(), 1, cmat->elem()+1, 2);
		}
		rval = cmat;
		delete real_part;
	} else {
		rval = real_part;
	}
	// give the matrix the parameterization from the existing one
	rval->pz = old_pz;

	return rval;
}

// Matlab exporter
static
int32_t
matrix_nbytes(Matrix* mp) {
// compute the number of bytes a Matrix will take in a MAT file.
// This includes all subelements and their tags
	Trace trc(1,"matrix_nbytes");
	int32_t rval{0};
	int32_t name_size = mp->mid().size();
	int32_t nr = mp->rsize();
	int32_t nc = mp->csize();
	// array flags
	rval += 8; // tag
	rval += 8; // data
	// dimensions array
	rval += 8; // tag
	rval += 8; // 2 dimensions
	// array name
	rval += 8; // tag
	rval += align64(name_size);
	// real part: assume double data
	rval += 8; // tag
	rval += 8*nr*nc;
	// imag part
	if (mp->is_complex()) {
		rval += 8; // tag
		rval += 8*nr*nc;
	}
	trc.dprint("returning ",rval," bytes");
	return rval;
}

template<typename MT, typename T>
bool
writedt(int fd, vector<T>& data, int32_t nbytes) {
// template to convert "data" from data type T to data type MT,
// then write to fd; nbytes is only used to verify the value
// written to the tag.
	int nitems = nbytes/sizeof(MT);
	assert(nitems*sizeof(MT) == (size_t)nbytes);
	if (data.size()*sizeof(T) != (size_t)nbytes) {
		throw runtime_error(vastr(data.size()*sizeof(T)," bytes of data, but ",
					" writedt called with nbytes ",nbytes));
	}

	// if nbytes is not a multiple of 8 we must use
	// a buffer with 8-byte alignment
	int32_t nb = align64(nbytes);

	// if type MT and T are the same and no padding necessary, write directly
	if (nb == nbytes && typeid(MT) == typeid(T)) {
		write(fd, &data[0], nbytes);
	} else {
		vector<MT> buf(nb/sizeof(MT));
		for (int i=0; i<nitems; i++)
			buf[i] = data[i];
		write(fd, &buf[0], nb);
	}
	return true;
}

template<typename T>
bool
writedtype (int fd, int32_t data_type, int32_t nbytes, vector<T>& data) {
// template to write nbytes of data to file descriptor fd.
// The data is type T, and it is to be written as data type "data_type".
// Matlab data types are int32_t documented in table 1-1, reference 1,
// page 1-6.
	switch (data_type) {
		case 1:
			return writedt<int8_t,T>(fd, data, nbytes);
		case 2:
			return writedt<uint8_t,T>(fd, data, nbytes);
		case 3:
			return writedt<int16_t,T>(fd, data, nbytes);
		case 4:
			return writedt<uint16_t,T>(fd, data, nbytes);
		case 5:
			return writedt<int32_t,T>(fd, data, nbytes);
		case 6:
			return writedt<uint32_t,T>(fd, data, nbytes);
		case 7:
			return writedt<float,T>(fd, data, nbytes);
		case 8:
			throw runtime_error("data type 8 is reserved");
		case 9:
			return writedt<double,T>(fd, data, nbytes);
		case 10:
			throw runtime_error("data type 10 is reserved");
		case 11:
			throw runtime_error("data type 11 is reserved");
		case 12:
			return writedt<int64_t,T>(fd, data, nbytes);
		case 13:
			return writedt<uint64_t,T>(fd, data, nbytes);
		case 14:
			throw runtime_error("cannot import MATLAB arrays");
		case 15:
			throw runtime_error("cannot import compressed MATLAB data yet");
		case 16:
		case 17:
		case 18:
			throw runtime_error("cannot import UTF encoded MATLAB data");
	}
	return 0.0;
}

bool
write_header(int fd) {
// Write a level-5 MAT-file header (128 bytes)

	vector<char> header(128, ' ');
	string text{vastr("MATLAB 5.0 MAT-file, Platform: ",get_uname("nodename"),
		", Created on ", stringTimeDate())};
	// bytes 117 to 124 contain an offset to subsystem-specific data
	// so only copy at most 116 char
	memcpy(&header[0], text.c_str(), std::min((size_t)116, text.size()));
	// the last 4 bytes are 2 16-bit flag fields: version & endian indicator
	if (bigendian()) {
		header[124] = 1;
		header[125] = 0;
		header[126] = 'M';
		header[127] = 'I';
	} else {
		header[124] = 0;
		header[125] = 1;
		header[126] = 'I';
		header[127] = 'M';
	}
	writedtype (fd, miINT8, 128, header);
	return true;
}

bool
Matlab::
exporter (const string& path, const vector<Matrix*>& ml, bool append) {
	Trace trc(1,"Matlab::exporter");
	bool rval{true};

	trc.dprint("path<",path,">, ",ml.size()," matrices");

	// open the file write-only, truncate if !append
	mode_t mode{S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH};
	int flags{O_WRONLY | O_CREAT};
	if (append)
		flags |= O_APPEND;
	else
		flags |= O_TRUNC;
	// use RAII class Fd to open, not necessary to close
	Fd pathfd(path, flags, mode);
	int fd = pathfd();

	// Level-5 MAT-file header 
	write_header(fd);

	// all tags should use this
	vector<int32_t> tag(2,0);

	for (auto mp : ml) {
		// replace periods in the mid with double-underscore; do this
		// BEFORE matrix_nbytes as it may change size of name
		string name{mp->mid()};
		period2underscore(name);
		// a Matrix has a Data Element tag followed by 4 or 5
		// subelements, each with a tag and data:
		// 1) Array flags
		// 2) Dimensions Array
		// 3) Array name
		// 4) Real part
		// 5) Imaginary part (if complex)
		int32_t data_type;
		int32_t nbytes;

		// Array Data Element tag
		tag[0] = miMATRIX;
		tag[1] = matrix_nbytes(mp);
		writedtype (fd, miINT32, 8, tag);

		// 1) Array Flags subelement tag
		nbytes = 8;
		tag[0] = miUINT32;
		tag[1] = 2*sizeof(uint32_t);
		writedtype (fd, miINT32, 8, tag);
		// ... and the Array Flags data: 2 uint32_t:
		// 1) flags, class in the 2 low-order bytes
		// 2) unused
		// write as 8 unint8_t so we can just set flags & class
		vector<uint8_t> flags(8,0);
		uint8_t af_class{6};
		uint8_t af_flags{0};
		if (mp->is_complex())
			af_flags = 8;    // 1 in bit 4 => complex
		if (bigendian()) {
			flags[3] = af_class;
			flags[2] = af_flags;
		} else {
			flags[0] = af_class;
			flags[1] = af_flags;
		}
		writedtype (fd, miUINT8, nbytes, flags);

		// Dimensions Array subelement tag...
		int ndim = 2;   // always for Matrix
		data_type = miINT32;
		nbytes = ndim*sizeof(int32_t);
		tag[0] = data_type;
		tag[1] = nbytes;
		writedtype (fd, miINT32, 8, tag);
		// ... and the dimensions
		vector<int32_t> dim(2);
		size_t nr = mp->rsize();
		size_t nc = mp->csize();
		dim[0] = nr;
		dim[1] = nc;
		writedtype (fd, data_type, nbytes, dim);

		// Array Name subelement tag with the actual name size.
		// The ref says these are miINT8 so we cast to char
		data_type = miINT8;
		nbytes = name.size();
		tag[0] = data_type;
		tag[1] = nbytes;
		writedtype (fd, miINT32, 8, tag);
		vector<char> charvec(nbytes);
		for (int32_t i=0; i<nbytes; i++)
			charvec[i] = name[i];
		writedtype (fd, data_type, nbytes, charvec);

		// real part subelement tag...
		data_type = miDOUBLE;
		nbytes = nr*nc*sizeof(double);
		tag[0] = data_type;
		tag[1] = nbytes;
		writedtype (fd, miINT32, 8, tag);
		// split into real & imag parts
		vector<double> real_data(nr*nc);
		vector<double> imag_data;
		if (mp->is_complex()) {
			blas_copy(nr*nc, mp->elem(), 2, &real_data[0], 1);
			imag_data = vector<double>(nr*nc);
			blas_copy(nr*nc, mp->elem()+1, 2, &imag_data[0], 1);
		} else {
			blas_copy(nr*nc, mp->elem(), 1, &real_data[0], 1);
		}

		// write the real data
		writedtype(fd, data_type, nbytes, real_data);

		// ... and the imaginary part
		if (!imag_data.empty()) {
			// Imag Part subelement tag: same as real part
			writedtype (fd, miINT32, 8, tag);
			// ... and the imag data
			writedtype(fd, data_type, nbytes, imag_data);
		}
	}

	return rval;
}

// } namespace Matlab

// Matrix Market: namespace MM
// namespace MM {

Matrix*
MM::
importer (const string& path, const string& output) {
// import a single matrix from a Matrix Market formatted file
	Trace trc(1,"MM::importer ",path);
	Matrix* rval{nullptr};

	// the default output matrix name is the path stripped of
	// extension & directories
	string mid{output};
	if (output.empty())
		mid = stringBasename(path);

	int nr;
	int nc;
	bool is_complex;
	vector<double> mat = MM::importer (path, nr, nc, is_complex);

	// create the return Matrix
	rval = new Matrix(mid, "MatrixMarket",nr,nc,is_complex);
	rval->data() = mat;
	trc.dprint("returning (",nr,",",nc,") ",is_complex?"complex":"real");
	return rval;
}

vector<double>
MM::
importer (const string& name, int& nr, int& nc, bool& is_complex) {
// import a single matrix from a Matrix Market formatted file
	Trace trc(1,"MM::importer(vector) ",name);
	int entries;
	int i, j;
	string line;
	vector<string> toks;

	// check access & get full path
	string actualpath = get_path(expenv(name), "mm");
	if (actualpath.empty())
		throw runtime_error(vastr("Matrix Market file \"",name,"\" is not available"));

	Ascii in(actualpath);

	// the first line looks like:
	//   %%MatrixMarket matrix coordinate real general
	in(line, "", false);
	// ECMAScript grammar, case-sensitive
	regex infore{R"(^%%MatrixMarket[ ]*matrix[ ]*(.*)[ ][ ]*([^ ]*) [ ]*([^ ]*))"};
	smatch mch;
	if (regex_match(line, mch, infore)) {
		string formatstr = mch[1];
		string dtypestr = mch[2];
		trc.dprint("got format <",formatstr,"> datatype <",dtypestr,">");
		if (formatstr.substr(0,5) != "coord") {
			string exc = vastr("cannot read ", formatstr,
					"-formatted Matrix Market files");
			trc.dprint("throwing exception: ",exc);
			throw runtime_error(exc);
		}
		// real or complex?
		if (compare(dtypestr, "com"))
			is_complex = true;
		else if (compare(dtypestr, "re"))
			is_complex = false;
		else {
			string exc = vastr("cannot read datatype \"", dtypestr,"\"");
			trc.dprint("throwing exception: ",exc);
			throw runtime_error(exc);
		}
	}

	// Next line is rows, columns, entries (toss comments)
	in(line, "%", false);
	toks = string2tok(line, " ");
	if (toks.size() != 3) {
		string exc = vastr("expecting 3 integers (number of rows, ",
				"columns, elements), got: ", line);
		trc.dprint("throwing exception: ",exc);
		throw runtime_error(exc);
	}
	str2int(toks[0], nr);
	str2int(toks[1], nc);
	str2int(toks[2], entries);
	trc.dprint("matrix is ",nr," by ",nc,", ",entries," entries");

	// Sanity check
	if (nr <= 0 || nr > 10000) {
		string exc = vastr("line ", in.line_number(), " of ", actualpath,
			": bad number of rows (", nr, ')');
		trc.dprint("throwing exception: ",exc);
		throw runtime_error(exc);
	}
	if (nc <= 0 || nc > 10000) {
		string exc = vastr("line ", in.line_number(), " of ", actualpath,
			": bad number of columns (", nc, ')');
		trc.dprint("throwing exception: ",exc);
		throw runtime_error(exc);
	}

	int nd{nr*nc};
	if (is_complex)
		nd *= 2;
	vector<double> rval(nd, 0.0);
	if (is_complex) {
		// a complex line has: i j real imag
		complex<double>* cp = (complex<double>*)&rval[0];
		double xr, xi;
		while(in(line, "%", false)) {
			toks = string2tok(line, " ");
			if (toks.empty())
				continue;
			if (toks.size() != 4) {
				string exc = vastr("expecting 4 items on line ", in.line_number(),
					", got: ", line);
				trc.dprint("throwing exception: ",exc);
				throw runtime_error(exc);
			}
			str2int(toks[0], i);
			str2int(toks[1], j);
			if (i < 1 || i > nr || j < 1 || j > nc) {
				string exc = vastr("bad indices on line ", in.line_number(),
					": (", i, ',', j, ")");
				trc.dprint("throwing exception: ",exc);
				throw runtime_error(exc);
			}
			str2double(toks[2], xr);
			str2double(toks[3], xi);
			cp[(i-1)+((j-1)*nr)] = complex<double>(xr, xi);
		}
		trc.dprintm(2*nr,nc,2*nr,cp,"read Complex Matrix Market");
	} else {
		// real data: i j x
		double* rp = &rval[0];
		while(in(line, "%", false)) {
			toks = string2tok(line, " ");
			if (toks.empty())
				continue;
			if (toks.size() != 3) {
				string exc = vastr("expecting 3 items on line ", in.line_number(),
					", got: ", line);
				trc.dprint("throwing exception: ",exc);
				throw runtime_error(exc);
			}
			str2int(toks[0], i);
			str2int(toks[1], j);
			if (i < 1 || i > nr || j < 1 || j > nc) {
				string exc = vastr("bad indices on line ", in.line_number(),
					": (", i, ',', j, ")");
				trc.dprint("throwing exception: ",exc);
				throw runtime_error(exc);
			}
			str2double(toks[2], rp[(i-1)+((j-1)*nr)]);
		}
		trc.dprintm(nr,nc,nr,rp,"read real Matrix Market");
	}
	trc.dprint("returning (",nr,",",nc,") ",is_complex?"complex":"real");
	return rval;
}


static string
fullpath (std::string const& path, std::string const& title);

bool
MM::
exporter (const string& path, const string& title,
		const double* rp, int nr, int nc) {
	Trace trc(1,"MM::exporter(double)");
	std::ofstream file;
	bool rval{true};

	trc.dprint("path <",path,"> nr ",nr,", nc ",nc);

	string mmpath = fullpath(path, title);
	if (rsubstr(mmpath,3) != ".mm")
		mmpath += ".mm";

	file.open(mmpath, std::ios::trunc);
	if (!file)
		throw std::system_error(errno, std::generic_category(),
			vastr("cannot open ", mmpath));

	// Count the number of non-zero elements for the header
	int nnz{0};
	for (int i=0; i<nr*nc; i++)
		if (std::abs(rp[i]) > 0.0)
			nnz++;

	file << "%%MatrixMarket matrix coordinate real general\n";
	file << nr << ' ' << nc << ' ' << nnz << std::endl;

	for (int j=0; j<nc; j++) {
		for (int i=0; i<nr; i++) {
			double x = rp[(i+j*nr)];
			if (std::abs(x) > 0.0)
				file << i+1 << ' ' << j+1 << ' ' << x << std::endl;
		}
	}

	return rval;
}

bool
MM::
exporter (const string& path, const string& title,
		const complex<double>* cp, int nr, int nc) {
	Trace trc(1,"MM::exporter(complex)");
	std::ofstream file;
	std::string mmpath;
	bool rval{true};

	trc.dprint("path <",path,"> nr ",nr,", nc ",nc);

	mmpath = fullpath(path, title);
	if (rsubstr(mmpath,3) != ".mm")
		mmpath += ".mm";

	file.open(mmpath, std::ios::trunc);
	if (!file)
		throw std::system_error(errno, std::generic_category(),
			vastr("cannot open ", mmpath));

	// Count the number of non-zero elements for the header
	int nnz = 0;
	for (int i=0; i<nr*nc; i++)
		if (std::abs(cp[i]) > 0.0)
			nnz++;

	file << "%%MatrixMarket matrix coordinate complex general\n";
	if (!title.empty()) {
		std::string::size_type idx = title.find('\n');
		if (idx == std::string::npos)
			idx = title.size();
		idx = std::min(idx, (std::string::size_type)60);
		file << "% " << title.substr(0,idx) << std::endl;
	}
	file << nr << ' ' << nc << ' ' << nnz << std::endl;

	for (int j=0; j<nc; j++) {
		for (int i=0; i<nr; i++) {
			complex<double> x = cp[(i+j*nr)];
			if (std::abs(x) > 0.0)
				file << i+1 << ' ' << j+1 << ' ' << x.real()
					<< " " << x.imag() << std::endl;
		}
	}

	return rval;
}

bool
MM::
exporter (const string& path, const Matrix* mat) {
/*------------------------------------------------------------------
 * export to a file (path) in Matrix Market format
 * Note: this writes the raw data_ array so it will not work if the
 *       matrix is a function of parameters
 * path can be "cout", "cerr" or a file name; if it is empty
 * it defaults to:
 *    cout       if the number of columns is small enough to fit
 *               within the page width as given by page_width(), or
 *    title.mm   if not
 *------------------------------------------------------------------*/
	Trace trc(1,"MM::exporter(Matrix)");

	trc.dprint("path <",path,"> mat ",mat->summary());

	string title = stringBasename(path);
	int nr = mat->rsize();
	int nc = mat->csize();
	if (mat->is_complex()) {
		return MM::exporter (path, title, mat->celem(), nr, nc);
	} else {
		return MM::exporter (path, title, mat->elem(), nr, nc);
	}
}

static
string
fullpath (string const& path, string const& title) {
// Returns a Matrix Market file name: the title with blanks
// replaced by underscores, slashes (/) replaced with
// underscores, and .mm appended if there is not one already
// If debug() is > 0 the path will be on the Flaps temporary
// directory if it is not a full path
	Trace trc(2,"MM::fullpath");
	string rval;

	trc.dprint("path ",path,", title ",title);

	if (path.empty() || path == "cout" || path == "cerr") {
		if (title.empty())
			rval = "noname";
		else
			rval = title;
	} else {
		rval = path;
	}

	// make sure the path ends in .mm
	if (rsubstr(rval, 3) != ".mm")
		rval += ".mm";

	// Replace blanks with underscores
	rval = replace_char(rval, ' ', '_');
	// If we were called from within a debug statement
	// and path is not a full path, put the file on Ftmp
//	if (debug() > 0 && rval[0] != '/') {
//		rval = vastr(getftmp(), '/', rval);
//	}
	trc.dprint("returning ",rval);
	return rval;
}

// } : namespace MM


// namespace Fig {

bool
Fig::
exporter (const string& file, vector<double>& nodes) {
// write an xfig file of the (x,y) coordinates in "nodes".
// nodes is a vector of (x,y,z) with length 3*nnodes
	ofstream stream(file);
	stream << "#FIG 3.2 produced by Fig::exporter\n";
	stream << "Landscape\n";
	stream << "Center\n";
	stream << "Inches\n";
	stream << "Letter\n";
	stream << "100\n";
	stream << "Single\n";
	stream << "-2\n";
	stream << "1200 2\n";

	int sub_type{3};				// 2)
	int line_style{0};			// 3)
	int thickness{4};				// 4)
	int pen_color{0};				// 5)
	int fill_color{0};			// 6)
	int depth{50};					// 7)
	int pen_style{1};				// 8) unused
	int area_fill{-1};			// 9) no fill
	double style_val{0.0};		// 10)
	int direction{1};				// 11)
	double angle{0.0};			// 12)
	int radius{4};					// 15-16)

	// compute a scale factor so that the largest x, y is 10000
	double xymax{0.0};
	double xmin{std::numeric_limits<double>::max()};
	double ymin{std::numeric_limits<double>::max()};
	size_t nnodes = nodes.size()/3;
	assert(nodes.size()%3 == 0);
	for (size_t i=0; i<nnodes; i++) {
		size_t k = 3*i;
		xmin = std::min(xmin, nodes[k]);
		ymin = std::min(ymin, nodes[k+1]);
		xymax = std::max(nodes[k], xymax);
		xymax = std::max(nodes[k+1], xymax);
	}
	double scale{10000.0/xymax};
	double xoffset{xmin};
	double yoffset{ymin};

	for (size_t i=0; i<nnodes; i++) {
		size_t k = 3*i;
		int x = (int)((nodes[k]-xoffset)*scale);
		int y = (int)((nodes[k+1]-yoffset)*scale);
		// stream << "1 3 0 3 0 0 50 1 -1 0.0 1 0.0 "
		stream << "1 " << sub_type << " " << line_style << " " << thickness
			<< " " << pen_color << " " << fill_color << " " << depth << " "
			<< pen_style << " " << area_fill << " " << style_val << " "
			<< direction << " " << angle << " "
			<< x << " " << y << " " << radius << " " << radius
			<< " " << x << " " << y << " " << x << " " << y << endl;
	}

	return true;
}

// } namespace Fig

// namespace Apf: ASCII plot files (.apf) {

bool
Apf::
exporter (string const& cid, vector<pair<string,string>> const& params,
		vector<double> const& data, const string& path, bool append) {
// write a single curve to a file (path) with data in a (nr,nc) array "data",
// each row corresponding to a parameter with name and description in "params"
	Trace trc(1,"Apf::exporter(vector)");

	size_t nr = params.size();
	size_t nc = data.size()/nr;
	if (data.size()%nr != 0)
		throw runtime_error(vastr("attempt to plot ",data.size()," data with ",nr," parameters"));

	// add a .apf suffix if there is not one
	string filename{path};
	if (rsubstr(path,4) != ".apf")
		filename += ".apf";

	std::ios_base::openmode om = std::ios::trunc;
	if (append)
		om = std::ios::app;
	ofstream ofs(filename, om);
	if (!ofs)
		throw runtime_error(vastr("cannot open ", filename));

	ofs << "$ " << cid << endl;
	trc.dprint("curve ",cid," has ",nr," parameters, ",nc," values");
	// write descriptions, e.g. +1 Velocity (m/s)
	int pno{1};
	for (auto& pj : params)
		ofs << "+" << pno++ << " " << get<1>(pj) << endl;

	// curve id
	ofs << cid << endl;

	// par ids
	string sep;
	for (auto& pj : params) {
		ofs << sep << get<0>(pj);
		sep = " ";
	}
	ofs << endl;

	// write the curve data
	for (size_t j=0; j<nc; j++) {
		const double* dp = &data[j*nr];
		sep = "";
		for (size_t i=0; i<nr; i++) {
			ofs << sep << setprecision(16) << *dp++;
			sep = " ";
		}
		ofs << endl;
	}
	ofs << "*EOF\n";
	return true;
}

bool
Apf::
exporter (vector<Curve*> const& curves, vector<string> const& toplot,
	const string& path, bool append) {
// write some Curves to a file (path)
	Trace trc(1,"Apf::exporter(curves)");

	// lambda to determine if a name is in toplot; always true if toplot is empty
	auto plot = [&](string const& name) {
		if (toplot.empty()) return true;
		if (std::find(toplot.begin(), toplot.end(), name) == toplot.end())
			return false;
		return true;
	};
	// add a .apf suffix if there is not one
	string filename{path};
	if (rsubstr(path,4) != ".apf")
		filename += ".apf";

	std::ios_base::openmode om = std::ios::trunc;
	if (append)
		om = std::ios::app;
	ofstream ofs(filename, om);
	if (!ofs) {
		throw runtime_error(vastr("cannot open ", filename));
	}

	for (auto cp : curves) {
		size_t nval{cp->params.nsolns()};
		// title
		//!! ofs << "$ " << cp->params.desc() << endl;
		ofs << "$ " << cp->cid() << endl;
		trc.dprint("curve ",cp->cid()," has ",cp->params.size()," parameters, ",nval," values");
		int pno{1};
		// parameter descriptions
		for (auto pj : cp->params.pmap()) {
			if (plot(pj.first)) {
				Par* pp{pj.second};
				ofs << "+" << pno++ << " " << pp->desc() << endl;
			}
		}
		// curve id XXX same as desc()?
		//!! ofs << cp->params.desc() << endl;
		ofs << cp->cid() << endl;

		// parameter names
		for (auto pj : cp->params.pmap()) {
			if (plot(pj.first)) {
				Par* pp{pj.second};
				ofs << " " << pp->name;
			}
		}
		ofs << endl;

		// write the curve data (solns), converting from EU to PU
		for (size_t k=0; k<nval; k++) {
			for (auto pj : cp->params.pmap()) {
				if (plot(pj.first)) {
					Par* pp{pj.second};
					ofs << " " << setprecision(16) << pp->solns[k]*pp->conv_factor();
				}
			}
			ofs << endl;
		}
		ofs << "*EOF\n";
	}
	return true;
}

static Curve*
read_curve(Ascii& in, string aid);

vector<Curve*>
Apf::
importer (const string& path, const string& output) {
// Read an ascii plot file (.apf)
	Trace trc(1,"Apf::importer ",path);

	// Attempt to open the apf file; if the path does not end in
	// ".apf" try first with that extension to avoid taking a file
	// like "pk" when there also exists a "pk.apf"
	string fullpath;
	if (rsubstr(path, 4) != ".apf") {
		string apfpath = path + ".apf";
		if (access(apfpath.c_str(), R_OK) != -1) {
			fullpath = apfpath;
		} else if (access(path.c_str(), R_OK) != -1) {
			fullpath = path;
		} else {
			throw runtime_error(vastr("cannot access ", path, " (also tried ",
				path, ".apf)"));
		}
	} else if (access(path.c_str(), R_OK) != -1) {
		fullpath = path;
	} else {
		if (path.empty())
			throw runtime_error("empty path passed to Apf::importer");
		throw runtime_error(vastr("cannot access ", path));
	}
	
	Ascii in(fullpath);

	// the analysis id for these curves is the base file name
	string aid{stringBasename(fullpath)};

	// read each curve
	vector<Curve*> rval;
	Curve* cp{nullptr};
	while((cp = read_curve(in, aid)) != nullptr)
		rval.push_back(cp);
	trc.dprint("returning ",rval.size()," curves");
	return rval;
}

static Curve*
read_curve(Ascii& in, string aid) {
// Read a curve from "in", return a pointer to the curve, giving
// the curve an analysis id "aid"

	// header:
	//   $ run title (up to 3 lines)
	//   +01 first parameter title
	//   +02 second parameter title...
	//   curve_id
	//   parname1 parname2 ...
	Trace trc(2,"read_curve");
	string line;
	string title;
	vector<string> parnames;
	map<int,string> pardesc;
	string curveid;
	while(in(line, "", false)) {
		if (line[0] == '$') {
			string::size_type idx = line.find_first_not_of(" \t", 1);
			if (idx != string::npos)
				title = line.substr(idx);
		} else if (line[0] == '+') {
			size_t idx = line.find_first_of(" ");
			string num{line.substr(1,idx-1)};
			int k;
			str2int(num, k);
			string desc{line.substr(idx+1)};
			pardesc.insert(make_pair(k, desc));
		} else {
			// neither: must be the curve id
			curveid = line;
			break;
		}
	}

	if (line.empty())		// EOF
		return nullptr;
			
	// Parameter names: one line of strings
	in(line, "", false);
	parnames = string2tok(line, " ");
	// create a Par for each
	vector<Par*> params;
	for (size_t i=0; i<parnames.size(); i++) {
		string lhs = parnames[i];
		// see if there was a desc given
		auto j = pardesc.find(i+1);
		if (j != pardesc.end())
			lhs += vastr('(',j->second,')');
		Par* pp = new Par(lhs,"0");
		pp->solns.clear();
		params.push_back(pp);
	}

	// data: one data point per line until *EOF
	while(in(line, "", false)) {
		if (line == "*EOF")
			break;
		vector<string> data = string2tok(line, " ");
		if (data.size() != params.size()) {
			string exc{vastr("line ",in.line_number(), " has ",data.size(),
					" data items, but there are ",params.size()," parameters")};
			throw runtime_error(exc);
		}
		for (size_t i=0; i<params.size(); i++) {
			double x;
			str2double(data[i],x);
			params[i]->solns.push_back(x);
		}
	}

	Curve* rval = new Curve(aid, curveid, "");
	rval->params.clear();
	for (auto& pp : params)
		rval->params.add(pp);

	rval->params.desc(title);
	trc.dprint("curve ",title," has ",rval->params.size()," parameter, ",rval->params.nsolns()," values");

	return rval;
}

bool
Apf::
plot(const string& path, const string& cid, const vector<vector<double>>& xs,
		const vector<string>& names, bool append) {
	Trace trc(1,"Apf::plot");

	trc.dprint("path ",path,", cid ",cid,", ",xs.size()," parameters");
	if (names.size() != xs.size()) {
		string exc{vastr("Apf::plot(",path,"): ",xs.size()," par in xs but ",
				names.size()," names: ",names)};
		throw runtime_error(exc);
	}
	// create a Curve for exporter with a parameter for each xs
	// the analysis id for these curves is the base file name
	string aid{stringBasename(path)};
	Curve* curve{new Curve(aid,cid,"")};
	curve->params.clear();
	for (size_t i=0; i<xs.size(); i++) {
		Par* pp{new Par(names[i], "0")};
		pp->clear_solns();
		for (size_t j=0; j<xs[i].size(); j++) {
			pp->solns.push_back(xs[i][j]);
		}
		curve->params.add(pp);
	}
	vector<Curve*> curves;
	curves.push_back(curve);
	vector<string> toplot;	// plot all
	return Apf::exporter (curves, toplot, path, append);
}

// } namespace Apf

// namespace LTI {
double
aeval(int elidx, vector<double> xb, vector<double> coef, double x) {
// evaluate the spline from simulink for an element of A:
//   rval = c[0]*d*d*d + c[1]*d*d + c[2]*d + c[3]
// where d = x - xb[i], i is row i of coef, c is row i, plane elidx
// of the (nb-1,4,nix) coef array

	int nb = xb.size();
	int nbm = nb - 1;
	double d{0.0};
	int idx{0};
	// check for x out of range: return boundary value
	if (x < xb[0])
		return coef[IJK(0,3,elidx,nbm,4)];
	if (x > xb[nbm]) {
		idx = nbm-1;
		d = xb[nbm] - xb[idx];
	}
	// look for interval containing x
	for (int i=0; i<nbm; i++) {
		if (x >= xb[i] && x <= xb[i+1]) {
			d = x - xb[i];
			idx = i;
			break;
		}
	}
	// double d = xb[nbm] - xb[nbm-1];
	//!! ap[elem[k].index(ns)] = ((coef[IJK(nbm-1,0,k,nbm,4)]*d +
	//!! 							       coef[IJK(nbm-1,1,k,nbm,4)])*d +
	//!! 						          coef[IJK(nbm-1,2,k,nbm,4)])*d +
	//!! 					             coef[IJK(nbm-1,3,k,nbm,4)];
	double rval = ((coef[IJK(idx,0,elidx,nbm,4)]*d +
				       coef[IJK(idx,1,elidx,nbm,4)])*d +
				       coef[IJK(idx,2,elidx,nbm,4)])*d +
				       coef[IJK(idx,3,elidx,nbm,4)];
	return rval;
}

bool
LTI::
importer (string const& path, Matrix& A, Matrix& B, Matrix& C, Matrix& D,
	vector<Itd>& Atd, vector<double>& itd, vector<double>& otd,
	vector<int>& Sidx, double rho) {
	Trace trc(1,"LTI::importer");

	// open the input file: throw exception if not available
	Ascii in(path);

	// a lambda to convert a string to a double, catching exceptions
	auto Stod = [&](string& s) {
		double rval;
		if (s.empty())
			throw runtime_error("expecting a double, got EOF");
		try {
			rval = stod(s);
		} catch (std::exception& t) {
			throw runtime_error(vastr("expecting a double on line ",
				in.line_number()," got \"",s,"\""));
		}
		return rval;
	};

	 // First non-comment line: 6 integers
	 //   ns:   order of A
	 //   ni:   number of input channels
	 //   no:   number of ouput channels
	 //   nb:   number of interpolation breakpoints
	 //         for the elements of A that vary
	 //   nix:  number of elements of A that vary
	 //   nin:  number of elements of A that are constant (this is redundant)
	string line;
	in(line,"#",false);
	vector<string> tok = string2tok(line, " ");
	if (tok.size() != 6)
		throw runtime_error(vastr("expecting 6 integers on line ",in.line_number()));

	int ns = stoi(tok[0]);
	int ni = stoi(tok[1]);
	int no = stoi(tok[2]);
	int nb = stoi(tok[3]);	// number of A breakpoints
	int nix = stoi(tok[4]);
	int nin = stoi(tok[5]);
	trc.dprint("ns ",ns,", ni ",ni,", no ",no,", nb ",nb,", nix ",nix,", nin ",nin);

	// Second (non-comment) line: interpolation variable name
	in(line, "#",false);
	tok = string2tok(line, " ");
	if (tok.size() != 1)
		throw runtime_error(vastr("expecting interpolation variable name on line ",
			in.line_number()));

	string paramA = tok[0];	// output
	// check that the parameter is in gpset
	// Replace paramA with the real name from gpset: allow different
	// case and longer name, e.g. "Altitude" instead of "alt"
	Par* parA{nullptr};
	// loop through all parameters in gpset, search for paramA
	// XXX need a method in gpset to do this?
	for (auto& elem : gpset::get().pmap()) {
		if (compare(paramA, elem.first)) {
			paramA = elem.first;
			parA = elem.second;
			break;
		}
	}
	if (parA == nullptr)
		throw runtime_error(vastr("the A-matrix interpolation parameter (",
			paramA,") from the LTI input file (",path,") has not been defined"));

	// nb double breakpoints for interpolation of A wrt parA
	// If parA has a conversion factor, don't apply it
	vector<double> xb;
	for (int i=0; i<nb; i++) {
		in(line, "#",false);
		xb.push_back(Stod(line));	// output
	}
	trc.dprint("interp A wrt ",paramA," breakpoints:");
	for (auto& xi : xb)
		cerr << setprecision(16) << xi << endl;
	
	// nix pairs of integers: row,col of elements of A that vary
	int row, col;
	vector<Elem> elem;
	//!! vector<pair<int,int>> aelem;
	for (int i=0; i<nix; i++) {
		in(line, "#",false);
		tok = string2tok(line, " ");
		if (tok.size() != 2)
			throw runtime_error(vastr("expected 2 integers, line ",
				in.line_number(),", got \"", line, "\""));

		str2int(tok[0],row);
		str2int(tok[1],col);
		//!! aelem.push_back({row-1,col-1});	// zero-based XXX just output interp coeff
		size_t r = row - 1; // 0b
		size_t c = col - 1; // 0b
		elem.push_back({r,c});
	}

	// indices of constant terms of A - ignore
	for (int i=0; i<nin; i++)
		in(line, "#",false);

	// Input, output time delays
	for (int i=0; i<ni; i++) {
		in(line, "#",false);
		itd.push_back(Stod(line));	// output
	}
	for (int i=0; i<no; i++) {
		in(line, "#",false);
		otd.push_back(Stod(line));	// output
	}

	// exponents for S: come as floats, convert to ints
	for (int i=0; i<ni; i++) {
		in(line, "#",false);
		double si = Stod(line);
		if (is_equal(si, 0.0, 5))
			Sidx.push_back(0);
		else if (is_equal(si, 1.0, 5))
			Sidx.push_back(1);
		else if (is_equal(si, 2.0, 5))
			Sidx.push_back(2);
		else
			throw runtime_error(vastr("illegal S-matrix exponent: ", si));
	}

	// A matrix: (ns,ns) complex by rows
	A = Matrix("A","LTI A",ns,ns,true);	// output
	complex<double>* ap = A.celem();
	for (int i=0; i<ns; i++) {
		for (int j=0; j<ns; j++) {
			in(line, "#",false);
			ap[IJ(i,j,ns)] = Stod(line);
		}
	}
	// B matrix: (ns,ni) double by rows
	B = Matrix("B","LTI B",ns,ni,false);	// output
	double* bp = B.elem();
	for (int i=0; i<ns; i++) {
		for (int j=0; j<ni; j++) {
			in(line, "#",false);
			bp[IJ(i,j,ns)] = Stod(line);
		}
	}
	// C matrix: (no,ns) double by rows
	C = Matrix("C","LTI C",no,ns,false);	// output
	double* cp = C.elem();
	for (int i=0; i<no; i++) {
		for (int j=0; j<ns; j++) {
			in(line, "#",false);
			//!! C.push_back(Stod(line));
			cp[IJ(i,j,no)] = Stod(line);
		}
	}
	// D matrix: (no,ni) double by rows
	D = Matrix("D","LTI D",no,ni,false);	// output
	double* dp = D.elem();
	for (int i=0; i<no; i++) {
		for (int j=0; j<ni; j++) {
			in(line, "#",false);
			dp[IJ(i,j,no)] = Stod(line);
		}
	}

	// cubic-spline coeff for the variable terms of A:
	// a (nb-1, 4, nix) array
	if (nb > 1 && nix > 0) {
		int nbm = nb - 1;
		vector<double> coef(nbm*4*nix);
		for (int k=0; k<nix; k++) {
			for (int i=0; i<nbm; i++) {
				for (int j=0; j<4; j++) {
					in(line, "#",false);
					coef[IJK(i,j,k,nbm,4)] = Stod(line);
				}
			}
		}
		// create a plotfile of element elidx for comparison with
		// the smoothed version below; av is (2,npt): x is first row, element second
		int npt{101};
		double del = (xb[nbm] - xb[0])/(npt-1);
		vector<double> av(2*npt);
		int elidx{0};
		for (int j=0; j<npt; j++) {
			double x = xb[0] + j*del;
			av[IJ(0,j,2)] = x;
			av[IJ(1,j,2)] = aeval(elidx, xb, coef, x);
		}
		//!! vector<tuple<string,string>> params{{paramA,"x"},
		//!! 	{vastr("A",aelem[elidx].first,aelem[elidx].second),"y"}};
		vector<pair<string,string>> params{{paramA,"x"},
			{"real", "x"}};
			//!! {vastr("A",elem[elidx]),"x"}};
		// av is (2,npt): first row is x, second is element
		Apf::exporter(vastr("A",elem[elidx]), params, av, "simulink.apf", false);
		
		// smooth the splines: create nb copies of A, replace the
		// interpolated terms, create an IntPz
		vector<Matrix*> matrices;
		matrices.push_back(&A);
		for (int i=0; i<nbm; i++)
			matrices.push_back(A.clone());
		// interpolate the variable terms
		for (int i=0; i<nbm; i++) {
			double* ap = matrices[i]->elem();
			for (int k=0; k<nix; k++)
				//!! ap[IJ(aelem[k].first,aelem[k].second,ns)] = coef[IJK(i,3,k,nbm,4)];
				ap[elem[k].index(ns)] = coef[IJK(i,3,k,nbm,4)];
		}
		// the last breakpoint must be interpolated
		double* ap = matrices[nbm]->elem();
		//!! double d = xb[nbm] - xb[nbm-1];
		for (int k=0; k<nix; k++)
			//!! ap[elem[k].index(ns)] = ((coef[IJK(nbm-1,0,k,nbm,4)]*d +
			//!! 									coef[IJK(nbm-1,1,k,nbm,4)])*d +
			//!! 									coef[IJK(nbm-1,2,k,nbm,4)])*d +
			//!! 									coef[IJK(nbm-1,3,k,nbm,4)];
			//!! ap[IJ(aelem[k].first,aelem[k].second,ns)] = aeval(k, xb, coef, xb[nbm]);
			ap[elem[k].index(ns)] = aeval(k,xb,coef,xb[nbm]);
		// give each A a Pz_const with the breakpoint value
		vector<Par*> par;
		par.push_back(parA);
		if (xb.size() != matrices.size())
			throw runtime_error(vastr(xb.size()," breakpoints, but ",matrices.size(),
				" matrices (LTI::importer)"));
		for (int i=0; i<nb; i++) {
			par[0]->value(xb[i]);
			matrices[i]->pz.push_back(new Pz_const(par, ns, ns));
		}
		// create an IntPz smoothed with rho
		IntPz* intpz = new IntPz(matrices, 0, rho);
		A.pz.clear();
		A.pz.push_back(intpz);
		vector<Elem> elem{{0,0}};	// XXX change aelem back to Elem?
		A.plot_apf(gpset::get(), "smoothed.apf",elem, 101, paramA);
	}

	// internal time delays (optional)
	int nitd{0};
	if (in(line, "#",false))
		nitd = stoi(line);

	for(int i=0; i<nitd; i++) {
		Itd itd;
		int ninda;
		in(line, "#",false);
		ninda = stoi(line);
		for (int j=0; j<ninda; j++) {
			in(line, "#",false);
			tok = string2tok(line, " ");
			if (tok.size() != 2)
				throw runtime_error(vastr("expected 2 integers, line ",in.line_number(),
					", got \"", line, '\"'));
			row = stoi(tok[0]);
			col = stoi(tok[1]);
			itd.elem.push_back(Elem(row-1,col-1));  // zero-based
		}
		Atd.push_back(itd);	// output
	}
	for(int i=0; i<nitd; i++) {
		in(line, "#",false);
		Atd[i].deltat = Stod(line);
	}
	return true;
}	// importer
// }	namespace LTI

#ifdef MAIN
#undef MAIN

#include "main.h"	// set debug, etc

bool
test_lti(const string& file) {
	Matrix A,B,C,D;
	vector<Itd> Atd;
	vector<double> itd, otd;
	vector<int> Sidx;
	double rho{2.0};
	LTI::importer (file, A, B, C, D, Atd, itd, otd, Sidx, rho);
	return true;
}

int
main (int argc, char **argv) {

	if (argc == 1) {
		cerr << "Usage: " << argv[0] << " file\n";
		cerr << "        the file extension determines the format\n";
		return 1;
	}

	string file(argv[1]);
	string basename = stringBasename(file);
	string::size_type idx = file.rfind(".");
	if (idx == string::npos) {
		cerr << "the input file name must have one of the extensions mat, mm, uf\n";
		return 1;
	}

	string exten = file.substr(idx);
	vector<Matrix*> ml;
	try {
		if (exten == ".mat") {
			ml = Matlab::importer(file);
			for (auto mp : ml)
				mp->mid(vastr(mp->mid(), "2"));
			Matlab::exporter(vastr(basename,"2.mat"), ml);
		} else if (exten == ".op4") {
			ml = Output4::importer(file);
		} else if (exten == ".mm") {
			ml.push_back(MM::importer(file));
		} else if (exten == ".uf") {
			ml = UF::importer (argv[1]);
		} else if (exten == ".lti") {
			test_lti (file);
		}
		// for visualizing with matview:
		for (auto mp : ml)
			MM::exporter(vastr(mp->mid(), ".mm"), mp);
	} catch (runtime_error& s) {
		std::cerr << "Error: " << s << std::endl;
		exit(1);
	}
}
#endif // MAIN
