//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

// Trace support is enabled in the build process by uncommenting
// CFLAG=-DTRACE in froot/Profile

// Class Trace: support for trace printing when debug turned on; debug is
// turned on to level "n" by:
// 1) settings{d=n} in a flaps script
// 2) Trace::debug(n) in a program
// 3) DEBUG=n in the environment
// Printing is done mainly with Trace::dprint* functions when the trace
// level (set with the Trace constructor) is >= debug level
//
#ifndef TRACE_H
#define TRACE_H

#include <cassert>
#include <complex>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <ctime>
#include <vector>

#include "Ad.h"
#include "message.h"
#include "text.h"
#include "vastr.h"

#ifdef TRACE
#  define T_(x) x
#else // NO TRACE
#  define T_(x)
#endif // TRACE

// XXX move Offset stuff to class Trace?
void
incrementOffset();

void
decrementOffset();

std::string DOffset(void);
void resetDOffset ();

// Functions for printing arrays not necessarily under debug
// XXX need new names?
namespace Dprint {
void
dprintm (int nr,int nc,int lda,const std::vector<Ad>& a,const std::string& title);
void
dprintm (int nr,int nc,int lda,const std::vector<std::complex<Ad>>& a,
		const std::string& title);
void
dprintm (int nr, int nc, int lda, const int* a, const std::string& title);
void
dprintm (int nr, int nc, int lda, const double* a, const std::string& title);
void
dprintm (int nr, int nc, int lda, const std::complex<double>* a, const std::string& title);
void
correction(int n, const double* x, const double* h, std::string const& title);
} // namespace Dprint

class Trace {
private:
	std::string name;
	int Level;
	bool add_newline{true};
	int olddebug{-1};
public:

	template<typename... Types>
	Trace(int lev, const std::string fcn, const Types&... args) : Level(lev) {
		// Trace constructor: print Enter message, increment offset
		if (debug() >= Level) {
			name = fcn + vastr(args...);
			*logfile() << DOffset() << "Enter " << name << " {\n";
			incrementOffset();
		}
	}
	template<typename... Types>
	Trace(int lev, int dbg, const std::string fcn) : name(fcn), Level(lev) {
		olddebug = debug();
		debug(dbg);
		if (dbg >= Level) {
			*logfile() << DOffset() << "Enter " << fcn << " {\n";
			incrementOffset();
		}
	}
	~Trace() {
		if (debug() >= Level) {
			decrementOffset();
			*logfile() << DOffset() << "Exit " << name << " }\n";
			if (olddebug != -1)
				debug(olddebug);
		}
	}

	//!! static std::ostream* logfile(const std::string& path = "");
	static std::ostream* logfile(const std::string& path="") {
	// open a new logfile or return a pointer to the current one
		static std::ostream* thestream{&std::cerr};
		std::ostream* rval = thestream;
		if (!path.empty())
			thestream = new std::ofstream(path);
		return rval;
	}

	// get/set global debug level
	static int debug(int d=-1);

	bool operator()() {return Level <= debug(); }

	template<typename T>
	int
	column_width(int nc, T* a) {
		std::string wid = strArray(1,nc,1,a);
		return (int)wid.size();
	}

	template<typename T, typename... Types>
	void
	dprint(const T& firstArg, const Types&... args) {
	// debug print with trailing newline
		if (debug() < Level)
			return;
		std::ostream* lf = logfile();
		// set float format to default
		lf->unsetf(std::ios::floatfield);
		lf->precision(6);
		if (add_newline)
			*lf << DOffset();
		*lf << vastr(firstArg, args...);
		*lf << "\n";
		lf->flush();
		add_newline = true;
	}

	template<typename T, typename... Types>
	void
	dprintn(const T& firstArg, const Types&... args) {
	// debug print without trailing newline
		std::ostream* lf = logfile();
		lf->unsetf(std::ios::floatfield);
		lf->precision(6);
		if (debug() >= Level) {
			if (add_newline)
				*lf << DOffset();
			*lf << firstArg;
			*lf << vastr(args...);
			lf->flush();
			add_newline = false;
		}
	}

	// dprintm: print an (nr,nc) matrix stored in an (lda,nc)
	// column-major array
#ifdef NEVER // duplicate
	template<typename... Types>
	void
	dprintm (int nr,int nc,int lda,const std::vector<Ad>& a,const Types&... args) {
	// vector<Ad> print: only the values, not derivatives
		if (debug() < Level)
			return;
		std::string title = vastr(args...);
		Dprint::dprintm(nr,nc,lda,a.vec(),title);
	}
#endif // NEVER // duplicate

	template<typename... Types>
	void
	dprintm (int nr,int nc,int lda,const std::vector<Ad>& a,const Types&... args) {
	// vector<Ad> print: only the values, not derivatives
		if (debug() < Level)
			return;
		std::string title = vastr(args...);
		Dprint::dprintm(nr,nc,lda,a,title);
	}

#ifdef NEVER // duplicate
	template<typename... Types>
	void
	dprintm (int nr,int nc,int lda,const std::vector<std::complex<Ad>>& a,const Types&... args) {
	// CAdvector print: only the values, not derivatives
		if (debug() < Level)
			return;
		std::string title = vastr(args...);
		Dprint::dprintm(nr,nc,lda,a.vec(),title);
	}
#endif // NEVER // duplicate

	template<typename... Types>
	void
	dprintm (int nr,int nc,int lda,const std::vector<std::complex<Ad>>& a,const Types&... args) {
	// CAdvector print: only the values, not derivatives
		if (debug() < Level)
			return;
		std::string title = vastr(args...);
		Dprint::dprintm(nr,nc,lda,a,title);
	}


	template<typename T, typename... Types>
	void
	dprintm (int nr, int nc, int lda, const T* a, const Types&... args) {
		if (debug() < Level)
			return;
		std::string title = vastr(args...);
		Dprint::dprintm(nr,nc,lda,a,title);
	}

	template<typename T, typename... Types>
	void
	dprintm (int nr, int nc, int lda, const std::vector<T>& a, const Types&... args) {
		if (debug() >= Level) {
			std::string title = vastr(args...);
			Dprint::dprintm(nr,nc,lda,&a[0],title);
		}
	}

	template<typename T, typename... Types>
	void
	dprintv (const std::vector<T>& a, const Types&... args) {
		if (debug() < Level)
			return;
		std::ostream* lf = logfile();
		*lf << DOffset();
		*lf << vastr(args...);
		if (a.empty()) {
			*lf << "  empty\n";
			return;
		}
		// print items ipl/line to keep line lengths < 80
		std::ostringstream os;
		os << a[0];
		int ipl = 80/(os.str().size()+2);
		if (ipl < 1) ipl = 1;
		*lf << "(" << a.size() << "): ";
		// if ipl>a.size put the vector on the same line as text
		if (ipl >= (int)a.size()) {
			for (auto& ai : a)
				*lf << "  " << ai;
		} else {
			for (size_t i=0; i<a.size(); i++) {
				if (i%ipl == 0) *lf << std::endl << DOffset();
				*lf << a[i];
				if (i != a.size()-1) *lf << ", ";
			}
		}
		*lf << std::endl;
	}

	template<typename T, typename... Types>
	void
	dprintvv (const std::vector<T>& a, const std::vector<T>& b,
			const Types&... args) {
		// print 2 vectors side-by-side
		if (debug() < Level)
			return;
		std::ostream* lf = logfile();
		size_t n = a.size();
		if (b.size() != n) {
			throw std::runtime_error(vastr("dprintvv called with vectors of length ",
						n, " and ", b.size()));
		}
		*lf << DOffset();
		*lf << vastr(args...);
		*lf << ": {\n";
		for (size_t i=0; i<n; i++) {
			*lf << DOffset() << "   " <<
				a[i] << "  " << b[i] << std::endl;
		}
		*lf << DOffset() << "}\n";
	}

	template<typename T, typename... Types>
	void
	dprintvvv (const std::vector<T>& a, const std::vector<T>& b,
			const std::vector<T>& c, const Types&... args) {
		// print 3 vectors side-by-side
		if (debug() < Level)
			return;
		std::ostream* lf = logfile();
		size_t n = a.size();
		if (b.size() != n || c.size() != n) {
			throw std::runtime_error(vastr("dprintvv called with vectors of length ",
						n, ", ", b.size(),", and ",c.size()));
		}
		*lf << DOffset();
		*lf << vastr(args...);
		*lf << ": {\n";
		for (size_t i=0; i<n; i++) {
			*lf << DOffset() << "   " <<
				a[i] << "  " << b[i] << "  " << c[i] << std::endl;
		}
		*lf << DOffset() << "}\n";
	}

	template<typename T, typename... Types>
	void
	dprintvp (const std::vector<T>& a, const Types&... args) {
		// print a vector of pointers to a type that has operator<<
		if (debug() >= Level) {
			std::ostream* lf = logfile();
			*lf << DOffset();
			*lf << vastr(args...);
			if (a.size() < 6) {
				*lf << ": ";
				for (auto& ai : a)
					*lf << "  " << *ai;
				*lf << std::endl;
			} else {
				*lf << ":\n";
				for (auto& ai : a)
					*lf << DOffset() << "   " << *ai << std::endl;
			}
		}
	}

	template<typename... Types>
	void
	correction (int n, const double* x, const double* h, const Types&... args) {
	// Special function for continuation correctors: print x, x+h, and x/h
		if (debug() < Level)
			return;
		std::string title = vastr(args...);
		Dprint::correction(n,x,h,title);
	}

	// return the name of the function being Trace'd
	std::string fcnname() { return name; }
};  // class Trace

class ObjTrace {
private:
	static int count;
	const char* name;
public:
	int ct;

	ObjTrace(char const* s): ct(++count) {
		if (Trace::debug()>0) {
			name = s;
			// dprint(1,name,' ',ct," constructed");
			std::cerr << (name?name:"Object") << ' ' << ct << " constructed\n";
		}
	}
	~ObjTrace() {
		if (Trace::debug()>0)
			// dprint(1,name,' ',ct," destroyed");
			std::cerr << name << ' ' << ct << " destroyed\n";
	}
	ObjTrace& operator=(const ObjTrace&) {
		return *this;
	}
};

#endif  // TRACE_H
