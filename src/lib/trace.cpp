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
 * trace: routines for debug tracing
 *------------------------------------------------------------------*/
#include <sys/stat.h>
#include <cstring>

#include "trace.h"
#include "exim.h"
#include "fio.h"
#include "matrix.h"

using namespace std;

int
Trace::
debug(int d) {
// set/get the global debug level
	static int thelevel{0};
	int oldlevel = thelevel;
	if (d >= 0)
		thelevel = d;
	return oldlevel;
}

// XXX not thread-safe
#define MAXLEVEL 200
static int Offlvl = 0;
//!! static char Offset[3*(MAXLEVEL+1)] = "";
void
incrementOffset() {
   Offlvl++;
	//!! memset (Offset, ' ', 3*Offlvl);
   //!! Offset[3*Offlvl] = '\0';
}
void
decrementOffset() {
   Offlvl--;
	if (Offlvl < 0) {
		cerr << "level less than zero\n";
		Offlvl = 0;
	}
   //!! Offset[3*Offlvl] = '\0';
}

void
resetDOffset () {
	Offlvl = 0;
   //!! Offset[0] = '\0';
}

string
DOffset() {
	//!! return Offset;
	return string(Offlvl, ' ');
}

// Trace implementation

namespace Dprint {
void
dprintm (int nr,int nc,int lda,const vector<Ad>& a,const string& title) {
// vector<Ad> print: print the values and derivatives
	std::vector<double> av(nr*nc);
	vector<string> toprint{""};
	vector<string> const& adparnames = Ad::adnames();
	for (auto& nm : adparnames)
		toprint.push_back(nm);
	for (auto& nm : toprint) {
		extract(a, nm, &av[0]);
		std::string width = strArray(nc,1,1,&av[0]);
		string ti{vastr(title," partial wrt ",nm)};
		if (nm.empty())
			ti = title;
		if ((int)width.size() > page_width()) {
			std::string path = vastr(getftmp(),'/',ti,".mm");
			MM::exporter(path, ti, &av[0], nr, nc);
		} else {
			std::ostream* lf = Trace::logfile();
			*lf << ti << endl << strArray (nr, nc, lda, &av[0]);
		}
	}
}

void
dprintm (int nr,int nc,int lda,const vector<complex<Ad>>& a,const string& title) {
// vector<complex<Ad>> print: print the values and derivatives
	std::vector<std::complex<double> > av(nr*nc);
	vector<string> toprint{""};
	vector<string> const& adparnames = Ad::adnames();
	for (auto& nm : adparnames)
		toprint.push_back(nm);
	for (auto& nm : toprint) {
		extract(a, nm, &av[0]);
		std::string width = strArray(nc,1,1,&av[0]);
		string ti{vastr(title," partial wrt ",nm)};
		if (nm.empty())
			ti = title;
		if ((int)width.size() > page_width()) {
			std::string path = vastr(getftmp(),'/',ti,".mm");
			MM::exporter(path, ti, &av[0], nr, nc);
		} else {
			std::ostream* lf = Trace::logfile();
			*lf << ti << endl << strArray (nr, nc, lda, &av[0]);
		}
	}
}

void
dprintm (int nr, int nc, int lda, const int* a, const string& title) {
	std::ostream* lf = Trace::logfile();
	*lf << title << endl << strArray (nr, nc, lda, a);
}

void
dprintm (int nr, int nc, int lda, const double* a, const string& title) {
	std::string width = strArray(nc,1,1,a);
	if ((int)width.size() > page_width()) {
		std::string path = vastr(getftmp(),'/',title,".mm");
		MM::exporter(path, title, a, nr, nc);
	} else {
	 	std::ostream* lf = Trace::logfile();
	 	*lf << title << endl << strArray (nr, nc, lda, a);
	}
}

void
dprintm (int nr, int nc, int lda, const complex<double>* a, const string& title) {
	std::string width = strArray(nc,1,1,a);
	if ((int)width.size() > page_width()) {
		std::string path = vastr(getftmp(),'/',title,".mm");
		MM::exporter(path, title, a, nr, nc);
	} else {
		std::ostream* lf = Trace::logfile();
		*lf << title << endl << strArray (nr, nc, lda, a);
	}
}

void
correction(int n, const double* x, const double* h, string const& title) {
// print on Trace::logfile() info about a correction
	std::ostream* lf = Trace::logfile();
	string offset(DOffset());
	double eps = std::numeric_limits<double>::epsilon();

	*lf << offset;
	*lf << title << " {\n";
	*lf << offset;
	*lf << "        x^i                h                  x^i+1           h/x\n";
	*lf << offset;
	*lf << "-------------------------------------------------------------------\n";

	for (int i=0; i<n; i++) {
		double ximh = x[i] - h[i];
		double t = std::max(std::abs(ximh), std::abs(x[i]));
		double del = 0.0;
		if (t > eps)
			del = safe_divide(std::abs(h[i]), t);
		*lf << offset << setw(3) << i+1;
		*lf << setw(16) << setprecision(8) << fixed
				<< x[i];
		*lf << setw(16) << setprecision(8) << fixed
				<< -h[i];
		*lf << setw(16) << setprecision(8) << fixed
				<< ximh;
		*lf << setw(16) << setprecision(8) << fixed
				<< del << endl;
	}
	*lf << "}\n";
}

} // namespace Dprint


#ifdef MAIN

#include "message.h"

class A {
	ObjTrace T;
	int a;
public:
	A() : T("A"), a(0) {}
	A(int y) : T("A"), a(y) {}
};

class B : public A {
	ObjTrace T;
public:
	int b;
	B(int x=-1) : T("B"), b(x) {}
};

ostream& operator<<(ostream& s, const B& t) {
	return s << "B:" << t.b;
}

string
fcn() {
	Debug dbg("fcn");
	Trace trc(2,"fcn", "1st arg", string("2nd arg"));
	B b(33);
	B c(34);
	A d(35);  // no operator<<
	// this line won't compile: d has no operator<<:
	// trc.dprintn("b = ",b, ", c = ",c, "obj with no operator<<:",d);
	trc.dprintn("b = ",b, ", c = ",c);
	return "fcn called";
}

int
main(int argc, char* argv[]) {
	Trace trc(1,__FILE__," argc ",argc);
	// test dprintv
	vector<uint64_t> ia = {
   70346511260311519, 70346545890586559, 70350978565267423, 70918291110453181, 70918291110485981, 
   70918291114913727, 70918291664396255, 70918360924942271, 71481275694161727, 71481275694182207, 
   71481275696521183, 71481275979657151, 71481311146794975, 71485743821483967, 72053126181150687 };
	trc.dprintv(ia, "a vector of uint64_t, length ",ia.size());
	vector<int> ib(4, 4);
	trc.dprintv(ib, "a vector of ints, length ",ib.size());

	// cout << "test of sprint: " << sprint("an int ",ia[0]," and float ",1.2)
	// 	<< endl;

	vector<double> d = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
	int m{3};
	int n{2};
	trc.dprintm(3,2,3,d, m," by ",n," double matrix");
}
#endif
