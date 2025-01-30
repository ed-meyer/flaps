//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

#include <iostream>
#include <sys/times.h>
#include <unistd.h>

#include "tyme.h"

using namespace std;

Tyme::
Tyme (const string& nm) : name(nm) {
	wcstart = ::times(&start);
}

Tyme::
~Tyme() {
	struct tms fini;
	clock_t wcticks = times(&fini) - wcstart;
	double clk_tck = (double)sysconf(_SC_CLK_TCK);  // ticks/sec APUE pg 36
	double wctime = wcticks/clk_tck;
	double usrtime = (double)(fini.tms_utime - start.tms_utime)/clk_tck;
	double systime = (double)(fini.tms_stime - start.tms_stime)/clk_tck;
	double childut = (double)(fini.tms_cutime - start.tms_cutime)/clk_tck;
	double childst = (double)(fini.tms_cstime - start.tms_cstime)/clk_tck;
	// cerr << name << " user time " << usrtime << ", system " << systime << endl;
	cerr << name << " times " << wctime << "(" << usrtime << ", " << systime << ")("
		<< childut << ", " << childst << ")\n";
}

#ifdef MAIN

#include <cmath>
#include "trace.h"
#include "vastr.h"

int
main() {
	T_(Trace trc(1,"tyme");)
	constexpr int niter{100000};
	constexpr int nsin{1000};
	constexpr double pi{atan(1.0)};
	constexpr double del{2.0*pi/nsin};
	Tyme tyme(vastr(niter*nsin," calls to sin()"));
	for (int j=0; j<niter; j++) {
		for (int i=0; i<nsin; i++) {
			double y = i*del;
			y = sin(y);
		}
	}
}
#endif // MAIN
