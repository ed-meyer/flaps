//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

#include <fftw3.h>
#include <vector>
#include <string>

#include "trace.h"
#include "hrpts.h"
using namespace std;

double
samples_per_stride(int first, int last, const double* accels);


vector<int>
hrpts (int first, int last, const vector<double>& accels) {
	int n = (int)accels.size();
	return hrpts(first, last, n, &accels[0]);
}

vector<int>
hrpts (int first, int last, int n, const double* accels) {
// Find the beginning and end of each stride in accels
// Assumptions:
//    accels[first] is a peak
//    first, last are 0-based
//    last is < n-1
	Trace trc(0,"hrpts");
	vector<int> rval;

	trc.dprint("first ",first,", last ",last," ",n," accels");

	if (first < 0) {
		string exc{vastr("hr first point < 0: ",first)};
		throw runtime_error(exc);
	}
	if (last < 0 || last >= n) {
		string exc{vastr("hr last point (",last,") out of the range 0-",n-1)};
		throw runtime_error(exc);
	}

	// sps is the average samples per stride: use it as a starting
	// guess for where the end of the first stride is
	int sps = samples_per_stride(first, last, accels);
	trc.dprint("initial samples per stride = ",sps);

	// double t{times[first]};
	// double tlast{times[last]};
	// trc.dprint("time interval: ",t,':',tlast);

	double big{std::numeric_limits<double>::max()};

	int begin{first};    // beginning of a stride
	rval.push_back(begin);
	int del = sps/10; // # samples in search interval x2
	while(true) {
		// search for a peak (end of a stride) between idx1 and idx2
		// ok if last is in the interval, quit if idx1 is past it
		int idx1 = begin + sps - del;
		int idx2 = begin + sps + del;
		if (idx2 > n-1)
			idx2 = n-1;   // ie last accel
		vector<int> cand;
		if (idx1 >= last) {
			trc.dprint("quitting: idx1 > last");
			break;
		}
		// get possible peaks in the interval; watch out for consecutive
		// points equal
		trc.dprint("search for peak in ",idx1,':',idx2);
		for (int i=idx1+1; i<idx2; i++) {
			double a1 = accels[i-1];
			double a = accels[i];
			double a2 = accels[i+1];
			int skip{0};
			if (a1 == a) {
				for (int j=i-2; j>=0; j--) {
					a1 = accels[j];
					if (a1 != a) {
						trc.dprint("moved index back from ",i-1," to ",j);
						break;
					}
				}
			}
			if (a == a2) {
				for (int j=i+2; j<=last; j++) {
					a2 = accels[j];
					if (a2 != a) {
						skip = j - (i+1);
						trc.dprint("moved index ahead from ",i+1," to ",j," skipping ",skip);
						break;
					}
				}
			}
			trc.dprint(a1," < ",a," > ",a2);
			if (a1 < a && a2 < a) {
				cand.push_back(i);
			}
			i += skip;
		}
		if (cand.empty()) {
			// rval.clear();
			trc.dprintv(rval,"no candidate peaks found, returning ");
			return rval;
		}
		// the largest peak is the end of the stride
		double peak{-big};
		int end;
		for (auto idx : cand) {
			if (accels[idx] > peak) {
				peak = accels[idx];
				end = idx;
			}
		}
		rval.push_back(end);
		trc.dprint("added ",end," to rval");
		// new sps...
		sps = end - begin;
		trc.dprint("new samples per stride = ",sps);
		// ... and the start of the next stride
		begin = end;
	}
	
	trc.dprintv(rval, "returning");
	return rval;
}

double
samples_per_stride(int first, int last, const double* accels) {
// Compute the average number of samples per stride over the
// interval between samples "first" and "last"
	Trace trc(0,"samples_per_stride");
	int n = last - first + 1;

	double* in = (double*)fftw_malloc(n*sizeof(double));
	int j{0};
	for (int i=first; i<=last; i++)
		in[j++] = accels[i];
	// vector<double> in(&accels[first], &accels[last+1]);
	// actually only need n/2+1
	int nc = n/2 + 1;
	if (n%2 == 1)
		nc++;
	trc.dprint(n," points, ",nc," complex coef");

	fftw_complex* out = (fftw_complex*)fftw_malloc(nc*sizeof(fftw_complex));
	// vector<complex<double> > out(nc, complex<double>(0.0));

	fftw_plan p = fftw_plan_dft_r2c_1d(n, in, out, FFTW_ESTIMATE);
	fftw_execute(p);

	// find the largest abs coef - ignore the DC (c0)
	double cmax{0.0};
	int imax{0};
	for (int i=1; i<nc-1; i++) {
		double absi = out[i][0]*out[i][0] + out[i][1]*out[i][1];
		trc.dprint("coef ",i," = ",absi);
		if (absi > cmax) {
			cmax = absi;
			imax = i;
		}
	}
	trc.dprint("largest coef = ",cmax,", no. ",imax);
	// The motion is assumed to be biphasic so the number of strides
	// is imax/2; imax may be odd in which case the number of strides
	// is (imax-1)/2 + 0.5.
	double nstride = double(imax)/2.0;

	// samples per stride (average) is the total number of samples/nstride
	int rval = (int)(double(n)/nstride);
	trc.dprint(nstride, " strides, samples per stride = ",rval);

	// cleanup
	fftw_destroy_plan(p);
	fftw_free(in);
	fftw_free(out);
	return rval;
}


#ifdef MAIN

vector<double>
readv(char* file) {
	ifstream ifs(file);
	string line;
	vector<double> rval;
	while(ifs.good()) {
		getline(ifs, line);
		if (!line.empty())
			rval.push_back(stod(line));
	}
	return rval;
}

int
main(int argc, char** argv) {
// Usage:
//   hrpts [first last acc_file]
// if called with no arguments a simple cosine wave is
// analyzed; otherwise it should have 3 arguments:
//   first, last   sample numbers delimiting the region of interest
//   acc_file      text file containing the accelerations, one per line
	vector<double> accels;
	vector<double> times;
	int first{0};
	int last{0};
	double pi{atan(1.0)*4.0};

	// hrpts first last acc_file
	if (argc == 4) {
		first = stoi(string(argv[1])) - 1;
		last = stoi(string(argv[2])) - 1;
		accels = readv(argv[3]);
		vector<int> pts = hrpts(first, last, accels);
	} else {
		int npts{65};
		double period{4.0*pi};  // stride period
		int nstride{2};  // should find 1 peak: end of first stride
		double del = (nstride*period)/(npts-1);
		int mpts = npts+1;  // 1 extra for finding the last peak
		ofstream ofs("/tmp/hrpts");
		for (int i=0; i<mpts; i++) {
			double wt = i*del;
			double x = cos(wt);
			accels.push_back(x);
			times.push_back(wt);
			ofs << x << endl;
		}
		system("vz -x points -y /tmp/hrpts&");
		last = npts - 1;  // 0b
		vector<int> pts = hrpts(first, last, accels);
		// what happens if we pick last as the next-to-last peak?
		last -= (npts-1)/(2*nstride);
		pts = hrpts(first, last, accels);
	}
	return 0;
}
#endif // MAIN
