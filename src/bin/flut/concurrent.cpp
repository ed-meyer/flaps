//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

#include <algorithm>
#include <thread>
#include <vector>

#include "concurrent.h"
#include "specs.h"
#include "trace.h"

using namespace std;

int
concurrently(int act) {
// set or return whether or not (0/1) to run concurrently
	static int crnt{1};  // default: concurrent
	int rval = crnt;
	if (act >= 0) {
		crnt = act;
	}
	return rval;
}

void
concurrent(Fstack* curves, void (*fcn)(Fstack*,int)) {
// Track curves concurrently by calling fcn if:
// - there are more than 1 hardware threads
// - there are more than 2 curves,
// - if concurrently() is true
// Otherwise track them sequentially, calling fcn one curve at a time
	Trace trc(1,"concurrent");

	trc.dprint(curves->size()," curves to trace");

	// leave one thread alone for the user
	const int hwthreads = std::thread::hardware_concurrency() - 1;


	if (curves->size() > 1 && concurrently() && hwthreads > 1) {
		trc.dprint(hwthreads," hardware threads");
		int ncon = std::min(static_cast<int>(curves->size()), hwthreads);
		flaps::info("tracking concurrently ", ncon, " at a time");
		size_t k{0};
		// do hwthreads curves at a time saving one for me
		int const max_threads = curves->size() - k;
		int const nthread = std::min(hwthreads, max_threads) - 1;
		trc.dprint("max_threads ",max_threads,", nthread ",nthread);
		vector<thread> threads;
		for (int i=0; i<nthread; i++) {
			threads.push_back(std::thread(fcn, curves, i+1));
		}
		fcn(curves,0);
		// ... wait until they are all finished
		std::for_each(threads.begin(), threads.end(),
				std::mem_fn(&std::thread::join));
	} else {
		flaps::info("tracking sequentially");
		fcn(curves,0);
	}
}
