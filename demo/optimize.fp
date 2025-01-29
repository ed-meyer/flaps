# Demonstration of using a matlab optimizer to adjust the betas
# in an RFA fit of DLM aero.

# Goland wing from references 1 and 2 section 9.3.
# References
#   1) @article{goland1945flutter,
#      title={The flutter of a uniform cantilever wing},
#      author={Goland, Martin},
#      year={1945},
#      publisher={American Society of Mechanical Engineers}
#    }
#   2) @book{bisplinghoff1955aeroelasticity,
#        title={Aeroelasticity},
#        author={Bisplinghoff, Raymond L and Ashley, Holt and Halfman, Robert L},
#        year={1955},
#        publisher={Addison-Wesley},
#        address={Reading, Mass.}
#    }
# Cantilevered straight wing: demonstrates beam modeling and doublet-lattice
# unsteady aerodynamics
   # using the default units: meters, Kg
   fem {
      nodes{1{0,0,0} : 13{0,6,0} }
      material{x=(1,0,0), gj=9.876e+5, eiz=3.0e+7, eix=9.773e+6, ea=2.0e+5}
      beam{1{0}, 2{3:5}, 3 : 13}
      mass{mass=18.14, moi=(1,4.395,0), cg=0.183, orient=x, node=(2:13)}
      dlm{panels=4, chord=1.83, rf=(0.001,0.1,0.3,2.0), mach=0.5, ac=-0.5}
		nmodes = 10
   }
	catalog{}

	# betas from gordon2008nonlinear: (0.0005, 0.005, 0.1, 0.5, 0.9)
 	pz { i=gaf[1234], o=rfa, beta=(0.0001,0.005,0.01) }

   flut {
      id=interp
      indep=(vtas[0:250],sigma,freq)
      alt=0
      mass=mass
      stif=stif
      gaf=gaf
      start{modes=(1,2)}
      target{sigma=0, vtas[10:250]}
   }
   flut {
      id=initial
      indep=(vtas[0:250],sigma,freq)
      alt=0
      mass=mass
      stif=stif
      gaf=rfa
      start{modes=(1,2)}
      target{sigma=0, vtas[10:250]}
   }
	# compute betas which result in the critical flutter speed
	# close to the speed with interpolated gaf
	call { optimize{beta=(0.01, 0.05, 0.1, 0.5, 1.0), maxfun=200} }
	# show the difference between interpolated, initial betas,
	# and optimized betas in the critical flutter mode
	vz {id=(interp,initial,optimize), x=vtas, y=sigma}
	vz {id=(interp,initial,optimize), x=sigma, y=freq}
	vz {y=f, betas.apf}
	vz {y=beta[12345], betas.apf}
end

fxing.cpp {{
#include <cstring>
#include <unistd.h>
#include <cstdio>
#include <fcntl.h>
#include <sys/mman.h>
int fxing(pset& targets) {
	Trace trc(1,"fxing");
	double target{124.851};
	cerr << "fxing called with " << targets.nsolns() << endl;
	Par* vtasp = targets.findp("vtas");
	cerr << "found " << vtasp->solns.size() << " vtas targets: " << vtasp->solns[0] << endl;
	double diff = abs(vtasp->solns[0] - target);
	cerr << "sending " << diff << endl;
	if (access("function", W_OK) != 0)
		cerr << "cannot open function: " << strerror(errno) << endl;
	int fd = open("function", O_RDWR);
	if (fd == -1)
		throw std::system_error(errno,std::generic_category(),
			vastr("cannot open function: ",strerror(errno)));

	// create a memory-map linked to the file
	size_t nbytes = 2*sizeof(double);
	double* ptr = (double*)mmap(0, nbytes, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
	if (ptr == MAP_FAILED) {
		string exc{vastr("cannot map function: ",strerror(errno))};
		cerr << "fxing throwing exception: " << exc << endl;
		throw std::system_error(errno,std::generic_category(),exc);
	}
	ptr[0] = 1.0;
	ptr[1] = diff;
	cerr << "sent function = " << diff << endl;
	return 0;
}
}}

optimize.cpp {{
#include <cstdio>
#include <cstring>
#include <fcntl.h>
#include <sys/mman.h>
#include "pac.h"
#include "run_program.h"
double*
open_octlab(string const& file, size_t n) {
// open shared memory file for communicating with octlab, with n doubles;
// returns a pointer to the shared memory with n+1 doubles - the first is
// a semaphore: if it is zero the data has been read and can be overwritten,
// otherwise it is the number of doubles following.
	Trace trc(1,"open_octlab");

	int fd = open(file.c_str(), O_RDWR | O_CREAT | O_TRUNC, S_IRWXU);
	if (fd == -1)
		throw std::system_error(errno,std::generic_category(),
			vastr("cannot open ",file,": ",strerror(errno)));

	// create the file with n+1 doubles: the first double is the length
	size_t nbytes = (n+1)*sizeof(double);
	vector<double> buf(n+1, 0.0);
	write(fd, buf.data(), nbytes);

	// create a memory-map linked to the file
	double* rval = (double*)mmap(0, nbytes, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
	if (rval == MAP_FAILED) {
		string exc{vastr("cannot map ",file,": ",strerror(errno))};
		trc.dprint("throwing exception: ",exc);
		throw std::system_error(errno,std::generic_category(),
			vastr("cannot map ",file,": ",strerror(errno)));
	}
	return rval;
}
	
string
script(const vector<double>& betas, int maxfunevals) {
// create a script to run in matlab/octave
// XXX run it in a different thread?

	string cmd("matlab");
	// string cmd("octave");
	ofstream ofs("optbeta.m", std::ios::trunc);
	ofs << "cd " << getftmp() << endl;
	ofs << "x = [";
	string sep;
	for (auto bi : betas) {
		ofs << sep << bi;
		sep = "; ";
	}
	ofs << "];\n";

	ofs << "opt = optimset(\'MaxFunEvals\'," << maxfunevals
		<< ",\'TolX\',0.1,\'TolFun\',0.1);\n";
	ofs << "beta = fminsearch(@optfcn,x,opt)\n";

	ofs << "enviar = memmapfile(\'betas\', \'Writable\',true,\'Format\', \'double\');\n";
	ofs << "len = length(beta);\n";
	ofs << "enviar.Data(2:len+1) = beta;\n";
	ofs << "enviar.Data(1) = -1;\n";
	ofs << "\'sent final betas\'\n";
	ofs << "quit\n";
	// create the command line to run: only matlab works - octave's memmapfile
	// doesn't work
	// unset LD_PRELOAD while running matlab (libstdc++ conflict)
	string cmdline;
	if (cmd == "matlab") {
		cmdline = vastr("LD_PRELOAD=\"\" ",cmd," -batch \"optbeta\"&");
	} else {
		cmdline = vastr(cmd, " optbeta.m&");
	}

	// run matlab/octave in the background
	int status = system(cmdline.c_str());

	return cmdline;
}

int
optimize(string const& pref) {
	Trace trc(1,"optimize");

	// change to temp directory using RAII class Chdir
	Chdir tmp(getftmp());

	// parse the options
	vector<double> initial_betas;
	int maxfunevals{60};
	vector<Pref*> unrec = Pref::parse(pref, {
		{"beta", [&](const Pref& p) { initial_betas = p.rvec; return true; }},
		{"maxfun", [&](const Pref& p) { maxfunevals = p.ivec[0]; return true; }}
	});
	trc.dprint("got betas: ",initial_betas);

	int nbeta = initial_betas.size();

	// open memory-mapped files to send/receive data to matlab/octave
	double* enviar = open_octlab("function", 1);
	double* recibir = open_octlab("betas", nbeta);

	// start matlab/octave with the initial betas
	string file = script(initial_betas, maxfunevals);

	vector<double> x(nbeta, 0.0);
	// vector<vector<double>> results;
	vector<double> results;
	vector<double> result(nbeta+1, 0.0);
	for (int iter=0; iter<maxfunevals; iter++) {
		// wait for a recibir message: if the length is -1 quit
		while (recibir[0] == 0.0)
			sleep(1);
		int len = recibir[0];
		cout << "got length = " << len << endl;
		if (len == -1)
			break;
		if (len != nbeta)
			throw runtime_error(vastr("expected ",nbeta," betas, got ",len));
		for (int i=0; i<len; i++) {
			x[i] = recibir[i+1];
			result[i] = x[i];
		}
		// results.push_back(x);
		// vappend(results, nbeta, x.data());
		results.insert(results.end(), x.cbegin(), x.cend());
		cout << "got x = " << x << endl;
		// signal that we got the message, then optfcn will wait for a reply
		recibir[0] = 0.0;
		trc.dprint("set recibir[0] to 0");

		// create new RFA aero
		ostringstream os;
		os << "o=rfa, i=gaf[1234],beta=(";
		string sep{""};
		for (auto bi : x) {
			os << sep << bi;
			sep = ", ";
		}
		os << ")";
		pid_t childpid;
		double cputime;
		run_program("pz", os.str(), true, childpid, cputime);

		// then flut: fxing will send the target back to optfcn
		string flutopt{"id=optimize, indep=(vtas[10:250], freq, sigma), alt=0,"
			"mass=mass, stif=stif, gaf=rfa, start{modes=(1,2)},"
			"target{sigma=0, vtas[10:250], process=fxing}"};
		flutopt += vastr(", plotfile=",tmp.from(),"/optimize.pf");
		run_program("flut", flutopt, true, childpid, cputime);

		// save the function value
		result[nbeta] = enviar[1];
		// vappend(results, nbeta+1, result.data());
		results.insert(results.end(), result.cbegin(), result.cend());
	}

	// get the final betas
	cerr << "waiting for final betas\n";
	vector<double> final_betas(nbeta,0.0);
	while (recibir[0] == 0.0)
		sleep(1);
	int len = recibir[0];
	cout << "got final betas, length = " << len << endl;
	for (int i=0; i<nbeta; i++) {
		final_betas[i] = recibir[i+1];
		result[i] = final_betas[i];
	}
	result[nbeta] = enviar[1];
	// vappend(results, nbeta+1, result.data());
	results.insert(results.end(), result.cbegin(), result.cend());
	cout << "got final betas = " << final_betas << endl;

	// plot results
	vector<pair<string,string>> params;
	for (int i=0; i<nbeta; i++) {
		string name = vastr("beta",i+1);
		params.push_back({name,name});
	}
	params.push_back({"f","function"});
	string plotfile = vastr(tmp.from(),"/betas.apf");
	Apf::exporter("betas", params, results, plotfile, false);

	return 0;
}
}}
