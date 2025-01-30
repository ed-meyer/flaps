//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

// Functions and classes for I/O in Flaps.
// Object I/O: Stroustrup 3rd Ed. sect 25.4.1 & 4th Ed. sect 22.2.4

#include "config.h"

#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <dirent.h>
#include <fcntl.h>
#include <filesystem>
#include <iostream>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h> // access

#include "conv.h"
#include "fio.h"
#include "fma.h"
#include "matrix.h"
#include "message.h"
#include "Regex.h"
#include "trace.h"

using namespace std;


// io_map maps object ids to the objects "get" function. It is
// created in register_class() below
map<string, FioGetFcn > *io_map{nullptr};


Fio*
Fio::
get (Receiver& s) {
// Get an object from a Receiver:
// - read the virtual id
// - look up the id in the map of registered objects
// - call the object's serializer
// Throws a runtime_error if the object's class has not been registered,
// e.g. if it was stored with fio::store(vector<Type>)
	T_(Trace trc(2,"Fio::get");)

	// read the object's id from the Receiver
	string id;
	s.serialize(id);
	
	// XXX just call io_map.at?
	map<string, FioGetFcn >::const_iterator p = io_map->find(id);
	if (p == io_map->end()) {
		ostringstream os;
		os << "reading " << s.file_ << ": ";
		os << "class (" << id << ") has not been registered (" <<
			io_map->size() << " classes have been)\n";
		for (map<string,FioGetFcn >::const_iterator q = io_map->begin();
			q != io_map->end(); q++)
			os << (*q).first << " is registered\n";

		T_(trc.dprint("throwing exception: ",os.str());)
		throw runtime_error(os.str());
	}

	FioGetFcn f = (*io_map)[id];
	T_(trc.dprint("calling ",id," get serializer");)

	Fio* rval = f(s);

	return rval;
}

void
Fio::
put (Sender& s, const Fio& m) {
// Write an object to a Sender: serialize it's vid, then
// call it's virtual put() function
// XXX need a way to ensure m.put is never called directly
	T_(Trace trc(2,"Fio::put ",m.vid());)
	s.serialize(m.vid());
	m.put(s);
}


bool
Fio::
register_class (string id, FioGetFcn pf) {
// add a class name (id) to the map of registered classes
	T_(Trace trc(3,"Fio::register_class ",id);)
	if (io_map == nullptr)
		io_map = new map<string, FioGetFcn >;
	(*io_map)[id] = pf;
	return true;
}

void
putVFio (Sender& s, const std::vector<Fio*>& t) {
// Serialize an STL vector of pointers to classes that are
// derived from Fio and have registered (register_class)
	T_(Trace trc(2,"putVFio ",t.size());)

	s.serialize(t.size());
	for (auto tp : t) {
		T_(trc.dprint("putting a ",tp->vid());)
		Fio::put(s, *tp);
	}
}

void
getVFio (Receiver& s, vector<Fio*>& t) {
// de-serialize a vector of Fio objects by calling Fio::get for
// each element of the vector (the first object de-serialized is
// the length of the vector)
	T_(Trace trc(2,"getVFio");)
	size_t n;
	s.serialize(n);
	T_(trc.dprint("reading ",n," objects");)
	if (n == 0)
		return;

	for (size_t i=0; i<n; i++) {
		Fio* tp = Fio::get(s);
		t.push_back(tp);
		T_(trc.dprint("got object ",i+1,": \"",tp->vid(),"\"");)
	}
	return;
}

// Sender constructor
Sender::
Sender(const string& file) {
	T_(Trace trc(2,"Sender constructor");)
	// create a file on datadir
	string path = vastr(get_datadir(), "/", file);
	T_(trc.dprint("opening ",path);)
	int flags = O_CREAT | O_WRONLY | O_TRUNC;;
	mode_t mode = S_IRUSR | S_IWUSR;
	ofd = open(path.c_str(), flags, mode);
	if (ofd < 0) {
		string exc = vastr("cannot open ",path);
		T_(trc.dprint("throwing exception: ",exc);)
		throw runtime_error(exc);
	}
}

Sender::
~Sender() {
	close(ofd);
}

void
Sender::
write (const char* s, size_t nbytes) {
// write nbytes to this->ofd, throw an exception if the
// write was not successfull
	errno = 0;
	ssize_t nw = ::write(ofd, s, nbytes);
	if (nw == -1)
		throw std::system_error(errno, std::generic_category(),"write failed");

	if (nw < (ssize_t)nbytes)
		throw std::system_error(errno, std::generic_category(),
			vastr("attempt to write ",nbytes,", but only ",nw," written"));
}

// Sender serializers
void
Sender::
serialize(int  t) {
	T_(Trace trc(3,"Sender::serialize(int) ",t);)
	// ofs << t << endl;
	this->write((char*)&t, sizeof(int));
}

void
Sender::
serialize(size_t  t) {
	T_(Trace trc(3,"Sender::serialize(size_t) ",t);)
	// ofs << t << endl;
	this->write((char*)&t, sizeof(size_t));
}

void
Sender::
serialize(double  t) {
	T_(Trace trc(3,"Sender::serialze(double) ",t);)
	// ofs << std::showpoint << std::setprecision(20) << t << endl;
	this->write((char*)&t, sizeof(double));
}

void
Sender::
serialize(const complex<double>& t) {
	T_(Trace trc(3,"Sender::serialze(complex) ",t);)
	this->serialize(t.real());
	this->serialize(t.imag());
}

void
Sender::
serialize(const string&  t) {
	T_(Trace trc(3,"Sender::serialize(string) ",t);)
	size_t m = t.size();
	this->serialize(m);
	if (m == 0) return;
	// ofs << t << "\n";
	this->write(t.c_str(), m);
}

void
Sender::
serialize(const vector<vector<double> >& t) {
	T_(Trace trc(3,"Sender::serialize(vector<vector<double>>");)
	size_t m = t.size();
	this->serialize(m);
	if (m == 0) return;
	for (size_t i=0; i<m; i++) {
		this->serialize(t[i]);
	}
}

void
Sender::
serialize (const std::optional<double>& t) {
// c++17 optional<double>
	int init{1};	// t has been initialized
	double v{0.0};
	if (!t)
		init = 0;
	else 
		v = *t;
	serialize(init);
	serialize(v);
}

// Receiver
Receiver::
Receiver(const string& file) {
// Receiver constructor: throws runtime_error exception if the
// file does not exist
	T_(Trace trc(2,"Receiver constructor");)
	// the file must be on datadir
	string datadir = get_datadir();
	string path = vastr(datadir,"/",file);
	int flags = 0;
	mode_t mode = O_RDONLY;
	ifd = open(path.c_str(), flags, mode);
	if (ifd < 0) {
		T_(trc.dprint("\"",path,"\" does not exist");)
		string err = file + " does not exist";
		throw runtime_error(err);
	}
	file_ = file;
}

Receiver::
~Receiver() {
	close(ifd);
}

// Receiver serializers
void
Receiver::
serialize(int&  t) {
	T_(Trace trc(3,"Receiver::serialize(int)");)
	// if (ifs.eof())
	// 	throw runtime_error("eof in serialize(int)");
	// string st;
	// std::getline(this->ifs, st);
	// T_(trc.dprint("read <",st,">");)
	// t = stoi(st, nullptr);
	// T_(trc.dprint("got ",st," => ",t);)
	read(ifd, &t, sizeof(int));
}

void
Receiver::
serialize(size_t&  t) {
	// int it;
	// this->serialize(it);
	// t = (size_t)it;
	read(ifd, &t, sizeof(size_t));
}

void
Receiver::
serialize(double&  t) {
	T_(Trace trc(3,"Receiver::serialize(double)");)
	// string st;
	// std::getline(this->ifs, st);
	// str2double(st, t);
	// T_(trc.dprint("got ",st," => ",t);)
	read(ifd, &t, sizeof(double));
}

void
Receiver::
serialize(complex<double>& t) {
	double tr{0.0};
	double ti{0.0};
	this->serialize(tr);
	this->serialize(ti);
	t = complex<double>(tr,ti);
}

void
Receiver::
serialize(string&  t) {
// a string might have spaces and newlinee so we must read
// the string length, then the characters.
	T_(Trace trc(3,"Receiver::serialize(string)");)
	size_t m;
	this->serialize(m);
	if (m == 0) return;
	char* buf = new char[m+1];
	// ifs.read(buf, streamsize(m+1));  // read terminating newline

	read(ifd, buf, m);
	t = string(buf, string::size_type(m));
	delete[] buf;
	T_(trc.dprint("got <",t,">");)
}


void
Receiver::
serialize(vector<vector<double> >& t) {
	T_(Trace trc(3,"Receiver::serialize(vector<vector<double>>");)
	size_t m;
	this->serialize(m);
	if (m == 0) return;
	for (size_t i=0; i<m; i++) {
		vector<double> ti;
		this->serialize(ti);
		t.push_back(ti);
	}
}

void
Receiver::
serialize (std::optional<double>& t) {
// c++17 optional<double>
	int init{1};	// t has been initialized?
	this->serialize(init);
	double v;
	this->serialize(v);
	if (init == 0)
		t = {};
	else
		t = v;
}

static
bool
strwnum(const std::string& a, const std::string& b) {
// compare two strings which may contain numbers,
// e.g. a="mode_10" and b="mode_9", is a < b? without 
// comparing 10 and 9 as ints a !< b
	size_t n = std::min(a.size(), b.size());
	for (size_t i=0; i<n; i++) {
		if (isdigit(a[i]) && isdigit(b[i])) {
			size_t anc{0};
			size_t bnc{0};
			int ia = stoi(&a[i], &anc);
			int ib = stoi(&b[i], &bnc);
			if (ia != ib) {
				return ia < ib;
			} else {
				// XXX what if anc != bnc - impossible?
				i += anc;
			}
		} else {
			if (a[i] != b[i]) {
				return a[i] < b[i];
			}
		}
	}
	return false;  // a == b so not <
}

vector<string>
fio::
catalog (const string& rxname, bool sum) {
// Returns a vector of strings which are the names (sum=false or not included)
// or a summary (name, size, datatype) (sum=true) of each item in the database
// whose name matches regular expression "rxname"
	T_(Trace trc(2,"fio::catalog ",rxname);)
	vector<string> entries = ls (get_datadir());
	vector<string> rval;

	// if no regex given return all entries
	if (rxname.empty()) {
		for (auto& ei : entries) {
			if (ei != "." && ei != "..")
				rval.push_back(ei);
		}
		T_(trc.dprintv(rval, "no name given: returning all entries");)
	} else {
		regex re = make_regex(rxname);
		for (auto& ei : entries) {
			T_(trc.dprintn("testing <",ei,">...");)
			if (regex_match(ei, re)) {
				rval.push_back(ei);
				T_(trc.dprint(" matches");)
			} else {
				T_(trc.dprint(" no match");)
			}
		}
	}
	if (rval.empty()) {
		T_(trc.dprint("returning empty");)
		return rval;
	}
	// sort the entries
	if (rval.size() > 1) {
		T_(trc.dprintv(rval,"before sorting");)
		sort(rval.begin(),rval.end(), strwnum);
	}

	// if "sum" we must read each object and return it's summary
	// instead of just the name
	if (sum) {
		for (auto& name : rval) {
			try {
				ostringstream os;
				Receiver rcv(name);
				Fio* op = Fio::get(rcv);
				Matrix* mp = dynamic_cast<Matrix*>(op);
				if (mp == nullptr) {
					Fma* fp = dynamic_cast<Fma*>(op);
					if (fp == nullptr)
						os << name << " (" << op->vid() << ")";
					else
						os << *fp;
				} else
					os << *mp;
				name = os.str();		// replace this element of rval
			} catch (runtime_error& s) {
				T_(trc.dprint("ignoring exception: ",s.what());)
			}
		}
	}

	T_(trc.dprintv(rval,"returning");)
	return rval;
}

vector<string>
fio::
save (const vector< string>& rxnames, const string& file) {
// Read the fid's, which may be regular-expressions, and
// tar them to "file" compressed
	T_(Trace trc(2,"fio::save");)
	ostringstream os;
	vector<string> entries;
	string wd = get_cwd();

	T_(trc.dprint("rxnames: ",rxnames);)
	// get mids of existing matrices that match the regex's
	for (auto& rxname : rxnames) {
		T_(trc.dprint("cataloging ",rxname);)
		vector<string> ent = catalog(rxname);
		for (auto& ei : ent)
			entries.push_back(ei);
	}

	T_(trc.dprint("entries: ",entries);)

	if (entries.empty())
		return entries;

	string path;
	if (file[0] == '/')
		path = file;
	else
		path = wd + '/' + file;

	// move to datadir to do the tar
	string datadir = get_datadir();
	Chdir dd(datadir);
	os << "tar zcf " << path;
	for (auto fp : entries)
		os << " " << fp;
	T_(trc.dprint("system(",os.str(),")");)
	int stat = system(os.str().c_str());
	if (stat != 0) {
		cerr << "failed to save matrices\n";
		entries.clear();
	}
	return entries;
}

vector<string>
fio::
restore(const string& file) {
// restore matrices saved on "file" (by fio::save()) to datadir
	T_(Trace trc(2,"fio::restore");)

	// create a full path if not already
	string fullpath{file};
	if (file[0] != '/') {
		fullpath = get_cwd() + '/' + file;
	}
	// change to datadir to untar
	string datadir{get_datadir()};
	Chdir wd(datadir);
	// run a shell command to untar
	string cmd{vastr("tar zvxf ",fullpath)};
	T_(trc.dprint("system(",cmd,")");)

	// open a pipe to the command...
	Fpipe pipe{cmd};

	// ... and read it's output: a list of matrices in the tar
	// file (v option to tar)
	FILE* fp = pipe();
	size_t bufsize{256};
	ssize_t n;
	char* buf{(char*)malloc(bufsize)};
	vector<string> rval;
	while ((n = getline(&buf, &bufsize, fp)) != -1) {
		char* last = &buf[n];
		while (*(last-1) == '\n')
			last--;
		string item{&buf[0], last};
		T_(trc.dprint("got ",n," char <",item,">");)
		rval.push_back(item);
	}

	// return a list of the objects restored
	return rval;
}

Fd::
Fd(const string& file, int flags, mode_t mode) {
// Fd constructor: RAII class to open "file" with options "flags"
// and optionally "mode". See "man open" for flags and mode
// Throws system_error if open fails
// The file is closed in the destructor when this object goes out of scope
	fd = open(file.c_str(), flags, mode);
	if (fd == -1)
		throw std::system_error(errno, std::generic_category(),
				vastr("cannot open ",file));
}

Fd::
~Fd() {
	if (close(fd) == -1) {
		// must not throw (Stroustrup sect 13.2)
		flaps::warning("Fd failed to close: ", strerror(errno));
	}
}

File::
File(const string& file, const string& mode) {
// File constructor: RAII class to open "file" with options "flags"
// and optionally "mode". See "man open" for flags and mode
// Throws system_error if open fails
// The file is closed in the destructor when this object goes out of scope
	fp = fopen(file.c_str(), mode.c_str());
	if (fp == nullptr)
		throw std::system_error(errno, std::generic_category(),
				vastr("cannot open ",file));
}

File::
~File() {
	if (fclose(fp) == EOF) {
		// must not throw (Stroustrup sect 13.2)
		flaps::warning("File failed to close: ", strerror(errno));
	}
}

Msgq::
Msgq(const string& name, int flags, mode_t mode, struct mq_attr* attr) {
	T_(Trace trc(1,"Msgq constructor");)
	name_ = name;
	if (attr == nullptr) {
		creator = false;
		des_ = mq_open(name_.c_str(), flags);
	} else {
		creator = true;
		des_ = mq_open(name_.c_str(), flags, mode, attr);
	}
	if (des_ == -1)
		throw std::system_error(errno, std::generic_category(),"mq_open");
	T_(trc.dprint("opened ",name,", creator? ",creator,", des ",des_);)
}

Msgq::
~Msgq() {
	T_(Trace trc(2,"Msgq destructor");)
	if (mq_close(des_) == -1)
		cerr << "mq_close failed, errno " << errno << endl;
	T_(trc.dprint("closed ",name_);)
	if (creator) {
		if (mq_unlink(name_.c_str()) == -1)
			cerr << "mq_close failed, errno " << errno << endl;
		T_(trc.dprint("unlinked ",name_);)
	}
}

void
Msgq::
send(const std::string& v) {
	int stat = mq_send(des_, v.c_str(), v.size(), 0);
	if (stat == -1)
		throw std::system_error(errno, std::generic_category(),"mq_send");
}

void
Msgq::
receive(std::string& v) {
// read a string from a message queue
	T_(Trace trc(2,"Msgq::receive(string)");)
	struct mq_attr attr;
	mq_getattr(des_, &attr);
	long msgsize = attr.mq_msgsize;
	std::vector<char> msg(msgsize);
	ssize_t nbytes{0};
	//!! struct timespec ts;
	//!! clock_gettime(CLOCK_REALTIME, &ts);
	//!! ts.tv_sec += 5;	// wait 5 sec
	// do {
	for (int i=0; i<5; i++) {
		errno = 0;
		nbytes = mq_receive(des_, msg.data(), msgsize, nullptr);
		T_(trc.dprint("errno ",errno,", nbytes ",nbytes);)
		if (errno == 0 && nbytes > 0) break;
		if (nbytes == -1 && errno != EAGAIN) {
			T_(trc.dprint("throwing system_error exception");)
			throw std::system_error(errno, std::generic_category(),"mq_receive");
		}
		sleep(1);
	}
	// } while (nbytes == -1 && clock_gettime(CLOCK_REALTIME, &ts) != 0);
	v.assign(msg.data(), nbytes);
}

void
Msgq::
notify(int signo) {
// notify this process that a message has been added to an
// empty queue by sending a SIGUSR1 signal
	struct sigevent sigev;
	sigev.sigev_notify = SIGEV_SIGNAL;	// send a signal, not a thread
	sigev.sigev_signo = signo;
	mq_notify(des_, &sigev);
}

Fdir::
Fdir(const string& path) {
// Fdir constructor: RAII method for opening a directory (DIR struct)
// Throws system_error if opendir fails
	dir = opendir(path.c_str());
	if (dir == nullptr)
		throw std::system_error(errno, std::generic_category(),
			vastr("opendir failed to open ",path));
}

Fdir::
~Fdir() {
	if (closedir(dir) == -1) {
		// must not throw (Stroustrup sect 13.2)
		flaps::warning("closedir failed to close: ",strerror(errno));
	}
}

Fpipe::
Fpipe(const string& cmd, const string& mode) {
// Fpipe constructor: RAII method for opening a pipe
// mode: must be either "r" or "w"
// Throws system_error if popen fails
	comm = cmd;
	fp = popen(cmd.c_str(), mode.c_str());
	if (fp == nullptr)
		throw std::system_error(errno, std::generic_category(),
			vastr("cannot open pipe to ",cmd));
}

Fpipe::
~Fpipe() {
	if (pclose(fp) == -1) {
		// must not throw (Stroustrup sect 13.2)
		flaps::warning("failed to close pipe to \"", comm, "\": ",
			strerror(errno));
	}
}


Chdir::
Chdir(const string& dir) : cwdfd(".", O_RDONLY), cwd(get_cwd()) {
// Chdir constructor: RAII method for changing directories so that
// the destructor changes back to the current directory.
// Method: open current directory, change to new dir returning then
// is simply a call to fchdir()
// Throws system_error if chdir fails
	if (chdir(dir.c_str()) == -1)
		throw std::system_error(errno, std::generic_category(),
			vastr("Chdir failed to change to ",dir));
}

Chdir::
~Chdir() {
// Chdir destructor: change back to where we were
	if (fchdir(cwdfd()) == -1) {
		// must not throw (Stroustrup sect 13.2)
		flaps::warning("Chdir destructor failed to change back: ",strerror(errno));
	}
}

string
getfroot() {
// returns the full path of the root directory of the flaps installation
	string froot = getEnv("FROOT");
	return froot;
}

string
getftmp () {
/*------------------------------------------------------------------
 * returns the full path of the temporary directory for this
 * Flaps job. This path will be the same as environment variable FTMP.
 *------------------------------------------------------------------*/
	T_(Trace trc(3,"getftmp");)
	string rval{getEnv("FTMP")};

#ifdef NEVER // return empty if none
	// if no temp directory has been created just return /tmp
	if (rval.empty())
		rval = "/tmp";
#endif // NEVER // return empty if none
	T_(trc.dprint("returning ",rval);)
	return rval;
}

string
get_datadir () {
/*------------------------------------------------------------------
 * returns the full path of the temporary directory containing
 * data for this Flaps job.
 *------------------------------------------------------------------*/
	T_(Trace trc(3,"get_datadir");)
	string rval{getftmp()};

	// data is kept in directory "data" under ftmp
	rval += "/data";
	T_(trc.dprint("returning ",rval);)
	return rval;
}

string
get_cwd() {
// return the current working directory as a string
	T_(Trace trc(3,"get_cwd");)
	string rval;
#ifdef _GNU_SOURCE
	// using GNU-specific replacement for getcwd()
	char* cwd = get_current_dir_name();
	T_(trc.dprint("get_current_dir_name returned ",cwd);)
	if (cwd == nullptr)
		throw std::system_error(errno, std::generic_category(),
			"cannot get current working directory: ");

	rval = string(cwd);
	T_(trc.dprint("returning ",rval);)
	free(cwd);
#else
	size_t len{PATH_MAX};
	char* buf = new char[len];
	char* wd = getcwd(buf, len);
	rval = string(wd);
	delete[] buf;
#endif
	return rval;
}

std::string
fio::
shortenpath(const std::string& path) {
    std::filesystem::path fullPath = std::filesystem::absolute(path);
    string cwd = std::filesystem::current_path().string();

	// first try removing cwd
	 string starts_with = fullPath.string().substr(0,cwd.size());
    if (starts_with == cwd)
	 	return fullPath.string().substr(cwd.size()+1);
	// then try removing home
	std::string homeDir = std::getenv("HOME");
	starts_with = fullPath.string().substr(0,homeDir.size());
	if (starts_with == homeDir)
		return "~" + fullPath.string().substr(homeDir.size());

#ifdef NEVER // needs work
    int commonPrefixLength = 0;
    while (commonPrefixLength < fullPath.string().length() &&
           commonPrefixLength < cwd.string().length() &&
           fullPath.string()[commonPrefixLength] == cwd.string()[commonPrefixLength]) {
        commonPrefixLength++;
    }

    std::string shortenedPath = "";
    for (int i = 0; i < cwd.string().length() - commonPrefixLength; i++) {
        shortenedPath += "../";
    }
    shortenedPath += fullPath.string().substr(commonPrefixLength);

    return shortenedPath;
#endif // NEVER // needs work
	return path;
}

string
getfwd () {
/*------------------------------------------------------------------
 * returns the full path of the working directory for this
 * flaps job.
 *------------------------------------------------------------------*/
	T_(Trace trc(1,"getfwd");)
	string rval{getEnv("FWD")};

	// if FWD is not already in the environment set it to the
	// current working directory
	if (rval.empty()) {
		rval = get_cwd();
		string fwd = vastr("FWD=",rval);
		putEnv(fwd);
	}
	T_(trc.dprint("returning ",rval);)
	return rval;
}

vector<string>
ls (string const& path) {
// returns a vector of the names of all files in "path"
	vector<string> rval;
	struct dirent* entry;
	ostringstream os;

	// open using RAII class Fdir...
	Fdir fdir(path);

	while ((entry = readdir(fdir())) != 0) {
		string name{entry->d_name};
		if (name != "." && name != "..")
			rval.push_back(name);
	}
	return rval;
}

#ifdef HAVE_NFTW
#ifdef HAVE_FTW_H
#include <ftw.h>
#endif // HAVE_FTW_H

static vector<string> RmErr;

static int
rmdirTreeFcn (char const* file, struct stat const* sb, int type, struct FTW* s) {
	T_(Trace trc(2,"rmdirTreeFcn ", file);)
	errno = 0;
	int status = remove(file);
	T_(trc.dprint("remove(",file,") returned status ",status,", errno ",errno);)
	T_(trc.dprint("cwd ",get_cwd());)
	if (status < 0) {
		string err;
		if (errno) {
			err = strerror(errno);
			RmErr.push_back(err);
		} else {
			err = vastr("remove returned ",status);
		}
		T_(trc.dprint("returning 1: ",err);)
		return 1;
	}
	return 0;
}
#endif /* HAVE_NFTW */

void
rmdirTree (string const& path) {
// Remove a directory and all the files in it. Does the
// equivalent of "rm -rf path" and in fact if the function
// nftw() is not available that is what is done
	T_(Trace trc(2,"rmdirTree ", path);)

#ifndef HAVE_NFTW
	ostringstream os;
	os << "rm -rf " << path;
	system(os.str().c_str());
#else // HAVE_NFTW
	int nopenfd = 20;
 	int flags = FTW_DEPTH;
	ostringstream os;

	RmErr.clear();

	// nftw (file tree walk) calls rmdirTreeFcn for each file
 	int status = nftw (path.c_str(), rmdirTreeFcn, nopenfd, flags);

	if (status != 0)
		throw std::system_error(errno, std::generic_category(),
			vastr("cannot remove ",path));

	if (!RmErr.empty()) {
		os.str("");
		for (size_t i=0; i<RmErr.size(); i++) {
			if (i > 0)
				os << "; ";
			os << RmErr[i];
		}
		throw runtime_error(os.str());
	}
#endif // HAVE_NFTW
}

string
removeHOME (string const& path) {
/*------------------------------------------------------------------
 * Given a path, returns a string containing that portion of
 * the path which is relative to the current "HOME" environment
 * variable if the path starts with it; otherwise the input path
 * is returned.
 *------------------------------------------------------------------*/
	T_(Trace trc(1,"removeHOME");)
	char const* home = getenv("HOME");
	// char* logname = getenv("LOGNAME");
	string rval;

	T_(trc.dprint("path<",path,">");)
	
	if (!home) {
		T_(trc.dprint("quick return: no HOME in env");)
		return path;
	}

	if (path == string(home)) {
		T_(trc.dprint("returning \".\": path==home");)
		return string(".");
	}

	string homestr(home);
	if (homestr[homestr.size()-1] != '/') {
		homestr += '/';
	}
	if (path.size() > homestr.size()) {
		if (path.substr(0,homestr.size()) == homestr) {
			rval = path.substr(homestr.size());
			T_(trc.dprint("returning \"",rval,"\"");)
			return rval;
		}
	}
	/*
	 * HOME not found - just return the input string
	 */
	T_(trc.dprint("returning input: \"",path,"\"");)
	return path;
}

static void
add2PATH (string& path) {
// add "path" to the PATH environment variable
	char* cp = getenv("PATH");
	string newpath = vastr(path,":",cp);
	putEnv(newpath);
}

Ftmpdir::
Ftmpdir(const string& path) {
// Ftmpdir constructor
// - put the Flaps working director in the environment if needed
// - create a new directory under /tmp: Ftmp+pid
// - create data directory under Ftmp
// - put FTMP=Ftmp+pid in the environment
	T_(Trace trc(1,"Ftmpdir constructor");)

	// put the current path in the environment: getfwd()
	// does this if FWD is not already in the environment
	string fwd = getfwd();

	// if FTMP is in the environment and it exists, do nothing
	// (leave _tmpdir empty)
	char* cp = getenv("FTMP");
	if (cp != nullptr) {
		int mode = R_OK | W_OK | X_OK;
		int stat = access(cp, mode);
		if (stat == 0) {
			return;
		} else {
			string msg{vastr("FTMP is in the environment (",cp,
				") but is not accessible: ",strerror(errno))};
			throw runtime_error(msg);
		}
	}

	ostringstream os;
	int status{0};
	// if "path" is empty create a new directory
	if (path.empty()) {
		// Create the temporary directory on /tmp so it does not
		// get saved in a cloud (Dropbox, Drive)
		os << "/tmp/Ftmp" << getpid();
		_tmpdir = os.str();
		status = mkdir(_tmpdir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
		if (status == -1)
			throw std::system_error(errno, std::generic_category(),
				vastr("cannot create temp directory (", _tmpdir));

		T_(trc.dprint("Flaps temporary directory: ", _tmpdir);)

	} else {
		_tmpdir = path;
	}

	// also create /tmp/Ftmp/data for storing matrices
	os.str("");
	os << _tmpdir << "/data";
	string datadir = os.str();
	status = mkdir(datadir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
	if (status == -1)
		throw std::system_error(errno, std::generic_category(),
			vastr("cannot create data directory (", datadir));

	// add this directory to the PATH
	add2PATH(_tmpdir);

	//Set env variable FTMP
	os.str("");
	os << "FTMP=" << _tmpdir;
	putEnv (os.str());
	// print a message if the temp dir is not to be deleted
	if (getenv("KEEPFTMP") != nullptr)
		flaps::info("the temporary directory ",_tmpdir," will not be deleted");
}

Ftmpdir::
~Ftmpdir() {
// Deletes the temporary directory unless KEEPFTMP is in the environment
// or if DEBUG is in the environment
	T_(Trace trc(1,"Ftmpdir destructor");)

	// sometimes multiple Ftmpdirs are created but only one
	// has an actual directory: quick return if we don't have one
	if (_tmpdir.empty()) {
		T_(trc.dprint("quick return: no ftmpdir");)
		return;
	}

	string msg;
	char* keep = getenv("KEEPFTMP");
	if (keep != nullptr)
		msg = "KEEPFTMP is in the environment";
	char* dbg = getenv("DEBUG");
	if (dbg != nullptr) {
		if (!msg.empty())
			msg += " and ";
		msg += "DEBUG is in the environment";
	}
	if (keep == nullptr && dbg == nullptr) {
		rmdirTree(_tmpdir);
	} else {
		flaps::info ("Keeping temporary directory ",_tmpdir,": ",msg);
	}
}

#ifdef MAIN

int
main (int argc, char** argv) {
	Ftmpdir ftmp;
	putEnv("KEEPFTMP=1");
	size_t i;
	// test Chdir
	cout << "Flaps temp directory: " << getftmp() << endl;
	cout << "current directory: " << get_cwd() << endl;
	{
		Chdir ftmp(getftmp());
		cout << "ftmp  directory: " << get_cwd() << endl;
	}
	cout << "should be back to original directory: " << get_cwd() << endl;
	// test ls()
	string path("/tmp");
	vector<string> entries = ls(path);
	cout << entries.size() << " entries in " << path << endl;
	for (auto& e : entries)
		cout << e << endl;
	string id("test");
	string altid("taste");
	vector<string> svec;
	svec.push_back("abc");
	svec.push_back("def");
	svec.push_back("ghi");
	try {
		fio::store (id, svec);
		fio::store (altid, svec);
		vector<string> newsvec;
		fio::fetch (id, newsvec);
		for (auto& ni : newsvec) {
			cout << ni << endl;
		}
		// doubles
		double pi{flaps::pi};
		vector<double> dbl;
		for (i=0; i<5; i++) {
			dbl.push_back(pi*(i+1));
		}
		string dblid("dbl");
		fio::store(dblid, dbl);
		// ... then read it
		vector<double> newdbl;
		fio::fetch(dblid, newdbl);
		for (i=0; i<newdbl.size(); i++) {
			cout << "read dbl " << newdbl[i];
			cout << ", diff " << newdbl[i]-dbl[i] << endl;
		}

		// catalog just the files "test" and "taste" with a regex
		vector<string> te = fio::catalog("t[ae].t.*");
		for (auto& ti : te)
			cout << "got file name \"" << ti << "\"\n";
		// catalog all files
		te = fio::catalog();
		for (auto& ti : te)
			cout << "got file name \"" << ti << "\"\n";
		// test shortenpath
		//!! string path = get_cwd();
		string path{"/home/eem2314/git/build/src/bin/flut"};
		cout << "shortenPath(" << path << " = " << shortenPath(path) << endl;
	} catch (runtime_error& s) {
		cout << "Exception: " << s << endl;
	}
}
#endif	// MAIN
