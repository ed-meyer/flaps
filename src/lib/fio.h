//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
#ifndef FIOb_h
#define FIOb_h

// Object I/O: Stroustrup 3rd Ed. sect 25.4.1 and 4th Ed. sect 22.2.4
// XXX need to convert from stdio to iostreams, i.e.:
//    std::ofstream ofs;
// instead of
//   	int ofd, ifd;
// this requires using iostreams in binary mode

// Flaps I/O functions
// I/O in Flaps consists of 4 layers:
//   fetch/store)
//     for the highest-level classes (Matrix, Curve, pset).
//     These classes inherit from class Fio and provide
//     functions store() and fetch() which open the relevant file
//     and call Fio functions Fio::get() or Fio::put() to serialize
//     objects to the file. These classes must register their get()
//     function (see register_class below) prior to calling Fio::get().
//   Fio::get()/Fio::put() classes that inherit from Fio and have registered
//     by calling register_class may call these functions to read/write
//     to an open file. Registration allows Fio::get() to read a file
//     without knowing what classes the objects on the file belong to;
//     if it is known for certain what classes are represented, the
//     class's get() function may be called directly.
//   get/put: Object I/O. Classes that use this inheirit from
//     class Fio and are required to provide functions get and
//     put which serialize members with calls to the lowest-level
//     serialization functions and Fio::get/Fio::put for Fio-derived
//     members.
//   serialize: classes Receiver and Sender represent files to be
//     read/written. These classes have several serialize() functions
//     which read/write low-level data like int, string, double, etc.
//     Classes get/put functions call these serialization functions.
//
// Any class "MyClass" that uses Fio must:
//
// 1) include the following member functions:
//      static bool regd;  // has been registered?
//      static std::string id() { return "MyClass"; }
//      std::string vid() const { return id(); }
//      void put (Sender& s) const;
//      static Fio* get (Receiver& s) { return new MyClass(s); }
//      MyClass(Receiver& s) {
//         // serialize basic type members...
//         // ... call Fio::get() for each Fio-derived member
//      }
//   In addition, 
//
// 2) register the class's "get" function with a statement like this
//    in the class's .c file:
//       bool MyClass::regd = Fio::register_class(MyClass::id(), MyClass::get);


class Fio;

#include <complex>
#include <deque>
#include <dirent.h>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <map>
#include <mqueue.h>
#include <optional>
#include <sstream>
#include <string>
#include <typeinfo>
#include <vector>

#include "csignals.h"
#include "message.h"

class Fpipe;

class Sender {
	// std::ofstream ofs;
	int ofd;
public:
	Sender(const std::string& file);
	~Sender();

	void serialize (int t);
	void serialize (size_t t);
	void serialize (double t);

	void serialize (const std::complex<double>& t);

	void serialize (const std::string& t);

	template<typename Type>
	void serialize (const std::vector<Type>& t);

	template<typename Type>
	void serialize (size_t n, const Type* t);

	template<typename Type>
	void serialize (const std::deque<Type>& t);
	// map<Type>
	template<typename Type>
	void serialize (const std::map<std::string,Type>& t);

	// matrix (vector<vector>)
	void serialize (const std::vector<std::vector<double> >& t);

	// c++17 optional<>
	void serialize (const std::optional<double>& t);

	// serializers write with my write():
	void write (const char* s, size_t nbytes);

	bool good() { return true; }
	bool eof() { return false; }
};  // Sender

class Receiver {
	// std::ifstream ifs;
	int ifd;
public:
	std::string file_;

	Receiver(const std::string& file);
	~Receiver();

	void serialize (int& t);
	void serialize (size_t& t);
	void serialize (double& t);

	void serialize (std::complex<double>& t);

	void serialize (std::string& t);

	template<typename Type>
	void serialize (std::vector<Type>& t);

	template<typename Type>
	void serialize (std::vector<std::vector<Type> >& t);

	template<typename Type>
	void serialize (size_t& n, Type*& t);

	template<typename Type>
	void serialize (std::deque<Type>& t);

	template<typename Type>
	void serialize (std::map<std::string,Type>& t);

	void serialize (std::vector<std::vector<double> >& t);

	// c++17 optional<>
	void serialize (std::optional<double>& t);

	// bool good() { return ifs.good(); }
	bool good() { return true; }
	bool eof() { return false; }
};  // Receiver

template<typename Type>
void
Sender::
serialize(size_t n, const Type* t) {
// An array of Types
	this->serialize(n);
	for (size_t i=0; i<n; i++)
		this->serialize(t[i]);
}

template<typename Type>
void
Sender::
serialize(const std::vector<Type>& t) {
// A vector of Types
	this->serialize(t.size(), &t[0]);
}

template<typename Type>
void
Sender::
serialize(const std::deque<Type>& t) {
	this->serialize(t.size());
	for (auto& ti : t) {
		this->serialize(ti);
	}
}

template<typename Type>
void
Sender::
serialize(const std::map<std::string,Type>& t) {
	this->serialize(t.size());
	for (auto& ti : t) {
		this->serialize(ti.first);
		this->serialize(ti.second);
	}
}

template<typename Type>
void
Receiver::
serialize(size_t& n, Type*& t) {
// An array of Types: deserialize n, then
// "t" is allocated and deserialized
	this->serialize(n);
	t = new Type[n];
	for (size_t i=0; i<n; i++)
		this->serialize(t[i]);
}

template<typename Type>
void
Receiver::
serialize(std::vector<Type>& t) {
// A vector of Types: deserialize n (t.size()), push_back
// each element of "t" (assumes "t" is empty)
	size_t n;
	this->serialize(n);
	Type ti;
	for (size_t i=0; i<n; i++) {
		this->serialize(ti);
		t.push_back(ti);
	}
}

template<typename Type>
void
Receiver::
serialize(std::deque<Type>& t) {
	size_t n;
	this->serialize(n);
	Type ti;
	for (size_t i=0; i<n; i++) {
		this->serialize(ti);
		t.push_back(ti);
	}
}


template<typename Type>
void
Receiver::
serialize(std::map<std::string,Type>& t) {
	size_t n;
	this->serialize(n);
	std::map<std::string,Type> rval;
	for (auto& ti : t) {
		Type tt;
		std::string key;
		this->serialize(key);
		this->serialize(tt);
		t.insert(make_pair(key,tt));
	}
}

using FioGetFcn = Fio* (*)(Receiver&);

class Fio {
public:
	Fio() {}

	// things which all Fio classes must have.
	// XXX The virtual destructors
	// must have the body implemented in the .c file (even if it's empty)
	// to avoid "undefined reference to 'vtable for ...'" loader errors
	virtual ~Fio() = default;
	virtual std::string vid() const { return std::string("Fio"); }
	virtual std::string summary() const { return vid(); }
	virtual void put (Sender&) const { throw std::runtime_error("Fio::put called"); }
	// get() is not virtual, but static because it must be
	// callable without and object, and it must be registered
	// by calling Fio::register_class
	// static Fio* get(Receiver& s)
	virtual Fio* clone() const { std::string msg(vid());
		msg += "::clone() called"; throw std::runtime_error(msg); }

	static bool register_class (std::string, Fio* (*)(Receiver&));

	// Special static Fio-only versions of get/put: get reads
	// an object id, looks up the registered get function for
	// that object, calls the get function, and returns an Fio*
	static Fio* get(Receiver& s);
	// put just writes m.vid, then calls m.put
	static void put(Sender& s, const Fio& m);
}; // Fio


// vectors of Fio*, e.g. vector<Pz*>
void getVFio (Receiver& s, std::vector<Fio*>& t);
void putVFio (Sender& s, const std::vector<Fio*>& t);


// fetch/store vector<Type> with Type a basic datatype (i.e. there is
// a Sender::serialize(Type))/Receiver::serialize(Type)
namespace fio {

template<typename Type>
void
fetch(const std::string& mid, std::vector<Type>& rval) {
	// the expected object name
	std::string exn{vastr("vector<",typeid(Type).name(),'>')};

	// open the file
	Receiver s(mid);
	// get the type name
	std::string name;
	s.serialize(name);
	if (name != exn)
		throw std::runtime_error(vastr("fetch(",mid,") expecting ",exn,", got ",name));
	// and serialize rval
	s.serialize(rval);
}

template<typename Type>
void
store(const std::string& mid, const std::vector<Type>& a) {
	// the object name
	std::string name{std::string("vector<")+typeid(Type).name()+std::string(">")};

	// open the file
	Sender s(mid);
	// put the type name
	s.serialize(name);
	// serialize the object
	s.serialize(a);
}

// list all the files in Ftmp matching regex "re"
std::vector<std::string>
catalog(const std::string& re="", bool sum=false);

// save some matrices in file "path" (compressed tar file)
std::vector<std::string> save (const std::vector<std::string>& rxname,
		const std::string& path);
// restore matrices from "path" to the data directory
std::vector<std::string> restore (const std::string& path);

template<typename Type>
bool
is_type(const std::string& name) {
	std::string id{name};
	Receiver s(id);
	s.serialize(id);
	if (id == Type::id())
		return true;
	return false;
}

// shorten a path by replacing HOME with ~, and ?
std::string shortenpath(const std::string& path);

} // namespace fio

// functions and Classes for various types of I/O
// The classes are for using "Resource acquisition is initialization"
// (RAII) Stroustrup V4 p 356
std::string get_cwd();
std::string getfroot();
std::string getftmp();
std::string get_datadir();
std::string getfwd();

char* errDir();
void rmdirTree (std::string const& path);
std::vector<std::string> ls (std::string const& path);

// RAII class to open a file using open()
class Fd {
private:
	int fd;
public:
	Fd(const std::string& file, int flags, mode_t mode=0);
	~Fd();
	int operator()() { return fd; }
};

// RAII class to open a file using fopen()
class File {
private:
	FILE* fp;
public:
	File(const std::string& file, const std::string& mode);
	~File();
	FILE* operator()() { return fp; }
};

// RAII class to open a POSIX message queue using mq_open()
// Reference:
// Stevens, W.R., "UNIX Network Programming Vol. 2, Interprocess Communication"
class Msgq {
private:
	std::string name_;
	mqd_t des_;
	bool creator{false};
public:
	// the 2 optional arguments (mode,attr) signify that the queue
	// is to be created and unlinked in the destructor.
	// "flags" is one of O_RDONLY, O_WRONLY or O_RDWR possibly bitwise ORed
	// with O_CREAT, O_EXCL, and O_NONBLOCK
	Msgq(const std::string& name, int flags,
			mode_t mode=0, struct mq_attr* attr = nullptr);
	~Msgq();
	mqd_t operator()() { return des_; }	// XXX unnecessary?

	// send vectors and strings
	template<typename T>
	void send(const std::vector<T>& v) {
		int stat = mq_send(des_, reinterpret_cast<const char*>(v.data()),
				v.size()*sizeof(T), 0);
		if (stat == -1)
			throw std::system_error(errno, std::generic_category(),"mq_send");
	}
	void send(const std::string& v);

	// receive vectors and strings
	template<typename T>
	void receive(std::vector<T>& v) {
		struct mq_attr attr;
		mq_getattr(des_, &attr);
		long msgsize = attr.mq_msgsize;
		std::vector<char> msg(msgsize);
		ssize_t nbytes = mq_receive(des_, msg.data(), msgsize, nullptr);
		if (nbytes == -1)
			throw std::system_error(errno, std::generic_category(),"mq_receive");
		if (nbytes%sizeof(T) != 0)
			throw std::runtime_error(vastr("received ",nbytes,
						" bytes, not a multiple of ", sizeof(T)));
		size_t len = nbytes/sizeof(T);
		v.resize(len);
		T* tp = reinterpret_cast<T*>(msg.data());
		for (size_t i=0; i<len; i++)
			v[i] = tp[i];
	}
	void receive(std::string& v);
	// notify me when a message is added to an empty queue
	void notify(int signo = SIGUSR1);
};

// RAII class to open a directory using opendir()
class Fdir {
private:
	DIR* dir;
public:
	Fdir(const std::string& path);
	~Fdir();
	DIR* operator()() { return dir; }
};

// RAII class to change directories while an instance is in scope
class Chdir {
	// open the starting directory to make returning easy
	Fd cwdfd;
	std::string cwd;		// where we changed from
public:
	Chdir(const std::string& dir);
	~Chdir();
	std::string from() { return cwd; }
};

// RAII class to open a pipe to a command
class Fpipe {
private:
	FILE* fp;
	std::string comm;
public:
	Fpipe(const std::string& cmd, const std::string& mode="r");
	~Fpipe();
	FILE* operator()() { return fp; }
};

// Ftmpdir creates a temporary directory that the fio functions
// use for storing data. The path to this directory is set in
// the environment as "FTMP"
// The only constructor takes one optional argument, the path
// to use as the temporary directory. If this argument is not
// included a new directory under /tmp is created
class Ftmpdir {
private:
	std::string _tmpdir;
public:
	Ftmpdir(const std::string& path="");
	~Ftmpdir();
}; // Ftmpdir

#endif	// FIOb_h
