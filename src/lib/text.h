//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

// utilities for manipulating text

#ifndef Text_h
#define Text_h

#include <algorithm>
#include <complex>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>


// get/put a string into the environment
std::string
getEnv (std::string const& var);
void
putEnv (std::string const& var);
// expand environment variables in a string, e.g. $HOME/local
std::string
expenv(const std::string& str);

// return a string of char that are in a but not b
std::string
string_diff(const std::string& a, const std::string& b);
// return a summary of the largest diffs between a & b
std::string
vector_diff(const std::vector<double>& a, const std::vector<double>& b);

// find the delim that matches the (opening or closing) delimeter at str[start]
std::string::size_type find_delim(std::string const& str,
		std::string::size_type start);

// returns the string delimited by "open" and "close" in string "str"
// where "open" is str[start]. Also returns "end" which is the
// index in str of character "close".
std::string delimitedString(std::string const& str, char open, char close,
		std::string::size_type start, std::string::size_type& end);


// reverse substr: last n char
std::string rsubstr(std::string const& s, std::string::size_type n);

// replace all instances of char 'c' with 't'
std::string replace_char(std::string const& s, char c, char t);

// stripquotes() works a little different than stripDelim():
// it only considers the first and last char when looking for quotes
void stripquotes(std::string& s);

// Adds '\' before each quote
std::string escapeQuotes(std::string const& opt);

// compare 2 strings up to length of b, ignore case & leading and trailing
// whitespace, try converting to doubles, comparing up to nsig significant
// figures
bool compare(std::string const& a, std::string const& b, int nsig=6);

// add an element to a vector if it is not already there
template<typename T>
void
add2vector(const T& s, std::vector<T>& v) {
	if (std::find(v.begin(), v.end(), s) == v.end())
		v.push_back(s);
}

// accumulate phrases in a message, e.g. "a=b & b!=c & c==d"
void add2msg(std::string& msg, const std::string& phrase);

// Tokenize a string
// 3 versions of string2tok: the first takes a single-character delimiter
// and no quotes; the second takes a string of delimiter characters, and
// the third has a third argument specifying what pairs of characters are
// to be considered "quotes". Quoted strings are considered a token without
// regard to "delim" characters inside the quotes.
std::vector<std::string>
string2tok (std::string const& str, char delim);

std::vector<std::string>
string2tok (std::string const& str, std::string const& delim);

std::vector<std::string>
string2tok (std::string const& str, std::string const& delim,
		std::vector<std::string> const& quotes);

std::vector<std::string>
cstr2tok(int nstr, char* cstr[], int start=0);

std::string stringLower(std::string const& a);
std::string stringUpper(std::string const& a);

std::string stripwhitespace(std::string const& s);

std::string roff (std::string const& s, unsigned int indent, unsigned int linelen);

std::string regex_double(); // regex string for detecting a double

// convert a string to int/double/complex
bool str2int(const std::string& s, int& x);
bool str2int(const std::string& s, int& x,
		std::string::size_type start, std::string::size_type& end);

bool str2double(const std::string& s, double& x);
bool str2double(const std::string& s, double& x, size_t start, size_t& end);
bool str2complex (const std::string& s, std::complex<double>& x);
bool str2complex (const std::string& s, std::complex<double>& x,
		std::string::size_type start, std::string::size_type& end);


std::string stringJustify (size_t width, std::string const& s, std::string const& just);
std::string stringCenter (size_t width, std::string const& s);
std::string stringIndent (size_t nchar, std::string const& s);

std::string stringTimeDate ();

// returns "path" stripped of the directory and extension
std::string stringBasename (std::string const& path);
// returns the directory portion of "path"
std::string stringDirname (std::string const& path);

int page_width(int newpw=-1);

std::string& separator();

void dos2Unix(std::string& str);

// get a piece of information about the kernel from uname(2)
// Valid options for "what":
//   sysname     Operating system name (e.g., "Linux")
//   nodename    Name within "some implementation-defined network"
//   release     Operating system release (e.g., "2.6.28") 
//   version     Operating system version
//   machine     Hardware identifier
//   domainname  NIS or YP domain name
std::string get_uname(const std::string& what="");

// get a string with a description of wait() status code, or an
// empty string if normal exit
std::string exitstatus(int status);

// an Ascii object opens a file and has an operator() which reads
// lines:
// 	Ascii file("/tmp/junk");
// 	string line;
// 	while(true) {
//			file(line, "//");
//			if (line.empty()) break;
//			   ...
//		}
class Ascii {
	std::ifstream ifs;
	int lineno{0};
	std::string path;
public:
	Ascii() {}
	Ascii(Ascii&& ap) = default;	// move constructor
	Ascii& operator=(Ascii&& rhs) = default;	// move assignment op
	Ascii(const std::string& name);	// constructor
	bool operator()(std::string& line, const std::string& comment="",
			bool blanks=true);
	bool good();
	int line_number() { return lineno; }
};

/*------------------------------------------------------------------
 * User-defined I/O manipulator for specifying formats for floating-
 * point numbers.
 * Taken from Stroustrop 4th Ed. 38.4.5.3
 * Usage:
 * Form gen4(4); // general format, precision 4
 *
 * void f(double d) {
 *    Form sci8 = gen4;
 *    sci8.scientific().precision(8);  // scientific format, precision 8
 *    cout << d << ' ' << gen4(d) << ' ' << sci8(d) << ' ' << d << '\n';
 * }
 *
 * A call f(1234.56789) writes
 *    1234.57 1235 1.23456789e+03 1234.57
 * Note how the use of a Form doesn't affect the state of the stream
 * so that the last output of d has the same default format as the first.
 *------------------------------------------------------------------*/

#include <iosfwd>
#include <iostream>

class Bound_form;

class Form {
	friend std::ostream& operator<<(std::ostream&, Bound_form const&);

	int prc;  // precision
	int wdt;  // width, 0 means as wide as necessary
	std::ios::fmtflags fmt;  // general, scientific, or fixed (21.4.3)
	char fillchar;
public:
	explicit Form(int p=6) : prc(p) {
		fmt = std::ios::scientific;  // general
		wdt = 0;  // wide as necessary
		fillchar = ' ';
	}
	Bound_form operator() (double d) const;   // make a Bound_form for *this and d

	Form& scientific() { fmt = std::ios::scientific; return *this; }
	Form& fixed() { fmt = std::ios::fixed; return *this; }
	Form& general() { fmt = std::ios::showpoint; return *this; }
	Form& point() { fmt |= std::ios::showpoint; return *this; }

	Form& uppercase();
	Form& lowercase();
	Form& precision(int p) {prc = p; return *this; }

	Form& width(int w) { wdt = w; return *this; }
	Form& fill(char t) { fillchar = t; return *this; }

	Form& plus(bool b = true);    // explicit plus sign
	Form& trailing_zeros(bool b = true);
};

struct Bound_form {
	Form const& f;
	double val;
	Bound_form(const Form& ff, double v) : f(ff), val(v) {}
};

std::ostream& operator<<(std::ostream&, Bound_form const&);

#endif // Text_h
