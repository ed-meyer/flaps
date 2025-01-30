//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-


#include <cstdlib>
#include <cstring>
#include <fstream>
#include <pwd.h>  // for HOME
#include <sstream>
#include <string>
#include <sys/stat.h> 		// for stat()
#include <sys/utsname.h>	// for uname()
#include <unistd.h>			// fileno

#include "csignals.h"
#include "fptype.h"
#include "matrix.h"
#include "Regex.h"
#include "text.h"
#include "trace.h"
#include "vastr.h"

using namespace std;

string
getEnv (string const& var) {
/*
 * get the value of environment variable "var"
 * if it exists; otherwise an empty string
 * is returned. This avoids the problem of initializing
 * a string from getenv() if it returns a null pointer.
 */
	char const* val = getenv(var.c_str());
	string empty;
	if (val != nullptr)
		return string(val);
	return empty;
}

void
putEnv (string const& var) {
// Just like std unix putenv() except that this one takes a
// constant string argument, and throws an exception if there
// is an error; more importantly the input string does not
// need to be static - it is copied to a new char array from the heap
	T_(Trace trc(1,"putEnv", ": ",var);)
	// putenv wants a pointer to permanent memory
	char* ev = new char[var.size()+1];
	strcpy (ev, var.c_str());
	if (putenv(ev) < 0)
		throw std::system_error(errno, std::generic_category(),
			vastr("could not put ", var, " into the environment"));
}

string
expenv(const string& str) {
/*------------------------------------------------------------------
 * Expand any environment variables (which are not enclosed
 * in single quotes) in a string. Environment variables must
 * be enclosed in curly braces if the following character
 * is a letter or digit; e.g.
 *     ${HOME}junk  ok - HOME is recognized
 *     $HOMEjunk    error - HOME is not recognized
 *     $HOME_junk   ok - HOME is recognized
 * Tilde is recognized as the current $HOME
 * The original string is not modified.
 * input:
 *    str      string with (possibly) environment variable(s)
 * returns string with environment variables expanded.
 *------------------------------------------------------------------*/
	T_(Trace trc(3,"expenv");)
	string rval;
	string::size_type i = 0;

   T_(trc.dprint("expanding <",str,">");)
	// quick return if str does not contain $ or ~
	string::size_type dollarpos = str.find('$');
	string::size_type tildepos = str.find('~');
	if (dollarpos == string::npos && tildepos == string::npos) {
		T_(trc.dprint("returning unchanged");)
		return str;
	}

	while (i<str.size()) {
		if (str[i] == '\'') {
			rval.push_back(str[i++]);
			while(i<str.size()) {
				rval.push_back(str[i]);
				if (str[i++] == '\'') {
					break;
				}
			}
			if (i >= str.size())
				break;
		}
		if (str[i] == '~') {
			if (str[i+1] == '/') {
				char const* home = getenv("HOME");
				while (*home != 0)
					rval.push_back(*home++);
				i++;
			} else if (isalpha(str[i+1])) {
				string::size_type idx = str.find('/', i);
				string::size_type last;
				string user;
				if (idx != string::npos) {
					user = str.substr(i+1, idx-i-1);
					last = idx;
				} else {
					user = str.substr(i+1);
					last = str.size();
				}
				struct passwd* pw = getpwnam(user.c_str());
				if (pw && pw->pw_dir) {
					char const* cp = pw->pw_dir;
					while(*cp != 0)
						rval.push_back(*cp++);
					i = last;
				} else {
					rval.push_back(str[i++]);
				}
			} else {
				rval.push_back(str[i++]);
			}
		} else if (str[i] != '$') {
			rval.push_back(str[i++]);
		} else {
			string ev;
			i++;
			if (str[i] == '{') {
				i++;
				while (i<str.size() && str[i] != '}')
					ev.push_back(str[i++]);
				i++;	// eat the closing brace
			} else {
				// move to the char which ends this env var: !(letter or digit)
				while (i<str.size() && std::isalnum(str[i]))
					ev.push_back(str[i++]);
			}
			T_(trc.dprint("got ev <",ev,">");)
			char* val = getenv(ev.c_str());
			// if this env var is undefined leave blank
			if (val) {
				rval += val;
				T_(trc.dprint("appended <",val,">");)
			}
		}
	}
	T_(trc.dprint("returning ",rval);)
	return rval;
}

string
string_diff(const string& a, const string& b) {
// returns a string containing the characters of "a" that differ from "b"
// Method:
//   "a" and "b" are tokenized using "./:[a-zA-Z][0-9]" as delimiters
//   then each token is compared
	string rval;
	for (size_t i=0; i<a.size(); i++) {
		if (a[i] != b[i])
			rval.push_back(a[i]);
	}
	return rval;
}

string
vector_diff(const vector<double>& a, const vector<double>& b) {
// returns a string summarizing the largest abs diffs between "a" and "b"
	size_t n = std::min(a.size(), b.size());
	vector<double> diff(n, 0.0);
	for (size_t i=0; i<n; i++)
		diff[i] = std::abs(a[i] - b[i]);
	return flaps::summarize<double>(diff, 80);
}


string::size_type
find_delim(string const& str, string::size_type start) {
// find the delimiter matching the character at position "start"
// in string "str", which may enclose 1 or more pairs of the
// same open/close delimiters; e.g. pz{gaf{rf=0},...}
// The character at "start" must be one of {}[]()<>, that is
// it may be either the open or close character of a pair; if it
// is a closer, searching will be backwards to find the opener.
	string openers("{([<");
	string closers("})]>");
	char startchar{str[start]};
	char endchar;
	bool findclose{false}; // search for opener or closer?
	ostringstream os;

	// if the start char is a closer, the end char
	// is the corresponding opener and vice-versa
	string::size_type pos = openers.find(startchar);
	if (pos == string::npos) {
		pos = closers.find(startchar);
		if (pos == string::npos) {
			throw runtime_error(vastr("find_delim called with \"",
					str, "\" and start = ", start));
		}
		endchar = openers[pos];
	} else {
		findclose = true;
		endchar = closers[pos];
	}

	string::size_type end = start;
	char next = str[end];
	char openclose[3];
	openclose[0] = startchar;
	openclose[1] = endchar;
	openclose[2] = 0;

	// search forwards or backwards depending on whether the
	// start char is an opener or a closer
	if (findclose) {
		while (end < str.size()-1 && next != endchar) {
			end = str.find_first_of(openclose, end+1);
			if (end == string::npos)
				return end;
			next = str[end];
			if (next == startchar) {
				end = find_delim(str, end);
				if (end == string::npos)
					return end;
			}
		}
	} else {
		while (end > 0 && next != endchar) {
			end = str.find_last_of(openclose, end-1);
			if (end == string::npos)
				return end;
			next = str[end];
			if (next == startchar) {
				end = find_delim(str, end);
				if (end == string::npos)
					return end;
			}
		}
	}
	return end;
}


string
delimitedString(string const& str, char open, char close,
		string::size_type start, string::size_type& end) {
// extract a string delimited by characters "open" and "close".
// Returns the enclosed string without the delimiters and:
//   end    index in "str" of the character "close"
	string rval;
	if (str[start] != open)
		return rval;
	string::size_type first = start + 1;
	if (first > str.size())
		return rval;
	end = str.find(close, first);
	if (end == string::npos)
		return rval;
	// watch out for embedded open-close
	string::size_type emb = str.find(open, first);
	if (emb != string::npos && emb < end)
		end = str.find(close, end+1);
	if (end == string::npos)
		throw runtime_error(vastr("unmatched ",open," in ",str));
	rval = str.substr(first, end-first);
	return rval;
}

string
rsubstr(string const& s, string::size_type len) {
// Returns the substring which is the last "len"
// characters of "s"
	if (len >= s.size())
		return s;
	if (len == 0)
		return "";
	return s.substr(s.size()-len);
}

string
replace_char(string const& s, char c, char t) {
// Replace all instances of character 'c' with 't'
// in string "s"
	string rval(s);
	// use a lamba for the action
	std::for_each(rval.begin(), rval.end(), [&](char& f){ if (f == c) f = t; });
	return rval;
}

string&
separator() {
// returns a (reference to static) string of dashes with width pagewidth
	static string rval;
	if (rval.empty())
		rval = string(page_width(), '-');
	return rval;
}

void
stripquotes(string& s) {
// Remove single or double quotes surrounding a
// string.
// NOTE: all quotes (inner and outer) will be removed
//       this is a change in behavior (7/7/08)
//       XXX changing back to removing ONLY outer quotes (6/21/14)
	string rval;
	string::size_type end;
	if (s[0] == '\'') 
		rval = delimitedString(s, '\'', '\'', 0, end);
	else if (s[0] == '\"')
		rval = delimitedString(s, '\"', '\"', 0, end);
	else
		rval = s;
	s = rval;
}

string
escapeQuotes(string const& str) {
/*------------------------------------------------------------------
 * Put an escape char (\) before each single or double quote
 * that is not already escaped
 *------------------------------------------------------------------*/
	string::size_type qc = str.find_first_of("\'\"");
	string::size_type current = 0;
	ostringstream os;

	if (qc == string::npos) {
		return str;
	}

	while (qc < str.size()) {
		if (qc > 0 && str[qc-1] == '\\') {
			current = qc + 1;
		} else {
			os << str.substr(current, qc-current) << '\\' << str[qc];
			current = qc + 1;
		}
		if (current >= str.size())
			break;
		qc = str.find_first_of("\'\"", current);
		if (qc == string::npos) {
			os << str.substr(current);
			break;
		}
	}
	return os.str();
}

bool
compare(const string& a, const string& b, int nsig) {
// Compare 2 strings, either of which may be a regular-expression, ignoring case.
// Strip leading & trailing blanks and only compare up to the
// number of characters in "b". 
// For example:
// output, out: true
// i, in:       false
// So if you want to check an option and it's full
// name is "input" but you want to allow "i", "in",
// "inp", etc, put the Tok in the second arg:
//     compare("input", pref->srhs)
	T_(Trace trc(3,"compare");)
	string aa = stripwhitespace(a);
	string bb = stripwhitespace(b);
	T_(trc.dprint("\"",aa,"\" =? \"",bb,"\"");)
	if (aa == bb)
		return true;
	// b may be shorter than a but not longer
	if (bb.size() < aa.size() && bb == aa.substr(0,bb.size()))
		return true;
	// one or the other may be a regex
	if (regex_match(aa, make_regex_ic(bb)))
		return true;
	if (regex_match(bb, make_regex_ic(aa)))
			return true;
	// check for both ints or doubles
	double x, y;
	if (str2double(aa, x) && str2double(bb, y) && is_equal(x, y, nsig))
		return true;
	return false;
}

void
add2msg(std::string& msg, const std::string& phrase) {
	if (!msg.empty()) {
		msg += " & ";
		msg += phrase;
	} else
		msg = phrase;
}

string
stringLower(string const& a) {
	string rval(a);
	for (size_t i=0; i<rval.size(); i++)
		rval[i] = tolower(rval[i]);
	return rval;
}

string
stringUpper(string const& a) {
	string rval(a);
	for (size_t i=0; i<rval.size(); i++)
		rval[i] = toupper(rval[i]);
	return rval;
}

string
stripwhitespace(string const& s) {
// Strip leading and trailing whitespace from a string
	T_(Trace trc(3,"stripwhitespace "," <",s,">");)
	string rval;
	string::size_type first = s.find_first_not_of(" \t\n");
	if (first != string::npos)
		rval = s.substr(first);
	else
		return rval;

	string::size_type last = rval.find_last_not_of(" \t\n");
	if (last != string::npos)
		rval = rval.substr(0, last+1);

	T_(trc.dprint("returning <",rval,">");)
	return rval;
}

string
stripDelim (string const& s, char open, char close) {
/*
 * Strip the leading and trailing delimiters from a string;
 * all char up to (and including) the "open" character will
 * be eliminated and all characters from (and including) the
 * "close" char to the end of the string will be eliminated
 * in the returned string
 */
	string::size_type first = s.find(open);
	string::size_type last = s.rfind(close);
	if (first == string::npos || last == string::npos)
		return s;
	if (first == last)
		return s;
	return s.substr(first+1, last-first-1);
}

vector<string>
string2tok (string const& str, char delim) {
// Simplest string2tok: single-character delimiter and no quotes
	T_(Trace trc(3,"string2tok(char)");)
	vector<string> rval;
	string::size_type first = 0, last;
	
	T_(trc.dprint("str<",str,"> delim<",delim,'>');)

	// is whitespace a delimiter?
	string white(" \t\n");
	bool skipwhite = (white.find_first_of(delim) == string::npos);

	size_t nch = str.size();
	while (first <= nch) {
		string tok;
		// skip whitespace
		if (skipwhite)
			first = str.find_first_not_of(white, first);
		if (first == string::npos)
			break;
		// advance "last" until a delim
		last = first;
		while (last < nch) {
			if (last >= str.size())
				break;
			// is this char a delim?
			if (str[last] == delim) {
				while (str[last+1] == delim) last++;
				break;
			}
			last++;
		}
		// strip leading & trailing whitespace
		//!! tok = stripwhitespace(str.substr(first, last-first));
		string::size_type lst = last;
		while (str[lst-1] == ' ' || str[lst-1] == '\t') lst--;
		tok = str.substr(first, lst-first);
		if (!tok.empty()) {
			T_(trc.dprint("got <",tok,">");)
			rval.push_back(tok);
		}
		last++;
		first = last;
	}
	return rval;
}

vector<string>
string2tok (string const& str, string const& delim, vector<string> const& quotes) {
// Break a string up into tokens delimited by any of the characters in
// "delim". Quoted strings are not tokenized; that is, "delim" characters
// are ignored inside quotes.
//    quotes     vector<string> each string is exactly two characters: the
//               opening and closing character.
	// use a simpler version if possible
	if (delim.size() == 1 && quotes.empty()) {
		char cdelim{delim[0]};
		return string2tok(str, cdelim);
	}
	T_(Trace trc(3,"string2tok(quotes)");)
	vector<string> rval;
	string::size_type first = 0, last;
	
	T_(trc.dprint("str<",str,"> delim<",delim,"> ",quotes.size()," quotes");)

	// is whitespace a delimiter?
	string white(" \t\n");
	bool skipwhite = (delim.find_first_of(white) == string::npos);

	size_t nch = str.size();
	while (first <= nch) {
		string tok;
		// skip whitespace
		if (skipwhite)
			first = str.find_first_not_of(white, first);
		if (first == string::npos)
			break;
		// advance "last" until a delim
		last = first;
		while (last < nch) {
			// check for a quote
			for (auto& quo : quotes) {
				if (str[last] == quo[0]) {
					last++;
					last = str.find(quo[1], last);
					if (last == string::npos) {
						throw runtime_error(vastr("no matching \"", quo[1],
								" in ", str));
					}
					last++;
					if (last == str.size())
						break;
				}
			}
			if (last >= str.size())
				break;
			// is this char a delim?
			if (delim.find(str[last]) != string::npos)
				break;
			last++;
		}
		// strip leading & trailing whitespace
		//!! tok = stripwhitespace(str.substr(first, last-first));
		string::size_type lst = last;
		while (str[lst-1] == ' ' || str[lst-1] == '\t') lst--;
		tok = str.substr(first, lst-first);
		if (!tok.empty()) {
			T_(trc.dprint("got <",tok,">");)
			rval.push_back(tok);
		}
		last++;
		first = last;
	}
	return rval;
}

vector<string>
string2tok (string const& str, string const& delim) {
// simpler version of string2tok: no quoted string handling

	if (delim.size() == 1) {
		char cdelim{delim[0]};
		return string2tok(str, cdelim);
	}
	vector<string> quotes;
	return string2tok(str, delim, quotes);
}

vector<string>
cstr2tok(int nstr, char* cstr[], int start) {
// Convert an array of C-strings to a vector of strings
	vector<string> rval;
	for (int i=start; i<nstr; i++) {
		rval.push_back(string(cstr[i]));
	}
	return rval;
}

string
roff (string const& s, unsigned int indent, unsigned int linelen) {
// do some simple formatting of a string:
//   - indent it by "indent" char
//   - restrict lines to no more than "linelen" char
			
	if (s.empty())
		return s;

	ostringstream os;
	string ind;
	if (indent > 0)
		ind = string(indent, ' ');

	string delim(" \t");
	vector<string> quotes;
	quotes.push_back("\'\'");
	quotes.push_back("\"\"");

	string::size_type start = 0;
	string theline;
	string::size_type idx = 0;
	// break the string into lines (theline), then format each line
	do {
		idx = s.find('\n', start);
		if (idx == string::npos)
			theline = s.substr(start);
		else {
			theline = s.substr(start, idx-start+1);
			start = idx+1;
		}

		// break the line up into "tokens": blank-separated words
		// or quote-delimited strings
		vector<string> tok = string2tok (theline, delim, quotes);
		string line;

		for (size_t i=0; i<tok.size(); i++) {
			if (line.size() + tok[i].size() + 1 >= linelen) {
				os << ind << line << std::endl;
				line.clear();
			} else {
				if (!line.empty())
					line += ' ';
			}
			line += tok[i];
		}
		if (!line.empty()) {
			os << ind << line;
		}
	} while (idx != string::npos && start < s.size());
	return os.str();
}

bool
str2int(const string& s, int& x, size_t start, size_t& end) {
// convert a portion of a string to int, starting at s[start].
// "end" is output: one past the last character in the int (0b).
// Returns false if s[start] does not begin an int, true othersize.
	static regex* re{nullptr};
	if (re == nullptr) {
		re = new regex("^[ \t]*[+-]?[0-9]+");
	}
	smatch mch;
	string sub;
	if (start > 0) {
		sub = s.substr(start);
		if (!regex_search(sub, mch, *re))
			return false;
	} else if (!regex_search(s, mch, *re)) {
		return false;
	}
	// convert the string to int with std::stod
	size_t nc;
	string dstr{mch[0].str()};
	x = std::stod(dstr, &nc);
	end = start + nc;
	return true;
}

bool
str2int(const string& s, int& x) {
	size_t start{0};
	size_t end{s.size()};
	return str2int(s, x, start, end);
}

string
regex_double() {
// returns a string which can be used in a regex to match a double
	return "[+-]?([0-9]+[.]?[0-9]*|[0-9]*[.]?[0-9]+)([dDeE][+-]?[0-9]+)?";
}

bool
str2double(const string& s, double& x, size_t start, size_t& end) {
// convert a string (s) starting at index "start" to a double (x),
// and set "end" to one past the last character of the double string.
// Returns true if the conversion was successful, false if
// the string does not represent a double and x is unchanged
	static regex* re{nullptr};
	if (re == nullptr) {
		re = new regex("^[ \t]*[+-]?([0-9]+[.]?[0-9]*|[0-9]*[.]?[0-9]+)"
			"([dDeE][+-]?[0-9]+)?[ \t\n]*");
	}
	smatch mch;
	string sub;
	if (start > 0) {
		sub = s.substr(start);
		if (!regex_search(sub, mch, *re))
			return false;
	} else if (!regex_search(s, mch, *re)) {
		return false;
	}

	size_t nc;
	string dstr{mch[0].str()};
	x = std::stod(dstr, &nc);
	if (nc == dstr.size())
		end = string::npos;
	else
		end = start + nc;
	return true;
}

bool
str2double(const string& s, double& x) {
// convert a string (s) to a double (x), return true if the conversion was
// successful, false if the string does not represent a double and x is unchanged
	size_t start{0};
	size_t end{0};
	return str2double(s, x, start, end);
}

bool
str2complex (const string& s, complex<double>& x) {
	string::size_type start{0};
	string::size_type end{string::npos};
	return str2complex(s, x, start, end);
}

bool
str2complex (string const& s, complex<double>& x, string::size_type start,
		string::size_type& end) {
/*------------------------------------------------------------------
 * Convert a string to a complex<double>. The string may use d, D, e,
 * or E for scientific notation (unlike strtod).
 * The string may have trailing non-float characters which
 * are ignored.
 * If the string is not a float or if there was an error
 * in converting the string a runtime_error exception is thrown.
 * Modifys "end": last character in the floating-point number
 * plus one.
 * The legal format for a complex scalar:
 *    (x + iy)
 * Note: even though a single float is technically also a
 * complex, we do not allow that here because there are times
 * (e.g. alge) when it is important to distinguish.
 *------------------------------------------------------------------*/
	T_(Trace trc(3,"str2complex");)

	T_(trc.dprint("<",s,"> start ",start);)

	if (start >= s.size()) {
		ostringstream os;
		os << "request to start at char " << start;
		if (s.empty())
			os << " in an empty string";
		else
			os << " in \"" << s << '\"';
		T_(trc.dprint("throwing exception: ",os.str());)
		throw runtime_error(os.str());
	}
	/*
	 * a complex scalar may begin with open paren '('
	 */
	if (s[start] == '(')
		start++;

	// end is either the first ) or npos
	end = s.find_first_of(")", start);

	// Grab the relevant part - enclosed in parentheses
	string str = s.substr(start, end-start);

	T_(trc.dprint("testing substring ",str);)

	smatch mch;
	if (regex_match(str, mch, make_regex(R"(\s*([^-+ ]*)\s*([-+])\s*[ij](.*))"))) {
		double re{0.0};
		double im{0.0};
		if (mch.size() < 4) {
			T_(trc.dprint("returning false: regex does not match");)
			return false;
		}
		if (str2double(mch[1],re) && str2double(mch[3],im)) {
			if (mch[2] == "-")
				im = -im;
			x = complex<double>(re,im);
			T_(trc.dprint("returning true: ",x);)
			return true;
		}
	}
	T_(trc.dprint("returning false");)
	return false;
}

string
stringJustify (size_t width, string const& s, string const& just) {
// create a string of width "width" containing string "s" placed
// according to "just":
//   "left"    left-justify
//   "right"   right-justify
//   "center"  centered
// XXX this could replace stringCenter()
	size_t i;
	size_t n = s.size();
	if (n > width)
		n = width;
	string rval(width, ' ');

	if (just == "left") {
		for (i=0; i<n; i++)
			rval[i] = s[i];
	} else if (just == "right") {
		for (i=0; i<n; i++)
			rval[width-1-i] = s[n-1-i];
	} else if (just == "center") {
		int npad = (width - s.size())/2;
		if (npad < 0) npad = 0;
		for (i=0; i<n; i++)
			rval[npad+i] = s[i];
	}
	return rval;
}

string
stringCenter (size_t width, string const& s) {
	if (width < s.size() + 2)
		return s;
	int npad = (width - s.size())/2;
	string pad(npad, ' ');
	string rval = pad + s;
	return rval;
}

string
stringIndent (size_t nchar, const string& s) {
	string rval;
	string indent(nchar, ' ');
	if (s[0] != '\n')
		rval = indent;
	for (size_t i=0; i<s.size(); i++) {
		rval += s[i];
		if (s[i] == '\n' && i != s.size()-1)
			rval += indent;
	}
	return rval;
}

string
stringTimeDate () {
/*------------------------------------------------------------------
 * Returns string containing the current local time
 * & date (without the newline added by asctime())
 *------------------------------------------------------------------*/
	time_t t = time(0);
	string rval(asctime(localtime(&t)));

	string::size_type idx = rval.find('\n');
	if (idx != string::npos)
		rval = rval.substr(0,idx);
	return rval;
}

string
stringBasename (string const& path) {
/*------------------------------------------------------------------
 * Returns the "basename" portion of a path: e.g. the
 * basename of /home/git/junk.trash is "junk"
 *------------------------------------------------------------------*/
	string rval;
	auto start = path.rfind('/');
	if (start == string::npos)
		start = 0;
	else
		start++;
	auto end = path.rfind('.');
	if (end == string::npos)
		rval = path.substr(start);
	else
		rval = path.substr(start, end-start);
	return rval;
}

string
stringDirname (string const& path) {
/*------------------------------------------------------------------
 * Returns the "dirname" portion of a path: e.g. the
 * dirname of /u/ba/eem2314/junk is "/u/ba/eem2314"
 *------------------------------------------------------------------*/
	string::size_type ipos = path.rfind('/');
	if (ipos == string::npos)
		return ".";
	return path.substr(0,ipos);
}

int
page_width(int newpw) {
// Returns the current desired width of text going to output.
// The pagewidth can be set either by setting environment variable
// PAGEWIDTH or calling this function with an integer argument, e.g.
//     page_width(380);
	T_(Trace trc(3,"page_width");)
	static int pagewidth{0};  // default
	int defaultpw{160};
	char *s;

	// check environment
	if (pagewidth == 0) {
		if ((s = getenv("PAGEWIDTH")) != nullptr) {
			pagewidth = atoi(s);
			if (pagewidth < 10 || pagewidth > 10000) {
				flaps::warning("illegal pagewidth specified with PAGEWIDTH "
					"environment variable (",s,")");
				pagewidth = defaultpw;
			}
			flaps::info("set page width to ",pagewidth," from the environment");
		} else {
			pagewidth = defaultpw;
			T_(trc.dprint("setting default page width ",pagewidth);)
		}
	}

	// arg included?
	int rval = pagewidth;
	if (newpw > 0) {
		pagewidth = newpw;
		T_(trc.dprint("re-setting page width to ",pagewidth);)
	}
	T_(trc.dprint("returning ",rval);)
	return rval;
}

void
dos2Unix(string& str) {
// Convert a string from dos-land to unix-land by removing
// all carriage-returns (\r). Dos (and Windows) end lines
// with the sequence "\r\n" whereas Unix ends lines with
// only "\n".
// The Mac ends lines with just "\r" so they must be replaced
// with "\n" if we ever start getting files from Macs (I don't
// know what convention OSX uses - maybe this only applies to
// OS9 and earlier)
	string::size_type idx = str.find('\r');
	while (idx != string::npos) {
		str.erase(idx, 1);
		idx = str.find('\r');
	}
}

// class Form implementation
Bound_form
Form::
operator() (double d) const { return Bound_form(*this,d); }

std::ostream& operator<<(std::ostream&  os, Bound_form const& bf) {
	std::ostringstream s;
	s.precision(bf.f.prc);
	s.setf(bf.f.fmt, std::ios_base::floatfield);
	s.width(bf.f.wdt);
	s.fill(bf.f.fillchar);
	s << bf.val;
	return os << s.str();
}

string
get_uname(const string& what) {
// get a piece of information about the kernel from uname(2)
// Valid options for "what":
//   sysname     Operating system name (e.g., "Linux")
//   nodename    Name within "some implementation-defined network"
//   release     Operating system release (e.g., "2.6.28") 
//   version     Operating system version
//   machine     Hardware identifier
//   domainname  NIS or YP domain name
	string rval;
	struct utsname buf;
	int status = uname(&buf);
	if (status == -1)
		throw std::system_error(errno,std::generic_category(),"uname error");
	
	// default (what empty):
	if (what.empty()) {
		rval = vastr(buf.sysname," ",buf.version," ",buf.machine);
		return rval;
	}

	if (what == "sysname")
		rval = buf.sysname;
	else if (what == "nodename")
			rval = buf.nodename;
	else if (what == "release")
		rval = buf.release;
	else if (what == "version")
		rval = buf.version;
	else if (what == "machine")
		rval = buf.machine;
	else if (what == "domainname")
		rval = buf.domainname;
	else
		throw runtime_error(vastr("unrecognized option to get_uname: ",what));
	return rval;
}

//---------------------------------------------------------------------
// Returns a string with a description of wait() status code
// or empty if normal exit
// Reference
// 1) this is known as APUI:
//   @book{stevens2008advanced,
//     title={Advanced programming in the UNIX environment},
//     author={Stevens, W Richard and Rago, Stephen A},
//     year={2008},
//     publisher={Addison-Wesley}
//   }
//
// Discussion:
//  the low-order 16 bits of status contain info regarding
//  the exit status of a child process.
//  Status differentiates between stopped and terminated
//  child processes.  If the child process has terminated,
//  status identifies the cause of termination and passes
//  useful info to the parent process. This is accomplished
//  in the following manner:
//  
//  If the child process has stopped, the high-order 8 bits
//  of status contain the number of the signal that caused
//  the process to stop and the low-order 8 bits are set
//  equal to 0177 (01111111).
//  
//  If the child process has terminated because of an exit
//  call, the low-order 8 bits of status are 0 and the
//  high-order 8 bits contain the low-order 8 bits of the
//  argument that the child process passed to exit(2).
//  
//  If the child process has terminated because of a
//  signal, the high-order 8 bits of status are 0 and the
//  low-order 8 bits contain the number of the signal that
//  caused the termination.  In addition, if the low-order
//  seventh bit (that is, bit 0200) is set, a core image
//  will have been produced; see man signal(2).
//  
//---------------------------------------------------------------------
string
exitstatus(int status) {
	string rval;

	if (status == 0)
		return rval;
	// I'm using the convention that an integer error code returned
	// by WEXITSTATUS is the number of errors encountered
	if ((WIFEXITED(status) != 0) && (WEXITSTATUS(status) != 0)) {
		int nerror = WEXITSTATUS(status);
		if (nerror == 1)
			rval = "one error";
		else
			rval = vastr(WEXITSTATUS(status)," errors");
	}
	
	if (WIFSIGNALED(status)) {
		int sig = WTERMSIG(status);
		rval = vastr(Signal::desc(sig), " (signal ", sig, ")");
	}
   return rval;
}

Ascii::
Ascii(const string& p) : ifs(p) {
// Ascii constructor: open path "p"
	lineno = 0;
	path = p;

	if (!ifs.good()) {
		if (p.empty()) {
			throw runtime_error("empty string passed to Ascii");
		}
		throw runtime_error(vastr("cannot access ", path));
	}
	// check that this is a regular file
	//!! int filedes = fileno(ifs.rdbuf()->stdiofile());
	struct stat fs;
	//!! fstat(filedes, &fs);
	stat(path.c_str(), &fs);
	if (!S_ISREG(fs.st_mode)) {
		throw runtime_error(vastr(path," is not a regular file"));
	}
}

#ifdef NEVER
Ascii&
Ascii::
operator=(Ascii& rhs) {
// assignment op
	return move(Ascii(rhs.path));
}
#endif // NEVER

bool
Ascii::
good() { return ifs.good(); }

bool
Ascii::
operator()(string& line, const string& comment, bool blanks) {
// read a line from "ifs", strip anything following "comment",
// skip empty lines, on eof, return false and empty line
	while(ifs) {
		std::getline(ifs, line);
		if (ifs.fail()) break;
		dos2Unix(line);		// remove \r from dos-land
		lineno++;
		// strip comments?
		if (!comment.empty()) {
			string::size_type idx = line.find(comment);
			if (idx == 0)
				continue;
			if (idx != string::npos)
				line = line.substr(0, idx);
		}
		// skip blank lines?
		if (!blanks) {
			string::size_type idx = line.find_first_not_of(" \t");
			if (idx != string::npos)
				return true;
		} else {
			return true;
		}
	}
	return false; // EOF
}


#ifdef MAIN

#include <iostream>

static void
usage(char const* prog) {
	cerr << "Usage: String option\n";
	cerr << "  -a string 'openclose':   test delimitedString\n";
	cerr << "  -b string:               test stringDirname, stringBasename\n";
	cerr << "  -c string string:        test compare\n";
	cerr << "  -d string:               test dos2Unix\n";
	cerr << "  -e string:               test expenv\n";
	cerr << "  -f number:               test str2complex\n";
	cerr << "  -F number:               test class Form\n";
	cerr << "  -g string:               test str2double\n";
	cerr << "  -i string:               test str2int\n";
	cerr << "  -j string width:         test stringJustify\n";
	cerr << "  -k string:               test get_uname\n";
	cerr << "  -m string:               test find_delim\n";
	cerr << "  -r string [pagewidth][indent]:   test roff\n";
	cerr << "  -s string:               test stripquotes\n";
	cerr << "  -t string delim [quotes]: test string2tok\n";
	cerr << "  -u string:               test flaps::info\n";
	cerr << "  -v:                      test Ascii\n";
	cerr << "  -w string:               test stripwhitespace\n";
	cerr << "  -z string 'f' 'c':       test replace_char\n";
	exit(1);
}

int
test_Form(int argc, char** argv) {
	if (argc != 2) {
		cerr << "Usage: " << argv[0] << " number\n";
		exit(1);
	}
	double d = strtod(argv[1], 0);
	Form sci12(12);
	Form fix2(3);
	Form gen20(20);
	sci12.scientific();
	fix2.fixed().precision(2);
	cout << "default: " << d << ", sci(12): " << sci12(d) <<
		", fix(2): " << fix2(d) << ", general(8): " << gen20(d) << '\n';
	return 0;
}

Ascii
open_ascii(const string& path) {
	return Ascii(path);	// should use move constructor
}
void
test_ascii(const string& path, const string& comment) {
	Ascii ifs = open_ascii(path);
	string line;
	// read, discard comment & blank lines
	while(ifs(line, comment, false))
		cout << "line " << ifs.line_number() << ": " << line << endl;
}

int
main (int argc, char** argv) {
	char* dbg = getenv("DEBUG");
	if (dbg != nullptr) Trace::debug(std::stoi(dbg));
	T_(Trace trc(1,argv[0]);)
	if (argc < 2 || argv[1][0] != '-') {
		usage(argv[0]);
	}
	char opt = argv[1][1];
	string argv2;
	if (argc > 2)
		argv2 = string(argv[2]);
	string argv3;
	if (argc > 3)
		argv3 = string(argv[3]);
	string argv4;
	if (argc > 4)
		argv4 = string(argv[4]);
	ostringstream os;

	try {
		if (opt == 'a') {
			char open{argv3[0]};
			char close{argv3[1]};
			string::size_type start{0};
			string::size_type end{0};
			cout << "testing delimitedString with open " << open <<
				", close " << close << endl;
			cout << delimitedString(argv2, open, close, start, end) <<
				", end " << end << endl;
		} else if (opt == 'b') {
			cout << "stringDirname(" << argv2 << ") = " << stringDirname(argv2) << endl;
			cout << "stringBasename(" << argv2 << ") = " << stringBasename(argv2) << endl;
		} else if (opt == 'c') {
			if (argc < 4)
				usage(argv[0]);
			string a(argv[2]);
			string b(argv[3]);
			bool result = compare(a,b);
			if (result) {
				cout << a << " equals " << b << endl;
			} else {
				cout << a << " does not equal " << b << endl;
			}
		} else if (opt == 'd') {
			cout << "diff between " << argv2 << " and " << argv3
				<< ": " << string_diff(argv2, argv3) << endl;
		} else if (opt == 'e') {
			cout << argv2 << " expands to " << expenv(argv2) << endl;
		} else if (opt == 'f') {
			size_t end;
			complex<double> xc;
			double xr;
			if (str2complex(argv2, xc, 0, end)) {
				cout << "got complex: " << xc << endl;
			} else if (str2double(argv2, xr, 0, end)) {
				cout << "got double: " << xr << endl;
			}
		} else if (opt == 'F') {
			test_Form(argc, argv);
		} else if (opt == 'g') {
			double x;
			size_t start = argv2.find_first_of("+-.0123456789");
			size_t end;
			if (str2double(argv2, x, start, end)) {
				cout << "str2double(" << argv2 << ") = " << x;
				if (end != string::npos)
					cout << ", trailing = " << argv2.substr(end) << endl;
				else
					cout << endl;
			} else {
				cout << argv2 << "is not a double" << endl;
			}
		} else if (opt == 'i') {
			size_t start = argv2.find_first_of("+-0123456789");
			size_t end;
			int k;
			if (str2int(argv2, k, start, end)) {
				cout << argv2 << " is an int\n";
				cout << "str2int(" << argv2 << ") = " << k;
				if (end != string::npos)
					cout << ", trailing = " << argv2.substr(end) << endl;
				else
					cout << endl;
			} else {
				cout << argv2 << " is not an int\n";
			}
		} else if (opt == 'j') {
			int width = 80;
			if (argc > 3)
				width = atoi(argv[3]);
			cout << "left-justified<" << stringJustify(width, argv2, "left") << ">\n";
			cout << "right-justified<" << stringJustify(width, argv2, "right") << ">\n";
			cout << "centered<" << stringJustify(width, argv2, "center") << ">\n";
		} else if (opt == 'k') {
			cout << "sysname: " << get_uname("sysname") << endl;
			cout << "nodename: " << get_uname("nodename") << endl;
			cout << "release: " << get_uname("release") << endl;
			cout << "version: " << get_uname("version") << endl;
			cout << "machine: " << get_uname("machine") << endl;
			cout << "domainname: " << get_uname("domainname") << endl;
		} else if (opt == 'm') {
			string::size_type close = argv2.find_last_not_of(" \t");
			if (close == string::npos) {
				cerr << "string is empty\n";
				exit(1);
			}
			string closers("}])>");
			if (closers.find(argv2[close]) == string::npos) {
				cerr << "last char is not a closing delimiter\n";
				exit(1);
			}
			string::size_type open = find_delim(argv2, close);
			if (open == string::npos) {
				cerr << "no opening " << argv2[close] << endl;
				exit(1);
			}
			cout << "opening " << argv2[close] << " is character " << open << endl;
			open = argv2.find_first_of(argv2[open]);
			close = find_delim(argv2, open);
			cout << "closing " << argv2[close] << " is character " << close << endl;
		} else if (opt == 'r') {
			unsigned int linelen = 80;
			unsigned int indent = 0;
			if (argc >= 4) {
				linelen = atoi(argv[3]);
				if (argc >= 5)
					indent = atoi(argv[4]);
			} else {
				cout << "Default line length (" << linelen << ") and indent ("
					<< indent << ")\n";
			}
			string out = roff(argv2, indent, linelen);
			cout << out << endl;
		} else if (opt == 's') {
			os << '\"' << argv2 << '\"';
			string arg(os.str());
			stripquotes(arg);
			cout << "stripquotes(" << arg << ") = <" << arg << ">\n";
		} else if (opt == 't') {
			vector<string> quotes;
			if (argc > 4) {
				for (int j=4; j<argc; j++)
					quotes.push_back(argv[j]);
			}
			vector<string> tok = string2tok(argv2, argv3, quotes);
			for (size_t i=0; i<tok.size(); i++)
				cout << i << " \"" << tok[i] << "\"\n";
		} else if (opt == 'u') {
			flaps::info("Test of flaps::info:\n", argv2);
		} else if (opt == 'v') {
			test_ascii(argv2,argv3);
		} else if (opt == 'w') {
			string newstr = stripwhitespace(argv2);
			cerr << "<" << argv2 << "> -> <" << newstr << ">\n";
		} else if (opt == 'z') {
			char f = argv3[0];
			char t = argv4[0];
			string str{argv2};
			string newstr = replace_char(str, f, t);
			cerr << "replaced " << f << " with " << t <<
				" in " << str << ": " << newstr << endl;
		}
	} catch (runtime_error& s) {
		cerr << "caught exception: " << s << endl;
	}
	return 0;
}
#endif // MAIN
