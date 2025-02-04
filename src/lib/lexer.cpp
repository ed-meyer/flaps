//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

// Preferences handling classes and functions.
// 1) The starting point is a string of preferences passed to Prefs::lexer
//    or if an empty string is passed lexer will take it from cin.
// 2) Prefs::lexer calls Tok::prelim to split the preferences into Pref*
// 3) Prefs::lexer creates an Iterator then calls Iterator::next() in a loop
//    to create each Pref*
// 4) the vector<Pref*> returned by lexer is then passed to
//    Prefs::lexer(vector<Pref*>, vector<Parser>) which loops through each
//    Pref looking for a match with the regex of each Parser. Prefs that do
//    not match any Parser are returned so the caller can try again with
//    different Parsers; this allows the caller to treat certain options first.
//
// The main lexical analysis goes on in:
//   Iterator::next ()
//         splits a preferences string into individual Pref
//         and returns the next in the list
//   Pref::Pref(lhs,rhs)
//         "Pref lexing constructor" - creates a Pref from the right- and left-hand sides


#include <cctype>
#include <fnmatch.h>
#include <string>
#include <vector>

#include "lexer.h"
#include "message.h"
#include "text.h"
#include "trace.h"
#include "vastr.h"

using namespace std;

static char open_paren{'('};
static char close_paren{')'};
static char single_quote{'\''};
static char double_quote{'\"'};
static char open_brace{'{'};
static char close_brace{'}'};
static char comma{','};
static char newline{'\n'};
static char blank{' '};
static char equals{'='};
static char pound_sign{'#'};

/*------------------------------------------------------------------
 * local (static) functions
 *------------------------------------------------------------------*/
static string get_thestring();
static vector<Tok*> prelim (string const& str);
static string rm_comments (string const& list);
static std::string::size_type find_close (std::string const& s,
		std::string::size_type from, char open, char close);
// expand a lhs implicit list
static vector<Tok*> expand_lhs(Tok* tok);
static vector<double> get_lhs_doubles(string const& opts);
static vector<string> tokenize(const string& str);


static
vector<int>
string2IntVec(vector<string> const& svec) {
// If all items are int convert and return
// XXX should go in text.cpp?
	vector<int> rval;
	int k;

	for (const auto& s : svec) {
		if (str2int(s, k)) {
			rval.push_back(k);
		} else {
			rval.clear();
			break;
		}
	}
	return rval;
}

static
vector<double>
string2RealVec(vector<string> const& svec) {
// If all items are double, convert and return,
// otherwise return empty
	vector<double> rval;
	double x;

	for (auto& s : svec) {
		if (str2double(s, x)) {
			rval.push_back(x);
		} else {
			rval.clear();
			break;
		}
	}
	return rval;
}

Tok::Iterator::
Iterator(string const& opt, bool quotes) {
// Iterator constructor
	prefstr = rm_comments(opt);
	current = 0;
}

Tok*
Tok::Iterator::
next () {
// Returns (a pointer to) the next preference in the preference string; this
// is the heart of the lexical analyzer lexer()
	T_(Trace trc(3,"Tok::Iterator::next");)
	Tok* rval{nullptr};

	T_(trc.dprint("current ",current);)

	// watch out for empty option list
	if (this->prefstr.empty()) {
		T_(trc.dprint("returning nullptr: empty option list");)
		return rval;
	}

	// skip whitespace
	while(1) {
		current = this->prefstr.find_first_not_of(" \t\n", current);
		if (current == string::npos) {
			T_(trc.dprint("returning nullptr: end of preferences");)
			return nullptr;
		}
		if (prefstr[current] != comma)
			break;
		else {
			current++;
			if (current == prefstr.size()) {
				T_(trc.dprint("returning nullptr: end of preferences");)
				return nullptr;
			}
		}
	}
	string::size_type start = current;
	string::size_type end = current;
	string::size_type lhs_term{0};

	// get the lhs; watch for quotes, parentheses, curly braces,
	// or square brackets: skip to the closing element if found
	// Also watch out for *=, /=, +=, -= (ops)
	char ch;
	string ops("*/+-");
	char op{0};
	while(current < prefstr.size()) {
		ch = prefstr[current];
		if (ch == single_quote || ch == double_quote) {
			current = find_close(prefstr, current+1, ch, ch);
		}
		if (ch == '[')
			current = find_close(prefstr, current+1, '[', ']');
		if (ch == '(')
			current = find_close(prefstr, current+1, '(', ')');
		if (ops.find(ch) != string::npos && prefstr[current+1] == '=') {
			op = ch;
			T_(trc.dprint("got operator <",op,">");)
			lhs_term = current;
			end = current+1;
			break;
		}
		// these guys end the lhs:
		if (ch == comma || ch == equals) {
			lhs_term = current;
			end = current;
			break;
		}
		// a newline *might* end the lhs, or preference-options
		// (opening curly brace) might be on the next line
		if (ch == newline) {
			string::size_type idx = prefstr.find_first_not_of(" \t\n", current);
			if (idx != string::npos && prefstr[idx] == open_brace) {
				current = idx;
				ch = prefstr[idx];
			} else {
				lhs_term = current;
				end = current;
				break;
			}
		}
		if (ch == open_brace)
			current = find_close(prefstr, current+1, open_brace, close_brace);
		current++;
	}
	// watch out for last lhs w/no rhs
	if (current >= prefstr.size()) {
		lhs_term = prefstr.size();
		end = prefstr.size();
	}
	T_(trc.dprint("finished w/lhs: end ",end,", lhs_term ",lhs_term,", current ",current);)
	string lhs = prefstr.substr(start, lhs_term-start);

	// get the rhs if the lhs ended with an =: might be a
	//   - list, e.g. (1 : 2 : 30)
	//   - operator, e.g. *= "eqn", then the eqn must be in quotes, so take
	//       up to the closing quote
	string rhs;

	// skip the lhs terminator
	current = end+1;
	// check for end of preferences
	if (prefstr[end] == equals && current < prefstr.size() && 
		(current = prefstr.find_first_not_of(" \t", current)) != string::npos) {
		start = current;
		while(current < prefstr.size()) {
			ch = prefstr[current];
			// a comma ends the rhs...
			if (ch == comma)
				break;
			// ... and a newline *might* end the rhs -
			// check for opening curly brace on the next line
			// XXX newline *always* ends a preference?
			if (ch == newline) {
				string::size_type idx = prefstr.find_first_not_of(" \t\n", current);
				if (idx != string::npos && prefstr[idx] == open_brace) {
					current = idx;
					ch = prefstr[idx];
				} else {
					break;
				}
			}
			if (ch == single_quote || ch == double_quote) {
				current = find_close(prefstr, current+1, ch, ch);
			}
			if (ch == '[')
				current = find_close(prefstr, current+1, '[', ']');
			if (ch == open_brace) {
				current =
					find_close(prefstr, current+1, open_brace, close_brace);
			}
			if (ch == open_paren) {
				current =
					find_close(prefstr, current+1, open_paren, close_paren);
			}
			current++;
		}
		// strip trailing whitespace
		rhs = stripwhitespace(prefstr.substr(start, current-start));

		// increment current to be ready for next call
		if (current+1 < prefstr.size())
			current++;
	}

	T_(trc.dprint("lhs<",lhs,"> op<",op,"> rhs<",rhs,">");)

	// Expand environment variables
	lhs = expenv(lhs);
	rhs = expenv(rhs);

	// Create the return Tok
	rval = new Tok(lhs, rhs, false, op);

	T_(trc.dprint("returning ",*rval);)
	return rval;
}

//------------------------------------------------------------------
// Parser implementation
//------------------------------------------------------------------

Parser::
Parser (string const& pat, std::function<bool(const Tok&)> f) :
	pattern(pat), rx(pat) {
// Main Parser constructor, the rx gets compiled above, when
// initializing rx with pat
	handler = f;
}

std::ostream&
operator<<(std::ostream& os, const Parser& t) {
	os << t.pattern;
	return os;
}

bool
flaps::
wildcard(const std::string& pattern, const std::string& str) {
// test a string against a wildcard pattern; see the manual for details
// on valid wildcard syntax
	T_(Trace trc(2,"wildcard");)
	bool rval{false};


	int res = fnmatch(pattern.c_str(), str.c_str(), FNM_EXTMATCH);
	T_(trc.dprint("testing ",str," against ",pattern,": ",res==0?"matches":"no match");)
	if (res == 0)
		rval = true;
	else if (res != FNM_NOMATCH)
		throw runtime_error(vastr("invalid wildcard pattern \"",pattern,"\""));
	return rval;
}

//------------------------------------------------------------------
// Tok implementation
//------------------------------------------------------------------

Tok::
Tok(string const& l, string const& r, bool leaveQuotes, char oper) {
/*------------------------------------------------------------------
 * Tok lexing constructor
 * Do some more lexical analysis and construct a Tok. May be all in the lhs, as in
 *     Tok("a=b","")
 * or split into rhs, lhs, as in
 *     Tok("a","b")
 * Multiple rhs are allowed, as in
 *     Tok("a","1,3 : 3 : 12,101)
 * but multiple lhs, e.g.
 *     Tok("(abcd : abcf)", "")
 * must be handled upstream since we can only return one Tok.
 * The right and left-hand sides may have options enclosed in braces,
 * as in
 *     Tok("a{x=y, z=1.0} = b{c=d, e=2}","")
 *------------------------------------------------------------------*/
	T_(Trace trc(3,"Tok lexing constructor");)

	T_(trc.dprint("lhs<",l,"> rhs<",r,"> leave quotes? ",leaveQuotes," op? ",oper);)

	op = oper;
	lhs = l;
	string rhs = r;

	// look for options following the lhs, e.g. modes{row=1}
	lopt = sub_opt(lhs);

	// ... and the rhs
	if (rhs.empty()) {
		dtype = Dtype::Invalid;
	} else {
		// look for multiple rhs, as in mode=(1,3:3:12)
#ifdef NEVER // embeddedlist fails with [1,1] = bilinear(0.01, 0.5, 1)
		svec = flaps::embeddedlist(rhs);
#else // NEVER // try embeddedlist
		svec = flaps::expand_list(rhs);
#endif // NEVER // try embeddedlist
		// check each rhs for trailing options: sub_opt removes the options
		// from s and returns them minus the curly-braces or an empty string
		// so that roptvec always is the same size as svec
		for (auto& s : svec)
			roptvec.push_back(sub_opt(s));
		// for convenience the first rhs opt is ropt, svec is srhs
		ropt = roptvec[0];
		srhs = svec[0];

		// classify the rhs according to the single rhs...
		// convert the string rhs to int, double, or complex
		dtype = Dtype::Char;
		ivec = string2IntVec(svec);
		rvec = string2RealVec(svec);
		if (!ivec.empty())
			dtype = Dtype::Int;
		if (!rvec.empty())
			dtype = Dtype::Real;
	}

	// strip leading, trailing whitespace from lhs, srhs, svec
	lhs = stripwhitespace(lhs);

	// strip leading, trailing quotes from lhs, srhs, svec
	if (!leaveQuotes) {
		stripquotes(lhs);
		stripquotes(srhs);
	}
	for (auto& si : svec) {
		si = stripwhitespace(si);
		stripquotes(si);
	}
	T_(trc.dprint("returning ",*this);)
}

/*------------------------------------------------------------------
 * Three versions of lexer():
 * 1) with just the option string lexer returns a vector of Tok*
 * 2) the string version takes a string containing all the options
 *    separated by commas or newlines, and a vector of Parsers
 * 3) the vector<Tok*> version takes a vector of Tok* and a vector of Parsers.
 * Both return a vector of Tok* which are the preferences that
 * were not identified by any of the handlers - this allows you to
 * handle only some options, then call lexer again with a different set of Parsers.
 *------------------------------------------------------------------*/

vector<Tok*>
flaps::
lexer (string const& setstr) {
// No handlers - just return a vector of Tok*
	vector<Parser> hl;
	return lexer(setstr, hl);
}

vector<Tok*>
flaps::
lexer (string const& setstr, const vector<Parser>& handlers) {
// Create a vector<Tok*> from "setstr", then call lexer(vector<Tok*>)
// "setstr" may be empty in which case it is read from cin
	T_(Trace trc(2,"flaps::lexer(string,vector<Parser>&)");)
	vector<Tok*> rval;

	string s = setstr;
	// If setstr is empty try to read it from cin
	if (s.empty())
		s = get_thestring();
	if (s.empty()) {
		T_(trc.dprint("no input string and nothing from cin");)
		return rval;
	}

	// lexer splits the preferences into Tok*
	vector<Tok*> setgs = prelim(s);

	// no handlers? just return setgs
	if (handlers.empty()) {
		T_(trc.dprint("returning ",setgs.size()," Tok: no handlers");)
		return setgs;
	}
	// lexer(vector<Tok*>) calls the appropriate handler for each Tok
	rval = flaps::lexer(setgs, handlers);
	T_(trc.dprint("returning ",rval.size()," unrecognized Tok");)
	return rval;
}

vector<Tok*>
flaps::
lexer (string const& setstr, vector<Parser>&& handlers) {
	T_(Trace trc(2,"flaps::lexer(string,vector<Parser>&&)");)
	vector<Parser> hl = handlers;
	return flaps::lexer(setstr, hl);
}

vector<Tok*>
flaps::
lexer (vector<Tok*> const& prefs, const vector<Parser>& handlers) {
//------------------------------------------------------------------
// For each Tok* in the input list check each of my
// Parsers. There are two possibilities if a Parser::pattern
// matches a Tok::lhs:
// 1) If the Parser::handler is a null pointer the Tok is
//    ignored
// 2) The Parser::handler function is called to parse the Tok;
//    if the function returns true we are finished with that Tok;
//    if the function returns false the search continues for a
//    Parser that matches. If none of the Parsers match the Tok*
//    is inserted into a vector which will be returned.
//------------------------------------------------------------------
	T_(Trace trc(2,"flaps::lexer(vector<Tok*>,vector<Parser>)");)
	vector<Tok*> rval;

	for (auto pref : prefs) {
		T_(trc.dprint("working on <",*pref,"> ");)
		T_(trc.dprint(pref->svec.size()," rhs, ",pref->roptvec.size()," rhs opts");)
		bool matched = false;
		// check the lhs of this Tok with the regex of each Parser until one matches
		for (auto& ap : handlers) {
			T_(trc.dprint("testing ",ap);)
			if (regex_match(pref->lhs, ap.rx)) {
				T_(trc.dprint("identified <",ap.pattern,">");)
				if (ap.handler == nullptr) {
					matched = true;
					T_(trc.dprint("ignoring and returning true: null handler");)
				} else {
					matched = ap.handler(*pref);
					T_(trc.dprint("handler returned ",matched?"true":"false");)
				}
			}
			if (matched)
				break;
		}
		// if none match put this Tok on the return vector
		if (!matched)
			rval.push_back(pref);
	}
	T_(trc.dprint("returning ",rval.size()," unrecognized Tok");)
	return rval;
}

static vector<Tok*>
prelim (string const& options) {
// Given a string of preferences separated by commas or newlines,
// break it into Tok pointers, return in a vector<Tok*>
	T_(Trace trc(3,"prelim");)
	vector<Tok*> rval;

	T_(trc.dprint("options string: <",options,'>');)

	// create an Iterator: breaks up "options" into Tok* which are
	// returned by Iterator::next()
	Tok::Iterator iterator(options);

	// if the lhs of an option has matlab-style lists, expand them and create
	// multiple options; Give each new lhs the lopt of the original, interpolating
	// if necessary
	Tok* op;
	while ((op = iterator.next())) {
		vector<Tok*> toby = expand_lhs(op);
		if (toby.size() <= 1) {
			rval.push_back(op);
		} else {
			for (auto ti : toby)
				rval.push_back(ti);
		}
	}

	// Special case: no rhs and lhs is an "rvalue", a number, e.g.
	//  lexer("x{0.1,0.2,0.3}") followed by lexer("0.1,0.2,0.3")
	vector<double> rv;
	vector<int> iv;
	bool special{true};
	for (auto pi : rval) {
		if (!pi->svec.empty() || !pi->lopt.empty() || !pi->ropt.empty()) {
			special = false;
			break;
		}
		double x;
		if (str2double(pi->lhs, x)) {
			rv.push_back(x);
			int i;
			if (str2int(pi->lhs, i))
				iv.push_back(i);
		} else {
			special = false;
			break;
		}
	}
	// ... delete all Toks, create a new one with just the
	// numbers in rvec and ivec
	if (special) {
		for (auto pi : rval)
			delete pi;
		rval.clear();
		Tok* solo = new Tok(options,"");
		solo->rvec = rv;
		if (iv.size() == rv.size())
			solo->ivec = iv;
		rval.push_back(solo);
		T_(trc.dprint("special case: put all lhs values in rvec: ",rv);)
	}

	T_(trc.dprint("returning ",rval.size()," Tok");)
	return rval;
}

bool
Tok::
compare(string const& pat) const {
// Compare my lhs with a ECMA extended regular-expression (pat) ignoring case
	return regex_match(lhs, regex(pat, regex_constants::icase));
}

bool
Tok::
rhscompare(string const& pat) const {
// Compare my rhs with "pat":
// - as 2 strings
// - as a regex and a string ignoring case
// - as a string and a regex ignoring case
// - as 2 doubles
	if (this->srhs == pat)
		return true;
	// Compare my rhs with a ECMA regular-expression (pat) ignoring case
	if (regex_match(this->srhs, regex(pat, regex_constants::icase)))
		return true;
	if (regex_match(pat, regex(this->srhs, regex_constants::icase)))
		return true;
	// check for both ints or doubles
	double a{0.0};
	double b{0.0};
	if (str2double(this->srhs, a) && str2double(pat, b) &&
			is_equal(a, b, 6))
		return true;

	return false;
}

std::ostream&
operator<<(std::ostream& os, const Tok& t) {

	// how many rhs?
	size_t n = t.svec.size();

	os << t.lhs;
	if (!t.lopt.empty())
		os << '{' << t.lopt << '}';

	// if op is not zero the '=' becomes op=
	string opeq{"="};
	if (t.op)
		opeq = string(1,t.op) + opeq;

	// the rhs may be a list of int, double, or strings
	if (n > 1) {
		os << " " << opeq << " (";
		for (size_t i=0; i<n; i++) {
			os << t.svec[i];
			if (!t.roptvec[i].empty())
				os << '{' << t.roptvec[i] << '}';
			if (i < n-1)
				os << ", ";
		}
		os << ")";
	} else if (!t.svec.empty()) {
		// quote the rhs if it contains space(s)
		string::size_type idx = t.srhs.find(blank);
		if (idx == string::npos) {
			os << ' ' << opeq << ' ' << t.srhs;
		} else {
			os << ' ' << opeq << " \"" << t.srhs << '\"';
		}
		if (!t.ropt.empty())
			os << '{' << t.ropt << '}';
	}
	return os;
}

static vector<string>
expand(const string& first, int incr, const string& last);

vector<string>
flaps::
implist(const string& str) {
// given an implicit list like 1:2:5, returns the expanded
// list (1,3,5)
	T_(Trace trc(2,"implist");)
	vector<string> rval;
	vector<string> toks = string2tok(str, ":");
	if (toks.size() < 2)
		return rval;
	int incr{1};
	string first = toks[0];
	string last;
	if (toks.size() == 2)
		last = toks[1];
	else if (toks.size() == 3) {
		if (!str2int(toks[1], incr))
			throw runtime_error(vastr("implicit-list increment (",toks[1],") must be an integer"));
		last = toks[2];
	}
	return expand(first, incr, last);
}

enum class Cc { letter, digit, blank, etc, end };
static Cc
charclass(char c) {
	if (isalpha(c))
		return Cc::letter;
	else if (isdigit(c))
		return Cc::digit;
	else if (isblank(c))
		return Cc::blank;
	return Cc::etc;
}

static vector<string>
tokenize(const string& str) {
// tokenize a string based on character classes, i.e. a group
// of letters forms a token, likewise a group of digits, blanks,
// or 'other' characters such as .,:#/[]{}, etc
	T_(Trace trc(1,"tokenize ",str);)
	vector<string> rval;
	size_t start{0};
	size_t n = str.size();
	Cc cc = charclass(str[0]);
	for (size_t i=1; i<=n; i++) {
		Cc ic = Cc::end;
		if (i < n)
			ic = charclass(str[i]);
		if (ic != cc) {
			rval.push_back(str.substr(start, i-start));
			cc = ic;
			start = i;
		}
	}
	T_(trc.dprint("returning ",rval);)
	return rval;
}

vector<string>
expand(const string& first, int incr, const string& last) {
// Given pieces of an implicit list, expand into the equivalent strings
	T_(Trace trc(2,"expand");)
	vector<string> rval;

	// split first & last into leader/var/trailer
	string leader;
	string trailer;
	string va;
	string vb;
	vector<string> firsttok = tokenize(first);
	vector<string> lasttok = tokenize(last);
	size_t n{firsttok.size()};
	if (n != lasttok.size())
		return rval;
		// throw runtime_error(vastr("could not parse \"",first,"\" : \"",last,"\""));
	size_t i{0};
	for (; i<n; i++) {
		if (firsttok[i] != lasttok[i])
			break;
		leader += firsttok[i];
	}
	// variables: add tokens until they match again
	va = firsttok[i];
	vb = lasttok[i];
	for (i++; i<n; i++) {
		if (firsttok[i] == lasttok[i])
			break;
		va += firsttok[i];
		vb += lasttok[i];
	}
	for (; i<n; i++)
		trailer += firsttok[i];
	T_(trc.dprint("leader<",leader,"> va<",va,"> vb<",vb,"> trailer<",trailer,">");)

	// now we are ready to create the list
	int v1;
	int v2;
	if (str2int(va, v1) && str2int(vb,v2)) {
		for (int v = v1; v<=v2; v += incr)
			rval.push_back(leader + std::to_string(v) + trailer);
	} else {
		// character var: find the first char (k) that differs in va, vb
		// that is the char that will be the variable
		int m = std::min(va.size(), vb.size());
		int k{0};
		for (int i=0; i<m; i++) {
			if (va[i] != vb[i]) {
				k = i;
				break;
			}
		}
		//!! if (va.size() > 1 || vb.size() > 1)
			//!! throw runtime_error("character variables in implicit lists must only be 1 character");
		for (string v = va; v[k] <= vb[k]; v[k] += incr)
			rval.push_back(leader + v + trailer);
	}
	T_(trc.dprint("returning ",rval.size()," items: ",rval);)
	return rval;
}

static vector<double>
get_lhs_doubles(string const& prefs) {
// parse "prefs" and get the doubles which are the lhs of
// each pref
	T_(Trace trc(1,"get_lhs_doubles");)
	vector<double> rval;
	// parse "prefs" with an rvalue vector of a lambda handler
	vector<Tok*> ol = flaps::lexer(prefs,
			vector<Parser>{{".*",[&](const Tok& p){
				T_(trc.dprint("got pref ",p);)
				rval = p.rvec;
				return true;}}
			});
	T_(trc.dprint("returning ",rval.size()," doubles: ",rval);)
	return rval;
}

static vector<Tok*>
expand_lhs(Tok* pref) {
// expand a lhs like "1{0,1,0} : 4{0,4,0}" into new Toks, e.g.
//  1{0,1,0}, 2{0,2,0}, 3{0,3,0}, 4{0,4,0}
// i.e. interpolate the lhs and any sub-opts they have
// Requirements:
// 1) no rhs
// 2) has 1 or 2 colons (:)
// 3) first and last both have either 1 or 0 lopt
// 4) only allow ints for stride?
// 5) first and last have the form prefix+double+suffix
	T_(Trace trc(1,"expand_lhs");)
	vector<Tok*> rval;
	constexpr auto nop = string::npos;

	// 1) rhs?
	if (!pref->svec.empty()) {
		rval.push_back(pref);	// return the input
		return rval;
	}

	string lhs = pref->lhs;

	// 2) split the lhs into first : stride : last
	auto idx = lhs.find(':');
	if (idx == nop || (lhs.rfind('[', idx) != nop && lhs.find(']',idx) != nop)) {
		rval.push_back(pref);	// return the input
		return rval;
	}
	auto i = idx-1;
	for (i=idx; lhs[i-1]==' '; i--);
	string first = lhs.substr(0,i);
	while(lhs[++idx] == ' ');
	string therest = lhs.substr(idx);
	T_(trc.dprint("first: ",first,", the rest: ",therest);)
	idx = therest.find(':');
	string last;
	string first_lopt;
	string last_lopt;
	string stride{"1"};
	if (idx == string::npos)
		last = therest;
	else {
		for (i=idx; therest[i-1]==' '; i--);
		stride = therest.substr(0,i);
		last = therest.substr(idx+1);
		T_(trc.dprint("stride: ",stride, ", last: ",last);)
	}
	// strip the first lopt if there and if input pref has
	// an lopt from 'last'
	last_lopt = pref->lopt;
	//!! if (!last_lopt.empty()) {
		T_(trc.dprint("last lopt: ",last_lopt);)
		idx = first.find('{');
		if (idx != string::npos) {
			first_lopt = first.substr(idx);
			string::size_type end;
			first_lopt = delimitedString(first_lopt, '{', '}', 0, end);
			first = first.substr(0,idx);
			T_(trc.dprint("first: ",first,", first lopt: ",first_lopt);)
		}
	//!! }

	// 3) first and last both have lopt
	if (!(last_lopt.empty() == first_lopt.empty())) {
		T_(trc.dprint("no expansion: first,last opt");)
		rval.push_back(pref);	// return the input
		return rval;
	}
	// 4) only allow ints for stride?
	int istride;
	if (!str2int(stride, istride)) {
		rval.push_back(pref);	// return the input
		return rval;
	}

	// 5) split first & last into prefix/number/suffix
	regex re{regex_double()};
	smatch mch;
	string fpre, fnum, fsuf;
	string lpre, lnum, lsuf;
	if (regex_search(first, mch, re)) {
		fpre = mch.prefix().str();
		fnum = mch[1].str();
		fsuf = mch.suffix().str();
		T_(trc.dprint("first: <",fpre,"><",fnum,"><",fsuf,">");)
	}
	if (regex_search(last, mch, re)) {
		lpre = mch.prefix().str();
		lnum = mch[1].str();
		lsuf = mch.suffix().str();
		T_(trc.dprint("last: <",lpre,"><",lnum,"><",lsuf,">");)
	}

	// 5) convert the first and last numbers to ints
	int a, b;
	if (!str2int(fnum,a) || !str2int(lnum,b)) {
		T_(trc.dprint("returning input Tok: first or last do not have a number");)
		rval.push_back(pref);
		return rval;
	}
	int n = (b - a)/istride + 1;
	// prepare for interpolation of the lopts: get doubles
	vector<double> first_val;
	vector<double> last_val;
	vector<double> del;
	if (!first_lopt.empty() && !last_lopt.empty()) {
		first_val = get_lhs_doubles(first_lopt);
		last_val = get_lhs_doubles(last_lopt);
		if (first_val.size() != last_val.size())
			throw runtime_error("lhs options are not the same size");
		for (size_t i=0; i<first_val.size(); i++)
			del.push_back((last_val[i] - first_val[i])/(n-1));
		T_(trc.dprint("dels: ",del);)
	}

	// create new Toks interpolating the lhs number and it's subopt
	for (int i=0; i<n; i++) {
		int v = a + i*istride;
		Tok* pi = new Tok(*pref);	// copy
		pi->lhs = vastr(fpre,v,fsuf); 
		if (!pi->rvec.empty())
			throw runtime_error(vastr("option \"",first,"\" has a rhs"));
		string sep{""};
		ostringstream os;
		for (size_t j=0; j<first_val.size(); j++) {
			double x = first_val[j] + i*del[j];
			os << sep << x;
			sep = ",";
			//!! pi->rvec.push_back(x);
		}
		pi->lopt = os.str();
		rval.push_back(pi);
	}
	T_(trc.dprint("returning ",rval.size()," new Tok");)
	return rval;
}

vector<string>
flaps::
embeddedlist (string const& str) {
// Given a string like
//   abcd(1,3,5,7:2:11)efg
// or equivalently
//   abcd(1:2:11)efg
// returns the expanded list in a vector of strings
	T_(Trace trc(2,"embeddedlist ",str);)
	vector<string> rval;

	// watch out for quoted string
	if (str[0] == '\"' || str[0] == '\'') {
		T_(trc.dprint("quick return: quoted string");)
		rval.push_back(str);
		return rval;
	}

	// first split into comma-separated tokens
	vector<string> quotes;
	vector<string> toks = string2tok(str, ",", quotes);
	T_(trc.dprint("comma-separated tokens: ",toks);)

	// split each token containing an implicit list into: head(list)tail
	for (size_t i=0; i<toks.size(); i++) {
		string toki = toks[i];
		string head;
		string list;
		string tail;
		string::size_type op = toki.find_first_of('(');
		// no parentheses? might still be a non-embedded list
		if (op == string::npos) {
			op = toki.find_first_of(':');
			if (op == string::npos) {
				rval.push_back(toki);
			} else {
				vector<string> tokitoks = expand_list('('+toki+')');
				for (auto tokj : tokitoks)
					rval.push_back(tokj);
			}
			continue;
		}
		if (op > 0)
			head = toki.substr(0,op);
		// find the matching paren
		string::size_type cp = find_delim(toki, op);
		if (cp == string::npos)
			throw runtime_error(vastr("no closing parenthesis in ",toki));
		list = toki.substr(op, cp-op+1);
		if (cp < toki.size())
			tail = toki.substr(cp+1);
		T_(trc.dprint("head \"",head,"\", list \"",list,"\", tail \"",tail,"\"");)
		// expand the list part
		vector<string> items = expand_list(list);
		// create output strings from the items
		for (auto i : items)
			rval.push_back(head + i + tail);
	}
	T_(trc.dprint("returning ",rval.size()," items: ",rval);)
	return rval;
}

vector<string>
flaps::
expand_list (string const& str) {
// Given a string like
//   (1,3,5,7:2:11)
// or equivalently
//   (1:2:11)
// returns the expanded list in a vector of strings
	T_(Trace trc(1,"expand_list ",str);)
	vector<string> rval;

	// watch out for quoted string
	if (str[0] == '\"' || str[0] == '\'') {
		T_(trc.dprint("quick return: quoted string");)
		rval.push_back(str);
		return rval;
	}

	// the string must start and end with parentheses; remove them
	string list;
	if (str[0] == '(' && str[str.size()-1] == ')') {
		list = str.substr(1, str.size()-2);
	} else {
#ifdef NEVER // do not require parentheses
		list = str;
#else // NEVER // do not require parentheses
		T_(trc.dprint("quick return: no parentheses");)
		rval.push_back(str);
		return rval;
#endif // NEVER // do not require parentheses
	}

	// first split into comma-separated tokens; watch out for
	// sub-options enclosed in {} XXX any others?
	vector<string> quotes{{"{}"}};
	vector<string> toks = string2tok(list, ",", quotes);
	T_(trc.dprint("comma-separated tokens: ",toks);)

	// check each one for colons (1 or 2)
	for (auto& tok : toks) {
		vector<string> bits = flaps::implist(tok);	
		if (bits.empty()) {
			rval.push_back(tok);
		} else {
			for (auto& bit : bits)
				rval.push_back(bit);
		}
	}

	T_(trc.dprint("returning ",rval.size()," items: ",rval);)
	return rval;
}

static string
get_thestring() {
// Read preferences from cin or return previously-read preferences, so a
// program can call this function multiple times and get the same string.
	T_(Trace trc(1,"get_thestring");)
	static string rval;
	if (rval.empty()) {
#ifdef NEVER // doesn't work with piped data
		if (cin.rdbuf()->in_avail() == 0) {
			T_(trc.dprint("nothing to read from cin: returning empty");)
			return rval;
		}
#endif // NEVER // doesn't work with piped data
		cin >> noskipws;
		char ch;
		while (cin >> ch) {
			rval.push_back(ch);
			if (cin.eof())
				break;
		}
	}
	T_(trc.dprint("got \"",rval,"\"");)
	return rval;
}

/*------------------------------------------------------------------
 * 
 *------------------------------------------------------------------*/

static
string
rm_comments (string const& list) {
// Remove comments from a multi-line string:
// the remainder of a line following a "hash" (pound_sign) character.
	T_(Trace trc(3,"rm_comments");)
	string rval;

	T_(trc.dprint("list: ",list);)
	
	if (list.empty()) {
		T_(trc.dprint("returning empty string");)
		return rval;
	}

	string::size_type first = 0;
	while(true) {
		string::size_type idx = list.find(pound_sign,first);
		if (idx == string::npos) {
			rval += list.substr(first);
			break;
		}
		while(list[idx-1] == blank) idx--;
		rval += list.substr(first, idx-first);
		first = list.find(newline, idx);
		if (first == string::npos)
			break;
	}
	T_(trc.dprint("returning ",rval);)
	return rval;
}

string::size_type
find_close (string const& s, string::size_type from,
		char open, char close) {
/*
 * Skip to the "close" character which matches "open",
 * watching for enclosed open/close pairs if "open" and
 * "close" are different.
 * s is assumed to point 1 past the opening char XXX it has
 * to be either this or point to the open char
 * Note: this function is recursive.
 */
	T_(Trace trc(3,"find_close ",open,close);)
	ostringstream os;
	string::size_type rval;

	// If open and close are the same all we have to do is look for close...
	if (open == close) {
		T_(trc.dprint("searching for <",close, "> in <",s.substr(from),"> starting at ",from);)
		rval = s.find(close, from);
		if (rval == string::npos) {
			throw runtime_error(vastr("no closing \'", close, "\' in \"", s, '\"'));
		}
	} else {

	// ... otherwise we must watch for open and skip enclosed open/close stuff
		rval = from;
		while (rval < s.size()) {
			if (s[rval] == open) {
				rval = find_close(s, rval+1, open, close);
				rval++;
			}
			if (s[rval] == close)
				break;
			rval++;
		}
		if (rval == s.size()) {
			throw runtime_error(vastr("no closing \'",close,"\' in \"",s,'\"'));
		}
	}
	T_(trc.dprint("found close at ",rval);)
	return rval;
}


string
sub_opt(string& opt) {
// Strip a trailing {} delimited set of options from a string, and return
// the delimited string (sans braces), and the original string without the
// options in "opt". Inner braces are allowed; double or single quotes
// protect the curly braces.
	T_(Trace trc(3,"sub_opt"," <",opt,'>');)
	string rval;

	// skip trailing whitespace, return if last
	// character is not close brace
	string::size_type close = opt.find_last_not_of(" \t\n");
	if (close == string::npos || opt[close] != close_brace) {
		T_(trc.dprint("returning empty: last char is not closing brace");)
		return rval;
	}
	string::size_type openb = find_delim(opt, close);
	if (openb == string::npos) {
		T_(trc.dprint("returning empty: no opening brace");)
		return rval;
	}

	if (openb+1 < opt.size() && close > 0) {
		string::size_type start = opt.find_first_not_of(" \t", openb+1);
		string::size_type end = opt.find_last_not_of(" \t", close-1);
		if (end >= start)
			rval = opt.substr(start, end-start+1);

		if (openb > 0) {
			end = opt.find_last_not_of(" \t", openb-1);
			// the opt may consist only of options
			if (end == string::npos)
				opt = "";
			else
				opt = opt.substr(0,end+1);
		} else {
			opt = "";
		}
	}
	T_(trc.dprint("returning \"",rval,"\"");)
	return rval;
}


#ifdef MAIN

#include "settings.h"
#include "main.h"

bool abcd(const Tok& p) {
	cerr << "abcd: got pref " << p << endl;
	return true;
}

bool
mass(const Tok& p) {
	cerr << "mass: got pref " << p << endl;
	return true;
}

void
test_handlers() {
	T_(Trace trc(1,"test_handlers");)
	// this is how to parse a preferences string using Parsers,
	// one with a lambda, the other a function (mass)
	vector<Parser> hl{
		{"a(bcd)?", [&](const Tok& p) {cerr << "abcd got Tok " << p << endl; return true; }},
		{"m[1-4]", mass}
	};
	string str{"abcd=junk, m2=a, m3=b, f=g"};
	vector<Tok*> vs = flaps::lexer(str, hl);
	T_(trc.dprint("lexer returned ",vs.size()," unrecognized preferences:",vs);)
	// try again with the unrecognized Tok and new Parsers
	vector<Parser> hl2{
		{"f(req)?", abcd},
		{"m[1-4]", mass}
	};
	vs = flaps::lexer(vs, hl2);

	T_(trc.dprint("lexer returned ",vs.size()," unrecognized preferences:",vs);)
}

void
test_settings() {
// Test Settings
	try {
		cout << "---------- testing Settings\n";
		double cf;
		if (Settings::defaults.get("minstepsize", cf))
			cout << "minstepsize: " << cf << endl;
		else
			cout << "minstepsize is not defined\n";
		bool b;
		if (Settings::defaults.get("constrained",b))
			cout << "constrained: " << b << endl;
		else
			cout << "constrained is not defined\n";
		Settings::defaults.set("constrained",true);
		Settings::defaults.get("constrained",b);
		cout << "new constrained: " << b << endl;
	} catch (runtime_error& s) {
		cerr << "caught exception: " << s << endl;
	}

	// test Settings::global
	Settings::global("d=2, env='path=/tmp'");

	try {
		// test exception: try setting a bool policy to int
		// XXX no longer throws - just resets
		int i{1};
		Settings::defaults.set("constrained",i);
	} catch (runtime_error& s) {
		cerr << "caught exception: " << s << endl;
	}
}

string
pref_details(Tok const& t) {
	ostringstream os;
	os << t << ":\n";
	os << "   lhs: <" << t.lhs << ">\n";
	os << "   " << t.svec.size() << " rhs strings:\n";
	for (size_t j=0; j<t.svec.size(); j++)
		os << "      <" << t.svec[j] << "> ropt<"
			<< t.roptvec[j] << ">\n";
	os << "   " << t.rvec.size() << " reals\n";
	os << "   " << t.ivec.size() << " ints\n";
	os << "   lopt<" << t.lopt << ">\n";
	os << "   ropt<" << t.ropt << ">\n";
	return os.str();
}

int
main(int argc, char **argv) {
// Usage:
// settings [-s][-i] ["op1,op2,..."]
//   with no arguments, test Parsers
//   -i:  interactive typing preferences,
//        otherwise the first arg is taken as a preferences string, e.g.
//          settings "a=b, c=d"
//   -s   test Settings
	T_(Trace trc(1,argv[0]);)

	// test flaps::lexer and prelim
	try {
		string preferences;
		if (argc == 1) {
			test_handlers();
		} else if (argc == 2 && string(argv[1]) == "-i") {
			cout << " preferences test: type an preferences list...end with ctrl-d\n";
			preferences = get_thestring();
			T_(trc.dprint("got preferences string <",preferences,">");)
		} else if (argc == 2 && string(argv[1]) == "-s") {
			test_settings();
		} else {
			// parse a set of prefs in quotes in argv[1] using an rvalue
			// vector of lambda handlers
			preferences = argv[1];
			vector<Tok*> ol = flaps::lexer(preferences, {
				{".*",[](const Tok& p){
					cout << "got " << pref_details(p) << endl;
					return true;}
				}
			});
		}
	} catch (runtime_error& s) {
		cerr << "caught exception: " << s << endl;
	}
	return 0;
}
#endif // MAIN
