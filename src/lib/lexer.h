//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
#ifndef LEXER_H
#define LEXER_H
/*
 * Split a string into comma- or newline-separated tokens and optionally
 * pass each token to a set of functions provided by the caller to
 * parse the tokens.
 *
 * User preferences always come to user-callable functions as the
 * one-and-only argument; programs such as flut get preferences from
 * cin (standard input stream).
 * Parsing the preferences is done by writing a set of Parsers which
 * take an individual preference and create a setting or some other
 * task. For example
 *    bool input(const Tok& p) { Settings::set("input",p.srhs); return true }
 *    bool output(const Tok& p) { Settings::set("output",p.srhs); return true }
 *    vector<Parser> hl {
 *       {"i(n|nput)?",input}, {"o(ut|utput)?",output}
 *    };
 *    vector<Tok> unrec = Settings::parse(prefstr, hl);
 * this will set 2 preferences: the names of the input and output.
 * If there were preferences in the prefstr that were not handled by
 * any of the Parsers in hl, they will be returned by parse and may
 * be treated by a subsequent call to parse with those and a different
 * set of Parsers:
 *    bool others(const Tok& p) { ... }
 *    vector<Parser> hl2 { {".*", others} };
 *    vector<Tok> unrec2 = Settings::parse(unrec, hl2);
 *
 * "Settings" is a namespace comprising several functions for storing
 * and retrieving data items.
 * Settings control the execution of Flaps programs and certain functions.
 * Settings are created from user-supplied preferences and from within
 * functions and programs at run time.
 */

#include <any>
// #include <boost/any.hpp>
#include <exception>
#include <functional>
#include <iostream>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "fptype.h"
#include "message.h"
#include "Regex.h"

enum class Dtype { Invalid, Char, Int, Real, Cmplx };
// std::ostream& operator<<(std::ostream& s, Dtype const& t);

class Tok;

class Parser {
public:
	std::string pattern;
	std::function<bool(const Tok&)> handler;
	std::regex rx;
	// constructors: the use of std::function allows for function objects,
	// lambdas and regular functions
	Parser (std::string const& pattern, std::function<bool(const Tok&)>);
	// default copy ctor, assignment op, dtor
	Parser (Parser const& from) noexcept = default;
	Parser (Parser&&) noexcept = default;
	Parser& operator=(Parser const&) = default;
	Parser& operator=(Parser&&) = default;
	~Parser() = default;
};


class Tok {
	Dtype dtype;
public:
	class Iterator {
		std::string prefstr;
		std::string::size_type current;
	public:
		Iterator(std::string const& opt, bool quotes = false);
		Tok* next();
	};  // Iterator
	std::string lhs;
	// possibly multiple rhs:
	std::vector<int> ivec;
	std::vector<double> rvec;
	std::vector<std::string> svec;
	std::vector<std::complex<double> > cvec;
	// srhs is svec[0] if !svec.empty()
	std::string srhs;

	// options on the lhs or rhs:
	std::string lopt;		// as in GSTIF{mass=...}
	std::string ropt;		// as in gc=2{mode2...}
	// If there are multiple rhs, each rhs may have its own trailing
	// option, e.g. i=(gm1{rf=0}, gm2{rf=0.004},...); roptvec[i] will be
	// an empty string if no trailing option
	std::vector<std::string> roptvec;
	// op is */+- in lhs *= rhs, etc
	char op{0};

	// constructors
	Tok() = default;
	// Primary lexing constructor
	Tok(std::string const& lhs, std::string const& rhs,
			bool leaveQuotes=false, char oper=0);
	// default copy/move constructors, assignment operators ok
	Tok(const Tok&) noexcept = default;
	Tok(Tok&&) noexcept = default;
	Tok& operator=(const Tok&) = default;
	Tok& operator=(Tok&&) = default;

	bool compare(std::string const& pat) const;
	bool rhscompare(std::string const& pat) const;

	// rhs datatype: use these instead of accessing dtype member:
	// if there is no rhs member svec will be empty and these return false
	bool is_char() const { return dtype == Dtype::Char; }
	bool is_int() const { return dtype == Dtype::Int; }
	bool is_real() const { return dtype == Dtype::Real; }
	bool is_complex() const { return dtype == Dtype::Cmplx; }
}; // class Tok


	/*------------------------------------------------------------------
	 * 4 versions of lexer():
	 * 1) (string) with just the option string (all the options separated
	 *    by commas or newlines) lexer returns a vector of Tok*
	 * 2) (string,vector<Parser>&) takes the option string and an lvalue
	 *    vector of Parsers to do the parsing
	 * 3) (string,vector<Parser>&&) takes the option string and an rvalue
	 *    vector of Parsers to do the parsing
	 * 4) the vector<Tok*> version takes a vector of Tok* and a vector of Parsers
	 * All return a vector of Tok* which are the preferences that
	 * were not identified by any of the handlers - this allow you to parse
	 * with one set of Parsers follow by parsing with a different set
	 *------------------------------------------------------------------*/
namespace flaps {
	// test a string against a wildcard pattern; see the manual for details
	// // on valid wildcard syntax
	bool wildcard(const std::string& pattern, const std::string& str);

	std::vector<Tok*> lexer (std::string const& setstr);
	std::vector<Tok*> lexer (std::string const& setstr,
			const std::vector<Parser>& handlers);		// lvalue handlers
	std::vector<Tok*> lexer (std::string const& setstr,
			std::vector<Parser>&& handlers);		// rvalue handlers, e.g. brace-
																// enclosed initializer list
	std::vector<Tok*> lexer (std::vector<Tok*> const& setgs,
			const std::vector<Parser>& handlers);

	// tokenize a string with commas and implicit lists
	std::vector<std::string> expand_list (std::string const& str);
	std::vector<std::string> embeddedlist (std::string const& str);
	// tokenize a string with an implicit list: first : incr : last
	std::vector<std::string> implist(const std::string& str);
}	// namespace flaps

// a Policies is a list of items (policies) that can be added to
// Each policy is one of:
//    bool
//    int
//    double
//    std::string
//    std::vector<std::string>
//    std::vector<int>
//
// Public interface:
//    bool is_defined(name)
//    These return value in second arg:
//      void operator()(string, bool)
//      void operator()(string, int)
//           ...
//    These set the value of a policy:
//      bool set(name, bool)       set policy "name", return existing value
//         ...
//
// a set of default policies is kept in default_policies; this means
// to get a default policy you must say default_policies()(name,val)
// An exception is thrown if the policy "name" and "val" have different
// datatypes.
// To check if a default policy is defined:
//    if (default_policies().is_defined(name)) ...



// split, e.g. mass{cg=5} into "mass" and "cg=5"
std::string
sub_opt(std::string& pref);

std::ostream& operator<<(std::ostream& s, const Parser& t);
std::ostream& operator<<(std::ostream& s, const Tok& t);

#endif  // LEXER_H
