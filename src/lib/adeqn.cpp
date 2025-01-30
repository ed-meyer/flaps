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
#include <sstream>
#include <string>

#include "adeqn.h"
#include "atmos.h"
#include "Ad.h"
#include "df.h"
#include "matrix.h"
#include "pset.h"
#include "trace.h"

using namespace std;

enum class Kind : char {
	number, text, end,
	plus='+', minus='-',mul='*',div='/',hat='^',lp='(',rp=')',comma=',',
	lbrace='{', rbrace='}'
};

struct Token {
	Kind kind;
	Ad number_value;
	string text_value;
};

class Adeqn {
public:
	Adeqn(const string& s, pset& pl, bool nexc=true) :
		ip{s}, params{pl}, noexc{nexc} {}
	Ad expr(bool get);
private:
	Token next(bool nofcn=false);	// read & return next token
	Ad term(bool get);
	Ad prim(bool get);
	bool adfcn(const string& name, Ad& result);
	bool matrix_element(const string& name, Ad& result);
	bool dfval(Ad& result);
	istringstream ip;
	pset& params;
	Token ct{Kind::end};  // current token
	bool noexc;				// catch exceptions?
};

Ad
Adeqn::
expr(bool get) {
// This is the main entry point to the adeqn parser
// handle add & subtract
	Ad left = term(get);
	for(;;) {
		switch (ct.kind) {
			case Kind::plus:
				left += term(true);
				break;
			case Kind::minus:
				left -= term(true);
				break;
			default:
				return left;
		}
	}
}

Ad
Adeqn::
term(bool get) {
// handle multiply & divide
	Ad left = prim(get);
	for (;;) {
		switch (ct.kind) {
			case Kind::mul:
				left *= prim(true);
				break;
			case Kind::div:
			{
				Ad d = prim(true);
				if (d.value() != 0.0) {
					left /= d;
				} else if (!noexc) {
#ifdef NEVER // check for 0/0
					throw runtime_error("divide by zero");
#else // NEVER // check for 0/0
					if (left != 0.0)
						left = sqrt(sqrt(std::numeric_limits<double>::max()));
#endif // NEVER // check for 0/0
				}
				break;
			}
			case Kind::hat:
			{
				Ad exp = prim(true);
				left = pow(left,exp);
				break;
			}
			default:
				return left;
		}
	}
}

Ad
Adeqn::
prim(bool get) {
// handle primaries: numbers, unary minus, expressions in parentheses
	if (get) next();

	switch (ct.kind) {
		case Kind::number:
		{
			Ad v = ct.number_value;
			next();   // pg 246 2nd paragraph
			return v;
		}
		case Kind::minus:		// unary minus
			return -prim(true);
		case Kind::lp:
		{
			auto e = expr(true);
			if (ct.kind != Kind::rp) throw runtime_error("expected )");
			next();		// eat rp
			return e;
		}
		default:
			throw runtime_error("primary expected");
	}
}

bool
Adeqn::
matrix_element(const string& name, Ad& result) {
	Matrix* mp = Matrix::find_mid(name);
	if (mp == nullptr)
		return false;
	char ch;
	ip.get(ch);
	if (ch != '[')
		throw runtime_error(vastr("missing [i,j] following ",name));
	int i, j;
	ip >> i; // 1b
	ip.get(ch);
	if (ch != ',')
		throw runtime_error(vastr("missing comma in [i,j] following ",name));
	ip >> j; // 1b
	ip.get(ch);
	if (ch != ']')
		throw runtime_error(vastr("missing ] following ",name));
	int nr = mp->rsize();
	int nc = mp->csize();
	if (i < 0 || i >= nr || j < 0 || j >= nc)
		throw runtime_error(vastr("illegal element (",i+1,',',j+1,") for ",name));
	vector<complex<Ad>> mat(nr*nc);
	mp->eval(params, mat);
	complex<Ad> crval = mat[i-1+nr*(j-1)]; // operator() 1b
	result = crval.real();
	return true;
}

bool
Adeqn::
dfval(Ad& result) {
// parse a call to dfval, declared as
//   dfval(string& dfid, int gcno, vector<double>& userpar)
	T_(Trace trc(2,"dfval");)

	char ch;
	next();
	if (ct.kind != Kind::lp) throw runtime_error("expected (");

	string dfid;
	int gcno;
	vector<double> userpar;
	// dfid
	next(true);
	dfid = ct.text_value;
	T_(trc.dprint("got dfid ",dfid);)

	next();
	if (ct.kind != Kind::comma)
		throw runtime_error("expected comma after dfid");

	// gcno
	next();
	if (ct.kind != Kind::number)
		throw runtime_error("expected gc number");
	gcno = static_cast<int>(ct.number_value.value());
	next();   // eat the comma
	if (ct.kind != Kind::comma)
		throw runtime_error("expected comma after gcno");
	T_(trc.dprint("got gcno ",gcno);)

	// vector<double> userpar or just a vector initializer:
	// {0.001, 2.0} XXX only allow initializer in an equation
	do { ip.get(ch); } while (isspace(ch));
	if (ch != '{')
		throw runtime_error("{double,...} expected");
	while(true) {
		userpar.push_back(expr(true).value());  // only take the value
		// no need to get next token: see pg 246 2nd paragraph
		if (ct.kind == Kind::rbrace) {
			break;
		} else if (ct.kind != Kind::comma) {
			throw runtime_error("expected comma");
		}
	}
	T_(trc.dprint("got ",userpar.size()," userpar");)
	result = ::dfval(params, dfid, gcno, userpar);
	T_(trc.dprint("returning dfval ",result);)
	return true;
}

bool
Adeqn::
adfcn(const string& name, Ad& result) {
// If "name" refers to a function taking 1 or more Ad
// arguments, get the arguments from expr(), evaluate the
// function, and return the result in "result"

	// dfval requires special handling
	if (name == "dfval")
		return dfval(result);

	using fcn1p = Ad(*)(const Ad&);
	using fcn2p = Ad(*)(const Ad&, const Ad&);
	using item1 = std::pair<string, fcn1p>;
	using item2 = std::pair<string, fcn2p>;
	static map<string, fcn1p> fcn1 {			// add more functions here...
			item1("abs",abs),
			item1("cos",cos),
			item1("sqrt",sqrt),
			item1("atmos::rho", atmos::rho),
			item1("atmos::press", atmos::press),
			item1("atmos::temp", atmos::temp),
			item1("atmos::vsonic", atmos::vsonic),
	};
	static map<string, fcn2p> fcn2 {			// ... or here
		item2("atan2",atan2), item2("pow",pow),
		item2("atmos::cdpress",atmos::cdpress),
		item2("atmos::vcas",atmos::vcas)
	};
	// builtin describing functions (df.h):
	// 2 arguments: freeplay(delta, gcno)
	using dffcn2p = Ad(*)(pset&, const Ad&, int);
	using dfitem2 = std::pair<string, dffcn2p>;
	static map<string, dffcn2p> dffcn2 {
		dfitem2("freeplay",df::freeplay)
	};
	// 3 arguments: bilinear(delta, ratio, gcno)
	using dffcn3p = Ad(*)(pset&, const Ad&, const Ad&, int);
	using dfitem3 = std::pair<string, dffcn3p>;
	static map<string, dffcn3p> dffcn3 {
		dfitem3("bilinear",df::bilinear)
	};

	// first check to see if it is a reference atmos property
	if (name == "rhoref") {
		result = atmos::rhoref();
		return true;
	}
	if (name == "spressref") {
		result = atmos::spressref();
		return true;
	}

	// check to see if it is a recognized function, i.e. check
	// fcn1, fcn2, dffcn2, dffcn3
	if (fcn1.find(name) == fcn1.end() &&   // add fcn3, etc here
			fcn2.find(name) == fcn2.end() &&
			dffcn2.find(name) == dffcn2.end() &&
			dffcn3.find(name) == dffcn3.end()) return false;

	// get arguments assuming the string now looks like
	//   (expr,...,expr)
	next();   // get the left parentheses (lp)
	if (ct.kind != Kind::lp)
		throw runtime_error("( expected");
	vector<Ad> args;
	while(true) {
		args.push_back(expr(true));
		// no need to get next token: see pg 246 2nd paragraph
		// next(); 	// get a comma or rp
		if (ct.kind == Kind::rp) {
			// next(); 	// eat rp
			break;
		} else if (ct.kind != Kind::comma) {
			throw runtime_error("expected comma");
		}
	}

	// search for function according to argument count
	if (args.size() == 1) {
		auto f = fcn1.find(name);
		if (f == fcn1.end())
			throw runtime_error(vastr(name," has not been implemented"));
		result = f->second(args[0]);
		return true;
	} else if (args.size() == 2) {
		auto f = fcn2.find(name);
		auto df = dffcn2.find(name);
		if (f != fcn2.end())
			result = f->second(args[0], args[1]);
		// 2-argument (not counting the pset) builtin df
		else if (df != dffcn2.end())
			result = df->second(params, args[0], static_cast<int>(args[1].value()));
		else
			throw runtime_error(vastr(name," has not been implemented"));
		return true;
	} else if (args.size() == 3) {
		// 3-argument (not counting the pset) builtin df
		auto f = dffcn3.find(name);
		if (f == dffcn3.end())
			throw runtime_error(vastr(name," has not been implemented"));
		result = f->second(params, args[0], args[1], static_cast<int>(args[2].value()));
		return true;
	} else {
		throw runtime_error(vastr(args.size(),"-argument functions not implemented"));
	}
}

Token
Adeqn::
next(bool nofcn) {
// get the next token from the input stream
	char ch;
	// skip whitespace except '\n'
	do {
		if (!ip.get(ch)) return ct = {Kind::end};
	} while (ch != '\n' && isspace(ch));

	switch (ch) {
		case '*': case '/': case '+': case '-': case '(':
		case ')': case '=': case ',': case '^' : case '{' : case '}':
			return ct = {static_cast<Kind>(ch)};
		case '0': case '1': case '2': case '3': case '4': case '5': case '6':
		case '7': case '8': case '9': case '.':
			ip.putback(ch);		// put the first digit back into the input stream...
			double x;
			ip >> x;		// ... then input the number
			ct.number_value = x;
			ct.kind = Kind::number;
			return ct;
		default:
			// get a name
			if (isalpha(ch)) {
				string name{ch};
				while(ip.get(ch))
					if (isalnum(ch) || ch == ':' || ch == '.')
						name += ch;
					else {
						ip.putback(ch);
						break;
					}
				// the name may be a function...
				if (nofcn) {
					ct.kind = Kind::text;
					ct.text_value = name;
					return ct;
				}
				if (adfcn(name, ct.number_value)) {
					ct.kind = Kind::number;
					return ct;
				}
				// ... a parameter, ...
				Par* pp = params.findp(name);
				if (pp != nullptr) {
					pp->eval(params);
					ct.number_value = pp->advalue();
					ct.kind = Kind::number;
					return ct;
				} else if (evparse(name) >= 0) {
					// ok if it is an eigenvector component - maybe not defined yet
					ct.number_value = Ad(0.0);
					ct.kind = Kind::number;
					return ct;
				}
				// ...or a matrix element
				Ad result;
				if (matrix_element(name, result)) {
					ct.number_value = result;
					ct.kind = Kind::number;
					return ct;
				}
				throw runtime_error(vastr("unrecognized name: ",name));
			}
			throw runtime_error(vastr("unrecognized character: ",ch));
	}
}

Ad
adeqn::
eval (pset& plt, const string& eqn, bool noexc) {
	T_(Trace trc(2,"adeqn::eval ",eqn);)
	Adeqn eq(eqn, plt, noexc);
	Ad rval;
	try {
		rval = eq.expr(true);
	} catch (runtime_error& s) {
		throw runtime_error(vastr("evaluating \"",eqn,"\": ",s.what()));
	}
	return rval;
}

std::vector<string>
adeqn::
dependson(pset& plt, string const& eqn, bool& isComplex) {
// Returns a vector of parameter names that "eqn" depends on.
// Method: turn on monitoring in "plt", call adeqn::eval with parse_only=true,
// then turn off monitoring in "plt" to retrieve the vector.
	T_(Trace trc(2,"adeqn::dependson");)

	T_(trc.dprint("evaluating \"",eqn,"\" to determine dependencies");)
	plt.monitor();
	bool noexc{true};		// do not throw exceptions
	adeqn::eval(plt, eqn, noexc);
	vector<string> rval = plt.rotinom();
	T_(trc.dprint("returning ",rval);)
	return rval;
}

#ifdef MAIN
#undef MAIN

int
main (int argc, char *argv[]) {
	vector<string> pnames{"p1", "p2"};
	gpset::get().add(new Par("p1(Test Par 1)","1"));
	gpset::get().add(new Par("p2(Test Par 2)","3"));
	Ad::initialize(pnames);
	gpset.get().realloc();

	string eqn{"3.0*sqrt(p1 + p2)*atan2(2.0,2.0)"};
	// take the first cmd line argument if given
	if (argc > 1)
		eqn = argv[1];
	// else {
	// 	cerr << "Usage: adeqn equation\n";
	// 	exit(1);
	// }


	Adeqn eq(eqn, gpset::get());
	Ad res = eq.expr(true);
	cout << eqn << " = " << res << endl;
}

#endif // MAIN

