//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
/*------------------------------------------------------------------
 * Create a gyro matrix comprising one or more spinning rigid bodies
 *
 *  gyro {
 *		modes=modes, modal
 *		name=gyro, o=gyro
 *    node = 212500 {
 *			inertia=1.0,
 *			orient=(1,0,0),
 *			spin=1.0,
 *		}...
 *	}
 *
 *	modes=mid   name of the modes matrix to be used in reducing
 *	            the gyro matrix to modal dof.
 *	            default: the output matrix is in nodal dof
 *
 *	o=mid       Name of the output matrix
 *	            Default: o=gyro
 *
 * node = node_number{ options }
 *	            Node number where the gyro forces are to act, followed
 *	            by options for this rigid body; the options are:
 *
 *    spin     spin rate possibly followed by units. If units
 *              are not included, rad/sec is assumed. For example,
 *                 spin = 2
 *                 spin = 2 rad/sec
 *              are equivalent.
 *              Default: 1 rad/sec
 *
 *    inertia  Moment of inertia in kg-m^2 (N-m-s^2), lbm-in^2,
 *    moi      lbf-in-s^2, or slug-ft^2. Units may follow the
 *             number: k, N, lbm, lbf, slug-ft^2, or slug-in^2 (default: k)
 *             Default: 1 k
 *
 *    orient   x, y, and z coordinates in the global reference
 *             frame, relative to the node and on the spin axis.
 *             This defines the orientation and direction of spin.
 *             For example
 *                   orient=(-1,0,0)
 *             defines a spin vector in the negative x direction,
 *             so by the right-hand rule the spin direction would
 *             be clockwise when viewed down the x-axis.
 *             It is not necessary to specify trailing zero coordinates,
 *             so this example is equivalent to
 *                orient=-1
 *
 *    for example:
 *       node = 212500 {inertia=1, orient=(1,0,0), spin=100rpm}...
 *
 * Method
 *   if a mass is spinning with an angular momentum vector h,
 *   and h passes through a node with angular displacements r,
 *   the moment acting on the node due to angular motion of h
 *   is given by
 *       m = dh/dt
 *   ie the rate of change of the angular momentum vector
 *       m = dr/dt x h = - h x dr/dt
 *   where x is the vector cross-product.
 *   Let H be the matrix form of the vector cross-product:
 *       h x = H
 *             |   0   -h_3    h_2 |
 *           = |  h_3    0    -h_1 |
 *             | -h_2   h_1     0  |
 *   then the moment vector is
 *             |   0    h_3   -h_2 |
 *         m = | -h_3    0     h_1 | dr/dt
 *             |  h_2  -h_1     0  |
 *   and the gyro matrix is H with units of I\omega:
 *      kg-m^2/sec = N-m-s or lbf-in-s
 *   The angular displacements are assumed to be in radians.
 *------------------------------------------------------------------*/


#include <assert.h>
#include <cmath>
#include <sstream>
#include <string>
#include <vector>

#include "config.h"
#include "blas.h"
#include "exim.h"
#include "lexer.h"
#include "matrix.h"
#include "trace.h"

using namespace std;

class GyroElement {
public:
	int node{0};              // node number where the gyro acts
	int dof[3]{0,0,0};        // 1b rows/cols where "elem" goes in the gross mtx
	double moi{0};            // Moment of inertia in kg-m^2 (N-m-s^2)
	vector<double> orient;    // spin orientation (x,y,z) relative to node
	vector<double> elem;      // 3x3 gyro element
	double spin{0};           // spin rate in rad/sec

	// default constructor
	GyroElement();
	// default copy constructor, assignment operator ok
};

class GSpecs {
public:
	string output;
	bool SI{true};
	Matrix* modes{nullptr};
	bool modal{false};
	vector<GyroElement> gelem;
	Matrix* free{nullptr};
};
/*------------------------------------------------------------------
 * scale factors
 * Reference: The International System of Units, Physical Constants
 *  and Conversion Factors, E.A. Mechtly, NASA SP-7012,
 *  Second Revision, 1973
 *------------------------------------------------------------------*/
constexpr double Massconv{32.174048556*12.0}; // lbm-in/lbf-s^2

ostream&
operator<<(ostream& s, const GyroElement& ge) {

	s << "  node ............ " <<  ge.node << endl;
	s << "  spin orientation relative to the node:\n";
	s << "    x ............  " << ge.orient[0] << endl;
	s << "    y ............  " << ge.orient[1] << endl;
	s << "    z ............  " << ge.orient[2] << endl;
	s << "  moment of inertia " << ge.moi << " kg-m^2 (N-m-s^2)" << endl;
	s << "  element gyro matrix:\n";
	int i;
	for (i=0; i<3; i++)
		s << "      " << ge.elem[i] << "  " << ge.elem[i+3] << "  "
			<< ge.elem[i+6] << endl;
	for (i=0; i<3; i++)
		if (ge.dof[i] > 0)
			break;
	if (i != 3) {
		vector<string> rotdof{"RX", "RY", "RZ"};
		s << "  element gyro matrix goes in rows/columns:\n";
		for (i=0; i<3; i++) {
			if (ge.dof[i] > 0) {
				s << "   " << rotdof[i] << " .............. " << ge.dof[i] << endl;
			}
		}
	}
	return s;
}

class Sparse {
public:
	vector<double> nz;
	vector<int> ri;
	vector<int> cs;

	Sparse() {}
	Matrix* to_full(const string& mid);
};

ostream&
operator<<(ostream& s, const Sparse& sp) {
	s << "row indices:\n" << sp.ri << endl;
	s << "column starts:\n" << sp.cs << endl;
	s << "non-zeros:\n" << sp.nz << endl;
	return s;
}

Matrix*
Sparse::
to_full(const string& mid) {
	int n = cs.size()-1;
	Matrix* rval = new Matrix(mid, "Gyro", n, n, false);
	double* el = rval->elem();
	for (int j=0; j<n; j++) {
		int first = cs[j];
		int last = cs[j+1];
		for (int k=first-1; k<last-1; k++) {
			int i = ri[k] - 1;
			el[i+j*n] = nz[k];
		}
	}
	return rval;
}

// Private functions
static bool parser (string const& optionList, GSpecs& specs);
static Sparse* gyroMatrix (GSpecs& specs);

bool
gyro (string const& optionlist) {
// main entry point for the Flaps "gyro" command
	GSpecs specs;

	try {
		if (!parser(optionlist, specs))
			return false;

		// Get the freedoms matrix: doubles like node.freedom
		specs.free = new Matrix("freedoms");

		// Create the gyro matrix in sparse format...
		Sparse* sparsegm = gyroMatrix(specs);
		if (sparsegm == nullptr)
			return false;

#ifdef NEVER // disable for now
		// ... and reduce it to modal dof if requested
		if (Modes != nullptr) {
			NumMat* tmp = sparsegm->mult("n", "n", Modes);
			Modes->transpose();
			result = Modes->mult("n", "n", tmp);
			result->mid = string(Name);
			flaps::info ("The output gyro matrix is in modal"
					" (generalized) coordinates associated with "
					"modes matrix \"",Modes->mid(),"\"");
			delete tmp;
			delete Modes;
			Modes = 0;
		} else {
			flaps::info ("The output gyro matrix is in nodal coordinates");
			result = sparsegm;
		}
#endif // NEVER : disable for now

		Matrix* result = sparsegm->to_full(specs.output);

		if (result != nullptr) {
			flaps::info ("Output: ", result->summary());
			MM::exporter(result->mid()+".mm", result);
			result->store();
			delete result;
		}
	} catch (runtime_error const& s) {
		flaps::error(s.what());
		return false;
	}

	return true;
}

// parsing functions
static bool node (const Tok& opt, GSpecs& sp);
static bool modes (const Tok& opt, GSpecs& sp);

static bool
parser (string const& options, GSpecs& specs) {
	T_(Trace trc(1,"parser");)
	bool rval = true;

	// Note: some options are parsed in node() since
	// they are associated with a node
	vector<Tok*> unrec = flaps::lexer(options, {
		{"^node(s)?", [&](const Tok& p) { return node(p,specs); }},
		{"^o(ut)?", [&](const Tok& p) { specs.output = p.srhs; return true; }},
		{"^output?", [&](const Tok& p) { specs.output = p.srhs; return true; }},
		{"^modes",[&](const Tok& p) { return modes(p,specs); }},
		{"^units",[&](const Tok& p) {
			if (p.srhs == "uscs") {
				specs.SI = false;
				return true;
			}
			return false;
		}}
	});

	if (!unrec.empty())
		flaps::warning(vastr(unrec.size()," specs were not recognized: ",unrec));

	// Defaults
	if (specs.output.empty()) {
		specs.output = "gyro";
		flaps::info ("The name of the output matrix is "
			"defaulted to \"",specs.output,"\"");
	}

	return rval;
}

static bool
modes (const Tok& opt, GSpecs& sp) {
	T_(Trace trc(1,"modes");)
	if (opt.srhs.empty()) {
		T_(trc.dprint("output gyro matrix will be in modal dof");)
	} else {
		sp.modes = new Matrix(opt.srhs);
		T_(trc.dprint("got Modes Matrix <",sp.modes->mid());)
	}
	sp.modal = true;
	return true;
}

double
moi_parse(const string& rhs, bool si) {
// parse a moment-of-inertia description: a float optionally
// followed by units. Valid units are:
//    lbm:        lbm-in^2
//    slug-ft^2   slug-ft^2 
//    slug        slug-in^2 
// Units of the return value are SI (kg-m^2 or N-m-s^2) if si==true,
// USCS otherwise (lbf-in-s^2)
	size_t end;
	double rval{0.0};
	if (!str2double(rhs, rval, 0, end))
		throw runtime_error(vastr("illegal moi spec: ",rhs));

	// followed by units?
	if (end == string::npos)
		return rval;

	string units = stripwhitespace(rhs.substr(end));
	if (units.empty())
		return rval;

	if (si) {
		constexpr double kg_lbm{0.45359237};
		constexpr double m_in{0.0254};
		constexpr double m_ft{0.3048};
		constexpr double kg_slug{14.594};
		if (units == "lbm") {
			rval *= kg_lbm*m_in*m_in;  // lbm-in^2 -> kg-m^2
		} else if (units == "slug-ft^2") {
			rval *= kg_slug*m_ft*m_ft;  // slug-ft^2 -> kg-m^2
		} else if (units == "slug-in^2") {
			rval *= kg_slug*m_in*m_in;  // slug-in^2 -> kg-m^2
		} else
			throw runtime_error(vastr("unrecognized moi units: ",units));
	} else {
		constexpr double lbm_lbf{386.0886};  // lbm-in/lbf-s^2
		constexpr double lbf_slug{1.0/12.0}; // (lbf-s^2/in)/slug
		if (units == "lbm") {
			rval /= lbm_lbf;
		} else if (units.substr(0,4) == "slug") {
			rval *= lbf_slug;    // slug-in^2 -> lbf-in-s^2
		} else
			throw runtime_error(vastr("unrecognized moi units: ",units));
	}
	return rval;
}

static bool
node (const Tok& opt, GSpecs& sp) {
/*
 * Get one or more node numbers and create a GyroElement for each node.
 * Each node number must be followed by a set of options describing the
 * gyro forces at that node. At least the orientation option must be given;
 * other options have defaults.
 * The general form of this option is:
 *    node = (node_number{ specs }, node_number{ specs }, ...)
 * for example:
 *    node = 212500 {
 *			inertia=1.0,
 *			orient=(1,0,0),
 *			spin=1.0,
 *		}...
 *	The specs for each node are:
 *    spin     spin rate possibly followed by units. If units
 *             are not included, rad/sec is assumed. For example,
 *                spin = 2
 *                spin = 2 rad/sec
 *                spin = 19.1 rpm
 *             are equivalent.
 *             Default: 1 rad/sec
 *
 *    inertia  Moment of inertia in lbf-in-sec**2
 *    moi
 *    moment
 *             Default: 1 lbf-in-sec**2
 *
 *    orient   x, y, and z coordinates in the global reference
 *             frame, relative to the node and on the spin axis.
 *             This defines the orientation and direction of spin.
 *             For example
 *                orient=(-1,0,0)
 *             defines a spin vector in the negative x direction,
 *             so by the right-hand rule the spin direction would
 *             be clockwise when viewed down the x-axis.
 *             It is not necessary to specify trailing zero coordinates,
 *             so this example is equivalent to
 *                orient=-1
 *
 */
	T_(Trace trc(1,"node");)

	assert(opt.ivec.size() == opt.roptvec.size());
	
	for (size_t i=0; i<opt.ivec.size(); i++) {
		GyroElement np;
		np.node = opt.ivec[i];
		// Parse the node options:
		//   (moment|inertia|moi)=1, orient=(x,y,z), spin = 1 rad/sec,
		// Specs (lambdas)
		// moi = float [lbm | lbf | slugs] default: kg-m^2
		// 1/386 lb_f-in-sec^2/lb_m-in^2
		// 12 lb_f-in-sec^2/slug-ft^2
		vector<Tok*> unrec = flaps::lexer(opt.roptvec[i], {
			{"moi", [&](const Tok& p) {
				np.moi = moi_parse(p.srhs, sp.SI);
				return true;
			}},
			{"orient", [&np](const Tok& p) { np.orient = p.rvec; return true; }}
		});

		if (!unrec.empty())
			throw runtime_error(vastr(unrec.size()," unrecognized specs: ",unrec));

		if (np.orient.empty())
			throw runtime_error("at least one orientation component must be included");

		T_(trc.dprint("got gyro element: ",np);)
		sp.gelem.push_back(np);
	}
	return true;
}

GyroElement::
GyroElement() {
// default constructor
	node = 0;
	moi = 0.0;
	orient = vector<double>(3, 0.0);
	elem = vector<double>(9, 0.0);
	spin = 0.0;
}

class Nz {
public:
	int i, j; // 1b row & column indices
	double x;

	Nz() { i = j = 0; x = 0.0; }
	Nz(int ii, int jj, double xx) : i(ii), j(jj), x(xx) {}
	bool operator==(Nz const& rhs) const { return (i == rhs.i && j == rhs.j); }
};

ostream&
operator<<(ostream& s, const Nz& nz) {
	s << "(" << nz.i << ", " << nz.j << "), nz = " << nz.x;
	return s;
}

static Sparse*
gyroMatrix (GSpecs& specs) {
/*------------------------------------------------------------------
 * Assemble all gyro elements into a Matrix
 * An element gyro matrix is formed by the relation
 *    .
 *    H = O x H
 * where O is the angular velocity and H is the moment of
 * inertia
 *   H = Io
 * where I is the inertia tensor and o is the spin vector.
 *------------------------------------------------------------------*/
	T_(Trace trc(0,"gyroMatrix");)
	double scale;
	double x[3];
	size_t i, j, k;
	vector<double> freedoms = specs.free->data();
	// the return matrix will be square with order n
	size_t n = freedoms.size();
	
	T_(trc.dprint("creating ",specs.output,", nelem: ",specs.gelem.size(),", freedoms:",freedoms);)

	// Allocate sparse structures large enough to hold all gyro elements
	// assuming all gyro elements are full 3x3 matrices
	vector<int> colnnz(n, 0);  // number of non-zeros in each column in rval
	vector<Nz> nzs;            // row/col (0b) of each non-zero in rval

	for (auto& ge : specs.gelem) {
		scale = 0.0;
		for (j=0; j<3; j++) {
			x[j] = ge.orient[j];
			scale += x[j]*x[j];
		}
		if (scale + 1.0 == 1.0)
			throw runtime_error(vastr("the orientation point for node ",
						ge.node, " is zero"));

		scale = sqrt(scale);

		if (ge.moi == 0.0) {
			flaps::info ("the moment of inertia for node ",
				ge.node, " is defaulted to 1.0 lbf-in-sec**2 (",
				Massconv, " lbm-in**2) (",1.0/12.0, "slugs-ft^2)");
			ge.moi = 1.0;
		}
		scale = ge.moi/scale;

		if (ge.spin == 0.0) {
			flaps::info("the spin rate for node ", ge.node,
				" is defaulted to 1.0 rad/sec");
			ge.spin = 1.0;
		}
		scale *= ge.spin;
			
		for (j=0; j<3; j++)
			x[j] *= scale;

		// form the 3,3 matrix of gyro terms for the 3 rotational dof
		// at this node
		ge.elem[0] = 0.0;			/* 1,1 */
		ge.elem[1] = -x[2];			/* 2,1 */
		ge.elem[2] = x[1];			/* 3,1 */
		ge.elem[3] = x[2];			/* 1,2 */
		ge.elem[4] = 0.0;			/* 2,2 */
		ge.elem[5] = -x[0];			/* 3,2 */
		ge.elem[6] = -x[1];			/* 1,3 */
		ge.elem[7] = x[0];			/* 2,3 */
		ge.elem[8] = 0.0;				/* 3,3 */

		// Determine which rows/columns get the element.
		for (j=0; j<freedoms.size(); j++) {
			double nn;
			double freedom = modf(freedoms[j], &nn);
			int inode = nn;
			int ifree = 10.0*freedom;
			if (inode == ge.node && (ifree > 3 && ifree <= 6)) {
				ge.dof[ifree-4] = j+1;
				T_(trc.dprint("freedom ",ifree,", node ",inode," is row ",j+1," 1b");)
			}
		}
		if ((ge.dof[0] + ge.dof[1] + ge.dof[2]) == 0)
			throw runtime_error(vastr("node ", ge.node,
				" has no retained rotational freedoms"));

		// create an NZ struct for each term of this ge that
		// corresponds to a retained dof and is non-zero
		for (j=0; j<3; j++) {
			if (ge.dof[j] == 0)
				continue;
			for (i=0; i<3; i++) {
				double xij = ge.elem[i+j*3];
				if (ge.dof[i] == 0 || xij == 0.0)
					continue;
				Nz nzij(ge.dof[i], ge.dof[j], xij);
				vector<Nz>::iterator p = find(nzs.begin(), nzs.end(), nzij);
				if (p == nzs.end()) {
					nzs.push_back(nzij);
					colnnz[ge.dof[j]-1]++;
					T_(trc.dprint("got nzs[",nzs.size(),"] = ",nzij);)
				} else {
					T_(trc.dprint("adding ",nzij.x," to ",*p);)
					p->x += nzij.x;
				}
			}
		}
		cout << ge << endl;
	}

	size_t nnz = nzs.size();

	T_(trc.dprint("nnz = ",nnz);)
	T_(trc.dprint("nnz in each column: ",colnnz);)

	vector<double> nz(nnz,0.0);
	vector<int> ri(nnz,0);
	vector<int> cs(n+1,0);

	// Now go through each column, count the number of non-zeros
	// and set the column start (1b) for the next column
	cs[0] = 1;
	for (j=0; j<n; j++)
		cs[j+1] = cs[j] + colnnz[j];

	// Insert elements in nz and ri by getting the column-start,
	// Let colnnz[j] keep track of how many elements we have inserted
	// into column j so far.
	for (j=0; j<n; j++)
		colnnz[j] = 0;

	for (k=0; k<nnz; k++) {
		int col = nzs[k].j - 1;	// 0-based
		int start = cs[col] - 1; // 0b
		nz[start+colnnz[col]] = nzs[k].x;
		ri[start+colnnz[col]] = nzs[k].i; // 1-based
		colnnz[col]++;
	}

	// Be sure the lower-right diagonal is included
	if (ri[nnz-1] != (int)n) {
		ri.push_back(n);
		nz.push_back(0.0);
		cs[n]++;
		nnz++;
	}

	// Create a Sparse and copy stuff
	string desc("Gyro Matrix");
	Sparse* rval = new Sparse();
	rval->ri = ri;
	rval->cs = cs;
	rval->nz = nz;

	T_(trc.dprint("returning: ",*rval);)

	return rval;
}
