//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// minimal beam model builder with gyroscopics, unsteady propeller
// aerodynamics, and mass matrices

#include "config.h"

#include <assert.h>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "config.h"
#include "Ad.h"
#include "blas.h"
#include "conv.h"
#include "exim.h"
#include "fem.h"
#include "lapack.h"
#include "lexer.h"
#include "main.h"  // setup signal, debugger handling
#include "matrix.h"
#include "Pz.h"
#include "Regex.h"
#include "settings.h"
#include "trace.h"
#include "version.cpp"

using namespace std;

bool parser (string const& progname, Fem& model);

int
main (int argc, char** argv) {
	string progname{argv[0]};
	Trace trc(1,progname, ":main");
	Fem& model = Fem::instance();
	try {
		parser(progname, model);

		// assemble the various elements into the output matrices
		model.assemble();
		// export it as a .mm file for checking with matview
		model.exportmm();
		// create the data needed for animated-mode visualization
		model.vzdata();
		// do modal reduction if requested, after vzdata so gct gets multiplied
		// BMA
		if (model.bma != nullptr) {
			model.bma->xform = model.bma->transform(model.retained());
			model.modal_reduction(model.bma->xform, "bma");
		}
		// CMS
		if (model.cms != nullptr) {
			model.cms->xform = model.cms->transform(model.retained());
			model.modal_reduction(model.cms->xform, "cms");
		}
		// free vibration of the gross model
		if (!model.modal().empty()) {
			// if a model id has not been specified give it "modal"
			string id{model.id()};
			if (id.empty())
				id = "modal";
			model.modal_reduction(model.modal(), id);
		}
	} catch (runtime_error& s) {
		flaps::error(s.what());
	}

	return flaps::nerrors();
}

// local functions
static double
norm2(const vector<double>& x) {
	return blas_snrm2(x.size(), &x[0], 1);
}

bool
parser (string const& progname, Fem& model) {
	Trace trc(1,"parser");
	Freev freev;

	vector<Tok*> unrec = flaps::lexer("", {
		{"beam(s)?", [&](const Tok& p) {
			vector<Beam*> beams = parse_beam(p,model.matl,freev);
			for (auto bp : beams)
				model.stif_elements.push_back(bp);
			for (auto fi : freev)
				model.retained().add(fi);
			return true;
		}},
		{"bma", [&](const Tok& p) { return model.bma_f(p); }},
		{"cms", [&](const Tok& p) { return model.cms_f(p); }},
		{"conmass", [&](const Tok& p) {
			return parse_conmass(p,model.retained(), model.stif_elements, model.mass_elements);
		}},
		{"dlm", [&](const Tok& p) {return model.dlm_f(p);}},
		{"gyro", [&](const Tok& p) {return model.gyro_f(p);}},
		{"id", [&](const Tok& p) { model.id(p.srhs); return true;}},
		{"nmodes", [&](const Tok& p) { model.freevib(true);
			return model.nmodes(p.ivec[0]); }},
		{"node(s)?", [&](const Tok& p) {return model.nodes_f(p);}},
		{"material", [&](const Tok& p) {model.matl = parse_material(p); return true;}},
		{"pgaf", [&](const Tok& p) {return model.pgaf_f(p);}},
		{"ss", [&](const Tok& p) { return model.ss_f(p); }},
		{"units", [&](const Tok& p) {
			if (compare(p.srhs, "si")) model.si_units = true;
			else if (compare(p.srhs, "uscs")) model.si_units = false;
			else throw runtime_error(vastr("unrecognized units: ",p.srhs));
			return true;
		}},
	});

	if (!unrec.empty())
		flaps::error("unrecognized specs: ",unrec);

	return true;
}

bool
par_or_value(const Tok& p, string& param, double& value) {
// Check a spec for either a parameter name or a value
// If it is a parameter set "param" and set "value" to zero,
// if is is a value set "value" and set "param" to empty
	if (p.is_char()) {
		Par* pp = gpset::find(p.srhs);
		if (pp == nullptr)
			throw runtime_error(vastr("parameter \"",p.srhs,"\" is undefined"));
		param = p.srhs;
		value = 0.0;
		return true;
	} else if (p.is_real()) {
		value = p.rvec[0];
		param.clear();
		return true;
	}
	return false;
}

bool
parse_orient(const Tok& p, vector<double>& orient, vector<Beam*> beams) {
// check a Tok for various ways of specifying an orientation:
//   1) orient = x | y | z
//   2) orient = (x,y,z)    3 doubles
//   3) orient = beam_id(x,y,z)   relative to the beam with id "beam_id"
// Output:
//   orient  normalized to 1.0
	Trace trc(2,"parse_orient");
	trc.dprint("spec<",p,">");
	// lambda to normalize
	auto normalize = [&](vector<double> v) {
		double nrm = blas_snrm2(3,&v[0],1);
		if (nrm > 0.0)
			blas_scal(3, 1.0/nrm, &v[0], 1);
	};
	// initialize
	orient = vector<double>(3, 0.0);

	// 1)
	string xyz{"xyz"};
	auto idx = xyz.find(p.srhs[0]);
	if (p.srhs.size() == 1 && idx != string::npos) {
		orient[idx] = 1.0;
		trc.dprint("returning: ",orient);
		return true;
	}

	// 2)
	if (p.rvec.size() == 3) {
		orient = p.rvec;
		normalize(orient);
		return true;
	}

	// 3) check for id(x,y,z) where "id" is the id of a beam
	string s = "([a-zA-Z][^(]*)\\((" + regex_double() + "),(" + regex_double()
			+ "),(" + regex_double() + ")\\)";
	regex re = make_regex(s);
	smatch mch;
	if (regex_match(p.srhs, mch, re)) {
		string id{mch[1].str()};
		double x, y, z;
		str2double(mch[2].str(), x);
		str2double(mch[5].str(), y);
		str2double(mch[8].str(), z);
		// find beam with "id" in the model
		for (auto bi : beams) {
			if (id == bi->id()) {
				vector<double> v{x,y,z};
				orient = vector<double>(3,0.0);
				blas_sgemv("t", 3, 3, 1.0, &bi->transform()[0], 3, &v[0],
					1, 0.0, &orient[0], 1);
				normalize(orient);
				trc.dprint("returning ",orient);
				return true;
			}
		}
		throw runtime_error(vastr("beam id \"",id,"\" has not been defined"));
	}
	trc.dprint("returning false: spec not recognized");
	return false;
}

// parsing functions: end in _f or begin with parse_

bool
Fem::
nodes_f (const Tok& p) {
// create a matrix of nodal data: (3*nnodes,1) real Matrix "coordinates"
// nodal data:
//   nodes{1=(x,y,z), 2=(x,y,z), ...}
//   nodes{1{x,y,z} : 21{x,y,z}, ...}
	Trace trc(2,"nodes_f");

	if (p.lopt.empty())
		throw runtime_error("bad nodes option, should be nodes{n1=(x,y,z),...}");

	vector<Tok*> unrec = flaps::lexer(p.lopt, {
		{".*",[&](const Tok& p) {
			int nodeno;
			vector<double> coord;
			if (!str2int(p.lhs, nodeno))
				return false;
			// check for n=(x,y,z) ...
			if (!p.rvec.empty()) {
				coord = p.rvec;
			// ... or n{x,y,z}
			} else {
				if (p.lopt.empty())
					return false;
				vector<Tok*> unrec = flaps::lexer(p.lopt, {
					{".*",[&](const Tok& q) {
						if (q.rvec.size() != 3)
							return false;
						coord = q.rvec;
						return true;
					}}
				});
				if (!unrec.empty())
					throw runtime_error(vastr("unrecognized coord option: ",unrec));
			}
			this->grid.add(Node(nodeno, coord));
			return true;
		}}
	});
	if (!unrec.empty())
		throw runtime_error(vastr("unrecognized node option: ",unrec));

	trc.dprint("structure nodes:\n",this->grid);
	// sort the nodes by node number?

	return true;	
}

Material
parse_material(const Tok& p) {
// material{x=(x,y,z), eix=d, eiz=d, gj=d, ea=d}
	Trace trc(2,"parse_material");
	Material rval;
	flaps::lexer(p.lopt, {
		// global coord of a point on the local x axis:
		{"x", [&](const Tok& p) { rval.x = p.rvec; return true; }},
		// bending about the local x axis:
		{"eix", [&](const Tok& p) { rval.eix = p.rvec[0]; return true; }},
		// bending about the local z axis:
		{"eiz", [&](const Tok& p) { rval.eiz = p.rvec[0]; return true; }},
		// torsion (rotation about local y):
		{"gj", [&](const Tok& p) { rval.gj = p.rvec[0]; return true; }},
		// extension in local y direction
		{"ea", [&](const Tok& p) { rval.ea = p.rvec[0]; return true; }},
	});
	trc.dprint("got ",rval);
	return rval;
}

ostream&
operator<<(ostream& s, const Material& t) {
	s << "pt on local x axis: (" << t.x[0] << ", " << t.x[1] << ", " << t.x[2] << endl;
	s << " eix: " << t.eix << endl;
	s << " eiz: " << t.eiz << endl;
	s << " gj:  " << t.gj << endl;
	s << " ea:  " << t.ea << endl;
	return s;
}

vector<Beam*>
parse_beam (const Tok& p, const Material& matl, Freev& freev) {
// Parse the defn of a set of beam elements, e.g.:
//  beam{ 1{0}, 2{3:5}, 3:13 }
// which is equivalent to
//  beam{ 1{0}, 2{3,4,5}, 3{3,4,5}, ... 13{3,4,5} }
	Trace trc(1,"parse_beam");
	Nodedof node_dof1;
	Nodedof node_dof2;
	int ordinal{0};
	vector<int> dof;
	vector<Beam*> rval;
	// call lexer with lambda handlers
	vector<Tok*> unrec = flaps::lexer(p.lopt, {
		{".*", [&](const Tok& p) {
			int node;
			if (!str2int(p.lhs, node)) return false;
			// dof: either as node=dof or node{dof}
			if (p.lopt.empty()) {
				if (!p.ivec.empty())
					dof = p.ivec;
			} else {
				// ...or special parse treatment of a set of numbers w/no rhs
				flaps::lexer(p.lopt, {
					{".*", [&](const Tok& q) { dof = q.ivec; return true;}}
				});
			}
			if (dof.size() == 1 && dof[0] == 0)
				dof.clear();
				
			if (ordinal++ == 0) {		// first time through...
				node_dof1 = Nodedof(node,dof);
			} else {
				node_dof2 = Nodedof(node,dof);	// ... all subsequent times
				// create a new Beam, put in rval
				Beam* bp = new Beam(node_dof1,node_dof2,matl,Fem::instance().grid);
				//!! stif_elements.push_back(bp);
				rval.push_back(bp);
				// add these retained freedoms to the Fem model XXX do this in assemble
				for (auto& fp : bp->freedoms())
					freev.add(fp);
				// ready for the next beam element
				node_dof1 = node_dof2;
			}
			return true;
		}}
	});
	if (!unrec.empty())
		throw runtime_error(vastr("unrecognized beam options: ",unrec));

	trc.dprint("material properties:\n",matl,"new stif elements: ",rval);
	return rval;
}

#ifdef NEVER // deprecate
bool
Fem::
beam_f (const Tok& p) {
// XXX deprecated: use parse_beam
// Parse the defn of a set of beam elements, e.g.:
//  beam{ 1{0}, 2{3:5}, 3:13 }
// which is equivalent to
//  beam{ 1{0}, 2{3,4,5}, 3{3,4,5}, ... 13{3,4,5} }
	Trace trc(1,"beam_f");
	Nodedof node_dof1;
	Nodedof node_dof2;
	int ordinal{0};
	vector<int> dof;
	vector<Beam*> new_elements;
	vector<Tok*> unrec = flaps::lexer(p.lopt, {
		{".*", [&](const Tok& p) {
			int node;
			if (!str2int(p.lhs, node)) return false;
			// dof: either as node=dof or node{dof}
			if (p.lopt.empty()) {
				if (!p.ivec.empty())
					dof = p.ivec;
			} else {
				// ...or special parse treatment of a set of numbers w/no rhs
				flaps::lexer(p.lopt, {
					{".*", [&](const Tok& q) { dof = q.ivec; return true;}}
				});
			}
			if (dof.size() == 1 && dof[0] == 0)
				dof.clear();
				
			if (ordinal++ == 0) {		// first time through...
				node_dof1 = Nodedof(node,dof);
			} else {
				node_dof2 = Nodedof(node,dof);	// ... all subsequent times
				// create a new Beam, put in stif_elements
				Beam* bp = new Beam(node_dof1,node_dof2,matl,this->grid);
				stif_elements.push_back(bp);
				new_elements.push_back(bp);
				// add these retained freedoms to the Fem model
				for (auto& fp : bp->freedoms())
					this->add(fp);
				// ready for the next beam element
				node_dof1 = node_dof2;
			}
			return true;
		}}
	});
	if (!unrec.empty())
		throw runtime_error(vastr("unrecognized beam options: ",unrec));

	trc.dprint("material properties:\n",matl,"new stif elements: ",new_elements);
	// cout << "material properties:\n" << matl << "new stif elements: \n" << new_elements;
	return true;
}
#endif // NEVER // deprecate

Beam::
Beam (Nodedof const& nodedof1, Nodedof const& nodedof2, Material const& matl, Grid const& grid) {
// Beam constructor
// create an (nf,nf) beam element where nf is the total number of retained
// freedoms at nodes 1 and 2. The (nf,nf) matrix element goes in the data_
// member; it is in the global coordinate system so it can be merged into
// the output matrix using the freedoms_ member as the row/column keys
// in blas::copy().
	Trace trc(2,"Beam constructor");
	static vector<int> lastseen; // last dof seen
	constexpr int n{12};	// size of the matrix before extracting retained dof:
								// all 6 dof at each node

	string id;
	node1_ = nodedof1.node;
	node2_ = nodedof2.node;
	vector<int> dof1 = nodedof1.dof;
	vector<int> dof2 = nodedof2.dof;
	trc.dprint("node 1: ",node1_,", dof ",dof1);
	trc.dprint("node 2: ",node2_,", dof ",dof2);

	// get the nodal coord for these 2 nodes from model.grid
	vector<double> x1 = grid.find(node1_)->coord();
	vector<double> x2 = grid.find(node2_)->coord();

	// compute a local coordinate system (xbar,...) with ybar along the
	// beam element.
	vector<double> orient{{x2[0]-x1[0]}, {x2[1]-x1[1]}, {x2[2]-x1[2]}};
	double l = blas_snrm2(3, &orient[0], 1); // the length of this beam element
	if (l < std::numeric_limits<double>::epsilon())
		throw runtime_error(vastr("nodes ",node1_," and ",node2_," are coincident"));

	trc.dprint("orientation: ",orient);

	// ybar (beam axis) is "orient" normalized to 1
	vector<double> ybar{orient};
	blas_scal(3, 1.0/l, &ybar[0], 1);

	// xbar is matl.x normalized to 1
	vector<double> xbar{matl.x};
	double xnorm = blas_snrm2(3, &xbar[0], 1);
	if (xnorm == 0.0)
		throw runtime_error("xbar is a zero vector");
	blas_scal(3, 1.0/xnorm, &xbar[0], 1);
	// zbar is xbar X ybar
	vector<double> zbar{0.0, 0.0, 1.0};
	xprod(&xbar[0], &ybar[0], &zbar[0]);
		
	// create the beam element in the local coord system with the beam oriented
	// along the local y axis; a beam element has 12 dof in the local coordinate
	// system: tx, ty, tz, rx, ry, rz at the 2 nodes.
	// Let EIz = Young's modulus (E) times the second moment of area about
	//     the z axis (determines bending stiffness in the x direction)
	//     EIx = E times the second moment of area about the x axis (determines
	//           bending stiffness in the z direction)
	// the units of EI are N-m^2 or lb_f-ft^2
	// Let r = EIz/EIx, t = EA/EIx, s = GJ/EIx
	double EIz = matl.eiz;
	double EIx = matl.eix;
	double GJ = matl.gj;
	double EA = matl.ea;
	double r = EIz/EIx;
	double s = GJ/EIx;
	double t = EA/EIx;
	// then K = EIx/l^3 *
	/*   tx    ty    tz    rx    ry    rz    tx    ty    tz    rx    ry    rz
	 tx  12r                          -6lr  -12r                          -6lr
	 ty        tl^2                               -tl^2
	 tz             12     6l                           -12    6l 
	 rx             6l    4l^2                          -6l   2l^2 
	 ry                         sl^2                                -sl^2
	 rz -6lr                         4l^2r   6lr                          2l^2r 
	 tx -12r                          6lr    12r                           6lr
	 ty       -tl^2                                tl^2
	 tz             -12  -6l                             12   -6l 
	 rx              6l  2l^2                           -6l   4l^2 
	 ry                        -sl^2                                 sl^2
	 rz -6lr                          2l^2r  6lr                          4l^2r 
	 */
	constexpr int tx1{0}, ty1{1}, tz1{2}, rx1{3}, ry1{4}, rz1{5};
	constexpr int tx2{6}, ty2{7}, tz2{8}, rx2{9}, ry2{10}, rz2{11};
	vector<double> el(n*n, 0.0);
	// insert all 12 dof into el (upper triangle only); rows/columns
	// not retained will be zeroed out before transforming to global
	// since dof1/2 are in the (xbar,ybar,zbar) coordinate system XXX ??
	// row 1: tx1
	el[IJ(tx1,tx1,n)] = 12.0*r;
	el[IJ(tx1,rz1,n)] = -6.0*l*r;
	el[IJ(tx1,tx2,n)] = -12.0*r;
	el[IJ(tx1,rz2,n)] = -6.0*l*r;
	// row 2: ty1
	el[IJ(ty1,ty1,n)] = t*l*l;
	el[IJ(ty1,ty2,n)] = -t*l*l;
	// row 3: tz1
	el[IJ(tz1,tz1,n)] = 12.0;
	el[IJ(tz1,rx1,n)] = 6.0*l;
	el[IJ(tz1,tz2,n)] = -12.0;
	el[IJ(tz1,rx2,n)] = 6.0*l;
	// row 4: rx1
	el[IJ(rx1,rx1,n)] = 4.0*l*l;
	el[IJ(rx1,tz2,n)] = -6.0*l;
	el[IJ(rx1,rx2,n)] = 2.0*l*l;
	// row 5: ry1
	el[IJ(ry1,ry1,n)] = s*l*l;
	el[IJ(ry1,ry2,n)] = -s*l*l;
	// row 6: rz1
	el[IJ(rz1,rz1,n)] = 4.0*l*l*r;
	el[IJ(rz1,tx2,n)] = 6.0*l*r;
	el[IJ(rz1,rz2,n)] = 2.0*l*l*r;
	// row 7: tx2
	el[IJ(tx2,tx2,n)] = 12.0*r;
	el[IJ(tx2,rz2,n)] = 6.0*l*r;
	// row 8: ty2
	el[IJ(ty2,ty2,n)] = t*l*l;
	// row 9: tz2
	el[IJ(tz2,tz2,n)] = 12.0;
	el[IJ(tz2,rx2,n)] = -6.0*l;
	// row 10: rx2
	el[IJ(rx2,rx2,n)] = 4.0*l*l;
	// row 11: ry2
	el[IJ(ry2,ry2,n)] = s*l*l;
	// row 12: rz2
	el[IJ(rz2,rz2,n)] = 4.0*l*l*r;

	// symmetrize
	for (int j=0; j<n; j++)
		for (int i=j+1; i<n; i++)
			el[IJ(i,j,n)] = el[IJ(j,i,n)];

	// scale the entire matrix by EIx/l^3
	blas_scal(n*n, EIx/(l*l*l), &el[0], 1);

	trc.dprintm(n,n,n,el,"beam element in local ",node1_,"-",node2_);

	// the beam is oriented along the local y axis (ybar); transform so that
	// it is in "global" coordinates; the transformation to global
	// coordinates is xbar, ybar and zbar as ROWS so we can use lapack::triprod
	tlg_ = vector<double>(9,0.0);
	blas_copy(3, &xbar[0], 1, &tlg_[0], 3);
	blas_copy(3, &ybar[0], 1, &tlg_[1], 3);
	blas_copy(3, &zbar[0], 1, &tlg_[2], 3);
	trc.dprintm(3,3,3,tlg_,"transformation from local to Fem");

	if (!is_identity(tlg_)) {
		// copy the (3,3) tlg_ into a (12,12) transformation matrix T
		// T is temporary, tlg_ is a member which may be needed later
		vector<double> T(n*n, 0.0);
		for (int j=0; j<3; j++) {
			for (int i=0; i<3; i++) {
				T[IJ(i,j,n)] = tlg_[IJ(i,j,3)];
				T[IJ(i+3,j+3,n)] = tlg_[IJ(i,j,3)];
				T[IJ(i+6,j+6,n)] = tlg_[IJ(i,j,3)];
				T[IJ(i+9,j+9,n)] = tlg_[IJ(i,j,3)];
			}
		}
		trc.dprintm(n,n,n,T,"beam",node1_,"-",node2_," transformation matrix");
		
		// transform the element in local to global: T'*el*T -> el
		vector<double> te(n*n,0.0);
		lapack::triprod(n, n, &T[0], &el[0], &te[0]);
		el = te;
	}

	trc.dprintm(n,n,n,el,"beam element in global ",node1_,"-",node2_);

	// the retained dof for this beam element: dof1 + dof2
	for (auto i: dof1)
		freedoms_.push_back(Freedom(node1_,i));
	for (auto i: dof2)
		freedoms_.push_back(Freedom(node2_,i));
	int mr = freedoms_.size();
	// all freedoms: the freedoms for el
	vector<Freedom> all;
	for (int j=1; j<7; j++)
		all.push_back(Freedom(node1_,j));
	for (int j=1; j<7; j++)
		all.push_back(Freedom(node2_,j));
	// extract the retained freedoms from el, put in data_
	data_ = vector<double>(mr*mr,0.0);
	blas::copy(all,all,1.0,el,freedoms_,freedoms_,0.0,data_);

	trc.dprint("all dof: ",all,", retain: ",freedoms_);
	trc.dprintm(mr,mr,mr,data_,"retained beam element");
}

bool
parse_conmass (const Tok& p, const Freev& ret,
	const vector<Beam*> beams, vector<Mass*>& masses) {
// parse concentrated-mass specs
// conmass{mass=w, moi=(Ixx,Iyy,Izz), cg=c, orient=(x,y,z)), node=(first:last),...}
// where w is the mass (Kg) and I is the moment of inertia (Kg-m^2)
	Trace trc(2,"parse_conmass");
	double totalmass{0.0};
	double mass{0.0};
	vector<double> totalmoi;
	vector<double> moi{3,0.0};		// default: zero
	string cgpar;
	double cg{0.0};
	vector<double> orient(3,0.0);
	vector<int> nodes;

	vector<Tok*> unrec = flaps::lexer(p.lopt, {
		{"totalmass", [&](const Tok& p) {totalmass = p.rvec[0]; return true;}},
		{"mass", [&](const Tok& p) {mass = p.rvec[0]; return true;}},
		{"totalmoi", [&](const Tok& p) { totalmoi = p.rvec; return true;}},
		{"moi", [&](const Tok& p) { moi = p.rvec; return true;}},
		{"cg", [&](const Tok& p) {return par_or_value(p,cgpar,cg); }},
		{"orient", [&](const Tok& p) { return parse_orient(p,orient,beams); }},
		{"node(s)?", [&](const Tok& p) { nodes = p.ivec; return true; }}
	});

	if (!unrec.empty())
		throw runtime_error(vastr("unrecognized specs: ",unrec));

	// create a Mass element for each node in "nodes" with "mass", "moi",
	// "cg", and "orient"

	trc.dprint(nodes.size()," nodes: ",nodes);
	if (nodes.empty())
		throw runtime_error("the nodes where conmasses are attached must be specified");

	// default moi: if only one given, apply to all 3
	if (moi.empty())
		moi = vector<double>{1.0,1.0,1.0};
	double last = moi[moi.size()-1];
	while(moi.size() < 3)
		moi.push_back(last);

	// XXX need retained's?
	if (ret.empty())
		throw runtime_error("retained freedoms have not been defined");
	trc.dprint("retained freedoms: ",ret);

	//  if totalmass specified instead of mass (per node), compute mass
	if (totalmass > 0.0) {
		mass = totalmass/nodes.size();
		flaps::info("mass per node = ",totalmass,"/",nodes.size()," = ",mass);
	}
	//  if totalmoi specified instead of moi (per node), compute moi
	if (!totalmoi.empty()) {
		moi = totalmoi;
		for(auto& mm : moi)
			mm /= nodes.size();
		flaps::info("moi per node = ",totalmoi,"/",nodes.size()," = ",moi);
	}

	// treat each node
	for (int node : nodes) {
		// compute the 1b row/col numbers (rc) where this element belongs
		vector<int> rc(6,0);
		for (int i=0; i<6; i++) {
			auto idx = std::find(ret.begin(), ret.end(), Freedom(node,i+1));
			if (idx != ret.end())
				rc[i] = idx - ret.begin() + 1; // 1b row/col
		}
		if (*std::max_element(rc.begin(),rc.end()) == 0) {
			//throw runtime_error(vastr("node ",node," is undefined"));
			flaps::warning("node ",node," has no retained freedoms");
			continue;
		}
		// create a new Mass & add it to masses
		Mass* me = new Mass(node,mass,moi,cgpar,cg,orient,rc);
		masses.push_back(me);
		trc.dprint("now have ",masses.size()," mass elements");
		trc.dprint(*me);
	}
	return true;
}

vector<double>
makeT(const vector<double>& orient) {
// given the orientation of the local x axis (xbar), and a point
// on the xbar-ybar plane, create a transformation T such that
//   vglobal = T' * vlocal
	Trace trc(2,"makeT");
	vector<double> rval(9,0.0);
	vector<double> xbar{orient};
	// normalize xbar
	double xn = blas_snrm2(3, &xbar[0], 1);
	if (xn == 0.0)
		throw runtime_error("zero-length orient");
	blas_scal(3, 1.0/xn, &xbar[0], 1);

	// svd xbar to get 2 vectors orthogonal to it
	vector<double> u(9,0.0), s(3,0.0), A{xbar};
	double vt;
	int info = lapack::dgesvd("a","a",3,1,&A[0],3,&s[0],&u[0],3,&vt,1);
	if (info != 0)
		throw runtime_error(vastr("makeT failed: svd error ",info));
	vector<double> vec1{u[3],u[4],u[5]};
	vector<double> vec2{u[6],u[7],u[8]};

	// tx is the matrix representation of the cross-product xbar X v
	vector<double> tx{0.0,xbar[2],-xbar[1],-xbar[2],0.0,xbar[0],xbar[1],-xbar[0],0.0};

	vector<double> txy(3,0.0);  // t X vec1
	blas_sgemv("n", 3,3,1.0,&tx[0],3,&vec1[0],1,0.0,&txy[0],1);

	// xbar,ybar,zbar are the rows of T
	blas_copy(3, &xbar[0], 1, &rval[0], 3);
	blas_copy(3, &vec1[0], 1, &rval[1], 3);
	blas_copy(3, &txy[0], 1, &rval[2], 3);
	return rval;
}

string
Fem::
mid(const string& name) {
	if (id().empty())
		return name;
	else
		return name + "." + id();
}

void
conmass(double mass, const vector<double>& moi, const vector<double>& cgv,
	vector<double>& elem) {
// Given a concentrated mass, moments, and the vector (in global)
// from the node where the mass is lumped to the cg, compute the (6,6)
// elemental mass matrix at the node

	assert(elem.size() == 36);
	for (int i=0; i<36; i++)
		elem[i] = 0.0;

	// diagonals: mass and moi
	for (int i=0; i<3; i++) {
		elem[IJ(i,i,6)] = mass;
		elem[IJ(i+3,i+3,6)] = moi[i];
	}
	// moments at the point
	if (!cgv.empty()) {
		assert(cgv.size() == 3);
		// rotation-translation coupling: the kinetic energy is
		// 1/2m(a X cgv)^2 where a is one of the coord axes.
		// each axis (rotational velocity) with the cg vector (cgv)
		for (int j=0; j<3; j++) {
			vector<double> prod(3,0.0);
			vector<double> axis(3,0.0);
			axis[j] = 1.0;
			xprod(&axis[0], &cgv[0], &prod[0]);
			blas_axpy(3, mass, &prod[0], 1, &elem[IJ(0,j+3,6)], 1);
			// lower-triangle
			for (int i=0; i<3; i++) {
				elem[IJ(j+3,i,6)] = elem[IJ(i,j+3,6)];
			}
		}
		// moments
		elem[IJ(3,3,6)] += mass*(cgv[1]*cgv[1] + cgv[2]*cgv[2]);
		elem[IJ(4,4,6)] += mass*(cgv[0]*cgv[0] + cgv[2]*cgv[2]);
		elem[IJ(5,5,6)] += mass*(cgv[0]*cgv[0] + cgv[1]*cgv[1]);
	}
}

bool
Fem::
dlm_f (const Tok& tok) {
// parse dlm options: dlm{ac=d,panels=d,...}
	Trace trc(1,"dlm_f");
	Dlm de;
	// parse the options with lambda handlers
	vector<Tok*> unrec = flaps::lexer(tok.lopt, {
		{"ac", [&](const Tok& p) {de.acs = p.rvec; return true;}},
		{"reflen", [&](const Tok& p) {de.reflen = p.rvec[0]; return true;}},
		{"panels", [&](const Tok& p) {de.panels = p.ivec; return true;}},
		{"semispan", [&](const Tok& p) {de.semi_span = p.rvec[0]; return true;}},
		{"ss", [&](const Tok& p) {de.ss = p.srhs; return true;}},
		{"chord", [&](const Tok& p) {de.chord = p.rvec[0]; return true;}},
		{"taper", [&](const Tok& p) {de.taper = p.rvec[0]; return true;}},
		{"sweep", [&](const Tok& p) {de.sweep = p.rvec[0]; return true;}},
		{"dihedral", [&](const Tok& p) {de.dihedral = p.rvec[0]; return true;}},
		{"rf", [&](const Tok& p) {de.rfs = p.rvec; return true;}},
		{"rsf", [&](const Tok& p) {de.rsfs = p.rvec; return true;}},
		{"mach", [&](const Tok& p) {de.machs = p.rvec; return true;}},
		{"nodes", [&](const Tok& p) {de.nodes = p.ivec; return true;}},
	});
	
	if (!unrec.empty())
		throw runtime_error(vastr(unrec.size()," unrecognized dlm specs: ", unrec));
	trc.dprint("got Dlm\n",de);

	// copy this Dlm to the Fem
	this->dlm = new Dlm(de);
	return true;
}

#ifdef NEVER // unused?
vector<double>
guyan(const vector<double>& elem) {
// Given a (12,12) stiffness element compute the Guyan transformation
// matrix to condense node 2 onto node1
//   T = |       I       |
//       | -k22^{-1} k21 |
// but only return the lower (6,6) portion:  -k22^{-1}k21
	Trace trc(1,"guyan");
	vector<double> k22(36, 0.0);
	vector<double> k21(36, 0.0);
	for (int i=0; i<6; i++) {
		for (int j=0; j<6; j++) {
			k22[IJ(i,j,6)] = elem[IJ(i+6,j+6,12)];
			k21[IJ(i,j,6)] = elem[IJ(i+6,j,12)];
		}
	}
	int n{6};
	int info = lapack::dposv("u",n,n,&k22[0],n,&k21[0],n);
	if (info != 0)
		throw runtime_error(vastr("dposv returned ",info));
	vector<double> rval(6*6, 0.0);
	for (int i=0; i<6; i++) {
		for (int j=0; j<6; j++) {
			rval[IJ(i,j,6)] = -k21[IJ(i,j,6)];
		}
	}
	trc.dprintm(6,6,6,&k22[0],"k22");
	trc.dprintm(6,6,6,&k21[0],"k21");
	trc.dprintm(6,6,6,&rval[0],"T");
	return rval;
}
#endif // NEVER // unused?

vector<Freedom>
Fem::
get_rotations(int node) {
// returns a vector of all rotational dof at "node"; if a dof is
// missing set it to Freedom(0,0) so the return vector is always
// 3 Freedoms
	Trace trc(1,"get_rotations");
	vector<Freedom> rval(3, Freedom());
	for (auto& fp : retained()) {
		if (fp.node() == node && fp.dof() > 3)
			rval[fp.dof()-4] = fp;
	}
	trc.dprint("returning ",rval);
	return rval;
}

void
Fem::
assemble() {
// Statically merge all substructures mass and stiffness matrices
// assemble all element matrices for the gyro, pgaf, and dlm,
// into single flaps matrices with mid "id_stif", "id_mass", etc
// Store them with tag '_nodal'
	Trace trc(1,"Fem::assemble");

	// assemble substructure stif matrices
	for (auto ss : substructures_)
		ss->assemble(substructures_);

	// collect dof from all the beams...
	for (auto bp : stif_elements)
		for (auto fp : bp->freedoms())
			this->retained_.add(fp);
	// ... and substructures
	for (auto ss : substructures_)
		for (auto fp : ss->retained())
			this->retained_.add(fp);

	size_t n = retained_.size();

	trc.dprint("gross matrices retained freedoms:\n",this->retained());
	cout << "Model nodes and freedoms\n" << this->retained() << endl;

	this->stif_ = new Matrix(mid("stif.nodal"), "beam fem", n, n, false);
	this->mass_ = new Matrix(mid("mass.nodal"), "fem mass", n, n, false);

	// XXX should we allow both substructures and a non-substructured Fem
	// OR if there are substructures make the non-ss a substructure?

	// merge substructure stif matrices into stif_nodal
	for (auto ss : substructures_) {
		blas::copy(ss->retained(), ss->retained(), 1.0, ss->stif(),
			this->retained(), this->retained(), 1.0, this->stif_->data());
	}

	// merge element stiffness matrices into the output stiffness matrix
	for (auto& ei : this->stif_elements)
		blas::copy(ei->freedoms(), ei->freedoms(), 1.0, ei->data(),
			this->retained(), this->retained(), 1.0, this->stif_->data());
	// store the nodal stif matrix
	this->stif_->store();
	trc.dprintm(n,n,n,this->stif_->data(),"nodal stif matrix");
	cout << "  Created nodal stiffness matrix \"" << *stif_ << "\"\n";

	// assemble the mass matrix
	//  - first merge substructure mass matrices into mass_nodal
	for (auto ss : substructures_) {
		blas::copy(ss->retained(), ss->retained(), 1.0, ss->mass(),
			this->retained(), this->retained(), 0.0, this->mass_->data());
	}
	//  ... then the fem model
	if (!this->mass_elements.empty()) {
		// if any of the mass elements has a cg parameter (cgpar),
		// create a CustomPz function
		bool needpz{false};
		for (auto& mi : this->mass_elements) {
			if (!mi->cgpar().empty()) {
				needpz = true;
				break;
			}
		}
		if (needpz) {
			// create a CustomPz function to evaluate the mass matrix...
			string entrypt;
			string path = make_massfcn(entrypt);
			// ... add it to custom.so ...
			Custom::add(path);
			// ... and add the CustomPz to the output matrix
			this->mass_->pz.push_back(new CustomPz(entrypt, n, n, 0));
		} else {
			// no cg parameter - assume Mass_elem::cg is the value, and
			// orient is the orientation; we don't have to worry about masses
			// being used multiple times like we do with substructures
			for (auto& mi : this->mass_elements) {
				// create the (6,6) elemental mass matrix...
				vector<double> elem(36, 0.0);
				vector<double> cgv{mi->orient()};
				for (auto& ci : cgv)
					ci *= mi->cg();
				conmass(mi->mass(), mi->moi(), cgv, elem);
				// ... and insert into the gross matrix
				double* m = this->mass_->elem();
				vector<int> rc = mi->rc();
				for (int j=0; j<6; j++) {
					if (rc[j] == 0) continue;
					for (int i=0; i<6; i++) {
						if (rc[i] == 0) continue;
						m[IJ(rc[i]-1,rc[j]-1,n)] = elem[IJ(i,j,6)];
					}
				}
			}
		}
	}
	this->mass_->store();
	trc.dprintm(n,n,n,this->mass_->data(),"nodal mass matrix");
	cout << "  Created nodal mass matrix \"" << *mass_ << "\"\n";

	// Gyro matrix XXX make this an assemble()?
	if (!this->gyro_elements.empty()) {
		bool needpz{false};
		this->gyro_ = new Matrix(mid("gyro.nodal"), "fem gyro", n, n, false);
		for (auto& ge : this->gyro_elements) {
			// set up the rc array: row/col numbers for each element
			ge.rc = vector<int>(3,0);
			int node = ge.node();
			for (int j=0; j<3; j++) {
				if (blas_snrm2(3, &ge.elem[3*j], 1) != 0.0) {
					auto idx = std::find(this->retained().begin(), this->retained().end(),
						Freedom(node,j+4));
					if (idx != this->retained().end())
						ge.rc[j] = (int)(idx-this->retained().begin()+1); // 1b row/col#
				}
			}
			if (!ge.amom_name().empty())
				needpz = true;
			double* el = this->gyro_->elem();
			int nr = this->gyro_->rsize();
			double amom{1.0};
			if (!needpz)
				amom = ge.amom_ad(gpset::get()).value();
			for (size_t j=0; j<ge.rc.size(); j++) {
				if (ge.rc[j] == 0) continue;
				for (size_t i=0; i<ge.rc.size(); i++) {
					if (ge.rc[i] == 0) continue;
					el[IJ(ge.rc[i]-1,ge.rc[j]-1,nr)] = amom*ge.elem[IJ(i,j,3)];
				}
			}
		}
		trc.dprintm(n,n,n,this->gyro_->data(),"gyro matrix");

		// create a CustomPz function to evaluate the gyro terms if
		// this any elements have amom_name
		if (needpz) {
			// create a function to evaluate the nodal gyro...
			string entrypt;
			string path = make_gyrofcn(entrypt,"nodal");
			// ... add it to custom.so ...
			Custom::add(path);
			// ... and add the CustomPz to the output gyro matrix
			this->gyro_->pz.push_back(new CustomPz(entrypt, n, n, 0));
		}
		this->gyro_->store();
		cout << "  Created nodal gyro matrix \"" << *gyro_ << "\"\n";
	}

	// create a dlm matrix
	if (this->dlm != nullptr)
		this->dlm->assemble(grid, retained());

	// propeller gaf (Pgaf) matrix, and also pgaffcn will be added
	// to a dlm matrix if one is created
	if (!this->pgaf_elements.empty()) {
		this->pgaf_ = new Matrix(mid("pgaf.nodal"), "propeller gaf", n, n, true);
		for (auto& gaf : this->pgaf_elements) {
			if (gaf.orient.empty()) {
				flaps::warning("defaulting pgaf orient to minus global x");
				gaf.orient = vector<double>{-1.0, 0.0, 0.0};
			}
			// orient is the local x axis (xbar) - compute ybar and zbar
			// to create a local to model transform: vglobal = T'*vlocal
			gaf.T = makeT(gaf.orient);
			trc.dprint("pgaf element:\n",gaf);
		}

		// create a cpp file which will evaluate the prop gaf...
		// XXX assuming pgaf_custom is only for modal
		string entrypt;
		string path = make_pgaffcn(entrypt, "nodal");
		// ... add it to custom.so ...
		Custom::add(path);
		// ... and include it in the pgaf parameterizations
		CustomPz* pz = new CustomPz(entrypt, n, n, 0);
		this->pgaf_->pz.push_back(pz);
		this->pgaf_->store();
		cout << "  Created nodal pgaf matrix \"" << *pgaf_ << "\"\n";
		// and the dlm matrix parameterizations
		if (this->dlm != nullptr) {
			this->dlm->gaf_nodal->pz.push_back(pz);
		}
	}

	// dlm gaf matrix
	if (this->dlm != nullptr) {
		this->dlm->gaf_nodal->store();
		cout << "  Created nodal gaf matrix \"" << *this->dlm->gaf_nodal << "\"\n";
	}

	// Create a modal transformation matrix if either Branch-modes
	// or free-vibration were requested
	if (bma != nullptr) {
		//!! assemble_branches();
		//!! if (cms != nullptr)
			//!! cms->reduce();
	} else if (this->freevib())
		this->modal(freevibration());
}

string
Fem::
make_massfcn(string& entrypt) {
// create a custom function to evaluate the mass matrix
// in a CustomPz

	// create the file on the cwd, return the full path
	string id{this->id()};
	entrypt = mid("massfcn");
	string rval{vastr(get_cwd(),"/",entrypt,".cpp")};
	ofstream file(rval, std::ios::out | std::ios::trunc);
	if (!file)
		throw runtime_error(vastr("cannot open ",rval));

	file << "#include \"fem.h\"\n";
	file << "int " << entrypt << " (pset& plt, int nr, int nc, vector<complex<Ad>>& ca) {\n";
	// create a vector of Mass's
	file << "   vector<Mass> masses {\n";
	string comma;
	for (auto me : mass_elements) {
		vector<double> orient = me->orient();
		if (!me->cgpar().empty() && (orient.size() != 3 || norm2(orient) == 0.0))
			throw runtime_error("the cg orientation was not specified");
		vector<double> moi = me->moi();
		vector<int> rc = me->rc();
		file << "      " << comma << "{" << me->node() << ", " << me->mass() << ", ";
		file << "{" << moi[0] << ", " << moi[1] << ", " << moi[2] << "},";
		file << "\"" << me->cgpar() << "\", ";
		file << me->cg() << ", ";
		file << "{" << orient[0] << ", " << orient[1] << ", " << orient[2] << "},";
		file << "{" << rc[0] << ", " << rc[1] << ", " << rc[2] << ", "
			  <<  rc[3] << ", " << rc[4] << ", " << rc[5] << "} }\n";
		comma = ",";
	}
	file << "   };\n";

	file << "	// create each (6,6) element mass matrix\n";
	file << "	for (auto& me : masses) {\n";
	file << "      vector<double> moi = me.moi();\n";
	file << "		// get the current value of the cg\n";
	file << "		Ad cg = me.cg();\n";
	file << "		if (!me.cgpar().empty())\n";
	file << "			cg = plt.parval(me.cgpar());\n";
	file << "		vector<Ad> elem(36, 0.0);\n";
	file << "		for (int i=0; i<3; i++) {\n";
	file << "			elem[IJ(i,i,6)] = me.mass();\n";
	file << "			elem[IJ(i+3,i+3,6)] = moi[i];\n";
	file << "		}\n";
	file << "		// rotation-translation coupling due to cg offset\n";
	file << "      vector<double> orient = me.orient();\n";
	file << "		for (int j=0; j<3; j++) {\n";
	file << "			vector<double> prod(3,0.0);\n";
	file << "			vector<double> axis(3,0.0);\n";
	file << "			axis[j] = 1.0;\n";
	file << "			xprod(&axis[0], &orient[0], &prod[0]);\n";
	file << "			for (int i=0; i<3; i++) {\n";
	file << "				elem[IJ(i,j+3,6)] += cg*me.mass()*prod[i];\n";
	file << "				elem[IJ(j+3,i,6)] = elem[IJ(i,j+3,6)];  // lower-tri\n";
	file << "			}\n";
	file << "		}\n";
	file << "		// moments\n";
	file << "		Ad t{cg*cg*me.mass()};\n";
	file << "		elem[IJ(3,3,6)] += t*(orient[1]*orient[1] + orient[2]*orient[2]);\n";
	file << "		elem[IJ(4,4,6)] += t*(orient[0]*orient[0] + orient[2]*orient[2]);\n";
	file << "		elem[IJ(5,5,6)] += t*(orient[0]*orient[0] + orient[1]*orient[1]);\n";
	file << "		// insert into ca\n";
	file << "		const vector<int>& rc = me.rc();\n";
	file << "		for (int j=0; j<6; j++) {\n";
	file << "			if (rc[j] == 0) continue;\n";
	file << "			for (int i=0; i<6; i++) {\n";
	file << "				if (rc[i] == 0) continue;\n";
	file << "				ca[rc[i]-1 + nr*(rc[j]-1)] += elem[i+j*6];\n";
	file << "			}\n";
	file << "		}\n";
	file << "	}\n";
	file << "   return 0;\n";
	file << "}\n";

	return rval;
}

string
Fem::
make_gyrofcn(string& entrypt, string const& prefix) {
// create a custom gyro function to evaluate the gyro matrix
// in a CustomPz. Create the file on the cwd, return both
// the full path and the entry point name "entrypt"
	Trace trc(1,"make_gyrofcn");
	if (prefix.empty())
		entrypt = mid("gyrofcn");
	else
		entrypt = mid(vastr(prefix, "gyrofcn"));
	string rval{vastr(get_cwd(),"/",entrypt,".cpp")};
	ofstream file(rval, std::ios::out | std::ios::trunc);
	if (!file)
		throw runtime_error(vastr("cannot open ",rval));

	file << "#include \"fem.h\"\n";

	// write a vector of Gyro_elem with initializer-list constructors
	file << "int " << entrypt << " (pset& plt, int nr, int nc, vector<complex<Ad>>& ca) {\n";
	file << "// multiply gyro elements by their amom\n";
	file << "   vector<Gyro_elem> gyros{";
	string comma;
	for (auto& ge : gyro_elements) {
		file << comma << "{";  // start ge
		comma = ",";
		file << ge.node() << ", \"" << ge.amom_name() << "\", " << ge.amom_v();
		file << ", ";
		// elem: 9 doubles
		string sep;
		file << "{";
		for (auto d : ge.elem) {
			file << sep << d;
			sep = ",";
		}
		file << "}";  // end elem
		// 1b row/column numbers in model matrix: error if not retained
		sep = "";
		file << comma << "{";  // start rc
		for (auto i : ge.rc) {
			file << sep << i;
			sep = ",";
		}
		file << "}";  // end rc
		// (3,nmodes) transform from rotations to gc
		sep = "";
		file << comma << "{";	// start modal
		for (auto i : ge.modal_) {
			file << sep << i;
			sep = ",";
		}
		file << "}";	// end modal
		file << "}\n";  // end ge
	}
	file << "};\n";
	file << "	// multiply each (3,3) gyro matrix by its amom\n";
	file << "	// if modal is (3,nr) => modal dof\n";
	file << "   // zero the output array\n";
	file << "	for (auto& ge : gyros) {\n";
	file << "      Ad amom = ge.amom_ad(plt);\n";
	file << "      if (ge.modal_.size() == (size_t)(3*nr)) {\n";
	file << "         // transform the (3,3) element to modal\n";
	file << "         vector<double> gm(nr*nr);\n";
	file << "         lapack::triprod(3,nr,&ge.modal_[0],&ge.elem[0],&gm[0]);\n";
	file << "         for (int i=0; i<nr; i++)\n";
	file << "            for (int j=0; j<nr; j++)\n";
	file << "               ca[IJ(i,j,nr)] += amom*gm[IJ(i,j,nr)];\n";
	file << "	   } else {\n";
	file << "         for (size_t i=0; i<ge.rc.size(); i++) {\n";
	file << "            if (ge.rc[i] == 0) continue;\n";
	file << "            for (size_t j=0; j<ge.rc.size(); j++) {\n";
	file << "               if (ge.rc[j] == 0) continue;\n";
	file << "               ca[IJ(ge.rc[i]-1,ge.rc[j]-1,nr)] = amom*ge.elem[IJ(i,j,3)];\n";
	file << "            }\n";
	file << "         }\n";
	file << "	   }\n";
	file << "	}\n";
	file << "   static int visit{0};\n";
	file << "   if (visit++ == 0)\n";
	file << "      Dprint::dprintm(nr,nc,nr,ca,\"gyro\");\n";
	file << "   return 0;\n";
	file << "}\n";
	return rval;
}

string
Fem::
make_pgaffcn(string& entrypt, string const& prefix) {
// create a custom function to evaluate the propeller
// unsteady aero matrix (pgaf) in a CustomPz

	// create the file on the cwd, return the full path
	if (prefix.empty())
		entrypt = mid("pgaffcn");
	else
		entrypt = mid(vastr(prefix,"pgaffcn"));
	string rval{vastr(get_cwd(),"/",entrypt,".cpp")};
	ofstream file(rval, std::ios::out | std::ios::trunc);
	if (!file)
		throw runtime_error(vastr("cannot open ",rval));

	file << "#include \"fem.h\"\n";

	file << "int " << entrypt << " (pset& plt, int nr, int nc, vector<complex<Ad>>& ca) {\n";
	file << "// unsteady aero for spinning props\n";
	file << "   Trace trc(1,\"" << entrypt << "\");\n";
	file << "   vector<Pgaf_elem> props{\n";
	string elsep{"   "};
	for (auto pg : pgaf_elements) {
		file << elsep << "{";
		elsep = ",";
		// node number, radius & value, lbar & value, beta & value members
		string membersep;
		file << membersep << pg.node();
		membersep = ",";
		file << membersep << "\"" << pg.radius_s << "\", " << pg.radius;
		file << membersep << "\"" << pg.lbar_s << "\", " << pg.lbar;
		file << membersep << "\"" << pg.beta_s << "\", " << pg.beta;
		// 1b row/column numbers for (rx,ry,rz) in model matrix
		file << membersep << "{";
		string rcsep;
		for (int i=4; i<7; i++) {
			auto idx = std::find(this->retained().begin(), this->retained().end(),
				Freedom(pg.node(),i));
			int rc{0};
			if (idx != this->retained().end())
				rc = (int)(idx-this->retained().begin()+1); // 1b row/col#
			file << rcsep << rc;
			rcsep = ",";
		}
		file << "}";
		// transformation to model (3,3) or modal (3,nmodes)
		file << membersep << "{";
		string tsep;
		for (auto i : pg.T) {
			file << tsep << i;
			tsep = ",";
		}
		file << "}";
		file << "}\n";
	}
file << "};\n";
file << "   Ad sigma = plt.parval(\"sigma\");\n";
//!! file << "   Ad freq = plt.parval(\"freq\");\n";
file << "   Ad vtas = plt.parval(\"vtas\", true);\n";
file << "   Ad rf = plt.parval(\"rf\", true);\n";
file << "   Ad rsf = plt.parval(\"rsf\", true);\n";
file << "   // watch out for zero vtas: return zero ca\n";
file << "   double eps{0.1};\n";
file << "   if (vtas.value() <= eps)\n";
file << "      return 0;\n";
file << "	// compute the (3,3) gaf matrix for each prop\n";
file << "	for (auto& pi : props) {\n";
file << "		Ad beta = pi.beta_v(plt);\n";
file << "		Ad cmq, czr, czpsi, cztheta, cmpsi;\n";
file << "		whirl_deriv(beta, cmq, czr, czpsi, cztheta, cmpsi);\n";
file << "		Ad lbar = pi.lbar_v(plt);\n";
file << "		Ad radius = pi.radius_v(plt);\n";
file << "      complex<Ad> sbar(rsf, rf);\n";
file << "		complex<Ad> t = flaps::pi*radius*radius*radius*(1.0 - lbar*sbar);\n";
file << "		// (3,3) prop gaf matrix with the shaft along the local x axis\n";
file << "		// so the only non-zeros are ry, rz\n";
file << "		vector<complex<Ad>> el(3*3);\n";
file << "		el[IJ(1,1,3)] = t*(2.0*cmq*sbar - lbar*cztheta);\n";
file << "		el[IJ(1,2,3)] = t*(2.0*cmpsi - lbar*(czpsi + czr*sbar));\n";
file << "		el[IJ(2,1,3)] = -el[IJ(1,2,3)];\n";
file << "		el[IJ(2,2,3)] = el[IJ(1,1,3)];\n";
file << "		// transform to modal ...\n";
file << "      if (pi.T.size() == (size_t)(3*nr)) {\n";
file << "         vector<complex<Ad>> telt(nr*nr);\n";
file << "         lapack::triprod(3,nr,&pi.T[0],&el[0],&telt[0]);\n";
file << "		   // insert this (nr,nr) into the return matrix ca\n";
file << "         for (int i=0; i<nr; i++)\n";
file << "            for (int j=0; j<nr; j++)\n";
file << "               ca[IJ(i,j,nr)] += telt[IJ(i,j,nr)];\n";
file << "      // ... or model coordinates\n";
file << "      } else if (pi.T.size() == 9) {\n";
file << "         vector<complex<Ad>> telt(9);\n";
file << "         lapack::triprod(3,3,&pi.T[0],&el[0],&telt[0]);\n";
file << "		   // insert this (3,3) into the return matrix ca\n";
file << "         for (int i=0; i<3; i++) {\n";
file << "            if (pi.rc[i] == 0) continue;\n";
file << "            for (int j=0; j<3; j++) {\n";
file << "               if (pi.rc[j] == 0) continue;\n";
file << "               ca[IJ(pi.rc[i]-1,pi.rc[j]-1,nr)] += telt[IJ(i,j,3)];\n";
file << "            }\n";
file << "         }\n";
file << "      } else {\n";
file << "         throw runtime_error(vastr(\"T is \",pi.T.size(),\n";
file << "            \", should be \",3*nr,\" or \",9));\n";
file << "	   }\n";
file << "	}\n";
file << "   static int visit{0};\n";
file << "   if (visit++ == 0)\n";
file << "      Dprint::dprintm(nr,nc,nr,ca,\"pgaf\");\n";
file << "   return 0;\n";
file << "}\n";
	return rval;
}

/*------------------------------------------------------------------
 * Create a gyro matrix comprising one or more spinning rigid bodies
 *
 *  gyro {
 *		modes=modes, modal
 *		name=gyro, o=gyro
 *    node = 212500 {
 *			inertia=1.0,
 *			orient=(1,0,0),
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
 *       node = 212500 {inertia=1, orient=(1,0,0), }..
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

/*------------------------------------------------------------------
 * scale factors
 * Reference: The International System of Units, Physical Constants
 *  and Conversion Factors, E.A. Mechtly, NASA SP-7012,
 *  Second Revision, 1973
 *------------------------------------------------------------------*/
constexpr double Massconv{32.174048556*12.0}; // lbm-in/lbf-s^2

bool
Fem::
gyro_f (const Tok& p) {
// gyro{node=n{amom=string|double, orient=(x,y,z)}}
// gyro{node=n{amom=string|double, orient=beam_id(x,y,z)}}

	// each call to gyro_node_f creates a gyro element in Model.gyro_elements
	vector<Tok*> unrec = flaps::lexer(p.lopt, {
		{"node", [&](const Tok& p) {return this->gyro_node_f(p);}},
	});

	if (!unrec.empty())
		flaps::warning(vastr(unrec.size()," specs were not recognized: ",unrec));

	return true;
}

bool
Fem::
gyro_node_f (const Tok& p) {
/*
 * Given a spec like
 *    node = node_number{ amom=string|double, orient=(x,y,z) }
 * The general form of this option is:
 *    node = node_number{ amom=string|double, orient=(x,y,z) }
 * for example:
 *    node = 21{amom=237.0, orient=(1,0,0)}
 *
 *	The specs for each node are:
 *    amom     Angular momentum in kg-m^2/s (N-m-s^2) or lbf-in-sec,
 *             or the name of an angular momentum parameter
 *                default: 1 kg-m^2/s
 *
 *    orient   x, y, and z (model) coordinates of a point on the spin axis,
 *             relative to the node and in the direction of the spin vector.
 *             This defines the orientation and direction of spin.
 *             For example
 *                orient=(-1,0,0)
 *             defines a spin vector in the negative x direction,
 *             so by the right-hand rule the spin direction would
 *             be clockwise when viewed down the x-axis.
 *             It is not necessary to specify trailing zero coordinates,
 *             so this example is equivalent to
 *                orient=-1
 *             default: (0, 1, 0)
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
 */
	Trace trc(1,"gyro_node_f");
	
	static double amom{1.0};
	static string amom_s;
	static vector<double> orient{1.0, 0.0, 0.0};

	// the node number is the rhs
	if (p.ivec.empty())
		throw runtime_error("a node number must be specified");
	int node = p.ivec[0];

	// Parse the node options:
	//   (amom=1, orient=(x,y,z)
	//   (amom=name, orient=(x,y,z)
	// Preference handlers (lambdas)
	if (!p.ropt.empty()) {
		vector<Tok*> unrec = flaps::lexer(p.ropt, {
			{"amom", [&](const Tok& p) { return par_or_value(p, amom_s, amom); }},
			{"orient", [&](const Tok& p) { return parse_orient(p,orient,this->stif_elements);}},
		});
		if (!unrec.empty())
			throw runtime_error(vastr("unrecognized gyro node option: ",*unrec[0]));
	}

	// sanity checks
	if (orient.empty())
		throw runtime_error("at least one orientation component must be included");
	if (orient.size() < 3)
		for (size_t i=orient.size(); i<3; i++)
			orient.push_back(0.0);

	// create the 3,3 gyro element: the angular momentum vector h is
	// the orientation vector normalized to 1 and scaled by the angular
	// momentum, unless the element has an amom parameter in which case it is 1
	vector<double> h{orient};
	double length = blas_snrm2(3, &orient[0], 1);
	if (length + 1.0 == 1.0)
		throw runtime_error(vastr("the orientation point for node ",
					node, " is zero"));

	double scale{length};

	amom = 1.0;

	for (int j=0; j<3; j++)
		h[j] *= scale;

	// Form the 3,3 matrix of gyro terms for the 3 rotational dof.
	// Because the orient vector (hence the momentum vector) is in global
	// coordinates, the gyro element matrix will also be in global.
	vector<double> el(9, 0.0);
	el[0] = 0.0;			// 1,1
	el[1] = -h[2];		// 2,1
	el[2] = h[1];		// 3,1
	el[3] = h[2];		// 1,2
	el[4] = 0.0;			// 2,2
	el[5] = -h[0];		// 3,2
	el[6] = -h[1];		// 1,3
	el[7] = h[0];		// 2,3
	el[8] = 0.0;			// 3,3

	// the freedoms for this gyro element are rx,ry,rz in model coordinates
	// but depending on the orientation one may not have any gyro
	// forces on it
	vector<int> rc(3,0);
#ifdef NEVER // move to Fem::assemble - no retained yet
	for (int j=0; j<3; j++) {
		if (blas_snrm2(3, &el[3*j], 1) != 0.0) {
			auto idx = std::find(this->retained().begin(), this->retained().end(),
				Freedom(node,j+4));
			if (idx != this->retained().end())
				rc[j] = (int)(idx-this->retained().begin()+1); // 1b row/col#
		}
	}
#endif // NEVER // move to Fem::assemble - no retained yet

	// Create the gyro element without modal transform
	vector<double> modal;
	Gyro_elem np(node, amom_s, amom, el, rc, modal);
	trc.dprint("got gyro element:\n",np);
	this->gyro_elements.push_back(np);
	return true;
}

bool
Fem::
pgaf_f (const Tok& p) {
// pgaf{node=n{ node options }, custom=string }
// node options:
//   radius=string|double, beta=string|double,
//   lbar=string|double, orient=(x,y,z), ydir=(x,y,z)
// all specs are optional except node number

	// each call to pgaf_node_f creates a pgaf element in Model.pgaf_elements
	vector<Tok*> unrec = flaps::lexer(p.lopt, {
		{"node", [&](const Tok& p) {return this->pgaf_node_f(p);}},
		{"custom",[&](const Tok& p) {this->pgaf_custom = p.srhs; return true;}}
	});

	if (!unrec.empty())
		flaps::warning(vastr(unrec.size()," specs were not recognized: ",unrec));

	return true;
}

bool
Fem::
pgaf_node_f (const Tok& p) {
// create a pgaf element from user specs
	Trace trc(1,"pgaf_node_f");
	Pgaf_elem np;
	// radius, beta, lbar, and custom are the only options
	// radius, beta, and lbar are parameter names, static so "persistent"
	static string radius_p;
	static double radius_v;
	static string beta_p;
	static double beta_v;
	static string lbar_p;
	static double lbar_v;
	static vector<double> orient{0.0,1.0,0.0};	// xbar default: -x
	static vector<double> ydir{-1.0,0.0,0.0};		// point in the xbar-ybar plane
																// default: -x
	// the node number is the rhs, the only required spec
	np.node_ = p.ivec[0];

	// Parse the node options:
	//   radius=1, beta=34, lbar=1, orient=(x,y,z)
	// Preference handlers (lambdas)
	if (!p.ropt.empty()) {
		vector<Tok*> unrec = flaps::lexer(p.ropt, {
			{"radius", [&](const Tok& p) {
				return par_or_value(p,radius_p,radius_v);}},
			{"beta", [&](const Tok& p) {
				return par_or_value(p,beta_p,beta_v);}},
			{"lbar", [&](const Tok& p) {
				return par_or_value(p,lbar_p,lbar_v);}},
			{"orient", [&](const Tok& p) { return parse_orient(p,orient,this->stif_elements); }}
		});

		if (!unrec.empty())
			throw runtime_error(vastr("unrecognized pgaf node option: ",*unrec[0]));
	}

	// set/create defaults
	np.radius_s = radius_p;
	if (np.radius_s.empty()) {
		if (radius_v != 0.0)
			np.radius = radius_v;
		else
			np.radius = 1.0;		// default radius
	} else {
		Par* pp = gpset::find(np.radius_s);
		np.radius = pp->value();
	}
	np.beta_s = beta_p;
	if (np.beta_s.empty()) {
		if (beta_v != 0.0)
			np.beta = beta_v;
		else
			np.beta = 34.0;		// default beta
	} else {
		Par* pp = gpset::find(np.beta_s);
		np.beta = pp->value();
	}
	np.lbar_s = lbar_p;
	if (np.lbar_s.empty()) {
		if (lbar_v != 0.0)
			np.lbar = lbar_v;
		else
			np.lbar = 1.0;		// default lbar
	} else {
		Par* pp = gpset::find(np.lbar_s);
		np.lbar = pp->value();
	}

	// put it on the vector of pgaf elements
	this->pgaf_elements.push_back(np);

	trc.dprint("created pgaf element:\n",np);

	return true;
}


// operator<< for most classes
ostream&
operator<<(ostream& s, const Node& t) {
	s << t.number_ << ": (" << t.coord_[0] << ", " << t.coord_[1]
		<< ", " << t.coord_[2] << ")\n";
	return s;
}

ostream&
operator<<(ostream& s, const Grid& t) {
	for (auto& ni : t.nodes_)
		s << ni << endl;
	return s;
}

ostream&
operator<<(ostream& s, const Freedom& t) {
	s << t.node() << '/' << t.dof();
	return s;
}

ostream&
operator<<(ostream& s, const vector<Freedom>& t) {
	if (t.empty())
		return s;
	int ord{1};  // ordinal number
	int last{t[0].node()};
	for (auto ti : t) {
		if (ti.node() != last) {
			last = ti.node();
			s << endl;
		}
		s << "  " << std::setw(5) << ord++ << ")   "
			<< ti.node() << '/' << ti.dof();
	}
	return s;
}

ostream&
operator<<(ostream& s, const Beam& t) {
	s << "beam";
	if (!t.id().empty())
		s << " " << t.id();
	s << " nodes (" <<  t.node1() << ", " << t.node2() << ")\n";
	s << "freedoms: ";
	string sep;
	for (auto fi : t.freedoms()) {
		s << sep << fi.node() << '/' << fi.dof();
		sep = " ";
	}
	s << endl;
	s << "transform from local to global:\n";
	const vector<double>& tlg = t.transform();
	s << "   " << tlg[0] << "  " << tlg[3] << "  " << tlg[6] << endl;
	s << "   " << tlg[1] << "  " << tlg[4] << "  " << tlg[7] << endl;
	s << "   " << tlg[2] << "  " << tlg[5] << "  " << tlg[8] << endl;
	return s;
}

ostream&
operator<<(ostream& s, const vector<Beam*>& t) {
	for (auto ni : t)
		s << *ni << endl;
	return s;
}
ostream&
operator<<(ostream& s, const Mass& t) {
	s << "node " << t.node() << " mass " << t.mass() <<  " (Kg)";
	vector<double> moi = t.moi();
	s << " moment-of-inertia: " << moi[0] << ", " << moi[1] << ", "
	  << moi[2] << " (Kg-m)\n";
	s << "cg: " << t.cg() << endl;
	vector<double> orient = t.orient();
	s << "orient: " << orient[0] << ", " << orient[1]
		<< ", " << orient[2] << endl;
	return s;
}

ostream&
operator<<(ostream& s, const Gyro_elem& ge) {

	s << "  node ............ " <<  ge.node() << " rows/cols: "
		<< ge.rc[0] << ", " << ge.rc[1] << ", " << ge.rc[2] << endl;
	string nm = ge.amom_name();
	if (nm.empty())
		s << "amom: " << ge.amom_v() << endl;
	else
		s << "amom parameter: " << ge.amom_name() << endl;
	s << "  element gyro matrix:\n";
	int i;
	if (!ge.elem.empty())
		for (i=0; i<3; i++)
			s << "      " << ge.elem[i] << "  " << ge.elem[i+3] << "  "
				<< ge.elem[i+6] << endl;
	return s;
}

ostream&
operator<<(ostream& s, const Pgaf_elem& pe) {

	s << "  node ............ " <<  pe.node() << endl;
	s << "  radius .......... " << pe.radius << " " << pe.radius_s << endl;
	s << "  lbar ............ " << pe.lbar << " " << pe.lbar_s << endl;
	s << "  beta ............ " << pe.beta << " " << pe.beta_s << endl;
	if (!pe.orient.empty())
		s << "  spin orientation relative to the node: " << pe.orient[0]
			<< ", " << pe.orient[1] << ", " << pe.orient[2] << endl;
	if (!pe.T.empty()) {
		s << "  transformation from local to model coord:\n";
		for (int i=0; i<3; i++)
			s << "      " << pe.T[i] << "  " << pe.T[i+3] << "  "
				<< pe.T[i+6] << endl;
	}
	if (!pe.custom.empty())
		s << "custom evaluation function: " << pe.custom << endl;
	return s;
}

