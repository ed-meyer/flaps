//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

#ifndef FEM_H
#define FEM_H 1

#include <complex>
#include <stdexcept>
#include <string>
#include <vector>

#include "exim.h"
#include "lapack.h"
#include "lexer.h"
#include "matrix.h"
#include "trace.h"

class Fem;

class Node {
	int number_;
	std::vector<double> coord_; // in the global coord system
public:

	Node(int n, std::vector<double> c) : number_(n), coord_(c) {}

	int number() const { return number_; }
	const std::vector<double> coord() const { return coord_; }
	double x() const { return coord_[0]; }
	double y() const { return coord_[1]; }
	double z() const { return coord_[2]; }

	friend std::ostream& operator<<(std::ostream& s, const Node& t);
};

// a Grid is a collection of Nodes
class Grid  {
	std::vector<Node> nodes_;
	double norm_{0.0};

public:
	Grid() {}

	void add(const Node& node) { nodes_.push_back(node); } // replace w/add2vector

	const std::vector<Node>& nodes() const { return nodes_; }

	const Node* find(int number) const {
		for (auto& ni : nodes_)
			if (ni.number() == number) return &ni;
		throw std::runtime_error(vastr("node ",number," has not been defined"));
	}

	double norm() {
		if (norm_ == 0.0)
			for (auto& ni : nodes_)
				norm_ = std::max(norm_, blas_snrm2(3,&ni.coord()[0],1));
		return norm_;
	}

	int nnodes() { return nodes_.size(); }

	friend std::ostream& operator<<(std::ostream& s, const Grid& t);
}; // Grid

// a Freedom is a node/dof combination
class Freedom {
	int node_;
	int dof_;  // 1-6

public:
	Freedom() : node_(0), dof_(0) {}
	Freedom(int node, int dof) : node_(node), dof_(dof) {}

	int node() const { return node_; }
	int dof() const { return dof_; }

	//!! double key() const { return (double(node_) + double(dof_)/10.0); }

	// operator== necessary for std::find
	bool operator==(const Freedom& t) const { return (t.node() == node_ && t.dof() == dof_); }
};

class Freev : public std::vector<Freedom> {
public:
	Freev() {}
	void add(const Freedom& f);	// add if not already here
};

// Material properties
struct Material {
	std::vector<double> x{0.0,0.0,0.0};
	double eix{0.0};
	double eiz{0.0};
	double gj{0.0};
	double ea{0.0};
};

std::ostream& operator<<(std::ostream& s, const Material& t);

struct Nodedof {
	int node;
	std::vector<int> dof;
	Nodedof() noexcept = default;
	Nodedof(int n, std::vector<int> const& d) : node(n), dof(d) {}
	Nodedof(Nodedof const&) noexcept = default;
	Nodedof& operator=(Nodedof const& rhs) = default;
};

// definition of a beam element
class Beam {
	std::string id_;
	int node1_;
	int node2_;
	std::vector<Freedom> freedoms_;	// nf retained freedoms
	std::vector<double> tlg_;			// (3,3) transform: global = tlg'*local
	Material matl;
	std::vector<double> data_;			// (nf,nf) beam element in global
public:
	std::vector<double> branchT;
	Beam() = default;
#ifdef NEVER // new ctor
	Beam(std::string ident, int n1, int n2, const std::vector<Freedom>& free,
		const std::vector<double>& trans, const std::vector<double>& elem);
		//!! : id_(ident), node1_(n1), node2_(n2), freedoms_(free),
			//!! tlg_(trans), data_(elem);
#else // NEVER // new ctor
	Beam (Nodedof const& n1, Nodedof const& n2, Material const& mtl, Grid const& grid);
#endif // NEVER // new ctor

	const std::string& id() const { return id_; }
	int node1() const { return node1_; }
	int node2() const { return node2_; }
	const std::vector<Freedom>& freedoms() const { return freedoms_; }
	const std::vector<double>& transform() const { return tlg_; }
	const std::vector<double>& data() const { return data_; }
};

std::ostream& operator<<(std::ostream& s, const Beam& t);
std::ostream& operator<<(std::ostream& s, const std::vector<Beam*>& t);

// concentrated masses have a mass, moi, orientation, cg
// Each Mass can only be used ONCE (see "used" member)
class Mass {
	int node_;
	double mass_;
	std::vector<double> moi_;
	std::string cgpar_;
	double cg_;
	std::vector<double> orient_;
	std::vector<int> rc_;		// 6 1b row/col numbers: 0 => not retained
	bool used_{false};
public:
	Mass(int node, double mass, const std::vector<double>& moi,
			const std::string& cgpar, double cg, const std::vector<double>& orient,
			const std::vector<int> rc) :
				node_(node), mass_(mass), moi_(moi), cgpar_(cgpar),
				cg_(cg), orient_(orient), rc_(rc) {}

	int node() const { return node_; }
	double mass() const { return mass_; }
	std::vector<double> moi() const { return moi_; }
	std::string cgpar() const { return cgpar_; }
	double cg() const { return cg_; }
	std::vector<double> orient() const { return orient_; }
	std::vector<int> rc() const { return rc_; }
	bool used() { return used_; }
};

class Gyro_elem {
	int node_{0};              // node number where the gyro acts
	std::string amom_s;			// amom parameter name
	double amom_{0};				// Angular Moment in kg-m^2/s (N-m-s^2)
public:
	std::vector<double> elem;	// 3x3 gyro element
	std::vector<int> rc;			// 1b row/column numbers in output matrix
	std::vector<double> modal_;	// transform from 3 rotations to modal gc

	// constructors
	Gyro_elem() {}
	// Gyro_elem(int node, const std::string& amomnm, double amom, 
	// 	const std::vector<double>& el, const std::vector<int>& r) :
	// 	node_(node), amom_s(amomnm), amom_(amom), elem(el), rc(r) {}
	Gyro_elem(int node, const std::string& amomnm, double amom,
		const std::vector<double>& el, const std::vector<int>& r,
		std::vector<double> modal={}) :
		node_(node), amom_s(amomnm), amom_(amom), elem(el), rc(r), modal_(modal) {}
	// default copy constructor, assignment operator ok

	int node() const { return node_; }
	std::vector<double>& modal() { return modal_; }

	// angular momentum parameter name
	std::string amom_name() const { return amom_s; }
	double amom_v() const { return amom_; }
	// return the current angular momentum
	Ad amom_ad(pset& plt) {
		if (amom_s.empty())
			return amom_;
		else
			return plt.parval(amom_s);
	}
};

// propeller gaf matrix
class Pgaf_elem {
public:
	int node_{0};              // node number where the pgaf acts
	std::string radius_s;		// parameter name
	double radius;					// value
	std::string lbar_s;			// l/radius parameter name
	double lbar;					// l/radius
	std::string beta_s;        // 3/4 blade angle parameter name
	double beta;					// 3/4 blade angle
	std::vector<int> rc;			// 1b row/col number in gross matrix XXX make 0b
	std::vector<double> T;		// transformation from local to model
	// these are only needed for vzpgaf?
	std::vector<double> orient;    // spin orientation (x,y,z) relative to node
	std::string custom;			// file or entry point for evaluating pgaf
	
	Pgaf_elem() {}
	Pgaf_elem(
			int node,
			const std::string& r, double rv,
			const std::string& l, double lv,
			const std::string& b, double bv,
			const std::vector<int> rows,
			const std::vector<double> t) :
		node_(node),
		radius_s(r), radius(rv),
		lbar_s(l), lbar(lv),
		beta_s(b), beta(bv),
		rc(rows), T(t) {};

	int node() const { return node_; }

	Ad radius_v(pset& plt) {
		if (radius_s.empty())
			return radius;
		else
			return plt.parval(radius_s);
	}

	Ad lbar_v(pset& plt) {
		if (lbar_s.empty())
			return lbar;
		else
			return plt.parval(lbar_s);
	}

	Ad beta_v(pset& plt) {
		if (beta_s.empty())
			return beta;
		else
			return plt.parval(beta_s);
	}
}; // Pgaf_elem

class Dlm {
public:
	std::string ss;				// substructure id, use fem if empty
	double reflen{0.0};					// reference length aka "b"
	std::vector<int> panels; // # chordwise, spanwise (m,n)
	double semi_span{0.0};
	double chord{0.0};
	double taper{1.0};
	double sweep{0.0};
	double dihedral{0.0};
	std::vector<int> nodes;  // fem grid node numbers for paneling
	std::vector<double> X;  // 3,m+1,n+1
	std::vector<int> conn;  // (4,m,n) 0b indexes into X
	// these may have multiple values; if so they are used in the interpolation
	std::vector<double> rfs;		// reduced frequencies
	std::vector<double> rsfs;		// reduced stability factors
	std::vector<double> machs;		// mach numbers
	std::vector<double> acs;			// x-coord of the aerodynamic center
	// matrices in nodal
	std::vector<Matrix*> nodal_gafs; // dlm gaf matrices at rf's
	Matrix* gaf_nodal{nullptr};   // interpolated DLM matrices in nodal
											//
	// transform from structural grid to aero grid: (nr,nc)
	// where nr = 3*(m+1)*(n+1), nc is the number of gc
	std::vector<double> T;
	// function to create T
	std::vector<double> transform(int m, int n, const std::vector<double>& X,
			const std::vector<Node>& grid, const std::vector<Freedom>& freedoms);

	// Create the output DLM matrices
	void assemble(const Grid& grid, const Freev& freev);

	// add aero panels to an amvz grid
	void vzaero(std::vector<int>& node_numbers, std::vector<double>& coords,
		std::vector<int>& segments, std::vector<double>& gct_nodal, int nc);

	Dlm() {}
};

class SS {
	std::string id_;
	Material matl;
	std::vector<Beam*> beams_;
	std::vector<Mass*> masses_;
	Freev retained_;	// all freedoms
	Freev attach_;		// attachments
	Freev interior_;
	std::vector<double> mass_;	// order retained.size() assembled matrix
	std::vector<double> stif_;	// order retained.size() assembled matrix
	std::vector<double> Gia_;	// Guyan transform interior-attachment
	std::vector<int> proj_;		// axis # (1-3) to project eigenvector
public:
	SS(const std::string& id, const std::vector<Beam*>& beams,
			const std::vector<Mass*>& masses, Freev& freev, const std::vector<int>& project);

	std::vector<double> kii;	// XXX computed in gia() used where?

	Grid grid;	// XXX unused?

	// assemble this substructure's mass and stiffness
	// allss is needed to determine attachment dof
	void assemble(const std::vector<SS*> allss);

	const Freev& retained() const { return retained_; }
	const Freev& attach() const { return attach_; }
	const Freev& interior() const { return interior_; }

	std::string id() const { return id_; }

	const std::vector<Beam*> beams() const { return beams_; }
	std::vector<double>& mass() { return mass_; }
	std::vector<double>& stif() { return stif_; }

	const std::vector<double>& Gia() { return Gia_; }

	// compute G_{ia} = -K_{ii}^{-1}K_{ia}
	std::vector<double> gia();

	// condense the interior dof to the attachment: maa, kaa
	std::vector<double> lump(std::vector<double>& kaa);

	// project each column of V onto the corresponding axis in proj_
	const std::vector<int>& proj() const { return proj_; }
	void projection(std::vector<double>& v);
}; // class SS


#ifdef NEVER // deprecate
class Branch : public SS {
	int nmodes_;
	std::vector<int> proj_;	// axis # (1-3) to project eigenvector
	//!! std::vector<Freedom> attach_;
	//!! std::vector<Freedom> interior_;
	std::vector<double> Gia_;
	std::vector<double> gia();	// compute Gia_
	std::vector<double> lumped_mass_;	// interior mass on attach
	std::vector<double> lumped_stif_;	// interior stif on attach
	std::vector<double> guyan();	// compute lumped_mass_
	std::vector<double> lump(std::vector<double>&);	// compute lumped_mass_, lumped_stif_
public:
	Grid* grid;

	// one-and-only constructor - no body
	Branch(const std::string& id, const std::vector<Beam*>& be,
		const std::vector<Mass*>& me,
		const std::vector<Freedom>& att,
		const std::vector<Freedom>& interior,
		int nm, const std::vector<int> project, Grid* g) :
		SS(id,be,me), nmodes_(nm), proj_(project),
		attach_(att), interior_(interior), grid(g) {}

	std::string id() const { return SS::id(); }
	const std::vector<Beam*> beams() const { return SS::beams(); }
	//!! const std::vector<Freedom>& attach() const { return attach_; }
	//!! const std::vector<Freedom>& interior() const { return interior_; }
	int nmodes() const { return nmodes_; }
	const std::vector<int>& proj() const { return proj_; }
	void assemble();
	const std::vector<double>& Gia() { return Gia_; }
	const std::vector<double>& lumped_stif() { return lumped_stif_; }
	const std::vector<double>& lumped_mass() { return lumped_mass_; }
	void projection(std::vector<double>& v);
	std::vector<double> kii;		// interior stif
}; // Branch

class Root : public SS {
	int nmodes_{0};
public:

	Root(const std::string& id, const std::vector<Beam*> beams,
		const std::vector<Mass*> masses, int nmodes) :
		SS(id,beams,masses), nmodes_(nmodes) {}

	std::string id() const { return SS::id(); }
	const std::vector<Beam*> beams() const { return SS::beams(); }
	int nmodes() const { return nmodes_; }
	void assemble();
};
#else // NEVER // deprecate
class Comp {
public:
	// Comp(const std::string& ssid, const std::string& nmstr)

	SS* ss{nullptr};
	int nmodes{0};
};
std::ostream& operator<<(std::ostream& s, const Comp& t);

class BMA {
public:
	Comp root;
	std::vector<Comp> branches;
	std::vector<double> xform;		// transform to BMA gc
	std::vector<double> transform(const Freev& gross);
};
std::ostream& operator<<(std::ostream& s, const BMA& t);

class CMS {
public:
	std::vector<Comp> ss;
	std::vector<double> xform;		// transform to CMS gc
	std::vector<double> transform(const Freev& gross);
};
#endif // NEVER // deprecate

class Fem {
	std::string id_;  // no default
	// retained dof are determined in beam_f as beam elements are added
	Freev retained_;
	// output matrices
	Matrix* stif_{nullptr};  // assembled model stiffness matrix
	Matrix* mass_{nullptr};  // assembled model mass matrix
	Matrix* gyro_{nullptr};  // assembled model gyro matrix
	Matrix* pgaf_{nullptr};  // assembled propeller gaf matrix
	std::string pgaf_custom;	// custom eval function for pgaf

	bool freevib_{false};
	int nmodes_{0};	// free-vibration modal reduction
	std::vector<double> modal_;	// transform from modal gc to nodal
	std::vector<int> nodes_;		// node numbers: nnodes
	std::vector<double> coords_;	// coordinates in same order: (3,nnodes)
	std::vector<int> conn_;			// connectivity in penlift format
	std::vector<double> gct_nodal_;	// transform from nodal to (tx,ty,tz)
	std::vector<double> gct_;			// transform from modal to (tx,ty,tz)
	size_t nvz_{0};	// number of visualization freedoms: rows of gct_nodal & gct_

	// substructuring stuff
	std::vector<SS*> substructures_;
public:
	BMA* bma{nullptr};
	CMS* cms{nullptr};

	bool si_units{true};
	Grid grid;
	Material matl;
	std::vector<Beam*> stif_elements;
	std::vector<Mass*> mass_elements;
	std::vector<Gyro_elem> gyro_elements;
	std::vector<Pgaf_elem> pgaf_elements;
	Dlm* dlm{nullptr};
	std::vector<double> bmt;	// branch-modes XXX unused?

	// constructor
	Fem() {}

	// Fem is an automatic singleton: the one-and-only instance is
	// in static member function "instance"
	static Fem& instance() {
		static Fem instance_;
		return instance_;
	}

	// member access functions
	Freev& retained() { return retained_; }
	//!! void add(const Freedom& f);  // { retained_.push_back(f); }
	// find node number "number"
	const Node* find(int number) const { return grid.find(number); }
	// get all rotational dof for "node"
	std::vector<Freedom> get_rotations(int node);


	// including the nmodes spec means use free-vibration modes to reduce
	bool freevib() { return freevib_; }
	void freevib(bool f) { freevib_ = f; }
	int nmodes() { return nmodes_; }
	bool nmodes(int n) { nmodes_ = n; return true; }
	// transformation from modal gc to nodal
	const std::vector<double>& modal() const { return modal_; }
	bool modal(const std::vector<double>& m) { modal_ = m; return true; }

	// id is used in naming the output matrices (see mid())
	std::string id() { return id_; }
	std::string id(const std::string& newid) {
		std::string oldid = id_;
		id_ = newid;
		return oldid;
	}
	std::string mid(const std::string& name);

	// parsing functions
	bool nodes_f (const Tok& p);
#ifdef NEVER // deprecate
	bool material_f (const Tok& p);
	bool beam_f (const Tok& p);
	bool mass_f (const Tok& p);
	bool mass_node_f(double mass, std::vector<double>& moi, const std::string& cg_p,
		double cg, std::vector<double>& orient, const std::vector<int>& nodes);
#endif // NEVER // deprecate
	bool gyro_f (const Tok& p);
	bool gyro_node_f (const Tok& p);
	bool pgaf_f (const Tok& p);
	bool pgaf_node_f (const Tok& p);
	int pgaffcn (pset& plt, int nr, int nc, std::vector<std::complex<Ad>>& ca);
	bool dlm_f (const Tok& p);
	bool ss_f (const Tok& p);
	//!! bool root_f (const Tok& p);
	//!! bool branch_f (const Tok& p);
	bool bma_f(const Tok& p);
	bool cms_f (const Tok& p);

	//!! bool parse_orient(const Tok& p, std::vector<double>& orient);

	SS* find_ss(const std::string& id);

	// assemble all output matrices
	void assemble();

	// these 2 return a transformation from gc to nodal dof
	std::vector<double> branch_transform();
	std::vector<double> cms_transform();
	// free-vibration with the nodal mass, stiffness and gyro
	std::vector<double> freevibration();
	// reduce the nodal matrices to modal
	void modal_reduction(const std::vector<double>& T, const std::string& id);
	//!! void cms();

	// functions for creating Custom functions
	std::string make_massfcn(std::string& entrypt);
	std::string make_gyrofcn(std::string& entrypt, std::string const& prefix="");
	std::string make_pgaffcn(std::string& entrypt, std::string const& prefix="");

	// create data & matrices necessary for amvz
	size_t nvz() { return nvz_; }		// # freedoms in vz grid
	void nvz(size_t n) { nvz_ = n; }	// i.e. rows of gct/gct_nodal
	void vzdata();
	void vzmatrices(const std::string& vzid, std::vector<double>& gct);
	void vzmatrices(const std::string& vzid, std::vector<std::complex<double>>& gct);

	void vzstif(std::vector<int>& node_numbers, std::vector<double>& coords,
		std::vector<int>& segments, std::vector<double>& gct_nodal);
	void vzpgaf(std::vector<int>& node_numbers, std::vector<double>& coords,
		std::vector<int>& segments, std::vector<double>& gct_nodal);

	//!! void store();
	void exportmm() {
		if (stif_ != nullptr)
			MM::exporter(stif_->mid(), stif_);
		if (mass_ != nullptr)
			MM::exporter(mass_->mid(), mass_);
		if (gyro_ != nullptr)
			MM::exporter(gyro_->mid(), gyro_);
	}

	friend std::ostream& operator<<(std::ostream& s, const std::vector<Freedom>& t);

}; // class Fem

std::ostream& operator<<(std::ostream& s, const std::vector<Freedom>& t);
std::ostream& operator<<(std::ostream& s, const std::vector<SS*>& t);
std::ostream& operator<<(std::ostream& s, const SS& t);

inline
void xprod(const double* a, const double* b, double* c) {
// vector cross product: c = a x b
	c[0] = a[1]*b[2] - a[2]*b[1];
	c[1] = a[2]*b[0] - a[0]*b[2];
	c[2] = a[0]*b[1] - a[1]*b[0];
}

inline
bool
is_identity(std::vector<double> t) {
// test if a (3,3) double matrix is the identity
	return (is_equal(t[0], 1.0, 8) && is_equal(t[4], 1.0, 8) &&
		is_equal(t[8], 1.0, 8));
}

template<typename T>
void append(const std::vector<T>& from, std::vector<T>& to) {
// append vector "from" to vector "to"
	Trace trc(1,"append");
	trc.dprint("appending ",from.size()," elements to ",to.size(),"-vector");
	std::vector<T> tmp{to};
	size_t newsize = tmp.size() + from.size();
	to = std::vector<T>(newsize, 0.0);
	blas_copy(tmp.size(), &tmp[0], 1, &to[0], 1);
	blas_copy(from.size(), &from[0], 1, &to[tmp.size()], 1);
	trc.dprint("returning ",to.size(),"-vector");
}

Material parse_material(const Tok& p);
bool parse_orient(const Tok& p, std::vector<double>& orient);
std::vector<Beam*> parse_beam(const Tok& p, const Material& matl, Freev& freev);
bool parse_conmass(const Tok& p, const Freev& ret,
		const std::vector<Beam*> beams, std::vector<Mass*>& masses);
void conmass(double mass, const std::vector<double>& moi,
		const std::vector<double>& cgv, std::vector<double>& elem);

#endif // FEM_H
