//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-


#include "config.h"

#include <assert.h>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "Ad.h"
#include "blas.h"
#include "conv.h"
#include "exim.h"
#include "fem.h"
#include "fma.h"
#include "lapack.h"
#include "matrix.h"
#include "Pz.h"
#include "trace.h"

using namespace std;

void
Fem::
vzdata() {
// create 4 member arrays needed for animated-mode visualization:
//   nodes_     node numbers
//   coords_    (x,y,z) at each node in the same order as "nodes_"
//   conn_      connectivity in penlift format
//   gct_nodal_ transform from retained Freedoms to translations at visualization grid
	Trace trc(1,"vzdata");

	size_t n = retained().size();	// # Freedom in Fem

	// stiffness and gyro: this will allocate gct_nodal
	if (this->stif_ != nullptr)
		vzstif(nodes_, coords_, conn_, gct_nodal_);

	// pgaf elements
	this->vzpgaf(nodes_, coords_, conn_, gct_nodal_);

	// doublet-lattice grid
	if (this->dlm != nullptr)
		this->dlm->vzaero(nodes_, coords_, conn_, gct_nodal_, n);

	// check sizes
	int nnodes = nodes_.size();
	if (nnodes == 0)
		throw runtime_error("no node numbers in vzdata");
	this->nvz(3*nnodes);	// set # rows in gct_nodal: # dof in vz grid
	if (coords_.size() != nvz())
		throw runtime_error(vastr("vzdata: there are ",nnodes," but ",
			coords_.size()," coordinates: should be 3*",nnodes));
	if (gct_nodal_.size()%3 != 0)
		throw runtime_error(vastr("gct_nodal is the wrong size (",gct_nodal_.size(),
			"): should be a multiple of 3"));

	// save this data as 4 matrices for use by amvz with tag "_nodal"
	vzmatrices("nodal", gct_nodal_);
}

// functions for creating a picture of the model for amviz
void
Fem::
vzstif(vector<int>& node_numbers, vector<double>& coords,
	vector<int>& segments, vector<double>& gct_nodal) {
// add all stiffness elements to the gc transformation matrix (gct_nodal):
// (tx,ty,tz), node_numbers, coordinates, and connectivities (segments)
	Trace trc(1,"vzstif");

	// create the "coordinates" matrix of the amvz grid
	int maxnode{0};
	for (auto& ni : this->grid.nodes()) {
		node_numbers.push_back(ni.number());
		coords.push_back(ni.coord()[0]);
		coords.push_back(ni.coord()[1]);
		coords.push_back(ni.coord()[2]);
		maxnode = std::max(maxnode, ni.number());
	}

	// segments in penlift format: 2 node numbers followed by 0 (penlift)
	// stiffness elements...
	for (auto& bi : stif_elements) {
		segments.push_back(bi->node1());
		segments.push_back(bi->node2());
		segments.push_back(0);
	}

	// Transformation from the model dof to a visualization grid: gct_nodal.
	// Each column of gct_nodal corresponds to (and in the same order as)
	// a dof in this->retained(). The rows are tx,ty and tz for each node
	// in the "coordinates" matrix. gct_nodal is thus
	// (3*nnodes,this->retained().size())
	// and transforms from Model dof to the (tx,ty,tz) dof expected by amvz.
	// Note that some nodes may not have translational dof; if they have
	// rotational dof create an interpolation with the other node of
	// the beam element

	// Assumption: gct_nodal has not been allocated yet
	int nr = 3*node_numbers.size();
	int nc = retained().size();
	gct_nodal = vector<double>(nr*nc, 0.0);
	flaps::info("structural dof are rows 1 to ",nr," in gct_nodal");

	auto begin = retained().begin();
	auto end = retained().end();
	for (size_t k=0; k<node_numbers.size(); k++) {
		int node = node_numbers[k];
		// for each retained translational freedom for node:
		int ntrans{0};
		for (int dof=1; dof < 4; dof++) {
			auto idx = std::find(begin,end,Freedom(node,dof));
			if (idx != end) {
				int i = 3*k + dof-1;
				int j = idx - begin;
				gct_nodal[IJ(i,j,nr)] = 1.0;
				ntrans++;
			}
		}
		if (ntrans > 0) continue;
		// if no translations retained attempt to use rotations on a beam
		// but use the freedoms in retained() to get the attached structure
		for (auto& bi : stif_elements) {
			// XXX need to do this also for node == bi.node1
			if (node == bi->node2()) {
				//!! vector<int> ret = this->node_retained(bi->node2());
				// get the nodal coord for these 2 nodes from model.grid
				vector<double> x1 = this->grid.find(bi->node1())->coord();
				vector<double> x2 = this->grid.find(bi->node2())->coord();
				vector<double> orient{{x2[0]-x1[0]}, {x2[1]-x1[1]}, {x2[2]-x1[2]}};
				for (int roti=4; roti<7; roti++) {
					// find the translation columns for node1
					int trans = roti - 3;
					auto tidx = std::find(begin,end,Freedom(bi->node1(),trans));
					if (tidx != end) {
						int i = trans-1 + 3*k;
						int j = tidx - begin;
						gct_nodal[IJ(i,j,nr)] = 1.0;
					}
					// search for rotations on node1
					auto idx = std::find(begin,end,Freedom(bi->node1(),roti));
					if (idx != end) {
						int j = idx - begin;
						vector<double> disp(3,0.0);
						vector<double> qp(3,0.0);
						qp[roti-4] = 1.0;
						xprod(&qp[0], &orient[0], &disp[0]);
						for (int ti=0; ti<3; ti++) {
							int i = ti + 3*k;
							gct_nodal[IJ(i,j,nr)] = disp[ti];
						}
					}
				}
			}
		}
	}

	trc.dprintm(nr, nc, nr, gct_nodal,"gct_nodal");
}

void
Fem::
vzpgaf (vector<int>& node_numbers, vector<double>& coords,
		vector<int>& segments, vector<double>& gct_nodal) {
// given a set of pgaf elements, create 4 nodes representing a propeller,
// add these to "nodes", "coords", "conn", and "gct_nodal" for visualization
// Input:
// Output:
//   node_numbers  has had the propeller nodes added
//   coords        has had the propeller node coordinates added
//   segments
//   gct_nodal     on output gct_nodal has been allocated nr*nc and the propellor
//                 terms inserted, where nr = coords.size() and nc = retained.size()
	Trace trc(1,"vzpgaf");

	if (this->pgaf_elements.empty())
		return;

	// nc is the number of dof in the Fem grid, and the number of
	// columns in gct_nodal
	int nc = retained().size();
	// input gct_nodal is (nr,nc): compute nr - should be 3*node_numbers.size()
	int nr = gct_nodal.size()/nc;
	if (nr != (int)(3*node_numbers.size()))
		throw runtime_error(vastr("gct_nodal size mis-match: nr=",nr,", ",
			node_numbers.size()," node numbers"));
	if (nr != (int)coords.size())
		throw runtime_error(vastr("coords size ",coords.size(),", nr ",nr));

	// new nodes start numbering at maxnode+1
	int maxnode{0};
	for (auto i : node_numbers)
		maxnode = std::max(maxnode, i);

	// append each propeller
	int iprop{0};
	for (const auto& pg : this->pgaf_elements) {
		// attachment is the node number of the structure where this prop attaches
		int attach{pg.node()};

		// vector<double> T = makeT(pg.orient);
		assert(pg.T.size() == 9);
		vector<double> vec1{pg.T[3],pg.T[4],pg.T[5]};
		vector<double> vec2{pg.T[6],pg.T[7],pg.T[8]};

		// get base coordinates XXX assuming node_numbers includes grid?
		auto idx = std::find(node_numbers.begin(), node_numbers.end(), attach);
		if (idx == node_numbers.end())
			throw runtime_error(vastr("pgaf node ",attach," has not been defined"));

		// base coordinates: where the prop shaft attaches to
		// the model (3 coords per node)
		// attach_idx is the ordinal node number
		int attach_idx = idx - node_numbers.begin();
		vector<double> base{coords[3*attach_idx], coords[3*attach_idx+1],
				coords[3*attach_idx+2]};

		// the hub is attached rigidly at a distance lbar*radius along "orient"
		double t{pg.lbar*pg.radius};
		vector<double> rigid{t*pg.orient[0], t*pg.orient[1], t*pg.orient[2]};

		// create 5 new nodes: propeller blades + hub, put coords on new_coords,
		// node numbers on new_node_numbers starting with maxnode+1
		// hub:
		vector<double> new_coords;
		vector<int> new_node_numbers;
		vector<double> hub{base[0]+rigid[0], base[1]+rigid[1], base[2]+rigid[2]};
		new_coords.push_back(hub[0]);
		new_coords.push_back(hub[1]);
		new_coords.push_back(hub[2]);
		int hub_node{++maxnode};
		new_node_numbers.push_back(hub_node);
		// blade 1:
		new_coords.push_back(hub[0]+pg.radius*vec1[0]);
		new_coords.push_back(hub[1]+pg.radius*vec1[1]);
		new_coords.push_back(hub[2]+pg.radius*vec1[2]);
		new_node_numbers.push_back(++maxnode);
		// blade 2:
		new_coords.push_back(hub[0]-pg.radius*vec1[0]);
		new_coords.push_back(hub[1]-pg.radius*vec1[1]);
		new_coords.push_back(hub[2]-pg.radius*vec1[2]);
		new_node_numbers.push_back(++maxnode);
		// blade 3:
		new_coords.push_back(hub[0]+pg.radius*vec2[0]);
		new_coords.push_back(hub[1]+pg.radius*vec2[1]);
		new_coords.push_back(hub[2]+pg.radius*vec2[2]);
		new_node_numbers.push_back(++maxnode);
		// blade 4:
		new_coords.push_back(hub[0]-pg.radius*vec2[0]);
		new_coords.push_back(hub[1]-pg.radius*vec2[1]);
		new_coords.push_back(hub[2]-pg.radius*vec2[2]);
		new_node_numbers.push_back(++maxnode);

		// connectivity
		segments.push_back(attach);
		segments.push_back(hub_node);
		segments.push_back(0);
		for (int i=1; i<=4; i++) {
			segments.push_back(hub_node);
			segments.push_back(hub_node+i);
			segments.push_back(0);
		}
		
		// add the motion of the prop to "gct_nodal": the motion at the
		// prop nodes is the translation at "attach" plus r X p where
		// r is the rotation vector at "attach" and p is the vector
		// from "attach" to prop node; use the fact that
		// r X p = -p X r and the matrix form of -p X is
		//       |   0    p_z   -p_y |
		// -pX = | -p_z    0     p_x |
		//       |  p_y  -p_1     0  |

		int tx{0};
		int ty{1};
		int tz{2};
		int rx{0};
		int ry{1};
		int rz{2};
		// treat all retained dof on our attachment node	
		vector<int> columns; // 0b column numbers of attach freedoms
		int ntrans{0}; // attach has ntrans translational freedoms
		const vector<Freedom>& ret = retained();  // retained freedoms in Fem
		for (int k=0; k<nc; k++) {
			if (ret[k].node() == attach) {
				columns.push_back(k);
				if (ret[k].dof() < 4)
					ntrans++;
			}
		}

		// new_gct_nodal is (3,5,nc)
		vector<double> new_gct_nodal(15*nc);

		if (ntrans == 0) {
		// copy translational rows of attach in gct_nodal to new_gct_nodal
			for (int i=0; i<3; i++) {
				int ia = i + 3*attach_idx;
				for (int j=0; j<5; j++) {
					for (int k=0; k<nc; k++) {
						new_gct_nodal[IJK(i,j,k,3,5)] = gct_nodal[IJ(ia,k,nr)];
					}
				}
			}
		}
		// XXX use columns
		for (int k=0; k<nc; k++) {
			int nodek{retained()[k].node()};
			int dofk{retained()[k].dof()};
			if (nodek == attach) {
				for (int j=0; j<5; j++) {
					// new_coords is (3,5)
					vector<double> p(3,0.0);
					for(int i=0; i<3; i++)
						p[i] = new_coords[IJ(i,j,3)] - base[i];
					int i = dofk-1;
					switch (dofk) {
						case 1:
						case 2:
						case 3:
							new_gct_nodal[IJK(i,j,k,3,5)] = 1.0;
							break;
						case 4:
							new_gct_nodal[IJK(ty,j,k,3,5)] = -p[rz];
							new_gct_nodal[IJK(tz,j,k,3,5)] = p[ry];
							break;
						case 5:
							new_gct_nodal[IJK(tx,j,k,3,5)] = p[rz];
							new_gct_nodal[IJK(tz,j,k,3,5)] = -p[rx];
							break;
						case 6:
							new_gct_nodal[IJK(tx,j,k,3,5)] = -p[ry];
							new_gct_nodal[IJK(ty,j,k,3,5)] = p[rx];
							break;
					}
				}
			}
		}
		trc.dprint("new node numbers ",new_node_numbers);
		trc.dprintm(3,5,3,&new_coords[0],"new coords");
		trc.dprintm(15,nc,15,&new_gct_nodal[0],"new gct_nodal");

		// add to the input vectors (append is in fem.h)
		append(new_node_numbers, node_numbers);
		append(new_coords, coords);

		// the new gct_nodal is (nr+15,nc)
		int newnr = nr+15;
		flaps::info("propeller ",++iprop," is rows ",nr+1," to ",newnr," in gct_nodal");
		vector<double> tmp(newnr*nc);
		for (int j=0; j<nc; j++)
			blas_copy(nr, &gct_nodal[j*nr], 1, &tmp[j*newnr], 1);
		for (int j=0; j<nc; j++)
			blas_copy(15, &new_gct_nodal[j*15], 1, &tmp[IJ(nr,j,newnr)], 1);
		gct_nodal = tmp;
		nr = newnr;
	}
}

void
Fem::
vzmatrices(const string& vzid, vector<double>& gct) {
// create and store a Fem structure tagged with "vzid" and
// descriptions nodes|coords|conn|gct
	Trace trc(1,"vzmatrices");

	if (gct.empty())
		throw runtime_error(vastr("vzmatrices (",vzid,") called with empty gct"));

	// store with extension ".vzid" XXX should we require a vzid?
	string ext;
	if (!vzid.empty()) {
		ext = "." + vzid;
		trc.dprint("extension: ",ext);
	}
	Fma fma;
	fma.nodes = nodes_;
	fma.coords = coords_;
	fma.segidx = UF::penlift2segidx(conn_, nodes_);
	// gct is a (nvz,nc) real vector - convert to complex
	size_t nr = this->nvz();
	size_t nc = gct.size()/nr;
	if (nr*nc != gct.size())
		throw runtime_error(vastr("the size of gct (",gct.size(),") should be "
			"a multiple of ",nr));
	fma.gct = vector<complex<double>>(gct.size(), complex<double>(0.0));
	blas_copy(gct.size(), gct.data(), 1, (double*)fma.gct.data(), 2);
	fma.vzid = vzid;
	fma.store();
	flaps::info("Created an Fma:\n",fma);

	// create a UF file for command-line visualization
	string path{"fma"};
	if (!ext.empty())
		path += ext;
	UF::exporter(path, nodes_, coords_, conn_, gct);

	// write gct as a .mm file locally
	MM::exporter(path+".mm", "gct", &gct[0], nr, nc);
}
