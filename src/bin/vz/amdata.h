//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

#ifndef AMDATA_H
#define AMDATA_H

// Data for 3D animation of structures using OpenGL
// disp_:  complex vector of displacements at each node (gctransform*gc)

#include <complex>
#include <mutex>
#include <string>
#include <thread>
#include <vector>
#include <wx/arrstr.h> // for wxArrayString

#include "matrix.h"

class Amdata {
	// nodal_disp_: a collection of picked complex displacements
	// (dx,dy,dz) at each node: 3*nnode-complex vectors
	std::vector<std::vector<std::complex<double>>> nodal_disp_;
	wxArrayString picked_names_;
	int picked_point_;  // 0b index in nodal_disp_
	int nnodes_;
	int ngc_;
	std::vector<std::complex<double>> gct_; // (3*nnode, ngc) transformation
	std::vector<double> coords_;       // nodal coordinates: (3*nnode,1) double
	std::vector<int> nodenumbers_;	// (nnode) node numbers in the same order
	std::vector<int> segidx_;	// pairs of indices into nodenumbers_
	double coord_norm_{0.0};
	double amplitude_;   // displacement scale: slider-value*nominal_scale
	double nominal_scale_;
	int period;          // for MyFrame::m_timer
	// the one-and-only instance
	static Amdata* instance_;
	// for locking while allocating instance_
	static std::mutex mutex_;
	std::string colormode_;
	int colorskew_;
	int drawundeformed_;
	// for drawing the coord axes:
	std::vector<float> origin_{0.0,0.0,0.0};
	int drawcs_{1};
protected:
	// private constructors
	Amdata(const std::string& ufname);
	Amdata(const std::string& nodes_name,
			const std::string& coords_name,
			const std::string& conn_name,
			const std::string& gct_name);
	void dostuff(std::vector<int>& nodenumbers_, std::vector<double>& coords_,
		std::vector<int>& penlift, std::vector<std::complex<double>>& gct_);

	double normalize_coord();
public:
	// delete the copy constructor & assignment operator
	Amdata(Amdata&) = delete;
	void operator=(const Amdata&) = delete;
	// access the one-and-only Amdata:
	static Amdata* instance(const std::string& path="");

	// initialize amvz: must be called first
	static void initialize(const std::string& ufname);
	static void initialize(const std::string& node_numbers,
			const std::string& coords,
			const std::string& segments, const std::string& gct);

	// returns the (0b) disp_coord_ index of this ev
	size_t add(const std::vector<std::complex<double>>& ev);
	// compute the coordinates of the displaced structure at time omegat:
	// XXX make this private, called from Amdata::render
	std::vector<double> coord_displacements(int omegat, int phase,
			std::vector<double>& colorvalue);
	int nnodes() { return nnodes_; }
	int ngc() { return ngc_; }
	static std::vector<double>& coordinates() { return instance()->coords_; }
	// compute the nodal displacements at time omegat+phase and feed
	// the segments to GL
	void render(int omegat, int phase);
	static std::string colormode(const std::string& newmode="");
	static int colorskew(int newskew=-1);
	static int npicked();
	static int picked_point(int newpoint=-1);
	void add_picked_name(const std::string& nm);
	//!! static void make_picked_names();
	wxArrayString picked_names() { return picked_names_; }
	static double amplitude(int newamp=-1);
	static int draw_undeformed(int newund=-1);
	int draw_coord_sys(int newcs=-1);
};


#endif // AMDATA_H
