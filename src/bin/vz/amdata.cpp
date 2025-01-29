//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

#include <wx/defs.h>
#include <wx/app.h>
#include <wx/menu.h>
#include <wx/dcclient.h>
#include <wx/wfstream.h>
#if wxUSE_ZLIB
#include <wx/zstream.h>
#endif

#include <wx/glcanvas.h>
#include <wx/wx.h>
#include <wx/gtk/print.h>
#include <GL/glu.h>

#undef Complex    // X11 macro

#include "trackball.h"
#include "amdata.h"
#include "amvz.h"
#include "exim.h"
#include "plotcurve.h"

#if !wxUSE_GLCANVAS
    #error "OpenGL required: set wxUSE_GLCANVAS to 1 and rebuild the library"
#endif

using namespace std;

mutex Amdata::mutex_;

// new stuff for handling mode menu
//!! wxArrayString Amdata::picked_names_;
// XXX maybe Frame should be a member of Amdata?
extern MyFrame* Frame;

// Only one Amdata, accesssed with Amdata::instance()
// XXX for backward compatability
// XXX change all public functions to static which just call instance()->real_version
Amdata* Amdata::instance_ = nullptr;
Amdata*
Amdata::
instance(const string& ufname) {
	if (instance_ == nullptr)
		throw runtime_error("Amdata::initialize() has not been called");
	return instance_;
}

void
Amdata::
initialize(const std::string& ufname) {
	if (Amdata::instance_ != nullptr)
		throw runtime_error("Amdata::initialize called more than once");
	instance_ = new Amdata(ufname);
}

void
Amdata::
initialize(const string& nodes, const string& coords,
		const string& conn, const string& gct) {
	if (Amdata::instance_ != nullptr)
		throw runtime_error("Amdata::initialize called more than once");
	instance_ = new Amdata(nodes, coords, conn, gct);
}


vector<int>
penlift2segments(const vector<int>& penlift) {
// convert connectivities in penlift format to segment format
	Trace trc(2,"penlift2segments");
	vector<int> rval;
	int prev{0};
   for (size_t i=0; i<penlift.size(); i++) {
      if (penlift[i] == 0) {
         prev = 0;
         continue;
      }
      if (prev > 0) {
         rval.push_back(prev);
         rval.push_back(penlift[i]);
         prev = penlift[i];
      } else if (penlift[i+1] > 0) {
         rval.push_back(penlift[i++]);
         rval.push_back(penlift[i]);
			prev = penlift[i];
		}
	}
	trc.dprint("connectivities in segment format",rval);
	return rval;
}

void
Amdata::
dostuff(vector<int>& nodenumbers_, vector<double>& coords_,
	vector<int>& penlift, vector<complex<double>>& gct_) {
// do some common setup for both constructors
	Trace trc(1,"dostuff");

	trc.dprint("coordinates: ",coords_);
	trc.dprint("connectivity: ",penlift);
	trc.dprint("node numbers: ",nodenumbers_);

	nnodes_ = nodenumbers_.size();
	ngc_ = gct_.size()/(3*nnodes_);

	// gct must be dimensioned (3*nnodes*ngc) where ngc is the
	// number of generalized coordinates; coords are (3*nnodes)
	// and must be in the same order as gct.
	if (gct_.size() != (size_t)3*nnodes_*ngc_)
		throw runtime_error(vastr("dimensions do not agree: gct_ size ",
			gct_.size(),", ",nnodes_));

	// normalize the coordinates and shift the origin to the bounding-box center
	coord_norm_ = normalize_coord();
	trc.dprint("normalized & shifted coords: ",coords_);
	trc.dprint("coord norm ",coord_norm_);

	// the connectivities are in penlift format - change them to 0b indices
	// into the nodenumbers array; then the start of start of ordinal
	// node j is 3*segidx[j].
	vector<int> seg = penlift2segments(penlift);
	segidx_.clear();
	for (auto ci : seg) {
		auto idx = find(nodenumbers_.begin(), nodenumbers_.end(), (int)ci);
		if (idx == nodenumbers_.end())
			throw runtime_error(vastr("node ",(int)ci," in the connectivity"
				" has not been defined"));
		segidx_.push_back((int)(idx - nodenumbers_.begin()));
	}
	trc.dprint("segment indices (segidx)",segidx_);

	colormode_ = "abs value";
	colorskew_ = 5;  // no skewing
	// nodal displacements will be scaled by amplitude_, which is
	// changed in MyFrame::OnAmpSlider (amvz.cpp)
	nominal_scale_ = 0.005;
	amplitude_ = nominal_scale_*pow(2.0, double(Init_amp));
	drawundeformed_ = 0;
	// set initial period for m_timer XXX necessary?
	period = Default_period/Init_freq;
}

Amdata::
Amdata(const string& ufname) {
// Amdata constructor (private) importing data from a Universal file.
// 4 arrays must be in the Universal file:
// 1) nodenumbers_  (nnodes int) node numbers
// 2) coords_       (3*nnodes double) coordinates of each node in, and
//                  in the same order as "nodenumbers_"
// 3) penlift       sequences of node numbers with zeros representing
//                  penlifts
// 4) gct           (3*nnodes*ngc) transformation from the ngc generalized
//                  (or physical) coordinates to (x,y,z) displacements at
//                  "coords_"
	Trace trc(1,"Amdata uf constructor");

	// import the uf file
	vector<int> penlift;
	UF::importer(ufname, nodenumbers_, coords_, penlift, gct_);
	
	// checks
	if (gct_.empty())
		throw runtime_error(vastr(ufname," does not have gct"));
	if (coords_.empty())
		throw runtime_error(vastr(ufname," does not have coordinates"));
	if (penlift.empty())
		throw runtime_error(vastr(ufname," does not have connectivity"));
	if (nodenumbers_.empty())
		throw runtime_error(vastr(ufname," does not have node numbers"));


#ifdef STANDALONE
	// visualizing the modes in a UF file: the gct columns *are* the modes
	// put each column of gctransform into nodal_disp_ in reverse
	// order: amvz orders by the most recent first
	// XXX try normal ordering
	size_t nr = 3*nodenumbers_.size();
	size_t nc = gct_.size()/nr;
	vector<complex<double>> disp(nr);
	complex<double>* cele = &gct_[0];
	for (size_t j=0; j<nc; j++) {
		//!! blas_copy(nr, &cele[IJ(0,nc-j-1,nr)], 1, &disp[0], 1);
		blas_copy(nr, &cele[IJ(0,j,nr)], 1, &disp[0], 1);
		double disp_norm = blas_scnrm2(nr, &disp[0], 1);
		if (disp_norm == 0.0)
			flaps::warning("mode ",nc-j," is zero");
		nodal_disp_.push_back(disp);
		picked_names_.Add(wxString(vastr("mode ",j+1)));
	}
	// start visualizing on the first one
	picked_point_ = 0;
#endif // STANDALONE
	// finish construction
	dostuff(nodenumbers_, coords_, penlift, gct_);
}

Amdata::
Amdata(const string& nodes_name, const string& coords_name,
		const string& conn_name, const string& gct_name) {
// Amdata constructor (private) reading data from the Flaps database.
// 4 matrices are necessary:
// 1) nodes.vzid    (nnodes,1) node numbers
// 2) coords.vzid   (3,nnodes) coordinates of each node in "nodenumbers"
// 3) conn.vzid     list of segments: pairs of node numbers separated by
//                  zero (penlift), the format required by OpenGL
// 4) gct.vzid      (3*nnodes,ngc) transformation from the ngc generalized
//                  (or physical) coordinates to (x,y,z) displacements at
//                  "coordinates"
// The data in these matrices is stored in vectors of doubles, so the int
// arrays nodenumbers and penlift must be cast to ints.
	Trace trc(1,"Amdata (nodes...) constructor");

	Matrix* nodes = new Matrix(nodes_name);
	Matrix* coords = new Matrix(coords_name);
	Matrix* conn = new Matrix(conn_name);
	Matrix* gct = new Matrix(gct_name);

	// coords
	coords_ = coords->data();	// copy

	// conn and nodenumbers are int arrays stored as doubles
	vector<int> penlift;
	for (auto d : conn->data())
		penlift.push_back(d);
	for (auto d : nodes->data())
		nodenumbers_.push_back(d);
	// gct_ is a vector<complex> but gct may be real
	size_t nel = gct->rsize()*gct->csize();
	if (gct->is_complex()) {
		complex<double>* cp = gct->celem();
		for (size_t i=0; i<nel; i++)
			gct_.push_back(cp[i]);
	} else {
		double* rp = gct->elem();
		for (size_t i=0; i<nel; i++)
			gct_.push_back(rp[i]);
	}

	// finish construction
	dostuff(nodenumbers_, coords_, penlift, gct_);
}

size_t
Amdata::
add(const vector<complex<double>>& ev) {
// add gctransform*ev to the list of nodal displacements, return
// the 0b index of the result in the nodal_disp_ array.
	Trace trc(1,"Amdata::add");
	complex<double> alpha{1.0};
	complex<double> beta{0.0};
	int m{3*nnodes()};
	int n{ngc()};
	// ev may be larger than the number of columns but not smaller - e.g.
	// if there are extra controls equations
	if ((int)ev.size() < n) {
		string exc{vastr("mismatch: ",ev.size()," gc but ",n," transform columns")};
		trc.dprint("throwing exception: ",exc);
		throw runtime_error(exc);
	}
	complex<double>* gcxf = &gct_[0];
	vector<complex<double>> disp(m);
	trc.dprint("eigenvector:",ev);
	blas_cgemv ("n", m, n, alpha, gcxf, m, &ev[0], 1, beta, &disp[0], 1);
	size_t rval = nodal_disp_.size();
	picked_point_ = rval;     // default picked point
	nodal_disp_.push_back(disp);
	return rval;
}

void
Amdata::
add_picked_name(const std::string& nm) {
// add a name to picked_names_, presumably the name of the most recently
// added disp (Amdata::add) XXX probably should include this fcn in add()
	picked_names_.Add(wxString(nm));
	// delete the existing button, create new one
	// XXX ModeMenu dups one in amvz: doesn't need to be gobal?
	if (Frame != nullptr) {
		wxChoice* ModeMenu = new wxChoice(Frame->GetToolBar(), ID_MODE_NUMBER,
            wxDefaultPosition, wxDefaultSize, picked_names_);
		Frame->GetToolBar()->RemoveTool(ID_MODE_NUMBER);
		Frame->GetToolBar()->AddControl(ModeMenu);
	}
}

vector<double>
Amdata::
coord_displacements(int omegat, int phase, vector<double>& colorvalue) {
// compute the (real) displacement at each node at time omegat
//   x[k] = Real(cos(omegat)+i*sin(omegat)*nodal_disp[k] + coord[k]
//       k = 0:n3-1
	Trace trc(1,"Amdata::coord_displacements");
	constexpr double deg2rad = atan(1.0)/45.0;
	constexpr double pi = 4.0*atan(1.0);
	static int visit{0};
	int n3 = 3*nnodes();

	vector<double> rval(n3, 0.0);
	Amdata* amd = instance();

	if ((int)colorvalue.size() != nnodes())
		throw runtime_error(vastr("colorvalue wrong size (",colorvalue.size(),
				") should be ",nnodes()));

	// copy the current picked nodal_disp_
	vector<complex<double>> disp = nodal_disp_[picked_point_];

	if (visit == 0)
		trc.dprint("omegat ",omegat,", disp = ",disp);
	visit++;

	// normalize the complex nodal displacements to amplitude_*coord_norm_
	// use the infinity norm and normalize colorvalues with it below
	double dispnorm = abs(disp[blas_icamax(n3, &disp[0], 1)-1]);
	double eps{1.0e-6};
	// return zero if tiny displacements
	if (dispnorm*amplitude_ <= eps) {
		return rval;
	}
	complex<double> ct = amplitude_*coord_norm_/dispnorm;
	blas_scal(n3, ct, &disp[0], 1);

	double omtp{(omegat + phase)*deg2rad};
	double scr{cos(omtp)};
	double sci{sin(omtp)};
	complex<double> csc(scr, sci);

	trc.dprint("omegat+phase ",omtp," rad, (cos,sin) ",csc);

	vector<complex<double>> cd(3);
	for (int i=0; i<nnodes(); i++) {
		int k = 3*i;
		cd[0] = disp[k]*csc;
		cd[1] = disp[k+1]*csc;
		cd[2] = disp[k+2]*csc;
		double x = cd[0].real();
		double y = cd[1].real();
		double z = cd[2].real();
		rval[k] = x + coords_[k];
		rval[k+1] = y + coords_[k+1];
		rval[k+2] = z + coords_[k+2];
		// double normd = blas_scnrm2(3,&cd[0],1);
		double normd = sqrt(x*x+y*y+z*z);
		int idx = blas_icamax(3,&cd[0],1)-1;
		// set a color value for this node
		if (amd->colormode_ == "abs value") {
			colorvalue[i] = normd;
		} else if (amd->colormode_ == "phase") {
			colorvalue[i] = atan2(cd[idx].imag(), cd[idx].real())/(2.0*pi) + 0.5;
		} else if (Amdata::colormode() == "abs phase") {
			colorvalue[i] = abs(atan2(cd[idx].imag(), cd[idx].real())/pi);
		}
		// XXX trying to find node 199 1b: set all colors to 0, 199 to 1
		// if (i == 32)
		// 	colorvalue[i] = 1.0;
		// else
		// 	colorvalue[i] = 0.0;
	}

	// scale colorvalues so that the largest colorvalue over an
	// oscillation cycle is 1.0. The largest (normalized) displacement 
	// over a cycle is amplitude_*coord_norm_
	blas_scal(nnodes(), 1.0/(amplitude_*coord_norm_), &colorvalue[0], 1);

	// skewed color variation
	if (amd->colorskew_ != 5) {
		double cs = ((double)amd->colorskew_)/5.0;
		for (int i=0; i<nnodes(); i++) {
			double ci = colorvalue[i];
			if (ci > 0.0)
				ci = pow(abs(10.0*ci), cs)/10.0;
			else if (ci < 0.0)
				ci = -pow(abs(10.0*ci), cs)/10.0;
			colorvalue[i] = ci;
		}
	}

	trc.dprint("nodal disp:",rval);
	return rval;
}

void
color_values (double t, float& red, float& green, float& blue) {
// given a number (t) between 0 and 1, set the rgb values
// so that they vary linearly with t=0 blue, t=0.5 green,
// and t=1 red

	if (t > 1.0)
		t = 1.0;
	if (t < -1.0)
		t = -1.0;

	if (t < 0.5) {
		blue = 1.0 - 2.0*t;
		green = 2.0*t;
		red = 0.0;
	} else {
		blue = 0.0;
		green = -2.0*t + 2.0;
		red = 2.0*t - 1.0;
	}
}

void
Amdata::
render(int omegat, int phase) {

	// if colormode == "none" use the reverse of the bg color
	float red = 1.0 - BGRed;
	float green = 1.0 - BGGreen;
	float blue = 1.0 - BGBlue;
	//!! Amdata* amd = Amdata::instance();

    /* clear color and depth buffers */
	glClearColor(BGRed, BGGreen, BGBlue, BGAlpha);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	
	// antialiasing
	glEnable(GL_LINE_SMOOTH);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	// glBlendFunc(1.0, 0.0);
	glHint(GL_LINE_SMOOTH_HINT, GL_DONT_CARE);

	// colorvalue is a vector of color values ranging from 0-1 to
	// be used by color_values to give rgb values for GL
	//!! vector<double> colorvalue(amd->nnodes_, 0.0);
	vector<double> colorvalue(this->nnodes_, 0.0);

	// compute the coordinates of the nodes at time omegat+phase
	vector<double> dispnodes = this->coord_displacements(omegat,
			phase, colorvalue);

	//! static double* vertexEnergy = 0;

	// we always have segidx_ which are pairs of 0b indices into the
	// nodenumbers_ array, i.e. the ordinal numbers of nodes.
	// The start of (x,y,z) data for ordinal node j in coordinates,
	// gct, and dispnodes arrays is thus 3*j
	// This is a more useful format for connectivity than either Universal
	// File's penlift format or GL's segment format.
	glLineWidth(2.0);
	glBegin(GL_LINES);
	for (size_t i=0; i<this->segidx_.size(); i++) {
		int k = 3*this->segidx_[i]; // start of this node
		// the Nodes vector is assumed to be in the same order as dispnodes
		float x = dispnodes[k];
		float y = dispnodes[k+1];
		float z = dispnodes[k+2];

		// set the color for this node and the next
		color_values(colorvalue[k/3], red, green, blue);

		glColor3f(red, green, blue);
		glVertex3f(x, y, z);
		i++;
		k = 3*this->segidx_[i];
		x = dispnodes[k];
		y = dispnodes[k+1];
		z = dispnodes[k+2];
		glVertex3f(x, y, z);
	}
	// draw the coord system if requested
	if (this->draw_coord_sys()) {
		vector<float> x;
		for (size_t i=0; i<3; i++) {
			x = origin_;
			color_values(1.0-i*0.5, red, green, blue);
			glColor3f(red, green, blue);
			glVertex3f(x[0], x[1], x[2]);
			x[i] += 0.1;
			glVertex3f(x[0], x[1], x[2]);
		}
	}
	glEnd();
	// draw the undeformed grid if requested
	if (this->drawundeformed_) {
		double* coord = &this->coords_[0];
		// red = green = blue = 0.5;  // gray
		red = 1.0 - BGRed;
		green = 1.0 - BGGreen;
		blue = 1.0 - BGBlue;
		glColor3f(red, green, blue);
		glLineWidth(0.2);
		glBegin(GL_LINES);
		for (size_t i=0; i<this->segidx_.size(); i++) {
			int k = 3*this->segidx_[i];
			float x = coord[k];
			float y = coord[k+1];
			float z = coord[k+2];
			glVertex3f(x, y, z);
			i++;
			k = 3*this->segidx_[i];
			x = coord[k];
			y = coord[k+1];
			z = coord[k+2];
			glVertex3f(x, y, z);
		}
		glEnd();
		// glDisable(GL_LINE_STIPPLE);
	}

	// glFlush();
	// wxGLCanvas::SwapBuffers();
} // Amdata::render()

string
Amdata::
colormode(const string& newmode) {
	string rval = instance()->colormode_;
	if (!newmode.empty())
		instance()->colormode_ = newmode;
	return rval;
}

int
Amdata::
colorskew(int newskew) {
	int rval = instance()->colorskew_;
	if (newskew >= 0)
		instance()->colorskew_ = newskew;
	return rval;
}

int
Amdata::
npicked() {
	return instance()->nodal_disp_.size();
}

int
Amdata::
picked_point(int newpoint) {
	int rval = instance()->picked_point_;
	if (newpoint >= 0)
		instance()->picked_point_ = newpoint;
	return rval;
}

double
Amdata::
amplitude(int newamp) {
	Amdata* amd = instance();
	double rval = amd->amplitude_;
	if (newamp >= 0)
		amd->amplitude_ = amd->nominal_scale_*pow(2.0, double(newamp));
	return rval;
}

int
Amdata::
draw_undeformed(int newund) {
// change drawundeformed_ setting, return current
	int rval = instance()->drawundeformed_;
	if (newund >= 0)
		instance()->drawundeformed_ = newund;
	return rval;
}

int
Amdata::
draw_coord_sys(int newcs) {
// change drawcs_ setting, return current
	int rval = instance()->drawcs_;
	if (newcs >= 0)
		instance()->drawcs_ = newcs;
	return rval;
}

double
Amdata::
normalize_coord() {
// Normalize the nodal coordinates by scaling so that the largest
// dimension of the enclosing box is 1.0 and the enclosing box
// is centered on the origin
	//!! vector<double>& coord = coord_->data();
	double rval{1.0};

	double big = std::numeric_limits<double>::max();
	double xmax=-big, xmin=big;
	double ymax=-big, ymin=big;
	double zmax=-big, zmin=big;
	for (int i=0; i<nnodes(); i++) {
		int k = 3*i;
		xmin = std::min(xmin, coords_[k]);
		ymin = std::min(ymin, coords_[k+1]);
		zmin = std::min(zmin, coords_[k+2]);
		xmax = std::max(xmax, coords_[k]);
		ymax = std::max(ymax, coords_[k+1]);
		zmax = std::max(zmax, coords_[k+2]);
	}
	double scale = abs(xmax - xmin);
	scale = std::max(scale, abs(ymax - ymin));
	scale = std::max(scale, abs(zmax - zmin));
	scale *= 2.0;
	double xc = (xmax + xmin)/2.0;
	double yc = (ymax + ymin)/2.0;
	double zc = (zmax + zmin)/2.0;
	for (int i=0; i<nnodes(); i++) {
		int k = 3*i;
		coords_[k] -= xc;
		coords_[k+1] -= yc;
		coords_[k+2] -= zc;
		coords_[k] /= scale;
		coords_[k+1] /= scale;
		coords_[k+2] /= scale;
	}
	// save the origin for drawing the coordinate system
	origin_.clear();
	origin_.push_back(-xc/scale);
	origin_.push_back(-yc/scale);
	origin_.push_back(-zc/scale);
	return rval;
}
