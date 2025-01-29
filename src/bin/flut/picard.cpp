//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

// For n-dimensional space evenly split into m intervals:
// Number of Simplex's in one interval = n!
// Total number of Simplexes = n!m^n
// Number of boundary faces = 2m^{n-1}n!
// So for example a 6th order problem with 4 intervals
// has 2*4^5*6! = 1,474,560 faces
//
// References
// Bowyer, A., "Computing Dirichlet tessellations", The Computer Journal,
//    24(2), pp. 162-166, 1981
// Watson, D.F., "Computing the n-dimensional Delaunay tesselation with
//    application to Voronoi polytopes", The Computer Journal,
//    24(2), pp. 167-172, 1981
// Stenger, F., "Computing the topological degree of a mapping in R^n",
//    Numerische Mathematik, 25(1), pp. 23-38, 1975
// Picard, E., "Sur le nombre des racines communes a plusieurs
//    equations simultanees", J. Math Pures Appl., 4(8), pp 5-24, 1892

#include "config.h"
#include "exim.h"
#include "lapack.h"
#include "picard.h"
#include "ridders.h"

using namespace std;

static bool Check_jacobian{false};

static int Nfcn = 0;   // number function calls
static int Nperm = 0;  // number of boundary faces permuted

void
callfcn(const vector<double>& x, vector<double>& f, vector<double>& jac,
	int (*fcn)(vector<double> const& x, vector<double>& f, vector<double>& jac));

static
vector<Cube>
triangulate (vector<vector<double> >& pts);

static
int
sgn (double a);

static double
angleReal (int n, const double *a, const double *b);

static bool
indices(int idx, vector<int> const& npts, vector<int>& pt);

double
check_jac(const vector<double>& x, Userfcn fcn, size_t& maxi, size_t& maxj);

vector<Face>
regionFaces (vector<Cube>& cubes, int idx, double idxval);

#ifdef NEVER // unused
static
void
face_coord(Face const& face,
		vector<double> const& xmin, vector<double> const& xmax,
		vector<int>& xlower, vector<int>& xupper) {
// determine which coord is constant
	Trace trc(1,"face_coord");
	size_t n = face.size();
	vector<int> xl(n,0);
	vector<int> xu(n,0);
	size_t i;

	assert(xmin.size() == n);
	assert(xmax.size() == n);
	assert(xlower.size() == n);
	assert(xupper.size() == n);

	for (i=0; i<n; i++) {      // each Point
		for (size_t j=0; j<n; j++) {   // each coord
			if (is_equal(face[i][j], xmin[j], 8))
				xl[j]++;
			if (is_equal(face[i][j], xmax[j], 8))
				xu[j]++;
		}
	}
	for (i=0; i<n; i++) {
		if (xl[i] == (int)n)
			xlower[i]++;
		if (xu[i] == (int)n)
			xupper[i]++;
	}
}
#endif // NEVER : unused

#ifdef NEVER // unused
static
void
printPicardFace(Face const& face) {
// if debug is on and this face has all Points at x[n]==+/-1
// print it
	Trace trc(1,"printPicardFace");
	if (trc() < 0) return;
	size_t n = face.size();
	size_t nm1 = n -1;
	int lim{0};
	for (size_t i=0; i<n; i++) {      // each Point
		// coord n-1 must be either +1 or -1; if it is not the
		// same as previous coord quit
		if (lim == 0)
			lim = (int)face[i][nm1];
		if (lim != (int)face[i][nm1])
			return;
	}
	trc.dprint("Picard (",lim,") boundary face: ",face);
}
#endif // NEVER : unused

static
int
countpts(vector<vector<double> > const& pts) {
	int rval{0};
	for (auto& p : pts)
		rval += p.size();
	return rval;
}

static
vector<int>
coordpts (vector<vector<double> > const& pts) {
	vector<int> rval;
	for (auto& p : pts)
		rval.push_back(p.size());
	return rval;
}

bool
Point::
operator==(Point const& p) const {
// two Points are equal if each coord is equal
	assert(size() == p.size());
	for (size_t i=0; i<size(); i++) {
		if (operator[](i) != p[i]) {
			return false;
		}
	}
	return true;
}

string
Point::
toString() const {
	ostringstream os;
	size_t i;
	size_t n = size();
	for (i=0; i<n; i++) {
		os << "  " << operator[](i);
	}
	return os.str();
}

static double
angleReal (int n, const double *a, const double *b) {
// Compute the angle (radians) between two vectors
	double rval = 0.0;
	double tcos, tsin, t;
	double abnorm = blas_snrm2(n,a,1)*blas_snrm2(n,b,1);
	double eps = sqrt(std::numeric_limits<double>::epsilon());

	if (abnorm <= eps)
		return 0.0;

	blas_dot (n, a, 1, b, 1, tcos);
	tcos /= abnorm;
	if (tcos > 1.0) tcos = 1.0;
	if (tcos < -1.0) tcos = -1.0;

	t = 1.0 - tcos*tcos;
	if (t < 0.0)
		tsin = 0.0;
	else
		tsin = sqrt(t);
	if (tsin != 0.0 || tcos != 0.0) {
		rval = atan2 (tsin, tcos);
	}
	rval = fabs(rval);
	return rval;
}

static
double
getangle (vector<double> const& a, vector<double> const& b) {
// compute the angle (in degrees) between 2 vectors
	Trace trc(1,"getangle");
	double dpr = 45.0/atan(1.0);
	double rval{0.0};
	trc.dprint("angle between ",a," and ",b);
	rval = angleReal(a.size(), &a[0], &b[0]);
	rval = dpr*abs(rval);
	trc.dprint("returning ",rval," degrees");
	return rval;
}

static
void
insertpt(vector<double> const& a, vector<double> const& b,
		vector<double> const& mid, vector<vector<double> >& pts) {
// insert "mid" between "a" and "b" in "pts"
	Trace trc(1,"insertpt");
	size_t n = a.size();
	double eps = sqrt(numeric_limits<double>::epsilon());

	trc.dprintvv(a,b,"between");
	// trc.dprintv(a,"between");
	// trc.dprintv(b,"and");

	for (size_t i=0; i<n; i++) {
		// if coord i changes:
		if (abs(b[i]-a[i]) > eps) {
			// find j such that pts[i][j] < mid[i] < pts[i][j+1]
			// (coord i's axis points), which are assumed
			// to be ordered increasing
			for (size_t j=1; j<pts[i].size(); j++) {
				// check that is not already there
				if (is_equal(mid[i], pts[i][j], 8)) {
					trc.dprint(mid[i]," is already in coord ",i);
					break;
				}
				if (is_greaterthan(pts[i][j], mid[i], 8)) {
					vector<double>::iterator p = pts[i].begin() + j;
					pts[i].insert(p, mid[i]);
					trc.dprintv(pts[i],"added ",mid[i]," to coord ",i);
					break;
				}
			}
		}
	}
}

static
int
addpts(vector<double> const& a, vector<double> const& b,
	vector<double> const& fa, vector<double> const& fb,
	Userfcn fcn, vector<vector<double> >& pts, int depth) {
// XXX args should be: (Point& a, Point& b, fcn, pts, depth)
// add points between a and b to the coord-axis points "pts"
// to keep the angle between function values at adjacent points
// below 90 degrees.
// add a mid-point, recursively add points to either
// side if the angle is still too large
// Returns the number of Points added
	Trace trc(1,"addpts depth ",depth);
	size_t n = a.size();
	assert(b.size() == n);
	assert(fa.size() == n);
	assert(fb.size() == n);
	vector<double> mid(n);
	vector<double> fmid(n);
	vector<double> jac(n*n);
	// int sgndet;

	vector<int> ncoordpts = coordpts(pts);
	int npts = countpts(pts);

	for (size_t i=0; i<n; i++) {
		mid[i] = (b[i] + a[i])/2.0;
	}
	// insert mid into the axis points
		insertpt(a, b, mid, pts);

	// evaluate f(mid)
	callfcn(mid, fmid, jac, fcn);

	// check angle between fa and f(mid)
	double angle = getangle(fa, fmid);
	if (angle > 90.0) {
		addpts(a, mid, fa, fmid, fcn, pts, depth+1);
	}
	// check angle between f(mid) and fb
	angle = getangle(fmid, fb);
	if (angle > 90.0) {
		addpts(mid, b, fmid, fb, fcn, pts, depth+1);
	}

	int rval = countpts(pts) - npts;
	vector<int> addedpts;
	for (size_t i=0; i<pts.size(); i++) {
		addedpts.push_back(pts[i].size() - ncoordpts[i]);
	}
	trc.dprintv(addedpts, "# added points in each coord:");
	trc.dprint("added ",rval," Points");
	return rval;
}


#ifdef NEVER // unused
int
intDet (size_t n, int* A) {
// returns the determinant of integer matrix A
	Trace trc(2,"intDet");
	int rval =  0;
	size_t i, j, kj;
	size_t nm1 = n - 1;

	static int depth = 0;
	if (depth == 0)
	depth++;

	if (n == 1) {
		rval = A[0];
	} else if (n == 2) {
		rval = (A[0]*A[3] - A[1]*A[2]);
	} else if (n == 3) {
		rval = A[0]*(A[4]*A[8] - A[5]*A[7]) - A[3]*(A[1]*A[8] - A[2]*A[7])
			+ A[6]*(A[1]*A[5] - A[2]*A[4]);
	} else if (n == 4) {
	 	int d31 = A[5]*(A[10]*A[15] - A[11]*A[14]) - A[9]*(A[6]*A[15] - A[7]*A[14])
	 		+ A[13]*(A[6]*A[11] - A[7]*A[10]);
	 	int d32 = A[1]*(A[10]*A[15] - A[11]*A[14]) - A[9]*(A[2]*A[15] - A[3]*A[14])
	 		+ A[13]*(A[2]*A[11] - A[3]*A[10]);
	 	int d33 = A[1]*(A[6]*A[15] - A[7]*A[14]) - A[5]*(A[2]*A[15] - A[3]*A[14])
	 		+ A[13]*(A[2]*A[7] - A[3]*A[6]);
	 	int d34 = A[1]*(A[6]*A[11] - A[7]*A[10]) - A[5]*(A[2]*A[11] - A[3]*A[10])
	 		+ A[9]*(A[2]*A[7] - A[3]*A[6]);
		rval = A[0]*d31 - A[4]*d32 + A[8]*d33 - A[12]*d34;
	} else if (n == 5) {
		int B[16];
		int scale = 1;
		for (j=0; j<5; j++) {
			size_t k = 0;
			for (kj=0; kj<n; kj++) {
				if (kj != j) {
					for (i=1; i<n; i++) {
						B[IJ(i-1,k,nm1)] = A[IJ(i,kj,n)];
					}
					k++;
				}
			}
			rval += scale*A[IJ(0,j,n)]*intDet(4, &B[0]);
			scale = -scale;
		}
	} else if (n == 6) {
		int B[25];
		int scale = 1;
		for (j=0; j<n; j++) {
			size_t k = 0;
			for (kj=0; kj<n; kj++) {
				if (kj != j) {
					for (i=1; i<n; i++) {
						B[IJ(i-1,k,nm1)] = A[IJ(i,kj,n)];
					}
					k++;
				}
			}
			rval += scale*A[IJ(0,j,n)]*intDet(nm1, &B[0]);
			scale = -scale;
		}
	} else {
		int scale = 1;
		vector<int> M(nm1*nm1);
		for (j=0; j<n; j++) {
			size_t k = 0;
			for (kj=0; kj<n; kj++) {
				if (kj != j) {
					for (i=1; i<n; i++) {
						M[IJ(i-1,k,nm1)] = A[IJ(i,kj,n)];
					}
					k++;
				}
			}
			rval += scale*A[IJ(0,j,n)]*intDet(nm1, &M[0]);
			scale = -scale;
		}
	}
	depth--;

	if (depth == 0) {
		delete time;
	}
	trc.dprint("returning ",rval);
	return rval;
}
#endif // NEVER : unused

double
Simplex::
orientation(int& expon) const {
// Returns the determinant of the (n,n) matrix whos
// columns are the elements of each Point in the Face
	Trace trc(2,"Simplex::orientation");
	size_t np1 = this->size();
	size_t n = np1 - 1;
	size_t i, j;

	assert(n == (*this)[0].size());

	vector<double> A(np1*np1);
	for (j=0; j<np1; j++) {     // each Point
		A[IJ(0,j,np1)] = 1.0;
		for (i=0; i<n; i++) {    // each coord
			A[IJ(i+1,j,np1)] = (*this)[j][i];
		}
	}
	// compute the determinant using template function determ:
	double rval = determ(np1, &A[0], expon);

	trc.dprint("returning ",rval,"*10^",expon);
	return rval;
}

int
Simplex::
orient (int dir) {
// Permute the points in this Simplex so that the orientation,
// given by det[Z_0 ... Z_n] where Z_i = [1 P_i]^t and P_i is
// the ith Point
// Returns: the number of permutations (0 or 1)
	int rval = 0;
	int expon;
	double det = orientation(expon);
	if (dir < 0 && det > 0.0) {
		permute();
		rval++;
	} else if (dir > 0 && det < 0.0) {
		permute();
		rval++;
	}
	return rval;
}


#ifdef NEVER // unused
double
Face::
orientation(int& expon) const {
// Orientation of a boundary Face: n Points forming
// an (n-1)-dimensional Simplex.
// Returns the determinant of the (n,n) matrix whos
// columns are the elements of each Point in this Face
// except for the coord which is constant across
// the Face. In addition, each column has a 1 in the
// first row
	Trace trc(2,"Face::orientation");
	size_t n = this->size();
	size_t i, j;
	size_t iconst{n+1};

	// determine which coord is constant
	for (i=0; i<n; i++) {    // each coord
		for (j=1; j<n; j++) {     // each Point
			assert(n == (*this)[j].size());
			if (!is_equal((*this)[j][i], (*this)[j-1][i], 8)) {
				break;
			}
		}
		if (j == n) {
			assert(iconst == n+1);
			iconst = i;
		}
	}
	trc.dprint("coord ",iconst," is constant");

	vector<double> A(n*n);
	for (j=0; j<n; j++) {     // each Point
		A[IJ(0,j,n)] = 1.0;
		size_t k{1};
		for (i=0; i<n; i++) {   // each coord except iconst
			if (i != iconst)
				A[IJ(k++,j,n)] = (*this)[j][i];
		}
	}
	trc.dprintm(n,n,n,&A[0],"computing determinant of:");
	double rval = determ(n, &A[0], expon);
	trc.dprint("returning ",rval,"*10^",expon);
	return rval;
}

int
Face::
orient (int dir) {
// put this Face-s n Points into a matrix A
	Trace trc(1,"Face::orient");
	int rval = 0;
	int expon;
	trc.dprint("orienting ",*this);

	double det = this->orientation(expon);

	trc.dprint("dir ",dir,", det ",det);
	if (dir < 0 && det > 0.0) {
		permute();
		rval++;
	} else if (dir > 0 && det < 0.0) {
		permute();
		rval++;
	}
	return rval;
}
#endif // NEVER : unused

bool
Face::
operator==(const Face& f) const {
// two Faces are equal if they have the same points, not
// necessarily in the same order
	assert(size() == f.size());
	assert(operator[](0).size() == f[0].size());
	size_t i, j;
	for (i=0; i<size(); i++) {
		bool found = false;
		for (j=0; j<f.size(); j++) {
			if (operator[](i) == f[j]) {
				found = true;
				break;
			}
		}
		if (!found) {
			return false;
		}
	}
	return true;
}

void
Face::
permute() {
// do a single permutation of the points
	Trace trc(2,"Face::permute");
	trc.dprintv(*this,"input face");
	Point tmp = operator[](0);
	operator[](0) = operator[](1);
	operator[](1) = tmp;
	trc.dprintv(*this,"output face");
#ifdef NEVER
	size_t n = operator[](0).size();
	for (size_t i=0; i<n; i++) {
		int tmp = operator[](0)[i];
		operator[](0)[i] = operator[](1)[i];
		operator[](1)[i] = tmp;
	}
#endif
}

string
Face::
toString() const {
	ostringstream os;
	size_t i, j;
	size_t n = size();
	for (i=0; i<n; i++) {
		for (j=0; j<n; j++) {
			os << "  " << operator[](j)[i];
		}
		os << endl;
	}
	return os.str();
}

void
Face::
plotgnu(int faceno) const {
// Write a Face to a gnuplot plot file named "face"faceno
// A Face has n n-dimensional Points
	Trace trc(1,"Face::plotgnu");
	size_t n = operator[](0).size();
	ostringstream os;
	assert(size() == n);

	if (n != 3) {
		trc.dprint("quick return: only treat 3D");
		return;
	}
	std::ios_base::openmode om = std::ios::trunc;

	string file = vastr("face",faceno);
	ofstream path(file, om);
	if (!path) {
		string err = vastr("cannot open ",file);
		trc.dprint("throwing exception: ",err);
		throw runtime_error(err);
	}

	// following stackoverflow.com/questions/22311919/gnuplot-solid-filled-triangle-in-3d-splot-function
	// first point
	path << (*this)[0][0] << " " << (*this)[0][1] << " "
			<< (*this)[0][2] << "\n\n";
	// the rest
	path << (*this)[1][0] << " " << (*this)[1][1] << " "
			<< (*this)[1][2] << "\n";
	path << (*this)[2][0] << " " << (*this)[2][1] << " "
			<< (*this)[2][2] << "\n";
}

void
plotgnu(vector<Face>& faces, int faceno) {
// Write a vector of Faces to a gnuplot plot file named
// "face"faceno
// A Face has n n-dimensional Points
	Trace trc(1,"plotgnu");
	size_t n = faces[0][0].size();
	ostringstream os;

	if (n != 3) {
		trc.dprint("quick return: only treat 3D");
		return;
	}
	std::ios_base::openmode om = std::ios::trunc;

	string file = vastr("face",faceno);
	ofstream path(file, om);
	if (!path) {
		string err = vastr("cannot open ",file);
		trc.dprint("throwing exception: ",err);
		throw runtime_error(err);
	}

	// following stackoverflow.com/questions/22311919/gnuplot-solid-filled-triangle-in-3d-splot-function
	for (size_t i=0; i<faces.size(); i++) {
		// first point
		path << faces[i][0][0] << " " << faces[i][0][1] << " "
				<< faces[i][0][2] << "\n\n";
		// the rest
		path << faces[i][1][0] << " " << faces[i][1][1] << " "
				<< faces[i][1][2] << "\n";
		path << faces[i][2][0] << " " << faces[i][2][1] << " "
				<< faces[i][2][2] << "\n\n";
	}
}

void
Face::
makegnucmd(int nbf) {
// make a gnuplot command file to plot "nbf" Faces
	string file("faces");
	std::ios_base::openmode om = std::ios::trunc;
	ofstream path(file, om);
	path << "set pm3d ftriangles\n";
	path << "set key off\n";
	path << "splot ";
	for (int i=0; i<nbf-1; i++)
		path << "\"face" << i << "\" with pm3d, ";
	path << "\"face" << nbf-1 << "\" with pm3d\n";
	path << "pause -1 \"press q to quit\"\n";
}

void
Face::
plot(int faceno, string const& file, bool append) const {
// Write a Face to an Apf plot file
// A Face has n n-dimensional Points
	size_t n = operator[](0).size();
	size_t i, j;
	ostringstream os;
	assert(size() == n);

	vector<EsaRun2> runs;
	// n parameters named xi
	vector<EsaPar> par;
	for (i=0; i<n; i++) {
		os.str("");
		os << "x" << i+1;
		par.push_back(EsaPar(os.str(), os.str()));
	}
	// original points: n+1 (n-1)-simplices, n(n+1)/2 lines
	EsaRun2 ptsrun("pts", "Nuclei");
	for (j=0; j<n; j++) {
		ptsrun.par.push_back(par[j]);
	}

	for (i=0; i<n; i++) {
		for (j=0; j<n; j++) {
			// ptsrun.par[j].values.clear();
			ptsrun.par[j].solns.push_back(operator[](i)[j]);
			// add a line back to the start point if n>2
			if (i == n-1) {
				if (n > 2)
					ptsrun.par[j].solns.push_back(operator[](0)[j]);
			} else {
				ptsrun.par[j].solns.push_back(operator[](i+1)[j]);
			}
		}
	}
	os.str("");
	os << "face" << faceno;
	ptsrun.id = os.str();
	runs.push_back(ptsrun);
	writeEsa2(file, runs, append);
}

void
plotFaces (vector<Face> const& faces, string const& file) {
	bool append = false;

	for (size_t i=0; i<faces.size(); i++) {
		faces[i].plot(i, file, append);
		append = true;
	}
}

Simplex::
Simplex (vector<Point> const& points) {
// Create a Simplex: n+1 n-dimensional Points
	Trace trc(2,"Simplex constructor");
	size_t n = points[0].size();
	size_t npts = points.size();
	assert (npts == n+1);

	for (size_t i=0; i<npts; i++) {
		trc.dprintv(points[i],"point ",i);
		assert(points[i].size() == n);
		this->push_back(points[i]);
	}
}

vector<Face>
Simplex::
faces() const {
// A simplex has n+1 Points, a Face has n of these points,
// so there are n+1 faces (n+1 things combined n at a time)
// Returns Faces ordered as in Stenger eq. 2.5, with odd faces
// permuted instead of the minus sign):
//   b[X_0...X_n] = \sum_{i=0}^n (-1)^i [X_0...X_{i-1} X_{i+1}...X_n]
	vector<Face> rval;
	size_t n = (*this)[0].size();
	size_t npts = this->size();
	size_t i, j;

	if (npts != n+1) {
		ostringstream os;
		os << "Simplex::faces: Simplex has " << npts << " Points, should have "
			<<  n+1;
		throw runtime_error(os.str());
	}

	for (i=0; i<npts; i++) {
		Face face;
		for (j=0; j<npts; j++) {
			if (j != i)
				face.push_back((*this)[j]);
		}
		// permute the odd faces
		if (i%2 == 1)
			face.permute();
		rval.push_back(face);
	}
	assert(rval.size() == npts);
	return rval;
}

vector<Face>
Simplex::
cube_bndry() {
// get the min and max indices for this simplex,
// then find all faces which have all points
// for some coord at the min or max
	Trace trc(2,"cube_bndry");
	vector<Face> all = faces();
	size_t n = this->size() - 1;
	int big{std::numeric_limits<int>::max()};
	vector<int> idxmin(n, big);
	vector<int> idxmax(n, -big);
	vector<Face> rval;
	size_t i, j, k;

	for (i=0; i<all.size(); i++) {         // each face
		for (j=0; j<all[i].size(); j++) {   // each point
			for (k=0; k<n; k++) {            // coord index
				if (all[i][j][k] < idxmin[k]) {
					idxmin[k] = all[i][j][k];
				}
				if (all[i][j][k] > idxmax[k]) {
					idxmax[k] = all[i][j][k];
				}
			}
		}
	}
	trc.dprintv(idxmin,"min indices");
	trc.dprintv(idxmax,"max indices");
	// now find all faces having a coord at min or max
	for (i=0; i<all.size(); i++) {   // face
		for (k=0; k<n; k++) {         // coord
			bool lower = false;
			bool upper = false;
			if (all[i][0][k] == idxmin[k]) {
				lower = true;
			} else if (all[i][0][k] == idxmax[k]) {
				upper = true;
			} else {
				continue;   // next coord
			}
			for (j=1; j<all[i].size(); j++) {  // point
				if ((lower && all[i][j][k] != idxmin[k]) ||
					(upper && all[i][j][k] != idxmax[k])) {
					break;
				}
			}
			if (j == all[i].size()) {
				rval.push_back(all[i]);
				trc.dprint("added face ",i);
			}
		}
	}
	// there should be 2n boundary faces
	// if (rval.size() != 2*n) {
	// 	cerr << "wrong number of boundary faces for n=" << n
	// 		<< ": " << rval.size();
	// }
	trc.dprint("returning ",rval.size()," boundary faces");
	return rval;
}

void
Simplex::
permute() {
// do a single permutation of the points
	Point tmp = (*this)[0];
	(*this)[0] = (*this)[1];
	(*this)[1] = tmp;
}

void
Simplex::
plot(int simpno, string const& file, bool append) const {
// Write a Simplex to an ESA plot file
// XXX just get this Simplex's Faces and call plotfaces?
	size_t n = (*this)[0].size();
	size_t i, j, k;
	ostringstream os;
	// 
	vector<EsaRun2> runs;
	// n parameter named xi
	vector<EsaPar> par;
	for (i=0; i<n; i++) {
		os.str("");
		os << "x" << i+1;
		par.push_back(EsaPar(os.str(), os.str()));
	}
	// original points: n+1 (n-1)-simplices, n(n+1)/2 lines
	EsaRun2 ptsrun("pts", "Nuclei");
	for (j=0; j<n; j++) {
		ptsrun.par.push_back(par[j]);
	}
	for (k=0; k<n; k++) {
		ptsrun.par[k].solns.clear();
	}
	// plot the points in vis "penlift" mode: pairs of points
	// define line segments
	for (i=0; i<n; i++) {
		for (j=i+1; j<n+1; j++) {
			for (k=0; k<n; k++) {
				ptsrun.par[k].solns.push_back((*this)[i][k]);
				ptsrun.par[k].solns.push_back((*this)[j][k]);
			}
		}
	}
	os.str("");
	os << "simplex" << simpno;
	ptsrun.id = os.str();
	runs.push_back(ptsrun);
	writeEsa2(file, runs, append);
}

void
plotSimplices (vector<Simplex> simplices, string const& file) {
	Trace trc(2,"plotSimplices");
	bool append = false;

	trc.dprint(simplices.size()," Simplex, to file ",file);
	for (size_t i=0; i<simplices.size(); i++) {
		simplices[i].plot(i, file, append);
		append = true;
	}
}

static
double
factorial (int n) {
// returns a double to avoid overflowing the
// range of 32-bit int's (2.147*10^9)
	double rval = 1.0;
	for (int i=n; i>1; i--)
		rval *= i;
	return rval;
}


static
int
sgn (double a) {
	if (a > 0.0)
		return 1;
	else if (a == 0.0)
		return 0;
	return -1;
}

bool
bigendian () {
  /* Are we little or big endian?  From Harbison&Steele.  */
  union {
    long l;
    char c[sizeof (long)];
  } u{1};
  return (u.c[sizeof (long) - 1] == 1);
}

static bool Be{bigendian()};


class Fx {
public:
	vector<vector<double> > fs;
	vector<vector<double> > jacs;
	// vector<int> sgndets;
	int nb;  // number of bytes from each coord

	// constructor
	Fx(int nbytes=4) : nb(nbytes) {}

	// fless is the comparator for fmap
	class fless {
	public:
		bool operator()(vector<char> const& a, vector<char> const& b) const {
			Trace trc(3,"fless");
			auto p = a.begin();
			auto q = b.begin();

			// find the first char where a, b are different
			while (p != a.end() && q != b.end() && *p == *q) {
				p++;
				q++;
			}
			// if a and b are the same => not enough sig figs?
			if (p == a.end())
				return false;
			return *p < *q;
		}
	};

	// fmap maps a vector<char> key to a size_t index,
	// an indexinto arrays fs and jacs
	map<vector<char>,size_t,fless> fmap;

	vector<char> makekey (vector<double> x) {
		Trace trc(3,"makekey");

		size_t n = x.size();
		union {
			double x;
			char c[8];
		} u;
		vector<char> rval;

		// take the most significant nb bytes,
		// pack them into a vector<char>
		if (Be) {
			throw runtime_error("big-endian arch not implemented");
		} else {
			for (size_t i=0; i<n; i++) {
				u.x = x[i];
				for (int j=0; j<nb; j++) {
					rval.push_back(u.c[7-j]);
				}
			}
		}
		trc.dprintv(rval,"returning");
		return rval;   // need move
	}

	bool find (const vector<double>& x, vector<double>& f, vector<double>& jac) {
		// search member xs for x, set "f" to corresponding "fs"
		Trace trc(3,"Fx::find");
		trc.dprintv(x,"x");
		bool found{false};

		auto k = fmap.find(makekey(x));
		if (k != fmap.end()) {
			size_t idx = k->second;
			f = fs[idx];
			jac = jacs[idx];
			// sgndet = sgndets[idx];
			found = true;
		}
		trc.dprint("returning ",found?"true":"false");
		return found;
	}

	void add (const vector<double>& x,
			const vector<double>& f, vector<double>& jac) {
		vector<char> key = makekey(x);
		if (fmap.find(key) != fmap.end()) {
			throw runtime_error(vastr("grid too fine: attempt to add existing point ",x));
		}
		size_t idx = fs.size();
		fmap[makekey(x)] = idx;
		fs.push_back(f);
		jacs.push_back(jac);
		// sgndets.push_back(sd);
	}
};

void
callfcn(Point& pt, int (*fcn)(vector<double> const& x,
			vector<double>& f, vector<double>& jac)) {
// transition to this interface (Point,fcn)
	callfcn(pt, pt.f, pt.jac, fcn);
}

double FnormMin{std::numeric_limits<double>::max()};
double FnormMax{-std::numeric_limits<double>::max()};

void
callfcn(const vector<double>& x, vector<double>& f,
		vector<double>& jac, Userfcn fcn) {
	Trace trc(1,"callfcn");
	// f need not be the correct size - resize if necessary
	size_t n = x.size();
	assert(n > 0);
	if (f.size() != n) {
		f = vector<double>(n, 0.0);
	}
	if (jac.size() != n*n) {
		jac = vector<double>(n*n, 0.0);
	}
	// previously computed f(x) are saved in fx
	static Fx fx;
	if (fx.find(x, f, jac)) {
		trc.dprintv(f,"found existing f(x): ");
		trc.dprintv(x,"at x: ");
		// trc.dprint("sgndet = ",sgndet);
		return;
	}

	// call the user's function
	fcn(x, f, jac);

	Nfcn++;

	// check for (near) zero f
	double fnorm = blas_snrm2(n, &f[0], 1);
	FnormMin = std::min(fnorm, FnormMin);
	FnormMax = std::max(fnorm, FnormMax);
	if (is_equal(fnorm, 0.0, 12)) {
		cerr << "zero function encountered at (";
		cerr << x[0];
		for (size_t i=1; i<n; i++)
			cerr << ", " << x[i];
		cerr << ")\n";
	}
	// save for future calls
	fx.add(x, f, jac);
	trc.dprintv(f,"add f(x)");
	trc.dprint("now have ",fx.fs.size());

	// check the Jacobian
	if (Check_jacobian) {
		size_t maxi, maxj;
		double maxre = check_jac(x, fcn, maxi, maxj);
		if (maxre > 1.0e-6) {
			cerr << vastr("max rel error in Jacobian at (",x,
					"): ",maxre," in the (", maxi,',',maxj,") term") << endl;
		}
	}

	return;
}

void
picardfcn(const vector<double>& px, vector<double>& pf, Userfcn fcn) {
	Trace trc(2,"picardfcn");
	size_t np1 = px.size();
	size_t n = np1 - 1;
	assert(n > 0);
	trc.dprint("n = ",n," pf.size = ",pf.size());
	assert(pf.size() == np1);
	// callfcn computes the n-dimensional f(x)
	vector<double> x(n, 0.0);
	vector<double> f(n, 0.0);
	vector<double> jac(n*n, 0.0);
	for (size_t i=0; i<n; i++) {
		x[i] = px[i];
	}
	// int sgndet;
	callfcn(x, f, jac, fcn);
	for (size_t i=0; i<n; i++) {
		pf[i] = f[i];
	}
	// last eqn in (n+1)-dimensional Picard problem: x[n]*det(jac)
	int expon;
	double det = determ(n, &jac[0], expon);
	pf[n] = px[n]*det;
	trc.dprint("f[n] = det(J))*x[n] = ",det,"*",px[n]);

	return;
}

vector<Face>
Cube::
faces() {
// returns all the n+1 Faces of all the
// n! Simplex's in this Cube
	vector<Face> rval;
	for (size_t i=0; i<this->size(); i++) {
		vector<Face> fi = (*this)[i].faces();
		for (size_t j=0; j<fi.size(); j++) {
			rval.push_back(fi[j]);
		}
	}
	return rval;
}

vector<Point>
Cube::
limits() {
// get the min and max coordinates for this Cube
// Returns 2 Points: the minimums and the maximums
	Trace trc(2,"Cube::limits");
	size_t n = this->order();
	double big{std::numeric_limits<double>::max()};
	vector<double> xmin(n, big);
	vector<double> xmax(n, -big);

	for (auto sp : *this) {              // each Simplex
		for (auto& pt : sp) {        // each Point
			assert(pt.size() == n);
			for (size_t k=0; k<n; k++) {   // coord
				if (pt[k] < xmin[k]) {
					xmin[k] = pt[k];
				}
				if (pt[k] > xmax[k]) {
					xmax[k] = pt[k];
				}
			}
		}
	}

	vector<Point> rval;
	rval.push_back(Point(xmin));
	rval.push_back(Point(xmax));
	trc.dprintv(rval,"returning min/max ");
	return rval;
}

vector<Face>
Cube::
boundary() {
// gather all faces which make up the boundary of this Cube
// get the min and max indices for this Cube,
// then find all faces which have all points
// for some coord at the min or max
	Trace trc(2,"Cube::boundary");
	vector<Face> all = faces();
	vector<Point> lim = this->limits();
	size_t n = this->order();
	vector<double> xmin(n);
	vector<double> xmax(n);
	vector<Face> rval;
	size_t i, j, k;

	for (i=0; i<n; i++) {
		xmin[i] = lim[0][i];
		xmax[i] = lim[1][i];
	}
	trc.dprintv(xmin,"min coords");
	trc.dprintv(xmax,"max coords");

	// now find all faces having a coord at min or max
	for (i=0; i<all.size(); i++) {   // face
		for (k=0; k<n; k++) {         // coord
			bool lower = false;
			bool upper = false;
			if (all[i][0][k] == xmin[k]) {
				lower = true;
			} else if (all[i][0][k] == xmax[k]) {
				upper = true;
			} else {
				continue;   // next coord
			}
			for (j=1; j<all[i].size(); j++) {  // point
				if ((lower && all[i][j][k] != xmin[k]) ||
					(upper && all[i][j][k] != xmax[k])) {
					break;
				}
			}
			if (j == all[i].size()) {
				rval.push_back(all[i]);
				trc.dprint("added face ",i);
			}
		}
	}
	// there should be 2n! boundary faces
	if (rval.size() != 2*factorial(n)) {
		ostringstream os;
	 	os << "wrong number of boundary faces for n=" << n
	 		<< ": " << rval.size() << ", should be "
			<< 2*factorial(n) << endl;
		throw runtime_error(os.str());
	}
	trc.dprint("returning ",rval.size()," boundary faces");
	return rval;
}

int
bisect(vector<vector<double> > pts, vector<double> a, vector<double> b) {
// add a mid-point to each coord between a and b
	Trace trc(1,"bisect");
	int nadded{0};
	size_t n = pts.size();
	size_t i, j;
	for (i=0; i<n; i++) {
		double ai = a[i];
		double bi = b[i];
		double mid = (ai + bi)/2.0;
		size_t m = pts[i].size();
		for (j=1; j<m; j++) {
			if (is_equal(pts[i][j], mid, 8)) {
				trc.dprint(mid," is already in coord ",i);
				break;
			}
			if (is_greaterthan(pts[i][j], mid, 8)) {
				auto q = pts[i].begin() + j;
				pts[i].insert(q, mid);
				trc.dprintv(pts[i], "added ",mid," to coord ",i);
				nadded++;
				break;
			}
		}
	}
	trc.dprint("added ",nadded," midpoints");
	return nadded;
}

#ifdef NEVER // unused
double
Cube::
td(vector<vector<double> > pts,
	int (*fcn)(vector<double> const&, vector<double>&, int&), int& nadded) {
// compute the topological degree for a Cube by
// summing
//    deg(f,P_n) = \frac{1}{2^n n!} \sum_k=1^m \delta sgn \Matrix{B}^k
// that is, for each of the m = 2n boundary faces compute
//    \Matrix{B}^k = \[ \Vector{f}(x_1^k) ...\]
// and sum the determinant of the sgn(\Matrix{B}^k)
// A Cube is the triangulation of one interval in each
// coordinate consisting of n! simplices
	Trace trc(2,"Cube::td");
	vector<Face> faces = this->boundary();
	size_t n = (*this)[0][0].size();
	size_t i, j, k;

	// sanity check
	assert(n == faces[0][0].size());

	nadded = 0;
	int sum{0};
	vector<int> Bk(n*n);
	for (k=0; k<faces.size(); k++) {  // Face k
		assert(faces[k].size() == n);
		for (j=0; j<n; j++) {          // Point j
			// allocate f and evaluate
			if (faces[k][j].f.empty()) {
				faces[k][j].f = vector<double>(n,0.0);
				callfcn(faces[k][j], faces[k][j].f, faces[k][j].jac, fcn);
			}
			for (i=0; i<n; i++) {    // coord i
				Bk[IJ(i,j,n)] = sgn(faces[k][j].f[i]);
			}
			// check angle between this and previous points
			if (j > 0) {
				double angle = getangle(faces[k][j-1].f, faces[k][j].f);
				if (angle > 90.0) {
					nadded += addpts(faces[k][j-1], faces[k][j],
							faces[k][j-1].f, faces[k][j].f, fcn, pts, 0);
				}
			}
		}
		trc.dprintm(n,n,n,&Bk[0],"computing determinant of:");
		// sum += intDet(n, &Bk[0]);
		int expon;
		sum += determ(n, &Bk[0], expon);
		trc.dprint("face ",k,": summation now = ",sum);
	}

	trc.dprint("added ",nadded," Points");

	// the td is \fac{1}{2^n n!} times sum
	int denom = pow(2, n)*factorial(n);
	if (sum % denom != 0) {
		trc.dprint("grid too course: sum = ",sum,", denom = ",denom);
		vector<Point> lim = this->limits();
		// callfcn(lim[0], fcn);
		// callfcn(lim[1], fcn);
		// int na = addpts(lim[0], lim[1], lim[0].f, lim[1].f, fcn, pts, 0);
		// trc.dprint("added ",na," Points");
		bisect(pts,  lim[0], lim[1]);
	}

	this->topdeg = (double)sum/(double)denom;
	trc.dprint("returning topological degree ",this->topdeg);
	return this->topdeg;
}  // Cube::td

double
Cube::
tdpicard(vector<double> const& xmin, vector<double> const& xmax,
		int (*fcn)(vector<double> const& x, vector<double>& f, int& sd),
		vector<int> nintvl) {
// compute the topological degree for a n-Cube which is
// one face of an (n+1)-dimensional Picard problem.
// by summing
//    deg(f,P_n) = \frac{1}{2^n n!} \sum_k=1^m \delta sgn \Matrix{B}^k
// that is, for each of the m = 2n boundary faces compute
//    \Matrix{B}^k = \[ \Vector{f}(x_1^k) ...\]
// and sum the determinant of the sgn(\Matrix{B}^k)
// A Cube is the triangulation of one interval in each
// coordinate consisting of n! simplices
	Trace trc(2,"Cube::tdpicard");
	size_t n = nintvl.size();
	size_t np1 = n + 1;
	size_t i, j, k;

	// a Cube "is a" vector<Simplex>
	int sum{0};
	vector<int> Bk(np1*np1);
	for (k=0; k<this->size(); k++) {         // Simplex k
		Simplex& sp = (*this)[k];
		assert(sp.size() == np1);
		for (j=0; j<np1; j++) {   // Point j
			if (sp->pts[j].x.empty()) {
				sp->pts[j].x = vector<double>(n, 0.0);
				sp->pts[j].f = vector<double>(n, 0.0);
				for (i=0; i<n; i++) {       // coord i
					double dx = xmax[i] - xmin[i];
					sp->pts[j].x[i] = xmin[i] + sp->pts[j][i]*dx/(double)nintvl[i];
				}
				callfcn(sp->pts[j].x, sp->pts[j].f, sp->pts[j].jac, fcn);
			}
			assert(sp->pts[j].f.size() == n);
			for (i=0; i<n; i++) {    // coord i
				Bk[IJ(i,j,np1)] = sgn(sp->pts[j].f[i]);
			}
			Bk[IJ(n,j,np1)] = sp->pts[j].sgndet;
		}
		// if (!ignsign) {
		// 	int nsign = signChanges(&Bk[0],np1);
		// 	if (nsign > 1) {
		// 		cerr << "Face " << k << " has " << nsign
		// 			<< " sign changes\n";
				// return false;
		// 	}
		// }
		trc.dprintm(np1,np1,np1,&Bk[0],"computing determinant of:");
		// sum += intDet(np1, &Bk[0]);
		int expon;
		sum += determ(np1, &Bk[0], expon);
		trc.dprint("Simplex ",k,": summation now = ",sum);
	}

	// the td is \fac{1}{2^n n!} times sum
	// int denom = pow(2, n)*factorial(n);
	// if (sum % denom != 0) {
	// 	trc.dprint("grid too course: sum = ",sum,", denom = ",denom);
	// }

	// this->topdeg = (double)sum/(double)denom;
	// trc.dprint("returning ",this->topdeg);
	// return this->topdeg;
	trc.dprint("returning sum ",sum);
	return sum;
}
#endif // NEVER : unused

void
Cube::
plot (string const& file) {
	Trace trc(2,"Cube::plot");
	plotSimplices (*this, file);
}

int
checkAngles (vector<vector<double> >& pts, Userfcn fcn) {
// for each coordinate (i):
//   - for each combination of values ot all other
//     coordinates:
//     - hold all other coordinates constant at this
//       combination of values
//     - for each value (k) of coordinate i compute
//       the angle between f(k) and f(k-1)
//     - if the angle is greater than maxangle add
//       points between x(k-1) and x(k)
	Trace trc(1,"checkAngles");
	int rval{0};
	size_t i, j, k;
	size_t n = pts.size();
	size_t nm1 = n - 1;
	double maxangle{90};

	for (i=0; i<n; i++) {
		int naddi{0};
		vector<int> npts;
		for (j=0; j<n; j++) {
			if (j != i)
				npts.push_back(pts[j].size());
		}
		trc.dprint("working on coord ",i,", npts = ",npts);
		int idx{0};
		vector<int> ind(nm1);
		vector<double> xk(n);
		vector<double> xkm1(n);
		vector<double> fk(n);
		vector<double> fkm1(n);
		vector<double> jac(n*n);
		while (indices(idx++, npts, ind)) {
			// set the value of all coordinates except i
			for (k=j=0; j<n; j++) {
				if (j != i) {
					xk[j] = pts[j][ind[k++]];
					xkm1[j] = xk[j];
				}
			}
			for (k=0; k<pts[i].size(); k++) {		// each value of coord i
				xk[i] = pts[i][k];
				callfcn(xk, fk, jac, fcn);
				// for (j=0; j<n; j++)
				// 	fk[j] = (double)sgn(fk[j]);
				if (k > 0) {
					double angle = getangle(fk, fkm1);
					if (angle > maxangle) {
						int nadd =
							addpts(xk, xkm1, fk, fkm1, fcn, pts, 0);
						rval += nadd;
						naddi += nadd;
						trc.dprint("function angle ",angle, " > ",maxangle,", added ",nadd," points");
					}
				}
				fkm1 = fk;
				xkm1[i] = xk[i];
			}
		}
		trc.dprint("finished working on coord ",i,", added ",naddi);
	}
	trc.dprint("returning ",rval, " added points, new # pts: ",coordpts(pts));
	return rval;
}

static
int
topdeg (vector<vector<double> >& pts, Userfcn fcn) {
// Compute the topological degree without the Picard
// extension of an nD region defined by points in n-space
// "pts". pts[i][j] is the jth point on the ith
// coordinate axis.
// Start with "pts", add points as needed to keep the
// angle between adjacent function values below 90 degrees,
// then re-triangulate with the new points. Repeat this
// process until no new points are added.
	Trace trc(1,"topdeg");
	int rval{0};
	size_t n = pts.size();
	size_t i, j, k;
	size_t idx;

	trc.dprintv(coordpts(pts),"begining # pts: ");

	// get xmin, xmax from pts
	vector<double> xmin(n);
	vector<double> xmax(n);
	for (i=0; i<n; i++) {
		xmin[i] = pts[i][0];
		xmax[i] = pts[i][pts[i].size()-1];
	}

	int sum;
	bool newgrid;
	vector<Cube> cubes;
	do {
		newgrid = false;
		sum = 0;
		size_t nbf{0};  // number of boundary faces
		Nperm = 0;      // number of permutations
		// triangulate the nD region with pts
		// Add to these points if necessary (see getangle())
		cubes = triangulate (pts);

		// compute the nD TD
		// with extra points inserted as needed to keep
		// the angle between function vectors down
		trc.dprint(cubes.size()," Cube");

		for (size_t ic=0; ic<cubes.size(); ic++) {
			for (idx=0; idx<cubes[ic].size(); idx++) {
				Simplex& sp = cubes[ic][idx];
				// get the boundary faces for this Simplex
				vector<Face> faces =
					sp.boundaryFaces(xmin, xmax);
				nbf += faces.size();

				// for each boundary Face: eval the function at
				// each of the n points, put sgn into Bk
				vector<double> Ak(n*n);
				vector<int> Bk(n*n);
				for (k=0; k<faces.size(); k++) {  // Face k
					assert(faces[k].size() == n);
					for (j=0; j<n; j++) {          // Point j
						if (faces[k][j].f.empty()) {
							faces[k][j].f = vector<double>(n, 0.0);
							trc.dprint("calling fcn for Cube ",ic, " Face ",k,", Point ",j,": ",faces[k][j]);
							callfcn(faces[k][j], faces[k][j].f, faces[k][j].jac, fcn);
						}
						for (i=0; i<n; i++) {
							Ak[IJ(i,j,n)] = faces[k][j].f[i];
							Bk[IJ(i,j,n)] = sgn(faces[k][j].f[i]);
						}
						// check angle between this and the previous point
						if (j > 0) {
							double angle = getangle(faces[k][j-1].f, faces[k][j].f);
							if (angle > 90.0) {
								addpts(faces[k][j-1], faces[k][j], faces[k][j-1].f,
										faces[k][j].f, fcn, pts, 0);
								newgrid = true;
							}
						}
					}
					trc.dprintm(n,n,n,&Bk[0],"computing det of:");
					int expon;
					int detk = determ(n, &Bk[0], expon);
					sum += detk;
					trc.dprint(n,"D cube/simplex/face ",ic,'/',idx,'/',k, " det ",detk," summation now = ",sum);
				}
			}
		}
		trc.dprint("processed ",nbf," boundary faces, ", Nperm," permutations");

		trc.dprint("sum of ",n,"D td:",sum);
		int denom = pow(2, n)*factorial(n);
		if (sum % denom != 0) {
			trc.dprint("grid too course: sum = ", sum,", denom = ",denom);
		}
		rval = sum/denom;
		trc.dprint("estimated TD: ",rval);

			// Compute the TD for each Cube in the final grid
#ifdef NEVER
		if (!newgrid) {
			trc.dprint("final grid has ",cubes.size(), " Cubes, computing TD for each Cube");
			int npts = countpts(pts);
			int nnewpts{0};
			for (i=0; i<cubes.size(); i++) {
				int nadded{0};
				trc.dprint("working on Cube ",i,"/",cubes.size());
				int td = cubes[i].td(pts, fcn, nadded);	
				nnewpts += nadded;
				trc.dprint("Cube ",i," td = ",td," added ",nadded," Points");
				if (td > 0) {
					vector<Point> lim = cubes[i].limits();
					cout << "possible root between " << lim[0]
						<< " and " << lim[1] << endl;
				} else if (nadded > 0) {
					newgrid = true;
				}
			}
		}
#endif // NEVER
		if (newgrid) {
			trc.dprint("added points: re-triangulate");
			trc.dprintv(coordpts(pts),"new # pts: ");
		}
	} while (newgrid);

	trc.dprint("returning TD = ",rval);
	return rval;
}  // topdeg

int
picard (vector<double> const& xmin, vector<double> const& xmax,
	Userfcn fcn, int inint, int& nDtd) {
// Given
//   xmin   n-vector of minimums for each independent variable
//   xmax   n-vector of maximums for each independent variable
// Returns  (int) TD with the Picard extension, also:
//   nDtd   (int) TD without the Picard extension

	Trace trc(1,"picard");
	size_t n = xmin.size(); // number of eqns in f(x)
	size_t i, j;
	ostringstream os;
	vector<int> nintvl(n,inint);  // initial intervals

	// if CHECKJAC is in the env call check_jac
	const char* chk = getenv("CHECKJAC");
	if (chk != nullptr)
		Check_jacobian = true;

	trc.dprintv(nintvl,"intervals");
	trc.dprintv(xmin,"minimums");
	trc.dprintv(xmax,"maximums");

	// sanity checks
	if (xmax.size() != n) {
		os << xmax.size() << " upper-limits given, should be " << n;
		throw runtime_error(os.str());
	}
	for (size_t i=0; i<n; i++) {
		if (xmin[i] >= xmax[i]) {
			os << "illegal limits on dimension " << i << "(0b): "
				<< xmin[i] << ':' << xmax[i];
			throw runtime_error(os.str());
		}
	}

	// create a set of evenly-spaced points along
	// each of the n coord axes
	vector<vector<double> > pts(n);
	for (i=0; i<n; i++) {
		double dx = (xmax[i] - xmin[i])/((double)nintvl[i]);
		size_t npt = nintvl[i] + 1;
		for (j=0; j<npt; j++) {
			pts[i].push_back(xmin[i] + j*dx);
		}
	}
	// compute the TD without the Picard extension (nDtd arg)
	// adding points to "pts" as necessary to keep the
	// angle between adjacent function values below 90 deg
	nDtd = topdeg (pts, fcn);
	trc.dprint("TD without Picard: ",nDtd);

	// check angle between function values, add pts if necessary
	int nadd{0};
	// do {
		nadd = checkAngles(pts, fcn);
	// } while (nadd > 0);
		trc.dprint(nadd," points added to ",n,"D grid");

	// compute the TD with the Picard extension
	// add a dimension to pts: the Picard extension
	size_t np1 = n + 1;
	vector<double> paxis(2);
	paxis[0] = -1.0;
	paxis[1] = 1.0;
	pts.push_back(paxis);

	vector<double> pxmin{xmin};
	pxmin.push_back(-1.0);
	vector<double> pxmax{xmax};
	pxmax.push_back(1.0);

	// triangulate the n+1 space
	vector<Cube> pcubes = triangulate(pts);

	// compute the determinant of Bk for the boundary Faces
	// of each Simplex
	int psum{0};
	size_t nbf{0};  // total number of boundary faces
	int faceno{0};
	size_t k;
	nbf = 2*np1;
	vector<int> Bk(np1*np1);
	vector<int> sumsl(np1, 0);
	vector<int> sumsu(np1,0);
	for (size_t idx=0; idx<np1; idx++) {     // each coord lower
		trc.dprint("working on coord ",idx," lower");
		// coord idx lower
		vector<Face> faces = regionFaces(pcubes, idx, pxmin[idx]);
		for (i=0; i<faces.size(); i++) {    // each facet of lower
			assert(faces[i].size() == np1);
			for (j=0; j<np1; j++) {         // each Point of facet
				if (faces[i][j].f.empty()) {
					faces[i][j].f = vector<double>(np1, 0.0);
					picardfcn(faces[i][j], faces[i][j].f, fcn);
				}
				for (k=0; k<np1; k++) {
					Bk[IJ(j,k,np1)] = sgn(faces[i][j].f[k]);
				}
			}
			trc.dprintm(np1,np1,np1,&Bk[0],"compute det of:");
			// int sumk = intDet(np1, &Bk[0]);
			int expon;
			int deti = determ(np1, &Bk[0], expon);
			psum += deti;
			sumsl[idx] += deti;
			trc.dprint(np1,"D region face ",idx, " lower face ",i,", det = ",deti, " summation now = ",psum);
		}
		plotgnu(faces, faceno++);

		trc.dprint("working on coord ",idx," upper");
		// coord idx upper
		faces = regionFaces(pcubes, idx, pxmax[idx]);
		for (i=0; i<faces.size(); i++) {
			for (j=0; j<np1; j++) {        // each Point of face
				if (faces[i][j].f.empty()) {
					faces[i][j].f = vector<double>(np1, 0.0);
					picardfcn(faces[i][j], faces[i][j].f, fcn);
				}
				for (k=0; k<np1; k++) {
					Bk[IJ(k,j,np1)] = sgn(faces[i][j].f[k]);
				}
			}
			trc.dprintm(np1,np1,np1,&Bk[0],"compute det of:");
			// int sumk = intDet(np1, &Bk[0]);
			int expon;
			int deti = determ(np1, &Bk[0], expon);
			psum += deti;
			sumsu[idx] += deti;
			trc.dprint(np1,"D region face ",idx, " upper face ",i,", det = ",deti, " summation now = ",psum);
		}
		plotgnu(faces, faceno++);
	}
	trc.dprint("lower sums: ",sumsl);
	trc.dprint("upper sums: ",sumsu);

	Face::makegnucmd(nbf);

	int denom = pow(2, np1)*factorial(np1);
	if (psum % denom != 0) {
		trc.dprint("grid too course: sum = ",psum,", denom = ",denom);
	}

	int rval = psum/denom;

	trc.dprint("topological degree with Picard extension: ",rval);
	trc.dprint("processed ",nbf," Picard boundary faces");
	// compute the number of calls to picardfcn and the
	// final grid size
	int sum{0};
	os.str("");
	for (i=0; i<np1; i++) {
		int prod{1};
		for (j=0; j<np1; j++) {
			if (j != i)
				prod *= pts[j].size()-1;
		}
		sum += prod;
		os << pts[i].size();
		if (i != n)
			os << " x ";
	}
	int ncall = 2*factorial(n)*np1*sum;
	trc.dprint("there should be ",ncall," calls to picardfcn");
	cerr << "final grid: " << os.str() << " points\n";
	cerr << "min, max f(x) = " << FnormMin
			<< ", " << FnormMax << endl;

	return rval;
}

static
bool
indices(int idx, vector<int> const& nitem, vector<int>& rval) {
// given an index (idx) into an m-dimensional array with "nitem"
// items in each coord, set the index of each coord in "rval"
	Trace trc(2,"indices ",idx);
	size_t m = nitem.size();
	size_t i;

	trc.dprintv(nitem,"number of items in each coord:");

	// initialize and check for idx out of range
	int maxidx = 1;
	for (i=0; i<m; i++) {
		rval[i] = 0;
		maxidx *= nitem[i];
	}
	if (idx >= maxidx) {
		trc.dprint("returning false: ",idx," out of the range 0-",maxidx-1);
		return false;
	}

	for (i=0; i<m; i++) {
		if (idx == 0)
			break;
		rval[i] = idx%nitem[i];
		idx /= nitem[i];
	}
	trc.dprintv(rval,"returning");
	return true;
}

vector<Face>
Simplex::
boundaryFaces (vector<double> const& lower,
		vector<double> const& upper) {
// Given n+1 n-dimensional Points, look for combinations of n
// Points in which a coord is at it's lower bound at each
// point, or a coord is at it's upper bound at each point.
// These n Points form a boundary Face, an (n-1)-dimensional
// hypersurface.
	Trace trc(2,"boundaryFaces");
	vector<Face> rval;
	size_t i, j;

	size_t n = (*this)[0].size();
	assert (this->size() == n+1);
	assert(lower.size() == n);
	assert(upper.size() == n);

	trc.dprint("searching Simplex ",*this);
	trc.dprintv(lower,"lower");
	trc.dprintv(upper,"upper");
	if (trc()) {
		int expon;
		double det = this->orientation(expon);
		assert(det > 0.0);
	}

	// check each coord...
	for (i=0; i<n; i++) {
		size_t nll = 0;
		size_t nul = 0;
		// ...at each point...
		for (j=0; j<this->size(); j++) {
			// ...to count the number of points at upper or lower
			if (is_equal((*this)[j][i], lower[i], 8)) {
				nll++;
			} else if (is_equal((*this)[j][i], upper[i],8)) {
				nul++;
			}
		}
		// if a coord is at a limit for n of the n+1 points,
		// we have found a boundary face
		double limit;
		// int dir{1};
		if (nll == n) {
			limit = lower[i];
			// dir = -1;
		} else if (nul == n) {
			limit = upper[i];
		}
		if (nll == n || nul == n) {
			Face face;
			bool needPermute = false;
			// test each point to see if coord i is at the limit
			// if so, add it to the face. It may be necessary to
			// permute the points to keep the orientation correct
			// see stenger1975computing.pdf eqn. 2.4
			for (j=0; j<(*this).size(); j++) {
				if (is_equal((*this)[j][i], limit, 8)) {
					face.push_back((*this)[j]);
				} else {
					// there will be only one point not added,
					// permute if j odd effectively mult by (-1)^j
					if (j%2 != 0)
						needPermute = true;
				}
			}
			if (needPermute)
			   face.permute();
			assert(face.size() == n);
			// check the orientation of this face
			// Nperm += face.orient(dir);
			// add it to the list of Face's to return
			rval.push_back(face);
		}
	}
	trc.dprintv(rval,"returning ",rval.size()," faces");
	return rval;
}

vector<Face>
regionFaces (vector<Cube>& cubes, int idx, double idxval) {
// Given a set of Cubes from the triangulation of an
// nD region, find all Faces ((n-1)D Simplexes) with 
// all n Points having coordinate "idx" = idxval.
// These Faces comprise a face of the nD region.
	Trace trc(2,"regionFaces");
	vector<Face> rval;
	size_t i, j, k;

	// n is the size of the first Point of the first
	// Simplex of the first Cube
	size_t n = cubes[0][0][0].size();
	assert(idx >= 0 && idx < (int)n);

	size_t nfact = factorial(n);
	size_t np1 = n+1;

	// for each Cube check each Simplex for Faces
	// with each Point[idx] == idxval
	for (i=0; i<cubes.size(); i++) {
		// each Cube must have n! Simplex
		assert(cubes[i].size() == nfact);
		// ... each Simplex...
		for (j=0; j<nfact; j++) {
			// each Simplex must have n+1 Points
			assert(cubes[i][j].size() == np1);
			Face face;
			bool needPermute = false;
			// ... count the number of points at idxval
			for (k=0; k<np1; k++) {
				// each Point must have n coord
				assert(cubes[i][j][k].size() == n);
				if (is_equal(cubes[i][j][k][idx], idxval, 8)) {
					face.push_back(cubes[i][j][k]);
				} else if (k%2 != 0) {
					// there will be only one point not added,
					// permute if k odd effectively mult by (-1)^k
					needPermute = true;
				}
			}
			// n Points found with idx = idxval?
			// may need to be permuted
			if (face.size() == n) {
				if (needPermute)
					face.permute();
				rval.push_back(face);
			}
		}
	}
	trc.dprintv(rval,"returning ",rval.size()," faces");
	return rval;
}

static
vector<Cube>
triangulate (vector<vector<double> >& pts) {
// input:  n-vector of vector<double>: pts in each coord
// compute the delaunay triangulation of a set of
// points in n dimensions.
// Returns: prod(nintvl) Cube, each Cube contains n! Simplex
	size_t n = pts.size();
	Trace trc(2,"triangulate",n,"D");
	size_t np1 = n + 1;
	size_t i, j, k;
	vector<double> pt(n);
	vector<int> nintvl(n);
	vector<Cube> rval;

	assert(n > 1);

	// nintvl are the number of intervals in each coord
	for (i=0; i<n; i++) {
		nintvl[i] = pts[i].size() - 1;
		trc.dprintv(pts[i], "coord ",i," pts: ");
	}
	trc.dprintv(nintvl,"number of intervals in each coord: ");

	// compute ncube = prod(nintvl)
	int prod = 1;
	for (i=0; i<n; i++) {
		prod *= nintvl[i];
	}
	trc.dprint("number of Cube: ",prod);
	if (prod == 0) {
		trc.dprint("quick return: too few pts in a coord");
		return rval;
	}

	// if ncube > 1 triangulate only the
	// first n-dimensional cube...
	vector<vector<double> > pts0(n);
	vector<double> x0(n);
	vector<double> x1(n);
	for (i=0; i<n; i++) {
		x0[i] = pts[i][0];
		x1[i] = pts[i][1];
	}
	Cube cube0 = proto(n);

	// ...then dup it
	size_t nsimp = cube0.size();  // # Simplex/Cube
	assert(nsimp == factorial(n));
	assert(cube0[0].size() == np1);
	size_t idx{0};
	vector<int> intvl(n,0);
	while (true) {
		if (!indices(idx++, nintvl, intvl)) {
			trc.dprint("finished duping Cubes");
			break;
		}
		// get the lower and upper values for each coord
		for (i=0; i<n; i++) {
			x0[i] = pts[i][intvl[i]];
			x1[i] = pts[i][intvl[i]+1];
		}
		// set coords to x0 or x1
		vector<Simplex> simplices;
		for (i=0; i<nsimp; i++) {  // each Simplex
			Simplex& simp = cube0[i];
			vector<Point> points;
			for (j=0; j<simp.size(); j++) {  // each pt on Simplex
				vector<double> xs;
				for (k=0; k<n; k++) {    // each coord
					if (simp[j][k] == 0.0)
						xs.push_back(x0[k]);
					else
						xs.push_back(x1[k]);
				}
				points.push_back(Point(xs));
			}
			simplices.push_back(Simplex(points));
		}
	 	// checkSimp(simplices);
		rval.push_back(Cube(simplices));
	}
	// plot the Cubes
	// if (trc() > 0) {
		for (i=0; i<rval.size(); i++) {
			ostringstream os;
			os << n << "DCube" << i;
			rval[i].plot(os.str());
		}
	// 	plotSimplices(simplices, "triangulation");
	// 	checkSimp(simplices);
	// }
	trc.dprint("returning ",rval.size()," Cube");
	return rval;
} // triangulate

ostream&
operator<<(ostream& s, const Point& t) {
	for (const auto& ti : t)
		s << " " << ti;
	return s;
}

ostream&
operator<<(ostream& s, const Face& t) {
	// for (const auto& ti : t)
	for (size_t i=0; i<t.size(); i++) {
		s << " (" << t[i] << ")";
		if (i != t.size()-1)
			s << "->";
	}
	return s;
}

ostream&
operator<<(ostream& s, const Simplex& t) {
// see also: Simplex::print()
	for (size_t i=0; i<t.size(); i++)
		s << " [" << t[i] << "]";
	return s;
}

ostream&
operator<<(ostream& s, const Cube& t) {
	for (size_t i=0; i<t.size(); i++)
		s << t[i] << " ";
	return s;
}

vector<double> CJx;
size_t CJidx;
Userfcn CJfcn;

void
cjfcn(int n, double xi, double* y) {
// compute y(x) for richardson extrapolation
	CJx[CJidx] = xi;
	vector<double> f(n);
	vector<double> jac(n*n);
	CJfcn (CJx, f, jac);
	for (int i=0; i<n; i++)
		y[i] = f[i];
}

double
check_jac(const vector<double>& x, Userfcn fcn, size_t& maxi, size_t& maxj) {
	Trace trc(1,"check_jac");
	double maxre{0.0};

#ifdef NEVER // needs work: convert to ridders
	size_t n = x.size();
	// get the exact jacobian
	vector<double> f(n);
	vector<double> exjac(n*n);
	fcn (x, f, exjac);

	vector<double> appjac(n*n, 0.0);

	CJfcn = fcn;

	double scale{1.0};
	for (CJidx=0; CJidx<n; CJidx++) {
		CJx = x;
		double xj = x[CJidx];
		double dx = std::max(xj*0.01, 0.01);
		double xmin = xj - dx;
		double xmax = xj + dx;
		trc.dprint("working on column ",CJidx," x[",CJidx,"] = ",xj);

		richardson (n,xj,xmin, xmax, scale,
				&appjac[IJ(0,CJidx,n)], cjfcn);
	}
	trc.dprintm(n,n,n,appjac,"approx jacobian");
	int expon;
	trc.dprint("determinant = ",determ(n,&appjac[0],expon));

	trc.dprintm(n,n,n,exjac,"exact jacobian");
	trc.dprint("determinant = ",determ(n,&exjac[0],expon));

	maxi = 0;
	maxj = 0;
	for (size_t i=0; i<n; i++) {
		for (size_t j=0; j<n; j++) {
			double e = abs(appjac[IJ(i,j,n)] - exjac[IJ(i,j,n)]);
			double re = e/std::max(0.1, abs(exjac[IJ(i,j,n)]));
			if (re > maxre) {
				maxi = i;
				maxj = j;
				maxre = re;
			}
		}
	}
	trc.dprint("max rel error in Jacobian: ",maxre," in the (", maxi,',',maxj,") term");

#endif // NEVER // needs work: convert to ridders
	return maxre;
}

#ifdef MAIN
#undef MAIN

#ifdef NEVER
vector<double>
makeEqn(size_t n, vector<double> roots, vector<vector<double> >& ps) {
// create a matrix A and a vector of vectors p such that
//   Ar_i - p_i = 0
// where r_i is a vector with all elements = roots[i]
// So if x = r_i the equation f(x) = Ax - p_i = 0

	// form a random (n,n) matrix
	vector<double> A(n*n, 0.0);
	for (i=0; i<n; i++) {
		for (j=0; j<n; j++) {
			A[IJ(i,j,n)] = xrand();
		}
	}
	// for each root mult A
	for (k=0; k<roots.size(); k++) {
		vector<double> r(n, roots[k]);
		vector<double> p(n, 0.0);
		for (i=0; i<n; i++) {
			for (j=0; j<n; j++) {
				p[i] += A[IJ(i,j,n)]*r[j];
			}
		}
		ps.push_back(p);
	}
	return A;
}
#endif // NEVER

int
linfcn (vector<double> const& x, vector<double>& f, vector<double>& jac) {
// f(i) = x_{i+1} - x[0]  (i=0, n-2)
// f(n-1) = (x_0 - r_0)(x_0 - r_1)...(x_0 - r_{m-1})
//   f'(n-1) = (x_0 - r_1)...(x_0 - r_{m-1}) + 
// m roots at x = [r_i r_i ... r_i]^t   (i=1, m)
	Trace trc(1,"linfcn");
	vector<double> roots;
	// double r = 0.62;
	// double q = 0.12;
	size_t n = x.size();
	size_t i, j, k;

	// set the root(s)
	roots.push_back(0.62);
	// roots.push_back(0.72);

	// f_i = prod(x_i - roots[i])
	double fac{1.0};
	for (i=0; i<n; i++)
		f[i] = 1.0;
	for (j=0; j<roots.size(); j++) {
		for (i=0; i<n; i++) {
			f[i] *= fac*(x[i] - roots[j]);
		}
		fac *= -1.0;
	}
	trc.dprint("f(",x,") = ",f);

	// df/dx = d/dx(x-r_0)*(r1-x)*(x-r2) ...
	//       = (r1-x)*(x-r2) ... + (x-r_0)*(-1)*... +
	for (i=0; i<n; i++) {
		double dfdx = 1.0;
		for (j=0; j<roots.size(); j++) {
			for (k=0; k<roots.size(); k++) {
				if (k != j) {
					dfdx *= (x[i] - roots[k]);
				}
			}
		}
		if (roots.size()%2 == 0)
			dfdx *= -1.0;
		jac[IJ(i,i,n)] = dfdx;
	}
	vector<double> det(2,0.0);
	// determinant(n, &jac[0], &det[0]);
	trc.dprintm(n,n,n,&jac[0],"Jacobian:");
	// int expon;
	// double d =  determ(n, &jac[0], expon);
	// if (det[0] < 0.0)
	// if (d < 0.0)
	// 	sd = -1;
	// else
	// 	sd = 1;

	return 0;
}

void
linear(size_t n, int inint) {
	vector<double> xmin(n, 0.0);
	vector<double> xmax(n, 1.0);
	int nDtd;
	try {
		double td = picard (xmin, xmax, linfcn, inint, nDtd);
		cerr << vastr("topological degree without/with Picard: ",
				nDtd, '/',td) << endl;
	} catch (runtime_error& s) {
		cerr << "caught exception: " << s << endl;
	}
}

int
parfcn (vector<double> const& x, vector<double>& f, vector<double>& jac) {
// f(0) = x_0^2 - x_1
// f[1] = x_0^2 + x_1 - 8
// f' = 2x_0  -1
//      2x_0   1
// det(f') = 4x_0
// sgn(det(f')) = sgn(x_0)
// roots: (-2,4) (2,4)
	Trace trc(1,"parfcn");
	size_t n = x.size();

	assert(n%2 == 0);  // only allow even orders

	trc.dprintv(x,"parfcn x");
	// 2m eqns:
	//   x_i^2 - x_{i+1} = x_i^2 + x_{i+1} - 8 = 0
	for (size_t i=0; i<n; i+=2) {
		f[i] = x[i]*x[i] - x[i+1];
		f[i+1] = x[i]*x[i] + x[i+1] - 8.0;
		for (size_t j=0; j<n; j+=2) {
			jac[IJ(i,j,n)] = 2.0*x[i];
			jac[IJ(i,j+1,n)] = -1.0;
			jac[IJ(i+1,j,n)] = 2.0*x[i];
			jac[IJ(i+1,j+1,n)] = 1.0;
		}
	}
	// double j00, j01, j10, j11;
	// j00 = 2.0*x[0]; j01 = -1.0;
	// j10 = 2.0*x[0]; j11 = 1.0;
	// double det = j00*j11 - j10*j01;
	// if (det < 0.0) sd = -1;
	// else sd = 1;
	// if (x[0] < 0.0)
		// sd = -1;
	// else if (x[0] == 0.0)
	// 	sd = 0;
	// else
	// 	sd = 1;
	return 0;
}

void
parabola(size_t n, int inint) {
	if (n%2 != 0) {
		cerr << "only even-orders allowed for problem=par\n";
		n++;
		cerr << "changing to n = " << n << endl;
	}
	try {
		vector<double> xmin(n);
		vector<double> xmax(n);
		for (size_t i=0; i<n; i+=2) {
			xmin[i] = -3.2;
			xmin[i+1] = 1.0;
			xmax[i] = 3.0;
			xmax[i+1] = 5.0;
		}
		int nDtd;
		double td = picard (xmin, xmax, parfcn, inint, nDtd);
		cerr << vastr("topological degree without/with Picard: ",
				nDtd, '/',td) << endl;
	} catch (runtime_error& s) {
		cerr << "caught exception: " << s << endl;
	}
}

static
void
determinant(size_t n, double* a, double* det) {
// Input:
//   n   order of matrix "a"
//   a   (n,n) real matrix
//   det 2-vector
// Compute the determinant of the (n,n) real matrix "a",
// return it as 2 numbers: d = det[0]*10^det[1]
	Trace trc(2,"determinant (double)");

	trc.dprintm(n,n,n,a,"A");

	// factorize the matrix
	vector<int> pivots(n);
	int info;
	dgetrf_(&n, &n, a, &n, &pivots[0], &info);
	if (info != 0) {
		ostringstream os;
		os << "determinant failed: dgetrf returned " << info;
		trc.dprint("returning 0: ",os.str());
		// throw runtime_error(os.str());
		det[0] = 0.0;
		det[1] = 0.0;
		return;
	}

	// compute the determinant as the product of the diagonals
	// of the upper-triangle. Scale so that 1 < det[0] < 10
	det[0] = 1.0;
	det[1] = 0.0;
	for (size_t i=0; i<n; i++) {
		if (pivots[i] != (int)(i+1))
			det[0] = -det[0];
		det[0] = a[IJ(i,i,n)]*det[0];
		if (abs(det[0]) == 0.0) break;
		while (abs(det[0]) < 1.0) {
			det[0] = 10.0*det[0];
			det[1] -= 1.0;
		}
		while (abs(det[0]) > 10.0) {
			det[0] = det[0]/10.0;
			det[1] += 1.0;
		}
	}
	trc.dprint("returning ",det[0]," 10^",det[1]);
}


void
testdeterm(int n) {
	Trace trc(1,"testdeterm");

	cout << "testing template<> determ"
		<< " with n = " << n << endl;
	// integer
	int expon;
	int i;
	vector<int> Ai(n*n);
	for (i=0; i<n*n; i++)
		Ai[i] = (int)(10.0*flaps::xrand());
	trc.dprintm(n,n,n,Ai,"calling determ with");
	int di = determ(n, &Ai[0], expon);
	trc.dprint("integer determ: ",di,"*10^",expon);
	// repeat with determinant(double)
	vector<double> Ad(n*n);
	for (i=0; i<n*n; i++)
		Ad[i] = (double)Ai[i];
	vector<double> det(2);
	determinant(n, &Ad[0], &det[0]);
	trc.dprint("using determinant(): ",det[0],"*10^",det[1]);

	// double
	for (i=0; i<n*n; i++)
		Ad[i] = 10.0*flaps::xrand();
	trc.dprintm(n,n,n,Ad,"calling determ with");
	double dd = determ(n, &Ad[0], expon);
	trc.dprint("double determ: ",dd,"*10^",expon);

	// repeat with determinant
	determinant(n, &Ad[0], &det[0]);
	trc.dprint("using determinant(): ",det[0],"*10^",det[1]);

	cout << "determ(double): " << dd << "*10^" << expon << endl;
	cout << "determinant(): " << det[0] << "*10^" << det[1] << endl;
}

int
main(int argc, char** argv) {
// Usage:
//   picard [lin|par] n
//    n  -        number of nonlinear equations
//                default: 2
//  for example
//     topdeg par 4
//  will do the parabola problem with n=4
	size_t n = 2;   // default order
	ostringstream os;
	string prob("lin");
	int inint{2};  // default # initial intervals

// topdeg n
	if (argc > 1) {
		prob = (argv[1]);
	}
	if (argc > 2) {
		n = atoi(argv[2]);
	}

	// test determ
	if (prob == "determ") {
		testdeterm(n);
		exit(0);
	}

	cout << "topological degree problem \""
		<< prob << " with n = " << n << endl;
	cout << (Be ? "Big":"Little") << "-endian arch\n";

	// test indices
#ifdef NEVER
	vector<int> ptidx(n);
	size_t idx = 0;
	while (true) {
		if (!indices(idx, npts, ptidx)) {
			break;
		}
		T("index %d:",idx);
		for (size_t j=0; j<n; j++) {
			TN("  %d", ptidx[j]);
		}
		T("\n");
		idx++;
	}
exit(0);
#endif // NEVER

// test delaunay: triangulate the region
#ifdef TEST_TRIANGULATE
	try {
		vector<Cube> cubes = triangulate (nintvl);
	} catch (runtime_error& s) {
		cerr << "caught exception: " << s << endl;
	}
#endif // TEST_TRIANGULATE

// compute the topological degree for one of the
// problems lin or par
	if (prob == "lin") {
		linear(n, inint);
	} else if (prob == "par") {
		parabola(n, inint);
	} else {
		cerr << "unrecognized problem: \"" << prob << endl;
	}

	return 0;
}

#endif // MAIN
