
# controls: generic 2-engine airplane with a nonlinear control law
#
# 51 dof:
#  1-40   airplane flexible
#  41     rudder/geared tab rotation
#  42     left elevator rotation
#  43     right elevator rotation
#  44     left aileron rotation
#  45     right aileron rotation
#  nacelle 3 dof/side 46-51

   import { ${FROOT}/demo/lco.op4, ${FROOT}/demo/lco.uf }

   # use matlab/octave to extract rows of gct, return ztipl, ztipr
   octlab { i=(extract.m,gct.lco), o=(ztipl, ztipr, ldiff, rdiff) }
   display {ztipl, ztipr}

   # experimenting with a range of gains
   parameters{
		clfac(Control-law factor)[0:10.0] = 0.01
		gain(Gain)[0.1:10] = 1,
		phase(Phase)[-120:120]<57.3 Deg/Rad> = 0
	}


   # create a T-matrix with a simple control equation
   pz {
      code=controls{ size=gstif, extra = 2 },
      o = Cont
		plot{indep=freq, ([44,44],[45,45])}
   }
   display{Cont{freq=4, vtas=50} }

   # interpolate gaf matrices wrt rf
   pz {
      i=(gaf0{rf=0}, gaf0008{rf=0.031496}, gaf0015{rf=0.059055}
         gaf002{rf=0.07874}, gaf003{rf=0.11811}, gaf01{rf=0.3937}, gaf05{rf=1.9685})
      o = pgaf
      plot{diag, indep=rf}
   }

   # open-loop (no control-law): including the nonlinear
   # controls equations but gcnorm=0
   # using rf instead of vtas - that way all solutions are
   # at valid rf with no extrapolation
   flut {
      id = vso0
      mass=mass,
      stif=gstif,
      gaf = pgaf,
      controls=Cont,
      indep=(rf, freq, sigma),
      gcnorm = 0.0,
      rho = 1.6
      start{freq[2.9:3]}
      target{vtas[20:500], sigma=0, sort=vtas}
		vzid = lco
   }

   # closed-loop vso at gcnorm=0, with the "linearize" option
   flut {
      id = vso
      linearize
      mass=mass
      stif=gstif
      gaf = pgaf
      controls=Cont
      indep=(rf, freq, sigma)
      gcnorm = 0.0
      rho = 1.6
      start{freq[2.9:3]}
      target{vtas[20:500], sigma=0, sort=vtas}
      plot=homotopy
      print=matrices
		vzid = lco
   }
   viz{id=(vso,homotopy.vso,vso0), color, x=sigma, y=freq}
   viz{id=(vso,vso0), color, x=vtas, y=sigma}

   # vsoe at a few vtas searching for sigma=0 (LCO)
   # allowing gain and phase to vary
   flut {
      id=vsoe, source=vso0
      indep=(vtas, sigma[-100:0.1], freq, gcnorm[0:1], gain, phase)
      optimize{sigma=1.0}
      start{vtas=(50,100,150,200,250,300)}
      target{sigma=0}
		vzid = lco
   }
   viz{id=vsoe, x=sigma, y=gcnorm}

   # voe: lco starting from vsoe
   flut{
      id=voe, source=vsoe,
      indep=(vtas[20:500],freq,gcnorm),
      sigma=0,
      start{vtas[20:500]},
   }
   viz{id=voe, x=vtas, y=gcnorm}

   # phase variation
   flut {
      id=phase, source=vso,
      linearize
      indep=(vtas[0:500], freq, phase)
      sigma=0
   }
   viz {id=phase, x=phase, y=vtas}

   # gain variation
   flut {
      id=gain, source=vso,
      linearize
      indep=(vtas[0:500], freq, gain),
      sigma=0,
   }
   viz{id=gain, x=gain, y=vtas}

   # gain-phase optimization for contour start points
   flut {
      id=opt, source=vso,
      linearize
      indep=(gain, phase, freq, vtas),
      optimize = vtas
      sigma=0
   }
   viz {id=opt, x=gain, y=vtas}
   viz {id=opt, x=phase, y=vtas}

   # gain-phase contours at specific velocities
   flut {
      id=contour, source=opt,
      linearize
      indep=(gain, phase, freq),
      vtas=(360,370,380,390,400,420,440)
   }
   viz {id=(contour,opt), x=gain, y=phase}
end

controls.cpp {{

complex<Ad> cgain(pset& plt, vector<complex<Ad>>& sensor);
Ad gfactor (Ad& absfb);
Ad pfactor (Ad& absfb);

inline complex<Ad>
numer3(double a, double b, double c, complex<Ad> s) { return a + s*(b + s*c);}
inline complex<Ad>
denom3(double a, double b, double c, complex<Ad> s) { return a + s*(b + s*c);}
inline complex<Ad>
denom2(double a, double b, complex<Ad> s) { return a + s*b; }
inline complex<Ad>
denom4(double a, double b, double c, double d, complex<Ad> s) {
   return (a + s*(b + s*(c + s*d))); }
vector<complex<Ad>> getsensor(pset& plt, const string& name);

int
controls (pset& plt, int nr, int nc, vector<complex<Ad>>& ca) {
   Trace trc(1,"controls");
   int fbl1b{44};  // 1b left aileron: feedback dof
   int fbl0b{43};  // 1b left aileron: feedback dof
   int fbr1b{45};  // 1b right aileron: feedback dof
   int fbr0b{44};  // 1b right aileron: feedback dof
   static int visit{0};

   if ((int)ca.size() != nr*nc) {
      throw runtime_error(vastr("controls: (",nr,',',nc,") got lenc ",ca.size()));
   }

   int nce = 2;  // number of extra equations
   int n = nr - nce; // number of normal modes
   trc.dprint("n = ", n, ", nr = ", nr);

   // get some parameter values
   Ad sigma = plt.parval("sigma");
   Ad freq = plt.parval("freq");
   complex<Ad> s(sigma, freq);
   Ad clfac = plt.parval("clfac");
   Ad gcnorm = plt.parval("gcnorm");

   // get the sensors as AD row vectors
   vector<complex<Ad>> lsens = getsensor(plt, "ztipl");
   vector<complex<Ad>> rsens = getsensor(plt, "ztipr");

   // insert the sensors into the last 2 rows; they will be multiplied
   // by the eigenvector in the flutter solution, but we want them
   // multiplied by the gc, so multiply the sensors by gcnorm
   bool linearize{false};
   Settings::defaults.get("linearize", linearize);
	// Specs& sp = flutspecs();
   if (visit == 0 && linearize)
      cerr << "linearizing the controls equations\n";
   for (int j=0; j<(int)lsens.size(); j++) {
      if (linearize) {
         ca[IJ(n,j,nr)] = lsens[j];
         ca[IJ(n+1,j,nr)] = rsens[j];
      } else {
         ca[IJ(n,j,nr)] = gcnorm*lsens[j];
         ca[IJ(n+1,j,nr)] = gcnorm*rsens[j];
      }
   }

   complex<Ad> cgainl = cgain(plt, lsens);
   complex<Ad> cgainr = cgain(plt, rsens);

   plt.setpar(cgainl.real(), "cgainl.real");
   plt.setpar(cgainl.imag(), "cgainl.imag");
   plt.setpar(cgainr.real(), "cgainr.real");
   plt.setpar(cgainr.imag(), "cgainr.imag");

   complex<Ad> claw = numer3(1.0, 1.0, 174.1,s)/denom3(1.0, 4.5, 174.1, s);

   // control law
   // ca(n+1,n+1) = 0.01*cgainl*claw;   // 1b
   // ca(n+2,n+2) = 0.01*cgainr*claw;
   ca[IJ(n,n,nr)] = -claw;   // 1b
   ca[IJ(n+1,n+1,nr)] = -claw;

   // feedback
   complex<Ad> fbl = matvij(plt, "stif", fbl0b, fbl0b);
   ca[IJ(fbl0b,n,nr)] = cgainl*clfac*fbl;
   complex<Ad> fbr = matvij(plt, "stif", fbr0b, fbr0b);
   ca[IJ(fbr0b,n+1,nr)] = cgainl*clfac*fbr;

   // nacelle damping
   // complex<Ad> cj = s*freq*freq*0.1;
   // trc.dprint("set diags 46-51 to ", cj, " = s*freq^2*0.1, s=", s);
   // for (int j=46; j<=51; j++)
     //  ca(j,j) = cj;

#ifdef NEVER
   if (visit == 0) {
      cerr << "T-matrix\n" << ca << endl;
   }
#endif // NEVER
   visit++;
   return 0;
}

vector<complex<Ad>>
getsensor(pset& plt, const string& name) {
// get a sensor matrix, convert it to vector<complex<Ad>> by calling
// it's eval() function; this also allows for parameterized
// sensor matrices.
   Matrix* sens = Matrix::find_mid(name);
   if (sens == nullptr) {
      vector<Matrix*> mp = Matrix::fetch(name);
		if (mp.size() == 1)
			sens = mp[0];
      if (sens == nullptr)
         throw runtime_error(vastr(name," is not available"));
      // put it in the global list of matrices for this analysis
      Matrix::insert(sens);
   }
   vector<complex<Ad>> rval(sens->rsize()*sens->csize());
   sens->eval(plt, rval);
   return rval;
}

complex<Ad>
cgain (pset& plt, vector<complex<Ad>>& sensor) {
// Computes the complex gain (gain*exp(i*phase)) as a
// function of the sensor amplitude
   complex<Ad> rval{0.0};

   // get the eigenvector as an AD vector
   vector<complex<Ad>> ev;
   plt.getadev(ev);
   // eigenvector parameters may not have been defined yet: return 0
   if (ev.empty())
      return complex<Ad>(0.0);

   // sensor should be (1,n), i.e. a row matrix
   complex<Ad> t;
   for (size_t j = 0; j<sensor.size(); j++)
      t += sensor[j]*ev[j];
   
   Ad gcnorm = plt.parval("gcnorm");
   Ad absgc = abs(gcnorm*t);
   plt.setpar(absgc, "ztip");

   // Ad absgc = plt.parval(vastr("absgc",dof1b));
   Ad gfac = gfactor(absgc);
   plt.setpar(gfac, "gfactor");

   Ad gain = plt.parval("gain");
   Ad pfac = pfactor(absgc);
   Ad phase = plt.parval("phase");
   //!! rval = gain*gfac*exp(complex<Ad>(0.0, phase*pfac));
   rval = gain*exp(complex<Ad>(0.0, phase));
   return rval;
}

Ad
gfactor (Ad& absfb) {
//------------------------------------------------------
// compute gain factor as function of abs(fbdof)
//------------------------------------------------------
   Ad rval{1.0};
   static Interp* interp{nullptr};
	double rvmin{0.0};
	double rvmax{2.0};
	if (absfb.value() < rvmin)
		absfb = rvmin;
	if (absfb.value() > rvmax)
		absfb = rvmax;

   if (interp == nullptr) {
      double width=0.01;
      double nom{1.0};   // nominal factor
      vector<double> rv;
      vector<double> rks;

      rv.push_back(rvmin);
      rks.push_back(nom);

      rv.push_back(0.01);
      // rks.push_back(2.0*nom);
      // rks.push_back(0.5*nom);
      rks.push_back(10.0*nom);

      // rv.push_back(0.03);
      // rks.push_back(3.0*nom);
      // rks.push_back(0.25*nom);
      rv.push_back(0.05);
      rks.push_back(20.0*nom);

      // rv.push_back(0.06);
      // rks.push_back(6.0*nom);
      // rks.push_back(0.125*nom);
      // rks.push_back(40.0*nom);
      rv.push_back(0.07);
      rks.push_back(20.0*nom);

      rv.push_back(rvmax);
      // rks.push_back(6.0);
      // rks.push_back(0.0625);
      // rks.push_back(rks[rks.size()-1]);
      rks.push_back(nom);

      interp = plsmooth(rv,rks,width,"gfactor","absfb","gfactor");
      // plot the interpolant
      interp->plot("gfactor", 100);
   }

   // evaluate the interpolant
   rval = interp->eval(absfb);

   return rval;
}

Ad
pfactor (Ad& absfb) {
//------------------------------------------------------
// compute phase factor as function of abs(fbdof)
//------------------------------------------------------
   Ad rval{1.0};
   static Interp* interp{nullptr};
	double rvmin{0.0};
	double rvmax{2.0};
	if (absfb.value() < rvmin)
		absfb = rvmin;
	if (absfb.value() > rvmax)
		absfb = rvmax;

   if (interp == nullptr) {
      double width=0.005;
      double nom{1.0};

      vector<double> rv;
      vector<double> rks;

      rv.push_back(0.0);
      rks.push_back(nom);

      rv.push_back(0.01);
      rks.push_back(2.0*nom);

      rv.push_back(0.1);
      rks.push_back(5.0*nom);

      rv.push_back(2.0);
      rks.push_back(rks[rks.size()-1]);

      interp = plsmooth(rv,rks,width,"pfactor","absfb","pfactor");
      // plot the interpolant
      interp->plot("pfactor", 100);
   }

   // evaluate the interpolant
   rval = interp->eval(absfb);

   return rval;
}

}}

extract.m {{
% gct is dimensioned (3*nnodes,nmodes) so to get the
% dz displacement for nodes 133 and 33:
l = gct__lco(3*133,:);
r = gct__lco(3*33,:);
lsign = sign(l);
rsign = sign(r);
mid = (abs(l)+abs(r))/2.0;
ztipl = lsign.*mid;
ztipr = rsign.*mid;
ldiff = ztipl - l
rdiff = ztipr - r
}}

gaffcn.cpp {{
int
gaffcn (pset& plt, int nr, int nc, vector<complex<Ad>>& ca) {
// bilinear factor applied to the ailerons
   Trace trc(2,"gaffcn");

   if (ca.size() != nr*nc)
      throw runtime_error(vastr("gaffcn: (",nr,',',nc,") got ca.size() ",ca.size()));

//  44     1b left aileron rotation
//  45     1b right aileron rotation

   Ad lfac = df::bilinear(plt, 0.0001, 0.2, 44);
   for (int i=0; i<nr; i++)
      ca[IJ(i,43,nr)] *= lfac;
   Ad absgcl = plt.absgcv(44);
   plt.setpar(lfac, "ailfacl");

   Ad rfac = df::bilinear(plt, 0.0001, 0.2, 45);
   for (int i=1; i<=nr; i++)
      ca[IJ(i,44,nr)] *= rfac;
   Ad absgcr = plt.absgcv(45);
   plt.setpar(rfac, "ailfacr");

   return 0;
}
}}
