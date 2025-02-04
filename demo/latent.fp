
#------------------------------------------------------------------
# generic 2-engine airplane
# freeplay on elevators, rudder, ailerons
# bilinear stiffness on nacelles
# bilinear factor on gaf
# Shows a latent lco
#------------------------------------------------------------------

# nonlinear dofs implemented in stiffcn and gaffcn:
#  41     rudder rotation
#  42     left elevator rotation
#  43     right elevator rotation
#  44     left aileron rotation
#  45     right aileron rotation
#  46-51  nacelles 3 dof/side

   import { ${FROOT}/demo/lco.op4}
   import { ${FROOT}/demo/lco.uf}
   catalog{}

   # give the stiffness matrix an evaluation function with freeplay and bilinear
   pz {
      i=gstif
      code=stiffcn
      o=stif
   }

   # gaf matrices are interpolated wrt rf, then passed
   # to gaffcn for nonlinear mods
   pz {
      i=(gaf0{rf=0}, gaf0008{rf=0.031496}, gaf0015{rf=0.059055}
         gaf002{rf=0.07874}, gaf003{rf=0.11811}, gaf01{rf=0.3937}, gaf05{rf=1.9685})
      o = pgaf
   }

   # VSO at gcnorm=0
   flut {
      id = vso
      mass=mass
      stif=stif
      gaf = pgaf
      indep=(vtas[5:70], freq[0:50], sigma)
      rho = 1.6
      gcnorm = 0.0
      start{vtas=5, modes=14}
      target{sigma=0,vtas[10:70]}
      freevib
      vzid = lco
   }
   viz {id=vso, x[-1:1]=sigma, y=freq}
   viz {id=vso, x=vtas, y=sigma}

   # vsoe starting at 50 m/s
   flut {
      id=vsoe, source=vso
      indep=(vtas[0:100], sigma[-100:0.001], freq, gcnorm[0:1])
      start{vtas=50}
      optimize{sigma=0.99, gcnorm=0.14106736}
      target{sigma=0}
   }
   viz{id=vsoe, x=sigma, y=gcnorm}

   # voe starting from vsoe
   flut {
      id=voe2, source=vsoe
      indep=(vtas[10:500], freq, gcnorm)
      sigma=0
      checklooping
   }
   viz {id=voe2, x=vtas, y=gcnorm, w=lcostab}
   viz {id=(voe2,vsoe), x=vtas, y=gcnorm}
   # these are the only gc that go beyond delta
   viz {id=voe2, x=vtas, y='fac4@([1-5])'}
   # check the interpolations of the 3 df's
   viz {x=q, y=f, bilinear@0.001@0.5.apf}
   viz {x=q, y=f, bilinear@0.01@1.2.apf}
   viz {x=q, y=f, gap@0.001.apf}
end

stiffcn.cpp {{

int
stiffcn (pset& plt, int nr, int nc, vector<complex<Ad>>& ca) {
   Trace trc(2,"stiffcn");

   if (ca.size() != nr*nc)
      throw runtime_error(vastr("stiffcn: (",nr,',',nc,") got ca.size() ",ca.size()));

// "DOF=(7 to 46,112,118,124,129,134,135 to 140)"
//  41     rudder rotation
//  42     left elevator rotation
//  43     right elevator rotation
//  44     left aileron rotation
//  45     right aileron rotation
//  46-51  nacelles 3 dof/side

   // Ad fpx0 = plt.parval("fpx0");  // rudder tab
   // Builtin bi;
   vector<int> fpdof;
   fpdof.push_back(41);
   fpdof.push_back(42);
   fpdof.push_back(43);
   fpdof.push_back(44);
   fpdof.push_back(45);

   size_t j;
   for (j=0; j<fpdof.size(); j++) {
      int k1b = fpdof[j];
      int k0b = k1b - 1;
      // Ad fac = bi.freeplay(plt, fpx0, k1b);
//     df = gap{dfid="gap", gap=0.001}
      Ad fac = dfval(plt, "gap", k1b, {0.001});
      plt.setpar(fac, vastr("fac",k1b));
      complex<Ad> t = ca[IJ(k0b,k0b,nr)];
      ca[IJ(k0b,k0b,nr)] *= fac;
      trc.dprint("stiffcn: scaled diag ",k1b," = ",ca[IJ(k0b,k0b,nr)],
         " = ",t,'*',fac);
   }

// bilinear stiffness dof
//  nacelle 3 dof/side 46-51

   vector<int> bldof;
   // Ad blx0 = plt.parval("blx0");
   // Ad k2k1 = plt.parval("k2k1");

   for (int i=0; i<6; i++)
      bldof.push_back(46+i);

   for (j=0; j<bldof.size(); j++) {
      int k1b = bldof[j];
      int k0b = k1b - 1;
      // Ad fac = bi.bilinear(plt, blx0, k2k1, k1b);
//     df = bilinear{ dfid="bl12", delta = 0.01, ratio = 1.2}
      Ad fac = dfval(plt, "bilinear", k1b, {0.01, 1.2});
      plt.setpar(fac, vastr("fac",k1b));
      complex<Ad> t = ca[IJ(k0b,k0b,nr)];  // 1b
      ca[IJ(k0b,k0b,nr)] *= fac;
      trc.dprint("stiffcn: scaled diag ",k1b," = ",ca[IJ(k0b,k0b,nr)],
         " = ",t,'*',fac);
   }

   return 0;
}
}}

gaffcn.cpp {{
int
gaffcn (pset& plt, int nr, int nc, vector<complex<Ad>>& ca) {
// bilinear factor applied to some of the low-frequency modes
   Trace trc(2,"gaffcn");

   if (ca.size() != nr*nc)
      throw runtime_error(vastr("gaffcn: (",nr,',',nc,") got ca.size() ",ca.size()));

//  41     rudder rotation
//  42     left elevator rotation
//  43     right elevator rotation
//  44     left aileron rotation
//  45     right aileron rotation
//  46-51  nacelles 3 dof/side

   vector<int> dof;
   // first few low-freq modes
   for (int k=1; k<5; k++)
      dof.push_back(k);
   //  nacelle 3 dof/side 46-51
   // for (int k=46; k<52; k++)
     //  dof.push_back(k);
   // control-surfaces
   // dof.push_back(41);
   // dof.push_back(42);
   // dof.push_back(43);
   // dof.push_back(44);
   // dof.push_back(45);
   // Ad gafx0 = plt.parval("gafx0");
   // Ad gafratio = plt.parval("gafratio");
   //!! complex<Ad> vfac = plt.parval("dpress");
   complex<Ad> vfac{1.0};

   for (size_t j=0; j<dof.size(); j++) {
      int k1b = dof[j];
      int k0b = k1b - 1;
      // Ad fac = bi.bilinear(plt, gafx0, gafratio, k1b);
//  df = bilinear{ dfid="blgaf", delta = 0.001, ratio = 0.5}
      Ad fac = dfval(plt, "bilinear", k1b, {0.001, 0.5});
      plt.setpar(fac, vastr("fac",k1b));
      // scale the entire column
      for (size_t i=1; i<=nr; i++) {
         ca[IJ(i,k0b,nr)] *= fac*vfac;
      }
   }
   return 0;
}
}}

gap {{
// functions for describing a freeplay describing function
// the pz option must have a value for the gap, e.g.
//   pz { df=gap{dfid="gap", gap=0.001}, ...}

vector<double>
df_qs(const vector<double>& params) {
   vector<double> rval;
   int nstep{101};
   assert(params.size() == 1);
   double gap = params[0];
   // put a bunch of points between gap and 1...
   //!! double qmax{1.0};
   double qmax{0.1};
   double del{(qmax-gap)/(double)(nstep-1)};
   for (int i=0; i<nstep; i++)
      rval.push_back(gap+i*del);
   // ... and a single point at 3
   // rval.push_back(3.0);
   return rval;
}

double
df_fq(double y, const vector<double>& params) {
// returns the factor which multiplies K to get force at displacement 'y'
   Trace trc(1,"df_fs");
   double gap = params[0];
   double rval{0.0};
   if (abs(y) > gap) {
      if (y > 0.0)
         rval = y-gap;
      else
         rval = y+gap;
   }
   trc.dprint("returning ",rval);
   return rval;
}

Interp*
df_interp(const string& iid, const vector<double>& params,
   const vector<double>& qs, const vector<double>& fs) {

   double gap = params[0];
   // fit a spline to qs/fs unsmoothed
   spline* sp = new spline(qs, fs);
   // f(q[0:gap]) = 0: fit a piecewise linear
   vector<double> q0{{0.0}, {gap}};
   vector<double> f0{{0.0}, {0.0}};
   plinear* pl = new plinear(q0,f0);
   // smoothed transition from the piecewise linear to the spline
   double width{0.1};
   double x0{gap*(1.0-width)};
   double y0{0.0};
   double yp0{0.0};
   double x1{gap*(1.0+width)};
   double y1 = sp->eval(x1);
   double yp1 = sp->deriv(x1);
   corner* cp = new corner(x0, y0, yp0, x1, y1, yp1);
   // first interpolant must be the corner
   Interp* coef = new Interp(iid, "q", "f");
   coef->add(cp);
   coef->add(pl);
   coef->add(sp);
   coef->plot(iid,2001);
   return coef;
}
}}

bilinear {{
// functions for describing a bilinear describing function.
// Two parameters (delta,ratio) are passed through dfval either as
// a vector<string> or an initializer-list.

vector<double>
df_qs(const vector<double>& params) {
   vector<double> rval;
   int nstep{101};
   double delta = params[0];
   double ratio = params[1];
   // put a bunch of points between gap and 1...
   double qmax{1.0};
   // double qmax{0.1};
   double del{(qmax-delta)/(double)(nstep-1)};
   for (int i=0; i<nstep; i++)
      rval.push_back(delta+i*del);
   // ... and a single point at 3
   // rval.push_back(3.0);
   return rval;
}

double
df_fq(double y, const vector<double>& params) {
// returns the factor which multiplies K to get force
// at dimensionless time "wt"
   Trace trc(1,"df_fq");
   double rval = y;
   double delta = params[0];
   double ratio = params[1];
   if (abs(y) > delta) {
      if (y > 0.0)
         rval = delta + ratio*(y-delta);
      else
         rval = -delta + ratio*(y+delta);
   }
   trc.dprint("returning ",rval);
   return rval;
}

Interp*
df_interp(const string& iid, const vector<double>& params,
      const vector<double>& qs, const vector<double>& fs) {
   // fit a spline to qs/fs unsmoothed
   double delta = params[0];
   spline* sp = new spline(qs, fs);
   // f(q[0:delta]) = 1: fit a piecewise linear
   vector<double> q0{{0.0}, {delta}};
   vector<double> f0{{1.0}, {1.0}};
   plinear* pl = new plinear(q0,f0);
   // smoothed transition from the piecewise linear to the spline
   double width{0.1};
   double x0{delta*(1.0-width)};
   double y0{1.0};
   double yp0{0.0};
   double x1{delta*(1.0+width)};
   double y1 = sp->eval(x1);
   double yp1 = sp->deriv(x1);
   corner* cp = new corner(x0, y0, yp0, x1, y1, yp1);
   // first interpolant must be the corner
   Interp* coef = new Interp(iid, "q", "f");
   coef->add(cp);
   coef->add(pl);
   coef->add(sp);
   coef->plot(iid,2001);
   return coef;
}
}}
