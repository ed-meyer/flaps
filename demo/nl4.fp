
# small nonlinear problem based on HA145B
# Demonstrates the use of the "df" option in pz to generate
# describing functions based on arbitrary descriptions of
# nonlinear behavior in the time domain.
#
# This is the fourth of 4 files demonstrating 4 different
# ways of implementing describing functions
# This one creates 3 different bilinear describing functions,
# one at each ratio in nl1.fp, nl2.fp and nl3.fp, using the C++ code in the
# block "bilinear" in this file. The code is compiled automatically before
# any commands are executed and made available to all flaps commands.
# The dfs are created automatically using the 3 functions in blocks "bilinear"
# and "freeplay" and evaluated with Flaps function dfval in equations
# defined in pz
#
# An example of latent LCO is shown in voe2

   import{${FROOT}/demo/ha145b.op4}

	# add a copy of QHHL7 at a high rf to extend the rf range
	octlab{i=(copy.m,QHHL7), o=QHHL8}

   pz { i = KHH, o=Stif,
      [1,1] *= dfval(bilinear, 1, {0.01, 0.5})
      [2,2] *= dfval(bilinear, 2, {0.01, 0.5})
      [3,3] *= dfval(bilinear, 3, {0.01, 2.0})
      [4,4] *= dfval(bilinear, 4, {0.01, 3.5})
   }

   # the mass matrix has freeplay on dof 10
   pz { i=MHH, o=Mass, [10,10] *= dfval(freeplay, 10, {0.001}) }

   pz {
      units=uscs
      out=GAF,
      i=(QHHL1{rf=0.00000001524}, QHHL2{rf=0.00001524}, QHHL3{rf=0.00076200},
         QHHL4{rf=0.001524}, QHHL5{rf=0.003048}, QHHL6{rf=0.00762},
         QHHL7{rf=0.01524}, QHHL8{rf=0.03}),
      beta = (0.0007, 0.007)
   }


# VSO - start points for SOE, VOE analyses
   flut {
      id=vso,
      indep=(vtas[0:500], freq[0.1:100], sigma),
      start{vtas=10, freq[0:7]},
      gcnorm = 0,
      alt=-2000.0,
      mass=Mass,
      stif=Stif,
      gaf=GAF,
      target{vtas[10:500], sigma=0},
   }
   vz {id=vso, x=vtas, y=sigma}

# VOE process starting from eta=0, sigma=0 solutions from VSO
  flut {
      id=voe0, source=vso,
      indep=(vtas[0:1000], freq, gcnorm[0:0.3]),
      sigma=0,
      start{vtas[20:500]},
   }
   vz {id=voe0, x=vtas, y=gcnorm, w=lcostab}

# SOE processes at a few velocities starting from vso
   flut {
      id=soe, source=vso,
      indep=(gcnorm[0:3], freq, sigma),
      vtas = (100,200,400),
      target{vtas[10:500], sigma=0},
   }
   vz {id=soe, x=sigma, y=gcnorm}

# VOE process starting from soe with sigma=0 
  flut {
      id=voe1, source=soe,
      indep=(vtas[0:1000], freq, gcnorm[0:0.3]),
      start{ordinal=1},
      sigma=0,
   }
   vz {id=voe1, x=vtas, y=gcnorm, w=lcostab}

# VSOE process starting from a few velocities
   flut {
      id=vsoe, source=vso,
      indep=(vtas[10:500], sigma[-100:0.05], freq, gcnorm[0:3]),
      constrained,
      start{vtas=(50,100,200)},
      optimize{sigma=1.0, gcnorm=0.5},
      target{vtas[10:500], sigma=0}
   }
   vz {id=vsoe, x=sigma, y=gcnorm}
   vz {id=vsoe, x=vtas, y=sigma}
   vz {id=vsoe, x=vtas, y=gcnorm}

# VOE process starting from vsoe with sigma=0 
# mode 1 shows a latent LCO
  flut {
      id=voe2, source=vsoe,
      indep=(vtas[0:1000], freq, gcnorm[0:3]),
      start{vtas[20:1000], ordinal=1},
      sigma=0,
      checklooping
   }
   vz {id=(vsoe,voe2), color, c=1, x=vtas, y=gcnorm}
end

bilinear {{
vector<double>
df_qs(const vector<double>& params) {
// return a vector of the abs values of the gc where we want
// the time-domain force function to be evaluated
// Input:
//  params   vector of the value of the parameters passed to
//           dfval as the 3rd argument
   Trace trc(1,"bilinear::df_qs");
   vector<double> rval;
   int nstep{201};
   double qmax{3.0};
   double delta = params[0];
   double del{(qmax-delta)/(double)(nstep-1)};
   // from 0 to delta the df is 1 so start at delta
   for (int i=0; i<nstep; i++)
      rval.push_back(delta+i*del);
   trc.dprint("returning ",rval);
   return rval;
}

double
df_fq(double y, const vector<double>& params) {
// returns the factor which multiplies K to get force at displacement 'y'
// with 2 parameters (params): delta and ratio. The params vector
// is passed by the user to dfval().
   Trace trc(1,"df_fq");
   assert(params.size() == 2);
   double delta = params[0];
   double ratio = params[1];
   double rval = y;
   if (abs(y) > delta) {
      if (y > 0.0)
         rval = delta + ratio*(y-delta);
      else
         rval = -delta + ratio*(y+delta);
   }
   trc.dprint("returning f(",y,") = ",rval);
   return rval;
}

Interp*
df_interp(const string& iid, const vector<double>& params,
   const vector<double>& qs, const vector<double>& fs) {
// construct a set of interpolation coefficients:
// 1) piecewise linear between 0 and delta: 1.0
// 2) a smoothing corner between delta*(1-width):(1+width)
// 3) a cubic spline for fs(qs)

   double delta = params[0];

   // fit a spline to qs/fs unsmoothed
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

freeplay {{
vector<double>
df_qs(const vector<double>& params) {
   vector<double> rval;
   int nstep{1001};  // use many steps to take significant time
   assert(params.size() == 1);
   double gap = params[0];
   // put a bunch of points between gap and 100*gap...
   double qmax{100.0*gap};
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
   //!! coef->plot(iid,2001);
   return coef;
}
}}

copy.m {{
QHHL8 = QHHL7;
}}
