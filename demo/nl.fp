
# small nonlinear problem based on HA145B
# Demonstrates the use of the "df" option in pz to generate
# describing functions based on arbitrary descriptions of
# nonlinear behavior in the time domain.
# custom function:
# - stiffcn:  bilinear stiffness (0.05/2.0) and freeplay on
#             first 4 dof (diagonals)
# 2 describing functions:
#   bldf:  bilinear factor with delta, ratio specified in pz
#   fpdf:  freeplay with gap specified in pz
# GAF: rational-function approximation to a set of unsteady
#      aero matrices at several (real) reduced frequencies
#
# An example of latent LCO is shown in voe2

   import{${FROOT}/demo/ha145b.op4}

   pz {
      i = KHH,
      o=Stif,
      code=stiffcn,
      df = fpdf{dfid="gap05", gap=0.05},
      df = bldf{ dfid="bl1", delta = 0.01, ratio = 0.5},
      df = bldf{ dfid="bl2", delta = 0.01, ratio = 2.0},
      df = bldf{ dfid="bl3", delta = 0.01, ratio = 3.5},
   }

   pz {
      out=GAF,
      i=(QHHL1{rf=0.00000001524}, QHHL2{rf=0.00001524}, QHHL3{rf=0.00076200},
         QHHL4{rf=0.001524}, QHHL5{rf=0.003048}, QHHL6{rf=0.00762},
         QHHL7{rf=0.01524}, QHHL7{rf=0.03}),
      beta = (0.0007, 0.007)
#     code = gaffcn,
#     df = bldf{dfid="bl.5", delta=0.01, ratio=0.5},
   }

# create a control-system matrix
#  pz {
#     code = controls{ size = KHH, extra = 2},
#     o = Cont {
#        gain1(gain 1)[-5:5] = 1,
#        gain2(gain 2)[-5:5] = 1,
#        phase1(phase 1)[-180:180] = 0,
#        phase2(phase 2)[-180:180] = 0,
#     }
#  }

# neutral-stability without controls
# output{debugger=valgrind}
#  flut {
#     id=vso0,
#     indep=(vtas[0:1000], freq[0:100], sigma),
#     start{vtas=10, freq[0:40]},
#     gcnorm = 0.0,
#     alt=-10000.0,
#     mass=MHH,
#     stif=Stif,
#     gaf=GAF,
#     plot = stepsize,
#     target{vtas[20:1000], sigma=0},
#  }
#  vz {id=vso0, x=vtas, y=sigma}

# VSO - start points for SOE, VOE analyses
   flut {
      id=vso,
      indep=(vtas[0:1000], freq[0:100], sigma),
      start{vtas=10, freq[0:40]},
#     start{vtas=10, mode_2},
      gcnorm = 0,
      alt=-10000.0,
      mass=MHH,
      stif=Stif,
      gaf=GAF,
#     controls=Cont,
      plot = stepsize,
      target{vtas[10:1000], sigma=0},
   }
   vz {id=vso, x=vtas, y=sigma}
# VOE process starting from eta=0, sigma=0 solutions from VSO
  flut {
      id=voe0, source=vso,
      indep=(vtas[0:3000], freq, gcnorm[0:3]),
      sigma=0,
      start{vtas[20:1000]},
   }
   vz {id=voe0, x=vtas, y=gcnorm, w=lcostability}

# SOE processes at a few velocities starting from vso
   flut {
      id=soe, source=vso,
      indep=(gcnorm[0:3], freq, sigma),
      vtas = (100,200,400),
      target{vtas[10:1000], sigma=0},
   }
   vz {id=soe, x=sigma, y=gcnorm}

# VOE process starting from soe with sigma=0 
  flut {
      id=voe1, source=soe,
      indep=(vtas[0:3000], freq, gcnorm[0:3]),
      start{ordinal=1},
      sigma=0,
   }
   vz {id=voe1, x=vtas, y=gcnorm, w=lcoStability}

# VSOE process starting from a few velocities
   flut {
      id=vsoe, source=vso,
      indep=(vtas, sigma[-100:0.05], freq, gcnorm[0:3]),
      constrained,
      start{vtas=100},
      start{vtas=300},
      start{vtas=500},
      optimize{sigma=1.0, gcnorm=0.5},
      target{vtas[10:1000], sigma=0}
   }
   vz {id=vsoe, x=sigma, y=gcnorm}
   vz {id=vsoe, x=vtas, y=sigma}
   vz {id=vsoe, x=vtas, y=gcnorm}

# VOE process starting from vsoe with sigma=0 
# mode 1 shows a latent LCO
  flut {
      id=voe2, source=vsoe,
      indep=(vtas[0:3000], freq, gcnorm[0:3]),
      start{vtas[20:3000], ordinal=1},
      sigma=0,
      checklooping
   }
   vz {id=(vsoe,voe2), color, c='mode_1\..*', x=vtas, y=gcnorm}
end



nastran {{
assign output4='flut10.op4',unit=12,status=unknown,form=formatted
$ DEC/CMS REPLACEMENT HISTORY, Element HA145B.DAT
$ *2     6-JUL-1994 14:08:16 A_BOYADJIAN "68 PLUS/G/ CHANGE DBSDIR: TO TPLDIR: FOR INCLUDE CARDS"
$ *1     5-JUL-1994 16:54:46 A_BOYADJIAN "68 PLUS/G/ NEW FOR V68 AERO_SS BOOK"
$ DEC/CMS REPLACEMENT HISTORY, Element HA145B.DAT
ID MSC, HA145B $ E_JOHNSON V68  5-JUL-1994
$ID MSC, HA145B
$$$$$$$$    HANDBOOK FOR AEROELASTIC ANALYSIS EXAMPLE HA145B    $$$$$$$$
$                                                                      $
$       MODEL DESCRIPTION       BAH JET TRANSPORT WING EXAMPLE         $
$                               CANTILEVERED WING WITH TEN BEAM        $
$                               ELEMENTS AND DUMBBELL MASSES           $
$                                                                      $
$       SOLUTION                PK FLUTTER ANALYSIS METHOD USING       $
$                               DOUBLET-LATTICE AERODYNAMICS AT        $
$                               MACH NO. 0.0                           $
$                                                                      $
$                                                                      $
$       OUTPUT                  TABULATED MODAL DEFLECTIONS PLUS       $
$                               X-Y PLOTS OF V-G FLUTTER DATA          $
$                                                                      $
$$$$$$$$                                                        $$$$$$$$
TIME 10 $ TIME IN CPU MINUTES
SOL 145  $ FLUTTER ANALYSIS
$diag=4,14
compile SUBDMAP=FLUTTER list
alter 'DO WHILE ( FLOOP\>=0 )' $
output4 KHH,MHH,QHHL,,//-3/12/0 $
CEND
TITLE = EXAMPLE HA145B: BAH JET TRANSPORT WING FLUTTER ANALYSI HA145B
LABEL = PK FLUTTER METHOD
  ECHO    = BOTH
  SPC     = 1    $ FUSELAGE CONSTRAINT
  SDAMP   = 2000 $ STRUCTURAL DAMPING
  METHOD  = 10   $ MODIFIED GIVENS FOR VIBRATION ANALYSIS
  SVEC    = ALL  $ PRINT VIBRATION MODES
SUBTI = CANTILEVERED, DOUBLET-LATTICE AERODYNAMICS AT MACH NO. 0.0
  FMETHOD = 40   $ PK-FLUTTER METHOD
  DISP    = ALL  $ PRINT FLUTTER MODES
OUTPUT(XYOUT)
  CSCALE 2.0
  PLOTTER NASTRAN
  CURVELINESYMBOL = -6
  YTTITLE = DAMPING  G
  YBTITLE =FREQUENCY  F Hz
  XTITLE  = VELOCITY  V (KEAS)
  XMIN  =    0.
  XMAX  = 2500.
  YTMIN =  -.6
  YTMAX =  +.1
  YBMIN =   0.
  YBMAX =  15.
  XTGRID LINES = YES
  XBGRID LINES = YES
  YTGRID LINES = YES
  YBGRID LINES = YES
  UPPER TICS  = -1
  TRIGHT TICS = -1
  BRIGHT TICS = -1
  XYPLOT VG / 1(G,F) 2(G,F) 3(G,F) 4(G,F) 5(G,F)
  CURVELINESYMBOL = -6
  YTTITLE = DAMPING  G
  YBTITLE =FREQUENCY  F Hz
  XTITLE  = VELOCITY  V (ft/sec)
  XMIN  =    0.
  XMAX  = 1600.
  YTMIN =  -1.7
  YTMAX =   +.1
  YBMIN =   0.
  YBMAX =   2.5
  XTGRID LINES = YES
  XBGRID LINES = YES
  YTGRID LINES = YES
  YBGRID LINES = YES
  UPPER TICS   = -1
  TRIGHT TICS  = -1
  BRIGHT TICS  = -1
  XYPLOT VG / 1(G,F)
BEGIN BULK
$*** * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * ***$
$                                                                       $
$        THE ANNOTATIONS IN THIS INPUT DECK ARE INTENDED TO             $
$        EXPLAIN THE DATA ON THE CARD IMAGES FOR THIS SPECIFIC          $
$        EXAMPLE WITHOUT REFERENCE TO THE VARIOUS MANUALS WHERE         $
$        MORE GENERAL DESCRIPTIONS WILL BE FOUND.                       $
$                                                                       $
$*** * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * ***$
INCLUDE 'TPLDIR:bahstru.dat'
$
INCLUDE 'TPLDIR:bahmass.dat'
$                                                                       $
$                  * * STRUCTURAL CONSTRAINTS * *                       $
$                                                                       $
$        THE SPC1 ENTRY CONSTRAINS DOFS OF THE LISTED GRID POINTS.      $
$                                                                       $
$       SID     C       G1      G2      ETC.
SPC1    1       345     11
$                                                                       $
$                   * * STRUCTURAL DAMPING * *                          $
$                                                                       $
$        THE PARAMETER KDAMP DETERMINES THE MANNER OF INCLUSION         $
$        OF STRUCTURAL DAMPING IN EQUATIONS OF MOTION (SEE HANDBOOK     $
$        FOR DYNAMIC ANALYSIS, SECT. 3.2.2).  IF SET TO -1, MODAL       $
$        DAMPING IS PUT INTO COMPLEX STIFFNESS MATRIX AS STRUCTURAL     $
$        DAMPING.                                                       $
$                                                                       $
$       N       V1      V2
PARAM   KDAMP   +1
$                                                                       $
$        THE TABDMP1 ENTRY DEFINES MODAL DAMPING AS A TABULAR           $
$        FUNCTION OF FREQUENCY.  THE DAMPING LEVELS ARE LINEAR          $
$        BETWEEN THE FREQUENCY AND DAMPING PAIRS AND ARE EXTRAP-        $
$        OLATED OUTSIDE THE TABULATED FREQUENCY RANGE.                  $
$                                                                       $
$       ID                                                              +TDP
TABDMP1 2000                                                            +T2000
$       F1      G1      F2      G2      ETC             ENDT
+T2000  0.0     0.0     10.0    0.0     ENDT
$                                                                       $
$                                                                       $
$                      * * * AERODYNAMIC DATA * * *                     $
$                                                                       $
$                          (SNAIL-IN-SEC SYSTEM)                        $
$                                                                       $
$                        * * ELEMENT GEOMETRY * *                       $
$                                                                       $
$        THE AERO ENTRY SPECIFIES THE AERO COORDINATE SYSTEM, THE       $
$        VELOCITY (USED FOR DATA RECOVERY), THE REFERENCE CHORD         $
$        AND FLUID DENSITY, PLUS SYMMETRY KEYS.  SYMXZ=1 INDICATES      $
$        THAT THE MODEL IS MOUNTED WITH A ROOT REFLECTION PLANE;        $
$        SYMXY = 0 INDICATES THAT THE MODEL IS MOUNTED FAR ENOUGH       $
$        FROM THE FLOOR SO THAT REFLECTION EFFECTS ARE NEGLIGIBLE.      $
$                                                                       $
$       ACSID   VELOCITY REFC   RHOREF  SYMXZ   SYMXY
AERO    1               131.232 1.1468-7  1
$                                                                       $
INCLUDE 'TPLDIR:aero5.dat'
$                                                                       $
$                  * * * SOLUTION SPECIFICATIONS * * *                  $
$                                                                       $
$                   * VIBRATION SOLUTION PARAMETERS *                   $
$                                                                       $
$        THE EIGR ENTRY SPECIFIES THE METHOD OF EXTRACTING THE EIGEN-   $
$        SOLUTIONS OF THE STRUCTURE IN A VACUUM, IN THIS CASE THE       $
$        MODIFIED GIVENS METHOD.  TEN MODES ARE DESIRED, NORMALIZED     $
$        ON THE MAXIMUM DISPLACEMENTS.                                  $
$                                                                       $
$       SID     METHOD  F1      F2              ND                      $
EIGR    10      MGIV                            10                      +EIGR
$       NORM    G       C                                               $
+EIGR   MAX
$                                                                       $
$                      * AERODYNAMIC CONDITIONS *                       $
$                                                                       $
$        ALL COMBINATIONS OF MACH NUMBER AND REDUCED FREQUENCY LISTED   $
$        ON THE MKAERO1 ENTRY AND ITS CONTINUATION CARD WILL BE USED    $
$        TO GENERATE GENERALIZED AERO FORCE MATRICES.  IF MORE THAN     $
$        EIGHT MACH NO.S OR REDUCED FREQUENCIES ARE REQUIRED A SECOND   $
$        MKAERO1 ENTRY IS NECESSARY.                                    $
$                                                                       $
$        M1     M2      M3      ETC
MKAERO1 0.                                                              +MK
$        K1     K2      K3      K4      K5      ETC
+MK     0.0000010.001   0.05    0.10    0.20    0.50    1.0
$                                                                       $
$                     *FLUTTER SOLUTION PARAMETERS *                    $
$                                                                       $
$        THE FLUTTER ENTRY DEFINES THE METHOD OF SOLUTION, IDENTIFIES   $
$        THE FLFACT ENTRIES THAT FOLL0W, SPECIFIES THE INTERPOLATION    $
$        METHOD, THE NUMBER OF ROOTS DESIRED IN THE OUTPUT AND THE      $
$        CRITERION FOR CONVERGENCE (DEFAULT IS 10-3).                   $
$                                                                       $
$       SID     METHOD  DENS    MACH    VEL     IMETH   NVALUE  EPS     $
FLUTTER 40      PK      1       2       4       L       10
$                                                                       $
$        FLFACT ENTRIES ARE USED TO SPECIFY DENSITY RATIOS, MACH NO.S   $
$        AND REDUCED FREQUENCIES/VELOCITIES FOR FLUTTER ANALYSES.       $
$        NEGATIVE VELOCITIES ARE SIGNALS TO COMPUTE AND PRINT EIGEN-    $
$        VECTORS.                                                       $
$                                                                       $
$       SID     F1      F2      F3      F4      F5      F6      F7      $
FLFACT  1       1.                                                      DENSITY
FLFACT  2       .0                                                      MACH NO
FLFACT  4       4800.   6000.   7200.   8400.   9600.   10800.  -12000. +FLF4
+FLF4   -13200. 14400.  15600.  16800.  16920.  17040.  17100.  17112.  +FLF4A
+FLF4A  17124.  17136.  17148.  17160.  18000.  -19200. -20400. 21600.  +FLF4B
+FLF4B  22800.  24000.  25200.                                          VELOCITY
$                                                                       $
$        THE PARAM,LMODES,N ENTRY SPECIFIES THAT N MODES ARE TO BE      $
$        USED IN THE FLUTTER ANALYSIS.                                  $
$                                                                       $
PARAM   LMODES  10
$                                                                       $
$        THE PARAM,VREF,C ENTRY SPECIFIES A CONVERSION FACTOR TO BE     $
$        USED TO CONVERT THE DIMENSIONS OF THE OUTPUT VELOCITIES BY     $
$        DIVIDING BY C, IN THIS CASE BY 12.0 IN/FT TO PRINT VEL-        $
$        OCITIES IN FT/SEC RATHER THAN IN/SEC.                          $
$                                                                       $
PARAM   VREF    20.253
$PARAM   VREF    12.0
$
ENDDATA
}}

stiffcn {{

int
stiffcn (pset& plt, int nr, int nc, CMatrix& ca) {
   Trace trc(1,"stiffcn");
   Timer tim("stiffcn");

#ifdef NEVER
   vector<int> dof;
   int ndof = nr;

   for (int i=0; i<ndof; i++)
      dof.push_back(i+1);

   for (size_t j=0; j<dof.size(); j++) {
      int k1b = dof[j];
      int k0b = k1b - 1;
      Real fac = dfval(plt, "bl2", k1b);
      plt.setpar(fac, "bl2");
      ca(k0b,k0b) *= fac;
   }
#else
   // dof 1 & 2 get bl1
   int k1b = 1;
   int k0b = k1b - 1;
   Real fac = dfval(plt, "bl1", k1b);
   ca(k1b,k1b) *= fac;
   k1b = 2;
   k0b = k1b - 1;
   fac = dfval(plt, "bl1", k1b);
   ca(k1b,k1b) *= fac;
   // dof 3 gets bl2
   k1b = 3;
   k0b = k1b - 1;
   fac = dfval(plt, "bl2", k1b);
   ca(k1b,k1b) *= fac;
   // dof 4 gets bl3
   k1b = 4;
   k0b = k1b - 1;
   fac = dfval(plt, "bl3", k1b);
   ca(k1b,k1b) *= fac;
#endif
   return 0;
}
}}

controls.cpp {{

int controls (pset& plt, int nr, int nc, CMatrix& ca) {

   double degprad(45.0/atan(1.0));
   vector<Real> gain;
   vector<Real> phase;
   gain.push_back(plt.parval("gain1"));
   gain.push_back(plt.parval("gain2"));
   phase.push_back(plt.parval("phase1")/degprad);
   phase.push_back(plt.parval("phase2")/degprad);
   size_t m = gain.size();
   Real one(1.0);
   Real k2k1(2.0);
   Real x0(0.01);
   Builtin bi;
   for (size_t j=1; j<=m; j++) {  // 1b
      // size_t k0b = nr - m + j;
      size_t k1b = nr - m + j;  // last m row/col 0b
      ca(k1b,j) = Complex(one);
      ca(k1b,k1b) = Complex(-one);
      Complex stifjj(matvij(plt, "stif", j, j)); // 1b
      int j1b = j+1; // 1b dof
      Real fac = bi.bilinear(plt, x0, k2k1, j);
      Complex ct(fac*gain[j-1]*cos(phase[j-1]), gain[j-1]*sin(phase[j-1]));
      // feedback to dof 2,3 0b:
      ca(j+m,k1b) = ct*stifjj;
   }
   return 0;
}
}}

gaffcn.cpp {{

int
gaffcn (pset& plt, int nr, int nc, CMatrix& ca) {
   Trace trc(1,"gaffcn");
   Timer tim("gaffcn");
   ostringstream os;

   Real dpress = plt.parval("dpress");
   vector<int> dof;  // 1b columns
   int ndof = nr;
   for (int i=0; i<ndof; i++)
      dof.push_back(i+1);

   for (size_t j=0; j<dof.size(); j++) {
      int k1b = dof[j];
      int k0b = k1b - 1;
      Real fac = dfval(plt, "bl.5", k1b);

      for (size_t i=1; i<=nr; i++) {
         ca(i,k1b) *= fac*dpress;
      }
   }
   return 0;
}
}}

fpdf {{

vector<double>
df_qs() {
   vector<double> rval;
   int nstep{101};
   double qmax{3.0};
   // gap is made available automatically from the df
   // option in pz
   double del{(qmax-gap)/(double)(nstep-1)};
   for (int i=0; i<nstep; i++)
      rval.push_back(gap+i*del);
   return rval;
}

double
df_fq(double wt) {
// returns the factor which multiplies K to get force
// at dimensionless time "wt"
   Trace trc(1,"df_fs");
   double sinwt = sin(wt);
   double y = q*sinwt;
   double rval{0.0};
   if (abs(y) > gap) {
      if (y > 0.0)
         rval = y-gap;
      else
         rval = y+gap;
   }
   // if (abs(y) > gap) return (y-gap)*sinwt;
   // return 0.0;
   rval *= sinwt;
   trc.dprint("returning ",rval);
   return rval;
}

void
df_interp(const vector<double> qs, const vector<double> fs) {
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
   coef.add(cp);
   coef.add(pl);
   coef.add(sp);
   coef.plot(dfid,"Gap DF",2001,false);
}
}}

bldf {{

vector<double>
df_qs() {
   vector<double> rval;
   int nstep{101};
   double qmax{3.0};
   double del{(qmax-delta)/(double)(nstep-1)};
   for (int i=0; i<nstep; i++)
      rval.push_back(delta+i*del);
   return rval;
}

double
df_fq(double wt) {
// returns the factor which multiplies K to get force
// at dimensionless time "wt"
   Trace trc(1,"df_fq");
   double sinwt = sin(wt);
   double y = q*sinwt;
   double rval = y;
   if (abs(y) > delta) {
      if (y > 0.0)
         rval = delta + ratio*(y-delta);
      else
         rval = -delta + ratio*(y+delta);
   }
   rval *= sinwt;
   trc.dprint("returning ",rval);
   return rval;
}

void
df_interp(const vector<double> qs, const vector<double> fs) {
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
   coef.add(cp);
   coef.add(pl);
   coef.add(sp);
   coef.plot(dfid,"Bilinear DF",2001,false);
}
}}
