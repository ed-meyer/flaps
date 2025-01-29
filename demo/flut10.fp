#-------------------------------------------------------------   
# basic pk flutter analyses using NASTRAN Aeroelastic Handbook
# example ha145b.dat found in the MSC/Nastran Test Problem Library
# page 350-?
#
# Description:
# Checks:
#   import of a NASTRAN OUTPUT4 file
#
# The flutter results should compare with the NASTRAN
# flutter results with one exception: mode 1 diverges
# at about 850 ktas in NASTRAN but diverges at about
# 1300 ktas in Apex flut. The reason for this is the
# difference in formulation of the flutter equation
# between Apex and NASTRAN (see Ref 1 pg 58).
#
# NASTRAN splits the aerodynamic matrix into real and imaginary
# parts, then modifies the imaginary part to appear
# as a viscous damping term in the equations, resulting
# in a real, quadratic eigenvalue problem. This is equivalent
# to adding sigma/omega times the imaginary part to the
# real part of the aero matrix, where the characteristic
# exponent s = sigma + i*omega. This does not change
# the solution much except when sigma is large or omega
# is small, as it is when the wing diverges.
# According to BAH this wing should diverge between 890
# and 1100 ktas, depending on the aero method used.
# Neither method actually predicts divergence however: the frequency
# goes to zero but sigma does not. Divergence must be computed
# using a different method.
#
# However, note what happens when we simulate a v-g type
# solution in flut (making sdamp indep instead of sigma):
# the diverging modes behave just as they should, with mode
# 1 diverging at 976 knots.
#
# The NASTRAN manual gives 1651 ft/s = 978 knots as the
# divergence speed.
#-------------------------------------------------------------

# [1] Read the ASCII file from NASTRAN
#     or run nastran using data in block "nastran"
#     and create an OUTPUT4 file

   import{${FROOT}/demo/flut10.op4}

# [1] Cubic-spline interpolate the aero wrt k-value
#     Note: Nastran takes ref chord as input (131.232)
#     but uses semi-chord in the definition of k-value,
#     so I set rf = k/131.232/2 where k is the Nastran
#     k-value

   param {
      out=GAF,
      i=QHHL1{rf=0.00000001524},
      i=QHHL2{rf=0.00001524},
      i=QHHL3{rf=0.00076200},
      i=QHHL4{rf=0.001524},
      i=QHHL5{rf=0.003048},
      i=QHHL6{rf=0.00762},
      i=QHHL7{rf=0.01524},
      plot = diag
   }

# [2] Neutral-stability CMCD (pk) run

   flut {
      id=pk,
      startregion{freq[0:20]},
      indep=(vtas[0:2500],freq,sigma),
      g = sigma,
      alt=0.0,
      mass=MHH,stif=KHH, gaf=GAF,
      title="EXAMPLE HA145B: BAH JET TRANSPORT WING FLUTTER",
      print=matrices,
      target{vtas=1000}
   }
   vis {id=pk, x=veas, y=sigma }
   vis {id=pk, x=sigma, y=freq }
end
# [2a] repeat using NASTRAN's aero modification

   flut {
      id=pk_nastran
      aeromod=nastran
      startregion{freq[0:20]}
      indep=(vtas[0:2500],freq,growth),
      alt=0.0,
      mass=MHH,stif=KHH, gaf=GAF,
      title="EXAMPLE HA145B: BAH JET TRANSPORT WING FLUTTER"
      target{vtas=1000}
   }

# [2b] ... and again using g-method
#   XXX mode 3 has convergence trouble above 2000 kts - maybe
#       because we are not using the second derivative when computing
#       partials of Q - see lib/Matrix.c:addGmethodDeriv()

   flut {
      id=pk_gmethod
      aeromod=gmethod
      startregion{freq[0:20]}
      indep=(vtas[0:2500],freq,growth),
      alt=0.0,
      mass=MHH,stif=KHH, gaf=GAF,
      title="EXAMPLE HA145B: BAH JET TRANSPORT WING FLUTTER"
      print=matrices
      target{vtas=1000}
   }
   vis {id=(pk, pk_nastran, pk_gmethod), color, x=veas, y=sigma }
   vis {id=(pk, pk_nastran, pk_gmethod), color, x=sigma, y=freq }

# [3] Interpolate aero wrt rf_nastran = rf*65.616
#     see comments above about the Nastran definition of k-value

   param {
      out=gaf,
      i=QHHL1{rf_nastran=0.000001}
      i=QHHL2{rf_nastran=0.001}
      i=QHHL3{rf_nastran=0.05}
      i=QHHL4{rf_nastran=0.1}
      i=QHHL5{rf_nastran=0.2}
      i=QHHL6{rf_nastran=0.5}
      i=QHHL7{rf_nastran=1.0}
   }

   catalog{}

# [4] Neutral-stability CMCD (pk) run using
#     rf_nastran interpolated aero.
#     Use rf as an indep and set vtas limits

   flut {
      id=pkb,
      startregion{freq[0:20]}
      indep=(rf[0.1524e-7:0.01524],freq,growth),
      rf_nastran = "rf*65.616"
      alt=0.0,
      vtas[0:2500]
      mass=MHH,stif=KHH, gaf=gaf,
      title="EXAMPLE HA145B: BAH JET TRANSPORT WING FLUTTER"
      print=matrices
      target{vtas=1000}
   }
   vis {id=(pk,pkb), color, x=veas, y=freq, ymin=0., ymax=16.0 }
   vis {id=(pk,pkb), color, x=veas, y=sigma }

# [5] Simulate a v-g solution: indep are vtas, freq,
#     and structural damping with zero growth rate.
   flut {
      id=vg,
      startregion{freq[0:20]}
      indep=(vtas[0:2500],freq,sdamp)
      growth=0
      g = "freq*sdamp*PI"
      alt=0.0
      mass=MHH,stif=KHH, gaf=GAF
      title="EXAMPLE HA145B: BAH JET TRANSPORT WING FLUTTER"
      print=matrices
      target{vtas=900}
   }

# [6] compare pk with the NASTRAN aero modification, the
#     g-method aero modification, and the v-g solution
#     Note that the pk solution does not diverge until way after the
#     true divergence speeds (976 and 2254)

   vis {id=(pk, pk_nastran, pk_gmethod, vg), color, x=veas, y=freq, ymin=0., ymax=16.0 }
   vis {id=(pk, pk_nastran, pk_gmethod, vg), color, x=veas, y=sdamp, y=growth, ymin=-1, ymax=1 }
   vis {id=(pk, pk_nastran, pk_gmethod, vg), color, x=g, y=freq, xmin=-5, xmax=5 }

# Check targets

   bool ok = true;
   vector<Real> values;
   values.push_back(1.76559);
   values.push_back(2.72031);
   values.push_back(7.0825);
   values.push_back(11.7408);
   values.push_back(11.1434);
   checkTarget("pk", "vtas", "freq", values, 3) && ok;

   values.clear();
   values.push_back(2.74505);
   values.push_back(7.03514);
   values.push_back(11.7411);
   values.push_back(11.0842);
   checkTarget("pk_nastran", "vtas", "freq", values, 3) && ok;

   values.clear();
   values.push_back(2.73565);
   values.push_back(7.03455);
   values.push_back(11.7411);
   values.push_back(11.0938);
   checkTarget("pk_gmethod", "vtas", "freq", values, 3) && ok;

   values.clear();
   values.push_back(1.76557);
   values.push_back(2.72031);
   values.push_back(7.08251);
   values.push_back(11.7408);
   values.push_back(11.1434);
   checkTarget("pkb", "vtas", "freq", values, 2) && ok;

   values.clear();
   values.push_back(1.22176);
   values.push_back(2.81357);
   values.push_back(7.13156);
   values.push_back(11.9078);
   values.push_back(11.6782);
   checkTarget("vg", "vtas", "freq", values, 2) && ok;

   if (!ok)
      exit(1);

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
