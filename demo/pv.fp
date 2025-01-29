# for details about this model, see flut10.ax

   import{${FROOT}/demo/ha145b.op4}

# [1] Cubic-spline interpolate the aero wrt k-value
#     Note: Nastran takes ref chord as input (131.232)
#     but uses semi-chord in the definition of k-value,
#     so I set rf = k/131.232/2 where k is the Nastran
#     k-value

# output{d=3}
   pz {
      out=GAF,
      i=QHHL1{rf=0.00000001524},
      i=QHHL2{rf=0.00001524},
      i=QHHL3{rf=0.00076200},
      i=QHHL4{rf=0.001524},
      i=QHHL5{rf=0.003048},
      i=QHHL6{rf=0.00762},
      i=QHHL7{rf=0.01524},
      plot = diag,
#     beta = (0.0005, 0.01)
   }

   print {KHH}
   print {QHHL7}

   new Par{k1fac[0.1:6] = 1.0}
   new Par{k2fac[0.1:4] = 1.0}
   pz {
      i=KHH,
      o=Stif,
      [1,1] = k1fac*KHH[1,1],
      [2,2] = k2fac*KHH[2,2]
   }

# [2] Neutral-stability CMCD (vso) run

# output{d=2}
# output{valgrind}
# output{profile}
# output{wait}
# output{malloc=malloc.out}
   flut {
      id=vso,
   plot=(std,stepsize),
      startregion{freq[3.4:3.6]},  # mode 2
      indep=(vtas[0:2500],freq,growth),
      alt=0.0,
      mass=MHH,
      stif=Stif,
      gaf=GAF,
      title="EXAMPLE HA145B: BAH JET TRANSPORT WING FLUTTER",
      target{sigma=0}
   }
   vz {id=vso, x=veas, y=sigma }

# vso run with multiple values of k1fac & k2fac
   flut {
      id=vso1,
#  plot=(std,stepsize),
      startregion{freq[2.5:5.1]},  # mode 2
      indep=(vtas[0:2500],freq,growth),
      k1fac = (0.8, 0.9, 1.1),
      k2fac = (0.8, 0.9, 1.1),
      alt=0.0,
      mass=MHH,
      stif=Stif,
      gaf=GAF,
      title="HA145B discrete parameter variations",
      target{sigma=0}
   }
   vz {id=vso1, x=sigma, y=freq }

# k1fac variation starting from each sigma=0 point in vso1
   flut {
      id=pv, source = vso1,
      indep=(k1fac, freq, vtas),
      startregion{vtas[500:650]},
      sigma=0,
      curvaturefactor=6,
      plot=(std,stepsize)
   }
   vz {id=pv, x=k1fac, y=vtas}

# compute an optimization curve varying k1fac & k2fac
# output{d=2}
   flut {
      id=opt, source=vso,
      indep=(k1fac, k2fac, freq, vtas[50:1000]),
      optimize=vtas,
      startregion{vtas[500:650]},
      sigma=0
   }
# end
   vz{id=opt, x=vtas, y=k1fac, y=k2fac}

# countours of sigma=0 varying k1fac & k2fac, starting from the
# optimization curve
   flut {
      id=contour, source = opt,
      indep=(k1fac, k2fac, freq),
      vtas=(300, 400, 500, 600, 650, 700, 750, 800),
      curvaturefactor=6,
      plot=(std,stepsize)
   }
   vz {id=(contour,opt), x=k1fac, y=k2fac}
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
