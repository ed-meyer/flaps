# XXX there is something wrong with this model; e.g. the
# free-vibration solution has mode 7 with freq=0.6 but
# normalizing on gc 47, a fore-aft wing tip mode
# Also the nacelles (and others) are in the connectivity
# but there are no nodes for them, so the nacells only show as tiny boxes
#
# 4-engine airplane with closed-loop
# Shows very strange behavior: there is a flutter mode which does
# not extend below about 100 knots eas

import{${FROOT}/demo/flut13.op4}
import{${FROOT}/demo/flut13.uf}

print{stif}
print{mass}
# end
# Create generalized airforce matrix that is a function
# of reduced frequency (K-Value)
param {o=gaf,
	i=gaf12{rf=0.000000},
	i=gaf11{rf=0.000100},
	i=gaf10{rf=0.000200},
	i=gaf9{rf=0.000500},
	i=gaf8{rf=0.001000},
	i=gaf7{rf=0.002000},
	i=gaf6{rf=0.003000},
	i=gaf5{rf=0.005000},
	i=gaf4{rf=0.007000},
	i=gaf3{rf=0.010000},
	i=gaf2{rf=0.020000},
	i=gaf1{rf=0.050000}
}


# Add branch-mode frequency parameters to the stiffness matrix
# param {
#    o=pstif
# 	AEKHH{mass=XXXX,
#       gc=118  {olnsb (OB LHS nacelle side bending      freq)=2.920}
#       gc=119  {olnvb (OB LHS nacelle vertical bending  freq)=5.211}
#       gc=120  {olnts (OB LHS nacelle torsion           freq)=6.955}
#       gc=121  {ilnsb (IB LHS nacelle side bending      freq)=2.562}
#       gc=122  {ilnvb (IB LHS nacelle vertical bending  freq)=4.683}
#       gc=123  {ilnts (IB LHS nacelle torsion           freq)=6.445}
#       gc=124  {irnsb (IB RHS nacelle side bending      freq)=2.562}
#       gc=125  {irnvb (IB RHS nacelle vertical bending  freq)=4.683}
#       gc=126  {irnts (IB RHS nacelle torsion           freq)=6.445}
#       gc=127  {ornsb (OB RHS nacelle side bending      freq)=2.920}
#       gc=128  {ornvb (OB RHS nacelle vertical bending  freq)=5.211}
#       gc=129  {ornts (OB RHS nacelle torsion           freq)=6.955}
#       gc=166  {olail (OB LHS aileron rotation           freq)=20.210}
#       gc=167  {orail (OB RHS aileron rotation           freq)=20.210}
#       gc=168  {ilail (IB LHS aileron rotation           freq)=27.570}
#       gc=169  {irail (IB RHS aileron rotation           freq)=27.570}
#       gc=170  {olelv (OB LHS elevator rotation          freq)=21.440}
#       gc=171  {orelv (OB RHS elevator rotation          freq)=21.440}
#       gc=172  {ilelv (IB LHS elevator rotation          freq)=27.780}
#       gc=173  {irelv (IB RHS elevator rotation          freq)=27.780}
#       gc=174  {uprud (Upper rudder rotation             freq)=18.930}
#       gc=175  {lwtab (tab rotation                      freq)=59.930}
#       gc=176  {lwrud (Lower rudder with tab geared ratio 1.0 freq)=8.160}
#       gc=177  {slab  (lower rudder with tab             freq)=99.0}
#       gc=178  {stabpt (stab pitch                       freq)=34.7}
#       gc=179  {olspl (OB LHS spoiler rotation           freq)=200.0}
#       gc=180  {orspl (OB RHS spoiler rotation           freq)=200.0}
#       gc=181  {mlspl (MB LHS spoiler rotation           freq)=200.0}
#       gc=182  {mrspl (MB RHS spoiler rotation           freq)=200.0}
#       gc=183  {ilspl (IB LHS spoiler rotation           freq)=200.0}
#       gc=184  {irspl (IB RHS spoiler rotation           freq)=200.0}
#    }
# }


#   SIGNALS FOR CONTROL LAWS
# -----------------------------------------------------------------------------------------
#  extract {
#   ADIRUTZ='modes,set=2'{rows=(3)},
#   ADIRURY='modes,set=2'{rows=(5)},
#}


# Create a control-law matrix (Cont) which is evaluated by
# subroutine controlLaw
# param {
#    code=controlLaw
#    size=mass
#    extra=6 
#    o = Cont{
#       gain(Gain)[0:5.0] = 1
#       phase(Phase)[-180:180] = 0
#    }
# }

output{ env = "fuelbif=81.7486" }
output{ env = "fuelhi=84" }
output{ env = "fuello=80" }


# vso soln at gcnorm=0
flut {
   id=vso,
   start{vtas=0, freq[0.1:2.0]},
   gaf=gaf,stif=stif,
   mass=mass,
   target{growth=0},
   indep=(vtas, sigma, freq[0.0:5.0]),
	alt=0,
   title="generic 4-engine"
}
   vis {id=vso, x=vtas, y=sigma}
end

# fuel variation starting with the pk crossing
flut {
   id=pv,source=pk
   sigma = 0
   indep=(alt[0:],freq,fuel[0:100]),
   gain=1.0,
}
end

# work backwards starting with the pv curve at 80
stab {
   id=pk80, source=pv
   indep=(alt,sigma,freq)
   fuel=${fuello}
}
end

controlLaw {{
C     ------------------------------------------------------------------
C     Note: This file is designed to run with the apex program for closed     
C           loop automated control law analyses on linux clusters.                      
C     ------------------------------------------------------------------
C     THE PROGRAM WAS ORIGINALLY WRITTEN BY ARUN NADKARNI FOR THE 777
C     MODELS IT WAS ADAPTED FOR USE WITH APEX/STAB ON LINUX CLUSTERS,              
C     CLEANED UP AND COMMENTS WERE INTRODUCED BY SHIMON LOWY IN OCT. 2005.  
C     THIS PROGRAM IS FOR FLUTTER GROUP USE ONLY.                           
C                                                                            
C     The program was modified in January 2009 For The 747-8 freighter.                              
C     - Scott Gunther 
C     ------------------------------------------------------------------
C     Control Law: Linearized MLA
C     Source: Memo BE326-C08-003 "747-8 Freighter Preliminary Control 
C     Law Release"
C     ------------------------------------------------------------------

      SUBROUTINE controlLaw (NDOF,C) 

C     NAME: Apex Variable Name
C     DTYPE: DataType of the Apex Variable        
C     ------------------------------------------------------------------
      CHARACTER*8 NAME, DTYPE

C     FILTERS: Switch for filters
C     ------------------------------------------------------------------
      LOGICAL FILTERS

C     OPENLOOP: Switch to decouple control law from the structure 
C     ------------------------------------------------------------------
      LOGICAL OPENLOOP

C     N: The number of structural modes
C     ------------------------------------------------------------------
      PARAMETER (N=184) 

C     Control Law Sensor Inputs 
C     ------------------------------------------------------------------
      REAL RY(N)
      REAL TZ(N)

C     C Matrix and characteristic exponent
C     ------------------------------------------------------------------
      COMPLEX C(NDOF,NDOF), P 

C     Transfer Functions
C     ------------------------------------------------------------------
      COMPLEX TF1,TF2,TF3,TF4,TD1,TD2
      COMPLEX NFLT1,NFLT2,NFLT3,NTCHFLTS
      COMPLEX TFPCU1,TFPCU2,TFPCU3
      COMPLEX CM1,CM2,CM3
      REAL G_FLAPS_IB,G_FLAPS_MID,G_MLA2

C     Actuator stiffness
C     ------------------------------------------------------------------
      COMPLEX CRK179,CRK180,CRK181,CRK182,CRK183,CRK184

      complex denom2, denom3, denom4

      REAL VQ(1)
      REAL GAIN,PHASE
      INTEGER FLAPS,SB

C     Gains at connections between CLAW and Struct Model
C     use OPENLOOP to connect or disconnect the control law
C-----------------------------------------------------------------------
      REAL KOPEN

      INTEGER IPMODE
      INTEGER NCE, NDOF
      DATA IPMODE/0/
      SAVE IPMODE
      SAVE RY TZ

C     Check to see if the total number of dofs is correct
C     ------------------------------------------------------------------
      NCE = NDOF - N
      IF (6 .NE. NCE) THEN
        PRINT *,'ERROR: The wrong an incorrect number of control'
        PRINT *,'equations was specified.'
      STOP
      ENDIF

C     Unit Conversions and constants 
C     ------------------------------------------------------------------
      RAD2DEG = 45.0/atan(1.0)
      DEG2RAD = 1.0/RAD2DEG

C     CONTROL LAW INPUT DATA FROM A/C STRUCTURE MEASURED BY ON-BOARD 
C     INSTRUMENTATION
C     ------------------------------------------------------------------
C     Note: fetch reads an Apex array (and name) defined by the Apex
c     "extract" function
C           scopy - copy to a new array with a different name

      IF (IPMODE .EQ. 0) THEN 
C     Airplane motions are defined at the ADIRU (Air Data Inertial 
C     Reference Unit) in terms of vertical translation (DOF 3) or
C     its derivatives
C     ------------------------------------------------------------------
C     --------------------------WARNING---------------------------------
        NAME = 'ADIRUTZ' 
        IPRINT = 0 
        CALL FETCH (NAME,NR,NC,VQ,IPMODE,DTYPE,IPRINT,IRR) 
        CALL SCOPY (N,VQ(IPMODE),NR,TZ, 1) 
        NAME = 'ADIRURY' 
        IPRINT = 0 
        CALL FETCH (NAME,NR,NC,VQ,IPMODE,DTYPE,IPRINT,IRR) 
        CALL SCOPY (N,VQ(IPMODE),NR,RY, 1) 
    
      ENDIF 
   
C     Creating the complex variable "P" and picking current values
C     of parameters
C     Note: Parval in Apex means current value
C     ------------------------------------------------------------------

      p = cmplx(parval('sigma'), parval('freq')) 
C     Ascertain that the run continues even if dividing by a "p" which
C     is zero
      if (cabs(p) .eq. 0.0) p = cmplx(0.001, 0.001) 

C     Get additional parameters from Apex
C     ------------------------------------------------------------------
      GAIN = parval('gain') 
      PHASE = parval('phase') / RAD2DEG 

C     ------------------------------------------------------------------
      FILTERS = .FALSE.

C     Variable used to disconnect the control law from the structure
C-----------------------------------------------------------------------
      OPENLOOP = .FALSE.

      IF (OPENLOOP) THEN
        KOPEN=0.
      ELSE
        KOPEN=1.
      ENDIF

C     MLA Control Law 747-8 Freighter
C     ------------------------------------------------------------------
c     TF1 = 404/(P*P+28.43*P+404)
      TF1 = 404/denom3(1.0, 28.43, 404.0, p)
c     TF2 = (0.3827*P*P+2.7385*P+2265.3)/(P*P+61.008*P+2265.3)
      TF2 = (0.3827*P*P+2.7385*P+2265.3)/denom3(1.0, 61.008, 2265.3, p)
      TF3 = 1/(0.05*P+1)
c     TF3 = 1/(0.05*P+1)
      TF4 = P/(0.05*P+1)
      TD1 = CEXP(-.0517*P)
      TD2 = CEXP(-.0467*P)
      TFPCU1 = 25.0/(P+25.0)
      TFPCU2 = 22.0/(P+22.0)
      TFPCU3 = 25.0/(P+25.0)
      FLAPS = 0
      SB = 1

C     Notch Filters
C     ------------------------------------------------------------------
c     NFLT1= (P*P+2.031*P+208.84)/(P*P+8.086*P+208.84)
      NFLT1= (P*P+2.031*P+208.84)/denom3(1.0, 8.086, 208.84, p)
c     NFLT2= (P*P+3.345*P+799.438)/(P*P+13.32*P+799.438)
      NFLT2= (P*P+3.345*P+799.438)/denom3(1.0, 13.32, 799.438, p)
c     NFLT3= (P*P+4.6536*P+1934.4)/(P*P+18.5263*P+1934.4)
      NFLT3= (P*P+4.6536*P+1934.4)/denom3(1.0, 18.5263, 1934.4, p)

      IF (FILTERS) THEN
        NTCHFLTS = NFLT1*NFLT2*NFLT3
      ELSE
        NTCHFLTS = 1.0
      ENDIF

C     Defining rotional stiffness of the spoiler commands
C     ------------------------------------------------------------------
      call matvij('stif', 179, 179, CRK179) 
      call matvij('stif', 180, 180, CRK180) 
      call matvij('stif', 181, 181, CRK181) 
      call matvij('stif', 182, 182, CRK182) 
      call matvij('stif', 183, 183, CRK183) 
      call matvij('stif', 184, 184, CRK184) 
C     modal DOFs 179,180,181,182,183,184
      C(179,N+1) = (-CRK179/RAD2DEG)*KOPEN
      C(180,N+1) = (-CRK180/RAD2DEG)*KOPEN
      C(181,N+3) = (-CRK181/RAD2DEG)*KOPEN
      C(182,N+3) = (-CRK182/RAD2DEG)*KOPEN
      C(183,N+2) = (-CRK183/RAD2DEG)*KOPEN
      C(184,N+2) = (-CRK184/RAD2DEG)*KOPEN

C     External inputs into the control law - diagonal terms
C     ------------------------------------------------------------------
C     OUTBOARD SPOILERS
      C(N+1,N+1) = -1.0

C     INBOARD SPOILERS
      C(N+2,N+2) = -1.0

C     MID SPOILERS
      C(N+3,N+3) = -1.0

C     ADIRU acceleration in the vertical direction (g)
      C(N+4,N+4) = -386.4

C     ADIRU Pitch Rate
      C(N+5,N+5) = -DEG2RAD

C     Diagonal Terms for control law equations
C     --------------------------------------------------------------

      C(N+6,N+6) = -1.0

C     Control Law Equations
C     ------------------------------------------------------------------
      C(N+6,N+4) = TF1*TD1*TF3
      C(N+6,N+5) = DEG2RAD*2.74*TF1*TD2*TF2*TF4
      CM1 = CEXP(CMPLX(0.,PHASE))
      C(N+1,N+6) =-20*NTCHFLTS*G_MLA2(FLAPS)*TFPCU1*GAIN*CM1
      CM2 = CEXP(CMPLX(0.,PHASE))
      C(N+2,N+6) =20*NTCHFLTS*G_FLAPS_IB(FLAPS,SB)*TFPCU2*GAIN*CM2
      CM3 = CEXP(CMPLX(0.,PHASE))
      C(N+3,N+6) =45*NTCHFLTS*G_FLAPS_MID(FLAPS,SB)*TFPCU3*GAIN*CM3
c
      c64 = cabs(c(n+6,n+4))
      c65 = cabs(c(n+6,n+5))
      c16 = cabs(c(n+1,n+6))
      c26 = cabs(c(n+2,n+6))
      c36 = cabs(c(n+3,n+6))
      call setpar(c64, 'c64')
      call setpar(c65, 'c65')
      call setpar(c16, 'c16')
      call setpar(c26, 'c26')
      call setpar(c36, 'c36')
c
C     Summing aircraft structure input signal contribution from each
C     mode and converting to various dynamic parameters (as needed).
C     ------------------------------------------------------------------
      DO 50 J = 1, N
C     Airplane Pitch Rate at ADIRU (Air Data Inertial Reference Unit) (deg/s)
      C(N+5,J) =  RY(J)*P*KOPEN
C     Airplane Vertical acceleration at ADIRU (Air Data Inertial Reference Unit) (g)
      C(N+4,J) =  TZ(J)*P*P*KOPEN
  50  CONTINUE
 
C     Scale the controls variables 
C     ------------------------------------------------------------------
      scale = 1.0e+6
      DO 500 j = 1, NCE
         DO 400 i = 1, NDOF
            C(i,NDOF-j+1) = scale * C(i,NDOF-j+1)
  400    CONTINUE
  500 CONTINUE
      RETURN
      END

      FUNCTION G_FLAPS_IB(FLAPS,SB)
C-----------------------------------------------------------------------
C G_FLAPS: COMPUTES GAIN VALUES FOR INBOARD AND MID FLAPS
C-----------------------------------------------------------------------
      REAL G_FL_UP_POS_IB,G_FL_UP_NEG,G_FL_DWN
      INTEGER FLAPS,SB

      G_FL_UP_NEG = (0.0-1.0)/(0.1+0.9)
      G_FL_UP_POS_IB = (0.8-0.0)/(1.9-1.2)
      G_FL_DWN = (1.0-0.0)/(1.9-1.3)
      
      IF ( FLAPS .EQ. 0 ) THEN
        G_FLAPS_IB = SB*G_FL_UP_POS_IB
      ELSE
        G_FLAPS_IB = SB*G_FL_DWN
      ENDIF

      RETURN
      END
C-----------------------------------------------------------------------

      FUNCTION G_FLAPS_MID(FLAPS,SB)
C-----------------------------------------------------------------------
C G_FLAPS: COMPUTES GAIN VALUES FOR INBOARD AND MID FLAPS
C-----------------------------------------------------------------------
      REAL G_FL_UP_POS_MID,G_FL_UP_NEG,G_FL_DWN
      INTEGER FLAPS,SB

      G_FL_UP_NEG = (0.0-1.0)/(0.1+0.9)
      G_FL_UP_POS_MID = (0.8-0.0)/(1.9-1.2)
      G_FL_DWN = (1.0-0.0)/(1.9-1.3)

      IF ( FLAPS .EQ. 0 ) THEN
        G_FLAPS_MID = SB*G_FL_UP_POS_MID
      ELSE
        G_FLAPS_MID = SB*G_FL_DWN
      ENDIF

      RETURN
      END
C-----------------------------------------------------------------------

      FUNCTION G_MLA2(FLAPS)
C-----------------------------------------------------------------------
C G_MLA2: COMPUTES MLA2 GAIN VALUES FOR OUTBOARD FLAPS
C-----------------------------------------------------------------------
      INTEGER FLAPS

      IF ( FLAPS .EQ. 0) THEN
        G_MLA2 = (0.5-0.0)/(1.6-1.3)
      ELSE
        G_MLA2 = 0
      ENDIF

      RETURN
      END
C-----------------------------------------------------------------------
}}
