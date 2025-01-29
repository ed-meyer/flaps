#-------------------------------------------------------------   
# simple LTI control-law using the Goland wing
#
# Checks:
#   1) ABCD parameterization
#   2) pk flutter solutions at several gain values
#   3) gain variation - one gain for both input channels (igain)
#-------------------------------------------------------------

   # using the default units: meters, Kg
   fem {
      nodes{1{0,0,0} : 13{0,6,0} }
      material{x=(1,0,0), gj=9.876e+5, eiz=3.0e+7, eix=9.773e+6, ea=2.0e+5}
      beam{1{0}, 2{3:5}, 3 : 13}
      conmass{mass=18.14, moi=(1,4.395,0), cg=0.183, orient=x, node=(2:13)}
      dlm{panels=4, chord=1.83, rf=(0.001,0.1,0.3,2.0), mach=0.5, ac=-0.152}
      nmodes = 10
   }
   catalog{}
   display{modes.modal}
#------------------------------------------------------------------
# create a psi matrix: extract rows from the modes matrix
# corresponding to node 13 tz & ry
#------------------------------------------------------------------
   octlab {octave, i=extractpsi.m, i=modes.modal, o=tip}
   display {tip}

   # gain parameters
   parameters{
      gain(input gain)[0:5]=1
      phase(input phase)[-3.14:3.14]=0
   }
   # create T
   lti{file=lti.dat,psi=tip, stif=stif.modal, e=(1{1},2{-1},3{0.4}),
      igain=(gain,gain), iphase=(phase,phase)
      ogain=(1,1,1), ophase=(0,0,0)
      smooth = 0.01
      o=abcd
   }
   catalog{}

   display{abcd{alt=100000,sigma=1,freq=1,derivs=(sigma,freq,alt)}}
   display{B,C,D}

   # VSO flutter analysis
   flut {
      id=vso
      indep=(rf,sigma,freq)
      alt=0
      mass=mass.modal
      stif=stif.modal
      gaf=gaf.modal
      controls = abcd
#     start{modes=(1,2)}
      start{modes=1}
      target{sigma=0, vtas[10:250]}
   plot=stepsize
   pac{maxsteps=40,curvaturefac=6}
   }
   viz{id=vso, x=vtas, y=sigma, w=rfrange}
   viz{id=vso, x=vtas, y=freq, w=rfrange}
   viz{id=vso, x=sigma, y=freq, w=rfrange}

   # gain variation
   flut {
      id=gain, source=vso
      indep=(gain,freq,vtas)
      sigma=0
      start{vtas[20:500]}
   }
   viz{id=gain, x=gain, y=vtas}
   # phase variation
   flut {
      id=phase, source=vso
      indep=(phase,freq,vtas)
      sigma=0
      start{vtas[20:500]}
   }
   viz{id=phase, x=phase, y=vtas}
   end


lti.dat {{
# ns   ni   no   nb   nix   nin
   3    2    3    3    1    8
# Independent Variable Points (alt Breakpoints)
altitude
-20000.0
100000.0
250000.0
# interpolated elements of A:
1 1
# constant elements of A:
1 2
1 3
2 1
2 2
2 3
3 1
3 2
3 3
# input time delays
0.1
0.0
# output time delays
0.2
0.0
0.0
# S:
0.0
2.0
# A: (ns,ns) = (3,3)
11.0
12.0
13.0
21.0
22.0
23.0
31.0
32.0
33.0
# B: (ns,ni) = (3,2)
0.11
0.12
0.21
0.22
0.31
0.32
# C: (no,ns) = (3,3)
0.01
0.0
0.0
0.0
0.02
0.0
0.0
0.0
0.03
# D: (no,ni) = (3,2)
0.0
0.0
0.0
0.0
0.0
0.0
# A coeff: (nb-1)*4*nix = 8: A[1,1] = 1000 + (alt+20000)
# (-20000:100000):
0.0
0.0
1.0
1000.0
# (100000:250000): 121000 - alt
0.0
0.0
-1.0
121000.0
# number of internal time delays:
1
# number of elements multiplied by this internal td
2
# elements that get this td
1 1
2 2
# the td
0.02
}}

extractpsi.m {{
% get the tz, ry displacements: rows 34, 36 (1b)
tip = [modes__modal(34,:); modes__modal(36,:)]
}}

