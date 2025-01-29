
# Demonstration of whirl flutter: a rigid propeller and
# shaft, aligned with the x axis, supported at node 1 but
# free to rotate about the y and z axes, spinning at node 2.
# Aero is treated in a function (gaffcn.cpp) that computes the
# component of propeller thrust in the y-z plane; the
# x-axis component is along the spin axis and is irrelevant.
# Reference:
#  [1] Reed III, Wilmer H and Bland, Samuel R, "An analytical
#      treatment of aircraft propeller precession instability",
#      NASA-TN D-659, 1961
# 
# There are two modes in this problem: a forward procession mode
# at around -8 Hz and a backward-precession mode at around 6 Hz.

   # change frequency units to rpm for precession
#  parameters{freq(Precession rate)[0:1000]<radps2rpm>=0}

   parameters{amom(angular momentum)[0:100000] = 2484.666 }

   # create some new parameters for quasi-steady aero
   parameters{
      beta(3/4 blade angle)[34:58]<degrees>=34.0,
      radius(Prop radius)=2.0574,
      lbar(l/r)[0.1:2]=0.377797,
      rf(Prop rf)[0:10]=freq*radius/vtas
      rsf(Prop rsf)=sigma*radius/vtas
      omegaz(get desc from reed)<0.1591549 Hz/(r/s)> = 1.23,
      lambda(Freq ratio)=freq/omegaz
   }
   # create a gyro matrix at unit angular momentum and parameterized by "amom"
   # using data from [1] pg. 30
   #  175 slug-ft^2 = 237.2681325 kg-m^2
   fem{
      nodes{1=(0,0,0),2=(0,1,0)}
      material{x=(1,0,0), gj=4.0e+5, eix=1.0e+5, eiz=1.0e+5, ea=2.0e+6, 1=0, 2=(5:6)}
      beam{1=0, 2=(5,6)}
      conmass{mass=500, moi=200.0, orient=(1,0,0)), node=2}
      gyro{node=2{orient=(-1,0,0), amom=amom}}
      pgaf{node=2{beta=beta, radius=radius, lbar=lbar, orient=(-1,0,0)}}
   }

   # display the gaf matrix and derivatives
   display {gyro.nodal}
   display{
      matview
      alt=0, vtas=10, freq=2, sigma=0, deriv=(0,vtas,freq,sigma)
      pgaf.nodal
   }
   # ratio of z/y stiffness
   parameters{kratio(z/y stiffness factor)[0.1:8.0]=1}
   pz {o=Stif, i=stif.nodal, [2,2] *= kratio }

   # varying velocity at a constant amom
   flut {
      id=vso,
      indep=(vtas[0:100], freq[-50:50], sigma)
      amom = 2484.666
      alt=0
      mass=mass.nodal, stif=Stif
      gyro=gyro.nodal
      gaf=pgaf.nodal
      gcnorm = 0
      target{sigma=0, vtas[20:100]}
      print=(vtas,freq,sigma,ev1.real,ev1.imag,ev2.real,ev2.imag)
      plot=homotopy
      vzid=nodal
   }
   viz{id=vso, x=vtas, y=sigma, w=rfrange}
   viz{id=vso, x=sigma, y=freq, w=rfrange}

   # amom variation
   flut {
      id=amom, source=vso
      indep=(amom[0:100000], freq, vtas)
      start{vtas[50:100]}
      sigma=0
      print=(vtas,freq,sigma,ev1.real,ev1.imag,ev2.real,ev2.imag)
   }
   viz {id=amom, x=amom, y=vtas}

   # 3/4 blade angle variation
   flut {
      id=beta, source=vso
      indep=(beta, freq, vtas)
      sigma=0
      start{vtas[50:100]}
      print=(beta,freq,vtas,ev1.real,ev1.imag,ev2.real,ev2.imag)
   }
   viz {id=beta, x=beta, y=vtas}

   # l/r (lbar) variation: confirms "Gimbal axes location" (p. 21 Reed)
   # increasing the shaft length is stabilizing, but leads to divergence
   flut {
      id=lbar, source=vso
      indep=(lbar[0.01:3], freq, vtas)
      sigma=0
      beta=34
      start{vtas[50:100]}
      print=(lbar,freq,vtas,ev1.real,ev1.imag,ev2.real,ev2.imag)
   }
   viz {id=lbar, x=lbar, y=vtas}
   viz {id=lbar, x=lbar, y=freq}

   # ratio of pitch to yaw stiffness variation
   flut {
      id=kratio, source=vso
      indep=(kratio, freq, vtas[10:100])
      sigma=0
      beta=34
      lbar=0.377797
      start{vtas[50:100]}
      print=(vtas,freq,kratio,ev1.real,ev1.imag,ev2.real,ev2.imag)
   }
   viz {id=kratio, x=kratio, y=vtas}

   end
