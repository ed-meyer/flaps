
# Demonstration of whirl flutter: 2 rigid propellers with
# shafts aligned with the z axis, attached to nodes 3 and 6.
# The shafts are rigid so the 2 dof that are relevant to the gyro
# and gaf terms are rx and ry at nodes 3 and 6.
# Aero is treated in a function (evtol_gaffcn) that computes the
# component of propeller thrust in the y-z plane; the
# x-axis component is along the spin axis and is irrelevant.
# Reference:
#  [1] Reed III, Wilmer H and Bland, Samuel R, "An analytical
#      treatment of aircraft propeller precession instability",
#      NASA-TN D-659, 1961

   # create some new parameters for quasi-steady aero
   parameters{
		spin1(Prop 1 spin rate)[0:2000]<radps2rpm> = 500
		spin2(Prop 2 spin rate)[0:2000]<radps2rpm> = 500
      beta1(3/4 blade angle)[34:58]<degrees>=34.0
      beta2(3/4 blade angle)[34:58]<degrees>=34.0
      radius(Prop radius)=1.0
      lbar1(l/r)[0.1:2]=1.0
      lbar2(l/r)[0.1:2]=1.0
#     ar(Advance ratio)=vtas/(2.0*spin*radius),
#     omegaz(get desc from reed)<0.1591549 Hz/(r/s)> = 1.23,
#     lambda(Freq ratio)=freq/omegaz
      rf(prop reduced freq) = freq*radius/vtas
   }

  # create a beam model with 2 propellers.
  # The beam is along the global x axis; local coordinates for beam
  # elements are the y-axis along the beam and local x (xbar) specified by
  # a point in the global x-positive y plane
# settings{wait}
  fem{
#   nodes{1=(0,0,0),2=(.57735,.57735,.57735),3=(1.1547,1.1547,1.1547),4=(1.73205,1.73205,1.73205),5=(2.3094,2.3094,2.3094),6=(2.886751,2.886751,2.886751),7=(3.4641,3.4641,3.4641) }
    nodes{1=(0,0,0),2=(1,0,0),3=(2,0,0),4=(3,0,0),5=(4,0,0),6=(5,0,0),7=(6,0,0) }
    beam{1=0, 2=(1:6), xbar=(0,-1,0), gj=1.0e+6, eiy=8.0e+6, eiz=2.0e+5, ea=16.0e+6}
    beam{2=(1:6), 3, id=prop1}, beam{3,4}, beam{4,5}, beam{5,6, id=prop2}, beam{6,7}
    mass{mass=100, moi=100, orient=y, node=(2:7)}
    gyro{node=3{orient=prop1(0,1,0),moi=237.0, spin=spin1},
	      node=6{orient=prop2(0,1,0), spin=spin2}}
    pgaf{node=3{beta=beta1, radius=radius, lbar=lbar1, ydir=(0,0,1)},
        node=6{beta=beta2, radius=radius, lbar=lbar2, ydir=(0,0,1)}}
  }
# end

   # create a matrix "gaf" which is populated in function "my_gaffcn"
# settings{d=2}
  # pz {i=evtol_gaf, code=evtol_gaffcn}

	# gyro matrix
	# pz {i=evtol_gyro, code=gyrofcn}

   # create a viscous damping matrix
#  parameters{vdfac(Viscous damping factor)[0:0.01]=0.001}
#  pz{i=evtol_stif, o=vdamp, [1,1] *= vdfac, [2,2] = 0}
#  display{vdamp}

   # varying velocity at a constant spin
# settings{wait}
# settings{d=2}
   flut {
      id=vso,
      indep=(vtas[0:400], freq[0:1000], sigma)
      alt=0
   start{vtas=0, mode=2}
#  start{vtas=0, freq[0:10]}
      mass=mass, stif=stif, gyro=gyro, gaf=pgaf
#     vdamp=vdamp
      target{sigma=0.00001, vtas[40:400]}
   }
   vz{id=vso, x=vtas, y=sigma}
   vz{id=vso, x=sigma, y=freq}

#ifdef NEVER
   # spin1 variation
   flut {
      id=spin1, source=vso
      indep=(spin1, freq, vtas)
      sigma=0
      start{vtas[50:100]}
   }

   # spin2 variation
   flut {
      id=spin2, source=vso
      indep=(spin2, freq, vtas)
      sigma=0
      start{vtas[50:100]}
   }
   vz{id=(spin1,spin2), x=vtas,y=(spin1,spin2)}

	# both the same
   flut {
      id=spin, source=vso
      indep=(spin2, freq, vtas)
		spin1=spin2
      sigma=0
      start{vtas[50:100]}
   }
   vz{id=spin, x=vtas,y=spin2}
#endif // NEVER

   # spin1-spin2 vtas optimization: start points for spin_contours
   flut {
      id=spinopt, source=vso
      indep=(spin1, spin2, freq, vtas)
      optimize{vtas=1}, sigma=0
      start{vtas[40:400], ordinal=1}
   }
   vz {id=spinopt, x=vtas, y=spin1, y=spin2}

   # spin1-spin2 vtas contours
   flut {
		id=spin_contours, source=spinopt
		indep=(spin1,spin2,freq)
		sigma=0
		vtas=(60,70,80,90,100)
   }
   vz {id=(spinopt,spin_contours), x=spin1, y=spin2}

   # lbar1-lbar2 vtas optimization: start points for lbar_contours
   flut {
      id=lbaropt, source=vso
      indep=(lbar1, lbar2, freq, vtas)
      optimize{vtas=1}, sigma=0
      start{vtas[40:400], ordinal=1}
   }
   vz {id=lbaropt, x=vtas, y=lbar1, y=lbar2}

   # lbar1-lbar2 vtas contours
   # l/r (lbar) variation: confirms "Gimbal axes location" (p. 21 Reed)
   # increasing the shaft length is stabilizing, but leads to divergence
   flut {
		id=lbar_contours, source=lbaropt
		indep=(lbar1,lbar2,freq)
		sigma=0
		vtas=(50,100,200,300,350)
   }
   vz {id=(lbaropt,lbar_contours), x=lbar1, y=lbar2}
end

   # 3/4 blade angle variation
   flut {
      id=beta, source=vso
      indep=(beta1, freq, vtas)
		beta2 = beta1
      sigma=0
      start{vtas[50:700]}
   }
   vz {id=beta, x=beta1, y=vtas}
end

   # beta1-beta2 vtas optimization: start points for beta_contours
   flut {
      id=betaopt, source=vso
      indep=(beta1, beta2, freq, vtas)
      optimize{vtas=1}, sigma=0
      start{vtas[40:500], ordinal=1}
   }
   vz {id=betaopt, x=vtas, y=beta1, y=beta2}

   # beta1-beta2 vtas contours
   flut {
		id=beta_contours, source=betaopt
		indep=(beta1,beta2,freq)
		sigma=0
		vtas=(130, 150, 180, 220)
   }
   vz {id=(betaopt,beta_contours), x=beta1, y=beta2}
end
   # l/r (lbar) variation: confirms "Gimbal axes location" (p. 21 Reed)
   # increasing the shaft length is stabilizing, but leads to divergence
   flut {
      id=lbar, source=vso
      indep=(lbar[0.01:3], freq, vtas)
      sigma=0
      beta=34
      start{vtas[50:700]}
   }
   vz {id=lbar, x=lbar, y=vtas}
   vz {id=lbar, x=lbar, y=freq}

#ifdef NEVER
   # stiffness variation
   flut {
      id=opt, source=vso
      indep=(sfx, sfy, freq, vtas)
      optimize=vtas
      sigma=0
      beta=34
      lbar=0.377797
   }
   vz {id=opt, x=sfy, y=vtas}
   # velocity contours varying sfx and sfy
   flut {
      id=contour, source=opt
      indep=(sfx,sfy,freq)
#   vtas=(40,50,60,70,80,90,100)
      vtas=(50 : 50 : 700)
      checklooping
   }
   vz {id=(contour,opt), x=sfx, y=sfy}

   # freeplay on the gimbals
   flut {
      id=voe, source=vso
      indep=(gcnorm[0:1], freq, vtas)
      sigma=0
      sfx = 1
      sfy = 1
      beta=34
      lbar=0.377797
   }
   vz {id=voe, x=vtas, y=gcnorm, w=lcostab}

   # soe at 50 m/s, 100 rpm
   flut {
      id=soe, source=vso
      indep=(sigma, freq, gcnorm[0:1])
      vtas=50
      spin = 100
      sfx = 1
      sfy = 1
      beta=34
      lbar=0.377797
   }
   vz {id=soe, x=gcnorm, y=sigma}
#endif // NEVER

   # vsoe at a few vtas
#ifdef NEVER
   flut {
      id=vsoe, source=vso
      indep=(beta, freq, gcnorm[0:1], lbar, spin, sigma, sfx, sfy, vtas, gapx, gapy)
      start{vtas=(20,30,40,50,60)}
      optimize{sigma=0.9,gcnorm=.435}
   }
   vz {id=vsoe, x=gcnorm, y=sigma}
#endif // NEVER

   end
