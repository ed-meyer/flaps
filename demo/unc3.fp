
#------------------------------------------------------------------
# generic 2-engine airplane
# freeplay on elevators, rudder, ailerons
# bilinear stiffness on nacelles
# bilinear factor on gaf
# Shows a latent lco
#------------------------------------------------------------------
   settings{units=metric}

#  settings{pagewidth=300}

# nonlinear dofs implemented in stiffcn and gaffcn:
#  41     rudder rotation
#  42     left elevator rotation
#  43     right elevator rotation
#  44     left aileron rotation
#  45     right aileron rotation
#  46-51  nacelles 3 dof/side

   import { ${FROOT}/demo/lco.op4}
   import { ${FROOT}/demo/lco.uf}

#  display {matview, gstif, mass}

# give the stiffness matrix an evaluation function with freeplay and bilinear
   new Par{gap41(mode 41 gap)=0.001}
   new Par{gap42(mode 42 gap)=0.001}
   new Par{gap43(mode 43 gap)=0.001}
   new Par{gap44(mode 44 gap)=0.001}
   new Par{gap45(mode 45 gap)=0.001}
   pz {
      i = gstif,
      o=Stif,
      [41,41] *= freeplay(gap41,41)
      [42,42] *= freeplay(gap42,42)
      [43,43] *= freeplay(gap43,43)
      [44,44] *= freeplay(gap44,44)
      [45,45] *= freeplay(gap45,45)
   }

# gaf matrices are interpolated wrt rf
   pz {
      i=(gaf0{rf=0}, gaf0008{rf=0.0008}, gaf0015{rf=0.0015}
         gaf002{rf=0.002}, gaf003{rf=0.003}, gaf01{rf=0.01}, gaf05{rf=0.05})
      o = pgaf
   }

# VSO at gcnorm=0
# settings{valgrind}
   flut {
      id = vso
      mass=mass
      stif=Stif
      gaf = pgaf
      indep=(vtas[10:300], freq[0:50], sigma)
      alt = 0
      gcnorm = 0.0
#     start{vtas=10, freq[0:10]}
      start{vtas=10, modes=14}
      target{sigma=0,vtas[10:300]}
      freevib
   }
   vz {id=vso, x[-1:1]=sigma, y=freq}
   vz {id=vso, x=vtas, y=sigma}

# voe starting from vso
   flut {
      id=voe0, source=vso
   maxstep=150
   uncert=(gap41{7%}, gap42{7%}, gap43{7%}, gap44{7%}, gap45{7%})
   debug=uncert
      indep=(gcnorm[0:1], freq, vtas[10:500])
      sigma=0
      start{vtas[20:500], ordinal=1},  # first xing only
      checklooping
   }
   vz {id=voe0, x=vtas, y=gcnorm}
end
# soe starting from vso: only mode 8
   flut {
      id=soe, source=vso
      indep=(sigma, freq, gcnorm[0:1])
      vtas=(50,100)
      start{mode=8}
   }
   vz{id=soe, x[0:0.02]=gcnorm, y=sigma}

# lco analysis starting from soe
   flut {
      id=voe1, source=soe
      indep=(gcnorm[0:1], freq, vtas[0:1000])
      sigma=0
      start{ordinal=1},  # first xing only
      checklooping
   }
   vz {id=voe1, x=vtas, y=gcnorm, w=lcostab}
# vsoe starting at 100
   flut {
      id=vsoe, source=vso
      indep=(vtas, sigma[-100:0.001], freq, gcnorm[0:1])
      start{vtas=100, mode=8}
      optimize{sigma=0.9, gcnorm=0.4359}
      target{sigma=0}
   }
   vz{id=vsoe, x=sigma, y=gcnorm}

# voe starting from vsoe
   flut {
      id=voe2, source=vsoe
      indep=(vtas[10:1000], freq, gcnorm)
      sigma=0
      checklooping
   }
   vz {id=(voe2,vsoe), x=vtas, y=gcnorm}
# these are the only gc that go beyond delta
   vz {id=voe2, x=vtas, y='fac4[12345]'}
   vz {id=voe2, x=vtas, y='fac[1234]'}
# check the interpolations of the 3 df's
   vz {x=q, y=f, bilinear@0.001@0.5.apf}
   vz {x=q, y=f, bilinear@0.01@1.2.apf}
   vz {x=q, y=f, gap@0.001.apf}
end


