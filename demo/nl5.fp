
#-------------------------------------------------------------
# Demonstrates:
# 1) pk flutter solution
# 2) nacelle side-bending frequency variation
# 3) right-click on v-g curve brings up amvis
#
# Modes:
#  1-6    rigid-body
#  7-60   elastic                 -> 1-54
#  61-66  6 nacelle modes         -> 55-60
#  67-71  rudder modes            -> 61-65
#  72-76  right elevator          -> 66-70
#  77-81  left elevator           -> 71-75
#  82-83  stabilizer pitch/roll   -> 76-77
#  84-88  APU
#
#-------------------------------------------------------------

   import { ${FROOT}/demo/nl5.op4}
   import { ${FROOT}/demo/nl5.uf}

# Create gaf interpolation coef
   pz {
      units=uscs
      i = gaf0{rf=0},
      i = gaf001{rf=0.001},
      i = gaf005{rf=0.005},
      i = gaf01{rf=0.01},
      o = gaf,
   }

# vso with linear stiffness
#ifdef NEVER
   flut {
      id=vso0
      start{vtas=10, freq[0:5]},
      mass=mass, stif=stif, gaf=gaf,
      indep=(vtas[10:350],sigma,freq),
      alt=0,
      target{growth=0}
#     plot=stepsize
   }
   vz {id=vso0, x=vtas, y=sigma}
   vz {id=vso0, x=sigma, y=freq}
#endif

# parameterize the stiffness matrix: nacelle sB frequencies
   parameters{
      gap61(Rudder gap) = 0.01,
      gap66(Rudder gap) = 0.01,
      gap71(Rudder gap) = 0.01
   }

   pz {
      o=stif,
      i=stif,
      [61,61] *= freeplay(gap61, 61)
      [66,66] *= freeplay(gap66, 66)
      [71,71] *= freeplay(gap71, 71)
   }

   display{stif}
   display{mass}

# vso at gcnorm=0
   flut {
      id=vso,
#     start{vtas=10, freq[0:50]},
      start{vtas=10, mode=13}
      mass=mass, stif=stif, gaf=gaf,
      indep=(vtas[10:350],sigma,freq),
      alt=0,
      target{sigma=0}
      freevib
   }
# end
   vz {id=vso, x=vtas, y=sigma}
   vz {id=vso, x=sigma, y=freq}

#ifdef NEVER
# voe starting from vso
   flut {
      id=voe0, source=vso
      indep=(gcnorm[0:1], freq, vtas[10:350])
      sigma=0
      start{vtas[20:350], ordinal=1}
   }
   vz {id=voe0, x=vtas, y=gcnorm, w=lcostab}
# uncertainty voe starting from vso
   flut {
      id=voehilo, source=vso
      indep=(gcnorm[0:1], freq, vtas[10:350])
      sigma=0
      start{vtas[20:350], ordinal=1}
      gap61=(0.0090,0.0110), gap66=(0.0090,0.0110), gap71=(0.0090,0.0110)
   }
   vz {id=voehilo, x=vtas, y=gcnorm}
#endif // NEVER

# soe starting from vso
   flut {
      id=soe, source=vso
      indep=(sigma, freq, gcnorm[0:1])
      target{sigma=0}
#     vtas=(100,200,300)
      vtas=50
   }
   vz {id=soe, x=gcnorm, y=sigma}
# voe starting from soe
   flut {
      id=voe1, source=soe
      indep=(vtas[10:350], freq, gcnorm[0:1])
      sigma=0
      start{ordinal=1}
      checklooping
   }
   vz{id=voe1, x=vtas, y=gcnorm, w=lcostab}
# gap uncertainty
   flut {
      id=voehilo, source=soe
      indep=(vtas[10:350], freq, gcnorm[0:1])
      sigma=0
      start{ordinal=1}
      checklooping
      gap61=(0.0099,0.0101), gap66=(0.0099,0.0101), gap71=(0.0099,0.0101)
   }
   vz {id=voehilo, x=vtas, y=gcnorm}
end
