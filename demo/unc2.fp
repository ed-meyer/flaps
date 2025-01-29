
# small nonlinear problem based on HA145B
# Demonstrates the use of built-in functions for describing
# function approximation of nonlinearities.
#
# This is the first of 4 files demonstrating 4 different
# ways of implementing describing functions. It scales the first
# 4 diagonals of the stiffness matrix by a built-in function
# (bilinear) in a "matrix element" parameterization.
#
# An example of latent LCO is shown in voe2

   import{${FROOT}/demo/ha145b.op4}

# Give the first four diagonals of the stiffness matrix
# bilinear stiffness; args are (delta, ratio, gcno)
   new Par{gap6(mode 6 gap)=0.001}
   new Par{gap7(mode 7 gap)=0.001}
   new Par{gap8(mode 8 gap)=0.001}
   new Par{gap10(mode 10 gap)=0.001}
#   new Par{fpsmooth(freeplay smoothing)=0.1}
   pz {
      i = KHH,
      o=Stif,
#     [3,3] *= freeplay(0.001, 3)
#     [4,4] *= freeplay(0.0001, 4)
#     [5,5] *= freeplay(0.001, 5)
      [6,6] *= freeplay(gap6,6)
      [7,7] *= freeplay(gap7,7)
      [8,8] *= freeplay(gap8,8)
#     [9,9] *= freeplay(0.001,9)
      [10,10] *= freeplay(gap10,10)
   }

   pz { i=MHH, o=Mass,
#     [2,2] *= freeplay(0.001,2)
#     [3,3] *= freeplay(0.0001,3)
#     [4,4] *= freeplay(0.0001,4)
#     [5,5] *= freeplay(0.0001,5)
#     [6,6] *= freeplay(0.0001,6)
#     [7,7] *= freeplay(0.0001,7)
#     [8,8] *= freeplay(0.0001,8)
#     [9,9] *= freeplay(0.0001,9)
      [10,10] *= freeplay(0.0001,10)
   }

   pz {
      out=GAF,
      i=(QHHL1{rf=0.00000001524}, QHHL2{rf=0.00001524}, QHHL3{rf=0.00076200},
         QHHL4{rf=0.001524}, QHHL5{rf=0.003048}, QHHL6{rf=0.00762},
         QHHL7{rf=0.01524}, QHHL7{rf=0.03}),
      beta = (0.0007, 0.007)
   }

# VSO - start points for SOE, VOE analyses
   flut {
      id=vso,
      indep=(vtas[0:1000], freq[0:100], sigma),
#     start{vtas=10, freq[0.1:40]},
    start{mode=7}
      gcnorm = 0,
      alt=-10000.0,
#     mass=Mass,
      mass=MHH,
      stif=Stif,
#     stif=KHH,
      gaf=GAF,
      target{vtas[10:1000], sigma=0},
   }
   vz {id=vso, x=vtas, y=sigma}
# end
# VOE process starting from eta=0, sigma=0 solutions from VSO
# Note: mode 5 takes a different path than in nl2.fp
# settings{d=2}
  flut {
      id=voe0, source=vso,
      indep=(vtas[0:3000], freq, gcnorm[0:1]),
      sigma=0,
      start{vtas[20:1000]},
      uncert=(gap6{1%},gap7{1%},gap8{1%},gap10{1%})
#     uncert=(gap8{0.01%})
#     plot=stepsize
      maxstep = 100
#     maxstep = 160
#     curvaturefac = 7
      debug=uncert
   }
   vz {id=voe0, x=vtas, y=gcnorm}
end

# VOE process at all combinations of the limits of each
# uncertain parameter
   flut {
      id=vsohilo
      indep=(vtas[0:1000], freq[0:100], sigma),
#     start{vtas=10, freq[0.1:40]},
    start{vtas=10, mode=7}
      gcnorm = 0,
      alt=-10000.0,
      gap6=(0.00099,0.00101), gap7=(0.00099,0.00101), gap8=(0.00099,0.00101), gap10=(0.00099,0.00101)
#     mass=Mass,
      mass=MHH,
      stif=Stif,
#     stif=KHH,
      gaf=GAF,
      target{vtas[10:1000], sigma=0},
   }
   vz {id=vsohilo, x=vtas, y=sigma}

# settings{wait}
   flut {
      id=voehilo, source=vsohilo
#  debug=start
      indep=(vtas[0:3000], freq, gcnorm[0:1])
      sigma=0
      gap6=(0.00099,0.00101), gap7=(0.00099,0.00101), gap8=(0.00099,0.00101), gap10=(0.00099,0.00101)
      start{vtas[20:1000]},
      maxstep = 160, curvaturefac = 7
#     plot=stepsize
   }
end
   vz {id=voehilo, x=vtas, y=gcnorm}
end
