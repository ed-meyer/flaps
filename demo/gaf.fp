#-------------------------------------------------------------   
# Demonstration of interpolation wrt 3 parameters: gaf(rsf,rf,mach)
#
# Description:
#  straight unswept beam model wing with cg variable between
#  -4 inches and 24 inches
#
# Checks:
#   1) matrix parameterizations:
#      - interpolation of a set of 45 gaf matrices with respect
#        to 1, 2 and 3 parameters
#   3) vso with each type of interpolation

# Beam model:                                       ------- y
#                o - lumped masses (cg=8)           |
#                + - structural nodes               |
#                                                   x
#        |
#        |
#   0    |---------+---------+---------+---------+---------+---------+
#        |
#   8    |         o         o         o         o         o         o
#        0        100       200       300       400       500       600
#
# Dublat grid:
#                                                                     1/4 chord at -31
#        |                                                                 /
# -36    |-----------------------------------------------------------|    /
#        |     |     |     |     |     |     |     |     |     |     | <-/
# -16    |-----------------------------------------------------------|
#        |     |     |     |     |     |     |     |     |     |     | <- and -11
#   4    |-----------------------------------------------------------|
#        |
#        0     60   120   180   240   300   360   420   480   540   600
#-------------------------------------------------------------
   # use USCS units for presentation
   settings{units=uscs}

   import{${FROOT}/demo/contour.op4}
   import{${FROOT}/demo/contour.uf}
   import{${FROOT}/demo/gaf.op4}
   catalog{all, summary}
#------------------------------------------------------------------
# parameterize the unsteady generalized airforce matrices
# 45 gaf matrices dependent on (rsf,rf,mach):
#   rsf:  3 values (-0.005, 0, 0.01)
#   rf:  5 values (0, 0.001, 0.002, 0.005, 0.01)
#   mach: 3 values (0, 051, 0.99)
#------------------------------------------------------------------

   # rf interpolation
   pz{
      i = gaf17{rsf=0,rf=0,mach=0.51}
      i = gaf20{rsf=0,rf=0.001,mach=0.51}
      i = gaf23{rsf=0,rf=0.002,mach=0.51}
      i = gaf26{rsf=0,rf=0.005,mach=0.51}
      i = gaf29{rsf=0,rf=0.01,mach=0.51}
      o=gf1
      plot{diagonals, indep=rf}
      units=uscs
   }
   # Display gf1, it's derivative with respect to rf, and the original
   # gaf matrices above and below rf=0.0015/.0254 = 0.059055 (SI) to
   # check the derivative
   display{gaf20, gf1{rf=0.059055}, gaf23}
   display{gf1{rf=0.059055, deriv=rf}}

   # rf-rsf interpolation at mach=0.51
   pz{
      i = gaf16{rsf=-0.005,rf=0,mach=0.51}
      i = gaf17{rsf=0,rf=0,mach=0.51}
      i = gaf18{rsf=0.01,rf=0,mach=0.51}
      i = gaf19{rsf=-0.005,rf=0.001,mach=0.51}
      i = gaf20{rsf=0,rf=0.001,mach=0.51}
      i = gaf21{rsf=0.01,rf=0.001,mach=0.51}
      i = gaf22{rsf=-0.005,rf=0.002,mach=0.51}
      i = gaf23{rsf=0,rf=0.002,mach=0.51}
      i = gaf24{rsf=0.01,rf=0.002,mach=0.51}
      i = gaf25{rsf=-0.005,rf=0.005,mach=0.51}
      i = gaf26{rsf=0,rf=0.005,mach=0.51}
      i = gaf27{rsf=0.01,rf=0.005,mach=0.51}
      i = gaf28{rsf=-0.005,rf=0.01,mach=0.51}
      i = gaf29{rsf=0,rf=0.01,mach=0.51}
      i = gaf30{rsf=0.01,rf=0.01,mach=0.51}
      o=gf2a
      plot{diag, indep=rf, rsf=0}
      units=uscs
   }

# rf-mach interpolation
   pz{
      i = gaf2{rsf=0,rf=0,mach=0}
      i = gaf5{rsf=0,rf=0.001,mach=0}
      i = gaf8{rsf=0,rf=0.002,mach=0}
      i = gaf11{rsf=0,rf=0.005,mach=0}
      i = gaf14{rsf=0,rf=0.01,mach=0}
      i = gaf17{rsf=0,rf=0,mach=0.51}
      i = gaf20{rsf=0,rf=0.001,mach=0.51}
      i = gaf23{rsf=0,rf=0.002,mach=0.51}
      i = gaf26{rsf=0,rf=0.005,mach=0.51}
      i = gaf29{rsf=0,rf=0.01,mach=0.51}
      i = gaf32{rsf=0,rf=0,mach=0.99}
      i = gaf35{rsf=0,rf=0.001,mach=0.99}
      i = gaf38{rsf=0,rf=0.002,mach=0.99}
      i = gaf41{rsf=0,rf=0.005,mach=0.99}
      i = gaf44{rsf=0,rf=0.01,mach=0.99}
      o=gf2b
      plot{diag, indep=rf, mach=0.25}
      units=uscs
   }

   # rsf-rf-mach interpolation
   pz{
      i = gaf1{rsf=-0.005,rf=0,mach=0}
      i = gaf2{rsf=0,rf=0,mach=0}
      i = gaf3{rsf=0.01,rf=0,mach=0}
      i = gaf4{rsf=-0.005,rf=0.001,mach=0}
      i = gaf5{rsf=0,rf=0.001,mach=0}
      i = gaf6{rsf=0.01,rf=0.001,mach=0}
      i = gaf7{rsf=-0.005,rf=0.002,mach=0}
      i = gaf8{rsf=0,rf=0.002,mach=0}
      i = gaf9{rsf=0.01,rf=0.002,mach=0}
      i = gaf10{rsf=-0.005,rf=0.005,mach=0}
      i = gaf11{rsf=0,rf=0.005,mach=0}
      i = gaf12{rsf=0.01,rf=0.005,mach=0}
      i = gaf13{rsf=-0.005,rf=0.01,mach=0}
      i = gaf14{rsf=0,rf=0.01,mach=0}
      i = gaf15{rsf=0.01,rf=0.01,mach=0}
      i = gaf16{rsf=-0.005,rf=0,mach=0.51}
      i = gaf17{rsf=0,rf=0,mach=0.51}
      i = gaf18{rsf=0.01,rf=0,mach=0.51}
      i = gaf19{rsf=-0.005,rf=0.001,mach=0.51}
      i = gaf20{rsf=0,rf=0.001,mach=0.51}
      i = gaf21{rsf=0.01,rf=0.001,mach=0.51}
      i = gaf22{rsf=-0.005,rf=0.002,mach=0.51}
      i = gaf23{rsf=0,rf=0.002,mach=0.51}
      i = gaf24{rsf=0.01,rf=0.002,mach=0.51}
      i = gaf25{rsf=-0.005,rf=0.005,mach=0.51}
      i = gaf26{rsf=0,rf=0.005,mach=0.51}
      i = gaf27{rsf=0.01,rf=0.005,mach=0.51}
      i = gaf28{rsf=-0.005,rf=0.01,mach=0.51}
      i = gaf29{rsf=0,rf=0.01,mach=0.51}
      i = gaf30{rsf=0.01,rf=0.01,mach=0.51}
      i = gaf31{rsf=-0.005,rf=0,mach=0.99}
      i = gaf32{rsf=0,rf=0,mach=0.99}
      i = gaf33{rsf=0.01,rf=0,mach=0.99}
      i = gaf34{rsf=-0.005,rf=0.001,mach=0.99}
      i = gaf35{rsf=0,rf=0.001,mach=0.99}
      i = gaf36{rsf=0.01,rf=0.001,mach=0.99}
      i = gaf37{rsf=-0.005,rf=0.002,mach=0.99}
      i = gaf38{rsf=0,rf=0.002,mach=0.99}
      i = gaf39{rsf=0.01,rf=0.002,mach=0.99}
      i = gaf40{rsf=-0.005,rf=0.005,mach=0.99}
      i = gaf41{rsf=0,rf=0.005,mach=0.99}
      i = gaf42{rsf=0.01,rf=0.005,mach=0.99}
      i = gaf43{rsf=-0.005,rf=0.01,mach=0.99}
      i = gaf44{rsf=0,rf=0.01,mach=0.99}
      i = gaf45 {rsf=0.01,rf=0.01,mach=0.99}
      o=gf3
      plot{diag, indep=rf, mach=0.25}
      units=uscs
   }
   # Display gf3 at rf=0.0015/.0254 = 0.059055 (SI), rsf=0, mach=0.6
   # the derivatives with respect to rf, rsf, and mach
   display{gf3{rsf=0, rf=0.059055, mach=0.6, deriv=(0,rsf,rf,mach)}}

   # 1) Flutter with rf interpolation
   flut {
      id=vso1,
      mass=gm2
      stif=stif
      gaf=gf1
      indep=(vtas[0:500],freq,sigma),
      start{mode=1}
      alt=300
      target{dpress=200},
      target{vtas[20:500], sigma=0},
      freevib
   }
   vz {id=vso1, x=vtas, y=sigma}
   vz {id=vso1, x=sigma, y=freq}

   # 2a) again with rf-rsf interpolation
   flut {
      id=vso2a
      mass=gm2
      stif=stif
      gaf=gf2a
      indep=(vtas[0:500],freq,sigma),
      start{mode=1}
      alt=300
      target{dpress=200},
      target{vtas[20:500], sigma=0},
      freevib
   }
   vz {id=vso2a, x=vtas, y=sigma}
   vz {id=vso2a, x=sigma, y=freq}

   # 2b) again with rf-mach interpolation
   flut {
      id=vso2b
      mass=gm2
      stif=stif
      gaf=gf2b
      indep=(vtas[0:500],freq,sigma),
      start{mode=1}
      alt=300
      target{dpress=200},
      target{vtas[20:500], sigma=0},
      freevib
   }
   vz {id=vso2b, x=vtas, y=sigma}
   vz {id=vso2b, x=sigma, y=freq}

   # 3) and again with 3d interpolation
   flut {
      id=vso3
      mass=gm2
      stif=stif
      gaf=gf3
      indep=(vtas[0:500],freq,sigma),
      start{mode=1}
      alt=300
      target{dpress=200},
      target{vtas[20:500], sigma=0},
   }
   vz {id=vso3, x=vtas, y=sigma}
   vz {id=vso3, x=sigma, y=freq}

   # finally, compare all 4 vso's
   vz {id=(vso1,vso2a,vso2b,vso3), x=vtas, y=sigma}
   vz {id=(vso1,vso2a,vso2b,vso3), x=sigma, y=freq}
end
