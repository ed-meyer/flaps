
#-------------------------------------------------------------   
# Demonstration of parameter contours
#
# Description:
#  straight unswept beam model wing with cg variable between
#  -4 inches and 24 inches
#
# Checks:
#   1) matrix parameterizations:
#      - interpolation of 4 mass matrices with cg locations
#      - frequency of mode 2 set by an equation
#   3) simple pk flutter
#   4) mode 2 frequency parameter variation
#   5) optimization wrt mode 2 frequency and cg
#   6) mode 2 freq - cg contours

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

   import{${FROOT}/demo/contour.op4}
   import{${FROOT}/demo/contour.uf}
   catalog{all, summary}
# create a file for xfig of the amvz aero grid
    export{o=contour.fig}
#------------------------------------------------------------------
# interpolate mass matrices at 4 cg locations
#------------------------------------------------------------------
   pz {
      o=Mass,
      i=(gm1{cg(Center Of Gravity)=-4}, gm2{cg=8}, gm3{cg=16}, gm4{cg=24}),
   }
#------------------------------------------------------------------
# Parameterize the stiffness matrix:
#  freq2:  new parameter: the frequency of dof 2
#------------------------------------------------------------------
   new Par{freq2(Second mode freq)[0:10]<HZPRS> = 2}
   pz {
      i=stif,
      o=Stif,
      [2,2] = "Mass[2,2]*(freq2^2)"
   }

   display {stif, gm1}
#  display { Stif{cg=8, freq2=2}, matview }
#------------------------------------------------------------------
# parameterize the unsteady generalized airforce matrices
#------------------------------------------------------------------
   pz{
      i = gaf0{rf=0},
      i = gaf001{rf=0.001},
      i = gaf002{rf=0.002},
      i = gaf005{rf=0.005},
      i = gaf01{rf=0.01},
      o=gaf
   }

#------------------------------------------------------------------
# Neutral-stability (vso) first mode only
#------------------------------------------------------------------
settings{valgrind}
   flut {
      id=vso
      mass=Mass,cg=12, freq2=2
      stif=Stif
      gaf=gaf
      indep=(vtas[0:500],freq,sigma)
      unc = (freq2{1%}, cg{2%})
      start{mode=1}
      alt=300
      target{dpress=200}
      target{vtas[10:500], sigma=0}
      freevib
   }
   vz {id=vso, x=vtas, y[-0.1:0.1]=growth}
   vz {id=vso, x=sigma, y=freq}
#------------------------------------------------------------------
# freq2 variation starting from 1st flutter xing (ordinal=1)
# Use start{vtas[10:100]} to get only the first xing. Also could
# do the same with start{ordinal=1}
#------------------------------------------------------------------
settings{wait}
   flut {
      id=freq2, source=vso
      indep=(freq2,vtas[30:300],freq)
      sigma=0
      start{ordinal=1}
      unc = (freq2{1%}, cg{2%})
   }
end
   vz {id=freq2, y=vtas, x=freq2}
end
#------------------------------------------------------------------
# Optimize vtas wrt mode 2 frequency and cg
#------------------------------------------------------------------
   flut {
      id=opt, source=vso
      indep=(cg,freq2,freq,vtas[40:200])
      optimize=vtas,
      sigma=0
   }
   vz {id=opt, x=cg, y=freq2}
#------------------------------------------------------------------
#  cg-Mode 2 contours
#------------------------------------------------------------------
   flut {
      id=contour, source=opt
      indep=(cg,freq2,freq)
      checklooping
      sigma=0
      vtas=(40 : 10 : 150)
   }
# plot the optimization results together with the contours...
   vz {id=(opt,contour), color, x=cg, y=freq2}

   end
