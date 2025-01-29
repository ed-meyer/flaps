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
   # use USCS units for presentation
   settings{units=uscs}

   import{${FROOT}/demo/contour.op4, ${FROOT}/demo/contour.uf}

   catalog{all, summary}
   # create a file for xfig of the amviz aero grid
   export{vzid=contour, o=contour.fig}
   #------------------------------------------------------------------
   # interpolate mass matrices at 4 cg locations
   #------------------------------------------------------------------
   pz {
      o=Mass
      i=(gm1{cg(Center Of Gravity)=-4}, gm2{cg=8}, gm3{cg=16}, gm4{cg=24})
   }
   #------------------------------------------------------------------
   # Parameterize the stiffness matrix:
   #  freq2:  new parameter: the frequency of dof 2
   #------------------------------------------------------------------
   parameters{freq2(Second mode freq)[0:10]<radps2Hz> = 2}

   pz {i=stif, o=Stif, [2,2] = "Mass[2,2]*(freq2^2)" }

   display {stif, gaf0}
   display { Stif{cg=8, freq2=2}, matview }
   #------------------------------------------------------------------
   # parameterize the unsteady generalized airforce matrices
   # the "units=uscs" option scales the input matrices and rf values
   # to effectively convert the entire flutter equation to SI units
   #------------------------------------------------------------------
   pz{
      i = gaf0{rf=0}
      i = gaf001{rf=0.001}
      i = gaf002{rf=0.002}
      i = gaf005{rf=0.005}
      i = gaf01{rf=0.01}
      o=gaf
      units = uscs
   }

   #------------------------------------------------------------------
   # Neutral-stability (vso) first mode only, using reduced-frequency
   # instead of vtas; also test various targets
   #------------------------------------------------------------------
   flut {
      id=vso
      mass=Mass,cg=12, freq2=2
      stif=Stif
      gaf=gaf
      indep=(rf,freq,growth[-0.025:0.3])  # USCS
      start{mode=1}
      alt=300    # ft
      target{dpress=100}
      target{vtas[10:500], growth=0}
      freevib
      plot=homotopy
      print=(rf,freq,growth,vtas)
		vzid=fma
   }
   # plot the homotopy process used to find start points
   # along with the vso solution; show 2 x-y combinations
   viz {id=(vso,homotopy.vso), v=(vtas:growth, growth:freq)}
#   viz {id=(vso,homotopy.vso), x=vtas, y[-0.1:0.1]=growth}
#   viz {id=(vso,homotopy.vso), x=growth, y=freq}
   #------------------------------------------------------------------
   # freq2 variation starting from 1st flutter xing (ordinal=1)
   #------------------------------------------------------------------
   flut {
      id=freq2, source=vso
      indep=(freq2,rf,freq)
      growth=0
      start{ordinal=1}
   }
   viz {id=freq2, y=vtas, x=freq2}
   #------------------------------------------------------------------
   # Optimize dpress wrt freq2 and cg to use as start points for contours
   #------------------------------------------------------------------
   flut {
      id=opt, source=vso
      indep=(cg,freq2,freq,dpress[10:100])
      optimize=dpress
      growth=0
      start{ordinal=1}
   }
   viz {id=opt, x=cg, y=freq2}
   #------------------------------------------------------------------
   #  cg-freq2 contours of constant dpress
   #------------------------------------------------------------------
   flut {
      id=contour, source=opt
      indep=(cg,freq2,freq)
      checklooping
      growth=0
      dpress=(10 : 10 : 90)
   }
   # 3D plot the optimization results together with the contours
   viz {id=(opt,contour), x=cg, y=freq2, z=dpress}

   end
