
# bif.fp: Tests bifucation.
# Similar to the case used in the paper
# Edward E. Meyer, "Continuation and Bifurcation in Linear Flutter Equations",
# AIAA Journal, V. 53, No. 10, 2015
#
# Demonstrates the use of the determinant of the Jacobian matrix augmented
# with the tangent to the curve to monitor bifurcation. The "bif(urcation)?"
# option is required to do this monitoring.
# If the sign of the determinant changes the default behavior is to retry
# the step with a smaller stepsize to try to avoid bifurcation, provided
# the stepsize is not too "small"

# Goland wing from references 1 and 2 section 9.3.
# References
#   1) @article{goland1945flutter,
#      title={The flutter of a uniform cantilever wing},
#      author={Goland, Martin},
#      year={1945},
#      publisher={American Society of Mechanical Engineers}
#    }
#   2) @book{bisplinghoff1955aeroelasticity,
#        title={Aeroelasticity},
#        author={Bisplinghoff, Raymond L and Ashley, Holt and Halfman, Robert L},
#        year={1955},
#        publisher={Addison-Wesley},
#        address={Reading, Mass.}
#    }
# Cantilevered straight wing: demonstrates beam modeling and doublet-lattice
# unsteady aerodynamics
   parameters{
      cg(CG offset)[-0.5:0.5]=0.183
      ac(Aero center)[-0.5:0.5]=-0.152
   }

   # using the default units: meters, Kg
   # The ac is adjusted from the value in goland.fp (-0.152) to give bifurcation in one
   # or both of modes 1 and 2.
   # The abs value of ac needs to be just small enough that one or both modes smoothly
   # return to  vtas=0; too large and mode 2 flutters, mode 1 reaches 250 (stable),
   # too small and mode 1 flutters, mode 2 is stable. To speed up this trial-and-error
   # approach, with "rcond:" printed in src/bin/flut/bifurcation.cpp:bif_tan(), use a script like
   #    flaps bif >out 2>err
   #    grep here bif.fp >>rcond
   #    grep rcond: err >>rcond
   #    cp rcond /tmp
   #    cp rcond rcond.bak
   #    cat rcond
   # to get bifurcation in either (or both) mode 1 or 2
   fem {
      nodes{1{0,0,0} : 13{0,6,0}}
      material{x=(1,0,0), gj=9.876e+5, eiz=3.0e+7, eix=9.773e+6, ea=2.0e+5}
      beam{1{0}, 2{3:5}, 3:13}
      conmass{mass=18.14, moi=(1,4.43300,0), cg=0.183, orient=x, node=(2:13)}
  dlm{panels=4, chord=1.83, rf=(0.001,0.1,0.3,2.0), mach=0.5, ac=-0.1462388559165075} # adjust ac here
      nmodes=10
   }
	catalog{}

   # VSO flutter analysis using shorthand "vvca" and setting Mach limits
   # "vvca=x" is shorthand for "indep=(rf,freq,sigma), alt=x, which
   # tracks in the limits of rf at altitude "x"
   flut {
      id=vso
      vvca=0, mach[0:0.8]
      mass=mass.modal, stif=stif.modal, gaf=gaf.modal,
      start{modes=(1,2)}
      bifurcation
      target{sigma=0, vtas[10:1000]}
		vzid = modal
   }
   # show 3 x-y combinations ("views")
   viz{id=vso, v=(sigma:freq, vtas:sigma, vtas:freq), uf=fma.modal.uf}
end
