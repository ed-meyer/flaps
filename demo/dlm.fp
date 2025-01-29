# Demonstration of some features of the doublet-lattice implementation
# using the Goland wing. 5 different interpolations of the unsteady aerodynamic
# matrices are used in simple VSO flutter calculations and a parameter variation:
# 1) rf interpolation
# 2) rf-rsf interpolation
# 3) rf-mach interpolation
# 4) rf-rsf-mach interpolation
# 5) rf-ac interpolation with an ac-variation showing the effect of ac
#    location on flutter
# Cantilevered straight wing: demonstrates beam modeling and doublet-lattice
# unsteady aerodynamics
# Reference
#   1) @article{goland1945flutter,
#      title={The flutter of a uniform cantilever wing},
#      author={Goland, Martin},
#      year={1945},
#      publisher={American Society of Mechanical Engineers}
#    }
# using the default units: meters, Kg

   # 1) rf interpolation
   fem {
      nodes{1{0,0,0} : 13{0,6,0} }
      material{x=(1,0,0), gj=9.876e+5, eiz=3.0e+7, eix=9.773e+6, ea=2.0e+5}
      beam{1{0}, 2{3:5}, 3 : 13}
      conmass{mass=18.14, moi=(1,4,0), cg=0.183, orient=x, node=(2:13)}
      dlm{panels=4, chord=1.83, rf=(0.001,0.1,0.3,2.0), rsf=0, mach=0.5, ac=-0.152}
      nmodes = 5
   }

   catalog{}

   # Interpolate gaf wrt rf
   pz {i= (gaf1.modal, gaf2.modal, gaf3.modal, gaf4.modal)
      o=gaf,
      plot{row=(1:3),col=(1:3)}
   }

   # VSO flutter analysis using rf interpolation
   flut {
      id=rf
      indep=(vtas[0:250],sigma,freq)
      alt=0
      mass=mass.modal
      stif=stif.modal
      gaf=gaf
      start{modes=(1:2)}
      target{sigma=0, vtas[10:1000]}
      vzid = modal
   }
   viz{id=rf, x=sigma, y=freq, w=rfrange}

   # 2) repeat with rf-rsf
   fem {
      nodes{1{0,0,0} : 13{0,6,0} }
      material{x=(1,0,0), gj=9.876e+5, eiz=3.0e+7, eix=9.773e+6, ea=2.0e+5}
      beam{1{0}, 2{3:5}, 3 : 13}
      conmass{mass=18.14, moi=(1,4,0), cg=0.183, orient=x, node=(2:13)}
      dlm{panels=4, chord=1.83, rf=(0.001,0.1,0.3,2.0), rsf=(-1,0,1), mach=0.5, ac=-0.152}
      nmodes = 5
   }

   flut {
      id=rf-rsf
      indep=(vtas[0:250],sigma,freq)
      alt=0
      mass=mass.modal
      stif=stif.modal
      gaf=gaf.modal
      start{modes=(1:2)}
      target{sigma=0, vtas[10:1000]}
      vzid = modal
   }
   viz{id=rf-rsf, x=sigma, y=freq, w=rfrange}

   # 3) repeat with rf-mach
   fem {
      nodes{1{0,0,0} : 13{0,6,0} }
      material{x=(1,0,0), gj=9.876e+5, eiz=3.0e+7, eix=9.773e+6, ea=2.0e+5}
      beam{1{0}, 2{3:5}, 3 : 13}
      conmass{mass=18.14, moi=(1,4,0), cg=0.183, orient=x, node=(2:13)}
      dlm{panels=4, chord=1.83, rf=(0.001,0.1,0.3,2.0), rsf=0, mach=(0,0.5,0.9), ac=-0.152}
      nmodes = 5
   }

   flut {
      id=rf-mach
      indep=(vtas[0:250],sigma,freq)
      alt=0
      mass=mass.modal
      stif=stif.modal
      gaf=gaf.modal
      start{modes=(1:2)}
      target{sigma=0, vtas[10:1000]}
      vzid = modal
   }
   viz{id=rf-mach, x=sigma, y=freq, w=rfrange}

   # 4) repeat with rf-rsf-mach
   fem {
      nodes{1{0,0,0} : 13{0,6,0} }
      material{x=(1,0,0), gj=9.876e+5, eiz=3.0e+7, eix=9.773e+6, ea=2.0e+5}
      beam{1{0}, 2{3:5}, 3 : 13}
      conmass{mass=18.14, moi=(1,4,0), cg=0.183, orient=x, node=(2:13)}
      dlm{panels=4, chord=1.83, rf=(0.001,0.1,0.3,2.0), rsf=(-1,0,1), mach=(0,0.5,0.9), ac=-0.152}
      nmodes = 5
   }

   flut {
      id=rf-rsf-mach
      indep=(vtas[0:250],sigma,freq)
      alt=0
      mass=mass.modal
      stif=stif.modal
      gaf=gaf.modal
      start{modes=(1:2)}
      target{sigma=0, vtas[10:1000]}
      vzid = modal
   }
   viz{id=(rf,rf-rsf,rf-mach,rf-rsf-mach), x=sigma, y=freq, w=rfrange}

end
