# Demonstration of a variable alt-constant-Mach (vacm) flutter solution.
# Using the Goland wing built with 24 beam elements and 96 doublet-lattice panels,
# and reduced with a free-vibration solution in fem.

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
   # A refined version of the Goland wing
   fem {
      nodes{ 1{0,0,0} : 49{0,6,0} }
      material{x=(1,0,0), gj=9.876e+5, eiz=3.0e+7, eix=9.773e+6, ea=2.0e+5}
      beam{1{0}, 2{3:5}, 3:49}
      conmass{mass=4.53, moi=(0.5,1.1,0), cg=0.183, orient=x, node=(2:49)}
      dlm{panels=4, chord=1.83, rf=(0.001,0.1,0.3,2.0), mach=0.5, ac=-0.152}
      nmodes = 10
   }

   catalog{summary}

   # create a viscous damping matrix
   octlab {octave, i=vdamp.m, i=stif.modal, o=vdamp}
   display{vdamp}

   # vacm VSO flutter analysis in generalized coordinate dof
   flut {
      id=vacm
      vacm=0.5		 # equivalent to indep=(alt,sigma,freq[0:150]), mach=0.5
      mass=mass.modal, stif=stif.modal, gaf=gaf.modal, vdamp=vdamp
      start{modes=(1:10)}
      target{sigma=0, veas[10:1000]}
      vzid = modal
   }
   viz{id=vacm, v=(veas:sigma, sigma:freq)}
end

vdamp.m {{
vdfac = 0.0005;
% note stif.modal has been replaced by stif__modal: necessary
% because Matlab/octave would take stif.modal as a struct. This
% replacement is done in Matlab::exporter and undone in Matlab::importer
vdamp = stif__modal*vdfac;
}}
