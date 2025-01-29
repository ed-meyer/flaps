# Demonstration of model reduction to a set of free-vibration modes
# using octlab to compute the modes. Starting with a model of the
# Goland wing built with 24 beam elements and 96 doublet-lattice panels,
# a flutter solution is done in physical dof. Then the mass, stiffness,
# and aero matrices are sent to matlab/octave where the free-vibration
# eigensolution is done, and the matrices reduced by triple-products with
# the eigenvectors. Using the reduced matrices the flutter solution is
# repeated and the results compared with the physical solution.

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
      nodes{ 1{0,0,0} : 49{0,6.096,0} }
      material{x=(1,0,0), gj=9.876e+5, eiz=3.0e+7, eix=9.773e+6, ea=2.0e+5}
      beam{1{0}, 2{3:5}, 3:49}
      conmass{totalmass=217.74, totalmoi=(1,52.68,0), cg=0.183, orient=x, node=(1:49)}
      dlm{panels=4, chord=1.83, rf=(0.001,0.1,0.3,2.0), mach=0.5, ac=-0.152}
      nmodes = 10
   }

   catalog{summary}


   # check first 2 natural frequencies: should be 50 & 87 r/s = 7.96 & 13.84 Hz
   octlab{i=(roots.m, mass.nodal{cg=0}, stif.nodal), o=roots}
   display{roots}

   # VSO flutter analysis: should flutter at 175.7 m/s, 10.5 Hz
	# Note the bifurction option is necessary
   # flutter solution with the model in nodal dof
   flut {
      id=nodal
      indep=(vtas[0:250],sigma,freq[0:500])
      alt=0
      mass=mass.nodal
      stif=stif.nodal
      gaf=gaf.nodal
      start{modes=(1:30)}
      target{sigma=0, vtas[10:1000]}
      vzid = nodal
   }

   # specify 3 x-y combos with the "v" option
   viz{id=nodal, v=(vtas:sigma,vtas:freq,sigma:freq), title=Nodal dof}

   # VSO flutter analysis in generalized coordinate dof
   flut {
      id=gc
      indep=(vtas[0:250],sigma,freq[0:150])
      alt=0
      mass=mass.modal
      stif=stif.modal
      gaf=gaf.modal
      start{modes=(1:10)}
      target{sigma=0, vtas[10:1000]}
		vzid = modal
   }
   viz{id=gc, v=(vtas:sigma,vtas:freq,sigma:freq), title=Generalized Coord}

   # reduce to free-vibration gc in octlab to compare with the
   # reduction done in fem
   octlab {octave, i=reduce.m, i=(.*.nodal, gct.nodal)
      o=(gmass,gstif,gct__octlab,gaf1,gaf2,gaf3,gaf4)}

   # interpolate the reduced gaf matrices
   pz { i =(gaf1{rf=0.001}, gaf2{rf=0.1}, gaf3{rf=0.3}, gaf4{rf=2}), o=dlm }

   # VSO flutter analysis in generalized coordinate dof
   flut {
      id=octlab
      indep=(vtas[0:250],sigma,freq[0:150])
      alt=0
      mass=gmass
      stif=gstif
      gaf=dlm
      start{modes=(1:10)}
      target{sigma=0, vtas[10:1000]}
      vzid = octlab
   }
   viz{id=octlab, x=vtas, y=sigma}
   # compare nodal, gc and octlab
   viz{id=(gc,nodal,octlab), x=vtas, y=sigma, w=rfrange,
      title='gc-nodal-octlab comparison'}
   viz{id=(gc,nodal,octlab), x=vtas, y=freq, w=rfrange,
      title='gc-nodal-octlab comparison'}
   viz{id=(gc,nodal,octlab), x=sigma, y=freq, w=rfrange,
      title='gc-nodal-octlab comparison'}

end

reduce.m {{
nmodes = 10;
stif = sparse(stif__nodal);
mass = sparse(mass__nodal);
[v,D] = eigs(stif,mass,nmodes,0);
gmass = v'*mass*v;
gstif = v'*stif*v;
gct__octlab = gct__nodal*v;
gaf1 = v'*gaf1__nodal*v;
gaf2 = v'*gaf2__nodal*v;
gaf3 = v'*gaf3__nodal*v;
gaf4 = v'*gaf4__nodal*v;
}}

roots.m {{
stif = stif__nodal;
mass = mass__nodal;
[v,d] = eigs(stif,mass,2,0);
roots = sqrt(diag(d))
}}
