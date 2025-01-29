# time-domain bifurcation

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
   # using the default units: meters, Kg
   fem {
      nodes{1{0,0,0} : 13{0,6,0} }
      material{x=(1,0,0), gj=9.876e+5, eiz=3.0e+7, eix=9.773e+6, ea=2.0e+5}
      beam{1{0}, 2{3:5}, 3 : 13}
      mass{mass=18.14, moi=(1,4.395,0), cg=0.183, orient=x, node=(2:13)}
      dlm{panels=4, chord=1.83, rf=(0.001,0.1,0.3,2.0), mach=0.5, ac=-0.5}
		nmodes = 10
   }
	catalog{}

	# using the optimal betas computed in optimize.fp
   # got final betas = -0.00021754, 0.00435947, 0.119997, 1.13081, 1.67805
 	# pz { i=gaf[1234], o=rfa, beta=(-0.00021754, 0.00435947, 0.119997, 1.13081, 1.67805) }
 	pz { i=gaf[1234], o=rfa, beta=(0.00435947, 0.119997, 1.13081, 1.67805) }

   flut {
      id=interp
      indep=(vtas[0:250],sigma,freq)
      alt=0
      mass=mass
      stif=stif
      gaf=gaf
      start{modes=(1,2)}
      target{sigma=0, vtas[10:250]}
   }
   flut {
      id=initial
      indep=(vtas[0:250],sigma,freq)
      alt=0
      mass=mass
      stif=stif
      gaf=rfa
      start{modes=(1,2)}
      target{sigma=0, vtas[10:250]}
   }
#  vz{id=(interp,rfa), x=vtas, y=sigma, w=rfrange}

	td {
      id=vso
      indep=(vtas[100:150], freq[0:100], sigma)
		maxstep = 21
      gcnorm = 0
      alt=0
      mass=mass, stif=stif, gaf=rfa
	}
end
