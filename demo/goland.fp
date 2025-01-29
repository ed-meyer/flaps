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
# unsteady aerodynamics with 2D interpolation (rf-ac)
   parameters{
      cg(CG offset)[-3:3]=0.183
      ac(Aero center)[-2:2]=-0.152
   }

   # using the default units: meters, Kg
   fem {
      nodes{1{0,0,0} : 13{0,6.096,0} }
      material{x=(1,0,0), gj=9.876e+5, eiz=3.0e+7, eix=9.773e+6, ea=2.0e+5}
      beam{1{0}, 2{3:5}, 3 : 13}
      conmass{totalmass=217.74, totalmoi=(1,52.68,0), cg=cg, orient=x, node=(1:13)}
      dlm{panels=4, chord=1.83, rf=(0.001,0.1,0.3,2.0), mach=0.5, ac=(-2,-0.152,2)}
   }

   # check first 2 natural frequencies: should be 50 & 87 r/s = 7.96 & 13.84 Hz
   octlab{i=(freevib.m, mass.nodal{cg=0}, stif.nodal), o=roots}
   display{roots}

   # VSO flutter analysis: should flutter at 175.7 m/s, 10.5 Hz
   flut {
      id=vso
      indep=(vtas[0:250],sigma,freq)
		bifurcation
      alt=0
      ac=-0.152
      cg = 0.183
      mass=mass.nodal
      stif=stif.nodal
      gaf=gaf.nodal
      freevibration
      start{modes=(1,2)}
      target{sigma=0, vtas[10:250]}
      vzid = nodal
   }
#  vz{id=vso, x=vtas, y=sigma, w=rfrange}
#  vz{id=vso, x=vtas, y=freq, w=rfrange}
#  vz{id=vso, x=sigma, y=freq, w=rfrange}
   # specify 3 x-y combos with the "v" option
   viz{id=vso, v=(vtas:sigma, vtas:freq, sigma:freq)}

   # mass cg variation
   flut {
      id=cg, source=vso
      indep=(cg, vtas[10:800], freq)
      start{vtas[10:500]}
      sigma=0
      vzid = nodal
   }
   viz{id=cg, x=cg, y=vtas, w=rfrange}

   # ac variation
   flut {
      id=ac, source=vso
      indep=(vtas,freq,ac)
      sigma=0
      start{vtas[10:250]}
      vzid = nodal
   }
   viz{id=ac, x=ac, y=vtas, w=rfrange}

   # ac-cg optimization: start pts for contours
   flut {
      id=optim, source=vso
      indep=(ac, cg, vtas[65:500], freq)
      sigma=0
      optimize=vtas
      start{vtas[65:300]}
      vzid = nodal
   }
#  viz {id=optim, x=cg, y=vtas, w=rfrange}

   # ac-cg contours
   flut {
      id=contour, source=optim
      indep=(ac, cg, freq)
      vtas=(125,140,160,180,200,210,220,230,240,250)
      sigma=0
      checklooping
      vzid = nodal
   }
   viz {id=(contour,optim), x=cg, y=ac, w=rfrange}

end

freevib.m {{
stif = stif__nodal;
mass = mass__nodal;
[v,d] = eigs(stif,mass,2,0);
roots = sqrt(diag(d));
}}
