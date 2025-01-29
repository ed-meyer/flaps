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
   fem {
      nodes{1=(0,0,0),2=(0,0.25,0),3=(0,0.5,0),4=(0,0.75,0), 5=(0,1,0),
			6=(0,1.25,0), 7=(0,1.5,0), 8=(0,1.75,0), 9=(0,2,0),
			10=(0,2.25,0), 11=(0,2.5,0), 12=(0,2.75,0), 13=(0,3,0),
			14=(0,3.25,0), 15=(0,3.5,0), 16=(0,3.75,0), 17=(0,4,0),
			18=(0,4.25,0), 19=(0,4.5,0), 20=(0,4.75,0), 21=(0,5,0),
			22=(0,5.25,0), 23=(0,5.5,0), 24=(0,5.75,0), 25=(0,6,0)}
      beam{1=0, 2=(3:5), x=(1,0,0), gj=9.88e+5, eiz=3.0e+7, eix=9.77e+6, ea=2.0e+5}
      beam{2=(3:5), 3}, beam{3,4}, beam{4,5}, beam{5,6}, beam{6,7},
      beam{7,8}, beam{8,9}, beam{9,10}, beam{10,11}, beam{11,12},
		beam{12,13}, beam{13,14}, beam{14,15}, beam{15,16}, beam{16,17},
		beam{17,18}, beam{18,19}, beam{19,20}, beam{20,21}, beam{21,22},
		beam{22,23}, beam{23,24}, beam{24,25},
      mass{mass=9.07, moi=(1,2.2,0), cg=0.183, orient=x, node=(2:25)}
      dlm{panels=4, chord=1.83, rf=(0.001,0.1,0.3,2.0), mach=0.5, ac=-0.152}
   }

	# flutter solution with the full model in physical dof
   flut {
      id=full
      indep=(vtas[0:250],sigma,freq)
      alt=0
      mass=mass
      stif=stif
      gaf=gaf
      start{modes=(1,2)}
      target{sigma=0, vtas[10:1000]}
   }
   vz{id=full, x=vtas, y=sigma, w=rfrange}
   vz{id=full, x=vtas, y=freq, w=rfrange}
   vz{id=full, x=sigma, y=freq, w=rfrange}

	# reduce to free-vibration gc in octlab
	octlab { i=(reduce.m, mass, stif, gct, gaf1,gaf2,gaf3,gaf4),
		o=(gmass,gstif,gct,gaf1,gaf2,gaf3,gaf4)}

	# interpolate the reduced gaf matrices
	pz {
		i = (gaf1{rf=0.001}, gaf2{rf=0.1}, gaf3{rf=0.3}, gaf4{rf=2})
		o=dlm
	}

# settings{d=2}
   # VSO flutter analysis
   flut {
      id=gc
      indep=(vtas[0:250],sigma,freq)
      alt=0
      mass=gmass
      stif=gstif
      gaf=dlm
      start{modes=(1,2)}
      target{sigma=0, vtas[10:1000]}
   }
# end
   vz{id=gc, x=vtas, y=sigma, w=rfrange}
   vz{id=gc, x=vtas, y=freq, w=rfrange}
   vz{id=gc, x=sigma, y=freq, w=rfrange}

	# compare gc and physical
   vz{id=(gc,full), x=vtas, y=sigma, w=rfrange}
   vz{id=(gc,full), x=vtas, y=freq, w=rfrange}
   vz{id=(gc,full), x=sigma, y=freq, w=rfrange}

end
   # mass cg variation
   flut {
      id=cg, source=vso
      indep=(cg, vtas[10:800], freq)
      start{vtas[10:500]}
      sigma=0
   }
   vz{id=cg, x=cg, y=vtas, w=rfrange}

   # ac variation
   flut {
      id=ac, source=vso
      indep=(vtas,freq,ac)
      sigma=0
      start{vtas[10:1000]}
   }
   vz{id=ac, x=ac, y=vtas, w=rfrange}

   # ac-cg optimization: start pts for contours
   flut {
      id=optim, source=vso
      indep=(ac, cg, vtas[65:300], freq)
      sigma=0
      optimize{vtas=1}
      start{vtas[50:300]}
   }
#  vz {id=optim, x=cg, y=vtas, w=rfrange}

   # ac-cg contours
   flut {
      id=contour, source=optim
      indep=(ac, cg, freq)
      vtas=(120,140,160,180,200)
      sigma=0
      checklooping
   }
   vz {id=(contour,optim), x=cg, y=ac, w=rfrange}

end

reduce.m {{
ext = [1:10];
[v,D] = eig(stif,mass)
t = v(:,ext);
gmass = t'*mass*t;
gstif = t'*stif*t;
gct = gct*t;
gct = gct(:,ext);
gaf1 = t'*gaf1*t;
gaf2 = t'*gaf2*t;
gaf3 = t'*gaf3*t;
gaf4 = t'*gaf4*t;
}}
