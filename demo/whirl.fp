# Whirl flutter and variation of some parameters.
# The model is the Goland wing with 2 propellers added

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

#  parameters{cg(CG offset)[-1:1]=0.1, ryscale(torsion scale)[0.01:10]=1.0}

   # create some new parameters for quasi-steady propeller aero
   parameters{
      amom1(Prop 1 amom)[0:5000] = 2000
      amom2(Prop 2 amom)[0:5000] = amom1
      beta1(3/4 blade angle)[34:58]<deg>=45.0
      beta2(3/4 blade angle)[34:58]<deg>=45.0
      radius(Prop radius)=0.5
      lbar1(l/r)[0.1:2]=0.5
      lbar2(l/r)[0.1:2]=0.5
   }

   # A refined version of the Goland wing with 24 beam elements,
   # (4,24) dlm panels. Using the default units: meters, Kg
   # Branch-mode formulation with the wing as the root substructure,
   # and 2 branches: the 2 nacelles and propellers. Nonlinearities are freeplay
   # on nacelle pitch (dof 22,24)
   # 1-20  wing vibration modes with propellers rigidly attached
   # 21,22 propeller 1 yaw and pitch
   # 23,24 propeller 2 yaw and pitch
   fem {
      nodes{1{0,0,0} : 25{0,6,0}, 109{-0.5,2,0}, 119{-0.5,4.5,0}}
      ss{
         id=wing
         material{x=(1,0,0), gj=9.88e+5, eiz=3.0e+7, eix=9.77e+6, ea=8.0e+7}
         beam{1{0}, 2{1:6}, 3 : 25}
         conmass{mass=9.07, moi=(1,2.2,1), cg=0.183, orient=x, node=(2:25)}
      }
      ss{
         id=prop1
         material{x=(0,1,0), eiz=8e+4, eix=8e+4, gj=1.0e+7, ea=8.0e+7}
         beam{9{1:6},109}
         conmass{mass=10, moi=(5,5,5), cg=0, orient=y, node=109}
         conmass{mass=9.07, moi=(1,2.2,1), cg=0.183, orient=x, node=9}
         project=(2,3)
      }
      ss{
         id=prop2
         material{x=(0,1,0), eiz=8e+4, eix=8e+4, gj=1.0e+7, ea=8.0e+7}
         beam{19{1:6},119}
         conmass{mass=10, moi=(5,5,5), cg=0, orient=y, node=119}
         conmass{mass=9.07, moi=(1,2.2,1), cg=0.183, orient=x, node=19}
         project=(2,3)
      }
      bma{root=wing{20}, branches=(prop1{2}, prop2{2})}
      cms{ss=(wing{20}, prop1{2}, prop2{2}) }
      gyro{node=109{orient=(-1,0,0), amom=amom1},
           node=119{orient=(-1,0,0), amom=amom2}}
      pgaf{node=109{beta=beta1, radius=radius, lbar=lbar1, orient=(-1,0,0)}
           node=119{beta=beta2, radius=radius, lbar=lbar2, orient=(-1,0,0)}}
      dlm{nodes=(1:25), panels=4, chord=1.83, rf=(0.001,0.1,0.3,2.0,5.0), mach=0.5, ac=-0.152}
   }
# end
   catalog{summary}

   # VSO flutter analysis on the full fem model
   flut {
      id=vsop
      indep=(rf, sigma, freq[0:150])
      mach[:1], print=(rf,vtas,sigma,freq,mach)
      alt=0
      amom1 = 2000
      mass=mass.nodal, stif=stif.nodal, gyro=gyro.nodal, gaf=gaf.nodal
      start{modes=(1:8)}
      target{sigma=0, vtas[10:300]}
      vzid=nodal
   }

   # frequencies of the nacelles, set to match the branch freq from fem
   # for comparison with vsop (see output from fem)
   parameters{ratio(Pitch/yaw freq ratio)[0.1:10] = 1}
   parameters { freq21(Prop 1 yaw)[1:80] = 26.7853
                freq22(Prop 1 pitch)[1:80] = ratio*freq21
                freq23(Prop 2 yaw)[1:80] = freq21
                freq24(Prop 2 pitch)[1:80] = freq22 }
   pz {i=stif.bma, code=stiffcn}

   # vso with the branch-mode reduced model
   flut {
      id=vso
      indep=(rf, sigma, freq[0:150])
      mach[:1], print=(rf,vtas,sigma,freq,mach)
      alt=0
      amom1 = 2000
      mass=mass.bma, stif=stif.bma, gaf=gaf.bma, gyro=gyro.bma
      start{modes=(1:8)}
      target{sigma=0, vtas[10:300]}
      vzid=bma
   }
   viz{id=vso, x=vtas, y=sigma}

   # compare nodal v branch modal
   viz{id=(vsop,vso), x=vtas, y[-0.5:0.5]=sigma, w=rfrange}
   viz{id=(vsop,vso), color, x=sigma, y=freq}

   # amom variation
   flut {
      id=amom, source=vso
      indep=(amom1, vtas, freq)
      sigma=0
      start{ordinal=1, vtas[10:300]}
   }
   viz{id=amom, x=amom1, y=vtas, w=rfrange}

   flut {
      id=ratio, source=vso
      indep=(ratio,vtas,freq)
      sigma=0
      start{ordinal=1, vtas[10:300]}
   }
   viz{id=ratio, x=ratio, y=vtas}

   # amom1-amom2 vtas optimization: start points for vtas contours
   flut {
      id=optim, source=vso
      indep=(amom1, amom2, vtas, freq)
      sigma=0
      optimize=vtas
      start{ordinal=1, vtas[10:300]}
   }
   viz {id=optim, x=amom1, y=vtas, w=rfrange}

   # amom1-amom2 contours
   flut {
      id=contour, source=optim
      indep=(amom1, amom2, freq)
      vtas=(54,80,100,120,150,200,300)
      sigma=0
      checklooping
   }
   viz {id=(contour,optim), c=2, x=amom1, y=amom2, w=rfrange}
   viz {id=(contour,optim), c=3, x=amom1, y=amom2, w=rfrange}

   # look for vtas-amom combinations causing flutter
   flut {
      id=vsopt, source=vso
      indep=(amom1,vtas,freq,sigma)
      optimize=sigma
      start{ordinal=1,vtas=(100,200)}
   }
   viz{id=vsopt, x=amom1, y=sigma}

   # track the flutter crossings found in vamom
   flut {
      id=vamom, source=vsopt
      indep=(amom1,vtas,freq)
      sigma=0
      start{vtas[20:300]}
   }
   viz {id=vamom, x=amom1, y=vtas}
end

stiffcn.cpp {{

int
stiffcn (pset& plt, int nr, int nc, vector<complex<Ad>>& ca) {
   Trace trc(2,"stiffcn");

   if (ca.size() != nr*nc)
      throw runtime_error(vastr("gkfcn: (",nr,',',nc,") got ca.size() ",ca.size()));

   // 4 propeller dof:
   // 21)   propeller 1 yaw
   // 22)   propeller 1 pitch
   // 23)   propeller 2 yaw
   // 24)   propeller 2 pitch

   // freeplay is on the pitch branch modes for the props: 22, 24
   Ad f21 = plt.parval("freq21");
   Ad f22 = plt.parval("freq22");
   Ad f23 = plt.parval("freq23");
   Ad f24 = plt.parval("freq24");
   // Ad gap1 = plt.parval("p1gap");
   // Ad gap2 = plt.parval("p2gap");
   // Ad fp21 = df::freeplay(plt, gap1, 21);
   // Ad fp21{1.0};
   // Ad fp22 = df::freeplay(plt, gap1, 22);
   // Ad fp23 = df::freeplay(plt, gap2, 23);
   // Ad fp23{1.0};
   // Ad fp24 = df::freeplay(plt, gap2, 24);
   // Ad eqf21 = sqrt(fp21)*f1;
   // Ad eqf22 = sqrt(fp22)*f22;
   // Ad eqf23 = sqrt(fp23)*f2;
   // Ad eqf24 = sqrt(fp24)*f24;

   // plt.setpar(fp21,"fp21");
   // plt.setpar(fp22,"fp22");
   // plt.setpar(fp23,"fp23");
   // plt.setpar(fp24,"fp24");

   // plt.setpar(eqf21,"eqf21");
   // plt.setpar(eqf22,"eqf22");
   // plt.setpar(eqf23,"eqf23");
   // plt.setpar(eqf24,"eqf24");

   // convert f21, f22, f23, f24 to radians
   f21 *= 2.0*flaps::pi;
   f22 *= 2.0*flaps::pi;
   f23 *= 2.0*flaps::pi;
   f24 *= 2.0*flaps::pi;

   // branches are mass-normalized to 1.0
   // ca indices are 0b so (21,21) -> [20*(nr+1)]
   // ca[20*(nr+1)] = f21*f21*fp21;
   // ca[21*(nr+1)] = f22*f22*fp22;
   // ca[22*(nr+1)] = f23*f23*fp23;
   // ca[23*(nr+1)] = f24*f24*fp24;
   ca[20*(nr+1)] = f21*f21;
   ca[21*(nr+1)] = f22*f22;
   ca[22*(nr+1)] = f23*f23;
   ca[23*(nr+1)] = f24*f24;
   trc.dprintm(nr,nc,nr,ca,"stiffcn ca");

   return 0;
}
}}

