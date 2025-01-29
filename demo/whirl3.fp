
   parameters{cg(CG offset)[-1:1]=0.5, ryscale(torsion scale)[0.01:10]=0.5}

#   parameters{rfrange=1.0-rf}

   # create some new parameters for quasi-steady aero
   parameters{
#     spin1(Prop 1 spin rate)[0:2000]<radps2rpm> = 500
#     spin2(Prop 2 spin rate)[0:2000]<radps2rpm> = 500
      spin1(Prop 1 spin rate)[0:2000]<radps2rpm> = 5
      spin2(Prop 2 spin rate)[0:2000]<radps2rpm> = spin1
      beta1(3/4 blade angle)[34:58]<degrees>=34.0
      beta2(3/4 blade angle)[34:58]<degrees>=34.0
      radius(Prop radius)=0.5
      lbar1(l/r)[0.1:2]=0.3
      lbar2(l/r)[0.1:2]=0.3
#     rf(prop reduced freq) = freq*radius/vtas
   }

   # using the default units: meters, Kg
   fem {
      nodes{1=(0,0,0),2=(0,1,0),3=(0,2,0),4=(0,3,0), 5=(0,4,0), 6=(0,5,0), 7=(0,6,0) }
      beam{1=0, 2=(3:6), xbar=(1,0,0), gj=4.0e+5, eiy=1.0e+5, eiz=2.0e+6, ea=2.0e+5}
      beam{2=(3:6), 3, id=prop1}, beam{3,4}, beam{4,5}, beam{5,6, id=prop2}, beam{6,7}
      mass{mass=5, moi=2.0, cg=0, orient=(1,0,0)), node=(2:7)}
      gyro{node=3{orient=(-1,0,0),moi=237.0, spin=spin1},
         node=6{orient=(-1,0,0), spin=spin2}}
      pgaf{node=3{beta=beta1, radius=radius, lbar=lbar1, orient=(-1,0,0)},
         node=6{beta=beta2, radius=radius, lbar=lbar2, orient=(-1,0,0)}}
      dlm{panels=2, chord=0.4, kval=(0.001,0.1,0.3,5.0), mach=0.5, ac=0}
   }

#  pz{i=gaf.*, o=gaf, code=mypgaf}

   # give the stiffness matrix an evaluation function to change the
   # torsional stiffness
   pz {
      i=stif
      code=stiffcn
      o=stif
   }

   # VSO flutter analysis
# settings{d=2}
# settings{wait}
   flut {
      id=vso
      indep=(vtas[0:200],sigma,freq)
      alt=0
      mass=mass, stif=stif, gyro=gyro, gaf=gaf
      start{modes=(1,2,3,4,5)}
#     start{modes=2}
      target{sigma=0, vtas[10:1000]}
   }
   vz{id=vso, x=vtas, y=sigma, w=rfrange}
   vz{id=vso, x=sigma, y=freq, w=rfrange}

   # search for flutter with vsoe
#  flut {
#     id=vsoe, source=vso
#     indep=(vtas,sigma,freq,spin1,spin2,beta)
#     target{sigma=0, vtas[20:1000]}
#     start{vtas=(40,80,120,200)}
#     optimize{sigma=0.9,spin1=0.1,spin2=0.1}
#  }
#  vz{id=vsoe, x=spin1, y=sigma, w=rfrange}
# end

   # spin variation
   flut {
      id=spin, source=vso
      indep=(spin1[5:500], vtas, freq)
      sigma=0
      start{ordinal=1, vtas[10:1000]}
   }
   vz{id=spin, x=spin1, y=vtas, w=rfrange}

   # mass cg variation
   flut {
      id=cg, source=vso
      indep=(cg, vtas, freq)
      start{ordinal=1, vtas[10:1000]}
      sigma=0
		checklooping
   }
   vz{id=cg, x=cg, y=vtas, w=rfrange}
   
   # torsion variation
   flut {
      id=torsion, source=vso
      indep=(ryscale, vtas, freq)
      start{ordinal=1, vtas[10:1000]}
      sigma=0
   }
   vz{id=torsion, x=ryscale, y=vtas, w=rfrange}

   # spin-cg optimization: start pts for contours
   flut {
      id=optim, source=vso
      indep=(spin1, cg, vtas, freq)
      sigma=0
      optimize=vtas
      start{ordinal=1, vtas[10:1000]}
   }
#  vz {id=optim, x=cg, y=vtas, w=rfrange}

   # cg-torsion contours
   flut {
      id=contour, source=optim
      indep=(spin1, cg, freq)
      vtas=(30,40,50,60,70)
      sigma=0
		checklooping
   }
   vz {id=(contour,optim), x=cg, y=spin1, w=rfrange}
end

   # torsion-cg optimization: start pts for contours
   flut {
      id=optim, source=vso
      indep=(ryscale, cg, vtas, freq)
      sigma=0
      optimize=vtas
      start{ordinal=1, vtas[10:1000]}
   }
#  vz {id=optim, x=cg, y=vtas, w=rfrange}

   # cg-torsion contours
   flut {
      id=contour, source=optim
      indep=(ryscale, cg, freq)
      vtas=(30,50,80,100,120)
      sigma=0
		checklooping
   }
   vz {id=(contour,optim), x=cg, y=ryscale, w=rfrange}
end

stiffcn.cpp {{

int
stiffcn (pset& plt, int nr, int nc, CMatrix& ca) {
   Trace trc(2,"stiffcn");

   if (ca.nel() != nr*nc)
      throw runtime_error(vastr("stiffcn: (",nr,',',nc,") got ca.nel() ",ca.nel()));

   Real ryscale = plt.parval("ryscale");
	if (ryscale.value() != 1.0) {
		for (int j=3; j<24; j += 4)
			for (int i=1; i<=nr; i++)
				ca(i,j) *= ryscale;
	}

   return 0;
}
}}

mypgaf.cpp {{
class myprop {
public:
   vector<int> rc;  // 1-based row,col numbers
   string beta;
   string lbar;
   string radius;
   vector<double> T;
   myprop(vector<int> rows, string be , string lb, string ra, vector<double> tr)
 : rc(rows),beta(be),lbar(lb),radius(ra), T(tr) {};
};
int mypgaf (pset& plt, int nr, int nc, CMatrix& ca) {
// unsteady aero for spinning props
   vector<myprop> props{{{10,11,12},"beta1","lbar1","radius",{-1,0,0,0,0,1,0,1,-0}
},{{28,29,30},"beta2","lbar2","radius",{-1,0,0,0,0,1,0,1,-0}
}};
   Real sigma = plt.parval("sigma");
   Real freq = plt.parval("freq");
   Real vtas = plt.parval("vtas");
   // grab rfrange so that it gets into the dependson for dmatrix
   Real rfrange = plt.parval("rfrange");
   // watch out for zero vtas: return zero ca
   if (vtas.value() <= 0.0)
      return 0;
   static int visit{0};
   if (visit == 0)
      Dprint::dprintm(nr,nc,nr,ca,"pgaf_before");
   // compute the (3,3) gaf matrix for each prop
   for (auto& pi : props) {
      Real beta = plt.parval(pi.beta);
      Real cmq, czr, czpsi, cztheta, cmpsi;
      whirl_deriv(beta, cmq, czr, czpsi, cztheta, cmpsi);
      Real lbar = plt.parval(pi.lbar);
      Real radius = plt.parval(pi.radius);
      Complex sbar(sigma*radius/vtas, freq*radius/vtas);
      Complex t = flaps::pi*radius*radius*radius*(1.0 - lbar*sbar);
      // (3,3) prop gaf matrix with the shaft along the local x axis
      // so the only non-zeros are ry, rz
      cadarray el(3,3);
      el(2,2) = t*(2.0*cmq*sbar - lbar*cztheta);
      el(2,3) = t*(2.0*cmpsi - lbar*(czpsi + czr*sbar));
      el(3,2) = -el(2,3);
      el(3,3) = el(2,2);
      // transform to model coordinates
      //!! ad_triple_product(3,&pi.T[0],el);
      // insert this (3,3) into the return matrix ca
      for (int i=0; i<3; i++)
         for (int j=0; j<3; j++)
            ca(pi.rc[i],pi.rc[j]) += el(i+1,j+1);
   }
   if (visit++ == 0)
      Dprint::dprintm(nr,nc,nr,ca,"pgaf");
   return 0;
}
}}
