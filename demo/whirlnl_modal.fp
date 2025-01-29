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

   # create some new parameters for quasi-steady aero
   parameters{
      spin1(Prop 1 spin rate)[0:3000]<radps2rpm> = 100
      spin2(Prop 2 spin rate)[0:3000]<radps2rpm> = spin1
      beta1(3/4 blade angle)[34:58]<degrees>=34.0
      beta2(3/4 blade angle)[34:58]<degrees>=34.0
      radius(Prop radius)=0.5
#     lbar1(l/r)[0.1:2]=0.1
#     lbar2(l/r)[0.1:2]=0.1
      lbar1(l/r)[0.1:2]=0.5
      lbar2(l/r)[0.1:2]=0.5
   }

   # using the default units: meters, Kg
# settings{d=2}
   fem {
      nodes{1=(0,0,0),2=(0,0.5,0),3=(0,1,0),4=(0,1.5,0), 5=(0,2,0),
         6=(0,2.5,0), 7=(0,3,0), 8=(0,3.5,0), 9=(0,4,0), 10=(0,4.5,0),
         11=(0,5,0), 12=(0,5.5,0), 13=(0,6,0),
         14=(-0.5,2,0), 15=(-0.5,4.5,0)}  # gimbal nodes
      beam{1=0, 2=(1:6), x=(1,0,0), gj=9.876e+5, eiz=3.0e+8, eix=9.773e+6, ea=8.0e+7}
      beam{2,3}, beam{3,4}, beam{4,5}, beam{5,6}, beam{6,7},
      beam{7,8}, beam{8,9}, beam{9,10}, beam{10,11}, beam{11,12}, beam{12,13},
      beam{5=(5,6),14, x=(0,1,0), eiz=2e+6, eix=2e+6, nomerge}, beam{10,15,nomerge}, # nacelles
      mass{mass=18.14, moi=(1,4.395,1), cg=0.183, orient=x, node=(2:13),
         mass=10, moi=(5,5,5), cg=0, orient=y, node=(14,15)}
      gyro{node=14{orient=(-1,0,0),moi=1.0, spin=spin1},
         node=15{orient=(-1,0,0), spin=spin2}}
      pgaf{node=14{beta=beta1, radius=radius, lbar=lbar1, orient=(-1,0,0)},
         node=15{beta=beta2, radius=radius, lbar=lbar2, orient=(-1,0,0)}}
      dlm{panels=4, nodes=(1:13), chord=1.83, rf=(0.001,0.1,0.3,2.0),
         mach=0.5, ac=-0.152}
   }
# end
   # reduce to free-vibration modes + assumed prop modes
# settings{d=2}
   octlab {i=(reduce.m, mass, stif, gyro, gaf[1234], gct),
      o=(gm,gg,gk,gct,gaf1,gaf2,gaf3,gaf4,gimbals,nac1,nac2,a1,a2), octave }
#  display{gm, gg, gk, a1, a2}
   catalog{summary}
   # create a new universal file with the reduced gct
   export{o=gct.uf}
   # give the stiffness matrix a custom evaluation function to
   # add freeplay to the gimbals
	# XXX using 2.35e+4 for gimbals instead of 4e+5 as in whirlnl.fp: factor of 17
   parameters{
      gap(Gimbal freeplay)[0:0.1]=0.001
      gimbal1(Prop 1 stiff) = 2.35e+4
      gimbal2(Prop 2 stiff) = gimbal1
   }
   pz { i=gk, code=stifred, o=gk }

   pz { i=gaf[1234], o=gaf, code=pgafred}

   pz { i=gg, code=gyrored, o=gg }

   # VSO flutter analysis
   flut {
      id=vso
      indep=(vtas[0:350],sigma,freq[0:150])
#   debug=start
      alt=0
      spin1 = 100
      gcnorm = 0.4
      mass=gm, stif=gk, gyro=gg, gaf=gaf
#     start{modes=(1:10)}
      start{modes=(1:5)}
      target{sigma=0, vtas[10:350]}
   }
   vz{id=vso, x=vtas, y=sigma, w=rfrange}
   vz{id=vso, x=sigma, y=freq, w=rfrange}

   # alternative way to do voe at various spin rates:
   # spin-vtas variation, then voe's at various spin
   flut {
      id=spin, source=vso
      indep=(vtas[10:350], freq, spin1[0:3000])
      sigma=0
      start{vtas[10:350]}
   }
   vz {id=spin, x=spin1, y=vtas}
   flut{
      id=voe, source=spin
      indep=(vtas[10:350], freq, gcnorm[0.001:1])
      spin1=(20,200,400,600,800,1000)
      sigma=0
   }
   vz{id=voe, x=vtas, y=gcnorm}
end
   # opt process for spin-gcnorm contours
   flut {
      id=opt, source=vso
      indep=(vtas[61:350],freq,spin1[10:1000],gcnorm[0.002:0.5])
      sigma=0
      optimize{vtas=0.8,gcnorm=0.1}
      start{vtas[20:350]}
   }
   vz{id=opt, x=vtas, y=gcnorm}
   vz{id=opt, x=vtas, y=spin1}

   # spin-gcnorm contours
   flut {
      id=contours, source=opt
      indep=(vtas,gcnorm,freq)
      sigma=0
      spin1=(70,100,120,140)
      start{ordinal=1}
   }
   vz{id=contours, x=vtas, y=gcnorm}
end
   # VOE process
   flut {
      id=voe, source=vso
      indep=(vtas[10:350], freq, gcnorm[0:0.5])
      sigma=0
      start{vtas[20:350]}
   }
   vz{id=voe, x=vtas, y=gcnorm, w=lcostab}

   # spin lco process
   flut {
      id=spinlco, source=vso
      indep=(spin1[50:300], freq, gcnorm[0:0.5])
      sigma=0
      start{vtas[20:350]}
   }
   vz{id=spinlco, x=spin1, y=gcnorm, w=lcostab}
end


reduce.m {{
ext = [1:14];
gimbal1y=73;
gimbal1z=74;
gimbal2y=75;
gimbal2z=76;
nac1y=23;
nac1z=24;
nac2y=53;
nac2z=54;
[v,D] = eig(stif,mass);
T = v(:,ext);
gm = T'*mass*T;
gk = T'*stif*T;
gg = T'*gyro*T;
gaf1 = T'*gaf1*T;
gaf2 = T'*gaf2*T;
gaf3 = T'*gaf3*T;
gaf4 = T'*gaf4*T;
gct = gct*T;
gimbals = T(73:76,:);
nac1 = T(23:24,:);
nac2 = T(53:54,:);
a = zeros(76);
a(gimbal1y,gimbal1y) = 1.0;
a(gimbal1y,nac1y) = -1.0;
a(nac1y,gimbal1y) = -1.0;
a(nac1y,nac1y) = 1.0;

a(gimbal1z,gimbal1z) = 1.0;
a(gimbal1z,nac1z) = -1.0;
a(nac1z,gimbal1z) = -1.0;
a(nac1z,nac1z) = 1.0;

a(gimbal2y,gimbal2y) = 1.0;
a(gimbal2y,nac2y) = -1.0;
a(nac2y,gimbal2y) = -1.0;
a(nac2y,nac2y) = 1.0;

a(gimbal2z,gimbal2z) = 1.0;
a(gimbal2z,nac2z) = -1.0;
a(nac2z,gimbal2z) = -1.0;
a(nac2z,nac2z) = 1.0;
xterms = T'*a*T;
a1 = xterms(1:2,:);
a2 = xterms(3:4,:);
}}

pgafred.cpp {{
#include "fem.h"
int pgafred (pset& plt, int nr, int nc, CMatrix& ca) {
// unsteady aero for spinning props
   Trace trc(1,"pgaffcn");
   vector<Pgaf_elem> props{
   {"radius", 0.5,"lbar1", 0.5,"beta1", 34,{0,1,2},{-1,0,0,0,1,0,0,0,-1}
},{"radius", 0.5,"lbar2", 0.5,"beta2", 34,{0,3,4},{-1,0,0,0,1,0,0,0,-1}
}};
   Real sigma = plt.parval("sigma");
   Real freq = plt.parval("freq");
   Real vtas = plt.parval("vtas", true);
   Real rf = plt.parval("rf", true);
   Real rsf = plt.parval("rsf", true);
   // watch out for zero vtas: return zero ca
   double eps{0.1};
   if (vtas.value() <= eps)
      return 0;
   // compute the (3,3) gaf matrix for each prop
   for (auto& pi : props) {
      Real beta = pi.beta_v(plt);
      Real cmq, czr, czpsi, cztheta, cmpsi;
      whirl_deriv(beta, cmq, czr, czpsi, cztheta, cmpsi);
      Real lbar = pi.lbar_v(plt);
      Real radius = pi.radius_v(plt);
      Complex sbar(rsf, rf);
      Complex t = flaps::pi*radius*radius*radius*(1.0 - lbar*sbar);
      // (3,3) prop gaf matrix with the shaft along the local x axis
      // so the only non-zeros are ry, rz
      cadarray el(3,3);
      el(2,2) = t*(2.0*cmq*sbar - lbar*cztheta);
      el(2,3) = t*(2.0*cmpsi - lbar*(czpsi + czr*sbar));
      el(3,2) = -el(2,3);
      el(3,3) = el(2,2);
      // transform to model coordinates
      ad_triple_product(3,&pi.T[0],el);
      // insert this (3,3) into the return matrix ca
      for (int i=0; i<3; i++) {
         if (pi.rc[i] == 0) continue;
         for (int j=0; j<3; j++) {
            if (pi.rc[j] == 0) continue;
            ca(pi.rc[i],pi.rc[j]) += el(i+1,j+1);
          }
       }
   }
   static int visit{0};
   if (visit++ == 0)
      Dprint::dprintm(nr,nc,nr,ca,"pgaf");
   return 0;
}
}}

gyrored.cpp {{
#include "fem.h"
int gyrored (pset& plt, int nr, int nc, CMatrix& ca) {
// multiply gyro elements by their spin
   vector<Gyro_elem> gyros{{14, "", 200, "spin1", 0, {0,-0,0,0,0,200,-0,-200,0},{0,1,2}}
,{15, "", 200, "spin2", 0, {0,-0,0,0,0,200,-0,-200,0},{0,3,4}}
};
   // multiply each (3,3) gyro matrix by its spin
   for (auto& ge : gyros) {
      Real spin;
      if (ge.spin_p.empty())
         spin = ge.spin;
      else
         spin = plt.parval(ge.spin_p);
      // if there is no moi parameter set it to 1: already applied
      Real moi;
      if (ge.moi_p.empty())
         moi = 1.0;
      else
         moi = plt.parval(ge.moi_p);
      for (int i=0; i<ge.rc.size(); i++) {
         if (ge.rc[i] == 0) continue;
         for (int j=0; j<ge.rc.size(); j++) {
            if (ge.rc[j] == 0) continue;
            ca(ge.rc[i],ge.rc[j]) *= spin*moi;
         }
      }
   }
   static int visit{0};
   if (visit++ == 0)
      Dprint::dprintm(nr,nc,nr,ca,"gyro");
   return 0;
}
}}

stifred.cpp {{
int
stifred (pset& plt, int nr, int nc, CMatrix& ca) {
   Trace trc(2,"stifred");

   if (ca.nel() != nr*nc)
      throw runtime_error(vastr("stifred: (",nr,',',nc,") got ca.nel() ",ca.nel()));
   // 4 extra dof:
   // 73)   14/5  ry node 14 real & complex  <-> 23) 5/5
   // 74)   14/6  rz node 14 real & complex  <-> 24) 5/6
   // 75)   15/5  ry node 15 real & complex  <-> 53) 10/5
   // 76)   15/6  rz node 15 real & complex  <-> 54) 10/6
   
   // freeplay is on the (y,z) rotational dof of nodes 14 (73,74) and
   // 15 (75,76). There is no stiffness between these dof so we add some
   Real k1 = plt.parval("gimbal1");
   Real k2 = plt.parval("gimbal2");
   Real gap = plt.parval("gap");
   vector<int> nacelle1{23, 24};
   vector<int> nacelle2{53, 54};
   vector<int> gim1{1, 2};
   vector<int> gim2{3, 4};
   //!! double smooth{0.0};  // no smoothing
   double smooth{0.05};
   Matrix nac1("a1");
   Matrix nac2("a2");
   double* n1 = nac1.elem();
   double* n2 = nac2.elem();

   assert(nc = (int)nac1.csize());
   assert(nac1.rsize() == 2);

   // apply freeplay
   for (int i=0; i<2; i++) {
      // Real absgc1 = plt.absgcp(gim1[i], nac1[i]);
      Real absgc1 = plt.absgcv(gim1[i]);
      plt.setpar(absgc1,vastr("absgc1.",i));
      Real fp1 = df::freeplay(gap, absgc1, smooth);
      plt.setpar(fp1, vastr("fp1.",i));
      Real absgc2 = plt.absgcv(gim2[i]);
      plt.setpar(absgc2, vastr("absgc2.",i));
      Real fp2 = df::freeplay(gap, absgc2, smooth);
      plt.setpar(fp2, vastr("fp2.",i));
      // trc.dprint("prop1 dof ",gim1[i],", ",nac1[i],", gc=",absgc1.value(),", fp = ",fp1);
      // trc.dprint("prop2 dof ",gim2[i],", ",nac2[i],", gc=",absgc2.value(),", fp = ",fp2);
      // prop1
      for (int j=0; j<nc; j++) {
         ca(gim1[i],j+1) = fp1*k1*n1[IJ(i,j,2)];
         ca(j+1,gim1[i]) = fp1*k1*n1[IJ(i,j,2)];
      }
      // ca(gim1[i],gim1[i]) = fp1*k1;

      // prop2
      for (int j=0; j<nc; j++) {
         ca(gim2[i],j+1) = fp2*k2*n2[IJ(i,j,2)];
         ca(j+1,gim2[i]) = fp2*k2*n2[IJ(i,j,2)];
      }
      // ca(gim2[i],gim2[i]) = fp2*k2;
   }
   trc.dprintm(nr,nc,nr,ca,"stifred ca");

   return 0;
}
}}
