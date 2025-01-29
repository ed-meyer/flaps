
# small nonlinear problem based on HA145B
# Demonstrates the use of built-in functions for describing
# function approximation of nonlinearities.
#
# This is the first of 4 files demonstrating 4 different
# ways of implementing describing functions. It scales the first
# 4 diagonals of the stiffness matrix by a built-in function
# (bilinear) in a "matrix element" parameterization.
#
# An example of latent LCO is shown in voe2

   import{${FROOT}/demo/ha145b.op4}

# Give the first four diagonals of the stiffness matrix
# bilinear stiffness; args are (delta, ratio, gcno)
   new Par{massfac=1}
   new Par{massfac2=0.001}
   new Par{massfac3=0.001}
   new Par{massfac4=0.001}

   pz { i=MHH, o=Mass,
      code = massfcn
   }

   new Par{stiffac=1}
   pz {
      i = KHH, o=Stif,
      code = stiffcn
   }
   display{Stif}

   pz {
      out=GAF,
      i=(QHHL1{rf=0.00000001524}, QHHL2{rf=0.00001524}, QHHL3{rf=0.00076200},
         QHHL4{rf=0.001524}, QHHL5{rf=0.003048}, QHHL6{rf=0.00762},
         QHHL7{rf=0.01524}, QHHL7{rf=0.03}),
      beta = (0.0007, 0.007)
   }

# VSO - start points for SOE, VOE analyses
   flut {
      id=vso,
      indep=(vtas[0:1000], freq[0:100], sigma),
#     start{vtas=10, freq[0.1:40]},
      start{vtas=10, mode=4}
      gcnorm = 0,
      alt=-10000.0,
      mass=Mass, stif=Stif, gaf=GAF,
      target{vtas[10:1000], sigma=0},
   }
   vz {id=vso, x=vtas, y=sigma}
# end
# VOE process starting from eta=0, sigma=0 solutions from VSO
# Note: mode 5 takes a different path than in nl2.fp
# settings{wait}
# settings{d=2}
  flut {
  sequential
      id=pv, source=vso,
      indep=(vtas[0:3000], freq, massfac[0.5:2]),
      sigma=0,
      start{vtas[20:1000]},
      uncert=(stiffac{5%})
      plot=stepsize
      debug=uncert
   }
   vz {id=pv, x=massfac, y=vtas}

  flut {
      id=pvhilo, source=vso,
      indep=(vtas[0:3000], freq, massfac[0.5:2]),
      sigma=0,
      start{vtas[20:1000]},
      stiffac=(0.95,1.05)
      plot=stepsize
#     debug=uncert
   }
   vz {id=(pv,pvhilo), x=massfac, y=vtas}
end

stiffcn.cpp {{
int
stiffcn (pset& plt, int nr, int nc, CMatrix& ca) {
   Trace trc(2,"stiffcn");

   if (ca.nel() != nr*nc)
      throw runtime_error(vastr("stiffcn: (",nr,',',nc,") got ca.nel() ",ca.nel()));

   Real fac = plt.parval("stiffac");
   for (int i=1; i<=4; i++)
      ca(i,i) *= fac;

   return 0;
}
}}

massfcn.cpp {{
int
massfcn (pset& plt, int nr, int nc, CMatrix& ca) {
   Trace trc(2,"massfcn");

   if (ca.nel() != nr*nc)
      throw runtime_error(vastr("stiffcn: (",nr,',',nc,") got ca.nel() ",ca.nel()));

   Real sdamp = plt.parval("sdamp");
   Complex damp{1.0, sdamp};
   int k1b = 41;
   int k0b = k1b - 1;
   Real fac = plt.parval("massfac");
   for (int i=1; i<=4; i++)
      ca(i,i) *= fac;
   return 0;
}
}}
