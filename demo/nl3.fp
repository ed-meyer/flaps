# small nonlinear problem based on HA145B
# Demonstrates the use of a custom describing function
# called from custom functions (stiffcn & massfcn)
#
# This is the third of four files demonstrating 4 different
# ways of implementing describing functions
#
# An example of latent LCO is shown in voe2

   import{${FROOT}/demo/ha145b.op4}

	# add a copy of QHHL7 at a high rf to extend the rf range
	octlab{i=(copy.m,QHHL7), o=QHHL8}

   pz { i = KHH, o=Stif, code=stiffcn }

   # the mass matrix has freeplay on dof 10
   pz { i = MHH, o=Mass, code=massfcn }

   pz {
      units=uscs
      out=GAF,
      i=(QHHL1{rf=0.00000001524}, QHHL2{rf=0.00001524}, QHHL3{rf=0.00076200},
         QHHL4{rf=0.001524}, QHHL5{rf=0.003048}, QHHL6{rf=0.00762},
         QHHL7{rf=0.01524}, QHHL8{rf=0.03}),
      beta = (0.0007, 0.007)
   }


# VSO - start points for SOE, VOE analyses
   flut {
      id=vso,
      indep=(vtas[0:500], freq[0:100], sigma),
      start{vtas=10, freq[0.1:7]},
      gcnorm = 0,
      alt=-2000.0,
      mass=Mass,
      stif=Stif,
      gaf=GAF,
      target{vtas[10:500], sigma=0},
   }
   vz {id=vso, x=vtas, y=sigma}

# VOE process starting from eta=0, sigma=0 solutions from VSO
  flut {
      id=voe0, source=vso,
      indep=(vtas[0:1000], freq, gcnorm[0:0.3]),
      sigma=0,
      start{vtas[20:500]},
   }
   vz {id=voe0, x=vtas, y=gcnorm, w=lcostab}

# SOE processes at a few velocities starting from vso
   flut {
      id=soe, source=vso,
      indep=(gcnorm[0:3], freq, sigma),
      vtas = (100,200,400),
      target{vtas[10:500], sigma=0},
   }
   vz {id=soe, x=sigma, y=gcnorm}

# VOE process starting from soe with sigma=0 
  flut {
      id=voe1, source=soe,
      indep=(vtas[0:1000], freq, gcnorm[0:0.3]),
      start{ordinal=1},
      sigma=0,
   }
   vz {id=voe1, x=vtas, y=gcnorm, w=lcostab}

# VSOE process starting from a few velocities
   flut {
      id=vsoe, source=vso,
      indep=(vtas[10:500], sigma[-100:0.05], freq, gcnorm[0:3]),
      constrained,
      start{vtas=(50,100,200)}
      optimize{sigma=1.0, gcnorm=0.5},
      target{vtas[10:500], sigma=0}
   }
   vz {id=vsoe, x=sigma, y=gcnorm}
   vz {id=vsoe, x=vtas, y=sigma}
   vz {id=vsoe, x=vtas, y=gcnorm}

# VOE process starting from vsoe with sigma=0 
# mode 1 shows a latent LCO
  flut {
      id=voe2, source=vsoe,
      indep=(vtas[0:1000], freq, gcnorm[0:3]),
      start{vtas[20:1000], ordinal=1},
      sigma=0,
      checklooping
   }
   vz {id=(vsoe,voe2), color, c=1, x=vtas, y=gcnorm}
end

stiffcn.cpp {{

Ad bilinear(pset& plt,const Ad& x0,const Ad& ratio,int gcno);

int
stiffcn (pset& plt, int nr, int nc, vector<complex<Ad>>& ca) {
// scale the first 4 diagonals by a bilinear describing function
   Ad fac = df::bilinear(plt, 0.01, 0.5, 1);
   ca[0] *= fac;
   // update a new parameter for the bilinear function on
   // dof 1: plot it versus absgc1
   plt.setpar(fac, "bl1");
   ca[1+nr] *= bilinear(plt, 0.01, 0.5, 2);
   ca[2*(1+nr)] *= bilinear(plt, 0.01, 0.5, 3);
   ca[3*(1+nr)] *= bilinear(plt, 0.01, 0.5, 4);
   return 0;
}

Ad
bilinear(pset& plt, Ad const& x0, Ad const& ratio, int gcno) {
// Describing function for a bilinear stiffness - a number between 1 and "ratio"
// which multiplies the nominal stiffness to give an equivalent stiffness.
// Input:
//    x0    (Ad) transition between k1 and k2
//    ratio ratio of the regions: k2/k1
//    gcno  (int) g.c. number (1b) which has the bilinear stiffness
// If x0 is zero:       return 1 regardless of abs(gcno)
// if abs(gcno) <= x0:  return 1
//    abs(gcno) > x0:   approaches "ratio" asymtotically

   Ad rval{1.0};
   if (x0.value() <= 0.0)
      return rval;
   Ad absgc = plt.absgcv(gcno);
   Ad gamma = absgc/x0;

   if (gamma.value() > 1.0) {
      Ad pi{4.0*atan(1.0)};
      Ad one(1.0);
      Ad two(2.0);
      Ad xi = asin(one/gamma);
      Ad gamsq = gamma*gamma;
      rval = ratio + (two/pi)*(one - ratio)*(xi + sqrt(gamsq-one)/gamsq);
   }
   return rval;
}
}}

massfcn.cpp {{

int
massfcn (pset& plt, int nr, int nc, vector<complex<Ad>>& ca) {
// scale diagonal 10 by a freeplay describing function
// Note: indices for the matrix are zero-based, but the index
// argument for freeplay is 1-based
   Ad fac = df::freeplay(plt, 0.001, 10);
   ca[9*(1+nr)] *= fac;
   // update a new parameter for the freeplay function on
   // dof 10: plot it versus absgc10
   plt.setpar(fac, "fac10");
   return 0;
}

Ad
freeplay(pset& plt, Ad const& gap, int gcno) {
// Describing function for freeplay - a number between 0 and 1
// which multiplies the nominal stiffness to give an equivalent stiffness.
// Input:
//    gap   (Ad) dead zone
//    gcno  (int) g.c. number (1b) which has the gap
// If gap is zero:       return 1 regardless of abs(gcno)
// if abs(gcno) <= gap:  return 0
//    abs(gcno) > gap:   approaches 1 asymtotically

   Ad rval{0.0};
   if (gap.value() <= 0.0)
      return rval;
   Ad absgc = plt.absgcv(gcno);

   if (absgc.value() > gap.value()) {
      Ad pi{4.0*atan(1.0)};
      Ad xi = 2.0*asin(gap/absgc);
      rval = 1.0 - (xi + sin(xi))/pi;
   }
   return rval;
}
}}

copy.m {{
QHHL8 = QHHL7;
}}
