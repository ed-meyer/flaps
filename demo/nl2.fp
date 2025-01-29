
# small nonlinear problem based on HA145B
# Demonstrates the use of built-in describing functions
# in custom functions (stiffcn, massfcn)
#
# This is the second of 4 files demonstrating 4 different
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
   viz {id=vso, x=vtas, y=sigma}

# VOE process starting from eta=0, sigma=0 solutions from VSO
  flut {
      id=voe0, source=vso,
      indep=(vtas[0:1000], freq, gcnorm[0:0.3]),
      sigma=0,
      start{vtas[20:500]},
   }
   viz {id=voe0, x=vtas, y=gcnorm, w=lcostab}

# SOE processes at a few velocities starting from vso
   flut {
      id=soe, source=vso,
      indep=(gcnorm[0:3], freq, sigma),
      vtas = (100,200,400),
      target{vtas[10:500], sigma=0},
   }
   viz {id=soe, x=sigma, y=gcnorm}

# VOE process starting from soe with sigma=0 
  flut {
      id=voe1, source=soe,
      indep=(vtas[0:1000], freq, gcnorm[0:0.3]),
      start{ordinal=1},
      sigma=0,
   }
   viz {id=voe1, x=vtas, y=gcnorm, w=lcostab}

# VSOE process starting from a few velocities
   flut {
      id=vsoe, source=vso,
      indep=(vtas[10:500], sigma[-100:0.05], freq, gcnorm[0:3]),
      constrained,
      start{vtas=(50,100,200)}
      optimize{sigma=1.0, gcnorm=0.5},
      target{vtas[10:500], sigma=0}
   }
   viz {id=vsoe, x=sigma, y=gcnorm}
   viz {id=vsoe, x=vtas, y=sigma}
   viz {id=vsoe, x=vtas, y=gcnorm}

# VOE process starting from vsoe with sigma=0 
# mode 1 shows a latent LCO
  flut {
      id=voe2, source=vsoe,
      indep=(vtas[0:1000], freq, gcnorm[0:3]),
      start{vtas[20:1000], ordinal=1},
      sigma=0,
      checklooping
   }
   viz {id=(vsoe,voe2), color, c=1, x=vtas, y=gcnorm}
end

stiffcn.cpp {{

int
stiffcn (pset& plt, int nr, int nc, vector<complex<Ad>>& ca) {
// scale the first 4 diagonals by a bilinear describing function
// Note: indices for the matrix are zero-based, but the index
// argument for bilinear is 1-based
   Ad fac = df::bilinear(plt, 0.01, 0.5, 1);
   ca[0] *= fac;
   // update a new parameter for the bilinear function on
   // dof 1: plot it versus absgc1
   plt.setpar(fac, "bl1");
   ca[1+nr] *= df::bilinear(plt, 0.01, 0.5, 2);
   ca[2*(1+nr)] *= df::bilinear(plt, 0.01, 0.5, 3);
   ca[3*(1+nr)] *= df::bilinear(plt, 0.01, 0.5, 4);
   return 0;
}
}}

massfcn.cpp {{

int
massfcn (pset& plt, int nr, int nc, vector<complex<Ad>>& ca) {
// scale diagonal 10 by a freeplay describing function
   Ad fac = df::freeplay(plt, 0.001, 10);
   ca[9*(1+nr)] *= fac;
   // update a new parameter for the freeplay function on
   // dof 10: plot it versus absgc10
   plt.setpar(fac, "fac10");
   return 0;
}
}}

copy.m {{
QHHL8 = QHHL7;
}}
