#------------------------------------------------------------------
# generic 2-engine airplane
# freeplay on elevators, rudder, ailerons
# bilinear stiffness on nacelles
# bilinear factor on gaf
# Shows 3 looping lco curves which do not extend to
# zero eta (gcnorm)
#------------------------------------------------------------------

# freeplay (gap):
output{env="fpx0=0.001"}
output{env="blx0=0.01"}
output{env="k2k1=2.0"}

output{env="gafx0=0.001"}
output{env="gafratio=0.5"}

output{env="SplineRho=2"}

#  1-40   airplane flexible
#  41-47  rudder flexible
#  48     tab rotation
#  49     rudder/geared tab rotation
#  50-54  5 left elevator flexible
#  55     left elevator rotation
#  56-60  5 right elevator flexible
#  61     right elevator rotation
#  62-65  4 left aileron flexible
#  66     left aileron rotation
#  67-70  4 right aileron flexible
#  71     right aileron rotation
#  72-74  3 left nacelle flexible
#  75-77  3 right nacelle flexible


   import { ${FROOT}/demo/flut12.op4}
   import { ${FROOT}/demo/flut12.uf}

   param {
      i=gstif,
      code=stiffcn,
      o=stif,
      fpsmooth=0.1,
      blsmooth=0.1,
      fpx0=$fpx0,
      blx0=$blx0,
      k2k1=$k2k1,
   }

   param {
      i = gaf0{rf=0},
      i = gaf0008{rf=0.0008},
      i = gaf0015{rf=0.0015},
      i = gaf002{rf=0.002},
      i = gaf003{rf=0.003},
      i = gaf01{rf=0.01},
      i = gaf05{rf=0.05},
      o = pgaf,
      code=gaffcn,
      gafx0 = $gafx0,
      gafratio = $gafratio
   }

# VSO at gcnorm=0
   flut {
      id = vso,
      mass=mass,
      stif=stif,
      gaf = pgaf,
      indep=(vtas[0:600], freq[0:50], sigma),
      alt = 0,
      gcnorm = 0.0,
#     start{vtas=10, mode_10,mode_15,mode_16},
      start{vtas=10, freq[0:25]},
      target{sigma=0,vtas[10:1000]},
   }
   vis {id=vso, x=sigma, y=freq}
   vis {id=vso, x=vtas, y=sigma}

# soe at vtas= 100, 200
   flut {
      id=soe, source=vso,
      indep=(sigma,freq,gcnorm[0:1]),
      vtas=(100,200),
       target{sigma=0},
   }
#  vis {id=soe, y[0:0.3]=gcnorm, x=sigma, r=mode_10.a, r=mode_15.a}
#  vis {id=soe, y[0:0.3]=gcnorm, x=sigma, r=mode_16.b}

# lco analysis starting from soe
   flut {
      id=voe, source=soe,
      indep=(gcnorm[0:1], freq, vtas[0:400]),
      sigma=0,
      start{ordinal=1},  # first xing only
      checklooping,
   }
   vis {id=voe, x=vtas, y=gcnorm, w=lcostab}

# lco analysis starting from vso
   flut {
      id=voe0, source=vso,
      indep=(gcnorm[0:1], freq, vtas[0:400]),
      sigma=0,
      start{ordinal=1},  # first xing only
      checklooping,
   }
   vis {id=voe0, x=vtas, y=gcnorm, w=lcostab}
end

stiffcn.cpp {{

int
stiffcn (const pset& plt, int nr, int nc, CMatrix& ca) {
   Trace trc(2,"stiffcn");
   Timer tim("stiffcn");

   if (ca.nel() != nr*nc) {
      ostringstream os;
      os << "stiffcn: (" << nr << ',' << nc << ") got ca.nel() " << ca.nel();
      throw Err(os.str());
   }
//  48     tab rotation
//  49     rudder/geared tab rotation
//  55     left elevator rotation
//  61     right elevator rotation
//  66     left aileron rotation
//  71     right aileron rotation

   Real fpx0 = plt.parval("fpx0");  // rudder tab
   Builtin bi;
   vector<int> fpdof;
   fpdof.push_back(49);
   fpdof.push_back(55);
   fpdof.push_back(61);
   fpdof.push_back(66);
   fpdof.push_back(71);

   Real tmax{0.0};
   size_t j;
   for (j=0; j<fpdof.size(); j++) {
      int k1b = fpdof[j];
      int k0b = k1b - 1;
      Real fac = bi.freeplay(plt, fpx0, k1b);
      Complex t = ca(k1b,k1b);
      ca(k1b,k1b) *= fac;
      trc.dprint("stiffcn: scaled diag ",k1b," = ",ca(k1b,k1b),
         " = ",t,'*',fac);
      Real tk = plt.absgcv(k1b)/fpx0;
      if (tk.value() > tmax.value())
         tmax = tk;
   }
   plt.setpar(tmax, "fpmax");

// bilinear stiffness dof
//  72-74  3 left nacelle flexible
//  75-77  3 right nacelle flexible

   vector<int> bldof;
   Real blx0 = plt.parval("blx0");
   Real k2k1 = plt.parval("k2k1");

   for (int i=0; i<6; i++)
      bldof.push_back(72+i);

   Real blmax{0.0};
   for (j=0; j<bldof.size(); j++) {
      int k1b = bldof[j];
      int k0b = k1b - 1;
      Real fac = bi.bilinear(plt, blx0, k2k1, k1b);
      Complex t = ca(k1b,k1b);  // zero-based indexes
      ca(k1b,k1b) *= fac;
      trc.dprint("stiffcn: scaled diag ",k1b," = ",ca(k1b,k1b),
         " = ",t,'*',fac);
      Real tk = plt.absgcv(k1b)/blx0;
      if (tk.value() > blmax.value())
         blmax = tk;
   }
   plt.setpar(blmax, "blmax");
   return 0;
}
}}

gaffcn.cpp {{
int
gaffcn (const pset& plt, int nr, int nc, CMatrix& ca) {
   Trace trc(2,"gaffcn");
   Timer tim("gaffcn");

   if (ca.nel() != nr*nc) {
      ostringstream os;
      os << "gaffcn: (" << nr << ',' << nc << ") got ca.nel() " << ca.nel();
      throw Err(os.str());
   }

   Builtin bi;

//  48     tab rotation
//  49     rudder/geared tab rotation
//  55     left elevator rotation
//  61     right elevator rotation
//  66     left aileron rotation
//  71     right aileron rotation
   vector<int> dof;
   // nacelles
   // for (int k=72; k<78; k++)
   //    dof.push_back(k);
   // control-surfaces
   dof.push_back(48);
   dof.push_back(49);
   dof.push_back(55);
   dof.push_back(61);
   dof.push_back(66);
   dof.push_back(71);
   Real gafx0 = plt.parval("gafx0");
   Real gafratio = plt.parval("gafratio");
   Complex vfac = plt.parval("dpress");

   Real blmax{0.0};
   for (size_t j=0; j<dof.size(); j++) {
      int k1b = dof[j];
      int k0b = k1b - 1;
      Real fac = bi.bilinear(plt, gafx0, gafratio, k1b);
      // scale the entire column
      for (size_t i=1; i<=nr; i++) {
         ca(i,k1b) *= fac*vfac;
      }
      Real ti = plt.absgcv(k1b)/gafx0;
      if (ti.value() > blmax.value())
         blmax = ti;
   }
   plt.setpar(blmax, "blgafmax");
   return 0;
}
}}
