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

# nonlinear dofs implemented in stiffcn and gaffcn:
#  41     rudder rotation
#  42     left elevator rotation
#  43     right elevator rotation
#  44     left aileron rotation
#  45     right aileron rotation
#  46-51  nacelles 3 dof/side

   import { ${FROOT}/demo/lco.op4}
   import { ${FROOT}/demo/lco.uf}

#  print {matview, gstif, mass}

   pz {
      i=gstif,
      code=stiffcn,
      o=stif,
      fpsmooth=0.1,
      blsmooth=0.1,
      fpx0=$fpx0,
      blx0=$blx0,
      k2k1=$k2k1,
   }

# gaf matrices are interpolated wrt rf, then passed
# to gaffcn for nonlinear mods
   pz {
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
      indep=(vtas[0:1000], freq[0:50], sigma),
      alt = 0,
      gcnorm = 0.0,
      start{vtas=10, freq[0:10]},
      target{sigma=0,vtas[10:1000]},
      freevib,
   }
   vis {id=vso, x[-1:1]=sigma, y=freq}
   vis {id=vso, x=vtas, y[-1:1]=sigma}

# lco analysis starting from vso
   flut {
      id=voe0, source=vso,
      indep=(gcnorm[0:1], freq, vtas[0:1000]),
      sigma=0,
      start{vtas[20:1000], ordinal=1},  # first xing only
      checklooping,
   }
   vis {id=voe0, x=vtas, y=gcnorm, w=lcostab}

# soe at various vtas
   flut {
      id=soe, source=vso,
      indep=(sigma,freq,gcnorm[0:1]),
      vtas=(50,80,100,200),
       target{sigma=0},
   }
#  vis {id=soe, y[0:0.3]=gcnorm, x=sigma, r=mode_10.a, r=mode_15.a}
#  vis {id=soe, y[0:0.3]=gcnorm, x=sigma, r=mode_16.b}

# lco analysis starting from soe
   flut {
      id=voe1, source=soe,
      indep=(gcnorm[0:1], freq, vtas[0:1000]),
      sigma=0,
      start{ordinal=1},  # first xing only
      checklooping,
   }
   vis {id=voe1, x=vtas, y=gcnorm, w=lcostab}
end

stiffcn.cpp {{

int
stiffcn (const pset& plt, int nr, int nc, CMatrix& ca) {
   Trace trc(2,"stiffcn");
   Timer tim("stiffcn");

   if (ca.nel() != nr*nc) {
      throw runtime_error(vastr("stiffcn: (",nr,',',nc,") got ca.nel() ",ca.nel()));
   }

// "DOF=(7 to 46,112,118,124,129,134,135 to 140)"
//  48     tab rotation
//  49     rudder/geared tab rotation -> 41
//  55     left elevator rotation -> 42
//  61     right elevator rotation -> 43
//  66     left aileron rotation -> 44
//  71     right aileron rotation -> 45
//  nacelle 3 dof/side 46-51

   Real fpx0 = plt.parval("fpx0");  // rudder tab
   Builtin bi;
   vector<int> fpdof;
   fpdof.push_back(41);
   fpdof.push_back(42);
   fpdof.push_back(43);
   fpdof.push_back(44);
   fpdof.push_back(45);

   Real tmax{0.0};
   size_t j;
   for (j=0; j<fpdof.size(); j++) {
      int k1b = fpdof[j];
      int k0b = k1b - 1;
      Real fac = bi.freeplay(plt, fpx0, k1b);
      Complex t = ca(k0b,k0b);
      ca(k0b,k0b) *= fac;
      trc.dprint("stiffcn: scaled diag ",k1b," = ",ca(k0b,k0b),
         " = ",t,'*',fac);
      Real tk = plt.absgcv(k1b)/fpx0;
      if (tk.value() > tmax.value())
         tmax = tk;
   }
   plt.setpar(tmax, "fpmax");

// bilinear stiffness dof
//  nacelle 3 dof/side 46-51

   vector<int> bldof;
   Real blx0 = plt.parval("blx0");
   Real k2k1 = plt.parval("k2k1");

   for (int i=0; i<6; i++)
      bldof.push_back(46+i);

   Real blmax{0.0};
   for (j=0; j<bldof.size(); j++) {
      int k1b = bldof[j];
      int k0b = k1b - 1;
      Real fac = bi.bilinear(plt, blx0, k2k1, k1b);
      Complex t = ca(k0b,k0b);  // zero-based indexes
      ca(k0b,k0b) *= fac;
      trc.dprint("stiffcn: scaled diag ",k1b," = ",ca(k0b,k0b),
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
// bilinear factor applied to some of the low-frequency modes
   Trace trc(2,"gaffcn");
   Timer tim("gaffcn");

   if (ca.nel() != nr*nc) {
      throw runtime_error(vastr("gaffcn: (",nr,',',nc,") got ca.nel() ",ca.nel()));
   }

   Builtin bi;

//  48     tab rotation
//  49     rudder/geared tab rotation
//  55     left elevator rotation
//  61     right elevator rotation
//  66     left aileron rotation
//  71     right aileron rotation
   vector<int> dof;
   // first few low-freq modes
   for (int k=1; k<5; k++)
      dof.push_back(k);
   //  nacelle 3 dof/side 46-51
   // for (int k=46; k<52; k++)
     //  dof.push_back(k);
   // control-surfaces
   // dof.push_back(41);
   // dof.push_back(42);
   // dof.push_back(43);
   // dof.push_back(44);
   // dof.push_back(45);
   Real gafx0 = plt.parval("gafx0");
   Real gafratio = plt.parval("gafratio");
   Complex vfac = plt.parval("dpress");

   Real blmax{0.0};
   for (size_t j=0; j<dof.size(); j++) {
      int k1b = dof[j];
      int k0b = k1b - 1;
      Real fac = bi.bilinear(plt, gafx0, gafratio, k1b);
      // scale the entire column
      for (size_t i=0; i<nr; i++) {
         ca(i,k0b) *= fac*vfac;
      }
      Real ti = plt.absgcv(k1b)/gafx0;
      if (ti.value() > blmax.value())
         blmax = ti;
   }
   plt.setpar(blmax, "blgafmax");
   return 0;
}
}}
