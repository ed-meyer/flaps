
#------------------------------------------------------------------
# generic 2-engine airplane
# 1) Nacelle frequency variations: contours of right/left
#    side-bending frequency
# 2) Using octlab to reduce the size of the matrices
#------------------------------------------------------------------

   import { ${FROOT}/demo/lco.op4}
   import { ${FROOT}/demo/lco.uf}
	catalog{summary}

   # run matlab/octave to reduce the matrices
   octlab{
		i=(reduce.m, gct.lco, mass, gstif),
		i=(gaf0, gaf0008, gaf0015, gaf002, gaf003, gaf01, gaf05)
	}

	display{mass}

   # the nacelles are now dof 41-46, 3 dof per side
   # parameterize the stiffness matrix: nacelle sB frequencies
   parameters{rsbfreq(Right Nac SB Freq)[1:6]<radps2Hz> = 4,
      lsbfreq(Left Nac SB Freq)[1:6]<0.1591549 Hz/(r/s)> = 4}

   pz {
      o=stif,
      i=gstif
      code = stiffcn
   }

   # gaf matrices are interpolated wrt rf
   pz {
      i=(gaf0{rf=0}, gaf0008{rf=0.031496}, gaf0015{rf=0.059055}
         gaf002{rf=0.07874}, gaf003{rf=0.11811}, gaf01{rf=0.3937}, gaf05{rf=1.9685})
      o = pgaf
   }

   # vso
   flut {
      id = vso
      mass=mass
      stif=stif
      gaf = pgaf
      indep=(vtas[0:600], freq[0:50], sigma)
      sdamp = 0.01
      alt = 0
      start{vtas=10, mode=6}
      target{sigma=0,vtas[10:1000]}
      freevib
		vzid = lco
   }
   viz {id=vso, x[-1:1]=sigma, y=freq}
   viz {id=vso, x=vtas, y[-1:1]=sigma}

   # optimize vtas to use as contour start points
   flut {
      id=opt, source=vso
      start{vtas[40:1000]}
      optimize{vtas=1}
      sigma = 0
      indep=(lsbfreq[1:4.6],rsbfreq[1:4.6],freq,vtas[40:420])
   }
   viz {id=opt, x=lsbfreq, y=vtas}

   # lsbfreq-rsbfreq contours at several vtas
   flut {
      id=contour, source=opt
      indep=(lsbfreq,rsbfreq,freq)
      sigma=0
      vtas=(383:2:400)
      start{ordinal=1}
      checklooping
   }
   viz {id=contour, x=lsbfreq, y=rsbfreq}
   viz {id=(opt,contour), color, x=lsbfreq, y=rsbfreq}

end

reduce.m {{
whos
ext = [1:40 46:51];
gct__lco = gct__lco(:,ext);
gstif = gstif(ext,ext);
mass = mass(ext,ext);
gaf0 = gaf0(ext,ext);
gaf0008 = gaf0008(ext,ext);
gaf0015 = gaf0015(ext,ext);
gaf002 = gaf002(ext,ext);
gaf003 = gaf003(ext,ext);
gaf01 = gaf01(ext,ext);
gaf05 = gaf05(ext,ext);
whos
}}

stiffcn.cpp {{
int
stiffcn (pset& plt, int nr, int nc, vector<complex<Ad>>& ca) {
   Trace trc(2,"stiffcn");

   if (ca.size() != nr*nc)
      throw runtime_error(vastr("stiffcn: (",nr,',',nc,") got ca.size() ",ca.size()));

//      [41,41] = "mass[41,41]*rsbfreq^2",
//      [44,44] = "mass[44,44]*lsbfreq^2"
   Ad sdamp = plt.parval("sdamp");
   complex<Ad> damp{1.0, sdamp};
   int k1b = 41;
   int k0b = k1b - 1;
   // get the side-bending freqs in rad/s
   Ad rsbfreq = plt.parval("rsbfreq");
   Ad lsbfreq = plt.parval("lsbfreq");
   // complex<Ad> rmass = matvij(plt,"mass", k1b, k1b);
	double rmass{0.9};
   ca[IJ(k0b,k0b,nr)] = damp*rsbfreq*rsbfreq*rmass;
   k1b = 44;
   k0b = k1b - 1;
   // complex<Ad> lmass = matvij(plt,"mass", k1b, k1b);
	double lmass{0.9};
   ca[IJ(k0b,k0b,nr)] = damp*lsbfreq*lsbfreq*lmass;

   return 0;
}
}}
