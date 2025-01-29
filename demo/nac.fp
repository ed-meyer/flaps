
#-------------------------------------------------------------
# Demonstrates:
# 1) running matlab/octave from a Flaps script
# 2) nacelle side-bending frequency variation
# 3) right-click on v-g curve brings up amvz
#-------------------------------------------------------------

# import the full model and export matrices to a matlab/octave .mat file
   import { ${FROOT}/demo/nac.op4}
   import { ${FROOT}/demo/nac.uf}
   export{o=nac.mat, mass, stif, gaf0, gaf001, gaf005, gaf01, gctransform}
# run matlab/octave to reduce the matrices using reduce.m below; note it is on
# the flaps temporary directory FTMP
!  matlab -r "run('$FTMP/reduce.m')"
# import the result
   import { reduced.mat }

# Create interpolation coef
   pz {
      i = gf0{rf=0},
      i = gf001{rf=0.001},
      i = gf005{rf=0.005},
      i = gf01{rf=0.01},
      o = gaf,
   }

# parameterize the stiffness matrix: nacelle sB frequencies
   new Par{rsbfreq(Right Nac SB Freq)[0:30]<0.1591549 Hz/(r/s)> = 4}
   new Par{lsbfreq(Left Nac SB Freq)[0:30]<0.1591549 Hz/(r/s)> = 4}

   pz {
      o=stif,
      i=kred,
      [31,31] = "mred[31,31]*(2*PI*rsbfreq)^2",
      [34,34] = "mred[34,34]*(2*PI*lsbfreq)^2"
   }

# vso at gcnorm=0
   flut {
      id=vso,
      mass=mred, stif=stif, gaf=gaf,
      indep=(vtas[5:1000],sigma,freq),
      rsbfreq = 4, lsbfreq = 4
      alt=0,
      start{vtas=10, freq[0:10]},
      target{growth=0}
		freevib
   }
   vz {id=vso, x=vtas, y=sigma}
   vz {id=vso, x=sigma, y=freq}

# nacelle sb frequency variation
#  flut {
#     id=sbpv, source=vso,
#     indep=(vtas, freq, sbfreq[1:8]),
#     sigma=0
#  }
#  vz{id=sbpv, x=sbfreq, y=vtas}
   
# optimize vtas to use as contour start points
   flut {
      id=opt, source=vso
#     start{mode_5}
      optimize = vtas
      sigma = 0
      indep=(lsbfreq,rsbfreq,freq,vtas)
   }
   vz {id=opt, x=vtas, y=lsbfreq, y=rsbfreq}

   flut {
      id=contour, source=opt
      indep=(lsbfreq,rsbfreq,freq)
      sigma=0
      vtas=(740:40:860)
   }
   vz {id=(opt,contour), color, x=lsbfreq, y=rsbfreq}

   end

reduce.m {{
cd $FWD
load nac.mat
ext = [1:36];
gf0 = gaf0(ext,ext);
gf001 = gaf001(ext,ext);
gf005 = gaf005(ext,ext);
gf01 = gaf01(ext,ext);
gctransform = gctransform(:,ext);
mred = mass(ext,ext);
kred = stif(ext,ext);
save('reduced.mat','mred','kred','gctransform','gf0','gf001','gf005','gf01','-nocompression')
quit
}}
