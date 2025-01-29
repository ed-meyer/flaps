
# extracted (apex-6.0/demo/stab12.extract.ax):
#  row=(7 to 46,103 to 110,112 to 140),
# 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
#             1 2 3  4  5  6  7  8  9 10 11 12 13 14
# 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40
# 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34
# 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60
# 35 36 37 38 39 40
# 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80
#
# 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100
#
# 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115 116 117 118 119 120
#          41  42  43  44  45  46  47  48      49  50  51  52  53  54  55  56  57
# 121 122 123 124 125 126 127 128 129 130 131 132 133 134 135 136 137 138 139 140
#  58  59  60  61  62  63  64  65  66  67  68  69  70  71  72  73  74  75  76  77
# so now
#   129 -> 66
#   134 -> 71
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

# 777-300 ER modal-suppression control-law
# This model demonstrates root-tracking modes which
# bump into control-law poles: the 4.66 Hz mode and the 6.856 Hz mode
#
# Results:
#      id             desc
#      --             ----
#    openloop     no control-law, no flutter
#    pk           closed-loop tracking 15 modes in the range 0.1-7 Hz
#------------------------------------------------------------------
# Original job from A. Nadkarni:
#    D3F002PLO1M100C02080vg57RUDgvpv180yd0ms1_NFon_pk_apexbeta
#
#================================================================
#       SELECTED BRANCH      SELECTED MODES            AIRPLANE  
#          MODAL DOF       MODAL DESCRIPTION          MODAL DOF  
#----------------------------------------------------------------
#           1 TO   6      6 RIGID BODY
#----------------------------------------------------------------
#           7 TO  100    94 AIRPLANE FLEXIBLE
#----------------------------------------------------------------
#                 101       FLAPERON ROTATION  LEFT
#                 102       FLAPERON ROTATION RIGHT
#----------------------------------------------------------------
#         103 TO  109     7 RUDDER FLEXIBLE
#                 110       TAB ROTATION   
#                 111*      RUDDER/FAIRED TAB ROTATION
#                 112       RUDDER/GEARED TAB ROTATION
#----------------------------------------------------------------
#         113 TO  117     5 ELEVATOR FLEXIBLE  LEFT
#                 118       ELEVATOR ROTATION  LEFT
#         119 TO  123     5 ELEVATOR FLEXIBLE RIGHT
#                 124       ELEVATOR ROTATION RIGHT
#----------------------------------------------------------------
#         125 TO  128     4 AILERON FLEXIBLE  LEFT
#                 129       AILERON ROTATION  LEFT
#         130 TO  133     4 AILERON FLEXIBLE RIGHT
#                 134       AILERON ROTATION RIGHT
#----------------------------------------------------------------
#         135 TO  137     3 NACELLE FLEXIBLE  LEFT
#         138 TO  140     3 NACELLE FLEXIBLE RIGHT
#         
#                    *(USED IN DUBLAT ONLY NOT FLUTTER=> ESET=111)
#---------------------------------------------------------------- 
#---------------------------------------------------------------- 
# TOTAL AVAILABLE MODES = 140 >>>>>> TOTAL SELECTED MODES = 140   
#================================================================ 

#--------------------------------------------------
# DOF's in the monset CLAW in ELFINI Models
# 
#--------------------------------------------------
#
#================================================================
#  Count        DOF       ID              Location
#  
#----------------------------------------------------------------
#  1 to 6        6        CL1             ADIRU E3-3
#  7 to 12       6        CL2             SAARU E2-7
#  13 to 18      6        CL3             ACE L1 E1-5
#  19 to 24      6        CL4             ACE L2 E1-1
#  25 to 30      6        CL5             ACE C  E2-5
#  31 to 36      6        CL6             ACE R  E5-2
#  37 to 39      3        CL7             Gust Suppression delta p sensor (10% chord on fin at Z=484. (fin MAC)
#  40 to 45      6        CL8             Modal Suppression aft body lateral accelerometer
#  46 to 51      6        CL9             Pilot Seat
#  52 to 54      3        CL10            LH Wingtip front spar, wing closure rib
#  55 to 57      3        CL11            RH Wingtip front spar, wing closure rib
#  58 to 60      3        CL12            LH horiz stab tip front spar
#  61 to 63      3        CL13            RH horiz stab tip front spar
#  64 to 66      3        CL14            Fin tip front spar
#  67 to 72      6        CL15            Aft Pressure Bulkhead at floor, LH side
#  73 to 78      6        CL16            Wing front spar a center line ( pitch rate sensor)

####
#  79 to 84      6        CL17            GS Antenna            |-\
#  85 to 90      6        CL18            Radio Altimeter          - New Nodes inserted on May  2002.
#  91 to 96      6        CL19            Loc. Transmitter      |-/         (SEE BELOW)
#  97 to 102     6        CL20            LH Nacelle reference point         === ===== 
#  103 to 108    6        CL21            RH Nacelle reference point

######################################################################################################


   import { ${FROOT}/demo/flut12.op4}
   import { ${FROOT}/demo/flut12.uf}

#  print {mass, gstif, gaf0015}


#  extract { rudder = "stif,set=F002"{col=112} }

#  export {
#     o=stab12.op4
#     "stif,set=F002"
#     PLC020
#     "gaf,..."
#     rudder
#  }
#  amv {
#     o=stab12.uf
#     iset=1
#   conn="subset,id=ON70,iset=1"
#   modes="modes,iset=1"
#  }
#     "nodes,iset=1"
#     "freedoms,iset=1"
# end


# output{d=2}
   param {
      i = gaf0{rf=0},
      i = gaf0008{rf=0.0008},
      i = gaf0015{rf=0.0015},
      i = gaf002{rf=0.002},
      i = gaf003{rf=0.003},
      i = gaf01{rf=0.01},
#     i = gaf05{rf=0.05},
      o = pgaf,
#     beta = (0.001, 0.005),
      plot{diag, rf}
   }
# end

# output{debugger=valgrind}
# output{d=2}
# output{env="WAITFORDEBUG=1"}
# output{timer=2}
# output{perf}
   flut {
      id = vso0,
      mass=mass,
      stif=gstif,
      gaf = pgaf,
      indep=(vtas, freq, sigma),
      gcnorm = 0,
      alt = 0,
      start{freq[0:5]},
      target{sigma=0},
#     plot=(std,stepsize),
#     curvaturefactor = 6
   }
# end
# output{debugger=valgrind}
   vis {id=vso0, x=vtas, y=sigma}
end

   flut {
      id=soe1, source=vso0
      constrained
      indep=(gcnorm[0:3], sigma, freq)
      vtas=650
      curvaturefactor=7
   }
   vis {id=soe1, x=gcnorm, y=sigma, iset=1, conn="subset,id=ON70"}
end
# output{d=2}
   vis {id=vso0, x=vtas, y=sigma, iset=1, conn="subset,id=ON70"}

   flut {
      id=voe1, source=vso0
      indep=(vtas, freq, gcnorm[0:3])
      sigma=0
   }
   vis {id=voe1, x=vtas, y=gcnorm, iset=1, conn="subset,id=ON70"}
end

   flut {
      id = ssit
      mass=PLC020
      stif=Stif
      controls = Cont
      indep=(alt, freq, sigma)
      mach = 0.8
      eset = 111
      startregion{freq[0.1:7]}
      xkmqc(km) = 0
   }
   vis {id=ssit, x=veas, y=growth, iset=1, conn=ON70}

   flut {
      id = pk
      mass=PLC020
      stif=Stif
      gaf = pgaf
      controls = Cont
      indep=(alt[-5000:], freq, sigma)
      mach = 0.80
      eset = 111
      startregion{freq[0.1:7]}
      print = matrices
      plot=(boxes,alt,vcas,cdpress,veas,sigma,freq)
      target{growth=0}
   }

   vis {id=pk, x=veas, y=growth, iset=1, conn=ON70}
# Plot the pk solution along with the control-law poles
# and the boxes used to find start points. This shows how
# you can combine esa files and flut results in one plot.
   vis {
      id=pk, esa=(boxes,poles), color
      x=sigma, y=freq, xmin=-10, xmax=1, ymin=0, ymax=8
      iset=1, conn=ON70
   }
# the 4.66 Hz mode hits a control-law pole (see usrtf) as shown here:
   vis {
      id=pk, esa=(boxes,poles), color
      x=sigma, y=freq, xmin=-3, xmax=0, ymin=4.5, ymax=5.2
      iset=1, conn=ON70
   }


# Gain variation
   flut {
      source=pk, id=gain
      indep=(gain[0:3], freq, veas), sigma=0
      startregion{ veas[500:600] }
      target{gain=2.5}
   }
   vis { id=gain, x=gain, y=veas, iset=1, conn=ON70 }
end
