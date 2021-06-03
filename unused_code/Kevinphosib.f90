!
!=======================================================================
      subroutine Kevinphosib (test,test10,Tc,Psur, &
       Vmax0,Trop,SHTI,HHTI,SLTI,HLTI,TRDA, TRDM,Atheta,Btheta, &
       effcon,respcp,Tran,Ref,bee,PhiSat,Phi_half, &
       fPAR,Green,gmudmu, &
       C3,C4, pO2m, pCO2i, WWW, PAR, Tice)
!=======================================================================      
! calculates the zenith angle for an array of land points
!-----------------------------------------------------------------------
!
      implicit none
!
! input variables
      real Tc     ! (K) canopy temperature
      real Psur  ! (mb) surface total air pressure 
      real Vmax0  ! (mol/m2/s) Maximum Rubisco catalytic capacity at canopy top
      real Trop   ! (K) Reference Temperature for Q10 response exponent
      real SHTI   ! (1/K) Slope high temperature photosynthesis inhibition function
      real HHTI   ! (K) half point High temp photosynthesis inhibition function
      real SLTI   ! (1/K) Slope low temperature photosynthesis inhibition function
      real HLTI   ! (K) half point low temp photosynthesis inhibition function
      real TRDA   ! (1/K) half point Respiration high Temp inhibition function
      real TRDM   ! (K) half point Respiration high Temp inhibition function
      real Atheta ! (-) Coupling factor between Rubisco-light assimilation rates
      real Btheta ! (-) Coupling factor between Rubisco-light-export assim rates
      real effcon ! (mol/mol) quantum efficiency of C3, C4 photosynthesis
      real respcp    ! autotrophic leaf respiration Fraction of Vmax
      real Tran(2,2) ! Leaf transmittance for green/brown plants
      real Ref(2,2)  ! Leaf reflectance for green/brown plants
!                      For LTran and LRef:
!                      (1,1)=shortwave, green plants
!                      (2,1)=longwave, green plants
!                      (1,2)=shortwave, brown plants
!                      (2,2)=longwave, brown plants
      real fPAR   ! (-) Canopy absorbed fraction of PAR
      real Green  ! (-) Canopy greeness fraction of LAI
      real gmudmu ! (-) Time-mean leaf projection
      real C3     ! (-) C3 plant area cover fraction
      real C4     ! (-) C3 plant area cover fraction
      real pO2m   ! (Pa) partial pressure of O2 in mixed layer  
      real pCO2i  ! (Pa) partial pressure of CO2 in leaf interior
      real PAR    ! (W/m2/s) vis direct and diffuse flux incident on canopy
      real Tice   ! (K) freezing temperature of ice
      real WWW    ! (-) fraction of water saturation in soil
      real bee    ! (-) soil wetness exponent
      real PhiSat ! (m) Soil tension at saturation
      real Phi_half ! (m) 1/2 Critical leaf water potential limit
!
! output variables
      real test, test10(10)   ! test variables
!
! internal variables
      real Qt     ! (-) temperature scaling exponent
      real Kc     ! (Pa) Michaelis-Menten CO2 constant (SE(1996) eqn C1)
      real Ko     ! (Pa) inhibition constant for O2 (SE(1996) eqn C1)
      real Spfy   ! (-) Rubisco specificity CO2 vs. O2 (Pa)(SE(1996) eqn C1)
      real gammas ! (Pa) CO2 compensation point (SE(1996) eqn C1)
      real TempH  ! (-) high temp inhibition factor
      real TempL  ! (-) low temp inhibition factor
      real ScaleTem   ! (-) temperature scaling factor for photosynthesis
      real ScaleH2O   ! (-) soil water scaling factor for photosynthesis
      real ScaleRub   ! (-) Rubisco/CO2 absorption factor for rubisco limited rate
      real ScalePAR   ! (-) PAR/CO2 absorption factor for PAR limited rate
      real ScaleExp   ! (-) Export/CO2 factor for export limited rate
      real ScaleRsp   ! (-) autotrophic leaf respiration scaling factor
      real Vm     ! (mol/m2/s) temperature scaled rubisco catalytic capacity
      real omegaC ! (mol/m2/s) Rubisco limted assimilation rate
      real omegaE ! (mol/m2/s) Light limited assimilation rate
      real omegaS ! (mol/m2/s) export (C3) or PEP (C4) limited assimilation rate
      real assim  ! (mol/m2/s) smoothed minimum rubisco, light, & export rates
      real RespC  ! (mol/m2/s) canopy autotrophic leaf respiration rate
      real aparkk ! (-) PAR use parameter (SE(96) eqn C15)
      real APARKKmax     ! (-) Theoretical maximum PAR use parameter
      real scatp  ! (-) leaf scattering coef for PAR (SE(96) eqn C14)
      real scatg  ! (-) green leaf absorption coef for PAR (SE(96) eqn C4)
      real TPrcor ! (K m3/mol) CO2 conductance conversion factor mol/m2/s to m/s
!
!-----------------------------------------------------------------------
! Calculate scaling factors
!-----------------------------------------------------------------------
! Temperature scaling exponent
      Qt=0.1*(TC-Trop)
!
! Michaelis-Menten constant for CO2 (SE(1992a) Tbl 2; SE(1996), eqn C1/C17)
      Kc=30.*2.1**Qt
!
! inhibition constant for O2 (SE(1992a) Tbl 2; SE(1996) eqn C1/C17)
      Ko=30000.*1.2**Qt
!
! Rubisco specificity CO2 vs. O2 (S)(SE(1992a) Tbl 2; SE(1996) eqn C1/C17)
      Spfy=2600.*0.57**Qt
!
! CO2 compensation point (assumes po2i=pO2m) (SE(1992a) Tbl 2; SE(1996) eqn C1)
      gammas=.5*pO2m/Spfy
!
! Temperature scaling factor
      TempH=1+exp(shti*(Tc-HhTI))
      TempL=1+exp(slti*(HLTI-Tc))
      ScaleTem=2.**Qt*(C3/TempH+C4/TempL/TempH)
!
! Rubisco/CO2 absorption scaling factor (assumes po2i=pO2m)
      ScaleRub=C3*((pCO2i-gammas)/(pCO2i+Kc*(1+pO2m/Ko)))+C4
!
! PAR/CO2 absorption scaling factor
      ScalePAR=C3*(pCO2i-gammas)/(pCO2i+2*gammas)+C4
!
! Export/CO2 absorption scaling factor
      ScaleExp=C3/2.+C4*2.e4*pCO2i/Psur/100.
!
! Autotrophic leaf respiration scaling factor
      ScaleRsp=2.**Qt/(1.+EXP(TRDA*(Tc-TRDM)))*(C3*respcp+C4*0.025)
!
! soil water scaling factor
      scaleH2O=1./(1.+EXP(0.02*(Phi_half-PhiSat*WWW**(-BEE))))
!
! green leaf absorption coefficient for PAR
      SCATg=1.-TRAN(1,1)-REF(1,1)
!
! canopy extinction coefficient for PAR
      SCATP=green*(1.-TRAN(1,1))+(1.-green)*(1.-TRAN(1,2))
!
! PAR use parameter (SE (1992a) eqn 31; SE(96) eqn C15)
      APARKK=green*fPAR/(SCATP*GMUDMU)
      APARKKmax=0.95/(1.-TRAN(1,1))/GMUDMU
!
! Rubisco catalytic capacity of canopy
      Vm=Vmax0*ScaleTem*ScaleH2O*APARKK
!
!-----------------------------------------------------------------------
! Calculate assimilation and respiration rates
!-----------------------------------------------------------------------
! Rubisco limted assimilation rate (SE(1996) eqn C1)
      omegaC=Vmax0*ScaleTem*ScaleH2O*APARKK*ScaleRub
!
! PAR limted assimilation rate (SE(1996) eqn C3)
      omegaE=PAR*4.6E-6*SCATg*GMUDMU*APARKK*EFFCON*ScalePAR
!
! export/PEP limited assimilation rate (SE(1992a) eqn 13; SE(1996) eqn C5)
      omegaS=Vmax0*ScaleTem*ScaleH2O*APARKK*ScaleExp
!
! smoothed minimum between light, Rubisco, and export rates
      call smoothmin (Atheta,Btheta,omegaC,omegaE,omegaS,assim)
!
! canopy autotrophic or leaf respiration rate (SE(1996) eqn C8)
      RespC=Vmax0*ScaleH2O*APARKK*ScaleRsp
      test=RespC
      test10(1)=RespC
!
!-----------------------------------------------------------------------
! Calculate conductances
!-----------------------------------------------------------------------
! Calculate CO2 conductance conversion factor from mol/m2/s to m/s
      TPRCOR=44.6*(Tice+25.0)*Psur*100./1.013e5
!
      return
      end
