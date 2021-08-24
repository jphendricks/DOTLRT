!
!==================SUBROUTINE PHOSIB===================================
      SUBROUTINE PHOSIB(pco2m,pco2ap,po2m,vmax0,tf,psur,green &
      ,          tran,ref,gmudmu,zlt,cas_cap_co2,tc,ta,trop,trda &
      ,          trdm,slti,shti,hltii,hhti,radn,etc &
      ,          ea,rb,ra,tm &
      ,          effcon,rstfac,binter,gradm,assimn &
      ,          rst,atheta,btheta,tgs,respcp &
      ,          aparkk,len,nsib &
      ,     omepot,assimpot,assimci,antemp,assimnp &
      ,     wsfws,wsfht,wsflt,wci,whs &
      ,     wags,wegs,aparc,pfd,assim,td &
      ,     tg &
      ,     drst,pdamp,qdamp,ecmass,dtt,bintc,tprcor,ansqr &
      ,     nsoil, respg, pco2c, pco2i &
      ,     pco2s,co2cap)
!=======================================================================
!     CALCULATION OF CANOPY CONDUCTANCE USING THE INTEGRATED   
!     MODEL RELATING ASSIMILATION AND STOMATAL CONDUCTANCE.
!     UNITS ARE CONVERTED FROM MKS TO BIOLOGICAL UNITS IN THIS ROUTINE.
!     BASE REFERENCE IS SE-92A
!
! UNITS
!      PCO2M, PCO2A, PCO2Ap, PCO2I, PO2M        : PASCALS
!      CO2A, CO2S, CO2I, H2OA, H2OS, H2OA       : MOL MOL-1
!      VMAX0, RESPN, ASSIM, GS, GB, GA, PFD     : MOL M-2 S-1
!      EFFCON                                   : MOL CO2 MOL QUANTA-1
!      GCAN, 1/RB, 1/RA, 1/RST                  : M S-1
!      EVAPKG                                   : KG M-2 S-1
!      Q                                        : KG KG-1
!
! CONVERSIONS
!      1 MOL H2O           = 0.018 KG
!      1 MOL CO2           = 0.044 KG
!      H2O (MOL MOL-1)     = EA / PSUR ( MB MB-1 )
!      H2O (MOL MOL-1)     = Q*MM/(Q*MM + 1)
!      GS  (CO2)           = GS (H2O) * 1./1.6 (applies to Ci-Cs pathway)
!pl    44.6 moles of air per cubic meter
!      GS  (MOL M-2 S-1 )  = GS (M S-1) * 44.6*TF/T*P/PO
!      PAR (MOL M-2 S-1 )  = PAR(W M-2) * 4.6*1.E-6
!      MM  (MOLAIR/MOLH2O) = 1.611
!
! OUTPUT
!       ASSIMN         CARBON ASSIMILATION FLUX (MOL M-2 S-1) 
!       RST            CANOPY RESISTANCE (S M-1)
!       1/RST          CANOPY CONDUCTANCE (m s-1)
!       RSTFAC(4)      CANOPY RESISTANCE STRESS FACTORS
!         RSTFAC(1) ( F(H-S) )               : EQUATION (17,18), SE-92A
!         RSTFAC(2) ( F(SOIL) )              : EQUATION (12 mod), SE-89
!         RSTFAC(3) ( F(TEMP) )              : EQUATION (5b)   , CO-92
!         RSTFAC(4) ( F(H-S)*F(SOIL)*F(TEMP)) 
!       RESPC          CANOPY RESPIRATION (MOL M-2 S-1)
!       RESPG          GROUND RESPIRATION (MOL M-2 S-1)
!       PCO2I          CANOPY INTERNAL CO2 CONCENTRATION (MOL MOL-1)
!       GSH2O          CANOPY CONDUCTANCE (MOL M-2 S-1)
!       H2OS           CANOPY SURFACE H2O CONCENTRATION (MOL MOL-1)
!       PCO2I          INTERNAL CO2 CONCENTRATION
! Modifications:
!  C.Zhang & Joe Berry added new diagnostics (10/19/92)
!  gs stomatal conduct. reduced for freezing soil by Jim Collatz (02/21/95)
!  DD Modified for multitasking; add gather/scatter indices (12/06/95)
!  Ian Baker Added pco2c (chloroplast partial co2) for isotope stuff (Sep 99)
!  Kevin Schaefer added some variable definitions to comments (11/16/00)
!  Kevin Schaefer added back green fraction to PAR use parameter (2/8/01)
!  Kevin Schaefer deleted 8 unused variables, added more definitions (5/4/01)
!  Kevin Schaefer added maximum PAR use parameter (5/14/01)
!  Kevin Schaefer deleted 25 unused variables (6/28/01)
!=======================================================================
!
      implicit none
!
! Input variables
      integer  &
        len    ! number of SiB points in chunk &
      , nsib   ! total number of SiB points &
      , nsoil  ! total number of soil layers
!
      REALX8 &
        vmax0(len)     ! Max Rubisco catalytic capacity at canopy top (mol/m2/s) &
      , respcp(len)    ! respiration fraction of VMAX &
      , green(len)     ! greeness fraction of canopy &
      , aparc(len)     ! absorbed fraction of incident PAR &
      , gmudmu(len)    ! time mean leaf projection divided by cosine sunangle &
      , zlt(len)       ! leaf area index &
      , trop(len)      ! Reference Temp for Q10 response exponent (K) &
      , trda(len)      ! slope temp canopy resp inhibition function (1/K) &
      , trdm(len)      ! 1/2 pt temp canopy resp inhibition function (K) &
      , slti(len)      ! Slope of low temp inhibition function (K) &
      , shti(len)      ! Slope high temp inhibition function (1/K) &
      , hltii(len)     ! 1/2 pt low temp photosynthesis inhibition function (K) &
      , hhti(len)      ! 1/2 pt High temp photosynthesis inhibition function (K) &
      , effcon(len)    ! quantum efficiency of C3, C4 photosynthesis (mol/mol) &
      , binter(len)    ! Ball-Berry Conductance Intercept (mol/m2/s) &
      , gradm(len)     ! Ball-Berry Conductance Photosynthesis Slope (none) &
      , atheta(len)    ! Coupling factor for Rubisco-light assim rates &
      , btheta(len)    ! Coupling factor for Rubisco-light-export assim rates &
      , tran(nsib,2,2) ! leaf transmittance &
      , ref(nsib,2,2)  ! leaf reflectance &
      , ra(len)        ! resistance from canopy air space to PBL &
      , rb(len)        ! resistance from canopy to canopy air space &
      , tc(len)        ! canopy temperature (K) &
      , ta(len)        ! canopy air space temperature (K) &
      , tm(len)        !  &
      , tg(len)        ! &
      , tgs(len)       ! ground surface temperature (K) &
      , tf             !  &
      , td(nsib,nsoil) !  &
      , etc(len)       ! canopy airspace saturation vapor pressure  &
      , rst(len)       ! &
      , pdamp          ! &
      , qdamp          ! &
      , dtt            ! time step (s) &
      , po2m           ! partial pressure of O2 in leaf interior (Pa) &
      , ea(len)        ! Canopy airspace water vapor pressure &
      , respg(len)     ! ground respiration (mole m-2 m-1) &
      , psur(len)      ! surface pressure (mb) &
      , radn(len,2,2)  ! vis/nir direct/diffuse flux at canopy top &
      , ecmass(len)    ! (mm) canopy evapotranspiration
!
! Output variables
      REALX8  &
        assimn(len)   ! (mol/m2/s) net canopy assimilation &
      , assim(len)    ! (mol/m2/s) gross canopy assimilation &
      , antemp(len)   ! (mol/m2/s) net assimilation Bottom stopped at zero &
      , ansqr(len)    ! antemp**2 &
      , omepot(len)   ! (mol/m2/s) potential PAR limited rate (SE1996) eqn C3) &
      , assimpot(len) ! (mol/m2/s) Potential assimilation &
      , assimci(len)  ! (mol/m2/s) stress limited assimilation &
      , assimnp(len)  ! (mol/m2/s) Potential Net assimilation &
      , rstfac(len,4) ! (-) water, humidity, temp, & total stress factors &
      , drst(len)     ! change in stomatal resistance &
      , bintc(len)    ! Ball-Berry canopy stomatal conductance, term 2 &
      , xgah2o(len)   ! bottom stopped water conductance from CAS to PBL &
      , whs(len)      ! Intermediate assim weighted rel humidty stress factor &
      , wsfws(len)    ! potential assimilation weighted Water stress &
      , wsfht(len)    ! potential assimilation weighted High temperature stress &
      , wsflt(len)    ! potential assimilation weighted Low temperature stress &
      , wci(len)      ! Intermediate assimilation weighted Ci &
      , wags(len)     ! Intermediate assimilation weighted stomatal conductance &
      , wegs(len)     ! Intermediate evaporation weighted stomatal conductance &
      , pfd(len)      ! photon flux density direct/diffuse PAR normal to leaf &
      , tprcor(len)   ! 
!     
! Internal variables:
      integer i &
      , ic1 &
      , ic &
      , l &
      , icconv(len) &
      , igath(len)
!
      REALX8  &
        assimny(len,6) ! (mol/m2/s) iterated net canopy assimilation &
      , assimy(len,6)  ! (mol/m2/s) iterated gross canopy assimilation &
      , c3(len)      ! c3 plant fraction per grid cell &
      , c4(len)      ! c4 plant fraction per grid cell &
      , range(len)   ! range of assimilation (compensation pt to max rate) &
      , gammas(len)  ! CO2 compensation point (SE (1996) eqn C1) &
      , aparkk(len)  ! PAR use parameter (SE (1996) eqn C15) &
      , gsh2o(len)   ! stomatal conductance for water &
      , gah2o(len)   ! water conductance from canopy air space to PBL &
      , gbh2o(len)   ! water conductance from canopy to canopy air space &
      , par(len)     ! PAR available for photosynthesis (SE (1996) eqn C13) &
      , rrkk(len)    ! combined Michaelis-Menten constants (SE (1996) eqn C1) &
      , omss(len)    ! Temp stressed Vm for omspot (SE (1996) eqn C5) &
      , vm(len)      ! stressed Rusbisco catalytic capacity (SE (1996) eqn C17) &
      , templ(len)   ! denominator Vm temp stress factor (SE (1996) eqn C17) &
      , temph(len)   ! numerator Vm temp stress factor (SE (1996) eqn C17) &
      , qt(len)      ! Q10 coeff for temp stress factor (SE (1996) eqn C17) &
      , scatp(len)   ! leaf scattering coef for PAR (SE (1996) eqn C14) &
      , scatg(len)   ! green leaf scattering coef for PAR &
      , park(len)    ! time-mean, rad weighted PAR ext coef (SE (1996) eqn C14) &
      , respn(len)   ! leaf respiration rate (SE (1996) eqn C8) &
      , respc(len)   ! canopy respiration rate &
      , zkc(len)     ! temp adj Michaelis-Menten CO2 const (SE (1996) eqn C1) &
      , zko(len)     ! inhibition constant for O2 (SE (1996) eqn C1) &
      , spfy(len)    ! temp stress Rubisco specificity CO2/O2 (SE (1996) eqn C1) &
      , cas_cap_co2(len) ! canopy air space cap on CO2 &
      , co2a(len)    ! (mol/mol) Canopy Air Space CO2 concentration &
      , pco2a        ! (Pa) Canopy Air Space CO2 partial pressure &
      , pco2ap(len)  ! (Pa) previous canopy air space CO2 partial pressure &
      , pco2c(len)   ! (Pa) chloroplast co2 partial pressure &
      , pco2i(len)   ! (Pa) stomate co2 partial pressure &
      , pco2in       ! (Pa) stomate co2 partial pressure from conductance &
      , PCO2Y(len,6) ! (Pa) Stomate CO2 partial pressure from rate eqn &
      , EYY(len,6)   ! (Pa) rate-conductance delta Stomate CO2 partial pressure &
      , pco2ipot     ! (Pa) Potential stomate CO2 partial pressure &
      , pco2s(len)   ! (Pa) Leaf surface CO2 partial pressure &
      , co2s(len)    ! (mol/mol) Leaf surface CO2 concentration &
      , co2m(len)    ! (mol/mol) Boundary Layer CO2 concentration &
      , pco2m(len)   ! (Pa) Boundary layer CO2 partial pressure &
      , co2cap(len)  ! (mol/mol) maximum CO2 concentration in air space &
      , soilfrz(len) ! soil freezing point &
      , cwsflt       ! low temperature stress coefficient &
      , cwsfht       ! high temperature stress coefficient &
      , cwsfws       ! water stress coefficient &
      , ccoms        ! Sink or peptide control coefficient &
      , ccomc        ! Rubisco control coefficient &
      , ascitemp     ! Bottom stopped stress limited top leaf photosynthesis &
      , dompdomc     ! intermediate values of control coefficient &
      , omsci        ! canopy carbo export/CO2 limit(SE (1996) eqn C21) &
      , ompci        ! canopy smoothed minimum omcpot & omepot &
      , omcci        ! canopy Rubisco limited rate (SE (1996) eqn C19) &
      , omcpot       ! Rubisco limited rate of photosynthesis (SE (1996) eqn C1) &
      , omppot       ! smoothed minimum of omcpot & omepot (SE (1996) eqn C6) &
      , sqrtin       ! root portion of solution to quadratic formula (a**2-4ac) &
      , omspot       ! carbo export (c3)/CO2 limit (c4) rate(SE (1996) eqn C5) &
      , gsh2oinf     ! Ball-Berry stomatal conductance &
      , h2osrh       ! relative humidity stress factor &
      , h2os         ! (mol/mol) leaf surface water mixing ratio &
      , ecmole       ! (mol/s) evapotranspiration rate in moles &
      , h2oa         ! (mol/mol) canopy air space water mixing ratio &
      , h2oi         ! (mol/mol) stomate water mixing ratio (assumed saturated) &
      , dtti         ! inverse of time step size &
      , soilfrztd    ! soil dew point &
      , soilfrztg    ! soil freezing point &
      , xgco2m(len)  ! scaling factor for chloroplast pCO2 &
      , APARKKmax    ! Theoretical maximum PAR use parameter
!
!pl introduce a co2 capacity 
!pl this will basically be the mass of air under the top of the canopy (in
!pl this case (CHEAS-RAMS) O(10-30m), that is, ground to displacemnt height.
!pl all the carbon fluxes are expresse as Mol C / m2 s and resistances for
!pl carbon are in m2 s / mol air
!pl one mole of gas occupies 22.4 cubic dm
!pl 1 cubic meter contains therefore 1000./22.4  = 44.6 moles of gas
!pl the units for the carbon capacity are mol air /m2. 
!pl (e.g. here 893 moles if thickness of the layer is 20m)
!pl this means that the units for pc02ap should be mol co2 / mol air, but
!pl it is also possible to keep just co2 pressure and convert
!
! First Loop over gridpoints
      do i = 1,len
!
        TPRCOR(i) = (TF+25.0)*PSUR(i)*100./1.013 E5
        co2cap(i) = cas_cap_co2(i) * 44.6 * tprcor(i)/ta(i)     ! moles air/m2
!
!pl this needs to be modified as in sibslv3 to automatically use the 
!pl thickness of the canopy air space. 
!
!       set c3/c4 plant fraction
        IF(EFFCON(i).GT.0.07) then
          C3(i)=1.
        else
          C3(i)=0.
        endif
         C4(i)=1.-C3(i)
!
!-----------------------------------------------------------------------
! Canopy and Leaf optical properties
!-----------------------------------------------------------------------
!       Average leaf extinction coefficient for PAR
        SCATP(I)=GREEN(i)*(TRAN(i,1,1)+REF(i,1,1)) &
           +(1.-GREEN(i))*(TRAN(i,1,2)+REF(i,1,2))
!
!       Green leaf extinction coefficient for PAR
        SCATG(i)=TRAN(i,1,1)+REF(i,1,1)
!
!       Time mean canopy extinction coefficient for PAR (SE (1996) eqn C14)
        PARK(i)=SQRT(1.-SCATP(i))*GMUDMU(i)
!
!       PAR use parameter (SE (1992a) eqn 31; SE (1996) eqn C15)
        APARKK(i)=APARC(i)/PARK(i)*Green(i)
        APARKKmax=0.95/SQRT(1.-TRAN(i,1,1)-REF(i,1,1))/GMUDMU(i)
!
!       Photon Flux Density: direct/diffuse PAR normal to projected leaf area
        PFD(i)=4.6E-6*GMUDMU(i)*(RADN(i,1,1)+RADN(i,1,2))
!
!       PAR available for photosynthesis  (SE (1996) eqn C13)
        PAR(i)=pfd(i)*EFFCON(i)*(1.-SCATG(i))
!
!-----------------------------------------------------------------------
! Preliminary photosynthesis calculations
!-----------------------------------------------------------------------
!       Q10 coefficient for temp stress (SE (1992a) Table 2; SE (1996) eqn C17)
        qt(i)=0.1*(TC(i)-TROP(i))
!
!       Leaf respiration (SE (1996) eqn C8)
        RESPN(i)=RESPCP(i)*VMAX0(i)*RSTFAC(i,2)
!
!itb    patch to prevent underflow if temp is too cool...
        if(TC(i) >= TRDM(i))then
          RESPC(i)=RESPN(i)*2.0**qt(i) &
              /(1.+EXP(TRDA(i)*(TC(i)-TRDM(i))))
        else
          RESPC(i) = RESPN(i) * 2.0**qt(i)
        endif
!
!       Adjust max catalytic capacity of Rusbisco (SE (1996) eqn C17)
        VM(i) = VMAX0(i) * 2.1**qt(i)
        TEMPL(i) = 1. + EXP(SLTI(i)*(HLTIi(i)-TC(i)))
        TEMPH(i) = 1. + EXP(SHTI(i)*(TC(i)-HHTI(i)))
        RSTFAC(i,3) = 1./( TEMPL(i)*TEMPH(i))
        VM(i)    = VM(i)/TEMPH(i) * RSTFAC(i,2)*C3(i) &
            + VM(i) * RSTFAC(i,2)*RSTFAC(i,3) * C4(i)
        RSTFAC(i,3)=2.1**qt(i)*(C3(i)/TempH(i)+C4(i)*TempL(i)/TempH(i))
!
!       Michaelis-Menten constant for CO2 (Kc) 
!       (SE (1992a) Table 2; SE (1996) eqn C1)
        ZKC(i) = 30. * 2.1**qt(i)
!
!       inhibition constant for O2 (Ko) (SE (1992a) Table 2; SE (1996) eqn C1/C17)
        ZKO(i) = 30000. * 1.2**qt(i)
!
!       temp stress Rubisco specificity CO2/O2 (S) 
!       (SE (92a) Tbl 2; SE (1996) eqn C1)
        SPFY(i) = 2600. * 0.57**qt(i)
!
!       CO2 compensation point (GAMMA-*) (SE (1992a) Table 2; SE (1996) eqn C1)
        GAMMAS(i) = 0.5 * PO2M/SPFY(i) * C3(i)
!
!pl     Convert from m/s to being mol/ (m2 sec)
        GSH2O(i)  = 1.0/RST(i) * 44.6*TPRCOR(i)/TC(i)
        GBH2O(i)  = 0.5/RB(i) * 44.6*TPRCOR(i)/TC(i)
        GAH2O(i)  = 1.0/RA(i) * 44.6*TPRCOR(i)/TM(i)
!
        xgah2o(i) = max(0.466 EEEXP 0, gah2o(i) )
        xgco2m(i) = 4000.0 * vmax0(i)*APARKK(i)
!
!       combined Michaelis-Menten constants (SE (1996) eqn C1)
        RRKK(i)=ZKC(i)*(1.+PO2M/ZKO(i))*C3(i) &
                     +VMAX0(i)/5.*(1.8**qt(i))*C4(i)
!
!       Soil freezing and dew points
        soilfrztg=1.+exp(-1.5*(max(270. EEEXP 0,tgs(i))-273.16))
        soilfrztd=1.+exp(-1.5*(max(270. EEEXP 0,td(i,nsoil))-273.16))
        soilfrz(i)=max(1./soilfrztg, 1./soilfrztd)
        soilfrz(i)=max(soilfrz(i),0.05 EEEXP 0)
!
!       Ball-Berry stomatal conductance, term 2 (SE (1992a) eqn 35)
        BINTC(i)=BINTER(i)*ZLT(i)*GREEN(i)*RSTFAC(i,2)*soilfrz(i)
!
!       export limited photosynthesi rate (SE (1992a) eqn 13, 37c)
        OMSS(i)=(VMAX0(i)/2.0)*(1.8**qt(i)) &
                        /TEMPL(i)*RSTFAC(i,2)*C3(i) &
                        + RRKK(i)*RSTFAC(i,2)*C4(i)
!
!       Possible range of assimilation: compensation pt to max assim rate
        RANGE(i)=PCO2M(i)*(1.-1.6/GRADM(i))-GAMMAS(i)
	icconv(i)=1
      enddo ! end first grid point loop
!
! Clear out some arrays
      DO IC = 1, 6
        do i = 1,len
          PCO2Y(i,IC) = 0.
          EYY(i,IC) = 0.
	enddo
      enddo
!
!pl beginning of PL's setup
      DO IC = 1, 6
!
!       estimate next PCO2
        CALL SORTIN( EYY, PCO2Y, RANGE, GAMMAS, ic,len )
!
!       Calculate assimilation rate
        CALL CYCALC(APARKK, VM, ATHETA, BTHETA,par, &
                    GAMMAS, RESPC, RRKK, OMSS, C3, C4, &
                    PCO2Y(1,ic), assimny(1,ic), assimy(1,ic),len)
!
        do i = 1,len
!pl       diagnose current CO2 flux in mol/(m2 s)
!pl       this is a modified ra that will get the right units
!pl       in the conservation equation. Its units are m2 s/(mol_air)
	  gah2o(i)=1./MAX(0.446 EEEXP 0,GAH2O(i))
!
!pl       prognose new CAS CO2 according to flux divergence
!pl       This is in mol C/mol air (same as PaC/PaAir)
          CO2A(i)=PCO2Ap(i)/(PSUR(i)*100.)
          co2m(i)=pco2m(i)/(PSUR(i)*100.)   
          CO2A(i)=(CO2A(i)+(dtt/co2cap(i))*   &
          (respg(i)-assimny(i,ic)+co2m(i)*gah2o(i))) &
                  /(1+dtt*gah2o(i)/co2cap(i)) 
!
          pco2a=co2a(i)*psur(i)*100.
          PCO2S(i)=PCO2A-(1.4/GBH2O(i)*ASSIMNy(i,ic)*PSUR(i)*100.)
! new variable name for intercellular CO2
          PCO2INT=PCO2S(i)-(1.6/GSH2O(i)*ASSIMNy(i,ic)*PSUR(i)*100.)
! PCO2IN is now defined as chloroplast CO2 and is used as intercellular
!  co2 was use in the past.
          PCO2IN=PCO2INT(i)-(1./xgco2m(i)*ASSIMNy(i,ic)*PSUR(i)*100.)
          EYY(i,IC)=PCO2Y(i,IC)-PCO2IN
        enddo
!
        if(ic.ge.2) then
	  ic1 = ic-1
          do i = 1,len
            if(abs(eyy(i,ic1)).ge.0.1)then
	      icconv(i) = ic
	    else
	      eyy(i,ic) = eyy(i,ic1)
	      pco2y(i,ic) = pco2y(i,ic1)
	    endif
	  enddo
	endif
!
      enddo
!pl end of PL's setup
!
! Loop over grid points
      do i = 1,len
        icconv(i) = min(icconv(i),6)
        igath(i) = i+(icconv(i)-1)*len
      enddo
!
! Loop over grid points
      do i = 1,len
        pco2i(i) = pco2y(igath(i),1)
        assimn(i) = assimny(igath(i),1)
        assim(i)  = assimy(igath(i),1)
        pco2c(i) = pco2i(i) - assimn(i)/xgco2m(i)*psur(i)*100.0
!
!pl     real C_A forecast with the iterated fluxes.
        CO2A(i)    = PCO2Ap(i) /   (PSUR(i)*100.)
        co2m(i)    = pco2m(i)  /   (PSUR(i)*100.)   
        CO2A(i) = (CO2A(i) + (dtt/co2cap(i)) * &
                    (respg(i) - assimn(i)  &
                    +co2m(i)*gah2o(i) ) ) &
                  / (1+dtt*gah2o(i)/co2cap(i))
!
!pl     convert from molC/mol air to Pascals 
        pco2ap(i) = co2a(i) * psur(i) * 100.
      enddo
!
      dtti = 1./dtt
!
! Loop over grid points
      do i = 1,len
!
        H2OI   = ETC(i)/PSUR(i)
        H2OA   =  EA(i)/PSUR(i)
        ECMOLE = 55.56 * ECMASS(i) * dtti
        H2OS = H2OA + ECMOLE / GBH2O(i)
        H2OS  = MIN( H2OS, H2OI )
        H2OS  = MAX( H2OS, 1. EEEXP -7)
        H2OSRH = H2OS/H2OI
!
!       qdamp and pdamp calculated and passed to here
!pl     This condition relaxed to 1/10 of previous. Old way made
!pl     CO2 on top of leaves always at least 1/2 of value at reference level.
        CO2S(i)=MAX(PCO2S(I),PCO2M(i)*0.05)/(PSUR(i)*100.)
!
!       Ball-Berry stomatal conductance (SE (1992a) eqn 35)     
        GSH2OINF = (GRADM(i)*MAX(1. EEEXP -12,ASSIMN(i))
     >            *H2OSRH*soilfrz(i)/CO2S(i))+BINTC(i)
!
!       change in stomatal resistance
        DRST(i) = RST(i) * QDAMP * ((GSH2O(i)-GSH2OINF)/ &
                  (PDAMP*GSH2O(i)+QDAMP*GSH2OINF))
!
!       change in stomatal resistance without damping factor.
        RSTFAC(i,1) = H2OS/H2OI
        RSTFAC(i,4) = RSTFAC(i,1)*RSTFAC(i,2)*RSTFAC(i,3)
      enddo
!
!-----------------------------------------------------------------------
! CARNEGIE new diagnostics (c.zhang & joe berry, 10/19/92)
!-----------------------------------------------------------------------
!
! Loop over grid points
      do i = 1,len
!-----------------------------------------------------------------------
! Potential assimilation rates
!-----------------------------------------------------------------------
!       Net assimilation at top of canopy
        ASSIMNp(i) = ASSIMN(i) / APARKK(i)
!
!       Bottom stopped assim.
        ANTEMP(i) = MAX(0. EEEXP 0,ASSIMNp(i))
!
!       Potential intercellular co2.
        PCO2IPOT = PSUR(i)*100.*(co2s(i)-(1.6*ASSIMNp(i)/ &
             ((GRADM(i)*ANTEMP(i)/co2s(i))+BINTC(i))))
!
!       Potential rubisco limitation.
        OMCPOT = VMAX0(i)*2.1**qt(i)*((PCO2IPOT-GAMMAS(i))/ &
             (PCO2IPOT+RRKK(i))*C3(i) + C4(i))
!
!       Potential light limitation.
        OMEPOT(i) = PAR(i)*((PCO2IPOT-GAMMAS(i))/ &
                 (PCO2IPOT+2.*GAMMAS(i))*C3(i) + C4(i))
!
!       root portion (a^2-4ac) Quadratic eqn for OMPPOT (SE (1996), eqn C6a)
        SQRTIN = MAX(0. EEEXP 0,((OMEPOT(i)+OMCPOT)**2- &
                     4.*ATHETA(i)*OMEPOT(i)*OMCPOT))
!
!       Quad 1. Intermediate  top leaf photosynthesis.
        OMPPOT = ((OMEPOT(i)+OMCPOT)-SQRT(SQRTIN))/(2.*ATHETA(i))
!
!       Solve Quadratic eqn for AssimPot (SE (1996), eqn C6b)
!       Potential sink or peptide limitation.
        OMSPOT = (VMAX0(i)/2.0)*(1.8**qt(i))*C3(i)  &
                      + RRKK(i)*PCO2IPOT*C4(i)
!
!       Quad 2.
        SQRTIN=MAX(0. EEEXP 0,((OMPPOT+OMSPOT)**2-4.*BTHETA(i)* &
              OMPPOT*OMSPOT))
!
!       Quad 2. Final Potential top leaf photosynthesis.
        ASSIMPOT(i) = ((OMSPOT+OMPPOT)-SQRT(SQRTIN))/(2.*BTHETA(i))
!
!-----------------------------------------------------------------------
! stress factor limited assimilation
!-----------------------------------------------------------------------
!       Stressed rubisco limitation.
        OMCCI=VM(i)*((PCO2IPOT-GAMMAS(i)) &
        /(PCO2IPOT+RRKK(i))*C3(i)+C4(i))
!
!       Quad 1.
        SQRTIN=MAX(0. EEEXP 0,(OMEPOT(i)+OMCCI)**2- &
                      4.*ATHETA(i)*OMEPOT(i)*OMCCI)
!
!       Quad 1. Intermediate stress limited top leaf photosynthesis.
        OMPCI = ((OMEPOT(i)+OMCCI)-SQRT(SQRTIN))/(2.*ATHETA(i))
!
!       Stressed sink or pep limitation.
        OMSCI = OMSS(i)*(C3(i) + PCO2IPOT*C4(i))
!
!       Quad 2.
        SQRTIN=MAX(0. EEEXP 0,(OMPCI+OMSCI)**2-4. &
        *BTHETA(i)*OMPCI*OMSCI)
! 
!       Quad 2. Final stress limited top leaf photosynthesis.
        ASSIMCI(i) = ((OMSCI+OMPCI)-SQRT(SQRTIN))/(2.*BTHETA(i))
!
!-----------------------------------------------------------------------
! CALCULATION OF CONTROL COEFFICIENTS
!-----------------------------------------------------------------------
!       Intermediate
        DOMPDOMC = (OMPCI-OMEPOT(i))/ &
                (2.*ATHETA(i)*OMPCI-OMCCI-OMEPOT(i))
!
!       Bottom stopped final stress limited top leaf photosynthesis
        ASCITEMP = MAX(ASSIMCI(i),1. EEEXP -12)
!
!       Rubisco control coefficient.
        CCOMC = (DOMPDOMC*(ASSIMCI(i)-OMSCI)/ &
            (2.*BTHETA(i)*ASSIMCI(i)-OMPCI-OMSCI))*OMCCI/ASCITEMP
!
!       Sink or pep control coefficient
        CCOMS = ((ASSIMCI(i)-OMPCI)/ &
             (2.*BTHETA(i)*ASSIMCI(i)-OMPCI-OMSCI))*OMSCI/ASCITEMP
!
!-----------------------------------------------------------------------
! Potential assimilation rates
!-----------------------------------------------------------------------
!       Canopy values (overwrites top leaf).
        OMEPOT(i) = OMEPOT(i)*APARKK(i)
        ASSIMPOT(i) = ASSIMPOT(i)*APARKK(i)
        ASSIMCI(i) = ASSIMCI(i)*APARKK(i)
        ASSIM(i) = ASSIM(i)*APARKK(i)
        ANTEMP(i) = ANTEMP(i)*APARKK(i)
        ANSQR(i) = ANTEMP(i)*ANTEMP(i)
        ASSIMNp(i) = ASSIMNp(i)*APARKK(i)
!
!-----------------------------------------------------------------------
! weighted stress factors
!-----------------------------------------------------------------------
!       Water stress
        WSFWS(i) = ASSIMPOT(i)*(1.-RSTFAC(i,2))*(CCOMC+CCOMS)
!
!       High temperature stress
        WSFHT(i) = ASSIMPOT(i)*(1.-1./TEMPH(i))*CCOMC
!
!       Low temperature stress
        WSFLT(i) = ASSIMPOT(i)*(1.-1./TEMPL(i)) &
        *(CCOMS*C3(i)+CCOMC*C4(i))
!
!       protection for wsfws, wsfht, and wsflt from <0 or >>xxx(2/24/93)
        cwsfws = (1.-RSTFAC(i,2))*(CCOMC+CCOMS)
        if(cwsfws.gt.1. .or. cwsfws.lt.0.) wsfws(i)=0.
        cwsfht = (1.-1./TEMPH(i))*CCOMC
        if(cwsfht.gt.1. .or. cwsfht.lt.0.) wsfht(i)=0.
        cwsflt = (1.-1./TEMPL(i))*(CCOMS*C3(i)+CCOMC*C4(i))
        if(cwsflt.gt.1. .or. cwsflt.lt.0.) wsflt(i)=0.
!
!       Intermediate assimilation weighted Ci
        WCI(i) = ANTEMP(i)*PCO2I(i)
!
!       Intermediate assimilation weighted relative humidty stress factor
        WHS(i) = ANTEMP(i)*RSTFAC(i,1)
!
!       Intermediate assimilation weighted stomatal conductance
        WAGS(i) = GSH2O(i)*ANTEMP(i)
!
!       Intermediate evaporation weighted stomatal conductance
!       (Step 1Step 2 after subroutine updat2)
        WEGS(i) = GSH2O(i)
!
      enddo
!
      RETURN
      END
!

