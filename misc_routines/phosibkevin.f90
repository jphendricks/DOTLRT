!
!==================SUBROUTINE PHOSIB===================================
      SUBROUTINE PHOSIB_kevin_old(len,nsib,nsoil,pCO2m,pCO2ap, &
        CAS_cap_CO2,pO2m,ea,etc,psur, &
        tc,ta,tm,td,tgs, &
        Tice,green,aparc,gmudmu,zlt, &
        vmax0,tran,ref,atheta, &
        btheta,Respcp,effcon,binter,gradm, &
        trop,trda,trdm,slti,shti, &
        hltii,hhti,rb,ra,rst, &
        radn,respg,pdamp,qdamp,ecmass, &
        dtt,assimn,assim, &
        Respc,NEE,pfd,rstfac,aparkk, &
        drst,bintc,tprcor,rhs,pCO2s, &
        pCO2i,pCO2c,CO2cap, &
        ScaleTem, &
        ScaleTex,ScaleH2O,ScaleRub,ScalePAR,ScaleExp, &
        ScaleRsp,ome,omc,oms,test, &
        test10)
!=======================================================================
! calculates canopy assimilation, respiration, and CO2 concentration.
!
! CONVERSIONS
!      1 mol H2O           = 0.018 KG
!      1 mol CO2           = 0.044 KG
!      H2O (mol mol-1)     = EA/PSUR ( MB MB-1 )
!      H2O (mol mol-1)     = Q*MM/(Q*MM + 1)
!      GS  (CO2)           = GS (H2O) * 1./1.6 (applies to Ci-Cs pathway)
!pl    44.6 moles of air per cubic meter
!      GS  (mol M-2 S-1 )  = GS (M S-1) * 44.6*Tice/T*P/PO
!      PAR (mol M-2 S-1 )  = PAR(W M-2) * 4.6*1.E-6
!      MM  (mol AIR/mol H2O) = 1.611
!
! Modifications:
!  C.Zhang & Joe Berry added new diagnostics (10/19/92)
!  gs stomatal conduct. reduced for freezing soil by Jim Collatz (02/21/95)
!  DD Modified for multitasking; add gather/scatter indices (12/06/95)
!  Ian Baker Added pCO2c (chloroplast partial CO2) for isotope stuff (Sep 99)
!  Kevin Schaefer added some variable definitions to comments (11/16/00)
!  Kevin Schaefer added back green fraction to PAR use parameter (2/8/01)
!  Kevin Schaefer deleted 8 unused variables, added more definitions (5/4/01)
!  Kevin Schaefer added maximum PAR use parameter (5/14/01)
!  Kevin Schaefer deleted 25 unused variables (6/28/01)
!  Kevin Schaefer replaced code to locate minimum pCO2i error (7/18/01)
!  Kevin Schaefer corrected mixed pCO2s, co2s, and pCO2m if statement (7/18/01)
!  Kevin Schaefer removed C3 from gammas calculation (7/22/01)
!  Kevin Schaefer corrected underflow patch for Respc (7/22/01)
!  Kevin Schaefer Changed Tf to Tice to match rest of SiB (7/22/01)
!  Kevin Schaefer canopy stomatal gsH2O for humid potential rate (7/23/01)
!  Kevin Schaefer changed xgCO2m to xgCO2c to match convention (7/23/01)
!  Kevin Schaefer replaced all diagnostics (7/23/01)
!  Kevin Schaefer added ScaleTex for export temperature scaling (7/29/01)
!  Kevin Schaefer moved pCO2i calculations from cycalc to phosib (7/29/01)
!  Kevin Schaefer added NEE, rhs, and ScaleTex to argument list (7/30/01)
!  Kevin Schaefer absolute value for Acon diagnostics (8/16/01)
!  Kevin Schaefer weighted AconPAR by assim so it goes to 0 @ night (8/16/01)
!  Kevin Schaefer included Vcover in arg list & aparkkmax calc (8/19/01)
!  Kevin Schaefer added soil water stress to zgco2c calculation (8/30/01)
!  Kevin Schaefer switched to pCO2c from pCO2i iteration (8/30/01)
!  Kevin Schaefer added optional switch to iterate pCO2i or pCO2c (10/17/01)
!  Kevin Schaefer corrected gammas calc (gammas=0 for C4 plant) (10/19/01)
!  Kevin Schaefer corrected error in TRPcor calc (Tice+25 to Tice) (11/6/01)
!  Kevin Schaefer corrected scatg calc (1.-t+r to 1-t-r) (11/6/01)
!  Kevin Schaefer changed Q10 base for scaletem from 2.0 to 2.1 (11/8/01)
!  Kevin Schaefer used correct assimn in pCO2c iteration (11/11/01)
!  Kevin Schaefer recalculated pCO2i and pCO2s after iteration loop (11/12/01)
!  Kevin Schaefer Acon stuff from assimn to assim; weight by assim (11/27/01)
!  Kevin Schaefer added RConLAI for respc (1/21/02)
!  Kevin Schaefer changed pCO2s to pCO2c for Aconhum (2/4/02)
!=======================================================================
!
      implicit none
!
! Input dimension variables
      integer len    ! number of SiB points in chunk
      integer nsib   ! total number of SiB points
      integer nsoil  ! total number of soil layers
!
! input time-dependent vegetation variables
      real green(len)     ! (-) greeness fraction of canopy
      real aparc(len)     ! (-) absorbed fraction of incident PAR
      real gmudmu(len)    ! (-) time mean leaf projection divided by cosz
      real zlt(len)       ! (m2/m2) leaf area index
!
! input time-independent vegetation variables
      real vmax0(len)     ! (mol/m2/s) Max Rubisco catalytic capacity at canopy top
      real Respcp(len)    ! (-) autotrophic leaf respiration fraction of VMAX0
      real trop(len)      ! (K) Reference Temp for Q10 response exponent
      real trda(len)      ! (1/K) slope temp canopy resp inhibition function
      real trdm(len)      ! (K) 1/2 pt temp canopy resp inhibition function
      real slti(len)      ! (K) Slope of low temp inhibition function
      real shti(len)      ! (1/K) Slope high temp inhibition function
      real hltii(len)     ! (K) 1/2 pt low temp photosynthesis inhibition function
      real hhti(len)      ! (K) 1/2 pt High temp photosynthesis inhibition function
      real effcon(len)    ! (mol CO2/mol quanta) photosynthesis quantum efficiency
      real binter(len)    ! (mol/m2/s) Ball-Berry Conductance Intercept
      real gradm(len)     ! (TBD) Ball-Berry Conductance Photosynthesis Slope
      real atheta(len)    ! (-) Coupling factor for Rubisco-light assim rates
      real btheta(len)    ! (-) Coupling factor for Rubisco-light-export assim rate
      real tran(nsib,2,2) ! (-) leaf transmittance
      real ref(nsib,2,2)  ! (-) leaf reflectance
!
! input state, prognostic, and diagnostic variables
      real ra(len)        ! (m/s) resistance from canopy air space to PBL
      real rb(len)        ! (m/s) resistance from canopy to canopy air space
      real tc(len)        ! (K) canopy temperature
      real ta(len)        ! (K) canopy air space temperature
      real tm(len)        ! (K)  
      real tgs(len)       ! (K) ground surface temperature
      real Tice           ! (K) freezing temperature
      real td(nsib,nsoil) ! (K) 
      real etc(len)       ! (mb) canopy airspace saturation vapor pressure
      real ea(len)        ! (mb) Canopy airspace water vapor pressure
      real rst(len)       ! (m/s) Stomatal resistance
      real pdamp          ! (TBD)
      real qdamp          ! (TBD)
      real dtt            ! (s) time step
      real pO2m           ! (Pa) partial pressure of O2 in leaf interior
      real respg(len)     ! (mole m-2 m-1) ground respiration
      real psur(len)      ! (mb) surface pressure
      real radn(len,2,2)  ! (W/m2/s) vis/nir direct/diffuse flux at canopy top
      real ecmass(len)    ! (mm) canopy evapotranspiration
      real c3(len)        ! (-) c3 plant fraction per grid cell
      real c4(len)        ! (-) c4 plant fraction per grid cell
!
! Output variables
      real assimn(len)    ! (mol/m2/s) net canopy assimilation
      real assim(len)     ! (mol/m2/s) gross canopy assimilation
      real Respc(len)     ! (mol/m2/s) canopy autotrophic respiration rate
      real NEE(len)       ! (mol/m2/s) Net Ecosystem Exchange
      real bintc(len)     ! (mol/m2/s) min Ball-Berry canopy stomatal conductance
      real pfd(len)       ! (mol/m2/s) leaf normal photon flux density dir/diff PAR
      real drst(len)      ! (m/s) change in stomatal resistance
      real tprcor(len)    ! (K) Temperature correction from STD Temp and Pressure
      real rstfac(len,4)  ! (-) water, humidity, temp, & total scaling factors
      real par(len)       ! (mol/m2/s) PAR used for photosynthesis
      real aparkk(len)    ! (-) PAR use parameter
      real ScaleTem(len)  ! (-) temperature scaling factor for rubisco rate
      real ScaleTex(len)  ! (-) temperature scaling factor for export rate
      real ScaleH2O(len)  ! (-) soil water scaling factor for photosynthesis
      real ScaleRub(len)  ! (-) Rubisco/CO2 absorption factor rubisco limited rate
      real ScalePAR(len)  ! (-) PAR/CO2 absorption factor for PAR limited rate
      real ScaleExp(len)  ! (-) Export/CO2 factor for export limited rate
      real ScaleRsp(len)  ! (-) autotrophic leaf respiration scaling factor
      real ome(len)       ! (mol/m2/s) Light limited assimilation rate
      real omc(len)       ! (mol/m2/s) Rubisco limited assimilation rate
      real oms(len)       ! (mol/m2/s) sink/export limited assimilation rate
!     
! Internal index variables:
      integer i        ! SiB point index
      integer ic             ! Iterative search index
      logical Iter_pCO2i ! iterate pCO2i (true) or pCO2c (false)
!
! internal assimilation variables
      real assimny(len,6) ! (mol/m2/s) iterated net canopy assimilation
      real assimy(len,6)  ! (mol/m2/s) iterated gross canopy assimilation
      real gsh2o(len)     ! (mol/m2/s) water stomatal conductance
      real gsh2oinf       ! (mol/m2/s) Ball-Berry stomatal conductance
      real gaH2O(len)     ! (mol/m2/s) water conductance canopy air space to PBL
      real gbh2o(len)     ! (mol/m2/s) water conductance canopy to canopy air space
      real xgaH2O(len)    ! (TBD) bottom stopped water conductance from CAS to PBL
      real vm(len)        ! (mol/m2/s) Canopy Rusbisco catalytic capacity
      real qt(len)       ! (-) Q10 coefficient for temperature scaling
      real Kc(len)       ! (TBD) Michaelis-Menten constant for CO2 
      real Ko(len)       ! (TBD) O2 inhibition constant
      real spfy(len)     ! (TBD) temp adjusted Rubisco specificity CO2/O2
      real gammas(len)   ! (Pa) CO2 compensation point
      real templ         ! (-) Low temp Vm stress factor
      real temph         ! (-) High temp Vm stress factor
      real scatp(len)     ! (-) leaf scattering coef for PAR
      real scatg(len)     ! (-) green leaf scattering coef for PAR
!
! CO2 Variables
      real pCO2m(len)     ! (Pa) Boundary layer CO2 partial pressure
      real CO2m(len)      ! (mol/mol) Boundary Layer CO2 concentration
      real CAS_cap_CO2(len) ! (mol/m2)canopy air space cap on CO2
      real CO2cap(len)    ! (mol/m2) maximum CO2 concentration in air space
      real pCO2a          ! (Pa) Canopy Air Space CO2 partial pressure
      real CO2a(len)      ! (mol/mol) Canopy Air Space CO2 concentration
      real pCO2ap(len)    ! (Pa) previous canopy air space CO2 partial pressure
      real CO2ap(len)     ! (mol/mol) previous canopy air space CO2 concentration
      real pCO2s(len)     ! (Pa) Leaf surface CO2 partial pressure
      real CO2s(len)      ! (mol/mol) Leaf surface CO2 concentration
      real pCO2i(len)     ! (Pa) stomate CO2 partial pressure
      real range(len)     ! (PA) range of pCO2i (compensation pt to pCO2m)
      real PCO2Y(len,6)   ! (Pa) Stomate CO2 partial pressure from rate eqn
      real EYY(len,6)     ! (Pa) Error between rate-conductance pCO2i 
      real pCO2c(len)     ! (Pa) chloroplast CO2 partial pressure
      real xgCO2c(len)    ! (mol/m2/s) conductance scaling factor chloroplast pCO2
!
! Water related variables
      real soilfrz(len)   ! (-) soil frozen (=0) or thawed (=1) flag
      real ecmole         ! (mol) evapotranspiration in moles
      real h2oa           ! (mol/mol) canopy air space water mixing ratio
      real h2os           ! (mol/mol) leaf surface water mixing ratio
      real rhs(len)       ! (-) relative humidity at leaf surface
      real h2oi           ! (mol/mol) stomate water mixing ratio (saturated)
      real soilfrztd      ! (K) soil dew point
      real soilfrztg      ! (K) soil freezing point
!
      real test, test10(10)   ! test variables
!
! set flag to iterate either pCO2i or pCO2c
      Iter_pCO2i=.false.
!
      do i = 1,len
        ScaleH2O(i)=rstfac(i,2)
!
!       set c3/c4 plant fraction
        IF(EFFCON(i).GT.0.07) then
          C3(i)=1.
        else
          C3(i)=0.
        endif
        C4(i)=1.-C3(i)
!
!       max CO2 concentration in canopy (moles air/m2)
        TPRCOR(i)=(Tice)*PSUR(i)*100./1.013e5
        CO2cap(i)=CAS_cap_CO2(i)*44.6*tprcor(i)/ta(i)
!
!       Average leaf extinction coefficient for PAR
        SCATP(I)=GREEN(i)*(TRAN(i,1,1)+REF(i,1,1)) &
           +(1.-GREEN(i))*(TRAN(i,1,2)+REF(i,1,2))
!
!       Green leaf absorption coefficient for PAR
        SCATG(i)=1.-TRAN(i,1,1)-REF(i,1,1)
!
!       Photon Flux Density: direct/diffuse PAR (visible) normal to LAI
        PFD(i)=4.6E-6*GMUDMU(i)*(RADN(i,1,1)+RADN(i,1,2))
!
!       PAR (visible) available for photosynthesis  (SE(96) eqn C13)
        PAR(i)=pfd(i)*EFFCON(i)*SCATG(i)
!
!       PAR use parameter (SE (1992a) eqn 31; SE(96) eqn C14-C15)
        APARKK(i)=Green(i)*APARC(i)/SQRT(1.-SCATP(i))/GMUDMU(i)
!
!       Q10 coefficient for temp scaling (SE(92a) Table 2; SE(96) eqn C17)
        qt(i)=0.1*(TC(i)-TROP(i))
!
!       Rubisco catalytic capacity of canopy (SE(96) eqn C17)
        TEMPL=1.+EXP(SLTI(i)*(HLTIi(i)-TC(i)))
        TEMPH=1.+EXP(SHTI(i)*(TC(i)-HHTI(i)))
        ScaleTem(i)=2.1**Qt(i)*(C3(i)/TEMPH+C4(i)/TEMPL/TEMPH)
        ScaleTex(i)=1.8**Qt(i)*(C3(i)/TEMPH+C4(i)/TEMPL/TEMPH)
!
        Vm(i)=VMAX0(i)*ScaleTem(i)*ScaleH2O(i)
!
!       Michaelis-Menten constant for CO2 (Kc)(SE(92a) Tbl 2; SE(96) eqn C1)
        Kc(i)=30.*2.1**qt(i)
!
!       O2 inhibition constant (Ko) (SE(92a) Table 2; SE(96) eqn C1/C17)
        Ko(i)=30000.*1.2**qt(i)
!
!       Rubisco specificity CO2/O2 (S)(SE(92a) Tbl 2; SE(96) eqn C1)
        SPFY(i)=2600.*0.57**qt(i)
!
!       CO2 compensation point (GAMMA*) (SE(92a) Tbl 2; SE(96) eqn C1)
        GAMMAS(i)=0.5*pO2m/SPFY(i)*C3(i)
!
!       Possible range of pCO2i: compensation pt to pCO2m
        RANGE(i)=pCO2m(i)*(1.-1.6/GRADM(i))-GAMMAS(i)
!
!       Autotrophic Leaf respiration (SE(96) eqn C8, C17)
!itb    patch to prevent underflow if Tc is too cool...
        if(TC(i) >= TRDM(i)-5.)then
          ScaleRsp(i)=2.**Qt(i)/(1.+EXP(TRDA(i)*(Tc(i)-TRDM(i))))
        else
          ScaleRsp(i)=2.**Qt(i)
        endif
        ScaleRsp(i)=ScaleRsp(i)
        Respc(i)=Vmax0(i)*ScaleH2O(i)*ScaleRsp(i)*APARKK(i) &
        *(Respcp(i)*C3(i)+0.025*C4(i))
!
!       Mixed layer CO2 conversion Pa to ppm
        co2m(i)=pCO2m(i)/(PSUR(i)*100.) 
!
!       Previous CAS CO2 conversion Pa to ppm
        CO2ap(i)=pCO2ap(i)/(PSUR(i)*100.)
!
!       Canopy to CAS Conductance, converting from m/s to mol/m2/s
        GBH2O(i)=0.5/RB(i)*44.6*TPRCOR(i)/TC(i)
!
!       CAS to boundary layer conductance, converting from m/s to mol/m2/s
!pl     modified to get the right units in conservation eqn (m2 s/mol_air)
        gaH2O(i)=1.0/RA(i)*44.6*TPRCOR(i)/TM(i)
!        gaH2O(i)=1./MAX(0.446,gaH2O(i))
!
        xgaH2O(i)=max(0.466,gaH2O(i))
        xgCO2c(i)=4000.0*vmax0(i)*APARKK(i)*RSTFAC(i,2)
!
!       Stomatal Conductance, converting from m/s to mol/m2/s
        GSH2O(i)=1.0/RST(i)*44.6*TPRCOR(i)/TC(i)
!
!       Soil freezing and dew points
        soilfrztg=1.+exp(-1.5*(max(270.,tgs(i))-273.16))
        soilfrztd=1.+exp(-1.5*(max(270.,td(i,nsoil))-273.16))
        soilfrz(i)=max(1./soilfrztg, 1./soilfrztd)
        soilfrz(i)=max(soilfrz(i),0.05)
      enddo ! end first grid point loop
!
! Clear out arrays
      PCO2Y=0.
      EYY=0.
!
! iterate internal/stomate CO2 to match conductance and flux based assimn
      DO IC = 1, 6
!       estimate internal/stomate CO2 (pCO2y)
        CALL SORTIN(EYY,PCO2Y,RANGE,GAMMAS,ic,len)
!
        do i = 1,len
!         Rubisco limited assim rate (OMEGA-C) SE(92a) eqn 11
!         Rubisco/CO2 absorption scaling factor assumes po2i=pO2m
          ScaleRub(i)=C3(i)*((PCO2Y(i,ic)-gammas(i))/ &
          (PCO2Y(i,ic)+Kc(i)*(1+pO2m/Ko(i))))+C4(i)
          OMC(i)=VMAX0(i)*ScaleTem(i)*ScaleH2O(i)*Aparkk(i)*Scalerub(i)
!
!         light limited assim rate (OMEGA-E) SE(92a) eqn 12
          ScalePAR(i)=C3(i)*(PCO2Y(i,ic)-gammas(i))/ &
          (PCO2Y(i,ic)+2*gammas(i))+C4(i)
          OME(i)=PAR(i)*ScalePAR(i)*APARKK(i)
!
!         export limited assim rate (OMEGA-S) SE(92a) eqn 13, 37c
          ScaleExp(i)=C3(i)/2.+C4(i)*2.e4*PCO2Y(i,ic)/Psur(i)/100.
          OMS(i)=Vmax0(i)*Scaletex(i)*ScaleH2O(i)*ScaleExp(i)*APARKK(i)
!
!         smooth between light, Rubisco, and export rates
          call CYCALC(assimy(i,ic),atheta(i),btheta(i), &
          ome(i), omc(i), oms(i))
!
!         net canopy assim rate (An) SE(92a) eqn 14-15
          assimny(i,ic)=assimy(i,ic)-Respc(i)
!
!pl       new CAS CO2 from assimilation rate
          CO2A(i)=(CO2Ap(i)+(dtt/CO2cap(i))*   &
          (respg(i)-assimny(i,ic)+CO2m(i)*gaH2O(i))) &
                  /(1+dtt*gaH2O(i)/CO2cap(i)) 
!
!         convert CAS CO2 to mb
          pCO2a=CO2a(i)*psur(i)*100.
!
!         leaf surface CO2 from CAS CO2
          pCO2s(i)=PCO2A-(1.4/GBH2O(i)*ASSIMNy(i,ic)*PSUR(i)*100.)
          CO2S(i)=pCO2s(i)/(PSUR(i)*100.)
          CO2S(i)=MAX(CO2S(I),CO2M(i)*0.05)
!
!         internal/stomate CO2 from leaf surface CO2
          pCO2i(i)=pCO2s(i)-(1.6/GSH2O(i)*ASSIMNy(i,ic)*PSUR(i)*100.)
!
!         chloroplast CO2 from internal/stomate CO2
          pCO2c(i)=pCO2i(i)-ASSIMNy(i,ic)/xgCO2c(i)*psur(i)*100.0
!
!         error between estimated (pCO2y) and calculated pCO2
          if(Iter_pCO2i) then
            EYY(i,IC)=PCO2Y(i,IC)-pCO2i(i)
          else
            EYY(i,IC)=PCO2Y(i,IC)-pCO2c(i)
          endif
        enddo
      enddo
!
      do i = 1,len
!       assume the last iteration has the lowest error
        if(Iter_pCO2i) then
          pCO2i(i)=pCO2y(i,6)
        else
          pCO2c(i)=pCO2y(i,6)
          pCO2i(i)=pCO2c(i)+ASSIMN(i)/xgCO2c(i)*psur(i)*100.0
        endif
        if(pCO2c(i).lt.0) print*, 'negative CO2',i
        assimn(i)=assimny(i,6)
        assim(i)=assimy(i,6)
        test10(1)=assim(i)*1.e6
        test10(2)=respc(i)*1.e6
        test10(3)=aparkk(i)
        test=assim(i)*1.e6
!
!       recalculate some variables based on converged value of pCO2c/pCO2i
        pCO2s(i)=pCO2i(i)+(1.6/GSH2O(i)*ASSIMN(i)*PSUR(i)*100.)
        ScaleRub(i)=C3(i)*((pCO2c(i)-gammas(i))/ &
        (pCO2c(i)+Kc(i)*(1+pO2m/Ko(i))))+C4(i)
        ScalePAR(i)=C3(i)*(pCO2c(i)-gammas(i))/ &
        (pCO2c(i)+2*gammas(i))+C4(i)
        ScaleExp(i)=C3(i)/2.+C4(i)*2.e4*pCO2c(i)/Psur(i)/100.
!
!       CAS CO2 from fluxes  
        CO2A(i)=(CO2Ap(i)+(dtt/CO2cap(i))*(respg(i)-assimn(i)  &
        +CO2m(i)*gaH2O(i)))/(1+dtt*gaH2O(i)/CO2cap(i))
!
!       Switch current and previous CAS CO2
        pCO2ap(i)=CO2a(i)*psur(i)*100.
!
!       Net Ecosystem Exchange
        NEE(i)=respg(i)-assimn(i)
!
!       internal/stomate water mixing ratio (assumed saturation) 
        H2OI=ETC(i)/PSUR(i)
!
!       CAS water mixing ratio
        H2OA=EA(i)/PSUR(i)
!
!       Water vapor into CAS during time step
        ECMOLE=55.56*ECMASS(i)/dtt
!
!       Leaf surface water mixing ratio
        H2OS=H2OA+ECMOLE/GBH2O(i)
        H2OS=MIN(H2OS,H2OI)
        H2OS=MAX(H2OS,1.e-7)

!
!       Leaf surface relative humidity
        rhs(i)=H2OS/H2OI
!
!       Ball-Berry stomatal conductance (SE(92a) eqn 35)
        BINTC(i)=BINTER(i)*ZLT(i)*GREEN(i)*RSTFAC(i,2)*soilfrz(i)     
        GSH2OINF=(GRADM(i)*MAX(1.e-12,ASSIMN(i)) &
                  *rhs(i)*soilfrz(i)/CO2S(i))+BINTC(i)
!
!       change in stomatal resistance
        DRST(i)=RST(i)*QDAMP*((GSH2O(i)-GSH2OINF)/ &
                  (PDAMP*GSH2O(i)+QDAMP*GSH2OINF))
      enddo
!
      RETURN
      END
!
!
!===================SUBROUTINE short_SORTIN===================================
!
      SUBROUTINE short_SORTIN( EYY, PCO2Y, RANGE, GAMMAS, IC)
!
!=======================================================================
!
!     ARRANGES SUCCESSIVE PCO2/ERROR PAIRS IN ORDER OF INCREASING PCO2.
!       ESTIMATES NEXT GUESS FOR PCO2 USING COMBINATION OF LINEAR AND
!       QUADRATIC FITS.
!
! modifications:
!  Kevin Schaefer deleted one variable (8/22/03)
!  Kevin Schaefer split ic<4 loop to speed up code (8/22/03)
!
!=======================================================================
!
      implicit none
                               
      integer ic                     
      REAL EYY(6), PCO2Y(6), range,gammas

!     work arrays
      REAL eyyi1,eyyi2,eyyi3,eyyis
      REAL eyyisp,pco2yis,pco2yisp
      REAL pco2b
      REAL pco2yi1,pco2yi2,pco2yi3
      REAL aterm, bterm, cterm, pco2yq, cc1, cc2, bc1, bc2, ac1, ac2
      REAL pco2yl, a, b, pmin, emin
      integer is
      logical bitx
      integer i, ix, i1, i2, i3, isp, n, j

!
      IF(IC<4) then
        PCO2Y(1)=GAMMAS+0.5*RANGE
!
        PCO2Y(2)=GAMMAS+RANGE*(0.5-0.3*SIGN(1.,EYY(1)))
!
        PCO2Y(3)=PCO2Y(1)-(PCO2Y(1)-PCO2Y(2)) &
                       /(EYY(1)-EYY(2)+1.E-10)*EYY(1)
        PMIN=MIN(PCO2Y(1),PCO2Y(2))
        EMIN=MIN(EYY(1),EYY(2))
        IF(EMIN.GT.0..AND.PCO2Y(3).GT.PMIN) PCO2Y(3)=GAMMAS
      else
!
      N = IC - 1
      bitx= abs(eyy(n)).gt.0.1
      if(.not.bitx) pco2y(ic) = pco2y(n)

      if(bitx) then
        DO 1000 J = 2, N
          A = EYY(J)
          B = PCO2Y(J)
          DO 2000 I = J-1,1,-1
            IF(EYY(I) .le. A ) go to 100
            EYY(I+1) = EYY(I)
            PCO2Y(I+1) = PCO2Y(I)
2000      CONTINUE
          i = 0
 100      continue
          EYY(I+1) = A
          PCO2Y(I+1) = B
1000    CONTINUE
      endif
!
!-----------------------------------------------------------------------
!
      if(bitx) then
        PCO2B = 0.
        IS    = 1
      endif

      DO 3000 IX = 1, N
          if(bitx) then
            IF( EYY(IX) .LT. 0. )  then
              PCO2B = PCO2Y(IX)
              IS = IX
            endif
          endif
3000    CONTINUE
      if(bitx) then
        I1 = IS-1
        I1 = MAX(1, I1)
        I1 = MIN(N-2, I1)
        I2 = I1 + 1
        I3 = I1 + 2
        ISP   = IS + 1
        ISP = MIN0( ISP, N )
        IS = ISP - 1
        eyyisp = eyy(isp)
        eyyis  = eyy(is)
        eyyi1 = eyy(i1)
        eyyi2 = eyy(i2)
        eyyi3 = eyy(i3)
        pco2yisp = pco2y(isp)
        pco2yis = pco2y(is)
        pco2yi1 = pco2y(i1)
        pco2yi2 = pco2y(i2)
        pco2yi3 = pco2y(i3)
      endif
!
      if(bitx) then

!itb...Neil Suits' patch to check for zero in the denominator...
          if(EYYis /= EYYisp)then
             PCO2YL=PCO2Yis-(PCO2Yis-PCO2Yisp)/(EYYis-EYYisp)*EYYis
          else
             PCO2YL = PCO2Yis * 1.01
          endif
!
!   METHOD USING A QUADRATIC FIT
!
          AC1 = EYYi1*EYYi1 - EYYi2*EYYi2
          AC2 = EYYi2*EYYi2 - EYYi3*EYYi3
          BC1 = EYYi1 - EYYi2
          BC2 = EYYi2 - EYYi3
          CC1 = PCO2Yi1 - PCO2Yi2
          CC2 = PCO2Yi2 - PCO2Yi3

!itb...Neil Suits' patch to prevent zero in denominator...
          if(BC1*AC2-AC1*BC2 /= 0.0 .and. AC1 /= 0.0)then
             BTERM = (CC1*AC2-CC2*AC1)/(BC1*AC2-AC1*BC2)
             ATERM = (CC1-BC1*BTERM)/AC1
             CTERM = PCO2Yi2-ATERM*EYYi2*EYYi2-BTERM*EYYi2
             PCO2YQ= CTERM
             PCO2YQ= MAX( PCO2YQ, PCO2B )
             PCO2Y(IC) = ( PCO2YL+PCO2YQ)/2.
          else
             PCO2Y(IC) = PCO2Y(IC) * 1.01
          endif

        endif
!
      endif
!
      pco2y(ic) = max(pco2y(ic),0.01)
!
      RETURN
      END
!
!
!==================SUBROUTINE PHOSIB===================================
      SUBROUTINE testPHOSIB(len,nsib,nsoil,pCO2m,pCO2ap, &
        CAS_cap_CO2,pO2m,ea,etc,psur, &
        tc,ta,tm,td,tgs, &
        Tice,green,aparc,gmudmu,zlt, &
        vmax0,tran,ref,atheta, &
        btheta,Respcp,effcon,binter,gradm, &
        trop,trda,trdm,slti,shti, &
        hltii,hhti,rb,ra,rst, &
        radn,respg,pdamp,qdamp,ecmass, &
        dtt,assimn,assim, &
        Respc,NEE,pfd,rstfac,aparkk, &
        drst,bintc,tprcor,rhs,pCO2s, &
        pCO2i,pCO2c,CO2cap, &
        ScaleTem, &
        ScaleTex,ScaleRub,ScalePAR,ScaleExp, &
        ScaleRsp,ome,omc,oms,test,test10)
!=======================================================================
! calculates canopy assimilation, respiration, and CO2 concentration.
!
! CONVERSIONS
!      1 mol H2O           = 0.018 KG
!      1 mol CO2           = 0.044 KG
!      H2O (mol mol-1)     = EA/PSUR ( MB MB-1 )
!      H2O (mol mol-1)     = Q*MM/(Q*MM + 1)
!      GS  (CO2)           = GS (H2O) * 1./1.6 (applies to Ci-Cs pathway)
!pl    44.6 moles of air per cubic meter
!      GS  (mol M-2 S-1 )  = GS (M S-1) * 44.6*Tice/T*P/PO
!      PAR (mol M-2 S-1 )  = PAR(W M-2) * 4.6*1.E-6
!      MM  (mol AIR/mol H2O) = 1.611
!
! Modifications:
!  C.Zhang & Joe Berry added new diagnostics (10/19/92)
!  gs stomatal conduct. reduced for freezing soil by Jim Collatz (02/21/95)
!  DD Modified for multitasking; add gather/scatter indices (12/06/95)
!  Ian Baker Added pCO2c (chloroplast partial CO2) for isotope stuff (Sep 99)
!  Kevin Schaefer added some variable definitions to comments (11/16/00)
!  Kevin Schaefer added back green fraction to PAR use parameter (2/8/01)
!  Kevin Schaefer deleted 8 unused variables, added more definitions (5/4/01)
!  Kevin Schaefer added maximum PAR use parameter (5/14/01)
!  Kevin Schaefer deleted 25 unused variables (6/28/01)
!  Kevin Schaefer replaced code to locate minimum pCO2i error (7/18/01)
!  Kevin Schaefer corrected mixed pCO2s, co2s, and pCO2m if statement (7/18/01)
!  Kevin Schaefer removed C3 from gammas calculation (7/22/01)
!  Kevin Schaefer corrected underflow patch for Respc (7/22/01)
!  Kevin Schaefer Changed Tf to Tice to match rest of SiB (7/22/01)
!  Kevin Schaefer canopy stomatal gsH2O for humid potential rate (7/23/01)
!  Kevin Schaefer changed xgCO2m to xgCO2c to match convention (7/23/01)
!  Kevin Schaefer replaced all diagnostics (7/23/01)
!  Kevin Schaefer added ScaleTex for export temperature scaling (7/29/01)
!  Kevin Schaefer moved pCO2i calculations from cycalc to phosib (7/29/01)
!  Kevin Schaefer added NEE, rhs, and ScaleTex to argument list (7/30/01)
!  Kevin Schaefer absolute value for Acon diagnostics (8/16/01)
!  Kevin Schaefer weighted AconPAR by assim so it goes to 0 @ night (8/16/01)
!  Kevin Schaefer included Vcover in arg list & aparkkmax calc (8/19/01)
!  Kevin Schaefer added soil water stress to zgco2c calculation (8/30/01)
!  Kevin Schaefer switched to pCO2c from pCO2i iteration (8/30/01)
!  Kevin Schaefer added optional switch to iterate pCO2i or pCO2c (10/17/01)
!  Kevin Schaefer corrected gammas calc (gammas=0 for C4 plant) (10/19/01)
!  Kevin Schaefer corrected error in TRPcor calc (Tice+25 to Tice) (11/6/01)
!  Kevin Schaefer corrected scatg calc (1.-t+r to 1-t-r) (11/6/01)
!  Kevin Schaefer changed Q10 base for scaletem from 2.0 to 2.1 (11/8/01)
!  Kevin Schaefer used correct assimn in pCO2c iteration (11/11/01)
!  Kevin Schaefer recalculated pCO2i and pCO2s after iteration loop (11/12/01)
!  Kevin Schaefer Acon stuff from assimn to assim; weight by assim (11/27/01)
!  Kevin Schaefer added RConLAI for respc (1/21/02)
!  Kevin Schaefer changed pCO2s to pCO2c for Aconhum (2/4/02)
!=======================================================================
!
      implicit none
!
! Input dimension variables
      integer len    ! number of SiB points in chunk
      integer nsib   ! total number of SiB points
      integer nsoil  ! total number of soil layers
!
! input time-dependent vegetation variables
      real green(len)     ! (-) greeness fraction of canopy
      real aparc(len)     ! (-) absorbed fraction of incident PAR
      real gmudmu(len)    ! (-) time mean leaf projection divided by cosz
      real zlt(len)       ! (m2/m2) leaf area index
!
! input time-independent vegetation variables
      real vmax0(len)     ! (mol/m2/s) Max Rubisco catalytic capacity at canopy top
      real Respcp(len)    ! (-) autotrophic leaf respiration fraction of VMAX0
      real trop(len)      ! (K) Reference Temp for Q10 response exponent
      real trda(len)      ! (1/K) slope temp canopy resp inhibition function
      real trdm(len)      ! (K) 1/2 pt temp canopy resp inhibition function
      real slti(len)      ! (K) Slope of low temp inhibition function
      real shti(len)      ! (1/K) Slope high temp inhibition function
      real hltii(len)     ! (K) 1/2 pt low temp photosynthesis inhibition function
      real hhti(len)      ! (K) 1/2 pt High temp photosynthesis inhibition function
      real effcon(len)    ! (mol CO2/mol quanta) photosynthesis quantum efficiency
      real binter(len)    ! (mol/m2/s) Ball-Berry Conductance Intercept
      real gradm(len)     ! (TBD) Ball-Berry Conductance Photosynthesis Slope
      real atheta(len)    ! (-) Coupling factor for Rubisco-light assim rates
      real btheta(len)    ! (-) Coupling factor for Rubisco-light-export assim rate
      real tran(nsib,2,2) ! (-) leaf transmittance
      real ref(nsib,2,2)  ! (-) leaf reflectance
!
! input state, prognostic, and diagnostic variables
      real ra(len)        ! (m/s) resistance from canopy air space to PBL
      real rb(len)        ! (m/s) resistance from canopy to canopy air space
      real tc(len)        ! (K) canopy temperature
      real ta(len)        ! (K) canopy air space temperature
      real tm(len)        ! (K)  
      real tgs(len)       ! (K) ground surface temperature
      real Tice           ! (K) freezing temperature
      real td(nsib,nsoil) ! (K) 
      real etc(len)       ! (mb) canopy airspace saturation vapor pressure
      real ea(len)        ! (mb) Canopy airspace water vapor pressure
      real rst(len)       ! (m/s) Stomatal resistance
      real pdamp          ! (TBD)
      real qdamp          ! (TBD)
      real dtt            ! (s) time step
      real pO2m           ! (Pa) partial pressure of O2 in leaf interior
      real respg(len)     ! (mole m-2 m-1) ground respiration
      real psur(len)      ! (mb) surface pressure
      real psur_pa(len)   ! (Pa) surface pressure in pascals
      real radn(len,2,2)  ! (W/m2/s) vis/nir direct/diffuse flux at canopy top
      real ecmass(len)    ! (mm) canopy evapotranspiration
!
! Output variables
      real assimn(len)    ! (mol/m2/s) net canopy assimilation
      real assim(len)     ! (mol/m2/s) gross canopy assimilation
      real Respc(len)     ! (mol/m2/s) canopy autotrophic respiration rate
      real NEE(len)       ! (mol/m2/s) Net Ecosystem Exchange
      real bintc(len)     ! (mol/m2/s) min Ball-Berry canopy stomatal conductance
      real pfd(len)       ! (mol/m2/s) leaf normal photon flux density dir/diff PAR
      real drst(len)      ! (m/s) change in stomatal resistance
      real tprcor(len)    ! (K) Temperature correction from STD Temp and Pressure
      real rstfac(len,4)  ! (-) water, humidity, temp, & total scaling factors
      real par(len)       ! (mol/m2/s) PAR used for photosynthesis
      real aparkk(len)    ! (-) PAR use parameter
      real ScaleTem(len)  ! (-) temperature scaling factor for rubisco rate
      real ScaleTex(len)  ! (-) temperature scaling factor for export rate
      real ScaleRub(len)  ! (-) Rubisco/CO2 absorption factor rubisco limited rate
      real ScalePAR(len)  ! (-) PAR/CO2 absorption factor for PAR limited rate
      real ScaleExp(len)  ! (-) Export/CO2 factor for export limited rate
      real ScaleRsp(len)  ! (-) autotrophic leaf respiration scaling factor
      real ome(len)       ! (mol/m2/s) Light limited assimilation rate
      real omc(len)       ! (mol/m2/s) Rubisco limited assimilation rate
      real oms(len)       ! (mol/m2/s) sink/export limited assimilation rate
      real omp            ! (mol/m2/s) combo ome/omc assimilation rate
!     
! Internal index variables:
      integer i        ! SiB point index
      integer ic             ! Iterative search index
      logical Iter_pCO2i ! iterate pCO2i (true) or pCO2c (false)
!
! internal assimilation variables
      real gsh2o(len)     ! (mol/m2/s) water stomatal conductance
      real gsh2oinf       ! (mol/m2/s) Ball-Berry stomatal conductance
      real gaH2O(len)     ! (mol/m2/s) water conductance canopy air space to PBL
      real gbh2o(len)     ! (mol/m2/s) water conductance canopy to canopy air space
      real xgaH2O(len)    ! (TBD) bottom stopped water conductance from CAS to PBL
      real vm(len)        ! (mol/m2/s) Canopy Rusbisco catalytic capacity
      real qt(len)       ! (-) Q10 coefficient for temperature scaling
      real Kc(len)       ! (TBD) Michaelis-Menten constant for CO2 
      real Ko(len)       ! (TBD) O2 inhibition constant
      real spfy(len)     ! (TBD) temp adjusted Rubisco specificity CO2/O2
      real gammas(len)   ! (Pa) CO2 compensation point
      real templ         ! (-) Low temp Vm stress factor
      real temph         ! (-) High temp Vm stress factor
      real temp
      real scatp(len)     ! (-) leaf scattering coef for PAR
      real scatg(len)     ! (-) green leaf scattering coef for PAR
!
! CO2 Variables
      real pCO2m(len)     ! (Pa) Boundary layer CO2 partial pressure
      real CO2m(len)      ! (mol/mol) Boundary Layer CO2 concentration
      real CAS_cap_CO2(len) ! (mol/m2)canopy air space cap on CO2
      real CO2cap(len)    ! (mol/m2) maximum CO2 concentration in air space
      real pCO2a          ! (Pa) Canopy Air Space CO2 partial pressure
      real CO2a(len)      ! (mol/mol) Canopy Air Space CO2 concentration
      real pCO2ap(len)    ! (Pa) previous canopy air space CO2 partial pressure
      real CO2ap(len)     ! (mol/mol) previous canopy air space CO2 concentration
      real pCO2s(len)     ! (Pa) Leaf surface CO2 partial pressure
      real CO2s(len)      ! (mol/mol) Leaf surface CO2 concentration
      real pCO2i(len)     ! (Pa) stomate CO2 partial pressure
      real range(len)     ! (PA) range of pCO2i (compensation pt to pCO2m)
      real PCO2Y(6)       ! (Pa) Stomate CO2 partial pressure from rate eqn
      real EYY(6)         ! (Pa) Error between rate-conductance pCO2i 
      real pCO2c(len)     ! (Pa) chloroplast CO2 partial pressure
      real xgCO2c(len)    ! (mol/m2/s) conductance scaling factor chloroplast pCO2
      real testpCO2(len)
      real sqrtin         ! (-) root portion of quadratic   
!
! Water related variables
      real soilfrz(len)   ! (-) soil frozen (=0) or thawed (=1) flag
      real ecmole         ! (mol) evapotranspiration in moles
      real h2oa           ! (mol/mol) canopy air space water mixing ratio
      real h2os           ! (mol/mol) leaf surface water mixing ratio
      real rhs(len)       ! (-) relative humidity at leaf surface
      real h2oi           ! (mol/mol) stomate water mixing ratio (saturated)
      real soilfrztd      ! (K) soil dew point
      real soilfrztg      ! (K) soil freezing point
!
      real test, test10(10)   ! test variables
!
! set flag to iterate either pCO2i or pCO2c
      Iter_pCO2i=.false.
!
! loop through sib points
      do i = 1,len
!
! convert pressure to pascals
        psur_pa(i)=psur(i)*100.
!
! max CO2 concentration in canopy (moles air/m2)
        TPRCOR(i)=(Tice)*PSUR(i)*100./1.013e5
        CO2cap(i)=CAS_cap_CO2(i)*44.6*tprcor(i)/ta(i)
!
! Mixed layer CO2 conversion Pa to ppm
        co2m(i)=pCO2m(i)/psur_pa(i) 
!
! Previous CAS CO2 conversion Pa to ppm
        CO2ap(i)=pCO2ap(i)/psur_pa(i)
!
! Canopy to CAS Conductance, converting from m/s to mol/m2/s
        GBH2O(i)=0.5/RB(i)*44.6*TPRCOR(i)/TC(i)
!
!
! Stomatal Conductance, converting from m/s to mol/m2/s
        GSH2O(i)=1.0/RST(i)*44.6*TPRCOR(i)/TC(i)
!
! Soil freezing and dew points
        soilfrztg=1.+exp(-1.5*(max(270.,tgs(i))-273.16))
        soilfrztd=1.+exp(-1.5*(max(270.,td(i,nsoil))-273.16))
        soilfrz(i)=max(1./soilfrztg, 1./soilfrztd)
        soilfrz(i)=max(soilfrz(i),0.05)
!
! Average leaf extinction coefficient for PAR
        SCATP(I)=GREEN(i)*(TRAN(i,1,1)+REF(i,1,1)) &
           +(1.-GREEN(i))*(TRAN(i,1,2)+REF(i,1,2))
!
! Green leaf absorption coefficient for PAR
        SCATG(i)=1.-TRAN(i,1,1)-REF(i,1,1)
!
! Photon Flux Density: direct/diffuse PAR (visible) normal to LAI
        PFD(i)=4.6E-6*GMUDMU(i)*(RADN(i,1,1)+RADN(i,1,2))
!
! PAR (visible) available for photosynthesis  (SE(96) eqn C13)
        PAR(i)=pfd(i)*EFFCON(i)*SCATG(i)
!
! PAR use parameter (SE (1992a) eqn 31; SE(96) eqn C14-C15)
        APARKK(i)=Green(i)*APARC(i)/SQRT(1.-SCATP(i))/GMUDMU(i)

! CAS to boundary layer conductance, converting from m/s to mol/m2/s
!pl     modified to get the right units in conservation eqn (m2 s/mol_air)
        gaH2O(i)=1.0/RA(i)*44.6*TPRCOR(i)/TM(i)
! gaH2O(i)=1./MAX(0.446,gaH2O(i))
!
        xgaH2O(i)=max(0.466,gaH2O(i))
        xgCO2c(i)=4000.0*vmax0(i)*APARKK(i)*RSTFAC(i,2)
!
! Q10 coefficient for temp scaling (SE(92a) Table 2; SE(96) eqn C17)
        qt(i)=0.1*(TC(i)-TROP(i))
!
! Michaelis-Menten constant for CO2 (Kc)(SE(92a) Tbl 2; SE(96) eqn C1)
        Kc(i)=30.*2.1**qt(i)
!
! O2 inhibition constant (Ko) (SE(92a) Table 2; SE(96) eqn C1/C17)
        Ko(i)=30000.*1.2**qt(i)
!
! Rubisco specificity CO2/O2 (S)(SE(92a) Tbl 2; SE(96) eqn C1)
        SPFY(i)=2600.*0.57**qt(i)
!
! Autotrophic Leaf respiration temperature scaling factor (SE(96) eqn C8, C17)
!itb    patch to prevent underflow if Tc is too cool...
        if(TC(i) >= TRDM(i)-5.)then
          temp=1.+EXP(TRDA(i)*(Tc(i)-TRDM(i)))
        else
          temp=1.
        endif
        ScaleRsp(i)=2.**Qt(i)/temp
!
!------------------------------------
! c3 plant type eqns
!------------------------------------
        IF(EFFCON(i).GT.0.07) then ! c3 plant type
!
! Rubisco catalytic capacity of canopy (SE(96) eqn C17)
          TEMPH=1.+EXP(SHTI(i)*(TC(i)-HHTI(i)))
          ScaleTem(i)=2.1**Qt(i)/TEMPH
          ScaleTex(i)=1.8**Qt(i)/TEMPH
!
! Rubisco capacity
          Vm(i)=VMAX0(i)*ScaleTem(i)*rstfac(i,2)
!
! CO2 compensation point (GAMMA*) (SE(92a) Tbl 2; SE(96) eqn C1)
          GAMMAS(i)=0.5*pO2m/SPFY(i)
!
! Possible range of pCO2i: compensation pt to pCO2m
          RANGE(i)=pCO2m(i)*(1.-1.6/GRADM(i))-GAMMAS(i)
!
! Autotrophic Leaf respiration (SE(96) eqn C8, C17)
          Respc(i)=Vmax0(i)*rstfac(i,2)*ScaleRsp(i)*APARKK(i)*Respcp(i)
!
! iterate pco2 to match flux and Ball-Berry conductance eqns 
          DO IC = 1, 6
!
! estimate internal/stomate CO2 (pCO2y)
            CALL short_SORTIN(EYY,PCO2Y,RANGE(i),GAMMAS(i),ic)
            testpco2(i)=PCO2Y(ic)
!
! Rubisco limited assim rate (OMEGA-C) SE(92a) eqn 11
! Rubisco/CO2 absorption scaling factor assumes po2i=pO2m
            ScaleRub(i)=(testpco2(i)-gammas(i))/ &
             (testpco2(i)+Kc(i)*(1+pO2m/Ko(i)))
            OMC(i)=VMAX0(i)*ScaleTem(i)*rstfac(i,2)*Scalerub(i) &
             *Aparkk(i)
!
! light limited assim rate (OMEGA-E) SE(92a) eqn 12
            ScalePAR(i)=(testpco2(i)-gammas(i))/ &
             (testpco2(i)+2*gammas(i))
            OME(i)=PAR(i)*ScalePAR(i)*APARKK(i)
!
! export limited assim rate (OMEGA-S) SE(92a) eqn 13, 37c
            ScaleExp(i)=0.5
            OMS(i)=Vmax0(i)*Scaletex(i)*rstfac(i,2)*ScaleExp(i) &
             *APARKK(i)
!
! smooth between light and Rubisco rates
            sqrtin=(OME(i)+OMC(i))**2-4.*ATHETA(i)*OME(i)*OMC(i)
            sqrtin=MAX(0.,sqrtin)
            OMP=((OME(i)+OMC(i))-SQRT(SQRTIN))/(2.*ATHETA(i))
!
! smooth between light, Rubisco, and export rates
            sqrtin=(OMP+OMS(i))**2-4.*BTHETA(i)*OMP*OMS(i)
            sqrtin=MAX(0.,sqrtin)
            assim(i)=((OMS(i)+OMP)-SQRT(SQRTIN))/(2.*BTHETA(i))
!
! net canopy assim rate (An) SE(92a) eqn 14-15
            assimn(i)=assim(i)-Respc(i)
!
!pl new CAS CO2 from assimilation rate
            CO2A(i)=(CO2Ap(i)+(dtt/CO2cap(i))*   &
            (respg(i)-assimn(i)+CO2m(i)*gaH2O(i))) &
                  /(1+dtt*gaH2O(i)/CO2cap(i)) 
!
! convert CAS CO2 to mb
            pCO2a=CO2a(i)*psur_pa(i)
!
! leaf surface CO2 from CAS CO2
            pCO2s(i)=PCO2A-(1.4/GBH2O(i)*ASSIMN(i)*psur_pa(i))
            CO2S(i)=pCO2s(i)/psur_pa(i)
            CO2S(i)=MAX(CO2S(I),CO2M(i)*0.05)
!
! internal/stomate CO2 from leaf surface CO2
            pCO2i(i)=pCO2s(i)-(1.6/GSH2O(i)*ASSIMN(i)*psur_pa(i))
!
! chloroplast CO2 from internal/stomate CO2
            pCO2c(i)=pCO2i(i)-ASSIMN(i)/xgCO2c(i)*psur_pa(i)
!
! error between estimated (pCO2y) and calculated pCO2
            if(Iter_pCO2i) then
              eyy(ic)=testpco2(i)-pCO2i(i)
            else
              eyy(ic)=testpco2(i)-pCO2c(i)
            endif
!
          enddo ! c3 pco2 iteration loop
!
!------------------------------------
! c4 plant type eqns
!------------------------------------
        else ! c4 plant type
!
! Rubisco catalytic capacity of canopy (SE(96) eqn C17)
          TEMPL=1.+EXP(SLTI(i)*(HLTIi(i)-TC(i)))
          TEMPH=1.+EXP(SHTI(i)*(TC(i)-HHTI(i)))
          ScaleTem(i)=2.1**Qt(i)/TEMPL/TEMPH
          ScaleTex(i)=1.8**Qt(i)/TEMPL/TEMPH
!
! Rubisco capacity
          Vm(i)=VMAX0(i)*ScaleTem(i)*rstfac(i,2)
!
! CO2 compensation point (GAMMA*) (SE(92a) Tbl 2; SE(96) eqn C1)
          GAMMAS(i)=0.
!
! Possible range of pCO2i: compensation pt to pCO2m
          RANGE(i)=pCO2m(i)*(1.-1.6/GRADM(i))-GAMMAS(i)
!
! Autotrophic Leaf respiration (SE(96) eqn C8, C17)
          Respc(i)=Vmax0(i)*rstfac(i,2)*ScaleRsp(i)*APARKK(i)*0.025
!
! iterate pco2 to match flux and Ball-Berry conductance eqns 
          DO IC = 1, 6
!
! estimate internal/stomate CO2 (pCO2y)
            CALL short_SORTIN(EYY,PCO2Y,RANGE(i),GAMMAS(i),ic)
            testpco2(i)=PCO2Y(ic)
!
! Rubisco limited assim rate (OMEGA-C) SE(92a) eqn 11
            OMC(i)=VMAX0(i)*ScaleTem(i)*rstfac(i,2)*Aparkk(i)
!
! light limited assim rate (OMEGA-E) SE(92a) eqn 12
            OME(i)=PAR(i)*APARKK(i)
!
! export limited assim rate (OMEGA-S) SE(92a) eqn 13, 37c
            ScaleExp(i)=2.e4*testpco2(i)/psur_pa(i)
            OMS(i)=Vmax0(i)*Scaletex(i)*rstfac(i,2)*ScaleExp(i) &
             *APARKK(i)
!
! smooth between light and Rubisco rates
            sqrtin=(OME(i)+OMC(i))**2-4.*ATHETA(i)*OME(i)*OMC(i)
            sqrtin=MAX(0.,sqrtin)
            OMP=((OME(i)+OMC(i))-SQRT(SQRTIN))/(2.*ATHETA(i))
!
! smooth between light, Rubisco, and export rates
            sqrtin=(OMP+OMS(i))**2-4.*BTHETA(i)*OMP*OMS(i)
            sqrtin=MAX(0.,sqrtin)
            assim(i)=((OMS(i)+OMP)-SQRT(SQRTIN))/(2.*BTHETA(i))
!
! net canopy assim rate (An) SE(92a) eqn 14-15
            assimn(i)=assim(i)-Respc(i)
!
!pl new CAS CO2 from assimilation rate
            CO2A(i)=(CO2Ap(i)+(dtt/CO2cap(i))*   &
            (respg(i)-assimn(i)+CO2m(i)*gaH2O(i))) &
                  /(1+dtt*gaH2O(i)/CO2cap(i)) 
!
! convert CAS CO2 to mb
            pCO2a=CO2a(i)*psur_pa(i)
!
! leaf surface CO2 from CAS CO2
            pCO2s(i)=PCO2A-(1.4/GBH2O(i)*ASSIMN(i)*psur_pa(i))
!
! internal/stomate CO2 from leaf surface CO2
            pCO2i(i)=pCO2s(i)-(1.6/GSH2O(i)*ASSIMN(i)*psur_pa(i))
!
! chloroplast CO2 from internal/stomate CO2
            pCO2c(i)=pCO2i(i)-ASSIMN(i)/xgCO2c(i)*psur_pa(i)
!
! error between estimated (pCO2y) and calculated pCO2
            if(Iter_pCO2i) then
              eyy(ic)=testpco2(i)-pCO2i(i)
            else
              eyy(ic)=testpco2(i)-pCO2c(i)
            endif
!
          enddo ! c4 pco2 iteration loop
!
        endif ! c3/c4 plant type if
!------------------------------------
!
! assume the last iteration has the lowest error
        if(Iter_pCO2i) then
          pCO2i(i)=pCO2y(6)
        else
          pCO2c(i)=pCO2y(6)
          pCO2i(i)=pCO2c(i)+ASSIMN(i)/xgCO2c(i)*psur_pa(i)
        endif
!
! recalculate leaf surface CO2 based on converged value of pCO2c/pCO2i
        pCO2s(i)=pCO2i(i)+(1.6/GSH2O(i)*ASSIMN(i)*psur_pa(i))
        CO2S(i)=pCO2s(i)/psur_pa(i)
        CO2S(i)=MAX(CO2S(I),CO2M(i)*0.05)
!
! CAS CO2 from fluxes  
        CO2A(i)=(CO2Ap(i)+(dtt/CO2cap(i))*(respg(i)-assimn(i)  &
        +CO2m(i)*gaH2O(i)))/(1+dtt*gaH2O(i)/CO2cap(i))
!
! Switch current and previous CAS CO2
        pCO2ap(i)=CO2a(i)*psur_pa(i)
!
! Net Ecosystem Exchange
        NEE(i)=respg(i)-assimn(i)
!
! internal/stomate water mixing ratio (assumed saturation) 
        H2OI=ETC(i)/PSUR(i)
!
! CAS water mixing ratio
        H2OA=EA(i)/PSUR(i)
!
! Water vapor into CAS during time step
        ECMOLE=55.56*ECMASS(i)/dtt
!
! Leaf surface water mixing ratio
        H2OS=H2OA+ECMOLE/GBH2O(i)
        H2OS=MIN(H2OS,H2OI)
        H2OS=MAX(H2OS,1.e-7)
!
! Leaf surface relative humidity
        rhs(i)=H2OS/H2OI
!
! Ball-Berry stomatal conductance (SE(92a) eqn 35)
        BINTC(i)=BINTER(i)*ZLT(i)*GREEN(i)*RSTFAC(i,2)*soilfrz(i)     
        GSH2OINF=(GRADM(i)*MAX(1.e-12,ASSIMN(i)) &
                  *rhs(i)*soilfrz(i)/CO2S(i))+BINTC(i)
!
! change in stomatal resistance
        DRST(i)=RST(i)*QDAMP*((GSH2O(i)-GSH2OINF)/ &
                  (PDAMP*GSH2O(i)+QDAMP*GSH2OINF))
!
        test = omp
        test10(1) = omp
      enddo ! sib point loop
!
      RETURN
      END
