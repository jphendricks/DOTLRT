!
!==================SUBROUTINE PHOSIB_pCO2c========================
      SUBROUTINE PHOSIB_pCO2c(len,nsib,nsoil,pCO2m,pCO2ap &
      , CAS_cap_CO2,pO2m,ea,etc,psur &
      , tc,ta,tm,td,tgs &
      , Tice,green,aparc,gmudmu,zlt &
      , vcover,vmax0,tran,ref,atheta &
      , btheta,respcp,effcon,binter,gradm &
      , trop,trda,trdm,slti,shti &
      , hltii,hhti,rb,ra,rst &
      , radn,respg,pdamp,qdamp,ecmass &
      , dtt,assimn,assim,antemp,assimnp &
      , respc,NEE,pfd,rstfac,aparkk &
      , drst,bintc,tprcor,rhs,pCO2s &
      , pCO2i,pCO2c,CO2cap,AconHum,AconLAI &
      , AconTem,AconWat,AconPAR,ScaleTem,ScaleTex &
      , ScaleH2O,ScaleRub,ScalePAR,ScaleExp,ScaleRsp &
      , ome,omc,oms,test,test10)
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
!  Kevin Schaefer corrected underflow patch for respc (7/22/01)
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
!=======================================================================
!
      implicit none
!
! Input dimension variables
      integer  &
        len    ! number of SiB points in chunk &
      , nsib   ! total number of SiB points &
      , nsoil  ! total number of soil layers
!
! input time-dependent vegetation variables
      real &
        green(len)     ! (-) greeness fraction of canopy &
      , aparc(len)     ! (-) absorbed fraction of incident PAR &
      , gmudmu(len)    ! (-) time mean leaf projection divided by cosz &
      , zlt(len)       ! (m2/m2) leaf area index
!
! input time-independent vegetation variables
      real &
        vcover(len)   ! vegetation cover fraction &
      , vmax0(len)     ! (mol/m2/s) Max Rubisco catalytic capacity at canopy top &
      , respcp(len)    ! (-) autotrophic leaf respiration fraction of VMAX0 &
      , trop(len)      ! (K) Reference Temp for Q10 response exponent &
      , trda(len)      ! (1/K) slope temp canopy resp inhibition function &
      , trdm(len)      ! (K) 1/2 pt temp canopy resp inhibition function &
      , slti(len)      ! (K) Slope of low temp inhibition function &
      , shti(len)      ! (1/K) Slope high temp inhibition function &
      , hltii(len)     ! (K) 1/2 pt low temp photosynthesis inhibition function &
      , hhti(len)      ! (K) 1/2 pt High temp photosynthesis inhibition function &
      , effcon(len)    ! (mol CO2/mol quanta) photosynthesis quantum efficiency &
      , binter(len)    ! (mol/m2/s) Ball-Berry Conductance Intercept &
      , gradm(len)     ! (TBD) Ball-Berry Conductance Photosynthesis Slope &
      , atheta(len)    ! (-) Coupling factor for Rubisco-light assim rates &
      , btheta(len)    ! (-) Coupling factor for Rubisco-light-export assim rate &
      , tran(nsib,2,2) ! (-) leaf transmittance &
      , ref(nsib,2,2)  ! (-) leaf reflectance
!
! input state, prognostic, and diagnostic variables
      real &
        ra(len)        ! (m/s) resistance from canopy air space to PBL &
      , rb(len)        ! (m/s) resistance from canopy to canopy air space &
      , tc(len)        ! (K) canopy temperature &
      , ta(len)        ! (K) canopy air space temperature &
      , tm(len)        ! (K)   &
      , tgs(len)       ! (K) ground surface temperature &
      , Tice           ! (K) freezing temperature &
      , td(nsib,nsoil) ! (K)  &
      , etc(len)       ! (mb) canopy airspace saturation vapor pressure &
      , ea(len)        ! (mb) Canopy airspace water vapor pressure &
      , rst(len)       ! (m/s) Stomatal resistance &
      , pdamp          ! (TBD) &
      , qdamp          ! (TBD) &
      , dtt            ! (s) time step &
      , pO2m           ! (Pa) partial pressure of O2 in leaf interior &
      , respg(len)     ! (mole m-2 m-1) ground respiration &
      , psur(len)      ! (mb) surface pressure &
      , radn(len,2,2)  ! (W/m2/s) vis/nir direct/diffuse flux at canopy top &
      , ecmass(len)    ! (mm) canopy evapotranspiration &
      , c3(len)        ! (-) c3 plant fraction per grid cell &
      , c4(len)        ! (-) c4 plant fraction per grid cell
!
! Output variables
      real  &
        assimn(len)    ! (mol/m2/s) net canopy assimilation &
      , assim(len)     ! (mol/m2/s) gross canopy assimilation &
      , antemp(len)    ! (mol/m2/s) net assimilation bottom stopped at zero &
      , assimnp(len)   ! (mol/m2/s) Canopy top Net assimilation &
      , respc(len)     ! (mol/m2/s) canopy autotrophic respiration rate &
      , NEE(len)       ! (mol/m2/s) Net Ecosystem Exchange &
      , AconHum(len)   ! (mol/m2/s) humidity effect on net assimilation &
      , AconLAI(len)   ! (mol/m2/s) LAI effect on net assimilation &
      , AconTem(len)   ! (mol/m2/s) Temperature effect on net assimilation &
      , AconWat(len)   ! (mol/m2/s) soil water effect on net assimilation &
      , AconPAR(len)   ! (mol/m2/s) available light effect on net assimilation &
      , bintc(len)     ! (mol/m2/s) min Ball-Berry canopy stomatal conductance &
      , pfd(len)       ! (mol/m2/s) leaf normal photon flux density dir/diff PAR &
      , drst(len)      ! (m/s) change in stomatal resistance &
      , tprcor(len)    ! (K) Temperature correction from STD Temp and Pressure &
      , rstfac(len,4)  ! (-) water, humidity, temp, & total scaling factors &
      , par(len)       ! (mol/m2/s) PAR used for photosynthesis &
      , aparkk(len)    ! (-) PAR use parameter &
      , APARKKmax(len) ! (-) Theoretical maximum PAR use parameter &
      , ScaleTem(len)  ! (-) temperature scaling factor for rubisco rate &
      , ScaleTex(len)  ! (-) temperature scaling factor for export rate &
      , ScaleH2O(len)  ! (-) soil water scaling factor for photosynthesis &
      , ScaleRub(len)  ! (-) Rubisco/CO2 absorption factor rubisco limited rate &
      , ScalePAR(len)  ! (-) PAR/CO2 absorption factor for PAR limited rate &
      , ScaleExp(len)  ! (-) Export/CO2 factor for export limited rate &
      , ScaleRsp(len)  ! (-) autotrophic leaf respiration scaling factor &
      , ome(len)       ! (mol/m2/s) Light limited assimilation rate &
      , omc(len)       ! (mol/m2/s) Rubisco limited assimilation rate &
      , oms(len)       ! (mol/m2/s) sink/export limited assimilation rate &
      , omp            ! (mol/m2/s) combo ome/omc assimilation rate
!     
! Internal index variables:
      integer i        ! SiB point index &
      , ic             ! Iterative search index
      logical Iter_pCO2i ! iterate pCO2i (true) or pCO2c (false)
!
! internal assimilation variables
      real  &
        assimny(len,6) ! (mol/m2/s) iterated net canopy assimilation &
      , assimy(len,6)  ! (mol/m2/s) iterated gross canopy assimilation &
      , gsh2o(len)     ! (mol/m2/s) water stomatal conductance &
      , gsh2oinf       ! (mol/m2/s) Ball-Berry stomatal conductance &
      , gaH2O(len)     ! (mol/m2/s) water conductance canopy air space to PBL &
      , gbh2o(len)     ! (mol/m2/s) water conductance canopy to canopy air space &
      , xgaH2O(len)    ! (TBD) bottom stopped water conductance from CAS to PBL &
      , vm(len)        ! (mol/m2/s) Canopy Rusbisco catalytic capacity &
      , qt(len)       ! (-) Q10 coefficient for temperature scaling &
      , Kc(len)       ! (TBD) Michaelis-Menten constant for CO2  &
      , Ko(len)       ! (TBD) O2 inhibition constant &
      , spfy(len)     ! (TBD) temp adjusted Rubisco specificity CO2/O2 &
      , gammas(len)   ! (Pa) CO2 compensation point &
      , templ         ! (-) Low temp Vm stress factor &
      , temph         ! (-) High temp Vm stress factor &
      , scatp(len)     ! (-) leaf scattering coef for PAR &
      , scatg(len)     ! (-) green leaf scattering coef for PAR
!
! temporary junk variables used for calculating potential rates
      real &
        Tassimn    ! (mol/m2/s) temporary assimn &
      , Tassim     ! (mol/m2/s) Temporary assim &
      , Trespc(len)     ! (mol/m2/s) temporary Respc &
      , TScaleTem(len)  ! (-) temperature scaling factor for rubisco rate &
      , TScaleTex  ! (-) temperature scaling factor for export rate &
      , TScaleH2O  ! (-) soil water scaling factor for photosynthesis &
      , TScaleRub  ! (-) Rubisco/CO2 absorption factor rubisco limited rate &
      , TScalePAR  ! (-) PAR/CO2 absorption factor for PAR limited rate &
      , TScaleExp  ! (-) Export/CO2 factor for export limited rate &
      , TScaleRsp  ! (-) autotrophic leaf respiration scaling factor &
      , Tome       ! (mol/m2/s) Light limited assimilation rate &
      , Tomc       ! (mol/m2/s) Rubisco limited assimilation rate &
      , Toms       ! (mol/m2/s) sink/export limited assimilation rate &
      , Tpar       ! (mol/m2/s) PAR used for photosynthesis &
      , Tgammas(len)    ! (Pa) CO2 compensation point &
      , TpCO2s     ! (Pa) Leaf surface CO2 partial pressure
!
! CO2 Variables
      real  &
        pCO2m(len)     ! (Pa) Boundary layer CO2 partial pressure &
      , CO2m(len)      ! (mol/mol) Boundary Layer CO2 concentration &
      , CAS_cap_CO2(len) ! (mol/m2)canopy air space cap on CO2 &
      , CO2cap(len)    ! (mol/m2) maximum CO2 concentration in air space &
      , pCO2a          ! (Pa) Canopy Air Space CO2 partial pressure &
      , CO2a(len)      ! (mol/mol) Canopy Air Space CO2 concentration &
      , pCO2ap(len)    ! (Pa) previous canopy air space CO2 partial pressure &
      , CO2ap(len)     ! (mol/mol) previous canopy air space CO2 concentration &
      , pCO2s(len)     ! (Pa) Leaf surface CO2 partial pressure &
      , CO2s(len)      ! (mol/mol) Leaf surface CO2 concentration &
      , pCO2i(len)     ! (Pa) stomate CO2 partial pressure &
      , range(len)     ! (PA) range of pCO2i (compensation pt to pCO2m) &
      , pCO2in         ! (Pa) stomate CO2 partial pressure from conductance &
      , pCO2cn         ! (Pa) chloroplast CO2 partial pressure from conductance &
      , PCO2Y(len,6)   ! (Pa) Stomate CO2 partial pressure from rate eqn &
      , EYY(len,6)     ! (Pa) Error between rate-conductance pCO2i  &
      , pCO2ipot(len)  ! (Pa) Potential stomate CO2 partial pressure &
      , pCO2c(len)     ! (Pa) chloroplast CO2 partial pressure &
      , xgCO2c(len)    ! (mol/m2/s) conductance scaling factor chloroplast pCO2
!
! Water related variables
      real  &
        soilfrz(len)   ! (-) soil frozen (=0) or thawed (=1) flag &
      , ecmole         ! (mol/s) evapotranspiration rate in moles &
      , h2oa           ! (mol/mol) canopy air space water mixing ratio &
      , h2os           ! (mol/mol) leaf surface water mixing ratio &
      , rhs(len)       ! (-) relative humidity at leaf surface &
      , h2oi           ! (mol/mol) stomate water mixing ratio (saturated) &
      , soilfrztd      ! (K) soil dew point &
      , soilfrztg      ! (K) soil freezing point
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
        TPRCOR(i)=(Tice+25.0)*PSUR(i)*100./1.013 E5
        CO2cap(i)=CAS_cap_CO2(i)*44.6*tprcor(i)/ta(i)
!
!       Average leaf extinction coefficient for PAR
        SCATP(I)=GREEN(i)*(TRAN(i,1,1)+REF(i,1,1)) &
           +(1.-GREEN(i))*(TRAN(i,1,2)+REF(i,1,2))
!
!       Green leaf absorption coefficient for PAR
        SCATG(i)=1.-TRAN(i,1,1)+REF(i,1,1)
!
!       Photon Flux Density: direct/diffuse PAR normal to projected leaf area
        PFD(i)=4.6E-6*GMUDMU(i)*(RADN(i,1,1)+RADN(i,1,2))
!
!       PAR available for photosynthesis  (SE(96) eqn C13)
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
        ScaleTem(i)=2.**Qt(i)*(C3(i)/TEMPH+C4(i)/TEMPL/TEMPH)
        ScaleTex(i)=1.8**Qt(i)*(C3(i)/TEMPH+C4(i)/TEMPL/TEMPH)
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
          ScaleRsp(i)=2.**Qt(i)/(1.+EXP(TRDA(i)*(Tc(i)-TRDM(i)))) &
          *(C3(i)*respcp(i)+C4(i)*0.025)
        else
          ScaleRsp(i)=2.**Qt(i)*(C3(i)*respcp(i)+C4(i)*0.025)
        endif
        RESPC(i)=VMAX0(i)*ScaleH2O(i)*ScaleRsp(i)*APARKK(i)
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
        gaH2O(i)=1./MAX(0.446,gaH2O(i))
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
          assimny(i,ic)=assimy(i,ic)-RESPC(i)
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
          pCO2c(i)=pCO2i(i)-assimn(i)/xgCO2c(i)*psur(i)*100.0
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
        endif
        if(pCO2c(i).lt.0.) print*, 'negative pCO2c', &
         '--Yall in a heap o trouble boy'
        assimn(i)=assimny(i,6)
        assim(i)=assimy(i,6)
!
!       Net Ecosystem Exchange
        NEE(i)=respg(i)-assimn(i)
!
!       CAS CO2 from fluxes  
        CO2A(i)=(CO2Ap(i)+(dtt/CO2cap(i))*(respg(i)-assimn(i)  &
        +CO2m(i)*gaH2O(i)))/(1+dtt*gaH2O(i)/CO2cap(i))
!
!       Switch current and previous CAS CO2
        pCO2ap(i)=CO2a(i)*psur(i)*100.
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
! Potential net Assimilation Rates
! 1) Change 1 variable to a max/reference value, keep all others the same
! 2) recalculate net assimilation rate
! 3) subtract potential-actual
!
!-----------------------------------------------------------------------
! Effect of Humidity on net assimilation
! set humidity=1, get potential pCO2i from conductance, 
! recalculate variables that depend on pCO2
      do i = 1,len
!       Bottom stopped assimn
        ANTEMP(i)=MAX(0.,ASSIMN(i))
!
!       Potential pCO2i for relative humidity=1.
        pCO2iPOT(i)=pCO2s(i)-1.6*ANTEMP(i)*PSUR(i)*100./ &
        (GRADM(i)*ANTEMP(i)*PSUR(i)*100./pCO2s(i)+BINTc(i))
!
!       Rubisco limited assim rate (OMEGA-C) SE(92a) eqn 11
!       Rubisco/CO2 absorption scaling factor assumes po2i=pO2m
        TScaleRub=C3(i)*((pCO2iPOT(i)-gammas(i))/ &
        (pCO2iPOT(i)+Kc(i)*(1+pO2m/Ko(i))))+C4(i)
        TOMC=VMAX0(i)*ScaleTem(i)*ScaleH2O(i)*Aparkk(i)*TScalerub
!
!       light limited assim rate (OMEGA-E) SE(92a) eqn 12
        TScalePAR=C3(i)*(pCO2iPOT(i)-gammas(i))/ &
        (pCO2iPOT(i)+2*gammas(i))+C4(i)
        TOME=PAR(i)*TScalePAR*APARKK(i)
!
!       export limited assim rate (OMEGA-S) SE(92a) eqn 13, 37c
        TScaleExp=C3(i)/2.+C4(i)*2.e4*pCO2iPOT(i)/Psur(i)/100.
        TOMS=Vmax0(i)*Scaletex(i)*ScaleH2O(i)*TScaleExp*APARKK(i)
!
!       smooth between light, Rubisco, and export rates
        call CYCALC(Tassim,atheta(i),btheta(i), &
        Tome, Tomc, Toms)
!
!       net potential canopy assim rate (An) SE(92a) eqn 14-15
        Tassimn=Tassim-RESPC(i)
!
!       Effect of Humidity on net assimilation
        AconHum(i)=abs(Tassimn-assimn(i))
      enddo
!
!-----------------------------------------------------------------------
! Effect of leaf mass on net assimilation
! Calculate Potential assimilation rate for max PAR use parameter
      do i = 1,len
!       maximum possible PAR use parameter
        APARKKmax(i)=vcover(i)*0.95/SQRT(1.-TRAN(i,1,1)-REF(i,1,1)) &
        /GMUDMU(i)
!
!       Effect of leaf mass on net assimilation
        AconLAI(i)=abs(assimn(i)*(APARKKmax(i)/APARKK(i)-1.))
      enddo
!
!-----------------------------------------------------------------------
! Effect of Temperature on net assimilation
!  Set Tc=Trop, simplify equations, and iterate for new matching pCO2c
!  Set Michaelis-Menten constant for CO2 Kc=30.
!  Set O2 inhibition constant Ko=30000.
!  Set Rubisco specificity CO2/O2 SPFY=2600.
!
      do i = 1,len
!       CO2 compensation point (GAMMA*) (SE(92a) Tbl 2; SE(96) eqn C1)
        Tgammas(i)=0.5*pO2m/2600.*C3(i)
!
!       Rubisco catalytic capacity of canopy (SE(96) eqn C17)
        TempL=1.+EXP(SLTI(i)*(HLTIi(i)-TROP(i)))
        TempH=1.+EXP(SHTI(i)*(TROP(i)-HHTI(i)))
        TScaleTem(i)=C3(i)/TempH+C4(i)/TempL/TempH
!
!       Autotrophic Leaf respiration (SE(96) eqn C8, C17) (no patch required)
        TScaleRsp=C3(i)*respcp(i)+C4(i)*0.025
        Trespc(i)=VMAX0(i)*ScaleH2O(i)*TScaleRsp*APARKK(i)
      enddo
!
! iterate chloroplast CO2 to match conductance and flux based assimn
! iterate only 4 times to save computer time
      DO IC = 1, 4
!       estimate internal/stomate CO2 (pCO2y) (assume range is same)
        CALL SORTIN(EYY,PCO2Y,RANGE,Tgammas,ic,len)
!
        do i = 1,len
!         Rubisco limited assim rate (OMEGA-C) SE(92a) eqn 11
!         Rubisco/CO2 absorption scaling factor assumes po2i=pO2m
          TScaleRub=C3(i)*((PCO2Y(i,ic)-Tgammas(i))/ &
          (PCO2Y(i,ic)+30.*(1+pO2m/30000.)))+C4(i)
          Tomc=VMAX0(i)*TScaleTem(i)*ScaleH2O(i)*Aparkk(i)*TScaleRub
!
!         light limited assim rate (OMEGA-E) SE(92a) eqn 12
          TScalePAR=C3(i)*(PCO2Y(i,ic)-Tgammas(i))/ &
          (PCO2Y(i,ic)+2*Tgammas(i))+C4(i)
          Tome=PAR(i)*TScalePAR*APARKK(i)
!
!         export limited assim rate (OMEGA-S) SE(92a) eqn 13, 37c
          TScaleExp=C3(i)/2.+C4(i)*2.e4*PCO2Y(i,ic)/Psur(i)/100.
          Toms=Vmax0(i)*TScaleTem(i)*ScaleH2O(i)*TScaleExp*APARKK(i)
!
!         smooth between light, Rubisco, and export rates
          call CYCALC(assimy(i,ic),atheta(i),btheta(i), &
          Tome, Tomc, Toms)
!
!         net canopy assim rate (An) SE(92a) eqn 14-15
          assimny(i,ic)=assimy(i,ic)-Trespc(i)
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
          TpCO2s=PCO2A-(1.4/GBH2O(i)*ASSIMNy(i,ic)*PSUR(i)*100.)
!
!         internal/stomate CO2 from leaf surface CO2
          pCO2iN=TpCO2s-(1.6/GSH2O(i)*ASSIMNy(i,ic)*PSUR(i)*100.)
!
!         chloroplast CO2 from internal/stomate CO2
          pCO2cN=pCO2iN-assimn(i)/xgCO2c(i)*psur(i)*100.0
!
!         error between estimated (pCO2y) and calculated pCO2
          if(Iter_pCO2i) then
            EYY(i,IC)=PCO2Y(i,IC)-pCO2iN
          else
            EYY(i,IC)=PCO2Y(i,IC)-pCO2cN
          endif
        enddo
      enddo
!
      do i = 1,len
!         Effect of Temperature on net assimilation
          AconTem(i)=abs(ASSIMNy(i,4)-assimn(i))
      enddo
!
!-----------------------------------------------------------------------
! Effect of soil water on net assimilation
! set ScaleH2O=1 and recalculate assimilation
      do i = 1,len
!       Autotrophic Leaf respiration (SE(96) eqn C8, C17)
        Trespc(i)=VMAX0(i)*ScaleRsp(i)*APARKK(i)
!
!       Rubisco limited assim rate (OMEGA-C) SE(92a) eqn 11
        Tomc=VMAX0(i)*ScaleTem(i)*Aparkk(i)*Scalerub(i)
!
!       export limited assim rate (OMEGA-S) SE(92a) eqn 13, 37c
        ScaleExp(i)=C3(i)/2.+C4(i)*2.e4*pCO2i(i)/Psur(i)/100.
        Toms=Vmax0(i)*Scaletex(i)*ScaleExp(i)*APARKK(i)
!
!       smooth between light, Rubisco, and export rates
        call CYCALC(Tassim,atheta(i),btheta(i), &
        ome(i), Tomc, Toms)
!
!       net canopy assim rate (An) SE(92a) eqn 14-15
        Tassimn=Tassim-Trespc(i)
!
!       Effect of soil water on net assimilation
        AconWat(i)=abs(Tassimn-assimn(i))
      enddo
!
!-----------------------------------------------------------------------
! Effect of light on net assimilation
! set available light to 200. W/m2/s and recalculate assimilation
      do i = 1,len
!       PAR available for photosynthesis  (SE(96) eqn C13)
        TPAR=4.6E-6*GMUDMU(i)*200.*EFFCON(i)*SCATG(i)
!
!       light limited assim rate (OMEGA-E) SE(92a) eqn 12
        Tome=Tpar*ScalePAR(i)*APARKK(i)
!
!       smooth between light, Rubisco, and export rates
        call CYCALC(Tassim,atheta(i),btheta(i), &
        Tome, omc(i), oms(i))
!
!       net canopy assim rate (An) SE(92a) eqn 14-15
        Tassimn=Tassim-RESPC(i)
!
!       Effect of light on net assimilation (weight by assim)
        AconPAR(i)=assim(i)*abs(Tassimn-assimn(i))
      enddo
!
!-------------------------------------------------------------
! Other Diagnostics
!-------------------------------------------------------------
      do i = 1,len
!       Net assimilation at canopy top
        ASSIMNp(i)=ASSIMN(i)/APARKK(i)
      enddo
!
      RETURN
      END
!
