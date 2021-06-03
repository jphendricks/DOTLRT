C==================SUBROUTINE PHOSIB===================================
C
        SUBROUTINE PHOSIB(pco2m,po2m,vmax0,tf,psur,green
       ^,          tran,ref,gmudmu,zlt,tc,trop,trda
       ^,          trdm,slti,shti,hltii,hhti,radn,etc,etgs,wc
       ^,          ea,em,rb,ra,tm
       ^,          effcon,rstfac,binter,gradm,assimn
       ^,          rst,atheta,btheta,tgs,respcp
       ^,          aparkk,len,nsib,
       >     omepot,assimpot,assimci,antemp,assimnp,
       >     wsfws,wsfht,wsflt,wci,whs,
       >     wags,wegs,aparc,pfd,assim,td,www,
       >     wopt,zm,wsat,tg,soilscale,zmstscale,zltrscale,zmlscale,
       >     drst,pdamp,qdamp,ecmass,dtt,bintc,tprcor,soilq10,ansqr,
       >     soilscaleold, nsoil, forcerestore
       .,    nsibpn,nsibgn,nsibtn,xpco2i,xc13o2,xc14o2
       .,    assim12, assim13, assim14, Delta_A, gammas, Vm,pco2i
       .,           del_assim,retroflux, pco2c
       .)
 
        implicit none
 
 c
 C
 C=======================================================================
 C
 C     CALCULATION OF CANOPY CONDUCTANCE USING THE INTEGRATED  
 C     MODEL RELATING ASSIMILATION AND STOMATAL CONDUCTANCE.
 C     UNITS ARE CONVERTED FROM MKS TO BIOLOGICAL UNITS IN THIS ROUTINE.
 C     BASE REFERENCE IS SE-92A
 C
 C                          UNITS
 C                         -------
 C
 C      PCO2M, PCO2A, PCO2I, PO2M                : PASCALS
 C      CO2A, CO2S, CO2I, H2OA, H2OS, H2OA       : MOL MOL-1
 C      VMAX0, RESPN, ASSIM, GS, GB, GA, PFD     : MOL M-2 S-1
 C      EFFCON                                   : MOL CO2 MOL QUANTA-1
 C      GCAN, 1/RB, 1/RA, 1/RST                  : M S-1
 C      EVAPKG                                   : KG M-2 S-1
 C      Q                                        : KG KG-1
 C
 C                       CONVERSIONS
 C                      -------------
 C
 C      1 MOL H2O           = 0.018 KG
 C      1 MOL CO2           = 0.044 KG
 C      H2O (MOL MOL-1)     = EA / PSUR ( MB MB-1 )
 C      H2O (MOL MOL-1)     = Q*MM/(Q*MM + 1)
 C      GS  (CO2)           = GS (H2O) * 1./1.6
 C      GS  (MOL M-2 S-1 )  = GS (M S-1) * 44.6*TF/T*P/PO
 C      PAR (MOL M-2 S-1 )  = PAR(W M-2) * 4.6*1.E-6
 C      MM  (MOLAIR/MOLH2O) = 1.611
 C
 C
 C                         OUTPUT
 C                      -------------
 C
 C      ASSIMN              = CANOPY NET ASSIMILATION RATE
 C      EA                  = CANOPY AIR SPACE VAPOR PRESSURE
 C      1/RST               = CANOPY CONDUCTANCE
 C      PCO2I               = INTERNAL CO2 CONCENTRATION
 C      RESPC               = CANOPY RESPIRATION
 C      RESPG               = GROUND RESPIRATION
 C
 C----------------------------------------------------------------------
 C
 C         RSTFAC(1) ( F(H-S) )               : EQUATION (17,18), SE-92A
 C         RSTFAC(2) ( F(SOIL) )              : EQUATION (12 mod), SE-89
 C         RSTFAC(3) ( F(TEMP) )              : EQUATION (5b)   , CO-92
 C         RSTFAC(4) ( F(H-S)*F(SOIL)*F(TEMP))
 C
 C-----------------------------------------------------------------------
 C
 
 c++++++++++++++++++++++++++++++OUTPUT+++++++++++++++++++++++++++++++++++
 c
 C       ASSIMN         CARBON ASSIMILATION FLUX (MOL M-2 S-1)
 C       RST            CANOPY RESISTANCE (S M-1)
 C       RSTFAC(4)      CANOPY RESISTANCE STRESS FACTORS
 C
 c++++++++++++++++++++++++++DIAGNOSTICS++++++++++++++++++++++++++++++++++
 C
 C       RESPC          CANOPY RESPIRATION (MOL M-2 S-1)
 C       RESPG          GROUND RESPIRATION (MOL M-2 S-1)
 C       PCO2I          CANOPY INTERNAL CO2 CONCENTRATION (MOL MOL-1)
 C       GSH2O          CANOPY CONDUCTANCE (MOL M-2 S-1)
 C       H2OS           CANOPY SURFACE H2O CONCENTRATION (MOL MOL-1)
 c
 c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
 c     Modifications:
 c       - gs (stomatal conductance reduced for freezing soils per Jim
 Collatz
 c         dd 950221     
 c
 c      Modified for multitasking - introduced gather/scatter indices
 c          - DD 951206
 
 C
 c     input arrays:
        integer len, nsib, nsoil, nsibpn,nsibgn,nsibtn
        logical forcerestore
 
        REALX8
       *    vmax0(len),psur(len),green(len),gmudmu(len),
       *    zlt(len),tc(len),trop(len),trda(len),
       *    slti(len),shti(len),hltii(len),hhti(len),
       *    etgs(len),wc(len),em(len),ra(len),rb(len),
       *    cog1(len),cog2(len),tm(len),effcon(len),
       *    binter(len),gradm(len),atheta(len),btheta(len),
       *    tgs(len),respcp(len),tran(nsib,2,2),ref(nsib,2,2),
       *    radn(len,2,2),ecmass(len),trdm(len),etc(len),
       *    aparc(len),www(len,2),tg(len),rst(len)
       *,    xc13o2(len)   ! canopy air C13o2 concentration
       *,    xc14o2(len)   ! canopy air C14o2 concentration
 
        REALX8 pdamp, qdamp, dtt, pco2m, tf, po2m
      
 c     output arrays:
 
 
        REALX8
       *     assimn(len),ea(len),rstfac(len,4),
       *     respc(len),respg(len),drst(len)
 czz new diagostics 10/14/92
 c
 c output arrays
        REALX8
       .     omepot(len),assimpot(len),assimci(len),
       >     assimnp(len),whs(len),antemp(len),
       >     wsfws(len),wsfht(len),wsflt(len),wci(len),
       >     wags(len),wegs(len),pfd(len),
       >     td(nsib,nsoil),zmlscale(len),assim(len),
       >     wopt(len),zm(len),wsat(len),
       >     soilscale(len,nsoil),zmstscale(len,2),zltrscale(len),
       >     tprcor(len),bintc(len),soilq10(len,nsoil),
       >     ansqr(len),soilscaleold(len)
        REALX8 , dimension(len,nsibpn) ::
       .    assim12  ! C12 assimilation
       .,   assim13  ! C13 assimilation
       .,   assim14  ! C14 assimilation
       .,   Delta_A  !
       .,   gammas   ! Gamma star - CO2 compensation point
       .,   Vm       ! current maximal RUBISCO capacity
       .,   del_assim  ! del of assimilate
       .,   retroflux  ! retrodiffused flux of co2 (for 18O)
       .,   pco2c    ! chloroplast CO2 concentration
 
      
 c Input/output arrays:
        REALX8 , dimension(nsib) :: pco2i
        REALX8 , dimension(nsib,nsibpn) ::
       .    xpco2i   ! stomate internal CO2 partial pressure
 
 
 c     work arrays:
 
        REALX8 PCO2Y(len,6), EYY(len,6),assimny(len,6),
       *     assimy(len,6)
        REALX8 c3(len),c4(len),range(len),
       *     aparkk(len),gah2o(len),
       *     gbh2o(len),
       *     par(len),rrkk(len),
       *     omss(len),gsh2o(len),pco2s(len),
       *     templ(len),temph(len),
       *     qt(len),co2s(len),scatp(len),scatg(len),
       *     park(len),respn(len),zkc(len),
       *     zko(len),spfy(len),xgah2o(len),xgsh2o(len)
        integer icconv(len),igath(len)
 
        REALX8 soilfrz(len)
        REALX8
       .     bl, tresp, soilq10td, b(2), woptzm, cwsflt, cwsfht, cwsfws,
       *     ccoms, ccomc, ascitemp, dompdomc, omsci,ompci,omcci,omeci,
       *     omcpot, omppot, sqrtin, omspot, pco2ipot, ohtp2, sttp2,
       *     gsh2oinf, h2osrh, h2os, ecmole, h2oa, h2oi, dtti,
       *     pco2in, pco2a, soilfrztd, soilfrztg, soilq10tg,xs
 
 
 ! Fractionation constants
        REALX8 ,parameter :: xab=2.9    ! canopy air space to leaf
 boundary layer
        REALX8 ,parameter :: xa=4.4     ! leaf boundary layer to stomatal
 cavity
        REALX8 ,parameter :: xas=0.7    ! liquid phase fractionation
        REALX8 ,parameter :: xes=1.1    ! dissolution
 ! C3
        REALX8 ,parameter :: xb=28.2    ! rubisco
 
 ! xr equal to 0.01115 for a delta value of -7.8 per mil
        REALX8 ,parameter :: xr=0.01124*(1.-0.0078)
        REALX8 ,parameter :: xrxx=xr/(1.+xr)
       
 
        integer i, ic1, ic, l
        integer,parameter :: ipn=1
 C
 C----------------------------------------------------------------------
 C
        do i = 1,len                !  LOOP OVER GRIDPOINT
          RESPG(i) = 0.
 C
 C----------------------------------------------------------------------
 C
          IF( EFFCON(i) .GT. 0.07 ) then
            C3(i) = 1.
          else
            c3(i) = 0.
          endif
          C4(i)     = 1. - C3(i)
 C
 C-----------------------------------------------------------------------
 C
 C     CALCULATION OF CANOPY PAR USE PARAMETER.
 C
 C      APARKK      (PI)     : EQUATION (31) , SE-92A
 C-----------------------------------------------------------------------
 C 
          SCATP(I) =     GREEN(i)   * ( TRAN(i,1,1) + REF(i,1,1) )
       &           +( 1.-GREEN(i) ) * ( TRAN(i,1,2) + REF(i,1,2) )
          SCATG(i) = TRAN(i,1,1) + REF(i,1,1)
          PARK(i) = SQRT(1.-SCATP(i)) * GMUDMU(i)
 
 c
 c Collatz-Bounoua commented the calculation of  aparc
 c replaced it with theone calculated in new_mapper.
 c
 cb        APARC(i) = 1. - EXP ( -PARK(i)*ZLT(i) )   ! lahouari
 c
          APARKK(i)   = APARC(i) / PARK(i) * GREEN(i)
 C-----------------------------------------------------------------------
 C
 C     Q-10 AND STRESS TEMPERATURE EFFECTS
 C
 C      QT          (QT)    : TABLE (2)     , SE-92A
 C-----------------------------------------------------------------------
 C
          qt(i) = 0.1*( TC(i) - TROP(i) )
          RESPN(i) = RESPCP(i) * VMAX0(i) * RSTFAC(i,2)
          RESPC(i) = RESPN(i) * 2.0**qt(i)
       *        /( 1. + EXP( TRDA(i)*(TC(i)-TRDM(i))))
          vm(i,ipn) = VMAX0(i) * 2.1**qt(i)
 
          TEMPL(i) = 1. + EXP(SLTI(i)*(HLTIi(i)-TC(i)))
          TEMPH(i) = 1. + EXP(SHTI(i)*(TC(i)-HHTI(i)))
          RSTFAC(i,3) = 1./( TEMPL(i)*TEMPH(i))
          vm(i,ipn)    = vm(i,ipn)/TEMPH(i) * RSTFAC(i,2)*C3(i)
       &      + vm(i,ipn) * RSTFAC(i,2)*RSTFAC(i,3) * C4(i)
 C
 C-----------------------------------------------------------------------
 C
 C     MICHAELIS-MENTEN CONSTANTS FOR CO2 AND O2, CO2/O2 SPECIFICITY,
 C     COMPENSATION POINT      
 C
 C      ZKC          (KC)     : TABLE (2)     , SE-92A
 C      ZKO          (KO)     : TABLE (2)     , SE-92A
 C      SPFY         (S)      : TABLE (2)     , SE-92A
 C      GAMMAS       (GAMMA-*): TABLE (2)     , SE-92A
 C      OMSS         (OMEGA-S): EQUATION (13) , SE-92A
 C      BINTC        (B*ZLT)  : EQUATION (35) , SE-92A
 C-----------------------------------------------------------------------
 C
          ZKC(i) = 30. * 2.1**qt(i)
 
 ! New isotope version - iterates for chloroplast CO2 concentration
 ! and not for stomate Co2 concentration
 ! km changed to 25 pa to compensate for concentration gradient from
 ! stomatal cavity to chloroplast now used explicitly
 !        ZKC(i) = 25. * 2.1**qt(i)
 
          ZKO(i) = 30000. * 1.2**qt(i)
          SPFY(i) = 2600. * 0.57**qt(i)
          gammas(i,ipn) = 0.5 * PO2M/SPFY(i) * C3(i)
          PFD(i)    = 4.6E-6 * GMUDMU(i)*
       *             ( RADN(i,1,1)+RADN(i,1,2) )
 C
          TPRCOR(i) = TF*PSUR(i)*100./1.013 E5
 
          GSH2O(i)  = 1.0/RST(i) * 44.6*TPRCOR(i)/TC(i)
          GBH2O(i)  = 0.5/RB(i) * 44.6*TPRCOR(i)/TC(i)
          GAH2O(i)  = 1.0/RA(i) * 44.6*TPRCOR(i)/TM(i)
          xgah2o(i)  = max(0.446 EEEXP 0 ,GAH2O(i))
 C
          RRKK(i)   = ZKC(i)*( 1. + PO2M/ZKO(i) ) * C3(i)
       &               + VMAX0(i)/5.* ( 1.8**qt(i)) * C4(i)
          PAR(i)    = pfd(i)*EFFCON(i)*( 1.-SCATG(i) )
          soilfrztg = 1.+exp(-1.5 * (max(270.0 EEEXP 0,tgs(i))-273.16) )
          soilfrztd = 1.+exp(-1.5 * (max(270.0 EEEXP 0,td
 (i,nsoil))-273.16) )
          soilfrz(i) = max(1./soilfrztg, 1./soilfrztd)
          soilfrz(i) = max( soilfrz(i), 0.05 EEEXP 0)
          BINTC(i)  = BINTER(i)*ZLT(i)*GREEN(i)*
       &                 RSTFAC(i,2) * soilfrz(i)
          OMSS(i)   = ( VMAX0(i)/2.0 ) * ( 1.8**qt(i) )
       &                  /TEMPL(i) * RSTFAC(i,2) * C3(i)
       &                  + RRKK(i) * RSTFAC(i,2) * C4(i)
 C
 C-----------------------------------------------------------------------
 C
 C     FIRST GUESS IS MIDWAY BETWEEN COMPENSATION POINT AND MAXIMUM
 C     ASSIMILATION RATE.
 C
 C-----------------------------------------------------------------------
          RANGE(i)    = PCO2M * ( 1. - 1.6/GRADM(i) ) -gammas(i,ipn)
        icconv(i) = 1
 
        enddo
       
 C
        DO 1000 IC = 1, 6
          do i = 1,len        ! LOOP OVER GRIDPOINT
            PCO2Y(i,IC) = 0.
            EYY(i,IC) = 0.
        enddo
 1000  CONTINUE
 C
 
        DO 2000 IC = 1, 6
 C
 !DIR$ INLINE
          CALL       SORTIN( EYY, PCO2Y, RANGE, GAMMAS, IC,len )
 C
 
 c        CALL       CYCALC( APARKK, VM, ATHETA, BTHETA,par,
 c     &                     GAMMAS, RESPC, RRKK, OMSS, C3, C4,
 c     &                     PCO2Y(1,IC), assimny(1,ic), assimy(1,ic),
 c     &                     len  )
 
        do i = 1,len
          if(PAR(i).gt.0.) then
             xs=(PCO2Y(i,ic)-gammas(i,ipn))*C3(i)
 
             OMCci=vm(i,ipn) * ( xs/(PCO2Y(i,ic)+RRKK(i) )     + C4(i))
             OMEci=PAR(i)* ( xs/(PCO2Y(i,ic)+2.*gammas(i,ipn)) + C4(i))
 
             xs= (OMEci+OMCci)/2.
             OMPci = ( xs - SQRT(
       .      MAX( 0. EEEXP 0,xs*xs-ATHETA(i)*OMEci*OMCci )))/ATHETA(i)
 
 
             OMSci  = OMSS(i) * ( C3(i) + PCO2Y(i,ic) * C4(i) )
 
             xs=(OMPci+OMSci)/2.
             ASSIMy(i,ic) = ( xs - SQRT(MAX( 0. EEEXP 0 ,
       .               xs*xs-BTHETA(i)*OMPci*OMSci)))/BTHETA(i)
          else
             ASSIMy(i,ic) = 0.
          endif
 
 
          ASSIMNy(i,ic)= ( ASSIMy(i,ic) - RESPC(i)) * APARKK(i)
 
 
 C
          PCO2A = PCO2M - (1.4/MAX(0.446 EEEXP 0,GAH2O(i)) *
       >             (ASSIMNy(i,ic) - RESPG(i))* PSUR(i)*100.)
          PCO2S(i) = PCO2A - (1.4/GBH2O(i) * ASSIMNy(i,ic)
       >            * PSUR(i)*100.)
          PCO2IN = PCO2S(i) - (1.6/GSH2O(i) * ASSIMNy(i,ic)
       >            * PSUR(i)*100.)
 
 
 ! New isotope version - iterates for chloroplast CO2 concentration
 ! and not for stomate Co2 concentration
 
 ! gradient due to resistence from mixed layer to canopy air space
 !        PCO2A(i) = PCO2M(i) - (1.4/MAX(0.446 EEEXP 0,GAH2O(i)) *
 !     >             (ASSIMNy(i,ic) - RESPG(i))* PSUR(i)*100.)
 
 ! gradient due to resistence from canopy air space to leaf surface
 !        PCO2S(i) = PCO2A(i) - (1.4/GBH2O(i) * ASSIMNy(i,ic)
 !     >            * PSUR(i)*100.)
 
 ! gradient due to resistence from leaf surface to stomatal cavity
 !        PCO2I(i) = PCO2S(i) - (1.6/GSH2O(i) * ASSIMNy(i,ic)
 !     >            * PSUR(i)*100.)
 
 ! gradient due to resistence from stomatal cavity to chloroplast
 !        PCO2IN = PCO2I(i) - (1.6/GSH2O(i) * ASSIMNy(i,ic)
 !     >            * PSUR(i)*100.)
 
 
          EYY(i,IC) = PCO2Y(i,IC) - PCO2IN
 
        enddo
 c
          if(ic.ge.2) then
        ic1 = ic-1
            do i = 1,len        ! LOOP OVER GRIDPOINT
              if(abs(eyy(i,ic1)).ge.0.1)then
              icconv(i) = ic
            else
              eyy(i,ic) = eyy(i,ic1)
              pco2y(i,ic) = pco2y(i,ic1)
            endif
          enddo
        endif
 
 C
 2000  CONTINUE
 C
        do i = 1,len        ! LOOP OVER GRIDPOINT
          icconv(i) = min(icconv(i),6)
          igath(i) = i+(icconv(i)-1)*len
        enddo
       
        do i = 1,len         ! LOOP OVER GRIDPOINT
          pco2i(i) = pco2y(igath(i),1)
          assimn(i) = assimny(igath(i),1)
        assim(i) = assimy(igath(i),1)
          xpco2i(i,1)=pco2i(i)
        enddo
 
 
 
 !      ic=1
 !      do i = 1,len
 !        if(PAR(i).gt.0.) then
 !           PCO2IN=xpco2i(i,1)
 !           xs=(PCO2IN-gammas(i,ipn))*C3(i)
 
 !           OMCci=vm(i,ipn) * ( xs/(PCO2IN+RRKK(i) )     + C4(i))
 !           OMEci=PAR(i)* ( xs/(PCO2IN+2.*gammas(i,ipn)) + C4(i))
 
 !           xs= (OMEci+OMCci)/2.
 !           OMPci = ( xs - SQRT(
 !     .      MAX( 0. EEEXP 0,xs*xs-ATHETA(i)*OMEci*OMCci )))/ATHETA(i)
 
 
 !           OMSci  = OMSS(i) * ( C3(i) + PCO2IN * C4(i) )
 
 !           xs=(OMPci+OMSci)/2.
 !           ASSIMy(i,ic) = ( xs - SQRT(MAX( 0. EEEXP 0 ,
 !     .               xs*xs-BTHETA(i)*OMPci*OMSci)))/BTHETA(i)
 !           xgsh2o(i)=GSH2O(i)
 !        else
 
 !           GSH2OINF = (GRADM(i) * MAX(1. EEEXP -12,ASSIMN(i))
 !     >              * H2OSRH * soilfrz(i) / CO2S(i)) + BINTC(i)
 !           BINTC(i)  = BINTER(i)*ZLT(i)*GREEN(i)*
 !     &                 RSTFAC(i,2) * soilfrz(i)
 !           GSH2O(i)  = 1.0/RST(i) * 44.6*TPRCOR(i)/TC(i)
 !           DRST(i) = RST(i) * QDAMP * ((GSH2O(i)-GSH2OINF)/
 !     &              (PDAMP*GSH2O(i)+QDAMP*GSH2OINF))
 
 !           xgsh2o(i)=BINTC(i)
 !           ASSIMy(i,ic) = 0.
 !        endif
 
 !        ASSIMNy(i,ic)= ( ASSIMy(i,ic) - RESPC(i)) * APARKK(i)
 
 !        PCO2S(i) = PCO2M - 1.4 * (
 !     >             (ASSIMNy(i,ic) - RESPG(i))/xgah2o(i)
 !     >             +ASSIMNy(i,ic)            /GBH2O(i)
 !     >                           )* PSUR(i)*100.
 
 !        xpco2i(i,1)= PCO2S(i) - (1.6/GSH2O(i) * ASSIMNy(i,ic)
 !     >            * PSUR(i)*100.)
 !         xpco2i(i,1)=max(0.,xpco2i(i,1))
 
 !         print*,
 !     . i,pco2i(i),xpco2i(i,1),pco2i(i)-xpco2i(i,1),icconv(i),' cccccc'
 !      enddo
 
 !        PCO2A = PCO2M - (1.4/MAX(0.446 EEEXP 0,GAH2O(i)) *
 !     >             (ASSIMNy(i,ic) - RESPG(i))* PSUR(i)*100.)
 !        PCO2S(i) = PCO2A - (1.4/GBH2O(i) * ASSIMNy(i,ic)
 !     >            * PSUR(i)*100.)
 !        PCO2IN = PCO2S(i) - (1.6/GSH2O(i) * ASSIMNy(i,ic)
 !     >            * PSUR(i)*100.)
 
 
 
 
        dtti = 1./dtt
        do i = 1,len        ! LOOP OVER GRIDPOINT
 cczzggrst5 - new code
          H2OI   = ETC(i)/PSUR(i)
          H2OA   =  EA(i)/PSUR(i)
          ECMOLE = 55.56 * ECMASS(i) * dtti  ! ecmass must be computed and
 passed in
          H2OS = H2OA + ECMOLE / GBH2O(i)
          H2OS  = MIN( H2OS, H2OI )
          H2OS  = MAX( H2OS, 1. EEEXP -7)
          H2OSRH = H2OS/H2OI
 C  need qdamp and pdamp calculated and passed to here!
          CO2S(i) = MAX(PCO2S(I),PCO2M*0.5) / (PSUR(i)*100.)
          GSH2OINF = (GRADM(i) * MAX(1. EEEXP -12,ASSIMN(i))
       >            * H2OSRH * soilfrz(i) / CO2S(i)) + BINTC(i)
          DRST(i) = RST(i) * QDAMP * ((GSH2O(i)-GSH2OINF)/
       &            (PDAMP*GSH2O(i)+QDAMP*GSH2OINF))
 c
          RSTFAC(i,1) = H2OS/H2OI
          RSTFAC(i,4) = RSTFAC(i,1)*RSTFAC(i,2)* RSTFAC(i,3)
        enddo
 C
 CZ CARNEGIE new diagnostics----start!!!(c.zhang&joe berry, 10/19/92)
 c-----------------------------------------------------------------------
 C  INPUTS:
 PSUR(i),CO2S,ASSIMN(i),GRADM(i),BINTC(i),VMAX0(i),RRKK(i),C3(i),
 C  
 C4(i),PAR(i),ATHETA(i),BTHETA(i),APARKK(i),OMSS(i),RSTFAC(i,2),TEMPH,
 C  
 TEMPL,RSTFAC(i,1),vm(i,ipn),ASSIM,GSH20(i),EFFCON(i),QT,gammas(i,ipn),
 C    PFD(i)
 C
 
          sttp2 = 73.**0.2
          ohtp2 = 100.**0.2
          do i = 1,len
 C-----------------------------------------------------------------------
 C CALCULATION OF POTENTIAL ASSIMILATION
 C-----------------------------------------------------------------------
 C
 C Make assimn a top leaf, not the canopy.
        ASSIMNp(i) = ASSIMN(i) / APARKK(i)
 C
 C Bottom stopped assim.
        ANTEMP(i) = MAX(0. EEEXP 0,ASSIMNp(i))
 
 
 C
 C Potential intercellular co2.
        PCO2IPOT = PSUR(i)*100.*(co2s(i)-(1.6*ASSIMNp(i)/
       &       ((GRADM(i)*ANTEMP(i)/co2s(i))+BINTC(i))))
 C
 C Potential rubisco limitation.
        OMCPOT = VMAX0(i)*2.1**qt(i)*((PCO2IPOT-gammas(i,ipn))/
       &       (PCO2IPOT+RRKK(i))*C3(i) + C4(i))
 C
 C Potential light limitation.
        OMEPOT(i) = PAR(i)*((PCO2IPOT-gammas(i,ipn))/
       &           (PCO2IPOT+2.*gammas(i,ipn))*C3(i) + C4(i))
 C
 
        xs=(OMEPOT(i)+OMCPOT)/2.
 C Quad 1.
 c      SQRTIN = MAX(0. EEEXP 0,((OMEPOT(i)+OMCPOT)**2-
 c     &               4.*ATHETA(i)*OMEPOT(i)*OMCPOT))
 C
 C Quad 1. Intermediate  top leaf photosynthesis.
 c      OMPPOT = ((OMEPOT(i)+OMCPOT)-SQRT(SQRTIN))/(2.*ATHETA(i))
 
        OMPPOT = (xs-SQRT(MAX(0. EEEXP 0,
       .                  xs*xs-ATHETA(i)*OMEPOT(i)*OMCPOT)))/ATHETA(i)
 
 
 C
 C Potential sink or pep limitation.
        OMSPOT = (VMAX0(i)/2.0)*(1.8**qt(i))*C3(i)
       &                + RRKK(i)*PCO2IPOT*C4(i)
 C
 
        xs=(OMPPOT+OMSPOT)/2.
 C Quad 2.
 c      SQRTIN=MAX(0. EEEXP 0,((OMPPOT+OMSPOT)**2-4.*BTHETA(i)*
 c     &        OMPPOT*OMSPOT))
 C
 C Quad 2. Final Potential top leaf photosynthesis.
 c      ASSIMPOT(i) = ((OMSPOT+OMPPOT)-SQRT(SQRTIN))/(2.*BTHETA(i))
 
        ASSIMPOT(i) = (xs-SQRT(MAX(0. EEEXP 0,
       .                  xs*xs-BTHETA(i)*OMPPOT*OMSPOT)))/BTHETA(i)
 
 
 C
 C-----------------------------------------------------------------------
 C CALCULATION OF STRESS FACTOR LIMITED ASSIMILATION
 C-----------------------------------------------------------------------
 C
 C Stressed rubisco limitation.
        OMCCI = vm(i,ipn)*((PCO2IPOT-gammas(i,ipn))/
       .               (PCO2IPOT+RRKK(i))*C3(i)+ C4(i))
 C
 
        xs=(OMEPOT(i)+OMCCI)/2.
 C Quad 1.
 c      SQRTIN = MAX(0. EEEXP 0,(OMEPOT(i)+OMCCI)**2-
 c     &                4.*ATHETA(i)*OMEPOT(i)*OMCCI)
 C
 C Quad 1. Intermediate stress limited top leaf photosynthesis.
        OMPCI = (xs-SQRT(MAX(0. EEEXP 0,
       .                xs*xs-ATHETA(i)*OMEPOT(i)*OMCCI)))/ATHETA(i)
 
 
 C
 C Stressed sink or pep limitation.
        OMSCI = OMSS(i)*(C3(i) + PCO2IPOT*C4(i))
 C
 
        xs=(OMSCI+OMPCI)/2.
 C Quad 2.
 c      SQRTIN = MAX(0. EEEXP
 0,(OMPCI+OMSCI)**2-4.*BTHETA(i)*OMPCI*OMSCI)
 C
 C Quad 2. Final stress limited top leaf photosynthesis.
 c      ASSIMCI(i) = ((OMSCI+OMPCI)-SQRT(SQRTIN))/(2.*BTHETA(i))
 
        ASSIMCI(i) = (xs-SQRT(MAX(0. EEEXP 0,
       .                  xs*xs-BTHETA(i)*OMPCI*OMSCI)))/BTHETA(i)
 
 C
 C-----------------------------------------------------------------------
 C CALCULATION OF CONTROL COEFFICIENTS
 C-----------------------------------------------------------------------
 C
 C Intermediate.
        DOMPDOMC = (OMPCI-OMEPOT(i))/
       &          (2.*ATHETA(i)*OMPCI-OMCCI-OMEPOT(i))
 C
 C Bottom stopped final stress limited top leaf photosynthesis.
        ASCITEMP = MAX(ASSIMCI(i),1. EEEXP -12)
 C
 C Rubisco control coefficient.
        CCOMC = (DOMPDOMC*(ASSIMCI(i)-OMSCI)/
       &      (2.*BTHETA(i)*ASSIMCI(i)-OMPCI-OMSCI))*OMCCI/ASCITEMP
 C
 C Sink or pep control coefficient.
        CCOMS = ((ASSIMCI(i)-OMPCI)/
       &       (2.*BTHETA(i)*ASSIMCI(i)-OMPCI-OMSCI))*OMSCI/ASCITEMP
 C
 C-----------------------------------------------------------------------
 C  OUTPUT:  POTENTIAL ASSIMILATION RATES TO BE SUMMED
 C-----------------------------------------------------------------------
 C Canopy values (overwrites top leaf).
 C
        OMEPOT(i) = OMEPOT(i)*APARKK(i)
        ASSIMPOT(i) = ASSIMPOT(i)*APARKK(i)
        ASSIMCI(i) = ASSIMCI(i)*APARKK(i)
        ASSIM(i) = ASSIM(i)*APARKK(i)
        ANTEMP(i) = ANTEMP(i)*APARKK(i)
        ANSQR(i) = ANTEMP(i)*ANTEMP(i)
        ASSIMNp(i) = ASSIMNp(i)*APARKK(i)
 
 
 !** C13/C12 discrimination for C3 plants
 !      Delta_A(i)=(xab*pco2a(i)+(xa-xab)*pco2s(i)+(xes+xas-xa)*pco2i(i)+
 !     .            (xb-xes-xas)*pco2c(i))/pco2a(i)
 !
 !      ASSIM13(i) = ASSIM(i)*(1-Delta_A(i)/1000)*xrxx
 !      ASSIM12(i) = ASSIM(i)-ASSIM13(i)
 
 
 !note: 1-delta_A/1000 is the fractionation factor for the discrimination
 !need to add a fractionation equation for C4 photosynthesis
 
 
 
 
 
 C-----------------------------------------------------------------------
 C OUTPUT:  WEIGHTED STRESS FACTORS AND OTHER DIAGNOSTIC OUTPUTS TO BE
 SUMMED
 C-----------------------------------------------------------------------
 C
 C Water stress.
        WSFWS(i) = ASSIMPOT(i)*(1.-RSTFAC(i,2))*(CCOMC+CCOMS)
 C
 C High temperature stress.
        WSFHT(i) = ASSIMPOT(i)*(1.-1./TEMPH(i))*CCOMC
 C
 C Low temperature stress.
        WSFLT(i) = ASSIMPOT(i)*(1.-1./TEMPL(i))*(CCOMS*C3(i)+CCOMC*C4(i))
 c
 c  protection for wsfws, wsfht, and wsflt from <0 or >>xxx(2/24/93)
        cwsfws = (1.-RSTFAC(i,2))*(CCOMC+CCOMS)
        if(cwsfws.gt.1. .or. cwsfws.lt.0.) wsfws(i)=0.
        cwsfht = (1.-1./TEMPH(i))*CCOMC
        if(cwsfht.gt.1. .or. cwsfht.lt.0.) wsfht(i)=0.
        cwsflt = (1.-1./TEMPL(i))*(CCOMS*C3(i)+CCOMC*C4(i))
        if(cwsflt.gt.1. .or. cwsflt.lt.0.) wsflt(i)=0.
 
 C
 C Intermediate assimilation weighted Ci.
        WCI(i) = ANTEMP(i)*PCO2I(i)
 C
 C Intermediate assimilation weighted relative humidty stress factor.
        WHS(i) = ANTEMP(i)*RSTFAC(i,1)
 C
 C Intermediate assimilation weighted stomatal conductance.
        WAGS(i) = GSH2O(i)*ANTEMP(i)
 C
 C Intermediate evaporation weighted stomatal conductance.(Step 1.
 C   Step 2 after subroutine updat2)
        WEGS(i) = GSH2O(i)
 c
 c      bl = (((100.* www(i,1))**0.2 - sttp2)/(sttp2 - ohtp2))**2
 c      bl = min(bl,10. EEEXP 0)
 c      zmlscale(i) = 0.8*0.75**bl + 0.2
 c
 c      zltrscale(i) = (exp(0.0693*(tgs(i)-298.15)))
 c      if (zltrscale(i) .lt. 0.17 ) zltrscale(i)= 0.
 c      zltrscale(i) = zltrscale(i) * zmlscale(i)
 
 c      make these zero until canned from the qp diagnostic, then scrap
        soilscaleold(i) = 0.0
         zmlscale(i) = 0.0
         zltrscale(i) = 0.0
        enddo
 C
 c    respiration diagnostics (changes made by cz, according J. Collatz)
 
        if(.not.forcerestore) then
 
        do i = 1,len
 
 ! new version
        woptzm = (wopt(i)/100.)**zm(i)
 
        b(1) = ( (www(i,1)**zm(i)-woptzm )
       //        (woptzm - 1.))**2.
        b(2) = ( (www(i,2)**zm(i)-woptzm )
       //        (woptzm - 1.))**2.
 
 !      print*,b(1),b(2)
 
 !! old version
 !      woptzm = wopt(i)**zm(i)
 !      print*,woptzm
 !      b(1) = (((100.*www(i,1))**zm(i)-woptzm)/
 !     >     (woptzm - 100.**zm(i)))**2
 !      b(2) = (((100.*www(i,2))**zm(i)-woptzm)/
 !     >     (woptzm - 100.**zm(i)))**2
 !      print*,b(1),b(2)
 
 
        b(1) = min(b(1),10. EEEXP 0)
        b(2) = min(b(2),10. EEEXP 0)
        zmstscale(i,1) = 0.8*wsat(i)**b(1) +0.2
        zmstscale(i,2) = 0.8*wsat(i)**b(2) +0.2
 c
 c     Changed Q10 value for soil respiration from 2.0 to 2.4
 c     following Raich and Schelsinger (1992, Tellus 44B, 81-89),
 c     Scott Denning, 9/14/95
 
        soilQ10(i,nsoil) = exp(0.087547 * (tg(i) - 298.15))
        soilscale(i,nsoil) = soilQ10(i,nsoil) * zmstscale(i,1)
        soilscale(i,1) = 0.0
        soilq10(i,1) = 0.0
        enddo
 
        do l = 3,nsoil
           do i = 1,len
              soilQ10(i,l-1) = exp(0.087547 * (td(i,l) - 298.15))
              soilscale(i,l-1) = soilQ10(i,l-1) * zmstscale(i,2)
           enddo
        enddo
 
        else
 
        do i = 1,len
        woptzm = wopt(i)**zm(i)
        b(2) = (((100.*www(i,2))**zm(i)-woptzm)/
       >     (woptzm - 100.**zm(i)))**2
        b(2) = min(b(2),10. EEEXP 0)
        zmstscale(i,1) = 0.8*wsat(i)**b(2) +0.2
 c
 c     Changed Q10 value for soil respiration from 2.0 to 2.4
 c     following Raich and Schelsinger (1992, Tellus 44B, 81-89),
 c     Scott Denning, 9/14/95
 
        soilQ10(i,1) = exp(0.087547 * (td(i,nsoil) - 298.15))
        soilscale(i,1) = soilQ10(i,1) * zmstscale(i,1)
        enddo
 
        endif
 CZZ carnegie new diagnostics end!!!
 
        RETURN
        END
