C==================SUBROUTINE PHOSIB===================================
C
      SUBROUTINE PHOSIB(pco2m,pco2ap,po2m,vmax0,tf,psur,green
     ^,          tran,ref,gmudmu,zlt,cas_cap_co2,tc,ta,trop,trda
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
     >     soilscaleold, nsoil, forcerestore, respg, pco2c, pco2i,
     >     pco2s,co2cap,cflux)

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
C      PCO2M, PCO2A, PCO2Ap, PCO2I, PO2M        : PASCALS
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
cpl the next line applies to the Ci to Cs pathway
C      GS  (CO2)           = GS (H2O) * 1./1.6
cpl 44.6 is the number of moles of air per cubic meter
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
c       - gs (stomatal conductance reduced for freezing soils per Jim Collatz
c         dd 950221      
c
c      Modified for multitasking - introduced gather/scatter indices
c          - DD 951206
c
citb   Added in pco2c (chloroplast partial co2) for neil's fractionation
citb   calculations
citb       - IB Sep99

C
c     input arrays:
      integer len, nsib, nsoil
      logical forcerestore

      REALX8 vmax0(len),psur(len),green(len),gmudmu(len),
     *    zlt(len),cas_cap_co2(len),tc(len),ta(len),
     *    trop(len),trda(len),
     *    slti(len),shti(len),hltii(len),hhti(len),
     *    etgs(len),wc(len),em(len),ra(len),rb(len),
     *    cog1(len),cog2(len),tm(len),effcon(len),
     *    binter(len),gradm(len),atheta(len),btheta(len),
     *    tgs(len),respcp(len),tran(nsib,2,2),ref(nsib,2,2),
     *    radn(len,2,2),ecmass(len),trdm(len),etc(len),
     *    aparc(len),www(len,2),tg(len),rst(len), cflux(len)
      REALX8 pdamp, qdamp, dtt, pco2m(len), pco2ap(len), tf, po2m
     
c     output arrays:


      REALX8 assimn(len),ea(len),rstfac(len,4),
     *     pco2i(len),respc(len),respg(len),drst(len)
czz new diagostics 10/14/92 
c
c output arrays
      REALX8 omepot(len),assimpot(len),assimci(len),
     >     assimnp(len),whs(len),antemp(len),
     >     wsfws(len),wsfht(len),wsflt(len),wci(len),
     >     wags(len),wegs(len),pfd(len),
     >     td(nsib,nsoil),zmlscale(len),assim(len),
     >     wopt(len),zm(len),wsat(len),
     >     soilscale(len,nsoil+1),zmstscale(len,2),zltrscale(len),
     >     tprcor(len),bintc(len),soilq10(len,nsoil+1),
     >     ansqr(len),soilscaleold(len)
     *     ,pco2c(len) !chloroplast pco2
     *     ,xgah2o(len) 
     *     ,xgco2m(len)
     
c     work arrays:

      REALX8 PCO2Y(len,6), EYY(len,6),assimny(len,6),
     *     assimy(len,6)
      REALX8 c3(len),c4(len),range(len),gammas(len),
     *     aparkk(len),gah2o(len),
     *     gbh2o(len),
     *     par(len),rrkk(len),
     *     omss(len),vm(len),gsh2o(len),pco2s(len),
     *     templ(len),temph(len),
     *     qt(len),co2s(len),scatp(len),scatg(len),
     *     park(len),respn(len),zkc(len),resa(len),
     *     zko(len),spfy(len), co2a(len), co2m(len),co2cap(len)

      integer icconv(len),igath(len)

      REALX8 soilfrz(len)
      REALX8 bl, tresp, soilq10td, b(2), woptzm, cwsflt, cwsfht, cwsfws,
     *     ccoms, ccomc, ascitemp, dompdomc, omsci, ompci, omcci,
     *     omcpot, omppot, sqrtin, omspot, pco2ipot, ohtp2, sttp2,
     *     gsh2oinf, h2osrh, h2os, ecmole, h2oa, h2oi, dtti,
     *     pco2in, pco2a, soilfrztd, soilfrztg, soilq10tg, 
     *     rstnew

      integer i, ic1, ic, l

      
      
cpl introduce a co2 capacity 
cpl this will basically be the mass of air under the top of the canopy (in
cpl this case (CHEAS-RAMS) O(10-30m), that is, ground to displacemnt height.

cpl all the carbon fluxes are expresse as Mol C / m2 s and resistances for
cpl carbon are in m2 s / mol air

cpl one mole of gas occupies 22.4 cubic dm
cpl 1 cubic meter contains therefore 1000./22.4  = 44.6 moles of gas
cpl the units for the carbon capacity are mol air /m2. 
cpl (e.g. here 893 moles if thickness of the layer is 20m)
cpl this means that the units for pc02ap should be mol co2 / mol air, but
cpl it is also possible to keep just co2 pressure and convert


      do i = 1,len                !  LOOP OVER GRIDPOINT


        TPRCOR(i) = TF*PSUR(i)*100./1.013 E5
        co2cap(i) = cas_cap_co2(i) * 44.6 * tprcor(i)/ta(i)     ! moles air / m2

cpl this needs to be modified as in sibslv3 to automatically use the 
cpl thickness of the canopy air space. 
      
C
C----------------------------------------------------------------------
C
cpl        RESPG(i) = 0. EEEXP -6 ! fixed respiration at 5 micromoles
cpl   no longe rused since we now have respsib
C
C----------------------------------------------------------------------
C
        IF( EFFCON(i) .GT. 0.07 ) then
          C3(i) = 1.
        else
          C3(i) = 0.
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
        SCATP(I) =     GREEN(i)   * 
     &           ( TRAN(i,1,1) + REF(i,1,1) )
     &           +( 1.-GREEN(i) ) *
     &           ( TRAN(i,1,2) + REF(i,1,2) )
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

citb...patch to prevent underflow if temp is too cool...
        if(TC(i) >= TRDM(i))then
            RESPC(i) = RESPN(i) * 2.0**qt(i)
     *        /( 1. + EXP( TRDA(i)*(TC(i)-TRDM(i))))
        else
            RESPC(i) = RESPN(i) * 2.0**qt(i)
        endif

        VM(i) = VMAX0(i) * 2.1**qt(i)
        TEMPL(i) = 1. + EXP(SLTI(i)*(HLTIi(i)-TC(i)))
        TEMPH(i) = 1. + EXP(SHTI(i)*(TC(i)-HHTI(i)))
        RSTFAC(i,3) = 1./( TEMPL(i)*TEMPH(i))
        VM(i)    = VM(i)/TEMPH(i) * RSTFAC(i,2)*C3(i)
     &      + VM(i) * RSTFAC(i,2)*RSTFAC(i,3) * C4(i)
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
        ZKO(i) = 30000. * 1.2**qt(i)
        SPFY(i) = 2600. * 0.57**qt(i)
        GAMMAS(i) = 0.5 * PO2M/SPFY(i) * C3(i)
        PFD(i)    = 4.6E-6 * GMUDMU(i)*
     *             ( RADN(i,1,1)+RADN(i,1,2) )
C
cpl these here all go from being m/s to being mol/ (m2 sec)
        GSH2O(i)  = 1.0/RST(i) * 44.6*TPRCOR(i)/TC(i)
        GBH2O(i)  = 0.5/RB(i) * 44.6*TPRCOR(i)/TC(i)
        GAH2O(i)  = 1.0/RA(i) * 44.6*TPRCOR(i)/TM(i)

        xgah2o(i) = max(0.466 EEEXP 0, gah2o(i) )
        xgco2m(i) = 4000.0 * vmax0(i)
C
citb...this is changed slightly from older version of code...
        RRKK(i)   = ZKC(i)*( 1. + PO2M/ZKO(i) ) * C3(i)
     &               + VMAX0(i)/5.* ( 1.8**qt(i)) * C4(i)
c        RRKK(i)   = ZKC(i)*( 1. + PO2M/ZKO(i) ) * C3(i)
c     &               + VMAX0(i)*200.0/psur(i)* ( 1.8**qt(i)) * C4(i)
        PAR(i)    = pfd(i)*EFFCON(i)*( 1.-SCATG(i) )
        soilfrztg = 1.+exp(-1.5 * (max(270.0 EEEXP 0,tgs(i))-273.16) )
        soilfrztd = 1.+exp(-1.5 * (max(270.0 EEEXP 0,td (i,nsoil))-273.16) )
        soilfrz(i) = max(1./soilfrztg, 1./soilfrztd)
        soilfrz(i) = max( soilfrz(i), 0.05 EEEXP 0)
        BINTC(i)  = BINTER(i)*ZLT(i)*GREEN(i)*
     &                 RSTFAC(i,2) * soilfrz(i)
citb...this is changed slightly from older version of code...
        OMSS(i)   = ( VMAX0(i)/2.0 ) * ( 1.8**qt(i) )
     &                  /TEMPL(i) * RSTFAC(i,2) * C3(i)
     &                  + RRKK(i) *RSTFAC(i,2) * C4(i)
c        OMSS(i)   = ( VMAX0(i)/2.0 ) * ( 1.8**qt(i) )
c     &                  /TEMPH(i) * RSTFAC(i,2) * C3(i)
c     &                  + RRKK(i) * RSTFAC(i,3)*RSTFAC(i,2) * C4(i)
C
C-----------------------------------------------------------------------
C
C     FIRST GUESS IS MIDWAY BETWEEN COMPENSATION POINT AND MAXIMUM
C     ASSIMILATION RATE.
C
C-----------------------------------------------------------------------
 
       
        RANGE(i)    = PCO2M(i) * ( 1. - 1.6/GRADM(i) ) - GAMMAS(i)
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

cpl beginning of PL's setup

c      do i=1,len
c       gah2o(i) =  1. / MAX(0.446 EEEXP 0,GAH2O(i))
c      enddo

cBio...HERE IS THE ITERATION LOOP.
cBio...
cBio...We iterate on PCO2C-sortin makes a 'first guess' at
cBio...then orders PCO2C/Error pairs on increasing error size,
cBio...then uses a combination of linear/quadratic fit to obtain 
cBio...the 'next best guess' as iteration count increases.
cBio...CYCALC uses that value of PCO2C to get the next value 
cBio...of ASSIMN. CO2A and CO2S follow.

      DO 2000 IC = 1, 6
C
        CALL       SORTIN( EYY, PCO2Y, RANGE, GAMMAS, ic,len )
        
        CALL       CYCALC( APARKK, VM, ATHETA, BTHETA,par,
     &                     GAMMAS, RESPC, RRKK, OMSS, C3, C4,
     &                     PCO2Y(1,ic), assimny(1,ic), assimy(1,ic),
     &                     len  )
C

      do i = 1,len

cpl first diagnose the current CO2 flux in mol / (m2 s)
cpl this is a modified ra that will get us the right units
cpl in the conservation equation. its units are m2 s / (mol_air)

c	resa(i) =   1. / MAX(0.446 EEEXP 0,GAH2O(i))

c	resb(i) =   1.4/GBH2O(i)

c	resc(i) =   1.6/GSH2O(i) 
     
cpl now prognose the new CAS CO2 according to flux divergence
cpl we are going to do this in mol C / mol air (same as PaC/PaAir)

        CO2A(i)    = PCO2Ap(i) /   (PSUR(i)*100.)
        co2m(i)    = pco2m(i)  /   (PSUR(i)*100.)   

        CO2A(i)   = (  CO2A(i) + (dtt/co2cap(i)) *  
     &              (respg(i) - assimny(i,ic)
     &              +co2m(i)*gah2o(i)        ) )
     &            / (1+dtt*gah2o(i)/ co2cap(i) ) 

        
c        PCO2A   =  PCO2Ap(i) + (dtt/co2cap) 
c     &            * (  (PSUR(i) *100. * respg(i))
c     &            +    (PCO2Y(i,IC)/(resb+resc))
c     &            +    (pco2m(i)   /       resa)  )
c     &            / (1.+ (dtt/co2cap(i)) * (1./resa + 1./(resb+resc)) )
      

cpl        PCO2A = PCO2M(i) - (1.4/MAX(0.446 EEEXP 0,GAH2O(i)) * 
cpl     >             (ASSIMNy(i,ic) - RESPG(i))* PSUR(i)*100.)


        pco2a = co2a(i) * psur(i) * 100.

        PCO2S(i) = PCO2A - (1.4/GBH2O(i) * ASSIMNy(i,ic)
     >            * PSUR(i)*100.)
c        PCO2IN   = PCO2S(i) - (1.6/GSH2O(i) * ASSIMNy(i,ic)
c     >            * PSUR(i)*100.)
changing to iterate on pco2c instead of pco2i
        PCO2IN = PCO2S(i) - ASSIMNy(i,ic)*psur(i)*100.0*
     &                          (1.6/gsh2o(i) +1.0/xgco2m(i))
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
C
      
cpl end of PL's setup

      do i = 1,len        ! LOOP OVER GRIDPOINT
        icconv(i) = min(icconv(i),6)
        igath(i) = i+(icconv(i)-1)*len
      enddo


      do i = 1,len         ! LOOP OVER GRIDPOINT

        pco2c(i) = pco2y(igath(i),1)
        assimn(i) = assimny(igath(i),1)
        assim(i)  = assimy(igath(i),1)


	
c        pco2i(i)  = pco2y  (i,icconv(i))
c        assimn(i) = assimny(i,icconv(i))
c	assim(i)  = assimy (i,icconv(i))

         pco2i(i) = pco2c(i) + assimn(i)/xgco2m(i)*psur(i)*100.0
         pco2s(i) = pco2i(i) + assimn(i)/gsh2o(i)*psur(i)*100.0


c        pco2c(i) = pco2i(i) - assimn(i)/xgco2m(i)*psur(i)*100.0

cpl now do the real C_A forecast with the iterated fluxes.

        CO2A(i)    = PCO2Ap(i) /   (PSUR(i)*100.)
        co2m(i)    = pco2m(i)  /   (PSUR(i)*100.)   
   
        CO2A(i) = (CO2A(i) + (dtt/co2cap(i)) *
     &              (respg(i) - assimn(i) 
     &              +co2m(i)*gah2o(i) ) )
     &            / (1+dtt*gah2o(i)/co2cap(i))
cpl go back from molC / mol air to Pascals 

        pco2ap(i) = co2a(i) * psur(i) * 100.

citb...carbon flux between CAS and reference level
         cflux(i) = gah2o(i)*(co2a(i)-co2m(i))

c        PCO2Ap(i) =  PCO2Ap(i) + (dtt/co2cap) 
c     &            * (  (PSUR(i) *100. * respg(i))
c     &            +    (PCO2i(i)   /(resb+resc))
c     &            +    (pco2m(i)   /       resa)  )
c     &            / (1.+ (dtt/co2cap) * (1./resa + 1./(resb+resc)) )
      enddo

C
      dtti = 1./dtt
      do i = 1,len        ! LOOP OVER GRIDPOINT
cczzggrst5 - new code
        H2OI   = ETC(i)/PSUR(i)
        H2OA   =  EA(i)/PSUR(i)
        ECMOLE = 55.56 * ECMASS(i) * dtti  ! ecmass must be computed and passed in
        H2OS = H2OA + ECMOLE / GBH2O(i)
        H2OS  = MIN( H2OS, H2OI )
        H2OS  = MAX( H2OS, 1. EEEXP -7)
        H2OSRH = H2OS/H2OI
C  need qdamp and pdamp calculated and passed to here!
cpl        CO2S(i) = MAX(PCO2S(I),PCO2M*0.5) / (PSUR(i)*100.)

cpl I have relaxed this condition to 1/10 of previous. The old way made
cpl the CO2 on top of the leaves always at least 1/2 of the value at the
cpl reference level.

        CO2S(i) = MAX(PCO2S(I),PCO2M(i)*0.05) / (PSUR(i)*100.)
        
cpl Ball-Berry equation right here !        

        GSH2OINF = (GRADM(i) * MAX(1. EEEXP -12,ASSIMN(i))
     >            * H2OSRH * soilfrz(i) / CO2S(i)) + BINTC(i)
     
cpl this is the change in stomatal resistance

        DRST(i) = RST(i) * QDAMP * ((GSH2O(i)-GSH2OINF)/
     &            (PDAMP*GSH2O(i)+QDAMP*GSH2OINF))


cpl this is the 'would be change' if we did not use the damping factor.

c        rstnew = (1./gsh2oinf) * 44.6 * tprcor(i)/tc(i)
c        DRST(i) = rstnew - RST(i)

c
        RSTFAC(i,1) = H2OS/H2OI
        RSTFAC(i,4) = RSTFAC(i,1)*RSTFAC(i,2)* RSTFAC(i,3)
      enddo
C
CZ CARNEGIE new diagnostics----start!!!(c.zhang&joe berry, 10/19/92)
c-----------------------------------------------------------------------
C  INPUTS: PSUR(i),CO2S,ASSIMN(i),GRADM(i),BINTC(i),VMAX0(i),RRKK(i),C3(i),
C    C4(i),PAR(i),ATHETA(i),BTHETA(i),APARKK(i),OMSS(i),RSTFAC(i,2),TEMPH,
C    TEMPL,RSTFAC(i,1),VM(i),ASSIM,GSH20(i),EFFCON(i),QT,GAMMAS(i),
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
      OMCPOT = VMAX0(i)*2.1**qt(i)*((PCO2IPOT-GAMMAS(i))/
     &       (PCO2IPOT+RRKK(i))*C3(i) + C4(i))
C
C Potential light limitation.
      OMEPOT(i) = PAR(i)*((PCO2IPOT-GAMMAS(i))/
     &           (PCO2IPOT+2.*GAMMAS(i))*C3(i) + C4(i))
C
C Quad 1.
      SQRTIN = MAX(0. EEEXP 0,((OMEPOT(i)+OMCPOT)**2-
     &               4.*ATHETA(i)*OMEPOT(i)*OMCPOT))
C
C Quad 1. Intermediate  top leaf photosynthesis.
      OMPPOT = ((OMEPOT(i)+OMCPOT)-SQRT(SQRTIN))/(2.*ATHETA(i))
C
C Potential sink or pep limitation.
      OMSPOT = (VMAX0(i)/2.0)*(1.8**qt(i))*C3(i) 
     &                + RRKK(i)*PCO2IPOT*C4(i)
C
C Quad 2.
      SQRTIN=MAX(0. EEEXP 0,((OMPPOT+OMSPOT)**2-4.*BTHETA(i)*
     &        OMPPOT*OMSPOT))
C
C Quad 2. Final Potential top leaf photosynthesis.
      ASSIMPOT(i) = ((OMSPOT+OMPPOT)-SQRT(SQRTIN))/(2.*BTHETA(i))
C
C-----------------------------------------------------------------------
C CALCULATION OF STRESS FACTOR LIMITED ASSIMILATION
C-----------------------------------------------------------------------
C
C Stressed rubisco limitation.
      OMCCI = VM(i)*((PCO2IPOT-GAMMAS(i))/(PCO2IPOT+RRKK(i))*C3(i) 
     &             + C4(i))
C
C Quad 1.
      SQRTIN = MAX(0. EEEXP 0,(OMEPOT(i)+OMCCI)**2-
     &                4.*ATHETA(i)*OMEPOT(i)*OMCCI)
C
C Quad 1. Intermediate stress limited top leaf photosynthesis.
      OMPCI = ((OMEPOT(i)+OMCCI)-SQRT(SQRTIN))/(2.*ATHETA(i))
C
C Stressed sink or pep limitation.
      OMSCI = OMSS(i)*(C3(i) + PCO2IPOT*C4(i))
C
C Quad 2.
      SQRTIN = MAX(0. EEEXP 0,(OMPCI+OMSCI)**2-4.*BTHETA(i)*OMPCI*OMSCI)
C 
C Quad 2. Final stress limited top leaf photosynthesis.
      ASSIMCI(i) = ((OMSCI+OMPCI)-SQRT(SQRTIN))/(2.*BTHETA(i))
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
C
C-----------------------------------------------------------------------
C OUTPUT:  WEIGHTED STRESS FACTORS AND OTHER DIAGNOSTIC OUTPUTS TO BE SUMMED
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

      RETURN
      END

