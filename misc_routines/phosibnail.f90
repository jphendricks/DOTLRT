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
     >     pco2s,co2cap,cflux,c4fract, Physnum, Phystype, PHYSMAX)

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
C      RESPC               = CANOPY RESPIRATION (NPH: per unit leaf area)
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
c
c NPH Added array structure and code for multiple physiological types
c NPH New subroutine arguments Physnum, Phystype and PHYSMAX. All physiol
c NPH variables now have extra array dimension (PHYSMAX). Output and calc
c NPH variables have array dimension (PHYSMAX+1) for separate results and
c NPH weighted total.
c NPH C3/C4 code based on stand alone of 9/1999. Finalized in SiB2.5 May 2001

C
c     input arrays:
      integer len, nsib, nsoil, PHYSMAX
      logical forcerestore

      REALX8 vmax0(len,PHYSMAX),psur(len),green(len),gmudmu(len),
     *    zlt(len),cas_cap_co2(len),tc(len),ta(len),
     *    trop(len,PHYSMAX),trda(len,PHYSMAX),
     *    slti(len,PHYSMAX),shti(len,PHYSMAX),
     *    hltii(len,PHYSMAX),hhti(len,PHYSMAX),
     *    etgs(len),wc(len),em(len),ra(len),rb(len),
     *    cog1(len),cog2(len),tm(len),effcon(len,PHYSMAX),
     *    binter(len,PHYSMAX),gradm(len,PHYSMAX),
     *    atheta(len,PHYSMAX),btheta(len,PHYSMAX),
     *    tgs(len),respcp(len,PHYSMAX),tran(nsib,2,2),ref(nsib,2,2),
     *    radn(len,2,2),ecmass(len),trdm(len,PHYSMAX),etc(len),
     *    aparc(len),www(len,2),tg(len),rst(len,PHYSMAX+1), cflux(len),
     *    c4fract(len)
      REALX8 pdamp, qdamp, dtt, pco2m(len), pco2ap(len), tf, po2m
      integer Physnum(len), Phystype(len,PHYSMAX), ilen
     
c     output arrays:


      REALX8 assimn(len,PHYSMAX+1),ea(len),rstfac(len,4),
     *     pco2i(len,PHYSMAX+1),respc(len,PHYSMAX+1),
     *     respg(len),drst(len,PHYSMAX+1)
czz new diagostics 10/14/92 
c
c output arrays
      REALX8 omepot(len,PHYSMAX+1),assimpot(len,PHYSMAX+1),
     >     assimci(len,PHYSMAX+1),
     >     assimnp(len,PHYSMAX+1),whs(len,PHYSMAX+1),
     >     antemp(len,PHYSMAX+1),
     >     wsfws(len,PHYSMAX+1),wsfht(len,PHYSMAX+1),
     >     wsflt(len,PHYSMAX+1),wci(len,PHYSMAX+1),
     >     wags(len,PHYSMAX+1),wegs(len,PHYSMAX+1),pfd(len),
     >     td(nsib,nsoil),zmlscale(len),assim(len,PHYSMAX+1),
     >     wopt(len),zm(len),wsat(len),
     >     soilscale(len,nsoil+1),zmstscale(len,2),zltrscale(len),
     >     tprcor(len),bintc(len,PHYSMAX+1),soilq10(len,nsoil+1),
     >     ansqr(len,PHYSMAX+1),soilscaleold(len)
     *     ,pco2c(len,PHYSMAX+1) !chloroplast pco2
     *     ,xgah2o(len) 
     *     ,xgco2m(len)
     
c     work arrays:

      REALX8 PCO2Y(len,6), EYY(len,6),assimny(len,6),
     *     assimy(len,6)
      REALX8 c3(len),c4(len),range(len),gammas(len),
     *     aparkk(len),gah2o(len),
     *     gbh2o(len),
     *     par(len),rrkk(len),
     *     omss(len),vm(len),gsh2o(len),pco2s(len,PHYSMAX+1),
     *     templ(len),temph(len),
     *     qt(len),co2s(len),scatp(len),scatg(len),
     *     park(len),respn(len),zkc(len),resa(len),
     *     zko(len),spfy(len), co2a(len), co2m(len),co2cap(len),
c NPH Add co2ap - working array for CAS [CO2] and separation of "green"
c NPH between physiological types (C3 and C4) 
     *     Pgreen(PHYSMAX), CO2Ap(len)

      integer icconv(len),igath(len)

      REALX8 soilfrz(len)
      REALX8 bl, tresp, soilq10td, b(2), woptzm, cwsflt, cwsfht, cwsfws,
     *     ccoms, ccomc, ascitemp, dompdomc, omsci, ompci, omcci,
     *     omcpot, omppot, sqrtin, omspot, pco2ipot, ohtp2, sttp2,
     *     gsh2oinf, h2osrh, h2os, ecmole, h2oa, h2oi, dtti,
     *     pco2in, pco2a, soilfrztd, soilfrztg, soilq10tg, 
     *     rstnew

      integer i, ic1, ic, l, Pindx

      
      
cpl introduce a co2 capacity 
cpl this will basically be the mass of air under the top of the canopy (in
cpl this case (CHEAS-RAMS) O(10-30m), that is, ground to displacemnt height.

cpl all the carbon fluxes are expresse as Mol C / m2 s and resistances for
cpl carbon are in m2 s / mol air

cpl one mole of gas occupies 22.4 cubic dm
cpl 1 cubic meter contains therefore 1000./22.4  = 44.6 moles of gas
cpl the units for the carbon capacity are mol air /m2. cpl (e.g. here 893 moles if thickness of the layer is 20m) cpl this means that the units for pc02ap should be mol co2 / mol air, but
cpl it is also possible to keep just co2 pressure and convert


c NPH Change to have one main gridpoint loop for better flow
c NPH Numerous old loops now commented out
c NPH Could be reinstated later if needed for e.g. parallel code


c START OF PRIMARY PHOSIB LOOP OVER GRIDPOINTS
      do i = 1,len                !  LOOP OVER GRIDPOINT

       TPRCOR(i) = TF*PSUR(i)*100./1.013 E5
       co2cap(i) = cas_cap_co2(i) * 44.6 * tprcor(i)/ta(i)     ! moles air / m2
       ilen = i

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


c NPH START PHYSIOLOGY LOOP
       do Pindx = 1, Physnum(i)   !  LOOP OVER VEG/PHYSIOL type

	if(Phystype(i,Pindx) .eq. 4) then
	 Pgreen(Pindx) = c4fract(i)
	 C3(i) = 0.
        else
	 Pgreen(Pindx) = (1. - c4fract(i))
	 C3(i) = 1.
        endif
	C4(i) = 1. - C3(i)

        APARKK(i) = (APARC(i)/PARK(i)) * (GREEN(i)*Pgreen(Pindx))

C-----------------------------------------------------------------------
C
C     Q-10 AND STRESS TEMPERATURE EFFECTS
C
C      QT          (QT)    : TABLE (2)     , SE-92A
C-----------------------------------------------------------------------
C
        qt(i) = 0.1*( TC(i) - TROP(i,Pindx) )
        RESPN(i) = RESPCP(i,Pindx) * VMAX0(i,Pindx) * RSTFAC(i,2)

citb...patch to prevent underflow if temp is too cool...
        if(TC(i) >= TRDM(i,Pindx))then
            RESPC(i,Pindx) = RESPN(i) * 2.0**qt(i)
     *        /( 1. + EXP( TRDA(i,Pindx)*(TC(i)-TRDM(i,Pindx))))
        else
            RESPC(i,Pindx) = RESPN(i) * 2.0**qt(i)
        endif

        VM(i) = VMAX0(i,Pindx) * 2.1**qt(i)
        TEMPL(i) = 1. + EXP(SLTI(i,Pindx)*(HLTIi(i,Pindx)-TC(i)))
        TEMPH(i) = 1. + EXP(SHTI(i,Pindx)*(TC(i)-HHTI(i,Pindx)))
        RSTFAC(i,3) = 1./( TEMPL(i)*TEMPH(i))
c NPH Note that separate rstfac3 not saved for different physiological
c NPH types, although calcs within phosib are separate. If the rstfac 
c NPH diagnostic output needs to be either separate or weighted average
c NPH we'd need to add code to that effect
c NPH Also note that rstfac1 can't easily be separated. This currently
c NPH depends on CAS humidity and leaf temperatures which are similar
c NPH between physiological types in this implementation (since energy bal
c NPH is not separated between phys types).
c NPH Futher note that RSTFAC3 for C3 uses only temph. Calculated/output 
c NPH RSTFAC is used for C4 P/S only.

        VM(i)    = VM(i)/TEMPH(i) * RSTFAC(i,2)*C3(i)
     &      + VM(i) * RSTFAC(i,2)*RSTFAC(i,3) * C4(i)

c	print*,Pindx, Vmax0(i,Pindx),tc(i),templ(i),temph(i)
c	print*,RSTFAC(i,2),www(i,1),www(i,2),www(i,3),VM(i)
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
        GSH2O(i)  = 1.0/RST(i,Pindx) * 44.6*TPRCOR(i)/TC(i)
        GBH2O(i)  = 0.5/RB(i) * 44.6*TPRCOR(i)/TC(i)
        GAH2O(i)  = 1.0/RA(i) * 44.6*TPRCOR(i)/TM(i)

        xgah2o(i) = max(0.466 EEEXP 0, gah2o(i) )
        xgco2m(i) = 4000.0 * vmax0(i,Pindx)
c NPH ISOTOPE NOTE - xgco2m may need to be redimensioned for phys type...
C
        RRKK(i)   = ZKC(i)*( 1. + PO2M/ZKO(i) ) * C3(i)
     &               + VMAX0(i,Pindx)/5.* ( 1.8**qt(i)) * C4(i)
        PAR(i)    = pfd(i)*EFFCON(i,Pindx)*( 1.-SCATG(i) )
        soilfrztg = 1.+exp(-1.5 * (max(270.0 EEEXP 0,tgs(i))-273.16) )
        soilfrztd = 1.+exp(-1.5 * (max(270.0 EEEXP 0,td (i,nsoil))-273.16) )
        soilfrz(i) = max(1./soilfrztg, 1./soilfrztd)
        soilfrz(i) = max( soilfrz(i), 0.05 EEEXP 0)
c NPH BINTC is in mol/m2/s (same units as BINTER)
        BINTC(i,Pindx)  = BINTER(i,Pindx)*ZLT(i)*GREEN(i)*Pgreen(Pindx)*
     &                 RSTFAC(i,2) * soilfrz(i)
        OMSS(i)   = ( VMAX0(i,Pindx)/2.0 ) * ( 1.8**qt(i) )
     &                  /TEMPL(i) * RSTFAC(i,2) * C3(i)
     &                  + RRKK(i) * RSTFAC(i,2) * C4(i)
C
C-----------------------------------------------------------------------
C
C     FIRST GUESS IS MIDWAY BETWEEN COMPENSATION POINT AND MAXIMUM
C     ASSIMILATION RATE.
C
C-----------------------------------------------------------------------
       
        RANGE(i)    = PCO2M(i) * ( 1. - 1.6/GRADM(i,Pindx) ) - GAMMAS(i)
	    icconv(i) = 1
	

      DO 1000 IC = 1, 6
          PCO2Y(i,IC) = 0.
          EYY(i,IC) = 0.
	  assimny(i,ic) = 0.
	  assimy(i,ic) = 0.
1000  CONTINUE


c NPH Move resa and [CO2] conversions out of SORTIN loop to avoid
c NPH repeat calcs and inverting resa each time through loop!
      resa(i) =   1. / MAX(0.446 EEEXP 0,GAH2O(i))
c       resb(i) =   1.4/GBH2O(i)
c       resc(i) =   1.6/GSH2O(i)
      CO2Ap(i) = PCO2Ap(i) /   (PSUR(i)*100.)
      co2m(i)  = pco2m(i)  /   (PSUR(i)*100.)


      DO 2000 IC = 1, 6
C
        CALL       SORTIN( EYY, PCO2Y, RANGE, GAMMAS, ic,len,ilen )
        
        CALL       CYCALC( APARKK, VM, ATHETA(1,Pindx), 
     &	                   BTHETA(1,Pindx), par, GAMMAS,
     &                     RESPC(1,Pindx), RRKK, OMSS, C3, C4,
     &                     PCO2Y(1,ic), assimny(1,ic), assimy(1,ic),
     &                     len, ilen )

c	print*,'1',APARKK(i),VM(i),atheta(i,Pindx),btheta(i,Pindx)
c	print*,'2',par(i),GAMMAS(i),respc(i,Pindx),rrkk(i),c3(i),c4(i)
c	print*,'3',pco2y(i,ic),assimny(i,ic),assimy(i,ic)

cpl now prognose the new CAS CO2 according to flux divergence
cpl we are going to do this in mol C / mol air (same as PaC/PaAir)

c NPH For mixed physiology, use previous time-step CAS pCO2. This saves 
c NPH iteration to solve for two assim fluxes. Is also a physically and 
c NPH physiologically reasonable solution to the problem.

        PCO2S(i,Pindx) = PCO2Ap(i) - (1.4/GBH2O(i) * ASSIMNy(i,ic)
     >            * PSUR(i)*100.)
        PCO2IN   = PCO2S(i,Pindx) - (1.6/GSH2O(i) * ASSIMNy(i,ic)
     >            * PSUR(i)*100.)
        EYY(i,IC) = PCO2Y(i,IC) - PCO2IN

        if(ic.ge.2) then
	ic1 = ic-1
            if(abs(eyy(i,ic1)).ge.0.1)then
	      icconv(i) = ic
	    else
	      eyy(i,ic) = eyy(i,ic1)
	      pco2y(i,ic) = pco2y(i,ic1)
	    endif
	endif

2000  CONTINUE

        icconv(i) = min(icconv(i),6)
        igath(i) = i+(icconv(i)-1)*len
        dtti = 1./dtt

c NPH Get final (converged) values of assimilation, CO2, etc.
	pco2i(i,Pindx) = pco2y(igath(i),1)
	assimn(i,Pindx) = assimny(igath(i),1)
	assim(i,Pindx) = assimy(igath(i),1)
        PCO2S(i,Pindx) = PCO2Ap(i) - (1.4/GBH2O(i) * ASSIMN(i,Pindx)
     >     *  PSUR(i)*100.)
        pco2c(i,Pindx) = pco2i(i,Pindx) - assimn(i,Pindx)/
     >      xgco2m(i)*psur(i)*100.0
c NPH ISOTOPE NOTE - PCO2S redimensioned for phys type..
c NPH ISOTOPE NOTE - pco2c redimensioned for phys type..

cczzggrst5 - new code
        H2OI   = ETC(i)/PSUR(i)
        H2OA   =  EA(i)/PSUR(i)
        ECMOLE = 55.56 * ECMASS(i) * dtti  ! ecmass must be computed and passed in
        H2OS = H2OA + ECMOLE / GBH2O(i)
        H2OS  = MIN( H2OS, H2OI )
        H2OS  = MAX( H2OS, 1. EEEXP -7)
        H2OSRH = H2OS/H2OI

        CO2S(i) = MAX(PCO2S(I,Pindx),PCO2M(i)*0.05) / (PSUR(i)*100.)
        
cpl Ball-Berry equation right here !        

        GSH2OINF = (GRADM(i,Pindx) * MAX(1. EEEXP -12,ASSIMN(i,Pindx))
     >            * H2OSRH * soilfrz(i) / CO2S(i)) + BINTC(i,Pindx)
     
cpl this is the change in stomatal resistance

c        DRST(i,Pindx) = RST(i,Pindx) * QDAMP * ((GSH2O(i)-GSH2OINF)/
c     &            (PDAMP*GSH2O(i)+QDAMP*GSH2OINF))
c NPH Do the "damping" factors convert to m/s or is DRST still in mol??
c NPH Given the form of PDAMP/QDAMP there does not appear to be any 
c NPH conversion so DRTS should still be is m.s/mol. Is it converted
c NPH somewhere else before adding to RST ???!!
c NPH In this new version the RST updates are done in phosib.
c NPH In case the above is wrong I'm using the direct method. I think
c NPH it makes more sense anyway (the damping is overkill!)

cpl this is the 'would be change' if we did not use the damping factor.

        rstnew = (1./gsh2oinf) * (44.6 * tprcor(i)/TC(i))
        DRST(i,Pindx) = rstnew - RST(i,Pindx)
c NPH RST now updated here for each physiological type and at end of
c NPH phosib for the canopy-weighted total. Units converted back to m/s
	RST(i,Pindx) = rstnew

c
        RSTFAC(i,1) = H2OS/H2OI
        RSTFAC(i,4) = RSTFAC(i,1)*RSTFAC(i,2)* RSTFAC(i,3)

c        print*,'Pindx,An,pi,pap,fC4,Pgreen,green,pfd,PAR',
c     &   Pindx,assimn(i,Pindx),pco2i(i,Pindx),pco2ap(i),
c     &   c4fract(i),Pgreen(Pindx),green(i),pfd(i),PAR(i)


CZ CARNEGIE new diagnostics----start!!!(c.zhang&joe berry, 10/19/92)
c-----------------------------------------------------------------------
C  INPUTS: PSUR(i),CO2S,ASSIMN(i),GRADM(i),BINTC(i),VMAX0(i),RRKK(i),C3(i),
C    C4(i),PAR(i),ATHETA(i),BTHETA(i),APARKK(i),OMSS(i),RSTFAC(i,2),TEMPH,
C    TEMPL,RSTFAC(i,1),VM(i),ASSIM,GSH20(i),EFFCON(i),QT,GAMMAS(i),
C    PFD(i)
C

        sttp2 = 73.**0.2
        ohtp2 = 100.**0.2
c NPH        do i = 1,len
C-----------------------------------------------------------------------
C CALCULATION OF POTENTIAL ASSIMILATION
C-----------------------------------------------------------------------
C
C Make assimn a top leaf, not the canopy.
      ASSIMNp(i,Pindx) = ASSIMN(i,Pindx) / APARKK(i)
C
C Bottom stopped assim.
      ANTEMP(i,Pindx) = MAX(0. EEEXP 0,ASSIMNp(i,Pindx))
C
C Potential intercellular co2.
      PCO2IPOT = PSUR(i)*100.*(co2s(i)-(1.6*ASSIMNp(i,Pindx)/
     &  ((GRADM(i,Pindx)*ANTEMP(i,Pindx)/co2s(i))+BINTC(i,Pindx))))
C
C Potential rubisco limitation.
      OMCPOT = VMAX0(i,Pindx)*2.1**qt(i)*((PCO2IPOT-GAMMAS(i))/
     &       (PCO2IPOT+RRKK(i))*C3(i) + C4(i))
C
C Potential light limitation.
      OMEPOT(i,Pindx) = PAR(i)*((PCO2IPOT-GAMMAS(i))/
     &       (PCO2IPOT+2.*GAMMAS(i))*C3(i) + C4(i))
C
C Quad 1.
      SQRTIN = MAX(0. EEEXP 0,((OMEPOT(i,Pindx)+OMCPOT)**2-
     &       4.*ATHETA(i,Pindx)*OMEPOT(i,Pindx)*OMCPOT))
C
C Quad 1. Intermediate  top leaf photosynthesis.
      OMPPOT = ((OMEPOT(i,Pindx)+OMCPOT)-SQRT(SQRTIN))/
     & (2.*ATHETA(i,Pindx))
C
C Potential sink or pep limitation.
      OMSPOT = (VMAX0(i,Pindx)/2.0)*(1.8**qt(i))*C3(i) 
     &       + RRKK(i)*PCO2IPOT*C4(i)
C
C Quad 2.
      SQRTIN=MAX(0. EEEXP 0,((OMPPOT+OMSPOT)**2-4.*BTHETA(i,Pindx)*
     &        OMPPOT*OMSPOT))
C
C Quad 2. Final Potential top leaf photosynthesis.
      ASSIMPOT(i,Pindx) = ((OMSPOT+OMPPOT)-SQRT(SQRTIN))/
     &        (2.*BTHETA(i,Pindx))
C
C-----------------------------------------------------------------------
C CALCULATION OF STRESS FACTOR LIMITED ASSIMILATION
C-----------------------------------------------------------------------
C
C Stressed rubisco limitation.
      OMCCI = VM(i)*((PCO2IPOT-GAMMAS(i))/(PCO2IPOT+RRKK(i))*C3(i) 
     &        + C4(i))
C
C Quad 1.
      SQRTIN = MAX(0. EEEXP 0,(OMEPOT(i,Pindx)+OMCCI)**2-
     &        4.*ATHETA(i,Pindx)*OMEPOT(i,Pindx)*OMCCI)
C
C Quad 1. Intermediate stress limited top leaf photosynthesis.
      OMPCI = ((OMEPOT(i,Pindx)+OMCCI)-SQRT(SQRTIN))/
     &        (2.*ATHETA(i,Pindx))
C
C Stressed sink or pep limitation.
      OMSCI = OMSS(i)*(C3(i) + PCO2IPOT*C4(i))
C
C Quad 2.
      SQRTIN = MAX(0. EEEXP 0,(OMPCI+OMSCI)**2-4.
     &        *BTHETA(i,Pindx)*OMPCI*OMSCI)
C 
C Quad 2. Final stress limited top leaf photosynthesis.
      ASSIMCI(i,Pindx) = ((OMSCI+OMPCI)-SQRT(SQRTIN))/
     &        (2.*BTHETA(i,Pindx))
C
C-----------------------------------------------------------------------
C CALCULATION OF CONTROL COEFFICIENTS
C-----------------------------------------------------------------------
C
C Intermediate.
      DOMPDOMC = (OMPCI-OMEPOT(i,Pindx))/
     &       (2.*ATHETA(i,Pindx)*OMPCI-OMCCI-OMEPOT(i,Pindx))
C
C Bottom stopped final stress limited top leaf photosynthesis.
      ASCITEMP = MAX(ASSIMCI(i,Pindx),1. EEEXP -12)
C
C Rubisco control coefficient.
      CCOMC = (DOMPDOMC*(ASSIMCI(i,Pindx)-OMSCI)/(2.*BTHETA(i,Pindx)
     &      *ASSIMCI(i,Pindx)-OMPCI-OMSCI))*OMCCI/ASCITEMP
C
C Sink or pep control coefficient.
      CCOMS = ((ASSIMCI(i,Pindx)-OMPCI)/ (2.*BTHETA(i,Pindx)
     &       *ASSIMCI(i,Pindx)-OMPCI-OMSCI))*OMSCI/ASCITEMP
C
C-----------------------------------------------------------------------
C  OUTPUT:  POTENTIAL ASSIMILATION RATES TO BE SUMMED
C-----------------------------------------------------------------------
C Canopy values (overwrites top leaf).
C 
      OMEPOT(i,Pindx) = OMEPOT(i,Pindx)*APARKK(i)
      ASSIMPOT(i,Pindx) = ASSIMPOT(i,Pindx)*APARKK(i)
      ASSIMCI(i,Pindx) = ASSIMCI(i,Pindx)*APARKK(i)
      ASSIM(i,Pindx) = ASSIM(i,Pindx)*APARKK(i)
      ANTEMP(i,Pindx) = ANTEMP(i,Pindx)*APARKK(i)
      ANSQR(i,Pindx) = ANTEMP(i,Pindx)*ANTEMP(i,Pindx)
      ASSIMNp(i,Pindx) = ASSIMNp(i,Pindx)*APARKK(i)
C
C-----------------------------------------------------------------------
C OUTPUT:  WEIGHTED STRESS FACTORS AND OTHER DIAGNOSTIC OUTPUTS TO BE SUMMED
C-----------------------------------------------------------------------
C
C Water stress.
      WSFWS(i,Pindx) = ASSIMPOT(i,Pindx)*(1.-RSTFAC(i,2))*(CCOMC+CCOMS)
C
C High temperature stress.
      WSFHT(i,Pindx) = ASSIMPOT(i,Pindx)*(1.-1./TEMPH(i))*CCOMC
C
C Low temperature stress.
      WSFLT(i,Pindx) = ASSIMPOT(i,Pindx)*(1.-1./TEMPL(i))*(CCOMS*
     &      C3(i)+CCOMC*C4(i))
c
c  protection for wsfws, wsfht, and wsflt from <0 or >>xxx(2/24/93)
      cwsfws = (1.-RSTFAC(i,2))*(CCOMC+CCOMS)
      if(cwsfws.gt.1. .or. cwsfws.lt.0.) wsfws(i,Pindx)=0.
      cwsfht = (1.-1./TEMPH(i))*CCOMC
      if(cwsfht.gt.1. .or. cwsfht.lt.0.) wsfht(i,Pindx)=0.
      cwsflt = (1.-1./TEMPL(i))*(CCOMS*C3(i)+CCOMC*C4(i))
      if(cwsflt.gt.1. .or. cwsflt.lt.0.) wsflt(i,Pindx)=0.

C
C Intermediate assimilation weighted Ci.
      WCI(i,Pindx) = ANTEMP(i,Pindx)*PCO2I(i,Pindx)
C
C Intermediate assimilation weighted relative humidty stress factor.
      WHS(i,Pindx) = ANTEMP(i,Pindx)*RSTFAC(i,1)
C
C Intermediate assimilation weighted stomatal conductance.
      WAGS(i,Pindx) = GSH2O(i)*ANTEMP(i,Pindx)
C
C Intermediate evaporation weighted stomatal conductance.(Step 1.
C   Step 2 after subroutine updat2)
      WEGS(i,Pindx) = GSH2O(i)

       enddo		! END PHYSIOLOGICAL-TYPE LOOP


c NPH  Sum output arrays
c NPH  NOTE: Leaf-level values calculated as weighted
c NPH  average of individual values, whereas canopy values are simple
c NPH  sums of the C3 and C4 components. All canopy values already
c NPH  calculated with PIfactor*green(i)*Pgreen(Pindx)
c
c Zero last element of output arrays
        assim(i,Physnum(i)+1)=0.
        assimn(i,Physnum(i)+1)=0.
        bintc(i,Physnum(i)+1)=0.
        respc(i,Physnum(i)+1)=0.
        RST(i,Physnum(i)+1)= RST(i,1)
	pco2i(i,Physnum(i)+1)=0.
	pco2c(i,Physnum(i)+1)=0.
	pco2s(i,Physnum(i)+1)=0.

	omepot(i,Physnum(i)+1)=0.
	assimpot(i,Physnum(i)+1)=0.
	assimci(i,Physnum(i)+1)=0.
	antemp(i,Physnum(i)+1)=0.
	assimnp(i,Physnum(i)+1)=0.
	wsfws(i,Physnum(i)+1)=0.
	wsfht(i,Physnum(i)+1)=0.
	wsflt(i,Physnum(i)+1)=0.
	wci(i,Physnum(i)+1)=0.
	whs(i,Physnum(i)+1)=0.
	wags(i,Physnum(i)+1)=0.
	wegs(i,Physnum(i)+1)=0.
	ansqr(i,Physnum(i)+1)=0.

C-----------------------------------------------------------------------
c NPH Loop through physnum(i) to calculate weighted total fluxes
C-----------------------------------------------------------------------
       do Pindx = 1, Physnum(i)         !  LOOP OVER VEG/PHYSIOL type

c Main phosib results
c Use total green fraction as weighting denominator (leaf-level)
	 assim(i,Physnum(i)+1)=assim(i,Physnum(i)+1) +
     &      assim(i,Pindx)*Pgreen(Pindx)/GREEN(i)
         assimn(i,Physnum(i)+1)=assimn(i,Physnum(i)+1) +
     &      assimn(i,Pindx)
         bintc(i,Physnum(i)+1)=bintc(i,Physnum(i)+1) +
     &      bintc(i,Pindx)
	 respc(i,Physnum(i)+1)=respc(i,Physnum(i)+1) +
     &      respc(i,Pindx)*Pgreen(Pindx)/GREEN(i)
         pco2i(i,Physnum(i)+1)=pco2i(i,Physnum(i)+1) +
     &       pco2i(i,Pindx)*Pgreen(Pindx)/GREEN(i)
         pco2s(i,Physnum(i)+1)=pco2s(i,Physnum(i)+1) +
     &       pco2s(i,Pindx)*Pgreen(Pindx)/GREEN(i)
         pco2c(i,Physnum(i)+1)=pco2c(i,Physnum(i)+1) +
     &       pco2c(i,Pindx)*Pgreen(Pindx)/GREEN(i)
c Carnegie diagnostics
c (May need to check that these still logical/correctly weighted)
         omepot(i,Physnum(i)+1)=omepot(i,Physnum(i)+1) +
     &      omepot(i,Pindx)
	 assimpot(i,Physnum(i)+1)=assimpot(i,Physnum(i)+1) +
     &      assimpot(i,Pindx)
	 assimci(i,Physnum(i)+1)=assimci(i,Physnum(i)+1) +
     &      assimci(i,Pindx)
         antemp(i,Physnum(i)+1)=antemp(i,Physnum(i)+1) +
     &      antemp(i,Pindx)
         assimnp(i,Physnum(i)+1)=assimnp(i,Physnum(i)+1) +
     &      assimnp(i,Pindx)
         wsfws(i,Physnum(i)+1)=wsfws(i,Physnum(i)+1) +
     &      wsfws(i,Pindx)*Pgreen(Pindx)/GREEN(i)
	 wsfht(i,Physnum(i)+1)=wsfht(i,Physnum(i)+1) +
     &      wsfht(i,Pindx)*Pgreen(Pindx)/GREEN(i)
         wsflt(i,Physnum(i)+1)=wsflt(i,Physnum(i)+1) +
     &       wsflt(i,Pindx)*Pgreen(Pindx)/GREEN(i)
	 wci(i,Physnum(i)+1)=wci(i,Physnum(i)+1) +
     &       wci(i,Pindx)*Pgreen(Pindx)/GREEN(i)
         whs(i,Physnum(i)+1)=whs(i,Physnum(i)+1) +
     &       whs(i,Pindx)*Pgreen(Pindx)/GREEN(i)
         wags(i,Physnum(i)+1)=wags(i,Physnum(i)+1) +
     &       wags(i,Pindx)*Pgreen(Pindx)/GREEN(i)
         wegs(i,Physnum(i)+1)=wegs(i,Physnum(i)+1) +
     &       wegs(i,Pindx)*Pgreen(Pindx)/GREEN(i)
         ansqr(i,Physnum(i)+1)=ansqr(i,Physnum(i)+1) +
     &       ansqr(i,Pindx)
       enddo                            ! End Physiol loop


c Sum resistance
c NPH Note that Rst now updated here (no longer in addinc.F)
	if (Physnum(i) .ge. 2) then
          do Pindx = 2, Physnum(i)         !  LOOP OVER VEG/PHYSIOL type
            RST(i,Physnum(i)+1)= 1.0 / ( (1.0/RST(i,Physnum(i)+1)) +
     &         (1.0/RST(i,Pindx)))
          enddo
        endif

c NPH Moved following RST limits from addinc
c bintc(i)- smallest canopy stomatal conductance needs to be passed in here.
c ---- c.zhang, 2/3/93
	 rst(i,Physnum(i)+1) = MIN( 1./(bintc(i,Physnum(1)+1)*tc(i)/
     &        (44.6 * tprcor(i))), rst(i,Physnum(i)+1) )
	 rst(i,Physnum(i)+1) = MAX( 10., rst(i,Physnum(i)+1) )

c NPH FINALLY use total assimn and resa to calculate new CO2A

	CO2A(i) = ( CO2Ap(i) + (dtt/co2cap(i))
     &     *  ( respg(i) - assimn(i,Physnum(i)+1)
     &     +  co2m(i)*gah2o(i) ) )
     &     /  (1+dtt*gah2o(i)/co2cap(i) ) 


c NPH Modified CO2A based on previous value and current fluxes
c NPH The following does not work! Why not?!
c	CO2A(i)   = CO2Ap(i) +                          ! previous CAS [CO2]
c     &     (dtt/co2cap(i)) *                           ! time/CAS
c     &     ( respg(i) - assimn(i,Physnum(i)+1) +       ! lower fluxes
c     &     (co2m(i)-CO2Ap(i))/(resa(i)) )    ! upper flux

c NPH Calculate PCO2A to carry to next time step (PCO2Ap)
       PCO2Ap(i) = co2a(i) * psur(i) * 100.

c NPH Trap extreme values of PCO2Ap assuming no CAS capacitance
       if(PCO2Ap(i) .ge. 300.) then
 	 PCO2Ap(i) = PCO2m(i) + ((respg(i)-ASSIMN(i,Physnum(i)+1))*
     >        resa(i) ) * PSUR(i)*100.
         print*,'correcting PCO2Ap down'
         CO2A(i) = pco2ap(i) / (psur(i) * 100.)
       endif
       if(PCO2Ap(i) .le. 0.) then
 	 PCO2Ap(i) = PCO2m(i) + ((respg(i)-ASSIMN(i,Physnum(i)+1))*
     >        resa(i) ) * PSUR(i)*100.
         print*,'correcting PCO2Ap up'
         CO2A(i) = pco2ap(i) / (psur(i) * 100.)
       endif


citb...carbon flux between CAS and reference level (in kg m-2 s-1)
	 cflux(i) = gah2o(i)*(co2a(i)-co2m(i))*0.012


c       print*,'An1,An2,Atot,Rst1,Rst2,Rstot,NBE', assimn(i,1),
c     &    assimn(i,2),assimn(i,3),RST(i,1),RST(i,2),RST(i,3),
c     &    (respg(i)-assimn(i,Physnum(i)+1))

c       print*,'gah2o,ca,cm,cflux', gah2o(i),co2a(i),
c     &    co2m(i),cflux/0.012
c       print*

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
c NPH      enddo

c-----------------------------------------------------------------------
      enddo             ! END PRIMARY LOOP OVER GRIDPOINTS

      RETURN
      END

