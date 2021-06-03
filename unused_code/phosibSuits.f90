!
!==================SUBROUTINE PHOSIB===================================
!
      SUBROUTINE PHOSIB_suits(pco2m,pco2ap,po2m,vmax0,tf,psur,green &
      ,          tran,ref,gmudmu,zlt,cas_cap_co2,tc,ta,trop,trda &
      ,          trdm,slti,shti,hltii,hhti,radn,etc,etgs,wc &
      ,          ea,em,rb,ra,tm &
      ,          effcon,rstfac,binter,gradm,assimn &
      ,          rst,atheta,btheta,tgs,respcp &
      ,          aparkk,len,nsib, &
          omepot,assimpot,assimci,antemp,assimnp, &
          wsfws,wsfht,wsflt,wci,whs, &
          wags,wegs,aparc,pfd,assim,td,www, &
          wopt,zm,wsat,tg,soilscale,zmstscale,zltrscale,zmlscale, &
          drst,pdamp,qdamp,ecmass,dtt,bintc,tprcor,soilq10,ansqr, &
          soilscaleold, nsoil, forcerestore, respg, pco2c, pco2i, &
          pco2s,co2cap,cflux,test, test10)

      implicit none

!
!
!=======================================================================
!
!     CALCULATION OF CANOPY CONDUCTANCE USING THE INTEGRATED   
!     MODEL RELATING ASSIMILATION AND STOMATAL CONDUCTANCE.
!     UNITS ARE CONVERTED FROM MKS TO BIOLOGICAL UNITS IN THIS ROUTINE.
!     BASE REFERENCE IS SE-92A
!
!                          UNITS
!                         -------
!
!      PCO2M, PCO2A, PCO2Ap, PCO2I, PO2M        : PASCALS
!      CO2A, CO2S, CO2I, H2OA, H2OS, H2OA       : MOL MOL-1
!      VMAX0, RESPN, ASSIM, GS, GB, GA, PFD     : MOL M-2 S-1
!      EFFCON                                   : MOL CO2 MOL QUANTA-1
!      GCAN, 1/RB, 1/RA, 1/RST                  : M S-1
!      EVAPKG                                   : KG M-2 S-1
!      Q                                        : KG KG-1
!
!                       CONVERSIONS
!                      -------------
!
!      1 MOL H2O           = 0.018 KG
!      1 MOL CO2           = 0.044 KG
!      H2O (MOL MOL-1)     = EA / PSUR ( MB MB-1 )
!      H2O (MOL MOL-1)     = Q*MM/(Q*MM + 1)
!pl the next line applies to the Ci to Cs pathway
!      GS  (CO2)           = GS (H2O) * 1./1.6
!pl 44.6 is the number of moles of air per cubic meter
!      GS  (MOL M-2 S-1 )  = GS (M S-1) * 44.6*TF/T*P/PO
!      PAR (MOL M-2 S-1 )  = PAR(W M-2) * 4.6*1.E-6
!      MM  (MOLAIR/MOLH2O) = 1.611
!
!
!                         OUTPUT
!                      -------------
!
!      ASSIMN              = CANOPY NET ASSIMILATION RATE
!      EA                  = CANOPY AIR SPACE VAPOR PRESSURE
!      1/RST               = CANOPY CONDUCTANCE
!      PCO2I               = INTERNAL CO2 CONCENTRATION
!      RESPC               = CANOPY RESPIRATION
!      RESPG               = GROUND RESPIRATION
!
!----------------------------------------------------------------------
!
!         RSTFAC(1) ( F(H-S) )               : EQUATION (17,18), SE-92A
!         RSTFAC(2) ( F(SOIL) )              : EQUATION (12 mod), SE-89
!         RSTFAC(3) ( F(TEMP) )              : EQUATION (5b)   , CO-92
!         RSTFAC(4) ( F(H-S)*F(SOIL)*F(TEMP))
!
!-----------------------------------------------------------------------
!

!++++++++++++++++++++++++++++++OUTPUT+++++++++++++++++++++++++++++++++++
!
!       ASSIMN         CARBON ASSIMILATION FLUX (MOL M-2 S-1) 
!       RST            CANOPY RESISTANCE (S M-1)
!       RSTFAC(4)      CANOPY RESISTANCE STRESS FACTORS 
!
!++++++++++++++++++++++++++DIAGNOSTICS++++++++++++++++++++++++++++++++++
!
!       RESPC          CANOPY RESPIRATION (MOL M-2 S-1)
!       RESPG          GROUND RESPIRATION (MOL M-2 S-1)
!       PCO2I          CANOPY INTERNAL CO2 CONCENTRATION (MOL MOL-1)
!       GSH2O          CANOPY CONDUCTANCE (MOL M-2 S-1)
!       H2OS           CANOPY SURFACE H2O CONCENTRATION (MOL MOL-1)
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!     Modifications:
!       - gs (stomatal conductance reduced for freezing soils per Jim Collatz
!         dd 950221
!      kevin Schaefer changed templ to temph in omss (11/8/01)
!      Kevin Schaefer added surf pressure to RRKK (11/8/01)
!      Kevin Schaefer for C3 OMS changed from templ to tempH (11/8/01)
!      Kevin Schaefer For C4 OMS added templ and tempH (11/8/01)      
!      Modified for multitasking - introduced gather/scatter indices
!          - DD 951206
!
!itb   Added in pco2c (chloroplast partial co2) for neil's fractionation
!itb   calculations
!itb       - IB Sep99

!
!     input arrays:
      integer len, nsib, nsoil
      logical forcerestore

      real vmax0(len),psur(len),green(len),gmudmu(len), &
         zlt(len),cas_cap_co2(len),tc(len),ta(len), &
         trop(len),trda(len), &
         slti(len),shti(len),hltii(len),hhti(len), &
         etgs(len),wc(len),em(len),ra(len),rb(len), &
         cog1(len),cog2(len),tm(len),effcon(len), &
         binter(len),gradm(len),atheta(len),btheta(len), &
         tgs(len),respcp(len),tran(nsib,2,2),ref(nsib,2,2), &
         radn(len,2,2),ecmass(len),trdm(len),etc(len), &
         aparc(len),www(len,2),tg(len),rst(len), cflux(len)
      real pdamp, qdamp, dtt, pco2m(len), pco2ap(len), tf, po2m
     
!     output arrays:


      real assimn(len),ea(len),rstfac(len,4), &
          pco2i(len),respc(len),respg(len),drst(len)
!zz new diagostics 10/14/92 
!
! output arrays
      real omepot(len),assimpot(len),assimci(len), &
          assimnp(len),whs(len),antemp(len), &
          wsfws(len),wsfht(len),wsflt(len),wci(len), &
          wags(len),wegs(len),pfd(len), &
          td(nsib,nsoil),zmlscale(len),assim(len), &
          wopt(len),zm(len),wsat(len), &
          soilscale(len,nsoil+1),zmstscale(len,2),zltrscale(len), &
          tprcor(len),bintc(len),soilq10(len,nsoil+1), &
          ansqr(len),soilscaleold(len), &
          pco2c(len), &
          xgah2o(len),  &
          xgco2m(len)
     
!     work arrays:

      real PCO2Y(len,6), EYY(len,6),assimny(len,6), &
          assimy(len,6)
      real c3(len),c4(len),range(len),gammas(len), &
          aparkk(len),gah2o(len), &
          gbh2o(len), &
          par(len),rrkk(len), &
          omss(len),vm(len),gsh2o(len),pco2s(len), &
          templ(len),temph(len), &
          qt(len),co2s(len),scatp(len),scatg(len), &
          park(len),respn(len),zkc(len),resa(len), &
          zko(len),spfy(len), co2a(len), co2m(len),co2cap(len)

      integer icconv(len),igath(len)

      real soilfrz(len)
      real bl, tresp, soilq10td, b(2), woptzm, cwsflt, cwsfht, cwsfws, &
          ccoms, ccomc, ascitemp, dompdomc, omsci, ompci, omcci, &
          omcpot, omppot, sqrtin, omspot, pco2ipot, ohtp2, sttp2, &
          gsh2oinf, h2osrh, h2os, ecmole, h2oa, h2oi, dtti, &
          pco2in, pco2a, soilfrztd, soilfrztg, soilq10tg,  &
          rstnew

      integer i, ic1, ic, l
!
      real test, test10(10)   ! test variables

      
      
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


      do i = 1,len                !  LOOP OVER GRIDPOINT


        TPRCOR(i) = TF*PSUR(i)*100./1.013e5
        co2cap(i) = cas_cap_co2(i) * 44.6 * tprcor(i)/ta(i)     ! moles air / m2

!pl this needs to be modified as in sibslv3 to automatically use the 
!pl thickness of the canopy air space. 
      
!
!----------------------------------------------------------------------
!
!pl        RESPG(i) = 0. E -6 ! fixed respiration at 5 micromoles
!pl   no longe rused since we now have respsib
!
!----------------------------------------------------------------------
!
        IF( EFFCON(i) .GT. 0.07 ) then
          C3(i) = 1.
        else
          C3(i) = 0.
        endif
         C4(i)     = 1. - C3(i)

!
!-----------------------------------------------------------------------
!
!     CALCULATION OF CANOPY PAR USE PARAMETER.
!
!      APARKK      (PI)     : EQUATION (31) , SE-92A
!-----------------------------------------------------------------------
!  
        SCATP(I) =     GREEN(i)   * &
                 ( TRAN(i,1,1) + REF(i,1,1) ) &
                +( 1.-GREEN(i) ) * &
                ( TRAN(i,1,2) + REF(i,1,2) )
        SCATG(i) = TRAN(i,1,1) + REF(i,1,1)
        PARK(i) = SQRT(1.-SCATP(i)) * GMUDMU(i)
!
! Collatz-Bounoua commented the calculation of  aparc
! replaced it with theone calculated in new_mapper.
!
!b        APARC(i) = 1. - EXP ( -PARK(i)*ZLT(i) )   ! lahouari
!
        APARKK(i)   = APARC(i) / PARK(i) * GREEN(i)
!-----------------------------------------------------------------------
!
!     Q-10 AND STRESS TEMPERATURE EFFECTS
!
!      QT          (QT)    : TABLE (2)     , SE-92A
!-----------------------------------------------------------------------
!
        qt(i) = 0.1*( TC(i) - TROP(i) )
        RESPN(i) = RESPCP(i) * VMAX0(i) * RSTFAC(i,2)
	print*, RSTFAC(i,2)


!itb...patch to prevent underflow if temp is too cool...
	if(TC(i) >= TRDM(i))then
            RESPC(i) = RESPN(i) * 2.0**qt(i) &
             /( 1. + EXP( TRDA(i)*(TC(i)-TRDM(i))))
        else
            RESPC(i) = RESPN(i) * 2.0**qt(i)
        endif

        VM(i) = VMAX0(i) * 2.1**qt(i)
        TEMPL(i) = 1. + EXP(SLTI(i)*(HLTIi(i)-TC(i)))
        TEMPH(i) = 1. + EXP(SHTI(i)*(TC(i)-HHTI(i)))
        RSTFAC(i,3) = 1./( TEMPL(i)*TEMPH(i))
        VM(i)    = VM(i)/TEMPH(i) * RSTFAC(i,2)*C3(i) &
           + VM(i) * RSTFAC(i,2)*RSTFAC(i,3) * C4(i)
!
!-----------------------------------------------------------------------
!
!     MICHAELIS-MENTEN CONSTANTS FOR CO2 AND O2, CO2/O2 SPECIFICITY,
!     COMPENSATION POINT       
!
!      ZKC          (KC)     : TABLE (2)     , SE-92A
!      ZKO          (KO)     : TABLE (2)     , SE-92A
!      SPFY         (S)      : TABLE (2)     , SE-92A
!      GAMMAS       (GAMMA-*): TABLE (2)     , SE-92A
!      OMSS         (OMEGA-S): EQUATION (13) , SE-92A
!      BINTC        (B*ZLT)  : EQUATION (35) , SE-92A
!-----------------------------------------------------------------------
!
        ZKC(i) = 30. * 2.1**qt(i)
        ZKO(i) = 30000. * 1.2**qt(i)
        SPFY(i) = 2600. * 0.57**qt(i)
        GAMMAS(i) = 0.5 * PO2M/SPFY(i) * C3(i)
        PFD(i)    = 4.6E-6 * GMUDMU(i)* &
                  ( RADN(i,1,1)+RADN(i,1,2) )
!
!pl these here all go from being m/s to being mol/ (m2 sec)
        GSH2O(i)  = 1.0/RST(i) * 44.6*TPRCOR(i)/TC(i)
        GBH2O(i)  = 0.5/RB(i) * 44.6*TPRCOR(i)/TC(i)
        GAH2O(i)  = 1.0/RA(i) * 44.6*TPRCOR(i)/TM(i)

        xgah2o(i) = max(0.466, gah2o(i) )
        xgco2m(i) = 4000.0 * vmax0(i) * APARKK(i) * RSTFAC(i,2)

!
! Kevin Schaefer added surf pressure to rrkk (11/8/01)
        RRKK(i)   = ZKC(i)*(1.+ PO2M/ZKO(i) ) * C3(i) &
                    + VMAX0(i)*(1.8**qt(i))*200./Psur(i)*C4(i)
        PAR(i)    = pfd(i)*EFFCON(i)*( 1.-SCATG(i) )
        soilfrztg = 1.+exp(-1.5 * (max(270.,tgs(i))-273.16) )
        soilfrztd = 1.+exp(-1.5 * (max(270.,td (i,nsoil))-273.16) )
        soilfrz(i) = max(1./soilfrztg, 1./soilfrztd)
        soilfrz(i) = max( soilfrz(i), 0.05)
        BINTC(i)  = BINTER(i)*ZLT(i)*GREEN(i)* &
                      RSTFAC(i,2) * soilfrz(i)
!
! Kevin Schaefer for C3 changed from templ to tempH (11/8/01)
! Kevin Schaefer For C4 added templ and tempH (11/8/01)
        OMSS(i)=(VMAX0(i)/2.0)*(1.8**qt(i)) &
                       /TEMPH(i)*RSTFAC(i,2)*C3(i) &
                       + RRKK(i)*RSTFAC(i,3)*RSTFAC(i,2)*C4(i)
!
!-----------------------------------------------------------------------
!
!     FIRST GUESS IS MIDWAY BETWEEN COMPENSATION POINT AND MAXIMUM
!     ASSIMILATION RATE.
!
!-----------------------------------------------------------------------
 
       
        RANGE(i)    = PCO2M(i) * ( 1. - 1.6/GRADM(i) ) - GAMMAS(i)
	    icconv(i) = 1
	
      enddo
      
!
      DO 1000 IC = 1, 6
        do i = 1,len        ! LOOP OVER GRIDPOINT
          PCO2Y(i,IC) = 0.
          EYY(i,IC) = 0.
	enddo
1000  CONTINUE
!

!pl beginning of PL's setup

!      do i=1,len
!       gah2o(i) =  1. / MAX(0.446 E 0,GAH2O(i))
!      enddo

      DO 2000 IC = 1, 6
!
        CALL       SORTIN( EYY, PCO2Y, RANGE, GAMMAS, ic,len )
        
        CALL       CYCALC_Suits( APARKK, VM, ATHETA, BTHETA,par, &
                          GAMMAS, RESPC, RRKK, OMSS, C3, C4, &
                          PCO2Y(1,ic), assimny(1,ic), assimy(1,ic), &
                          len,test, test10  )
!

      do i = 1,len

!pl first diagnose the current CO2 flux in mol / (m2 s)
!pl this is a modified ra that will get us the right units
!pl in the conservation equation. its units are m2 s / (mol_air)

!	resa(i) =   1. / MAX(0.446 E 0,GAH2O(i))

!	resb(i) =   1.4/GBH2O(i)

!	resc(i) =   1.6/GSH2O(i) 
     
!pl now prognose the new CAS CO2 according to flux divergence
!pl we are going to do this in mol C / mol air (same as PaC/PaAir)

        CO2A(i)    = PCO2Ap(i) /   (PSUR(i)*100.)
        co2m(i)    = pco2m(i)  /   (PSUR(i)*100.)   

        CO2A(i)   = (  CO2A(i) + (dtt/co2cap(i)) *  &
                    (respg(i) - assimny(i,ic) &
                   +co2m(i)*gah2o(i)        ) ) &
                 / (1+dtt*gah2o(i)/ co2cap(i) ) 

        pco2a = co2a(i) * psur(i) * 100.

        PCO2S(i) = PCO2A - (1.4/GBH2O(i) * ASSIMNy(i,ic) &
                 * PSUR(i)*100.)
     
!  NSS (7/18/01); J Collatz suggested that we iterate for pco2c 
!  instead of pco2i.  I have instituted that code here, although I 
!  have done so using 'arguments' already defined in our simulations.
    
        PCO2I(i) =PCO2S(i)-(1.6/GSH2O(i)*ASSIMNy(i,ic)*PSUR(i)*100.)

        PCO2c(i)=PCO2I(i)-(1./xgco2m(i)*ASSIMNy(i,ic)*PSUR(i)*100.)
        EYY(i,IC)=PCO2Y(i,IC)-PCO2c(i)
     
!        PCO2IN   = PCO2S(i) - (1.6/GSH2O(i) * ASSIMNy(i,ic)
!     >            * PSUR(i)*100.)
!        EYY(i,IC) = PCO2Y(i,IC) - PCO2IN

      enddo
!
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
!
2000  CONTINUE
!   
!pl end of PL's setup

      do i = 1,len        ! LOOP OVER GRIDPOINT
        icconv(i) = min(icconv(i),6)
        igath(i) = i+(icconv(i)-1)*len
      enddo


      do i = 1,len         ! LOOP OVER GRIDPOINT

        pco2c(i) = pco2y(igath(i),1)
        assimn(i) = assimny(igath(i),1)
        assim(i)  = assimy(igath(i),1)


	
!        pco2c(i) = pco2i(i) - assimn(i)/xgco2m(i)*psur(i)*100.0

!pl now do the real C_A forecast with the iterated fluxes.

        CO2A(i)    = PCO2Ap(i) /   (PSUR(i)*100.)
        co2m(i)    = pco2m(i)  /   (PSUR(i)*100.)   
   
        CO2A(i) = (CO2A(i) + (dtt/co2cap(i)) * &
                   (respg(i) - assimn(i) &
                    +co2m(i)*gah2o(i) ) ) &
                 / (1+dtt*gah2o(i)/co2cap(i))
!pl go back from molC / mol air to Pascals 

        pco2ap(i) = co2a(i) * psur(i) * 100.

!itb...carbon flux between CAS and reference level
         cflux(i) = gah2o(i)*(co2a(i)-co2m(i))*0.012

!        PCO2Ap(i) =  PCO2Ap(i) + (dtt/co2cap) 
!     &            * (  (PSUR(i) *100. * respg(i))
!     &            +    (PCO2i(i)   /(resb+resc))
!     &            +    (pco2m(i)   /       resa)  )
!     &            / (1.+ (dtt/co2cap) * (1./resa + 1./(resb+resc)) )
      enddo
!
      dtti = 1./dtt
      do i = 1,len        ! LOOP OVER GRIDPOINT
!czzggrst5 - new code
        H2OI   = ETC(i)/PSUR(i)
        H2OA   =  EA(i)/PSUR(i)
        ECMOLE = 55.56 * ECMASS(i) * dtti  ! ecmass must be computed and passed in
        H2OS = H2OA + ECMOLE / GBH2O(i)
        H2OS  = MIN( H2OS, H2OI )
        H2OS  = MAX( H2OS, 1.e-7)
        H2OSRH = H2OS/H2OI
!  need qdamp and pdamp calculated and passed to here!
!pl        CO2S(i) = MAX(PCO2S(I),PCO2M*0.5) / (PSUR(i)*100.)

!pl I have relaxed this condition to 1/10 of previous. The old way made
!pl the CO2 on top of the leaves always at least 1/2 of the value at the
!pl reference level.

        CO2S(i) = MAX(PCO2S(I),PCO2M(i)*0.05) / (PSUR(i)*100.)
        
!pl Ball-Berry equation right here !        

        GSH2OINF = (GRADM(i) * MAX(1.e-12,ASSIMN(i)) &
                 * H2OSRH * soilfrz(i) / CO2S(i)) + BINTC(i)
     
!pl this is the change in stomatal resistance

        DRST(i) = RST(i) * QDAMP * ((GSH2O(i)-GSH2OINF)/ &
                 (PDAMP*GSH2O(i)+QDAMP*GSH2OINF))


!pl this is the 'would be change' if we did not use the damping factor.

!        rstnew = (1./gsh2oinf) * 44.6 * tprcor(i)/tc(i)
!        DRST(i) = rstnew - RST(i)

!
        RSTFAC(i,1) = H2OS/H2OI
        RSTFAC(i,4) = RSTFAC(i,1)*RSTFAC(i,2)* RSTFAC(i,3)
      enddo
!
!Z CARNEGIE new diagnostics----start!!!(c.zhang&joe berry, 10/19/92)
!-----------------------------------------------------------------------
!  INPUTS: PSUR(i),CO2S,ASSIMN(i),GRADM(i),BINTC(i),VMAX0(i),RRKK(i),C3(i),
!    C4(i),PAR(i),ATHETA(i),BTHETA(i),APARKK(i),OMSS(i),RSTFAC(i,2),TEMPH,
!    TEMPL,RSTFAC(i,1),VM(i),ASSIM,GSH20(i),EFFCON(i),QT,GAMMAS(i),
!    PFD(i)
!

        sttp2 = 73.**0.2
        ohtp2 = 100.**0.2
        do i = 1,len
!-----------------------------------------------------------------------
! CALCULATION OF POTENTIAL ASSIMILATION
!-----------------------------------------------------------------------
!
! Make assimn a top leaf, not the canopy.
      ASSIMNp(i) = ASSIMN(i) / APARKK(i)
!
! Bottom stopped assim.
      ANTEMP(i) = MAX(0.,ASSIMNp(i))
!
! Potential intercellular co2.
      PCO2IPOT = PSUR(i)*100.*(co2s(i)-(1.6*ASSIMNp(i)/ &
            ((GRADM(i)*ANTEMP(i)/co2s(i))+BINTC(i))))
!
! Potential rubisco limitation.
      OMCPOT = VMAX0(i)*2.1**qt(i)*((PCO2IPOT-GAMMAS(i))/ &
            (PCO2IPOT+RRKK(i))*C3(i) + C4(i))
!
! Potential light limitation.
      OMEPOT(i) = PAR(i)*((PCO2IPOT-GAMMAS(i))/ &
                (PCO2IPOT+2.*GAMMAS(i))*C3(i) + C4(i))
!
! Quad 1.
      SQRTIN = MAX(0.,((OMEPOT(i)+OMCPOT)**2- &
                    4.*ATHETA(i)*OMEPOT(i)*OMCPOT))
!
! Quad 1. Intermediate  top leaf photosynthesis.
      OMPPOT = ((OMEPOT(i)+OMCPOT)-SQRT(SQRTIN))/(2.*ATHETA(i))
!
! Potential sink or pep limitation.
      OMSPOT = (VMAX0(i)/2.0)*(1.8**qt(i))*C3(i) &
                      + RRKK(i)*PCO2IPOT*C4(i)
!
! Quad 2.
      SQRTIN=MAX(0.,((OMPPOT+OMSPOT)**2-4.*BTHETA(i)* &
             OMPPOT*OMSPOT))
!
! Quad 2. Final Potential top leaf photosynthesis.
      ASSIMPOT(i) = ((OMSPOT+OMPPOT)-SQRT(SQRTIN))/(2.*BTHETA(i))
!
!-----------------------------------------------------------------------
! CALCULATION OF STRESS FACTOR LIMITED ASSIMILATION
!-----------------------------------------------------------------------
!
! Stressed rubisco limitation.
      OMCCI = VM(i)*((PCO2IPOT-GAMMAS(i))/(PCO2IPOT+RRKK(i))*C3(i) &
                   + C4(i))
!
! Quad 1.
      SQRTIN = MAX(0.,(OMEPOT(i)+OMCCI)**2- &
                     4.*ATHETA(i)*OMEPOT(i)*OMCCI)
!
! Quad 1. Intermediate stress limited top leaf photosynthesis.
      OMPCI = ((OMEPOT(i)+OMCCI)-SQRT(SQRTIN))/(2.*ATHETA(i))
!
! Stressed sink or pep limitation.
      OMSCI = OMSS(i)*(C3(i) + PCO2IPOT*C4(i))
!
! Quad 2.
      SQRTIN = MAX(0.,(OMPCI+OMSCI)**2-4.*BTHETA(i)*OMPCI*OMSCI)
! 
! Quad 2. Final stress limited top leaf photosynthesis.
      ASSIMCI(i) = ((OMSCI+OMPCI)-SQRT(SQRTIN))/(2.*BTHETA(i))
!
!-----------------------------------------------------------------------
! CALCULATION OF CONTROL COEFFICIENTS
!-----------------------------------------------------------------------
!
! Intermediate.
      DOMPDOMC = (OMPCI-OMEPOT(i))/ &
               (2.*ATHETA(i)*OMPCI-OMCCI-OMEPOT(i))
!
! Bottom stopped final stress limited top leaf photosynthesis.
      ASCITEMP = MAX(ASSIMCI(i),1.e-12)
!
! Rubisco control coefficient.
      CCOMC = (DOMPDOMC*(ASSIMCI(i)-OMSCI)/ &
           (2.*BTHETA(i)*ASSIMCI(i)-OMPCI-OMSCI))*OMCCI/ASCITEMP
!
! Sink or pep control coefficient.
      CCOMS = ((ASSIMCI(i)-OMPCI)/ &
            (2.*BTHETA(i)*ASSIMCI(i)-OMPCI-OMSCI))*OMSCI/ASCITEMP
!
!-----------------------------------------------------------------------
!  OUTPUT:  POTENTIAL ASSIMILATION RATES TO BE SUMMED
!-----------------------------------------------------------------------
! Canopy values (overwrites top leaf).
! 
      OMEPOT(i) = OMEPOT(i)*APARKK(i)
      ASSIMPOT(i) = ASSIMPOT(i)*APARKK(i)
      ASSIMCI(i) = ASSIMCI(i)*APARKK(i)
      ASSIM(i) = ASSIM(i)*APARKK(i)
      ANTEMP(i) = ANTEMP(i)*APARKK(i)
      ANSQR(i) = ANTEMP(i)*ANTEMP(i)
      ASSIMNp(i) = ASSIMNp(i)*APARKK(i)
!
!-----------------------------------------------------------------------
! OUTPUT:  WEIGHTED STRESS FACTORS AND OTHER DIAGNOSTIC OUTPUTS TO BE SUMMED
!-----------------------------------------------------------------------
!
! Water stress.
      WSFWS(i) = ASSIMPOT(i)*(1.-RSTFAC(i,2))*(CCOMC+CCOMS)
!
! High temperature stress.
      WSFHT(i) = ASSIMPOT(i)*(1.-1./TEMPH(i))*CCOMC
!
! Low temperature stress.
      WSFLT(i) = ASSIMPOT(i)*(1.-1./TEMPL(i))*(CCOMS*C3(i)+CCOMC*C4(i))
!
!  protection for wsfws, wsfht, and wsflt from <0 or >>xxx(2/24/93)
      cwsfws = (1.-RSTFAC(i,2))*(CCOMC+CCOMS)
      if(cwsfws.gt.1. .or. cwsfws.lt.0.) wsfws(i)=0.
      cwsfht = (1.-1./TEMPH(i))*CCOMC
      if(cwsfht.gt.1. .or. cwsfht.lt.0.) wsfht(i)=0.
      cwsflt = (1.-1./TEMPL(i))*(CCOMS*C3(i)+CCOMC*C4(i))
      if(cwsflt.gt.1. .or. cwsflt.lt.0.) wsflt(i)=0.

!
! Intermediate assimilation weighted Ci.
      WCI(i) = ANTEMP(i)*PCO2I(i)
!
! Intermediate assimilation weighted relative humidty stress factor.
      WHS(i) = ANTEMP(i)*RSTFAC(i,1)
!
! Intermediate assimilation weighted stomatal conductance.
      WAGS(i) = GSH2O(i)*ANTEMP(i)
!
! Intermediate evaporation weighted stomatal conductance.(Step 1.
!   Step 2 after subroutine updat2)
      WEGS(i) = GSH2O(i)
!
!      bl = (((100.* www(i,1))**0.2 - sttp2)/(sttp2 - ohtp2))**2
!      bl = min(bl,10. E 0)
!      zmlscale(i) = 0.8*0.75**bl + 0.2
!
!      zltrscale(i) = (exp(0.0693*(tgs(i)-298.15)))
!      if (zltrscale(i) .lt. 0.17 ) zltrscale(i)= 0.
!      zltrscale(i) = zltrscale(i) * zmlscale(i)

!      make these zero until canned from the qp diagnostic, then scrap
      soilscaleold(i) = 0.0
       zmlscale(i) = 0.0
       zltrscale(i) = 0.0
      enddo

	  test10(1)=RESPC(1)*APARKK(1)*1.e6
          test10(2)=assim(1)*1.e6


      RETURN
      END
