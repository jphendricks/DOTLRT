!
!=========================================================================
subroutine phosib(sib,sib_loc)
!=========================================================================
! CALCULATION OF CANOPY CONDUCTANCE USING THE INTEGRATED   
! MODEL RELATING ASSIMILATION AND STOMATAL CONDUCTANCE.
! UNITS ARE CONVERTED FROM MKS TO BIOLOGICAL UNITS IN THIS ROUTINE.
! BASE REFERENCE IS SE-92A
!
! UNITS
!  PCO2M, PCO2A, PCO2Ap, PCO2I, PO2M        : PASCALS
!  CO2A, CO2S, CO2I, H2OA, H2OS, H2OA       : MOL MOL-1
!  VMAX0, RESPN, ASSIM, GS, GB, GA, PFD     : MOL M-2 S-1
!  EFFCON                                   : MOL CO2 MOL QUANTA-1
!  GCAN, 1/RB, 1/RA, 1/RST                  : M S-1
!  EVAPKG                                   : KG M-2 S-1
!  Q                                        : KG KG-1
!
! CONVERSIONS
!  1 MOL H2O           = 0.018 KG
!  1 MOL CO2           = 0.044 KG
!  H2O (MOL MOL-1)     = EA / PSUR ( MB MB-1 )
!  H2O (MOL MOL-1)     = Q*MM/(Q*MM + 1)
!  GS  (CO2)           = GS (H2O) * 1./1.6 (pl applies to the Ci to Cs pathway)
!  GS  (MOL M-2 S-1 )  = GS (M S-1) * 44.6*TF/T*P/PO (pl 44.6 is number of moles air per cubic meter)
!  PAR (MOL M-2 S-1 )  = PAR(W M-2) * 4.6*1.E-6
!  MM  (MOLAIR/MOLH2O) = 1.611
!
! OUTPUT and diagnostics
!  EA                  = CANOPY AIR SPACE VAPOR PRESSURE
!  RSTFAC(1) ( F(H-S) )               : EQUATION (17,18), SE-92A
!  RSTFAC(2) ( F(SOIL) )              : EQUATION (12 mod), SE-89
!  RSTFAC(3) ( F(TEMP) )              : EQUATION (5b)   , CO-92
!  RSTFAC(4) ( F(H-S)*F(SOIL)*F(TEMP))
!  ASSIMN         CARBON ASSIMILATION FLUX (MOL M-2 S-1) 
!  RST            CANOPY RESISTANCE (S M-1)
!  RSTFAC(4)      CANOPY RESISTANCE STRESS FACTORS 
!  resp_can       CANOPY Autotrophic RESPIRATION (MOL M-2 S-1)
!  resp_grnd      GROUND RESPIRATION (MOL M-2 S-1)
!  PCO2I          CANOPY INTERNAL CO2 CONCENTRATION (MOL MOL-1)
!  GSH2O          CANOPY CONDUCTANCE (MOL M-2 S-1)
!  H2OS           CANOPY SURFACE H2O CONCENTRATION (MOL MOL-1)
!
! Modifications:
!  gs (stomatal conductance reduced for freezing soils per Jim Collatz (dd 950221)      
!  Modified for multitasking - introduced gather/scatter indices (dd 951206)
!  Ian Baker Added in pco2c (chloroplast partial co2) for neil's fractionation
!    calculations (Sep 1999)
!  Kevin Schaefer resp_can/assim multiplied by aparkk when calculated (5/31/05)
!  Kevin Schaefer changed convergence check so that it actually worked (5/22/13)
!  Kevin Schaefer corrected calculation of rrkk and omss for C4 (5/23/13)
!  Kevin Schaefer applied frost_stress, which was deleted somehow (5/23/13)
!  Kevin Schaefer deleted rrkk completely (5/23/13)
!  Kevin Schaefer deleted all potential assim diagnostics (wegs, wags, whs, wci,
!    wsfht, wsflt, wsfws, ansqr, antemp, omepot, assimnp, assimci, assimpot) (5/23/13)
!  Kevin Schaefer moved potential rate calculations from cycalc (5/24/13)
!  Kevin Schaefer corrected frost recovery from 4 to 2 C/day (5/24/13)
!-------------------------------------------------------------------------

use kinds
use sibtype
use cfrax
use sib_const_module, only: &
    po2m, &
    pdb,  &
    dtt,  &
    dti
use physical_parameters, only: &
    p0 => p0_sfc, &
    tice

implicit none
!
! input/output variables
type(sib_t), intent(inout) :: sib
type(sib_local_vars), intent(inout) :: sib_loc ! variables local to SiB
!
! local variables
integer icconv    !
integer ic  ! iteration loop
integer i   ! loop variable
real(kind=dbl_kind) co2cap    ! conversion factor from ppm to mole/m2 in the CAS
real(kind=dbl_kind) c3        ! C3 flag
real(kind=dbl_kind) c4        ! C4 flag
real(kind=dbl_kind) scatp     ! 
real(kind=dbl_kind) scatg     !
real(kind=dbl_kind) park      !
real(kind=dbl_kind) qt        !
real(kind=dbl_kind) respn     !
real(kind=dbl_kind) vm        !
real(kind=dbl_kind) templ     !
real(kind=dbl_kind) temph     !
real(kind=dbl_kind) zkc       !
real(kind=dbl_kind) zko       !
real(kind=dbl_kind) spfy      !
real(kind=dbl_kind) gammas    !
real(kind=dbl_kind) gah2o     !
real(kind=dbl_kind) gbh2o     !
real(kind=dbl_kind) gsh2o     !
real(kind=dbl_kind) xgco2m    !
real(kind=dbl_kind) soilfrz   !
real(kind=dbl_kind) soilfrztg !
real(kind=dbl_kind) soilfrztd !
real(kind=dbl_kind) bintc     !
real(kind=dbl_kind) omss      !
real(kind=dbl_kind) range1    !
real(kind=dbl_kind) par       !
real(kind=dbl_kind) co2a      ! CAS CO2 concentration (mol C/mol air)
real(kind=dbl_kind) co2m      ! reference level CO2 concentration (mol C/mol air)
real(kind=dbl_kind) co2s      ! leaf surface CO2 concentration (mol C/mol air)
real(kind=dbl_kind) pco2a     ! intermediate CAS CO2 concentration (Pa)
real(kind=dbl_kind) pco2s     ! intermediate leaf surface CO2 partial pressure (Pa)
real(kind=dbl_kind) pco2i     ! intermediate stomatal (internal) CO2 partial pressure (Pa)
real(kind=dbl_kind) pco2c     ! intermediate leaf chloroplast CO2 partial pressure (Pa)
real(kind=dbl_kind) h2oi      !
real(kind=dbl_kind) h2oa      !
real(kind=dbl_kind) h2os      !
real(kind=dbl_kind) h2osrh    !
real(kind=dbl_kind) ecmole    !
real(kind=dbl_kind) pco2y(6)  !
real(kind=dbl_kind) eyy(6)    !
real(kind=dbl_kind) assimny(6)!
real(kind=dbl_kind) assimy(6) !
real(kind=dbl_kind) gsh2oinf  !
real(kind=dbl_kind) drst(5)   ! delta of stomatal resistance (sec/m)
real(kind=dbl_kind) rstfac3(6)! intermediate temperature stress factor
real(kind=dbl_kind) zln2      ! used in calculating pdamp,qdamp
real(kind=dbl_kind) ghalf     ! used in calculating pdamp,qdamp
real(kind=dbl_kind) dttin     ! used in calculating pdamp,qdamp
real(kind=dbl_kind) dmin      ! used in calculating pdamp,qdamp
real(kind=dbl_kind) pdamp
real(kind=dbl_kind) qdamp
real(kind=dbl_kind) tprcor    ! temperature correction (K)
real(kind=dbl_kind) rstar     ! universal gas constant (N m mole^-1 K^-1)
real(kind=dbl_kind) ome       ! (mol/m2/s) light limited assimilation
real(kind=dbl_kind) omc       ! (mol/m2/s) Rubisco limited assimilation
real(kind=dbl_kind) omp       ! (mol/m2/s) weighted omc/ome limited assimilation
real(kind=dbl_kind) oms       ! (mol/m2/s) sink/export limited assimilation
real(kind=dbl_kind) sqrtin    ! quadratic weighting factor
!
! assign universal gas constant
    rstar = 8.3143
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
! calculate damping factors
    zln2 = 6.9314718e-1
    ghalf = 1.0257068e1
    dttin = 3.6e3 
    dmin = 6.0e1
    pdamp = exp (-1.0 * zln2*(dtt*dmin)/(dttin*ghalf))
    qdamp = 1.0 - pdamp
    tprcor = tice*sib%drvr%ps*100.0/p0
    co2cap = sib%diag%cas_cap_co2 * 44.6 * tprcor/sib%prog%ta  ! moles air / m2
    co2cap = sib%diag%cas_cap_co2 * sib%drvr%ps*100.0 /rstar/sib%prog%ta
!
! canopy PAR use factor
    scatp = sib%param%green * (sib%param%tran(1,1)+sib%param%ref(1,1))   &
     + (1.-sib%param%green) * (sib%param%tran(1,2)+sib%param%ref(1,2))
    scatg =  sib%param%tran(1,1) +  sib%param%ref(1,1)
    park = sqrt(1.-scatp)* sib%param%gmudmu
!
!itb...Niall integrates physfrac into aparkk. I'm not sure I like
!itb...doing it that way--SO I WON'T, FOR NOW...
    sib%diag%aparkk   = sib%param%aparc / park * sib%param%green
!
!itb...start PHYSIOLOGY LOOP here...
!itb...potentially 5 different physiologies can share the same
!itb...soil and CAS (making it different from a normal tile).
!itb...loop will cycle out of an unused physiology type.
!
! zero out physiology-specific values
    sib%prog%rst(6)     = 0.0
    sib%diag%pco2c(6)   = 0.0
    sib%diag%pco2i(6)   = 0.0
    sib%diag%pco2s(6)   = 0.0
    sib%diag%assimn(6)  = 0.0
    sib%diag%assim(6)   = 0.0
    sib%diag%ggl(6)     = 0.0
    sib%diag%resp_can(6)= 0.0
    rstfac3(6)          = 0.0
!
! Physiology loop
    phys_loop : do i=1,5
      if ( sib%param%physfrac(i) == 0.0 ) cycle phys_loop
!
! set physiology type
        if( sib%param%phystype(i) == 3) then
            c3 = 1._dbl_kind
            c4 = 0._dbl_kind
        elseif(sib%param%phystype(i) == 4) then
            c3 = 0._dbl_kind
            c4 = 1._dbl_kind
        else
            print*,'loop index=',i,' phystype=',sib%param%phystype(i)
            stop'ERROR:UNKNOWN PHYSIOLOGY TYPE IN PHOSIB'
        endif
!
!-------------------------------------------------------------------------
! Q-10 AND STRESS TEMPERATURE EFFECTS
!      QT          (QT)    : TABLE (2)     , SE-92A
!-------------------------------------------------------------------------
!
! Temperature response quocient
        qt = 0.1 * ( sib%prog%tc - sib%param%trop(i) )
        respn = sib%param%respcp(i) * sib%param%vmax0(i) * sib%diag%rstfac(2)
!
! canopy autotrophic respiration, with patch to prevent underflow if temp is too cool...
        if(sib%prog%tc >= sib%param%trdm(i) )then
            sib%diag%resp_can(i) = respn * 2.0**qt/( 1. + EXP( sib%param%trda(i)*(sib%prog%tc-sib%param%trdm(i) )))
        else
            sib%diag%resp_can(i) = respn * 2.0**qt
        endif
        sib%diag%resp_can(i)=sib%diag%resp_can(i)*sib%diag%aparkk
!
! Low and high temperature inhibition functions
        templ = 1. + EXP(sib%param%slti(i) * (sib%param%hltii(i) - sib%prog%tc))
        temph = 1. + EXP(sib%param%shti(i) * (sib%prog%tc - sib%param%hhti(i) ))
!
! frost stress: bottom-stopped min allowed canopy temperature at -20C
        if(sib%prog%tc < sib%param%tfrost) sib%prog%tcmin = sib%prog%tc
        sib%prog%tcmin = MAX(sib%prog%tcmin, 253.15_dbl_kind)
        sib%diag%frost_stress=(sib%prog%tcmin-tice)/(sib%param%tfrost-tice)
        sib%diag%frost_stress=min(sib%diag%frost_stress,1._dbl_kind)
        sib%diag%frost_stress=max(sib%diag%frost_stress,0._dbl_kind)
!
! frost stress: min canopy temp recovery at tcmin_rr (C/day)
        if(sib%prog%tc > sib%param%tfrost) sib%prog%tcmin = sib%prog%tcmin + (sib%param%tcmin_rr *dtt/86400.0_dbl_kind)
      sib%prog%tcmin=min(sib%prog%tcmin,sib%prog%tc)
!
! overall temperature scaling factor
        rstfac3(i) = 1./( templ*temph)
!
! temperature scaled Vmax        
        vm = sib%param%vmax0(i)*sib%diag%rstfac(2)*sib%diag%frost_stress * 2.1**qt
        vm = vm/temph*c3 + vm*rstfac3(i)*c4
!
!-------------------------------------------------------------------------
! MICHAELIS-MENTEN CONSTANTS FOR CO2 AND O2, CO2/O2 SPECIFICITY, COMPENSATION POINT       
!      ZKC          (KC)     : TABLE (2)     , SE-92A
!      ZKO          (KO)     : TABLE (2)     , SE-92A
!      SPFY         (S)      : TABLE (2)     , SE-92A
!      GAMMAS       (GAMMA-*): TABLE (2)     , SE-92A
!      OMSS         (OMEGA-S): EQUATION (13) , SE-92A
!      BINTC        (B*ZLT)  : EQUATION (35) , SE-92A
!-------------------------------------------------------------------------
        zkc     = 30. * 2.1**qt
        zko     = 30000. * 1.2**qt
        spfy    = 2600. * 0.57**qt
        gammas  = 0.5 * po2m/spfy * c3
!
! photon flux density with underflow check
      if((sib%drvr%radvbc+sib%drvr%radvdc)<1.E-6) then
            sib%diag%pfd=0.0
        else
            sib%diag%pfd=4.6E-6*sib%param%gmudmu*(sib%drvr%radvbc+sib%drvr%radvdc)
        endif
!
!...convert resistance to conductance  to  mol/ (m2 sec)
!...(44.6 mol m^-3 conversion factor)
        if ( sib%prog%rst(i) == 0. ) sib%prog%rst(i) = sib%prog%rst(1)
        gsh2o  = 1.0/sib%prog%rst(i) * 44.032476*tprcor/sib%prog%tc
        gbh2o  = 0.5/sib%diag%rb     * 44.032476*tprcor/sib%prog%tc
        gah2o  = 1.0/sib%diag%ra     * 44.032476*tprcor/sib%drvr%tm

        xgco2m = 4000.0*sib%param%vmax0(i)*sib%diag%aparkk*sib%diag%rstfac(2)
        par    = sib%diag%pfd*sib%param%effcon(i)*(1.-scatg)

        soilfrztg = 1.+exp(-1.5 * (max(270.0_dbl_kind,sib%prog%td(1))-273.16))
        soilfrztd = 1.+exp(-1.5 * (max(270.0_dbl_kind,sib%prog%td(2))-273.16))
        soilfrz   = max(1./soilfrztg, 1./soilfrztd)
        soilfrz   = max( soilfrz, 0.05_dbl_kind)

        bintc  = sib%param%binter(i) * sib%param%zlt * sib%param%green * sib%diag%rstfac(2) * soilfrz
        omss   = sib%param%vmax0(i)*sib%diag%rstfac(2)*sib%diag%frost_stress*1.8**qt* &
          (c3/2._dbl_kind/temph+c4/sib%drvr%ps/100.*1.e4/templ/temph)
!
!-------------------------------------------------------------------------
! iteration loop
!-------------------------------------------------------------------------
!Bio We iterate on PCO2C-sortin makes a 'first guess' at
!Bio then orders PCO2C/Error pairs on increasing error size,
!Bio then uses a combination of linear/quadratic fit to obtain 
!Bio the 'next best guess' as iteration count increases.
!Bio CYCALC uses that value of PCO2C to get the next value 
!Bio of ASSIMN. CO2A and CO2S follow.
!
! first guess midway between compensation point and max assimilation rate
        range1  = sib%drvr%pco2m * ( 1. - 1.6/sib%param%gradm(i) ) - gammas
        pco2y = 0.
        eyy   = 0.
!
! iteration loop
        do ic = 1, 6
            icconv = ic
!
! first guess at pco2c
            call sortin( eyy, pco2y, range1, gammas, ic)
!
! minimum potential rates
! ome (mol/m2/s) Light limited assimilation
! omc (mol/m2/s) Rubisco limited assimilation
! oms (mol/m2/s) sink/export limited assimilation
            omc = vm *((pco2y(ic)-gammas)/(pco2y(ic) + zkc*(1.+po2m/zko))*c3    +  c4)
            ome = par*(pco2y(ic)-gammas)/(pco2y(ic)+2.*gammas)*c3 + par * c4
            sqrtin= MAX( 0.0_dbl_kind, ( (ome+omc)**2._dbl_kind - 4._dbl_kind*sib%param%atheta(i)*ome*omc ) )
            omp  = ( ( ome+omc ) - SQRT( sqrtin ) ) / ( 2.*sib%param%atheta(i) )
            oms  = omss * c3 + omss*pco2y(ic) * c4
            sqrtin= MAX( 0.0_dbl_kind, ( (omp+oms)**2 - 4.*sib%param%btheta(i)*omp*oms ) ) 
            assimy(ic) = (( ( oms+omp ) - SQRT( sqrtin ) ) / ( 2.*sib%param%btheta(i) ))*sib%diag%aparkk
            assimny(ic) = assimy(ic) - sib%diag%resp_can(i) 
!
!pl prognose the new CAS CO2 according to flux divergence
!pl we are going to do this in mol C / mol air (same as PaC/PaAir)
            co2a    = sib%prog%pco2ap /   (sib%drvr%ps*100.)

            co2m    = sib%drvr%pco2m  /   (sib%drvr%ps*100.) 

            co2a=(co2a+(dtt/co2cap)*(sib%diag%resp_grnd-assimny(ic)+co2m*gah2o))/(1+dtt*gah2o/co2cap)
            
          pco2a = co2a * sib%drvr%ps * 100.
!

!itb...intermediate leaf surface CO2
            pco2s = pco2a - (1.4/gbh2o * assimny(ic) * sib%drvr%ps*100.)
!
!itb...intermediate leaf internal CO2
            pco2i = pco2s - assimny(ic) * sib%drvr%ps * 100.0 * 1.6/gsh2o
!
!itb...intermediate leaf chloroplast CO2-this is what we iterate on 
            pco2c = pco2i - assimny(ic) * sib%drvr%ps * 100.0 * 1.0/xgco2m*c3
          sib%diag%testvar1=omc*1.e6
          sib%diag%testvar2=ome*1.e6
          sib%diag%testvar3=oms*1.e6
!
! error convergence
            eyy(ic) = pco2y(ic) - pco2c
          if(abs(eyy(ic))<0.1) exit
        enddo ! iteration loop
!
! save physiology-specific values
        sib%diag%pco2c(i)  = pco2y(icconv)
        sib%diag%assimn(i) = assimny(icconv)
        sib%diag%assim(i)  = assimy(icconv)
        sib%diag%pco2i(i)  = sib%diag%pco2c(i) + sib%diag%assimn(i)/xgco2m*sib%drvr%ps*100.0*c3
        sib%diag%pco2s(i)  = sib%diag%pco2i(i) + sib%diag%assimn(i)/gsh2o *sib%drvr%ps*100.0
!
!  update stomatal resistance...
        h2oi   = sib_loc%etc / sib%drvr%ps
        h2oa   =  sib%prog%ea / sib%drvr%ps
        ecmole = 55.56 * sib%diag%ecmass * dti 
        h2os = h2oa + ecmole / gbh2o
        h2os  = min( h2os, h2oi )
        h2os  = max( h2os, 1.0e-7_dbl_kind)

        h2osrh = h2os / h2oi
        sib%diag%rstfac(1) = h2os/h2oi

        !Bio relaxed this condition to 1/10 of previous (.05 vs .5). The old way made
        !Bio the CO2 on top of the leaves always at least 1/2 of the value at the
        !Bio reference level.
        co2s = MAX(sib%diag%pco2s(i),sib%drvr%pco2m*0.05) / (sib%drvr%ps*100.)
!
!Bio Ball-Berry stromatal conductance equation     
        gsh2oinf = (sib%param%gradm(i) *  MAX(1.0e-12_dbl_kind,sib%diag%assimn(i))  &
            * h2osrh * soilfrz / co2s) + bintc

        !Bio this is the change in stomatal resistance
        !itb...this has been brought here from ADDINC 

        drst(i) = sib%prog%rst(i) * qdamp * ((gsh2o-gsh2oinf)/(pdamp*gsh2o+qdamp*gsh2oinf))

        bintc = bintc * sib%prog%tc / ( 44.032476 * tprcor)

        sib%prog%rst(i) = sib%prog%rst(i) + drst(i)

        ! bintc(i)- smallest canopy stomatal conductance needs to be passed in here.
        ! ---- c.zhang, 2/3/93

        sib%prog%rst(i)=MIN( 1./bintc, sib%prog%rst(i) )

        !...leaf conductance...
        sib%diag%ggl(i) = 1.0 / (sib%prog%rst(i)* sib%diag%rc)

        !
        !itb...determine weighted mean values for the gridcell
        !
        sib%diag%pco2c(6)   = sib%diag%pco2c(6)  + sib%diag%pco2c(i)  *  sib%param%physfrac(i)
        sib%diag%pco2i(6)   = sib%diag%pco2i(6)  + sib%diag%pco2i(i)  *  sib%param%physfrac(i)
        sib%diag%pco2s(6)   = sib%diag%pco2s(6)  + sib%diag%pco2s(i)  *  sib%param%physfrac(i)
        sib%diag%assimn(6)  = sib%diag%assimn(6) + sib%diag%assimn(i) *  sib%param%physfrac(i)
        sib%diag%assim(6)   = sib%diag%assim(6)  + sib%diag%assim(i)  *  sib%param%physfrac(i)
        rstfac3(6)     = rstfac3(6)    + rstfac3(i)    * sib%param%physfrac(i) 
        sib%prog%rst(6)     = sib%prog%rst(6)    + sib%prog%rst(i)    *  sib%param%physfrac(i)
        sib%diag%ggl(6)     = sib%diag%ggl(6)    + sib%diag%ggl(i)    *  sib%param%physfrac(i)
        sib%diag%resp_can(6)= sib%diag%resp_can(6)+ sib%diag%resp_can(i)*  sib%param%physfrac(i)

        !...CFRAX...
        !...at the end of each physiology loop, call the CFRAX code
        call cfrax_physloop(sib,i,c3)
    enddo phys_loop   ! PHYSIOLOGY LOOP

    !pl now do the real C_A forecast with the iterated fluxes.

    co2a    = sib%prog%pco2ap /   (sib%drvr%ps*100.)
    co2m    = sib%drvr%pco2m  /   (sib%drvr%ps*100.) 

    !itb...carbon flux between CAS and reference level (mol C m^-2 sec^-1)
    sib%diag%cflux = gah2o*(co2a-co2m)

    sib%prog%expand=sib%prog%pco2ap*sib%diag%cas_cap_co2/rstar/sib%prog%ta- sib%prog%cas_old
  
! original semi-implicit time differencing for CAS CO2
    co2a=(co2a+(dtt/co2cap)*(sib%diag%resp_grnd-sib%diag%assimn(6)+co2m*gah2o))/(1+dtt*gah2o/co2cap)
    sib%prog%pco2ap = co2a * sib%drvr%ps * 100.
    sib%diag%cflux = gah2o*(co2a-co2m)
!
! moles per m2 CO2 in canopy air space
    sib%prog%cas=sib%prog%pco2ap*sib%diag%cas_cap_co2/rstar/sib%prog%ta

    sib%diag%rstfac(3) = rstfac3(6)
    sib%diag%rstfac(4) = sib%diag%rstfac(1) * sib%diag%rstfac(2) *  sib%diag%rstfac(3)

    !...CFRAX...
    !...one last call to calculate canopy-mean discrimination values
    call cfrax_final(sib)
end subroutine phosib
