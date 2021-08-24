!
!=====================SUBROUTINE CYCALC=================================
      SUBROUTINE CYCALC_Suits( APARKK, VM, ATHETA, BTHETA, par, &
                        GAMMAS, RESPC, RRKK, OMSS, C3, C4,&
                        PCO2I, ASSIMN, assim, len,test, test10)
!=======================================================================
!     CALCULATION EQUIVALENT TO STEPS IN FIGURE 4 OF SE-92A
!     C4 CALCULATION BASED ON CO-92.
! OUTPUT
!       PCO2I          CANOPY INTERNAL CO2 CONCENTRATION (MOL MOL-1)
!       GSH2O          CANOPY CONDUCTANCE (MOL M-2 S-1)
!       H2OS           CANOPY SURFACE H2O CONCENTRATION (MOL MOL-1)
! DIAGNOSTICS
!       OMC            RUBISCO LIMITED ASSIMILATION (MOL M-2 S-1)
!       OME            LIGHT LIMITED ASSIMILATION (MOL M-2 S-1)
!       OMS            SINK LIMITED ASSIMILATION (MOL M-2 S-1)
!       CO2S           CANOPY SURFACE CO2 CONCENTRATION (MOL MOL-1)
!-----------------------------------------------------------------------
      implicit none

      integer len
      real aparkk(len),vm(len),atheta(len), &
          btheta(len),gammas(len),par(len), &
          respc(len),rrkk(len),omss(len), &
          c3(len),c4(len),pco2i(len), &
          assimn(len),assim(len)

!    local variables
      real ome, omc, omp, oms, sqrtin
      integer i
      real test, test10(10)   ! test variables
 
!-----------------------------------------------------------------------
!     CALCULATE ASSIMILATION RATE
!      OMC         (OMEGA-C): EQUATION (11) , SE-92A
!      OME         (OMEGA-E): EQUATION (12) , SE-92A
!      OMS         (OMEGA-S): EQUATION (13) , SE-92A
!      ASSIMN      (A-N)    : EQUATION (14,15), SE-92A
!-----------------------------------------------------------------------

      do i = 1,len
        OMC = VM(i) *(PCO2I(i)-GAMMAS(i))/(PCO2I(i) + RRKK(i))*C3(i) + VM(i) * C4(i)
        OME = PAR(i)*(PCO2I(i)-GAMMAS(i))/(PCO2I(i)+2.*GAMMAS(i))*C3(i) + PAR(i) * C4(i)
        SQRTIN= MAX(0., ( (OME+OMC)**2 - 4.*ATHETA(i)*OME*OMC ) )
        OMP  = ( ( OME+OMC ) - SQRT( SQRTIN ) ) / ( 2.*ATHETA(i) )
        OMS  = OMSS(i) * C3(i) + OMSS(i)*PCO2I(i) * C4(i)
        SQRTIN= MAX(0., ( (OMP+OMS)**2 - 4.*BTHETA(i)*OMP*OMS ) )
        ASSIM(i) = ( ( OMS+OMP ) - SQRT( SQRTIN ) ) /( 2.*BTHETA(i) )
        ASSIMN(i)= ( ASSIM(i) - RESPC(i)) * APARKK(i)
      enddo
      test=OMP
      test10(1)=OMP
!
      RETURN
      END
