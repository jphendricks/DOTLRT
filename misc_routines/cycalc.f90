!
!=======================================================================
      SUBROUTINE CYCALC(ASSIM,atheta,btheta,ome, omc, oms)
!=======================================================================
! CYCALC calculates the smoothed minimum of light, rubisco, and export
! assimilation rates
!
! References
!   Sellers et al. (1992a),(1996a)
!   Collatz et al. (1992) for C4 calculations
!
! Modifications
!  Kevin Schaefer commented code (5/17/01)
!  Kevin Schaefer included APARKK in assim calculation (7/19/01)
!  Kevin Schaefer used the C4 respcp value for C4 respc (7/21/01)
!  Kevin Schaefer moved all pCO2i & Temp calculations from Phosib (7/22/01)
!  Kevin Schaefer included C4 high/low Temp stress in export rate (7/22/01)
!-----------------------------------------------------------------------
!
      implicit none
!
! input variables
      real atheta   ! (-) OMC OME coupling parameter (bc input)
      real btheta   ! (-) OMC/OME WS coupling parameter (bc input)
      real ome      ! (mol/m2/s) Light limited assimilation rate
      real omc      ! (mol/m2/s) Rubisco limited assimilation rate
      real oms      ! (mol/m2/s) sink/export limited assimilation rate
      real omp      ! (mol/m2/s) combo ome/omc assimilation rate
!
! output variables
      real assim ! (mol/m2/s) gross canopy assimilation
!
! Internal variables
      real sqrtin        ! (-) root portion of quadratic
!
      real test, test10(10)   ! test variables
!
! smooth between light and Rubisco rates
      SQRTIN=MAX(0.,((OME+OMC)**2-4.*ATHETA*OME*OMC))
      OMP=((OME+OMC)-SQRT(SQRTIN))/(2.*ATHETA)
!
! smooth between light, Rubisco, and export rates
      SQRTIN=MAX(0., ((OMP+OMS)**2-4.*BTHETA*OMP*OMS))
      ASSIM=((OMS+OMP)-SQRT(SQRTIN))/(2.*BTHETA)
      
      test = omp
      test10(1)= omp
!
      RETURN
      END
!
