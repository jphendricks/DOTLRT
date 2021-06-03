!
!=======================================================================
      subroutine smoothmin (Atheta,Btheta, &
        omegaC,omegaE,omegaS,assim)
!=======================================================================      
! calculates smoothed minimum of PAR, Rubisco, & export assimilation rates
! or any associated scaling factor or variable.
!-----------------------------------------------------------------------
!
      implicit none
!
! input variables
      real Atheta ! (-) Coupling factor between Rubisco-light assimilation rates
      real Btheta ! (-) Coupling factor between Rubisco-light-export assim rates
      real omegaC ! (mol/m2/s) Rubisco limted assimilation rate
      real omegaE ! (mol/m2/s) Light limited assimilation rate
      real omegaS ! (mol/m2/s) export (C3) or PEP (C4) limited assimilation rate
      real omegaP ! (mol/m2/s) smoothed minimum rubisco & light rates
!
! output variables
      real assim  ! (mol/m2/s) smoothed minimum rubisco, light, & export rates
!
! internal variables
      real sqrtin ! root portion of smoothing quadratic
!
! Interpolate between light and Rubisco rates
      sqrtin=MAX(0.,((omegaE+omegaC)**2.-4.*Atheta*omegaE*omegaC))
      omegaP=((omegaE+omegaC)-SQRT(sqrtin))/(2.*Atheta)
!
! Interpolate between light, Rubisco, and export rates
      sqrtin=MAX(0.,((omegaP+omegaS)**2.-4.*Btheta*omegaP*omegaS))
      assim=((omegaS+omegaP)-SQRT(sqrtin))/(2.*Btheta)
!
      return
      end
