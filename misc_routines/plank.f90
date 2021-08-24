!
!=======================================================================
      subroutine Plank_wavelength
!=======================================================================      
! Calculates plancks function at specified wavelength & wavelength interval
!----------------------------------------------------------------------
!
      use sibscan_Variables
!
      implicit none
!
      real demon ! denominator of planck's function
!
! calculate denominator of planck's function
      demon=exp(0.014403/Tplank/wave)-1.
!
! calculate planck's function
      black=pi*1.19e-16/wave**5./demon*Dwave
!
      return
      end
!
!=======================================================================
      subroutine Plank_frequency
!=======================================================================      
! Calculates plancks function at specified frequency
!----------------------------------------------------------------------
!
      use sibscan_Variables
!
      implicit none
!
      real temp ! temporary variable
      real MegaHz ! input frequency converted to megahertz
!
      MegaHz=freq*1.e-9
!
      temp=exp(plank*freq/Boltz/Tplank)-1.

      black=1.473e-23*MegaHz*MegaHz*MegaHz/temp
!
      return
      end
!
