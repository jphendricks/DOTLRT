!===================================================================
subroutine calcprofile_d()
!===================================================================
! calculates vertical profile of particle size distribution parameters 
! p and q are constant for each phase and are parameters of gamma distribution of hydrometeor size
! k0 and a0 depend only on phase density and not temperature or frequency.
!
! History:
! TBD       Gail Skofronick Jackson wrote Original code 
! 11/1/1998 Marian Klein Modified code for MRT
! 4/1/2003  Ronald Richter Converted code from Marian Klein and Albin Gasiewski from PASCAL to FORTRAN
! 8/1/2003  Reg Hill tested Differentiated and derivatives 
! 8/28/2003 Reg Hill inserted "a0_const" to save computer time in subroutine HYDROMETEOR_MASTER_5PH_D
! 9/26/2020 Kevin Schaefer deleted unused variables and tabs 
! 2/3/2021  Kevin Schaefer cleaned up code, deleted all unused code
!-------------------------------------------------------------------
use dotlrt_variables

! internal variables
  integer ilev   ! (-) level index
  real(8) dk0_dw ! (?) derivative of k0 wrt water density

! derivatives of k0 and a0 with respect to cloud water density (each phase)
! is multipled by cloud water density (each phase) to remove infinities
! when cloud water density (each phase) goes to zero.

! Cloud liquid water
  do ilev= 1, nlev
    if( atm(ilev)%clw%dens .ne. 0.0d0 ) then
      ! Marshall Palmer distribution
      atm(ilev)%clw%p  = 0.0d0
      atm(ilev)%clw%q  = 1.0d0
      ! Constant a0, k0 varies w/r/t M
      dk0_dw = 2.0d0 * 1989436788.65d0
      atm(ilev)%clw%dk0_dw = dk0_dw
      atm(ilev)%clw%k0 = dk0_dw * atm(ilev)%clw%dens
      atm(ilev)%clw%a0 = 0.01d0
      atm(ilev)%clw%da0_dw = 0.0d0
      atm(ilev)%clw%a0_const = .true.
    else
      atm(ilev)%clw%p  = 0.0d0
      atm(ilev)%clw%q  = 0.0d0
      atm(ilev)%clw%k0 = 0.0d0
      atm(ilev)%clw%a0 = 0.0d0
      atm(ilev)%clw%dk0_dw = 0.0d0
      atm(ilev)%clw%da0_dw = 0.0d0
      atm(ilev)%clw%a0_const = .true.
    end if

! rain
    if( atm(ilev)%rain%dens .ne.  0.0d0 ) then
      ! Marshall Palmer distribution
      atm(ilev)%rain%p  = 0.0d0
      atm(ilev)%rain%q  = 1.0d0
      atm(ilev)%rain%k0 = 16000.0d0 
      atm(ilev)%rain%a0 = 0.22331d0 * ( atm(ilev)%rain%dens**0.250d0 )
      if( atm(ilev)%rain%dens < 1.d0-8 ) then
        atm(ilev)%rain%da0_dw = 0.0d0
      else
        atm(ilev)%rain%da0_dw = 0.250d0 * (atm(ilev)%rain%a0)/(atm(ilev)%rain%dens)
      end if
      atm(ilev)%rain%dk0_dw = 0.0d0
      atm(ilev)%rain%a0_const = .false.
    else
      atm(ilev)%rain%p  = 0.0d0
      atm(ilev)%rain%q  = 0.0d0
      atm(ilev)%rain%k0 = 0.0d0
      atm(ilev)%rain%a0 = 0.0d0
      atm(ilev)%rain%dk0_dw = 0.0d0
      atm(ilev)%rain%da0_dw = 0.0d0
      atm(ilev)%rain%a0_const = .true.
    end if

! ice
    if( atm(ilev)%ice%dens .ne. 0.0d0 ) then
      ! Sekhon-Sriv. Distribution, for two phase ice
      atm(ilev)%ice%p  = 0.0d0
      atm(ilev)%ice%q  = 1.0d0
      ! Constant a0, k0 varies w/r/t M, assumes no size change when changing phase
      ! (see Bauer and Schluessel) - allows for polydispersive sizes and for mie
      ! and Rayleigh scattering
      dk0_dw = 2.0d0 * 1989436788.65d0
      atm(ilev)%ice%dk0_dw = dk0_dw
      atm(ilev)%ice%k0 = dk0_dw * atm(ilev)%ice%dens
      atm(ilev)%ice%a0 =  0.01d0
      atm(ilev)%ice%da0_dw = 0.d0
      atm(ilev)%ice%a0_const = .true.
    else
      atm(ilev)%ice%p  = 0.0d0
      atm(ilev)%ice%q  = 0.0d0
      atm(ilev)%ice%k0 = 0.0d0
      atm(ilev)%ice%a0 = 0.0d0
      atm(ilev)%ice%dk0_dw = 0.0d0
      atm(ilev)%ice%da0_dw = 0.0d0
      atm(ilev)%ice%a0_const = .true.
    end if

! snow
    if( atm(ilev)%snow%dens .ne. 0.0d0 ) then
      ! Tao, Prasad, Alder snow size distributions converted to a0 and k0 form
      ! (see ntbk#3 pg 111)
      atm(ilev)%snow%p  = 0.0d0
      atm(ilev)%snow%q  = 1.0d0
      ! If following Rutledge and Hobbs (1983), then parameterization: cloudsnowk0 = 20000.0
      ! Instead, use Rutledge and Hobbs (1984)
      atm(ilev)%snow%k0 = 8000.0d0
      atm(ilev)%snow%dk0_dw = 0.d0
      atm(ilev)%snow%a0 = 0.3757d0 * ( atm(ilev)%snow%dens**0.25d0 )
      atm(ilev)%snow%da0_dw=0.25d0*(atm(ilev)%snow%a0)/(atm(ilev)%snow%dens)
      atm(ilev)%snow%a0_const = .false.
    else
      atm(ilev)%snow%p  = 0.0d0
      atm(ilev)%snow%q  = 0.0d0
      atm(ilev)%snow%k0 = 0.0d0
      atm(ilev)%snow%a0 = 0.0d0
      atm(ilev)%snow%dk0_dw = 0.0d0
      atm(ilev)%snow%da0_dw = 0.0d0
      atm(ilev)%snow%a0_const = .true.
    end if

! graupel
    if( atm(ilev)%grpl%dens .ne. 0.0d0 ) then
      ! Tao, Prasad, Alder graupel size distributions converted to a0 and k0 form (see ntbk#3 pg 111)
      atm(ilev)%grpl%p  = 0.0d0
      atm(ilev)%grpl%q  = 1.0d0
      atm(ilev)%grpl%k0 = 8000.0d0 
      atm(ilev)%grpl%dk0_dw = 0.0d0
      atm(ilev)%grpl%a0 = 0.3340d0 * ( atm(ilev)%grpl%dens**0.25d0 )
      atm(ilev)%grpl%da0_dw=0.25d0*(atm(ilev)%grpl%a0)/(atm(ilev)%grpl%dens)
      atm(ilev)%grpl%a0_const = .false.
    else
      atm(ilev)%grpl%p  = 0.0d0
      atm(ilev)%grpl%q  = 0.0d0
      atm(ilev)%grpl%k0 = 0.0d0
      atm(ilev)%grpl%a0 = 0.0d0
      atm(ilev)%grpl%dk0_dw = 0.0d0
      atm(ilev)%grpl%da0_dw = 0.0d0
      atm(ilev)%grpl%a0_const = .true.
    end if
  end do 

return
end subroutine calcprofile_d
