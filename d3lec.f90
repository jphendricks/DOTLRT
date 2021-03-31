!
!===============================================================================
subroutine dielectric_const_water( freq )
!===============================================================================
! calculates the complex dielectric constant of water
! as a function of temperature, salinity, and frequency
! complex dielectric constant of water from Ray, Applied Optics Vol 11, No 8
!
! History:
!   4/1/1997 M. Klein Rewrote to Delphi Pascal
!   10/5/20 Kevin Schaefer commented code
!  10/16/20 Kevin Schaefer switched to variable module and replaced duplicate variables
!------------------------------------------------------------------------------
Use dotlrt_variables
implicit none

! input variables
real(8) freq         ! (GHz) frequency < c/1000 microns
!
! really poorly named internal variables
real(8) L
real(8) LS
real(8) N
real(8) T
real(8) SIGMA
real(8) EI
real(8) ES
real(8) ALPHA
real(8) S   ! (?) salinity fraction in different units
real(8) DEL
real(8) SIG
real(8) D
real(8) ep
real(8) epp

! the calculations start assuming freshwater
  T = surf_inp%surf_temp - 273.2d0
  SIGMA = 12.5664d8
  EI = 5.27137d0 + T * ( 0.0216474d0 - 0.00131198d0 * T )
  T = T - 25.0d0
  ES = 78.54d0 * (1.0d0 + T * ( -4.579d-3 + T * ( 1.19d-5 - 2.8d-8 * T ) ) )
  ALPHA = 0.0609265d0 - 16.8129d0 / surf_inp%surf_temp
  LS = 0.33836d-5 * dexp( 2513.98d0 / surf_inp%surf_temp)
  

  if( trim(surf_inp%diel_mod)/='Fresh') then ! salt water
    S = surf_inp%salinity * 1000.0d0
    N = S * (1.707d-2 + S * (1.205d-5 + 4.058d-9 * S))
    DEL = -T
    SIG = 0.0d0

    if( trim(surf_inp%diel_mod)/='Salt') then ! seawater
        SIG = S * ( 0.182521d0 + S * ( -1.46192d-3        &
            + S * ( 2.09324d-5 - 1.28205d-7 * S ) ) )     &
            * dexp( - DEL * ( 2.033d-2 + DEL * ( 1.266d-4 &
            + 2.464d-6 * DEL ) - S * ( 1.849d-5           &
            + DEL * ( -2.551d-7 + 2.551d-8 * DEL ) ) ) )
    endif

    if( trim(surf_inp%diel_mod)=='NaCl' ) then ! NaCl water
        SIG = N * ( 10.394d0 + N * ( -2.3776d0 &
            + N * ( 0.68258d0 + N * ( -0.13538d0 &
            + N * 1.0086d-2 ) ) ) )              &
            * ( 1.0d0 + DEL * ( -1.962d-2 + DEL &
            * 8.08d-5 ) - DEL * N * ( 3.020d-5 &
            + 3.922d-5 * DEL + N * ( 1.721d-5 - 6.584d-6 * DEL ) ) )
    endif
    SIGMA = SIG / 8.854d-12
    ES = ES * ( 1.0d0 + N * ( -0.2551d0 + N * ( 5.151d-2 - 6.889d-3 * N ) ) )
    LS = LS * ( 1.0d0 + N * ( T * 0.1463d-2 - 0.04896d0 + N * ( -0.02967d0 + N * 5.644d-3 ) ) )
  end if

  S = dsin(ALPHA*90.0d0*pi/180.0d0)  !dSINd(ALPHA*90.0d0)
  L = 2.99776d-1/FREQ
  T = (LS/L)**(1.0d0-ALPHA)
  D = 1.0d0+T*(2.0d0*S+T)
  EP = EI+(ES-EI)*(1.0d0+T*S)/D
  EPP = (ES-EI)*T*dcos(ALPHA*90.0d0*pi/180.0d0)/D+L*SIGMA/18.8496d8 !(ES-EI)*T*dcosd(ALPHA*90.0d0)/D+L*SIGMA/18.8496d8
  surf_inp%dielectric = dcmplx(ep,-epp)

end subroutine dielectric_const_water
