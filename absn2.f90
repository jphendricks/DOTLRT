! Real function absn2
! INPUTS:
! temperature = TEMPERATURE (K)
! pressure  = PRESSURE (MB)
! frequency  = FREQUENCY (GHZ)
! OUTPUTS:
!    absn2   = ABSORPTION COEFFICIENT DUE TO NITROGEN IN AIR  (NEPER/KM)
!   dabsn2_t = DERIVATIVE OF absn2 WITH RESPECT TO TEMPERATURE
!   dabsn2_p = DERIVATIVE OF absn2 WITH RESPECT TO PRESSURE
! The derivatives were programmed and tested by R. Hill Aug. 2003
real(8) function absn2(temperature, pressure, frequency, dabsn2_t, dabsn2_p)
  implicit none
  real(8) temperature, pressure, frequency, dabsn2_t, dabsn2_p
  real(8) Th

  Th = 300.0d0 / temperature
  absn2 = 5.87d-14 * pressure * pressure * frequency * frequency * (Th**4.5d0)
  dabsn2_t = -absn2*4.5d0/temperature
  dabsn2_p = absn2*2.d0/pressure

  return
end function absn2
