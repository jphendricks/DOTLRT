!
!==================================================================
subroutine fresnel_refl( iang )
!==================================================================
! Computes the power reflection coefficients (rv and rh) for a
! specular dielectric surface.
! Theta is in degrees measured with respect to the surface normal.
! Note: kappa,mu are assumed negative for lossy media.
!
! History:
!  10/5/20 Kevin Schaefer commented the code
!  10/16/2020 Kevin Schaefer switched to variables module
!-----------------------------------------------------------------
use dotlrt_variables
implicit none
!
! inputs
  integer(4) iang   ! angle index
!
! internal
  real(8) costheta ! cosine of theta
  double complex mu
  double complex gamv
  double complex gamh
  double complex term1
  double complex term2

! the calculations
  mu = dcmplx(1.0d0, 0.0d0)
  costheta = dcos(surf_inp%theta(iang)*pi/180.0d0) 
  term1 = cdsqrt(surf_inp%dielectric*mu+costheta*costheta-1)
  term2 = surf_inp%dielectric*costheta
  gamv = (term2-term1)/(term2+term1)
  term2 = mu*costheta
  gamh = (term2-term1)/(term2+term1)
  surf_inp%vr(iang) = gamv * dconjg(gamv)
  surf_inp%hr(iang) = gamh * dconjg(gamh)

end subroutine fresnel_refl
