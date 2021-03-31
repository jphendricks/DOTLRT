!=========================================================================
subroutine ocean_refl(freq, temp, salinity, theta, ref_vert, ref_hor)
!=========================================================================
! calculates ocean surface reflectivity in the Fresnel approximation.
! Assumes specular surface with reflectivity determined from Fresnel coefficients.
!
! history:
!   10/5/20 Kevin Schaefer added implicit none
!   10/5/20 Kevin Schaefer commented code, removed dead code
!   10/5/20 Kevin Schaefer removed 'external d3lec, fresnel_refl'
!------------------------------------------------------------------------
implicit none
!
! input vriables
real(8), intent(in) :: freq      ! (Ghz) frequency
real(8), intent(in) :: temp      ! (K) surface temperature
real(8), intent(in) :: salinity  ! (?) salinity fraction (avg = 0.035)
real(8), intent(in) :: theta     ! (deg) zenith angle, degree from normal incidence
!
! output variables
real(8), intent(out) :: ref_vert ! (-) reflectance vertical polarized
real(8), intent(out) :: ref_hor  ! (-) reflectance horizontal polarized
!
! local variables
integer(4) :: type = 1           ! (-) water type flag 0=NaCl solution 1=Seawater 2=pure water
double complex kappa             ! (-) complex dielectric constant of water
double complex mu                ! (-) real component of complex number

! real component
mu = dcmplx(1.0d0, 0.0d0)

! calculate dielectric constant of water
call d3lec(kappa, freq, temp, salinity, type)

! calculate reflectivity
call fresnel_refl(ref_vert, ref_hor, kappa, mu, theta)

return
end subroutine ocean_refl
