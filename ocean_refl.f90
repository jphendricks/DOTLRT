! Subroutine ocean_refl calculates ocean surface reflectivity in the Fresnel approximation.
! Assumes specular surface with reflectivity determined from Fresnel coefficients.
!     freq : GHz
!     temp : K
!     sal  : salinity fraction (avg = 0.035)
!     theta: deg from normal incidence
subroutine ocean_refl(freq, temp, sal, theta, rv, rh)
    real(8), intent(in)  :: freq, temp, sal, theta
    real(8), intent(out) :: rv, rh
    integer(4) :: i = 1
    double complex kappa, mu
    external d3lec, fresnel_refl
    mu = dcmplx(1.0d0, 0.0d0)
    call d3lec(kappa,freq,temp,sal,i)
    ! dilec2(kappa,freq,temp,sal,i)
    ! call fresnel_refl(rv,rh,kappa,cmplx(1.0d0,0.0d0),theta)
    call fresnel_refl(rv,rh,kappa,mu,theta)
    ! ref(rh,rv,kappa,theta)
    return
end subroutine ocean_refl
