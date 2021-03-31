! Computes the power reflection coefficients (rv and rh) for a
! specular dielectric surface.
! Theta is in degrees measured with respect to the surface normal.
! Note: eps_i,mu_i are assumed negative for lossy media.

subroutine fresnel_refl( rv, rh, eps, mu, theta )
    implicit none
    real(8) :: pi = 3.1415926536d0
    real(8) theta
    real(8) rv, rh
    real(8) c
    double complex eps, mu
    double complex gamv, gamh, term1, term2
    c = dcos(theta*pi/180.0d0) !dcosd(theta)
    term1 = cdsqrt(eps*mu+c*c-1)
    term2 = eps*c
    gamv = (term2-term1)/(term2+term1)
    term2 = mu*c
    gamh = (term2-term1)/(term2+term1)
    rv = gamv * dconjg(gamv)
    rh = gamh * dconjg(gamh)
end subroutine fresnel_refl
