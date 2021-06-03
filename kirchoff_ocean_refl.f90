!=======================================================================================
subroutine kirchoff_ocean_refl(iang, freq, slope_var)
!=======================================================================================
! This subroutine calculates the vertical and horizontal reflectivities (rv and rh) of
! the ocean surface using the Kirchoff approximation.
! Assumes geometrically optical facets with Gaussian distributed slopes and
! isotropic azimuthal orientations. If the slope variance is less than zero, 
! a specular ocean surface is assumed. Shadowing is neglected.
!     freq      : GHz
!     temp      : K
!     sal       : salinity fraction (avg = 0.035)
!     theta_i   : degrees from average normal incidence
!     slope_var : variance of facet slopes
!
! history:
!  10/16/20 Kevin Schaefer switched to variable module and replaced duplicate variables
!---------------------------------------------------------------------------------------
use dotlrt_variables
    implicit none
    real(8), intent(in)  :: freq, slope_var
  integer(4) iang   ! angle index
    integer(4) :: az_max = 25 ! Number of azimuthal panels in integration
    integer(4) :: el_max = 20 ! Number of elevation panels in integration
    integer(4) :: slope_stds = 10 ! Number of slope standard deviations to integrate over
    integer(4) :: j, k
    real(8)    :: rv_k,rh_k, phi_k, theta_k, d_phi_k, d_theta_k, sin_theta_i,  &
                 cos_theta_i, sin_theta_k, cos_theta_k, sin_phi_k, cos_phi_k, &
                 cos_theta_d_k, theta_k_max, norm,sec_exp,two_pi, int_v,int_h, &
                 int_d, v_dot_q_s, v_dot_p_s
    double complex :: gamv, gamh, term1, term2
    surf%vref(iang) = 0.0d0
    surf%href(iang) = 0.0d0                                                                           
    if( freq      .ge.   1.0d0 .and. freq      .le. 1000.0d0 .and. &
        surf%temp      .ge. 200.0d0 .and. surf%temp      .le.  350.0d0 .and. &
        surf%sal       .ge.   0.0d0 .and. surf%sal       .le.    1.0d0 .and. &
        surf%theta(iang)   .ge.   0.0d0 .and. surf%theta(iang)   .le.   90.0d0 .and. &
        slope_var .ge.   0.0d0 .and. slope_var .le.    0.5d0 ) then
        if( slope_var .ge. 1.0d-4 ) then
            ! facetted ocean surface
            two_pi = 2.0d0 * pi
            sin_theta_i = dsin(surf%theta(iang)*pi/180.0d0)
            cos_theta_i = dcos(surf%theta(iang)*pi/180.0d0)
            d_phi_k = two_pi/dble(az_max)
            theta_k_max = datan(slope_stds*dsqrt(slope_var))
            d_theta_k = theta_k_max/dble(el_max)
            surf%vref(iang) = 0.0d0
            surf%href(iang) = 0.0d0                                                                           
            norm = 0.0d0
            ! Integrate over theta: use trapezoidal rule out to "slope_stds" slope standard deviations}
            do j = 0, el_max
                theta_k = j*d_theta_k
                sin_theta_k = dsin(theta_k)
                cos_theta_k = dcos(theta_k)
                ! Integrate over phi: use trapezoidal rule and periodicity of integrand}
                int_v = 0.0d0
                int_h = 0.0d0
                int_d = 0.0d0
                do k = 0, (az_max-1)
                    phi_k = k*d_phi_k
                    sin_phi_k = dsin(phi_k)
                    cos_phi_k = dcos(phi_k)
                    cos_theta_d_k = sin_theta_i * sin_theta_k * cos_phi_k + cos_theta_i * cos_theta_k
                    if( cos_theta_d_k .gt. 0.0d0 ) then
                        term1 = surf%diel + cos_theta_d_k * cos_theta_d_k -1.0d0
                        term1 = cdsqrt(term1)
                        term2 = cos_theta_d_k
                        gamh = (term2-term1) / (term2+term1)
                        term2 = surf%diel * term2
                        gamv = (term2-term1) / (term2+term1)
                        rv_k = dble(gamv * conjg(gamv))
                        rh_k = dble(gamh * conjg(gamh))
                        if( theta_k .gt. 0.0d0 .and. phi_k .gt. 0.0d0 .and. phi_k .lt. two_pi ) then
                             v_dot_q_s = 1.0d0/(1.0d0+((sin_theta_i*cos_theta_k-cos_theta_i*sin_theta_k*cos_phi_k) &
                                         /(sin_theta_k*sin_phi_k))**2)
                        else
                             v_dot_q_s = 0.0d0
                        end if
                        v_dot_p_s = 1.0d0-v_dot_q_s
                        int_v = int_v+(v_dot_p_s*rv_k+v_dot_q_s*rh_k)*cos_theta_d_k
                        int_h = int_h+(v_dot_q_s*rv_k+v_dot_p_s*rh_k)*cos_theta_d_k
                        int_d = int_d+cos_theta_d_k
                    end if
                end do
                sec_exp = dexp(-(sin_theta_k/cos_theta_k)**2/(2.0d0*slope_var))/(cos_theta_k)**2
                if( j .eq. 0 .or. j .eq. el_max ) THEN
                    surf%vref(iang) = surf%vref(iang)+int_v*sec_exp/2.0d0
                    surf%href(iang) = surf%href(iang)+int_h*sec_exp/2.0d0
                    norm = norm+int_d*sec_exp/2.0d0
                else
                    surf%vref(iang) = surf%vref(iang)+int_v*sec_exp
                    surf%href(iang) = surf%href(iang)+int_h*sec_exp
                    norm = norm+int_d*sec_exp
                end if
            end do
            surf%vref(iang) = dble(surf%vref(iang) / norm)
            surf%href(iang) = dble(surf%href(iang) / norm)
        else
            ! specular ocean surface}
            cos_theta_i = dcos(surf%theta(iang)*pi/180.0d0)
            term1 = cdsqrt(surf%diel + cos_theta_i * cos_theta_i - 1.0d0)
            term2 = cos_theta_i
            gamh = (term2-term1) / (term2+term1)
            term2 = surf%diel * term2
            gamv = (term2-term1) / (term2+term1)
            surf%vref(iang) = dble(gamv * conjg(gamv))
            surf%href(iang) = dble(gamh * conjg(gamh))
        end if
    end if
end subroutine kirchoff_ocean_refl
