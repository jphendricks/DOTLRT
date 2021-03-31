! calc_cristoffel_weights.f90
!
! Purpose: Calculates Cristoffel weights and angles for Gauss-Lobotto 
!          quadrature over theinterval[-1..1]

! Converted from original Pascal by:
!
! Bill Otto NOAA/ETL SET June 2003  william.d.otto@noaa.gov
! Modified by Ron Richter 7/14/03  - removed num_quad_angles from arguments because
! it is in variables_unit,  pi (PI) is in variables_unit also so removed local definition

subroutine calc_cristoffel_weights(num_quad_angles) ! cris_quad_wghts)
    use variables
    real(8) p1, p2, p3, pp, z, z1, sum
!    real(8) cris_quad_wghts(100)
    integer  i, j
!    sum = 0.0d0
    do i = 1, (num_quad_angles+1)/2
        z = dcos(pi*(i-0.25d0)/(num_quad_angles+0.5d0))
        z1 = 2.0d0
        do while( dabs(z-z1) .gt. 3.0d-14)
            p1 =1.0d0
            p2 =0.0d0
            do j = 1, num_quad_angles 
                p3 = p2
                p2 = p1
                p1 = ((2*j-1)*z*p2-(j-1)*p3)/dble(j)
            end do
            pp = num_quad_angles*(z*p1-p2)/(z*z-1.0d0)
            z1 = z
            z = z1-p1/pp
        end do
        quad_angle_array(i) = (180.0d0 / pi) * dacos(z)
        quad_angle_array(num_quad_angles+1-i) = 180.0d0-quad_angle_array(i)
        cris_quad_wghts(i) = 2.0d0/((1.0d0-z*z)*pp*pp)
        cris_quad_wghts(num_quad_angles+1-i) = cris_quad_wghts(i)
!        sum = sum + cris_quad_wghts(i) + cris_quad_wghts(num_quad_angles+1-i)
    end do
    return
end ! calc_cristoffel_weights


!subroutine calc_cristoffel_weights()
! Calculates Cristoffel weights and angles for Gauss-Lobotto quadrature over the interval[-1..1]

!double precision p1, p2, p3, pp, z, z1
!integer i,j

!  do i = 1, trunc((num_quad_angles+1)/2)
!    z = dcos(pi*(i-0.25)/(num_quad_angles+0.5))
!    z1 = 2.0d0
!    do while( dabs(z-z1) > 3.0d-14 )
!      p1 = 1.0d0
!      p2 = 0.0d0
!      do j = 1, num_quad_angles
!        p3 = p2
!        p2 = p1
!        p1 = ( ( 2.0d0 * dble(j) - 1.0d0 ) * z * p2 - ( dble(j) - 1.0d0 ) * p3 ) / dble(j)
!      end do
!      pp = num_quad_angles * ( z * p1 - p2 ) / ( z * z - 1.0d0 )
!      z1 = z
!      z = z1 - p1 / pp
!    end do while
!    quad_angle_array(i)                   = ( 180.0d0 / pi ) * dacos(z)
!    quad_angle_array(num_quad_angles+1-i) = 180.0d0 - quad_angle_array(i)
!    cris_quad_wghts(i)                    = 2.0d0 / ( ( 1.0d0 - z * z ) * pp * pp )
!    cris_quad_wghts(num_quad_angles+1-i)  = cris_quad_wghts(i)
!  end do ! i
!end subroutine calc_cristoffel_weights
