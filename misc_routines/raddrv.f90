!
!---------------------------------------------------------------------
subroutine raddrv(nsib,swdown,sunang,radvbc,radvdc,radnbc,radndc,test,test10)
!---------------------------------------------------------------------
! radiation radive code to use the downward sw at bottom 
! and the formulation to estimate radvbc,radvdc, radndc, radndc
! Modifications:
!  Kevin Schaefer removed minimum value stemp (swdwn) (5/22/13)
!  Kevin Schaefer removed stemp (local swdown) variable (5/22/13)
!  Kevin Schaefer changed cloud/difrat calc from sunang to localcosz (5/22/13)
!  Kevin Schaefer set cloud to zero if stemp is zero (5/22/13)
!---------------------------------------------------------------------

use kinds
use sib_const_module, only:  &
    cosz_min
implicit none

integer(kind=int_kind) :: nsib, i
real(kind=real_kind) ::  &
    cloud,          &  ! cloud cover fraction
    difrat,         &  ! diffuse fraction
    vnrat              ! visible/nir fraction
real(kind=dbl_kind) ::  &
    swdown(nsib),   &  ! shortwave downwelling radiation
    sunang(nsib),   &  ! cosine of solar zenith angle 
    localcosz          ! local version cosine of solar zenith angle
real(kind=dbl_kind) ::  &
    radvbc(nsib),   &  ! visible direct radiation
    radvdc(nsib),   &  ! visible diffuse radiation
    radnbc(nsib),   &  ! NIR direct radiation
    radndc(nsib)       ! NIR diffuse radiation

real(kind=dbl_kind),parameter :: c1 = 580.
real(kind=dbl_kind),parameter :: c2 = 464.
real(kind=dbl_kind),parameter :: c3 = 499.
real(kind=dbl_kind),parameter :: c4 = 963.
real(kind=dbl_kind),parameter :: c5 = 1160.
real test, test10(10)     ! test variables

    do i=1,nsib
        localcosz = max( 0.001_dbl_kind, sunang(i) )

        cloud = (c5 * localcosz - swdown(i)) / (c4 * localcosz)                   
        cloud = max(cloud,0.)                                                
        cloud = min(cloud,1.)
        if(swdown(i)==0.) cloud =0.

        difrat = 0.0604 / ( localcosz-0.001 + 1.0e-10 ) + 0.0683
        if ( difrat < 0. ) difrat = 0.
        if ( difrat > 1. ) difrat = 1.
        difrat = difrat + ( 1. - difrat ) * cloud

        vnrat = ( c1 - cloud*c2 ) / ( ( c1 - cloud*c3 ) + ( c1 - cloud*c2 ) )

        radvbc(i) = (1.-difrat)*vnrat     *swdown(i)
        radvdc(i) =      difrat*vnrat     *swdown(i)
        radnbc(i) = (1.-difrat)*(1.-vnrat)*swdown(i)
        radndc(i) =      difrat*(1.-vnrat)*swdown(i)
    enddo
    test=swdown(1)
    test10(1)=swdown(1)

end subroutine raddrv
