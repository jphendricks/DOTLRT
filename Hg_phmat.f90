!=======================================================
SUBROUTINE HG_phmat()
!=======================================================
! The Henyey-Greenstein (HG) scattering phase matrix is abreviated as HGph.
! This subroutine calculates HG phase matrix for a set of quadrature angles,
!       HGph(k,i,j)=(1-g^2)/(2*pi)*INT(from x=0 to x=pi) dx
!       [1+g^2+2*g*cos(t1)*cos(t2)+2*g*sin(t1)*sin(t2)*cos(x)]^(-3/2)
! where t1=quad_ang(i)*pi/180, t2=quad_ang(j)*pi/180, and g=HGg(k)
! It also calculates a derivative of the phase matrix with respect 
! to parameter _g:
!               dHGph(k,i,j)= d[HGph(k,i,j)]/dg
!
! -------------------   INPUT ------------------------
!    nhg          ( I)    number of HG parameters to be calculated
!  HGg(ng)       (DP)    array of values of parameter _g
! -------------------  OUTPUT ------------------------
!  HGph(ng,nstream,nstream)  (DP)  array of HG phase matrix
! dHGph(ng,nstream,nstream)  (DP)  array of the appropriate derivatives
!                              of HG phase matrix
!-------------------------------------------------------
! History
!   9/26/202 Kevin Schaefer deleted unused variables
!  10/17/2020 Kevin Schaefer moved quadrature angle assignment to configure
!  12/12/2020 Kevin Schaefer moved nhg assignment to dotlrt_variables
!  12/13/2020 Kevin Schaefer moved sin/cosine quad angles to configure
!-------------------------------------------------------

use dotlrt_variables
IMPLICIT NONE

integer i
integer j
integer k
real(8) g
real(8) t1
real(8) t2
real(8) t3
real(8) el
real(8) kl
real(8) kappa
real(8) kappa_sq
real(8) dkappa_sq_dg
real(8) xrf
real(8) yrf
real(8) zrf
real(8) rf
real(8) rd
external rf, rd

  do k = 1, nhg
    HGg(k)= (2.0d0*k-nhg-1.0d0)/dble(nhg)
  end do

DO k=1,nhg
 g=HGg(k)

 IF( DABS(g) >= 1) THEN
  WRITE(*,*) 'HG_phmat: k= ',k,' , HGg(k) = ', g, ' is not allowed'
  STOP
 END IF

 DO i=1,nstream
  DO j=i,nstream

   IF( g >= 0.d0 ) THEN

    t1=1+g*g+2*g*(cos_ang(i)*cos_ang(j)+sin_ang(i)*sin_ang(j))
    t2=1+g*g+2*g*(cos_ang(i)*cos_ang(j)-sin_ang(i)*sin_ang(j))
    kappa_sq=4*g*sin_ang(i)*sin_ang(j)/t1
    kappa = DSQRT(kappa_sq)
    xrf = 0.0d0
    yrf = 1.0d0 - kappa_sq
    zrf = 1.0d0
    el = rf(xrf,yrf,zrf) - (kappa_sq/3.0d0) * rd(xrf,yrf,zrf)

    ! "The Adding-doubling method, Scott Prahl" P.114
    HGph(k,i,j)=1.0d0/pi*(1-g*g)/dsqrt(t1)/t2*el 
    HGph(k,j,i)=HGph(k,i,j)

    dkappa_sq_dg=4*sin_ang(i)*sin_ang(j)*(t1-2*g*(g+cos_ang(i)*cos_ang(j)+sin_ang(i)*sin_ang(j)))/t1/t1

    IF( kappa > 0.01d0) THEN

     xrf = 0.0d0
     yrf = 1.0d0 - kappa_sq
     zrf = 1.0d0
     kl = rf(xrf,yrf,zrf)

     t3=(el-kl)/kappa_sq
    ELSE
     t3=-pi/512*( (((25*kappa_sq+30)*kappa_sq+48)*kappa_sq+128))
    END IF

   dHGph(k,i,j) = (-2*g/(1-g*g)                                        &
                   -  (g+cos_ang(i)*cos_ang(j)+sin_ang(i)*sin_ang(j))/t1                   &
                   -2*(g+cos_ang(i)*cos_ang(j)-sin_ang(i)*sin_ang(j))/t2                   &
                   +t3/el*dkappa_sq_dg/2             )*HGph(k,i,j)

   dHGph(k,j,i)=dHGph(k,i,j)

!----------------------   
    ELSE   ! g < 0
!----------------------   

    t1=1+g*g+2*g*(cos_ang(i)*cos_ang(j)-sin_ang(i)*sin_ang(j))
    t2=1+g*g+2*g*(cos_ang(i)*cos_ang(j)+sin_ang(i)*sin_ang(j))
    kappa_sq=4*DABS(g)*sin_ang(i)*sin_ang(j)/t1
    kappa   =DSQRT(kappa_sq)
    xrf = 0.0d0
    yrf = 1.0d0 - kappa_sq
    zrf = 1.0d0
    el = rf(xrf,yrf,zrf) - (kappa_sq/3.0d0) * rd(xrf,yrf,zrf)

    HGph(k,i,j)=1.0d0/pi*(1-g*g)/dsqrt(t1)/t2*el
    HGph(k,j,i)=HGph(k,i,j)

    dkappa_sq_dg=4*sin_ang(i)*sin_ang(j)*(t1-2*g*(g+cos_ang(i)*cos_ang(j)-sin_ang(i)*sin_ang(j)))/t1/t1

    IF( kappa > 0.01d0) THEN

     xrf = 0.0d0
     yrf = 1.0d0 - kappa_sq
     zrf = 1.0d0
     kl = rf(xrf,yrf,zrf)

     t3=(el-kl)/kappa_sq
    ELSE
     t3=-pi/512*( (((25*kappa_sq+30)*kappa_sq+48)*kappa_sq+128))
    END IF

   dHGph(k,i,j) =-( 2*g/(1-g*g)                                        &
                   +  (g+cos_ang(i)*cos_ang(j)-sin_ang(i)*sin_ang(j))/t1                   &
                   +2*(g+cos_ang(i)*cos_ang(j)+sin_ang(i)*sin_ang(j))/t2                   &
                   +t3/el*dkappa_sq_dg/2             )*HGph(k,i,j)

   dHGph(k,j,i)=dHGph(k,i,j)

    END IF ! g

  END DO ! j
 END DO ! i
END DO ! k

END SUBROUTINE HG_phmat
