
! SUBROUTINE HG_phmat( nteta, teta , ng, HGg, HGph , dHGph)
!-------------------------------------------------------
!   The Henyey-Greenstein (HG) scattering phase matrix is abreviated as HGph.
!   This subroutine calculates HG phase matrix for a set of angles:
!   teta(i), teta(j) (which are in degrees:  0 <= teta(i) <= 180,
!   i=1,..,nteta)  and for _ng values of asymmetry parameter _g:
!   g(1),...,g(ng)  ( |g| < 1 ).
!                 
!
!       HGph(k,i,j)=(1-g^2)/(2*pi)*INT(from x=0 to x=pi) dx
!
!       [1+g^2+2*g*cos(t1)*cos(t2)+2*g*sin(t1)*sin(t2)*cos(x)]^(-3/2)
!
!       where t1=teta(i)*pi/180, t2=teta(j)*pi/180, and g=HGg(k)
!
!   It also calculates a derivative of the phase matrix with respect 
!   to parameter _g:
!
!               dHGph(k,i,j)= d[HGph(k,i,j)]/dg
!
! -------------------   INPUT ------------------------
!
!   nteta        ( I)    number of angles
!
!   teta(nteta)  (DP)    array of angles (in deg) for which phase
!                        matrix and its derivative will be calculated
!
!    ng          ( I)    number of HG parameters to be calculated
!
!  HGg(ng)       (DP)    array of values of parameter _g
!
! -------------------  OUTPUT ------------------------
!
!  HGph(ng,nteta,nteta)  (DP)  array of HG phase matrix
!
! dHGph(ng,nteta,nteta)  (DP)  array of the appropriate derivatives
!                              of HG phase matrix
!
!-------------------------------------------------------
!-------------------------------------------------------

SUBROUTINE HG_phmat()
use variables

 IMPLICIT NONE

!INTEGER                                    , INTENT ( IN) :: nteta !, ng
! DOUBLE PRECISION, DIMENSION(nteta)         , INTENT ( IN) :: teta
! DOUBLE PRECISION, DIMENSION(ng   )         , INTENT ( IN) :: HGg
!DOUBLE PRECISION, DIMENSION(ng,nteta,nteta), INTENT (OUT) :: HGph,dHGph

INTEGER          :: i,j,k, nteta
real(8) ::   g, t1, t2, t3, el, eli, kl, kli, kappa, kappa_sq        &
                  , dkappa_sq_dg ! pi,

real(8),DIMENSION(nang) :: ct , st
real(8) xrf, yrf, zrf, rf, rd
external rf, rd

!-------------------------------------------------------

! pi=4*DATAN(1.d0)
nteta = nang
      do i = 1, nang
!        if (i == 1) then
!          teta(i) = 0.0
!        else if (i == nang) then
!          teta(i) = 180.0
!        else
!          teta(i) = quad_angle_array(i-1)
!        endif
          teta(i) = quad_angle_array(i)
      enddo


 DO i=1,nteta
  ct(i)=DCOS(teta(i)*pi/180.0d0)
  st(i)=DSIN(teta(i)*pi/180.0d0)
 END DO

  ng = 101 ! 10001

  do k = 1, ng
    HGg(k)= (2.0d0*k-ng-1.0d0)/dble(ng)
  end do


DO k=1,ng
 g=HGg(k)

 IF( DABS(g) >= 1) THEN
  WRITE(*,*) 'HG_phmat: k= ',k,' , HGg(k) = ', g, ' is not allowed'
  STOP
 END IF

 DO i=1,nteta
  DO j=i,nteta

   IF( g >= 0.d0 ) THEN

    t1=1+g*g+2*g*(ct(i)*ct(j)+st(i)*st(j))
    t2=1+g*g+2*g*(ct(i)*ct(j)-st(i)*st(j))
    kappa_sq=4*g*st(i)*st(j)/t1
    kappa = DSQRT(kappa_sq)
!    eli = DELE(kappa_sq)
    xrf = 0.0d0
    yrf = 1.0d0 - kappa_sq
    zrf = 1.0d0
    el = rf(xrf,yrf,zrf) - (kappa_sq/3.0d0) * rd(xrf,yrf,zrf)

    ! "The Adding-doubling method, Scott Prahl" P.114
    HGph(k,i,j)=1.0d0/pi*(1-g*g)/dsqrt(t1)/t2*el 
    HGph(k,j,i)=HGph(k,i,j)

    dkappa_sq_dg=4*st(i)*st(j)*(t1-2*g*(g+ct(i)*ct(j)+st(i)*st(j)))/t1/t1

    IF( kappa > 0.01d0) THEN

!     kli = DELK(kappa_sq)
     xrf = 0.0d0
     yrf = 1.0d0 - kappa_sq
     zrf = 1.0d0
     kl = rf(xrf,yrf,zrf)

     t3=(el-kl)/kappa_sq
    ELSE
     t3=-pi/512*( (((25*kappa_sq+30)*kappa_sq+48)*kappa_sq+128))
    END IF

   dHGph(k,i,j) = (-2*g/(1-g*g)                                        &
                   -  (g+ct(i)*ct(j)+st(i)*st(j))/t1                   &
                   -2*(g+ct(i)*ct(j)-st(i)*st(j))/t2                   &
                   +t3/el*dkappa_sq_dg/2             )*HGph(k,i,j)

   dHGph(k,j,i)=dHGph(k,i,j)

!----------------------   
    ELSE   ! g < 0
!----------------------   

    t1=1+g*g+2*g*(ct(i)*ct(j)-st(i)*st(j))
    t2=1+g*g+2*g*(ct(i)*ct(j)+st(i)*st(j))
    kappa_sq=4*DABS(g)*st(i)*st(j)/t1
    kappa   =DSQRT(kappa_sq)
!    eli = DELE(kappa_sq)
    xrf = 0.0d0
    yrf = 1.0d0 - kappa_sq
    zrf = 1.0d0
    el = rf(xrf,yrf,zrf) - (kappa_sq/3.0d0) * rd(xrf,yrf,zrf)

    HGph(k,i,j)=1.0d0/pi*(1-g*g)/dsqrt(t1)/t2*el
    HGph(k,j,i)=HGph(k,i,j)

    dkappa_sq_dg=4*st(i)*st(j)*(t1-2*g*(g+ct(i)*ct(j)-st(i)*st(j)))/t1/t1

    IF( kappa > 0.01d0) THEN

!     kli = DELK(kappa_sq)
     xrf = 0.0d0
     yrf = 1.0d0 - kappa_sq
     zrf = 1.0d0
     kl = rf(xrf,yrf,zrf)

     t3=(el-kl)/kappa_sq
    ELSE
     t3=-pi/512*( (((25*kappa_sq+30)*kappa_sq+48)*kappa_sq+128))
    END IF

   dHGph(k,i,j) =-( 2*g/(1-g*g)                                        &
                   +  (g+ct(i)*ct(j)-st(i)*st(j))/t1                   &
                   +2*(g+ct(i)*ct(j)+st(i)*st(j))/t2                   &
                   +t3/el*dkappa_sq_dg/2             )*HGph(k,i,j)

   dHGph(k,j,i)=dHGph(k,i,j)

    END IF ! g

  END DO ! j
 END DO ! i
END DO ! k

END SUBROUTINE HG_phmat
