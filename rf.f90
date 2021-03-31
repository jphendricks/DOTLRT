!====================================================================================
real(8) FUNCTION rf(x,y,z)
!====================================================================================
!  This function computes Carlson's elliptic integral of the first kind, denoted Rf(x,y,z).
!   x, y, and z must be nonnegative, and at most one can be zero.
!   TINY must be at least 5 times the machine underflow limit.
!   BIG must be at most one-fifth of the machine overflow limit.
!
! History
!  1986-92 copyuright Numerical Recipes Software "1"7%A2.
!  9/26/2020 Kevin Schaefer changed pause to stop in error statement
!  12/12/2020 Kevin Schaefer brought code to 21st century by removing goto statements
!------------------------------------------------------------------------------------
implicit none
  real(8) x,y,z,ERRTOL,TINY,BIG,THIRD,C1,C2,C3,C4
  PARAMETER (ERRTOL=.08d0, TINY=1.5d-38, BIG=3.d37,THIRD=1.0d0/3.0d0, &
            C1=1.0d0/24.0d0,C2=.1d0,C3=3.0d0/44.0d0,C4=1.0d0/14.0d0)
  REAL(8) alamb,ave,delx,dely,delz,e2,e3,sqrtx,sqrty,sqrtz,xt,yt,zt
  integer iter  ! (-) iteration index
  integer,parameter :: max_iter=10000 ! (-) max allowed number of iterations

  if( min(x,y,z)<0.0d0 .or. min(x+y,x+z,y+z)<TINY .or. max(x,y,z)>BIG) stop 'invalid arguments in rf'

  xt=x
  yt=y
  zt=z
  do iter = 1, max_iter
    sqrtx=dsqrt(xt)
    sqrty=dsqrt(yt)
    sqrtz=dsqrt(zt)
    alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
    xt = 0.25d0 * (xt+alamb)
    yt = 0.25d0 * (yt+alamb)
    zt = 0.25d0 * (zt+alamb)
    ave=THIRD*(xt+yt+zt)
    delx=(ave-xt)/ave
    dely=(ave-yt)/ave
    delz=(ave-zt)/ave
    if(max(dabs(delx),dabs(dely),dabs(delz))<ERRTOL) exit
  enddo
  e2=delx*dely-delz**2
  e3=delx*dely*delz
  rf=(1.0d0+(C1*e2-C2-C3*e3)*e2+C4*e3)/dsqrt(ave)

return
END function rf
