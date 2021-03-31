!  (C) Copr. 1986-92 Numerical Recipes Software "1"7%A2.
!  This function computes Carlson's elliptic integral of the first kind, denoted Rf(x,y,z).
!   x, y, and z must be nonnegative, and at most one can be zero.
!   TINY must be at least 5 times the machine underflow limit.
!   BIG must be at most one-fifth of the machine overflow limit.
      real(8) FUNCTION rf(x,y,z)
      implicit none
      real(8) x,y,z,ERRTOL,TINY,BIG,THIRD,C1,C2,C3,C4
      PARAMETER (ERRTOL=.08d0, TINY=1.5d-38, BIG=3.d37,THIRD=1.0d0/3.0d0, &
                C1=1.0d0/24.0d0,C2=.1d0,C3=3.0d0/44.0d0,C4=1.0d0/14.0d0)
      REAL(8) alamb,ave,delx,dely,delz,e2,e3,sqrtx,sqrty,sqrtz,xt,yt,zt
      if( min(x,y,z) .lt. 0.0d0 .or. min(x+y,x+z,y+z) .lt. TINY .or. &
          max(x,y,z) .gt. BIG) pause 'invalid arguments in rf'
      xt=x
      yt=y
      zt=z
1     continue
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
      if(max(dabs(delx),dabs(dely),dabs(delz)).gt.ERRTOL)goto 1
      e2=delx*dely-delz**2
      e3=delx*dely*delz
      rf=(1.0d0+(C1*e2-C2-C3*e3)*e2+C4*e3)/dsqrt(ave)
      return
      END function rf
