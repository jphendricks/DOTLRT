!  (C) Copr. 1986-92 Numerical Recipes Software "1"7%A2.
!  This function computes Carlson's elliptic integral of the second kind, denoted Rd(x,y,z).
!   x and y must be nonnegative, and at most one can be zero; z must be positive.
!   TINY must be at least twice the negative 2/3 power of the machine underflow limit.
!   BIG must be at most 0.1xERRTOL times the negative 2/3 power of the machine overflow limit.
! 9/26/2020 Kevin Schaefer changed pause to stop in error statement
      real(8) function rd(x,y,z)
      implicit none
      REAL(8) x,y,z,ERRTOL,TINY,BIG,C1,C2,C3,C4,C5,C6
      PARAMETER (ERRTOL=.05d0,TINY=1.0d-25,BIG=4.5d21,C1=3.0d0/14.0d0,C2=1.0d0/6.0d0, &
               C3=9.0d0/22.0d0,C4=3.0d0/26.0d0,C5=0.25d0*C3,C6=1.5d0*C4)
      REAL(8) alamb,ave,delx,dely,delz,ea,eb,ec,ed,ee,fac,sqrtx,sqrty, &
             sqrtz,sum,xt,yt,zt
      if( min(x,y) .lt. 0.0d0 .or. min(x+y,z) .lt. TINY .or. &
          max(x,y,z) .gt. BIG) stop 'invalid arguments in rd'
      xt=x
      yt=y
      zt=z
      sum=0.0d0
      fac=1.0d0
1     continue
        sqrtx=dsqrt(xt)
        sqrty=dsqrt(yt)
        sqrtz=dsqrt(zt)
        alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
        sum=sum+fac/(sqrtz*(zt+alamb))
        fac=0.25d0*fac
        xt=0.25d0*(xt+alamb)
        yt=0.25d0*(yt+alamb)
        zt=0.25d0*(zt+alamb)
        ave=0.2d0*(xt+yt+3.0d0*zt)
        delx=(ave-xt)/ave
        dely=(ave-yt)/ave
        delz=(ave-zt)/ave
      if(max(abs(delx),dabs(dely),dabs(delz)).gt.ERRTOL)goto 1
      ea=delx*dely
      eb=delz*delz
      ec=ea-eb
      ed=ea-6.0d0*eb
      ee=ed+ec+ec
      rd=3.0d0*sum+fac*(1.0d0+ed*(-C1+C5*ed-C6*delz*ee)+delz*(C2*ee+delz*(-C3* &
          ec+delz*C4*ea)))/(ave*dsqrt(ave))
      return
      END function rd
