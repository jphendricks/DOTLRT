!====================================================================================
subroutine jacobi(a,n,np,d,v,nrot)
!====================================================================================
! This subroutine computes all eigenvalues and eigenvectors of a real symmetric NxN matrix called 'a'.
! On output, elements of 'a' above the diagonal are destroyed.
!  'd' is a vector of length N that returns the eigenvalues of 'a'.
!  'v' is an NxN matrix whose columns contian, on output, the normalised eigenvectors of 'a'.
!  'nrot' returns the number of Jacobi rotations that were required in the computation.
!
! History
!  1986-92 copyright Numerical Recipes Software "1"7%A2.
!  12/12/2020 Kevin Schaefer brought code into 21st century by removing continue statements
!------------------------------------------------------------------------------------
implicit none

! Internal
! output
  integer n,np,nrot,NMAX
  real(8) a(np,np),d(np),v(np,np)
  PARAMETER (NMAX=500)
  integer i,ip,iq,j
  real(8) c,g,h,s,sm,t,tau,theta,tresh,b(NMAX),z(NMAX)

  do ip=1,n
    do iq=1,n
      v(ip,iq)=0.0d0
    enddo
    v(ip,ip)=1.0d0
  enddo

  do ip=1,n
    b(ip)=a(ip,ip)
    d(ip)=b(ip)
    z(ip)=0.0d0
  enddo

  nrot=0
  do i=1,50
    sm=0.0d0
    do ip=1,n-1
      do iq=ip+1,n
        sm=sm+dabs(a(ip,iq))
      enddo
    enddo

    if(sm.eq.0.d0)return
    if(i.lt.4)then
      tresh=0.2d0*sm/dble(n)**2
    else
      tresh=0.0d0
    endif

    do ip=1,n-1
      do iq=ip+1,n
        g=100.0d0*dabs(a(ip,iq))
        if((i.gt.4).and.(dabs(d(ip))+g==dabs(d(ip))).and.(dabs(d(iq))+g==dabs(d(iq))))then
          a(ip,iq)=0.0d0
        else if(dabs(a(ip,iq))>tresh)then
          h=d(iq)-d(ip)
          if(dabs(h)+g.eq.dabs(h))then
            t=a(ip,iq)/h
          else
            theta=0.5d0*h/a(ip,iq)
            t=1.0d0/(dabs(theta)+dsqrt(1.0d0+theta**2))
            if(theta.lt.0.0d0)t=-t
          endif
          c=1.0d0/dsqrt(1.0d0+t**2)
          s=t*c
          tau=s/(1.0d0+c)
          h=t*a(ip,iq)
          z(ip)=z(ip)-h
          z(iq)=z(iq)+h
          d(ip)=d(ip)-h
          d(iq)=d(iq)+h
          a(ip,iq)=0.0d0

          do j=1,ip-1
            g=a(j,ip)
            h=a(j,iq)
            a(j,ip)=g-s*(h+g*tau)
            a(j,iq)=h+s*(g-h*tau)
          enddo

          do j=ip+1,iq-1
            g=a(ip,j)
            h=a(j,iq)
            a(ip,j)=g-s*(h+g*tau)
            a(j,iq)=h+s*(g-h*tau)
          enddo

          do j=iq+1,n
            g=a(ip,j)
            h=a(iq,j)
            a(ip,j)=g-s*(h+g*tau)
            a(iq,j)=h+s*(g-h*tau)
          enddo

          do j=1,n
            g=v(j,ip)
            h=v(j,iq)
            v(j,ip)=g-s*(h+g*tau)
            v(j,iq)=h+s*(g-h*tau)
          enddo
          nrot=nrot+1
        endif
      enddo
    enddo

    do ip=1,n
      b(ip)=b(ip)+z(ip)
      d(ip)=b(ip)
      z(ip)=0.0d0
    enddo
  enddo

return
END subroutine jacobi
