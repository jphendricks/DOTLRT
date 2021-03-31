!  (C) Copr. 1986-92 Numerical Recipes Software "1"7%A2.
!  This subroutine computes all eigenvalues and eigenvectors of a real symmetric NxN matrix called 'a'.
!  On output, elements of 'a' above the diagonal are destroyed.
!  'd' is a vector of length N that returns the eigenvalues of 'a'.
!  'v' is an NxN matrix whose columns contian, on output, the normalised eigenvectors of 'a'.
!  'nrot' returns the number of Jacobi rotations that were required in the computation.
subroutine jacobi(a,n,np,d,v,nrot)
    implicit none
      integer n,np,nrot,NMAX
      real(8) a(np,np),d(np),v(np,np)
      PARAMETER (NMAX=500)
      integer i,ip,iq,j
      real(8) c,g,h,s,sm,t,tau,theta,tresh,b(NMAX),z(NMAX)
      do 12 ip=1,n
        do 11 iq=1,n
          v(ip,iq)=0.0d0
11      continue
        v(ip,ip)=1.0d0
12    continue
      do 13 ip=1,n
        b(ip)=a(ip,ip)
        d(ip)=b(ip)
        z(ip)=0.0d0
13    continue
      nrot=0
      do 24 i=1,50
        sm=0.0d0
        do 15 ip=1,n-1
          do 14 iq=ip+1,n
            sm=sm+dabs(a(ip,iq))
14        continue
15      continue
        if(sm.eq.0.d0)return
        if(i.lt.4)then
          tresh=0.2d0*sm/dble(n)**2
        else
          tresh=0.0d0
        endif
        do 22 ip=1,n-1
          do 21 iq=ip+1,n
            g=100.0d0*dabs(a(ip,iq))
            if((i.gt.4).and.(dabs(d(ip))+ &
     g.eq.dabs(d(ip))).and.(dabs(d(iq))+g.eq.dabs(d(iq))))then
              a(ip,iq)=0.0d0
            else if(dabs(a(ip,iq)).gt.tresh)then
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
              do 16 j=1,ip-1
                g=a(j,ip)
                h=a(j,iq)
                a(j,ip)=g-s*(h+g*tau)
                a(j,iq)=h+s*(g-h*tau)
16            continue
              do 17 j=ip+1,iq-1
                g=a(ip,j)
                h=a(j,iq)
                a(ip,j)=g-s*(h+g*tau)
                a(j,iq)=h+s*(g-h*tau)
17            continue
              do 18 j=iq+1,n
                g=a(ip,j)
                h=a(iq,j)
                a(ip,j)=g-s*(h+g*tau)
                a(iq,j)=h+s*(g-h*tau)
18            continue
              do 19 j=1,n
                g=v(j,ip)
                h=v(j,iq)
                v(j,ip)=g-s*(h+g*tau)
                v(j,iq)=h+s*(g-h*tau)
19            continue
              nrot=nrot+1
            endif
21        continue
22      continue
        do 23 ip=1,n
          b(ip)=b(ip)+z(ip)
          d(ip)=b(ip)
          z(ip)=0.0d0
23      continue
24    continue
      !pause 'too many iterations in jacobi'
      return
END subroutine jacobi
