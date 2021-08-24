!
!--------------------------------------------------------------------------------
subroutine equ_coeff(nl,dl,k0,k,kx,ky,z0,z,a)
!--------------------------------------------------------------------------------
! calculates the elements of the coefficient matrix in ax=b based on formulation
! in the paper and saves them in matrix a.

implicit none
!
! define input variables
integer,         intent(in)  :: nl
real(kind=8),    intent(in)  :: dl(nl-1)
real(kind=8),    intent(in)  :: k0
real(kind=8),    intent(in)  :: kx
real(kind=8),    intent(in)  :: ky
real(kind=8),    intent(in)  :: z0
complex(kind=8), intent(in)  :: k(nl) ! (?) Wavenumber per layer
complex(kind=8), intent(in)  :: z(nl) ! Intrinsic impedance per layer
!
! define output variable
complex(kind=8), intent(out) :: a(4*nl,4*nl)
!
! define local variables
real(kind=8)    :: kro
complex(kind=8) :: i
complex(kind=8) :: kz0
complex(kind=8) :: kz(nl)
integer :: ii
integer :: jj
real(kind=8) exp_01
real(kind=8) exp_02
real(kind=8) exp_03
real(kind=8) exp_04


        i=(0.d0,1.d0)

        kro=sqrt(kx**2+ky**2)
        kz0=sqrt(k0**2-kro**2)
        do ii=1,nl
            kz(ii)=sqrt(k(ii)**2-kro**2)
        end do

        a(1,1)=kx/kro 
        a(1,2)=kz0*ky/k0/kro 
        a(1,3)=-kx/kro 
        a(1,4)=-kz(1)*ky/k(1)/kro 
        a(1,5)=-kx/kro
        a(1,6)=ky*kz(1)/k(1)/kro
        
        a(2,1)=ky/kro 
        a(2,2)=-kz0*kx/k0/kro 
        a(2,3)=-ky/kro 
        a(2,4)=kz(1)*kx/k(1)/kro 
        a(2,5)=-ky/kro 
        a(2,6)=-kx*kz(1)/k(1)/kro
        
        a(3,1)=-kz0*ky/k0/kro/z0
        a(3,2)=kx/kro/z0
        a(3,3)=kz(1)*ky/k(1)/kro/z(1)
        a(3,4)=-kx/kro/z(1)
        a(3,5)=-ky*kz(1)/k(1)/kro/z(1)
        a(3,6)=-kx/kro/z(1)
        
        a(4,1)=kz0*kx/k0/kro/z0
        a(4,2)=ky/kro/z0
        a(4,3)=-kz(1)*kx/k(1)/kro/z(1)
        a(4,4)=-ky/kro/z(1)
        a(4,5)=kx*kz(1)/k(1)/kro/z(1)
        a(4,6)=-ky/kro/z(1)
        
        do ii=1,4
            do jj=7,4*nl
                a(ii,jj)=(0.d0,0.d0)
            enddo
        end do

        do ii=2,nl-1

            do jj=1,4*ii-6
                a(4*ii-3,jj)=(0.d0,0.d0)
                a(4*ii-2,jj)=(0.d0,0.d0)
                a(4*ii-1,jj)=(0.d0,0.d0)
                a(4*ii,jj)=(0.d0,0.d0)
            enddo

            exp_01=exp(-i*kz(ii-1)*dl(ii-1))
	    exp_02=exp(i*kz(ii-1)*dl(ii-1))
	    exp_03=exp(-i*kz(ii)*dl(ii-1))
	    exp_04=exp(i*kz(ii)*dl(ii-1))
            a(4*ii-3,4*ii-5)=kx/kro                          *exp_01
            a(4*ii-3,4*ii-4)=kz(ii-1)*ky/k(ii-1)/kro         *exp_01
            a(4*ii-3,4*ii-3)=kx/kro                          *exp_02
            a(4*ii-3,4*ii-2)=-ky*kz(ii-1)/k(ii-1)/kro        *exp_02
            a(4*ii-3,4*ii-1)=-kx/kro                         *exp_03
            a(4*ii-3,4*ii)=-kz(ii)*ky/k(ii)/kro              *exp_03
            a(4*ii-3,4*ii+1)=-kx/kro                         *exp_04
            a(4*ii-3,4*ii+2)=ky*kz(ii)/k(ii)/kro             *exp_04

            a(4*ii-2,4*ii-5)=ky/kro                          *exp_01
            a(4*ii-2,4*ii-4)=-kz(ii-1)*kx/k(ii-1)/kro        *exp_01
            a(4*ii-2,4*ii-3)=ky/kro                          *exp_02
            a(4*ii-2,4*ii-2)=kx*kz(ii-1)/k(ii-1)/kro         *exp_02
            a(4*ii-2,4*ii-1)=-ky/kro                         *exp_03
            a(4*ii-2,4*ii)=kz(ii)*kx/k(ii)/kro               *exp_03
            a(4*ii-2,4*ii+1)=-ky/kro                         *exp_04
            a(4*ii-2,4*ii+2)=-kx*kz(ii)/k(ii)/kro            *exp_04
                
            a(4*ii-1,4*ii-5)=-ky*kz(ii-1)/k(ii-1)/kro/z(ii-1)*exp_01
            a(4*ii-1,4*ii-4)=kx/kro/z(ii-1)                  *exp_01
            a(4*ii-1,4*ii-3)=ky*kz(ii-1)/k(ii-1)/kro/z(ii-1) *exp_02
            a(4*ii-1,4*ii-2)=kx/kro/z(ii-1)                  *exp_02
            a(4*ii-1,4*ii-1)=ky*kz(ii)/k(ii)/kro/z(ii)       *exp_03
            a(4*ii-1,4*ii)=-kx/kro/z(ii)                     *exp_03
            a(4*ii-1,4*ii+1)=-ky*kz(ii)/k(ii)/kro/z(ii)      *exp_04
            a(4*ii-1,4*ii+2)=-kx/kro/z(ii)                   *exp_04      

            a(4*ii,4*ii-5)=kx*kz(ii-1)/k(ii-1)/kro/z(ii-1)   *exp_01
            a(4*ii,4*ii-4)=ky/kro/z(ii-1)                    *exp_01
            a(4*ii,4*ii-3)=-kx*kz(ii-1)/k(ii-1)/kro/z(ii-1)  *exp_02
            a(4*ii,4*ii-2)=ky/kro/z(ii-1)                    *exp_02
            a(4*ii,4*ii-1)=-kx*kz(ii)/k(ii)/kro/z(ii)        *exp_03
            a(4*ii,4*ii)=-ky/kro/z(ii)                       *exp_03
            a(4*ii,4*ii+1)=kx*kz(ii)/k(ii)/kro/z(ii)         *exp_04
            a(4*ii,4*ii+2)=-ky/kro/z(ii)                     *exp_04

            do jj=4*ii+3,4*nl
                a(4*ii-3,jj)=(0.d0,0.d0)
                a(4*ii-2,jj)=(0.d0,0.d0)
                a(4*ii-1,jj)=(0.d0,0.d0)
                a(4*ii,jj)=(0.d0,0.d0)
            enddo
        enddo

        do ii=4*nl-3,4*nl
            do jj=1,4*nl-6
                a(ii,jj)=(0.d0,0.d0)
            enddo
        enddo


        exp_01=exp(-i*kz(nl-1)*dl(nl-1))
	exp_02=exp(i*kz(nl-1)*dl(nl-1))
	exp_03=exp(i*kz(nl)*dl(nl-1))

        a(4*nl-3,4*nl-5)=kx/kro                          *exp_01
        a(4*nl-3,4*nl-4)=kz(nl-1)*ky/k(nl-1)/kro         *exp_01
        a(4*nl-3,4*nl-3)=kx/kro                          *exp_02
        a(4*nl-3,4*nl-2)=-ky*kz(nl-1)/k(nl-1)/kro        *exp_02
        a(4*nl-3,4*nl-1)=-kx/kro                         *exp_03
        a(4*nl-3,4*nl)=ky*kz(nl)/k(nl)/kro               *exp_03

        a(4*nl-2,4*nl-5)=ky/kro                          *exp_01
        a(4*nl-2,4*nl-4)=-kz(nl-1)*kx/k(nl-1)/kro        *exp_01
        a(4*nl-2,4*nl-3)=ky/kro                          *exp_02
        a(4*nl-2,4*nl-2)=kx*kz(nl-1)/k(nl-1)/kro         *exp_02
        a(4*nl-2,4*nl-1)=-ky/kro                         *exp_03
        a(4*nl-2,4*nl)=-kx*kz(nl)/k(nl)/kro              *exp_03
        
        a(4*nl-1,4*nl-5)=-ky*kz(nl-1)/k(nl-1)/kro/z(nl-1)*exp_01
        a(4*nl-1,4*nl-4)=kx/kro/z(nl-1)                  *exp_01
        a(4*nl-1,4*nl-3)=ky*kz(nl-1)/k(nl-1)/kro/z(nl-1) *exp_02
        a(4*nl-1,4*nl-2)=kx/kro/z(nl-1)                  *exp_02
        a(4*nl-1,4*nl-1)=-ky*kz(nl)/k(nl)/kro/z(nl)      *exp_03
        a(4*nl-1,4*nl)=-kx/kro/z(nl)                     *exp_03

        a(4*nl,4*nl-5)=kx*kz(nl-1)/k(nl-1)/kro/z(nl-1)   *exp_01
        a(4*nl,4*nl-4)=ky/kro/z(nl-1)                    *exp_01
        a(4*nl,4*nl-3)=-kx*kz(nl-1)/k(nl-1)/kro/z(nl-1)  *exp_02
        a(4*nl,4*nl-2)=ky/kro/z(nl-1)                    *exp_02
        a(4*nl,4*nl-1)=kx*kz(nl)/k(nl)/kro/z(nl)         *exp_03
        a(4*nl,4*nl)=-ky/kro/z(nl)                       *exp_03

return
end subroutine equ_coeff
