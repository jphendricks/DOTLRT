
! Subroutines md_inv32.f90 and dlinrg.f90 are at the bottom of this file.

 SUBROUTINE core95 (nlr,nang,h,a0,b0,f,r,u,v,da0,db0,df,dr,du0,dv0     &
                   ,obs_lev,nvar,LP)
!  This subroutine solves radiative transfer equation,
!       specifics are below the list of modifications.
!   Written and modified by Alex Voronovich, Alexander.Voronovich@noaa.gov
!---------------------------------------------------------------------
!
!           MODIFICATION  05/15/03
!
!   Subroutine CORE95.F90 is renamed subroutine DTB32.F90 from ACCEL6
!
!---------------------------------------------------------------------
!
!           MODIFICATION  05/08/03
!
!   Symmetries of matrices _t, _r (i.e., p1 and p2), _rup, _rdn and
!   the relations:
!
!            S11^T=S22  , S12^T=S12 , S21^T=S21
!
!   are taken into account
!
!---------------------------------------------------------------------
!
!           MODIFICATION  04/28/03
!
!  In the first loop over _ilr acceleration due to symmetry of some matrices
!  is taken into account
!
!  Subroutine md_inv9 => md_inv31    - includes argument lama12    
!
!---------------------------------------------------------------------
!
!           MODIFICATION  3/4/03
!
!  It is taken into account if _rup, _rdn operators are diagonal  
!
!---------------------------------------------------------------------
!
!           MODIFICATION  2/14/03
!
!  Diagonal perturbations are treated differently  
!
!---------------------------------------------------------------------
!
!           MODIFICATION  1/31/03
!
!  Logical matrix input argument LP(0:nlr,nvar+3)  (layers' profile)
!  is added. Its _ilr row contains the following information:
!
!  To calculate _u,_v  ;  To calculate _du,_dv  ; Is it diagonal layer;
!  Is it diagonal layer with respect to _ivar-th parameter variation
!
!  Diagonal layers are treated differently
!
!
!----------------------------------------------------------------------
!
!  This subroutine solves radiative transfer equation taken in the form:
!
!  boundary condition at level n=N :    v_N (z=H) =  f_(N+1)
!
!  withiin layer n=N :
!
!  du_N/dz = -A0_N * u_N - B0_N * v_N + f_N
!  dv_N/dz =  B0_N * u_N + A0_N * v_N - f_N
!
!     ......................
!
!  within layer n=k :
!
!  du_k/dz = -A0_k * u_k - B0_k * v_k + f_k
!  dv_k/dz =  B0_k * u_k + A0_k * v_k - f_k
!
!     ......................
!
!  within layer n=1 :
!
!  du_1/dz = -A0_1 * u_1 - B0_1 * v_1 + f_1
!  dv_1/dz =  B0_1 * u_1 + A0_1 * v_1 - f_1
!
!  boundary condition at level n=0 :    u_0 (z=0) = R*v_0 (z=0) + f_0
!
!  where R(i), i=1,..,nang is reflection coefficient.
!
!  The solution is calculated by adding layer technique.
!
!
!  Output values _du0(ilr,n,nvar),_dv0(ilr,n,nvar) are derivatives of 
!  the up- and downwelling  brightness temperatures at the observation
!  level _obs_lev  corresponding to variations of the level properties
!  at level _ilr. Here _obs_lev could be equal to  any number within 
!  interval 
!
!                         0 <= obs_lev <= nrl
!
!  and 1 <= ilr <= nlr. For future use _du0, _dv0 are defined in fact
!  for ilr=0 also (when reflectivity properties will be varied).
!  However, in the present version of the code _du0(0,n,ivar) and 
! _dv0(0,n,ivar) do not have meaning.
!
!  Index _ivar corresponds to vector of different variations.
!  
!---------------------------------------------------------------------
!---------------------------------------------------------------------

!! USE MSIMSL

IMPLICIT NONE
integer*4 nrot
INTEGER                     , INTENT( IN) :: nlr,nang,obs_lev,nvar

DOUBLE PRECISION, DIMENSION(  nlr  )          , INTENT( IN) ::   h
DOUBLE PRECISION, DIMENSION(0:nlr+1,nang)     , INTENT( IN) ::   f
DOUBLE PRECISION, DIMENSION(0:nlr+1,nang,nvar), INTENT( IN) ::   df
DOUBLE PRECISION, DIMENSION(  nlr  ,nang,nang), INTENT( IN) ::  a0, b0 
DOUBLE PRECISION, DIMENSION(  nlr  ,nang,nang,nvar)                   &
                                              , INTENT( IN) :: da0,db0 

DOUBLE PRECISION, DIMENSION(nang,nang)        , INTENT( IN) ::   r,dr
LOGICAL         , DIMENSION(0:nlr,nvar+3)     , INTENT( IN) :: LP

DOUBLE PRECISION, DIMENSION(0:nlr,nang)       , INTENT(OUT) :: u,v
DOUBLE PRECISION, DIMENSION(0:nlr,nang,nvar)  , INTENT(OUT) :: du0,dv0

!DOUBLE PRECISION, DIMENSION(  nlr,nang,nang) ::  a, b ,am1, p1, p2, q  
!DOUBLE PRECISION, DIMENSION(  nlr,nang,nang,nvar)                     &
!                                             :: da,db,dp1,dp2

DOUBLE PRECISION, allocatable ::  a(:,:,:), b(:,:,:), am1(:,:,:), p1(:,:,:), p2(:,:,:), q(:,:,:)  
DOUBLE PRECISION, allocatable :: da(:,:,:,:), db(:,:,:,:), dp1(:,:,:,:), dp2(:,:,:,:)

DOUBLE PRECISION, DIMENSION(0:nlr,nang,nang) ::rup,rdn,s11,s12,s21,s22

DOUBLE PRECISION, DIMENSION(0:nlr,nang)      :: ust,vst
DOUBLE PRECISION, DIMENSION(0:nlr,nang,nvar) :: dz2,dz3

DOUBLE PRECISION, DIMENSION(  nlr,nang)      :: unh
DOUBLE PRECISION, DIMENSION(  nlr,nang,nvar) :: dunh

DOUBLE PRECISION, DIMENSION( nang  )         :: y0,y1,y2,dy0,b1,b2    

DOUBLE PRECISION, DIMENSION(nang,nang) ::  t0, t1, t2, w0, w1, ma, mab &
                                         , a11, a12, a21, a22, a11m1   &
                                         , a22m1,w2,t3,w3,lama12nm

DOUBLE PRECISION, DIMENSION(nang,nang,nvar) ::  dt0,dt1,dt2,dma,dmab   &
                                               ,da11,da12,da21,da22    &
                                               ,da11m1,da22m1,dw0,dw1

DOUBLE PRECISION, DIMENSION(nang)      :: lama,lamab,lama12
DOUBLE PRECISION, DIMENSION(nang,nvar) :: dlama,dlamab,db1,db2

DOUBLE PRECISION, DIMENSION(nang,nvar) :: dy0v,dy1v,dy2v

INTEGER                   :: ilr,n,m,k, oblv, ivar
DOUBLE PRECISION          :: s1,s2
LOGICAL, DIMENSION(  nlr) :: q_d
LOGICAL, DIMENSION(0:nlr) :: rup_d,rdn_d,s_d
LOGICAL                   :: a_d
integer alloc_err

INTEGER :: i,j
character*120 debugout

EXTERNAL md_inv32


allocate(a(nlr,nang,nang),b(nlr,nang,nang),am1(nlr,nang,nang), &
         p1(nlr,nang,nang),p2(nlr,nang,nang),q(nlr,nang,nang))
allocate(da(nlr,nang,nang,nvar),db(nlr,nang,nang,nvar), &
         dp1(nlr,nang,nang,nvar),dp2(nlr,nang,nang,nvar))

!-------------------------------------------------------------
 
!  OPEN(115, file='chk.tst',form='formatted')           

oblv=obs_lev

u(  0,:)=0.d0
v(  0,:)=0.d0
u(nlr,:)=0.d0
v(nlr,:)=0.d0

!kz, Symmetrization of matrices 
DO ilr=1,nlr
    IF( .NOT. LP(ilr,3) ) THEN
        a(ilr,:,:) = a0(ilr,:,:)+b0(ilr,:,:)
        b(ilr,:,:) = a0(ilr,:,:)-b0(ilr,:,:)

        ! symmetrization of matrices _a, _b below - principally not needed

        a(ilr,:,:)=(a(ilr,:,:)+TRANSPOSE(a(ilr,:,:)))/2   
        b(ilr,:,:)=(b(ilr,:,:)+TRANSPOSE(b(ilr,:,:)))/2   
    END IF

    IF(       LP(ilr,3) ) THEN
        ! blw
        a(ilr,:,:) = 0.0d0
        b(ilr,:,:) = 0.0d0
        ! blw
        DO n=1,nang
            a(ilr,n,n) = a0(ilr,n,n)+b0(ilr,n,n) !kz, Eqn. 30
            b(ilr,n,n) = a0(ilr,n,n)-b0(ilr,n,n)
        END DO
    END IF

    DO ivar=1,nvar
        IF( .NOT. LP(ilr,3+ivar) ) THEN
            da(ilr,:,:,ivar) = da0(ilr,:,:,ivar)+db0(ilr,:,:,ivar)
            db(ilr,:,:,ivar) = da0(ilr,:,:,ivar)-db0(ilr,:,:,ivar)

            ! symmetrization of matrices _da, _db below - principally not needed

            da(ilr,:,:,ivar)=(da(ilr,:,:,ivar)+TRANSPOSE(da(ilr,:,:,ivar)))/2   
            db(ilr,:,:,ivar)=(db(ilr,:,:,ivar)+TRANSPOSE(db(ilr,:,:,ivar)))/2 
        END IF

        IF(       LP(ilr,3+ivar) ) THEN
            ! blw
            da(ilr,:,:,ivar) = 0.0d0
            db(ilr,:,:,ivar) = 0.0d0
            ! blw
            DO n=1,nang
                da(ilr,n,n,ivar) = da0(ilr,n,n,ivar)+db0(ilr,n,n,ivar)
                db(ilr,n,n,ivar) = da0(ilr,n,n,ivar)-db0(ilr,n,n,ivar)
            END DO
        END IF
    END DO
END DO

rup(0,:,:) = r(  :,:)
ust(0,:  ) = f(0,:)

rup_d(0)=LP(0,3)     ! rup_d(0) is _TRUE if reflection matrix is diagonal

!------------------------------------------
!------------------------------------------
DO ilr=1,nlr   ! the first loop over _ilr
!------------------------------------------
!------------------------------------------

    rup_d(ilr)=rup_d(ilr-1) .AND. LP(ilr,3)
    q_d(ilr)=.FALSE.

    !----------------------------------------------------------
    IF( .NOT. LP(ilr,3)) THEN   !  non-diagonal case: part 0
    !----------------------------------------------------------
    !!  CALL DEVCSF(nang,a(ilr,:,:),nang,lama,ma,nang)
        call jacobi(a(ilr,:,:),nang,nang,lama,ma,nrot)
    !----------------------------------------------------------
    END IF                ! end of non-diagonal case: part 0
    !----------------------------------------------------------

    !----------------------------------------------------------
    IF(LP(ilr,3)) THEN   ! diagonal case: part 0
    !----------------------------------------------------------
       DO n=1,nang
            lama(n)=a(ilr,n,n)
       END DO
    !---------------------------------------------------------
    END IF                     ! end of diagonal case: part 0
    !---------------------------------------------------------

    lama12=DSQRT(lama) !kz: square root of eigenvalues of A

!    write(debugout,*) "a="
!    call mexPrintf(debugout//achar(10))

!    do i=1,nang
!       do j=1,nang
!          write(debugout,*) a(ilr,i,j)
!          call mexPrintf(debugout)
!       end do
!       call mexPrintf(achar(10))
!    end do
    

    DO ivar=0,nvar                    ! sic - ivar=0,..
        IF( .NOT. LP(ilr,ivar+3)) THEN
        DO n=1,nang
            lama12nm(n,n)=lama(n)
            DO m=n+1,nang
                lama12nm(n,m)=lama12(n)*lama12(m)
                lama12nm(m,n)=lama12nm(n,m)
            END DO
        END DO
        EXIT
        END IF
    END DO

    !---------------------------------------------------------
    IF(.NOT. LP(ilr,3) ) THEN  !    non-diagonal case: part 1
    !---------------------------------------------------------
        DO ivar=1,nvar

            IF(.NOT. LP(ilr,ivar+3) ) t3=MATMUL(da(ilr,:,:,ivar),ma)
 
            DO n=1,nang
                DO m=n,nang

                    !------------------------------------------------------------------
                    IF( .NOT. LP(ilr,ivar+3) ) THEN   ! non-diagonal perturbation
                    !------------------------------------------------------------------
                        s1=DOT_PRODUCT(ma(:,n),t3(:,m))
                    !------------------------------------------------------------------
                    ELSE                              !     diagonal perturbation
                    !------------------------------------------------------------------
                        s1=0.d0
                        DO k=1,nang
                            s1=s1+ma(k,n)*da(ilr,k,k,ivar)*ma(k,m)        
                        END DO
                    !------------------------------------------------------------------
                    END IF     ! end of diagonal/non-diagonal perturbation
                    !------------------------------------------------------------------

                    IF( n == m ) THEN
                        dlama(n,ivar)=s1
                        t2(n,n)=0.d0
                    ELSE
                        t2(n,m)=s1/(lama(m)-lama(n))
                        t2(m,n)=-t2(n,m)
                    END IF
                 END DO  ! m
             END DO  ! n

             dma(:,:,ivar)=MATMUL(ma,t2)  !kz, Eqn. 94

        END DO  ! ivar

        t3=MATMUL(TRANSPOSE(ma),b(ilr,:,:))

		DO n=1,nang
		    DO m=n,nang
				s1=0.d0
				s2=0.d0
				DO k=1,nang
					s1=s1+ma(n,k)/lama(k)*ma(m,k)
					s2=s2+t3(n,k)        *ma(k,m)
				END DO
				am1(ilr,n,m)=s1
				am1(ilr,m,n)=s1
				t1(n,m)=s2
				t1(m,n)=s2
			END DO
		 END DO


		 !   t1=MATMUL(t3,ma)=MATMUL(MATMUL(TRANSPOSE(ma),b(ilr,:,:)),ma)

		 t0=t1*lama12nm
		 w2=t1/lama12nm/2

		 DO ivar=1,nvar

			DO n=1,nang
				DO m=n,nang
					w3(n,m)=dlama(n,ivar)*lama(m)+lama(n)*dlama(m,ivar)
					w3(m,n)=w3(n,m)
				END DO
			END DO

			!---------------------------------------------------------------
			IF( .NOT. LP(ilr,ivar+3) ) THEN ! non-diagonal perturbation
			!---------------------------------------------------------------
				t2= MATMUL(db(ilr,:,:,ivar),ma)
				DO n=1,nang
					DO m=n,nang
						s1=0.d0
						DO k=1,nang
							s1=s1+ma(k,n)*t2(k,m)
						END DO
						t1(n,m)=s1
						t1(m,n)=s1
					END DO
				END DO
			!     t1= MATMUL(MATMUL(TRANSPOSE(ma),db(ilr,:,:,ivar)),ma)
			!---------------------------------------------------------------
			ELSE                            !     diagonal perturbation
			!---------------------------------------------------------------
				DO n=1,nang
					DO m=n,nang
						s1=0.d0
						DO k=1,nang
							s1=s1+ma(k,n)*db(ilr,k,k,ivar)*ma(k,m)
						END DO
						t1(n,m)=s1
						t1(m,n)=s1
					END DO
				END DO
			!---------------------------------------------------------------
			END IF         ! end of diagonal/non-diagonal perturbation
			!---------------------------------------------------------------

			t2=MATMUL(t3,dma(:,:,ivar))

			dt0(:,:,ivar)=lama12nm*(t2+TRANSPOSE(t2)+t1)+w3*w2             

		 END DO  ! ivar
	 !---------------------------------------------------------
	 END IF                 ! end of non-diagonal case: part 1
	 !---------------------------------------------------------

	 !---------------------------------------------------------
	 IF( LP(ilr,3) ) THEN                 !   diag case: part 1
     !---------------------------------------------------------
		DO ivar=1,nvar
			!-----------------------------------------------------------------
			IF( .NOT. LP(ilr,ivar+3) ) THEN  ! non-diagonal perturbation	
			!-----------------------------------------------------------------
				DO n=1,nang   
					dlama(n,ivar)=da(ilr,n,n,ivar)
					dma(n,n,ivar)=0.d0
					DO m=n+1,nang
						dma(n,m,ivar)= da(ilr,n,m,ivar)/(lama(m)-lama(n))
						dma(m,n,ivar)=   -dma(n,m,ivar)
					END DO
				END DO
			!-----------------------------------------------------------------
			ELSE                             !     diagonal perturbation
			!-----------------------------------------------------------------
				DO n=1,nang
					dlama(n,ivar)=da(ilr,n,n,ivar)
				END DO
			!-----------------------------------------------------------------
			END IF           ! end of diagonal/non-diagonal perturbation
			!-----------------------------------------------------------------
		END DO   ! ivar

		! blw
		t0 = 0.0d0
		! blw
		DO n=1,nang
			t0(n,n)=b(ilr,n,n)*lama(n)
		END DO

		DO ivar=1,nvar
			!-----------------------------------------------------------------
			IF( .NOT. LP(ilr,ivar+3) ) THEN  ! non-diagonal perturbation
			!-----------------------------------------------------------------
				DO n=1,nang
					DO m=n,nang
						dt0(n,m,ivar)=( dma(m,n,ivar)*b(ilr,m,m)+b(ilr,n,n)*dma(n,m,ivar) &
							+db(ilr,n,m,ivar))*lama12nm(n,m)
						dt0(m,n,ivar)=dt0(n,m,ivar)
					END DO
					dt0(n,n,ivar)=dt0(n,n,ivar)+dlama(n,ivar)*b(ilr,n,n)
				END DO
			!-----------------------------------------------------------------
			ELSE                             !     diagonal perturbation
			!-----------------------------------------------------------------     
				! blw
				dt0(:,:,ivar) = 0.0d0
				! blw
				DO n=1,nang
					dt0(n,n,ivar)=db(ilr,n,n,ivar)*lama(n)+dlama(n,ivar)*b(ilr,n,n)
				END DO
			!-----------------------------------------------------------------
			END IF           ! end of diagonal/non-diagonal perturbation
			!-----------------------------------------------------------------
		END DO   ! ivar
		!--------------------------------------------------------

	END IF                     ! end of diagonal case: part 1
	!--------------------------------------------------------

	!--------------------------------------------------------
	IF(LP(ilr,3)) THEN  !                diagonal case
	!--------------------------------------------------------
		DO n=1,nang
			lamab(n)=t0(n,n)        !kz, Eqn. 46??
		END DO
	!--------------------------------------------------------
	ELSE                !            non-diagonal case
	!--------------------------------------------------------
	!!  CALL DEVCSF(nang,t0,nang,lamab,mab,nang)
		call jacobi(t0,nang,nang,lamab,mab,nrot)
	!--------------------------------------------------------
	END IF              ! diagonal / non-diagonal case
	!--------------------------------------------------------

	!--------------------------------------------------------
	IF( .NOT. LP(ilr,3) ) THEN  ! non-giagonal case: part 2
	!--------------------------------------------------------
   DO ivar=1,nvar
    t1=MATMUL(dt0(:,:,ivar),mab)
    DO n=1,nang
     DO m=n,nang
      s1=DOT_PRODUCT(mab(:,n),t1(:,m))
      IF( n == m ) THEN
       dlamab(n,ivar)=s1
       t2(n,n)=0.d0
      ELSE
       t2(n,m)=s1/(lamab(m)-lamab(n))
       t2(m,n)=-t2(n,m)
      END IF
     END DO
    END DO
    dmab(:,:,ivar)=MATMUL(mab,t2)
   END DO  ! ivar
!------------------------------------------------------
 END IF              ! end of non-diagonal case: part 2
!------------------------------------------------------

!------------------------------------------------------
 IF( LP(ilr,3) ) THEN          !  diagonal case: part 2
!------------------------------------------------------
  DO ivar=1,nvar
   DO n=1,nang
    dlamab(n,ivar)=dt0(n,n,ivar)
   END DO
!-----------------------------------------------------------------
    IF( .NOT. LP(ilr,ivar+3) ) THEN  ! non-diagonal perturbation
!-----------------------------------------------------------------
     DO n=1,nang
      dmab(n,n,ivar)=0.d0
      DO m=n+1,nang
       dmab(n,m,ivar)= dt0(n,m,ivar)/(lamab(m)-lamab(n))
       dmab(m,n,ivar)=-dmab(n,m,ivar)
      END DO
     END DO
!-----------------------------------------------------------------
    END IF                             ! end of diagonal perturbation
!-----------------------------------------------------------------
  END DO
!------------------------------------------------------
 END IF                      ! end of diag case: part 2
!------------------------------------------------------

!----------------   calculation of _p1, _p2   (T, R operators) ----------

CALL md_inv32 (nang, h(ilr), lama, ma, lamab, mab, lama12, t1, t2             &
                         ,dlama,dma,dlamab,dmab,dt1,dt2,nvar,LP(ilr,:))

!write(debugout,*) "t"
!call mexPrintf(debugout//achar(10))

!do i=1,nang
!   do j=1,nang
!      write(debugout,*) t1(i,j)
!      call mexPrintf(debugout)
!   end do
!   call mexPrintf(achar(10))
!end do

!write(debugout,*) "r"
!call mexPrintf(debugout//achar(10))

!do i=1,nang
!   do j=1,nang
!      write(debugout,*) t2(i,j)
!      call mexPrintf(debugout)
!   end do
!   call mexPrintf(achar(10))
!end do



!------------------------------------------------------------------------

!------------------------------------------------------
 IF(.NOT. LP(ilr,3) ) THEN !  non-diagonal case: part 3
!------------------------------------------------------
   p1(ilr,:,:)=t1
   p2(ilr,:,:)=t2
!------------------------------------------------------
 ELSE                      !    diagonal case: part 3
!------------------------------------------------------
! blw
  p1(ilr,:,:) = 0.0d0
  p2(ilr,:,:) = 0.0d0
! blw
   DO n=1,nang
     p1(ilr,n,n  )= t1(n,n)
     p2(ilr,n,n  )= t2(n,n)
   END DO
!------------------------------------------------------
 END IF                    ! end of diagonal case: part 3
!------------------------------------------------------

 DO ivar=1,nvar
!---------------------------------------------------------------------
    IF( LP(ilr,3) .AND. LP(ilr,ivar+3) ) THEN  ! diagonal layer and
!                                              ! diagonal perturbation
!---------------------------------------------------------------------
! blw
  dp1(ilr,:,:,ivar) = 0.0d0
  dp2(ilr,:,:,ivar) = 0.0d0
! blw
   DO n=1,nang
    dp1(ilr,n,n,ivar)=dt1(n,n,ivar)
    dp2(ilr,n,n,ivar)=dt2(n,n,ivar)
   END DO
!---------------------------------------------------------------------
    ELSE                         !  non-diagonal layer or perturbation
!---------------------------------------------------------------------
   dp1(ilr,:,:,ivar)=dt1(:,:,ivar)
   dp2(ilr,:,:,ivar)=dt2(:,:,ivar)
!---------------------------------------------------------------------
    END IF           ! end of non-diagonal/diagonal layer/perturbation
!---------------------------------------------------------------------
 END DO    ! ivar

!                          calculation of matrix RUP

!-------------------------------------------------------
 IF( .NOT. LP(ilr,3) ) THEN ! non-diagonal case: part 4
!-------------------------------------------------------

!-----------------------------------------------------------
  IF( .NOT. rup_d(ilr-1) ) THEN   ! non-diagonal _rup(ilr-1)
!-----------------------------------------------------------
   t0=-MATMUL(p2(ilr,:,:),rup(ilr-1,:,:)) 
   t2= MATMUL(p1(ilr,:,:),rup(ilr-1,:,:))
!-----------------------------------------------------------
  ELSE                            !     diagonal _rup(ilr-1)
!-----------------------------------------------------------
   DO n=1,nang
    DO m=1,nang
     t0(n,m)=-p2(ilr,n,m)*rup(ilr-1,m,m)
     t2(n,m)= p1(ilr,n,m)*rup(ilr-1,m,m)
    END DO
   END DO
!-----------------------------------------------------------
  END IF                 ! non-diagonal/diagonal _rup(ilr-1)
!-----------------------------------------------------------

   DO n=1,nang
    t0(n,n)= t0(n,n)+1.d0
   END DO
   CALL dlinrg(nang,t0,nang,t1,nang)
     q(ilr,:,:)=t1                         ! q(ilr) = (1-p2*t)**(-1)

!   rup(ilr,:,:)=p2(ilr,:,:)+MATMUL(t2,MATMUL(q(ilr,:,:),p1(ilr,:,:)) )

   w0=MATMUL(q(ilr,:,:),p1(ilr,:,:))
   DO n=1,nang
    DO m=n,nang
     s1=p2(ilr,n,m)
     DO k=1,nang
      s1=s1+t2(n,k)*w0(k,m)
     END DO
     rup(ilr,n,m)=s1
     rup(ilr,m,n)=s1
    END DO
   END DO

!write(debugout,*) "Rup"
!call mexPrintf(debugout//achar(10))

!do i=1,nang
!   do j=1,nang
!      write(debugout,*) rup(1,i,j)
!      call mexPrintf(debugout)
!   end do
!   call mexPrintf(achar(10))
!end do

   
!-------------------------------------------------------
 END IF              ! end of non-diagonal case: part 4
!-------------------------------------------------------

!-------------------------------------------------------
 IF( LP(ilr,3) ) THEN           ! diagonal case: part 4
!-------------------------------------------------------
!------------------------------------------------------------
  IF( .NOT. rup_d(ilr-1) ) THEN    !  non-diagonal rup(ilr-1)
!------------------------------------------------------------
   DO n=1,nang
     s1=p2(ilr,n,n)
    DO m=1,nang
     t0(n,m)=-s1*rup(ilr-1,n,m)
    END DO
    t0(n,n)=t0(n,n)+1.d0
   END DO
   CALL dlinrg(nang,t0,nang,t1,nang)
   q(ilr,:,:)=t1
   t2=MATMUL(rup(ilr-1,:,:),t1)
   DO n=1,nang
    DO m=n,nang
     rup(ilr,n,m)=p1(ilr,n,n)*t2(n,m)*p1(ilr,m,m)
     rup(ilr,m,n)=rup(ilr,n,m)
    END DO
    rup(ilr,n,n)=rup(ilr,n,n)+p2(ilr,n,n)
   END DO
!-----------------------------------------------------------
  ELSE                      ! diagonal rup(ilr-1)
!-----------------------------------------------------------
!   rup(ilr,:,:)=0.d0            ! remove  after debugging
!     q(ilr,:,:)=0.d0            ! remove  after debugging
! blw 911
  q(ilr,:,:) = 0.0d0
  rup(ilr,:,:) = 0.0d0
! blw 911
   DO n=1,nang
    q  (ilr,n,n)=1/(1-p2(ilr,n,n)*rup(ilr-1,n,n))
    rup(ilr,n,n)= p1(ilr,n,n)*rup(ilr-1,n,n)*q(ilr,n,n)*p1(ilr,n,n)     &
                 +p2(ilr,n,n)
   END DO
   q_d(ilr)=.TRUE.
!-----------------------------------------------------------
 END IF              ! non-diagonal/diagonal rup(ilr-1)
!-----------------------------------------------------------

!-------------------------------------------------------
 END IF                   ! end of diagonal case: part 4
!-------------------------------------------------------

!--------------------------------------------------------
 IF( .NOT. LP(ilr,3) ) THEN  ! non-diagonal case: part 5
!--------------------------------------------------------
  y0=MATMUL(am1(ilr,:,:),f(ilr,:))
  t1=p1(ilr,:,:)+p2(ilr,:,:)
  unh(ilr,:)= y0-MATMUL(t1,y0)
  DO ivar=1,nvar
    dy0=MATMUL(am1(ilr,:,:),df(ilr,:,ivar)-MATMUL(da(ilr,:,:,ivar),y0))
   dunh(ilr,:,ivar)=dy0-MATMUL(dp1(ilr,:,:,ivar)+dp2(ilr,:,:,ivar),y0)  &
                       -MATMUL(t1,dy0) 
  END DO

!  write(debugout,*) "am1"
!  call mexPrintf(debugout//achar(10))
!  do i=1,nang
!     do j=1,nang
!        write(debugout,*) am1(ilr,i,j)
!        call mexPrintf(debugout)
!     end do
!     call mexPrintf(achar(10))
!  end do

!  write(debugout,*) "f"
!  call mexPrintf(debugout//achar(10))
!  do i=1,nang
!     write(debugout,*) f(ilr,i)
!     call mexPrintf(debugout)
!  end do
!  call mexPrintf(achar(10))
  
!  write(debugout,*) "y0"
!  call mexPrintf(debugout//achar(10))
!  do i=1,nang
!     write(debugout,*) y0(i)
!     call mexPrintf(debugout)
!  end do
!  call mexPrintf(achar(10))
  
!  write(debugout,*) "unh="
!  call mexPrintf(debugout//achar(10))
!  do i=1,nang
!     write(debugout,*) unh(ilr,i)
!     call mexPrintf(debugout)
!  end do
!  call mexPrintf(achar(10))
  
!---------------------------------------------------------
  IF( .NOT. rup_d(ilr-1)) THEN   ! non-diagonal rup(ilr-1)
!---------------------------------------------------------
   y0=ust(ilr-1,:)+MATMUL(rup(ilr-1,:,:),unh(ilr,:))
   ust(ilr,:)=unh(ilr,:)+MATMUL(p1(ilr,:,:),y0+MATMUL(rup(ilr-1,:,:)    &
        ,MATMUL(q (ilr,:,:),   MATMUL(p2(ilr,:,:),y0))))
   
!--------------------------------------------------------
  ELSE  ! diagonal rup(ilr-1)
!--------------------------------------------------------
   DO n=1,nang
    y0(n)=ust(ilr-1,n)+rup(ilr-1,n,n)*unh(ilr,n)
   END DO
   y1=MATMUL(q(ilr,:,:),MATMUL(p2(ilr,:,:),y0))
   DO n=1,nang
    y2(n)=y0(n)+rup(ilr-1,n,n)*y1(n)
   END DO
   ust(ilr,:)=unh(ilr,:)+MATMUL(p1(ilr,:,:),y2)
!--------------------------------------------------------
  END IF  ! non-diagonal/diagonal rup(ilr-1)
!--------------------------------------------------------     
!-------------------------------------------------------
 END IF               ! end of non-diagonal case: part 5
!-------------------------------------------------------

!-------------------------------------------------------
  IF( LP(ilr,3) ) THEN           ! diagonal case: part 5
!-------------------------------------------------------
   DO n=1,nang
    y0(n)=f(ilr,n)/a(ilr,n,n)
    unh(ilr,n)=(1-p1(ilr,n,n)-p2(ilr,n,n))*y0(n)
   END DO
    DO ivar=1,nvar
!-----------------------------------------------------------------
     IF( .NOT. LP(ilr,ivar+3) ) THEN  ! non-diagonal perturbation
!-----------------------------------------------------------------
      y1=df(ilr,:,ivar)-MATMUL(da(ilr,:,:,ivar),y0)
!-----------------------------------------------------------------
     ELSE                             !     diagonal perturbation
!-----------------------------------------------------------------
      DO n=1,nang
       y1(n)=df(ilr,n,ivar)-da(ilr,n,n,ivar)*y0(n)
      END DO
!-----------------------------------------------------------------
     END IF           ! end of diagonal/non-diagonal perturbation
!-----------------------------------------------------------------
     DO n=1,nang
      dy0(n)=y1(n)/a(ilr,n,n)
     END DO
!-----------------------------------------------------------------
    IF( .NOT. LP(ilr,ivar+3) ) THEN  ! non-diagonal perturbation
!-----------------------------------------------------------------
     dunh(ilr,:,ivar)=dy0-MATMUL(dp1(ilr,:,:,ivar)+dp2(ilr,:,:,ivar),y0)
!-----------------------------------------------------------------
    ELSE                             !     diagonal perturbation
!-----------------------------------------------------------------
     DO n=1,nang
      dunh(ilr,n,ivar)=dy0(n)-(dp1(ilr,n,n,ivar)+dp2(ilr,n,n,ivar))*y0(n)
     END DO
!-----------------------------------------------------------------
    END IF           ! end of diagonal/non-diagonal perturbation
!-----------------------------------------------------------------

     DO n=1,nang
      dunh(ilr,n,ivar)=dunh(ilr,n,ivar)-(p1(ilr,n,n)+p2(ilr,n,n))*dy0(n)
     END DO

!---------------
    END DO    !   ivar

!----------------------------------------------------------
  IF( .NOT. rup_d(ilr-1) ) THEN   ! non-diagonal rup(ilr-1)
!---------------------------------------------------------- 
   y0=ust(ilr-1,:)+MATMUL(rup(ilr-1,:,:),unh(ilr,:))
   DO n=1,nang
    y1(n)=p2(ilr,n,n)*y0(n)
   END DO
   y2=y0+MATMUL(rup(ilr-1,:,:),MATMUL(q(ilr,:,:),y1))
!---------------------------------------------------------
  ELSE                           !     diagonal rup(ilr-1)
!---------------------------------------------------------
   DO n=1,nang
    y0(n)=ust(ilr-1,n)+rup(ilr-1,n,n)*unh(ilr,n)
    y1(n)=p2(ilr,n,n)*y0(n)    
   END DO
   b1=MATMUL(q(ilr,:,:),y1)
   DO n=1,nang
    y2(n)=y0(n)+rup(ilr-1,n,n)*b1(n)    
   END DO
!--------------------------------------------------------
  END IF               ! non-diagonal/diagonal rup(ilr-1)
!--------------------------------------------------------

   DO n=1,nang
    ust(ilr,n)=unh(ilr,n)+p1(ilr,n,n)*y2(n)
   END DO

!------------------------------------------------------
  END IF                 ! end of diagonal case: part 5
!------------------------------------------------------

!---------------------------------------------------------
!---------------------------------------------------------
END DO         !  end of the first _ilr loop
!---------------------------------------------------------
!---------------------------------------------------------


!-----------------------------------------------------------------
!  The second loop over _ilr calculates vectors _u and _v which
!  are up- and downwelling brightness temperatures. They are also
!  calculated in the 5-th loop over _ilr as a by-product of
!  calculations of derivatives a little bit faster. 
!  This mode of calculation should be used when variation of the
!  _ilr layer parameters is not performed, however brightness
!  temperature at the top OR at the bottom of the layer still has
!  to be calculated
!------------------------------------------------------------------

!------------------------------------------------------------------
!           Loop 2 is not used in this version ! 
!           It could be "uncommented" for future use
!------------------------------------------------------------------
!
!  v(nlr,:)=f(nlr+1,:)
!  u(nlr,:)=ust(nlr,:)+MATMUL(rup(nlr,:,:),v(nlr,:))
!
!!-------------------------------------------------
!!-------------------------------------------------
! DO ilr=nlr,1,-1   !         the second _ilr loop
!!-------------------------------------------------
!!-------------------------------------------------
!
!!------------------------------------------------------------------
!  IF( .NOT. LP(ilr,3) ) THEN   ! non-diagonal layer
!!------------------------------------------------------------------
!   v(ilr-1,:)=MATMUL(q(ilr,:,:),unh(ilr,:)+MATMUL(p1(ilr,:,:),v  (ilr  ,:))  &
!                                          +MATMUL(p2(ilr,:,:),ust(ilr-1,:)))
!!------------------------------------------------------------------
!  ELSE                         !     diagonal layer
!!------------------------------------------------------------------
!
!   DO n=1,nang
!    y1(n)=unh(ilr,n)+p1(ilr,n,n)*v(ilr,n)+p2(ilr,n,n)*ust(ilr-1,n)
!   END DO
!
!!-----------------------------------------------------
!  IF( .NOT. q_d(ilr) ) THEN    !  non-diagonal _q(ilr)
!!-----------------------------------------------------
!   v(ilr-1,:)=MATMUL(q(ilr,:,:),y1)
!!-----------------------------------------------------
!  ELSE                         !      diagonal _q(ilr)
!!-----------------------------------------------------
!   DO n=1,nang
!    v(ilr-1,n)=q(ilr,n,n)*y1(n)
!   END DO
!!-----------------------------------------------------
!  END IF               ! non-diagonal/diagonal _q(ilr)
!!-----------------------------------------------------
!
!!------------------------------------------------------------------
!  END IF                       ! end of non-diagonal/diagonal layer
!!------------------------------------------------------------------
!
!
!!-----------------------------------------------------------
!  IF( .NOT. rup_d(ilr-1) ) THEN   ! non-diagonal rup(ilr-1)
!!-----------------------------------------------------------
!   u(ilr-1,:)=ust(ilr-1,:)+MATMUL(rup(ilr-1,:,:),v(ilr-1,:)) 
!!----------------------------------------------------------- 
!  ELSE                           !      diagonal rup(ilr-1)
!!-----------------------------------------------------------
!   DO n=1,nang
!    u(ilr-1,n)=ust(ilr-1,n)+rup(ilr-1,n,n)*v(ilr-1,n)
!   END DO
!!-----------------------------------------------------------
!  END IF               ! non-diagonal/diagonal rup(ilr-1)
!!-----------------------------------------------------------
!
!!-------------------------------------------------
!!-------------------------------------------------
! END DO            !  end of the second _ilr loop
!!-------------------------------------------------
!!-------------------------------------------------

 vst(nlr  ,:)=f(nlr+1,:)
 vst(nlr-1,:)=unh(nlr,:)+MATMUL(p1(nlr,:,:),vst(nlr,:))

 rdn=0.d0                        ! could be deleted after debugging

 rdn(nlr  ,:,:)=0.d0
 rdn(nlr-1,:,:)=p2(nlr,:,:)

 rdn_d(nlr  )=.TRUE.
 rdn_d(nlr-1)=LP(nlr,3)

!------------------------------------------------------------------
 DO ilr=nlr-1,1,-1     !  the third loop: calculation of _rdn, _vst
!------------------------------------------------------------------

   rdn_d(ilr-1)=rdn_d(ilr).AND.LP(ilr,3)

!----------------------------------------------------------
 IF( .NOT. LP(ilr,3) ) THEN ! non-diagonal case: part 6
!----------------------------------------------------------

!------------------------------------------------------
 IF( .NOT. rdn_d(ilr) ) THEN  ! rdn(ilr) - non-diagonal 
!------------------------------------------------------
  t0=-MATMUL(rdn(ilr,:,:),p2(ilr,:,:))
  DO n=1,nang
!!    t0(n,n)=1+t0(n,n)
    t0(n,n) = 1.0d0 + t0(n,n)
  END DO
  CALL dlinrg(nang,t0,nang,t1,nang)
  t2=MATMUL(p1(ilr,:,:),t1)

!  rdn(ilr-1,:,:)=p2(ilr,:,:)+MATMUL(MATMUL(t2,rdn(ilr,:,:)),p1(ilr,:,:))

   w0=MATMUL(rdn(ilr,:,:),p1(ilr,:,:))
   DO n=1,nang
    DO m=n,nang
     s1=p2(ilr,n,m)
     DO k=1,nang
      s1=s1+t2(n,k)*w0(k,m)
     END DO
     rdn(ilr-1,n,m)=s1
     rdn(ilr-1,m,n)=s1
    END DO
   END DO

  vst(ilr-1,:)=unh(ilr,:)+MATMUL(t2,vst(ilr,:)+MATMUL(rdn(ilr,:,:)      &
                                                     ,unh(ilr,:  ) ))
!-----------------------------
 ELSE   ! rdn(ilr) - diagonal
!-----------------------------
  DO n=1,nang
   DO m=1,nang
    t0(n,m)=-rdn(ilr,n,n)*p2(ilr,n,m)
   END DO
!!   t0(n,n)=1+t0(n,n)
   t0(n,n) = 1.0d0 + t0(n,n)
  END DO
  CALL dlinrg(nang,t0,nang,t1,nang)
  t2=MATMUL(p1(ilr,:,:),t1)
  DO n=1,nang
   DO m=n,nang
    s1=p2(ilr,n,m)
    DO k=1,nang
     s1=s1+t2(n,k)*rdn(ilr,k,k)*p1(ilr,k,m)
    END DO
    rdn(ilr-1,n,m)=s1
    rdn(ilr-1,m,n)=s1
   END DO
   y0(n)=vst(ilr,n)+rdn(ilr,n,n)*unh(ilr,n)
  END DO
  
  vst(ilr-1,:)=unh(ilr,:)+MATMUL(t2,y0)
!-----------------------------------------
 END IF ! rdn(ilr) - non-diagonal/diagonal
!-----------------------------------------
!------------------------------------------------------
 END IF              ! end of non-diagonal case: part 6
!------------------------------------------------------

!------------------------------------------------------
 IF( LP(ilr,3) ) THEN           ! diagonal case: part 6
!------------------------------------------------------
!------------------------------------------------------
 IF( .NOT. rdn_d(ilr) ) THEN  ! rdn(ilr) - non-diagonal 
!------------------------------------------------------
  DO n=1,nang
   DO m=1,nang
    t0(n,m)=-rdn(ilr,n,m)*p2(ilr,m,m)
   END DO
!!   t0(n,n)=1+t0(n,n)
   t0(n,n) = 1.0d0 + t0(n,n)
  END DO
  CALL dlinrg(nang,t0,nang,t1,nang)
  DO n=1,nang
   DO m=1,nang
    t2(n,m)=p1(ilr,n,n)*t1(n,m)
   END DO
  END DO
  t1=MATMUL(t2,rdn(ilr,:,:))
  DO n=1,nang
   DO m=n,nang
    rdn(ilr-1,n,m)=t1(n,m)*p1(ilr,m,m)
    rdn(ilr-1,m,n)=rdn(ilr-1,n,m)
   END DO
    rdn(ilr-1,n,n)=rdn(ilr-1,n,n)+p2(ilr,n,n)
  END DO

  vst(ilr-1,:)=unh(ilr,:)+MATMUL(t2,vst(ilr,:)+MATMUL(rdn(ilr,:,:)      &
                                                     ,unh(ilr,:  ) ))
!-----------------------------
 ELSE   ! rdn(ilr) - diagonal
!-----------------------------
! blw 911
  rdn(ilr-1,:,:) = 0.0d0
! blw 911
  DO n=1,nang
   s1=p1(ilr,n,n)/(1-rdn(ilr,n,n)*p2(ilr,n,n))
   rdn(ilr-1,n,n)=p2(ilr,n,n)+p1(ilr,n,n)*s1*rdn(ilr,n,n)
   vst(ilr-1,n  )=unh(ilr,n)+s1*(vst(ilr,n)+rdn(ilr,n,n)*unh(ilr,n))
  END DO
!-----------------------------------------
 END IF ! rdn(ilr) - non-diagonal/diagonal
!-----------------------------------------
!------------------------------------------------------
 END IF                         ! diagonal case: part 6
!------------------------------------------------------

!--------------------------------------------------------------------
END DO   !  end of the third loop over _ilr: calculation of _rdn,_vst
!--------------------------------------------------------------------

 oblv=MAX0(0  ,oblv)
 oblv=MIN0(nlr,oblv)

!write(debugout,*) "oblv=",oblv
!call mexPrintf(debugout//achar(10))

!-----------   the forth _ilr loop: calculation of s11,s12,s21,s22

 s11=0.d0                      !  to remove after debugging
 s12=0.d0                      !  to remove after debugging
 s21=0.d0                      !  to remove after debugging
 s22=0.d0                      !  to remove after debugging

 DO n=1,nang
  s11(oblv,n,n)=1.d0
  s22(oblv,n,n)=1.d0
 END DO
 s_d(oblv)=.TRUE.

 IF(oblv < nlr) THEN
  IF(.NOT.LP(oblv+1,3)) THEN
   s11(oblv+1,:,:)=p1(oblv+1,:,:)
   s22(oblv+1,:,:)=p1(oblv+1,:,:)
   s12(oblv+1,:,:)=p2(oblv+1,:,:)
   s21(oblv+1,:,:)=p2(oblv+1,:,:)
  ELSE
   DO n=1,nang
    s11(oblv+1,n,n)=p1(oblv+1,n,n)
    s22(oblv+1,n,n)=p1(oblv+1,n,n)
    s12(oblv+1,n,n)=p2(oblv+1,n,n)
    s21(oblv+1,n,n)=p2(oblv+1,n,n)
   END DO
  END IF
  s_d(oblv+1)=LP(oblv+1,3)
 END IF

!--------------------------------------------------------------------- 
DO ilr=oblv+2,nlr  ! the first half of the forth _ilr loop: top layers
!---------------------------------------------------------------------

  s_d(ilr)=s_d(ilr-1).AND.LP(ilr,3)

!------------------------------------------------------------
  IF( .NOT. LP(ilr,3) ) THEN    !   non-diagonal case: part 7
!------------------------------------------------------------

!------------------------------------------------------------
   IF( .NOT. s_d(ilr-1) ) THEN    !   non-diagonal s_d(ilr-1)
!------------------------------------------------------------
 t0=-MATMUL(s12(ilr-1,:,:), p2(ilr  ,:,:))
 w0=-MATMUL( p2(ilr  ,:,:),s12(ilr-1,:,:))
 DO n=1,nang
!!  t0(n,n)=1+t0(n,n)
!!  w0(n,n)=1+w0(n,n)
  t0(n,n) = 1.0d0 + t0(n,n)
  w0(n,n) = 1.0d0 + w0(n,n)
 END DO
 CALL dlinrg(nang,t0,nang,t1,nang)   ! t1= (1-s12*r)**(-1)
 CALL dlinrg(nang,w0,nang,w1,nang)   ! w1= (1-r*s12)**(-1)
  t2=MATMUL(t1,s11(ilr-1,:,:))
  w2=MATMUL(w1, p1(ilr  ,:,:))

!  s11(ilr,:,:)=MATMUL( p1(ilr  ,:,:),t2)
!  s22(ilr,:,:)=MATMUL(s22(ilr-1,:,:),w2)
!  s12(ilr,:,:)=p2(ilr   ,:,:)+MATMUL( MATMUL(p1(ilr,:,:),s12(ilr-1,:,:)),w2)
!  s21(ilr,:,:)=s21(ilr-1,:,:)+MATMUL( MATMUL(s22(ilr-1,:,:),p2(ilr,:,:)),t2)


   s11(ilr,:,:)=MATMUL( p1(ilr  ,:,:),t2)
   s22(ilr,:,:)=TRANSPOSE(s11(ilr,:,:))

   w3=MATMUL(p1(ilr,:,:),s12(ilr-1,:,:))
   t3=MATMUL(s22(ilr-1,:,:),p2(ilr,:,:))

   DO n=1,nang
    DO m=n,nang
     s1=p2(ilr,n,m)
     s2=s21(ilr-1,n,m)
     DO k=1,nang
      s1=s1+w3(n,k)*w2(k,m)
      s2=s2+t3(n,k)*t2(k,m)
     END DO
     s12(ilr,n,m)=s1
     s21(ilr,n,m)=s2
     s12(ilr,m,n)=s1
     s21(ilr,m,n)=s2      
    END DO
   END DO

!----------------------------------
   ELSE     !  diagonal s_d(ilr-1)
!----------------------------------
    DO n=1,nang
     DO m=1,nang
      t0(n,m)=-s12(ilr-1,n,n)*p2(ilr,n,m)
      w0(n,m)=-s12(ilr-1,m,m)*p2(ilr,n,m)
     END DO
!!     t0(n,n)=1+t0(n,n)
!!     w0(n,n)=1+w0(n,n)
     t0(n,n) = 1.0d0 + t0(n,n)
     w0(n,n) = 1.0d0 + w0(n,n)
    END DO
    CALL dlinrg(nang,t0,nang,t1,nang)   ! t1= (1-s12*r)**(-1)
    CALL dlinrg(nang,w0,nang,w1,nang)   ! w1= (1-r*s12)**(-1)
    t2=MATMUL(p1(ilr,:,:),t1)
    w2=MATMUL(w1,p1(ilr,:,:))
    t3=MATMUL(p2(ilr,:,:),t1)

!    DO n=1,nang
!     DO m=1,nang
!      s11(ilr,n,m)=t2(n,m)*s11(ilr-1,m,m)
!      s22(ilr,n,m)=        s22(ilr-1,n,n)*w2(n,m)       
!      s21(ilr,n,m)=s22(ilr-1,n,n)*t3(n,m)*s11(ilr-1,m,m)
!       s1=p2(ilr,n,m)
!       DO k=1,nang
!        s1=s1+p1(ilr,n,k)*s12(ilr-1,k,k)*w2(k,m)
!       END DO
!      s12(ilr,n,m)=s1       
!     END DO
!     s21(ilr,n,n)=s21(ilr,n,n)+s21(ilr-1,n,n)
!    END DO

    DO n=1,nang
     DO m=1,nang
      s11(ilr,n,m)=t2(n,m)*s11(ilr-1,m,m)
     END DO
    END DO

    s22(ilr,:,:)=TRANSPOSE(s11(ilr,:,:))

    DO n=1,nang
     DO m=n,nang       
      s21(ilr,n,m)=s22(ilr-1,n,n)*t3(n,m)*s11(ilr-1,m,m)
      s21(ilr,m,n)=s21(ilr,n,m)
       s1=p2(ilr,n,m)
       DO k=1,nang
        s1=s1+p1(ilr,n,k)*s12(ilr-1,k,k)*w2(k,m)
       END DO
      s12(ilr,n,m)=s1
      s12(ilr,m,n)=s1       
     END DO
     s21(ilr,n,n)=s21(ilr,n,n)+s21(ilr-1,n,n)
    END DO

!------------------------------------------------------------
  END IF    ! end of the non-diagonal/diagonal s_d(ilr-1)
!------------------------------------------------------------

!------------------------------------------------------------
  END IF               ! end of the non-diagonal case: part 7
!------------------------------------------------------------
!------------------------------------------------------------
  IF( LP(ilr,3) ) THEN           !      diagonal case: part 7
!------------------------------------------------------------

!------------------------------------------------------------
   IF( .NOT. s_d(ilr-1) ) THEN    !   non-diagonal s_d(ilr-1)
!------------------------------------------------------------
   DO n=1,nang
    DO m=1,nang
     t0(n,m)=-s12(ilr-1,n,m)*p2(ilr,m,m)
     w0(n,m)=-s12(ilr-1,n,m)*p2(ilr,n,n)
     t3(n,m)= s22(ilr-1,n,m)*p2(ilr,m,m)
    END DO
!!    t0(n,n)=1+t0(n,n)
!!    w0(n,n)=1+w0(n,n)
    t0(n,n) = 1.0d0 + t0(n,n)
    w0(n,n) = 1.0d0 + w0(n,n)
   END DO
   CALL dlinrg(nang,t0,nang,t1,nang)   ! t1= (1-s12*r)**(-1)
   CALL dlinrg(nang,w0,nang,w1,nang)   ! w1= (1-r*s12)**(-1)
   t2=MATMUL(t1,s11(ilr-1,:,:))
   w2=MATMUL(s22(ilr-1,:,:),w1)
   w3=MATMUL(s12(ilr-1,:,:),w1)


!---
!   DO n=1,nang
!    DO m=1,nang
!     s11(ilr,n,m)=p1(ilr,n,n)*t2(n,m)
!     s22(ilr,n,m)=p1(ilr,m,m)*w2(n,m)
!     s12(ilr,n,m)=p1(ilr,n,n)*w3(n,m)*p1(ilr,m,m)
!    END DO
!     s12(ilr,n,n)=s12(ilr,n,n)+p2(ilr,n,n)
!   END DO
!   s21(ilr,:,:)=s21(ilr-1,:,:)+MATMUL(t3,t2)
!---

   DO n=1,nang
    DO m=1,nang
     s11(ilr,n,m)=p1(ilr,n,n)*t2(n,m)
    END DO
   END DO

    s22(ilr,:,:)=TRANSPOSE(s11(ilr,:,:))

   DO n=1,nang
    DO m=n,nang
     s12(ilr,n,m)=p1(ilr,n,n)*w3(n,m)*p1(ilr,m,m)
     s12(ilr,m,n)=s12(ilr,n,m)
     s1=s21(ilr-1,n,m)
     DO k=1,nang
      s1=s1+t3(n,k)*t2(k,m)
     END DO
     s21(ilr,n,m)=s1
     s21(ilr,m,n)=s1
    END DO
     s12(ilr,n,n)=s12(ilr,n,n)+p2(ilr,n,n)
   END DO

!----------------------------------
   ELSE     !  diagonal s_d(ilr-1)
!----------------------------------
! blw 911
   s11(ilr,:,:) = 0.0d0
   s22(ilr,:,:) = 0.0d0
   s12(ilr,:,:) = 0.0d0
   s21(ilr,:,:) = 0.0d0
! blw 911
    DO n=1,nang
     s2=1-s12(ilr-1,n,n)*p2(ilr,n,n)
     s1=p1(ilr,n,n)/s2
     s11(ilr,n,n)=s1*s11(ilr-1,n,n)
     s22(ilr,n,n)=s1*s22(ilr-1,n,n)
     s12(ilr,n,n)=s1*s12(ilr-1,n,n)*p1(ilr,n,n)+p2(ilr,n,n)
     s21(ilr,n,n)=s21(ilr-1,n,n)+p2(ilr,n,n)/s2*s11(ilr-1,n,n)      &
                                               *s22(ilr-1,n,n)
    END DO
!------------------------------------------------------------
  END IF    ! end of the non-diagonal/diagonal s_d(ilr-1)
!------------------------------------------------------------
!------------------------------------------------------------
  END IF                   ! end of the diagonal case: part 7
!------------------------------------------------------------

!---------------------------------------------------------------------
END DO                ! end of the  first half of the forth _ilr loop
!---------------------------------------------------------------------

 IF(oblv > 0) THEN
  IF(.NOT.LP(oblv,3)) THEN
   s11(oblv-1,:,:)=p1(oblv,:,:)
   s22(oblv-1,:,:)=p1(oblv,:,:)
   s12(oblv-1,:,:)=p2(oblv,:,:)
   s21(oblv-1,:,:)=p2(oblv,:,:)
  ELSE
! blw 911
  s11(oblv-1,:,:) = 0.0d0
  s12(oblv-1,:,:) = 0.0d0
  s21(oblv-1,:,:) = 0.0d0
  s22(oblv-1,:,:) = 0.0d0
! blw 911
   DO n=1,nang
    s11(oblv-1,n,n)=p1(oblv,n,n)
    s22(oblv-1,n,n)=p1(oblv,n,n)
    s12(oblv-1,n,n)=p2(oblv,n,n)
    s21(oblv-1,n,n)=p2(oblv,n,n)
   END DO
  END IF
  s_d(oblv-1)=LP(oblv,3)
 END IF

!-----------------------------------------------------------------------
DO ilr=oblv-1,1,-1      !        the second half of the forth _ilr loop
!-----------------------------------------------------------------------

  s_d(ilr-1)=s_d(ilr).AND.LP(ilr,3)

!------------------------------------------------------------
  IF( .NOT.LP(ilr,3) ) THEN    !    non-diagonal case: part 8
!------------------------------------------------------------

!------------------------------------------------------------
   IF( .NOT. s_d(ilr) ) THEN    !   non-diagonal s_d(ilr)
!------------------------------------------------------------
 t0=-MATMUL(s21(ilr,:,:), p2(ilr,:,:))
 w0=-MATMUL( p2(ilr,:,:),s21(ilr,:,:))
 DO n=1,nang
  t0(n,n)=t0(n,n)+1.d0
  w0(n,n)=w0(n,n)+1.d0
 END DO
 CALL dlinrg(nang,t0,nang,t1,nang)   ! t1= (1-s21*r)**(-1)
 CALL dlinrg(nang,w0,nang,w1,nang)   ! w1= (1-r*s21)**(-1)
 t2=MATMUL(t1,s22(ilr,:,:))
 w2=MATMUL(w1, p1(ilr,:,:))

!---
!  s11(ilr-1,:,:)=MATMUL(s11(ilr,:,:),w2)
!  s22(ilr-1,:,:)=MATMUL( p1(ilr,:,:),t2)
!  s12(ilr-1,:,:)=s12(ilr,:,:)+MATMUL( MATMUL(s11(ilr,:,:), p2(ilr,:,:)),t2)
!  s21(ilr-1,:,:)= p2(ilr,:,:)+MATMUL( MATMUL( p1(ilr,:,:),s21(ilr,:,:)),w2)
!---

   s11(ilr-1,:,:)=   MATMUL(s11(ilr  ,:,:),w2)
   s22(ilr-1,:,:)=TRANSPOSE(s11(ilr-1,:,:))

   t3=MATMUL(s11(ilr,:,:), p2(ilr,:,:))
   w3=MATMUL( p1(ilr,:,:),s21(ilr,:,:))

   DO n=1,nang
    DO m=n,nang
     s1=s12(ilr,n,m)
     s2= p2(ilr,n,m)
     DO k=1,nang
      s1=s1+t3(n,k)*t2(k,m)
      s2=s2+w3(n,k)*w2(k,m)
     END DO
     s12(ilr-1,n,m)=s1
     s21(ilr-1,n,m)=s2
     s12(ilr-1,m,n)=s1
     s21(ilr-1,m,n)=s2      
    END DO
   END DO


!----------------------------------
   ELSE     !  diagonal s_d(ilr)
!----------------------------------
    DO n=1,nang
     DO m=1,nang
      t0(n,m)=-s21(ilr,n,n)*p2(ilr,n,m)
      w0(n,m)=-s21(ilr,m,m)*p2(ilr,n,m)
     END DO
!!     t0(n,n)=t0(n,n)+1
!!     w0(n,n)=w0(n,n)+1
     t0(n,n) = t0(n,n) + 1.0d0
     w0(n,n) = w0(n,n) + 1.0d0
    END DO
   CALL dlinrg(nang,t0,nang,t1,nang)   ! t1= (1-s21*r)**(-1)
   CALL dlinrg(nang,w0,nang,w1,nang)   ! w1= (1-r*s21)**(-1)
   w2=MATMUL(w1,p1(ilr,:,:))
   t2=MATMUL(p1(ilr,:,:),t1)
   t3=MATMUL(p2(ilr,:,:),t1)

!----
!   DO n=1,nang
!    DO m=1,nang
!     s11(ilr-1,n,m)=s11(ilr,n,n)*w2(n,m)
!     s22(ilr-1,n,m)=             t2(n,m)*s22(ilr,m,m)
!     s12(ilr-1,n,m)=s11(ilr,n,n)*t3(n,m)*s22(ilr,m,m)
!
!     s1=p2(ilr,n,m)
!     DO k=1,nang
!      s1=s1+p1(ilr,n,k)*s21(ilr,k,k)*w2(k,m)
!     END DO
!     s21(ilr-1,n,m)=s1
!    END DO
!    s12(ilr-1,n,n)=s12(ilr-1,n,n)+s12(ilr,n,n)
!   END DO
!----

   DO n=1,nang
    DO m=1,nang
     s11(ilr-1,n,m)=s11(ilr,n,n)*w2(n,m)    
    END DO
   END DO

   s22(ilr-1,:,:)=TRANSPOSE(s11(ilr-1,:,:))

   DO n=1,nang
    DO m=n,nang
     s12(ilr-1,n,m)=s11(ilr,n,n)*t3(n,m)*s22(ilr,m,m)
     s12(ilr-1,m,n)=s12(ilr-1,n,m)
     s1=p2(ilr,n,m)
     DO k=1,nang
      s1=s1+p1(ilr,n,k)*s21(ilr,k,k)*w2(k,m)
     END DO
     s21(ilr-1,n,m)=s1
     s21(ilr-1,m,n)=s1
    END DO
    s12(ilr-1,n,n)=s12(ilr-1,n,n)+s12(ilr,n,n)
   END DO

!------------------------------------------------------------
  END IF    ! end of the non-diagonal/diagonal s_d(ilr)
!------------------------------------------------------------
!------------------------------------------------------------
  END IF               ! end of the non-diagonal case: part 8
!------------------------------------------------------------
!------------------------------------------------------------
  IF( LP(ilr,3) ) THEN           !      diagonal case: part 8
!------------------------------------------------------------
!------------------------------------------------------------
   IF( .NOT. s_d(ilr) ) THEN    !   non-diagonal s_d(ilr)
!------------------------------------------------------------
 DO n=1,nang
    DO m=1,nang
     t0(n,m)=-s21(ilr,n,m)*p2(ilr,m,m)
     w0(n,m)=-s21(ilr,n,m)*p2(ilr,n,n)
     t3(n,m)= s11(ilr,n,m)*p2(ilr,m,m)
    END DO
!!    t0(n,n)=1+t0(n,n)
!!    w0(n,n)=1+w0(n,n)
    t0(n,n) = 1.0d0 + t0(n,n)
    w0(n,n) = 1.0d0 + w0(n,n)
   END DO
   CALL dlinrg(nang,t0,nang,t1,nang)   ! t1= (1-s21*r)**(-1)
   CALL dlinrg(nang,w0,nang,w1,nang)   ! w1= (1-r*s21)**(-1)
   t2=MATMUL(t1,s22(ilr,:,:))
   w2=MATMUL(s11(ilr,:,:),w1)
   w3=MATMUL(s21(ilr,:,:),w1)

!--
!   DO n=1,nang
!    DO m=1,nang
!     s11(ilr-1,n,m)=p1(ilr,m,m)*w2(n,m)
!     s22(ilr-1,n,m)=p1(ilr,n,n)*t2(n,m)
!     s21(ilr-1,n,m)=p1(ilr,n,n)*w3(n,m)*p1(ilr,m,m)
!    END DO
!     s21(ilr,n,n)=s21(ilr,n,n)+p2(ilr,n,n)
!   END DO
!   s12(ilr-1,:,:)=s12(ilr,:,:)+MATMUL(t3,t2)
!--

   DO n=1,nang
    DO m=1,nang
     s11(ilr-1,n,m)=p1(ilr,m,m)*w2(n,m)
    END DO
   END DO

   s22(ilr-1,:,:)=TRANSPOSE(s11(ilr-1,:,:))

   DO n=1,nang
    DO m=n,nang
     s21(ilr-1,n,m)=p1(ilr,n,n)*w3(n,m)*p1(ilr,m,m)
     s21(ilr-1,m,n)=s21(ilr-1,n,m)
     s1=s12(ilr,n,m)
     DO k=1,nang
      s1=s1+t3(n,k)*t2(k,m)
     END DO
     s12(ilr-1,n,m)=s1
     s12(ilr-1,m,n)=s1
    END DO
     s21(ilr,n,n)=s21(ilr,n,n)+p2(ilr,n,n)
   END DO

!----------------------------------
   ELSE     !  diagonal s_d(ilr)
!----------------------------------
! blw 911
s11(ilr-1,:,:) = 0.0d0
s12(ilr-1,:,:) = 0.0d0
s21(ilr-1,:,:) = 0.0d0
s22(ilr-1,:,:) = 0.0d0
! blw 911
    DO n=1,nang
     s2=1-s21(ilr,n,n)*p2(ilr,n,n)
     s1=p1(ilr,n,n)/s2
     s11(ilr-1,n,n)=s1*s11(ilr,n,n)
     s22(ilr-1,n,n)=s1*s22(ilr,n,n)
     s21(ilr-1,n,n)=s1*s21(ilr,n,n)*p1(ilr,n,n)+p2(ilr,n,n)
     s12(ilr-1,n,n)=s12(ilr,n,n)+p2(ilr,n,n)/s2*s11(ilr,n,n)        &
                                               *s22(ilr,n,n)
    END DO
!------------------------------------------------------------
  END IF    ! end of the non-diagonal/diagonal s_d(ilr)
!------------------------------------------------------------
!-----------------------------------------------------------
  END IF                  ! end of the diagonal case: part 8
!-----------------------------------------------------------

!---------------------------------------------------------------------
END DO                ! end of the second half of the forth _ilr loop
!---------------------------------------------------------------------

! OPEN(103,file='dtb31.tst',FORM='FORMATTED')

!----------------------------------------------------------------------
 DO ilr=nlr,1,-1   !      the fifth _ilr loop
!----------------------------------------------------------------------
   
   a11=0.d0   
   a12=0.d0
   a21=0.d0
   a22=0.d0
  da11=0.d0
  da12=0.d0
  da21=0.d0
  da22=0.d0

   a_d=LP(ilr,3).AND.rup_d(ilr-1).AND.rdn_d(ilr)

!------------ calculation of  a11, a12, a21, a22, b1, b2                !kz, eqn. (88)
!------------            and da11,da12,da21,da22,db1,db2 
!-----------------------------------------------------------
  IF( .NOT. LP(ilr,3) ) THEN     ! non-diagonal case: part 9
!-----------------------------------------------------------

  IF ( rdn_d(ilr  ) ) THEN
    DO n=1,nang
     DO m=1,nang
      a11(n,m)=-p2(ilr,n,m)*rdn(ilr  ,m,m)
      a21(n,m)=-p1(ilr,n,m)*rdn(ilr  ,m,m)
     END DO
    END DO
  ELSE
    a11=-MATMUL(p2(ilr,:,:),rdn(ilr  ,:,:))
    a21=-MATMUL(p1(ilr,:,:),rdn(ilr  ,:,:))
  END IF

  IF ( rup_d(ilr-1) ) THEN
    DO n=1,nang
     DO m=1,nang
      a12(n,m)=-p1(ilr,n,m)*rup(ilr-1,m,m)
      a22(n,m)=-p2(ilr,n,m)*rup(ilr-1,m,m)
     END DO
    END DO
  ELSE
    a12=-MATMUL(p1(ilr,:,:),rup(ilr-1,:,:))
    a22=-MATMUL(p2(ilr,:,:),rup(ilr-1,:,:))
  END IF

! a11=-MATMUL(p2(ilr,:,:),rdn(ilr  ,:,:))
! a12=-MATMUL(p1(ilr,:,:),rup(ilr-1,:,:))
! a21=-MATMUL(p1(ilr,:,:),rdn(ilr  ,:,:))
! a22=-MATMUL(p2(ilr,:,:),rup(ilr-1,:,:))

!-----------------------------------------------------------
  END IF              ! end of the non-diagonal case: part 9
!-----------------------------------------------------------
!-----------------------------------------------------------
  IF( LP(ilr,3) ) THEN          !      diagonal case: part 9
!-----------------------------------------------------------

  IF ( rdn_d(ilr  ) ) THEN
! blw 911
  a11 = 0.0d0
  a21 = 0.0d0
! blw 911
    DO n=1,nang
      a11(n,n)=-p2(ilr,n,n)*rdn(ilr  ,n,n)
      a21(n,n)=-p1(ilr,n,n)*rdn(ilr  ,n,n)
    END DO
  ELSE
    DO n=1,nang
     DO m=1,nang
      a11(n,m)=-p2(ilr,n,n)*rdn(ilr  ,n,m)
      a21(n,m)=-p1(ilr,n,n)*rdn(ilr  ,n,m)
     END DO
    END DO
  END IF


  IF ( rup_d(ilr-1) ) THEN
! blw 91
  a12 = 0.0d0
  a22 = 0.0d0
! blw 911
    DO n=1,nang
      a12(n,n)=-p1(ilr,n,n)*rup(ilr-1,n,n)
      a22(n,n)=-p2(ilr,n,n)*rup(ilr-1,n,n)
    END DO
  ELSE
    DO n=1,nang
     DO m=1,nang
      a12(n,m)=-p1(ilr,n,n)*rup(ilr-1,n,m)
      a22(n,m)=-p2(ilr,n,n)*rup(ilr-1,n,m)
     END DO
    END DO
  END IF

!-----------------------------------------------------------
  END IF                  ! end of the diagonal case: part 9
!-----------------------------------------------------------

   DO n=1,nang
    a11(n,n)=a11(n,n)+1.d0
    a22(n,n)=a22(n,n)+1.d0
   END DO

   DO ivar=1,nvar
!-----------------------------------------------------------------------
    IF( LP(ilr,3) .AND. LP(ilr,ivar+3) ) THEN  ! diagonal layer and
                                               ! diagonal perturbation
!-----------------------------------------------------------------------

     IF( rdn_d(ilr).AND.rup_d(ilr-1) ) THEN
! blw 911
  da11(:,:,ivar) = 0.0d0
  da12(:,:,ivar) = 0.0d0
  da21(:,:,ivar) = 0.0d0
  da22(:,:,ivar) = 0.0d0
! blw 911
       DO n=1,nang
        da11(n,n,ivar)=-dp2(ilr,n,n,ivar)*rdn(ilr  ,n,n)
        da12(n,n,ivar)=-dp1(ilr,n,n,ivar)*rup(ilr-1,n,n)
        da21(n,n,ivar)=-dp1(ilr,n,n,ivar)*rdn(ilr  ,n,n)
        da22(n,n,ivar)=-dp2(ilr,n,n,ivar)*rup(ilr-1,n,n)
       END DO
     ELSE
       DO n=1,nang
        DO m=1,nang
        da11(n,m,ivar)=-dp2(ilr,n,n,ivar)*rdn(ilr  ,n,m)
        da12(n,m,ivar)=-dp1(ilr,n,n,ivar)*rup(ilr-1,n,m)
        da21(n,m,ivar)=-dp1(ilr,n,n,ivar)*rdn(ilr  ,n,m)
        da22(n,m,ivar)=-dp2(ilr,n,n,ivar)*rup(ilr-1,n,m)
        END DO
       END DO
     END IF

!-----------------------------------------------------------------------
    ELSE                         !  non-diagonal layer or perturbation
!-----------------------------------------------------------------------

     IF( rdn_d(ilr).AND.rup_d(ilr-1) ) THEN
       DO n=1,nang
        DO m=1,nang
        da11(n,m,ivar)=-dp2(ilr,n,m,ivar)*rdn(ilr  ,m,m)
        da12(n,m,ivar)=-dp1(ilr,n,m,ivar)*rup(ilr-1,m,m)
        da21(n,m,ivar)=-dp1(ilr,n,m,ivar)*rdn(ilr  ,m,m)
        da22(n,m,ivar)=-dp2(ilr,n,m,ivar)*rup(ilr-1,m,m)
        END DO
       END DO
     ELSE
      da11(:,:,ivar)=-MATMUL(dp2(ilr,:,:,ivar),rdn(ilr  ,:,:))
      da12(:,:,ivar)=-MATMUL(dp1(ilr,:,:,ivar),rup(ilr-1,:,:))
      da21(:,:,ivar)=-MATMUL(dp1(ilr,:,:,ivar),rdn(ilr  ,:,:))
      da22(:,:,ivar)=-MATMUL(dp2(ilr,:,:,ivar),rup(ilr-1,:,:))
     END IF
    
!-----------------------------------------------------------------------
    END IF           ! end of non-diagonal/diagonal layer/perturbation
!-----------------------------------------------------------------------
   END DO

    IF( rdn_d(ilr).AND.rup_d(ilr-1) ) THEN
     DO n=1,nang
      y1(n)= ust(ilr-1,n) + rup(ilr-1,n,n)*unh(ilr,n)
      y2(n)= vst(ilr  ,n) + rdn(ilr  ,n,n)*unh(ilr,n)
       DO ivar=1,nvar
        dy1v(n,ivar)= rup(ilr-1,n,n)*dunh(ilr,n,ivar)
        dy2v(n,ivar)= rdn(ilr  ,n,n)*dunh(ilr,n,ivar)
       END DO
      END DO
    ELSE
     y1= ust(ilr-1,:) + MATMUL(rup(ilr-1,:,:), unh(ilr,:))
     y2= vst(ilr  ,:) + MATMUL(rdn(ilr  ,:,:), unh(ilr,:))
      DO ivar=1,nvar
       dy1v(:,ivar)= MATMUL(rup(ilr-1,:,:), dunh(ilr,:,ivar))
       dy2v(:,ivar)= MATMUL(rdn(ilr  ,:,:), dunh(ilr,:,ivar))
      END DO
    END IF

!-----------------------------------------------------------
  IF( .NOT. LP(ilr,3) ) THEN    ! non-diagonal case: part 10
!-----------------------------------------------------------
   b1=MATMUL(p1(ilr,:,:),y1) + MATMUL(p2(ilr,:,:),y2)
   b2=MATMUL(p2(ilr,:,:),y1) + MATMUL(p1(ilr,:,:),y2)
!-----------------------------------------------------------
  END IF           !   end of the non-diagonal case: part 10
!-----------------------------------------------------------
!-----------------------------------------------------------
  IF( LP(ilr,3) ) THEN              ! diagonal case: part 10
!-----------------------------------------------------------
   DO n=1,nang
    b1(n)=p1(ilr,n,n)*y1(n)+p2(ilr,n,n)*y2(n)
    b2(n)=p2(ilr,n,n)*y1(n)+p1(ilr,n,n)*y2(n)
   END DO
!-----------------------------------------------------------
  END IF            !      end of the diagonal case: part 10
!-----------------------------------------------------------

!-------------------------------------------------------------
  IF( .NOT. LP(ilr,3) ) THEN    ! non-diagonal case: part 11
!-------------------------------------------------------------
   DO ivar=1,nvar
    db1(:,ivar)=  MATMUL(dp1(ilr,:,:,ivar), y1    )                    &
                + MATMUL(dp2(ilr,:,:,ivar), y2    )                    &
                + MATMUL( p1(ilr,:,:),dy1v(:,ivar))                    &
                + MATMUL( p2(ilr,:,:),dy2v(:,ivar))
    db2(:,ivar)=  MATMUL(dp2(ilr,:,:,ivar), y1    )                    &
                + MATMUL(dp1(ilr,:,:,ivar), y2    )                    &
                + MATMUL( p2(ilr,:,:),dy1v(:,ivar))                    &
                + MATMUL( p1(ilr,:,:),dy2v(:,ivar))
   END DO

!-----------------------------------------------------------
  END IF            ! end of the non-diagonal case: part 11
!-----------------------------------------------------------
!-----------------------------------------------------------
  IF( LP(ilr,3) ) THEN             ! diagonal case: part 11
!-----------------------------------------------------------
  DO ivar=1,nvar
!-----------------------------------------------------------------------
    IF( LP(ilr,3) .AND. LP(ilr,ivar+3) ) THEN  ! diagonal layer and
                                               ! diagonal perturbation
!-----------------------------------------------------------------------
      DO n=1,nang
       db1(n,ivar)=dp1(ilr,n,n,ivar)*y1(n)+dp2(ilr,n,n,ivar)*y2(n)
       db2(n,ivar)=dp2(ilr,n,n,ivar)*y1(n)+dp1(ilr,n,n,ivar)*y2(n)
      END DO
!-----------------------------------------------------------------------
    ELSE                         !  non-diagonal layer or perturbation
!-----------------------------------------------------------------------
    db1(:,ivar)=  MATMUL(dp1(ilr,:,:,ivar), y1    )                    &
                + MATMUL(dp2(ilr,:,:,ivar), y2    )                   
               
    db2(:,ivar)=  MATMUL(dp2(ilr,:,:,ivar), y1    )                    &
                + MATMUL(dp1(ilr,:,:,ivar), y2    )   
!-----------------------------------------------------------------------
    END IF           ! end of non-diagonal/diagonal layer/perturbation
!-----------------------------------------------------------------------                 
    DO n=1,nang
     db1(n,ivar)=db1(n,ivar)+p1(ilr,n,n)*dy1v(n,ivar)+p2(ilr,n,n)*dy2v(n,ivar)
     db2(n,ivar)=db2(n,ivar)+p2(ilr,n,n)*dy1v(n,ivar)+p1(ilr,n,n)*dy2v(n,ivar)
    END DO
   END DO
!-----------------------------------------------------------
  END IF               ! end of the diagonal case: part 11
!-----------------------------------------------------------

  IF( a_d ) THEN
   a11m1=0.d0
   a22m1=0.d0
      w0=0.d0
      w1=0.d0 
   DO n=1,nang
    a11m1(n,n)=1/a11(n,n)
    a22m1(n,n)=1/a22(n,n)
       w0(n,n)=a21(n,n)*a11m1(n,n)
       w1(n,n)=a12(n,n)*a22m1(n,n)
         y0(n)=b2(n)-w0(n,n)*b1(n)
         y1(n)=b1(n)-w1(n,n)*b2(n)
   END DO
  ELSE
   CALL dlinrg(nang,a11,nang,a11m1,nang)
   CALL dlinrg(nang,a22,nang,a22m1,nang)
   w0=MATMUL(a21,a11m1)
   w1=MATMUL(a12,a22m1)
   t0=a22-MATMUL(w0,a12)
   t1=a11-MATMUL(w1,a21)
   y0=b2 -MATMUL(w0, b1)
   y1=b1 -MATMUL(w1, b2)
  END IF

   DO ivar=1,nvar

    dt0(:,:,ivar)=0.d0
    dt1(:,:,ivar)=0.d0

    IF( a_d.AND.LP(ilr,ivar+3) ) THEN  ! all diag

     DO n=1,nang
  s1=da21(n,n,ivar)*a11m1(n,n)-a21(n,n)*da11(n,n,ivar)*a11m1(n,n)**2 ! =dw0
  s2=da12(n,n,ivar)*a22m1(n,n)-a12(n,n)*da22(n,n,ivar)*a22m1(n,n)**2 ! =dw1

      dy0v(n,ivar)=db2(n,ivar)-s1*b1(n)-w0(n,n)*db1(n,ivar)
      dy1v(n,ivar)=db1(n,ivar)-s2*b2(n)-w1(n,n)*db2(n,ivar)
       
      dt0(n,n,ivar)=da22(n,n,ivar)-s1*a12(n,n)-w0(n,n)*da12(n,n,ivar)
      dt1(n,n,ivar)=da11(n,n,ivar)-s2*a21(n,n)-w1(n,n)*da21(n,n,ivar)
     END DO

    ELSE    ! a_d .OR. LP(ilr,ivar+3)  - non diag

     da11m1(:,:,ivar)=-MATMUL( MATMUL(a11m1,da11(:,:,ivar)), a11m1 )
     da22m1(:,:,ivar)=-MATMUL( MATMUL(a22m1,da22(:,:,ivar)), a22m1 )

     dw0(:,:,ivar) = MATMUL(da21(:,:,ivar), a11m1          )              &
                    +MATMUL( a21          ,da11m1(:,:,ivar))
     dw1(:,:,ivar) = MATMUL(da12(:,:,ivar), a22m1          )              &
                    +MATMUL( a12          ,da22m1(:,:,ivar))

     dy0v(:,ivar)=db2(:,ivar) - MATMUL(dw0(:,:,ivar), b1        )         &
                              - MATMUL( w0          ,db1(:,ivar)) 
     dy1v(:,ivar)=db1(:,ivar) - MATMUL(dw1(:,:,ivar), b2        )         &
                              - MATMUL( w1          ,db2(:,ivar))
     dt0 (:,:,ivar)=da22(:,:,ivar) - MATMUL(dw0(:,:,ivar), a12(:,:     )) &
                                   - MATMUL( w0          ,da12(:,:,ivar))
     dt1 (:,:,ivar)=da11(:,:,ivar) - MATMUL(dw1(:,:,ivar), a21(:,:     )) &
                                   - MATMUL( w1          ,da21(:,:,ivar))
 
    END IF   ! non-diagonal/diagonal a_d and perturbation

   END DO

  IF( a_d ) THEN
   w0=0.d0
   w1=0.d0
   DO n=1,nang
!!    w0(n,n)=1/(a22(n,n)-a21(n,n)*a12(n,n)/a11(n,n))
!!    w1(n,n)=1/(a11(n,n)-a21(n,n)*a12(n,n)/a22(n,n))
    w0(n,n) = 1.0d0 /(a22(n,n)-a21(n,n)*a12(n,n)/a11(n,n))
    w1(n,n) = 1.0d0 /(a11(n,n)-a21(n,n)*a12(n,n)/a22(n,n))
   END DO
  ELSE
   CALL dlinrg(nang,t0,nang,w0,nang)  ! w0=(a22-a21*a11**(-1)*a12)**(-1)
   CALL dlinrg(nang,t1,nang,w1,nang)  ! w1=(a11-a12*a22**(-1)*a21)**(-1)
  END IF

!---------------------------------------------------------------------
!  This mode of calculation of brightness temperature can be used
!  if variation of the properties of the layer takes place and
!  matrices _w0,_w1 and vectors _y0, _y1 are known anyway. Otherwise,
!  the calculation should proceed according to the second _ilr loop.
!---------------------------------------------------------------------

 IF( ilr==nlr ) v(nlr  ,:)=f(nlr+1,:)

 IF(.NOT. a_d) THEN      !   non-diagonal _a => non-diagonal _w0, _w1
      u(ilr  ,:)=unh(ilr,:)+MATMUL(w1,y1)
      v(ilr-1,:)=unh(ilr,:)+MATMUL(w0,y0)
 ELSE                    !      diagonal _a =>      diagonal _w0, _w1
     DO n=1,nang
        u(ilr  ,n)=unh(ilr,n)+w1(n,n)*y1(n)
        v(ilr-1,n)=unh(ilr,n)+w0(n,n)*y0(n)
     END DO
 END IF

 IF( ilr==1   ) THEN
     IF(.NOT. rup_d(0) ) THEN      !   non-diagonal  _rup_d(0)
         u(0,:)=ust(0,:)+MATMUL(rup(0,:,:),v(0,:))
     ELSE                          !       diagonal  _rup_d(0)
        DO n=1,nang
         u(0,n)=ust(0,n)+rup(0,n,n)*v(0,n)
        END DO
     END IF
  END IF

!-------------------

DO ivar=1,nvar
!-------------------------------------------------
!   it should be y3(n)=u(ilr,n), dy3(n)=du(ilr,n)
!-------------------------------------------------

  IF( .NOT. a_d ) THEN             ! non-diagonal _w0, _w1

!   dz2(ilr-1,:,ivar)=dunh(ilr,:,ivar) + MATMUL(w0,dy0v(:,ivar))         &
!                        -MATMUL(MATMUL( MATMUL(w0,dt0(:,:,ivar)),w0),y0) 
!   dz3(ilr  ,:,ivar)=dunh(ilr,:,ivar) + MATMUL(w1,dy1v(:,ivar))         &
!                        -MATMUL(MATMUL( MATMUL(w1,dt1(:,:,ivar)),w1),y1) 

   dz2(ilr-1,:,ivar)=dunh(ilr,:,ivar)                                   &
            +MATMUL(w0,dy0v(:,ivar)-MATMUL(dt0(:,:,ivar),MATMUL(w0,y0))) 

   dz3(ilr  ,:,ivar)=dunh(ilr,:,ivar)                                   &
            +MATMUL(w1,dy1v(:,ivar)-MATMUL(dt1(:,:,ivar),MATMUL(w1,y1))) 

  ELSE                             !    diagonal _w0, _w1

   IF ( .NOT. LP(ilr,ivar+3) ) THEN   ! non-diagonal perturbation

     DO n=1,nang
      s1=dy0v(n,ivar)
      s2=dy1v(n,ivar)
      DO m=1,nang
       s1=s1-dt0(n,m,ivar)*w0(m,m)*y0(m)
       s2=s2-dt1(n,m,ivar)*w1(m,m)*y1(m)
      END DO
      dz2(ilr-1,n,ivar)=dunh(ilr,n,ivar)+w0(n,n)*s1
      dz3(ilr  ,n,ivar)=dunh(ilr,n,ivar)+w1(n,n)*s2
     END DO

   ELSE                               !     diagonal perturbation

     DO n=1,nang
      dz2(ilr-1,n,ivar)= dunh(ilr,n,ivar)                                &
                    +w0(n,n)*(dy0v(n,ivar)-dt0(n,n,ivar)*w0(n,n)*y0(n))
      dz3(ilr  ,n,ivar)= dunh(ilr,n,ivar)                                &
                    +w1(n,n)*(dy1v(n,ivar)-dt1(n,n,ivar)*w1(n,n)*y1(n))
     END DO

   END IF                             ! non-diagonal/diagonal perturbation

  END IF                              ! non-diagonal/diagonal   _w0, _w1

END DO   ! ivar

!-------------------------------------------------------------------------------
 END DO            !  end of the fifth _ilr loop
!-------------------------------------------------------------------------------


  DO ilr=1,oblv    !  lower part

   IF( .NOT. s_d(ilr) ) THEN    ! non-diagonal matrices s11,..,s22

    IF( .NOT. rdn_d(oblv) ) THEN    ! non-diagonal matrix rdn(oblv)

     t0=-MATMUL(s12(ilr,:,:),rdn(oblv,:,:))
     DO n=1,nang
      t0(n,n)=t0(n,n)+1.d0
     END DO

     CALL dlinrg(nang,t0,nang,t1,nang)  ! t1=(1-S12*R)**(-1)

     DO ivar=1,nvar
      du0(ilr,:,ivar)=MATMUL(t1,MATMUL(s11(ilr,:,:),dz3(ilr,:,ivar)))
      dv0(ilr,:,ivar)=MATMUL(rdn(oblv,:,:),du0(ilr,:,ivar))
     END DO   ! ivar

    ELSE                         !     diagonal matrix rdn(oblv)

     DO n=1,nang
      DO m=1,nang
       t0(n,m)=-s12(ilr,n,m)*rdn(oblv,m,m)
      END DO
      t0(n,n)=t0(n,n)+1.d0
     END DO

     CALL dlinrg(nang,t0,nang,t1,nang)  ! t1=(1-S12*R)**(-1)

     DO ivar=1,nvar
      du0(ilr,:,ivar)=MATMUL(t1,MATMUL(s11(ilr,:,:),dz3(ilr,:,ivar)))
!      dv0(ilr,:,ivar)=MATMUL(rdn(oblv,:,:),du0(ilr,:,ivar))
      DO n=1,nang
        dv0(ilr,n,ivar)=rdn(oblv,n,n)*du0(ilr,n,ivar)
      END DO
     END DO   ! ivar

    END IF                      ! non-diagonal/diagonal matrix rdn(oblv)

   ELSE                         !     diagonal matrices s11,..,s22

    IF( .NOT. rdn_d(oblv) ) THEN    ! non-diagonal matrix rdn_d(oblv)

    DO n=1,nang
     DO m=1,nang
      t0(n,m)=-s12(ilr,n,n)*rdn(oblv,n,m)
     END DO
     t0(n,n)=t0(n,n)+1.d0
    END DO

    CALL dlinrg(nang,t0,nang,t1,nang)  ! t1=(1-S12*R)**(-1)

    DO ivar=1,nvar
     DO n=1,nang
      s1=0.d0
      DO m=1,nang
       s1=s1+t1(n,m)*s11(ilr,m,m)*dz3(ilr,m,ivar)
      END DO
      du0(ilr,n,ivar)=s1
     END DO
     dv0(ilr,:,ivar)=MATMUL(rdn(oblv,:,:),du0(ilr,:,ivar))     
    END DO

    ELSE     !              diagonal matrix rdn_d(oblv)

     DO ivar=1,nvar
      DO n=1,nang
       du0(ilr,n,ivar)=s11(ilr,n,n)/(1-s12(ilr,n,n)*rdn(oblv,n,n))*dz3(ilr,n,ivar)
       dv0(ilr,n,ivar)=rdn(oblv,n,n)*du0(ilr,n,ivar) 
      END DO    
     END DO

    END IF   ! non-diagonal/diagonal matrix rdn_d(oblv)

   END IF          !     non-diagonal/diagonal matrices s11,..,s22


  END DO   ! end of the lower part


!--------------------------------------------------------------------------


!      du0(ilr,n,ivar)=ds1  !   du0 is variation of upwelling radiation at
!                           !  level _obs_lev in respond to variation at 
!                           !  level _ilr


!--------------------------------------------------------------------------

  DO ilr=oblv+1,nlr    !  upper part

   IF( .NOT. s_d(ilr-1) ) THEN    ! non-diagonal matrices s11,..,s22

    IF( .NOT. rup_d(oblv) ) THEN    ! non-diagonal matrix rup(oblv)

     w0=-MATMUL(s21(ilr-1,:,:),rup(oblv,:,:))
     DO n=1,nang
      w0(n,n)=w0(n,n)+1.d0
     END DO

     CALL dlinrg(nang,w0,nang,w1,nang)  ! w1=(1-S21*R)**(-1)

     DO ivar=1,nvar
       dv0(ilr,:,ivar)=MATMUL(w1,MATMUL(s22(ilr-1,:,:),dz2(ilr-1,:,ivar)))
       du0(ilr,:,ivar)=MATMUL(rup(oblv,:,:),dv0(ilr,:,ivar))
     END DO   !  ivar

    ELSE                           !     diagonal matrix rup(oblv)

     DO n=1,nang
      DO m=1,nang
       w0(n,m)=-s21(ilr-1,n,m)*rup(oblv,m,m)
      END DO
      w0(n,n)=w0(n,n)+1.d0
     END DO

     CALL dlinrg(nang,w0,nang,w1,nang)  ! w1=(1-S21*R)**(-1)

     DO ivar=1,nvar
      dv0(ilr,:,ivar)=MATMUL(w1,MATMUL(s22(ilr-1,:,:),dz2(ilr-1,:,ivar)))
      DO n=1,nang
        du0(ilr,n,ivar)=rup(oblv,n,n)*dv0(ilr,n,ivar)
      END DO
     END DO   ! ivar

    END IF                ! non-diagonal/diagonal matrix rup(oblv)

   ELSE                         !     diagonal matrices s11,..,s22

    IF( .NOT. rup_d(oblv) ) THEN    ! non-diagonal matrix rup_d(oblv)

    DO n=1,nang
     DO m=1,nang
      w0(n,m)=-s21(ilr-1,n,n)*rup(oblv,n,m)
     END DO
     w0(n,n)=w0(n,n)+1.d0
    END DO

    CALL dlinrg(nang,w0,nang,w1,nang)  ! w1=(1-S21*R)**(-1)

    DO ivar=1,nvar
     DO n=1,nang
      s1=0.d0
      DO m=1,nang
       s1=s1+w1(n,m)*s22(ilr-1,m,m)*dz2(ilr-1,m,ivar)
      END DO
      dv0(ilr,n,ivar)=s1
     END DO
     du0(ilr,:,ivar)=MATMUL(rup(oblv,:,:),dv0(ilr,:,ivar))     
    END DO

    ELSE     !              diagonal matrix rup(oblv)

     DO ivar=1,nvar
      DO n=1,nang
       dv0(ilr,n,ivar)= s22(ilr-1,n,n)/(1-s21(ilr-1,n,n)       &
                       *rup(oblv,n,n))*dz2(ilr-1,n,ivar)
       du0(ilr,n,ivar)=rup(oblv,n,n)*dv0(ilr,n,ivar) 
      END DO    
     END DO

    END IF   ! non-diagonal/diagonal matrix rup_d(oblv)

   END IF          !     non-diagonal/diagonal matrices s11,..,s22

  END DO   ! end of the upper part

!----------------------------------------------------------

deallocate(a,b,am1,p1,p2,q,stat=alloc_err)
deallocate(da,db,dp1,dp2,stat=alloc_err)

END SUBROUTINE core95




SUBROUTINE md_inv32 (nang,h, lama, ma, lamab, mab, lama12d, p1, p2               &
                         ,dlama,dma,dlamab,dmab,dp1,dp2,nvar,L)

!---------------------------------------------------------------------
!
!  This subroutine inverts a matrix of a special form:
!
!   t(n,m)= SUM [i,j=1 to i,j=nang] ma(n,i)*w(i,j)*ma(m,j)                      !kz, eqn. (51)
!
!   where
!
!   w(n,m)=SUM [k=1 to k=nang] a(n,k)*b(m,k)/dzeta(k)-2*nang*delta(n,m)
!
!   and 
!
!   dzeta(k)= SQRT(lamab(k))/sinh(SQRT(lamab(k)*h))                              !kz, eqn. (54)
!
!   a(n,k)=[SQRT(lama(n))+SQRT(lamab(k))/tanh(SQRT(lamab(k)*h/2))      &         !kz, eqn. (52)
!                                       /SQRT(lama(n))] * mab(n,k)
!
!   b(m,k)=[SQRT(lama(m))+SQRT(lamab(k))*tanh(SQRT(lamab(k)*h/2))      &         !kz, eqn. (53)
!                                       /SQRT(lama(m))] * mab(m,k)
!
!   and both _ma(n,k) and _mab(n,k) are orthogonal matrices.
!
!   It is assumed h > 0.
!
!   Parameter _nvar correspond to a vector of variations: 
!   the "derivatives" vectors have in fact one extra dimension.
!
!  Logical array L characterizes "diagonal properties" and
!  corresponds to the _ilr-th row of matrix LP
!----------------------------------------------------------------------
!---------------------------------------------------------------------


IMPLICIT NONE

INTEGER                                    , INTENT( IN) :: nang,nvar

DOUBLE PRECISION, DIMENSION(nang,nang)     , INTENT( IN) ::  ma,  mab   
DOUBLE PRECISION, DIMENSION(nang,nang,nvar), INTENT( IN) :: dma, dmab 
DOUBLE PRECISION, DIMENSION(nang     ), INTENT( IN) ::   lama,  lamab, lama12d  
DOUBLE PRECISION, DIMENSION(nang,nvar), INTENT( IN) ::  dlama, dlamab   

DOUBLE PRECISION                      , INTENT( IN) :: h
LOGICAL         , DIMENSION(nvar+3   ), INTENT( IN) :: L

DOUBLE PRECISION, DIMENSION(nang,nang)     , INTENT(OUT) ::  p1, p2
DOUBLE PRECISION, DIMENSION(nang,nang,nvar), INTENT(OUT) :: dp1,dp2

INTEGER          :: k,n,m,ivar
DOUBLE PRECISION :: g,x,em,t,lama12,u1,u2,u3,v1,v2,v3
DOUBLE PRECISION, DIMENSION(nvar) :: dg,dx,dem,dt

DOUBLE PRECISION, DIMENSION(nang)      ::  dzeta, ad,b1d,b2d,c1d,c2d
DOUBLE PRECISION, DIMENSION(nang,nvar) :: ddzeta
DOUBLE PRECISION, DIMENSION(nang,nang) :: a, b1, b2, am1, bm1, c1, c2, ac1, ac2  
DOUBLE PRECISION, DIMENSION(nang,nang,nvar) :: da,db1,db2,dam1,dbm1   &
                                              ,dc1,dc2
DOUBLE PRECISION :: s1,s2
!-------------------------------------------------------------

!-----------------------------
 DO  k=1,nang   ! loop over k
!-----------------------------

  g=DSQRT(lamab(k))
  x= g*h
  DO ivar=1,nvar
   dg(ivar)=dlamab(k,ivar)/2/g
   dx(ivar)=dg(ivar)*h
  END DO

  IF( x < 1.d-3) THEN
   dzeta(k)=1/h/(1+x**2/6*(1+x**2/20))
   DO ivar=1,nvar
    ddzeta(k,ivar)= -h*dzeta(k)**2*x/3*(1+x**2/10)*dx(ivar)
   END DO

          x= x/2
         dx=dx/2
          t=2/h*(1+x**2/3*(1-  x**2/15))
      DO ivar=1,nvar
         dt(ivar)=2/h* 2*x   /3*(1-2*x**2/15) *dx(ivar)
      END DO
   
  ELSE
         em=DEXP(-x)
       DO ivar=1,nvar
        dem(ivar)=-em*dx(ivar)
       END DO

   dzeta(k)=2*g*      em /(1-em*em)
           t= g*(1+em)/(1-em)  

   DO ivar=1,nvar
  ddzeta(k,ivar)=2/(1-em*em)*(g*(1+em*em)/(1-em*em)*dem(ivar)+dg(ivar)*em)
        dt(ivar)=1/(1-em)*(dg(ivar)*(1+em)+2*g/(1-em)*dem(ivar))
   END DO
           
  END IF


!-----------------------------------------------------------------
    IF( .NOT. L(3) ) THEN  ! non-diagonal layer
!-----------------------------------------------------------------
     DO n=1,nang
      lama12=lama12d(n)
      x=1/t/lama12
      u1=lama12+t/lama12
      u2=lama12+lamab(k)*x
      u3=lama12-lamab(k)*x
       a(n,k)=u1*mab(n,k)
      b1(n,k)=u2*mab(n,k)
      b2(n,k)=u3*mab(n,k)

      DO ivar=1,nvar
       v1=lamab(k)/lama(n)
       v2=dlamab(k,ivar)-lamab(k)*dt(ivar)/t
       v3=x*mab(n,k)
       da (n,k,ivar)=u1*dmab(n,k,ivar)                                  &
               +((1-t/lama(n))*dlama(n,ivar)/2+dt(ivar))/lama12*mab(n,k)
       db1(n,k,ivar)=u2*dmab(n,k,ivar)+((t-v1)*dlama(n,ivar)/2+v2)*v3
       db2(n,k,ivar)=u3*dmab(n,k,ivar)+((t+v1)*dlama(n,ivar)/2-v2)*v3 
      END DO   !  ivar
     END DO  !  (n)
!-----------------------------------------------------------------
   ELSE                             !     diagonal layer
!-----------------------------------------------------------------

     m=0
     DO ivar=1,nvar
      IF( .NOT. L(ivar+3) ) THEN
       m=1
       EXIT
      END IF
     END DO          ! if there are non-diagonal perturbations, m==1


!---  do this block only if there are non-diagonal perturbations ( m > 0 )

  IF( m > 0 ) THEN
!------------------       
     DO n=1,nang
      lama12=lama12d(n)
      x =1/t/lama12
      u1=lama12+t/lama12
      u2=lama12+lamab(k)*x
      u3=lama12-lamab(k)*x
       DO ivar=1,nvar
!-------------------------------------------------------------------
        IF( .NOT. L(ivar+3) ) THEN       ! non-diagonal perturbation
!-------------------------------------------------------------------
         da (n,k,ivar)=u1*dmab(n,k,ivar)                                    
         db1(n,k,ivar)=u2*dmab(n,k,ivar)
         db2(n,k,ivar)=u3*dmab(n,k,ivar)
!--------------------
        END IF
!--------------------
       END DO ! ivar
     END DO  ! n
!-----------------------------------------------------------
  END IF   ! block for non-diagonal perturbations  ( m > 0 )
!-----------------------------------------------------------

     lama12=lama12d(k)
         x =1/t/lama12
     ad (k)=lama12+t/lama12
     b1d(k)=lama12+lamab(k)*x
     b2d(k)=lama12-lamab(k)*x

     v1=lamab(k)/lama(k)

       DO ivar=1,nvar
        v2=dlamab(k,ivar)-lamab(k)*dt(ivar)/t
!-------------------------------------------------------------------
        IF( .NOT. L(ivar+3) ) THEN       ! non-diagonal perturbation
!-------------------------------------------------------------------
     da (k,k,ivar)=da (k,k,ivar)                                       &
                  +((1-t/lama(k))*dlama(k,ivar)/2+dt(ivar))/lama12
     db1(k,k,ivar)=db1(k,k,ivar)+((t-v1)*dlama(k,ivar)/2+v2)*x
     db2(k,k,ivar)=db2(k,k,ivar)+((t+v1)*dlama(k,ivar)/2-v2)*x
!--------------------------------------------------------------------
        ELSE                             !     diagonal perturbation
!--------------------------------------------------------------------
     da (k,k,ivar)=((1-t/lama(k))*dlama(k,ivar)/2+dt(ivar))/lama12
     db1(k,k,ivar)=((t-v1)*dlama(k,ivar)/2+v2)*x
     db2(k,k,ivar)=((t+v1)*dlama(k,ivar)/2-v2)*x
!-------------------------------------------------------------------
        END IF           ! end of diagonal/non-diagonal perturbation
!------------------------------------------------------------------- 
      END DO   !  ivar
!-----------------------------------------------------------------
   END IF           ! end of diagonal/non-diagonal layer
!-----------------------------------------------------------------

!---------------------------
 END DO  ! (k)
!---------------------------

!----------------------------------
 IF ( L(3) ) THEN  ! diagonal layer
!----------------------------------
  p1=0.d0          ! remove after debugging
  p2=0.d0          ! remove after debugging
  DO n=1,nang
   c1d(n)  = dzeta(n)/b1d(n)/ad(n)
   c2d(n)  =  -b2d(n)/b1d(n) 
   p1 (n,n)= 2*c1d(n)
   p2 (n,n)=   c2d(n)+p1(n,n)
  END DO
  DO ivar=1,nvar
!-----------------------------------------------------------------
    IF( .NOT. L(ivar+3) ) THEN  ! non-diagonal perturbation
!-----------------------------------------------------------------

    DO n=1,nang
     DO m=1,nang
      dam1(n,m,ivar)=- da(n,m,ivar)/ ad(n)/ ad(m)
      dbm1(n,m,ivar)=-db1(n,m,ivar)/b1d(n)/b1d(m)
     END DO
    END DO

    DO n=1,nang
     DO m=n,nang
      dp1(n,m,ivar)= dma(n,m,ivar)*c1d(m)+c1d(n)*dma(m,n,ivar)          &
                    +dbm1(m,n,ivar)*dzeta(m)/ad (m)                     &
                    +dam1(n,m,ivar)*dzeta(n)/b1d(n)
      dp2(n,m,ivar)= dma (n,m,ivar)*c2d(m)+c2d(n)*dma(m,n,ivar)         &
                    -dbm1(m,n,ivar)*b2d(m)-db2(m,n,ivar)/b1d(n)

      dp1(m,n,ivar)=dp1(n,m,ivar)
      dp2(m,n,ivar)=dp2(n,m,ivar)

     END DO
      dp1(n,n,ivar)=dp1(n,n,ivar)+ddzeta(n,ivar)/ad(n)/b1d(n)
    END DO

!-----------------------------------------------------------------
    ELSE                             !     diagonal perturbation
!-----------------------------------------------------------------
    DO n=1,nang
      s1=- da(n,n,ivar)/ ad(n)/ ad(n)
      s2=-db1(n,n,ivar)/b1d(n)/b1d(n)
      dp1(n,n,ivar)= ddzeta(n,ivar)/b1d(n)/ad(n)                        &
                     +dzeta(n)* (s1/b1d(n)+s2/ad(n))   
      dp2(n,n,ivar)=-s2*b2d(n)-db2(n,n,ivar)/b1d(n)
    END DO
!-----------------------------------------------------------------
    END IF           ! end of diagonal/non-diagonal perturbation
!-----------------------------------------------------------------
  END DO     ! (ivar)

  dp1=  2*dp1
  dp2=dp2+dp1

  RETURN

!-------------------------------
 END IF  ! end of diagonal layer
!-------------------------------


  CALL dlinrg(nang, a,nang,am1,nang)
  CALL dlinrg(nang,b1,nang,bm1,nang)


   DO ivar=1,nvar
    dam1(:,:,ivar)=-MATMUL(MATMUL(am1,da (:,:,ivar)),am1)
    dbm1(:,:,ivar)=-MATMUL(MATMUL(bm1,db1(:,:,ivar)),bm1)
   END DO
     DO n=1,nang
      DO m=n,nang
       s1=0.d0
       s2=0.d0
       DO k=1,nang
        s1=s1+bm1(k,n)*dzeta(k)*am1(k,m)
        s2=s2-bm1(k,n)         * b2(m,k)
       END DO
       c1(n,m)=s1
       c2(n,m)=s2
       c1(m,n)=s1
       c2(m,n)=s2
      END DO
     END DO

     DO ivar=1,nvar
      DO n=1,nang
       DO m=n,nang
        s1=0.d0
        s2=0.d0
        DO k=1,nang
   s1=s1+(dbm1(k,n,ivar)*dzeta(k)+bm1(k,n)*ddzeta(k,ivar))*am1(k,m)     &
        +  bm1(k,n)*dzeta(k)*dam1(k,m,ivar)      
   s2=s2-dbm1(k,n,ivar)*b2(m,k)-bm1(k,n)*db2(m,k,ivar)
        END DO
        dc1(n,m,ivar)=s1
        dc2(n,m,ivar)=s2
        dc1(m,n,ivar)=s1
        dc2(m,n,ivar)=s2
       END DO
      END DO
     END DO


    b1=MATMUL(ma,c1)
    b2=MATMUL(ma,c2)
    DO n=1,nang
     DO m=n,nang
      s1=0.d0
      s2=0.d0
      DO k=1,nang
       s1=s1+ b1(n,k)*ma(m,k)
       s2=s2+ b2(n,k)*ma(m,k)
      END DO
      p1(n,m)=2*s1
      p2(n,m)=  s2
      p1(m,n)=2*s1
      p2(m,n)=  s2
     END DO
    END DO

!   p1=2*MATMUL(MATMUL(ma,c1),TRANSPOSE(ma))
!   p2=  MATMUL(MATMUL(ma,c2),TRANSPOSE(ma))

   DO ivar=1,nvar

    ac1=MATMUL(dma(:,:,ivar), c1)
    ac2=MATMUL(dma(:,:,ivar), c2)
    am1=MATMUL( ma,dc1(:,:,ivar))
    bm1=MATMUL( ma,dc2(:,:,ivar))

    DO n=1,nang
     DO m=n,nang
      s1=0.d0
      s2=0.d0
      DO k=1,nang
       s1=s1+ ac1(n,k)*ma(m,k)+b1(n,k)*dma(m,k,ivar)+am1(n,k)*ma(m,k)
       s2=s2+ ac2(n,k)*ma(m,k)+b2(n,k)*dma(m,k,ivar)+bm1(n,k)*ma(m,k)
      END DO
      dp1(n,m,ivar)=2*s1
      dp2(n,m,ivar)=  s2
      dp1(m,n,ivar)=2*s1
      dp2(m,n,ivar)=  s2
     END DO
    END DO

!    dp1(:,:,ivar)=2*( MATMUL( MATMUL(dma(:,:,ivar),c1)                 &
!                             +MATMUL(ma,dc1(:,:,ivar)),TRANSPOSE(ma))  &
!                     +MATMUL( MATMUL(ma,c1),TRANSPOSE(dma(:,:,ivar))) )
!    dp2(:,:,ivar)=    MATMUL( MATMUL(dma(:,:,ivar),c2)                 &
!                             +MATMUL(ma,dc2(:,:,ivar)),TRANSPOSE(ma))  &
!                     +MATMUL( MATMUL(ma,c2),TRANSPOSE(dma(:,:,ivar))) 
   END DO

   p2= p2+ p1
  dp2=dp2+dp1

!-------------------------------------------------------------
END SUBROUTINE  md_inv32

subroutine dlinrg(n,a,nn,ainv,nnn)
! Organizes input to gaussj, calls gaussj
implicit none
integer*4 n, nn, nnn, np, m, mp, j
real*8, dimension(n,n) :: a, ainv, aa
real*8, dimension(n)   :: b
    np = n
    m = 1
    mp = 1
    aa = a
    do j = 1, n
        b(j) = 0.0d0
    end do
    call gaussj(aa,n,np,b,m,mp)
    ainv = aa
end subroutine dlinrg

