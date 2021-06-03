!====================================================================
 subroutine core95 (nlev,nang,thick,a0,b0,f,r,u,v, &
                    da0,db0,df,dtb_up_dp,dtb_dn_dp,obs_lev,nvar,diagonal)
!====================================================================
! This subroutine solves radiative transfer equation,
! Written and modified by Alex Voronovich,Alexander.Voronovich@noaa.gov
!
! Histoy:
! 1/31/2003 Logical matrix input argument diagonal(0:nlev,nvar+3)  (layers' profile)
! is added. Its _ilev row contains the following information:
! To calculate _u,_v  ;  To calculate _du,_dv  ; Is it diagonal layer;
! Is it diagonal layer with respect to _ivar-th parameter variation
! Diagonal layers are treated differently
! 2/14/2003 Diagonal perturbations are treated differently
! 3/4/2003 if _rup,_rdn operators are diagonal
! 04/28/2003 In the first loop over _ilev acceleration due to symmetry of some matrices
! is taken into account Subroutine md_inv9 => md_inv31    - includes argument lama12  
! 05/08/2003 Symmetries of matrices _t,_r (i.e.,p1 and p2),_rup,_rdn and
! the relations: S11^T=S22,S12^T=S12,S21^T=S21 are taken into account
! 5/15/2003 Subroutine CORE95.F90 is renamed subroutine DTB32.F90 from ACCEL6
! 9/26/2020 Kevin Schaefer deleted unused variables, arguments, and tabs
! 1/27/2021 Kevin Schaefer deleted commented and unused code
! 5/14/2021 Kevin Schaefer cleaned up variable definitions and code
! 5/15/2021 Kevin Schaefer changed Jacobian variable names
! 5/16/2021 Kevin Schaefer changed logical LP to diagonal
!
! This subroutine solves radiative transfer equation taken in the form:
! boundary condition at level n=N :    v_N (z=H) =  f_(N+1)
! withiin layer n=N :
!
! du_N/dz = -A0_N * u_N - B0_N * v_N + f_N
! dv_N/dz =  B0_N * u_N + A0_N * v_N - f_N
!
! within layer n=k :
!
! du_k/dz = -A0_k * u_k - B0_k * v_k + f_k
! dv_k/dz =  B0_k * u_k + A0_k * v_k - f_k
!
! within layer n=1 :
!
! du_1/dz = -A0_1 * u_1 - B0_1 * v_1 + f_1
! dv_1/dz =  B0_1 * u_1 + A0_1 * v_1 - f_1
!
! boundary condition at level n=0 :    u_0 (z=0) = R*v_0 (z=0) + f_0
! where R(i),i=1,..,nang is reflection coefficient.
!
! The solution is calculated by adding layer technique.
! Output values dtb_up_dp(ilev,n,nvar), dtb_dn_dp(ilev,n,nvar) are derivatives of
! the up- and downwelling  brightness temperatures at the observation
! level _obs_lev  corresponding to variations of the level properties
! at level _ilev. Here _obs_lev could be equal to  any number within
! interval
!                      0 <= obs_lev <= nrl
! and 1 <= ilev <= nlev. For future use dtb_up_dp, dtb_dn_dp are defined in fact
! for ilev=0 also (when reflectivity properties will be varied).
! However,in the present version of the code dtb_up_dp(0,n,ivar) and
! dtb_dn_dp(0,n,ivar) do not have meaning.
! Index _ivar corresponds to vector of different variations.
!----------------------------------------------------------------------
! use dotlrt_variables, only: testval ! for debugging only

IMPLICIT NONE

! Input variables
integer,intent(in) :: nlev         ! (-) number of layers
integer,intent(in) :: nang         ! (-) number quadrature angles
integer,intent(in) :: obs_lev      ! (-) observation layer
integer,intent(in) :: nvar         ! (-) number of variables
real(8),intent(in) :: thick(nlev)  ! (km) layer thickness
real(8),intent(in) :: f(0:nlev+1,nang)
real(8),intent(in) :: r(nang,nang)       ! (1/km) reflectivity
real(8),intent(in) :: df(0:nlev+1,nang,nvar)
real(8),intent(in) :: a0(nlev,nang,nang)
real(8),intent(in) :: b0(nlev,nang,nang)
real(8),intent(in) :: da0(nlev,nang,nang,nvar)
real(8),intent(in) :: db0(nlev,nang,nang,nvar)
logical,intent(in) :: diagonal(0:nlev,nvar+3)  ! (-) true if diagonal case

! output variables
real(8),intent(out) :: u(0:nlev,nang)        ! (K) upwelling brightness temperature at observation height
real(8),intent(out) :: v(0:nlev,nang)        ! (K) upwelling brightness temperature at observation height
real(8),intent(out) :: dtb_up_dp(0:nlev,nang,nvar) ! (K param) jacobian upwelling brightness temperature wrt parameter
real(8),intent(out) :: dtb_dn_dp(0:nlev,nang,nvar) ! (K param) jacobian downwelling brightness temperature wrt parameter

! local variables
integer*4 nrot
real(8),allocatable :: a(:,:,:)
real(8),allocatable :: b(:,:,:)
real(8),allocatable :: am1(:,:,:)
real(8),allocatable :: p1(:,:,:)
real(8),allocatable :: p2(:,:,:)
real(8),allocatable :: q(:,:,:)
real(8),allocatable :: da(:,:,:,:)
real(8),allocatable :: db(:,:,:,:)
real(8),allocatable :: dp1(:,:,:,:)
real(8),allocatable :: dp2(:,:,:,:)
real(8) rup(0:nlev,nang,nang)
real(8) rdn(0:nlev,nang,nang)
real(8) s11(0:nlev,nang,nang)
real(8) s12(0:nlev,nang,nang)
real(8) s21(0:nlev,nang,nang)
real(8) s22(0:nlev,nang,nang)
real(8) ust(0:nlev,nang)
real(8) vst(0:nlev,nang)
real(8) dz2(0:nlev,nang,nvar)
real(8) dz3(0:nlev,nang,nvar)
real(8) unh(nlev,nang)
real(8) dunh(nlev,nang,nvar)
real(8) y0(nang)
real(8) y1(nang)
real(8) y2(nang)
real(8) dy0(nang)
real(8) b1(nang)
real(8) b2(nang) 
real(8) t0(nang,nang)
real(8) t1(nang,nang)
real(8) t2(nang,nang)
real(8) w0(nang,nang)
real(8) w1(nang,nang)
real(8) ma(nang,nang)
real(8) mab(nang,nang)
real(8) a11(nang,nang)
real(8) a12(nang,nang)
real(8) a21(nang,nang)
real(8) a22(nang,nang)
real(8) a11m1(nang,nang)
real(8) a22m1(nang,nang)
real(8) w2(nang,nang)
real(8) t3(nang,nang)
real(8) w3(nang,nang)
real(8) lama12nm(nang,nang)
real(8) dt0(nang,nang,nvar)
real(8) dt1(nang,nang,nvar)
real(8) dt2(nang,nang,nvar)
real(8) dma(nang,nang,nvar)
real(8) dmab(nang,nang,nvar)
real(8) da11(nang,nang,nvar)
real(8) da12(nang,nang,nvar)
real(8) da21(nang,nang,nvar)
real(8) da22(nang,nang,nvar)
real(8) da11m1(nang,nang,nvar)
real(8) da22m1(nang,nang,nvar)
real(8) dw0(nang,nang,nvar)
real(8) dw1(nang,nang,nvar)
real(8) lama(nang)
real(8) lamab(nang)
real(8) lama12(nang)
real(8) dlama(nang,nvar)
real(8) dlamab(nang,nvar)
real(8) db1(nang,nvar)
real(8) db2(nang,nvar)
real(8) dy0v(nang,nvar)
real(8) dy1v(nang,nvar)
real(8) dy2v(nang,nvar)
integer ilev  ! (-) layer index
integer n
integer m
integer k
integer oblv  ! (-) local version of obs_lev
integer ivar  ! (-) variable index
real(8) s1
real(8) s2
logical q_d(nlev)
logical rup_d(0:nlev)
logical rdn_d(0:nlev)
logical s_d(0:nlev)
logical a_d
integer alloc_err

! allocate variables
  allocate(a(nlev,nang,nang))
  allocate(b(nlev,nang,nang))
  allocate(am1(nlev,nang,nang))
  allocate(p1(nlev,nang,nang))
  allocate(p2(nlev,nang,nang))
  allocate(q(nlev,nang,nang))
  allocate(da(nlev,nang,nang,nvar))
  allocate(db(nlev,nang,nang,nvar))
  allocate(dp1(nlev,nang,nang,nvar))
  allocate(dp2(nlev,nang,nang,nvar))

oblv=obs_lev
u(0,:)=0.d0
v(0,:)=0.d0
u(nlev,:)=0.d0
v(nlev,:)=0.d0

!kz,Symmetrization of matrices
  DO ilev=1,nlev
    IF(.NOT. diagonal(ilev,3)) THEN
      a(ilev,:,:) = a0(ilev,:,:)+b0(ilev,:,:)
      b(ilev,:,:) = a0(ilev,:,:)-b0(ilev,:,:)
      ! symmetrization of matrices _a,_b below - principally not needed
      a(ilev,:,:)=(a(ilev,:,:)+TRANSPOSE(a(ilev,:,:)))/2 
      b(ilev,:,:)=(b(ilev,:,:)+TRANSPOSE(b(ilev,:,:)))/2 
    END IF

    IF(diagonal(ilev,3)) THEN
      a(ilev,:,:) = 0.0d0
      b(ilev,:,:) = 0.0d0
      DO n=1,nang
        a(ilev,n,n) = a0(ilev,n,n)+b0(ilev,n,n) !kz,Eqn. 30
        b(ilev,n,n) = a0(ilev,n,n)-b0(ilev,n,n)
      END DO
    END IF

    DO ivar=1,nvar
      IF(.NOT. diagonal(ilev,3+ivar)) THEN
        da(ilev,:,:,ivar) = da0(ilev,:,:,ivar)+db0(ilev,:,:,ivar)
        db(ilev,:,:,ivar) = da0(ilev,:,:,ivar)-db0(ilev,:,:,ivar)
        ! symmetrization of matrices _da,_db below - principally not needed
        da(ilev,:,:,ivar)=(da(ilev,:,:,ivar)+TRANSPOSE(da(ilev,:,:,ivar)))/2 
        db(ilev,:,:,ivar)=(db(ilev,:,:,ivar)+TRANSPOSE(db(ilev,:,:,ivar)))/2
      END IF
      IF(diagonal(ilev,3+ivar)) THEN
        da(ilev,:,:,ivar) = 0.0d0
        db(ilev,:,:,ivar) = 0.0d0
        DO n=1,nang
          da(ilev,n,n,ivar) = da0(ilev,n,n,ivar)+db0(ilev,n,n,ivar)
          db(ilev,n,n,ivar) = da0(ilev,n,n,ivar)-db0(ilev,n,n,ivar)
        END DO
      END IF
    END DO
  END DO

  rup(0,:,:) = r(:,:)
  ust(0,:) = f(0,:)
  rup_d(0)=diagonal(0,3) ! TRUE if reflection matrix is diagonal

  DO ilev=1,nlev ! the first loop over _ilev
    rup_d(ilev)=rup_d(ilev-1) .AND. diagonal(ilev,3)
    q_d(ilev)=.FALSE.
    IF(.NOT. diagonal(ilev,3)) THEN ! non-diagonal case: part 0
      call jacobi(a(ilev,:,:),nang,nang,lama,ma,nrot)
    END IF ! end of non-diagonal case: part 0

    IF(diagonal(ilev,3)) THEN ! diagonal case: part 0
     DO n=1,nang
       lama(n)=a(ilev,n,n)
     END DO
    END IF ! end of diagonal case: part 0

    lama12=DSQRT(lama) !kz: square root of eigenvalues of A

    DO ivar=0,nvar
      IF(.NOT. diagonal(ilev,ivar+3)) THEN
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

!-------------------------------------------------------------
    IF(.NOT. diagonal(ilev,3)) THEN ! non-diagonal case: part 1
      DO ivar=1,nvar
        IF(.NOT. diagonal(ilev,ivar+3)) t3=MATMUL(da(ilev,:,:,ivar),ma)
        DO n=1,nang
          DO m=n,nang
            IF(.NOT. diagonal(ilev,ivar+3)) THEN ! non-diagonal perturbation
              s1=DOT_PRODUCT(ma(:,n),t3(:,m))
            ELSE ! diagonal perturbation
              s1=0.d0
              DO k=1,nang
                s1=s1+ma(k,n)*da(ilev,k,k,ivar)*ma(k,m)      
              END DO
            END IF ! end of diagonal/non-diagonal perturbation
            IF(n == m) THEN
              dlama(n,ivar)=s1
              t2(n,n)=0.d0
            ELSE
              t2(n,m)=s1/(lama(m)-lama(n))
              t2(m,n)=-t2(n,m)
            END IF
          END DO ! m
        END DO ! n
        dma(:,:,ivar)=MATMUL(ma,t2) !kz,Eqn. 94
      END DO ! ivar

      t3=MATMUL(TRANSPOSE(ma),b(ilev,:,:))

      DO n=1,nang
        DO m=n,nang
          s1=0.d0
          s2=0.d0
          DO k=1,nang
            s1=s1+ma(n,k)/lama(k)*ma(m,k)
            s2=s2+t3(n,k)*ma(k,m)
          END DO
          am1(ilev,n,m)=s1
          am1(ilev,m,n)=s1
          t1(n,m)=s2
          t1(m,n)=s2
        END DO
      END DO

      t0=t1*lama12nm
      w2=t1/lama12nm/2

      DO ivar=1,nvar
        DO n=1,nang
          DO m=n,nang
               w3(n,m)=dlama(n,ivar)*lama(m)+lama(n)*dlama(m,ivar)
               w3(m,n)=w3(n,m)
          END DO
        END DO

        IF(.NOT. diagonal(ilev,ivar+3)) THEN ! non-diagonal perturbation
          t2= MATMUL(db(ilev,:,:,ivar),ma)
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
        ELSE ! diagonal perturbation
          DO n=1,nang
            DO m=n,nang
              s1=0.d0
              DO k=1,nang
                s1=s1+ma(k,n)*db(ilev,k,k,ivar)*ma(k,m)
              END DO
              t1(n,m)=s1
              t1(m,n)=s1
            END DO
          END DO
        END IF ! end of diagonal/non-diagonal perturbation
        t2=MATMUL(t3,dma(:,:,ivar))
        dt0(:,:,ivar)=lama12nm*(t2+TRANSPOSE(t2)+t1)+w3*w2           
      END DO ! ivar
    END IF ! end of non-diagonal case: part 1


!-------------------------------------------------------------
    IF(diagonal(ilev,3)) THEN ! diag case: part 1
      DO ivar=1,nvar
         IF(.NOT. diagonal(ilev,ivar+3)) THEN ! non-diagonal perturbation 
            DO n=1,nang 
               dlama(n,ivar)=da(ilev,n,n,ivar)
               dma(n,n,ivar)=0.d0
               DO m=n+1,nang
                  dma(n,m,ivar)= da(ilev,n,m,ivar)/(lama(m)-lama(n))
                  dma(m,n,ivar)=   -dma(n,m,ivar)
               END DO
            END DO
         ELSE ! diagonal perturbation
            DO n=1,nang
               dlama(n,ivar)=da(ilev,n,n,ivar)
            END DO
         END IF ! end of diagonal/non-diagonal perturbation
      END DO ! ivar

      t0 = 0.0d0
      DO n=1,nang
         t0(n,n)=b(ilev,n,n)*lama(n)
      END DO

      DO ivar=1,nvar
         IF(.NOT. diagonal(ilev,ivar+3)) THEN ! non-diagonal perturbation
            DO n=1,nang
               DO m=n,nang
                  dt0(n,m,ivar)=(dma(m,n,ivar)*b(ilev,m,m)+b(ilev,n,n)*dma(n,m,ivar) &
                     +db(ilev,n,m,ivar))*lama12nm(n,m)
                  dt0(m,n,ivar)=dt0(n,m,ivar)
               END DO
               dt0(n,n,ivar)=dt0(n,n,ivar)+dlama(n,ivar)*b(ilev,n,n)
            END DO
                  ELSE ! diagonal perturbation
            dt0(:,:,ivar) = 0.0d0
            DO n=1,nang
               dt0(n,n,ivar)=db(ilev,n,n,ivar)*lama(n)+dlama(n,ivar)*b(ilev,n,n)
            END DO
                  END IF ! end of diagonal/non-diagonal perturbation
               END DO ! ivar

   END IF ! end of diagonal case: part 1

!-------------------------------------------------------------
   IF(diagonal(ilev,3)) THEN ! diagonal case
      DO n=1,nang
         lamab(n)=t0(n,n) !kz,Eqn. 46??
      END DO
   ELSE ! non-diagonal case
      call jacobi(t0,nang,nang,lamab,mab,nrot)
   END IF ! diagonal / non-diagonal case

   IF(.NOT. diagonal(ilev,3)) THEN ! non-giagonal case: part 2
   DO ivar=1,nvar
    t1=MATMUL(dt0(:,:,ivar),mab)
    DO n=1,nang
     DO m=n,nang
      s1=DOT_PRODUCT(mab(:,n),t1(:,m))
      IF(n == m) THEN
       dlamab(n,ivar)=s1
       t2(n,n)=0.d0
      ELSE
       t2(n,m)=s1/(lamab(m)-lamab(n))
       t2(m,n)=-t2(n,m)
      END IF
     END DO
    END DO
    dmab(:,:,ivar)=MATMUL(mab,t2)
   END DO ! ivar
 END IF ! end of non-diagonal case: part 2

!-------------------------------------------------------------
 IF(diagonal(ilev,3)) THEN ! diagonal case: part 2
  DO ivar=1,nvar
   DO n=1,nang
    dlamab(n,ivar)=dt0(n,n,ivar)
   END DO
    IF(.NOT. diagonal(ilev,ivar+3)) THEN ! non-diagonal perturbation
     DO n=1,nang
      dmab(n,n,ivar)=0.d0
      DO m=n+1,nang
       dmab(n,m,ivar)= dt0(n,m,ivar)/(lamab(m)-lamab(n))
       dmab(m,n,ivar)=-dmab(n,m,ivar)
      END DO
     END DO
    END IF ! end of diagonal perturbation
  END DO
 END IF ! end of diag case: part 2

! calculation of _p1,_p2   (T,R operators)
CALL md_inv32 (nang,thick(ilev),lama,ma,lamab,mab,lama12,t1,t2, &
              dlama,dma,dlamab,dmab,dt1,dt2,nvar,diagonal(ilev,:))

!-------------------------------------------------------------
 IF(.NOT. diagonal(ilev,3)) THEN ! non-diagonal case: part 3
   p1(ilev,:,:)=t1
   p2(ilev,:,:)=t2
 ELSE ! diagonal case: part 3
  p1(ilev,:,:) = 0.0d0
  p2(ilev,:,:) = 0.0d0
   DO n=1,nang
     p1(ilev,n,n)= t1(n,n)
     p2(ilev,n,n)= t2(n,n)
   END DO
 END IF ! end of diagonal case: part 3

 DO ivar=1,nvar
    IF(diagonal(ilev,3) .AND. diagonal(ilev,ivar+3)) THEN ! diagonal layer and diagonal perturbation
  dp1(ilev,:,:,ivar) = 0.0d0
  dp2(ilev,:,:,ivar) = 0.0d0
   DO n=1,nang
    dp1(ilev,n,n,ivar)=dt1(n,n,ivar)
    dp2(ilev,n,n,ivar)=dt2(n,n,ivar)
   END DO
    ELSE ! non-diagonal layer or perturbation
   dp1(ilev,:,:,ivar)=dt1(:,:,ivar)
   dp2(ilev,:,:,ivar)=dt2(:,:,ivar)
    END IF ! end of non-diagonal/diagonal layer/perturbation
 END DO ! ivar

!-------------------------------------------------------------
! calculation of matrix RUP
 IF(.NOT. diagonal(ilev,3)) THEN ! non-diagonal case: part 4

  IF(.NOT. rup_d(ilev-1)) THEN ! non-diagonal _rup(ilev-1)
   t0=-MATMUL(p2(ilev,:,:),rup(ilev-1,:,:))
   t2= MATMUL(p1(ilev,:,:),rup(ilev-1,:,:))
  ELSE ! diagonal _rup(ilev-1)
   DO n=1,nang
    DO m=1,nang
     t0(n,m)=-p2(ilev,n,m)*rup(ilev-1,m,m)
     t2(n,m)= p1(ilev,n,m)*rup(ilev-1,m,m)
    END DO
   END DO
  END IF ! non-diagonal/diagonal _rup(ilev-1)

   DO n=1,nang
    t0(n,n)= t0(n,n)+1.d0
   END DO
   CALL dlinrg(nang,t0,t1)
     q(ilev,:,:)=t1

! rup(ilev,:,:)=p2(ilev,:,:)+MATMUL(t2,MATMUL(q(ilev,:,:),p1(ilev,:,:)))

   w0=MATMUL(q(ilev,:,:),p1(ilev,:,:))
   DO n=1,nang
    DO m=n,nang
     s1=p2(ilev,n,m)
     DO k=1,nang
      s1=s1+t2(n,k)*w0(k,m)
     END DO
     rup(ilev,n,m)=s1
     rup(ilev,m,n)=s1
    END DO
   END DO
 
 END IF ! end of non-diagonal case: part 4

!-------------------------------------------------------------
 IF(diagonal(ilev,3)) THEN ! diagonal case: part 4
  IF(.NOT. rup_d(ilev-1)) THEN ! non-diagonal rup(ilev-1)
   DO n=1,nang
     s1=p2(ilev,n,n)
    DO m=1,nang
     t0(n,m)=-s1*rup(ilev-1,n,m)
    END DO
    t0(n,n)=t0(n,n)+1.d0
   END DO

   CALL dlinrg(nang,t0,t1)
   q(ilev,:,:)=t1
   t2=MATMUL(rup(ilev-1,:,:),t1)
   DO n=1,nang
    DO m=n,nang
     rup(ilev,n,m)=p1(ilev,n,n)*t2(n,m)*p1(ilev,m,m)
     rup(ilev,m,n)=rup(ilev,n,m)
    END DO
    rup(ilev,n,n)=rup(ilev,n,n)+p2(ilev,n,n)
   END DO
  ELSE ! diagonal rup(ilev-1)
  q(ilev,:,:) = 0.0d0
  rup(ilev,:,:) = 0.0d0
   DO n=1,nang
    q  (ilev,n,n)=1/(1-p2(ilev,n,n)*rup(ilev-1,n,n))
    rup(ilev,n,n)= p1(ilev,n,n)*rup(ilev-1,n,n)*q(ilev,n,n)*p1(ilev,n,n) &
                 +p2(ilev,n,n)
   END DO
   q_d(ilev)=.TRUE.
 END IF ! non-diagonal/diagonal rup(ilev-1)

 END IF ! end of diagonal case: part 4

!-------------------------------------------------------------
    IF(.NOT. diagonal(ilev,3)) THEN ! non-diagonal case: part 5
      !if(ilev == 7) print*, 'ok1 dunh', ilev, diagonal(ilev,3)
      y0=MATMUL(am1(ilev,:,:),f(ilev,:))
      t1=p1(ilev,:,:)+p2(ilev,:,:)
      unh(ilev,:)= y0-MATMUL(t1,y0)
      DO ivar=1,nvar
        dy0=MATMUL(am1(ilev,:,:),df(ilev,:,ivar)-MATMUL(da(ilev,:,:,ivar),y0))
        dunh(ilev,:,ivar)=dy0-MATMUL(dp1(ilev,:,:,ivar)+dp2(ilev,:,:,ivar),y0)-MATMUL(t1,dy0)
      END DO

      IF(.NOT. rup_d(ilev-1)) THEN ! non-diagonal rup(ilev-1)
        y0=ust(ilev-1,:)+MATMUL(rup(ilev-1,:,:),unh(ilev,:))
        ust(ilev,:)=unh(ilev,:)+MATMUL(p1(ilev,:,:),y0+MATMUL(rup(ilev-1,:,:),&
        MATMUL(q (ilev,:,:),MATMUL(p2(ilev,:,:),y0))))
      ELSE ! diagonal rup(ilev-1)
        DO n=1,nang
          y0(n)=ust(ilev-1,n)+rup(ilev-1,n,n)*unh(ilev,n)
        END DO
        y1=MATMUL(q(ilev,:,:),MATMUL(p2(ilev,:,:),y0))
        DO n=1,nang
          y2(n)=y0(n)+rup(ilev-1,n,n)*y1(n)
        END DO
        ust(ilev,:)=unh(ilev,:)+MATMUL(p1(ilev,:,:),y2)
      END IF ! non-diagonal/diagonal rup(ilev-1)
    END IF ! end of non-diagonal case: part 5

!-------------------------------------------------------------
    IF(diagonal(ilev,3)) THEN ! diagonal case: part 5
      !if(ilev == 7) print*, 'ok2 dunh', ilev, diagonal(ilev,3)
      DO n=1,nang
        y0(n)=f(ilev,n)/a(ilev,n,n)
        unh(ilev,n)=(1-p1(ilev,n,n)-p2(ilev,n,n))*y0(n)
      END DO
      DO ivar=1,nvar
        IF(.NOT. diagonal(ilev,ivar+3)) THEN ! non-diagonal perturbation
          y1=df(ilev,:,ivar)-MATMUL(da(ilev,:,:,ivar),y0)
        ELSE ! diagonal perturbation
          DO n=1,nang
            y1(n)=df(ilev,n,ivar)-da(ilev,n,n,ivar)*y0(n)
          END DO
        END IF ! end of diagonal/non-diagonal perturbation
        DO n=1,nang
          dy0(n)=y1(n)/a(ilev,n,n)
        END DO
        IF(.NOT. diagonal(ilev,ivar+3)) THEN ! non-diagonal perturbation
          dunh(ilev,:,ivar)=dy0-MATMUL(dp1(ilev,:,:,ivar)+dp2(ilev,:,:,ivar),y0)
        ELSE ! diagonal perturbation
          DO n=1,nang
            dunh(ilev,n,ivar)=dy0(n)-(dp1(ilev,n,n,ivar)+dp2(ilev,n,n,ivar))*y0(n)
          END DO
        END IF ! end of diagonal/non-diagonal perturbation
        DO n=1,nang
          dunh(ilev,n,ivar)=dunh(ilev,n,ivar)-(p1(ilev,n,n)+p2(ilev,n,n))*dy0(n)
        END DO
      END DO ! ivar

      IF(.NOT. rup_d(ilev-1)) THEN ! non-diagonal rup(ilev-1)
        y0=ust(ilev-1,:)+MATMUL(rup(ilev-1,:,:),unh(ilev,:))
        DO n=1,nang
          y1(n)=p2(ilev,n,n)*y0(n)
        END DO
        y2=y0+MATMUL(rup(ilev-1,:,:),MATMUL(q(ilev,:,:),y1))
      ELSE ! diagonal rup(ilev-1)
        DO n=1,nang
          y0(n)=ust(ilev-1,n)+rup(ilev-1,n,n)*unh(ilev,n)
          y1(n)=p2(ilev,n,n)*y0(n)  
        END DO
        b1=MATMUL(q(ilev,:,:),y1)
        DO n=1,nang
          y2(n)=y0(n)+rup(ilev-1,n,n)*b1(n)  
        END DO
      END IF ! non-diagonal/diagonal rup(ilev-1)

      DO n=1,nang
        ust(ilev,n)=unh(ilev,n)+p1(ilev,n,n)*y2(n)
      END DO
    END IF ! end of diagonal case: part 5
  END DO ! end of the first _ilev loop

!------------------------------------------------------------------
! The second loop over _ilev calculates vectors _u and _v which
! are up- and downwelling brightness temperatures. They are also
! calculated in the 5-th loop over _ilev as a by-product of
! calculations of derivatives a little bit faster.
! This mode of calculation should be used when variation of the
! _ilev layer parameters is not performed,however brightness
! temperature at the top OR at the bottom of the layer still has
! to be calculated
! Loop 2 is not used in this version and was deleted
!------------------------------------------------------------------

 vst(nlev,:)=f(nlev+1,:)
 vst(nlev-1,:)=unh(nlev,:)+MATMUL(p1(nlev,:,:),vst(nlev,:))

 rdn=0.d0 ! could be deleted after debugging

 rdn(nlev,:,:)=0.d0
 rdn(nlev-1,:,:)=p2(nlev,:,:)

 rdn_d(nlev)=.TRUE.
 rdn_d(nlev-1)=diagonal(nlev,3)

 DO ilev=nlev-1,1,-1 ! the third loop: calculation of _rdn,_vst

   rdn_d(ilev-1)=rdn_d(ilev).AND.diagonal(ilev,3)

 IF(.NOT. diagonal(ilev,3)) THEN ! non-diagonal case: part 6

 IF(.NOT. rdn_d(ilev)) THEN ! rdn(ilev) - non-diagonal
  t0=-MATMUL(rdn(ilev,:,:),p2(ilev,:,:))
  DO n=1,nang
    t0(n,n) = 1.0d0 + t0(n,n)
  END DO

  CALL dlinrg(nang,t0,t1)
  t2=MATMUL(p1(ilev,:,:),t1)

! rdn(ilev-1,:,:)=p2(ilev,:,:)+MATMUL(MATMUL(t2,rdn(ilev,:,:)),p1(ilev,:,:))

   w0=MATMUL(rdn(ilev,:,:),p1(ilev,:,:))
   DO n=1,nang
    DO m=n,nang
     s1=p2(ilev,n,m)
     DO k=1,nang
      s1=s1+t2(n,k)*w0(k,m)
     END DO
     rdn(ilev-1,n,m)=s1
     rdn(ilev-1,m,n)=s1
    END DO
   END DO

  vst(ilev-1,:)=unh(ilev,:)+MATMUL(t2,vst(ilev,:)+MATMUL(rdn(ilev,:,:),unh(ilev,:)))
 ELSE ! rdn(ilev) - diagonal
  DO n=1,nang
   DO m=1,nang
    t0(n,m)=-rdn(ilev,n,n)*p2(ilev,n,m)
   END DO
   t0(n,n) = 1.0d0 + t0(n,n)
  END DO
  CALL dlinrg(nang,t0,t1)
  t2=MATMUL(p1(ilev,:,:),t1)
  DO n=1,nang
   DO m=n,nang
    s1=p2(ilev,n,m)
    DO k=1,nang
     s1=s1+t2(n,k)*rdn(ilev,k,k)*p1(ilev,k,m)
    END DO
    rdn(ilev-1,n,m)=s1
    rdn(ilev-1,m,n)=s1
   END DO
   y0(n)=vst(ilev,n)+rdn(ilev,n,n)*unh(ilev,n)
  END DO

  vst(ilev-1,:)=unh(ilev,:)+MATMUL(t2,y0)
 END IF ! rdn(ilev) - non-diagonal/diagonal
 END IF ! end of non-diagonal case: part 6

 IF(diagonal(ilev,3)) THEN ! diagonal case: part 6
 IF(.NOT. rdn_d(ilev)) THEN ! rdn(ilev) - non-diagonal
  DO n=1,nang
   DO m=1,nang
    t0(n,m)=-rdn(ilev,n,m)*p2(ilev,m,m)
   END DO
   t0(n,n) = 1.0d0 + t0(n,n)
  END DO
  CALL dlinrg(nang,t0,t1)
  DO n=1,nang
   DO m=1,nang
    t2(n,m)=p1(ilev,n,n)*t1(n,m)
   END DO
  END DO
  t1=MATMUL(t2,rdn(ilev,:,:))
  DO n=1,nang
   DO m=n,nang
    rdn(ilev-1,n,m)=t1(n,m)*p1(ilev,m,m)
    rdn(ilev-1,m,n)=rdn(ilev-1,n,m)
   END DO
    rdn(ilev-1,n,n)=rdn(ilev-1,n,n)+p2(ilev,n,n)
  END DO

  vst(ilev-1,:)=unh(ilev,:)+MATMUL(t2,vst(ilev,:)+MATMUL(rdn(ilev,:,:),unh(ilev,:)))
 ELSE ! rdn(ilev) - diagonal
  rdn(ilev-1,:,:) = 0.0d0
  DO n=1,nang
   s1=p1(ilev,n,n)/(1-rdn(ilev,n,n)*p2(ilev,n,n))
   rdn(ilev-1,n,n)=p2(ilev,n,n)+p1(ilev,n,n)*s1*rdn(ilev,n,n)
   vst(ilev-1,n)=unh(ilev,n)+s1*(vst(ilev,n)+rdn(ilev,n,n)*unh(ilev,n))
  END DO
 END IF ! rdn(ilev) - non-diagonal/diagonal
 END IF ! diagonal case: part 6

END DO ! end of the third loop over _ilev: calculation of _rdn,_vst

 oblv=MAX0(0,oblv)
 oblv=MIN0(nlev,oblv)

! the forth _ilev loop: calculation of s11,s12,s21,s22

 s11=0.d0 ! to remove after debugging
 s12=0.d0 ! to remove after debugging
 s21=0.d0 ! to remove after debugging
 s22=0.d0 ! to remove after debugging

 DO n=1,nang
  s11(oblv,n,n)=1.d0
  s22(oblv,n,n)=1.d0
 END DO
 s_d(oblv)=.TRUE.

 IF(oblv < nlev) THEN
  IF(.NOT.diagonal(oblv+1,3)) THEN
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
  s_d(oblv+1)=diagonal(oblv+1,3)
 END IF

DO ilev=oblv+2,nlev ! the first half of the forth _ilev loop: top layers

  s_d(ilev)=s_d(ilev-1).AND.diagonal(ilev,3)

  IF(.NOT. diagonal(ilev,3)) THEN ! non-diagonal case: part 7
   IF(.NOT. s_d(ilev-1)) THEN ! non-diagonal s_d(ilev-1)
 t0=-MATMUL(s12(ilev-1,:,:),p2(ilev,:,:))
 w0=-MATMUL(p2(ilev,:,:),s12(ilev-1,:,:))
 DO n=1,nang
  t0(n,n) = 1.0d0 + t0(n,n)
  w0(n,n) = 1.0d0 + w0(n,n)
 END DO
 CALL dlinrg(nang,t0,t1) ! t1= (1-s12*r)**(-1)
 CALL dlinrg(nang,w0,w1) ! w1= (1-r*s12)**(-1)
  t2=MATMUL(t1,s11(ilev-1,:,:))
  w2=MATMUL(w1,p1(ilev,:,:))

   s11(ilev,:,:)=MATMUL(p1(ilev,:,:),t2)
   s22(ilev,:,:)=TRANSPOSE(s11(ilev,:,:))

   w3=MATMUL(p1(ilev,:,:),s12(ilev-1,:,:))
   t3=MATMUL(s22(ilev-1,:,:),p2(ilev,:,:))

   DO n=1,nang
    DO m=n,nang
     s1=p2(ilev,n,m)
     s2=s21(ilev-1,n,m)
     DO k=1,nang
      s1=s1+w3(n,k)*w2(k,m)
      s2=s2+t3(n,k)*t2(k,m)
     END DO
     s12(ilev,n,m)=s1
     s21(ilev,n,m)=s2
     s12(ilev,m,n)=s1
     s21(ilev,m,n)=s2    
    END DO
   END DO

   ELSE ! diagonal s_d(ilev-1)
    DO n=1,nang
     DO m=1,nang
      t0(n,m)=-s12(ilev-1,n,n)*p2(ilev,n,m)
      w0(n,m)=-s12(ilev-1,m,m)*p2(ilev,n,m)
     END DO
     t0(n,n) = 1.0d0 + t0(n,n)
     w0(n,n) = 1.0d0 + w0(n,n)
    END DO
    CALL dlinrg(nang,t0,t1) ! t1= (1-s12*r)**(-1)
    CALL dlinrg(nang,w0,w1) ! w1= (1-r*s12)**(-1)
    t2=MATMUL(p1(ilev,:,:),t1)
    w2=MATMUL(w1,p1(ilev,:,:))
    t3=MATMUL(p2(ilev,:,:),t1)

    DO n=1,nang
     DO m=1,nang
      s11(ilev,n,m)=t2(n,m)*s11(ilev-1,m,m)
     END DO
    END DO

    s22(ilev,:,:)=TRANSPOSE(s11(ilev,:,:))

    DO n=1,nang
     DO m=n,nang     
      s21(ilev,n,m)=s22(ilev-1,n,n)*t3(n,m)*s11(ilev-1,m,m)
      s21(ilev,m,n)=s21(ilev,n,m)
       s1=p2(ilev,n,m)
       DO k=1,nang
        s1=s1+p1(ilev,n,k)*s12(ilev-1,k,k)*w2(k,m)
       END DO
      s12(ilev,n,m)=s1
      s12(ilev,m,n)=s1     
     END DO
     s21(ilev,n,n)=s21(ilev,n,n)+s21(ilev-1,n,n)
    END DO

  END IF ! end of the non-diagonal/diagonal s_d(ilev-1)

  END IF ! end of the non-diagonal case: part 7
  IF(diagonal(ilev,3)) THEN ! diagonal case: part 7

   IF(.NOT. s_d(ilev-1)) THEN ! non-diagonal s_d(ilev-1)
   DO n=1,nang
    DO m=1,nang
     t0(n,m)=-s12(ilev-1,n,m)*p2(ilev,m,m)
     w0(n,m)=-s12(ilev-1,n,m)*p2(ilev,n,n)
     t3(n,m)= s22(ilev-1,n,m)*p2(ilev,m,m)
    END DO
    t0(n,n) = 1.0d0 + t0(n,n)
    w0(n,n) = 1.0d0 + w0(n,n)
   END DO
   CALL dlinrg(nang,t0,t1) ! t1= (1-s12*r)**(-1)
   CALL dlinrg(nang,w0,w1) ! w1= (1-r*s12)**(-1)
   t2=MATMUL(t1,s11(ilev-1,:,:))
   w2=MATMUL(s22(ilev-1,:,:),w1)
   w3=MATMUL(s12(ilev-1,:,:),w1)

   DO n=1,nang
    DO m=1,nang
     s11(ilev,n,m)=p1(ilev,n,n)*t2(n,m)
    END DO
   END DO

    s22(ilev,:,:)=TRANSPOSE(s11(ilev,:,:))

   DO n=1,nang
    DO m=n,nang
     s12(ilev,n,m)=p1(ilev,n,n)*w3(n,m)*p1(ilev,m,m)
     s12(ilev,m,n)=s12(ilev,n,m)
     s1=s21(ilev-1,n,m)
     DO k=1,nang
      s1=s1+t3(n,k)*t2(k,m)
     END DO
     s21(ilev,n,m)=s1
     s21(ilev,m,n)=s1
    END DO
     s12(ilev,n,n)=s12(ilev,n,n)+p2(ilev,n,n)
   END DO

   ELSE ! diagonal s_d(ilev-1)
   s11(ilev,:,:) = 0.0d0
   s22(ilev,:,:) = 0.0d0
   s12(ilev,:,:) = 0.0d0
   s21(ilev,:,:) = 0.0d0
    DO n=1,nang
     s2=1-s12(ilev-1,n,n)*p2(ilev,n,n)
     s1=p1(ilev,n,n)/s2
     s11(ilev,n,n)=s1*s11(ilev-1,n,n)
     s22(ilev,n,n)=s1*s22(ilev-1,n,n)
     s12(ilev,n,n)=s1*s12(ilev-1,n,n)*p1(ilev,n,n)+p2(ilev,n,n)
     s21(ilev,n,n)=s21(ilev-1,n,n)+p2(ilev,n,n)/s2*s11(ilev-1,n,n)*s22(ilev-1,n,n)
    END DO
  END IF ! end of the non-diagonal/diagonal s_d(ilev-1)
  END IF ! end of the diagonal case: part 7

END DO ! end of the  first half of the forth _ilev loop

 IF(oblv > 0) THEN
  IF(.NOT.diagonal(oblv,3)) THEN
   s11(oblv-1,:,:)=p1(oblv,:,:)
   s22(oblv-1,:,:)=p1(oblv,:,:)
   s12(oblv-1,:,:)=p2(oblv,:,:)
   s21(oblv-1,:,:)=p2(oblv,:,:)
  ELSE
  s11(oblv-1,:,:) = 0.0d0
  s12(oblv-1,:,:) = 0.0d0
  s21(oblv-1,:,:) = 0.0d0
  s22(oblv-1,:,:) = 0.0d0
   DO n=1,nang
    s11(oblv-1,n,n)=p1(oblv,n,n)
    s22(oblv-1,n,n)=p1(oblv,n,n)
    s12(oblv-1,n,n)=p2(oblv,n,n)
    s21(oblv-1,n,n)=p2(oblv,n,n)
   END DO
  END IF
  s_d(oblv-1)=diagonal(oblv,3)
 END IF

DO ilev=oblv-1,1,-1 ! the second half of the forth _ilev loop
  s_d(ilev-1)=s_d(ilev).AND.diagonal(ilev,3)

  IF(.NOT.diagonal(ilev,3)) THEN ! non-diagonal case: part 8

   IF(.NOT. s_d(ilev)) THEN ! non-diagonal s_d(ilev)
 t0=-MATMUL(s21(ilev,:,:),p2(ilev,:,:))
 w0=-MATMUL(p2(ilev,:,:),s21(ilev,:,:))
 DO n=1,nang
  t0(n,n)=t0(n,n)+1.d0
  w0(n,n)=w0(n,n)+1.d0
 END DO
 CALL dlinrg(nang,t0,t1) ! t1= (1-s21*r)**(-1)
 CALL dlinrg(nang,w0,w1) ! w1= (1-r*s21)**(-1)
 t2=MATMUL(t1,s22(ilev,:,:))
 w2=MATMUL(w1,p1(ilev,:,:))

   s11(ilev-1,:,:)=   MATMUL(s11(ilev,:,:),w2)
   s22(ilev-1,:,:)=TRANSPOSE(s11(ilev-1,:,:))

   t3=MATMUL(s11(ilev,:,:),p2(ilev,:,:))
   w3=MATMUL(p1(ilev,:,:),s21(ilev,:,:))

   DO n=1,nang
    DO m=n,nang
     s1=s12(ilev,n,m)
     s2= p2(ilev,n,m)
     DO k=1,nang
      s1=s1+t3(n,k)*t2(k,m)
      s2=s2+w3(n,k)*w2(k,m)
     END DO
     s12(ilev-1,n,m)=s1
     s21(ilev-1,n,m)=s2
     s12(ilev-1,m,n)=s1
     s21(ilev-1,m,n)=s2    
    END DO
   END DO


   ELSE ! diagonal s_d(ilev)
    DO n=1,nang
     DO m=1,nang
      t0(n,m)=-s21(ilev,n,n)*p2(ilev,n,m)
      w0(n,m)=-s21(ilev,m,m)*p2(ilev,n,m)
     END DO
     t0(n,n) = t0(n,n) + 1.0d0
     w0(n,n) = w0(n,n) + 1.0d0
    END DO
   CALL dlinrg(nang,t0,t1) ! t1= (1-s21*r)**(-1)
   CALL dlinrg(nang,w0,w1) ! w1= (1-r*s21)**(-1)
   w2=MATMUL(w1,p1(ilev,:,:))
   t2=MATMUL(p1(ilev,:,:),t1)
   t3=MATMUL(p2(ilev,:,:),t1)

   DO n=1,nang
    DO m=1,nang
     s11(ilev-1,n,m)=s11(ilev,n,n)*w2(n,m)  
    END DO
   END DO

   s22(ilev-1,:,:)=TRANSPOSE(s11(ilev-1,:,:))

   DO n=1,nang
    DO m=n,nang
     s12(ilev-1,n,m)=s11(ilev,n,n)*t3(n,m)*s22(ilev,m,m)
     s12(ilev-1,m,n)=s12(ilev-1,n,m)
     s1=p2(ilev,n,m)
     DO k=1,nang
      s1=s1+p1(ilev,n,k)*s21(ilev,k,k)*w2(k,m)
     END DO
     s21(ilev-1,n,m)=s1
     s21(ilev-1,m,n)=s1
    END DO
    s12(ilev-1,n,n)=s12(ilev-1,n,n)+s12(ilev,n,n)
   END DO

  END IF ! end of the non-diagonal/diagonal s_d(ilev)
  END IF ! end of the non-diagonal case: part 8
  IF(diagonal(ilev,3)) THEN ! diagonal case: part 8
   IF(.NOT. s_d(ilev)) THEN ! non-diagonal s_d(ilev)
 DO n=1,nang
    DO m=1,nang
     t0(n,m)=-s21(ilev,n,m)*p2(ilev,m,m)
     w0(n,m)=-s21(ilev,n,m)*p2(ilev,n,n)
     t3(n,m)= s11(ilev,n,m)*p2(ilev,m,m)
    END DO
    t0(n,n) = 1.0d0 + t0(n,n)
    w0(n,n) = 1.0d0 + w0(n,n)
   END DO
   CALL dlinrg(nang,t0,t1) ! t1= (1-s21*r)**(-1)
   CALL dlinrg(nang,w0,w1) ! w1= (1-r*s21)**(-1)
   t2=MATMUL(t1,s22(ilev,:,:))
   w2=MATMUL(s11(ilev,:,:),w1)
   w3=MATMUL(s21(ilev,:,:),w1)

   DO n=1,nang
    DO m=1,nang
     s11(ilev-1,n,m)=p1(ilev,m,m)*w2(n,m)
    END DO
   END DO

   s22(ilev-1,:,:)=TRANSPOSE(s11(ilev-1,:,:))

   DO n=1,nang
    DO m=n,nang
     s21(ilev-1,n,m)=p1(ilev,n,n)*w3(n,m)*p1(ilev,m,m)
     s21(ilev-1,m,n)=s21(ilev-1,n,m)
     s1=s12(ilev,n,m)
     DO k=1,nang
      s1=s1+t3(n,k)*t2(k,m)
     END DO
     s12(ilev-1,n,m)=s1
     s12(ilev-1,m,n)=s1
    END DO
     s21(ilev,n,n)=s21(ilev,n,n)+p2(ilev,n,n)
   END DO

   ELSE ! diagonal s_d(ilev)
s11(ilev-1,:,:) = 0.0d0
s12(ilev-1,:,:) = 0.0d0
s21(ilev-1,:,:) = 0.0d0
s22(ilev-1,:,:) = 0.0d0
    DO n=1,nang
     s2=1-s21(ilev,n,n)*p2(ilev,n,n)
     s1=p1(ilev,n,n)/s2
     s11(ilev-1,n,n)=s1*s11(ilev,n,n)
     s22(ilev-1,n,n)=s1*s22(ilev,n,n)
     s21(ilev-1,n,n)=s1*s21(ilev,n,n)*p1(ilev,n,n)+p2(ilev,n,n)
     s12(ilev-1,n,n)=s12(ilev,n,n)+p2(ilev,n,n)/s2*s11(ilev,n,n)*s22(ilev,n,n)
    END DO
  END IF ! end of the non-diagonal/diagonal s_d(ilev)
  END IF ! end of the diagonal case: part 8

END DO ! end of the second half of the forth _ilev loop

 DO ilev=nlev,1,-1 ! the fifth _ilev loop
 
   a11=0.d0 
   a12=0.d0
   a21=0.d0
   a22=0.d0
  da11=0.d0
  da12=0.d0
  da21=0.d0
  da22=0.d0

   a_d=diagonal(ilev,3).AND.rup_d(ilev-1).AND.rdn_d(ilev)

! calculation of  a11,a12,a21,a22,b1,b2 !kz,eqn. (88)
! and da11,da12,da21,da22,db1,db2
  IF(.NOT. diagonal(ilev,3)) THEN ! non-diagonal case: part 9

  IF (rdn_d(ilev)) THEN
    DO n=1,nang
     DO m=1,nang
      a11(n,m)=-p2(ilev,n,m)*rdn(ilev,m,m)
      a21(n,m)=-p1(ilev,n,m)*rdn(ilev,m,m)
     END DO
    END DO
  ELSE
    a11=-MATMUL(p2(ilev,:,:),rdn(ilev,:,:))
    a21=-MATMUL(p1(ilev,:,:),rdn(ilev,:,:))
  END IF

  IF (rup_d(ilev-1)) THEN
    DO n=1,nang
     DO m=1,nang
      a12(n,m)=-p1(ilev,n,m)*rup(ilev-1,m,m)
      a22(n,m)=-p2(ilev,n,m)*rup(ilev-1,m,m)
     END DO
    END DO
  ELSE
    a12=-MATMUL(p1(ilev,:,:),rup(ilev-1,:,:))
    a22=-MATMUL(p2(ilev,:,:),rup(ilev-1,:,:))
  END IF

  END IF ! end of the non-diagonal case: part 9
  IF(diagonal(ilev,3)) THEN ! diagonal case: part 9

  IF (rdn_d(ilev)) THEN
  a11 = 0.0d0
  a21 = 0.0d0
    DO n=1,nang
      a11(n,n)=-p2(ilev,n,n)*rdn(ilev,n,n)
      a21(n,n)=-p1(ilev,n,n)*rdn(ilev,n,n)
    END DO
  ELSE
    DO n=1,nang
     DO m=1,nang
      a11(n,m)=-p2(ilev,n,n)*rdn(ilev,n,m)
      a21(n,m)=-p1(ilev,n,n)*rdn(ilev,n,m)
     END DO
    END DO
  END IF


  IF (rup_d(ilev-1)) THEN
  a12 = 0.0d0
  a22 = 0.0d0
    DO n=1,nang
      a12(n,n)=-p1(ilev,n,n)*rup(ilev-1,n,n)
      a22(n,n)=-p2(ilev,n,n)*rup(ilev-1,n,n)
    END DO
  ELSE
    DO n=1,nang
     DO m=1,nang
      a12(n,m)=-p1(ilev,n,n)*rup(ilev-1,n,m)
      a22(n,m)=-p2(ilev,n,n)*rup(ilev-1,n,m)
     END DO
    END DO
  END IF

  END IF ! end of the diagonal case: part 9

   DO n=1,nang
    a11(n,n)=a11(n,n)+1.d0
    a22(n,n)=a22(n,n)+1.d0
   END DO

   DO ivar=1,nvar
    IF(diagonal(ilev,3) .AND. diagonal(ilev,ivar+3)) THEN ! diagonal layer and diagonal perturbation
     IF(rdn_d(ilev).AND.rup_d(ilev-1)) THEN
  da11(:,:,ivar) = 0.0d0
  da12(:,:,ivar) = 0.0d0
  da21(:,:,ivar) = 0.0d0
  da22(:,:,ivar) = 0.0d0
       DO n=1,nang
        da11(n,n,ivar)=-dp2(ilev,n,n,ivar)*rdn(ilev,n,n)
        da12(n,n,ivar)=-dp1(ilev,n,n,ivar)*rup(ilev-1,n,n)
        da21(n,n,ivar)=-dp1(ilev,n,n,ivar)*rdn(ilev,n,n)
        da22(n,n,ivar)=-dp2(ilev,n,n,ivar)*rup(ilev-1,n,n)
       END DO
     ELSE
       DO n=1,nang
        DO m=1,nang
        da11(n,m,ivar)=-dp2(ilev,n,n,ivar)*rdn(ilev,n,m)
        da12(n,m,ivar)=-dp1(ilev,n,n,ivar)*rup(ilev-1,n,m)
        da21(n,m,ivar)=-dp1(ilev,n,n,ivar)*rdn(ilev,n,m)
        da22(n,m,ivar)=-dp2(ilev,n,n,ivar)*rup(ilev-1,n,m)
        END DO
       END DO
     END IF

    ELSE ! non-diagonal layer or perturbation
     IF(rdn_d(ilev).AND.rup_d(ilev-1)) THEN
       DO n=1,nang
        DO m=1,nang
        da11(n,m,ivar)=-dp2(ilev,n,m,ivar)*rdn(ilev,m,m)
        da12(n,m,ivar)=-dp1(ilev,n,m,ivar)*rup(ilev-1,m,m)
        da21(n,m,ivar)=-dp1(ilev,n,m,ivar)*rdn(ilev,m,m)
        da22(n,m,ivar)=-dp2(ilev,n,m,ivar)*rup(ilev-1,m,m)
        END DO
       END DO
     ELSE
      da11(:,:,ivar)=-MATMUL(dp2(ilev,:,:,ivar),rdn(ilev,:,:))
      da12(:,:,ivar)=-MATMUL(dp1(ilev,:,:,ivar),rup(ilev-1,:,:))
      da21(:,:,ivar)=-MATMUL(dp1(ilev,:,:,ivar),rdn(ilev,:,:))
      da22(:,:,ivar)=-MATMUL(dp2(ilev,:,:,ivar),rup(ilev-1,:,:))
     END IF
  
    END IF ! end of non-diagonal/diagonal layer/perturbation
   END DO

    IF(rdn_d(ilev).AND.rup_d(ilev-1)) THEN
     DO n=1,nang
      y1(n)= ust(ilev-1,n) + rup(ilev-1,n,n)*unh(ilev,n)
      y2(n)= vst(ilev,n) + rdn(ilev,n,n)*unh(ilev,n)
       DO ivar=1,nvar
        dy1v(n,ivar)= rup(ilev-1,n,n)*dunh(ilev,n,ivar)
        dy2v(n,ivar)= rdn(ilev,n,n)*dunh(ilev,n,ivar)
       END DO
      END DO
    ELSE
     y1= ust(ilev-1,:) + MATMUL(rup(ilev-1,:,:),unh(ilev,:))
     y2= vst(ilev,:) + MATMUL(rdn(ilev,:,:),unh(ilev,:))
      DO ivar=1,nvar
       dy1v(:,ivar)= MATMUL(rup(ilev-1,:,:),dunh(ilev,:,ivar))
       dy2v(:,ivar)= MATMUL(rdn(ilev,:,:),dunh(ilev,:,ivar))
      END DO
    END IF

  IF(.NOT. diagonal(ilev,3)) THEN ! non-diagonal case: part 10
   b1=MATMUL(p1(ilev,:,:),y1) + MATMUL(p2(ilev,:,:),y2)
   b2=MATMUL(p2(ilev,:,:),y1) + MATMUL(p1(ilev,:,:),y2)
  END IF ! end of the non-diagonal case: part 10
  IF(diagonal(ilev,3)) THEN ! diagonal case: part 10
   DO n=1,nang
    b1(n)=p1(ilev,n,n)*y1(n)+p2(ilev,n,n)*y2(n)
    b2(n)=p2(ilev,n,n)*y1(n)+p1(ilev,n,n)*y2(n)
   END DO
  END IF ! end of the diagonal case: part 10

  IF(.NOT. diagonal(ilev,3)) THEN ! non-diagonal case: part 11
   DO ivar=1,nvar
    db1(:,ivar)=  MATMUL(dp1(ilev,:,:,ivar),y1) &
                + MATMUL(dp2(ilev,:,:,ivar),y2) &
                + MATMUL(p1(ilev,:,:),dy1v(:,ivar)) &
                + MATMUL(p2(ilev,:,:),dy2v(:,ivar))
    db2(:,ivar)=  MATMUL(dp2(ilev,:,:,ivar),y1) &
                + MATMUL(dp1(ilev,:,:,ivar),y2) &
                + MATMUL(p2(ilev,:,:),dy1v(:,ivar)) &
                + MATMUL(p1(ilev,:,:),dy2v(:,ivar))
   END DO

  END IF ! end of the non-diagonal case: part 11
  IF(diagonal(ilev,3)) THEN ! diagonal case: part 11
  DO ivar=1,nvar
    IF(diagonal(ilev,3) .AND. diagonal(ilev,ivar+3)) THEN ! diagonal layer and diagonal perturbation
      DO n=1,nang
       db1(n,ivar)=dp1(ilev,n,n,ivar)*y1(n)+dp2(ilev,n,n,ivar)*y2(n)
       db2(n,ivar)=dp2(ilev,n,n,ivar)*y1(n)+dp1(ilev,n,n,ivar)*y2(n)
      END DO
    ELSE ! non-diagonal layer or perturbation
    db1(:,ivar)=  MATMUL(dp1(ilev,:,:,ivar),y1) + MATMUL(dp2(ilev,:,:,ivar),y2)                 
    db2(:,ivar)=  MATMUL(dp2(ilev,:,:,ivar),y1) + MATMUL(dp1(ilev,:,:,ivar),y2) 
    END IF ! end of non-diagonal/diagonal layer/perturbation
    DO n=1,nang
     db1(n,ivar)=db1(n,ivar)+p1(ilev,n,n)*dy1v(n,ivar)+p2(ilev,n,n)*dy2v(n,ivar)
     db2(n,ivar)=db2(n,ivar)+p2(ilev,n,n)*dy1v(n,ivar)+p1(ilev,n,n)*dy2v(n,ivar)
    END DO
   END DO
  END IF! end of the diagonal case: part 11

  IF(a_d) THEN
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
   CALL dlinrg(nang,a11,a11m1)
   CALL dlinrg(nang,a22,a22m1)
   w0=MATMUL(a21,a11m1)
   w1=MATMUL(a12,a22m1)
   t0=a22-MATMUL(w0,a12)
   t1=a11-MATMUL(w1,a21)
   y0=b2 -MATMUL(w0,b1)
   y1=b1 -MATMUL(w1,b2)
  END IF

   DO ivar=1,nvar

    dt0(:,:,ivar)=0.d0
    dt1(:,:,ivar)=0.d0

    IF(a_d.AND.diagonal(ilev,ivar+3)) THEN ! all diag

     DO n=1,nang
  s1=da21(n,n,ivar)*a11m1(n,n)-a21(n,n)*da11(n,n,ivar)*a11m1(n,n)**2 ! =dw0
  s2=da12(n,n,ivar)*a22m1(n,n)-a12(n,n)*da22(n,n,ivar)*a22m1(n,n)**2 ! =dw1

      dy0v(n,ivar)=db2(n,ivar)-s1*b1(n)-w0(n,n)*db1(n,ivar)
      dy1v(n,ivar)=db1(n,ivar)-s2*b2(n)-w1(n,n)*db2(n,ivar)
     
      dt0(n,n,ivar)=da22(n,n,ivar)-s1*a12(n,n)-w0(n,n)*da12(n,n,ivar)
      dt1(n,n,ivar)=da11(n,n,ivar)-s2*a21(n,n)-w1(n,n)*da21(n,n,ivar)
     END DO

    ELSE ! a_d .OR. diagonal(ilev,ivar+3)  - non diag

     da11m1(:,:,ivar)=-MATMUL(MATMUL(a11m1,da11(:,:,ivar)),a11m1)
     da22m1(:,:,ivar)=-MATMUL(MATMUL(a22m1,da22(:,:,ivar)),a22m1)

     dw0(:,:,ivar) = MATMUL(da21(:,:,ivar),a11m1) + MATMUL(a21,da11m1(:,:,ivar))
     dw1(:,:,ivar) = MATMUL(da12(:,:,ivar),a22m1) + MATMUL(a12,da22m1(:,:,ivar))

     dy0v(:,ivar)=db2(:,ivar) - MATMUL(dw0(:,:,ivar),b1) - MATMUL(w0,db1(:,ivar))
     dy1v(:,ivar)=db1(:,ivar) - MATMUL(dw1(:,:,ivar),b2) - MATMUL(w1,db2(:,ivar))
     dt0(:,:,ivar)=da22(:,:,ivar) - MATMUL(dw0(:,:,ivar),a12(:,:)) - MATMUL(w0,da12(:,:,ivar))
     dt1(:,:,ivar)=da11(:,:,ivar) - MATMUL(dw1(:,:,ivar),a21(:,:)) - MATMUL(w1,da21(:,:,ivar))

    END IF ! non-diagonal/diagonal a_d and perturbation

   END DO

  IF(a_d) THEN
   w0=0.d0
   w1=0.d0
   DO n=1,nang
    w0(n,n) = 1.0d0 /(a22(n,n)-a21(n,n)*a12(n,n)/a11(n,n))
    w1(n,n) = 1.0d0 /(a11(n,n)-a21(n,n)*a12(n,n)/a22(n,n))
   END DO
  ELSE
   CALL dlinrg(nang,t0,w0) ! w0=(a22-a21*a11**(-1)*a12)**(-1)
   CALL dlinrg(nang,t1,w1) ! w1=(a11-a12*a22**(-1)*a21)**(-1)
  END IF

! This mode of calculation of brightness temperature can be used
! if variation of the properties of the layer takes place and
! matrices _w0,_w1 and vectors _y0,_y1 are known anyway. Otherwise,
! the calculation should proceed according to the second _ilev loop.

 IF(ilev==nlev) v(nlev,:)=f(nlev+1,:)

 IF(.NOT. a_d) THEN ! non-diagonal _a => non-diagonal _w0,_w1
      u(ilev,:)=unh(ilev,:)+MATMUL(w1,y1)
      v(ilev-1,:)=unh(ilev,:)+MATMUL(w0,y0)
 ELSE ! diagonal _a => diagonal _w0,_w1
     DO n=1,nang
        u(ilev,n)=unh(ilev,n)+w1(n,n)*y1(n)
        v(ilev-1,n)=unh(ilev,n)+w0(n,n)*y0(n)
     END DO
 END IF

 IF(ilev==1) THEN
     IF(.NOT. rup_d(0)) THEN ! non-diagonal  _rup_d(0)
         u(0,:)=ust(0,:)+MATMUL(rup(0,:,:),v(0,:))
     ELSE ! diagonal  _rup_d(0)
        DO n=1,nang
         u(0,n)=ust(0,n)+rup(0,n,n)*v(0,n)
        END DO
     END IF
  END IF

  !if(ilev == 7) print*, 'ok dz3',ilev
DO ivar=1,nvar

  IF(.NOT. a_d) THEN ! non-diagonal _w0,_w1

   dz2(ilev-1,:,ivar)=dunh(ilev,:,ivar) &
            +MATMUL(w0,dy0v(:,ivar)-MATMUL(dt0(:,:,ivar),MATMUL(w0,y0)))

   dz3(ilev,:,ivar)=dunh(ilev,:,ivar) &
            +MATMUL(w1,dy1v(:,ivar)-MATMUL(dt1(:,:,ivar),MATMUL(w1,y1)))

  ELSE ! diagonal _w0,_w1

   IF (.NOT. diagonal(ilev,ivar+3)) THEN ! non-diagonal perturbation

     DO n=1,nang
      s1=dy0v(n,ivar)
      s2=dy1v(n,ivar)
      DO m=1,nang
       s1=s1-dt0(n,m,ivar)*w0(m,m)*y0(m)
       s2=s2-dt1(n,m,ivar)*w1(m,m)*y1(m)
      END DO
      dz2(ilev-1,n,ivar)=dunh(ilev,n,ivar)+w0(n,n)*s1
      dz3(ilev,n,ivar)=dunh(ilev,n,ivar)+w1(n,n)*s2
     END DO

   ELSE ! diagonal perturbation

     DO n=1,nang
      dz2(ilev-1,n,ivar)= dunh(ilev,n,ivar) &
                    +w0(n,n)*(dy0v(n,ivar)-dt0(n,n,ivar)*w0(n,n)*y0(n))
      dz3(ilev,n,ivar)= dunh(ilev,n,ivar) &
                    +w1(n,n)*(dy1v(n,ivar)-dt1(n,n,ivar)*w1(n,n)*y1(n))
     END DO

   END IF ! non-diagonal/diagonal perturbation
  END IF ! non-diagonal/diagonal   _w0,_w1
END DO ! ivar
  END DO ! end of the fifth _ilev loop

  DO ilev=1,oblv ! lower part
    IF(.NOT. s_d(ilev)) THEN ! non-diagonal matrices s11,..,s22
      IF(.NOT. rdn_d(oblv)) THEN ! non-diagonal matrix rdn(oblv)
        t0=-MATMUL(s12(ilev,:,:),rdn(oblv,:,:))
        DO n=1,nang
          t0(n,n)=t0(n,n)+1.d0
        END DO
        CALL dlinrg(nang,t0,t1) ! t1=(1-S12*R)**(-1)
        DO ivar=1,nvar
          dtb_up_dp(ilev,:,ivar)=MATMUL(t1,MATMUL(s11(ilev,:,:),dz3(ilev,:,ivar)))
          dtb_dn_dp(ilev,:,ivar)=MATMUL(rdn(oblv,:,:),dtb_up_dp(ilev,:,ivar))
        END DO ! ivar
      ELSE ! diagonal matrix rdn(oblv)
        DO n=1,nang
          DO m=1,nang
            t0(n,m)=-s12(ilev,n,m)*rdn(oblv,m,m)
          END DO
          t0(n,n)=t0(n,n)+1.d0
        END DO
        CALL dlinrg(nang,t0,t1) ! t1=(1-S12*R)**(-1)
        !if(ilev == 7) print*, 'ok dtb_up_dp',ilev
        DO ivar=1,nvar
          dtb_up_dp(ilev,:,ivar)=MATMUL(t1,MATMUL(s11(ilev,:,:),dz3(ilev,:,ivar)))
          DO n=1,nang
            dtb_dn_dp(ilev,n,ivar)=rdn(oblv,n,n)*dtb_up_dp(ilev,n,ivar)
          END DO
        END DO ! ivar
      END IF ! non-diagonal/diagonal matrix rdn(oblv)

    ELSE ! diagonal matrices s11,..,s22
      IF(.NOT. rdn_d(oblv)) THEN ! non-diagonal matrix rdn_d(oblv)
        DO n=1,nang
          DO m=1,nang
            t0(n,m)=-s12(ilev,n,n)*rdn(oblv,n,m)
          END DO
          t0(n,n)=t0(n,n)+1.d0
        END DO

        CALL dlinrg(nang,t0,t1) ! t1=(1-S12*R)**(-1)
        DO ivar=1,nvar
          DO n=1,nang
            s1=0.d0
            DO m=1,nang
              s1=s1+t1(n,m)*s11(ilev,m,m)*dz3(ilev,m,ivar)
            END DO
            dtb_up_dp(ilev,n,ivar)=s1
          END DO
          dtb_dn_dp(ilev,:,ivar)=MATMUL(rdn(oblv,:,:),dtb_up_dp(ilev,:,ivar))   
        END DO
      ELSE ! diagonal matrix rdn_d(oblv)
        DO ivar=1,nvar
          DO n=1,nang
            dtb_up_dp(ilev,n,ivar)=s11(ilev,n,n)/(1-s12(ilev,n,n)*rdn(oblv,n,n))*dz3(ilev,n,ivar)
            dtb_dn_dp(ilev,n,ivar)=rdn(oblv,n,n)*dtb_up_dp(ilev,n,ivar)
          END DO  
        END DO
      END IF ! non-diagonal/diagonal matrix rdn_d(oblv)
    END IF ! non-diagonal/diagonal matrices s11,..,s22
  END DO ! end of the ilev loop

!--------------------------------------------------------------------------

  DO ilev=oblv+1,nlev ! upper part

   IF(.NOT. s_d(ilev-1)) THEN ! non-diagonal matrices s11,..,s22

    IF(.NOT. rup_d(oblv)) THEN ! non-diagonal matrix rup(oblv)

     w0=-MATMUL(s21(ilev-1,:,:),rup(oblv,:,:))
     DO n=1,nang
      w0(n,n)=w0(n,n)+1.d0
     END DO

     CALL dlinrg(nang,w0,w1) ! w1=(1-S21*R)**(-1)

     DO ivar=1,nvar
       dtb_dn_dp(ilev,:,ivar)=MATMUL(w1,MATMUL(s22(ilev-1,:,:),dz2(ilev-1,:,ivar)))
       dtb_up_dp(ilev,:,ivar)=MATMUL(rup(oblv,:,:),dtb_dn_dp(ilev,:,ivar))
     END DO ! ivar

    ELSE ! diagonal matrix rup(oblv)

     DO n=1,nang
      DO m=1,nang
       w0(n,m)=-s21(ilev-1,n,m)*rup(oblv,m,m)
      END DO
      w0(n,n)=w0(n,n)+1.d0
     END DO

     CALL dlinrg(nang,w0,w1) ! w1=(1-S21*R)**(-1)

     DO ivar=1,nvar
      dtb_dn_dp(ilev,:,ivar)=MATMUL(w1,MATMUL(s22(ilev-1,:,:),dz2(ilev-1,:,ivar)))
      DO n=1,nang
        dtb_up_dp(ilev,n,ivar)=rup(oblv,n,n)*dtb_dn_dp(ilev,n,ivar)
      END DO
     END DO ! ivar

    END IF ! non-diagonal/diagonal matrix rup(oblv)

   ELSE ! diagonal matrices s11,..,s22
    IF(.NOT. rup_d(oblv)) THEN ! non-diagonal matrix rup_d(oblv)
    DO n=1,nang
     DO m=1,nang
      w0(n,m)=-s21(ilev-1,n,n)*rup(oblv,n,m)
     END DO
     w0(n,n)=w0(n,n)+1.d0
    END DO
    CALL dlinrg(nang,w0,w1) ! w1=(1-S21*R)**(-1)
    DO ivar=1,nvar
     DO n=1,nang
      s1=0.d0
      DO m=1,nang
       s1=s1+w1(n,m)*s22(ilev-1,m,m)*dz2(ilev-1,m,ivar)
      END DO
      dtb_dn_dp(ilev,n,ivar)=s1
     END DO
     dtb_up_dp(ilev,:,ivar)=MATMUL(rup(oblv,:,:),dtb_dn_dp(ilev,:,ivar))   
    END DO
    ELSE ! diagonal matrix rup(oblv)
     DO ivar=1,nvar
      DO n=1,nang
       dtb_dn_dp(ilev,n,ivar)= s22(ilev-1,n,n)/(1-s21(ilev-1,n,n)*rup(oblv,n,n))*dz2(ilev-1,n,ivar)
       dtb_up_dp(ilev,n,ivar)=rup(oblv,n,n)*dtb_dn_dp(ilev,n,ivar)
      END DO  
     END DO
    END IF ! non-diagonal/diagonal matrix rup_d(oblv)
   END IF ! non-diagonal/diagonal matrices s11,..,s22
  END DO ! end of the upper part

! deallocate variables
deallocate(a,stat=alloc_err)
deallocate(b,stat=alloc_err)
deallocate(am1,stat=alloc_err)
deallocate(p1,stat=alloc_err)
deallocate(p2,stat=alloc_err)
deallocate(q,stat=alloc_err)
deallocate(da,stat=alloc_err)
deallocate(db,stat=alloc_err)
deallocate(dp1,stat=alloc_err)
deallocate(dp2,stat=alloc_err)

END SUBROUTINE core95

!=============================================================
SUBROUTINE md_inv32 (nang,h,lama,ma,lamab,mab,lama12d,p1,p2, &
                        dlama,dma,dlamab,dmab,dp1,dp2,nvar,L)
!=============================================================
! This subroutine inverts a matrix of a special form:
! t(n,m)= SUM [i,j=1 to i,j=nang] ma(n,i)*w(i,j)*ma(m,j) !kz,eqn. (51)
! where
! w(n,m)=SUM [k=1 to k=nang] a(n,k)*b(m,k)/dzeta(k)-2*nang*delta(n,m)
! and
! dzeta(k)= SQRT(lamab(k))/sinh(SQRT(lamab(k)*h)) !kz,eqn. (54)
! a(n,k)=[SQRT(lama(n))+SQRT(lamab(k))/tanh(SQRT(lamab(k)*h/2)) & !kz,eqn. (52)
!                                    /SQRT(lama(n))] * mab(n,k)
!
! b(m,k)=[SQRT(lama(m))+SQRT(lamab(k))*tanh(SQRT(lamab(k)*h/2)) & !kz,eqn. (53)
!                                    /SQRT(lama(m))] * mab(m,k)
! and both _ma(n,k) and _mab(n,k) are orthogonal matrices.
! It is assumed h > 0.
! Parameter _nvar correspond to a vector of variations:
! the "derivatives" vectors have in fact one extra dimension.
! Logical array L characterizes "diagonal properties" and
! corresponds to the _ilev-th row of matrix diagonal


IMPLICIT NONE

INTEGER,INTENT(IN) :: nang,nvar

DOUBLE PRECISION,DIMENSION(nang,nang),INTENT(IN) ::  ma,mab 
DOUBLE PRECISION,DIMENSION(nang,nang,nvar),INTENT(IN) :: dma,dmab
DOUBLE PRECISION,DIMENSION(nang),INTENT(IN) ::   lama,lamab,lama12d
DOUBLE PRECISION,DIMENSION(nang,nvar),INTENT(IN) ::  dlama,dlamab 

DOUBLE PRECISION,INTENT(IN) :: h
LOGICAL,DIMENSION(nvar+3),INTENT(IN) :: L

DOUBLE PRECISION,DIMENSION(nang,nang),INTENT(OUT) ::  p1,p2
DOUBLE PRECISION,DIMENSION(nang,nang,nvar),INTENT(OUT) :: dp1,dp2

INTEGER :: k,n,m,ivar
DOUBLE PRECISION :: g,x,em,t,lama12,u1,u2,u3,v1,v2,v3
DOUBLE PRECISION,DIMENSION(nvar) :: dg,dx,dem,dt

DOUBLE PRECISION,DIMENSION(nang) ::  dzeta,ad,b1d,b2d,c1d,c2d
DOUBLE PRECISION,DIMENSION(nang,nvar) :: ddzeta
DOUBLE PRECISION,DIMENSION(nang,nang) :: a,b1,b2,am1,bm1,c1,c2,ac1,ac2
DOUBLE PRECISION,DIMENSION(nang,nang,nvar) :: da,db1,db2,dam1,dbm1,dc1,dc2
DOUBLE PRECISION :: s1,s2

 DO  k=1,nang ! loop over k

  g=DSQRT(lamab(k))
  x= g*h
  DO ivar=1,nvar
   dg(ivar)=dlamab(k,ivar)/2/g
   dx(ivar)=dg(ivar)*h
  END DO

  IF(x < 1.d-3) THEN
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

    IF(.NOT. L(3)) THEN ! non-diagonal layer
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
       da (n,k,ivar)=u1*dmab(n,k,ivar) &
               +((1-t/lama(n))*dlama(n,ivar)/2+dt(ivar))/lama12*mab(n,k)
       db1(n,k,ivar)=u2*dmab(n,k,ivar)+((t-v1)*dlama(n,ivar)/2+v2)*v3
       db2(n,k,ivar)=u3*dmab(n,k,ivar)+((t+v1)*dlama(n,ivar)/2-v2)*v3
      END DO ! ivar
     END DO ! (n)
   ELSE ! diagonal layer

     m=0
     DO ivar=1,nvar
      IF(.NOT. L(ivar+3)) THEN
       m=1
       EXIT
      END IF
     END DO ! if there are non-diagonal perturbations,m==1

! do this block only if there are non-diagonal perturbations (m > 0)

  IF(m > 0) THEN
     DO n=1,nang
      lama12=lama12d(n)
      x =1/t/lama12
      u1=lama12+t/lama12
      u2=lama12+lamab(k)*x
      u3=lama12-lamab(k)*x
       DO ivar=1,nvar
        IF(.NOT. L(ivar+3)) THEN ! non-diagonal perturbation
         da (n,k,ivar)=u1*dmab(n,k,ivar)                                  
         db1(n,k,ivar)=u2*dmab(n,k,ivar)
         db2(n,k,ivar)=u3*dmab(n,k,ivar)
        END IF
       END DO ! ivar
     END DO ! n
  END IF ! block for non-diagonal perturbations  (m > 0)

     lama12=lama12d(k)
         x =1/t/lama12
     ad (k)=lama12+t/lama12
     b1d(k)=lama12+lamab(k)*x
     b2d(k)=lama12-lamab(k)*x

     v1=lamab(k)/lama(k)

       DO ivar=1,nvar
        v2=dlamab(k,ivar)-lamab(k)*dt(ivar)/t
        IF(.NOT. L(ivar+3)) THEN ! non-diagonal perturbation
     da (k,k,ivar)=da (k,k,ivar)+((1-t/lama(k))*dlama(k,ivar)/2+dt(ivar))/lama12
     db1(k,k,ivar)=db1(k,k,ivar)+((t-v1)*dlama(k,ivar)/2+v2)*x
     db2(k,k,ivar)=db2(k,k,ivar)+((t+v1)*dlama(k,ivar)/2-v2)*x
        ELSE ! diagonal perturbation
     da (k,k,ivar)=((1-t/lama(k))*dlama(k,ivar)/2+dt(ivar))/lama12
     db1(k,k,ivar)=((t-v1)*dlama(k,ivar)/2+v2)*x
     db2(k,k,ivar)=((t+v1)*dlama(k,ivar)/2-v2)*x
        END IF ! end of diagonal/non-diagonal perturbation
      END DO ! ivar
   END IF ! end of diagonal/non-diagonal layer

 END DO ! (k)

 IF (L(3)) THEN ! diagonal layer
  p1=0.d0 ! remove after debugging
  p2=0.d0 ! remove after debugging
  DO n=1,nang
   c1d(n)  = dzeta(n)/b1d(n)/ad(n)
   c2d(n)  =  -b2d(n)/b1d(n)
   p1 (n,n)= 2*c1d(n)
   p2 (n,n)=   c2d(n)+p1(n,n)
  END DO
  DO ivar=1,nvar
    IF(.NOT. L(ivar+3)) THEN ! non-diagonal perturbation

    DO n=1,nang
     DO m=1,nang
      dam1(n,m,ivar)=- da(n,m,ivar)/ ad(n)/ ad(m)
      dbm1(n,m,ivar)=-db1(n,m,ivar)/b1d(n)/b1d(m)
     END DO
    END DO

    DO n=1,nang
     DO m=n,nang
      dp1(n,m,ivar)= dma(n,m,ivar)*c1d(m)+c1d(n)*dma(m,n,ivar) &
                    +dbm1(m,n,ivar)*dzeta(m)/ad (m) &
                    +dam1(n,m,ivar)*dzeta(n)/b1d(n)
      dp2(n,m,ivar)= dma (n,m,ivar)*c2d(m)+c2d(n)*dma(m,n,ivar) &
                    -dbm1(m,n,ivar)*b2d(m)-db2(m,n,ivar)/b1d(n)

      dp1(m,n,ivar)=dp1(n,m,ivar)
      dp2(m,n,ivar)=dp2(n,m,ivar)

     END DO
      dp1(n,n,ivar)=dp1(n,n,ivar)+ddzeta(n,ivar)/ad(n)/b1d(n)
    END DO

    ELSE ! diagonal perturbation
    DO n=1,nang
      s1=- da(n,n,ivar)/ ad(n)/ ad(n)
      s2=-db1(n,n,ivar)/b1d(n)/b1d(n)
      dp1(n,n,ivar)= ddzeta(n,ivar)/b1d(n)/ad(n) &
                     +dzeta(n)* (s1/b1d(n)+s2/ad(n)) 
      dp2(n,n,ivar)=-s2*b2d(n)-db2(n,n,ivar)/b1d(n)
    END DO
    END IF ! end of diagonal/non-diagonal perturbation
  END DO ! (ivar)

  dp1=  2*dp1
  dp2=dp2+dp1

  RETURN
 END IF ! end of diagonal layer

  CALL dlinrg(nang,a,am1)
  CALL dlinrg(nang,b1,bm1)


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
   s1=s1+(dbm1(k,n,ivar)*dzeta(k)+bm1(k,n)*ddzeta(k,ivar))*am1(k,m) &
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

   DO ivar=1,nvar

    ac1=MATMUL(dma(:,:,ivar),c1)
    ac2=MATMUL(dma(:,:,ivar),c2)
    am1=MATMUL(ma,dc1(:,:,ivar))
    bm1=MATMUL(ma,dc2(:,:,ivar))

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
   END DO

   p2= p2+ p1
  dp2=dp2+dp1

END SUBROUTINE  md_inv32

!=============================================================
subroutine dlinrg(n,a,ainv)
!=============================================================
! Organizes input to gaussj,calls gaussj
!
! History:
! 9/26/2020 Kevin Schaefer deleted unused variables and arguments
implicit none

integer*4 n,np,m,mp,j
real*8,dimension(n,n) :: a,ainv,aa
real*8,dimension(n)   :: b
    np = n
    m = 1
    mp = 1
    aa = a
    do j = 1,n
        b(j) = 0.0d0
    end do
    call gaussj(aa,n,np,b,m,mp)
    ainv = aa
end subroutine dlinrg

