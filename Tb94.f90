!======================================================
SUBROUTINE  tb94()
!======================================================
! This subroutine prepares data for CORE95.F90 which numerically 
! solves RT equations represented in discrete form
!
! Scattering matrices _phaseff(k,i,j) and _phasefb(k,i,j) are supposed
! to be  symmetric with respect to indeces (i,j) (only upper half of
! those matrices is used for calculations). No specific normalization
! of the phase matrices is required.
!
! Within layers all environmental parameters are assumed to be constant. 
! "d"-values representing derivatives wrt appropriate parameters
!
! issues:
!  12/12/20 (Kevin Schaefer) dt_srf, dt_cb, and dsurf_reflec always zero
!
! History:
! 4/24/2003  blw eps_err introduced for max allowed relative error in brightness temperature.
! 11/11/2003 blw CORE95.F90 called rather then DTB27.F90 because it has special treatment for diagonal layers.
! 9/11/2004  blw eps_err deleted due to parrallelization errors
! 5/20/2005  blw fixed diagonal layer problem
! 9/26/2020  Kevin Schaefer deleted unused variables
! 9/26/2020  Kevin Schaefer deleted unused arguments for core95 call
! 12/9/2020  Kevin Schaefer deleted unused debugging code
! 12/12/2020 Kevin Schaefer deleted commented code, outdated equations
! 12/12/2020 Kevin Schaefer moved comment variable definitions to dotlrt_variables
!------------------------------------------------------
use dotlrt_variables
IMPLICIT NONE

real(8), DIMENSION(nlev) :: h
real(8), allocatable :: a0(:,:,:)
real(8), allocatable :: b0(:,:,:)
real(8), allocatable :: da0(:,:,:,:)
real(8), allocatable :: db0(:,:,:,:)
real(8) f(0:nlev+1,nangover2)
real(8) df(0:nlev+1,nangover2,nvar)
real(8) r(nangover2,nangover2)
real(8) dr(nangover2,nangover2)
real(8) u(0:nlev,nangover2)
real(8) v(0:nlev,nangover2)
real(8) du(0:nlev,nangover2,nvar)
real(8) dv(0:nlev,nangover2,nvar)
real(8) t1(nangover2)
real(8) t2(nangover2)
real(8) pt(nangover2)
real(8) pt_sc(nangover2)
real(8) dpt_g(nangover2,max_nphase)
real(8) dpt_sc(nangover2,max_nphase)
real(8) at(nangover2,nangover2)
real(8) bt(nangover2,nangover2)
real(8) at_sc(nangover2,nangover2)
real(8) bt_sc(nangover2,nangover2)
real(8) dat_g(nangover2,nangover2,max_nphase)
real(8) dbt_g(nangover2,nangover2,max_nphase)
real(8) dat_sc(nangover2,nangover2,max_nphase)
real(8) dbt_sc(nangover2,nangover2,max_nphase)
LOGICAL LP(0:nlev,nvar+3) ! (-) logical matrix to determine if layer will be treated as diagonal or not
integer i
integer j
integer ivar
integer iphase ! (-) hydrometeor phase index
integer ilev ! (-) layer or level index
integer jlev ! (-) secondary layer index
real(8) dt_srf
real(8) dt_cb
real(8) albedo
real(8) dsurf_reflec(nangover2)

allocate( a0(nlev,nangover2,nangover2), b0(nlev,nangover2,nangover2) )
allocate( da0(nlev,nangover2,nangover2,nvar), db0(nlev,nangover2,nangover2,nvar) )


         dt_srf = 0.d0
          dt_cb = 0.d0
   dsurf_reflec = 0.d0
   dscat_cloud1 = 1.0d0
       dal_gas1 = 1.0d0
    dabs_cloud1 = 1.0d0
  dtemperature1 = 1.0d0

  DO i=1,nangover2
    t1(i) = DSQRT(quad_wts(i)/cos_ang(i))
    t2(i) = DSQRT(quad_wts(i)*cos_ang(i))
  END DO
  
  DO ilev=1,nlev
    DO i=1,nangover2
      DO j=i,nangover2
        at(i,j) = t1(i) * phaseff1(ilev,i,j) * t1(j)
        bt(i,j) = t1(i) * phasefb1(ilev,i,j) * t1(j)
        at(j,i) = at(i,j)
        bt(j,i) = bt(i,j)
        at_sc(i,j) = t1(i) * phaseff1_sc(ilev,i,j) * t1(j)
        bt_sc(i,j) = t1(i) * phasefb1_sc(ilev,i,j) * t1(j)
        at_sc(j,i) = at_sc(i,j)
        bt_sc(j,i) = bt_sc(i,j)
        do iphase = 1, nphase
          dat_g(i,j,iphase)  = t1(i) * dphaseff1_g(ilev,i,j,iphase)  * t1(j)
          dbt_g(i,j,iphase)  = t1(i) * dphasefb1_g(ilev,i,j,iphase)  * t1(j)
          dat_sc(i,j,iphase) = t1(i) * dphaseff1_sc(ilev,i,j,iphase) * t1(j)
          dbt_sc(i,j,iphase) = t1(i) * dphasefb1_sc(ilev,i,j,iphase) * t1(j)
          dat_g(j,i,iphase)  = dat_g(i,j,iphase)
          dbt_g(j,i,iphase)  = dbt_g(i,j,iphase)
          dat_sc(j,i,iphase) = dat_sc(i,j,iphase)
          dbt_sc(j,i,iphase) = dbt_sc(i,j,iphase)
        end do ! iphase
      END DO
    END DO

    DO i = 1, nangover2
      pt(i)=0.d0
      pt_sc(i) = 0.0d0
      do iphase = 1, nphase
        dpt_g(i,iphase)  = 0.0d0
        dpt_sc(i,iphase) = 0.0d0
      end do ! iphase
      DO j = 1, nangover2
        pt(i) =  pt(i) + ( at(j,i)+ bt(j,i))*t2(j)
        pt_sc(i) =  pt_sc(i) + ( at_sc(j,i)+ bt_sc(j,i))*t2(j)
        do iphase = 1, nphase
          dpt_g(i,iphase) = dpt_g(i,iphase) + ( dat_g(j,i,iphase) + dbt_g(j,i,iphase) ) * t2(j)
          dpt_sc(i,iphase) = dpt_sc(i,iphase) + ( dat_sc(j,i,iphase)+ dbt_sc(j,i,iphase) ) * t2(j)
        end do ! iphase
      END DO
   END DO
   
   h(ilev)=altitude(ilev+1)-altitude(ilev)
   
    DO i=1,nangover2
      DO j=1,nangover2
        a0(ilev,i,j) = - at_sc(i,j)
        b0(ilev,i,j) = - bt_sc(i,j)
        da0(ilev,i,j,1) = 0.d0
        db0(ilev,i,j,1) = 0.d0
        da0(ilev,i,j,2) = 0.d0
        db0(ilev,i,j,2) = 0.d0
        do iphase = 1, nphase
          ivar = 1 + 2 * iphase
          da0(ilev,i,j,ivar) = - dat_sc(i,j,iphase)
          db0(ilev,i,j,ivar) = - dbt_sc(i,j,iphase)
          ivar = 2 + 2 * iphase
          da0(ilev,i,j,ivar) = - dat_g(i,j,iphase)
          db0(ilev,i,j,ivar) = - dbt_g(i,j,iphase)
        end do ! iphase
        IF( i == j) THEN
          a0(ilev,i,i) = a0(ilev,i,i) + abs_total1(ilev) / cos_ang(i) + pt_sc(i) / t2(i)
          da0(ilev,i,i,2) = da0(ilev,i,i,2) + ( dal_gas1(ilev) + dabs_cloud1(ilev) ) / cos_ang(i)
          do iphase = 1, nphase
            ivar = 1 + 2 * iphase
            da0(ilev,i,i,ivar) = da0(ilev,i,i,ivar)+ dscat_cloud1(ilev) * pt(i) / t2(i) 
            ! this is why we need phase11 without scat_cloud
            ivar = 2 + 2 * iphase
            da0(ilev,i,i,ivar) = da0(ilev,i,i,ivar) + dpt_g(i,iphase) / t2(i)
          end do ! iphase
        END IF    ! i == j

      END DO   ! loop over j

      f(ilev,i)= abs_total1(ilev) * temperature1(ilev) * t1(i)
      
      df(ilev,i,1) = abs_total1(ilev) * dtemperature1(ilev) * t1(i)
      df(ilev,i,2) = ( dal_gas1(ilev) + dabs_cloud1(ilev) ) *  temperature1(ilev) * t1(i)
      do iphase = 1, nphase
        ivar = 1 + 2 * iphase
        df(ilev,i,ivar) = 0.d0
        ivar = 2 + 2 * iphase
        df(ilev,i,ivar) = 0.d0
      end do ! iphase
    END DO   ! loop over i
  END DO  ! loop over layers
  
  r=0.d0
  dr=0.d0
  DO i=1,nangover2
     r(i,i) =  surf_reflec(i)
    dr(i,i) = dsurf_reflec(i)
  END DO
  DO i=1,nangover2
    f(0     ,i) = (1- r(i,i)) * surf%temp * t2(i)
    f(nlev+1,i) =t_cosmic  * t2(i)
    DO ivar=1,nvar
      df(0     ,i,ivar) = -dr(i,i) * surf%temp * t2(i)+ (1- r(i,i)) * dt_srf * t2(i)
      df(nlev+1,i,ivar) = dt_cb * t2(i)
    END DO
  END DO
 
! defining diagonal layers (04/24/03) 
  LP = .false. ! all layers/perturbations are considered non-diagonal
  DO ilev=1,nlev
    albedo = scat_cloud1(ilev) / ( scat_cloud1(ilev) + abs_total1(ilev) )
    IF( albedo .le. d_albedo ) THEN
      LP(ilev,3)=.true.   ! diagonal layer
      do ivar = 1, nvar
        LP(ilev,3+ivar)=.true.
      end do
    END IF
  END DO

  CALL core95 (nlev,nangover2,h,a0,b0,f,r,u,v,da0,db0,df,du,dv,obs_lev,nvar,LP)
  
  DO ilev = 0,nlev
    DO i=1,nangover2
      tb_pl(ilev,i)     = u(ilev,i)/t2(i)
      tb_mn(ilev,i)     = v(ilev,i)/t2(i)
      DO ivar = 1, nvar
        dtb_pl(ilev,i,ivar) = du(ilev,i,ivar) / t2(i)
        dtb_mn(ilev,i,ivar) = dv(ilev,i,ivar) / t2(i)
        ! correction for absorption derivative
        if( ivar == 2 ) then
          jlev = ilev
          if( ilev == 0 ) jlev = 1
          dtb_pl(ilev,i,ivar) = dtb_pl(ilev,i,ivar) * 0.25d0 * da0(jlev,i,i,ivar)
          dtb_mn(ilev,i,ivar) = dtb_mn(ilev,i,ivar) * 0.25d0 * db0(jlev,i,i,ivar)
        end if
      END DO
    END DO
  END DO

deallocate(a0)
deallocate(b0)
deallocate(da0)
deallocate(db0)

END SUBROUTINE tb94
