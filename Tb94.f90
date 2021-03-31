
SUBROUTINE  tb94()

!------------------------------------------------------
!                         Modification 04/24/03
!
!   Parameter _eps_err is introduced which determines max allowed 
!   relative error in brightness temperature. It is used to calculate
!   logical matrix LP which determines if the layer will be treated as
!   diagonal or not.
!
!   Subroutine CORE95.F90 is called here rather then DTB27.F90. The 
!   difference between two is that CORE95.F90 has special treatment for
!   diagonal layers.
!------------------------------------------------------
!
!   This subroutine prepares data for CORE95.F90 which numerically 
!   solves RT equations represented in the following discrete form:
!
!   c(i)* d/dz Tb_pl(k,i) =
!
!  -[al_gas(k) + abs_cloud(k) + scat_cloud(k)*SUM {from j=1 to j=nangover2}
!
!                           gamma(j)*(phaseff(k,j,i)+phasefb(k,j,i))]*
!  * tb_pl(k,i) +
!
!  +  scat_cloud(k)* 
!
!   * SUM {from j=1 to j=nangover2} gamma(j)*phaseff(k,i,j)*tb_pl(k,j) +
!
!  +  scat_cloud(k)* 
!
!   * SUM {from j=1 to j=nangover2} gamma(j)*phasefb(k,i,j)*tb_mn(k,j)
!
!  + ( al_gas(k,i) + abs_cloud(k))*temperature(k) ;
! 
!
!  -c(i)* d/dz Tb_mn(k,i) =
!
!  -[al_gas(k) + abs_cloud(k) + scat_cloud(k)*SUM {from j=1 to j=nangover2}
!
!                           gamma(j)*(phaseff(k,j,i)+phasefb(k,j,i))]*
!  * tb_mn(k,i) +
!
!  +  scat_cloud(k)* 
!
!   * SUM {from j=1 to j=nangover2} gamma(j)*phasefb(k,i,j)*tb_pl(k,j) +
!
!  +  scat_cloud(k)* 
!
!   * SUM {from j=1 to j=nangover2} gamma(j)*phaseff(k,i,j)*tb_mn(k,j)
!
!  + ( al_gas(k,i) + abs_cloud(k))*temperature(k) ;
!
!   Here:
!
!   i= 1,..,nangover2; k=1,2,...,nlr, and cs(i)=COS(teta(i))
!
!
!   Boundary conditions: 
!
!   k=0    : tb_pl(i)=surf_reflec(i)*tb_mn(i)+(1-surf_reflec(i))*t_surf
!
!   k=nlr+1: tb_mn(i)=t_cb
!
!   i= 1,..,nangover2
!
!   Here index _k corresponds to the number of the layer, index _i
!   corresponds to the number of the angle _teta(i), and
!
!   gamma(j)=cris_quad_angle(j) - weights associated with numerical
!   integratrion.
!
!   Scattering matrices _phaseff(k,i,j) and _phasefb(k,i,j) are supposed
!   to be  symmetric with respect to indeces (i,j) (only upper half of
!   those matrices is used for calculations). No specific normalization
!   of the phase matrices is required.
!
!   Within layers all environmental parameters (al_gas,.., phaseff,..
!   , etc. are assumed to be constant). However, the brightness
!   temperatures _tb_pl and _tb_mn , of course, vary within layers.
!                        
!
!                                INPUT:
!
!    nlr        ( I) -   number of layers. The first layer: k=1 is
!                        located at:
!
!                              _altitude(1  ) < z < altitude(2   ),
!
!                        and the last layer at:
!
!                              _altitude(nlr) < z < altitude(nlr+1)
!
!    nangover2       ( I) -   number of angles
!
!    nvar       ( I) -   number of variable parameters
!
!            cs(nangover2)      (DP) -  cosines of the incidence angles.
!
!  cris_quad_wghts(nangover2)   (DP) - weights, which could be used for 
!                                 numerical intagration
!
!
!  altitude(nlr+1)  (DP) - elevations (see above)
!
!  temperature(nlr) (DP) - thermodynamic temperature of the layer
!
!       al_gas(nlr) (DP) - gas absorbtion (1/length) at a layer
!
!   scat_cloud(nlr) (DP) - coefficient of scattering by clouds
!                          (1/length) at a layer. This parameter is
!                          equivalent to _n0 in the theoretical
!                          development
!
!    abs_cloud(nlr) (DP) - absorbtion by clouds (1/length) at a layer.
!                          This parameter is equaivalent to 
!                          _n0*alpha_ap in the theoretical development.
!
!    t_srf          (DP) - thermodynamic temperature of the surface
!
!    t_cb           (DP) - thermodynamic temperature of the cosmics
!                          background
!
!   surf_reflec(nangover2)      (DP) - surface reflectivity. The surface
!                                 is assmed to be plane in this
!                                 application. Basic subroutine
!                                 DTB23.F90 does not require this.
!
!   phaseff(nlr,nangover2,nangover2) (DP) - scattering matrix in the forward 
!                                 direction
! 
!   phasefb(nlr,nangover2,nangover2) (DP) - scattering matrix in the backward
!                                 direction
!
!   "d"-values representing derivatives of the have the corresponding
!       values with respect to appropriate parameters have the same
!       dimensions
!
!   eps_err                (DP) - max allowed relative error in brightness
!                                 temperature
!
!
!                           OUTPUT
!
!  tb_pl(0:nlr,nangover2)   (DP)  - upwelling brightness temperatures
!                              at a boundaries between layers at
!                              a particular angle ( k=0 corresponds
!                              to the surface, and k=nlr+1 to the
!                              upper boundary of the atmosphere)
!
!  tb_mn(0:nlr,nangover2)   (DP)  - downwelling brightness temperatures
!                              at a boundaries between layers.
!
! dtb_pl(0:nlr,nangover2,nvar) (DP)  - incremental btightness temperature
! dtb_mn(0:nlr,nangover2,nvar)         profiles, i.e. derivatives of the
!                                 Tb profiles with respect to a
!                                 parameter.
!
!------------------------------------------------------
use variables
IMPLICIT NONE

real(8), DIMENSION(  nlr1)               :: h
!real(8), DIMENSION(  nlr1  ,nangover2,nangover2)      ::  a0, b0
!real(8), DIMENSION(  nlr1  ,nangover2,nangover2,nvar) :: da0, db0
real(8), allocatable      ::  a0(:,:,:), b0(:,:,:)
real(8), allocatable :: da0(:,:,:,:), db0(:,:,:,:)

real(8), DIMENSION(0:nlr1+1,nangover2)        ::  f
real(8), DIMENSION(0:nlr1+1,nangover2,nvar)   :: df
real(8), DIMENSION(nangover2   ,nangover2)        :: r, dr
real(8), DIMENSION(0:nlr1  ,nangover2)        ::  u, v
real(8), DIMENSION(0:nlr1  ,nangover2,nvar)   :: du,dv
real(8), DIMENSION(nangover2)      :: t1,t2,pt,dpt, pt_sc
real(8), DIMENSION(nangover2,max_number_h2o_phases)      :: dpt_g, dpt_sc
real(8), DIMENSION(nangover2,nangover2) :: at, bt,dat,dbt, at_sc, bt_sc
real(8), DIMENSION(nangover2,nangover2,max_number_h2o_phases) :: dat_g, dbt_g, dat_sc, dbt_sc
LOGICAL         , DIMENSION(0:nlr1,nvar+3) :: LP
INTEGER :: i,j, ilr, ivar, hydrometeor_phase, jvar, jlr

double precision dt_srf, dt_cb, t_srf
real(8) albedo
real(8), DIMENSION(nangover2) :: dsurf_reflec

integer alloc_err

character*120 debugout

allocate( a0(nlr1,nangover2,nangover2), b0(nlr1,nangover2,nangover2) )
allocate( da0(nlr1,nangover2,nangover2,nvar), db0(nlr1,nangover2,nangover2,nvar) )

!f = 0.0d0
    DO i=1,nang
         cs(i)= DCOS(teta(i)/180.d0*pi)
    END DO
          t_srf = surf_inp%surf_temp
        dt_srf  = 0.d0
        dt_cb   = 0.d0
   dsurf_reflec = 0.d0

  ! linear or logarithmic derivative
   dscat_cloud1 = 1.0d0 ! scat_cloud1  ! change to dscat_cloud1 = 1.d0 after check !
       dal_gas1 = 1.0d0 !     al_gas1  ! change to     dal_gas1 = 1.d0 after check !
    dabs_cloud1 = 1.0d0 !  abs_cloud1  ! change to  dabs_cloud1 = 1.d0 after check !
  dtemperature1 = 1.0d0

  DO i=1,nangover2
    t1(i) = DSQRT(cris_quad_wghts(i)/cs(i))
    t2(i) = DSQRT(cris_quad_wghts(i)*cs(i))
  END DO
  
  DO ilr=1,nlr1
    !--------------  _at and _bt :
    DO i=1,nangover2
      DO j=i,nangover2
        at(i,j) = t1(i) * phaseff1(ilr,i,j) * t1(j)
        bt(i,j) = t1(i) * phasefb1(ilr,i,j) * t1(j)
        at(j,i) = at(i,j)
        bt(j,i) = bt(i,j)
        at_sc(i,j) = t1(i) * phaseff1_sc(ilr,i,j) * t1(j)
        bt_sc(i,j) = t1(i) * phasefb1_sc(ilr,i,j) * t1(j)
        at_sc(j,i) = at_sc(i,j)
        bt_sc(j,i) = bt_sc(i,j)
        do hydrometeor_phase = 1, number_h2o_phases
          dat_g(i,j,hydrometeor_phase)  = t1(i) * dphaseff1_g(ilr,i,j,hydrometeor_phase)  * t1(j)
          dbt_g(i,j,hydrometeor_phase)  = t1(i) * dphasefb1_g(ilr,i,j,hydrometeor_phase)  * t1(j)
          dat_sc(i,j,hydrometeor_phase) = t1(i) * dphaseff1_sc(ilr,i,j,hydrometeor_phase) * t1(j)
          dbt_sc(i,j,hydrometeor_phase) = t1(i) * dphasefb1_sc(ilr,i,j,hydrometeor_phase) * t1(j)
          dat_g(j,i,hydrometeor_phase)  = dat_g(i,j,hydrometeor_phase)
          dbt_g(j,i,hydrometeor_phase)  = dbt_g(i,j,hydrometeor_phase)
          dat_sc(j,i,hydrometeor_phase) = dat_sc(i,j,hydrometeor_phase)
          dbt_sc(j,i,hydrometeor_phase) = dbt_sc(i,j,hydrometeor_phase)
        end do ! hydrometeor_phase
      END DO
    END DO

    DO i = 1, nangover2
      pt(i)=0.d0
      pt_sc(i) = 0.0d0
      do hydrometeor_phase = 1, number_h2o_phases
        dpt_g(i,hydrometeor_phase)  = 0.0d0
        dpt_sc(i,hydrometeor_phase) = 0.0d0
      end do ! hydrometeor_phase
      DO j = 1, nangover2
        pt(i) =  pt(i) + ( at(j,i)+ bt(j,i))*t2(j)
        pt_sc(i) =  pt_sc(i) + ( at_sc(j,i)+ bt_sc(j,i))*t2(j)
        do hydrometeor_phase = 1, number_h2o_phases
          dpt_g(i,hydrometeor_phase)  =   dpt_g(i,hydrometeor_phase)              &
                                      + ( dat_g(j,i,hydrometeor_phase)            &
                                        + dbt_g(j,i,hydrometeor_phase) ) * t2(j)
          dpt_sc(i,hydrometeor_phase) =   dpt_sc(i,hydrometeor_phase)             &
                                      + ( dat_sc(j,i,hydrometeor_phase)           &
                                        + dbt_sc(j,i,hydrometeor_phase) ) * t2(j)
        end do ! hydrometeor_phase
      END DO
   END DO
   
   h(ilr)=altitude(ilr+1)-altitude(ilr)
   
    DO i=1,nangover2
      DO j=1,nangover2
        a0(ilr,i,j) = - at_sc(i,j)
        b0(ilr,i,j) = - bt_sc(i,j)
        !----new---------------------------------------------
        !   ivar=1         : variation of temperature 
        !   ivar=2         : variation of (gas_abs+cloud_abs)
        !   ivar=1+2*phase : variation of scattering by clouds (phase)
        !   ivar=2+2*phase : variation of HG assymetry parameter (phase)
        !----new---------------------------------------------
        !----old---------------------------------------------
        !   ivar=1 : variation of scattering by clouds
        !   ivar=2 : variation of (gas_abs+cloud_abs)
        !   ivar=3 : variation of temperature
        !   ivar=4 : variation of HG assymetry parameter
        !----old---------------------------------------------
        !----general case-------------------------------------------------
        !  da0(ilr,i,j)=-scat_cloud(ilr)*dat(i,j)-dscat_cloud(ilr)*at(i,j)
        !  db0(ilr,i,j)=-scat_cloud(ilr)*dbt(i,j)-dscat_cloud(ilr)*bt(i,j)
        !----general case-------------------------------------------------
        !----old----------------------------------------------------------
        ! da0(ilr,i,j,1) =                        - dscat_cloud1(ilr) * at(i,j)
        ! db0(ilr,i,j,1) =                        - dscat_cloud1(ilr) * bt(i,j)
        ! da0(ilr,i,j,2) = 0.d0      ! both dat=0,dbt=0 and dscat_cloud=0
        ! db0(ilr,i,j,2) = 0.d0
        ! da0(ilr,i,j,3) = 0.d0
        ! db0(ilr,i,j,3) = 0.d0
        ! da0(ilr,i,j,4) = - scat_cloud1(ilr) * dat(i,j)  !  dscat_cloud=0
        ! db0(ilr,i,j,4) = - scat_cloud1(ilr) * dbt(i,j)
        !----old----------------------------------------------------------
        da0(ilr,i,j,1) = 0.d0
        db0(ilr,i,j,1) = 0.d0
        da0(ilr,i,j,2) = 0.d0
        db0(ilr,i,j,2) = 0.d0
        do hydrometeor_phase = 1, number_h2o_phases
          ivar = 1 + 2 * hydrometeor_phase
          da0(ilr,i,j,ivar) = - dat_sc(i,j,hydrometeor_phase)
          db0(ilr,i,j,ivar) = - dbt_sc(i,j,hydrometeor_phase)
          ivar = 2 + 2 * hydrometeor_phase
          da0(ilr,i,j,ivar) = - dat_g(i,j,hydrometeor_phase)
          db0(ilr,i,j,ivar) = - dbt_g(i,j,hydrometeor_phase)
        end do ! hydrometeor_phase
        IF( i == j) THEN
          a0(ilr,i,i) = a0(ilr,i,i) + abs_total1(ilr) / cs(i) + pt_sc(i) / t2(i)
          !----general case---------------------------------------------------------------------
          ! da0(ilr,i,i)= da0(ilr,i,i) + (dal_gas(ilr)+dabs_cloud(ilr))/cs(i) &
          !                            +  scat_cloud(ilr)*dpt(i)/t2(i)        &
          !                            + dscat_cloud(ilr)* pt(i)/t2(i)      
          !----general case---------------------------------------------------------------------
          da0(ilr,i,i,2) =   da0(ilr,i,i,2)                            &
                         + ( dal_gas1(ilr) + dabs_cloud1(ilr) ) / cs(i)
          do hydrometeor_phase = 1, number_h2o_phases
            ivar = 1 + 2 * hydrometeor_phase
            da0(ilr,i,i,ivar) = da0(ilr,i,i,ivar)                &
                              + dscat_cloud1(ilr) * pt(i) / t2(i) ! this is why we need phase11
                                                                   ! without scat_cloud
            ivar = 2 + 2 * hydrometeor_phase
            da0(ilr,i,i,ivar) = da0(ilr,i,i,ivar) &
                              + dpt_g(i,hydrometeor_phase) / t2(i)
          end do ! hydrometeor_phase
        END IF    ! i == j

      END DO   !  loop over j

      f(ilr,i)= abs_total1(ilr) * temperature1(ilr) * t1(i)
      
      !----general case---------------------------------------------------
      !    df(ilr,i) = ( dal_gas(ilr) + dabs_cloud(ilr) ) *  temperature(ilr) * t1(i)  &
      !              + (  al_gas(ilr) +  abs_cloud(ilr) ) * dtemperature(ilr) * t1(i)
      !----general case---------------------------------------------------
      df(ilr,i,1) = abs_total1(ilr) * dtemperature1(ilr) * t1(i)
      df(ilr,i,2) = ( dal_gas1(ilr) + dabs_cloud1(ilr) ) *  temperature1(ilr) * t1(i)
      do hydrometeor_phase = 1, number_h2o_phases
        ivar = 1 + 2 * hydrometeor_phase
        df(ilr,i,ivar) = 0.d0
        ivar = 2 + 2 * hydrometeor_phase
        df(ilr,i,ivar) = 0.d0
      end do ! hydrometeor_phase
    END DO   !  loop over i
  END DO  ! loop over layers ilr
  
  r=0.d0
  dr=0.d0
  DO i=1,nangover2
     r(i,i) =  surf_reflec(i)
    dr(i,i) = dsurf_reflec(i)
  END DO
!  do i = 1, nangover2
!    do j = 1, nangover2
!      r(i,j) = surface_reflection(i,j)
!    end do
!  end do
  DO i=1,nangover2
    f(0     ,i) = (1- r(i,i)) * t_srf * t2(i)
    f(nlr1+1,i) =               t_cb  * t2(i)
    DO ivar=1,nvar
      df(0     ,i,ivar) =   -dr(i,i) *   t_srf * t2(i)   &
                        + (1- r(i,i)) * dt_srf * t2(i)
      df(nlr1+1,i,ivar) = dt_cb * t2(i)
    END DO
  END DO

!  do i = 1,nlr1
!     write(debugout,*) "absTotal, h = ", abs_total1(i),h(i)
!     call mexPrintf(debugout//achar(10))
!  end do

 
  !------------------- defining diagonal layers (04/24/03) ------------
  LP = .false. ! all layers/perturbations are considered non-diagonal
! blw 911     11 November 2004
! This apparently does not work properly in parallel processing
! activated 9 March 2005
! problem was in core95 - off-diagonal elements set to zero now
! when layer is diagonal
!! There is a problem with diagonalizing layers ! blw  20 May 2005
  DO ilr=1,nlr1
!!    IF( scat_cloud1(ilr) .le. eps_err*abs_total1(ilr) ) THEN
    albedo = scat_cloud1(ilr) / ( scat_cloud1(ilr) + abs_total1(ilr) )
    IF( albedo .le. d_albedo ) THEN
      LP(ilr,3)=.true.   ! diagonal layer
      do ivar = 1, nvar
        LP(ilr,3+ivar)=.true.
      end do
    END IF
  END DO
! blw 911

!write(debugout,*) "nlr1=",nlr1
!call mexPrintf(debugout//achar(10))

!write(debugout,*) "a0"
!call mexPrintf(debugout//achar(10))

!do i=1,nangover2
!   do j=1,nangover2
!      write(debugout,*) a0(1,i,j)
!      call mexPrintf(debugout)
!   end do
!   call mexPrintf(achar(10))
!end do

!write(debugout,*) "b0"
!call mexPrintf(debugout//achar(10))

!do i=1,nangover2
!   do j=1,nangover2
!      write(debugout,*) b0(1,i,j)
!      call mexPrintf(debugout)
!   end do
!   call mexPrintf(achar(10))
!end do

!write(debugout,*) "r="
!call mexPrintf(debugout//achar(10))
!do i=1,nangover2
!   write(debugout,*) r(i,i)
!   call mexPrintf(debugout)
!end do
!call mexPrintf(achar(10))

!write(debugout,*) "f="
!call mexPrintf(debugout//achar(10))
!do j=0,nlr1
!   do i=1,nangover2
!      write(debugout,*) f(j,i)
!      call mexPrintf(debugout)
!   end do
!   call mexPrintf(achar(10))
!end do

!-------------------------------------------------------------------- 
  CALL core95 (nlr1,nangover2,h,a0,b0,f,r,u,v,da0,db0,df,dr,du,dv,obs_lev1,nvar,LP)
!----------------------------------    

!write(debugout,*) "u="
!call mexPrintf(debugout//achar(10))
!do i=1,nang
!   write(debugout,*) u(1,i)
!   call mexPrintf(debugout)
!end do
!call mexPrintf(achar(10))

  
  DO ilr = 0,nlr1
    DO i=1,nangover2
      tb_pl(ilr,i)     = u(ilr,i)/t2(i)
      tb_mn(ilr,i)     = v(ilr,i)/t2(i)
      DO ivar = 1, nvar
        dtb_pl(ilr,i,ivar) = du(ilr,i,ivar) / t2(i)
        dtb_mn(ilr,i,ivar) = dv(ilr,i,ivar) / t2(i)
        ! correction for absorption derivative
        if( ivar == 2 ) then
          jlr = ilr
          if( ilr == 0 ) jlr = 1
          dtb_pl(ilr,i,ivar) = dtb_pl(ilr,i,ivar) * 0.25d0 * da0(jlr,i,i,ivar)
          dtb_mn(ilr,i,ivar) = dtb_mn(ilr,i,ivar) * 0.25d0 * db0(jlr,i,i,ivar) !?
        end if
      END DO
    END DO
  END DO

deallocate(a0,b0,stat=alloc_err)
deallocate(da0,db0,stat=alloc_err)

END SUBROUTINE tb94
