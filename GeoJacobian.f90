!=======================================================================
subroutine GeoJacobian()
!=======================================================================
! GeoJacobian.f90
!   Geophysical Jacobian:
!     1- total absorption (gaseous + hydrometeor absorption) temperature derivative:
!        dKab_dT(level)
!     2- total absorption (gaseous only) pressure derivative:
!        dKab_dp(level)
!     3- total absorption (gaseous only) water vapor derivative:
!        dKab_dq(level)
!     4- hydrometeor scatter temperature derivative by phase:
!        dKsc_dT(level,iphase)
!     5- hydrometeor asymmetry temperature derivative by phase:
!        dg_dT(level,iphase)
!     6- total absorption (hydrometeor only) water density derivative by phase:
!        dKab_dw(level,iphase)
!     7- hydrometeor scatter water density derivative by phase:
!        dKsc_dw(level,iphase)
!     8- hydrometeor asymmetry water density derivative by phase:
!        dg_dw(level,iphase)
!
!-----------------------------------------------------------------------
! History:
!   9/26/2020 Kevin Schaefer deleted unused variables
!  12/12/2020 Kevin Schaefer changed to standard indexing variables
!-----------------------------------------------------------------------

use dotlrt_variables
implicit none

! internal variables
  integer iphase ! (-) hydrometeor phase index
  integer ilev   ! (-) vertical level index

    ! geophysical Jacobian
    do ilev = 1, nlev
      ! gaseous + hydrometeor absorption - temperature derivative
      dKab_dT(ilev) =   gas_prof(ilev)%dabsn2_dt     &
                    +   gas_prof(ilev)%dabsh2o_dt    &
                    +   gas_prof(ilev)%do2abs_dt     &
                    + hydro_prof(ilev,1)%dcloudab_dt &
                    + hydro_prof(ilev,2)%dcloudab_dt &
                    + hydro_prof(ilev,3)%dcloudab_dt &
                    + hydro_prof(ilev,4)%dcloudab_dt &
                    + hydro_prof(ilev,5)%dcloudab_dt
      ! gaseous absorption - pressure derivative
      dKab_dp(ilev) = gas_prof(ilev)%dabsn2_dp  &
                    + gas_prof(ilev)%dabsh2o_dp &
                    + gas_prof(ilev)%do2abs_dp
      ! gaseous absorption - water vapor derivative
      dKab_dq(ilev) = gas_prof(ilev)%dabsh2o_dw &
                    + gas_prof(ilev)%do2abs_dw
      ! hydrometeor phases
      do iphase = 1, nphase
        ! hydrometeor scatter - temperature derivative
        dKsc_dT(ilev,iphase) = hydro_prof(ilev,iphase)%dcloudsc_dt
        ! hydrometeor asymmetry - temperature derivative
        dg_dT(ilev,iphase) = hydro_prof(ilev,iphase)%dcloudg_dt
      end do

      ! iphase = 1 (cloud liquid) a0 = constant
      dKab_dw(ilev,1) = atm(ilev)%dclw_k0_dw  * hydro_prof(ilev,1)%dcloudab_dk0
      dKsc_dw(ilev,1) = atm(ilev)%dclw_k0_dw  * hydro_prof(ilev,1)%dcloudsc_dk0
      dg_dw(ilev,1) = atm(ilev)%dclw_k0_dw    * hydro_prof(ilev,1)%dcloudg_dk0

      ! iphase = 2 (rain) k0 = constant
      dKab_dw(ilev,2) = atm(ilev)%drain_a0_dw  * hydro_prof(ilev,2)%dcloudab_da0
      dKsc_dw(ilev,2) = atm(ilev)%drain_a0_dw  * hydro_prof(ilev,2)%dcloudsc_da0
      dg_dw(ilev,2) = atm(ilev)%drain_a0_dw    * hydro_prof(ilev,2)%dcloudg_da0

      ! iphase = 3 (ice) a0 = constant
      dKab_dw(ilev,3) = atm(ilev)%dice_k0_dw  * hydro_prof(ilev,3)%dcloudab_dk0
      dKsc_dw(ilev,3) = atm(ilev)%dice_k0_dw  * hydro_prof(ilev,3)%dcloudsc_dk0
      dg_dw(ilev,3) = atm(ilev)%dice_k0_dw    * hydro_prof(ilev,3)%dcloudg_dk0

      ! iphase = 4 (snow) k0 = constant
      dKab_dw(ilev,4) = atm(ilev)%dsnow_a0_dw  * hydro_prof(ilev,4)%dcloudab_da0
      dKsc_dw(ilev,4) = atm(ilev)%dsnow_a0_dw  * hydro_prof(ilev,4)%dcloudsc_da0
      dg_dw(ilev,4) = atm(ilev)%dsnow_a0_dw    * hydro_prof(ilev,4)%dcloudg_da0

      ! iphase = 5 (graupel) k0 = constant
      dKab_dw(ilev,5) = atm(ilev)%dgrpl_a0_dw  * hydro_prof(ilev,5)%dcloudab_da0
      dKsc_dw(ilev,5) = atm(ilev)%dgrpl_a0_dw  * hydro_prof(ilev,5)%dcloudsc_da0 
      dg_dw(ilev,5) = atm(ilev)%dgrpl_a0_dw    * hydro_prof(ilev,5)%dcloudg_da0
    end do ! ilev

end subroutine GeoJacobian
