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
!-----------------------------------------------------------------------

    use dotlrt_variables
    implicit none
    integer ilr1, iphase

    ! geophysical Jacobian
    do ilr1 = 1, nlr1
      ! gaseous + hydrometeor absorption - temperature derivative
      dKab_dT(ilr1) =   gas_prof(ilr1)%dabsn2_dt     &
                    +   gas_prof(ilr1)%dabsh2o_dt    &
                    +   gas_prof(ilr1)%do2abs_dt     &
                    + hydro_prof(ilr1,1)%dcloudab_dt &
                    + hydro_prof(ilr1,2)%dcloudab_dt &
                    + hydro_prof(ilr1,3)%dcloudab_dt &
                    + hydro_prof(ilr1,4)%dcloudab_dt &
                    + hydro_prof(ilr1,5)%dcloudab_dt
      ! gaseous absorption - pressure derivative
      dKab_dp(ilr1) = gas_prof(ilr1)%dabsn2_dp  &
                    + gas_prof(ilr1)%dabsh2o_dp &
                    + gas_prof(ilr1)%do2abs_dp
      ! gaseous absorption - water vapor derivative
      dKab_dq(ilr1) = gas_prof(ilr1)%dabsh2o_dw &
                    + gas_prof(ilr1)%do2abs_dw
      ! hydrometeor phases
      do iphase = 1, nphase
        ! hydrometeor scatter - temperature derivative
        dKsc_dT(ilr1,iphase) = hydro_prof(ilr1,iphase)%dcloudsc_dt
        ! hydrometeor asymmetry - temperature derivative
          dg_dT(ilr1,iphase) = hydro_prof(ilr1,iphase)%dcloudg_dt
      end do
      ! iphase = 1 (cloud liquid) a0 = constant
      dKab_dw(ilr1,1) = atm(ilr1)%dclw_k0_dw  &
                      * hydro_prof(ilr1,1)%dcloudab_dk0
      dKsc_dw(ilr1,1) = atm(ilr1)%dclw_k0_dw  &
                      * hydro_prof(ilr1,1)%dcloudsc_dk0
        dg_dw(ilr1,1) = atm(ilr1)%dclw_k0_dw  &
                      * hydro_prof(ilr1,1)%dcloudg_dk0
      ! iphase = 2 (rain) k0 = constant
      dKab_dw(ilr1,2) = atm(ilr1)%drain_a0_dw  &
                      * hydro_prof(ilr1,2)%dcloudab_da0
      dKsc_dw(ilr1,2) = atm(ilr1)%drain_a0_dw  &
                      * hydro_prof(ilr1,2)%dcloudsc_da0
        dg_dw(ilr1,2) = atm(ilr1)%drain_a0_dw  &
                      * hydro_prof(ilr1,2)%dcloudg_da0
      ! iphase = 3 (ice) a0 = constant
      dKab_dw(ilr1,3) = atm(ilr1)%dice_k0_dw  &
                      * hydro_prof(ilr1,3)%dcloudab_dk0
      dKsc_dw(ilr1,3) = atm(ilr1)%dice_k0_dw  &
                      * hydro_prof(ilr1,3)%dcloudsc_dk0
        dg_dw(ilr1,3) = atm(ilr1)%dice_k0_dw  &
                      * hydro_prof(ilr1,3)%dcloudg_dk0
      ! iphase = 4 (snow) k0 = constant
      dKab_dw(ilr1,4) = atm(ilr1)%dsnow_a0_dw  &
                      * hydro_prof(ilr1,4)%dcloudab_da0
      dKsc_dw(ilr1,4) = atm(ilr1)%dsnow_a0_dw  &
                      * hydro_prof(ilr1,4)%dcloudsc_da0
        dg_dw(ilr1,4) = atm(ilr1)%dsnow_a0_dw  &
                      * hydro_prof(ilr1,4)%dcloudg_da0
      ! iphase = 5 (graupel) k0 = constant
      dKab_dw(ilr1,5) = atm(ilr1)%dgrpl_a0_dw  &
                      * hydro_prof(ilr1,5)%dcloudab_da0
      dKsc_dw(ilr1,5) = atm(ilr1)%dgrpl_a0_dw  &
                      * hydro_prof(ilr1,5)%dcloudsc_da0 
        dg_dw(ilr1,5) = atm(ilr1)%dgrpl_a0_dw  &
                      * hydro_prof(ilr1,5)%dcloudg_da0
    end do ! ilr1

end subroutine GeoJacobian
