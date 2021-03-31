! GeoJacobian.f90
!   Geophysical Jacobian:
!     1- total absorption (gaseous + hydrometeor absorption) temperature derivative:
!        dKab_dT(level)
!     2- total absorption (gaseous only) pressure derivative:
!        dKab_dp(level)
!     3- total absorption (gaseous only) water vapor derivative:
!        dKab_dq(level)
!     4- hydrometeor scatter temperature derivative by phase:
!        dKsc_dT(level,hydrometeor_phase)
!     5- hydrometeor asymmetry temperature derivative by phase:
!        dg_dT(level,hydrometeor_phase)
!     6- total absorption (hydrometeor only) water density derivative by phase:
!        dKab_dw(level,hydrometeor_phase)
!     7- hydrometeor scatter water density derivative by phase:
!        dKsc_dw(level,hydrometeor_phase)
!     8- hydrometeor asymmetry water density derivative by phase:
!        dg_dw(level,hydrometeor_phase)

subroutine GeoJacobian()
    use variables
!    use type_kinds
!    use readprofile
    implicit none
    integer ilr1, hydrometeor_phase, jvar_sc, jvar_g
    real(8) z

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
      do hydrometeor_phase = 1, number_h2o_phases
        ! hydrometeor scatter - temperature derivative
        dKsc_dT(ilr1,hydrometeor_phase) = hydro_prof(ilr1,hydrometeor_phase)%dcloudsc_dt
        ! hydrometeor asymmetry - temperature derivative
          dg_dT(ilr1,hydrometeor_phase) = hydro_prof(ilr1,hydrometeor_phase)%dcloudg_dt
      end do
      ! hydrometeor_phase = 1 (cloud liquid) a0 = constant
      dKab_dw(ilr1,1) = atm_inp%prof(ilr1)%dcloud_liq_k0_dw  &
                      * hydro_prof(ilr1,1)%dcloudab_dk0
      dKsc_dw(ilr1,1) = atm_inp%prof(ilr1)%dcloud_liq_k0_dw  &
                      * hydro_prof(ilr1,1)%dcloudsc_dk0
        dg_dw(ilr1,1) = atm_inp%prof(ilr1)%dcloud_liq_k0_dw  &
                      * hydro_prof(ilr1,1)%dcloudg_dk0
      ! hydrometeor_phase = 2 (rain) k0 = constant
      dKab_dw(ilr1,2) = atm_inp%prof(ilr1)%dcloud_rn_a0_dw  &
                      * hydro_prof(ilr1,2)%dcloudab_da0
      dKsc_dw(ilr1,2) = atm_inp%prof(ilr1)%dcloud_rn_a0_dw  &
                      * hydro_prof(ilr1,2)%dcloudsc_da0
        dg_dw(ilr1,2) = atm_inp%prof(ilr1)%dcloud_rn_a0_dw  &
                      * hydro_prof(ilr1,2)%dcloudg_da0
      ! hydrometeor_phase = 3 (ice) a0 = constant
      dKab_dw(ilr1,3) = atm_inp%prof(ilr1)%dcloud_ice_k0_dw  &
                      * hydro_prof(ilr1,3)%dcloudab_dk0
      dKsc_dw(ilr1,3) = atm_inp%prof(ilr1)%dcloud_ice_k0_dw  &
                      * hydro_prof(ilr1,3)%dcloudsc_dk0
        dg_dw(ilr1,3) = atm_inp%prof(ilr1)%dcloud_ice_k0_dw  &
                      * hydro_prof(ilr1,3)%dcloudg_dk0
      ! hydrometeor_phase = 4 (snow) k0 = constant
      dKab_dw(ilr1,4) = atm_inp%prof(ilr1)%dcloud_snow_a0_dw  &
                      * hydro_prof(ilr1,4)%dcloudab_da0
      dKsc_dw(ilr1,4) = atm_inp%prof(ilr1)%dcloud_snow_a0_dw  &
                      * hydro_prof(ilr1,4)%dcloudsc_da0
        dg_dw(ilr1,4) = atm_inp%prof(ilr1)%dcloud_snow_a0_dw  &
                      * hydro_prof(ilr1,4)%dcloudg_da0
      ! hydrometeor_phase = 5 (graupel) k0 = constant
      dKab_dw(ilr1,5) = atm_inp%prof(ilr1)%dcloud_grpl_a0_dw  &
                      * hydro_prof(ilr1,5)%dcloudab_da0
      dKsc_dw(ilr1,5) = atm_inp%prof(ilr1)%dcloud_grpl_a0_dw  &
                      * hydro_prof(ilr1,5)%dcloudsc_da0 
        dg_dw(ilr1,5) = atm_inp%prof(ilr1)%dcloud_grpl_a0_dw  &
                      * hydro_prof(ilr1,5)%dcloudg_da0
    end do ! ilr1

end subroutine GeoJacobian
