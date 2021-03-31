! Jacobian.f90
! The Total Jacobian is derived using both the Geophysical Jacobian
! and the Radiative Jacobian by means of the differentiation chain rule.

!   Total Jacobian:
!     1-temperature  dTb_dT(level,i0,jud)
!     2-pressure     dTb_dp(level,i0,jud)
!     3-water vapor  dTb_dq(level,i0,jud)
!     4-liquid       dTb_dw(level,i0,1,jud)    hydrometeor phase 1
!     5-rain         dTb_dw(level,i0,2,jud)    hydrometeor phase 2
!     6-ice          dTb_dw(level,i0,3,jud)    hydrometeor phase 3
!     7-snow         dTb_dw(level,i0,4,jud)    hydrometeor phase 4
!     8-graupel      dTb_dw(level,i0,5,jud)    hydrometeor phase 5
! jud = 1, 2 (up, down)

!   Radiative Jacobian:
!     dTb / dT   = dtb_pl(ilr1,i0,1)    variation by temperature
!     dTb / dKab = dtb_pl(ilr1,i0,2)    variation by absorption by hydrometeors
!     dTb / dKsc = dtb_pl(ilr1,i0,jsc)  variation by scattering by hydrometeors
!                                       jsc = 3, 5, 7, 9, 11
!     dTb / dg   = dtb_pl(ilr1,i0,jg)   variation by asymmetry in hydrometeor scatter
!                                       jg = 4, 6, 8, 10, 12

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

! subroutine Jacobian(geo_unit, atm_data)
subroutine Jacobian()
    use variables
!    use type_kinds
!    use readprofile
    implicit none
!    integer geo_unit
    integer ilr1, hydrometeor_phase, jvar_sc, jvar_g
!    real(8) z

    ! total Jacobian

  do i0 = 1, nangover2
    do ilr1 = 1, nlr1
      dTb_dT(ilr1,i0,1) = dtb_pl(ilr1,i0,1) + dKab_dT(ilr1) * dtb_pl(ilr1,i0,2)
      ! hydrometeor phases
      do hydrometeor_phase = 1, number_h2o_phases
        jvar_sc = 1 + 2 * hydrometeor_phase
        jvar_g  = 2 + 2 * hydrometeor_phase
        dTb_dT(ilr1,i0,1) = dTb_dT(ilr1,i0,1)                                              &
                     + dKsc_dT(ilr1,hydrometeor_phase) * dtb_pl(ilr1,i0,jvar_sc) &
                     +  dg_dT(ilr1,hydrometeor_phase) * dtb_pl(ilr1,i0,jvar_g)
      end do ! hydrometeor_phase
      dTb_dp(ilr1,i0,1) = dKab_dp(ilr1) * dtb_pl(ilr1,i0,2)
      dTb_dq(ilr1,i0,1) = dKab_dq(ilr1) * dtb_pl(ilr1,i0,2)
      ! hydrometeor phases
      do hydrometeor_phase = 1, number_h2o_phases
        jvar_sc = 1 + 2 * hydrometeor_phase
        jvar_g  = 2 + 2 * hydrometeor_phase
        dTb_dw(ilr1,i0,hydrometeor_phase,1) = dKab_dw(ilr1,hydrometeor_phase) * dtb_pl(ilr1,i0,2)       &
                                          + dKsc_dw(ilr1,hydrometeor_phase) * dtb_pl(ilr1,i0,jvar_sc) &
                                          +   dg_dw(ilr1,hydrometeor_phase) * dtb_pl(ilr1,i0,jvar_g)
      end do ! hydrometeor_phase
    end do ! ilr1

!    write(geo_unit,'(2(x,e16.8),i8,x,e16.8)') atm_data%lat, atm_data%lng, nlr1, passband_freq(nfreq)
!    z=0.d0
!    do ilr1=1,nlr1
!      z=z+h1(ilr1)
!      write(geo_unit,'(9(x,e16.8))') z                                , &
!        dTb_dT(ilr1,i0)                                        / h1(ilr1), &
!        dTb_dp(ilr1,i0)   * atm_inp%prof(ilr1)%pressure        / h1(ilr1), & ! derivative with respect to ln(pressure)
!        dTb_dq(ilr1,i0)   * atm_inp%prof(ilr1)%vapor_density   / h1(ilr1), & ! derivative with respect to ln(water vapor density)
!        dTb_dw(ilr1,i0,1) * atm_inp%prof(ilr1)%cloud_liq_dens  / h1(ilr1), & ! derivative with respect to ln(condensed  liquid water density)
!        dTb_dw(ilr1,i0,2) * atm_inp%prof(ilr1)%cloud_rn_dens   / h1(ilr1), & ! derivative with respect to ln(condensed    rain water density)
!        dTb_dw(ilr1,i0,3) * atm_inp%prof(ilr1)%cloud_ice_dens  / h1(ilr1), & ! derivative with respect to ln(condensed     ice water density)
!        dTb_dw(ilr1,i0,4) * atm_inp%prof(ilr1)%cloud_snow_dens / h1(ilr1), & ! derivative with respect to ln(condensed    snow water density)
!        dTb_dw(ilr1,i0,5) * atm_inp%prof(ilr1)%cloud_grpl_dens / h1(ilr1)    ! derivative with respect to ln(condesned graupel water density)
!    end do
  end do ! i0

end subroutine Jacobian
