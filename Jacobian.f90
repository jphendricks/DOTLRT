!=======================================================================================
subroutine Jacobian()
!=======================================================================================
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
! History:
!  10/17/2020 Kevin Schaefer removed all commented code
!----------------------------------------------------------------------------------------
    use dotlrt_variables
    implicit none
    integer ilr1, iphase, jvar_sc, jvar_g

    ! total Jacobian

  do i0 = 1, nangover2
    do ilr1 = 1, nlr1
      dTb_dT(ilr1,i0,1) = dtb_pl(ilr1,i0,1) + dKab_dT(ilr1) * dtb_pl(ilr1,i0,2)
      ! hydrometeor phases
      do iphase = 1, nphase
        jvar_sc = 1 + 2 * iphase
        jvar_g  = 2 + 2 * iphase
        dTb_dT(ilr1,i0,1) = dTb_dT(ilr1,i0,1)                                              &
                     + dKsc_dT(ilr1,iphase) * dtb_pl(ilr1,i0,jvar_sc) &
                     +  dg_dT(ilr1,iphase) * dtb_pl(ilr1,i0,jvar_g)
      end do ! iphase
      dTb_dp(ilr1,i0,1) = dKab_dp(ilr1) * dtb_pl(ilr1,i0,2)
      dTb_dq(ilr1,i0,1) = dKab_dq(ilr1) * dtb_pl(ilr1,i0,2)
      ! hydrometeor phases
      do iphase = 1, nphase
        jvar_sc = 1 + 2 * iphase
        jvar_g  = 2 + 2 * iphase
        dTb_dw(ilr1,i0,iphase,1) = dKab_dw(ilr1,iphase) * dtb_pl(ilr1,i0,2)       &
                                          + dKsc_dw(ilr1,iphase) * dtb_pl(ilr1,i0,jvar_sc) &
                                          +   dg_dw(ilr1,iphase) * dtb_pl(ilr1,i0,jvar_g)
      end do ! iphase
    end do ! ilr1

  end do ! i0

end subroutine Jacobian
