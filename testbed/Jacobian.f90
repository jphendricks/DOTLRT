!=======================================================================================
subroutine Jacobian()
!=======================================================================================
! The Total Jacobian is derived using both the Geophysical Jacobian
! and the Radiative Jacobian by means of the differentiation chain rule.

!   Total Jacobian:
!     1-temperature  dTb_dT(level,iang,jud)
!     2-pressure     dTb_dp(level,iang,jud)
!     3-water vapor  dTb_dq(level,iang,jud)
!     4-liquid       dTb_dw(level,iang,1,jud)    hydrometeor phase 1
!     5-rain         dTb_dw(level,iang,2,jud)    hydrometeor phase 2
!     6-ice          dTb_dw(level,iang,3,jud)    hydrometeor phase 3
!     7-snow         dTb_dw(level,iang,4,jud)    hydrometeor phase 4
!     8-graupel      dTb_dw(level,iang,5,jud)    hydrometeor phase 5
! jud = 1, 2 (up, down)

!   Radiative Jacobian:
!     dTb / dT   = dtb_pl(ilev,iang,1)    variation by temperature
!     dTb / dKab = dtb_pl(ilev,iang,2)    variation by absorption by hydrometeors
!     dTb / dKsc = dtb_pl(ilev,iang,jsc)  variation by scattering by hydrometeors
!                                       jsc = 3, 5, 7, 9, 11
!     dTb / dg   = dtb_pl(ilev,iang,jg)   variation by asymmetry in hydrometeor scatter
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
!  12/12/2020 Kevin Schaefer changed to standard indexing variables
!----------------------------------------------------------------------------------------
  use dotlrt_variables
  implicit none

! internal variables
  integer jvar_sc
  integer jvar_g
  integer ilev   ! (-) vertical level index
  integer iphase ! (-) hydrometeor phase index
  integer iang   ! (-) stream angle index

! total Jacobian
  do iang = 1, nangover2
    do ilev = 1, nlev
      dTb_dT(ilev,iang,1) = dtb_pl(ilev,iang,1) + dKab_dT(ilev) * dtb_pl(ilev,iang,2)

      
      do iphase = 1, nphase ! hydrometeor phases
        jvar_sc = 1 + 2 * iphase
        jvar_g  = 2 + 2 * iphase
        dTb_dT(ilev,iang,1) = dTb_dT(ilev,iang,1) + dKsc_dT(ilev,iphase) * dtb_pl(ilev,iang,jvar_sc) &
          +  dg_dT(ilev,iphase) * dtb_pl(ilev,iang,jvar_g)
      end do
      dTb_dp(ilev,iang,1) = dKab_dp(ilev) * dtb_pl(ilev,iang,2)
      dTb_dq(ilev,iang,1) = dKab_dq(ilev) * dtb_pl(ilev,iang,2)

      do iphase = 1, nphase ! hydrometeor phases
        jvar_sc = 1 + 2 * iphase
        jvar_g  = 2 + 2 * iphase
        dTb_dw(ilev,iang,iphase,1) = dKab_dw(ilev,iphase) * dtb_pl(ilev,iang,2) &
          + dKsc_dw(ilev,iphase) * dtb_pl(ilev,iang,jvar_sc) &
          + dg_dw(ilev,iphase) * dtb_pl(ilev,iang,jvar_g)
      end do
    end do
  end do

end subroutine Jacobian
