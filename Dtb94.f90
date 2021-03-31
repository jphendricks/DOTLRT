!=======================================================================
SUBROUTINE dtb94 ()
!=======================================================================
! This subroutine further prepares data for tb94 which calls CORE95.F90 which
! then calculates brightness temperature and the radiative Jacobian.
! Below, that Jacobian is called the incremental profiles. 
!
! Issues:
!  12/12/2020 (Kevin Schaefer) cannot remove psuedo levelsinterpolation
!
! History:
!  12/17/2002 blw updated routine from DTB90.F90 to DTB93.F90
!  12/17/2002 blw added interolation between layers
!  4/24/2003  blw updated routine from DTB93.F90 to DTB93.F90
!  4/24/2003  blw added eps_err
!  9/26/2020  Kevin Schaefer deleted unused variables
!  12/12/2020 Kevin schaefer removed all outdated comments and variable defs
!  12/12/2020 Kevin Schaefer removed external tb94 statement
!-----------------------------------------------------------------------
use dotlrt_variables
implicit none

integer i
integer j
integer k
integer k1
integer k2
integer iphase ! (-) hydrometeor phase index
integer ilev   ! (-) layer or level index
real(8) h(nlev+1)
real(8) k1m1
real(8) m1k1m1
real(8) k2m1
real(8) m1k2m1

  DO ilev=1,nlev
    h(ilev)=altitude(ilev+1)-altitude(ilev)
  enddo

  !----  the first layer: no interpolation 
  altitude1(1)=altitude(1)
  k=0
  DO k1=1,npseudo
    k = k+1
    k1m1 = dble(k1) / dble(npseudo)
    m1k1m1 = 1.0d0 - k1m1
    altitude1(k+1) = m1k1m1 * altitude(1) + k1m1 * altitude(2)
    h1(k) =  altitude1(k+1) - altitude1(k)
    temperature1(k) = temperature(1)
    abs_O2_1(k) = abs_O2(1)
    abs_H2O_1(k) = abs_H2O(1)
    abs_cloud1(k) = abs_cloud(1)
    scat_cloud1(k) = scat_cloud(1)
    al_gas1(k) = abs_O2_1(1) + abs_H2O_1(1)
    abs_total1(k) = al_gas1(k) + abs_cloud1(k)
    DO i=1,nangover2
      DO j=1,nangover2
        phaseff1(k,i,j) = phaseff(1,i,j)
        phasefb1(k,i,j) = phasefb(1,i,j)
        phaseff1_sc(k,i,j) = phaseff_sc(1,i,j)
        phasefb1_sc(k,i,j) = phasefb_sc(1,i,j)
        do iphase = 1, nphase
          dphasefb1_g(k,i,j,iphase)  = dphasefb_g(1,i,j,iphase)
          dphaseff1_g(k,i,j,iphase)  = dphaseff_g(1,i,j,iphase)
          dphasefb1_sc(k,i,j,iphase) = dphasefb_sc(1,i,j,iphase)
          dphaseff1_sc(k,i,j,iphase) = dphaseff_sc(1,i,j,iphase)
        enddo
      enddo
    enddo

   enddo

! higher layers: linear interpolation 
  DO ilev=2,nlev
    DO k1=1,npseudo
      k=k+1
      k1m1 = dble(k1) / dble(npseudo)
      m1k1m1 = 1.0d0 - k1m1
      altitude1(k+1) = m1k1m1 * altitude(ilev) + k1m1 * altitude(ilev+1)

      ! linear iterpolation: if one sets k2=npseudo, then all parameters
      ! at pseudo-levels will be the same as at the basic levels. This
      ! can be used for checking purposes. The choise k2=k1 leads to
      ! linear interpolation.
      k2=k1                          !   linear interpolation in effect
      h1(k) =   altitude1(k+1)-altitude1(k)

      k2m1 = dble(k2) / dble(npseudo)
      m1k2m1 = 1.0d0 - k2m1

      temperature1(k) = m1k2m1 * temperature(ilev-1) + k2m1 * temperature(ilev)
      abs_O2_1(k) = m1k2m1 *      abs_O2(ilev-1) + k2m1 *      abs_O2(ilev)
      abs_H2O_1(k) = m1k2m1 *     abs_H2O(ilev-1) + k2m1 *     abs_H2O(ilev)
      abs_cloud1(k) = m1k2m1 *   abs_cloud(ilev-1) + k2m1 *   abs_cloud(ilev)
      scat_cloud1(k) = m1k2m1 *  scat_cloud(ilev-1) + k2m1 *  scat_cloud(ilev)

      al_gas1(k) = abs_O2_1(k)+abs_H2O_1(k)
      abs_total1(k) = al_gas1(k) + abs_cloud1(k)
      DO i=1,nangover2
        DO j=1,nangover2
          phaseff1(k,i,j)    = m1k2m1 * phaseff(ilev-1,i,j)    + k2m1 * phaseff(ilev,i,j)
          phasefb1(k,i,j)    = m1k2m1 * phasefb(ilev-1,i,j)    + k2m1 * phasefb(ilev,i,j)
          phaseff1_sc(k,i,j) = m1k2m1 * phaseff_sc(ilev-1,i,j) + k2m1 * phaseff_sc(ilev,i,j)
          phasefb1_sc(k,i,j) = m1k2m1 * phasefb_sc(ilev-1,i,j) + k2m1 * phasefb_sc(ilev,i,j)
          do iphase = 1, nphase
            dphasefb1_g(k,i,j,iphase)  = m1k2m1 * dphasefb_g(ilev-1,i,j,iphase)  &
                                                  +   k2m1 * dphasefb_g(ilev,i,j,iphase)
            dphaseff1_g(k,i,j,iphase)  = m1k2m1 * dphaseff_g(ilev-1,i,j,iphase)  &
                                                  +   k2m1 * dphaseff_g(ilev,i,j,iphase)
            dphasefb1_sc(k,i,j,iphase) = m1k2m1 * dphasefb_sc(ilev-1,i,j,iphase) &
                                                  +   k2m1 * dphasefb_sc(ilev,i,j,iphase)
            dphaseff1_sc(k,i,j,iphase) = m1k2m1 * dphaseff_sc(ilev-1,i,j,iphase) &
                                                  +   k2m1 * dphaseff_sc(ilev,i,j,iphase)
          enddo ! iphase
        enddo
      enddo
    enddo ! npsuedo
  enddo ! ilev

!-------  end of forming interpolated values

 CALL  tb94()

END SUBROUTINE dtb94
