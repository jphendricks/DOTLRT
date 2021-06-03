!===================================================================================
SUBROUTINE do_tb_gvh94()
!===================================================================================
! this routine calculates the Henyey-Greenstein (HG) phase function using
! linear interpolation of a lookup table
!
! History:
!  4/24/2003  added eps_err max relative error in brightness temperature allowed
!  12/16/2002 pre-calculated HG phase matrix into lookup table
!  9/26/2020  Kevin Schaefer deleted unused variables
!  10/15/2020 Kevin Schaefer deleted commented code
!  12/12/2020 Kevin Schaefer removed external dtb94 statement
!  12/12/2020 Kevin schaefer removed all outdated comments and variable defs
!  5/19/2021  Kevin Schaefer changed left interp index to eliminate Jagged Jacobians
!------------------------------------------------------------------------------------
use dotlrt_variables
implicit none

! local variables
  integer iphase    ! (-) hydrometeor phase index
  integer i         ! (-) stream angle index
  integer j         ! (-) stream angle index
  integer i_left    ! (-) left interpolation index
  integer i_right   ! (-) right interpolation index
  integer ilev      ! (-) layer index
  real(8) g         ! (-) local version asymmetry
  real(8) HGgs(2)   ! (-) left and right value 
  real(8) HGphs(2)  ! (-) left and right value
  real(8) dHGphs(2) ! (-) left and right value
  real(8) HGs       ! (-) interpolated HG phase
  real(8) dHGs      ! (-) interpolated HG phase Jacobian

! Transfer to separate variables
  do i=1, nlev
    g_asymmetry(i) = atm(i)%asymmetry 
    altitude(i) = atm(i)%hgt_mid
    temperature(i) = atm(i)%temp
    abs_O2(i) = atm(i)%abs_O2
    abs_H2O(i) = atm(i)%abs_H2O
    abs_cloud(i) = atm(i)%abs_cloud
    scat_cloud(i) = atm(i)%scat_cloud
    abs_total(i) = abs_O2(i) + abs_H2O(i) + abs_cloud(i)
  end do
	
  altitude(nlev+1) = 2 * altitude(nlev) - altitude(nlev-1)

! initialize phase matrices
  phase11 = 0.0d0
  phase11_sc = 0.0d0

! loop over layers
  DO ilev = 1,nlev

! for each layr, sum  HG scattering matrix over hydrometeor phases
    do iphase = 1, nphase 

! transfer asymmetry to local variab;e
      g = hydro_prof(ilev,iphase)%cloudg

! calculate indices for linear interpolation
      i_left = min(nhg, max(1, int((nhg*(g + 1.0d0) + 1.0d0)/2.0d0)))
      i_right = min(nhg, (i_left + 1))

! interpolate full HG matrix
      DO j=1,nang
        DO i = 1,nang
          HGgs(1) = HGg(i_left)
          HGgs(2) = HGg(i_right)
          HGphs(1) = HGph(i_left,i,j)
          HGphs(2) = HGph(i_right,i,j)
          dHGphs(1) = dHGph(i_left,i,j)
          dHGphs(2) = dHGph(i_right,i,j)
          HGs = HGphs(1) + ( (HGphs(2)-HGphs(1)) / (HGgs(2)-HGgs(1)) ) * (g-HGgs(1))
          dHGs = dHGphs(1) + ( (dHGphs(2)-dHGphs(1)) / (HGgs(2)-HGgs(1)) ) * (g-HGgs(1))
          phase11(ilev,i,j) = HGs
          phase11_sc(ilev,i,j) = phase11_sc(ilev,i,j) + HGs * hydro_prof(ilev,iphase)%cloudsc
          d_phase11_g(ilev,i,j,iphase) = dHGs * hydro_prof(ilev,iphase)%cloudsc
          d_phase11_sc(ilev,i,j,iphase) =  HGs
        END DO  ! loop over i
      END DO  ! loop over j

! transfer to reduced HG phase matrix
      DO j=1,nangover2
        DO i=1,nangover2
          phasefb(ilev,j,i)= phase11(ilev,j,i)
          phaseff(ilev,j,i)= phase11(ilev,j,nang+1-i)
          phasefb_sc(ilev,j,i)= phase11_sc(ilev,j,i)
          phaseff_sc(ilev,j,i)= phase11_sc(ilev,j,nang+1-i)
          dphasefb_g(ilev,j,i,iphase)  = d_phase11_g(ilev,j,i,iphase)
          dphaseff_g(ilev,j,i,iphase)  = d_phase11_g(ilev,j,nang+1-i,iphase)
          dphasefb_sc(ilev,j,i,iphase) = d_phase11_sc(ilev,j,i,iphase)
          dphaseff_sc(ilev,j,i,iphase) = d_phase11_sc(ilev,j,nang+1-i,iphase)
        END DO
      END DO

    end do ! end iphase loop    
  END DO ! end ilev loop

! transfer surface reflectivity
  DO i=1,nangover2
    IF (ipol == 1) surf_reflec(i) = surf%vref(i)
    IF (ipol == 0) surf_reflec(i) = surf%href(i)
  END DO

! call the core95 setup routine
  CALL dtb94 ()

END subroutine do_tb_gvh94
