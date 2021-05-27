!===================================================================================
SUBROUTINE do_tb_gvh94()
!===================================================================================
!  This program reads input file and calls subroutines that further prepare the data for the
!  calculation of brightness temperature and incremental profiles by CORE95.F90. 
!  This subroutine reads input file and calls subroutines that further prepare the data for the
!  calculation of brightness temperature and incremental profiles by CORE95.F90.
!  Brightness temperature and incremental profiles are determined
!
! History:
!  4/24/2003  added eps_err max relative error in brightness temperature allowed
!  12/16/2002 pre-calculated HG phase matrix 
!  9/26/2020  Kevin Schaefer deleted unused variables
!  10/15/2020 Kevin Schaefer deleted commented code
!  12/12/2020 Kevin Schaefer removed external dtb94 statement
!  12/12/2020 Kevin schaefer removed all outdated comments and variable defs
!------------------------------------------------------------------------------------
use dotlrt_variables
implicit none

    INTEGER iphase, i, j, iHG
    integer ilev ! (-) layer or level index
    real(8) g
    real(8) HGgs(2), HGphs(2), dHGphs(2), HGs, dHGs
    integer jHGphs(2)

    do i=1, nlev
      g_asymmetry(i) = atm(i)%asymmetry 
         altitude(i) = atm(i)%hgt
      temperature(i) = atm(i)%temp
           abs_O2(i) = atm(i)%abs_O2
          abs_H2O(i) = atm(i)%abs_H2O
        abs_cloud(i) = atm(i)%abs_cloud
       scat_cloud(i) = atm(i)%scat_cloud
       abs_total(i) = abs_O2(i) + abs_H2O(i) + abs_cloud(i)
    end do
	
    altitude(nlev+1) = 2 * altitude(nlev) - altitude(nlev-1)

    ! CALCULATE REDUCED PHASE MATRIX
    phase11 = 0.0d0
    phase11_sc = 0.0d0
    DO ilev = 1,nlev
      do iphase = 1, nphase ! sum over hydrometeor phases HG scattering matrix
        g = hydro_prof(ilev,iphase)%cloudg
        iHG = min( nhg, max( 1, int(( nhg * ( g + 1.0d0 ) + 1.0d0 ) / 2.0d0 ) ) )
        jHGphs(1) = max(  1, (iHG - 1) )
        jHGphs(2) = min( nhg, (iHG + 1) )
        DO j=1,nang
          DO i = 1,nang
            HGgs(1) = HGg(jHGphs(1))
            HGgs(2) = HGg(jHGphs(2))
            HGphs(1) = HGph(jHGphs(1),i,j)
            HGphs(2) = HGph(jHGphs(2),i,j)
            dHGphs(1) = dHGph(jHGphs(1),i,j)
            dHGphs(2) = dHGph(jHGphs(2),i,j)
            HGs = HGphs(1) + ( (HGphs(2)-HGphs(1)) / (HGgs(2)-HGgs(1)) ) * (g-HGgs(1))
            dHGs = dHGphs(1) + ( (dHGphs(2)-dHGphs(1)) / (HGgs(2)-HGgs(1)) ) * (g-HGgs(1))
            phase11(ilev,i,j) = HGs
            phase11_sc(ilev,i,j) = phase11_sc(ilev,i,j) + HGs * hydro_prof(ilev,iphase)%cloudsc
            d_phase11_g(ilev,i,j,iphase) = dHGs * hydro_prof(ilev,iphase)%cloudsc
            d_phase11_sc(ilev,i,j,iphase) =  HGs
          END DO  ! loop over i
        END DO  ! loop over j
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
     end do ! iphase     
    END DO ! ilev

    DO i=1,nangover2
      IF (int(ipol) == 1) surf_reflec(i) = surf%vref(i)
      IF (int(ipol) == 0) surf_reflec(i) = surf%href(i)
    END DO

    CALL dtb94 ()

END subroutine do_tb_gvh94
