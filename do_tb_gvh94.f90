!===================================================================================
SUBROUTINE do_tb_gvh94()
!===================================================================================
!  This program reads input file and calls subroutines that further prepare the data for the
!  calculation of brightness temperature and incremental profiles by CORE95.F90. 
!  This subroutine reads input file and calls subroutines that further prepare the data for the
!  calculation of brightness temperature and incremental profiles by CORE95.F90.
!  Brightness temperature and incremental profiles are determined
!  as functions of altitude. This subroutine also allows to insert pseudo-layers
!  for more detailed account for z-dependency
!------------------------------------------------------------------------------------
!                        Modification of 04/24/03
!   This modification introdices new parameter _eps_err (max relative error in brightness
!   temperature allowed) which will be used to determine whether a particular matrix can
!   be considered diagonal or not.
!   It also calls subroutine DTB94.F90 which uses parameter _eps_err (the only use of this
!   parameter by DO_TB_GVH94.F90 is to pass this  parameter to DTB94.F90).
!
!------------------------------------------------------------------------------------
!                        Modification of 12/16/02
!
!   This is a modification of subroutine DO_TB_GVH9.F90 which uses pre-calculated
!   values of HG phase matrix along with its derivative for parameter _g running 
!   from (1-ng)/ng,(3-ng)/ng,...,(ng-1)/ng   [ _ng values of parameter _g altogether] 
!   and for angles _quad_angle_array which are extended by inclusion angles 0 and 180 deg.
!   This extension is consistent with the extension which is done in DO_TB_GVH93 below.
!
!   There are three new input parameters: 
!
!                             nHG   - ( I)
!
!   sc_mat(nHG, nang_I+2, nang_I+2) - (DP)
!
!  dsc_mat(nHG, nang_I+2, nang_I+2) - (DP)
!
!   Scattering matrix _sc_mat should be calculated before calling this subroutine
!   with the help of subroutine HG_phmat:
!
!   CALL HG_phmat (nang_I+2, qa_array, nHG, HGg, sc_mat, dsc_mat)

!   Here: _qa_array is coincident with the array which bears the same name
!          inside DO_TB_GVH91. It is as follows: 
!
!            qa_array(1       )=0; 
!            qa_array(k       )=quad_angle_array(k-1), k=2,..,nang_I+1;
!            qa_array(nang_I+2)=180.
!
!            HGg(k)= (2*k-nHG-1.)/nHG , k=1,2,...,nHG
!
!  Parameter _nHG is an extra input parameter which should be large enough
!  (say, nHG=101).
!
!  Parameter _nvar should be set to 4 rather then 3: _nvar=4. Variations with respect
!  to HG parameter _g will be also included into the output.
!------------------------------------------------------------------------------------
!  9/26/2020 Kevin Schaefer deleted unused variables
!  10/15/20 Kevin Schaefer deleted commented code
!------------------------------------------------------------------------------------
use dotlrt_variables
implicit none

    INTEGER iphase, ilr, i, j, iHG
    real(8) g
    real(8) HGgs(2), HGphs(2), dHGphs(2), HGs, dHGs
    integer jHGphs(2)
    
    EXTERNAL    dtb94   ! InterpSR

    do i=1, nlev
      g_asymmetry(i) = atm(i)%asymmetry 
         altitude(i) = atm(i)%hgt
      temperature(i) = atm(i)%temperature
           abs_O2(i) = atm(i)%abs_O2
          abs_H2O(i) = atm(i)%abs_H2O
        abs_cloud(i) = atm(i)%abs_cloud
       scat_cloud(i) = atm(i)%scat_cloud
       abs_total(i) = abs_O2(i) + abs_H2O(i) + abs_cloud(i)
    end do
	
    altitude(nlev+1) = 2 * altitude(nlev) &
                                       - altitude(nlev-1)


    ! CALCULATE REDUCED PHASE MATRIX
    phase11 = 0.0d0
    phase11_sc = 0.0d0
    DO ilr = 1,nlev
      do iphase = 1, nphase ! sum over hydrometeor phases HG scattering matrix
        g = hydro_prof(ilr,iphase)%cloudg
        !iHG = min( ng, max( 1, INT(ng*(g+1)/2+1) ) ) !? INT(ng*(g+1)/2)
        iHG = min( ng, max( 1, int(( ng * ( g + 1.0d0 ) + 1.0d0 ) / 2.0d0 ) ) )
        jHGphs(1) = max(  1, (iHG - 1) )
        jHGphs(2) = min( ng, (iHG + 1) )
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
            phase11(ilr,i,j) = HGs
            phase11_sc(ilr,i,j) = phase11_sc(ilr,i,j) + HGs * hydro_prof(ilr,iphase)%cloudsc
            d_phase11_g(ilr,i,j,iphase) = dHGs * hydro_prof(ilr,iphase)%cloudsc
            d_phase11_sc(ilr,i,j,iphase) =  HGs
          END DO  ! loop over i
        END DO  ! loop over j
        DO j=1,nangover2
          DO i=1,nangover2
            phasefb(ilr,j,i)= phase11(ilr,j,i)
            phaseff(ilr,j,i)= phase11(ilr,j,nang+1-i)
            phasefb_sc(ilr,j,i)= phase11_sc(ilr,j,i)
            phaseff_sc(ilr,j,i)= phase11_sc(ilr,j,nang+1-i)
            dphasefb_g(ilr,j,i,iphase)  = d_phase11_g(ilr,j,i,iphase)
            dphaseff_g(ilr,j,i,iphase)  = d_phase11_g(ilr,j,nang+1-i,iphase)
            dphasefb_sc(ilr,j,i,iphase) = d_phase11_sc(ilr,j,i,iphase)
            dphaseff_sc(ilr,j,i,iphase) = d_phase11_sc(ilr,j,nang+1-i,iphase)
          END DO
        END DO
     end do ! iphase     
    END DO ! ilr


    DO i=1,nangover2
      IF (ipol == 1) surf_reflec(i) = surf_inp%vr(i)
      IF (ipol == 0) surf_reflec(i) = surf_inp%hr(i)
    END DO

    CALL dtb94 ()

END subroutine do_tb_gvh94
