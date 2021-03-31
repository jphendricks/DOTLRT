
!  This program reads input file and calls subroutines that further prepare the data for the
!  calculation of brightness temperature and
!  incremental profiles by CORE95.F90, details are below the list of modifications. 
!------------------------------------------------------------------------------------
!                        Modification of 04/24/03
!
!   This modification introdices new parameter _eps_err (max relative error in brightness
!   temperature allowed) which will be used to determine whether a particular matrix can
!   be considered diagonal or not.
!
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
!    qa_array = teta

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
!   
!------------------------------------------------------------------------------------

!-----------------------------------------------------
!
!  This subroutine reads input file and calls subroutines that further prepare the data for the
!  calculation of brightness temperature and incremental profiles by CORE95.F90.
!  Brightness temperature and incremental profiles are determined
!  as functions of altitude. This subroutine also allows to insert pseudo-layers
!  for more detailed account for z-dependency
!
!-----------------------------------------------------
SUBROUTINE do_tb_gvh94()
use variables
!    use type_kinds
!    use readprofile
    implicit none

    INTEGER hydrometeor_phase, ilr, i, j, iHG
    real(8) g
    real(8) HGgs(2), HGphs(2), dHGphs(2), HGs, dHGs
    integer jHGphs(2)
    EXTERNAL    dtb94   ! InterpSR

!    do i = 1, surf_inp%num_angles
!      reflec_angle_array(i) = surf_inp%theta(i)
!      surfinp_reflecv(i) = surf_inp%vr(i)
!      surfinp_reflech(i) = surf_inp%hr(i)
!    end do

    do i=1, atm_inp%num_levels
      g_asymmetry(i) = atm_inp%prof(i)%asymmetry 
         altitude(i) = atm_inp%prof(i)%height
      temperature(i) = atm_inp%prof(i)%temperature
           abs_O2(i) = atm_inp%prof(i)%abs_O2
          abs_H2O(i) = atm_inp%prof(i)%abs_H2O
        abs_cloud(i) = atm_inp%prof(i)%abs_cloud
       scat_cloud(i) = atm_inp%prof(i)%scat_cloud
        abs_total(i) = abs_O2(i) + abs_H2O(i) + abs_cloud(i)
    end do
	altitude(atm_inp%num_levels+1) = 2 * altitude(atm_inp%num_levels) &
                                       - altitude(atm_inp%num_levels-1)

!    !retain reflectivity at zero degrees as follows
!	surf_reflecv(1) = surfinp_reflecv(1)
!	surf_reflech(1) = surfinp_reflech(1) 

!    !retain reflectivity at 90 degrees as follows
!	surf_reflecv(nangover2+1) = surfinp_reflecv(surf_inp%num_angles) 
!	surf_reflech(nangover2+1) = surfinp_reflech(surf_inp%num_angles)

!	CALL InterpSR(nangover2-1, surf_inp%num_angles, teta(2), reflec_angle_array, &
!                   surfinp_reflecv, surfinp_reflech, surf_reflecv, surf_reflech)

    ! CALCULATE REDUCED PHASE MATRIX
    phase11 = 0.0d0
    phase11_sc = 0.0d0
    DO ilr = 1,atm_inp%num_levels
      do hydrometeor_phase = 1, number_h2o_phases ! sum over hydrometeor phases HG scattering matrix
        g = hydro_prof(ilr,hydrometeor_phase)%cloudg
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
            phase11_sc(ilr,i,j) = phase11_sc(ilr,i,j) + HGs * hydro_prof(ilr,hydrometeor_phase)%cloudsc
            d_phase11_g(ilr,i,j,hydrometeor_phase) = dHGs * hydro_prof(ilr,hydrometeor_phase)%cloudsc
            d_phase11_sc(ilr,i,j,hydrometeor_phase) =  HGs
          END DO  ! loop over i
        END DO  ! loop over j
        DO j=1,nangover2
          DO i=1,nangover2
            phasefb(ilr,j,i)= phase11(ilr,j,i)
            phaseff(ilr,j,i)= phase11(ilr,j,nang+1-i)
            phasefb_sc(ilr,j,i)= phase11_sc(ilr,j,i)
            phaseff_sc(ilr,j,i)= phase11_sc(ilr,j,nang+1-i)
            dphasefb_g(ilr,j,i,hydrometeor_phase)  = d_phase11_g(ilr,j,i,hydrometeor_phase)
            dphaseff_g(ilr,j,i,hydrometeor_phase)  = d_phase11_g(ilr,j,nang+1-i,hydrometeor_phase)
            dphasefb_sc(ilr,j,i,hydrometeor_phase) = d_phase11_sc(ilr,j,i,hydrometeor_phase)
            dphaseff_sc(ilr,j,i,hydrometeor_phase) = d_phase11_sc(ilr,j,nang+1-i,hydrometeor_phase)
          END DO
        END DO
      end do ! hydrometeor_phase
    END DO ! ilr

!    DO i=1,nang
!         cs(i)= DCOS(teta(i)/180.d0*pi)
!    END DO

    DO i=1,nangover2
      ! inp_pol (1:vertical; 0:horizontal)
      !IF (inp_pol == 1) surf_reflec(i) = surf_reflecv(i)
      !IF (inp_pol == 0) surf_reflec(i) = surf_reflech(i)
      IF (inp_pol == 1) surf_reflec(i) = surf_inp%vr(i)
      IF (inp_pol == 0) surf_reflec(i) = surf_inp%hr(i)
    END DO

    CALL dtb94 ()

END subroutine do_tb_gvh94
