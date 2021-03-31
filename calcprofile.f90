! calcprofile.f90
! History:
!     1 - Original from Gail Skofronick Jackson
!     2 - Modified for MRT by Marian Klein, November 1998
!     3 - Provided in PASCAL by Marian Klein and Albin Gasiewski
!         NOAA ETL/ET1 Microwave System Development Division
!     4 - Converted April 2003 from PASCAL to FORTRAN by Ronald Richter
!         ron.richter@noaa.gov NOAA/ETL SET
!         FORTRAN 90 Portland Group Compiler on Red Hat Linux
!     5 - Differentiated and derivatives were tested by Reg Hill Aug. 2003
!     6 - The logical variable "a0_is_constant" was inserted by Reg Hill 8/28/03 to identify cases for
!         which a0 is constant because a great savings of computer time is obtained for that case.
!         The savings is in subroutine HYDROMETEOR_MASTER_5PH_D
!
! uses structure "atm_inp" in "variables_unit.f90"
!
! p and q vary depending on the phase of condensed water, but are constant given for each phase;
!         they are parameters of the gamma distribution of hydrometeor sizes.
! k0 and a0 depend only on water density in each of the phases;
!         they do not depend upon temperature or microwave frequency.

!c subroutine calcprofile_d(a0_is_constant, dk0_dw, da0_dw)
subroutine calcprofile_d()
    use variables
    integer index  
    integer i,j
!c	LOGICAL a0_is_constant(100,5) ! 100 atmospheric levels is arbitrary here; 5 phases
!c	DOUBLE PRECISION:: dk0_dw(100,5), da0_dw(100,5) ! k0 and a0 depend only on w,
	                               ! not temp or freq, and p and q are constants
    DOUBLE PRECISION dk0_dw ! w !w=water density, 
! the size of the arrays, i.e. 100, must be moved to variables.f90

   character*120 debugout

! tests
!   double precision cloud_w_dens(2), cloud_k0(2), cloud_a0(2), dkdw, dadw
! derivatives of k0 and a0 with respect to cloud water density (each phase)
! is multipled by cloud water density (each phase) to remove infinities
! when cloud water density (each phase) goes to zero.
    do index= 1, atm_inp%num_levels
        if( atm_inp%prof(index)%cloud_liq_dens .ne. 0.0d0 ) then
            ! Marshall Palmer distribution
            atm_inp%prof(index)%cloud_liq_p  = 0.0d0
            atm_inp%prof(index)%cloud_liq_q  = 1.0d0
            ! Constant a0, k0 varies w/r/t M
			! dk0_dw(index,1) = 2.0d0 * 1989436788.65d0
            dk0_dw = 2.0d0 * 1989436788.65d0
            atm_inp%prof(index)%dcloud_liq_k0_dw = dk0_dw
       !     atm_inp%prof(index)%dcloud_liq_k0_dw = atm_inp%prof(index)%dcloud_liq_k0_dw &
       !                                          * atm_inp%prof(index)%cloud_liq_dens
            ! atm_inp%prof(index)%cloud_liq_k0 = dk0_dw(index,1)*atm_inp%prof(index)%cloud_liq_dens
            atm_inp%prof(index)%cloud_liq_k0 = dk0_dw * atm_inp%prof(index)%cloud_liq_dens
            atm_inp%prof(index)%cloud_liq_a0 = 0.01d0
			! da0_dw(index,1) = 0.d0
            atm_inp%prof(index)%dcloud_liq_a0_dw = 0.0d0
            a0_is_constant(index,1) = .true.

            !write(debugout,*) "liq:a0,k0=", atm_inp%prof(index)%cloud_liq_a0, atm_inp%prof(index)%cloud_liq_k0
            !call mexPrintf(debugout//achar(10))

            !tests
            !cloud_w_dens(1) = (1.0d0-1.0d-6) * atm_inp%prof(index)%cloud_liq_dens
            !cloud_w_dens(2) = (1.0d0+1.0d-6) * atm_inp%prof(index)%cloud_liq_dens
            !do j = 1, 2
            !  cloud_k0(j) = dk0_dw * cloud_w_dens(j)
            !end do
            !dk_dw = ( cloud_k0(2) - cloud_k0(1) ) / ( cloud_w_dens(2) - cloud_w_dens(1) )
        else
            atm_inp%prof(index)%cloud_liq_p  = 0.0d0
            atm_inp%prof(index)%cloud_liq_q  = 0.0d0
            atm_inp%prof(index)%cloud_liq_k0 = 0.0d0
            atm_inp%prof(index)%cloud_liq_a0 = 0.0d0
            ! da0_dw(index,1) = 0.d0; dk0_dw(index,1) = 0.d0
            atm_inp%prof(index)%dcloud_liq_k0_dw = 0.0d0
            atm_inp%prof(index)%dcloud_liq_a0_dw = 0.0d0
			a0_is_constant(index,1) = .true.
        end if
        
        if( atm_inp%prof(index)%cloud_rn_dens .ne.  0.0d0 ) then
            ! Marshall Palmer distribution
            atm_inp%prof(index)%cloud_rn_p  = 0.0d0
            atm_inp%prof(index)%cloud_rn_q  = 1.0d0
            atm_inp%prof(index)%cloud_rn_k0 = 16000.0d0 
            atm_inp%prof(index)%cloud_rn_a0 = 0.22331d0 &
                   * ( atm_inp%prof(index)%cloud_rn_dens**0.250d0 )
			! da0_dw(index,2) = 0.250d0 * (atm_inp%prof(index)%cloud_rn_a0)/(atm_inp%prof(index)%cloud_rn_dens)
            if( atm_inp%prof(index)%cloud_rn_dens < 1.d0-8 ) then
              atm_inp%prof(index)%dcloud_rn_a0_dw = 0.0d0
            else
              atm_inp%prof(index)%dcloud_rn_a0_dw = 0.250d0 * (atm_inp%prof(index)%cloud_rn_a0)/(atm_inp%prof(index)%cloud_rn_dens)
            end if
      !      atm_inp%prof(index)%dcloud_rn_a0_dw = 0.250d0 * (atm_inp%prof(index)%cloud_rn_a0)
            ! dk0_dw(index,2) = 0.d0
            atm_inp%prof(index)%dcloud_rn_k0_dw = 0.0d0
            a0_is_constant(index,2) = .false.

            !write(debugout,*) "rain:a0,k0=", atm_inp%prof(index)%cloud_rn_a0, atm_inp%prof(index)%cloud_rn_k0
            !call mexPrintf(debugout//achar(10))

            !tests
            !cloud_w_dens(1) = (1.0d0-1.0d-6) * atm_inp%prof(index)%cloud_rn_dens
            !cloud_w_dens(2) = (1.0d0+1.0d-6) * atm_inp%prof(index)%cloud_rn_dens
            !do j = 1, 2
            !  cloud_a0(j) = 0.223d0 * ( cloud_w_dens(j)**0.250d0 )
            !end do
            !da_dw = ( cloud_a0(2) - cloud_a0(1) ) / ( cloud_w_dens(2) - cloud_w_dens(1) )
        else
            atm_inp%prof(index)%cloud_rn_p  = 0.0d0
            atm_inp%prof(index)%cloud_rn_q  = 0.0d0
            atm_inp%prof(index)%cloud_rn_k0 = 0.0d0
            atm_inp%prof(index)%cloud_rn_a0 = 0.0d0
			! da0_dw(index,2) = 0.d0  ;  dk0_dw (index,2) = 0.d0
            atm_inp%prof(index)%dcloud_rn_k0_dw = 0.0d0
            atm_inp%prof(index)%dcloud_rn_a0_dw = 0.0d0
			a0_is_constant(index,2) = .true.
        end if
        
        if( atm_inp%prof(index)%cloud_ice_dens .ne. 0.0d0 ) then
            ! Sekhon-Sriv. Distribution, for two phase ice
            atm_inp%prof(index)%cloud_ice_p  = 0.0d0
            atm_inp%prof(index)%cloud_ice_q  = 1.0d0
            ! Constant a0, k0 varies w/r/t M, assumes no size change when changing phase
            ! (see Bauer and Schluessel) - allows for polydispersive sizes and for mie
            ! and Rayleigh scattering
			! dk0_dw(index,3) = 2.0d0 * 1989436788.65d0
            dk0_dw = 2.0d0 * 1989436788.65d0
            atm_inp%prof(index)%dcloud_ice_k0_dw = dk0_dw
     !       atm_inp%prof(index)%dcloud_ice_k0_dw = atm_inp%prof(index)%dcloud_ice_k0_dw &
     !                                            * atm_inp%prof(index)%cloud_ice_dens
            ! atm_inp%prof(index)%cloud_ice_k0 = dk0_dw (index,3) * atm_inp%prof(index)%cloud_ice_dens 
            atm_inp%prof(index)%cloud_ice_k0 = dk0_dw * atm_inp%prof(index)%cloud_ice_dens
            atm_inp%prof(index)%cloud_ice_a0 =  0.01d0
			! da0_dw(index,3) = 0.d0
            atm_inp%prof(index)%dcloud_ice_a0_dw = 0.d0
			a0_is_constant(index,3) = .true.

            !tests
            !cloud_w_dens(1) = (1.0d0-1.0d-6) * atm_inp%prof(index)%cloud_ice_dens
            !cloud_w_dens(2) = (1.0d0+1.0d-6) * atm_inp%prof(index)%cloud_ice_dens
            !do j = 1, 2
            !  cloud_k0(j) = dk0_dw * cloud_w_dens(j)
            !end do
            !dk_dw = ( cloud_k0(2) - cloud_k0(1) ) / ( cloud_w_dens(2) - cloud_w_dens(1) )
        else
            atm_inp%prof(index)%cloud_ice_p  = 0.0d0
            atm_inp%prof(index)%cloud_ice_q  = 0.0d0
            atm_inp%prof(index)%cloud_ice_k0 = 0.0d0
            atm_inp%prof(index)%cloud_ice_a0 = 0.0d0
			! da0_dw(index,3) = 0.d0  ;  dk0_dw(index,3) = 0.d0
            atm_inp%prof(index)%dcloud_ice_k0_dw = 0.0d0
            atm_inp%prof(index)%dcloud_ice_a0_dw = 0.0d0
			a0_is_constant(index,3) = .true.
        end if
        
        if( atm_inp%prof(index)%cloud_snow_dens .ne. 0.0d0 ) then
            ! Tao, Prasad, Alder snow size distributions converted to a0 and k0 form
            ! (see ntbk#3 pg 111)
            atm_inp%prof(index)%cloud_snow_p  = 0.0d0
            atm_inp%prof(index)%cloud_snow_q  = 1.0d0
            ! If following Rutledge and Hobbs (1983), then parameterization:
            !                                         cloudsnowk0 = 20000.0
            ! Instead, use Rutledge and Hobbs (1984)
            atm_inp%prof(index)%cloud_snow_k0 = 8000.0d0
		    ! dk0_dw(index,4) = 0.d0
            atm_inp%prof(index)%dcloud_snow_k0_dw = 0.d0
            atm_inp%prof(index)%cloud_snow_a0 = 0.3757d0 &
                    * ( atm_inp%prof(index)%cloud_snow_dens**0.25d0 )
            ! da0_dw(index,4) = 0.25d0* (atm_inp%prof(index)%cloud_snow_a0)/(atm_inp%prof(index)%cloud_snow_dens)
            atm_inp%prof(index)%dcloud_snow_a0_dw=0.25d0*(atm_inp%prof(index)%cloud_snow_a0)/(atm_inp%prof(index)%cloud_snow_dens)
         !   atm_inp%prof(index)%dcloud_snow_a0_dw = 0.25d0* (atm_inp%prof(index)%cloud_snow_a0)
            a0_is_constant(index,4) = .false.

            !tests
            !cloud_w_dens(1) = (1.0d0-1.0d-6) * atm_inp%prof(index)%cloud_snow_dens
            !cloud_w_dens(2) = (1.0d0+1.0d-6) * atm_inp%prof(index)%cloud_snow_dens
            !do j = 1, 2
            !  cloud_a0(j) = 0.3757d0 * cloud_w_dens(j)**0.25d0
            !end do
            !da_dw = ( cloud_a0(2) - cloud_a0(1) ) / ( cloud_w_dens(2) - cloud_w_dens(1) )
        else
            atm_inp%prof(index)%cloud_snow_p  = 0.0d0
            atm_inp%prof(index)%cloud_snow_q  = 0.0d0
            atm_inp%prof(index)%cloud_snow_k0 = 0.0d0
            atm_inp%prof(index)%cloud_snow_a0 = 0.0d0
            ! da0_dw(index,4) = 0.d0  ;  dk0_dw(index,4) = 0.d0
            atm_inp%prof(index)%dcloud_snow_k0_dw = 0.0d0
            atm_inp%prof(index)%dcloud_snow_a0_dw = 0.0d0
	        a0_is_constant(index,4) = .true.
        end if
        
        if( atm_inp%prof(index)%cloud_grpl_dens .ne. 0.0d0 ) then
            ! Tao, Prasad, Alder graupel size distributions converted to a0 and k0 form
            ! (see ntbk#3 pg 111)
            atm_inp%prof(index)%cloud_grpl_p  = 0.0d0
            atm_inp%prof(index)%cloud_grpl_q  = 1.0d0
            atm_inp%prof(index)%cloud_grpl_k0 = 8000.0d0 
		    ! dk0_dw(index,5) = 0.d0
            atm_inp%prof(index)%dcloud_grpl_k0_dw = 0.0d0
            atm_inp%prof(index)%cloud_grpl_a0 = 0.3340d0 &
                    * ( atm_inp%prof(index)%cloud_grpl_dens**0.25d0 )
		    ! da0_dw(index,5) = 0.25d0 *(atm_inp%prof(index)%cloud_grpl_a0)/(atm_inp%prof(index)%cloud_grpl_dens)
            atm_inp%prof(index)%dcloud_grpl_a0_dw=0.25d0*(atm_inp%prof(index)%cloud_grpl_a0)/(atm_inp%prof(index)%cloud_grpl_dens)
!            atm_inp%prof(index)%dcloud_grpl_a0_dw = 0.25d0 *(atm_inp%prof(index)%cloud_grpl_a0)
			a0_is_constant(index,5) = .false.

            !write(debugout,*) "grp:a0,k0=", atm_inp%prof(index)%cloud_grpl_a0, atm_inp%prof(index)%cloud_grpl_k0
            !call mexPrintf(debugout//achar(10))

            !tests
            !cloud_w_dens(1) = (1.0d0-1.0d-6) * atm_inp%prof(index)%cloud_grpl_dens
            !cloud_w_dens(2) = (1.0d0+1.0d-6) * atm_inp%prof(index)%cloud_grpl_dens
            !do j = 1, 2
            !  cloud_a0(j) = 0.3340d0 * ( cloud_w_dens(j)**0.25d0 )
            !end do
            !da_dw = ( cloud_a0(2) - cloud_a0(1) ) / ( cloud_w_dens(2) - cloud_w_dens(1) ) 
        else
            atm_inp%prof(index)%cloud_grpl_p  = 0.0d0
            atm_inp%prof(index)%cloud_grpl_q  = 0.0d0
            atm_inp%prof(index)%cloud_grpl_k0 = 0.0d0
            atm_inp%prof(index)%cloud_grpl_a0 = 0.0d0
	        ! da0_dw(index,5) = 0.d0 ; dk0_dw(index,5) = 0.d0
            atm_inp%prof(index)%dcloud_grpl_k0_dw = 0.0d0
            atm_inp%prof(index)%dcloud_grpl_a0_dw = 0.0d0
	        a0_is_constant(index,5) = .true.
        end if
    end do 
    return
end subroutine calcprofile_d
