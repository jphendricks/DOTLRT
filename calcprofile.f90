!===================================================================
subroutine calcprofile_d()
!===================================================================
! p and q vary depending on the phase of condensed water, but are constant given for each phase;
!         they are parameters of the gamma distribution of hydrometeor sizes.
! k0 and a0 depend only on water density in each of the phases;
!         they do not depend upon temperature or microwave frequency.
! History:
!   9/26/2020 Kevin Schaefer deleted unused variables and tabs 
!     1 - Original from Gail Skofronick Jackson
!     2 - Modified for MRT by Marian Klein, November 1998
!     3 - Provided in PASCAL by Marian Klein and Albin Gasiewski
!         NOAA ETL/ET1 Microwave System Development Division
!     4 - Converted April 2003 from PASCAL to FORTRAN by Ronald Richter
!         ron.richter@noaa.gov NOAA/ETL SET
!         FORTRAN 90 Portland Group Compiler on Red Hat Linux
!     5 - Differentiated and derivatives were tested by Reg Hill Aug. 2003
!     6 - The logical variable "a0_const" was inserted by Reg Hill 8/28/03 to identify cases for
!         which a0 is constant because a great savings of computer time is obtained for that case.
!         The savings is in subroutine HYDROMETEOR_MASTER_5PH_D
!-------------------------------------------------------------------
    use dotlrt_variables
    integer index  
    DOUBLE PRECISION dk0_dw ! w !w=water density, 

! tests
!   double precision cloud_w_dens(2), cloud_k0(2), cloud_a0(2), dkdw, dadw
! derivatives of k0 and a0 with respect to cloud water density (each phase)
! is multipled by cloud water density (each phase) to remove infinities
! when cloud water density (each phase) goes to zero.
    do index= 1, nlev
        if( atm(index)%clw_dens .ne. 0.0d0 ) then
            ! Marshall Palmer distribution
            atm(index)%clw_p  = 0.0d0
            atm(index)%clw_q  = 1.0d0
            ! Constant a0, k0 varies w/r/t M
            dk0_dw = 2.0d0 * 1989436788.65d0
            atm(index)%dclw_k0_dw = dk0_dw

            atm(index)%clw_k0 = dk0_dw * atm(index)%clw_dens
            atm(index)%clw_a0 = 0.01d0
            atm(index)%dclw_a0_dw = 0.0d0
            a0_const(index,1) = .true.
        else
            atm(index)%clw_p  = 0.0d0
            atm(index)%clw_q  = 0.0d0
            atm(index)%clw_k0 = 0.0d0
            atm(index)%clw_a0 = 0.0d0
            atm(index)%dclw_k0_dw = 0.0d0
            atm(index)%dclw_a0_dw = 0.0d0
            a0_const(index,1) = .true.
        end if
        
        if( atm(index)%rain_dens .ne.  0.0d0 ) then
            ! Marshall Palmer distribution
            atm(index)%rain_p  = 0.0d0
            atm(index)%rain_q  = 1.0d0
            atm(index)%rain_k0 = 16000.0d0 
            atm(index)%rain_a0 = 0.22331d0 &
                   * ( atm(index)%rain_dens**0.250d0 )
         ! da0_dw(index,2) = 0.250d0 * (atm(index)%rain_a0)/(atm(index)%rain_dens)
            if( atm(index)%rain_dens < 1.d0-8 ) then
              atm(index)%drain_a0_dw = 0.0d0
            else
              atm(index)%drain_a0_dw = 0.250d0 * (atm(index)%rain_a0)/(atm(index)%rain_dens)
            end if
      !      atm(index)%drain_a0_dw = 0.250d0 * (atm(index)%rain_a0)
            ! dk0_dw(index,2) = 0.d0
            atm(index)%drain_k0_dw = 0.0d0
            a0_const(index,2) = .false.

            !write(debugout,*) "rain:a0,k0=", atm(index)%rain_a0, atm(index)%rain_k0
            !all mexPrintf(debugout//achar(10))

            !tests
            !loud_w_dens(1) = (1.0d0-1.0d-6) * atm(index)%rain_dens
            !loud_w_dens(2) = (1.0d0+1.0d-6) * atm(index)%rain_dens
            !do j = 1, 2
            !  cloud_a0(j) = 0.223d0 * ( cloud_w_dens(j)**0.250d0 )
            !end do
            !da_dw = ( cloud_a0(2) - cloud_a0(1) ) / ( cloud_w_dens(2) - cloud_w_dens(1) )
        else
            atm(index)%rain_p  = 0.0d0
            atm(index)%rain_q  = 0.0d0
            atm(index)%rain_k0 = 0.0d0
            atm(index)%rain_a0 = 0.0d0
			! da0_dw(index,2) = 0.d0  ;  dk0_dw (index,2) = 0.d0
            atm(index)%drain_k0_dw = 0.0d0
            atm(index)%drain_a0_dw = 0.0d0
            a0_const(index,2) = .true.
        end if
        
        if( atm(index)%ice_dens .ne. 0.0d0 ) then
            ! Sekhon-Sriv. Distribution, for two phase ice
            atm(index)%ice_p  = 0.0d0
            atm(index)%ice_q  = 1.0d0
            ! Constant a0, k0 varies w/r/t M, assumes no size change when changing phase
            ! (see Bauer and Schluessel) - allows for polydispersive sizes and for mie
            ! and Rayleigh scattering
         ! dk0_dw(index,3) = 2.0d0 * 1989436788.65d0
            dk0_dw = 2.0d0 * 1989436788.65d0
            atm(index)%dice_k0_dw = dk0_dw
     !       atm(index)%dice_k0_dw = atm(index)%dice_k0_dw &
     !                                            * atm(index)%ice_dens
            ! atm(index)%ice_k0 = dk0_dw (index,3) * atm(index)%ice_dens 
            atm(index)%ice_k0 = dk0_dw * atm(index)%ice_dens
            atm(index)%ice_a0 =  0.01d0
         ! da0_dw(index,3) = 0.d0
            atm(index)%dice_a0_dw = 0.d0
         a0_const(index,3) = .true.

            !tests
            !loud_w_dens(1) = (1.0d0-1.0d-6) * atm(index)%ice_dens
            !loud_w_dens(2) = (1.0d0+1.0d-6) * atm(index)%ice_dens
            !do j = 1, 2
            !  cloud_k0(j) = dk0_dw * cloud_w_dens(j)
            !end do
            !dk_dw = ( cloud_k0(2) - cloud_k0(1) ) / ( cloud_w_dens(2) - cloud_w_dens(1) )
        else
            atm(index)%ice_p  = 0.0d0
            atm(index)%ice_q  = 0.0d0
            atm(index)%ice_k0 = 0.0d0
            atm(index)%ice_a0 = 0.0d0
         ! da0_dw(index,3) = 0.d0  ;  dk0_dw(index,3) = 0.d0
            atm(index)%dice_k0_dw = 0.0d0
            atm(index)%dice_a0_dw = 0.0d0
            a0_const(index,3) = .true.
        end if
        
        if( atm(index)%snow_dens .ne. 0.0d0 ) then
            ! Tao, Prasad, Alder snow size distributions converted to a0 and k0 form
            ! (see ntbk#3 pg 111)
            atm(index)%snow_p  = 0.0d0
            atm(index)%snow_q  = 1.0d0
            ! If following Rutledge and Hobbs (1983), then parameterization:
            !                                         cloudsnowk0 = 20000.0
            ! Instead, use Rutledge and Hobbs (1984)
            atm(index)%snow_k0 = 8000.0d0
          ! dk0_dw(index,4) = 0.d0
            atm(index)%dsnow_k0_dw = 0.d0
            atm(index)%snow_a0 = 0.3757d0 &
                    * ( atm(index)%snow_dens**0.25d0 )
            ! da0_dw(index,4) = 0.25d0* (atm(index)%snow_a0)/(atm(index)%snow_dens)
            atm(index)%dsnow_a0_dw=0.25d0*(atm(index)%snow_a0)/(atm(index)%snow_dens)
         !   atm(index)%dsnow_a0_dw = 0.25d0* (atm(index)%snow_a0)
            a0_const(index,4) = .false.

            !tests
            !loud_w_dens(1) = (1.0d0-1.0d-6) * atm(index)%snow_dens
            !loud_w_dens(2) = (1.0d0+1.0d-6) * atm(index)%snow_dens
            !do j = 1, 2
            !  cloud_a0(j) = 0.3757d0 * cloud_w_dens(j)**0.25d0
            !end do
            !da_dw = ( cloud_a0(2) - cloud_a0(1) ) / ( cloud_w_dens(2) - cloud_w_dens(1) )
        else
            atm(index)%snow_p  = 0.0d0
            atm(index)%snow_q  = 0.0d0
            atm(index)%snow_k0 = 0.0d0
            atm(index)%snow_a0 = 0.0d0
            ! da0_dw(index,4) = 0.d0  ;  dk0_dw(index,4) = 0.d0
            atm(index)%dsnow_k0_dw = 0.0d0
            atm(index)%dsnow_a0_dw = 0.0d0
           a0_const(index,4) = .true.
        end if
        
        if( atm(index)%grpl_dens .ne. 0.0d0 ) then
            ! Tao, Prasad, Alder graupel size distributions converted to a0 and k0 form
            ! (see ntbk#3 pg 111)
            atm(index)%grpl_p  = 0.0d0
            atm(index)%grpl_q  = 1.0d0
            atm(index)%grpl_k0 = 8000.0d0 
          ! dk0_dw(index,5) = 0.d0
            atm(index)%dgrpl_k0_dw = 0.0d0
            atm(index)%grpl_a0 = 0.3340d0 &
                    * ( atm(index)%grpl_dens**0.25d0 )
          ! da0_dw(index,5) = 0.25d0 *(atm(index)%grpl_a0)/(atm(index)%grpl_dens)
            atm(index)%dgrpl_a0_dw=0.25d0*(atm(index)%grpl_a0)/(atm(index)%grpl_dens)
!            atm(index)%dgrpl_a0_dw = 0.25d0 *(atm(index)%grpl_a0)
         a0_const(index,5) = .false.

            !write(debugout,*) "grp:a0,k0=", atm(index)%grpl_a0, atm(index)%grpl_k0
            !all mexPrintf(debugout//achar(10))

            !tests
            !loud_w_dens(1) = (1.0d0-1.0d-6) * atm(index)%grpl_dens
            !loud_w_dens(2) = (1.0d0+1.0d-6) * atm(index)%grpl_dens
            !do j = 1, 2
            !  cloud_a0(j) = 0.3340d0 * ( cloud_w_dens(j)**0.25d0 )
            !end do
            !da_dw = ( cloud_a0(2) - cloud_a0(1) ) / ( cloud_w_dens(2) - cloud_w_dens(1) ) 
        else
            atm(index)%grpl_p  = 0.0d0
            atm(index)%grpl_q  = 0.0d0
            atm(index)%grpl_k0 = 0.0d0
            atm(index)%grpl_a0 = 0.0d0
           ! da0_dw(index,5) = 0.d0 ; dk0_dw(index,5) = 0.d0
            atm(index)%dgrpl_k0_dw = 0.0d0
            atm(index)%dgrpl_a0_dw = 0.0d0
           a0_const(index,5) = .true.
        end if
    end do 
    return
end subroutine calcprofile_d
