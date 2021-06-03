!=========================================================================
program dot
!=========================================================================
! driver program for DOTLRT
!
! History:
!  9/16/2020  Kevin Schaefer created program
!  9/29/2020  Kevin Schaefer added reading WRF output
!  10/5/2020  Kevin Schaefer removed all ancient, unused code from gaussq.f90
!  10/5/2020  Kevin Schaefer changed surfinp to prevent negative values
!  11/30/2020 Kevin Schaefer deleted call to surface characteristics in main loop
!  11/30/2020 Kevin Schaefer recalculate surface reflectances each time frequency changes
!-------------------------------------------------------------------------
!
  use dotlrt_variables
  use profiles
  use dotlrt_output
!
  implicit none

! internal variables
  integer ichan        ! (-) channel index
  integer ilon         ! (-) longitude index
  integer ilat         ! (-) latitude index

! variables for timing
  integer trate  ! conversion between count and clock time
  integer tstart ! time count at start
  integer tend   ! time count at end
  real tdelta    ! elapsed time

! collect command line input arguments
call get_command_argument(1, file_in)
call get_command_argument(2, out_path)

! start time
  call system_clock(tstart) ! start counting

! set up all inputs
  call setup_all_inputs( )

! set up all outputs
  call setup_all_outputs( )
!
! execution time
  call system_clock(tend,trate) ! stop counting
  tdelta = real(tend-tstart)/real(trate)
  print*, 'Execution time Setup (s):', tdelta

! start time
  call system_clock(tstart) ! start counting

! run MRT
  do ichan = chan_strt,chan_stop
    call extract_channel(ichan)
    if(flag_print_full) print*, 'Channel: ',ichan, channel%lo_freq

    do ilon = lon_strt,lon_stop
      do ilat = lat_strt,lat_stop
        if (trim(prof_src) == 'WRF') call construct_atmospheric_profile(ichan, ilon, ilat)
        if (trim(prof_src) == 'single') call construct_single_surf_ref()

        call mrt( )
        print*, Tbo_mat

        if(save_rad_file) call assign_to_output(ilon, ilat)
        if(save_sing_rad) call write_radiation_profile(ilon, ilat)
      enddo
    enddo
    if(save_rad_file) call write_tb_channel(ichan)
    if(save_jac_file) call create_jacobian_file(ichan)
  enddo

! execution time
  call system_clock(tend,trate) ! stop counting
  tdelta = real(tend-tstart)/real(trate)
  print*, 'Execution time DOTLRT (s):', tdelta

! print message
  print*, 'DOTLRT Run Complete'

end program dot
