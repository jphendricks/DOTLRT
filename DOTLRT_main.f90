!=========================================================================
program dot
!=========================================================================
! driver program for DOTLRT
!
! History:
!  9/16/20  Kevin Schaefer created program
!  9/29/20  Kevin Schaefer added reading WRF output
!  10/5/20  Kevin Schaefer removed all ancient, unused code from gaussq.f90
!  10/5/20  Kevin Schaefer changed surfinp to prevent negative values
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
  integer clock_rate   ! conversion between count and clock time
  integer clock_start  ! time count at stat
  integer clock_stop   ! time count at end
  real e_time          ! elapsed time

! start time
  call system_clock(clock_start) ! start counting

! set up all inputs
  call setup_all_inputs( )

! set up all outputs
  call setup_all_outputs( )
!
! execution time
  call system_clock(clock_stop,clock_rate) ! stop counting
  e_time = real(clock_stop-clock_start)/real(clock_rate)
  print*, 'Execution time Setup (s):', e_time

! start time
  call system_clock(clock_start) ! start counting

! run MRT
  do ichan = chan_strt,chan_stop
    call extract_channel(ichan)
    do ilon = lon_strt,lon_stop
      do ilat = lat_strt,lat_stop

        if (trim(prof_src) == 'WRF') then
          call construct_atmospheric_profile(ilon, ilat)
          call construct_surf_characteristics(ilon, ilat)
        endif

        call mrt( )
        print*, Tbo_mat

        if(save_rad_file) call assign_to_output(ilon, ilat)
        if(save_sing_rad) call write_radiation_profile(ilon, ilat)
      enddo
    enddo

    if(save_rad_file) call create_rad_file(ichan)
  enddo

! execution time
  call system_clock(clock_stop,clock_rate) ! stop counting
  e_time = real(clock_stop-clock_start)/real(clock_rate)
  print*, 'Execution time DOTLRT (s):', e_time

! print message
  print*, 'DOTLRT Run Complete'

end program dot
