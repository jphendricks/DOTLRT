!============================================================================
subroutine get_instr_spec( )
!============================================================================
! This routine:
!     1 - reads instrument specification file containing parameters:
!         a - number of channels
!         b - frequencies
!         c - bandwidths
!         d - noise parameters
!     2 - computes instrument noise covariance matrix assuming
!         channels are non-overlapping.
!
! history:
!  10/17/2020 Kevin schaefer switched to reading an instrument file
!----------------------------------------------------------------------------
  use dotlrt_variables
  implicit none

! internal variables
  character*250 junk ! junk variable for reading
  integer ichan   ! (-) channel index
  real(8) temp(6) ! (-) temporary read variable

! print message
  print*, 'Read Instrument Specifications'

! open instrument file
  open(unit=20, file=trim(file_instr), form='formatted', status='old')

! read number of channels
  read(20,*) nchan  ! number of channels

! read in all channel specs
  do ichan=1,nchan
    read(20,*) temp,junk
    instr_spec(ichan)%lo_freq   = temp(1)
    instr_spec(ichan)%if1_freq  = temp(2)
    instr_spec(ichan)%if2_freq  = temp(3)
    instr_spec(ichan)%bandwidth = temp(4)
    instr_spec(ichan)%dtrms     = temp(5)
    instr_spec(ichan)%desig     = temp(6)
    instr_spec(ichan)%num       = ichan
    instr_spec(ichan)%name      = trim(junk)

    instr_spec(ichan)%if1_freq   = instr_spec(ichan)%if1_freq  / 1000.0d0 
    instr_spec(ichan)%if2_freq   = instr_spec(ichan)%if2_freq  / 1000.0d0 
    instr_spec(ichan)%bandwidth  = instr_spec(ichan)%bandwidth / 1000.0d0 
  enddo
!
! close instrument file
  close(unit=20)

end subroutine get_instr_spec

!============================================================================
subroutine extract_channel(ichan )
!============================================================================
! extracts a single channel from the full instrument specification file
!
! history:
!  10/17/2020 Kevin schaefer created routine
!----------------------------------------------------------------------------
  use dotlrt_variables
  implicit none

! input
  integer, intent(in) :: ichan   ! (-) channel index

! copy variable tree
  channel%lo_freq   = instr_spec(ichan)%lo_freq
  channel%if1_freq  = instr_spec(ichan)%if1_freq
  channel%if2_freq  = instr_spec(ichan)%if2_freq
  channel%bandwidth = instr_spec(ichan)%bandwidth
  channel%dtrms     = instr_spec(ichan)%dtrms
  channel%desig     = instr_spec(ichan)%desig 
  channel%num       = instr_spec(ichan)%num 
  channel%name      = instr_spec(ichan)%name 
  
end subroutine extract_channel
