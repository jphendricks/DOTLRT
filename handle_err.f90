!==================================================
subroutine handle_err(status, subroutine, line)
!==================================================
! prints netcdf error messages
!
! history
!   9/30/20 Kevin Schaefer added routine and line number flags
!   9/30/20 Kevin Schaefer deleted if status /= nf90_noerr
!--------------------------------------------------
use netcdf
implicit none
            
integer, intent(in) :: status          ! netcdf operation status variable
character*50, intent(in) :: subroutine ! subroutine where error occured
integer, intent(in) :: line            ! line where error ocured
      
print*, trim(subroutine), ' at ', line
print*, trim(nf90_strerror(status))
stop 'Program'
      
end subroutine handle_err
