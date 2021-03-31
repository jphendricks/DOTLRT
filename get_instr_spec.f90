! get_instr_spec.f90
! This routine:
!     1 - reads instrument specification file containing parameters:
!         a - number of channels
!         b - frequencies
!         c - bandwidths
!         d - noise parameters
!     2 - computes instrument noise covariance matrix assuming
!         channels are non-overlapping.
subroutine get_instr_spec( instrspec )
use variables
  implicit none
  integer i, j
  real(8) instrspec(5)
    instr_spec%num_channels = 1
    do i = 1, instr_spec%num_channels
      instr_spec%chan(i)%lo_freq    = instrspec(1)
      instr_spec%chan(i)%if1_freq   = instrspec(2)
      instr_spec%chan(i)%if2_freq   = instrspec(3)
      instr_spec%chan(i)%bandwidth  = instrspec(4)
      instr_spec%chan(i)%dtrms      = instrspec(5)
      instr_spec%chan(i)%chan_desig = 1

      instr_spec%chan(i)%if1_freq   = instr_spec%chan(i)%if1_freq  / 1000.0d0 
      instr_spec%chan(i)%if2_freq   = instr_spec%chan(i)%if2_freq  / 1000.0d0 
      instr_spec%chan(i)%bandwidth  = instr_spec%chan(i)%bandwidth / 1000.0d0 
    end do
    instr_spec%text = ' '
    do i=1,instr_spec%num_channels
        do j=1,instr_spec%num_channels
            instr_noise_cov_matrix(i,j) = 0.0d0
        end do
    end do
    do i=1,instr_spec%num_channels
        instr_noise_cov_matrix(i,i) = instr_spec%chan(i)%dtrms &
                                    * instr_spec%chan(i)%dtrms
    end do
    instr_spec%inf=.true.
end subroutine get_instr_spec
