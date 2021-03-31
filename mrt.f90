!========================================================================================
subroutine mrt( )
!========================================================================================
! History:
!  5/27/2005  Bob Weber created program
!  9/26/2020  Kevin Schaefer deleted unused variables
!  10/12/2020 Kevin Schaefer cleaned up code
!  10/15/2020 Kevin Schaefer deleted all arguments duplicated in variables module 
!  10/15/2020 Kevin Schaefer deleted all arguments duplicated in data_assimilation module 
!  12/13/2020 Kevin Schaefer converted to readable indeces
! ---------------------------------------------------------------------------------------

use dotlrt_variables
use dotlrt_output

implicit none

! mex
real(8), dimension(2) :: tau

integer iup    ! (-) up down index
integer iang   ! (-) quad angle index index
integer jpol   ! (-) polarization index
integer ilev   ! (-) atmospheric level index
integer iphase ! (-) hydrometeor phase index

real(8), allocatable :: Tb_obs(:,:)
real(8), allocatable :: dTb_dT_obs(:,:)
real(8), allocatable :: dTb_dp_obs(:,:)
real(8), allocatable :: dTb_dq_obs(:,:)
real(8), allocatable :: dTb_dw_obs(:,:,:)
real(8), allocatable :: Tbo_streams(:)
real(8), allocatable :: dTb_dT_streams(:,:,:)
real(8), allocatable :: dTb_dp_streams(:,:,:)
real(8), allocatable :: dTb_dq_streams(:,:,:)
real(8), allocatable :: dTb_dw_streams(:,:,:,:)

! allocate variables
allocate(Tb_obs(0:nlev,updown))
allocate(dTb_dT_obs(nlev,updown))
allocate(dTb_dp_obs(nlev,updown))
allocate(dTb_dq_obs(nlev,updown))
allocate(dTb_dw_obs(nlev,nphase,updown))
allocate(Tbo_streams(nangover2))
allocate(dTb_dT_streams(nlev,nangover2,updown))
allocate(dTb_dp_streams(nlev,nangover2,updown))
allocate(dTb_dq_streams(nlev,nangover2,updown))
allocate(dTb_dw_streams(nlev,nangover2,nphase,updown))

! model hydrometeor distribution parameters
  call calcprofile_d()

do jpol = 1, npol ! 1 = horizontal, 2 = vertical polarization
  call calc_fbw_temp_weight_scat( Tbo_mat(jpol), tau, Tb_obs, dTb_dT_obs, dTb_dp_obs, dTb_dq_obs, dTb_dw_obs, &
    Tbo_streams, dTb_dT_streams, dTb_dp_streams, dTb_dq_streams, dTb_dw_streams)

  do iang = 1,nangover2
    Tbo_str_mat(iang,jpol) = Tbo_streams(iang)
  end do

  do ilev = 1, nlev
    do iang = 1, nangover2
      dTb_dT_str_mat(ilev,iang,jpol) = dTb_dT_streams(ilev,iang,1)
      dTb_dp_str_mat(ilev,iang,jpol) = dTb_dp_streams(ilev,iang,1)
      dTb_dq_str_mat(ilev,iang,jpol) = dTb_dq_streams(ilev,iang,1)
      do iphase = 1, nphase
        dTb_dw_str_mat(ilev,iang,iphase,jpol) = dTb_dw_streams(ilev,iang,iphase,1)
      end do
    end do
  end do

  do iup = 1, updown
    tau_mat(iup) = tau(iup)
  end do

  do ilev = 0, nlev
    Tb_obs_mat(ilev,jpol) = Tb_obs(ilev,1) 
  end do
                
  do ilev = 1, nlev
    dTb_dT_obs_mat(ilev,jpol) = dTb_dT_obs(ilev,1)
    dTb_dp_obs_mat(ilev,jpol) = dTb_dp_obs(ilev,1)
    dTb_dq_obs_mat(ilev,jpol) = dTb_dq_obs(ilev,1)
    do iphase = 1, nphase
      dTb_dw_obs_mat(ilev,iphase,jpol) = dTb_dw_obs(ilev,iphase,1)
    end do
  end do
end do ! jpol

! deallocate memory
deallocate(Tb_obs)
deallocate(dTb_dT_obs)
deallocate(dTb_dp_obs)
deallocate(dTb_dq_obs)
deallocate(dTb_dw_obs)
deallocate(Tbo_streams)
deallocate(dTb_dT_streams)
deallocate(dTb_dp_streams)
deallocate(dTb_dq_streams)
deallocate(dTb_dw_streams)

end ! program mrt

!========================================================================================
subroutine ex_time(num_segment, id_segment)
!========================================================================================
! History:
!  1/27/2021 Kevin Schaefer created routine
! ---------------------------------------------------------------------------------------

use dotlrt_variables
use dotlrt_output

implicit none

! input variables
integer num_segment
character*40 id_segment

! calculate execution time
  call system_clock(time_stop,time_rate) ! stop counting
  time_del= real(time_stop-time_start)/real(time_rate)
  time_seg(num_segment) = time_seg(num_segment)+time_del
  num_call(num_segment)=num_call(num_segment)+1
  seg_name(num_segment)= trim(id_segment)
  nseg=num_segment

end subroutine ex_time
