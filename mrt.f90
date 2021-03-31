!========================================================================================
subroutine mrt( )
!========================================================================================
! History:
!  5/272005 Bob Weber created program
!  9/26/2020 Kevin Schaefer deleted unused variables
!  10/12/2020 Kevin Schaefer cleaned up code
!  10/15/2020 Kevin Schaefer deleted all arguments duplicated in variables module 
!  10/15/2020 Kevin Schaefer deleted all arguments duplicated in data_assimilation module 
! ---------------------------------------------------------------------------------------

use dotlrt_variables
use dotlrt_output

implicit none

! mex
integer(4) jpol
real(8), dimension(2) :: tau


integer i, j, iphase
real(8) Tbo
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
allocate(Tb_obs(0:nlr1,2))
allocate(dTb_dT_obs(nlr1,2))
allocate(dTb_dp_obs(nlr1,2))
allocate(dTb_dq_obs(nlr1,2))
allocate(dTb_dw_obs(nlr1,5,2))
allocate(Tbo_streams(nangover2))
allocate(dTb_dT_streams(nlr1,nangover2,2))
allocate(dTb_dp_streams(nlr1,nangover2,2))
allocate(dTb_dq_streams(nlr1,nangover2,2))
allocate(dTb_dw_streams(nlr1,nangover2,5,2))

! model hydrometeor distribution parameters
call calcprofile_d()
      
do jpol = 1, 2 ! horizontal, vertical polarization
  ipol = jpol - 1  ! polarization number 1=vertical 0=horizontal
  call calc_fbw_temp_weight_scat( Tbo, tau, Tb_obs, dTb_dT_obs, dTb_dp_obs, &
                                        dTb_dq_obs, dTb_dw_obs, Tbo_streams, &
                                        dTb_dT_streams, dTb_dp_streams, dTb_dq_streams, &
                                        dTb_dw_streams)

  Tbo_mat(jpol) = Tbo
  do j = 1,nangover2
    Tbo_str_mat(j,jpol) = Tbo_streams(j)
  end do

  do i = 1, nlev
    do j = 1, nangover2
      dTb_dT_str_mat(i,j,jpol) = dTb_dT_streams(i,j,1)
      dTb_dp_str_mat(i,j,jpol) = dTb_dp_streams(i,j,1)
      dTb_dq_str_mat(i,j,jpol) = dTb_dq_streams(i,j,1)
      do iphase = 1, nphase
        dTb_dw_str_mat(i,j,iphase,jpol) = dTb_dw_streams(i,j,iphase,1)
      end do
    end do
  end do

  do j = 1, 2
    tau_mat(j) = tau(j)
  end do

  do j = 0, nlev
    ! Tb_obs(j,1) - upwelling; Tb_obs(j,2) - downwelling
    Tb_obs_mat(j,jpol) = Tb_obs(j,1) 
  end do
                
  do j = 1, nlev
    dTb_dT_obs_mat(j,jpol) = dTb_dT_obs(j,1)
    dTb_dp_obs_mat(j,jpol) = dTb_dp_obs(j,1)
    dTb_dq_obs_mat(j,jpol) = dTb_dq_obs(j,1)
    do iphase = 1, nphase
      dTb_dw_obs_mat(j,iphase,jpol) = dTb_dw_obs(j,iphase,1)
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
