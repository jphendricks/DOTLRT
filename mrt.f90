! mrt.f90
! DOTLRTv1.0
! 27 May 2005    Bob Weber
! program mrt
! mex
subroutine mrt( inpheight,              & ! observation height (AGL, km)
                inptheta,               & ! observation angle from nadir (degrees)
                numsbfreqs,             & ! number of frequency points per sideband
                instrspec,              &
                nstreams,               & ! number of streams (total, up and down)
                stream_angles,          & ! stream angles
                nsurfangles,            & ! number of angles for surface reflectivity
                surfinp,                & ! surface reflectivity (vertical, horizontal polarization)
                                          ! for specified angles of incidence
                atm_inp_num_levels,     & ! excluding surface level
                atminp,                 &
                Tb_obs_mat, dTb_dT_obs_mat,      &
                dTb_dp_obs_mat, dTb_dq_obs_mat, &
                dTb_dw_obs_mat, Tbo_mat, tau_mat, Tbo_streams_mat, &
                dTb_dT_streams_mat, dTb_dp_streams_mat, dTb_dq_streams_mat, &
                dTb_dw_streams_mat)
! mex
  use variables
  implicit none

  ! mex
  integer numsbfreqs
  real(8) inpheight, inptheta
  real(8) instrspec(5)
  integer(4) atm_inp_num_levels, nstreams, nsurfangles, jpol
  real(8) surfinp(nsurfangles,3)
  real(8) atminp(atm_inp_num_levels+1,9), stream_angles(nstreams)
  real(8), dimension(0:atm_inp_num_levels,2) :: Tb_obs_mat
  real(8), dimension(atm_inp_num_levels,2) :: dTb_dT_obs_mat
  real(8), dimension(atm_inp_num_levels,2) :: dTb_dp_obs_mat
  real(8), dimension(atm_inp_num_levels,2) :: dTb_dq_obs_mat
  real(8), dimension(atm_inp_num_levels,5,2) :: dTb_dw_obs_mat
  real(8), dimension(2) :: Tbo_mat, tau, tau_mat
  real(8), dimension(nstreams/2, 2) :: Tbo_streams_mat
  real(8), dimension(atm_inp_num_levels,nstreams/2,2) :: dTb_dT_streams_mat, dTb_dp_streams_mat, dTb_dq_streams_mat
  real(8), dimension(atm_inp_num_levels,nstreams/2,5,2) :: dTb_dw_streams_mat
  ! mex

  integer i, j, k, jlatlon, jchan
  real(8) Tbo
  real(8), dimension(:,:), allocatable :: Tb_obs, dTb_dT_obs, dTb_dp_obs, dTb_dq_obs
  real(8), dimension(:,:,:), allocatable :: dTb_dw_obs
  real(8), dimension(:), allocatable :: Tbo_streams
  real(8), dimension(:,:,:), allocatable :: dTb_dT_streams, dTb_dp_streams, dTb_dq_streams
  real(8), dimension(:,:,:,:), allocatable :: dTb_dw_streams

  character*120 debugout

  ! establish configuration parameters, alocate memory !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call configure( inpheight,          & ! observation height (AGL, km)
                    inptheta,           & ! observation angle from nadir (degrees)
                    numsbfreqs,         & ! number of frequency points per sideband
                    instrspec,          &
                    nstreams,           & ! number of streams (total, up and down)
                    stream_angles,      & ! stream angles
                    nsurfangles,        & ! number of angles for surface reflectivity
                    surfinp,            & ! surface reflectivity (vertical, horizontal polarization)
                    atm_inp_num_levels  )

       allocate( Tb_obs(0:nlr1,2), &
                 dTb_dT_obs(nlr1,2),   &
                 dTb_dp_obs(nlr1,2),   &
                 dTb_dq_obs(nlr1,2),   &
                 dTb_dw_obs(nlr1,5,2), &
                 Tbo_streams(nangover2), &
                 dTb_dT_streams(nlr1,nangover2,2), &
                 dTb_dp_streams(nlr1,nangover2,2), &
                 dTb_dq_streams(nlr1,nangover2,2), &
                 dTb_dw_streams(nlr1,nangover2,5,2))

       call getprofile( atminp )

      ! model hydrometeor distribution parameters
      call calcprofile_d()
      
      do jpol = 1, 2 ! horizontal, vertical
        inp_pol = jpol - 1
        call calc_fbw_temp_weight_scat( Tbo, tau, Tb_obs, dTb_dT_obs, dTb_dp_obs, &
                                        dTb_dq_obs, dTb_dw_obs, Tbo_streams, &
                                        dTb_dT_streams, dTb_dp_streams, dTb_dq_streams, &
                                        dTb_dw_streams)

        Tbo_mat(jpol) = Tbo
        do j = 1,nangover2
           Tbo_streams_mat(j,jpol) = Tbo_streams(j)
        end do

        do i = 1, atm_inp_num_levels
           do j = 1, nangover2
              dTb_dT_streams_mat(i,j,jpol) = dTb_dT_streams(i,j,1)
              dTb_dp_streams_mat(i,j,jpol) = dTb_dp_streams(i,j,1)
              dTb_dq_streams_mat(i,j,jpol) = dTb_dq_streams(i,j,1)
              do k = 1, number_h2o_phases
                 dTb_dw_streams_mat(i,j,k,jpol) = dTb_dw_streams(i,j,k,1)
              end do
           end do
        end do

        do j = 1, 2
          tau_mat(j) = tau(j)
        end do

        do j = 0, atm_inp_num_levels
          ! Tb_obs(j,1) - upwelling; Tb_obs(j,2) - downwelling
          Tb_obs_mat(j,jpol) = Tb_obs(j,1) 
        end do

        ! for debug
        !do j=0,atm_inp_num_levels
        !   write(debugout,*) Tb_obs(j,1)
        !   call mexPrintf(debugout)
        !end do
        !call mexPrintf(achar(10))
                
        do j = 1, atm_inp_num_levels
          dTb_dT_obs_mat(j,jpol) = dTb_dT_obs(j,1)
          dTb_dp_obs_mat(j,jpol) = dTb_dp_obs(j,1)
          dTb_dq_obs_mat(j,jpol) = dTb_dq_obs(j,1)
          do k = 1, number_h2o_phases
            dTb_dw_obs_mat(j,k,jpol) = dTb_dw_obs(j,k,1)
          end do
        end do
      end do ! jpol

      ! deallocate memory, close input / output files - start !!!!!!!!!!!!!!!!!!
        deallocate(     Tb_obs, &
                    dTb_dT_obs, &
                    dTb_dp_obs, &
                    dTb_dq_obs, &
                    dTb_dw_obs, &
                    Tbo_streams, &
                    dTb_dT_streams, &
                    dTb_dp_streams, &
                    dTb_dq_streams, &
                    dTb_dw_streams)
end ! program mrt
