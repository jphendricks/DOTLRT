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
!  12/18/2020 Agnes Lim Update code to handle output filename instead of path
!  12/18/2020 Agnes Lim Add MPI codeh
!-------------------------------------------------------------------------
!
  use mpi
  use dotlrt_variables
  use profiles
  use dotlrt_output
!
  implicit none

! internal variables
  integer ichan        ! (-) channel index
  integer ilon         ! (-) longitude index
  integer ilat         ! (-) latitude index

! Variables for in  
  integer :: input_count

! variables for MPI  
  integer, parameter :: master=0
  integer :: status(MPI_STATUS_SIZE)
  integer :: rank, numtasks, ierror
  integer :: dest, src, chunksize, extra
  integer :: tag1, tag2, tag3, tag4, tag5, tag6
  integer :: stride_start, stride_len, stride_end
  integer :: local_var
  real :: tstart, tend, tstart_io, tend_io, tstart_mrt, tend_mrt, tstart_write, tend_write
  real :: tstart_send, tend_send, tstart_recv, tend_recv
  integer, allocatable :: leftovers(:)

  ! Temporary holding arrays for DOTLRT on each processor  
  real(8), allocatable :: Tbo_temp(:,:), Tbo_streams_temp(:,:,:)
  real(8), allocatable :: dTb_dp_temp(:,:,:,:), dTb_dT_temp(:,:,:,:), dTb_dq_temp(:,:,:,:)
  real(8), allocatable :: dTb_dw_temp(:,:,:,:,:)

! These are for naster processors only  
  real(8), allocatable :: Tbo_mat_proc(:,:), Tbo_streams_mat_proc(:,:,:)
  real(8), allocatable :: dTb_dp_mat_proc(:,:,:,:), dTb_dT_mat_proc(:,:,:,:), dTb_dq_mat_proc(:,:,:,:)
  real(8), allocatable :: dTb_dw_mat_proc(:,:,:,:,:)

! variables 
  integer :: i, j, k, l, m, n
  integer :: total_grid_points

  tag1=1
  tag2=2
  tag3=3
  tag4=4
  tag5=5
  tag6=6

  call MPI_INIT(ierror)
  call CPU_TIME(tstart)
  if(ierror /= MPI_SUCCESS) then
     write(*,*) 'Bad NP'
     stop
  end if
  call MPI_COMM_SIZE(MPI_COMM_WORLD, numtasks, ierror) !N
  if(ierror /= MPI_SUCCESS) then
     write(*,*) 'Bad Size'
     stop
  end if
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror) !0,1,2, .. N-1
  if(ierror /= MPI_SUCCESS) then
    write(*,*) 'Bad Rank'
    stop
  end if

! collect command line input arguments
  if(rank == master) then
     input_count=iargc()
     call get_command_argument(1, file_in)
     call get_command_argument(2, file_out)
     call get_command_argument(3, file_jac_prefix)
     if(input_count==4) then
        call get_command_argument(4, out_path)
     end if
  end if

  call MPI_BCAST(file_in, 250, MPI_CHARACTER, master, MPI_COMM_WORLD, ierror)

! set up all inputs
  call CPU_time(tstart_io)
  call setup_all_inputs( )
  call CPU_time(tend_io)
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)

  total_grid_points=nlon*nlat

  chunksize=total_grid_points/numtasks

! set up all outputs
  if(rank == master) then
     call setup_all_outputs( )
!    This 2 lines for debugging only    
     extra=mod(total_grid_points,numtasks)
     write(*,*) 'chunksize=', chunksize, 'extra=', extra
!     allocate(Tbo_mat_proc(total_grid_points,npol))
!     allocate(Tbo_streams_mat_proc(total_grid_points,nstream/2,npol))
!     allocate(dTb_dp_mat_proc(total_grid_points,nlev,nstream/2,npol))
!     allocate(dTb_dT_mat_proc(total_grid_points,nlev,nstream/2,npol))
!     allocate(dTb_dq_mat_proc(total_grid_points,nlev,nstream/2,npol))
!     allocate(dTb_dw_mat_proc(total_grid_points,nlev,nstream/2,nphase,npol))
!     Tbo_mat_proc(:,:)=0.0
!     Tbo_streams_mat_proc(:,:,:)=0.0
!     dTb_dp_mat_proc(:,:,:,:)=0.0
!     dTb_dT_mat_proc(:,:,:,:)=0.0
!     dTb_dq_mat_proc(:,:,:,:)=0.0
!     dTb_dw_mat_proc(:,:,:,:,:)=0.0
  end if

  allocate(leftovers(numtasks))
  leftovers(:)=0
  extra=mod(total_grid_points,numtasks)
  leftovers(numtasks)=extra

!  allocate(Tbo_temp(chunksize+extra,npol))
!  allocate(Tbo_streams_temp(chunksize+extra,nstream/2,npol))
!  allocate(dTb_dp_temp(chunksize+extra,nlev,nstream/2,npol))
!  allocate(dTb_dT_temp(chunksize+extra,nlev,nstream/2,npol))
!  allocate(dTb_dq_temp(chunksize+extra,nlev,nstream/2,npol))
!  allocate(dTb_dw_temp(chunksize+extra,nlev,nstream/2,nphase,npol))

! run MRT
  do ichan = chan_strt,chan_stop
    call extract_channel(ichan)
    if(flag_print_full) print*, 'Channel: ',ichan, channel%lo_freq
    call CPU_TIME(tstart_mrt)
    do i = 1,chunksize+leftovers(rank+1)
      local_var = rank*chunksize+i
      ilon = (local_var-1) / nlat + 1
      ilat = local_var-(ilon-1)*nlat      
      if (trim(prof_src) == 'WRF') call construct_atmospheric_profile(ichan,ilon,ilat)
      call mrt()
      if(ilat==292 .and. ilon==397) then
        write(400+rank,*) 'Channel = ', ichan
        write(400+rank, *) ilat, ilon 
        do j=1, npol
           write(400+rank, '(8F10.5)') (Tbo_str_mat(k,j), k=1,nstream/2)
        end do
        do j=1, npol
           do k=1, nstream/2
                write(400+rank, '(74E20.12)') (dTb_dp_str_mat(l,k,j), l=1,nlev)
           end do
        end do
        do j=1, npol
           do k=1, nstream/2
                write(400+rank, '(74E20.12)') (dTb_dT_str_mat(l,k,j), l=1,nlev)
           end do
        end do
        do j=1, npol
           do k=1, nstream/2
                write(400+rank, '(74E20.12)') (dTb_dq_str_mat(l,k,j), l=1,nlev)
           end do
        end do
        call flush(400+rank)
        write(600+rank,*) 'Channel = ', ichan
        write(600+rank, *) ilat, ilon 
        do j=1,5
           do k=1,npol
              do l=1,nstream/2
                 write(600+rank,'(74E20.12)') (dTb_dw_str_mat(m,l,j,k), m=1,nlev)
              end do
           end do
        end do 
        call flush(500+rank)
      end if
!       write(300+rank,fmt='(I7,I12,2I5)') i, local_var, ilon, ilat
!      Tbo_temp(i,:)=Tbo_mat
!      Tbo_streams_temp(i,:,:)=Tbo_str_mat
!      dTb_dp_temp(i,:,:,:)=dTb_dp_str_mat
!      dTb_dT_temp(i,:,:,:)=dTb_dT_str_mat
!      dTb_dq_temp(i,:,:,:)=dTb_dq_str_mat
!      dTb_dw_temp(i,:,:,:,:)=dTb_dw_str_mat
    end do
    call CPU_TIME(tend_mrt)

    CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)

    call CPU_TIME(tstart_send)
    if(rank > master) then
!     Send BT to master
      n=chunksize+leftovers(rank+1)
!      call MPI_send(Tbo_temp(1:n,:), n*npol, MPI_DOUBLE_PRECISION, master, tag1, MPI_COMM_WORLD, ierror)
!      call MPI_send(Tbo_streams_temp(1:n,:,:), n*(nstream/2)*npol,MPI_DOUBLE_PRECISION, master, tag2, MPI_COMM_WORLD, ierror)
!      call MPI_send(dTb_dp_temp(1:n,:,:,:), n*nlev*(nstream/2)*npol,MPI_DOUBLE_PRECISION, master, tag3, MPI_COMM_WORLD, ierror)
!      call MPI_send(dTb_dT_temp(1:n,:,:,:), n*nlev*(nstream/2)*npol,MPI_DOUBLE_PRECISION, master, tag4, MPI_COMM_WORLD, ierror)
!      call MPI_send(dTb_dq_temp(1:n,:,:,:), n*nlev*(nstream/2)*npol,MPI_DOUBLE_PRECISION, master, tag5, MPI_COMM_WORLD, ierror)
!      call MPI_send(dTb_dw_temp(1:n,:,:,:,:), n*nlev*(nstream/2)*nphase*npol,MPI_DOUBLE_PRECISION, master, tag6, MPI_COMM_WORLD, ierror)
    end if
    call CPU_TIME(tend_send)

    call CPU_TIME(tstart_recv)
    if(rank == master) then
!     BT on rank 0 stays.
      stride_len=chunksize+leftovers(rank+1)
      stride_start=1
      stride_end=stride_start+stride_len-1
!      Tbo_mat_proc(stride_start:stride_end,:)=Tbo_temp
!      Tbo_streams_mat_proc(stride_start:stride_end,:,:)=Tbo_streams_temp
!      dTb_dp_mat_proc(stride_start:stride_end,:,:,:)=dTb_dp_temp
!      dTb_dT_mat_proc(stride_start:stride_end,:,:,:)=dTb_dT_temp
!      dTb_dq_mat_proc(stride_start:stride_end,:,:,:)=dTb_dq_temp
!      dTb_dw_mat_proc(stride_start:stride_end,:,:,:,:)=dTb_dw_temp

!    BT from rank > 0
      do src=1, numtasks-1
        n=chunksize+leftovers(src+1)
!        call MPI_recv(Tbo_temp(1:n,:),n*npol,MPI_DOUBLE_PRECISION, src, tag1,MPI_COMM_WORLD, status, ierror)
!        call MPI_recv(Tbo_streams_temp(1:n,:,:),n*(nstream/2)*npol,MPI_DOUBLE_PRECISION, src, tag2, MPI_COMM_WORLD, status, ierror)
!        call MPI_recv(dTb_dp_temp(1:n,:,:,:),n*nlev*(nstream/2)*npol,MPI_DOUBLE_PRECISION, src, tag3, MPI_COMM_WORLD, status, ierror)
!        call MPI_recv(dTb_dT_temp(1:n,:,:,:),n*nlev*(nstream/2)*npol,MPI_DOUBLE_PRECISION, src, tag4, MPI_COMM_WORLD, status, ierror)
!        call MPI_recv(dTb_dq_temp(1:n,:,:,:),n*nlev*(nstream/2)*npol,MPI_DOUBLE_PRECISION, src, tag5, MPI_COMM_WORLD, status, ierror)
!        call MPI_recv(dTb_dw_temp(1:n,:,:,:,:),n*nlev*(nstream/2)*nphase*npol,MPI_DOUBLE_PRECISION, src, tag6, MPI_COMM_WORLD, status, ierror)
        stride_len=chunksize+leftovers(src+1)
        stride_start=src*chunksize+1
        stride_end=stride_start+stride_len-1
!        Tbo_mat_proc(stride_start:stride_end,:)=Tbo_temp
!        Tbo_streams_mat_proc(stride_start:stride_end,:,:)=Tbo_streams_temp
!        dTb_dp_mat_proc(stride_start:stride_end,:,:,:)=dTb_dp_temp
!        dTb_dT_mat_proc(stride_start:stride_end,:,:,:)=dTb_dT_temp
!        dTb_dq_mat_proc(stride_start:stride_end,:,:,:)=dTb_dq_temp
!        dTb_dw_mat_proc(stride_start:stride_end,:,:,:,:)=dTb_dw_temp
      end do
    end if
    call CPU_TIME(tend_recv)

    CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)

!    if(rank == master) then
!       do i=1,total_grid_points
!         ilon = (i-1) / nlat + 1
!         ilat = i-(ilon-1)*nlat
!         Tbo_wrt(ilon,ilat,:)=Tbo_mat_proc(i,:)
!         Tbo_str_wrt(ilon,ilat,:,:)=Tbo_streams_mat_proc(i,:,:)
!         dTb_dp_str_wrt(ilon,ilat,:,:,:)=dTb_dp_mat_proc(i,:,:,:)
!         dTb_dT_str_wrt(ilon,ilat,:,:,:)=dTb_dT_mat_proc(i,:,:,:)
!         dTb_dq_str_wrt(ilon,ilat,:,:,:)=dTb_dq_mat_proc(i,:,:,:)
!         dTb_dc_str_wrt(ilon,ilat,:,:,:)=dTb_dw_mat_proc(i,:,:,1,:)
!         dTb_dr_str_wrt(ilon,ilat,:,:,:)=dTb_dw_mat_proc(i,:,:,2,:)
!         dTb_di_str_wrt(ilon,ilat,:,:,:)=dTb_dw_mat_proc(i,:,:,3,:)
!         dTb_ds_str_wrt(ilon,ilat,:,:,:)=dTb_dw_mat_proc(i,:,:,4,:)
!         dTb_dg_str_wrt(ilon,ilat,:,:,:)=dTb_dw_mat_proc(i,:,:,5,:)
!       end do
!       call CPU_TIME(tstart_write)
!       if(save_rad_file) call write_tb_channel(ichan)
!       call CPU_TIME(tend_write)
!       if(save_jac_file) call create_jacobian_file(ichan)
!    end if
!    write(*,fmt='(A20,I4,3F12.6)') 'Time of each proc:', rank, tend_mrt-tstart_mrt, tend_send-tstart_send, tend_recv-tstart_recv
  enddo

 

  call CPU_TIME(tend)
!  write(*,fmt='(A20,I4,2F12.6)') 'Time of each proc:', rank, tend-tstart, tend_io-tstart_io

  call MPI_FINALIZE(ierror)
  if(ierror /= MPI_SUCCESS) then
         write(*,*) 'Bad Finalize'
         stop
  end if
! print message
  print*, 'DOTLRT Run Complete'

end program dot
