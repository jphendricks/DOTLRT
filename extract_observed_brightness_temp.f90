!=======================================================================
  subroutine extract_observed_brightness_temp (point)
!=======================================================================
! extracts observed brightness temperatures from swath
!
! Modifications:
!  7/2/2021  Kevin Schaefer created subroutine
!----------------------------------------------------------------------
  use scan_Variables
  use dotlrt_variables
  use profiles
  use dotlrt_output

  implicit none

! input variables
  integer point ! (-) point number to extract

! local variables

! transfer from swath to observed brightness temperatures
! print*, Tbo_obs(:,1), swath_tb_chan(:,point)
  Tbo_obs(:,1) = swath_tb_chan(:,point)
  Tbo_obs(:,2) = Tbo_obs(:,1)

  return
  end subroutine extract_observed_brightness_temp

!!=======================================================================
!  subroutine read_observed_brightness_temp ()
!!=======================================================================
!! reads in observed brightness temperatures
!!
!! Modifications:
!!  4/8/2021  Kevin Schaefer created subroutine
!!  5/21/2021 Kevin Schaefer removed print stat ment
!!----------------------------------------------------------------------
!  use scan_Variables
!  use dotlrt_variables
!  use profiles
!  use dotlrt_output
!
!  implicit none
!
!! local variables
!  integer ivar    ! (-) variable index
!  integer ichan   ! (-) channel index
!  integer ifil    ! (-) file index
!  integer count   ! (-) channel index
!  integer status  ! (-) read status variable
!  integer num  ! (-) read status variable
!  real(8) val  ! (-) input read variable
!  real(8) obs(nchan) ! (k) full observed brightness temp
!  character*250 filename  ! file name
!
!! find observation file
!  do ifil=1,numfiles
!    if (trim(files(ifil)%type) == 'FY3_text')  filename = trim(files(ifil)%path)
!  enddo
!  print*, '    read Obs Tb: ', trim(filename)
!
!! check number of channels
!  open(unit=20, file=trim(filename), form='formatted', status='old')
!  read(20,*) junk ! read header
!  count=0
!  do ichan=1,100
!    read(20,*, iostat=status) num, val
!    if(status<0) exit
!    count=count+1
!  enddo
!  close(unit=20)
!  if(count/=nchan) then
!    print*, 'Error: wrong Tb observation file'
!    print*, trim(filename)
!    print*, count, nchan
!    stop
!  endif
!
!! read in observations
!  open(unit=20, file=trim(filename), form='formatted', status='old')
!  read(20,*) junk ! read header
!  do ichan=1,nchan
!    read(20,*) num, val
!    obs(ichan) = val
!  enddo
!  close(unit=20)
!
!! transfer to observed brightness temperatures
!  ivar=0
!  do ichan=minchan, maxchan
!    ivar=ivar+1
!    Tbo_obs(ivar,:) = obs(ichan)
!  enddo
!
!  return
!  end
!
!!============================================================================
!subroutine read_init_precip_lookup_table( )
!!============================================================================
!! This routine reads in an initial brightness temperature lookup table
!!
!! history:
!!  6/29/2021 Kevin schaefer created routine
!!----------------------------------------------------------------------------
!  use scan_Variables
!  use dotlrt_variables
!  implicit none
!
!! internal variables
!  character*250 filename ! (-) local filename
!  integer status  ! (-) read status
!  integer ifil    ! (-) file index
!  integer irec    ! (-) record index
!  integer max     ! (-) max allowed number of records
!
!! find init Tb filename
!  do ifil=1,numfiles
!    if (trim(files(ifil)%type) == 'init_precip')  filename = trim(files(ifil)%path)
!  enddo
!
!! count records
!  open(unit=20, file=trim(filename), form='formatted', status='old')
!  read(20,*, iostat=status) junk
!  init_nlook=0
!  max = 1000
!  do irec=1,max
!    read(20,*, iostat=status) junk
!    if(status<0) exit
!    init_nlook=init_nlook+1
!  enddo
!  close(unit=20)
!  print*, '    Number lookup Table records:', init_nlook
!
!! allocate lookup table
!  allocate(init_look(init_nlook,nchan))
!  allocate(init_dens(init_nlook))
!  allocate(init_precip(nchannel))
!
!! open init tb file
!  open(unit=20, file=trim(filename), form='formatted', status='old')
!  read(20,*, iostat=status) junk
!
!! read in records
!  do irec=1,init_nlook
!    read(20,*) init_dens(irec),init_look(irec,:)
!    if(dabs(init_dens(irec)) < 1.d-6) init_dens(irec) = 0.d0
!    !print*, irec, init_dens(irec), init_precip(irec,1)
!  enddo
!!
!! close lookup table file
!  close(unit=20)
!
!end subroutine read_init_precip_lookup_table
!
!!============================================================================
!subroutine calc_init_precip( )
!!============================================================================
!! This routine calculates initial total column precipitation from a lookup table
!!
!! history:
!!  6/29/2021 Kevin schaefer created routine
!!----------------------------------------------------------------------------
!  use scan_Variables
!  use dotlrt_variables
!  implicit none
!
!! internal variables
!  integer ichan    ! (-) file index
!  integer ivar    ! (-) record index
!  integer irec    ! (-) record index
!  real(8) val  ! (-) tb ratio
!  real(8) max  ! (K) max tb in lookup table
!
!! find nearest Tb
!  ivar=0
!  do ichan=minchan, maxchan
!    ivar=ivar+1
!    max = maxval(init_look(:,ichan))
!    bad_point = .false.
!    if(Tbo_obs(ivar,1) > max) then
!      !print*, 'bad point ', max, Tbo_obs(ivar,1)
!      bad_point = .true.
!    else
!      do irec=init_nlook, 1, -1
!        if(init_look(irec,ichan) >= Tbo_obs(ivar,1)) then
!          val=(Tbo_obs(ivar,1)-init_look(irec,ichan))/(init_look(irec+1,ichan)-init_look(irec,ichan))
!          init_precip(ivar) = init_dens(irec)+(init_dens(irec+1)-init_dens(irec))*val
!          !print*, irec
!          !print*, init_look(irec+1,ichan),Tbo_obs(ivar,1),init_look(irec,ichan)
!          !print*, init_dens(irec+1),init_precip(ivar), init_dens(irec)
!          exit
!        endif
!      enddo
!    endif
!  enddo
!
!end subroutine calc_init_precip
!!==============================================================================
!subroutine read_FY_swath()
!!==============================================================================
!! This program reads in FY data files
!!
!! History:
!!  6/30/2021 Kevin Schaefer created routine
!!------------------------------------------------------------------------------
!  use scan_Variables
!  use dotlrt_variables
!  use profiles
!  use netcdf
!
!implicit none
!!
!! Local variables
!  integer ncid       ! (-) netcdf file id
!  integer groupid    ! (-) netcdf group id
!  integer varid      ! (-) netcdf variable id
!  integer numdim     ! (-) number of dimensions
!  integer dimlen     ! (-) value of dimension
!  integer status     ! (-) netcdf status code
!  integer dimid(20)  ! (-) dimension IDs
!  integer parent     ! (-) track parent id
!  character(len = nf90_max_name) :: dimname
!  character*50 subroutine ! subroutine name
!  character*250 filename  ! file name
!  integer idim       ! (-) dimension index
!  integer ifil       ! (-) file index
!  integer ivar       ! (-) variable index
!  integer ichan      ! (-) channel index
!!
!! set subroutine name
!  subroutine='read_FY_netcdf'
!
!! find masked observation swath file
!  do ifil=1,numfiles
!    if (trim(files(ifil)%type) == 'FY3_vec')  filename = trim(files(ifil)%path)
!  enddo
!  print*, '    read Obs FY: ', trim(filename)
!
!! Open swath file
!  status = nf90_open(filename, NF90_NOWRITE, ncid)
!  if (status /= nf90_noerr) then
!    call handle_err(status, subroutine, 1)
!  endif
!
!! get Data group id
!  status = nf90_inq_ncid(ncid, 'aoi_data', groupid)
!  if (status /= nf90_noerr) call handle_err(status, subroutine, 2)
!
!! get dimension ids of tb data
!  parent=0
!  status = nf90_inq_dimids (groupid, numdim, dimid, parent)
!    if (status /= nf90_noerr) call handle_err(status, subroutine, 3)
!
!! get tb swath dimensions
!  do idim = 1,numdim
!    status = nf90_inquire_dimension(groupid, dimid(idim), dimname, dimlen)
!    if (status /= nf90_noerr) call handle_err(status, subroutine, 4)
!    if (trim(dimname) == 'phony_dim_0') then
!      swath_npts = dimlen
!      print*, '      Number Points: ', dimlen
!    endif
!    if (trim(dimname) == 'phony_dim_1') then
!      swath_nchan = dimlen
!      print*, '      Number Channels: ', dimlen
!    endif
!    if (trim(dimname) == 'phony_dim_2') then
!      swath_nlev = dimlen
!      print*, '      Number Levels: ', dimlen
!    endif
!  enddo
!
!! allocate swath data variables
!  allocate(swath_tb_full(swath_nchan, swath_npts))
!  allocate(swath_tb_chan(nchannel,swath_npts))
!  allocate(swath_angle(swath_npts))
!  allocate(swath_temp(swath_nlev,swath_npts))
!
!! allocate optimization metrics
!  allocate(swath_precip(swath_npts))
!  allocate(swath_point(swath_npts))
!  allocate(swath_iter(swath_npts))
!  allocate(swath_cost(swath_npts))
!  allocate(swath_tb_sim(swath_npts))
!  allocate(swath_res(swath_npts))
!  allocate(swath_ex(swath_npts))
!  allocate(swath_bad(swath_npts))
!
!! initialize metrics arrays to standard missing value
!  swath_precip(:) = -999.d0
!  swath_point(:)  = -999.d0
!  swath_iter(:)   = -999.d0
!  swath_cost(:)   = -999.d0
!  swath_tb_sim(:) = -999.d0
!  swath_res(:)    = -999.d0
!  swath_ex(:)     = -999.d0
!  swath_bad(:)    = -999.d0
!
!! get brightness temperatures
!  status = nf90_inq_varid(groupid, 'L1c_bt', varid)
!  if (status /= nf90_noerr) call handle_err(status, subroutine, 5)
!  status = nf90_get_var(groupid, varid, swath_tb_full)
!  if (status /= nf90_noerr) call handle_err(status, subroutine, 6)
!
!! get nadir angles
!  status = nf90_inq_varid(groupid, 'L1c_inc_ang', varid)
!  if (status /= nf90_noerr) call handle_err(status, subroutine, 7)
!  status = nf90_get_var(groupid, varid, swath_angle)
!  if (status /= nf90_noerr) call handle_err(status, subroutine, 8)
!
!! get vertical temperature profiles
!  status = nf90_inq_varid(groupid, 'temperature', varid)
!  if (status /= nf90_noerr) call handle_err(status, subroutine, 9)
!  status = nf90_get_var(groupid, varid, swath_temp)
!  if (status /= nf90_noerr) call handle_err(status, subroutine, 10)
!
!! Close File
!  status = nf90_close(ncid  )
!  if (status /= nf90_noerr) call handle_err(status, subroutine, 11)
!
!! subset brightness temperatures by channel
!  ivar=0
!  do ichan=minchan, maxchan
!    ivar=ivar+1
!    swath_tb_chan(ivar,:) = swath_tb_full(ichan,:)
!  enddo
!
!end subroutine read_FY_swath
!
