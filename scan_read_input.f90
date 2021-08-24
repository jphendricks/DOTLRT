!
!=======================================================================
  subroutine scan_read_input
!=========================================================================
! reads input parameters for the scan program
!
! Modifications:
!  3/34/2004 Kevin Schaefer split off from main program
!  1/23/2021 Kevin Schaefer cleaned up code for MRT
!--------------------------------------------------------------------------
  use scan_Variables
  use dotlrt_variables
  use profiles
  use dotlrt_output

  implicit none

! Internal variables
  integer ivar  ! (-) variable index
  integer ifil  ! (-) file index
  integer iman  ! (-) manipulation index
  integer ileg  ! (-) legend index
  integer num   ! (-) variable number
!
!--------------------------------------------------------------------------
! input control variables for various options
!--------------------------------------------------------------------------
! Open input file
  OPEN (9, File='Scan.in', Status='Unknown')
!
! read execution control variables
! junk=garbage variable for variable descriptors in input file
!
! read input file header line
  read (9,10) junk
!
! read input variables
  read (9,16) junk, scantype
!
! scanning options
  read (9,10) junk
  read (9,11) junk, flagDim
  read (9,11) junk, Xvar
  read (9,11) junk, YVar
  read (9,11) junk, nline
  read (9,11) junk, numpts
  read (9,11) junk, minchan
  read (9,11) junk, maxchan
  read (9,11) junk, save_2d
  read (9,11) junk, save_3d
  read (9,17) junk, src_obs
  read (9,17) junk, src_init
  read (9,12) junk, calc_psuedo
  read (9,11) junk, n_psuedo
  read (9,10) junk
  do ivar=1,n_psuedo
    read (9,*) num, psuedo_man(ivar)%doit, psuedo_man(ivar)%typ, &
    psuedo_man(ivar)%ind1, psuedo_man(ivar)%ind2, &
    psuedo_man(ivar)%val1, psuedo_man(ivar)%val2, psuedo_man(ivar)%val3
  enddo
  read (9,12) junk, calc_init
  read (9,11) junk, n_init
  read (9,10) junk
  do ivar=1,n_init
    read (9,*) num, init_man(ivar)%doit, init_man(ivar)%typ, &
    init_man(ivar)%ind1, init_man(ivar)%ind2, &
    init_man(ivar)%val1, init_man(ivar)%val2, init_man(ivar)%val3
  enddo
!
! plotting options
  read (9,10) junk
  read (9,11) junk, PlotFlag
  read (9,11) junk, color_con
  read (9,16) junk, LabY
  read (9,11) junk, Yscale
  read (9,15) junk, Ymin
  read (9,15) junk, Ymax
  read (9,17) junk, Title
  read (9,11) junk, LegFlag
  read (9,11) junk, nLeg
  read (9,10) junk
  do ileg=1,nLeg
    read (9,*) num, Legend(ileg)
  enddo
  read (9,11) junk, DevFlag
  read (9,17) junk, plotfile
  read (9,11) junk, Scale
  read (9,15) junk, base
  read (9,11) junk, NC
  read (9,15) junk, mincont
  read (9,15) junk, maxcont
  read (9,20) junk, barformat
!
! input files
  read (9,10) junk
  read (9,11) junk, numfiles
  allocate(files(numfiles))
  read (9,10) junk
  do ifil=1,numfiles
    read(9,*) files(ifil)%num, files(ifil)%type,files(ifil)%path
  enddo
!
! output manipulations
  read (9,10) junk
  read (9,11) junk,nout_man
  read (9,10) junk
  do iman = 1, nout_man
    read (9,*) out_man(iman)%doit, out_man(iman)%typ, out_man(iman)%ind1,out_man(iman)%ind2, out_man(iman)%val1, out_man(iman)%val2
  enddo
!
! close input file
  close (9)
!
! open data file with scan variable info
  open (unit=1, file='Var.dat', form='formatted')
!
! read number of variables to scan, allocate scan variables
  read(1,11) junk,NumScanVar
  allocate(Label(NumScanVar))
  allocate(Range(NumScanVar))
  allocate(Start(NumScanVar))
  allocate(Stop(NumScanVar))
  allocate(PlotMin(NumScanVar))
  allocate(PlotMax(NumScanVar))
  allocate(Typical(NumScanVar))
  allocate(Value(NumScanVar))
!
! read variable information
  read(1,19) junk
  do ivar=1,NumScanVar
  read(1,18) num, Label(ivar), Start(ivar), Stop(ivar), PlotMin(ivar), PlotMax(ivar), Typical(ivar)

! assign typical values
  Value(ivar) = Typical(ivar)

! calculate range
  range(ivar) = Stop(ivar)-Start(ivar)
  enddo
!
! close data file with scan variable info
  close(1)
!
! standard formats for input
10    Format (a45)
11    Format (a45, I8)
12    Format (a45, l2)
!13    Format (a45, E15.8)
!14    Format (a45, a40)
15    Format (a45, f6.4)
16    Format (a45, a25)
17    Format (a45, a100)
18    format (i2,2x,a25,6(f10.4))
19    Format (a100)
20    Format (a45, a8)

! assign sums
  Btotal=0.

! print message on what you are scanning
  If(scantype=='cloud')      print*, 'scan mrt generic cloud with ', numpts, ' pts'
  If(scantype=='mrt_opt')    print*, 'scan mrt optimization with ', numpts, ' pts'
  If(scantype=='mrt_func')   print*, 'scan mrt function with ', numpts, ' pts'
  If(scantype=='mrt_chan')   print*, 'scan mrt channels with ', numpts, ' pts'
  If(scantype=='cosz')       print*, 'scan cosz with ', numpts, ' pts'
  If(scantype=='ext_coef')   print*, 'scan ext coef', numpts, ' pts'
  if(scantype=='L_band_VWC') print*, 'scan L-band VWC retrieval with', numpts, ' pts'
  If(scantype=='plank')      print*, 'scan planks law', numpts, ' pts'
  If(scantype=='Sun')        print*, 'scan Sunrise/Sunset ', numpts, ' pts'

! check variables
  if(flagDim == 2 .and. yvar /= 0) then
    print*, 'Error: for 2D scan, yvar must equal zero'
    print*, 'Dim: ', flagDim, 'Yvar: ', yvar
    stop
  endif

! allocate plotting variables
  max_nline=100
  allocate(X(max_nline,numpts))
  allocate(Y(max_nline,numpts))
  allocate(zval(max_nline,numpts,numpts))
  allocate(X3(numpts))
  allocate(Y3(numpts))
  X=0.
  Y=0.
  X3=0.
  Y3=0.
  zval=0.

! allocate some airmoss variables
  allocate(diel(2))
  allocate(len_lay(2))
  allocate(std_lay(2))
  allocate(DLayer(1))

! setup for DOTLRT
  do ifil=1,numfiles
    if (trim(files(ifil)%type) == 'atm_txt')  atm_txt_file = trim(files(ifil)%path)
    if (trim(files(ifil)%type) == 'atm_wrf')  atm_wrf_file = trim(files(ifil)%path)
    if (trim(files(ifil)%type) == 'out_path')  out_path = trim(files(ifil)%path)
  enddo
  if (trim(prof_src) == 'single') file_in = trim(atm_txt_file)
  if (trim(prof_src) == 'WRF') file_in = trim(atm_wrf_file)
  file_in = trim(atm_wrf_file)

  call setup_all_inputs( )
  nlev_ref = nlev
  call setup_all_outputs( )
  call extract_channel(chan_strt)
  atm_temp=atm
  atm_psuedo=atm

  nchannel = maxchan - minchan + 1
! Set up for observed bightness temperatures
  if (calc_psuedo) then
    print*, 'Create Psuedo Data'
    allocate(Tbo_obs(nchannel,npol))
    allocate(Tbo_sim(nchannel,npol))
    allocate(jacob(nchannel))
    allocate(res(nchannel))
    allocate(cost_chan(nchannel))
    Print*, 'Read FY data'
    call read_FY_swath()
    call calculate_atm_profile (atm_psuedo, psuedo_lay, n_psuedo, psuedo_man)
  endif

! Set up initial atmospheric profile
  if (calc_init) then
    print*, 'Create Initial Profile'
    call calculate_atm_profile (atm_init, init_lay, n_init, init_man)
    atm = atm_init
    cur_lay = init_lay

    call read_init_precip_lookup_table( )
  endif

! check consistancy on 2-D vs. 3-D plots
  if (flagDim==2) then
    maxi=numpts
    maxj=1
  else if (flagDim==3) then
    maxi=numpts
    maxj=numpts
  else
    print*, 'error: specify 2-D or 3-D'
    stop
  endif

! calculate X and Y increments, Dx and Dy
  call Plot2D3D('DeltaX  ')
  if (flagDim==3) call Plot2D3D('DeltaY  ')

  allocate(myvals(12))
  allocate(myarray(5,nchannel,numpts,numpts,12))
  allocate(tarray(numpts))
  allocate(parray(numpts))

  return
  end

!=======================================================================
  subroutine read_observed_brightness_temp ()
!=======================================================================
! reads in observed brightness temperatures
!
! Modifications:
!  4/8/2021  Kevin Schaefer created subroutine
!  5/21/2021 Kevin Schaefer removed print stat ment
!----------------------------------------------------------------------
  use scan_Variables
  use dotlrt_variables
  use profiles
  use dotlrt_output

  implicit none

! local variables
  integer ivar    ! (-) variable index
  integer ichan   ! (-) channel index
  integer ifil    ! (-) file index
  integer count   ! (-) channel index
  integer status  ! (-) read status variable
  integer num  ! (-) read status variable
  real(8) val  ! (-) input read variable
  real(8) obs(nchan) ! (k) full observed brightness temp
  character*250 filename  ! file name

! find observation file
  do ifil=1,numfiles
    if (trim(files(ifil)%type) == 'FY3_text')  filename = trim(files(ifil)%path)
  enddo
  print*, '    read Obs Tb: ', trim(filename)

! check number of channels
  open(unit=20, file=trim(filename), form='formatted', status='old')
  read(20,*) junk ! read header
  count=0
  do ichan=1,100
    read(20,*, iostat=status) num, val
    if(status<0) exit
    count=count+1
  enddo
  close(unit=20)
  if(count/=nchan) then
    print*, 'Error: wrong Tb observation file'
    print*, trim(filename)
    print*, count, nchan
    stop
  endif

! read in observations
  open(unit=20, file=trim(filename), form='formatted', status='old')
  read(20,*) junk ! read header
  do ichan=1,nchan
    read(20,*) num, val
    obs(ichan) = val
  enddo
  close(unit=20)

! transfer to observed brightness temperatures
  ivar=0
  do ichan=minchan, maxchan
    ivar=ivar+1
    Tbo_obs(ivar,:) = obs(ichan)
  enddo

  return
  end

!============================================================================
subroutine read_init_precip_lookup_table( )
!============================================================================
! This routine reads in an initial brightness temperature lookup table
!
! history:
!  6/29/2021 Kevin schaefer created routine
!----------------------------------------------------------------------------
  use scan_Variables
  use dotlrt_variables
  implicit none

! internal variables
  character*250 filename ! (-) local filename
  integer status  ! (-) read status
  integer ifil    ! (-) file index
  integer irec    ! (-) record index
  integer max     ! (-) max allowed number of records

! find init Tb filename
  do ifil=1,numfiles
    if (trim(files(ifil)%type) == 'init_precip')  filename = trim(files(ifil)%path)
  enddo

! count records
  open(unit=20, file=trim(filename), form='formatted', status='old')
  read(20,*, iostat=status) junk
  init_nlook=0
  max = 1000
  do irec=1,max
    read(20,*, iostat=status) junk
    if(status<0) exit
    init_nlook=init_nlook+1
  enddo
  close(unit=20)
  print*, '    Number lookup Table records:', init_nlook

! allocate lookup table
  allocate(init_look(init_nlook,nchan))
  allocate(init_dens(init_nlook))
  allocate(init_precip(nchannel))

! open init tb file
  open(unit=20, file=trim(filename), form='formatted', status='old')
  read(20,*, iostat=status) junk

! read in records
  do irec=1,init_nlook
    read(20,*) init_dens(irec),init_look(irec,:)
    if(dabs(init_dens(irec)) < 1.d-6) init_dens(irec) = 0.d0
    !print*, irec, init_dens(irec), init_precip(irec,1)
  enddo
!
! close lookup table file
  close(unit=20)

end subroutine read_init_precip_lookup_table

!============================================================================
subroutine calc_init_precip( )
!============================================================================
! This routine calculates initial total column precipitation from a lookup table
!
! history:
!  6/29/2021 Kevin schaefer created routine
!----------------------------------------------------------------------------
  use scan_Variables
  use dotlrt_variables
  implicit none

! internal variables
  integer ichan    ! (-) file index
  integer ivar    ! (-) record index
  integer irec    ! (-) record index
  real(8) val  ! (-) tb ratio
  real(8) max  ! (K) max tb in lookup table

! find nearest Tb
  ivar=0
  do ichan=minchan, maxchan
    ivar=ivar+1
    max = maxval(init_look(:,ichan))
    bad_point = .false.
    if(Tbo_obs(ivar,1) > max) then
      !print*, 'bad point ', max, Tbo_obs(ivar,1)
      bad_point = .true.
    else
      do irec=init_nlook, 1, -1
        if(init_look(irec,ichan) >= Tbo_obs(ivar,1)) then
          val=(Tbo_obs(ivar,1)-init_look(irec,ichan))/(init_look(irec+1,ichan)-init_look(irec,ichan))
          init_precip(ivar) = init_dens(irec)+(init_dens(irec+1)-init_dens(irec))*val
          !print*, irec
          !print*, init_look(irec+1,ichan),Tbo_obs(ivar,1),init_look(irec,ichan)
          !print*, init_dens(irec+1),init_precip(ivar), init_dens(irec)
          exit
        endif
      enddo
    endif
  enddo

end subroutine calc_init_precip
!==============================================================================
subroutine read_FY_swath()
!==============================================================================
! This program reads in FY data files
!
! History:
!  6/30/2021 Kevin Schaefer created routine
!------------------------------------------------------------------------------
  use scan_Variables
  use dotlrt_variables
  use profiles
  use netcdf

implicit none
!
! Local variables
  integer ncid       ! (-) netcdf file id
  integer groupid    ! (-) netcdf group id
  integer varid      ! (-) netcdf variable id
  integer numdim     ! (-) number of dimensions
  integer dimlen     ! (-) value of dimension
  integer status     ! (-) netcdf status code
  integer dimid(20)  ! (-) dimension IDs
  integer parent     ! (-) track parent id
  character(len = nf90_max_name) :: dimname
  character*50 subroutine ! subroutine name
  character*250 filename  ! file name
  integer idim       ! (-) dimension index
  integer ifil       ! (-) file index
  integer ivar       ! (-) variable index
  integer ichan      ! (-) channel index
!
! set subroutine name
  subroutine='read_FY_netcdf'

! find masked observation swath file
  do ifil=1,numfiles
    if (trim(files(ifil)%type) == 'FY3_vec')  filename = trim(files(ifil)%path)
  enddo
  print*, '    read Obs FY: ', trim(filename)

! Open swath file
  status = nf90_open(filename, NF90_NOWRITE, ncid)
  if (status /= nf90_noerr) then
    call handle_err(status, subroutine, 1)
  endif

! get Data group id
  status = nf90_inq_ncid(ncid, 'aoi_data', groupid)
  if (status /= nf90_noerr) call handle_err(status, subroutine, 2)

! get dimension ids of tb data
  parent=0
  status = nf90_inq_dimids (groupid, numdim, dimid, parent)
    if (status /= nf90_noerr) call handle_err(status, subroutine, 3)

! get tb swath dimensions
  do idim = 1,numdim
    status = nf90_inquire_dimension(groupid, dimid(idim), dimname, dimlen)
    if (status /= nf90_noerr) call handle_err(status, subroutine, 4)
    if (trim(dimname) == 'phony_dim_0') then
      swath_npts = dimlen
      print*, '      Number Points: ', dimlen
    endif
    if (trim(dimname) == 'phony_dim_1') then
      swath_nchan = dimlen
      print*, '      Number Channels: ', dimlen
    endif
    if (trim(dimname) == 'phony_dim_2') then
      swath_nlev = dimlen
      print*, '      Number Levels: ', dimlen
    endif
  enddo

! allocate swath data variables
  allocate(swath_tb_full(swath_nchan, swath_npts))
  allocate(swath_tb_chan(nchannel,swath_npts))
  allocate(swath_angle(swath_npts))
  allocate(swath_temp(swath_nlev,swath_npts))

! allocate optimization metrics
  allocate(swath_precip(swath_npts))
  allocate(swath_point(swath_npts))
  allocate(swath_iter(swath_npts))
  allocate(swath_cost(swath_npts))
  allocate(swath_tb_sim(swath_npts))
  allocate(swath_res(swath_npts))
  allocate(swath_ex(swath_npts))
  allocate(swath_bad(swath_npts))

! initialize metrics arrays to standard missing value
  swath_precip(:) = -999.d0
  swath_point(:)  = -999.d0
  swath_iter(:)   = -999.d0
  swath_cost(:)   = -999.d0
  swath_tb_sim(:) = -999.d0
  swath_res(:)    = -999.d0
  swath_ex(:)     = -999.d0
  swath_bad(:)    = -999.d0

! get brightness temperatures
  status = nf90_inq_varid(groupid, 'L1c_bt', varid)
  if (status /= nf90_noerr) call handle_err(status, subroutine, 5)
  status = nf90_get_var(groupid, varid, swath_tb_full)
  if (status /= nf90_noerr) call handle_err(status, subroutine, 6)

! get nadir angles
  status = nf90_inq_varid(groupid, 'L1c_inc_ang', varid)
  if (status /= nf90_noerr) call handle_err(status, subroutine, 7)
  status = nf90_get_var(groupid, varid, swath_angle)
  if (status /= nf90_noerr) call handle_err(status, subroutine, 8)

! get vertical temperature profiles
  status = nf90_inq_varid(groupid, 'temperature', varid)
  if (status /= nf90_noerr) call handle_err(status, subroutine, 9)
  status = nf90_get_var(groupid, varid, swath_temp)
  if (status /= nf90_noerr) call handle_err(status, subroutine, 10)

! Close File
  status = nf90_close(ncid  )
  if (status /= nf90_noerr) call handle_err(status, subroutine, 11)

! subset brightness temperatures by channel
  ivar=0
  do ichan=minchan, maxchan
    ivar=ivar+1
    swath_tb_chan(ivar,:) = swath_tb_full(ichan,:)
  enddo

end subroutine read_FY_swath

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

