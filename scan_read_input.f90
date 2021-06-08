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

! print messages on what you are scanning
  If(scantype=='cloud')      print*, 'scan mrt generic cloud with ', numpts, ' pts'
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

! do setup for DOTLRT
  do ifil=1,numfiles
    if (trim(files(ifil)%type) == 'atm_prof')  file_in = trim(files(ifil)%path)
    if (trim(files(ifil)%type) == 'out_path')  out_path = trim(files(ifil)%path)
  enddo
  call setup_all_inputs( )
  nlev_ref = nlev
  call setup_all_outputs( )
  call extract_channel(chan_strt)
  temp_atm=atm
  psuedo_atm=atm
  n_up_man = 8

! Set up for observed bightness temperatures
  if (calc_psuedo) then
    print*, 'Create Psuedo Data'
    nchannel = maxchan - minchan + 1
    allocate(Tbo_obs(nchannel,npol))
    allocate(Tbo_sim(nchannel,npol))
    allocate(cost_chan(nchannel))
    call calculate_atm_profile (psuedo_atm, psuedo_lay, n_psuedo, psuedo_man)
  endif

! Set up initial atmospheric profile
  if (calc_init) then
    print*, 'Create Initial Profile'
    call calculate_atm_profile (init_atm, init_lay, n_init, init_man)
    atm = init_atm
    cur_lay = init_lay
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

! find observation file
  do ifil=1,numfiles
    if (trim(files(ifil)%type) == 'obs_path')  obs_path = trim(files(ifil)%path)
  enddo
  print*, '    read Obs Tb: ', trim(obs_path)

! check number of channels
  open(unit=20, file=trim(obs_path), form='formatted', status='old')
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
    print*, trim(obs_path)
    print*, count, nchan
    stop
  endif

! read in observations
  open(unit=20, file=trim(obs_path), form='formatted', status='old')
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
