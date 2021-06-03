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
  integer ivar ! (-) variable index
  integer ifil ! (-) file index
  integer iman ! (-) manipulation index
  integer num  ! (-) variable number
  integer i
  integer j
  integer k
  integer l
  integer m
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
  read (9,11) junk, scanFlag
  read (9,16) junk, scantype
!
! scanning options
  read (9,10) junk
  read (9,11) junk, flagDim
  read (9,11) junk, Xvar
  read (9,11) junk, YVar
  read (9,11) junk, NumY
  read (9,11) junk, numpts
  read (9,11) junk, minchan
  read (9,11) junk, maxchan
  read (9,11) junk, save_2d
  read (9,11) junk, save_3d
  read (9,11) junk, obs_calc
  read (9,11) junk, obs_num
  read (9,11) junk, obs_x
  read (9,11) junk, obs_y
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
  read (9,16) junk, Legend(1)
  read (9,16) junk, Legend(2)
  read (9,16) junk, Legend(3)
  read (9,16) junk, Legend(4)
  read (9,16) junk, Legend(5)
  read (9,16) junk, Legend(6)
  read (9,16) junk, Legend(7)
  read (9,16) junk, Legend(8)
  read (9,16) junk, Legend(9)
  read (9,16) junk, Legend(10)
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
  read(9,*),files(ifil)%num, files(ifil)%type,files(ifil)%path
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
!12    Format (a45, a2)
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

  If (scanflag==1) then

! print messages on what you are scanning
  If(scantype=='cloud')      print*, 'scan mrt generic cloud with ', numpts, ' pts'
  If(scantype=='mrt_func')   print*, 'scan mrt function with ', numpts, ' pts'
  If(scantype=='mrt_chan')   print*, 'scan mrt channels with ', numpts, ' pts'
  If(scantype=='cosz')       print*, 'scan cosz with ', numpts, ' pts'
  If(scantype=='ext_coef')   print*, 'scan ext coef', numpts, ' pts'
  if(scantype=='L_band_VWC') print*, 'scan L-band VWC retrieval with', numpts, ' pts'
  If(scantype=='plank')      print*, 'scan planks law', numpts, ' pts'
  If(scantype=='Sun')        print*, 'scan Sunrise/Sunset ', numpts, ' pts'

! allocate plotting variables
  max_NumY=10
  allocate(X(max_NumY,numpts))
  allocate(Y(max_NumY,numpts))
  allocate(zval(max_NumY,numpts,numpts))
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
  temp_atm=atm
  obs_atm=atm
  call setup_all_outputs( )
  call extract_channel(chan_strt)

! cost function setup
  nchannel = maxchan - minchan + 1
  allocate(Tbo_obs(nchannel,npol))
  allocate(Tbo_sim(nchannel,npol))
  allocate(cost_chan(nchannel))

  if (obs_calc == 1) call assign_obs_variables ()

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
  endif
  
! allocate index table variables
  allocate(myvals(12))
  !print*, '(nchannel, numpts) (', nchannel, ',', numpts, ')'
  allocate(myarray(5,nchannel,numpts,numpts,12))
  !do i = 1,5
  !do j = 1,nchannel
  !do k = 1,numpts
  !do l = 1,numpts
  !do m = 1,12
  !  myarray(i,j,k,l,m) = i+j+k+l+m
  !  print*, 'myarray(i,j,k,l,m)', i,j,k,l,m, myarray(i,j,k,l,m)
  !enddo
  !enddo
  !enddo
  !enddo
  !enddo
  allocate(tparray(numpts,numpts,12,5))
  allocate(tarray(numpts))
  allocate(parray(numpts))

  return
  end

