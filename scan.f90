!=========================================================================
     program scan
!=========================================================================
! calculates output for various subroutines by scanning through ranges of input variables
!
! Modifications:
! 2/1/2000  Kevin Schaefer created scan
! 7/1/2001  Kevin Schaefer added generic scanning
! 1/17/2021 Kevin Schaefer split off from scan and restructures for MRT
! 1/30/2021 Kevin Schaefer got rid of biome loop and if_scan
! 4/4/2021  Kevin Schaefer moved write code to separate routine
!--------------------------------------------------------------------------

  use scan_Variables
  use dotlrt_variables
  use profiles
  use dotlrt_output
  use netcdf

  IMPLICIT NONE

  integer ix   ! x value index
  integer iy   ! y value index
  integer iseg ! segment index

! variables for timing
  integer clock_rate   ! conversion between count and clock time
  integer clock_start  ! time count at stat
  integer clock_stop   ! time count at end
  real(8) e_time       ! elapsed time

! allocate sib variable tree
  nsib=1
  call clock

! read inputs
  call scan_read_input

! start measuring execution time
  call system_clock(clock_start) ! start counting

  !call write_netcdf('myfilename.nc')
  open (unit=444, file='results.txt', status='unknown')
  write(444,*) '['

! set start values X variable
  call Plot2D3D('StartX  ')
  do ix=1,maxi
!
! set start values Y variables
    if (flagDim==3) call Plot2D3D('StartY  ')
    do iy=1,maxj

! call scan control
      call scan_control(ix,iy)

! assign plotting variables
      call Plot2D3D('AssXval ')
      X(1,ix)=Xval
      if(Flagdim==2) then
        do kkk=1,Numy
          X(kkk,ix)=Xval
          Y(kkk,ix)=test10(kkk)
        enddo
        print*, trim(Label(xvar)), value(xvar)
      endif
      if(Flagdim==3) then
        zval(1,ix,iy)=test
        zval(2,ix,iy)=testvar1
        x3(ix)=value(xvar)
        y3(iy)=value(yvar)
        !print*, trim(Label(xvar)), value(xvar), " :: ", trim(Label(yvar)), value(yvar)
        tarray(ix) = value(xvar)
        parray(iy) = value(yvar)
      endif

! increment Y variable
      if (flagDim==3) call Plot2D3D('IncrY   ')
    enddo ! end y-loop

! Increment X variable
    call Plot2D3D('IncrX   ')
  enddo ! end x-loop
  write(444,*) ']'
  write(*,*) 'FINISHED Program'

! write Output
  !call write_output
  !print*, 'nchannel = ', nchannel
  !print*, 'tarray = ', tarray
  !print*, 'parray = ', parray
  !print*, 'myarray = ', myarray
  !call write_netcdf('index_table_FY_50.nc')

! plot output
!  if (PlotFlag==1) then
!    if (flagDim==2) then
!      if(Yscale==1) then
!        Ymin=minval(Y)
!        Ymax=maxval(Y)
!      endif
!      call output_manipulation
!      call MultiLine (Numpts,X,Y,Numpts,X,Y, Xmin,Xmax,Ymin,Ymax, DevFlag,LabX,LabY,title,nchan, 1,Numy, Plotfile, legFlag, Legend)
!    else if (flagDim==3) then
!      call output_manipulation
!      call Contour_color (zval,numpts,numpts, PlotMin(Xvar),PlotMax(Xvar),PlotMin(yvar),PlotMax(yvar), Dx, Dy,  &
!      DevFlag, mincont, maxcont, nc, nchannel, Minchan,Maxchan,  &
!      LabX,LabY,Title,Plotfile, barformat)
!    endif
!  end if ! plot data
!
!! print execution times
!  if(print_ex) then
!    do iseg = 1, nseg
!      print*, seg_name(iseg), time_seg(iseg), num_call(iseg)
!    enddo
!  endif
!  call system_clock(clock_stop,clock_rate) ! stop counting
!  e_time = real(clock_stop-clock_start)/real(clock_rate)
!  print*, 'Overall Execution time (s):', e_time

!--------------------------------------------------
! INTERNAL SUBROUTINES
!---------------------------------------------------
    contains
!
!---------------------------------------------------
    subroutine clock
!---------------------------------------------------
! Call system clock again to determine how long code took to run. 
! The values of Tstart and Tend are in internal system time units. 
! The variable specifier count_rate provides the conversion from 
! internal system time units to seconds. 
!value(xvar), trim(Label(yvar)), value(yvar)
    call system_clock(count_rate=clock2)
    call system_clock(count=Tend)
    if(Tend >= Tstart) then
               Time=float(Tend-Tstart)/float(clock2)
    else
               Time=0.0
    endif
!
    end subroutine clock

    end program scan
!
!=======================================================================
    subroutine output_manipulation
!=========================================================================
! manipulates output from scan subroutine
!
! Modifications:
!  Kevin Schaefer made routine (6/8/14))
!--------------------------------------------------------------------------
!
    use scan_Variables
!
  IMPLICIT NONE
!
    real var1, var2 ! local temporary variables
    integer iman ! manipulation index
    integer ilin ! line number index
    integer iptx  ! x point number index
    integer ipty  ! y point number index
    integer ibio ! biome index
!
! output manipulations
    do iman=1,nout_man
!
! 2-D value Integral 
      if(out_man(iman)%doit.and.flagDim==2.and.out_man(iman)%typ=='sum') then
        print*, '2-D Integral sum'
        var2=(xmax-xmin)/real(numpts)
     do ilin=1,numy
       var1=0.
       do iptx=1,Numpts
         var1=var1+Y(ilin,iptx)
       enddo
       var1=var1*var2
       print*, trim(legend(ilin)),': ',var1
     enddo
      endif
!
! 3-D value Integral nv, MinBiome,MaxBiome
      if(out_man(iman)%doit.and.flagDim==3.and.out_man(iman)%typ=='sum') then
        print*, '3-D Integral sum'
          var1=0.
          do ipty=1,Numpts
            do iptx=1,Numpts
              var1=var1+zval(ibio,iptx,ipty)
            enddo
          enddo
          var1=var1*dx*dy
          print*, 'biome ',ibio,' sum: ',var1
      endif
!
! end manipulation loop
    enddo
!
    return
    end

!
!=======================================================================
    subroutine write_output
!=========================================================================
! writes output from scan program
!
! Modifications:
!  04/04/2021 Kevin Schaefer created routine
!--------------------------------------------------------------------------
  use scan_Variables
  implicit none

! Internal variables
  integer ix   ! x value index
  integer iy   ! y value index
  integer ipt  ! point index
  integer i3d  ! 3d index
  Character*100 form       ! writing format
  Character*500 write_line ! writing line
  Character*20 var_char    ! character version of variable

!--------------------------------------------------------------------------
! save line data as formatted text file
!--------------------------------------------------------------------------
  if(save_2d==1) then
    open(87, file=trim(plotfile)//'_'//trim(scantype)//'.dat', form='formatted')

! header
    write(form,'(i2)') numy
    form='(a15,1x,'//trim(form)//'(a15,1x))'
    Write(87,form) labx,legend(1:numy)

! main body
    write(form,'(i2)') numy
    form='(e15.8,1x,'//trim(form)//'(e15.8,1x))'
    do ipt=1,Numpts
      Write(87,form) x(1,ipt),Y(1:numy,ipt)
    enddo
    close(unit=87)
  endif

!--------------------------------------------------------------------------
! Save line data as CSV file
!--------------------------------------------------------------------------
  if(save_2d==2) then
    open(87, file=trim(plotfile)//'_'//trim(scantype)//'.csv', form='formatted')
	
! header
    write_line=trim(labx)
    do iy=1,numy
      write_line=trim(adjustl(write_line))//','//trim(adjustl(legend(iy)))
    enddo
    Write(87,*) adjustl(write_line)

! main body
    form='(e15.8,1x)'
    do ipt=1,Numpts
      write(var_char, '(e15.8)') x(1,ipt)
      write_line=trim(adjustl(var_char))
      do iy=1,numy
        write(var_char, '(e15.8)') Y(iy,ipt)
        write_line=trim(adjustl(write_line))//','//trim(adjustl(var_char))
      enddo
      Write(87,*) adjustl(write_line)
    enddo

    close(unit=87)
  endif

!--------------------------------------------------------------------------
! Save contour data as CSV file
!--------------------------------------------------------------------------
  if(save_3d==2) then
    do i3d = 1,2
      if(i3d == 1) open(87, file=trim(plotfile)//'_'//trim(scantype)//'_3d_01.csv', form='formatted')
      if(i3d == 2) open(87, file=trim(plotfile)//'_'//trim(scantype)//'_3d_02.csv', form='formatted')

! header
      write(var_char,'(i2)') Numpts
      write_line=trim(adjustl(var_char))
      form='(e15.8,1x)'
      do ix=1,Numpts
        write(var_char, '(e15.8)') x3(ix)
        write_line=trim(adjustl(write_line))//','//trim(adjustl(var_char))
      enddo
      Write(87,*) adjustl(write_line)

! main body
      form='(e15.8,1x)'
      do iy=1,Numpts
        write(var_char, '(e15.8)') y3(iy)
        write_line=trim(adjustl(var_char))
        do ix=1,Numpts
          write(var_char, '(e15.8)') zval(i3d,ix,iy)
          write_line=trim(adjustl(write_line))//','//trim(adjustl(var_char))
        enddo
        Write(87,*) adjustl(write_line)
      enddo

      close(unit=87)
    enddo
  endif


  return
  end

!
!=======================================================================
    subroutine write_netcdf(filename)
!=========================================================================
! writes output from scan program
!
! Modifications:
!  04/04/2021 Kevin Schaefer created routine
!--------------------------------------------------------------------------
  use scan_Variables
  use dotlrt_variables
  use profiles
  use dotlrt_output
  use netcdf
  implicit none

  character(len=*), intent(in) :: filename
  !real(8),          intent(in) :: myarray(:,:,:,:,:)
  !real(8),          intent(in) :: tarray(:)
  !real(8),          intent(in) :: parray(:)
  !integer,          intent(in) :: nchannel
  !integer,          intent(in) :: npts
  
  character (len = *), parameter :: HYDR_DIM_NAME = "hydro_type"
  character (len = *), parameter :: CHAN_DIM_NAME = "channel"
  character (len = *), parameter :: TEMP_DIM_NAME = "temperature"
  character (len = *), parameter :: DENS_DIM_NAME = "density"
  character (len = *), parameter ::  HAB_DIM_NAME = "cloud_absorption_coef"
  character (len = *), parameter ::  HSC_DIM_NAME = "cloud_scattering_coef"
  character (len = *), parameter ::    G_DIM_NAME = "cloud_asymytry_factor"
  character (len = *), parameter :: DHAB_DIM_NAME = "dcloud_absorption_coef"
  character (len = *), parameter :: DHSC_DIM_NAME = "dcloud_scattering_coef"
  character (len = *), parameter ::   DG_DIM_NAME = "dcloud_asymytry_factor"
  !hab     ! (-) cloud absorption coefficient
  !hsc     ! (-) cloud scattering coefficient
  !g       ! (-) cloud asymytry factor
  !dhab(3) ! (?) derivative absorption wrt to temp, k0, a0
  !dhsc(3) ! (?) derivative scattering wrt to temp, k0, a0
  !dg(3)   ! (?) derivative asymytry factor wrt to temp, k0, a0
  character (len = *), parameter :: VALS_DIM_NAME = "values"
  character (len = *), parameter :: UNITS  = "units"
  character (len = *), parameter :: HYDR_UNITS  = "hydro"
  character (len = *), parameter :: CHAN_UNITS  = "channel"
  character (len = *), parameter :: TEMP_UNITS  = "kelven"
  character (len = *), parameter :: DENS_UNITS  = "g/m^3"
  character (len = *), parameter ::  HAB_UNITS  = "g/m^3"
  character (len = *), parameter ::  HSC_UNITS  = "g/m^3"
  character (len = *), parameter ::    G_UNITS  = "g/m^3"
  character (len = *), parameter :: VALS_UNITS  = "some-unit"

  integer, parameter :: NDIMS = 3
  integer, parameter :: NHYDROS = 5

  integer :: NCHANNS

  !real :: array(NCHANNS, NHYDROS, npts, npts)

  integer :: dimids(NDIMS)

  integer :: ichan, ihydro, ipts, jpts, vpts, indx
  integer :: ncid 
  !integer :: chan_dimid
  !integer :: hydr_dimid
  integer :: temp_dimid
  integer :: dens_dimid
  integer :: vals_dimid
  !integer ::  hab_dimid
  !integer ::  hsc_dimid
  !integer ::    g_dimid
  !integer :: dhab_dimid
  !integer :: dhsc_dimid
  !integer ::   dg_dimid

  !integer :: chan_varid
  !integer :: hydr_varid
  integer :: temp_varid
  integer :: dens_varid
  integer :: vals_varid
  !integer ::  hab_varid
  !integer ::  hsc_varid
  !integer ::    g_varid
  !integer :: dhab_varid
  !integer :: dhsc_varid
  !integer ::   dg_varid


  character(NF90_MAX_NAME) :: var_name
  !character(NF90_MAX_NAME) :: vars(12)

  integer :: vars(12)
  integer :: npts
  integer, allocatable :: chan_varid(:)
  real(8), allocatable :: vals(:,:,:)
  !real(8), allocatable :: array_full(:,:,:,:,:)
  !real(8), allocatable :: array_t(:)
  !real(8), allocatable :: array_p(:)

  vars(1)  = 1
  vars(2)  = 2 
  vars(3)  = 3 
  vars(4)  = 4 
  vars(5)  = 5 
  vars(6)  = 6 
  vars(7)  = 7 
  vars(8)  = 8 
  vars(9)  = 9 
  vars(10) = 10
  vars(11) = 11
  vars(12) = 12
  !vars(1)  = 'hab'
  !vars(2)  = 'hsc'
  !vars(3)  = 'g'
  !vars(4)  = 'dhab_t'
  !vars(5)  = 'dhab_k0'
  !vars(6)  = 'dhab_a0'
  !vars(7)  = 'dhsc_t'
  !vars(8)  = 'dhsc_k0'
  !vars(9)  = 'dhsc_a0'
  !vars(10) = 'dg_t'
  !vars(11) = 'dg_k0'
  !vars(12) = 'dg_a0'

  NCHANNS = nchannel
  npts = numpts

  allocate(chan_varid(NCHANNS*NHYDROS))
  allocate(vals(npts,npts,12))
  !allocate(array_full(NHYDROS,NCHANNS,npts,npts,12))
  !allocate(array_t(npts))
  !allocate(array_p(npts))
  
  !do ichan=1,NCHANNS
  !  do ihydro=1,NHYDROS
  !    write(var_name, '(A6,A,I2,A,A7)') "chann", "_", ichan, "_", get_hydro_name(ihydro)
  !    print*, ichan, ihydro, trim(var_name)
  !  end do
  !end do
  !call exit(0)

  !print*, 'size(marray) = ', size(myarray)
  !do ichan=1,NCHANNS
  !  do ihydro=1,NHYDROS
  !    do ipts=1,npts
  !      do jpts=1,npts
  !        do vpts=1,12
  !          array_full(ihydro,ichan,ipts,jpts,vpts) = 1.0/vpts
  !          print*, 'ihydro, ichan, ipts, ipts, vpts, myarray', ihydro, ichan, ipts, jpts, vpts
  !          print*, 'myarray', myarray(ihydro,ichan,ipts,jpts,vpts)
  !        end do
  !      end do
  !    end do
  !  end do
  !end do

  !do ipts=1,npts
  !  print*, 'ipts = ', ipts
  !  array_t(ipts) = ipts
  !  array_p(ipts) = ipts
  !end do

  !print*, 'myarray = ', myarray
  !print*, 'tarray  = ', tarray
  !print*, 'parray  = ', parray

  ! Create the file. 
  !print*, 'create'
  call check( nf90_create(filename, nf90_clobber, ncid) , 'nf90_create')
  !print*, 'end create'

  ! Define the dimensions.
  !call check( nf90_def_dim(ncid, HYDR_DIM_NAME, NHYDROS,  hydr_dimid) )
  !call check( nf90_def_dim(ncid, CHAN_DIM_NAME, NCHANNS,  chan_dimid) )
  call check( nf90_def_dim(ncid, TEMP_DIM_NAME, npts, temp_dimid), 'nf90_def_dim')
  call check( nf90_def_dim(ncid, DENS_DIM_NAME, npts, dens_dimid), 'nf90_def_dim')
  call check( nf90_def_dim(ncid, VALS_DIM_NAME, 12,   vals_dimid), 'nf90_def_dim')
  !call check( nf90_def_dim(ncid,  HAB_DIM_NAME,  1,    hab_dimid), 'nf90_def_dim')
  !call check( nf90_def_dim(ncid,  HSC_DIM_NAME,  1,    hsc_dimid), 'nf90_def_dim')
  !call check( nf90_def_dim(ncid,    G_DIM_NAME,  1,      g_dimid), 'nf90_def_dim')
  !call check( nf90_def_dim(ncid, DHAB_DIM_NAME,  1,   dhab_dimid), 'nf90_def_dim')
  !call check( nf90_def_dim(ncid, DHSC_DIM_NAME,  1,   dhsc_dimid), 'nf90_def_dim')
  !call check( nf90_def_dim(ncid,   DG_DIM_NAME,  1,     dg_dimid), 'nf90_def_dim')


  ! Define the coordinate variables. They will hold the coordinate
  ! information, that is, the latitudes and longitudes. A varid is
  ! returned for each.
  ! call check( nf90_def_var(ncid, HYDR_DIM_NAME, NF90_REAL, hydr_dimid, hydr_varid) )
  ! call check( nf90_def_var(ncid, CHAN_DIM_NAME, NF90_REAL, chan_dimid, chan_varid) )
  call check( nf90_def_var(ncid, TEMP_DIM_NAME, NF90_REAL, temp_dimid, temp_varid) ,'nf90_def_var')
  call check( nf90_def_var(ncid, DENS_DIM_NAME, NF90_REAL, dens_dimid, dens_varid) ,'nf90_def_var')
  call check( nf90_def_var(ncid, VALS_DIM_NAME, NF90_REAL, vals_dimid, vals_varid) ,'nf90_def_var')
  !call check( nf90_def_var(ncid,  HAB_DIM_NAME, NF90_CHAR,  hab_dimid,  hab_varid) ,'nf90_def_var')
  !call check( nf90_def_var(ncid,  HSC_DIM_NAME, NF90_CHAR,  hsc_dimid,  hsc_varid) ,'nf90_def_var')
  !call check( nf90_def_var(ncid,    G_DIM_NAME, NF90_CHAR,    g_dimid,    g_varid) ,'nf90_def_var')
  !call check( nf90_def_var(ncid, DHAB_DIM_NAME, NF90_CHAR, dhab_dimid, dhab_varid) ,'nf90_def_var')
  !call check( nf90_def_var(ncid, DHSC_DIM_NAME, NF90_CHAR, dhsc_dimid, dhsc_varid) ,'nf90_def_var')
  !call check( nf90_def_var(ncid,   DG_DIM_NAME, NF90_CHAR,   dg_dimid,   dg_varid) ,'nf90_def_var')

  ! Assign units attributes to coordinate var data. This attaches a
  ! text attribute to each of the coordinate variables, containing the
  ! units.
  ! call check( nf90_put_att(ncid, hydr_varid, UNITS, HYDR_UNITS) )
  ! call check( nf90_put_att(ncid, chan_varid, UNITS, CHAN_UNITS) )
  call check( nf90_put_att(ncid, temp_varid, UNITS, TEMP_UNITS) ,'nf90_put_att')
  call check( nf90_put_att(ncid, dens_varid, UNITS, DENS_UNITS) ,'nf90_put_att')
  call check( nf90_put_att(ncid, vals_varid, UNITS, VALS_UNITS) ,'nf90_put_att')
  !call check( nf90_put_att(ncid,  hab_varid, UNITS,  HAB_UNITS) ,'nf90_put_att')
  !call check( nf90_put_att(ncid,  hsc_varid, UNITS,  HSC_UNITS) ,'nf90_put_att')
  !call check( nf90_put_att(ncid,    g_varid, UNITS,    G_UNITS) ,'nf90_put_att')
  !call check( nf90_put_att(ncid, dhab_varid, UNITS,  HAB_UNITS) ,'nf90_put_att')
  !call check( nf90_put_att(ncid, dhsc_varid, UNITS,  HSC_UNITS) ,'nf90_put_att')
  !call check( nf90_put_att(ncid,   dg_varid, UNITS,    G_UNITS) ,'nf90_put_att')

  ! Define the netCDF variables. The dimids array is used to pass the
  ! dimids of the dimensions of the netCDF variables.
  dimids = (/ temp_dimid, dens_dimid, vals_dimid /)
  do ichan=1,NCHANNS
    do ihydro=1,NHYDROS
      write(var_name, '(A5,"_",I2.2,"_",A)') "chann", ichan, trim(get_hydro_name(ihydro))
      indx = ihydro + NHYDROS*(ichan-1)
      !print*, 'indx = ', indx, ' , var_name = ', var_name
      call check( nf90_def_var(ncid, trim(var_name), NF90_REAL, dimids, chan_varid(indx) ) , 'nf90_def_var')
    end do
  end do

  ! Assign units attributes to the pressure and temperature netCDF
  ! variables.
  do ichan=1,NCHANNS
    do ihydro=1,NHYDROS
      indx = ihydro + NHYDROS*(ichan-1)
      !print*, 'indx = ', indx, ' ihydro = ', ihydro, ' ichan = ', ichan
      call check( nf90_put_att(ncid, chan_varid(indx), UNITS, VALS_UNITS) , 'nf90_put_att')
    end do
  end do


  ! End define mode.
  call check( nf90_enddef(ncid) , 'nf90_enddef')

  !print*, 'enddef'

  ! Write the coordinate variable data. This will put the latitudes
  ! and longitudes of our data grid into the netCDF file.
  !print*, 'tarray', tarray(:)
  call check( nf90_put_var(ncid, temp_varid, tarray) , 'nf90_put_var')
  !print*, 'parray', parray(:)
  call check( nf90_put_var(ncid, dens_varid, parray) , 'nf90_put_var')
  !print*, 'vars', vars(:)
  call check( nf90_put_var(ncid, vals_varid, vars(1:12)) , 'nf90_put_var')
  !call check( nf90_put_var(ncid,  hab_varid, 'hab'),  'nf90_put_var')
  !call check( nf90_put_var(ncid,  hsc_varid, 'hsc'),  'nf90_put_var')
  !call check( nf90_put_var(ncid,    g_varid, 'hab'),  'nf90_put_var')
  !call check( nf90_put_var(ncid, dhab_varid, 'dhab'),  'nf90_put_var')
  !call check( nf90_put_var(ncid, dhsc_varid, 'dhsc'),  'nf90_put_var')
  !call check( nf90_put_var(ncid,   dg_varid, 'dg'),  'nf90_put_var')

  ! Write the pretend data. This will write our surface pressure and
  ! surface temperature data. The arrays of data are the same size as
  ! the netCDF variables we have defined.
  do ichan=1,NCHANNS
    do ihydro=1,NHYDROS
      indx = ihydro + NHYDROS*(ichan-1)
      vals = myarray(ihydro, ichan,:,:,:)
      !write(*,*) size(vals)
      !write(*,*) size(array_full(ihydro, ichan,:,:,:))
      indx = ihydro + NHYDROS*(ichan-1)
      !print*, 'put_var', indx
      call check( nf90_put_var(ncid, chan_varid(indx), vals) , 'nf90_put_var')
    end do
  end do

  ! Close the file.
  call check( nf90_close(ncid) , 'nf90_close')

contains
  subroutine check(status, context)
    integer, intent ( in) :: status
    character(len=*), intent ( in) :: context
    print*, 'running - ', context 
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
  end subroutine check  

  function get_hydro_name(ihydro)
    integer, intent ( in) :: ihydro
    character(7) :: get_hydro_name
      
    select case (ihydro)
        case (1) 
           get_hydro_name = 'clw'
        case (2)
           get_hydro_name = 'rain'
        case (3) 
           get_hydro_name = 'snow'
        case (4) 
           get_hydro_name = 'ice'
        case (5) 
           get_hydro_name = 'graupel'
        case default
           get_hydro_name = 'invalid'
     end select
  end function get_hydro_name

end subroutine write_netcdf

! !READGRID - read a netCDF gridfile
! !:========================================================================
! subroutine read_index_table(infile,xpos,ypos,idata,NX,NY)
!   use scan_Variables
!   use dotlrt_variables
!   use profiles
!   use dotlrt_output
!   use netcdf
!   implicit none
! 
!   character(len=*), intent(in) :: infile
!   real(4), dimension(NX), intent(OUT) :: xpos
!   real(4), dimension(NY), intent(OUT) :: ypos
!   real(4), dimension(NX,NY), intent(OUT) :: idata
!   integer(4), intent(IN) :: NX, NY
!   integer(4), dimension(2) :: dimids
!   integer :: ncid, xtype, ndims, varid
!   character(len=*) :: xname, yname, vname
!   !Open netCDF file
!   !:-------:-------:-------:-------:-------:-------:-------:-------:
!   call check(nf90_open(infile, nf90_nowrite, ncid))
!   !Get the values of the coordinates and put them in xpos & ypos
!   !netcdf index_table_25 {
!   !dimensions:
!   !	temperature = 75 ;
!   !	density = 75 ;
!   !	values = 12 ;
!   !variables:
!   !	float temperature(temperature) ;
!   !		temperature:units = "kelven" ;
!   !	float density(density) ;
!   !		density:units = "g/m^3" ;
!   !	float values(values) ;
!   !		values:units = "some-unit" ;
!   !	float chann_01_clw(values, density, temperature) ;
!   !		chann_01_clw:units = "some-unit" ;
!   !	float chann_01_rain(values, density, temperature) ;
!   !		chann_01_rain:units = "some-unit" ;
!   !	float chann_01_snow(values, density, temperature) ;
!   !		chann_01_snow:units = "some-unit" ;
!   !	float chann_01_ice(values, density, temperature) ;
!   !		chann_01_ice:units = "some-unit" ;
!   !	float chann_01_graupel(values, density, temperature) ;
!   !		chann_01_graupel:units = "some-unit" ;
! 
!   !:-------:-------:-------:-------:-------:-------:-------:-------:
!   call check(nf90_inquire_variable(ncid,1,vname,xtype,ndims,dimids))
!   call check(nf90_inq_varid(ncid,vname,varid))
!   call check(nf90_get_var(ncid,varid,xpos))
!   call check(nf90_inquire_variable(ncid,2,vname,xtype,ndims,dimids))
!   call check(nf90_inq_varid(ncid,vname,varid))
!   call check(nf90_get_var(ncid,varid,ypos))
!   !Get the values of the perturbations and put them in idata
!   !:-------:-------:-------:-------:-------:-------:-------:-------:
!   call check(nf90_inquire_variable(ncid,3,vname,xtype,ndims,dimids))
!   call check(nf90_inq_varid(ncid,vname,varid))
!   call check(nf90_get_var(ncid,varid,idata))
!   !Close netCDF file
!   !:-------:-------:-------:-------:-------:-------:-------:-------:
!   call check(nf90_close(ncid))
! 
! contains
!   subroutine check(status, context)
!     integer, intent ( in) :: status
!     character(len=*), intent ( in) :: context
!     print*, 'running - ', context 
!     if(status /= nf90_noerr) then 
!       print *, trim(nf90_strerror(status))
!       stop "Stopped"
!     end if
!   end subroutine check  
! end subroutine read_index_table
!:=========================================================================

