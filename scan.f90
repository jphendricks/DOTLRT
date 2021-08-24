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
  use netcdf_utilities_mod

  IMPLICIT NONE

  integer ix   ! x value index
  integer iy   ! y value index

  integer dsize
  integer tsize
  integer vsize

100 format(5(A,","))

  if (gen_index_table) then
     open(61,file='data-hab.txt',  action='write', status='replace')
     open(62,file='data-hsc.txt',  action='write', status='replace')
     open(63,file='data-g.txt',    action='write', status='replace')
     open(64,file='data-dhab.txt', action='write', status='replace')
     open(65,file='data-dhsc.txt', action='write', status='replace')
     open(66,file='data-dg.txt',   action='write', status='replace')

     write(61,100) 'error', 'vals', 'real', 'dens', 'temp'
     write(62,100) 'error', 'vals', 'real', 'dens', 'temp'
     write(63,100) 'error', 'vals', 'real', 'dens', 'temp'
     write(64,100) 'error', 'vals', 'real', 'dens', 'temp'
     write(65,100) 'error', 'vals', 'real', 'dens', 'temp'
     write(66,100) 'error', 'vals', 'real', 'dens', 'temp'
  endif

! allocate sib variable tree
  nsib=1
  call clock

! read inputs
  call scan_read_input

  print*, 'file_index_table = ', file_index_table
  ncltid = nc_open_file_readonly(file_index_table)
  print*, 'ncltid ' , ncltid

  if (use_index_table) then

     vsize =  nc_get_dimension_size(ncltid, 'values')
     dsize =  nc_get_dimension_size(ncltid, 'density')
     tsize =  nc_get_dimension_size(ncltid, 'temperature')

     if (dbg_index_table) then
         print*, 'tsize = ', tsize
         print*, 'dsize = ', dsize
         print*, 'vsize = ', vsize
     endif

     if (.not. allocated(dens_lt))      allocate(dens_lt(dsize))
     if (.not. allocated(temp_lt))      allocate(temp_lt(tsize))
     if (.not. allocated(lookup_table)) allocate(lookup_table(5,24,tsize,dsize,vsize))

     print*, "SIZE dens_lt", size(dens_lt)
     print*, "SIZE temp_lt", size(temp_lt)

     call nc_get_variable(ncltid, 'density',     dens_lt)
     call nc_get_variable(ncltid, 'temperature', temp_lt)

    print*, ncltid, lookup_table(1,1,1,1,1)
    call fill_index_table(ncltid, lookup_table)

  end if ! use_index_table

  if (.not. allocated(dens_lt))      allocate(dens_lt(numpts))
  if (.not. allocated(temp_lt))      allocate(temp_lt(numpts))
  if (.not. allocated(lookup_table)) allocate(lookup_table(5,24,numpts,numpts,12))

  !call write_netcdf('myfilename.nc')
  !open (unit=444, file='results.txt', status='unknown')
  !write(444,*) '['

! set start values X variable
  ! print*, '-----------------------------'
  ! print*, '  hab-     hydro_d(1),      hab,                      hydro_d(1)              (dens                       temp                  )'
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
      endif

! increment Y variable
      if (flagDim==3) call Plot2D3D('IncrY   ')
    enddo ! end y-loop

! Increment X variable
    call Plot2D3D('IncrX   ')
  enddo ! end x-loop
  !write(444,*) ']'

  if (gen_index_table) then
     close(61)
     close(62)
     close(63)
     close(64)
     close(65)
     close(66)
  endif

  write(*,*) 'FINISHED Program'

! write Output
  !call write_output
  !print*, 'nchannel = ', nchannel
  if (gen_index_table) then
     !call write_index_table('index_table_10_conv.nc')
     call write_index_table('index_table.nc')
  end if

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
    real(8) var1, var2 ! local temporary variables
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
     do ilin=1,nline
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
    write(form,'(i2)') nline
    form='(a15,1x,'//trim(form)//'(a15,1x))'
    Write(87,form) labx,legend(1:nline)

! main body
    write(form,'(i2)') nline
    form='(e15.8,1x,'//trim(form)//'(e15.8,1x))'
    do ipt=1,Numpts
      Write(87,form) x(1,ipt),Y(1:nline,ipt)
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
    do iy=1,nline
      write_line=trim(adjustl(write_line))//','//trim(adjustl(legend(iy)))
    enddo
    Write(87,*) adjustl(write_line)

! main body
    form='(e15.8,1x)'
    do ipt=1,Numpts
      write(var_char, '(e15.8)') x(1,ipt)
      write_line=trim(adjustl(var_char))
      do iy=1,nline
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
  if(flagDim==3 .and. save_3d==2) then
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
    subroutine write_index_table(filename)
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

  !character (len = *), parameter :: HYDR_DIM_NAME = "hydro_type"
  !character (len = *), parameter :: CHAN_DIM_NAME = "channel"
  character (len = *), parameter :: TEMP_DIM_NAME = "temperature"
  character (len = *), parameter :: DENS_DIM_NAME = "density"
  !character (len = *), parameter ::  HAB_DIM_NAME = "cloud_absorption_coef"
  !character (len = *), parameter ::  HSC_DIM_NAME = "cloud_scattering_coef"
  !character (len = *), parameter ::    G_DIM_NAME = "cloud_asymytry_factor"
  !character (len = *), parameter :: DHAB_DIM_NAME = "dcloud_absorption_coef"
  !character (len = *), parameter :: DHSC_DIM_NAME = "dcloud_scattering_coef"
  !character (len = *), parameter ::   DG_DIM_NAME = "dcloud_asymytry_factor"
  !hab     ! (-) cloud absorption coefficient
  !hsc     ! (-) cloud scattering coefficient
  !g       ! (-) cloud asymytry factor
  !dhab(3) ! (?) derivative absorption wrt to temp, k0, a0
  !dhsc(3) ! (?) derivative scattering wrt to temp, k0, a0
  !dg(3)   ! (?) derivative asymytry factor wrt to temp, k0, a0
  character (len = *), parameter :: VALS_DIM_NAME = "values"
  character (len = *), parameter :: UNITS  = "units"
  !character (len = *), parameter :: HYDR_UNITS  = "hydro"
  !character (len = *), parameter :: CHAN_UNITS  = "channel"
  character (len = *), parameter :: TEMP_UNITS  = "kelven"
  character (len = *), parameter :: DENS_UNITS  = "g/m^3"
  !character (len = *), parameter ::  HAB_UNITS  = "g/m^3"
  !character (len = *), parameter ::  HSC_UNITS  = "g/m^3"
  !character (len = *), parameter ::    G_UNITS  = "g/m^3"
  character (len = *), parameter :: VALS_UNITS  = "some-unit"

  integer, parameter :: NDIMS = 3
  integer, parameter :: NHYDROS = 5
  integer, parameter :: NVARIABLES = 12

  integer :: NCHANNS

  integer :: dimids(NDIMS)

  integer :: i, ichan, ihydro, indx
  integer :: myid

  integer :: temp_dimid
  integer :: dens_dimid
  integer :: vals_dimid

  integer :: temp_varid
  integer :: dens_varid
  integer :: vals_varid

  character(NF90_MAX_NAME) :: var_name

  integer, allocatable :: vars(:)
  integer :: npts
  integer, allocatable :: chan_varid(:)
  real(8), allocatable :: vals(:,:,:)

  NCHANNS = nchannel
  npts = numpts

  allocate(vars(NVARIABLES))
  do i = 1,NVARIABLES
     vars(i)  = i
  enddo

  allocate(chan_varid(NCHANNS*NHYDROS))
  allocate(vals(npts,npts,NHYDROS))

  ! Create the file.
  call check( nf90_create(filename, nf90_clobber, myid))

  ! Define the dimensions.
  call check( nf90_def_dim(myid, TEMP_DIM_NAME, npts,       temp_dimid))
  call check( nf90_def_dim(myid, DENS_DIM_NAME, npts,       dens_dimid))
  call check( nf90_def_dim(myid, VALS_DIM_NAME, NVARIABLES, vals_dimid))


  ! Define the coordinate variables. They will hold the coordinate
  ! information, that is, the latitudes and longitudes. A varid is
  ! returned for each.
  call check( nf90_def_var(myid, TEMP_DIM_NAME, NF90_DOUBLE, temp_dimid, temp_varid))
  call check( nf90_def_var(myid, DENS_DIM_NAME, NF90_DOUBLE, dens_dimid, dens_varid))
  call check( nf90_def_var(myid, VALS_DIM_NAME, NF90_DOUBLE, vals_dimid, vals_varid))

  ! Assign units attributes to coordinate var data. This attaches a
  ! text attribute to each of the coordinate variables, containing the
  ! units.
  call check( nf90_put_att(myid, temp_varid, UNITS, TEMP_UNITS))
  call check( nf90_put_att(myid, dens_varid, UNITS, DENS_UNITS))
  call check( nf90_put_att(myid, vals_varid, UNITS, VALS_UNITS))

  ! Define the netCDF variables. The dimids array is used to pass the
  ! dimids of the dimensions of the netCDF variables.
  dimids = (/ dens_dimid, temp_dimid, vals_dimid /)
  do ichan=1,NCHANNS
    do ihydro=1,NHYDROS
      write(var_name, '(A5,"_",I2.2,"_",A)') "chann", ichan, trim(get_hydro_name(ihydro))
      indx = ihydro + NHYDROS*(ichan-1)
      call check( nf90_def_var(myid, trim(var_name), NF90_DOUBLE, dimids, chan_varid(indx) ))
    end do
  end do

  ! Assign units attributes to the pressure and temperature netCDF
  ! variables.
  do ichan=1,NCHANNS
    do ihydro=1,NHYDROS
      indx = ihydro + NHYDROS*(ichan-1)
      call check( nf90_put_att(myid, chan_varid(indx), UNITS, VALS_UNITS))
    end do
  end do


  ! End define mode.
  call check( nf90_enddef(myid))

  ! Write the coordinate variable data. This will put the latitudes
  ! and longitudes of our data grid into the netCDF file.

  call check( nf90_put_var(myid, temp_varid, tarray))
  call check( nf90_put_var(myid, dens_varid, parray))
  call check( nf90_put_var(myid, vals_varid, vars(1:NVARIABLES)))

  ! Write the pretend data. This will write our surface pressure and
  ! surface temperature data. The arrays of data are the same size as
  ! the netCDF variables we have defined.
  do ichan=1,NCHANNS
    do ihydro=1,NHYDROS
      indx = ihydro + NHYDROS*(ichan-1)
      vals = myarray(ihydro, ichan,:,:,:)
      call check( nf90_put_var(myid, chan_varid(indx), vals))
    end do
  end do

  ! Close the file.
  call check( nf90_close(myid))

contains
  subroutine check(status, context)
    integer, intent ( in) :: status
    character(len=*), intent ( in), optional :: context
    if (present(context)) print*, 'running - ', context
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

end subroutine write_index_table
