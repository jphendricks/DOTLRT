!    This is part of the netCDF package.
!     Copyright 2006 University Corporation for Atmospheric Research/Unidata.
!     See COPYRIGHT file for conditions of use.

!     This is a very simple example which writes a 2D array of
!     sample data. To handle this in netCDF we create two shared
!     dimensions, "x" and "y", and a netCDF variable, called "data".

!     This example demonstrates the netCDF Fortran 90 API. This is part
!     of the netCDF tutorial, which can be found at:
!     http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-tutorial
      
!     Full documentation of the netCDF Fortran 90 API can be found at:
!     http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-f90

!     $Id: simple_read.f90,v 1.7 2006/12/09 18:44:58 russ Exp $

program simple_read
  use netcdf
  use netcdf_utilities_mod

  implicit none

  ! This is the name of the data file we will create.
  character (len = *), parameter :: FILE_NAME = "/home/johe6518/index_tables/index_table_25.nc"

  call read_index_table(FILE_NAME)

contains
!READGRID - read a netCDF gridfile
!:========================================================================
function nc_open(infile) result(ncid)
  use netcdf
  implicit none
  
  character(len=*), intent(in) :: infile

  integer :: ncid

  call check(nf90_open(infile, nf90_nowrite, ncid), 'nc_open')

end function nc_open

!:========================================================================
subroutine nc_close(ncid)
  use netcdf
  implicit none
  
  integer, intent(in) :: ncid

  call check(nf90_close(ncid), 'nc_close')

end subroutine nc_close

!:========================================================================
subroutine nc_get_var(ncid, varname, vals)

  integer,                   intent(in)  :: ncid
  character(len=*),          intent(in)  :: varname
  real(8), dimension(12,75,75), intent(out) :: vals

  character(NF90_MAX_NAME)               :: vname
  character(NF90_MAX_NAME), dimension(3) :: dimnames
  integer, dimension(3)                  :: varsize
  integer, dimension(NF90_MAX_VAR_DIMS)  :: vardimids
  integer                                :: ndims, varid, status
  real(8), dimension(:,:,:), allocatable :: varvals
  integer :: vsize, dsize, tsize

  vname = varname

  allocate(varvals(12,75,75))

  vsize =  nc_get_dimension_size(ncid, 'values')
  dsize =  nc_get_dimension_size(ncid, 'density')
  tsize =  nc_get_dimension_size(ncid, 'temperature')

  print *, 'vsize  = ', vsize
  print *, 'dsize  = ', dsize
  print *, 'tsize  = ', tsize

  call nc_get_variable_size(ncid, varname, varsize)

  print *, 'varsize  = ', varsize(1), varsize(2), varsize(3)
  
  call nc_get_variable_dimension_names(ncid, varname, dimnames)
  print *, 'dimnames  = ', dimnames(1), dimnames(2), dimnames(3)
  call nc_get_variable(ncid, vname, varvals)

  vals = varvals 

  deallocate(varvals)

end subroutine

!:========================================================================
subroutine read_index_table(infile)
  use netcdf
  implicit none

  character(len=*), intent(in) :: infile
  !real(4), dimension(NX), intent(OUT) :: xpos
  !real(4), dimension(NY), intent(OUT) :: ypos
  !real(4), dimension(NX,NY), intent(OUT) :: idata
  !integer(4), intent(IN) :: NX, NY
  !integer(4), dimension(2) :: dimids
  integer :: ncid, ndims, varid
  integer :: ivals, jdens, ktemp
  real(8),    dimension(:,:,:), allocatable :: vals
  character(NF90_MAX_NAME), dimension(3) :: dimnames
  integer, dimension(3)                  :: varsize
  integer :: vsize, dsize, tsize
  character(NF90_MAX_NAME) :: varname
  !character(len=*) :: xname, yname, vname
  !Open netCDF file
  !:-------:-------:-------:-------:-------:-------:-------:-------:
  allocate(vals(75,75,12))

  ncid = nc_open_file_readonly(infile)


  !Get the values of the coordinates and put them in xpos & ypos
  !netcdf index_table_25 {
  !dimensions:
  !	temperature = 75 ;
  !	density = 75 ;
  !	values = 12 ;
  !variables:
  !	float temperature(temperature) ;
  !		temperature:units = "kelven" ;
  !	float density(density) ;
  !		density:units = "g/m^3" ;
  !	float values(values) ;
  !		values:units = "some-unit" ;
  !	float chann_01_clw(values, density, temperature) ;
  !		chann_01_clw:units = "some-unit" ;
  !	float chann_01_rain(values, density, temperature) ;
  !		chann_01_rain:units = "some-unit" ;
  !	float chann_01_snow(values, density, temperature) ;
  !		chann_01_snow:units = "some-unit" ;
  !	float chann_01_ice(values, density, temperature) ;
  !		chann_01_ice:units = "some-unit" ;
  !	float chann_01_graupel(values, density, temperature) ;
  !		chann_01_graupel:units = "some-unit" ;
  vsize =  nc_get_dimension_size(ncid, 'values')
  dsize =  nc_get_dimension_size(ncid, 'density')
  tsize =  nc_get_dimension_size(ncid, 'temperature')

  print *, 'vsize  = ', vsize
  print *, 'dsize  = ', dsize
  print *, 'tsize  = ', tsize

  varname='chann_01_snow'
  call nc_get_variable_size(ncid, varname, varsize)

  print *, 'varsize  = ', varsize(1), varsize(2), varsize(3)
  !call exit(0) 
  call nc_get_variable_dimension_names(ncid, varname, dimnames)
  print *, 'dimnames  = ', dimnames(1), dimnames(2), dimnames(3)
  call  nc_get_variable(ncid, varname, vals)
  do ktemp=1,75
    do jdens=1,75
      do ivals=1,12
        if(vals(ktemp,jdens,ivals) .gt. 0.0) then
           print *, 'vals(ktemp,jdens,ivals) = ', vals(ktemp,jdens,ivals), ivals, jdens, ktemp
        endif
      enddo
    enddo
  enddo

  !:-------:-------:-------:-------:-------:-------:-------:-------:
  !call check(nf90_inquire_variable(ncid,1,vname,xtype,ndims,dimids))
  !call check(nf90_inq_varid(ncid,vname,varid))
  !call check(nf90_get_var(ncid,varid,xpos))
  !call check(nf90_inquire_variable(ncid,2,vname,xtype,ndims,dimids))
  !call check(nf90_inq_varid(ncid,vname,varid))
  !call check(nf90_get_var(ncid,varid,ypos))
  !Get the values of the perturbations and put them in idata
  !:-------:-------:-------:-------:-------:-------:-------:-------:
  !call check(nf90_inquire_variable(ncid,3,vname,xtype,ndims,dimids))
  !call check(nf90_inq_varid(ncid,vname,varid))
  !call check(nf90_get_var(ncid,varid,idata))
  !Close netCDF file
  !:-------:-------:-------:-------:-------:-------:-------:-------:
  call nc_close(ncid)
  deallocate(vals)
end subroutine read_index_table

subroutine check(status, context)
  integer, intent ( in) :: status
  character(len=*), intent ( in) :: context
  print*, 'running - ', context 
  if(status /= nf90_noerr) then 
    print *, 'status = ', status
    print *, trim(nf90_strerror(status))
    stop "Stopped"
  end if
end subroutine check  
!:=========================================================================

end program simple_read
