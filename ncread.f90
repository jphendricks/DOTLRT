!     This is a very simple example which writes a 2D array of
!     sample data. To handle this in netCDF we create two shared
!     dimensions, "x" and "y", and a netCDF variable, called "data".

program simple_read
  use netcdf_utilities_mod
  use netcdf

  implicit none

  ! This is the name of the data file we will create.
  character (len = *), parameter :: FILE_NAME = "/home/johe6518/index_tables/index_table_25.nc"
  
  character (NF90_MAX_NAME) :: varname
  integer :: ichan, jchan, vsize, dsize, tsize, ncid
  real(8), dimension(:,:,:), allocatable :: vals

  real(8), parameter :: mfreq = 123.2
  real(8), parameter :: mtair = 300.1
  real(8), parameter :: mdens = 5.3
  integer, parameter :: mphase = 1

  !chann_01_clw
  !chann_01_rain
  !chann_01_snow
  !chann_01_ice
  !chann_01_graupel

  ncid = nc_open_file_readonly(FILE_NAME)
  !call scan_read_input

  call get_var_index_table(ncid,mfreq,mtair,mdens,mphase)

  vsize =  nc_get_dimension_size(ncid, 'values')
  dsize =  nc_get_dimension_size(ncid, 'density')
  tsize =  nc_get_dimension_size(ncid, 'temperature')

  allocate(vals(tsize,dsize,vsize))

  !do ichan = 1,2
  !    write(varname,2000) ichan,'clw'
  !    call get_var_index_table(ncid,varname,vals)
  !    print *, trim(varname), ' ', maxval(vals)
  !    write(varname,2000) ichan,'rain'
  !    call get_var_index_table(ncid,varname,vals)
  !    print *, trim(varname), ' ', maxval(vals)
  !    write(varname,2000) ichan,'snow'
  !    call get_var_index_table(ncid,varname,vals)
  !    print *, trim(varname), ' ', maxval(vals)
  !    write(varname,2000) ichan,'ice'
  !    call get_var_index_table(ncid,varname,vals)
  !    print *, trim(varname), ' ', maxval(vals)

  !write(varname,2000) 1,'graupel'
  !call get_var_index_table(ncid,varname,vals)
  !print *, trim(varname), ' ', maxval(vals)

  !enddo

  deallocate(vals)

  call nc_close_file(ncid)

  2000 format('chann_',i0.2,'_',A)

contains



end program simple_read
