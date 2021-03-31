!
!===================================================================================
    subroutine create_qp2( time )
!===================================================================================
! creates a new qp2 diagnostic output file
!
! Modifications:
!   Kevin Schaefer changed to standard diagnostic variable nomenclature (10/22/08)
!   Kevin Schaefer added implicit none and use sib_io/const_module (10/22/08)
!   Kevin Schaefer removed long argument list of diagnostic control variables (10/22/08)
!   Kevin Schaefer changed to local instead of module netcdf id numbers (2/22/14)
!   Kevin Schaefer open and close file within routine instead of always open (2/22/14)
!   Kevin Schaefer added error handling checks (2/23/14)
!-----------------------------------------------------------------------------------
    #ifdef PGF
    use netcdf
    use typeSizes
    #endif
    use kinds
    use timetype
    use sib_const_module
    use sib_io_module
    
    implicit none
!
! input variables
    type(time_struct), intent(in) :: time
!
! local variables
    integer(kind=int_kind) :: status       ! return status of netcdf functions
    integer(kind=int_kind) did_lat      ! latitude dimension id #
    integer(kind=int_kind) did_lon      ! longitude dimension id #
    integer(kind=int_kind) did_time     ! time dimension id #
    integer(kind=int_kind) did_char     ! char_len dimension id #
    integer(kind=int_kind) did_subcount ! subcount variable id #
    integer(kind=int_kind) filid  ! netcdf file id number
    integer(kind=int_kind) varid  ! netcdf variable id number
    integer(kind=int_kind) ivar    ! variable index
    character(len=256) filename ! filename string
    character(len=20) subname    ! subroutine name for error checking
    character(len=40) units     ! variable units
    character(len=80) longname  ! variable description
    integer(kind=int_kind) unit_len ! lenth of units string
    integer(kind=int_kind) long_len ! length of long def string
!
! set subroutine name
    subname='create_qp2'
!
! create file
    if(qp_type=='mon') then
      if(flg_qpsubfold) then
        write( filename, '(a,i4.4,i2.2,a)' ) trim(out_path)//'qp2/hsib_', time%year, time%month, '.qp2.nc'
      else
        write( filename, '(a,i4.4,i2.2,a)' ) trim(out_path)//'hsib_', time%year, time%month, '.qp2.nc'
      endif
    elseif(qp_type=='ann') then
      if(flg_qpsubfold) then
        write( filename, '(a,i4.4,a)' ) trim(out_path)//'qp2/hsib_', time%year, '.qp2.nc'
      else
        write( filename, '(a,i4.4,a)' ) trim(out_path)//'hsib_', time%year, '.qp2.nc'
      endif
    else
      print*, 'Error: incorrect qp file format specified: ', trim(qp_type)
      stop
    endif
    status = nf90_create( trim(filename), nf90_clobber, filid)
    if(status/=nf90_noerr) then
      print*, trim(filename)
      call handle_err(status,subname,1)
    endif
!
! define global attributes
    call global_atts( filid, 'sib3', 'lat/lon', '1.0', drvr_type,  &
        biome_source, soil_source, soref_source, ndvi_source, c4_source,  &
        d13cresp_source )
!
! define dimensions
    status = nf90_def_dim( filid, 'time', nf90_unlimited, did_time )
    if(status/=nf90_noerr) call handle_err(status,subname,2)
    status = nf90_def_dim( filid, 'char_len', 10, did_char )
    if(status/=nf90_noerr) call handle_err(status,subname,3)
    status = nf90_def_dim( filid, 'latitude', jhr, did_lat )
    if(status/=nf90_noerr) call handle_err(status,subname,4)
    status = nf90_def_dim( filid, 'longitude', ihr, did_lon )
    if(status/=nf90_noerr) call handle_err(status,subname,5)
    status = nf90_def_dim( filid, 'subcount', subcount, did_subcount )
    if(status/=nf90_noerr) call handle_err(status,subname,6)
!
! define grid variables
    status = nf90_def_var( filid, 'latitude', nf90_float, (/did_lat/), varid )
    if(status/=nf90_noerr) call handle_err(status,subname,7)
    status = nf90_put_att( filid, varid, 'units', 'degrees_north' )
    if(status/=nf90_noerr) call handle_err(status,subname,8)
    status = nf90_put_att( filid, varid, 'quantity', 'latitude' )
    if(status/=nf90_noerr) call handle_err(status,subname,9)
    
    status = nf90_def_var( filid, 'longitude', nf90_float, (/did_lon/), varid )
    if(status/=nf90_noerr) call handle_err(status,subname,10)
    status = nf90_put_att( filid, varid, 'units', 'degrees_east' )
    if(status/=nf90_noerr) call handle_err(status,subname,11)
    status = nf90_put_att( filid, varid, 'quantity', 'longitude' )
    if(status/=nf90_noerr) call handle_err(status,subname,12)
    
    status = nf90_def_var( filid, 'lonindex', nf90_int, (/did_subcount/), varid )
    if(status/=nf90_noerr) call handle_err(status,subname,13)
    status = nf90_put_att( filid, varid, 'long_name', 'Longitude array index' )
    if(status/=nf90_noerr) call handle_err(status,subname,14)
    status = nf90_put_att( filid, varid, 'units', 'index-integer' )
    if(status/=nf90_noerr) call handle_err(status,subname,15)

    status = nf90_def_var( filid, 'latindex', nf90_int, (/did_subcount/), varid )
    if(status/=nf90_noerr) call handle_err(status,subname,16)
    status = nf90_put_att( filid, varid, 'long_name', 'Latitude array index' )
    if(status/=nf90_noerr) call handle_err(status,subname,17)
    status = nf90_put_att( filid, varid, 'units', 'index-integer' )
    if(status/=nf90_noerr) call handle_err(status,subname,18)
!
! define time variables
    status = nf90_def_var( filid, 'year', nf90_int, (/did_time/), varid )
    if(status/=nf90_noerr) call handle_err(status,subname,19)
    status = nf90_put_att( filid, varid, 'long_name', 'year at end of time averaging period' )
    if(status/=nf90_noerr) call handle_err(status,subname,20)
    status = nf90_put_att( filid, varid, 'units', 'year' )
    if(status/=nf90_noerr) call handle_err(status,subname,21)
    
    status = nf90_def_var( filid, 'month', nf90_int, (/did_time/), varid )
    if(status/=nf90_noerr) call handle_err(status,subname,22)
    status = nf90_put_att( filid, varid, 'long_name', 'month of year at end of time averaging period' )
    if(status/=nf90_noerr) call handle_err(status,subname,23)
    status = nf90_put_att( filid, varid, 'units', 'month' )
    if(status/=nf90_noerr) call handle_err(status,subname,24)
    
    status = nf90_def_var( filid, 'DOM', nf90_int, (/did_time/), varid )
    if(status/=nf90_noerr) call handle_err(status,subname,25)
    status = nf90_put_att( filid, varid, 'long_name', 'Day of Month (DOM) at end of time averaging period' )
    if(status/=nf90_noerr) call handle_err(status,subname,26)
    status = nf90_put_att( filid, varid, 'units', 'day' )
    if(status/=nf90_noerr) call handle_err(status,subname,27)
    
    status = nf90_def_var( filid, 'HOD', nf90_float, (/did_time/), varid )
    if(status/=nf90_noerr) call handle_err(status,subname,28)
    status = nf90_put_att( filid, varid, 'long_name', 'Hour of Day (HOD) at end of time averaging period' )
    if(status/=nf90_noerr) call handle_err(status,subname,29)
    status = nf90_put_att( filid, varid, 'units', 'hour' )
    if(status/=nf90_noerr) call handle_err(status,subname,30)

    status = nf90_def_var( filid, 'DOY', nf90_float, (/did_time/), varid )
    if(status/=nf90_noerr) call handle_err(status,subname,31)
    status = nf90_put_att( filid, varid, 'long_name', 'Day of Year (DOY) at end of time averaging period' )
    if(status/=nf90_noerr) call handle_err(status,subname,32)
    status = nf90_put_att( filid, varid, 'units', 'day' )
    if(status/=nf90_noerr) call handle_err(status,subname,33)
    
    status = nf90_def_var( filid, 'seconds', nf90_int, (/did_time/), varid )
    if(status/=nf90_noerr) call handle_err(status,subname,34)
    status = nf90_put_att( filid, varid, 'long_name', 'seconds since 0:00 GMT January first at end of time averaging period' )
    if(status/=nf90_noerr) call handle_err(status,subname,35)
    status = nf90_put_att( filid, varid, 'units', 'sec' )
    if(status/=nf90_noerr) call handle_err(status,subname,36)
!
! define diagnostic output variables
    do ivar = 1, nqp2var
        if ( doqp2(ivar) ) then
	    status = nf90_def_var( filid, trim(nameqp2(ivar)), nf90_float, (/did_subcount,did_time/), varid )
            if(status/=nf90_noerr) call handle_err(status,subname,37)
            call get_units( listqp2(ivar), longname, long_len, units, unit_len )
            status = nf90_put_att( filid, varid, 'long_name', trim(longname) )
            if(status/=nf90_noerr) call handle_err(status,subname,38)
            status = nf90_put_att( filid, varid, 'missing_value', 1.e36 )
            if(status/=nf90_noerr) call handle_err(status,subname,39)
        endif
    enddo
!
! switch from definition mode to data mode
    status = nf90_enddef( filid )
    if(status/=nf90_noerr) call handle_err(status,subname,40)
!
! write grid variables that do not vary with time
    status=nf90_inq_varid(filid, 'latitude', varid)
    if(status/=nf90_noerr) call handle_err(status,subname,41)
    status = nf90_put_var( filid, varid, latitude )
    if(status/=nf90_noerr) call handle_err(status,subname,42)

    status=nf90_inq_varid(filid, 'longitude', varid)
    if(status/=nf90_noerr) call handle_err(status,subname,43)
    status = nf90_put_var( filid, varid, longitude )
    if(status/=nf90_noerr) call handle_err(status,subname,44)

    status=nf90_inq_varid(filid, 'latindex', varid)
    if(status/=nf90_noerr) call handle_err(status,subname,45)
    status = nf90_put_var( filid, varid, sublat )
    if(status/=nf90_noerr) call handle_err(status,subname,46)

    status=nf90_inq_varid(filid, 'lonindex', varid)
    if(status/=nf90_noerr) call handle_err(status,subname,47)
    status = nf90_put_var( filid, varid, sublon )
    if(status/=nf90_noerr) call handle_err(status,subname,48)
!
! close file
    status = nf90_close(filid)
    if(status/=nf90_noerr) call handle_err(status,subname,49)

    end subroutine create_qp2
!
!===================================================================================
    subroutine write_qp2( time )
!===================================================================================
! writes time average output to an existing qp2 diagnostic output file
!
! Modifications:
!   Kevin Schaefer changed to standard diagnostic variable nomenclature (10/22/08)
!   Kevin Schaefer added implicit none and use sib_io/const_module (10/22/08)
!   Kevin Schaefer removed long argument list of diagnostic control variables (10/22/08)
!   Kevin Schaefer changed to local instead of module netcdf id numbers (2/22/14)
!   Kevin Schaefer open and close file within routine instead of always open (2/22/14)
!   Kevin Schaefer added error handling checks (2/23/14)
!-----------------------------------------------------------------------------------
    #ifdef PGF
    use netcdf
    use typeSizes
    #endif
    use kinds
    use timetype
    use sib_const_module
    use sib_io_module
    
    implicit none
!
! input variables
    type(time_struct), intent(in) :: time
!
! local variables
    integer status               ! return status of netcdf functions
    integer dimid                ! dimension id #
    integer(kind=int_kind) ivar  ! variable index
    integer(kind=int_kind) step  ! next time step in qp file
    character(len=10) name       ! temporary variable name
    character(len=20) subname    ! subroutine name for error checking
    character(len=256) filename  ! filename string
    integer(kind=int_kind) varid ! netcdf variable id number
    integer(kind=int_kind) filid ! netcdf file id number
!
! set subroutine name
    subname='write_qp2'
!
! open file
    if(qp_type=='mon') then
      if(flg_qpsubfold) then
        write( filename, '(a,i4.4,i2.2,a)' ) trim(out_path)//'qp2/hsib_', time%year, time%month, '.qp2.nc'
      else
        write( filename, '(a,i4.4,i2.2,a)' ) trim(out_path)//'hsib_', time%year, time%month, '.qp2.nc'
      endif
    elseif(qp_type=='ann') then
      if(flg_qpsubfold) then
        write( filename, '(a,i4.4,a)' ) trim(out_path)//'qp2/hsib_', time%year, '.qp2.nc'
      else
        write( filename, '(a,i4.4,a)' ) trim(out_path)//'hsib_', time%year, '.qp2.nc'
      endif
    endif
    status=nf90_open(trim(filename),nf90_write,filid)
    if(status/=nf90_noerr) call handle_err(status,trim(subname),1)
!
! find next time step
    status = nf90_inq_dimid( filid, 'time', dimid )
    if(status/=nf90_noerr) call handle_err(status,trim(subname),2)
    status = nf90_inquire_dimension( filid, dimid, name, step )
    if(status/=nf90_noerr) call handle_err(status,trim(subname),3)
    step = step + 1
!
! write time variables
    status=nf90_inq_varid(filid, 'year', varid)
    if(status/=nf90_noerr) call handle_err(status,trim(subname),4)
    status = nf90_put_var( filid, varid, time%year, (/step/) )
    if(status/=nf90_noerr) call handle_err(status,trim(subname),5)

    status=nf90_inq_varid(filid, 'month', varid)
    if(status/=nf90_noerr) call handle_err(status,trim(subname),6)
    status = nf90_put_var( filid, varid,  time%month, (/step/) )
    if(status/=nf90_noerr) call handle_err(status,trim(subname),7)

    status=nf90_inq_varid(filid, 'DOM', varid)
    if(status/=nf90_noerr) call handle_err(status,trim(subname),8)
    status = nf90_put_var( filid, varid,  time%day, (/step/) )
    if(status/=nf90_noerr) call handle_err(status,trim(subname),9)

    status=nf90_inq_varid(filid, 'DOY', varid)
    if(status/=nf90_noerr) call handle_err(status,trim(subname),10)
    status = nf90_put_var( filid, varid,  time%real_doy, (/step/) )
    if(status/=nf90_noerr) call handle_err(status,trim(subname),11)

    status=nf90_inq_varid(filid, 'HOD', varid)
    if(status/=nf90_noerr) call handle_err(status,trim(subname),12)
    status = nf90_put_var( filid, varid,  time%hour, (/step/) )
    if(status/=nf90_noerr) call handle_err(status,trim(subname),13)

    status=nf90_inq_varid(filid, 'seconds', varid)
    if(status/=nf90_noerr) call handle_err(status,trim(subname),14)
    status = nf90_put_var( filid, varid,  time%sec_year, (/step/) )
    if(status/=nf90_noerr) call handle_err(status,trim(subname),15)
!
! write diagnostiic variables
    do ivar = 1, nqp2var
      if ( doqp2(ivar) ) then
        status=nf90_inq_varid(filid, trim(nameqp2(ivar)), varid)
        if(status/=nf90_noerr) call handle_err(status,trim(subname),16)
        status = nf90_put_var(filid, varid, qp2sib(:,ivar),(/1,step/), (/subcount,1/) )
        if(status/=nf90_noerr) call handle_err(status,trim(subname),17)
      endif
    enddo
!
! close file
    status = nf90_close(filid)
    if(status/=nf90_noerr) call handle_err(status,trim(subname),18)

end subroutine write_qp2
!
!=======================================================================
    subroutine create_qp3( time )
!=======================================================================
! creates a new qp3 diagnostic output file for depth dependent time averages 
!
! Modifications:
!   Kevin Schaefer changed levels to actual soil depth (6/11/07)
!   Kevin Schaefer added implicit none and use sib_io_module (6/11/07)
!   Kevin Schaefer changed to standard diagnostic variable nomenclature (10/22/08)
!   Kevin Schaefer removed long argument list of diagnostic control variables (10/22/08)
!   Kevin Schaefer changed to local instead of module netcdf id numbers (2/22/14)
!   Kevin Schaefer open and close file within routine instead of always open (2/22/14)
!   Kevin Schaefer added error handling checks (2/23/14)
!-----------------------------------------------------------------------
!
    #ifdef PGF
    use netcdf
    use typeSizes
    #endif
    use kinds
    use timetype
    use sib_const_module
    use sib_io_module
!
implicit none
!
! input variables
    type(time_struct), intent(in) :: time
!
! local variables
    integer(kind=int_kind) status       ! return status of netcdf functions
    integer(kind=int_kind) did_lat      ! latitude dimension id #
    integer(kind=int_kind) did_lon      ! longitude dimension id #
    integer(kind=int_kind) did_time     ! time dimension id #
    integer(kind=int_kind) did_char     ! char_len dimension id #
    integer(kind=int_kind) did_nsoil    ! number of soil layer dimension id #
    integer(kind=int_kind) did_nsnow    ! number of soil layer dimension id #
    integer(kind=int_kind) did_npool    ! number of pools dimension id #
    integer(kind=int_kind) did_subcount ! subcount dimension id #
    integer(kind=int_kind) ivar  ! index variable
    character(len=40) units      ! variable units
    character(len=80) longname   ! variable description
    integer(kind=int_kind) unit_len ! lenth of units string
    integer(kind=int_kind) long_len ! length of long def string
    character(len=256) filename  ! file name
    character(len=20) subname    ! subroutine name for error checking
    integer(kind=int_kind) filid ! netcdf file id number
    integer(kind=int_kind) varid ! netcdf variable id number
!
! set subroutine name
    subname='create_qp3'
!
! create file
    if(qp_type=='mon') then
      if(flg_qpsubfold) then
        write( filename, '(a,i4.4,i2.2,a)' ) trim(out_path)//'qp3/hsib_',time%year, time%month, '.qp3.nc'
      else
        write( filename, '(a,i4.4,i2.2,a)' ) trim(out_path)//'hsib_',time%year, time%month, '.qp3.nc'
      endif
    elseif(qp_type=='ann') then
      if(flg_qpsubfold) then
        write( filename, '(a,i4.4,a)' ) trim(out_path)//'qp3/hsib_',time%year, '.qp3.nc'
      else
        write( filename, '(a,i4.4,a)' ) trim(out_path)//'hsib_',time%year, '.qp3.nc'
      endif
    else
      print*, 'Error: incorrect qp file format specified: ', trim(qp_type)
      stop
    endif
    status = nf90_create( trim(filename), nf90_clobber, filid)
    if(status/=nf90_noerr) call handle_err(status,subname,1)
!
! define global attributes
    call global_atts( filid, 'sib3', 'lat/lon', '1.0', drvr_type,  &
        biome_source, soil_source, soref_source, ndvi_source, c4_source,  &
        d13cresp_source )
!
! define dimensions
    status = nf90_def_dim( filid, 'time', nf90_unlimited, did_time )
    if(status/=nf90_noerr) call handle_err(status,subname,2)
    status = nf90_def_dim( filid, 'char_len', 10, did_char )
    if(status/=nf90_noerr) call handle_err(status,subname,3)
    status = nf90_def_dim( filid, 'latitude', jhr, did_lat )
    if(status/=nf90_noerr) call handle_err(status,subname,4)
    status = nf90_def_dim( filid, 'longitude', ihr, did_lon )
    if(status/=nf90_noerr) call handle_err(status,subname,5)
    status = nf90_def_dim( filid, 'nsoil', nsoil, did_nsoil )
    if(status/=nf90_noerr) call handle_err(status,subname,6)
    status = nf90_def_dim( filid, 'nsnow', nsnow, did_nsnow )
    if(status/=nf90_noerr) call handle_err(status,subname,7)
    status = nf90_def_dim( filid, 'npool', npool, did_npool )
    if(status/=nf90_noerr) call handle_err(status,subname,8)
    status = nf90_def_dim( filid, 'subcount', subcount, did_subcount )
    if(status/=nf90_noerr) call handle_err(status,subname,9)
!
! define grid variables    
    status = nf90_def_var( filid, 'latitude', nf90_float, (/did_lat/), varid )
    if(status/=nf90_noerr) call handle_err(status,subname,10)
    status = nf90_put_att( filid, varid, 'units', 'degrees_north' )
    if(status/=nf90_noerr) call handle_err(status,subname,11)
    status = nf90_put_att( filid, varid, 'quantity', 'latitude' )
    if(status/=nf90_noerr) call handle_err(status,subname,12)

    status = nf90_def_var( filid, 'longitude', nf90_float, (/did_lon/), varid )
    if(status/=nf90_noerr) call handle_err(status,subname,13)
    status = nf90_put_att( filid, varid, 'units', 'degrees_east' )
    if(status/=nf90_noerr) call handle_err(status,subname,14)
    status = nf90_put_att( filid, varid, 'quantity', 'longitude' )
    if(status/=nf90_noerr) call handle_err(status,subname,15)

    status = nf90_def_var( filid, 'lonindex', nf90_int, (/did_subcount/), varid )
    if(status/=nf90_noerr) call handle_err(status,subname,16)
    status = nf90_put_att( filid, varid, 'long_name',  'longitude index array' )
    if(status/=nf90_noerr) call handle_err(status,subname,17)
    status = nf90_put_att( filid, varid, 'units', 'index-integer' )
    if(status/=nf90_noerr) call handle_err(status,subname,18)

    status = nf90_def_var( filid, 'latindex', nf90_int, (/did_subcount/), varid )
    if(status/=nf90_noerr) call handle_err(status,subname,19)
    status = nf90_put_att( filid, varid, 'long_name', 'latitude index array' )
    if(status/=nf90_noerr) call handle_err(status,subname,20)
    status = nf90_put_att( filid, varid, 'units', 'index-integer' )
    if(status/=nf90_noerr) call handle_err(status,subname,21)
!
! define soil configuration
    status = nf90_def_var( filid, 'soil_layer_dz', nf90_float, (/did_nsoil/), varid )
    if(status/=nf90_noerr) call handle_err(status,subname,22)
    status = nf90_put_att( filid, varid, 'long_name', 'thickness of each soil layer' )
    if(status/=nf90_noerr) call handle_err(status,subname,23)
    status = nf90_put_att( filid, varid, 'units', 'm' )
    if(status/=nf90_noerr) call handle_err(status,subname,24)
    
    status = nf90_def_var( filid, 'soil_layer_top', nf90_float, (/did_nsoil/), varid )
    if(status/=nf90_noerr) call handle_err(status,subname,25)
    status = nf90_put_att( filid, varid, 'long_name', 'top of each soil layer' )
    if(status/=nf90_noerr) call handle_err(status,subname,26)
    status = nf90_put_att( filid, varid, 'units', 'm' )
    if(status/=nf90_noerr) call handle_err(status,subname,27)
    
    status = nf90_def_var( filid, 'soil_layer_mid', nf90_float, (/did_nsoil/), varid )
    if(status/=nf90_noerr) call handle_err(status,subname,28)
    status = nf90_put_att( filid, varid, 'long_name', 'middle of each soil layer' )
    if(status/=nf90_noerr) call handle_err(status,subname,29)
    status = nf90_put_att( filid, varid, 'units', 'm' )
    if(status/=nf90_noerr) call handle_err(status,subname,30)
    
    status = nf90_def_var( filid, 'soil_layer_bot', nf90_float, (/did_nsoil/), varid )
    if(status/=nf90_noerr) call handle_err(status,subname,31)
    status = nf90_put_att( filid, varid, 'long_name', 'bottom of each soil layer' )
    if(status/=nf90_noerr) call handle_err(status,subname,32)
    status = nf90_put_att( filid, varid, 'units', 'm' )
    if(status/=nf90_noerr) call handle_err(status,subname,33)

    status = nf90_def_var( filid, 'pool_name', nf90_char, (/did_char,did_npool/), varid )
    if(status/=nf90_noerr) call handle_err(status,subname,34)
    status = nf90_put_att( filid, varid, 'long_name', 'name of each biogeochemical pool' )
    if(status/=nf90_noerr) call handle_err(status,subname,35)
    status = nf90_put_att( filid, varid, 'units', '-' )
    if(status/=nf90_noerr) call handle_err(status,subname,36)
!
! define time variables
    status = nf90_def_var( filid, 'year', nf90_int, (/did_time/), varid )
    if(status/=nf90_noerr) call handle_err(status,subname,37)
    status = nf90_put_att( filid, varid, 'long_name', 'year at end of time averaging period' )
    if(status/=nf90_noerr) call handle_err(status,subname,38)
    status = nf90_put_att( filid, varid, 'units', 'year' )
    if(status/=nf90_noerr) call handle_err(status,subname,39)
    
    status = nf90_def_var( filid, 'month', nf90_int, (/did_time/), varid )
    if(status/=nf90_noerr) call handle_err(status,subname,40)
    status = nf90_put_att( filid, varid, 'long_name', 'month of year at end of time averaging period' )
    if(status/=nf90_noerr) call handle_err(status,subname,41)
    status = nf90_put_att( filid, varid, 'units', 'month' )
    if(status/=nf90_noerr) call handle_err(status,subname,42)
    
    status = nf90_def_var( filid, 'DOM', nf90_int, (/did_time/), varid )
    if(status/=nf90_noerr) call handle_err(status,subname,43)
    status = nf90_put_att( filid, varid, 'long_name', 'Day of Month (DOM) at end of time averaging period' )
    if(status/=nf90_noerr) call handle_err(status,subname,44)
    status = nf90_put_att( filid, varid, 'units', 'day' )
    if(status/=nf90_noerr) call handle_err(status,subname,45)
    
    status = nf90_def_var( filid, 'HOD', nf90_float, (/did_time/), varid )
    if(status/=nf90_noerr) call handle_err(status,subname,46)
    status = nf90_put_att( filid, varid, 'long_name', 'Hour of Day (HOD) at end of time averaging period' )
    if(status/=nf90_noerr) call handle_err(status,subname,47)
    status = nf90_put_att( filid, varid, 'units', 'hour' )
    if(status/=nf90_noerr) call handle_err(status,subname,48)

    status = nf90_def_var( filid, 'DOY', nf90_float, (/did_time/), varid )
    if(status/=nf90_noerr) call handle_err(status,subname,49)
    status = nf90_put_att( filid, varid, 'long_name', 'Day of Year (DOY) at end of time averaging period' )
    if(status/=nf90_noerr) call handle_err(status,subname,50)
    status = nf90_put_att( filid, varid, 'units', 'day' )
    if(status/=nf90_noerr) call handle_err(status,subname,51)
    
    status = nf90_def_var( filid, 'seconds', nf90_int, (/did_time/), varid )
    if(status/=nf90_noerr) call handle_err(status,subname,52)
    status = nf90_put_att( filid, varid, 'long_name', 'seconds since 0:00 GMT January first at end of time averaging period' )
    if(status/=nf90_noerr) call handle_err(status,subname,53)
    status = nf90_put_att( filid, varid, 'units', 'sec' )
    if(status/=nf90_noerr) call handle_err(status,subname,54)
!
! define diagnostic variables
    do ivar = 1, nqp3var
        if ( doqp3(ivar) ) then
	    if(dimqp3(ivar)==1) then
              status = nf90_def_var( filid, trim(nameqp3(ivar)),  nf90_float,  (/did_subcount,did_nsoil,did_time/), varid )
              if(status/=nf90_noerr) call handle_err(status,subname,55)

	    elseif(dimqp3(ivar)==2) then
              status = nf90_def_var( filid, trim(nameqp3(ivar)),  nf90_float,  (/did_subcount,did_nsnow,did_time/), varid )
              if(status/=nf90_noerr) call handle_err(status,subname,56)

	    elseif(dimqp3(ivar)==3) then
              status = nf90_def_var( filid, trim(nameqp3(ivar)),  nf90_float,  (/did_subcount,did_npool,did_time/), varid )
              if(status/=nf90_noerr) call handle_err(status,subname,57)
	    endif
            call get_units( listqp3(ivar), longname, long_len, units, unit_len )
            status = nf90_put_att( filid, varid, 'long_name', trim(longname) )
            if(status/=nf90_noerr) call handle_err(status,subname,58)
            status = nf90_put_att( filid, varid, 'units', trim(units) )
            if(status/=nf90_noerr) call handle_err(status,subname,59)
            status = nf90_put_att( filid, varid, 'missing_value', 1.e36 )
            if(status/=nf90_noerr) call handle_err(status,subname,60)
        endif
    enddo
!
! switch from definition mode to data mode
    status = nf90_enddef( filid )
    if(status/=nf90_noerr) call handle_err(status,subname,61)
!
! assign values to grid variables that do not vary with time
    status=nf90_inq_varid(filid, 'latitude', varid)
    if(status/=nf90_noerr) call handle_err(status,subname,62)
    status = nf90_put_var( filid, varid, latitude )
    if(status/=nf90_noerr) call handle_err(status,subname,63)

    status=nf90_inq_varid(filid, 'longitude', varid)
    if(status/=nf90_noerr) call handle_err(status,subname,64)
    status = nf90_put_var( filid, varid, longitude )
    if(status/=nf90_noerr) call handle_err(status,subname,65)

    status=nf90_inq_varid(filid, 'latindex', varid)
    if(status/=nf90_noerr) call handle_err(status,subname,66)
    status = nf90_put_var( filid, varid, sublat )
    if(status/=nf90_noerr) call handle_err(status,subname,67)

    status=nf90_inq_varid(filid, 'lonindex', varid)
    if(status/=nf90_noerr) call handle_err(status,subname,68)
    status = nf90_put_var( filid, varid, sublon )
    if(status/=nf90_noerr) call handle_err(status,subname,69)

    status=nf90_inq_varid(filid, 'soil_layer_dz', varid)
    if(status/=nf90_noerr) call handle_err(status,subname,71)
    status = nf90_put_var( filid, varid, soil_dz )
    if(status/=nf90_noerr) call handle_err(status,subname,72)

    status=nf90_inq_varid(filid, 'soil_layer_top', varid)
    if(status/=nf90_noerr) call handle_err(status,subname,73)
    status = nf90_put_var( filid, varid, soil_z_top )
    if(status/=nf90_noerr) call handle_err(status,subname,74)

    status=nf90_inq_varid(filid, 'soil_layer_mid', varid)
    if(status/=nf90_noerr) call handle_err(status,subname,75)
    status = nf90_put_var( filid, varid, soil_z_mid )
    if(status/=nf90_noerr) call handle_err(status,subname,76)

    status=nf90_inq_varid(filid, 'soil_layer_bot', varid)
    if(status/=nf90_noerr) call handle_err(status,subname,77)
    status = nf90_put_var( filid, varid, soil_z_bot )
    if(status/=nf90_noerr) call handle_err(status,subname,78)

    status=nf90_inq_varid(filid, 'pool_name', varid)
    if(status/=nf90_noerr) call handle_err(status,subname,79)
    status = nf90_put_var( filid, varid, pool_name )
    if(status/=nf90_noerr) call handle_err(status,subname,80)
!
! close qp3 file
    status = nf90_close(filid)
    if(status/=nf90_noerr) call handle_err(status,subname,81)

    end subroutine create_qp3
!
!===================================================================================
    subroutine write_qp3( time )
!===================================================================================
! writes depth dependent time averages to existing qp3 diagnostic output file
!
! Modifications:
!   Kevin Schaefer removed _date print statement (8/18/04)
!   Kevin Schaefer changed to standard diagnostic variable nomenclature (10/22/08)
!   Kevin Schaefer added implicit none ad use sib_io/const module(10/22/08)
!   Kevin Schaefer removed long argument list of diagnostic control variables (10/22/08)
!   Kevin Schaefer changed to local instead of module netcdf id numbers (2/22/14)
!   Kevin Schaefer open and close file within routine instead of always open (2/22/14)
!   Kevin Schaefer added error handling checks (2/23/14)
!-----------------------------------------------------------------------------------
    #ifdef PGF
    use netcdf
    use typeSizes
    #endif
    use kinds
    use timetype
    use sib_const_module
    use sib_io_module
    
    implicit none
!
! input variables
    type(time_struct), intent(in) :: time
!
! local variables
    integer(kind=int_kind) ivar   ! variable index
    integer(kind=int_kind) step   ! next time step in qp file
    character(len=10) name        ! temporary variable name
    integer(kind=int_kind) status ! return status of netcdf functions
    integer(kind=int_kind) dimid  ! dimension id #
    character(len=256) filename   ! filename string
    character(len=20) subname     ! subroutine name for error checking
    integer(kind=int_kind) varid  ! netcdf variable id number
    integer(kind=int_kind) filid  ! netcdf file id number
!
! set subroutine name
    subname='write_qp3'
!
! open file
    if(qp_type=='mon') then
      if(flg_qpsubfold) then
        write( filename, '(a,i4.4,i2.2,a)' ) trim(out_path)//'qp3/hsib_',time%year, time%month, '.qp3.nc'
      else
        write( filename, '(a,i4.4,i2.2,a)' ) trim(out_path)//'hsib_',time%year, time%month, '.qp3.nc'
      endif
    elseif(qp_type=='ann') then
      if(flg_qpsubfold) then
        write( filename, '(a,i4.4,a)' ) trim(out_path)//'qp3/hsib_',time%year, '.qp3.nc'
      else
        write( filename, '(a,i4.4,a)' ) trim(out_path)//'hsib_',time%year, '.qp3.nc'
      endif
    endif
    status=nf90_open(trim(filename),nf90_write,filid)
    if(status/=nf90_noerr) call handle_err(status,trim(subname),1)
!
! find next time step
    status = nf90_inq_dimid( filid, 'time', dimid )
    if(status/=nf90_noerr) call handle_err(status,trim(subname),2)
    status = nf90_inquire_dimension( filid, dimid, name, step )
    if(status/=nf90_noerr) call handle_err(status,trim(subname),3)
    step = step + 1
!
! write time variables
    status=nf90_inq_varid(filid, 'year', varid)
    if(status/=nf90_noerr) call handle_err(status,trim(subname),4)
    status = nf90_put_var( filid, varid, time%year, (/step/) )
    if(status/=nf90_noerr) call handle_err(status,trim(subname),5)

    status=nf90_inq_varid(filid, 'month', varid)
    if(status/=nf90_noerr) call handle_err(status,trim(subname),6)
    status = nf90_put_var( filid, varid,  time%month, (/step/) )
    if(status/=nf90_noerr) call handle_err(status,trim(subname),7)

    status=nf90_inq_varid(filid, 'DOM', varid)
    if(status/=nf90_noerr) call handle_err(status,trim(subname),8)
    status = nf90_put_var( filid, varid,  time%day, (/step/) )
    if(status/=nf90_noerr) call handle_err(status,trim(subname),9)

    status=nf90_inq_varid(filid, 'DOY', varid)
    if(status/=nf90_noerr) call handle_err(status,trim(subname),10)
    status = nf90_put_var( filid, varid,  time%real_doy, (/step/) )
    if(status/=nf90_noerr) call handle_err(status,trim(subname),11)

    status=nf90_inq_varid(filid, 'HOD', varid)
    if(status/=nf90_noerr) call handle_err(status,trim(subname),12)
    status = nf90_put_var( filid, varid,  time%hour, (/step/) )
    if(status/=nf90_noerr) call handle_err(status,trim(subname),13)

    status=nf90_inq_varid(filid, 'seconds', varid)
    if(status/=nf90_noerr) call handle_err(status,trim(subname),14)
    status = nf90_put_var( filid, varid,  time%sec_year, (/step/) )
    if(status/=nf90_noerr) call handle_err(status,trim(subname),15)
!
! write time mean diagnostics  trim(nameqp3(n))
    do ivar = 1, nqp3var
        if ( doqp3(ivar) ) then
            status=nf90_inq_varid(filid, trim(nameqp3(ivar)), varid)
            if(status/=nf90_noerr) call handle_err(status,trim(subname),16)

	    if(dimqp3(ivar)==1) then ! soil variables
              status = nf90_put_var( filid, varid,qp3_soil(:,:,ivar), (/1,1,step/),(/subcount,nsoil,1/) )
              if(status/=nf90_noerr) call handle_err(status,trim(subname),17)

	    elseif(dimqp3(ivar)==2) then ! snow variables
	      status = nf90_put_var( filid, varid,qp3_snow(:,:,ivar), (/1,1,step/),(/subcount,nsnow,1/) )
              if(status/=nf90_noerr) call handle_err(status,trim(subname),18)

	    elseif(dimqp3(ivar)==3) then ! pool variables
              status = nf90_put_var( filid, varid,qp3_pool(:,1:npool,ivar), (/1,1,step/),(/subcount,npool,1/) )
              if(status/=nf90_noerr) call handle_err(status,trim(subname),19)
	    endif
        endif
    enddo
!
! close qp3 file
    status = nf90_close(filid)
    if(status/=nf90_noerr) call handle_err(status,trim(subname),20)

end subroutine write_qp3
!
!===================================================================================
    subroutine global_atts (fileID, runname, grid, version, driver,  &
    biome_source, soil_source, soref_source, ndvi_source, c4_source,  &
    d13cresp_source )
!===================================================================================
! sets global attributes for all new diagnostic output files 
! netcdf file (fileID) must be in define mode.
!
! Modifications:
!   Kevin Schaefer added implicit none (10/22/08)
!-----------------------------------------------------------------------------------
    #ifdef PGF
    use netcdf
    use typeSizes
    #endif
    use sib_io_module
    
    implicit none
!
! input parameters
    integer, intent(in) :: fileID
    character(len=*), intent(in) :: runname
    character(len=*), intent(in) :: grid
    character(len=*), intent(in) :: version
    character(len=*), intent(in) :: driver
    character(len=100), intent(in) :: biome_source
    character(len=100), intent(in) :: soil_source
    character(len=100), intent(in) :: soref_source
    character(len=100), intent(in) :: ndvi_source
    character(len=100), intent(in) :: c4_source
    character(len=100), intent(in) :: d13cresp_source

! local variables
    integer :: status
    character(len=30) :: current_time
    character(len=8) :: t_date
    character(len=10) :: t_time
    character(len=5) :: zone
    integer, dimension(8) :: values
!
! set creation date and time
    call date_and_time(t_date, t_time, zone, values)
    current_time = t_date(5:6) // "/" // t_date(7:8) // "/" // t_date(1:4)   &
        //" at "// t_time(1:2) // ":" //t_time(3:4) // " "      &
        // zone // " GMT "
!
! add standard global attributes
    status = nf90_put_att ( fileID, nf90_global, 'calendar', 'noleap' )
    status = nf90_put_att ( fileID, nf90_global, 'institution', 'Colorado State University' )
    status = nf90_put_att ( fileID, nf90_global, 'history', 'Created: '//current_time )
    status = nf90_put_att ( fileID, nf90_global, 'run', runname )
    status = nf90_put_att ( fileID, nf90_global, 'grid', grid )
    status = nf90_put_att ( fileID, nf90_global, 'version', version )
    status = nf90_put_att ( fileID, nf90_global, 'Driver_Data', driver )
    status = nf90_put_att ( fileID, nf90_global, 'biome_source', trim(biome_source) )
    status = nf90_put_att ( fileID, nf90_global, 'soil_source', trim(soil_source) )
    status = nf90_put_att ( fileID, nf90_global, 'soref_source', trim(soref_source) )
    status = nf90_put_att ( fileID, nf90_global, 'ndvi_source', trim(ndvi_source) )
    status = nf90_put_att ( fileID, nf90_global, 'c4_source', trim(c4_source) )
    status = nf90_put_att ( fileID, nf90_global, 'd13cresp_source', trim(d13cresp_source) )
    status = nf90_put_att ( fileID, nf90_global, 'qp_type', trim(qp_type) )

    end subroutine global_atts
!
!===================================================================================
    subroutine get_units(description, longname, long_len, units, unit_len)
!===================================================================================
!  Purpose:
!   extracts a string enclosed within parentheses - used for 
!   units which are contained in a general description string.
!   returns the units string (units) and its length (unit_len),
!   the description string with the units removed (longname),
!   and its length (long_len).
!   note: embedded parentheses are ok
!
!  Variables:
!   Module parameters used:
!   MAXUNITS, MAXLONGN
!   
!  Bugs:
!   1) if the rightmost parenthesis is unmatched, units will be 
!      set to " " (one space) - this is to be interpreted as "none"
!   2) if a "(" is unmatched, it will be part of the returned
!      longname string
!   3) strings of only units (i.e. entire string enclosed in 
!      parentheses) do not work.
!
! Modifications:
!   Kevin Schaefer added implicit none (10/22/08)
!-----------------------------------------------------------------------------------
    implicit none
!
! input/output variables
    character(len=*), intent(in) :: description
    character(len=*), intent(out) :: units
    character(len=*), intent(out) :: longname
    integer, intent(out) :: unit_len, long_len
!
! local variables
    integer :: n, start_paren, end_paren, paren_count

    paren_count = 0
    start_paren = len_trim(description)
    end_paren = len_trim(description)

    do n = len(description), 1, -1
        if (description(n:n)==")") then
            if (paren_count == 0) then
                end_paren = n
            endif
            paren_count = paren_count + 1
        else if (description(n:n) == "(") then
            paren_count = paren_count - 1
            if (paren_count == 0) then
                start_paren = n
                exit
            endif
        end if
    end do

    !   in case of confusion, clear units and return unaltered description
    !   note: start_paren > end_paren should not be possible, but just in case...
    !         start_paren = end_paren occurs when there are no units
    !         start_paren = end_paren-1 occurs when units are "()"
    !   FIXME: n==1 is too limiting - what if I wanted only units?
    if (n == 1 .or. start_paren >= end_paren) then   ! no units
        units = " "
        unit_len = 1
        longname = trim(description)
        long_len = len_trim(longname)
    else if (start_paren == (end_paren-1)) then      ! "()" case
        units = " "
        unit_len = 1
        longname = trim(description(:start_paren-1))// &
            description(end_paren+1:)
        long_len = len_trim(longname)
    else                                             ! normal units
        units = description(start_paren+1:end_paren-1)
        unit_len = len_trim(units)
        longname = trim(description(:start_paren-1))// &
            description(end_paren+1:)
        long_len = len_trim(longname)
    end if

end subroutine get_units
