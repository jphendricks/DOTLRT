!
!==============================================================================
subroutine write_atm_profile(filid, varname, nlon, nlat, nlev_max, atm_profile)
!==============================================================================
! This program writes 3-dimension atmospheric profiles to netcdf files
!
! history:
!   12/27/2020 Kevin Schaefer created routine
!-------------------------------------------------------------------------------
use netcdf
implicit none

! input variables
integer,intent(in) :: filid        ! (-) netcdf file id
character*80,intent(in) :: varname ! (-) variable name to write
integer,intent(in) :: nlat         ! (-) number of latitude points
integer,intent(in) :: nlon         ! (-) number of longitude points
integer,intent(in) :: nlev_max     ! (-) max number of vertical levels
real(8),intent(in) :: atm_profile(nlon, nlat, nlev_max) ! (variable) profile to write

! internal variables
character*50 subname ! (-) subroutine name
integer varid        ! (-) netcdf variable ID
integer status       ! (-) netcdf operation status flag

! set subroutine name
subname='write_atm_profile'

! get variable id
  status = nf90_inq_varid(filid, trim(varname), varid)
  if(status /= nf90_NoErr) then
    print*, trim(varname)
    call handle_err(status, subname, 1)
  endif

! write variable
  status = nf90_put_var(filid, varid, atm_profile)
  if(status/=nf90_noerr) call handle_err(status,trim(subname),2)
!
return
end subroutine write_atm_profile
!
!===================================================================================
  subroutine create_profile_file()
!===================================================================================
! creates a netcdf atmospheric profile file
!
! History:
!  12/27/2020 Kevn Schaefer created routine
!-----------------------------------------------------------------------------------
  use dotlrt_variables
  use dotlrt_output
  use profiles
  use netcdf
    
  implicit none
!
! local variables
  integer filid            ! (-) netcdf file id number
  integer status           ! (-) return status of netcdf functions
  integer did_lat          ! (-) latitude dimension id number
  integer did_lon          ! (-) longitude dimension id number
  integer did_lev          ! (-) level dimension id number
  integer ivar             ! (-) variable index
  integer varid(n_out_var) ! (-) variable id number
  character(250) filename  ! (-) filename string
  character(250) varname   ! (-) variable name string
  character*50 subname     ! (-) subroutine name

! print message
  print*, '    Create atmospheric profile netcdf output file'

! set subroutine name
  subname='create_profile_file'
!
! create file
  filename = trim(out_path)//'Profile.nc'
  if(flag_print_full) print*, '    ', trim(filename)

  status = nf90_create( trim(filename), cmode=or(nf90_clobber,nf90_64bit_offset), ncid=filid)
  if(status/=nf90_noerr) then
    print*, trim(filename)
    call handle_err(status,subname,1)
  endif
!
! define dimensions
    status = nf90_def_dim( filid, 'latitude', nlat, did_lat )
    if(status/=nf90_noerr) call handle_err(status,subname,2)

    status = nf90_def_dim( filid, 'longitude', nlon, did_lon )
    if(status/=nf90_noerr) call handle_err(status,subname,3)

    status = nf90_def_dim( filid, 'level', nlev_max, did_lev )
    if(status/=nf90_noerr) call handle_err(status,subname,6)

! define variables common to all files
  do ivar = 1, n_out_var
    if (var(ivar)%file =='All') then
      status = nf90_def_var( filid, var(ivar)%name, nf90_float, (/did_lon, did_lat/), varid(ivar))
      if(status/=nf90_noerr) call handle_err(status,subname,7)

      status = nf90_put_att( filid, varid(ivar), 'units', trim(var(ivar)%units) )
      if(status/=nf90_noerr) call handle_err(status,subname,8)

      status = nf90_put_att( filid, varid(ivar), 'long_name', trim(var(ivar)%long_name) )
      if(status/=nf90_noerr) call handle_err(status,subname,9)
    
      status = nf90_put_att( filid, varid(ivar), 'missing_value', var(ivar)%missing )
      if(status/=nf90_noerr) call handle_err(status,subname,10)
    endif
  enddo

! Define variables for profile file
  do ivar = 1, n_out_var
    if (var(ivar)%file =='Profile') then
      varname = trim(var(ivar)%name)
      status = nf90_def_var( filid, trim(varname), nf90_float, (/did_lon, did_lat, did_lev/), varid(ivar))
      if(status/=nf90_noerr) call handle_err(status,subname,11)

      status = nf90_put_att( filid, varid(ivar), 'units', trim(var(ivar)%units) )
      if(status/=nf90_noerr) call handle_err(status,subname,13)

      status = nf90_put_att( filid, varid(ivar), 'long_name', trim(var(ivar)%long_name) )
      if(status/=nf90_noerr) call handle_err(status,subname,14)
    
      status = nf90_put_att( filid, varid(ivar), 'missing_value', var(ivar)%missing )
      if(status/=nf90_noerr) call handle_err(status,subname,15)
    endif
  enddo

! switch from definition mode to data mode
  status = nf90_enddef( filid )
  if(status/=nf90_noerr) call handle_err(status,subname,18)

! write latitude
  varname='XLAT' 
  status = nf90_put_var(filid, varid(1), xlat)
  if(status/=nf90_noerr) call handle_err(status,trim(subname),19)

! write longitude
  varname='XLONG' 
  status = nf90_put_var(filid, varid(2), xlong)
  if(status/=nf90_noerr) call handle_err(status,trim(subname),20)

! write variables
  varname='height' 
  call write_atm_profile(filid, varname, nlon, nlat, nlev_max, height_mid)

! close file
  status = nf90_close(filid)
  if(status/=nf90_noerr) call handle_err(status,subname,21)

  end subroutine create_profile_file
!===================================================================================
  subroutine write_tb_channel(ichan)
!===================================================================================
! writes brightness temperature into existing output file
!
! History:
!  12/8/2020 Kevin schaefer split routine off from create_tb_file
!-----------------------------------------------------------------------------------
  use dotlrt_variables
  use dotlrt_output
  use profiles
  use netcdf
    
  implicit none

! inputs
  integer ichan  ! (-) channel index
!
! local variables
  integer filid           ! (-) netcdf file id number
  integer status          ! (-) return status of netcdf functions
  integer var_id          ! (-) variable id number
  character(250) filename ! (-) filename string
  character(250) varname  ! (-) variable name string
  character*50 subname    ! (-) subroutine name
  character*2 chan_txt    ! (-) text version of channel number

! print message
!  if(flag_print_full) print*, 'Write to Tb netcdf output file'

! set subroutine name
  subname='write_tb_channel'

! create file name
!  filename = trim(out_path)//'Tb.nc'
!  if(flag_print_full) print*, '    ', trim(filename)
  filename = trim(file_out)
  if(flag_print_full) print*, '    ', trim(filename) 
 
! open file
  status=nf90_open(trim(filename),nf90_write,filid)
  if(status/=nf90_noerr) then
    print*, trim(filename)
    call handle_err(status,subname,1)
  endif

! create variable name
  write( chan_txt, '(i2.2)' ) ichan
  varname = 'Tbo_stream'//'_'//chan_txt

! get variable ID
  status=nf90_inq_varid(filid, trim(varname), var_id)
  if(status/=nf90_noerr) call handle_err(status,trim(subname),16)

! write variable
  status = nf90_put_var(filid, var_id, Tbo_str_wrt)
  if(status/=nf90_noerr) call handle_err(status,trim(subname),18)

! close file
  status = nf90_close(filid)
  if(status/=nf90_noerr) call handle_err(status,subname,21)

  end subroutine write_tb_channel
!
!===================================================================================
  subroutine create_tb_file()
!===================================================================================
! creates a new brightness temperature diagnostic output file
!
! History:
!  10/20/2020 Kevn Schaefer created routine from create_qp2 from SiBCASA
!  12/8/2020 Kevin schaefer changed routine to brightness temp for each frequency
!-----------------------------------------------------------------------------------
  use dotlrt_variables
  use dotlrt_output
  use profiles
  use netcdf
    
  implicit none

! inputs
  integer ichan  ! (-) channel index
!
! local variables
  integer filid        ! (-) netcdf file id number
  integer status       ! (-) return status of netcdf functions
  integer did_lat      ! (-) latitude dimension id number
  integer did_lon      ! (-) longitude dimension id number
  integer did_stream   ! (-) stream dimension id number
  integer did_pol      ! (-) polarization dimension id number
  integer did_lev      ! (-) level dimension id number
  integer var_id        ! (-) variable id number
  integer ivar         ! (-) variable index
  integer varid(n_out_var) ! variable id number
  character(250) filename  ! (-) filename string
  character(250) varname   ! (-) variable name string
  character(250) text      ! (-) generic text string
  character*50 subname     ! (-) subroutine name
  character*2 chan_txt     ! (-) text version of channel number
  logical flag_var         ! (-) flag specifying variable type

! print message
  print*, 'Creat radiation netcdf output file'

! set subroutine name
  subname='create_tb_file'
!
! create file
!  filename = trim(out_path)//'Tb.nc'
  filename = trim(file_out)
  if(flag_print_full) print*, '    ', trim(filename)

  status = nf90_create( trim(filename), cmode=or(nf90_clobber,nf90_64bit_offset), ncid=filid)
  if(status/=nf90_noerr) then
    print*, trim(filename)
    call handle_err(status,subname,1)
  endif
!
! define dimensions
    status = nf90_def_dim( filid, 'latitude', nlat, did_lat )
    if(status/=nf90_noerr) call handle_err(status,subname,2)

    status = nf90_def_dim( filid, 'longitude', nlon, did_lon )
    if(status/=nf90_noerr) call handle_err(status,subname,3)

    status = nf90_def_dim( filid, 'stream', nstream_surf, did_stream )
    if(status/=nf90_noerr) call handle_err(status,subname,4)

    status = nf90_def_dim( filid, 'polarization', npol, did_pol )
    if(status/=nf90_noerr) call handle_err(status,subname,5)

    status = nf90_def_dim( filid, 'level', nlev, did_lev )
    if(status/=nf90_noerr) call handle_err(status,subname,6)

! define variables common to all files
  do ivar = 1, n_out_var
    ! check if variable applies to Tb file
    flag_var = .false.
    if (var(ivar)%file =='All') flag_var = .true.

    if(flag_var) then ! create a variable
      if(trim(var(ivar)%dim)=='(nlat,nlon)') then
        status = nf90_def_var( filid, var(ivar)%name, nf90_float, (/did_lon, did_lat/), varid(ivar))
        if(status/=nf90_noerr) call handle_err(status,subname,7)
      endif

      status = nf90_put_att( filid, varid(ivar), 'units', trim(var(ivar)%units) )
      if(status/=nf90_noerr) call handle_err(status,subname,8)

      status = nf90_put_att( filid, varid(ivar), 'long_name', trim(var(ivar)%long_name) )
      if(status/=nf90_noerr) call handle_err(status,subname,9)
    
      status = nf90_put_att( filid, varid(ivar), 'missing_value', var(ivar)%missing )
      if(status/=nf90_noerr) call handle_err(status,subname,10)
    endif
  enddo

! Define Frequency specific variables
  do ivar = 1, n_out_var
    ! check if variable applies to Tb file
    flag_var = .false.
    if (var(ivar)%file =='Tb') flag_var = .true.

    if(flag_var) then ! create a variable
      do ichan = 1, nchan
        write( chan_txt, '(i2.2)' ) ichan
        varname = trim(var(ivar)%name)//'_'//chan_txt
        if(trim(var(ivar)%dim)=='(nlat,nlon,npol)') then
          print*, trim(varname), filid
          status = nf90_def_var( filid, trim(varname), nf90_float, (/did_lon, did_lat, did_pol/), var_id)
          if(status/=nf90_noerr) call handle_err(status,subname,11)

        elseif(trim(var(ivar)%dim)=='(nlat,nlon,nstream,npol)') then
          status = nf90_def_var( filid, trim(varname), nf90_float, (/did_lon, did_lat, did_stream, did_pol/), var_id)
          if(status/=nf90_noerr) call handle_err(status,subname,12)
        endif

        status = nf90_put_att( filid, var_id, 'units', trim(var(ivar)%units) )
        if(status/=nf90_noerr) call handle_err(status,subname,13)

        text = trim(var(ivar)%long_name)//' for channel '//trim(chan_txt)
        status = nf90_put_att( filid, var_id, 'long_name', trim(text) )
        if(status/=nf90_noerr) call handle_err(status,subname,14)
    
        status = nf90_put_att( filid, var_id, 'missing_value', var(ivar)%missing )
        if(status/=nf90_noerr) call handle_err(status,subname,15)
	
        status = nf90_put_att( filid, var_id, 'frequency', instr_spec(ichan)%lo_freq )
        if(status/=nf90_noerr) call handle_err(status,subname,16)

        status = nf90_put_att( filid, var_id, 'instrument', instr_spec(ichan)%name )
        if(status/=nf90_noerr) call handle_err(status,subname,17)

      enddo
    endif
  enddo
!
! switch from definition mode to data mode
  status = nf90_enddef( filid )
  if(status/=nf90_noerr) call handle_err(status,subname,18)

! write output
  status = nf90_put_var(filid, varid(1), xlat)
  if(status/=nf90_noerr) call handle_err(status,trim(subname),19)

  status = nf90_put_var(filid, varid(2), xlong)
  if(status/=nf90_noerr) call handle_err(status,trim(subname),20)

! close file
  status = nf90_close(filid)
  if(status/=nf90_noerr) call handle_err(status,subname,21)

  end subroutine create_tb_file
!
!===================================================================================
  subroutine assign_to_output(ilon, ilat)
!===================================================================================
! sets up all variable specifications fornetcdf output
!
! History:
!  10/21/2020 Kevn Schaefer created routine
!-----------------------------------------------------------------------------------
  use dotlrt_variables
  use dotlrt_output
  use netcdf
    
  implicit none

! local variables
  integer ilon         ! (-) longitude index
  integer ilat         ! (-) latitude index

  Tbo_wrt(ilon,ilat,:) = Tbo_mat
  Tbo_str_wrt(ilon,ilat,:,:) = Tbo_str_mat
  dTb_dp_str_wrt(ilon,ilat,:,:,:) = dTb_dp_str_mat
  dTb_dT_str_wrt(ilon,ilat,:,:,:) = dTb_dT_str_mat
  dTb_dq_str_wrt(ilon,ilat,:,:,:) = dTb_dq_str_mat
  dTb_dc_str_wrt(ilon,ilat,:,:,:) = dTb_dw_str_mat(:,:,1,:)
  dTb_dr_str_wrt(ilon,ilat,:,:,:) = dTb_dw_str_mat(:,:,2,:)
  dTb_di_str_wrt(ilon,ilat,:,:,:) = dTb_dw_str_mat(:,:,3,:)
  dTb_ds_str_wrt(ilon,ilat,:,:,:) = dTb_dw_str_mat(:,:,4,:)
  dTb_dg_str_wrt(ilon,ilat,:,:,:) = dTb_dw_str_mat(:,:,5,:)

  end subroutine assign_to_output
!
!===================================================================================
  subroutine define_variables
!===================================================================================
! sets up all variable specifications for netcdf output
!
! History:
!  10/21/2020 Kevn Schaefer created routine
!-----------------------------------------------------------------------------------
  use dotlrt_variables
  use dotlrt_output
  use netcdf
    
  implicit none

! local variables
  integer ivar    ! variable index
  integer status  ! read status
  character*250 junk ! junk variable for reading

! print message
  if(flag_print_full) print*, '    Read variable definitions'

! open variable definition file
  open(unit=20, file=trim(file_var), form='formatted', status='old')
  read(20,*) junk
  
! count variables
  n_out_var=0
  do ivar=1,max_out_var
    read(20,*, iostat=status) junk
    if(status<0) exit
    n_out_var=n_out_var+1
  enddo
  if(flag_print_full) print*, '    Number output variables:', n_out_var

! close variable definition file
  close(unit=20)

! re-open variable definition file
  open(unit=20, file=trim(file_var), form='formatted', status='old')
  read(20,*) junk
  
! read variable defitions
  do ivar=1,n_out_var
    read(20,*, iostat=status) var(ivar)%file,var(ivar)%name, var(ivar)%units, var(ivar)%ndim, var(ivar)%dim, var(ivar)%long_name
  enddo
!
! close variable definition file
  close(unit=20)

! set standard stuff
  do ivar = 1, n_out_var
    var(ivar)%missing=missing
  enddo

  end subroutine define_variables
!
!===================================================================================
  subroutine create_jacobian_file(ichan)
!===================================================================================
! creates a new radiation diagnostic output file
!
! History:
!  10/20/2020 Kevn Schaefer created routine from create_qp2 from SiBCASA
!  1/27/2021  Kevin Schaefer changed to Jacobian 
!-----------------------------------------------------------------------------------
  use dotlrt_variables
  use dotlrt_output
  use profiles
  use netcdf
    
  implicit none

! inputs
  integer ichan  !(-) channel index
!
! local variables
  integer filid      ! netcdf file id number
  integer status     ! return status of netcdf functions
  integer did_lat    ! latitude dimension id number
  integer did_lon    ! longitude dimension id number
  integer did_stream ! stream dimension id number
  integer did_pol    ! polarization dimension id number
  integer did_lev    ! level dimension id number
  integer ivar       ! variable index
  integer varid(n_out_var) ! variable id number
  character(250) filename  ! filename string
  character*50 subname     ! subroutine name
  logical flag_var         ! (-) flag specifying variable type

! print message
  print*, 'Write Jacobian netcdf output file'

! set subroutine name
  subname='create_jac_file'
!
! create file
!  write( filename, '(a,i2.2,a)' ) trim(out_path)//'Jacobian_', ichan, '.nc'
  write( filename, '(a,i2.2,a)' ) trim(file_jac_prefix)//'_', ichan, '.nc'
  if(flag_print_full) print*, '    ', trim(filename)

  status = nf90_create( trim(filename), cmode=or(nf90_clobber,nf90_64bit_offset), ncid=filid)
  if(status/=nf90_noerr) then
    print*, trim(filename)
    call handle_err(status,subname,1)
  endif
!
! define dimensions
    status = nf90_def_dim( filid, 'latitude', nlat, did_lat )
    if(status/=nf90_noerr) call handle_err(status,subname,2)

    status = nf90_def_dim( filid, 'longitude', nlon, did_lon )
    if(status/=nf90_noerr) call handle_err(status,subname,3)

    status = nf90_def_dim( filid, 'stream', nstream_surf, did_stream )
    if(status/=nf90_noerr) call handle_err(status,subname,4)

    status = nf90_def_dim( filid, 'polarization', npol, did_pol )
    if(status/=nf90_noerr) call handle_err(status,subname,5)

    status = nf90_def_dim( filid, 'level', nlev, did_lev )
    if(status/=nf90_noerr) call handle_err(status,subname,6)

! define variables common to all files
  do ivar = 1, n_out_var
    ! check if variable applies to Tb file
    flag_var = .false.
    if (var(ivar)%file =='All') flag_var = .true.

    if(flag_var) then ! create a variable
      if(trim(var(ivar)%dim)=='(nlat,nlon)') then
        status = nf90_def_var( filid, var(ivar)%name, nf90_float, (/did_lon, did_lat/), varid(ivar))
        if(status/=nf90_noerr) call handle_err(status,subname,7)
      endif

      status = nf90_put_att( filid, varid(ivar), 'units', trim(var(ivar)%units) )
      if(status/=nf90_noerr) call handle_err(status,subname,8)

      status = nf90_put_att( filid, varid(ivar), 'long_name', trim(var(ivar)%long_name) )
      if(status/=nf90_noerr) call handle_err(status,subname,9)
    
      status = nf90_put_att( filid, varid(ivar), 'missing_value', var(ivar)%missing )
      if(status/=nf90_noerr) call handle_err(status,subname,10)
    endif
  enddo

! define jacobian variables
  do ivar = 1, n_out_var
    ! check if variable applies to Tb file
    flag_var = .false.
    if (var(ivar)%file =='Jacobian') flag_var = .true.

    if (flag_var) then
      if(trim(var(ivar)%dim)=='(nlat,nlon)') then
        status = nf90_def_var( filid, var(ivar)%name, nf90_float, (/did_lon, did_lat/), varid(ivar))
        if(status/=nf90_noerr) call handle_err(status,subname,7)

      elseif(trim(var(ivar)%dim)=='(nlat,nlon,npol)') then
        status = nf90_def_var( filid, var(ivar)%name, nf90_float, (/did_lon, did_lat, did_pol/), varid(ivar))
        if(status/=nf90_noerr) call handle_err(status,subname,8)

      elseif(trim(var(ivar)%dim)=='(nlat,nlon,nstream,npol)') then
        status = nf90_def_var( filid, var(ivar)%name, nf90_float, (/did_lon, did_lat, did_stream, did_pol/), varid(ivar))
        if(status/=nf90_noerr) call handle_err(status,subname,9)

      elseif(trim(var(ivar)%dim)=='(nlat,nlon,nlev,nstream,npol)') then
        status = nf90_def_var( filid, var(ivar)%name, nf90_float, (/did_lon, did_lat, did_lev, did_stream, did_pol/), varid(ivar))
        if(status/=nf90_noerr) call handle_err(status,subname,10)
    endif

    status = nf90_put_att( filid, varid(ivar), 'units', trim(var(ivar)%units) )
    if(status/=nf90_noerr) call handle_err(status,subname,11)

    status = nf90_put_att( filid, varid(ivar), 'long_name', trim(var(ivar)%long_name) )
    if(status/=nf90_noerr) call handle_err(status,subname,12)
    
    status = nf90_put_att( filid, varid(ivar), 'missing_value', var(ivar)%missing )
    if(status/=nf90_noerr) call handle_err(status,subname,13)
    endif
  enddo
!
! switch from definition mode to data mode
 status = nf90_enddef( filid )
 if(status/=nf90_noerr) call handle_err(status,subname,14)

! write output
! beware: If you change the order in the var def file, you must change these indeces
  status = nf90_put_var(filid, varid(1), xlat)
  if(status/=nf90_noerr) call handle_err(status,trim(subname),15)

  status = nf90_put_var(filid, varid(2), xlong)
  if(status/=nf90_noerr) call handle_err(status,trim(subname),16)

  !status = nf90_put_var(filid, varid(5), dTb_dp_str_wrt)
  if(status/=nf90_noerr) call handle_err(status,trim(subname),19)

  !status = nf90_put_var(filid, varid(6), dTb_dT_str_wrt)
  if(status/=nf90_noerr) call handle_err(status,trim(subname),20)

  !status = nf90_put_var(filid, varid(7), dTb_dq_str_wrt)
  if(status/=nf90_noerr) call handle_err(status,trim(subname),21)

  !status = nf90_put_var(filid, varid(8), dTb_dc_str_wrt)
  if(status/=nf90_noerr) call handle_err(status,trim(subname),22)

  !status = nf90_put_var(filid, varid(8), dTb_dr_str_wrt)
  if(status/=nf90_noerr) call handle_err(status,trim(subname),23)

  !status = nf90_put_var(filid, varid(10), dTb_di_str_wrt)
  if(status/=nf90_noerr) call handle_err(status,trim(subname),24)

  !status = nf90_put_var(filid, varid(11), dTb_ds_str_wrt)
  if(status/=nf90_noerr) call handle_err(status,trim(subname),25)

  !status = nf90_put_var(filid, varid(12), dTb_dg_str_wrt)
  if(status/=nf90_noerr) call handle_err(status,trim(subname),26)

! close file
  status = nf90_close(filid)
  if(status/=nf90_noerr) call handle_err(status,subname,27)

  end subroutine create_jacobian_file
!
!====================================================================
  subroutine write_text_profile()
!====================================================================
! writes a single text profile from WRF run
!
! History:
!  10/20/2020  Kevin Schaefer created routine
!--------------------------------------------------------------------
  use dotlrt_variables
  use profiles
  implicit none
!
! internal variables  
  character*250 filename ! (-) output filename
  character*250 line     ! (-) output line
  character*250 fmt      ! (-) output format spec
  character*25 txt_val1  ! (-) text version of value
  character*25 txt_val2  ! (-) text version of value
  integer ilev           ! (-) level index
  integer ivar           ! (-) value index
  character*25 values(nvar_prof)! (varies) temporary write variable
!
! profile values
! values(1) (km) height
! values(2) (mb) atmospheric pressure
! values(3) (K) atmospheric temperature
! values(4) (g/m^3) water vapor density
! values(5) (g/m^3) cloud liquid water density
! values(6) (g/m^3) rain density
! values(7) (g/m^3) ice density
! values(8) (g/m^3) snow density
! values(9) (g/m^3) graupel density

! print message
  print*, 'Write Single Text Atmospheric Profile'

! open file
  filename=trim(out_path)//'Atm_prof'
  fmt='(i3.3)'
  write(txt_val1, fmt) save_ilon
  write(txt_val2, fmt) save_ilat
  filename=trim(filename)//'_'//trim(txt_val1)//'_'//trim(txt_val2)//'.txt'
  print*, '    ', trim(filename)
!
! open profile file
  open(unit=20, file=trim(filename), form='formatted')

! write number of levels
  fmt='(i3.0)'
  write(txt_val1, fmt) nlev
  fmt='(a3)'  
  write(20,fmt) trim(adjustl(txt_val1))

! convert level values to text and write to file
  do ilev = 1, nlev
    fmt = '(f15.8)'
    write(values(1), fmt) height_mid(save_ilon, save_ilat, ilev) ! (km)
    write(values(2), fmt) press(save_ilon, save_ilat, ilev)  ! (mb)
    write(values(3), fmt) temp(save_ilon, save_ilat, ilev)   ! (K)
    fmt = '(e15.8)'
    write(values(4), fmt) QVAPOR(save_ilon, save_ilat, ilev) ! (g/m^3)
    write(values(5), fmt) QCLOUD(save_ilon, save_ilat, ilev) ! (g/m^3)
    write(values(6), fmt) QRAIN(save_ilon, save_ilat, ilev)  ! (g/m^3)
    write(values(7), fmt) QICE(save_ilon, save_ilat, ilev)   ! (g/m^3)
    write(values(8), fmt) QSNOW(save_ilon, save_ilat, ilev)  ! (g/m^3)
    write(values(9), fmt) QGRAUP(save_ilon, save_ilat, ilev) ! (g/m^3)

    do ivar =1, nvar_prof
      if(ivar==1) then
        line=trim(values(ivar))
      else
        line=trim(line)//','//trim(values(ivar))
      endif
    enddo
    write(20,*) trim(adjustl(line))
  end do

! close profile file
  close(unit=20)

  end subroutine write_text_profile
!
!====================================================================
  subroutine write_radiation_profile(ilon, ilat)
!====================================================================
! writes a single text file with TOA radiation output
!
! History:
!  10/21/2020  Kevin Schaefer created routine
!--------------------------------------------------------------------
  use dotlrt_variables
  use dotlrt_output
  use profiles
  implicit none

! inputs
  integer ilon         ! (-) longitude index
  integer ilat         ! (-) latitude index

!
! internal variables  
  logical flag_write     ! (-) flag to write the file
  character*250 filename ! (-) output filename
  character*250 fmt      ! (-) output format spec
  character*25 txt_val1  ! (-) text version of value
  character*25 txt_val2  ! (-) text version of value
  character*25 txt_val3  ! (-) text version of value
  integer iang           ! (-) level index

! match grid point
  flag_write = .false.
  if (ilon==save_ilon .and. ilat==save_ilat) flag_write = .true.
  if (.not.flag_write) return

! print message
  print*, 'Write Single Text radiation Profile'

! open file
  filename=trim(out_path)//'rad_prof'
  fmt='(i3.3)'
  write(txt_val1, fmt) save_ilon
  write(txt_val2, fmt) save_ilat
  fmt='(i2.2)'
  write(txt_val3, fmt) nstream
  filename=trim(filename)//'_'//trim(txt_val1)//'_'//trim(txt_val2)//'_'//trim(txt_val3)//'.txt'
  print*, '    ', trim(filename)
!
! open profile file
  open(unit=20, file=trim(filename), form='formatted')
  fmt='(a15, 1x, a15, 1x, a15)'
  write(20, fmt) 'Angle (deg)', 'Tb hor (K)', 'Tb vert (K)'

! convert level values to text and write to file
    fmt = '(f15.8, 1x, f15.8, 1x, f15.8)'
  do iang = 1, nstream_surf
    write(20,fmt) quad_ang(iang), Tbo_str_mat(iang,1), Tbo_str_mat(iang,2)
  end do

! close profile file
  close(unit=20)

  end subroutine write_radiation_profile
!
!====================================================================
  subroutine write_atmospheric_profile(atm_loc, filename)
!====================================================================
! writes atmospheric profile to standard text file
!
! History:
!  3/18/2021 Kevin Schaefer created routine
!  5/23/2021 Kevin Schaefer switched to local atmospheric profile
!--------------------------------------------------------------------
  use dotlrt_variables
  use profiles
  implicit none

! input variables
  type(profile_type) atm_loc(max_nlev) ! (variable) local atmospheric profile
  character*250 filename

! internal variables  
  character*250 line     ! (-) output line
  character*250 fmt      ! (-) output format spec
  character*25 txt_val1  ! (-) text version of value
  integer ilev           ! (-) level index
  integer ivar           ! (-) value index
  character*25 values(nvar_prof)! (varies) temporary write variable

! profile values
! values(1) (km) height
! values(2) (mb) atmospheric pressure
! values(3) (K) atmospheric temperature
! values(4) (g/m^3) water vapor density
! values(5) (g/m^3) cloud liquid water density
! values(6) (g/m^3) rain density
! values(7) (g/m^3) ice density
! values(8) (g/m^3) snow density
! values(9) (g/m^3) graupel density

! print message
  print*, 'Write Single Text Atmospheric Profile'
  print*, '    ', trim(filename)
!
! open profile file
  open(unit=20, file=trim(filename), form='formatted')

! write number of levels
  fmt='(i3.0)'
  write(txt_val1, fmt) nlev
  fmt='(a3)'  
  write(20,fmt) trim(adjustl(txt_val1))

! convert level values to text and write to file
  do ilev = 1, nlev
    fmt = '(f15.8)'
    write(values(1), fmt) atm_loc(ilev)%hgt_mid   ! (km)
    write(values(2), fmt) atm_loc(ilev)%press     ! (mb)
    write(values(3), fmt) atm_loc(ilev)%temp      ! (K)
    fmt = '(e15.8)'
    write(values(4), fmt) atm_loc(ilev)%humid     ! (g/m^3)
    write(values(5), fmt) atm_loc(ilev)%clw%dens  ! (g/m^3)
    write(values(6), fmt) atm_loc(ilev)%rain%dens ! (g/m^3)
    write(values(7), fmt) atm_loc(ilev)%ice%dens  ! (g/m^3)
    write(values(8), fmt) atm_loc(ilev)%snow%dens ! (g/m^3)
    write(values(9), fmt) atm_loc(ilev)%grpl%dens ! (g/m^3)

    do ivar =1, nvar_prof
      if(ivar==1) then
        line=trim(values(ivar))
      else
        line=trim(line)//','//trim(values(ivar))
      endif
    enddo
    write(20,*) trim(adjustl(line))
  end do

! close profile file
  close(unit=20)

  end subroutine write_atmospheric_profile
