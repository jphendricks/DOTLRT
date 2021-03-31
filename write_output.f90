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
    read(20,*, iostat=status) var(ivar)%name, var(ivar)%units, var(ivar)%ndim, var(ivar)%dim, var(ivar)%long_name
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
  subroutine create_rad_file(ichan)
!===================================================================================
! creates a new radiation diagnostic output file
!
! History:
!  10/20/2020 Kevn Schaefer created routine from create_qp2 from SiBCASA
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
  integer filid  ! netcdf file id number
  integer status       ! return status of netcdf functions
  integer did_lat      ! latitude dimension id #
  integer did_lon      ! longitude dimension id #
  integer did_stream      ! stream dimension id #
  integer did_pol      ! polarization dimension id #
  integer did_lev      ! level dimension id #
  integer ivar    ! variable index
  integer varid(n_out_var)    ! variable id number
  character(250) filename ! filename string
  character*50 subname ! subroutine name

! print message
  print*, 'Write radiation netcdf output file'

! set subroutine name
  subname='create_rad_file'
!
! create file
  write( filename, '(a,i2.2,a)' ) trim(out_path)//'rad_', ichan, '.nc'
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

! define variables
  do ivar = 1, n_out_var
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
  enddo
!
! switch from definition mode to data mode
 status = nf90_enddef( filid )
 if(status/=nf90_noerr) call handle_err(status,subname,14)

! write output
  status = nf90_put_var(filid, varid(1), xlat)
  if(status/=nf90_noerr) call handle_err(status,trim(subname),15)

  status = nf90_put_var(filid, varid(2), xlong)
  if(status/=nf90_noerr) call handle_err(status,trim(subname),16)

  status = nf90_put_var(filid, varid(3), Tbo_wrt)
  if(status/=nf90_noerr) call handle_err(status,trim(subname),17)

  status = nf90_put_var(filid, varid(4), Tbo_str_wrt)
  if(status/=nf90_noerr) call handle_err(status,trim(subname),18)

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

  end subroutine create_rad_file
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
    write(values(1), fmt) height(save_ilon, save_ilat, ilev) ! (km)
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
    write(20,fmt) quad_angle(iang), Tbo_str_mat(iang,1), Tbo_str_mat(iang,2)
  end do

! close profile file
  close(unit=20)

  end subroutine write_radiation_profile
