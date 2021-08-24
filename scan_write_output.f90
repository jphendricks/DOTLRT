!
!===================================================================================
  subroutine create_retrieval_output_file()
!===================================================================================
! creates a netcdf retieval output file
!
! History:
!  7/4/2021 Kevn Schaefer created routine
!-----------------------------------------------------------------------------------
  use scan_Variables
  use dotlrt_variables
  use dotlrt_output
  use profiles
  use netcdf

  implicit none
!
! local variables
  integer filid            ! (-) netcdf file id number
  integer status           ! (-) return status of netcdf functions
  integer did_npts         ! (-) latitude dimension id number
  integer ivar             ! (-) variable index
  integer varid(n_out_var) ! (-) variable id number
  character(250) filename  ! (-) filename string
  integer ifil       ! (-) file index
  character*50 subname     ! (-) subroutine name
  type (variable_spec) ret_var(max_out_var) !(varies) retrieval output variable tree

! print message
  print*, '    Create netcdf retrieval output file'

! set subroutine name
  subname='create_retrieval_output_file'

! find masked observation swath file
  print*, '    Read variable definitions'
  do ifil=1,numfiles
    if (trim(files(ifil)%type) == 'FY3_var')  filename = trim(files(ifil)%path)
  enddo

! open variable definition file
  open(unit=20, file=trim(filename), form='formatted', status='old')
  read(20,*) junk

! count variables
  n_out_var=0
  do ivar=1,max_out_var
    read(20,*, iostat=status) junk
    if(status<0) exit
    n_out_var=n_out_var+1
  enddo
  print*, '    Number output variables:', n_out_var

! close variable definition file
  close(unit=20)

! re-open variable definition file
  open(unit=20, file=trim(filename), form='formatted', status='old')
  read(20,*) junk

! read variable defitions
  do ivar=1,n_out_var
    read(20,*, iostat=status) ret_var(ivar)%file,ret_var(ivar)%name, ret_var(ivar)%units,&
      ret_var(ivar)%ndim, ret_var(ivar)%dim, ret_var(ivar)%long_name
  enddo
!
! close variable definition file
  close(unit=20)

! set standard stuff
  do ivar = 1, n_out_var
    ret_var(ivar)%missing=missing
  enddo

! find masked observation swath file
  do ifil=1,numfiles
    if (trim(files(ifil)%type) == 'FY3_out')  filename = trim(files(ifil)%path)//'FY3_out.nc'
  enddo
  print*, '    write out FY: ', trim(filename)

  status = nf90_create( trim(filename), cmode=or(nf90_clobber,nf90_64bit_offset), ncid=filid)
  if(status/=nf90_noerr) then
    print*, trim(filename)
    call handle_err(status,subname,1)
  endif

! define dimensions
  status = nf90_def_dim( filid, 'swath_npts', swath_npts, did_npts )
  if(status/=nf90_noerr) call handle_err(status,subname,2)

! define variables common to all files
  do ivar = 1, n_out_var
    status = nf90_def_var( filid, ret_var(ivar)%name, nf90_float, (/did_npts/), varid(ivar))
    if(status/=nf90_noerr) call handle_err(status,subname,7)

    status = nf90_put_att( filid, varid(ivar), 'units', trim(ret_var(ivar)%units) )
    if(status/=nf90_noerr) call handle_err(status,subname,8)

    status = nf90_put_att( filid, varid(ivar), 'long_name', trim(ret_var(ivar)%long_name) )
    if(status/=nf90_noerr) call handle_err(status,subname,9)

    status = nf90_put_att( filid, varid(ivar), 'missing_value', ret_var(ivar)%missing )
    if(status/=nf90_noerr) call handle_err(status,subname,10)
  enddo

! switch from definition mode to data mode
  status = nf90_enddef( filid )
  if(status/=nf90_noerr) call handle_err(status,subname,18)

! write output
! varid index
! 1 swath_precip(swath_npts)
! 2 swath_point(swath_npts)
! 3 swath_iter(swath_npts)
! 4 swath_cost(swath_npts)
! 5 swath_tb_sim(swath_npts)
! 6 swath_res(swath_npts)
! 7 swath_ex(swath_npts)

! write estimated precipitation and metrics
  status = nf90_put_var(filid, varid(1), swath_precip)
  if(status/=nf90_noerr) call handle_err(status,trim(subname),19)

  status = nf90_put_var(filid, varid(2), swath_point)
  if(status/=nf90_noerr) call handle_err(status,trim(subname),20)

  status = nf90_put_var(filid, varid(3), swath_iter)
  if(status/=nf90_noerr) call handle_err(status,trim(subname),21)

  status = nf90_put_var(filid, varid(4), swath_cost)
  if(status/=nf90_noerr) call handle_err(status,trim(subname),22)

  status = nf90_put_var(filid, varid(5), swath_tb_sim)
  if(status/=nf90_noerr) call handle_err(status,trim(subname),23)

  status = nf90_put_var(filid, varid(6), swath_res)
  if(status/=nf90_noerr) call handle_err(status,trim(subname),24)

  status = nf90_put_var(filid, varid(7), swath_ex)
  if(status/=nf90_noerr) call handle_err(status,trim(subname),25)

  status = nf90_put_var(filid, varid(8), swath_bad)
  if(status/=nf90_noerr) call handle_err(status,trim(subname),26)

! close file
  status = nf90_close(filid)
  if(status/=nf90_noerr) call handle_err(status,subname,27)

  end subroutine create_retrieval_output_file


