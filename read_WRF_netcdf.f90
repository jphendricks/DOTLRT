
!==============================================================================
SUBROUTINE read_WRF_netcdf_dimensions()
!==============================================================================
! This program reads in WRF model output to extract input data for MRT
!
! Input
!  WRF Model Output file
!
! Output
!  U10, V10, T, P, PB, Q2, T2, PSFC, QVAPOR, QCLOUD, LANDUSE, HGT,
!	TSK, XLAT, XLONG
!
! History:
!  5/22/2009 Agness Lim created program (alim@ssec.wisc.edu)
!  9/29/2020 Kevin Schaefer adapted routine for DOTLRT
!  10/14/2020 Kevin Schaefer added use profiles module and allocations
!------------------------------------------------------------------------------

use dotlrt_variables
use profiles
use netcdf

IMPLICIT NONE
!
! Local variables
integer ncdfID     ! (-) netcdf file id
integer numdim     ! (-) number of dimensions
integer numvars    ! (-) number of variables
integer natt       ! (-) number of attributes
integer nunlimited ! (-) number unlimited dimensions
integer dimlen     ! (-) value of dimension
integer status     ! (-) netcdf operation status code
integer nlon_st    ! (-) longitude stagger
integer nlat_st    ! (-) latitude stagger
integer nlev_st    ! (-) level stagger
integer ndate      ! (-) number of dates
integer idim       ! (-) dimension index      
character(len = nf90_max_name) :: dimnam
character*50 subroutine ! subroutine name
!
! set subroutine name
  subroutine='read_WRF_netcdf_dim'
!
! Open file
  status = nf90_open(file_in, NF90_NOWRITE, ncdfID)
  if (status /= nf90_noerr) then
    print*, trim(file_in)
    call handle_err(status, subroutine, 1)
  endif
!
!------------------------------------------------------------------------------
! read dimension information
!------------------------------------------------------------------------------
  status = nf90_inquire( ncdfID,numdim,numvars,natt, nunlimited )
  if (status /= nf90_noerr) call handle_err(status, subroutine, 2)
  write(6,*) trim(file_in)
!
! loop through dimensions
  do idim = 1,numdim
    status = nf90_inquire_dimension( ncdfID,idim,dimnam,dimlen )
    if (status /= nf90_noerr) call handle_err(status, subroutine, 3)

    if( dimnam == 'Time' ) then
      ntime=dimlen
      write(6,*) 'Number of times in file:         ', ntime
    elseif( dimnam == 'DateStrLen' ) then
      ndate=dimlen
      if(flag_print_full) write(6,*) 'Length of date string:           ', ndate
    elseif( dimnam == 'west_east_stag') then
      nlon_st=dimlen
      if(flag_print_full) write(6,*) 'West-east staggered gridpoints:  ', nlon_st
    elseif( dimnam == 'west_east') then
      nlon=dimlen
      write(6,*) 'West-east gridpoints (nlon):     ', nlon
    elseif( dimnam == 'south_north') then
      nlat=dimlen
      write(6,*) 'South-north gridpoints (nlat):   ', nlat
    elseif( dimnam == 'bottom_top') then
      nlev_max=dimlen
      write(6,*) 'Number of half eta levels (nlev_max):', nlev_max
    elseif( dimnam == 'south_north_stag') then
      nlat_st=dimlen
      if(flag_print_full) write(6,*) 'South-north staggered gridpoints:', nlat_st
    elseif( dimnam == 'bottom_top_stag') then
      nlev_st=dimlen
      if(flag_print_full) write(6,*) 'Number of full eta levels:       ', nlev_st
    endif
  enddo 

  status = nf90_close(ncdfID  )
  if (status /= nf90_noerr) call handle_err(status, subroutine, 4)

!------------------------------------------------------------------------------
! allocate profile variables
!------------------------------------------------------------------------------
! temperature variables 
  allocate(T(nlon, nlat, nlev_max))
  allocate(temp(nlon, nlat, nlev_max))

! height variables
  allocate(PH(nlon, nlat, nlev_max))
  allocate(PHB(nlon, nlat, nlev_max))
  allocate(height_ter(nlon, nlat))
  allocate(height_mid(nlon, nlat, nlev_max))

! pressure variables
  allocate(P(nlon, nlat, nlev_max))
  allocate(PB(nlon, nlat, nlev_max))
  allocate(press(nlon, nlat, nlev_max))

! hydrometeor variables
  allocate(QVAPOR(nlon, nlat, nlev_max))
  allocate(QICE(nlon, nlat, nlev_max))
  allocate(QCLOUD(nlon, nlat, nlev_max))
  allocate(QSNOW(nlon, nlat, nlev_max))
  allocate(QGRAUP(nlon, nlat, nlev_max))
  allocate(QRAIN(nlon, nlat, nlev_max))

! surface reflectance variables
  allocate(landmask(nlon, nlat))
  allocate(lakemask(nlon, nlat))
  allocate(TSK(nlon, nlat))
  allocate(sref_hor(nlon, nlat, max_nang))
  allocate(sref_ver(nlon, nlat, max_nang))
  
! allocate variables required for output files
  if(save_rad_file) then
    allocate(XLAT(nlon, nlat))
    allocate(XLONG(nlon, nlat))
  endif
  
! optional ancillary data variables
  if (flag_read_anc) then
    allocate(PSFC(nlon, nlat))
    allocate(Q2(nlon, nlat))
    allocate(T2(nlon, nlat))
  endif

! allocate optional variables for ocean model
  if (trim(ocean_mod) == 'Wilheit') then
    allocate(u10m(nlon, nlat))
    allocate(v10m(nlon, nlat))
    allocate(wind(nlon, nlat))
  endif

end subroutine read_WRF_netcdf_dimensions
