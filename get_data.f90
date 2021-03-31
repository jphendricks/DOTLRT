!
!==============================================================================
subroutine get_wrf_data()
!==============================================================================
! This program reads wrf NC files
!
! history:
!  Jason Otkin Initial code  
!  5/26/2009 Agnes Lim adapted code (ahlim@ssec.wisc.edu)
!  9/29/2020 Kevin Schaefer switched to variables module
!  10/14/2020 Kevin Schaefer moved variable definitions to profiles module
!  10/14/2020 Kevin Schaefer changed wind variables to optional
!  10/16/2020 Kevin Schaefer added total pressure, temperature, height calculations
!  10/20/2020 Kevin Schaefer accounted for land altitude in height
!-------------------------------------------------------------------------------

use dotlrt_variables
use profiles
use netcdf
implicit none

! local variables
integer ncdfID       ! (-) netcdf file id
integer status       ! (-) netcdf operation status code
integer ilon         ! (-) longitude index
integer ilat         ! (-) latitude index
integer ilev         ! (-) level index
character*80 varname ! (-) name of variable to extract
real(8) air_density  ! (kg/m3) density of air
character*50 subroutine ! subroutine name
!
! set subroutine name
  subroutine='get_wrf_data'

!
! Open file
  status = nf90_open(file_wrf, NF90_NOWRITE, ncdfID)
  if (status /= nf90_noerr) call handle_err(status, subroutine, 1)

!-------------------------------------------------------------------------------
! All required WRF inputsd
!-------------------------------------------------------------------------------
write(*,*) "Reading required WRF output"

!  Surface Skin Temperature
  if(flag_print_full) write(*,*) "    Surface Skin Temperature"
  varname='TSK'
  call get_surface_map(ncdfID, varname, nlon, nlat, ntime, TSK)

!  Land Mask
  if(flag_print_full) write(*,*) "    Land Mask"
  varname='LANDMASK'
  call get_surface_map(ncdfID, varname, nlon, nlat, ntime, LANDMASK)

!  Lake Mask
  if(flag_print_full) write(*,*) "    Lake Mask"
  varname='LAKEMASK'
  call get_surface_map(ncdfID, varname, nlon, nlat, ntime, lakemask)

!  Potential Temperature Perturbation
  if(flag_print_full) write(*,*) "    Potential Temperature Perturbation"
  varname='T'
  call get_atm_profile(ncdfID, varname, nlon, nlat, nlev_max, ntime, T)

!  geopotential perturbation
  if(flag_print_full) write(*,*) "    Geopotential Perturbation"
  varname='PH'
  call get_atm_profile(ncdfID, varname, nlon, nlat, nlev_max, ntime, PH)

!  base state geopotential
  if(flag_print_full) write(*,*) "    Base State Geopotential"
  varname='PHB'
  call get_atm_profile(ncdfID, varname, nlon, nlat, nlev_max, ntime, PHB)

! Terrain Height
  if(flag_print_full) write(*,*) "    Terrian Height"
  varname='HGT'
  call get_surface_map(ncdfID, varname, nlon, nlat, ntime, HGT)

!  pressure perturbation
  if(flag_print_full) write(*,*) "    Pressure Perturbation"
  varname='P'
  call get_atm_profile(ncdfID, varname, nlon, nlat, nlev_max, ntime, P)

!  base state pressure
  if(flag_print_full) write(*,*) "    Base State Pressure"
  varname='PB'
  call get_atm_profile(ncdfID, varname, nlon, nlat, nlev_max, ntime, PB)

!  Water Vapor Mixing Ratio
  if(flag_print_full) write(*,*) "    Water Vapour Mixing Ratio"
  varname='QVAPOR'
  call get_atm_profile(ncdfID, varname, nlon, nlat, nlev_max, ntime, QVAPOR)

!  Cloud Liquid Water Mixing Ratio
  if(flag_print_full) write(*,*) "    Cloud Water Mixing Ratio"
  varname='QCLOUD'
  call get_atm_profile(ncdfID, varname, nlon, nlat, nlev_max, ntime, QCLOUD)

!  Rain Mixing Ratio
  if(flag_print_full) write(*,*) "    Rainwater Mixing Ratio"
  varname='QRAIN'
  call get_atm_profile(ncdfID, varname, nlon, nlat, nlev_max, ntime, QRAIN)

!  Ice Mixing Ratio
  if(flag_print_full) write(*,*) "    Ice Mixing Ratio"
  varname='QICE'
  call get_atm_profile(ncdfID, varname, nlon, nlat, nlev_max, ntime, QICE)

!  Snow Mixing Ratio
  if(flag_print_full) write(*,*) "    Snow Mixing Ratio"
  varname='QSNOW'
  call get_atm_profile(ncdfID, varname, nlon, nlat, nlev_max, ntime, QSNOW)

!  Graupel Mixing Ratio
  if(flag_print_full) write(*,*) "    Graupel Mixing Ratio"
  varname='QGRAUP'
  call get_atm_profile(ncdfID, varname, nlon, nlat, nlev_max, ntime, QGRAUP)

!-------------------------------------------------------------------------------
! get optional input for ocean model
!-------------------------------------------------------------------------------
  if (trim(ocean_mod) == 'Wilheit') then
    write(*,*) "Reading required Wilheit ocean model inputs"

! 10 meter U wind component
    if(flag_print_full) write(*,*) "    10m u-wind"
    varname='U'
    call get_surface_map(ncdfID, varname, nlon, nlat, ntime, u10m)
      
! 10 meter V wind component
    if(flag_print_full) write(*,*) "    10m v-wind"
    varname='V'
    call get_surface_map(ncdfID, varname, nlon, nlat, ntime, v10m)

! calculate total horizontal wind speed
    if(flag_print_full) write(*,*) "    Calculate Total Wind Speed"
    do ilon = 1,nlon
      do ilat = 1,nlat
        wind(ilon, ilat) = sqrt(u10m(ilon, ilat)*u10m(ilon, ilat) + v10m(ilon, ilat)*v10m(ilon, ilat))
      enddo
    enddo

! deallocate component winds
    deallocate(u10m)
    deallocate(v10m)
  endif

!-------------------------------------------------------------------------------
! get optional input required for output files
!-------------------------------------------------------------------------------
  if(save_rad_file) then
    write(*,*) "Reading stuff required for output files"

!  Latitude
    if(flag_print_full) write(*,*) "    Latitude"
    varname='XLAT'
    call get_surface_map(ncdfID, varname, nlon, nlat, ntime, XLAT)

!  Longitude
    if(flag_print_full) write(*,*) "    Longitude"
    varname='XLONG'
    call get_surface_map(ncdfID, varname, nlon, nlat, ntime, XLONG)
  endif
!-------------------------------------------------------------------------------
! get ancillary WRF output
!-------------------------------------------------------------------------------
  if (flag_read_anc) then
    write(*,*) "Reading ancillary WRF inputs"
      
!  Surface Pressure
    if(flag_print_full) write(*,*) "    Surface Pressure"
    varname='PSFC'
    call get_surface_map(ncdfID, varname, nlon, nlat, ntime, PSFC)

!  Humidity at 2 meters
    if(flag_print_full) write(*,*) "    Humidity at 2m"
    varname='Q2'
    call get_surface_map(ncdfID, varname, nlon, nlat, ntime, Q2)

!  Temperature at 2 meters
    if(flag_print_full) write(*,*) "    Temperature at 2m"  
    varname='T2'
    call get_surface_map(ncdfID, varname, nlon, nlat, ntime, T2)
  endif


  status = nf90_close(ncdfID  )
  if (status /= nf90_noerr) call handle_err(status, subroutine, 2)

!-------------------------------------------------------------------------------
! calculations and conversions
!-------------------------------------------------------------------------------
  write(*,*) "calculations and unit conversions"

! combined land mask differentiating between land, lake, and ocean
! original:
!   landmask 1 = land, 0 = water
!   lakemask 1 = lake, 0 = non-lake
! combined:
!   2 = land, 1 = lake, 0 = ocean
  if(flag_print_full) write(*,*) "    Create Land/Lake/Ocean Mask"
  landmask = landmask*2 + lakemask
  deallocate(lakemask)

! Calculate total pressure
! Convert from Pa to mb
  if(flag_print_full) write(*,*) "    Calculate Total Pressure"
  press=(P+PB)/100.d0
  deallocate(P)
  deallocate(PB)

! Calculate height above 1000 or 1013 mb from geopotential
! convert from m to km
  if(flag_print_full) write(*,*) "    Calculate Height from Geopotential"
  do ilev = 1,nlev_max
    do ilat = 1,nlat
      do ilon = 1,nlon
        height(ilon,ilat,ilev) = ((PH(ilon,ilat,ilev)+PHB(ilon,ilat,ilev))/grav-HGT(ilon,ilat))/1000.d0
      enddo
    enddo
  enddo
  deallocate(PH)
  deallocate(PHB)
  deallocate(HGT)

! calculate temperature from potential temperature
  if(flag_print_full) write(*,*) "    Calculate temperature from potential temperature"
  do ilev = 1,nlev_max
    do ilat = 1,nlat
      do ilon = 1,nlon
        temp(ilon,ilat,ilev) = (T(ilon,ilat,ilev)+base_pot_temp)*(press(ilon,ilat,ilev)*100/press_ref)**gamma
      enddo
    enddo
  enddo

! Convert hydrometeors from mixing ratio to densities
! get density using pressure, temperature, and ideal gas law
! convert from kg/m3 to g/m3
  if(flag_print_full) write(*,*) "    Convert Mixing Ratios to densities"
  do ilev = 1,nlev_max
    do ilat = 1,nlat
      do ilon = 1,nlon
      
        air_density=press(ilon,ilat,ilev)/temp(ilon,ilat,ilev)/rgas
	
        QVAPOR(ilon,ilat,ilev) = QVAPOR(ilon,ilat,ilev)*air_density*1000.d0
        QICE(ilon,ilat,ilev)   = QICE(ilon,ilat,ilev)  *air_density*1000.d0
        QCLOUD(ilon,ilat,ilev) = QCLOUD(ilon,ilat,ilev)*air_density*1000.d0
        QSNOW(ilon,ilat,ilev)  = QSNOW(ilon,ilat,ilev) *air_density*1000.d0
        QGRAUP(ilon,ilat,ilev) = QGRAUP(ilon,ilat,ilev)*air_density*1000.d0
        QRAIN(ilon,ilat,ilev)  = QRAIN(ilon,ilat,ilev) *air_density*1000.d0
      enddo
    enddo
  enddo

return
end subroutine get_wrf_data
!
!==============================================================================
subroutine get_atm_profile(ncdfID, varname, nlon, nlat, nlev_max, ntime, atm_profile)
!==============================================================================
! This program reads 3-dimension atmospheric profiles from wrf netcdf files
!
! history:
!   9/30/20 Kevin Schaefer created routine
!-------------------------------------------------------------------------------
use netcdf
implicit none
!
! input variables
integer,intent(in) :: ncdfID       ! netcdf file id
character*80,intent(in) :: varname ! (-) variable name to extract
integer,intent(in) :: nlat         ! (-) number of latitude points
integer,intent(in) :: nlon         ! (-) number of longitude points
integer,intent(in) :: nlev_max     ! (-) max number of vertical levels
integer,intent(in) :: ntime        ! (-) time index
!
! output variables
real(8),intent(out) :: atm_profile(nlon, nlat, nlev_max)
!
! internal variables
character*50 subroutine ! subroutine name
integer varID  ! netcdf variable ID
integer status ! netcdf operation status flag
integer lon_ll ! longitude start index lower left
integer lat_ll ! latitude start index lower left
!
! set subroutine name
subroutine='get_atm_profile'
!
! start indeces
lon_ll=1
lat_ll=1
!
! get variable id
status = nf90_inq_varid(ncdfID, trim(varname), varID)
if(status /= nf90_NoErr) then
  print*, trim(varname)
  call handle_err(status, subroutine, 1)
endif
!
! get variable
status = nf90_get_var(ncdfID,varID,atm_profile,start=(/lon_ll,lat_ll,1,ntime/), count=(/nlon,nlat,nlev_max,ntime/) )
if (status /= nf90_noerr) call handle_err(status, subroutine, 2)
!
return
end subroutine get_atm_profile
!
!==============================================================================
subroutine get_surface_map(ncdfID, varname, nlon, nlat, ntime, surf_map)
!==============================================================================
! This program reads a 2-dimension surface map from wrf netcdf files
!
! history:
!   9/30/20 Kevin Schaefer created routine
!-------------------------------------------------------------------------------
use netcdf
implicit none
!
! input variables
integer,intent(in) :: ncdfID       ! netcdf file id
character*80,intent(in) :: varname ! (-) variable name to extract
integer,intent(in) :: nlat         ! (-) number of latitude points
integer,intent(in) :: nlon         ! (-) number of longitude points
integer,intent(in) :: ntime        ! (-) time index
!
! output variables
real(8),intent(out) :: surf_map(nlon, nlat)
!
! internal variables
character*50 subroutine ! subroutine name
integer varID  ! netcdf variable ID
integer status ! netcdf operation status flag
integer lon_ll ! longitude start index
integer lat_ll ! latitude start index
!
! set subroutine name
subroutine='get_surface_map'
!!
! start indeces
lon_ll=1
lat_ll=1
!
! get variable id
status = nf90_inq_varid(ncdfID, trim(varname), varID)
if(status /= nf90_NoErr) then
  print*, trim(varname)
  call handle_err(status, subroutine, 1)
endif
!
! get variable
status = nf90_get_var(ncdfID,varID,surf_map,start=(/lon_ll,lat_ll,ntime/), count=(/nlon,nlat,ntime/) )
if (status /= nf90_noerr) call handle_err(status, subroutine, 2)
!
return
end subroutine get_surface_map
