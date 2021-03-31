!
!====================================================================
subroutine hydro_layer_geometry(profile)
!====================================================================
! calculates detailed statistics of hydrometeor layer
!
! History:
!  11/4/2020 Kevin Schaefer created routine
!  1/21/2021 Kevin Schaefer added htop, hbot, hmax, hdel, ave, max, min
!  1/24/2021 Kevin Schaefer changed name to hydro_layer_geometry
!--------------------------------------------------------------------
  use dotlrt_variables
  use profiles
  implicit none

! inputs
  real(8) profile(nlev) ! (g/m3) generic hydrometeor profile

! internal variables  
  integer ilev     ! (-) level index
  integer indx_top ! (-) level index
  integer indx_bot ! (-) level index
  real(8) var      ! (varies) generic variance
  real(8) nlayer   ! (-) number of layers
  real(8) tot      ! (varies) generic total

! find bottom of hydrometeor layer
  do ilev = 1, nlev
    indx_bot=ilev
    if(profile(ilev)/=0.d0) exit
  end do
  geom%h_bot=atm(indx_bot)%hgt_bot

! find top of hydrometeor layer
  do ilev = nlev,1,-1
    indx_top=ilev
    if(profile(ilev)/=0.d0) exit
  end do
  geom%h_top=atm(indx_top)%hgt_top
  
! calculate thickness of hydrometeor layer
  geom%h_del=geom%h_top-geom%h_bot

! hydrometeor min and max in layer
  geom%min=minval(profile(indx_bot:indx_top))
  geom%max=maxval(profile(indx_bot:indx_top))

! hydrometeor density average in layer
  tot=0.d0
  do ilev = indx_bot,indx_top
    tot = tot + profile(ilev)
  end do
  nlayer=dble(indx_top-indx_bot+1)
  geom%ave=tot/nlayer

! hydrometeor density standard deviation in layer
  var=0.d0
  do ilev = indx_bot,indx_top
    var=var+(profile(ilev)-geom%ave)*(profile(ilev)-geom%ave)
  end do
  var=var/nlayer
  geom%std=sqrt(var)

! calculate height of maximum of hydrometeor value
  geom%h_max=atm(indx_bot)%hgt
  do ilev = indx_bot,indx_top
    if (geom%max == profile(ilev)) geom%h_max=atm(ilev)%hgt
  end do

! calculate mean height of hydrometeor layer
  geom%h_ave=0.d0
  do ilev = indx_bot,indx_top
    geom%h_ave=geom%h_ave+atm(ilev)%hgt
  end do
  geom%h_ave=geom%h_ave/nlayer

! calculate variance and standard deviation of height in hydrometeor layer
  var=0.d0
  do ilev = indx_bot,indx_top
    var=var+(atm(ilev)%hgt-geom%h_ave)*(atm(ilev)%hgt-geom%h_ave)
  end do
  var=var/nlayer
  geom%h_std=sqrt(var)

end subroutine hydro_layer_geometry
!
!====================================================================
subroutine construct_reduced_dim_atm_profile()
!====================================================================
! constructs the reduced dimension atmospheric profile for hydrometeors
! reduced means no vertical dimension
! Used to match profile to cluster
! assumes the layer height is the exact middle of the layer
!
! History:
!  11/3/2020 Kevin Schaefer created routine
!--------------------------------------------------------------------
  use dotlrt_variables
  use profiles
  implicit none
!
! internal variables  
  integer ilev  ! (-) level index
  real(8) h_del ! (km) delta height layer bottom to middle
!
! variable key
! atm(ilev)%hgt         ! (km) height above sea level
! atm(ilev)%press       ! (mb) atmospheric pressure
! atm(ilev)%temp ! (K) atmospheric temperature
! atm(ilev)%humid  ! (g m-3) water vapor mixing ratio
! atm(ilev)%clw_dens    ! (g m-3) cloud liquid water mixing ratio
! atm(ilev)%rain_dens   ! (g m-3) rain mixing ratio
! atm(ilev)%ice_dens    ! (g m-3) ice mixing ratio
! atm(ilev)%snow_dens   ! (g m-3) snow mixing ratio
! atm(ilev)%grpl_dens   ! (g m-3) graupel mixing ratio

! calculate layer thicknesses
  do ilev = 1, nlev
    if (ilev==1) then
      atm(ilev)%hgt_bot=0.d0
    else
      atm(ilev)%hgt_bot=atm(ilev-1)%hgt_top
    endif
    h_del = atm(ilev)%hgt - atm(ilev)%hgt_bot
    atm(ilev)%hgt_top = atm(ilev)%hgt + h_del
    atm(ilev)%hgt_del = (atm(ilev)%hgt_top - atm(ilev)%hgt_bot) * 1000 ! convert to meters
    !print*, atm(ilev)%hgt_bot,atm(ilev)%hgt, atm(ilev)%hgt_top, atm(ilev)%hgt_del
  end do

! calculate total hydrometeor density
  do ilev = 1, nlev
    atm(ilev)%hydro_dens = atm(ilev)%clw_dens
    atm(ilev)%hydro_dens = atm(ilev)%hydro_dens + atm(ilev)%rain_dens
    atm(ilev)%hydro_dens = atm(ilev)%hydro_dens + atm(ilev)%ice_dens
    atm(ilev)%hydro_dens = atm(ilev)%hydro_dens + atm(ilev)%snow_dens
    atm(ilev)%hydro_dens = atm(ilev)%hydro_dens + atm(ilev)%grpl_dens
  end do

! calculate column total hydrometeors
  reduced%clw%tot   = 0.d0
  reduced%rain%tot  = 0.d0
  reduced%ice%tot   = 0.d0
  reduced%snow%tot  = 0.d0
  reduced%grpl%tot  = 0.d0
  reduced%hydro%tot = 0.d0
  do ilev = 1, nlev
    reduced%clw%tot   = reduced%clw%tot   + atm(ilev)%hgt_del*atm(ilev)%clw_dens
    reduced%rain%tot  = reduced%rain%tot  + atm(ilev)%hgt_del*atm(ilev)%rain_dens
    reduced%ice%tot   = reduced%ice%tot   + atm(ilev)%hgt_del*atm(ilev)%ice_dens
    reduced%snow%tot  = reduced%snow%tot  + atm(ilev)%hgt_del*atm(ilev)%snow_dens
    reduced%grpl%tot  = reduced%grpl%tot  + atm(ilev)%hgt_del*atm(ilev)%grpl_dens
    reduced%hydro%tot = reduced%hydro%tot + atm(ilev)%hgt_del*atm(ilev)%hydro_dens
  end do

! calculate the geometry of each hydrometeor layer
  if(reduced%clw%tot/=0.d0) then
    call hydro_layer_geometry(atm(:)%clw_dens)
    reduced%clw=geom
  endif
  if(reduced%rain%tot/=0.d0) then
    call hydro_layer_geometry(atm(:)%rain_dens)
    reduced%rain=geom
  endif
  if(reduced%ice%tot/=0.d0) then
    call hydro_layer_geometry(atm(:)%ice_dens)
    reduced%ice=geom
  endif
  if(reduced%snow%tot/=0.d0) then
    call hydro_layer_geometry(atm(:)%snow_dens)
    reduced%snow=geom
  endif
  if(reduced%grpl%tot/=0.d0) then
    call hydro_layer_geometry(atm(:)%grpl_dens)
    reduced%grpl=geom
  endif
  if(reduced%hydro%tot/=0.d0) then
    call hydro_layer_geometry(atm(:)%hydro_dens)
    reduced%hydro=geom
  endif

end subroutine construct_reduced_dim_atm_profile
!
!====================================================================
subroutine construct_atmospheric_profile(ichan,ilon,ilat)
!====================================================================
! extracts a single profile from WRF run
!
! History:
!  10/16/2020 Kevin Schaefer created routine
!  11/30/2020 Kevin Schaefer added surface inputs from lookup tables
!--------------------------------------------------------------------
  use dotlrt_variables
  use profiles
  implicit none
!
! input dotlrt_variables
  integer ichan ! (-) channel index
  integer ilon  ! (-) longitude index
  integer ilat  ! (-) latitude index
!
! internal variables  
  integer ilev  ! (-) level index
!
! variable key
! height(ilon, ilat, nlev) ! (km) height above sea level
! press(ilon, ilat, nlev)  ! (mb) atmospheric pressure
! temp(ilon, ilat, lev)    ! (K) atmospheric temperature
! QVAPOR(ilon, ilat, nlev) ! (g m-3) water vapor mixing ratio
! QCLOUD(ilon, ilat, nlev) ! (g m-3) cloud liquid water mixing ratio
! QRAIN(ilon, ilat, nlev)  ! (g m-3) rain mixing ratio
! QICE(ilon, ilat, nlev)   ! (g m-3) ice mixing ratio
! QSNOW(ilon, ilat, nlev)  ! (g m-3) snow mixing ratio
! QGRAUP(ilon, ilat, nlev) ! (g m-3) graupel mixing ratio
! landmask(ilon, ilat)     ! (-) land/lake/ocean mask
! TSK(ilon, ilat)          ! (K) surface skin temperature
! wind(ilon, ilat)         ! (m s-1) wind speed

! general stuff
  nlev=nlev

! vertical atmospheric profile
  do ilev = 1, nlev
    atm(ilev)%hgt       = height_mid(ilon, ilat, ilev) ! (km)
    atm(ilev)%press     = press(ilon, ilat, ilev)  ! (mb)
    atm(ilev)%temp      = temp(ilon, ilat, ilev)   ! (K)
    atm(ilev)%humid     = QVAPOR(ilon, ilat, ilev) ! (g/m^3)
    atm(ilev)%clw_dens  = QCLOUD(ilon, ilat, ilev) ! (g/m^3)
    atm(ilev)%rain_dens = QRAIN(ilon, ilat, ilev)  ! (g/m^3)
    atm(ilev)%ice_dens  = QICE(ilon, ilat, ilev)   ! (g/m^3)
    atm(ilev)%snow_dens = QSNOW(ilon, ilat, ilev)  ! (g/m^3)
    atm(ilev)%grpl_dens = QGRAUP(ilon, ilat, ilev) ! (g/m^3)
  end do

! surface inputs
  surf%href=sref_hor(ilon,ilat,ichan,:)
  surf%vref=sref_ver(ilon,ilat,ichan,:)
  surf%temp = TSK(ilon, ilat)

  call construct_reduced_dim_atm_profile()

end subroutine Construct_atmospheric_profile
!
!====================================================================
subroutine read_text_profile()
!====================================================================
! reads in single profile from WRF run
! Assumes csv format
! primarily for testing dotlrt
!
! History:
!  9/18/2020  Kevin Schaefer created routine
!  10/16/2020 Kevin Schaefer deleted all arguments duplicated in variables module
!  10/16/2020 Kevin Schaefer added assignment to profile variable tree
!  1/24/2021  Kevin Schaefer added reduced dimension call
!--------------------------------------------------------------------
  use dotlrt_variables
  use profiles
  implicit none
!
! internal variables  
  real(8), allocatable :: atminp(:,:)   ! (varies) atmospheric profile
  integer ilev     ! (-) level index
  integer ival     ! (-) value index
  real(8) loc_var(11) ! (varies) temporary read variable
!
! profile values
! loc_var(1) ! (km) height
! loc_var(2) ! (mb) atmospheric pressure
! loc_var(3) ! (K) atmospheric temperature
! loc_var(4) ! (g/m^3) water vapor density
! loc_var(5) ! (g/m^3) cloud liquid water density
! loc_var(6) ! (g/m^3) rain density
! loc_var(7) ! (g/m^3) ice density
! loc_var(8) ! (g/m^3) snow density
! loc_var(9) ! (g/m^3) graupel density
  print*, 'Read Single Text Atmospheric Profile'

! open profile file
  open(unit=20, file=trim(file_in), form='formatted', status='old')

! read number of levels
  read(20,*) nlev

! read profile
  allocate(atminp(nlev,9))
  do ilev = 1, nlev
    read(20,*) loc_var
    do ival=1,9
      atminp(ilev,ival)=loc_var(ival)
    enddo
  end do

! close profile file
  close(unit=20)

! assign profile to internal variables
  do ilev = 1, nlev
    atm(ilev)%hgt         = atminp(ilev,1) ! (km)
    atm(ilev)%press       = atminp(ilev,2) ! (mb)
    atm(ilev)%temp        = atminp(ilev,3) ! (K)
    atm(ilev)%humid       = atminp(ilev,4) ! (g/m^3)
    atm(ilev)%clw_dens    = atminp(ilev,5) ! (g/m^3)
    atm(ilev)%rain_dens   = atminp(ilev,6) ! (g/m^3)
    atm(ilev)%ice_dens    = atminp(ilev,7) ! (g/m^3)
    atm(ilev)%snow_dens   = atminp(ilev,8) ! (g/m^3)
    atm(ilev)%grpl_dens   = atminp(ilev,9) ! (g/m^3)
  end do

  call construct_reduced_dim_atm_profile()

! dellaocate
  deallocate(atminp)

end subroutine read_text_profile
