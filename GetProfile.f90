!
!====================================================================
subroutine hydro_layer_geometry(atm_loc, profile)
!====================================================================
! calculates detailed statistics of hydrometeor layer
!
! History:
!  11/4/2020 Kevin Schaefer created routine
!  1/21/2021 Kevin Schaefer added htop, hbot, hmax, hdel, ave, max, min
!  1/24/2021 Kevin Schaefer changed name to hydro_layer_geometry
!  7/8/2021  Kevin Schaefer switched to local atm profile variable tree
!  7/8/2021  Kevin Schaefer corrected frz_height calculation
!--------------------------------------------------------------------
  use dotlrt_variables
  use profiles
  implicit none

! input variable
  type(profile_type) atm_loc(max_nlev) ! (variable) local atmospheric profile
  real(8) profile(nlev) ! (g/m3) generic hydrometeor profile

! internal variables  
  integer ilev     ! (-) level index
  integer indx_top ! (-) level index
  integer indx_bot ! (-) level index
  real(8) var      ! (varies) generic variance
  real(8) nlayer   ! (-) number of layers
  real(8) slope    ! (km/K) slope of H vs temp
  real(8) t_max    ! (K) maximum temp in profile
  real(8) t_min    ! (K) minimum temp in profile
  integer frz_lev  ! (-) freezing level index
  integer max_lev  ! (-) peak temp level index

! calculate column total
  geom%tot = 0.d0
  do ilev = 1, nlev
    geom%tot    = geom%tot    + atm_loc(ilev)%hgt_del*profile(ilev)
  end do

! check if the column has any hydrometeors at all
  if(geom%tot == 0.d0) then
    geom%h_bot = 0.d0
    geom%h_top = 0.d0
    geom%h_del = 0.d0
    geom%min   = 0.d0
    geom%max   = 0.d0
    geom%ave   = 0.d0
    geom%std   = 0.d0
    geom%h_max = 0.d0
    geom%h_ave = 0.d0
    return
  endif

! find bottom of hydrometeor layer
  do ilev = 1, nlev
    indx_bot=ilev
    if(profile(ilev)/=0.d0) exit
  end do
  geom%h_bot=atm_loc(indx_bot)%hgt_bot

! find top of hydrometeor layer
  do ilev = nlev,1,-1
    indx_top=ilev
    if(profile(ilev)/=0.d0) exit
  end do
  geom%h_top=atm_loc(indx_top)%hgt_top
  
! calculate thickness of hydrometeor layer
  geom%h_del=geom%h_top-geom%h_bot

! hydrometeor min and max in layer
  geom%min=minval(profile(indx_bot:indx_top))
  geom%max=maxval(profile(indx_bot:indx_top))

! average hydrometeor density per layer
  nlayer=dble(indx_top-indx_bot+1)
  geom%ave=geom%tot/nlayer

! standard deviation hydrometeor density per layer
  var=0.d0
  do ilev = indx_bot,indx_top
    var=var+(profile(ilev)-geom%ave)*(profile(ilev)-geom%ave)
  end do
  var=var/nlayer
  geom%std=sqrt(var)

! calculate height of maximum of hydrometeor value
  geom%h_max=atm_loc(indx_bot)%hgt_mid
  do ilev = indx_bot,indx_top
    if (geom%max == profile(ilev)) geom%h_max=atm_loc(ilev)%hgt_mid
  end do

! calculate mean height of hydrometeor layer
  geom%h_ave=0.d0
  do ilev = indx_bot,indx_top
    geom%h_ave=geom%h_ave+atm_loc(ilev)%hgt_mid
  end do
  geom%h_ave=geom%h_ave/nlayer

! calculate variance and standard deviation of height in hydrometeor layer
  var=0.d0
  do ilev = indx_bot,indx_top
    var=var+(atm_loc(ilev)%hgt_mid-geom%h_ave)*(atm_loc(ilev)%hgt_mid-geom%h_ave)
  end do
  var=var/nlayer
  geom%h_std=sqrt(var)
  
! calculate freeze height
  t_min = minval(atm_loc(:)%temp)
  t_max = maxval(atm_loc(:)%temp)
  if(t_max < t_frz) geom%h_frz = atm_loc(1)%hgt_bot
  if(t_min > t_frz) geom%h_frz = atm_loc(nlev)%hgt_top
  if(t_max >= t_frz .and. t_min <= t_frz) then
    do ilev = 1, nlev
      if (atm_loc(ilev)%temp == t_max) max_lev = ilev
    enddo
    do ilev = max_lev, nlev
      if (atm_loc(ilev)%temp < t_frz) exit
      frz_lev = ilev
    enddo
    slope = (atm_loc(frz_lev)%hgt_mid-atm_loc(frz_lev+1)%hgt_mid)/&
      (atm_loc(frz_lev)%temp-atm_loc(frz_lev+1)%temp)
    geom%h_frz=atm_loc(frz_lev)%hgt_mid+slope*(t_frz-atm_loc(frz_lev)%temp)
  endif

end subroutine hydro_layer_geometry

!====================================================================
subroutine construct_atm_layer_geometry(atm_loc)
!====================================================================
! calculates atmospheric layer geometry
!
! History:
!  4/10/2021 Kevin Schaefer extracted routine from construct_reduced_dim_atm_profile
!  5/23/2021 Kevin Schaefer added local atmospheric profile
!  7/8/2021  Kevin Schaefer switched to local atm profile variable tree
!--------------------------------------------------------------------
  use dotlrt_variables
  use profiles
  implicit none

! input variables
  type(profile_type) atm_loc(max_nlev) ! (variable) local atmospheric profile
!
! internal variables  
  integer ilev  ! (-) level index

! calculate layer tops
  do ilev = 1, nlev-1
    atm_loc(ilev)%hgt_top = (atm_loc(ilev)%hgt_mid+atm_loc(ilev+1)%hgt_mid)/2
  end do

! calculate layer bottoms
  do ilev = 2, nlev
    atm_loc(ilev)%hgt_bot = (atm_loc(ilev)%hgt_mid+atm_loc(ilev-1)%hgt_mid)/2
  end do

! top of top layer and bottom of bottom layer
  atm_loc(nlev)%hgt_top = atm_loc(nlev)%hgt_mid + (atm_loc(nlev)%hgt_mid-atm_loc(nlev)%hgt_bot)
  atm_loc(1)%hgt_bot = atm_loc(1)%hgt_mid - (atm_loc(1)%hgt_top-atm_loc(1)%hgt_mid)

! calculate layer thicknesses
! convert  from km to meters
  do ilev = 1, nlev
    atm_loc(ilev)%hgt_del = (atm_loc(ilev)%hgt_top - atm_loc(ilev)%hgt_bot) * 1000
  end do

! print geometry
!  do ilev = 1, nlev
!    print*, atm_loc(ilev)%hgt_bot, atm_loc(ilev)%hgt_mid, atm_loc(ilev)%hgt_top, atm_loc(ilev)%hgt_del
!  end do

end subroutine construct_atm_layer_geometry

!====================================================================
subroutine construct_reduced_dim_atm_profile(atm_loc)
!====================================================================
! constructs the reduced dimension atmospheric profile for hydrometeors
! reduced means no vertical dimension
! Used to match profile to cluster
! assumes the layer height is the exact middle of the layer
!
! History:
!  11/3/2020 Kevin Schaefer created routine
!  2/22/2021 Kevin Schaefer corrected errors in fractions per layer
!  4/10/2021 Kevin Schaefer moved layer thickness calc to separate routine
!--------------------------------------------------------------------
  use dotlrt_variables
  use profiles
  implicit none

! input variable
  type(profile_type) atm_loc(max_nlev) ! (variable) local atmospheric profile

!
! internal variables  
  integer ilev  ! (-) level index
!
! variable key
! atm_loc(ilev)%hgt_mid     ! (km) height above sea level
! atm_loc(ilev)%press       ! (mb) atmospheric pressure
! atm_loc(ilev)%temp        ! (K) atmospheric temperature
! atm_loc(ilev)%humid       ! (g m-3) water vapor mixing ratio
! atm_loc(ilev)%clw%dens    ! (g m-3) cloud liquid water mixing ratio
! atm_loc(ilev)%rain%dens   ! (g m-3) rain mixing ratio
! atm_loc(ilev)%ice%dens    ! (g m-3) ice mixing ratio
! atm_loc(ilev)%snow%dens   ! (g m-3) snow mixing ratio
! atm_loc(ilev)%grpl%dens   ! (g m-3) graupel mixing ratio

! calculate total hydrometeor density
  do ilev = 1, nlev
    atm_loc(ilev)%hydro%dens = atm_loc(ilev)%clw%dens
    atm_loc(ilev)%hydro%dens = atm_loc(ilev)%hydro%dens + atm_loc(ilev)%rain%dens
    atm_loc(ilev)%hydro%dens = atm_loc(ilev)%hydro%dens + atm_loc(ilev)%ice%dens
    atm_loc(ilev)%hydro%dens = atm_loc(ilev)%hydro%dens + atm_loc(ilev)%snow%dens
    atm_loc(ilev)%hydro%dens = atm_loc(ilev)%hydro%dens + atm_loc(ilev)%grpl%dens
  end do

! calculate total precipitation density
  do ilev = 1, nlev
    atm_loc(ilev)%precip%dens = atm_loc(ilev)%rain%dens
    atm_loc(ilev)%precip%dens = atm_loc(ilev)%precip%dens + atm_loc(ilev)%snow%dens
    atm_loc(ilev)%precip%dens = atm_loc(ilev)%precip%dens + atm_loc(ilev)%grpl%dens
  end do

! calculate total cloud density
  do ilev = 1, nlev
    atm_loc(ilev)%cloud%dens = atm_loc(ilev)%clw%dens
    atm_loc(ilev)%cloud%dens = atm_loc(ilev)%cloud%dens + atm_loc(ilev)%ice%dens
  end do

! calculate the layer geometry of each hydrometeor
    call hydro_layer_geometry(atm_loc, atm_loc(:)%clw%dens)
    cld%clw=geom

    call hydro_layer_geometry(atm_loc, atm_loc(:)%rain%dens)
    cld%rain=geom

    call hydro_layer_geometry(atm_loc, atm_loc(:)%ice%dens)
    cld%ice=geom

    call hydro_layer_geometry(atm_loc, atm_loc(:)%snow%dens)
    cld%snow=geom

    call hydro_layer_geometry(atm_loc, atm_loc(:)%grpl%dens)
    cld%grpl=geom

    call hydro_layer_geometry(atm_loc, atm_loc(:)%cloud%dens)
    cld%cloud=geom

    call hydro_layer_geometry(atm_loc, atm_loc(:)%precip%dens)
    cld%precip=geom

    call hydro_layer_geometry(atm_loc, atm_loc(:)%hydro%dens)
    cld%hydro=geom

! calculate fraction of column total variables
  if(cld%hydro%tot > 0d0) then
    cld%clw%frac  = cld%clw%tot/cld%hydro%tot
    cld%rain%frac = cld%rain%tot/cld%hydro%tot
    cld%ice%frac  = cld%ice%tot/cld%hydro%tot
    cld%snow%frac = cld%snow%tot/cld%hydro%tot
    cld%grpl%frac = cld%grpl%tot/cld%hydro%tot
    cld%hydro%frac = 1.d0
  endif
  
! calculate fraction per layer of column total cloud
  if(cld%cloud%tot > 0d0) then
    do ilev = 1, nlev
      atm_loc(ilev)%clw%fcloud   = atm_loc(ilev)%hgt_del*atm_loc(ilev)%clw%dens/cld%cloud%tot
      atm_loc(ilev)%rain%fcloud  = 0.d0
      atm_loc(ilev)%ice%fcloud   = atm_loc(ilev)%hgt_del*atm_loc(ilev)%ice%dens/cld%cloud%tot
      atm_loc(ilev)%snow%fcloud  = 0.d0
      atm_loc(ilev)%grpl%fcloud  = 0.d0
      atm_loc(ilev)%cloud%fcloud = atm_loc(ilev)%hgt_del*atm_loc(ilev)%cloud%dens/cld%cloud%tot
      atm_loc(ilev)%precip%fcloud= 0.d0
      atm_loc(ilev)%hydro%fcloud = 0.d0
    end do
  endif

! calculate fraction per layer of column total precip
  if(cld%precip%tot > 0d0) then
    do ilev = 1, nlev
      atm_loc(ilev)%clw%fprecip   = 0.d0
      atm_loc(ilev)%rain%fprecip  = atm_loc(ilev)%hgt_del*atm_loc(ilev)%rain%dens/cld%precip%tot
      atm_loc(ilev)%ice%fprecip   = 0.d0
      atm_loc(ilev)%snow%fprecip  = atm_loc(ilev)%hgt_del*atm_loc(ilev)%snow%dens/cld%precip%tot
      atm_loc(ilev)%grpl%fprecip  = atm_loc(ilev)%hgt_del*atm_loc(ilev)%grpl%dens/cld%precip%tot
      atm_loc(ilev)%cloud%fprecip = 0.d0
      atm_loc(ilev)%precip%fprecip= atm_loc(ilev)%hgt_del*atm_loc(ilev)%precip%dens/cld%precip%tot
      atm_loc(ilev)%hydro%fprecip = 0.d0
    end do
  endif

! calculate fraction per layer of column total hydrometeor
  if(cld%hydro%tot > 0d0) then
    do ilev = 1, nlev
      atm_loc(ilev)%clw%fhydro   = atm_loc(ilev)%hgt_del*atm_loc(ilev)%clw%dens/cld%hydro%tot
      atm_loc(ilev)%rain%fhydro  = atm_loc(ilev)%hgt_del*atm_loc(ilev)%rain%dens/cld%hydro%tot
      atm_loc(ilev)%ice%fhydro   = atm_loc(ilev)%hgt_del*atm_loc(ilev)%ice%dens/cld%hydro%tot
      atm_loc(ilev)%snow%fhydro  = atm_loc(ilev)%hgt_del*atm_loc(ilev)%snow%dens/cld%hydro%tot
      atm_loc(ilev)%grpl%fhydro  = atm_loc(ilev)%hgt_del*atm_loc(ilev)%grpl%dens/cld%hydro%tot
      atm_loc(ilev)%cloud%fhydro = atm_loc(ilev)%hgt_del*atm_loc(ilev)%cloud%dens/cld%hydro%tot
      atm_loc(ilev)%precip%fhydro= atm_loc(ilev)%hgt_del*atm_loc(ilev)%precip%dens/cld%hydro%tot
      atm_loc(ilev)%hydro%fhydro = atm_loc(ilev)%hgt_del*atm_loc(ilev)%hydro%dens/cld%hydro%tot
    end do
  endif
  
! calculate fraction per layer of column total per species
  if(cld%clw%tot > 0d0) then
    do ilev = 1, nlev
      atm_loc(ilev)%clw%ftot = atm_loc(ilev)%hgt_del*atm_loc(ilev)%clw%dens/cld%clw%tot
    end do
  endif
  if(cld%rain%tot > 0d0) then
    do ilev = 1, nlev
      atm_loc(ilev)%rain%ftot = atm_loc(ilev)%hgt_del*atm_loc(ilev)%rain%dens/cld%rain%tot
    end do
  endif
  if(cld%ice%tot > 0d0) then
    do ilev = 1, nlev
      atm_loc(ilev)%ice%ftot = atm_loc(ilev)%hgt_del*atm_loc(ilev)%ice%dens/cld%ice%tot
    end do
  endif
  if(cld%snow%tot > 0d0) then
    do ilev = 1, nlev
      atm_loc(ilev)%snow%ftot = atm_loc(ilev)%hgt_del*atm_loc(ilev)%snow%dens/cld%snow%tot
    end do
  endif
  if(cld%grpl%tot > 0d0) then
    do ilev = 1, nlev
      atm_loc(ilev)%grpl%ftot = atm_loc(ilev)%hgt_del*atm_loc(ilev)%grpl%dens/cld%grpl%tot
    end do
  endif

end subroutine construct_reduced_dim_atm_profile
!
!====================================================================
subroutine construct_atmospheric_profile(atm_loc,ilon,ilat)
!====================================================================
! extracts a single profile from WRF run
!
! History:
!  10/16/2020 Kevin Schaefer created routine
!  11/30/2020 Kevin Schaefer added surface inputs from lookup tables
! 7/8/2021 Kevin Schaefer switched to local atmospheric variable tree
!--------------------------------------------------------------------
  use dotlrt_variables
  use profiles
  implicit none
!
! input dotlrt_variables
  type(profile_type) atm_loc(max_nlev) ! (variable) local atmospheric profile
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
! QVAPOR(ilon, ilat, nlev) ! (g m-3) water vapor density
! QCLOUD(ilon, ilat, nlev) ! (g m-3) cloud liquid water density
! QRAIN(ilon, ilat, nlev)  ! (g m-3) rain density
! QICE(ilon, ilat, nlev)   ! (g m-3) ice density
! QSNOW(ilon, ilat, nlev)  ! (g m-3) snow density
! QGRAUP(ilon, ilat, nlev) ! (g m-3) graupel density
! landmask(ilon, ilat)     ! (-) land/lake/ocean mask
! TSK(ilon, ilat)          ! (K) surface skin temperature
! wind(ilon, ilat)         ! (m s-1) wind speed

! general stuff
  nlev=nlev

! vertical atmospheric profile
  do ilev = 1, nlev
    atm_loc(ilev)%hgt_mid   = height_mid(ilon, ilat, ilev) ! (km)
    atm_loc(ilev)%press     = press(ilon, ilat, ilev)  ! (mb)
    atm_loc(ilev)%temp      = temp(ilon, ilat, ilev)   ! (K)
    atm_loc(ilev)%humid     = QVAPOR(ilon, ilat, ilev) ! (g/m^3)
    atm_loc(ilev)%clw%dens  = QCLOUD(ilon, ilat, ilev) ! (g/m^3)
    atm_loc(ilev)%rain%dens = QRAIN(ilon, ilat, ilev)  ! (g/m^3)
    atm_loc(ilev)%ice%dens  = QICE(ilon, ilat, ilev)   ! (g/m^3)
    atm_loc(ilev)%snow%dens = QSNOW(ilon, ilat, ilev)  ! (g/m^3)
    atm_loc(ilev)%grpl%dens = QGRAUP(ilon, ilat, ilev) ! (g/m^3)
  end do

! surface inputs
  surf%href=sref_hor(ilon,ilat,:)
  surf%vref=sref_ver(ilon,ilat,:)
  surf%temp = TSK(ilon, ilat)

  call construct_atm_layer_geometry(atm_loc)
  call construct_reduced_dim_atm_profile(atm_loc)

end subroutine Construct_atmospheric_profile
!
!====================================================================
subroutine read_text_profile(atm_loc, filename)
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
!  5/23/2021  Kevin Schaefer added local atm profile
!--------------------------------------------------------------------
  use dotlrt_variables
  use profiles
  implicit none

! input variables
  type(profile_type) atm_loc(max_nlev) ! (variable) local atmospheric profile
  character*250 filename

! internal variables  
  real(8) atm_read(max_nlev,9)   ! (varies) temporary atmospheric profile
  integer ilev     ! (-) level index
  integer ival     ! (-) value index
  real(8) loc_var(11) ! (varies) temporary read variable

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
  open(unit=20, file=trim(filename), form='formatted', status='old')

! read number of levels
  read(20,*) nlev

! read profile
  do ilev = 1, nlev
    read(20,*) loc_var
    do ival=1,9
      atm_read(ilev,ival)=loc_var(ival)
    enddo
  end do

! close profile file
  close(unit=20)

! assign profile to internal variables
  do ilev = 1, nlev
    atm_loc(ilev)%hgt_mid     = atm_read(ilev,1) ! (km)
    atm_loc(ilev)%press       = atm_read(ilev,2) ! (mb)
    atm_loc(ilev)%temp        = atm_read(ilev,3) ! (K)
    atm_loc(ilev)%humid       = atm_read(ilev,4) ! (g/m^3)
    atm_loc(ilev)%clw%dens    = atm_read(ilev,5) ! (g/m^3)
    atm_loc(ilev)%rain%dens   = atm_read(ilev,6) ! (g/m^3)
    atm_loc(ilev)%ice%dens    = atm_read(ilev,7) ! (g/m^3)
    atm_loc(ilev)%snow%dens   = atm_read(ilev,8) ! (g/m^3)
    atm_loc(ilev)%grpl%dens   = atm_read(ilev,9) ! (g/m^3)
  end do

! summed variables
  do ilev = 1, nlev
    atm_loc(ilev)%cloud%dens  = atm_loc(ilev)%clw%dens   + atm_loc(ilev)%ice%dens
    atm_loc(ilev)%precip%dens = atm_loc(ilev)%rain%dens  + atm_loc(ilev)%snow%dens + atm_loc(ilev)%grpl%dens
    atm_loc(ilev)%hydro%dens  = atm_loc(ilev)%cloud%dens + atm_loc(ilev)%precip%dens
  end do

! layer geometry
  call construct_atm_layer_geometry(atm_loc)
  call construct_reduced_dim_atm_profile(atm_loc)

end subroutine read_text_profile
