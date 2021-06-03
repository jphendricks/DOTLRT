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
  integer frz_lev  ! (-) freezing level index
  real(8) slope    ! (km/K) slope of H vs temp

! calculate column total
  geom%tot = 0.d0
  do ilev = 1, nlev
    geom%tot    = geom%tot    + atm(ilev)%hgt_del*profile(ilev)
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
  
! calculate freeze height
  do ilev = 1, nlev
    if (atm(ilev)%temp < t_frz) exit
    frz_lev = ilev
  enddo
  slope = (atm(frz_lev)%hgt-atm(frz_lev+1)%hgt)/(atm(frz_lev)%temp-atm(frz_lev+1)%temp)
  geom%h_frz=atm(frz_lev)%hgt+slope*(t_frz-atm(frz_lev)%temp)

end subroutine hydro_layer_geometry

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
!  2/22/2021 Kevin Schaefer corrected errors in fractions per layer
!--------------------------------------------------------------------
  use dotlrt_variables
  use profiles
  implicit none
!
! internal variables  
  integer ilev  ! (-) level index
!
! variable key
! atm(ilev)%hgt         ! (km) height above sea level
! atm(ilev)%press       ! (mb) atmospheric pressure
! atm(ilev)%temp        ! (K) atmospheric temperature
! atm(ilev)%humid       ! (g m-3) water vapor mixing ratio
! atm(ilev)%clw%dens    ! (g m-3) cloud liquid water mixing ratio
! atm(ilev)%rain%dens   ! (g m-3) rain mixing ratio
! atm(ilev)%ice%dens    ! (g m-3) ice mixing ratio
! atm(ilev)%snow%dens   ! (g m-3) snow mixing ratio
! atm(ilev)%grpl%dens   ! (g m-3) graupel mixing ratio

! calculate layer tops and bottoms
  do ilev = 1, nlev-1
    atm(ilev)%hgt_top = (atm(ilev)%hgt+atm(ilev+1)%hgt)/2
  end do
  do ilev = 2, nlev
    atm(ilev)%hgt_bot = (atm(ilev)%hgt+atm(ilev-1)%hgt)/2
  end do
  atm(nlev)%hgt_top = atm(nlev)%hgt + (atm(nlev)%hgt-atm(nlev)%hgt_bot)
  atm(1)%hgt_bot = atm(1)%hgt - (atm(1)%hgt_top-atm(1)%hgt)

! calculate layer thicknesses
  do ilev = 1, nlev
    atm(ilev)%hgt_del = (atm(ilev)%hgt_top - atm(ilev)%hgt_bot) * 1000 ! convert to meters
  end do

! calculate total hydrometeor density
  do ilev = 1, nlev
    atm(ilev)%hydro%dens = atm(ilev)%clw%dens
    atm(ilev)%hydro%dens = atm(ilev)%hydro%dens + atm(ilev)%rain%dens
    atm(ilev)%hydro%dens = atm(ilev)%hydro%dens + atm(ilev)%ice%dens
    atm(ilev)%hydro%dens = atm(ilev)%hydro%dens + atm(ilev)%snow%dens
    atm(ilev)%hydro%dens = atm(ilev)%hydro%dens + atm(ilev)%grpl%dens
  end do

! calculate total precipitation density
  do ilev = 1, nlev
    atm(ilev)%precip%dens = atm(ilev)%rain%dens
    atm(ilev)%precip%dens = atm(ilev)%precip%dens + atm(ilev)%snow%dens
    atm(ilev)%precip%dens = atm(ilev)%precip%dens + atm(ilev)%grpl%dens
  end do

! calculate total cloud density
  do ilev = 1, nlev
    atm(ilev)%cloud%dens = atm(ilev)%clw%dens
    atm(ilev)%cloud%dens = atm(ilev)%cloud%dens + atm(ilev)%ice%dens
  end do

! calculate the layer geometry of each hydrometeor
    call hydro_layer_geometry(atm(:)%clw%dens)
    cld%clw=geom

    call hydro_layer_geometry(atm(:)%rain%dens)
    cld%rain=geom

    call hydro_layer_geometry(atm(:)%ice%dens)
    cld%ice=geom

    call hydro_layer_geometry(atm(:)%snow%dens)
    cld%snow=geom

    call hydro_layer_geometry(atm(:)%grpl%dens)
    cld%grpl=geom

    call hydro_layer_geometry(atm(:)%cloud%dens)
    cld%cloud=geom

    call hydro_layer_geometry(atm(:)%precip%dens)
    cld%precip=geom

    call hydro_layer_geometry(atm(:)%hydro%dens)
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
      atm(ilev)%clw%fcloud   = atm(ilev)%hgt_del*atm(ilev)%clw%dens/cld%cloud%tot
      atm(ilev)%rain%fcloud  = 0.d0
      atm(ilev)%ice%fcloud   = atm(ilev)%hgt_del*atm(ilev)%ice%dens/cld%cloud%tot
      atm(ilev)%snow%fcloud  = 0.d0
      atm(ilev)%grpl%fcloud  = 0.d0
      atm(ilev)%cloud%fcloud = atm(ilev)%hgt_del*atm(ilev)%cloud%dens/cld%cloud%tot
      atm(ilev)%precip%fcloud= atm(ilev)%hgt_del*atm(ilev)%precip%dens/cld%cloud%tot
      atm(ilev)%hydro%fcloud = atm(ilev)%hgt_del*atm(ilev)%hydro%dens/cld%cloud%tot
      ! atm(ilev)%precip%fcloud= 0.d0
      ! atm(ilev)%hydro%fcloud = 0.d0
    end do
  endif

! calculate fraction per layer of column total precip
  if(cld%precip%tot > 0d0) then
    do ilev = 1, nlev
      atm(ilev)%clw%fprecip   = 0.d0
      atm(ilev)%rain%fprecip  = atm(ilev)%hgt_del*atm(ilev)%rain%dens/cld%precip%tot
      atm(ilev)%ice%fprecip   = 0.d0
      atm(ilev)%snow%fprecip  = atm(ilev)%hgt_del*atm(ilev)%snow%dens/cld%precip%tot
      atm(ilev)%grpl%fprecip  = atm(ilev)%hgt_del*atm(ilev)%grpl%dens/cld%precip%tot
      atm(ilev)%cloud%fprecip = atm(ilev)%hgt_del*atm(ilev)%cloud%dens/cld%precip%tot
      ! atm(ilev)%cloud%fprecip = 0.d0
      atm(ilev)%precip%fprecip= atm(ilev)%hgt_del*atm(ilev)%precip%dens/cld%precip%tot
      atm(ilev)%hydro%fprecip = atm(ilev)%hgt_del*atm(ilev)%hydro%dens/cld%precip%tot
      ! atm(ilev)%hydro%fprecip = 0.d0
    end do
  endif

! calculate fraction per layer of column total hydrometeor
  if(cld%hydro%tot > 0d0) then
    do ilev = 1, nlev
      atm(ilev)%clw%fhydro   = atm(ilev)%hgt_del*atm(ilev)%clw%dens/cld%hydro%tot
      atm(ilev)%rain%fhydro  = atm(ilev)%hgt_del*atm(ilev)%rain%dens/cld%hydro%tot
      atm(ilev)%ice%fhydro   = atm(ilev)%hgt_del*atm(ilev)%ice%dens/cld%hydro%tot
      atm(ilev)%snow%fhydro  = atm(ilev)%hgt_del*atm(ilev)%snow%dens/cld%hydro%tot
      atm(ilev)%grpl%fhydro  = atm(ilev)%hgt_del*atm(ilev)%grpl%dens/cld%hydro%tot
      atm(ilev)%cloud%fhydro = atm(ilev)%hgt_del*atm(ilev)%cloud%dens/cld%hydro%tot
      atm(ilev)%precip%fhydro= atm(ilev)%hgt_del*atm(ilev)%precip%dens/cld%hydro%tot
      atm(ilev)%hydro%fhydro = atm(ilev)%hgt_del*atm(ilev)%hydro%dens/cld%hydro%tot
    end do
  endif
  
! calculate fraction per layer of column total per species
  if(cld%clw%tot > 0d0) then
    do ilev = 1, nlev
      atm(ilev)%clw%ftot = atm(ilev)%hgt_del*atm(ilev)%clw%dens/cld%clw%tot
    end do
  endif
  if(cld%rain%tot > 0d0) then
    do ilev = 1, nlev
      atm(ilev)%rain%ftot = atm(ilev)%hgt_del*atm(ilev)%rain%dens/cld%rain%tot
    end do
  endif
  if(cld%ice%tot > 0d0) then
    do ilev = 1, nlev
      atm(ilev)%ice%ftot = atm(ilev)%hgt_del*atm(ilev)%ice%dens/cld%ice%tot
    end do
  endif
  if(cld%snow%tot > 0d0) then
    do ilev = 1, nlev
      atm(ilev)%snow%ftot = atm(ilev)%hgt_del*atm(ilev)%snow%dens/cld%snow%tot
    end do
  endif
  if(cld%grpl%tot > 0d0) then
    do ilev = 1, nlev
      atm(ilev)%grpl%ftot = atm(ilev)%hgt_del*atm(ilev)%grpl%dens/cld%grpl%tot
    end do
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
    atm(ilev)%clw%dens  = QCLOUD(ilon, ilat, ilev) ! (g/m^3)
    atm(ilev)%rain%dens = QRAIN(ilon, ilat, ilev)  ! (g/m^3)
    atm(ilev)%ice%dens  = QICE(ilon, ilat, ilev)   ! (g/m^3)
    atm(ilev)%snow%dens = QSNOW(ilon, ilat, ilev)  ! (g/m^3)
    atm(ilev)%grpl%dens = QGRAUP(ilon, ilat, ilev) ! (g/m^3)
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
    atm(ilev)%clw%dens    = atminp(ilev,5) ! (g/m^3)
    atm(ilev)%rain%dens   = atminp(ilev,6) ! (g/m^3)
    atm(ilev)%ice%dens    = atminp(ilev,7) ! (g/m^3)
    atm(ilev)%snow%dens   = atminp(ilev,8) ! (g/m^3)
    atm(ilev)%grpl%dens   = atminp(ilev,9) ! (g/m^3)
  end do

  call construct_reduced_dim_atm_profile()

! dellaocate
  deallocate(atminp)

end subroutine read_text_profile
