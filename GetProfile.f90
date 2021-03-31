!
!====================================================================
subroutine mean_height_hydro_layer(profile, h_ave, h_std)
!====================================================================
! calculates mean and standard deviation of hydrometeor layer
!
! History:
!  11/4/2020 Kevin Schaefer created routine
!--------------------------------------------------------------------
  use dotlrt_variables
  use profiles
  implicit none

! inputs
  real(8) profile(nlev) ! (g/m3) generic hydrometeor profile

! output
  real(8) h_ave    ! (km) mean height of hydrometeor layer
  real(8) h_std    ! (km) height standard deviation of hydrometeor layer

! internal variables  
  integer ilev     ! (-) level index
  integer indx_top ! (-) level index
  integer indx_bot ! (-) level index
  real(8) h_var    ! (km2) height variance of hydrometeor layer
  real(8) nlayer   ! (-) number of layers

! find bottom level index
  do ilev = 1, nlev
    indx_bot=ilev
    if(profile(ilev)/=0.d0) exit
  end do

! find top level index
  do ilev = nlev,1,-1
    indx_top=ilev
    if(profile(ilev)/=0.d0) exit
  end do
  nlayer=dble(indx_top-indx_bot+1)
!  print*, indx_bot, indx_top, nlayer

! calculate mean height  of hydrometeor layer
  h_ave=0.d0
  do ilev = indx_bot,indx_top
    h_ave=h_ave+atm(ilev)%hgt
  end do
  h_ave=h_ave/nlayer

! calculate variance and standard deviation  of hydrometeor layer
  h_var=0.d0
  do ilev = indx_bot,indx_top
    h_var=h_var+(atm(ilev)%hgt-h_ave)*(atm(ilev)%hgt-h_ave)
  end do
  h_var=h_var/nlayer
  h_std=sqrt(h_var)
!  print*, h_ave, h_var, h_std


end subroutine mean_height_hydro_layer
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

! layer thickness calculations
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

! column total calculations
  atm_reduced%clw_tot  = 0.d0
  atm_reduced%rain_tot = 0.d0
  atm_reduced%ice_tot  = 0.d0
  atm_reduced%snow_tot = 0.d0
  atm_reduced%grpl_tot = 0.d0
  do ilev = 1, nlev
    atm_reduced%clw_tot  = atm_reduced%clw_tot  + atm(ilev)%hgt_del*atm(ilev)%clw_dens
    atm_reduced%rain_tot = atm_reduced%rain_tot + atm(ilev)%hgt_del*atm(ilev)%rain_dens
    atm_reduced%ice_tot  = atm_reduced%ice_tot  + atm(ilev)%hgt_del*atm(ilev)%ice_dens
    atm_reduced%snow_tot = atm_reduced%snow_tot + atm(ilev)%hgt_del*atm(ilev)%snow_dens
    atm_reduced%grpl_tot = atm_reduced%grpl_tot + atm(ilev)%hgt_del*atm(ilev)%grpl_dens
  end do
  !print*, atm_reduced%clw_tot, atm_reduced%rain_tot, atm_reduced%ice_tot, atm_reduced%snow_tot, atm_reduced%grpl_tot

! Mean height and standard deviation of hydrometeor layers
  atm_reduced%clw_h_ave  = 0.d0
  atm_reduced%rain_h_ave = 0.d0
  atm_reduced%ice_h_ave  = 0.d0
  atm_reduced%snow_h_ave = 0.d0
  atm_reduced%grpl_h_ave = 0.d0

  atm_reduced%clw_h_std  = 0.d0
  atm_reduced%rain_h_std = 0.d0
  atm_reduced%ice_h_std  = 0.d0
  atm_reduced%snow_h_std = 0.d0
  atm_reduced%grpl_h_std = 0.d0

  if(atm_reduced%clw_tot/=0.d0)  call mean_height_hydro_layer(atm(:)%clw_dens,  atm_reduced%clw_h_ave,  atm_reduced%clw_h_std)
  if(atm_reduced%rain_tot/=0.d0) call mean_height_hydro_layer(atm(:)%rain_dens, atm_reduced%rain_h_ave, atm_reduced%rain_h_std)
  if(atm_reduced%ice_tot/=0.d0)  call mean_height_hydro_layer(atm(:)%ice_dens,  atm_reduced%ice_h_ave,  atm_reduced%ice_h_std)
  if(atm_reduced%snow_tot/=0.d0) call mean_height_hydro_layer(atm(:)%snow_dens, atm_reduced%snow_h_ave, atm_reduced%snow_h_std)
  if(atm_reduced%grpl_tot/=0.d0) call mean_height_hydro_layer(atm(:)%grpl_dens, atm_reduced%grpl_h_ave, atm_reduced%grpl_h_std)

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
  integer i
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

  if(ilat==292 .and. ilon==397) then
    write(300, *) ilat, ilon, nlev
    write(300, '(74F20.12)') (press(ilon, ilat, i), i=1,nlev)
    write(300, '(74F20.12)') (temp(ilon, ilat, i), i=1,nlev)
    write(300, '(74E20.12)') (QVAPOR(ilon, ilat, i), i=1,nlev)
    write(300, '(74E20.12)') (QCLOUD(ilon, ilat, i), i=1,nlev)
    write(300, '(74E20.12)') (QRAIN(ilon, ilat, i), i=1,nlev)
    write(300, '(74E20.12)') (QICE(ilon, ilat, i), i=1,nlev)
    write(300, '(74E20.12)') (QSNOW(ilon, ilat, i), i=1,nlev)
    write(300, '(74E20.12)') (QGRAUP(ilon, ilat, i), i=1,nlev)
    call flush(300)
  end if

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

! dellaocate
  deallocate(atminp)

end subroutine read_text_profile
