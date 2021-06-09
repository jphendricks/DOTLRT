!=======================================================================
  subroutine construct_precip_profile ()
!=======================================================================      
!  Constructs a cloud  
!
! Modifications:
!  3/16/2021 Kevin Schaefer created subroutine
!----------------------------------------------------------------------
  use scan_Variables
  use dotlrt_variables
  use profiles
  use dotlrt_output

  implicit none

! local variables
  integer ilev            ! (-) atm level index
  integer ngau            ! (-) number slices for Gaussian intergral
  integer igau            ! (-) Gaussian intergral index
  integer indx_top        ! (-) layer index top cloud/pecip layer
  integer indx_bot        ! (-) layer index bottom cloud/pecip layer
  real(8) del_h           ! (-) raw gaussian value top of layer
  real(8) loc_h           ! (-) raw gaussian value bottom of layer
  real(8) frac_gau(nlev)  ! (-) raw gaussian fraction per layer
  real(8) frac_pre(nlev)  ! (-) precip fraction of total column precip per layer
  real(8) frac_cld(nlev)  ! (-) cloud fraction of total column precip per layer
  real(8) frac_clw(nlev)  ! (-) clw fraction of total column cloud per layer
  real(8) frac_rain(nlev) ! (-) rain fraction of total column cloud per layer
  real(8) frac_ice(nlev)  ! (-) ice fraction of total column cloud per layer
  real(8) frac_snow(nlev) ! (-) snow fraction of total column cloud per layer
  real(8) frac_grpl(nlev) ! (-) graupel fraction of total column cloud per layer
  real(8) hgt_var_pre     ! (km) precip height var
  real(8) hgt_var_cld     ! (km) cloud height var
  real(8) total           ! (-) local total
  real(8) slope           ! (km/K) slope of H vs temp
  integer frz_lev         ! (-) freezing level index
  real(8) hgt_frz         ! (km) freezing height
  real(8) hgt_mix_top     ! (km) top of mixed ice/liquid layer
  real(8) hgt_mix_bot     ! (km) bottom of mixed ice/liquid layer

! calculate cloud and precip characteristics
  cloud%precip%top = cloud%precip%pk + 2.d0*cloud%precip%std
  cloud%precip%bot = cloud%precip%pk - 2.d0*cloud%precip%std
  cloud%precip%bot = max (0.d0, atm(1)%hgt_bot, cloud%precip%bot)
  hgt_var_pre = cloud%precip%std*cloud%precip%std

  cloud%cloud%top = cloud%cloud%pk + 2.d0*cloud%cloud%std
  cloud%cloud%bot = cloud%cloud%pk - 2.d0*cloud%cloud%std
  cloud%cloud%bot = max (0.d0, atm(1)%hgt_bot, cloud%cloud%bot)
  hgt_var_cld = cloud%cloud%std*cloud%cloud%std

! find top and bottom of precip layer
  do ilev=1,nlev
    indx_bot=ilev
    if (atm(ilev)%hgt_top > cloud%precip%bot) exit
  enddo
  do ilev=nlev,1,-1
    indx_top=ilev
    if (atm(ilev)%hgt_bot < cloud%precip%top) exit
  enddo

! calculate raw, unscaled precip Gaussian fraction per layer
  frac_gau(:) = 0.d0
  total=0.d0
  ngau=100
  do ilev=indx_bot,indx_top
    if (ilev == indx_bot) then
      del_h= (atm(ilev)%hgt_top - cloud%precip%bot)/dble(ngau)
      loc_h = cloud%precip%bot
    elseif (ilev == indx_top) then
      del_h= (cloud%precip%top - atm(ilev)%hgt_bot)/dble(ngau)
      loc_h = atm(ilev)%hgt_bot
    else
      del_h= (atm(ilev)%hgt_top - atm(ilev)%hgt_bot)/dble(ngau)
      loc_h = atm(ilev)%hgt_bot
    endif

    frac_gau(ilev) = 0.d0
    do igau = 1, ngau
      frac_gau(ilev) = frac_gau(ilev)+ del_h*exp(-(loc_h-cloud%precip%pk)**2.d0/hgt_var_pre)
      loc_h = loc_h + del_h
    enddo
    if (frac_gau(ilev) < 10.d-15) frac_gau(ilev) = 0.d0

    total = total + frac_gau(ilev)
    !print*, ilev, cloud%precip%pk, cloud%precip%std, frac_gau(ilev)
   enddo

! calculate precip fraction of total column precipitation per layer
  do ilev=1,nlev
    frac_pre(ilev) = frac_gau(ilev)/total
    !if (ilev <20) print*, ilev, atm(ilev)%hgt, frac_gau(ilev), frac_pre(ilev)
  enddo

! find top and bottom of cloud layer
  do ilev=1,nlev
    indx_bot=ilev
    if (atm(ilev)%hgt_top > cloud%cloud%bot) exit
  enddo
  do ilev=nlev,1,-1
    indx_top=ilev
    if (atm(ilev)%hgt_bot < cloud%cloud%top) exit
  enddo

! calculate cloud raw, unscaled cloud Gaussian fraction per layer
  frac_gau(:) = 0.d0
  total=0.d0
  ngau=100
  do ilev=indx_bot,indx_top
    if (ilev == indx_bot) then
      del_h= (atm(ilev)%hgt_top - cloud%cloud%bot)/dble(ngau)
      loc_h = cloud%precip%bot
    elseif (ilev == indx_top) then
      del_h= (cloud%cloud%top - atm(ilev)%hgt_bot)/dble(ngau)
      loc_h = atm(ilev)%hgt_bot
    else
      del_h= (atm(ilev)%hgt_top - atm(ilev)%hgt_bot)/dble(ngau)
      loc_h = atm(ilev)%hgt_bot
    endif

    frac_gau(ilev) = 0.d0
    do igau = 1, ngau
      frac_gau(ilev) = frac_gau(ilev)+ del_h*exp(-(loc_h-cloud%cloud%pk)**2.d0/hgt_var_cld)
      loc_h = loc_h + del_h
    enddo
    if (frac_gau(ilev) < 10.d-15) frac_gau(ilev) = 0.d0

    total = total + frac_gau(ilev)
    !print*, ilev, cloud%cloud%pk, cloud%cloud%std, frac_gau(ilev)
  enddo

! calculate cloud fraction of total column cloud per layer
  do ilev=1,nlev
    frac_cld(ilev) = frac_gau(ilev)/total
    !if(ilev<20) print*, ilev, atm(ilev)%hgt, frac_gau(ilev), frac_cld(ilev)
  enddo

! calculate freeze height
  do ilev = 1, nlev
    if (atm(ilev)%temp < t_frz) exit
    frz_lev = ilev
  enddo
  slope = (atm(frz_lev)%hgt-atm(frz_lev+1)%hgt)/(atm(frz_lev)%temp-atm(frz_lev+1)%temp)
  hgt_frz=atm(frz_lev)%hgt+slope*(t_frz-atm(frz_lev)%temp)
  hgt_mix_top = hgt_frz + 1.d0
  hgt_mix_bot = hgt_frz - 1.d0

! calculate liquid and ice fractions of total cloud per layer 
  do ilev=1,nlev
    if (atm(ilev)%hgt < hgt_mix_bot) then
      frac_clw(ilev) = 1.d0
      frac_ice(ilev) = 0.d0
    elseif (atm(ilev)%hgt > hgt_mix_top) then
      frac_clw(ilev) = 0.d0
      frac_ice(ilev) = 1.d0
    else
      frac_ice(ilev) = 0.5d0*(atm(ilev)%hgt-hgt_frz)+0.5d0
      frac_clw(ilev)  = 1.d0 - frac_ice(ilev)
    endif
  enddo

! calculate liquid and ice fractions of total precip per layer 
  do ilev=1,nlev
    if (atm(ilev)%hgt < hgt_mix_bot) then
      frac_rain(ilev) = 1.d0
      frac_snow(ilev) = 0.d0
      frac_grpl(ilev) = 0.d0
    elseif (atm(ilev)%hgt > hgt_mix_top) then
      frac_rain(ilev) = 0.d0
      frac_snow(ilev) = 0.5d0
      frac_grpl(ilev) = 0.5d0
    else
      frac_rain(ilev) = frac_clw(ilev)
      frac_snow(ilev) = 0.5d0*frac_ice(ilev)
      frac_grpl(ilev) = 0.5d0*frac_ice(ilev)
    endif
  enddo

! calculate fraction of total column cloud water per species per layer
  do ilev=1,nlev
    atm(ilev)%clw%fprecip  = 0.d0
    atm(ilev)%rain%fprecip = frac_pre(ilev)*frac_rain(ilev)
    !if(ilev<20) print*, ilev, atm(ilev)%hgt, frac_pre(ilev), frac_rain(ilev), atm(ilev)%rain%fprecip

    atm(ilev)%ice%fprecip  = 0.d0
    atm(ilev)%snow%fprecip = frac_pre(ilev)*frac_snow(ilev)
    atm(ilev)%grpl%fprecip = frac_pre(ilev)*frac_grpl(ilev)

    atm(ilev)%clw%fcloud  = frac_cld(ilev)*frac_clw(ilev)
    atm(ilev)%rain%fcloud = 0.d0
    atm(ilev)%ice%fcloud  = frac_cld(ilev)*frac_ice(ilev)
    atm(ilev)%snow%fcloud = 0.d0
    atm(ilev)%grpl%fcloud = 0.d0
  enddo

! construct vertical density profiles
  call construct_hydromet_profile (cld%cloud%tot,  atm(:)%clw%fcloud,   atm(:)%clw%dens)
  call construct_hydromet_profile (cld%precip%tot, atm(:)%rain%fprecip, atm(:)%rain%dens)
  call construct_hydromet_profile (cld%cloud%tot,  atm(:)%ice%fcloud,   atm(:)%ice%dens)
  call construct_hydromet_profile (cld%precip%tot, atm(:)%snow%fprecip, atm(:)%snow%dens)
  call construct_hydromet_profile (cld%precip%tot, atm(:)%grpl%fprecip, atm(:)%grpl%dens)

  return                                                                    
  end

!=======================================================================
  subroutine construct_cloud_profile (top, bot)
!=======================================================================      
!  Constructs a cloud  
!
! Modifications:
!  3/16/2021 Kevin Schaefer created subroutine
!----------------------------------------------------------------------
  use scan_Variables
  use dotlrt_variables
  use profiles
  use dotlrt_output

  implicit none

! input variables
  real(8) top             ! (km) cloud top height  
  real(8) bot             ! (km) cloud bottom height  

! local variables
  integer ilev            ! (-) atm level index
  real(8) gaus_frac(nlev) ! (-) raw gaussian fraction per layer
  real(8) cld_frac(nlev)  ! (-) cloud fraction of total column cloud per layer
  real(8) clw_frac(nlev)  ! (-) clw fraction of total column cloud per layer
  real(8) ice_frac(nlev)  ! (-) ice fraction of total column cloud per layer
  real(8) thickness       ! (km) cloud thickness
  real(8) peak            ! (km) cloud density peak
  real(8) hgt_std         ! (km) cloud height standard deviation
  real(8) hgt_var         ! (km) cloud height var
  real(8) total           ! (-) local total
  real(8) slope           ! (km/K) slope of H vs temp
  integer frz_lev         ! (-) freezing level index
  real(8) hgt_frz         ! (km) freezing height
  real(8) hgt_mix_top     ! (km) top of mixed ice/liquid layer
  real(8) hgt_mix_bot     ! (km) bottom of mixed ice/liquid layer

! cloud thickness and standard deviation
  thickness = top - bot
  peak = (top + bot)/2.d0
  if (thickness <=0.d0) return
  hgt_std = thickness/6
  hgt_var = hgt_std*hgt_std

! calculate raw, unscaled Gaussian fraction per layer
  total=0.d0
  do ilev=1,nlev
    if(atm(ilev)%hgt <top .and. atm(ilev)%hgt > bot) then
      gaus_frac(ilev)= exp(-(atm(ilev)%hgt-peak)**2.d0/hgt_var)
    else
      gaus_frac(ilev)= 0.d0
    endif
    total = total + gaus_frac(ilev)
  enddo

! calculate cloud fraction of total column cloud per layer
  do ilev=1,nlev
    cld_frac(ilev) = gaus_frac(ilev)/total
  enddo
  
! calculate freeze height
  do ilev = 1, nlev
    if (atm(ilev)%temp < t_frz) exit
    frz_lev = ilev
  enddo
  slope = (atm(frz_lev)%hgt-atm(frz_lev+1)%hgt)/(atm(frz_lev)%temp-atm(frz_lev+1)%temp)
  hgt_frz=atm(frz_lev)%hgt+slope*(t_frz-atm(frz_lev)%temp)
  hgt_mix_top = hgt_frz + 1.d0
  hgt_mix_bot = hgt_frz - 1.d0

! calculate clw and ice fraction of total cloud water per layer 
  do ilev=1,nlev
    if (atm(ilev)%hgt < hgt_mix_bot) then
      clw_frac(ilev) = 1.d0
      ice_frac(ilev) = 0.d0
    elseif (atm(ilev)%hgt > hgt_mix_top) then
      clw_frac(ilev) = 0.d0
      ice_frac(ilev) = 1.d0
    else
      ice_frac(ilev) = 0.5d0*(atm(ilev)%hgt-hgt_frz)+0.5d0
      clw_frac(ilev) = 1.d0 - ice_frac(ilev)
    endif
  enddo

! calculate fraction of total column cloud water per species per layer
  do ilev=1,nlev
    atm(ilev)%clw%fcloud  = cld_frac(ilev)*clw_frac(ilev)
    atm(ilev)%rain%fcloud = 0.d0
    atm(ilev)%ice%fcloud  = cld_frac(ilev)*ice_frac(ilev)
    atm(ilev)%snow%fcloud = 0.d0
    atm(ilev)%grpl%fcloud = 0.d0
  enddo

! construct vertical density profiles
  call construct_hydromet_profile (cld%cloud%tot,  atm(:)%clw%fcloud,  atm(:)%clw%dens)
  call construct_hydromet_profile (cld%cloud%tot, atm(:)%rain%fcloud, atm(:)%rain%dens)
  call construct_hydromet_profile (cld%cloud%tot,  atm(:)%ice%fcloud,  atm(:)%ice%dens)
  call construct_hydromet_profile (cld%cloud%tot, atm(:)%snow%fcloud, atm(:)%snow%dens)
  call construct_hydromet_profile (cld%cloud%tot, atm(:)%grpl%fcloud, atm(:)%grpl%dens)

  return                                                                    
  end

!=======================================================================
  subroutine construct_hydromet_profile (col_tot, frac_lay, prof_den)
!=======================================================================      
!  Constructs an atmospheric profile of single hydrometeor  
!
! Modifications:
!  2/28/2021 Kevin Schaefer created subroutine
!----------------------------------------------------------------------
  use scan_Variables
  use dotlrt_variables
  use profiles
  use dotlrt_output

  implicit none

! input variables
  real(8) frac_lay(nlev) ! (-) fractions per layer of total colum hydromet mass
  real(8) col_tot        ! (g/m2) column total mass  

! output variables
  real(8) prof_den(nlev) ! (g/m3) hydromet density profile

! local variables
  integer ilev      ! (-) atm level index
  integer ibad      ! (-) bad value index
  integer nbad      ! (-) total number of bad layers
  real(8) mass_max(nlev) ! (g/m2) maximum mass per level
  real(8) mass_lay(nlev) ! (g/m2) mass per level
  real(8) mass_ext(nlev) ! (g/m2) extra mass per level
  real(8) x_mass    ! (g/m2) extra hydromet mass
  real(8) max_den   ! (g/m3) maximum possible hydromet density
  real(8) frac_tot  ! (g/m3) maximum possible hydromet density

! maximum allowed density
  max_den = 1d4

! calculate initial mass per layer
  do ilev=1,nlev
    mass_lay(ilev) = frac_lay(ilev)*col_tot
    mass_max(ilev) = max_den*atm(ilev)%hgt_del
  enddo

! redistribute mass if required
  do ibad = 1,10  ! Iterate a max of 10 times
    nbad=0
    x_mass = 0.d0
    frac_tot=0.d0
    do ilev=1,nlev
      mass_ext(ilev) = 0.
      if(mass_lay(ilev) > mass_max(ilev)) then
        mass_ext(ilev) = mass_lay(ilev) - mass_max(ilev)
        mass_lay(ilev) = mass_max(ilev)
        x_mass = x_mass + mass_ext(ilev)
        nbad = nbad+1.d0
        frac_tot = frac_tot + frac_lay(ilev)
      endif
      !if (mass_lay(ilev) > 0.d0) print*, ilev, mass_lay(ilev), mass_max(ilev)
    enddo

    do ilev=1,nlev
      if(mass_lay(ilev) > mass_max(ilev)) then
        mass_lay(ilev) = mass_lay(ilev) + frac_lay(ilev)/frac_tot*x_mass
      endif
    enddo
    
    if (nbad == 0) exit
  enddo

! calculate density per layer
  do ilev=1,nlev
    prof_den(ilev) = mass_lay(ilev)/atm(ilev)%hgt_del
    !if (prof_den(ilev) > 0.d0) print*, ilev, prof_den(ilev), max_den
  enddo

  return                                                                    
  end
!
!=======================================================================
  subroutine assign_obs_variables ()
!=======================================================================      
! Constructs an atmospheric profile from generic cloud characteristics 
!
! Modifications:
!  2/13/2021 Kevin Schaefer created subroutine
!  2/21/2021 Kevin Schaefer added all hydrometeor species
!----------------------------------------------------------------------
  use scan_Variables
  use dotlrt_variables
  use profiles
  use dotlrt_output

  implicit none

! local variables
  integer ichan   ! (-) channel index
  integer ivar    ! (-) variable index
  real(8) col_tot ! (g/m2) local column total value
  logical do_prof ! (-) logical to construct vertical hydrometeor profile

! print Message
  print*, 'Create Psuedo Observation'
  print*, 'obs_x=', obs_x, 'obs_y=', obs_y

! calculate reduced dimension stuff
  call construct_reduced_dim_atm_profile()

  if (obs_x == 17 .or. obs_y == 17) then ! total column clw
    obs_den=typical(17)
    print*, typical(17)
    col_tot=10.d0**obs_den
    print*, 'clw', obs_den, col_tot
    call construct_hydromet_profile (col_tot, atm(:)%clw%ftot, atm(:)%clw%dens)
    call construct_reduced_dim_atm_profile()
  endif

  if (obs_x == 18 .or. obs_y == 18) then ! total column rain
    obs_den=typical(18)
    col_tot=10.d0**obs_den
    print*, 'rain', obs_den, col_tot
    call construct_hydromet_profile (col_tot, atm(:)%rain%ftot, atm(:)%rain%dens)
    call construct_reduced_dim_atm_profile()
  endif

  if (obs_x == 19 .or. obs_y == 19) then ! total column ice
    obs_den=typical(19)
    col_tot=10.d0**obs_den
    print*, 'ice', obs_den, col_tot
    call construct_hydromet_profile (col_tot, atm(:)%ice%ftot, atm(:)%ice%dens)
    call construct_reduced_dim_atm_profile()
  endif

  if (obs_x == 20 .or. obs_y == 20) then ! total column snow
    obs_den=typical(20)
    col_tot=10.d0**obs_den
    print*, 'snow', obs_den, col_tot
    call construct_hydromet_profile (col_tot, atm(:)%snow%ftot, atm(:)%snow%dens)
    call construct_reduced_dim_atm_profile()
  endif

  if (obs_x == 21 .or. obs_y == 21) then ! total column graupel
    obs_den=typical(21)
    col_tot=10.d0**obs_den
    print*, 'grpl', obs_den, col_tot
    call construct_hydromet_profile (col_tot, atm(:)%grpl%ftot, atm(:)%grpl%dens)
    call construct_reduced_dim_atm_profile()
  endif

  if (obs_x==22.or.obs_y==22) then  ! total column cloud
    obs_den=typical(22)
    col_tot=10.d0**obs_den
    print*, 'cloud', obs_den, col_tot
    call construct_hydromet_profile (col_tot, atm(:)%clw%fcloud,  atm(:)%clw%dens)
    call construct_hydromet_profile (col_tot, atm(:)%ice%fcloud,  atm(:)%ice%dens)
    call construct_reduced_dim_atm_profile()
  endif

  if (obs_x==23.or.obs_y==23) then ! total column precip
    obs_den=typical(23)
    col_tot=10.d0**obs_den
    print*, 'precip', obs_den, col_tot
    call construct_hydromet_profile (col_tot, atm(:)%rain%fprecip, atm(:)%rain%dens)
    call construct_hydromet_profile (col_tot, atm(:)%snow%fprecip, atm(:)%snow%dens)
    call construct_hydromet_profile (col_tot, atm(:)%grpl%fprecip, atm(:)%grpl%dens)
    call construct_reduced_dim_atm_profile()
  endif

  if (obs_x == 24 .or. obs_y == 24) then ! total column hydrometeor
    obs_den=typical(24)
    col_tot=10.d0**obs_den
    print*, 'hydro', obs_den, col_tot
    call construct_hydromet_profile (col_tot, atm(:)%clw%fhydro,  atm(:)%clw%dens)
    call construct_hydromet_profile (col_tot, atm(:)%rain%fhydro, atm(:)%rain%dens)
    call construct_hydromet_profile (col_tot, atm(:)%ice%fhydro,  atm(:)%ice%dens)
    call construct_hydromet_profile (col_tot, atm(:)%snow%fhydro, atm(:)%snow%dens)
    call construct_hydromet_profile (col_tot, atm(:)%grpl%fhydro, atm(:)%grpl%dens)
    call construct_hydromet_profile (col_tot, atm(:)%hydro%ftot,  atm(:)%hydro%dens)
    call construct_reduced_dim_atm_profile()
  endif
  
  do_prof = .false.
  if(obs_x==25) do_prof = .true.
  if(obs_y==25) do_prof = .true.
  if(obs_x==26) do_prof = .true.
  if(obs_y==26) do_prof = .true.
  if(obs_x==27) do_prof = .true.
  if(obs_y==27) do_prof = .true.
  if(obs_x==28) do_prof = .true.
  if(obs_y==28) do_prof = .true.
  if (do_prof) then
    cloud%cloud%pk = typical(25)
    cloud%cloud%std = typical(26)
    cloud%precip%pk = typical(27)
    cloud%precip%std = typical(28)
    print*, cloud%cloud%pk, cloud%cloud%std, cloud%precip%pk, cloud%precip%std
    call construct_precip_profile()
  endif

! recalculate reduced dimension stuff
  call construct_reduced_dim_atm_profile()
  call print_reduced_profile ()

! calculate observed brightness temperatures
  ivar=0
  do ichan=minchan, maxchan
    ivar=ivar+1
    call extract_channel(ichan)
    call construct_single_surf_ref()
    call mrt( )
    Tbo_obs(ivar,:) = Tbo_str_mat(1,:)
  enddo

  !call write_current_profile()
  !stop

  return                                                                    
  end
!
!=======================================================================
  subroutine print_reduced_profile ()
!=======================================================================      
! prints reduced profile characteristics 
!
! Modifications:
!  4/1/2021 Kevin Schaefer created subroutine
!----------------------------------------------------------------------
  use scan_Variables
  use dotlrt_variables
  use profiles
  use dotlrt_output

  implicit none

! local variables
  integer ilev  ! (-) level index
  character*100 fmt

! recalculate reduced dimension stuff
  call construct_reduced_dim_atm_profile()
  fmt = '(a6,2x,5(a10,2x))'
  print(fmt),'phase', 'Col_tot', 'log10(tot)', 'h_peak', 'h_ave', 'h_std'

  fmt = '(a6,2x,5(f10.3,2x))'
  print(fmt),'clw',  cld%clw%tot,  log10(cld%clw%tot), cld%clw%h_max, cld%clw%h_ave, cld%clw%h_std
  print(fmt),'rain', cld%rain%tot, log10(cld%rain%tot), cld%rain%h_max, cld%rain%h_ave, cld%rain%h_std
  print(fmt),'ice',  cld%ice%tot,  log10(cld%ice%tot), cld%ice%h_max, cld%ice%h_ave, cld%ice%h_std
  print(fmt),'snow', cld%snow%tot, log10(cld%snow%tot), cld%snow%h_max, cld%snow%h_ave, cld%snow%h_std
  print(fmt),'grpl', cld%grpl%tot, log10(cld%grpl%tot), cld%grpl%h_max, cld%grpl%h_ave, cld%grpl%h_std
  print(fmt),'cloud', cld%cloud%tot, log10(cld%cloud%tot), cld%cloud%h_max, cld%cloud%h_ave, cld%cloud%h_std
  print(fmt),'precip', cld%precip%tot, log10(cld%precip%tot), cld%precip%h_max, cld%precip%h_ave, cld%precip%h_std
  print(fmt),'hydro',cld%hydro%tot,log10(cld%hydro%tot), cld%hydro%h_max, cld%hydro%h_ave, cld%hydro%h_std

  print*, 'Hydrometeor densities (g/m3)'
  fmt = '(9(a10,1x))'
  print(fmt), 'ht (km)', 'clw', 'rain', 'ice', 'snow', 'grpl', 'cloud', 'precip', 'hydro'
  fmt = '(9(f10.3,1x))'
  do ilev = 1, 20
    print(fmt), atm(ilev)%hgt,atm(ilev)%clw%dens, atm(ilev)%rain%dens, atm(ilev)%ice%dens, atm(ilev)%snow%dens, &
    atm(ilev)%grpl%dens, atm(ilev)%cloud%dens, atm(ilev)%precip%dens, atm(ilev)%hydro%dens
  enddo

  return                                                                    
  end

!=======================================================================
  subroutine cost_function ()
!=======================================================================      
! calculates cost function
!
! Modifications:
!   2/16/2021 Kevin Schaefer created subroutine
!----------------------------------------------------------------------
  use scan_Variables
  use dotlrt_variables
  use profiles
  use dotlrt_output

  implicit none

! local variables
  integer ivar  ! (-) variable index
  integer ichan ! (-) channel index
  real(8) del1  ! (K) difference between observed simulated brightness temp hor pol
  real(8) del2  ! (K) difference between observed simulated brightness temp vert pol

! cost function
  ivar=0
  cost_tot=0.d0
  do ichan=minchan, maxchan
      ivar=ivar+1
      del1 = Tbo_obs(ivar,1)-Tbo_sim(ivar,1)
      del2 = Tbo_obs(ivar,2)-Tbo_sim(ivar,2)
      cost_chan(ivar) = del1*del1 + del2*del2
      cost_tot = cost_tot + cost_chan(ivar)
  enddo
  return                                                                    
  end
!
!=======================================================================
      subroutine ZenithAngle (len, sind, cosd, &
        latsib, lonsib, tofday, cosz, test, test10)
!=======================================================================      
! calculates the zenith angle for an array of lat/lon points.
! The cosH variable accounts for the difference between GMT and local time
!
! Modifications:
!   Kevin Schaefer created ZenithAngle subroutine (5/31/01)
!   Kevin Schaefer corrected hour angle calculation (5/31/01)
!   Kevin Schaefer removedminimum cosz limit (7/21/01)
!----------------------------------------------------------------------
!
      implicit none
!
! input variables
      integer len      ! number of points
      real(8) sind        ! sine of solar declination
      real(8) cosd        ! cosine of solar declination
      real(8) latsib(len) ! SiB point latitude
      real(8) lonsib(len) ! SiB point longitude
      real(8) tofday      ! time of day (GMT, in hours)
      real(8) test, test10(10)     ! test variables
!
! output variables
      real(8) cosz(len)   ! cosine of Solar zenith angle
!
! internal variables
      real(8) pid180      ! conversion from degrees to radians
      real(8) sinlat      ! sine of latitude
      real(8) coslat      ! cosine of latitude
      real(8) HrAng       ! hour angle; longitude of Sun from Greenwhich meridian
      real(8) cosH        ! cosine delta longitude between Sun & SiB point
      integer i        ! sib point index
!
! calculate conversion from degrees to radians
      pid180=3.14159/180.
!
! Calculate hour angle (longitude of Sun from Greenwhich meridian)
      HrAng=(12.-tofday)*360./24.
!
! Calculate cosine of solar zenith angle for each SiB point
      do i = 1,len
        cosH=cos(pid180*(HrAng-lonsib(i)))
        sinlat=sin(pid180*latsib(i))
        coslat=cos(pid180*latsib(i))
        cosz(i)=coslat*cosd*cosH+sinlat*sind
      enddo
      test = cosz(1)
      test10(1) = cosz(1)
!
      return                                                                    
      end
!
!=======================================================================
      subroutine SunriseSunset
!=======================================================================      
! Calculates sunrise, sunset, noon, and midnight
!----------------------------------------------------------------------
!
      use scan_Variables
      use dotlrt_variables
      use profiles
      use dotlrt_output
!
      implicit none
!
! local variables
      real(8) tempval
      real(8) nightlat ! latitude of 100% darkness
      real(8) daylat   ! latitude of 100% sunlight
      real(8) lenday   ! length of day
!
! calculate solar declination
      dec=23.5*sin(1.72e-2*(DOY-80))
      sind=sin(pi*dec/180.)
      cosd=cos(pi*dec/180.)
!
! Calculate noon and midnight
      noon=12.-24.*lon/360.
      midnight=12.-24.*(lon-180.)/360.
      if(midnight.gt.24.) midnight=midnight-24.
!
! Calculate sunrise and sunset
      nightlat=atan(-1./tan(dec*pi/180.))*180./pi
      daylat=atan(1./tan(dec*pi/180.))*180./pi
      if (dec*nightlat.lt.dec*lat.and.dec*lat.lt.dec*daylat) then
        tempval=-tan(lat*pi/180.)*tan(dec*pi/180.)
        sunrise=12.-24./360.*(lon-180./pi*acos(tempval))
        sunset=12.-24./360.*(lon+180./pi*acos(tempval))
      else if (dec*nightlat.ge.dec*lat) then
        tempval=1.
        sunrise=noon
        sunset=noon
      else if (dec*lat.ge.dec*daylat) then
        tempval=-1.
        sunrise=midnight
        sunset=midnight
      endif
      if(sunrise.gt.24.) sunrise=sunrise-24.
      if(sunset.lt.0.) sunset=sunset+24.
      lenday=24./pi*acos(tempval)
!
      return
      end
!
!=======================================================================
      SUBROUTINE ABSORB_L93(FREQY,TEMP,PRES,RH,ATM_EXT)
!=======================================================================
! Converts retrieval quantities to parameters needed by Liebe93
!
! Modifications:
!  Kevin Schaefer changed to implicit none (2/21/02)
!
      IMPLICIT NONE
!
      real(8) FREQY    ! (GHz) frequency of radiation
      real(8) TEMP     ! (K) air temperature
      real(8) PRES     ! (mb) air pressure
      real(8) RH       ! (%) relative humidity
      real(8) atm_ext  ! (-) extinction coefficient
      real(8) T        ! (K) air temperature
      real(8) TS       ! (K) reference air temperature (?)
      real(8) EWS
      real(8) ESLOG
      real(8) ES
      real(8) WS       ! (kPa) saturation water vapor pressure
      real(8) W        
      real(8) E        ! water vapor partial pressure [kPA]
      real(8) Tc       ! temperature [C]
      real(8) Pb       ! total pressure [KPa]
      real(8) R        ! cloud water content [g/m^3]
!
! calculate water vapor partial pressure [kPA]
      T=TEMP
      TS=373.16
      EWS=1013.246
      ESLOG=-7.90298*(TS/T-1)+5.02808*LOG10(TS/T) &
       -1.3816E-07*(10**(11.344*(1-T/TS))-1) &
       +8.1328E-03*(10**(-3.49149*(TS/T-1))-1)+LOG10(EWS)
      ES=10**ESLOG
      WS=0.622*1000*ES/PRES
      W=(RH/100.)*WS
      E=(W/(W+0.622e3)*0.1*PRES)
!
! convert temperature to centegrade
      Tc=TEMP-273.15
!
! calculate total atmospheric pressure
      Pb=0.1*PRES+E
!
! set cloud water content
      R=0.
!
! calculate extinction coefficient
      CALL MPM93 (FREQY, Pb, E, Tc, R, ATM_EXT)
!
      RETURN
      END
!
!=======================================================================
      SUBROUTINE MPM93(F, Pbkpa, Ekpa, Tc, W, ABSCOF)
!======================================================================= 
! Computes volume absorption coefficient for an atmospheric
! layer given the meteorological properties. The allowed frequency
! range is from 1 to 1000 GHz.  This routine is hacked from Liebe's
! GAS1 subroutine in his MPM93 model, taking out rain and dispersion
! computations.  Included is dry air attenuation, oxygen and "psuedo"
! water vapor line-continuum absorption, and Rayleigh cloud droplet
! absorption. 
!    Parameters:
!      F       frequency (GHz)
!      Pbkpa   total pressure (kPa)
!      Ekpa    water vapor pressure (kPa)
!      Tc      temperature (C)
!      W       cloud liquid water content (g/m^3)
!      ABSCOF  absorption coefficient (km^-1)
!--------------------------------------------------------------------- 
!  MPM93 - subroutines adapted by Jeff Haferman (NASA/GSFC 5/97)
!  from Liebe's MPM93 model.  His comments are included below.
!  I've based this adaptation on Frank Evans' MPM92 extraction.
!---------------------------------------------------------------------
!       ATMOSPHERIC ATTENUATION AND DELAY RATES UP TO 1000 GHz
!       June 1993
!       Hans J. Liebe     (303-497-3310)    
!       George A. Hufford (       -3457)
!       Michael G. Cotton (       -7346)
!       Institute for Telecommunication Sciences
!       NTIA/ITS.S3 
!       325 BROADWAY
!       Boulder, CO  80303,  USA
!
!       FAX   :  (303) 497-5993 (ITS), 497-3680 (ITS.S2)
!       E-Mail:  HLIEBE@NTIA.ITS.BLDRDOC.GOV
!
! COMMENTS:
!   The Millimeter-wave Propagation Model (MPM85) was reported in Ref.
! [1]. Molecular absorption by O2, H2O, and N2 is considered, as well as
! dielectric loss for haze and fog/cloud conditions (Rayleigh absorption
! approximation), and dielectric plus scatter losses (aR**b -
! approximation to Mie's theory) under rain conditions. The complex
! atmospheric refractivity N (or path-specific rates of attenuation A and
! delay B) were continued to be upgraded as discussed in [2] - [7].
! 
! Features of the current version, MPM93, are:
! - Haze model to predict the water droplet density for 
!       U = 80 to 99.95%RH , when a hygroscopic aerosol reference density
!       wa(80%RH) and a climatic code ('A, B, C, or D') are provided [2],[3]   
! - Improved model for the dielectric properties of liquid water to
!       calculate RAYLEIGH absorption and delay by suspended water droplets
!       for haze, fog, and cloud conditions [6],[7]
! - Rain attenuation model for Laws & Parsons drop-sizes by Olsen et al. 
!       [11], and associated dispersive delay, approximated from results 
!       reported by Zuffery [12]
! - New temperature-dependent linewidth data (b3 to b6) for the water
!       vapor lines below 1 THz, and a 5 percent increase in the 
!       strength b1 of the 22-GHz and 183-GHz lines [9]
! - New set of line mixing coefficients (a5, a6) for dry air, and 
!       their improved fit to the extensive 60-GHz lab. data [8],[9]
! - Revised water vapor saturation pressure equation [10] 
! - Approximation for Zeeman (O2) [4] and Doppler (H2O) line-broadening
!       to cover heights up to 100 km.
! - New pseudo-line water vapor continuum formulation [9]   
! - Detailed treatment of the anisotropic, mesospheric Zeeman effect
!   of O2 microwave lines [5]. The ZPM  code [9].
! 
! REFERENCES
!  [1] H. Liebe, "An updated model for millimeter-wave propagation in
!       moist air", Radio Science, vol. 20, no. 5, pp. 1069-1089, 1985.
!  [2] H. Liebe,"A contribution to modeling atmospheric mm-wave properties",
!       FREQUENZ, vol.41, no. 1/2, pp. 31-36, 1987.
!  [3] H. Liebe and D. Layton, "MM-wave Properties of the Atmosphere:
!       Laboratory Studies and Propagation Modeling",
!       NTIA Report 87-224, 80p., Oct. 1987 (NTIS Order No. PB88-164215/AF).
!  [4] H. Liebe,"MPM89 - An atmospheric mm-wave propagation model",
!       Int. J. IR & MM Waves, vol.10, no.6, pp. 631-650, June 1989.
!  [5] G. Hufford and H. Liebe, "MM-Wave Propagation in the Mesosphere",
!       NTIA Report 89-249, 67p., Sept. 1989 (NTIS Order No. PB90-119868/AS).
!  [6] H. Liebe, T. Manabe, and G. Hufford, "Mm-wave attenuation and delay
!       rates due to fog/cloud conditions", IEEE Trans. Ant. Prop.,
!       vol. 37, no. 12, pp. 1617-1623, Dec. 1989.
!  [7] H. Liebe, G. Hufford (ice), and T. Manabe, "A model for the complex
!       refractivity of water (ice) at frequencies below 1 THz",
!       Int. J. IR & MM Waves, vol. 12, no. 7, 659-682, 1991.
!  [8] H. Liebe, P. Rosenkranz, and G. Hufford, "Atmospheric 60-GHz   
!       oxygen spectrum: New laboratory measurements and line parameters", 
!       J. Quant. Spectr. Rad. Transf., vol. 48, no. 5/6, pp. 629-643, 1992.
!  [9] H. Liebe, G. Hufford, and M. Cotton, "Propagation modeling of moist air 
!       and suspended water/ice particles at frequencies below 1000 GHz", 
!       Proc. AGARD Conf. Paper 3/1-10, Palma De Mallorca, Spain, May 1993.
! [10] W. Boegel, "Neue Naeherungsgleichungen fuer den Saettigungsdruck des
!       Wasserdampfes, DFVLR Bericht DLR-FB 77-52, 1977.
! [11] R.L. Olsen, D.V. Rogers, and D.B. Hodge, "The aRb relation in the
!       calculation of rain attenuation",
!       IEEE Trans. Ant. Prop., vol. AP-26, no. 2, pp. 318-329, 1978.
! [12] C.H. Zuffery, "A study of rain effects on EM waves in the
!       1 to 600 GHz range", MS-THesis, Dept. Electrical Eng.,
!       University of Colorado, Boulder,  CO 80309, Feb., 1972.
!
      IMPLICIT NONE
      real(8)    F, Pbkpa, Ekpa, Tc, W, ABSCOF
      INTEGER IFIRST, I, ICE

      real(8) AT1, AT2, AT3, AT4
      real(8) GAMMA, S, DELTA, So, GAMMAo, Sn
      real(8) GAMH, GAMD2, DELH
      real(8) fD, fS, Eps, Epinf, Eopt
      real(8) Ai, Bi, fice
      real(8) V, P, Pb, E
      
      COMPLEX(8) ZN, ZNw, ZEp, ZF, ZFo, ZFn

! Common block for oxygen and water vapor lines
      real(8)    F0O2(44), A(6,44)
      real(8)    F0H2O(35), B(6,35)
      real(8)    A1(44),A2(44),A3(44),A4(44),A5(44),A6(44)
      real(8)    B1(35),B2(35),B3(35),B4(35),B5(35),B6(35)
      ! COMMON /MWLINES2/ F0O2,A, F0H2O,B

      DATA IFIRST/0/
      DATA ICE/0/     ! hardcoded by JLH

!
!     The following data was in mwlines.93.data      
      data  F0O2 /  50.474239,   50.987747,   51.503349,   52.021412, &
                   52.542393,   53.066906,   53.595749,   54.130001,    &
                   54.671158,   55.221367,   55.783802,   56.264774, &
                   56.363388,   56.968204,   57.612484,   58.323875, &
                   58.446590,   59.164207,   59.590984,   60.306061, &
                   60.434776,   61.150558,   61.800156,   62.411217, &
                   62.486259,   62.997978,   63.568520,   64.127769, &
                   64.678902,   65.224068,   65.764771,   66.302094, &
                   66.836830,   67.369598,   67.900864,   68.431007, &
                   68.960312,  118.750343,  368.498352,  424.763123, &
                  487.249359,  715.393127,  773.839661,  834.145325 /
      data A1 /      0.094,    0.246,    0.608,    1.414,    3.102, &
                    6.410,   12.470,   22.800,   39.180,   63.160, &
                   95.350,   54.890,  134.400,  176.300,  214.100, &
                  238.600,  145.700,  240.400,  211.200,  212.400, &
                  246.100,  250.400,  229.800,  193.300,  151.700, &
                  150.300,  108.700,   73.350,   46.350,   27.480, &
                   15.300,    8.009,    3.946,    1.832,    0.801, &
                    0.330,    0.128,   94.500,    6.790,   63.800, &
                   23.500,    9.960,   67.100,   18.000  /
      data A2 /      9.694,    8.694,    7.744,    6.844,    6.004, &
                    5.224,    4.484,    3.814,    3.194,    2.624, &
                    2.119,    0.015,    1.660,    1.260,    0.915, &
                    0.626,    0.084,    0.391,    0.212,    0.212, &
                    0.391,    0.626,    0.915,    1.260,    0.083, &
                    1.665,    2.115,    2.620,    3.195,    3.815, &
                    4.485,    5.225,    6.005,    6.845,    7.745, &
                    8.695,    9.695,    0.009,    0.049,    0.044, &
                    0.049,    0.145,    0.130,    0.147 /
      data A3 /      0.890,    0.910,    0.940,    0.970,    0.990, &
                    1.020,    1.050,    1.070,    1.100,    1.130, &
                    1.170,    1.730,    1.200,    1.240,    1.280, &
                    1.330,    1.520,    1.390,    1.430,    1.450, &
                    1.360,    1.310,    1.270,    1.230,    1.540, &
                    1.200,    1.170,    1.130,    1.100,    1.070, &
                    1.050,    1.020,    0.990,    0.970,    0.940, &
                    0.920,    0.900,    1.630,    1.920,    1.930, &
                    1.920,    1.810,    1.820,    1.810 /
      data A4 /      0.000,    0.000,    0.000,    0.000,    0.000, &
                    0.000,    0.000,    0.000,    0.000,    0.000, &
                    0.000,    0.000,    0.000,    0.000,    0.000, &
                    0.000,    0.000,    0.000,    0.000,    0.000, &
                    0.000,    0.000,    0.000,    0.000,    0.000, &
                    0.000,    0.000,    0.000,    0.000,    0.000, &
                    0.000,    0.000,    0.000,    0.000,    0.000, &
                    0.000,    0.000,    0.000,    0.600,    0.600, &
                    0.600,    0.600,    0.600,    0.600 /
      data A5 /      0.240,    0.220,    0.197,    0.166,    0.136, &
                    0.131,    0.230,    0.335,    0.374,    0.258, &
                   -0.166,    0.390,   -0.297,   -0.416,   -0.613, &
                   -0.205,    0.748,   -0.722,    0.765,   -0.705, &
                    0.697,    0.104,    0.570,    0.360,   -0.498, &
                    0.239,    0.108,   -0.311,   -0.421,   -0.375, &
                   -0.267,   -0.168,   -0.169,   -0.200,   -0.228, &
                   -0.240,   -0.250,   -0.036,    0.000,    0.000, &
                    0.000,    0.000,    0.000,    0.000 /
      data A6 /      0.790,    0.780,    0.774,    0.764,    0.751, &
                    0.714,    0.584,    0.431,    0.305,    0.339, &
                    0.705,   -0.113,    0.753,    0.742,    0.697, &
                    0.051,   -0.146,    0.266,   -0.090,    0.081, &
                   -0.324,   -0.067,   -0.761,   -0.777,    0.097, &
                   -0.768,   -0.706,   -0.332,   -0.298,   -0.423, &
                   -0.575,   -0.700,   -0.735,   -0.744,   -0.753, &
                   -0.760,   -0.765,    0.009,    0.000,    0.000, &
                    0.000,    0.000,    0.000,    0.000 /
      data F0H2O / 22.235081,   67.803963,  119.995941,  183.310089, &
                 321.225647,  325.152924,  336.222595,  380.197357, &
                 390.134521,  437.346680,  439.150818,  443.018280, &
                 448.001068,  470.888947,  474.689117,  488.491119, &
                 503.568542,  504.482697,  547.676453,  552.020935, &
                 556.935974,  620.700806,  645.866150,  658.005310, &
                 752.033203,  841.053955,  859.962341,  899.306702, &
                 902.616150,  906.207336,  916.171570,  923.118408, &
                 970.315002,  987.926758, 1780.000000 /
      data B1 /     0.01130,    0.00012,    0.00008,    0.24200,    &
                   0.00483,    0.14990,    0.00011,    1.15200, &
                   0.00046,    0.00650,    0.09218,    0.01976, &
                   1.03200,    0.03297,    0.12620,    0.02520, &
                   0.00390,    0.00130,    0.97010,    1.47700, &
                  48.74000,    0.50120,    0.00713,    0.03022, &
                  23.96000,    0.00140,    0.01472,    0.00605, &
                   0.00426,    0.01876,    0.83400,    0.00869, &
                   0.89720,   13.21000, 2230.00000 /
      data B2 /   2.143,    8.735,    8.356,    0.668,    6.181,     &
                 1.540,    9.829,    1.048,    7.350,    5.050, &
                 3.596,    5.050,    1.405,    3.599,    2.381, &
                 2.853,    6.733,    6.733,    0.114,    0.114, &
                 0.159,    2.200,    8.580,    7.820,    0.396, &
                 8.180,    7.989,    7.917,    8.432,    5.111, &
                 1.442,   10.220,    1.920,    0.258,    0.952 /
      data B3 /   2.811,    2.858,    2.948,    3.050,    2.303, &
                 2.783,    2.693,    2.873,    2.152,    1.845, &
                 2.100,    1.860,    2.632,    2.152,    2.355, &
                 2.602,    1.612,    1.612,    2.600,    2.600, &
                 3.210,    2.438,    1.800,    3.210,    3.060, &
                 1.590,    3.060,    2.985,    2.865,    2.408, &
                 2.670,    2.900,    2.550,    2.985,   17.620 /
      data B4 /   4.80,    4.93,    4.78,    5.30,    4.69,    4.85, &
                 4.74,    5.38,    4.81,    4.23,    4.29,    4.23, &
                 4.84,    4.57,    4.65,    5.04,    3.98,    4.01, &
                 4.50,    4.50,    4.11,    4.68,    4.00,    4.14, &
                 4.09,    5.76,    4.09,    4.53,    5.10,    4.70, &
                 4.78,    5.00,    4.94,    4.55,   30.50 /
      data B5 /   0.69,    0.69,    0.70,    0.64,    0.67,    0.68, &
                 0.69,    0.54,    0.63,    0.60,    0.63,    0.60, &
                 0.66,    0.66,    0.65,    0.69,    0.61,    0.61, &
                 0.70,    0.70,    0.69,    0.71,    0.60,    0.69, &
                 0.68,    0.33,    0.68,    0.68,    0.70,    0.70, &
                 0.70,    0.70,    0.64,    0.68,    2.00 /
      data B6 /   1.00,    0.82,    0.79,    0.85,    0.54,    0.74, &
                 0.61,    0.89,    0.55,    0.48,    0.52,    0.50, &
                 0.67,    0.65,    0.64,    0.72,    0.43,    0.45, &
                 1.00,    1.00,    1.00,    0.68,    0.50,    1.00, &
                 0.84,    0.45,    0.84,    0.90,    0.95,    0.53, &
                 0.78,    0.80,    0.67,    0.90,    5.00 /
!
!---------------------------------------------------------------------
!
      do i = 1, 44
        A(1,i) = A1(i)
        A(2,i) = A2(i)
        A(3,i) = A3(i)
        A(4,i) = A4(i)
        A(5,i) = A5(i)
        A(6,i) = A6(i)
      end do
      do i = 1, 35
        B(1,i) = B1(i)
        B(2,i) = B2(i)
        B(3,i) = B3(i)
        B(4,i) = B4(i)
        B(5,i) = B5(i)
        B(6,i) = B6(i)
      end do
         
! Only read in line data the first time called
!      IF (IFIRST.EQ.0) THEN
!        IFIRST = 1
!	CALL READLINES2
!      ENDIF

! Relative inverse temperature
      V=300./(Tc+273.15)
! This version inputs E.
! Note MPM93 has pressure in mb, whereas MPM92 uses kPA
      Pb = 10.*Pbkpa
      E  = 10.*Ekpa
      P=Pb-E
      IF(P.LT.0)THEN
        P=0.
        Pb=E
      ENDIF

! For OXYGEN
      ZN=CMPLX(0.,0.)
      DO 10 I=1,44
       GAMMA=0.
       S=A(1,I)*P*V**3*EXP(A(2,I)*(1.-V))*1.E-6
       GAMMA=A(3,I)*(P*V**(0.8-A(4,I))+1.1*E*V)*1.E-3
       GAMMA=(GAMMA**2+(25*0.6E-4)**2)**0.5
       DELTA=(A(5,I)+A(6,I)*V)*(P+E)*(V**0.8)*1.E-3
       ZF=F/F0O2(I)*(CMPLX(1.,-DELTA)/CMPLX(F0O2(I)-F,-GAMMA)- &
                    CMPLX(1.,DELTA)/CMPLX(F0O2(I)+F,GAMMA))
       ZN=ZN+S*ZF
10    CONTINUE

! OXYGEN LINE ABSORPTION  
! Cannot be less than 0.  
      AT1=.182*F*AIMAG(ZN)
      IF(AT1.LT.0.) AT1=0.
!
! DRY AIR CONTINUUM
      ZN=CMPLX(0.,0.)
      So=6.14E-5*P*V**2
      GAMMAo=0.56E-3*(P+E)*V**.8
      ZFo=-F/CMPLX(F,GAMMAo)
      Sn=1.40E-12*p**2*V**3.5
      ZFn=CMPLX(0.,F/(1.93E-5*F**1.5+1.))
      ZN=So*ZFo+Sn*ZFn

! NONRESONAT DRY AIR ABSORPTION 
      AT2=.182*F*AIMAG(ZN)
! 
! WATER VAPOR
      ZN=CMPLX(0.,0.)
      DO 20 I=1,35
       GAMH=0.
       S=B(1,I)*E*V**3.5*EXP(B(2,I)*(1.-V))  
! Doppler approximation.
       GAMH=B(3,I)*(P*V**B(5,I)+B(4,I)*E*V**B(6,I))*1.E-3
       GAMD2=1E-12/V*(1.46*F0H2O(I))**2
       GAMH=0.535*GAMH+(0.217*GAMH**2+GAMD2)**0.5
       DELH=0.
       ZF=F/F0H2O(I)*(CMPLX(1.,-DELH)/CMPLX(F0H2O(I)-F,-GAMH)- &
                   CMPLX(1.,DELH)/CMPLX(F0H2O(I)+F,GAMH))
       ZN=ZN+S*ZF
20           CONTINUE

! WATER VAPOR LINE ABSORPTION 
! SEE LIEBE'S COMMENT REGARDING "PSUEDO-LINE WATER VAPOR CONTINUUM" - JLH
      AT3=.182*F*AIMAG(ZN)

! 
! LIQUID WATER PERMITTIVITY [8]
! Use exponential form for gamma for T<0 extrapolation (a la Frank Evans)
      IF(ICE.EQ.0)THEN
!JLH    fD=20.20-146.4*(V-1)+316*(V-1)**2
        fD=20.1*exp(7.88*(1-V))
        fS=39.8*fD 
        Eps=103.3*(V-1)+77.66
        Epinf=0.0671*Eps
        Eopt=3.52
! Complex Permittivity of water (double-Debye model)
        ZEp=Eps-f*((Eps-Epinf)/CMPLX(f,fD)+(Epinf-Eopt)/CMPLX(f,fS))
!
! ICE PERMITTIVITY [8]
      ELSE
        Ai=(62.*V-11.6)*1.E-4*EXP(-22.1*(V-1.))
        Bi=.542E-6*(-24.17+116.79/V+(V/(V-.9927))**2)
        Eps=3.15
! Complex Permittivity of Ice 
        fice=f
        IF(f.LT..001)fice=.001
        ZEp=CMPLX(3.15,Ai/fice+Bi*fice)
      ENDIF
! SUSPENDED PARTICLE RAYLEIGH APPROXIMATION [6]
      ZNw=1.5*W*((ZEp-1.)/(ZEp+2.)-1.+3./(Eps+2.))
!

! SUSPENDED WATER DROPLET EXTINCTION 
      AT4=.182*F*AIMAG(ZNw)

      ABSCOF=0.23026*(AT1+AT2+AT3+AT4)
!
      RETURN
      END 
!
!
!=======================================================================
      subroutine Plank_wavelength
!=======================================================================      
! Calculates plancks function at specified wavelength & wavelength interval
!----------------------------------------------------------------------
!
      use scan_Variables
      use dotlrt_variables
      use profiles
      use dotlrt_output
!
      implicit none
!
      real(8) demon ! denominator of planck's function
!
! calculate denominator of planck's function
      demon=exp(0.014403/Tplank/wave)-1.
!
! calculate planck's function
      black=pi*1.19e-16/wave**5./demon*Dwave
!
      return
      end
!
!=======================================================================
      subroutine Plank_frequency
!=======================================================================      
! Calculates plancks function at specified frequency
!----------------------------------------------------------------------
!
      use scan_Variables
      use dotlrt_variables
      use profiles
      use dotlrt_output
!
      implicit none
!
      real(8) tempval ! temporary variable
      real(8) MegaHz ! input frequency converted to megahertz
!
      MegaHz=freq*1.e-9
!
      tempval=exp(plank*freq/Boltz/Tplank)-1.

      black=1.473e-23*MegaHz*MegaHz*MegaHz/tempval
!
      return
      end
!
!
!=======================================================================
      SUBROUTINE dielectric_constant(depth,satfrac, sandfrac,kroot,mass_org,&
      root_depth, org_depth,rho_om_max, poros_om,test, test10)
!=======================================================================
! Calculate volumetric water content, porosity, 
! and dielectric constant as a function of depth
!
! Modifications:
!  Kevin Schaefer made subroutine (11/29/18)
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
! Input variables
    real(8) depth      ! (m) depth from surface
    real(8) satfrac    ! (-) saturation fraction of soil porosity
    real(8) sandfrac   ! (%) sand fraction of soil texture
    real(8) kRoot      ! (1/m) Exp const for decrease root-zone organic matter w/ depth
    real(8) mass_org   ! (kg/m2) mass of organic matter in top 1 m of soil
    real(8) root_depth ! (m) maximum rooting depth
    real(8) org_depth  ! (m) thickness of organic layer
    real(8) rho_om_max ! (kg/m3) maximum organic matter density
    real(8) poros_om   ! (-) soil porosity for pure organic soil
!
! Output variables
    real(8) test, test10(10) ! test variables
!
! Local variables
    real(8) vwc_org_exp  ! (-) vwc mix of organic and mineral soil
    real(8) vwc_min      ! (-) vwc pure mineral soil
    real(8) vwc_org_lay  ! (-) vwc organic soil layer
    real(8) vwc_perc     ! (%) vwc in percent
    real(8) poros        ! (-) porosity for organic mineral soil mix
    real(8) poros_min    ! (-) mineral soil porosity
    real(8) apor_min     ! (-) ave porosity pure mineral soil
    real(8) apor_org_lay ! (-) ave porosity organic soil layer
    real(8) orgfrac      ! (-) organic soil fraction
    real(8) a, b, c, d   ! (-) dielectric curve fit coefficients
    real(8) di_org_exp   ! (-) dielectric constant mix of organic and mineral soil
    real(8) di_min       ! (-) dielectric constant pure mineral soil
    real(8) di_org_lay   ! (-) dielectric constant organic soil layer
!
! set curve fit constants
    a=-0.0001
    b=0.0208
    c=-0.2
    d=4.
!
! calculate mineral soil porosity
    poros_min=0.489-0.00126*Sandfrac
!
!-----------------------------------------------------------------------
! Di-electric constant for for pure mineral soil
!-----------------------------------------------------------------------
    vwc_min=satfrac*poros_min
    apor_min=poros_min
    vwc_perc=vwc_min*100.
    di_min=a*vwc_perc**3.+b*vwc_perc**2.+c*vwc_perc+d
!
!-----------------------------------------------------------------------
! Di-electric constant for exponential decrease in organic matter
!-----------------------------------------------------------------------
! organic soil fraction
    orgfrac=kRoot*mass_org/(1.-exp(-kRoot*root_depth))*exp(-kRoot*(depth))/rho_om_max
    if(orgfrac>1.) orgfrac=1.
    if(orgfrac<0.) orgfrac=0.
    if(depth>root_depth) orgfrac=0.
!
! soil porosity
    poros=(1.-orgfrac)*poros_min+orgfrac*poros_om
    if(poros>1.) poros=1.
!
! volume of water and pore space
    vwc_org_exp=poros*satfrac
    vwc_perc=vwc_org_exp*100.
    di_org_exp=(a*vwc_perc**3.+b*vwc_perc**2.+c*vwc_perc+d)
!
!-----------------------------------------------------------------------
! Di-electric constant for organic soil layer
!-----------------------------------------------------------------------
!
! Calculate total water volume and pore space
    if(depth<=org_depth) then
      vwc_org_lay=satfrac*poros_om
      apor_org_lay=poros_om
      vwc_perc=vwc_org_lay*100.
      di_org_lay=a*vwc_perc**3.+b*vwc_perc**2.+c*vwc_perc+d
    else
      vwc_org_lay=vwc_min
      apor_org_lay=poros_min
      vwc_perc=vwc_org_lay*100.
      di_org_lay=a*vwc_perc**3.+b*vwc_perc**2.+c*vwc_perc+d
    endif
!
! assign test variables
    test=di_org_lay
    test10(1)=di_org_exp
    test10(2)=di_org_lay
    test10(3)=di_min

    RETURN
    END
!
!=======================================================================
subroutine dtess_eau (len, pl, tl, ess, dtess)
!!=======================================================================
!eau_sat computes the saturation mixing ratio, vapor pressure, and saturation,
!and their derivatives with respect to temperature over water, ice and mixed-
!phase clouds. the parameterization is that used in the Community Climate Com-
!munity Climate Model at NCAR.
!Laura D. Fowler /slikrock (07-01-01).

!send comments to laura@atmos.colostate.edu.

!subroutines called:
!none.

!argument list variables:
!input arguments:
!----------------

use scan_Variables, only: hltm, rv
      
implicit none
      
integer, intent(in) :: len                !length of vector.
real(8), intent(in), dimension(len):: &
   pl,               &!pressure                                           (Pa).
   tl                 !temperature                                         (K).

!output arguments:
!-----------------

real(8), intent(out), dimension(len) :: &
   ess,              &!saturation vapor pressure                          (Pa).
   dtess              !derivative of es with respect to temperature     (Pa/K).

!local variables:

integer i

real(8), parameter:: twmin=173.16     !lowest allowed temperature boundary for water       (K).
real(8), parameter:: twmax=373.16     !highest allowed temperature boundary for water      (K).     
real(8), parameter:: timin=173.16     !lowest allowed temperature boundary for ice         (K).
real(8), parameter:: timax=273.16     !highest allowed temperature boundary for ice        (K).

real(8) tstl 
real(8) t0tl

real(8), dimension(len):: &
      esw ,     dtesw ,      esi ,     dtesi ,&
      esm ,     dtesm ,      tl0 ,            &
    wghtm 

!ccm parameterization of saturation vapor pressure over water and ice:

real(8), parameter:: &
   ps = 1013.246,                   &!reference pressure             (hPa).
   ts = 373.16,                     &!reference temperature            (K).
   t0 = 273.16,                     &!freezing temperature             (K)
   lsub = hltm+0.3336d+06, &!
   tbgmin   = 253.15d0,      &!
   tbgmax   = 273.15d0        !
  
real(8):: &
       e1 ,   e2 ,     f ,    f1 ,&
       f2 ,   f3 ,    f4 ,    f5 ,&
   lphase , term1 , term2 ,&
   term3     

!------------------------------------------------------------------------------

!initialization of different arrays:

tl0    = tl
esw    = 0.0d0
esi    = 0.0d0
esm    = 0.0d0
dtesw  = 0.0d0
dtesi  = 0.0d0
dtesm  = 0.0d0

ess    = 0.0d0
dtess  = 0.0d0

!saturation over water:

do i = 1, len

   tl0(i)    = max(twmin,tl0(i))
   tl0(i)    = min(twmax,tl0(i))
   tstl      = ts / tl0(i)
   e1        = 11.344*(1.0 - tl0(i)/ts)
   e2        = -3.49149*(tstl - 1.0)
   f1        = -7.90298*(tstl - 1.0)
   f2        = 5.02808*log10(tstl)
   f3        = -1.3816*(10.0**e1-1.0)/10000000.0
   f4        = 8.1328*(10.0**e2-1.0)/1000.0
   f5        = log10(ps)
   f         = f1 + f2 + f3 + f4 + f5

   esw(i)    = (10.0**f)*1.e+02
   esw(i)    = min(esw(i),pl(i)*0.9)
   dtesw(i)  = hltm*esw(i)/(rv*tl0(i)*tl0(i))

   ess(i)    = esw(i)
   dtess(i)  = dtesw(i)


!saturation over ice:

   if(tl0(i).lt.timax) then

      tl0(i)    = max(tl0(i),timin)
      t0tl      = t0 / tl0(i)
      term1     = 2.01889049/(t0tl)
      term2     = 3.56654*log(t0tl)
      term3     = 20.947031*(t0tl)

      esi(i)    = 575.185606e10*exp(-(term1 + term2 + term3))
      esi(i)    = min(esi(i),pl(i)*0.9)
      dtesi(i)  = lsub*esi(i)/(rv*tl0(i)*tl0(i))

      ess(i)    = esi(i)
      dtess(i)  = dtesi(i)

   endif

!interpolated saturation variables:

   if(tl0(i).lt.tbgmax .and. tl0(i).ge.tbgmin) then

      wghtm(i)  = (tl0(i)-tbgmin)/(tbgmax-tbgmin)
      lphase    = hltm*wghtm(i)+lsub*(1.-wghtm(i))
      esm(i)    = wghtm(i)*esw(i) + (1.-wghtm(i))*esi(i)
      esm(i)    = min(esm(i),pl(i)*0.9)
      dtesm(i)  = lphase*esm(i)/(rv*tl0(i)*tl0(i))

      ess(i)    = esm(i)
      dtess(i)  = dtesm(i)

   endif
enddo

end subroutine dtess_eau

