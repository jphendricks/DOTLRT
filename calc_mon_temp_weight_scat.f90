!========================================================================
SUBROUTINE calc_mon_temp_weight_scat( ifreq, Tb_inp, tau, tb_pl_inp, tb_mn_inp, dtb_pl_inp, dtb_mn_inp, Tbo_streams_inp) 
!========================================================================
! This routine:
!     1 - Calculates vertical and horizontal polarization monochromatic
!         temperature weighting vectors for absorbing and scattering
!         atmosphere using the discrete ordinate technique
!     2 - Iterative (perturbation) technique used originally
!     3 - Includes surface relection/emission and polarization coupling effects
!         (sea surface emissivity is modelled, but emissivity over land is not)
!     4 - Observation angle measured with respect to nadir in deg
!     5 - frequency in GHz
!     6 - geometric height in kilometers above surface level
!     7 - Assumes planar stratified model with specularly reflecting surface
!     8 - Stopping criterion determined by frac_error
!     9 - Extinction coefficient used to calculate opacity increments
! History:
!     1 - Original provided in PASCAL by Marian Klein and Albin Gasiewski
!         NOAA ETL/ET1 Microwave System Development Division
!     2 - Converted April 2003 from PASCAL to FORTRAN by Bill Otto
!         william.d.otto@noaa.gov NOAA/ETL SET
!         FORTRAN 90 Portland Group Compiler on Red Hat Linux
!     3 - Modified 31 July 2003 by Bob Weber
!         FORTRAN 90 COMPAQ Compiler on Windows 2000 / XP
!     9/26/2020 Kevin Schaefer deleted unused variables
!     1/27/2021 Kevin Schaefer cleaned up code
!-----------------------------------------------------------------------
  use dotlrt_variables
  implicit none

  integer ifreq  ! frequency index
  real(8), dimension(2) :: tau
  real(8) tb_pl_inp(0:nlev), tb_mn_inp(0:nlev)
  real(8) dtb_pl_inp(0:nlev,nvar), dtb_mn_inp(0:nlev,nvar)
  real(8) Tbo_streams_inp(nangover2)

  integer last_angle, next_angle, jud, k
  real(8) theta_obs, Tbplinp, Tbmninp, Tb_inp

  real(8) theta            ! observation view angle measured with respect to nadir
  REAL(8) cos_observ_angle ! secant of observation angle
  REAL(8) height           ! height of level in km at height
  REAL(8) pressure         ! pressure of level in mb at height}
  REAL(8) temp             ! temperature of level in K at height}
  REAL(8) humid            ! (g/m3) water vapor density of level at height}
  REAL(8) lp,lq,lk0,la0    ! hydrometeor liquid size distribution parameters}
  REAL(8) ip,iq,ik0,ia0    ! hydrometeor frozen size distribution parameters}
  REAL(8) cabs_o2          ! oxygen and nitrogen absorption in nepers/km at height}
  REAL(8) cabs_h2o         ! water vapor absorption in nepers/km at height}
  REAL(8) abs_cloud_liq    ! cloud liquid absorption in nepers/km at height}
  REAL(8) scat_cloud_liq   ! cloud liquid scattering in nepers/km at height}
  REAL(8) abs_cloud_ice    ! cloud ice absorption in nepers/km at height}
  REAL(8) scat_cloud_ice   ! cloud ice scattering in nepers/km at height}
  REAL(8) abs_cloud_rn     ! cloud rain absorption in nepers/km at height}
  REAL(8) scat_cloud_rn    ! cloud rain scattering in nepers/km at height}
  REAL(8) abs_cloud_snow   ! {cloud snow absorption in nepers/km at height}
  REAL(8) scat_cloud_snow  ! {cloud snow scattering in nepers/km at height}
  REAL(8) abs_cloud_grpl   ! cloud groupel absorption in nepers/km at height}
  REAL(8) scat_cloud_grpl  ! cloud groupel scattering in nepers/km at height}
  REAL(8) g_cloud_liq      ! cloud liquid particle scattering asymmetry at height}
  REAL(8) g_cloud_ice      ! cloud ice particle scattering asymmetry at height}
  REAL(8) g_cloud_rn       ! cloud rain particle scattering asymmetry at height}
  REAL(8) g_cloud_snow     ! cloud snow particle scattering asymmetry at height}
  REAL(8) g_cloud_grpl     ! cloud grpl particle scattering asymmetry at height}
  REAL(8) abs_tot          ! total atmospheric absorption in nepers/km at height}
  REAL(8) scat_tot         ! total atmospheric scattering in nepers/km at height}
  REAL(8) ext_tot          ! total atmospheric extinction in nepers/km at height}
  INTEGER next_level       ! level_type {next level of atmospheric input data}
  INTEGER last_level       ! level_type {last level of atmospheric input data}
  REAL(8) llw, nlw         ! level interpolation weights}
  REAL(8) law,naw          ! angle interpolation weights}
  real(8) frequency
  real(8) dabsn2_t, dabsn2_p
  real(8) dabh2o_t, dabh2o_v, dabh2o_p
  real(8) do2abs_t, do2abs_v,do2abs_p
  real(8) albedo
  integer phase, j
  real(8), DIMENSION(3) :: dhab, dhsc, dg
  REAL(8), external :: o2abs, absn2, abh2o

! Sanity Check
  if( nlev .le. 0 ) then
    print*,' number of levels = 0 '
    return
  end if

! check for valid direction of propagation
  theta = dabs(obs_theta)
  do while( theta .gt. 180.0d0 )
    theta = theta - 180.0d0
  end do
  if( dabs(theta - 90.0d0) .le. 0.1d0 ) then
    print*,' theta ~ 90 degrees '
    return
  end if

! linear interpolation coefficients for angle - start
  if( obs_theta <= 90.0d0 ) then      ! upwelling - looking down
    jud = 1
  else if( obs_theta > 90.0d0 ) then  ! downwelling - looking up
    jud = 2
  end if
  theta_obs = obs_theta

! Check for valid angle
  if( jud == 1 ) then
    if( theta_obs .lt. quad_ang(1) ) theta_obs = quad_ang(1)
    if( theta_obs .gt. quad_ang(nangover2) ) theta_obs = quad_ang(nangover2)
    ! Determine angle boundaries
    next_angle = 2
    do while( (quad_ang(next_angle) .lt. theta_obs ) .and. (next_angle .lt. nangover2) )
      next_angle = next_angle + 1
    end do
    last_angle = next_angle - 1
    law = (quad_ang(next_angle) - theta_obs) / (quad_ang(next_angle) - quad_ang(last_angle))
  else if ( jud == 2 ) then
    if( theta_obs .gt. quad_ang(nang) ) theta_obs = quad_ang(nang)
    if( theta_obs .lt. quad_ang(nangover2+1) ) theta_obs = quad_ang(nangover2+1)
    ! Determine level boundaries
    next_angle = 2
    do while( (quad_ang(nang - next_angle + 1) .gt. theta_obs ) .and. (next_angle .lt. nangover2) )
      next_angle = next_angle+1
    end do
    last_angle = next_angle - 1
    ! Interpolate over atmospheric information at 'height'
    law = (quad_ang(nang-next_angle+1) - theta_obs) / (quad_ang(nang-next_angle+1) - quad_ang(nang-last_angle+1))
  end if
  naw = 1.0d0 - law

! linear interpolation coefficients for angle - end
  frequency = passband_freq(ifreq)

! calculate absorption, scattering, extinction, and albedo profiles at levels
  call calc_tot_ext(frequency)

! check for valid observation height
  height = obs_height
  if( height < 0.0d0 ) height = 0.0d0
  if( height < ( atm(1)%hgt + 0.0001d0 ) ) then
    height = ( atm(1)%hgt + 0.0001d0 )
  end if
  if( height > ( atm(nlev)%hgt - 0.0001d0 ) ) then
    height = ( atm(nlev)%hgt - 0.0001d0 )
  end if

! determine level boundaries at height
  next_level = 2
  do while( ( atm(next_level)%hgt < height ) .and. ( next_level < nlev ) )
    next_level = next_level + 1
  end do
  last_level = next_level-1
  ! observation level at or above observation height
  !  0 <= obs_lev <= nlev
  obs_lev  = next_level

! calculate absorption, scattering, extinction, and albedo at height - start
! linear interpolation coefficients for height
  llw = (atm(next_level)%hgt - height) / (atm(next_level)%hgt - atm(last_level)%hgt)
  nlw = 1.0d0 - llw
      
! interpolate to height temperature, pressure, water vapor
  temp = llw * atm(last_level)%temp + nlw * atm(next_level)%temp
  pressure = exp( llw * log(atm(last_level)%press) + nlw * log(atm(next_level)%press) )
  humid = llw * atm(last_level)%humid + nlw * atm(next_level)%humid

! gaseous absorption at height
  cabs_o2 = o2abs( temp, pressure, humid, frequency, do2abs_t, do2abs_v,do2abs_p ) &
              + absn2( temp, pressure, frequency, dabsn2_t, dabsn2_p )
  cabs_h2o = abh2o( temp, pressure, humid, frequency, dabh2o_t, dabh2o_v, dabh2o_p )

! interpolate to height hydrometeor distribution parameters for all hydrometeor phases - start
! cloud liquid water------------------------------------------------------------
  phase = 1 ! cloud liquid
  lp  = 0.0d0
  lq  = 0.0d0
  lk0 = 0.0d0
  la0 = 0.0d0
  if( (atm(last_level)%clw_k0 .gt. 0.0d0) .and. (atm(next_level)%clw_k0 .gt. 0.0d0)) then
    lp  = llw * atm(last_level)%clw_p + nlw * atm(next_level)%clw_p
    lq  = llw * atm(last_level)%clw_q + nlw * atm(next_level)%clw_q
    lk0 = llw * atm(last_level)%clw_k0 + nlw * atm(next_level)%clw_k0
    la0 = llw * atm(last_level)%clw_a0 + nlw * atm(next_level)%clw_a0
  else
    if( atm(last_level)%clw_k0 .gt. 0.0d0 ) then
      lp  = atm(last_level)%clw_p
      lq  = atm(last_level)%clw_q
      lk0 = llw * atm(last_level)%clw_k0
      la0 = atm(last_level)%clw_a0
    end if
    if( atm(next_level)%clw_k0 .gt. 0.0d0 ) then
      lp  = atm(next_level)%clw_p
      lq  = atm(next_level)%clw_q
      lk0 = nlw * atm(next_level)%clw_k0
      la0 = atm(next_level)%clw_a0
    end if
  end if
  call hydrometeor_master_5ph_d( frequency, phase, temp, lp, lq, lk0, la0, &
                                     abs_cloud_liq, scat_cloud_liq, g_cloud_liq,      &
                                     dhab, dhsc, dg,                                  &
                                     a0_const(next_level,phase) )

! rain------------------------------------------------------------
  phase = 2 ! rain
  ip  = 0.0d0
  iq  = 0.0d0
  ik0 = 0.0d0
  ia0 = 0.0d0
  if( (atm(last_level)%rain_k0 .gt. 0.0d0) .and. (atm(next_level)%rain_k0 .gt. 0.0d0)) then
    ip  = llw * atm(last_level)%rain_p + nlw * atm(next_level)%rain_p
    iq  = llw * atm(last_level)%rain_q + nlw * atm(next_level)%rain_q
    ik0 = llw * atm(last_level)%rain_k0 + nlw * atm(next_level)%rain_k0
    ia0 = llw * atm(last_level)%rain_a0 + nlw * atm(next_level)%rain_a0
  else
    if( atm(last_level)%rain_k0 .gt. 0.0d0 ) then
      ip  = atm(last_level)%rain_p
      iq  = atm(last_level)%rain_q
      ik0 = llw * atm(last_level)%rain_k0
      ia0 = atm(last_level)%rain_a0
    end if
    if(atm(next_level)%rain_k0 .GT. 0) THEN
      ip  = atm(next_level)%rain_p
      iq  = atm(next_level)%rain_q
      ik0 = nlw * atm(next_level)%rain_k0
      ia0 = atm(next_level)%rain_a0
    end if
  end if
  call hydrometeor_master_5ph_d( frequency, phase, temp, ip, iq, ik0, ia0, &
                                     abs_cloud_rn, scat_cloud_rn, g_cloud_rn,         &
                                     dhab, dhsc, dg,                                  &
                                     a0_const(next_level,phase) )

! ice------------------------------------------------------------
  phase = 3 ! ice
  ip  = 0.0d0
  iq  = 0.0d0
  ik0 = 0.0d0
  ia0 = 0.0d0
  if( (atm(last_level)%ice_k0 .gt. 0.0d0 ) .and. (atm(next_level)%ice_k0 .gt. 0.0d0) ) then
    ip  = llw * atm(last_level)%ice_p + nlw * atm(next_level)%ice_p
    iq  = llw * atm(last_level)%ice_q + nlw * atm(next_level)%ice_q
    ik0 = llw * atm(last_level)%ice_k0 + nlw * atm(next_level)%ice_k0
    ia0 = llw * atm(last_level)%ice_a0 + nlw * atm(next_level)%ice_a0
  else
    if( atm(last_level)%ice_k0 .gt. 0.0d0 ) then
      ip  = atm(last_level)%ice_p
      iq  = atm(last_level)%ice_q
      ik0 = llw * atm(last_level)%ice_k0
      ia0 = atm(last_level)%ice_a0
    end if
    if( atm(next_level)%ice_k0 .gt. 0.0d0 ) then
      ip  = atm(next_level)%ice_p
      iq  = atm(next_level)%ice_q
      ik0 = nlw * atm(next_level)%ice_k0
      ia0 = atm(next_level)%ice_a0
    end if
  end if
  call hydrometeor_master_5ph_d( frequency, phase, temp, ip, iq, ik0, ia0, &
                                     abs_cloud_ice, scat_cloud_ice, g_cloud_ice,      &
                                     dhab, dhsc, dg,                                  &
                                     a0_const(next_level,phase) )

! snow------------------------------------------------------------
  phase = 4 ! snow
  ip  = 0.0d0
  iq  = 0.0d0
  ik0 = 0.0d0
  ia0 = 0.0d0
  if( (atm(last_level)%snow_k0 .gt. 0.0d0) .and. (atm(next_level)%snow_k0 .gt. 0.0d0)) then
    ip  = llw * atm(last_level)%snow_p + nlw * atm(next_level)%snow_p
    iq  = llw * atm(last_level)%snow_q + nlw * atm(next_level)%snow_q
    ik0 = llw * atm(last_level)%snow_k0 + nlw * atm(next_level)%snow_k0
    ia0 = llw * atm(last_level)%snow_a0 + nlw * atm(next_level)%snow_a0
  else
    IF (atm(last_level)%snow_k0 .gt. 0.0d0) THEN
      ip  = atm(last_level)%snow_p
      iq  = atm(last_level)%snow_q
      ik0 = llw * atm(last_level)%snow_k0
      ia0 = atm(last_level)%snow_a0
    END IF
    if( atm(next_level)%snow_k0 .gt. 0.0d0) THEN
      ip  = atm(next_level)%snow_p
      iq  = atm(next_level)%snow_q
      ik0 = nlw * atm(next_level)%snow_k0
      ia0 = atm(next_level)%snow_a0
    end if
  end if
  call hydrometeor_master_5ph_d( frequency, phase, temp, ip, iq, ik0, ia0, &
                                     abs_cloud_snow, scat_cloud_snow, g_cloud_snow,   &
                                     dhab, dhsc, dg,                                  &
                                     a0_const(next_level,phase) )

! graupel------------------------------------------------------------
  phase = 5 ! graupel
  ip  = 0.0d0
  iq  = 0.0d0
  ik0 = 0.0d0
  ia0 = 0.0d0
  if( (atm(last_level)%grpl_k0 .gt. 0.0d0) .and. (atm(next_level)%grpl_k0 .gt. 0.0d0)) THEN
    ip  = llw * atm(last_level)%grpl_p + nlw * atm(next_level)%grpl_p
    iq  = llw * atm(last_level)%grpl_q + nlw * atm(next_level)%grpl_q
    ik0 = llw * atm(last_level)%grpl_k0 + nlw * atm(next_level)%grpl_k0
    ia0 = llw * atm(last_level)%grpl_a0 + nlw * atm(next_level)%grpl_a0
  else
    if( atm(last_level)%grpl_k0 .gt. 0.0d0) then
      ip  = atm(last_level)%grpl_p
      iq  = atm(last_level)%grpl_q
      ik0 = llw * atm(last_level)%grpl_k0
      ia0 = atm(last_level)%grpl_a0
    end if
    if( atm(next_level)%grpl_k0 .gt. 0.0d0) then
      ip  = atm(next_level)%grpl_p
      iq  = atm(next_level)%grpl_q
      ik0 = nlw * atm(next_level)%grpl_k0
      ia0 = atm(next_level)%grpl_a0
    end if
  end if
  call hydrometeor_master_5ph_d( frequency, phase, temp, ip, iq, ik0, ia0, &
                                     abs_cloud_grpl, scat_cloud_grpl, g_cloud_grpl,   &
                                     dhab, dhsc, dg,                                  &
                                     a0_const(next_level,phase) )

! total extinction
  abs_tot = cabs_o2 + cabs_h2o + abs_cloud_liq + abs_cloud_ice &
              + abs_cloud_rn + abs_cloud_snow + abs_cloud_grpl
  scat_tot = scat_cloud_liq + scat_cloud_ice + scat_cloud_rn &
               + scat_cloud_snow + scat_cloud_grpl
  ext_tot = abs_tot + scat_tot
  albedo = scat_tot / ext_tot

! calculate absorption, scattering, extinction, and albedo at height - end
! calculate opacity - start
! tau(1) = opacity from obs_height downward to surface at obs_theta
! tau(2) = opacity from obs_height upward              at obs_theta
  tau(1) = atm(1)%ext_tot * atm(1)%hgt ! opacity of first layer
  do j = 2, last_level
    tau(1) = tau(1) + atm(j)%ext_tot * ( atm(j)%hgt- atm(j-1)%hgt )
  end do
  tau(2) = tau(1) ! opacity from surface up to last level
  tau(1) = tau(1) + ext_tot * ( obs_height - atm(last_level)%hgt )
  do j = next_level, nlev
    tau(2) = tau(2) + atm(j)%ext_tot * ( atm(j)%hgt- atm(j-1)%hgt )
  end do
  tau(2) = tau(2) - tau(1)
  cos_observ_angle = cos( obs_theta * pi / 180.0d0 )
  if( dabs(cos_observ_angle) > 0.0d0 ) then
    tau(1) = tau(1) / dabs(cos_observ_angle)
    tau(2) = tau(2) / dabs(cos_observ_angle)
  end if

! radiative transfer solution
  call do_tb_gvh94()

! interpolate to height and angle brightness temperatures
  Tbplinp = law * ( llw * tb_pl(last_level,last_angle)   &
                      + nlw * tb_pl(next_level,last_angle) ) &
              + naw * ( llw * tb_pl(last_level,next_angle)   &
                      + nlw * tb_pl(next_level,next_angle) )
  Tbmninp = law * ( llw * tb_mn(last_level,last_angle)   &
                      + nlw * tb_mn(next_level,last_angle) ) &
              + naw * ( llw * tb_mn(last_level,next_angle)   &
                      + nlw * tb_mn(next_level,next_angle) )
  if( obs_theta <= 90.0d0 ) then      ! upwelling - looking down
    Tb_inp = Tbplinp
  else if( obs_theta > 90.0d0 ) then  ! downwelling - looking up
    Tb_inp = Tbmninp
  end if

! interpolate to height brightness temperatures, 08/21/2015, by K.Zhang
  do j = 1, nangover2
    if(obs_theta <= 90.0d0) then ! upwelling - looking down
      Tbo_streams_inp(j) = llw * tb_pl(last_level,j) + nlw * tb_pl(next_level,j)
    else ! downwelling - looking up
      Tbo_streams_inp(j) = llw * tb_mn(last_level,j) + nlw * tb_mn(next_level,j)
    end if
  end do

! interpolate profiles to observation angle
  do j = 0, nlev
    tb_pl_inp(j) = law * tb_pl(j,last_angle) + naw * tb_pl(j,next_angle)
    tb_mn_inp(j) = law * tb_mn(j,last_angle) + naw * tb_mn(j,next_angle)
    do k = 1, nvar
      dtb_pl_inp(j,k) = law * dtb_pl(j,last_angle,k) + naw * dtb_pl(j,next_angle,k)
      dtb_mn_inp(j,k) = law * dtb_mn(j,last_angle,k) + naw * dtb_mn(j,next_angle,k)
    end do
  end do

return
end subroutine calc_mon_temp_weight_scat
