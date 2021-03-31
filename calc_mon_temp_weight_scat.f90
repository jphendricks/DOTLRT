! calc_mon_temp_weight_scat.f90
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
SUBROUTINE calc_mon_temp_weight_scat( Tb_inp, tau, tb_pl_inp, tb_mn_inp, dtb_pl_inp, dtb_mn_inp ) !  frequency,
  use variables
  implicit none
  real(8), dimension(2) :: tau
  real(8) tb_pl_inp(0:atm_inp%num_levels), tb_mn_inp(0:atm_inp%num_levels)
  real(8) dtb_pl_inp(0:atm_inp%num_levels,nvar), dtb_mn_inp(0:atm_inp%num_levels,nvar)

  integer last_angle, next_angle, jud, k
  real(8) theta_obs, Tbplinp, Tbmninp, Tb_inp

  real(8) theta          ! observation view angle measured with respect to nadir
  REAL(8) cos_observ_angle       ! secant of observation angle
  REAL(8) height                      ! height of level in km at height
  REAL(8) pressure                    ! pressure of level in mb at height}
  REAL(8) temp                 ! temperature of level in K at height}
  REAL(8) vapor_density               ! water vapor density of level in g/m**3 at height}
  REAL(8) lp,lq,lk0,la0               ! hydrometeor liquid size distribution parameters}
  REAL(8) ip,iq,ik0,ia0               ! hydrometeor frozen size distribution parameters}
  REAL(8) cabs_o2                      ! oxygen and nitrogen absorption in nepers/km at height}
  REAL(8) cabs_h2o                     ! water vapor absorption in nepers/km at height}
  REAL(8) abs_cloud_liq               ! cloud liquid absorption in nepers/km at height}
  REAL(8) scat_cloud_liq              ! cloud liquid scattering in nepers/km at height}
  REAL(8) abs_cloud_ice               ! cloud ice absorption in nepers/km at height}
  REAL(8) scat_cloud_ice              ! cloud ice scattering in nepers/km at height}
  REAL(8) abs_cloud_rn                ! cloud rain absorption in nepers/km at height}
  REAL(8) scat_cloud_rn               ! cloud rain scattering in nepers/km at height}
  REAL(8) abs_cloud_snow              ! {cloud snow absorption in nepers/km at height}
  REAL(8) scat_cloud_snow             ! {cloud snow scattering in nepers/km at height}
  REAL(8) abs_cloud_grpl              ! cloud groupel absorption in nepers/km at height}
  REAL(8) scat_cloud_grpl             ! cloud groupel scattering in nepers/km at height}
  REAL(8) g_cloud_liq                 ! cloud liquid particle scattering asymmetry at height}
  REAL(8) g_cloud_ice                 ! cloud ice particle scattering asymmetry at height}
  REAL(8) g_cloud_rn                  ! cloud rain particle scattering asymmetry at height}
  REAL(8) g_cloud_snow                ! cloud snow particle scattering asymmetry at height}
  REAL(8) g_cloud_grpl                ! cloud grpl particle scattering asymmetry at height}
  REAL(8) abs_tot                     ! total atmospheric absorption in nepers/km at height}
  REAL(8) scat_tot                    ! total atmospheric scattering in nepers/km at height}
  REAL(8) ext_tot                     ! total atmospheric extinction in nepers/km at height}
  INTEGER next_level                 ! level_type {next level of atmospheric input data}
  INTEGER last_level                 ! level_type {last level of atmospheric input data}
  REAL(8) llw, nlw                    ! level interpolation weights}
  REAL(8) vr                          ! vertical surface reflectivity}
  REAL(8) hr                          ! horizontal surface reflectivity}
  REAL(8) law,naw                     ! angle interpolation weights}
  REAL(8) tauhd                       ! atmospheric opacity from "height" increasing downwards}
  REAL(8) dtauhd                      ! atmospheric opacity from "height" increment}
  REAL(8) tausu                       ! atmospheric opacity from "surface" increasing upwards}
  REAL(8) dtausu                      ! atmospheric opacity from "surface" increment}
  REAL(8) trans_sh_v                  ! transmittance to "height" from "surface" including reflection loss for vertical}
  REAL(8) trans_sh_h                  ! transmittance to "height" from "surface" including reflection loss for horizontal}
  REAL(8) tau_diff                    ! opacity difference}
  REAL(8) accum                       ! General accumulator}
  INTEGER iter_order                  ! Smallint
  REAL(8) max_albedo                  ! Double
  INTEGER max_iter_order              ! LongInt
  real(8)                 :: lat 
  real(8)                 :: lng 
  real(8)                 :: frequency
  REAL(8) term11,term12,term21,term22 ! Double
  REAL(8) theta_index,theta_p_index   ! Double
  REAL(8) temp_exp                    ! Double
  REAL(8) max_weight_incr             ! Double
  REAL(8) weight_error_to_comp_est    ! Double
  real(8) dabsn2_t, dabsn2_p
  real(8) dabh2o_t, dabh2o_v, dabh2o_p
  real(8) do2abs_t, do2abs_v,do2abs_p
  real(8) albedo
  integer phase, j
  real(8), DIMENSION(3) :: dhab, dhsc, dg
  REAL(8), external :: o2abs, absn2, abh2o

  if( atm_inp%num_levels .le. 0 ) then
    print*,' number of levels = 0 '
    return
  end if

  ! check for valid direction of propagation
  theta = dabs(inp_theta)
  do while( theta .gt. 180.0d0 )
    theta = theta - 180.0d0
  end do
  if( dabs(theta - 90.0d0) .le. 0.1d0 ) then
    print*,' theta ~ 90 degrees '
    return
  end if
  ! linear interpolation coefficients for angle - start
    if( inp_theta <= 90.0d0 ) then      ! upwelling - looking down
      jud = 1
    else if( inp_theta > 90.0d0 ) then  ! downwelling - looking up
      jud = 2
    end if
    theta_obs = inp_theta
    ! Check for valid angle
    if( jud == 1 ) then
      if( theta_obs .lt. teta(1) ) theta_obs = teta(1)
      if( theta_obs .gt. teta(nangover2) ) theta_obs = teta(nangover2)
      ! Determine angle boundaries
      next_angle = 2
      do while( (teta(next_angle) .lt. theta_obs ) .and. &
                (next_angle .lt. nangover2) )
           next_angle = next_angle + 1
      end do
      last_angle = next_angle - 1
      law = (teta(next_angle) - theta_obs) &
          / (teta(next_angle) - teta(last_angle))
    else if ( jud == 2 ) then
      if( theta_obs .gt. teta(nang) ) theta_obs = teta(nang)
      if( theta_obs .lt. teta(nangover2+1) ) theta_obs = teta(nangover2+1)
      ! Determine level boundaries
      next_angle = 2
      do while( (teta(nang - next_angle + 1) .gt. theta_obs ) .and. &
                (next_angle .lt. nangover2) )
           next_angle = next_angle+1
      end do
      last_angle = next_angle - 1
      ! Interpolate over atmospheric information at 'height'
      law = (teta(nang-next_angle+1) - theta_obs) &
          / (teta(nang-next_angle+1) - teta(nang-last_angle+1))
    end if
    naw = 1.0d0 - law
  ! linear interpolation coefficients for angle - end

  frequency = passband_freq(nfreq)
  call cpu_time(times(1))
  ! calculate absorption, scattering, extinction, and albedo profiles at levels
    call calc_tot_ext(frequency)

  ! check for valid observation height
    height = inp_height
    if( height < 0.0d0 ) height = 0.0d0
    if( height < ( atm_inp%prof(1)%height + 0.0001d0 ) ) then
        height = ( atm_inp%prof(1)%height + 0.0001d0 )
    end if
    if( height > ( atm_inp%prof(atm_inp%num_levels)%height - 0.0001d0 ) ) then
        height = ( atm_inp%prof(atm_inp%num_levels)%height - 0.0001d0 )
    end if

  ! determine level boundaries at height
    next_level = 2
    do while( ( atm_inp%prof(next_level)%height < height ) .and. &
              ( next_level < atm_inp%num_levels ) )
                next_level = next_level + 1
    end do
    last_level = next_level-1
  ! observation level at or above observation height
  !  0 <= obs_lev <= atm_inp%num_levels
    obs_lev  = next_level
    obs_lev1 = obs_lev * m1

  ! calculate absorption, scattering, extinction, and albedo at height - start
    ! linear interpolation coefficients for height
      llw = (atm_inp%prof(next_level)%height - height)                          &
          / (atm_inp%prof(next_level)%height - atm_inp%prof(last_level)%height)
      nlw = 1.0d0 - llw
    ! interpolate to height temperature, pressure, water vapor
      temp = llw * atm_inp%prof(last_level)%temperature      &
           + nlw * atm_inp%prof(next_level)%temperature
      pressure = exp( llw * log(atm_inp%prof(last_level)%pressure)  &
                    + nlw * log(atm_inp%prof(next_level)%pressure) )
      vapor_density = llw * atm_inp%prof(last_level)%vapor_density  &
                    + nlw * atm_inp%prof(next_level)%vapor_density
    ! gaseous absorption at height
      cabs_o2 = o2abs( temp, pressure, vapor_density, frequency, do2abs_t, do2abs_v,do2abs_p ) &
              + absn2( temp, pressure, frequency, dabsn2_t, dabsn2_p )
      cabs_h2o = abh2o( temp, pressure, vapor_density, frequency, dabh2o_t, dabh2o_v, dabh2o_p )

    ! interpolate to height hydrometeor distribution parameters for all hydrometeor phases - start
      phase = 1 ! cloud liquid
      lp  = 0.0d0
      lq  = 0.0d0
      lk0 = 0.0d0
      la0 = 0.0d0
      if( (atm_inp%prof(last_level)%cloud_liq_k0 .gt. 0.0d0) .and. &
          (atm_inp%prof(next_level)%cloud_liq_k0 .gt. 0.0d0)) then
        lp  = llw * atm_inp%prof(last_level)%cloud_liq_p &
            + nlw * atm_inp%prof(next_level)%cloud_liq_p
        lq  = llw * atm_inp%prof(last_level)%cloud_liq_q &
            + nlw * atm_inp%prof(next_level)%cloud_liq_q
        lk0 = llw * atm_inp%prof(last_level)%cloud_liq_k0 &
            + nlw * atm_inp%prof(next_level)%cloud_liq_k0
        la0 = llw * atm_inp%prof(last_level)%cloud_liq_a0 &
            + nlw * atm_inp%prof(next_level)%cloud_liq_a0
      else
        if( atm_inp%prof(last_level)%cloud_liq_k0 .gt. 0.0d0 ) then
          lp  = atm_inp%prof(last_level)%cloud_liq_p
          lq  = atm_inp%prof(last_level)%cloud_liq_q
          lk0 = llw * atm_inp%prof(last_level)%cloud_liq_k0
          la0 = atm_inp%prof(last_level)%cloud_liq_a0
        end if
        if( atm_inp%prof(next_level)%cloud_liq_k0 .gt. 0.0d0 ) then
          lp  = atm_inp%prof(next_level)%cloud_liq_p
          lq  = atm_inp%prof(next_level)%cloud_liq_q
          lk0 = nlw * atm_inp%prof(next_level)%cloud_liq_k0
          la0 = atm_inp%prof(next_level)%cloud_liq_a0
        end if
      end if
      call hydrometeor_master_5ph_d( frequency, phase, temp, lp, lq, lk0, la0, &
                                     abs_cloud_liq, scat_cloud_liq, g_cloud_liq,      &
                                     dhab, dhsc, dg,                                  &
                                     a0_is_constant(next_level,phase) )
!        hydro_prof(next_level,phase)%cloudab = abs_cloud_liq
!        hydro_prof(next_level,phase)%cloudsc = scat_cloud_liq
!        hydro_prof(next_level,phase)%cloudg = g_cloud_liq
!        hydro_prof(next_level,phase)%dcloudab_dt = dhab(1)
!        hydro_prof(next_level,phase)%dcloudab_dk0 = dhab(2)
!        hydro_prof(next_level,phase)%dcloudab_da0 = dhab(3)
!        hydro_prof(next_level,phase)%dcloudsc_dt = dhsc(1)
!        hydro_prof(next_level,phase)%dcloudsc_dk0 = dhsc(2)
!        hydro_prof(next_level,phase)%dcloudsc_da0 = dhsc(3)
!        hydro_prof(next_level,phase)%dcloudg_dt = dg(1)
!        hydro_prof(next_level,phase)%dcloudg_dk0 = dg(2)
!        hydro_prof(next_level,phase)%dcloudg_da0 = dg(3)
      phase = 2 ! rain
      ip  = 0.0d0
      iq  = 0.0d0
      ik0 = 0.0d0
      ia0 = 0.0d0
      if( (atm_inp%prof(last_level)%cloud_rn_k0 .gt. 0.0d0) .and. &
          (atm_inp%prof(next_level)%cloud_rn_k0 .gt. 0.0d0)) then
        ip  = llw * atm_inp%prof(last_level)%cloud_rn_p &
            + nlw * atm_inp%prof(next_level)%cloud_rn_p
        iq  = llw * atm_inp%prof(last_level)%cloud_rn_q &
            + nlw * atm_inp%prof(next_level)%cloud_rn_q
        ik0 = llw * atm_inp%prof(last_level)%cloud_rn_k0 &
            + nlw * atm_inp%prof(next_level)%cloud_rn_k0
        ia0 = llw * atm_inp%prof(last_level)%cloud_rn_a0 &
            + nlw * atm_inp%prof(next_level)%cloud_rn_a0
      else
        if( atm_inp%prof(last_level)%cloud_rn_k0 .gt. 0.0d0 ) then
          ip  = atm_inp%prof(last_level)%cloud_rn_p
          iq  = atm_inp%prof(last_level)%cloud_rn_q
          ik0 = llw * atm_inp%prof(last_level)%cloud_rn_k0
          ia0 = atm_inp%prof(last_level)%cloud_rn_a0
        end if
        if(atm_inp%prof(next_level)%cloud_rn_k0 .GT. 0) THEN
          ip  = atm_inp%prof(next_level)%cloud_rn_p
          iq  = atm_inp%prof(next_level)%cloud_rn_q
          ik0 = nlw * atm_inp%prof(next_level)%cloud_rn_k0
          ia0 = atm_inp%prof(next_level)%cloud_rn_a0
        end if
      end if
      call hydrometeor_master_5ph_d( frequency, phase, temp, ip, iq, ik0, ia0, &
                                     abs_cloud_rn, scat_cloud_rn, g_cloud_rn,         &
                                     dhab, dhsc, dg,                                  &
                                     a0_is_constant(next_level,phase) )
!        hydro_prof(next_level,phase)%cloudab = abs_cloud_rn
!        hydro_prof(next_level,phase)%cloudsc = scat_cloud_rn
!        hydro_prof(next_level,phase)%cloudg = g_cloud_rn
!        hydro_prof(next_level,phase)%dcloudab_dt = dhab(1)
!        hydro_prof(next_level,phase)%dcloudab_dk0 = dhab(2)
!        hydro_prof(next_level,phase)%dcloudab_da0 = dhab(3)
!        hydro_prof(next_level,phase)%dcloudsc_dt = dhsc(1)
!        hydro_prof(next_level,phase)%dcloudsc_dk0 = dhsc(2)
!        hydro_prof(next_level,phase)%dcloudsc_da0 = dhsc(3)
!        hydro_prof(next_level,phase)%dcloudg_dt = dg(1)
!        hydro_prof(next_level,phase)%dcloudg_dk0 = dg(2)
!        hydro_prof(next_level,phase)%dcloudg_da0 = dg(3)
      phase = 3 ! ice
      ip  = 0.0d0
      iq  = 0.0d0
      ik0 = 0.0d0
      ia0 = 0.0d0
      if( (atm_inp%prof(last_level)%cloud_ice_k0 .gt. 0.0d0 ) .and. &
          (atm_inp%prof(next_level)%cloud_ice_k0 .gt. 0.0d0) ) then
        ip  = llw * atm_inp%prof(last_level)%cloud_ice_p &
            + nlw * atm_inp%prof(next_level)%cloud_ice_p
        iq  = llw * atm_inp%prof(last_level)%cloud_ice_q &
            + nlw * atm_inp%prof(next_level)%cloud_ice_q
        ik0 = llw * atm_inp%prof(last_level)%cloud_ice_k0 &
            + nlw * atm_inp%prof(next_level)%cloud_ice_k0
        ia0 = llw * atm_inp%prof(last_level)%cloud_ice_a0 &
            + nlw * atm_inp%prof(next_level)%cloud_ice_a0
      else
        if( atm_inp%prof(last_level)%cloud_ice_k0 .gt. 0.0d0 ) then
          ip  = atm_inp%prof(last_level)%cloud_ice_p
          iq  = atm_inp%prof(last_level)%cloud_ice_q
          ik0 = llw * atm_inp%prof(last_level)%cloud_ice_k0
          ia0 = atm_inp%prof(last_level)%cloud_ice_a0
        end if
        if( atm_inp%prof(next_level)%cloud_ice_k0 .gt. 0.0d0 ) then
          ip  = atm_inp%prof(next_level)%cloud_ice_p
          iq  = atm_inp%prof(next_level)%cloud_ice_q
          ik0 = nlw * atm_inp%prof(next_level)%cloud_ice_k0
          ia0 = atm_inp%prof(next_level)%cloud_ice_a0
        end if
      end if
      call hydrometeor_master_5ph_d( frequency, phase, temp, ip, iq, ik0, ia0, &
                                     abs_cloud_ice, scat_cloud_ice, g_cloud_ice,      &
                                     dhab, dhsc, dg,                                  &
                                     a0_is_constant(next_level,phase) )
!        hydro_prof(next_level,phase)%cloudab = abs_cloud_ice
!        hydro_prof(next_level,phase)%cloudsc = scat_cloud_ice
!        hydro_prof(next_level,phase)%cloudg = g_cloud_ice
!        hydro_prof(next_level,phase)%dcloudab_dt = dhab(1)
!        hydro_prof(next_level,phase)%dcloudab_dk0 = dhab(2)
!        hydro_prof(next_level,phase)%dcloudab_da0 = dhab(3)
!        hydro_prof(next_level,phase)%dcloudsc_dt = dhsc(1)
!        hydro_prof(next_level,phase)%dcloudsc_dk0 = dhsc(2)
!        hydro_prof(next_level,phase)%dcloudsc_da0 = dhsc(3)
!        hydro_prof(next_level,phase)%dcloudg_dt = dg(1)
!        hydro_prof(next_level,phase)%dcloudg_dk0 = dg(2)
!        hydro_prof(next_level,phase)%dcloudg_da0 = dg(3)
      phase = 4 ! snow
      ip  = 0.0d0
      iq  = 0.0d0
      ik0 = 0.0d0
      ia0 = 0.0d0
      if( (atm_inp%prof(last_level)%cloud_snow_k0 .gt. 0.0d0) .and. &
          (atm_inp%prof(next_level)%cloud_snow_k0 .gt. 0.0d0)) then
        ip  = llw * atm_inp%prof(last_level)%cloud_snow_p &
            + nlw * atm_inp%prof(next_level)%cloud_snow_p
        iq  = llw * atm_inp%prof(last_level)%cloud_snow_q &
            + nlw * atm_inp%prof(next_level)%cloud_snow_q
        ik0 = llw * atm_inp%prof(last_level)%cloud_snow_k0 &
            + nlw * atm_inp%prof(next_level)%cloud_snow_k0
        ia0 = llw * atm_inp%prof(last_level)%cloud_snow_a0 &
            + nlw * atm_inp%prof(next_level)%cloud_snow_a0
      else
        IF (atm_inp%prof(last_level)%cloud_snow_k0 .gt. 0.0d0) THEN
          ip  = atm_inp%prof(last_level)%cloud_snow_p
          iq  = atm_inp%prof(last_level)%cloud_snow_q
          ik0 = llw * atm_inp%prof(last_level)%cloud_snow_k0
          ia0 = atm_inp%prof(last_level)%cloud_snow_a0
        END IF
        if( atm_inp%prof(next_level)%cloud_snow_k0 .gt. 0.0d0) THEN
          ip  = atm_inp%prof(next_level)%cloud_snow_p
          iq  = atm_inp%prof(next_level)%cloud_snow_q
          ik0 = nlw * atm_inp%prof(next_level)%cloud_snow_k0
          ia0 = atm_inp%prof(next_level)%cloud_snow_a0
        end if
      end if
      call hydrometeor_master_5ph_d( frequency, phase, temp, ip, iq, ik0, ia0, &
                                     abs_cloud_snow, scat_cloud_snow, g_cloud_snow,   &
                                     dhab, dhsc, dg,                                  &
                                     a0_is_constant(next_level,phase) )
!        hydro_prof(next_level,phase)%cloudab = abs_cloud_snow
!        hydro_prof(next_level,phase)%cloudsc = scat_cloud_snow
!        hydro_prof(next_level,phase)%cloudg = g_cloud_snow
!        hydro_prof(next_level,phase)%dcloudab_dt = dhab(1)
!        hydro_prof(next_level,phase)%dcloudab_dk0 = dhab(2)
!        hydro_prof(next_level,phase)%dcloudab_da0 = dhab(3)
!        hydro_prof(next_level,phase)%dcloudsc_dt = dhsc(1)
!        hydro_prof(next_level,phase)%dcloudsc_dk0 = dhsc(2)
!        hydro_prof(next_level,phase)%dcloudsc_da0 = dhsc(3)
!        hydro_prof(next_level,phase)%dcloudg_dt = dg(1)
!        hydro_prof(next_level,phase)%dcloudg_dk0 = dg(2)
!        hydro_prof(next_level,phase)%dcloudg_da0 = dg(3)
      phase = 5 ! graupel
      ip  = 0.0d0
      iq  = 0.0d0
      ik0 = 0.0d0
      ia0 = 0.0d0
      if( (atm_inp%prof(last_level)%cloud_grpl_k0 .gt. 0.0d0) .and. &
          (atm_inp%prof(next_level)%cloud_grpl_k0 .gt. 0.0d0)) THEN
        ip  = llw * atm_inp%prof(last_level)%cloud_grpl_p &
            + nlw * atm_inp%prof(next_level)%cloud_grpl_p
        iq  = llw * atm_inp%prof(last_level)%cloud_grpl_q &
            + nlw * atm_inp%prof(next_level)%cloud_grpl_q
        ik0 = llw * atm_inp%prof(last_level)%cloud_grpl_k0 &
            + nlw * atm_inp%prof(next_level)%cloud_grpl_k0
        ia0 = llw * atm_inp%prof(last_level)%cloud_grpl_a0 &
            + nlw * atm_inp%prof(next_level)%cloud_grpl_a0
      else
        if( atm_inp%prof(last_level)%cloud_grpl_k0 .gt. 0.0d0) then
          ip  = atm_inp%prof(last_level)%cloud_grpl_p
          iq  = atm_inp%prof(last_level)%cloud_grpl_q
          ik0 = llw * atm_inp%prof(last_level)%cloud_grpl_k0
          ia0 = atm_inp%prof(last_level)%cloud_grpl_a0
        end if
        if( atm_inp%prof(next_level)%cloud_grpl_k0 .gt. 0.0d0) then
          ip  = atm_inp%prof(next_level)%cloud_grpl_p
          iq  = atm_inp%prof(next_level)%cloud_grpl_q
          ik0 = nlw * atm_inp%prof(next_level)%cloud_grpl_k0
          ia0 = atm_inp%prof(next_level)%cloud_grpl_a0
        end if
      end if
      call hydrometeor_master_5ph_d( frequency, phase, temp, ip, iq, ik0, ia0, &
                                     abs_cloud_grpl, scat_cloud_grpl, g_cloud_grpl,   &
                                     dhab, dhsc, dg,                                  &
                                     a0_is_constant(next_level,phase) )
!        hydro_prof(next_level,phase)%cloudab = abs_cloud_grpl
!        hydro_prof(next_level,phase)%cloudsc = scat_cloud_grpl
!        hydro_prof(next_level,phase)%cloudg = g_cloud_grpl
!        hydro_prof(next_level,phase)%dcloudab_dt = dhab(1)
!        hydro_prof(next_level,phase)%dcloudab_dk0 = dhab(2)
!        hydro_prof(next_level,phase)%dcloudab_da0 = dhab(3)
!        hydro_prof(next_level,phase)%dcloudsc_dt = dhsc(1)
!        hydro_prof(next_level,phase)%dcloudsc_dk0 = dhsc(2)
!        hydro_prof(next_level,phase)%dcloudsc_da0 = dhsc(3)
!        hydro_prof(next_level,phase)%dcloudg_dt = dg(1)
!        hydro_prof(next_level,phase)%dcloudg_dk0 = dg(2)
!        hydro_prof(next_level,phase)%dcloudg_da0 = dg(3)
    ! interpolate to height hydrometeor distribution parameters for all hydrometeor phases - end
    ! total extinction
      abs_tot = cabs_o2 + cabs_h2o + abs_cloud_liq + abs_cloud_ice &
              + abs_cloud_rn + abs_cloud_snow + abs_cloud_grpl
      scat_tot = scat_cloud_liq + scat_cloud_ice + scat_cloud_rn &
               + scat_cloud_snow + scat_cloud_grpl
      ext_tot = abs_tot + scat_tot
      albedo = scat_tot / ext_tot
  ! calculate absorption, scattering, extinction, and albedo at height - end
  ! calculate opacity - start
    ! tau(1) = opacity from inp_height downward to surface at inp_theta
    ! tau(2) = opacity from inp_height upward              at inp_theta
    tau(1) = atm_inp%prof(1)%ext_tot * atm_inp%prof(1)%height ! opacity of first layer
    do j = 2, last_level
      tau(1) = tau(1) + atm_inp%prof(j)%ext_tot &
                      * ( atm_inp%prof(j)%height- atm_inp%prof(j-1)%height )
    end do
    tau(2) = tau(1) ! opacity from surface up to last level
    tau(1) = tau(1) + ext_tot * ( inp_height - atm_inp%prof(last_level)%height )
    do j = next_level, atm_inp%num_levels
      tau(2) = tau(2) + atm_inp%prof(j)%ext_tot &
                      * ( atm_inp%prof(j)%height- atm_inp%prof(j-1)%height )
    end do
    tau(2) = tau(2) - tau(1)
    cos_observ_angle = cos( inp_theta * pi / 180.0d0 )
    if( dabs(cos_observ_angle) > 0.0d0 ) then
      tau(1) = tau(1) / dabs(cos_observ_angle)
      tau(2) = tau(2) / dabs(cos_observ_angle)
    end if
  ! calculate opacity - end

  call cpu_time(times(2))
  time_k = time_k + ( times(2) - times(1) )
  call cpu_time(times(1))
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
    if( inp_theta <= 90.0d0 ) then      ! upwelling - looking down
      Tb_inp = Tbplinp
    else if( inp_theta > 90.0d0 ) then  ! downwelling - looking up
      Tb_inp = Tbmninp
    end if
    ! interpolate profiles to observation angle
    do j = 0, atm_inp%num_levels
      tb_pl_inp(j) = law * tb_pl(j,last_angle) + naw * tb_pl(j,next_angle)
      tb_mn_inp(j) = law * tb_mn(j,last_angle) + naw * tb_mn(j,next_angle)
      do k = 1, nvar
        dtb_pl_inp(j,k) = law * dtb_pl(j,last_angle,k) + naw * dtb_pl(j,next_angle,k)
        dtb_mn_inp(j,k) = law * dtb_mn(j,last_angle,k) + naw * dtb_mn(j,next_angle,k)
      end do
    end do

  call cpu_time(times(2))
  time_rt = time_rt + ( times(2) - times(1) )
return
end subroutine calc_mon_temp_weight_scat
