!============================================================================================
subroutine calc_tot_ext(freq)
!============================================================================================
! Calculates atmospheric partial absorption and scattering and total extinction for each level
!
! history:
!   9/26/2020 Kevin Schaefer deleted unused variables
!  10/20/2020 Kevin Schaefer deleted unused code
!--------------------------------------------------------------------------------------------
  use dotlrt_variables
  integer level
  real(8) freq
  real(8) abs_cloud_liq
  real(8) scat_cloud_liq
  real(8) g_cloud_liq
  real(8) abs_cloud_rn
  real(8) scat_cloud_rn
  real(8) g_cloud_rn
  real(8) abs_cloud_ice
  real(8) scat_cloud_ice
  real(8) g_cloud_ice
  real(8) abs_cloud_snow
  real(8) scat_cloud_snow
  real(8) g_cloud_snow
  real(8) abs_cloud_grpl
  real(8) scat_cloud_grpl
  real(8) g_cloud_grpl
  integer phase
  real(8), DIMENSION(3) :: dhab, dhsc, dg
  real(8) dabsn2_t, dabsn2_p
  real(8) dabh2o_t, dabh2o_v, dabh2o_p
  real(8) do2abs_t, do2abs_v,do2abs_p
  real(8), external  :: o2abs, absn2, abh2o
  
  do level = 1, nlev
    gas_prof(level)%absn2 = absn2(atm(level)%temp, &
                                  atm(level)%press,freq, &
                                  dabsn2_t, dabsn2_p)
    gas_prof(level)%absh2o = abh2o(atm(level)%temp, &
                                   atm(level)%press, &
                                   atm(level)%humid,freq, &
                                   dabh2o_t, dabh2o_v, dabh2o_p)
    gas_prof(level)%o2abs = o2abs(atm(level)%temp, &
                                  atm(level)%press, &
                                  atm(level)%humid,freq, &
                                  do2abs_t, do2abs_v,do2abs_p)
    atm(level)%abs_o2 = gas_prof(level)%o2abs + gas_prof(level)%absn2
    atm(level)%abs_h2o = gas_prof(level)%absh2o
 
    gas_prof(level)%dabsn2_dt = dabsn2_t
    gas_prof(level)%dabsn2_dp = dabsn2_p

    gas_prof(level)%dabsh2o_dt = dabh2o_t
    gas_prof(level)%dabsh2o_dp = dabh2o_p
    gas_prof(level)%dabsh2o_dw = dabh2o_v

    gas_prof(level)%do2abs_dt = do2abs_t
    gas_prof(level)%do2abs_dp = do2abs_p
    gas_prof(level)%do2abs_dw = do2abs_v

    phase = 1 ! cloud liquid
    call hydrometeor_master_5ph_d( freq, phase, &
                                   atm(level)%temp, &
                                   atm(level)%clw_p, &
                                   atm(level)%clw_q, &
                                   atm(level)%clw_k0, &
                                   atm(level)%clw_a0, &
                                   abs_cloud_liq, scat_cloud_liq, g_cloud_liq, &
                                   dhab, dhsc, dg, &
                                   a0_const(level,phase) )
    hydro_prof(level,phase)%cloudab = abs_cloud_liq
    hydro_prof(level,phase)%cloudsc = scat_cloud_liq
    hydro_prof(level,phase)%cloudg = g_cloud_liq
    hydro_prof(level,phase)%dcloudab_dt = dhab(1)
    hydro_prof(level,phase)%dcloudab_dk0 = dhab(2)
    hydro_prof(level,phase)%dcloudab_da0 = dhab(3)
    hydro_prof(level,phase)%dcloudsc_dt = dhsc(1)
    hydro_prof(level,phase)%dcloudsc_dk0 = dhsc(2)
    hydro_prof(level,phase)%dcloudsc_da0 = dhsc(3)
    hydro_prof(level,phase)%dcloudg_dt = dg(1)
    hydro_prof(level,phase)%dcloudg_dk0 = dg(2)
    hydro_prof(level,phase)%dcloudg_da0 = dg(3)

    !write(debugout,*) "liquid: abs,sc,g=", abs_cloud_liq, scat_cloud_liq, g_cloud_liq
    !call mexPrintf(debugout//achar(10))

    phase = 2 ! rain
    call hydrometeor_master_5ph_d( freq, phase, &
                                   atm(level)%temp, &
                                   atm(level)%rain_p, &
                                   atm(level)%rain_q, &
                                   atm(level)%rain_k0, &
                                   atm(level)%rain_a0, &
                                   abs_cloud_rn, scat_cloud_rn, g_cloud_rn, &
                                   dhab, dhsc, dg, &
                                   a0_const(level,phase) )
    hydro_prof(level,phase)%cloudab = abs_cloud_rn
    hydro_prof(level,phase)%cloudsc = scat_cloud_rn
    hydro_prof(level,phase)%cloudg = g_cloud_rn
    hydro_prof(level,phase)%dcloudab_dt = dhab(1)
    hydro_prof(level,phase)%dcloudab_dk0 = dhab(2)
    hydro_prof(level,phase)%dcloudab_da0 = dhab(3)
    hydro_prof(level,phase)%dcloudsc_dt = dhsc(1)
    hydro_prof(level,phase)%dcloudsc_dk0 = dhsc(2)
    hydro_prof(level,phase)%dcloudsc_da0 = dhsc(3)
    hydro_prof(level,phase)%dcloudg_dt = dg(1)
    hydro_prof(level,phase)%dcloudg_dk0 = dg(2)
    hydro_prof(level,phase)%dcloudg_da0 = dg(3)

    !write(debugout,*) "liquid: abs,sc,g=", abs_cloud_rn, scat_cloud_rn, g_cloud_rn
    !call mexPrintf(debugout//achar(10))

    phase = 3 ! ice
    call hydrometeor_master_5ph_d( freq, phase, &
                                   atm(level)%temp, &
                                   atm(level)%ice_p, &
                                   atm(level)%ice_q, &
                                   atm(level)%ice_k0, &
                                   atm(level)%ice_a0, &
                                   abs_cloud_ice, scat_cloud_ice, g_cloud_ice, &
                                   dhab, dhsc, dg, &
                                   a0_const(level,phase) )
    hydro_prof(level,phase)%cloudab = abs_cloud_ice
    hydro_prof(level,phase)%cloudsc = scat_cloud_ice
    hydro_prof(level,phase)%cloudg = g_cloud_ice
    hydro_prof(level,phase)%dcloudab_dt = dhab(1)
    hydro_prof(level,phase)%dcloudab_dk0 = dhab(2)
    hydro_prof(level,phase)%dcloudab_da0 = dhab(3)
    hydro_prof(level,phase)%dcloudsc_dt = dhsc(1)
    hydro_prof(level,phase)%dcloudsc_dk0 = dhsc(2)
    hydro_prof(level,phase)%dcloudsc_da0 = dhsc(3)
    hydro_prof(level,phase)%dcloudg_dt = dg(1)
    hydro_prof(level,phase)%dcloudg_dk0 = dg(2)
    hydro_prof(level,phase)%dcloudg_da0 = dg(3)

    phase = 4 ! snow
    call hydrometeor_master_5ph_d( freq, phase, &
                                   atm(level)%temp, &
                                   atm(level)%snow_p, &
                                   atm(level)%snow_q, &
                                   atm(level)%snow_k0, &
                                   atm(level)%snow_a0, &
                                   abs_cloud_snow, scat_cloud_snow, g_cloud_snow, &
                                   dhab, dhsc, dg, &
                                   a0_const(level,phase) )
    hydro_prof(level,phase)%cloudab = abs_cloud_snow
    hydro_prof(level,phase)%cloudsc = scat_cloud_snow
    hydro_prof(level,phase)%cloudg = g_cloud_snow
    hydro_prof(level,phase)%dcloudab_dt = dhab(1)
    hydro_prof(level,phase)%dcloudab_dk0 = dhab(2)
    hydro_prof(level,phase)%dcloudab_da0 = dhab(3)
    hydro_prof(level,phase)%dcloudsc_dt = dhsc(1)
    hydro_prof(level,phase)%dcloudsc_dk0 = dhsc(2)
    hydro_prof(level,phase)%dcloudsc_da0 = dhsc(3)
    hydro_prof(level,phase)%dcloudg_dt = dg(1)
    hydro_prof(level,phase)%dcloudg_dk0 = dg(2)
    hydro_prof(level,phase)%dcloudg_da0 = dg(3)

    phase = 5 ! graupel
    call hydrometeor_master_5ph_d( freq, phase, &
                                   atm(level)%temp, &
                                   atm(level)%grpl_p, &
                                   atm(level)%grpl_q, &
                                   atm(level)%grpl_k0, &
                                   atm(level)%grpl_a0, &
                                   abs_cloud_grpl, scat_cloud_grpl, g_cloud_grpl, &
                                   dhab, dhsc, dg, &
                                   a0_const(level,phase) )
    hydro_prof(level,phase)%cloudab = abs_cloud_grpl
    hydro_prof(level,phase)%cloudsc = scat_cloud_grpl
    hydro_prof(level,phase)%cloudg = g_cloud_grpl
    hydro_prof(level,phase)%dcloudab_dt = dhab(1)
    hydro_prof(level,phase)%dcloudab_dk0 = dhab(2)
    hydro_prof(level,phase)%dcloudab_da0 = dhab(3)
    hydro_prof(level,phase)%dcloudsc_dt = dhsc(1)
    hydro_prof(level,phase)%dcloudsc_dk0 = dhsc(2)
    hydro_prof(level,phase)%dcloudsc_da0 = dhsc(3)
    hydro_prof(level,phase)%dcloudg_dt = dg(1)
    hydro_prof(level,phase)%dcloudg_dk0 = dg(2)
    hydro_prof(level,phase)%dcloudg_da0 = dg(3)

    atm(level)%abs_cloud = abs_cloud_liq + abs_cloud_rn + abs_cloud_ice &
                                  + abs_cloud_snow + abs_cloud_grpl
    atm(level)%scat_cloud = scat_cloud_liq + scat_cloud_rn + scat_cloud_ice &
                                       + scat_cloud_snow + scat_cloud_grpl
    if( (scat_cloud_liq + scat_cloud_rn + scat_cloud_ice + &
         scat_cloud_snow + scat_cloud_grpl) .ne. 0.0d0 ) then
         atm(level)%asymmetry = (scat_cloud_liq * g_cloud_liq &
                                       +  scat_cloud_rn * g_cloud_rn &
                                       +  scat_cloud_ice * g_cloud_ice &
                                       +  scat_cloud_snow * g_cloud_snow &
                                       +  scat_cloud_grpl * g_cloud_grpl) &
                                       / (scat_cloud_liq + scat_cloud_rn + scat_cloud_ice &
                                       +  scat_cloud_snow + scat_cloud_grpl) 
    else
        atm(level)%asymmetry = 0
    end if
    atm(level)%ext_tot = atm(level)%abs_o2 &
                                + atm(level)%abs_h2o &
                                + atm(level)%abs_cloud &
                                + atm(level)%scat_cloud
    atm(level)%albedo = atm(level)%scat_cloud &
                               / atm(level)%ext_tot
    atm(level)%bb_spec_int = atm(level)%temp
    
  end do
return
end subroutine calc_tot_ext
