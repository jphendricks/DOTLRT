!============================================================================================
subroutine calc_tot_ext(freq)
!============================================================================================
! Calculates atmospheric partial absorption and scattering and total extinction for each ilev
!
! history:
!   9/26/2020 Kevin Schaefer deleted unused variables
!  10/20/2020 Kevin Schaefer deleted unused code
!--------------------------------------------------------------------------------------------
  use dotlrt_variables
  integer ilev
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

  do ilev = 1, nlev
    gas_prof(ilev)%absn2 = absn2(atm(ilev)%temp, atm(ilev)%press,freq, dabsn2_t, dabsn2_p)
    gas_prof(ilev)%absh2o = abh2o(atm(ilev)%temp, atm(ilev)%press, atm(ilev)%humid,freq, &
                                   dabh2o_t, dabh2o_v, dabh2o_p)
    gas_prof(ilev)%o2abs = o2abs(atm(ilev)%temp, atm(ilev)%press, atm(ilev)%humid,freq, &
                                  do2abs_t, do2abs_v,do2abs_p)
    atm(ilev)%abs_o2 = gas_prof(ilev)%o2abs + gas_prof(ilev)%absn2
    atm(ilev)%abs_h2o = gas_prof(ilev)%absh2o
 
    gas_prof(ilev)%dabsn2_dt = dabsn2_t
    gas_prof(ilev)%dabsn2_dp = dabsn2_p

    gas_prof(ilev)%dabsh2o_dt = dabh2o_t
    gas_prof(ilev)%dabsh2o_dp = dabh2o_p
    gas_prof(ilev)%dabsh2o_dw = dabh2o_v

    gas_prof(ilev)%do2abs_dt = do2abs_t
    gas_prof(ilev)%do2abs_dp = do2abs_p
    gas_prof(ilev)%do2abs_dw = do2abs_v

    phase = 1 ! cloud liquid
    call system_clock(time_start)
    call hydrometeor_master_5ph_d(freq, phase, &
                                   atm(ilev)%temp, &
                                   atm(ilev)%clw%p, &
                                   atm(ilev)%clw%q, &
                                   atm(ilev)%clw%k0, &
                                   atm(ilev)%clw%a0, &
                                   abs_cloud_liq, scat_cloud_liq, g_cloud_liq, &
                                   dhab, dhsc, dg, &
                                   atm(ilev)%clw%a0_const)

    name='hydro'
    call ex_time(2, name)
    call system_clock(time_start)
    hydro_prof(ilev,phase)%cloudab = abs_cloud_liq
    hydro_prof(ilev,phase)%cloudsc = scat_cloud_liq
    hydro_prof(ilev,phase)%cloudg = g_cloud_liq
    hydro_prof(ilev,phase)%dcloudab_dt = dhab(1)
    hydro_prof(ilev,phase)%dcloudab_dk0 = dhab(2)
    hydro_prof(ilev,phase)%dcloudab_da0 = dhab(3)
    hydro_prof(ilev,phase)%dcloudsc_dt = dhsc(1)
    hydro_prof(ilev,phase)%dcloudsc_dk0 = dhsc(2)
    hydro_prof(ilev,phase)%dcloudsc_da0 = dhsc(3)
    hydro_prof(ilev,phase)%dcloudg_dt = dg(1)
    hydro_prof(ilev,phase)%dcloudg_dk0 = dg(2)
    hydro_prof(ilev,phase)%dcloudg_da0 = dg(3)

    phase = 2 ! rain
    call system_clock(time_start)
    call hydrometeor_master_5ph_d(freq, phase, &
                                   atm(ilev)%temp, &
                                   atm(ilev)%rain%p, &
                                   atm(ilev)%rain%q, &
                                   atm(ilev)%rain%k0, &
                                   atm(ilev)%rain%a0, &
                                   abs_cloud_rn, scat_cloud_rn, g_cloud_rn, &
                                   dhab, dhsc, dg, &
                                   atm(ilev)%rain%a0_const)
    name='hydro'
    call ex_time(2, name)
    call system_clock(time_start)
    hydro_prof(ilev,phase)%cloudab = abs_cloud_rn
    hydro_prof(ilev,phase)%cloudsc = scat_cloud_rn
    hydro_prof(ilev,phase)%cloudg = g_cloud_rn
    hydro_prof(ilev,phase)%dcloudab_dt = dhab(1)
    hydro_prof(ilev,phase)%dcloudab_dk0 = dhab(2)
    hydro_prof(ilev,phase)%dcloudab_da0 = dhab(3)
    hydro_prof(ilev,phase)%dcloudsc_dt = dhsc(1)
    hydro_prof(ilev,phase)%dcloudsc_dk0 = dhsc(2)
    hydro_prof(ilev,phase)%dcloudsc_da0 = dhsc(3)
    hydro_prof(ilev,phase)%dcloudg_dt = dg(1)
    hydro_prof(ilev,phase)%dcloudg_dk0 = dg(2)
    hydro_prof(ilev,phase)%dcloudg_da0 = dg(3)

    phase = 3 ! ice
    call system_clock(time_start)
    call hydrometeor_master_5ph_d(freq, phase, &
                                   atm(ilev)%temp, &
                                   atm(ilev)%ice%p, &
                                   atm(ilev)%ice%q, &
                                   atm(ilev)%ice%k0, &
                                   atm(ilev)%ice%a0, &
                                   abs_cloud_ice, scat_cloud_ice, g_cloud_ice, &
                                   dhab, dhsc, dg, &
                                   atm(ilev)%ice%a0_const)
    name='hydro'
    call ex_time(2, name)
    call system_clock(time_start)
    hydro_prof(ilev,phase)%cloudab = abs_cloud_ice
    hydro_prof(ilev,phase)%cloudsc = scat_cloud_ice
    hydro_prof(ilev,phase)%cloudg = g_cloud_ice
    hydro_prof(ilev,phase)%dcloudab_dt = dhab(1)
    hydro_prof(ilev,phase)%dcloudab_dk0 = dhab(2)
    hydro_prof(ilev,phase)%dcloudab_da0 = dhab(3)
    hydro_prof(ilev,phase)%dcloudsc_dt = dhsc(1)
    hydro_prof(ilev,phase)%dcloudsc_dk0 = dhsc(2)
    hydro_prof(ilev,phase)%dcloudsc_da0 = dhsc(3)
    hydro_prof(ilev,phase)%dcloudg_dt = dg(1)
    hydro_prof(ilev,phase)%dcloudg_dk0 = dg(2)
    hydro_prof(ilev,phase)%dcloudg_da0 = dg(3)

    phase = 4 ! snow
    call system_clock(time_start)
    call hydrometeor_master_5ph_d(freq, phase, &
                                   atm(ilev)%temp, &
                                   atm(ilev)%snow%p, &
                                   atm(ilev)%snow%q, &
                                   atm(ilev)%snow%k0, &
                                   atm(ilev)%snow%a0, &
                                   abs_cloud_snow, scat_cloud_snow, g_cloud_snow, &
                                   dhab, dhsc, dg, &
                                   atm(ilev)%snow%a0_const)
    name='hydro'
    call ex_time(2, name)
    call system_clock(time_start)
    hydro_prof(ilev,phase)%cloudab = abs_cloud_snow
    hydro_prof(ilev,phase)%cloudsc = scat_cloud_snow
    hydro_prof(ilev,phase)%cloudg = g_cloud_snow
    hydro_prof(ilev,phase)%dcloudab_dt = dhab(1)
    hydro_prof(ilev,phase)%dcloudab_dk0 = dhab(2)
    hydro_prof(ilev,phase)%dcloudab_da0 = dhab(3)
    hydro_prof(ilev,phase)%dcloudsc_dt = dhsc(1)
    hydro_prof(ilev,phase)%dcloudsc_dk0 = dhsc(2)
    hydro_prof(ilev,phase)%dcloudsc_da0 = dhsc(3)
    hydro_prof(ilev,phase)%dcloudg_dt = dg(1)
    hydro_prof(ilev,phase)%dcloudg_dk0 = dg(2)
    hydro_prof(ilev,phase)%dcloudg_da0 = dg(3)

    phase = 5 ! graupel
    call system_clock(time_start)
    call hydrometeor_master_5ph_d(freq, phase, &
                                   atm(ilev)%temp, &
                                   atm(ilev)%grpl%p, &
                                   atm(ilev)%grpl%q, &
                                   atm(ilev)%grpl%k0, &
                                   atm(ilev)%grpl%a0, &
                                   abs_cloud_grpl, scat_cloud_grpl, g_cloud_grpl, &
                                   dhab, dhsc, dg, &
                                   atm(ilev)%grpl%a0_const)
    name='hydro'
    call ex_time(2, name)
    hydro_prof(ilev,phase)%cloudab = abs_cloud_grpl
    hydro_prof(ilev,phase)%cloudsc = scat_cloud_grpl
    hydro_prof(ilev,phase)%cloudg = g_cloud_grpl
    hydro_prof(ilev,phase)%dcloudab_dt = dhab(1)
    hydro_prof(ilev,phase)%dcloudab_dk0 = dhab(2)
    hydro_prof(ilev,phase)%dcloudab_da0 = dhab(3)
    hydro_prof(ilev,phase)%dcloudsc_dt = dhsc(1)
    hydro_prof(ilev,phase)%dcloudsc_dk0 = dhsc(2)
    hydro_prof(ilev,phase)%dcloudsc_da0 = dhsc(3)
    hydro_prof(ilev,phase)%dcloudg_dt = dg(1)
    hydro_prof(ilev,phase)%dcloudg_dk0 = dg(2)
    hydro_prof(ilev,phase)%dcloudg_da0 = dg(3)

    atm(ilev)%abs_cloud = abs_cloud_liq + abs_cloud_rn + abs_cloud_ice &
                                  + abs_cloud_snow + abs_cloud_grpl
    atm(ilev)%scat_cloud = scat_cloud_liq + scat_cloud_rn + scat_cloud_ice &
                                       + scat_cloud_snow + scat_cloud_grpl
    if((scat_cloud_liq + scat_cloud_rn + scat_cloud_ice + &
         scat_cloud_snow + scat_cloud_grpl) .ne. 0.0d0 ) then
         atm(ilev)%asymmetry = (scat_cloud_liq * g_cloud_liq &
                                       +  scat_cloud_rn * g_cloud_rn &
                                       +  scat_cloud_ice * g_cloud_ice &
                                       +  scat_cloud_snow * g_cloud_snow &
                                       +  scat_cloud_grpl * g_cloud_grpl) &
                                       / (scat_cloud_liq + scat_cloud_rn + scat_cloud_ice &
                                       +  scat_cloud_snow + scat_cloud_grpl) 
    else
        atm(ilev)%asymmetry = 0
    end if
    atm(ilev)%ext_tot = atm(ilev)%abs_o2 &
                                + atm(ilev)%abs_h2o &
                                + atm(ilev)%abs_cloud &
                                + atm(ilev)%scat_cloud
    atm(ilev)%albedo = atm(ilev)%scat_cloud &
                               / atm(ilev)%ext_tot
    atm(ilev)%bb_spec_int = atm(ilev)%temp
    
  end do
return
end subroutine calc_tot_ext
