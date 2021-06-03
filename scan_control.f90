!
!=======================================================================
  subroutine scan_control()
!=========================================================================
! controls all calls to all scan subroutines
!
! Modifications:
!  5/27/2013 Kevin Schaefer split off from main program
!--------------------------------------------------------------------------
!
  use scan_Variables
  use dotlrt_variables
  use profiles
  use dotlrt_output
!
  IMPLICIT NONE
  
! local variables
  integer ichan ! channel index
  integer iph   ! phase index
  integer ivar  ! variable index index

! scasn MRT with generic cloud
  if (scantype=='cloud') then
  endif

! Dielectric calculations
  if (scantype=='dielectric') then
    do iph = 1, 5
      phase=iph
      if (phase == 1) call h2o_liquid_dielectric(freq, tair, epsil, depsil_dt)
      if (phase == 2) call h2o_liquid_dielectric(freq, tair, epsil, depsil_dt)
      if (phase == 3) call h2o_ice_dielectric(freq, dmin1(t_frz,tair), epsil, depsil_dt) 
      if (phase == 4) call h2o_mixed_dielectric(freq, tair, phase, epsil, depsil_dt) 
      if (phase == 5) call h2o_mixed_dielectric(freq, tair, phase, epsil, depsil_dt) 

      test10(iph)=dreal(depsil_dt)
    enddo
  endif

! hydrometeor scattering
  if (scantype=='hm_ph5') then
    atm(1)%clw%dens=gen_hm%dens
    atm(1)%rain%dens=gen_hm%dens
    atm(1)%ice%dens=gen_hm%dens
    atm(1)%snow%dens=gen_hm%dens
    atm(1)%grpl%dens=gen_hm%dens
    call calcprofile_d()
    do iph = 1, 5
      phase=iph
      if(phase == 1) gen_hm=atm(1)%clw
      if(phase == 2) gen_hm=atm(1)%rain
      if(phase == 3) gen_hm=atm(1)%ice
      if(phase == 4) gen_hm=atm(1)%snow
      if(phase == 5) gen_hm=atm(1)%grpl
      call hydrometeor_master_5ph_d(freq, phase, tair, gen_hm%p, gen_hm%q, gen_hm%k0, gen_hm%a0, &
                                   hab, hsc, g, dhab, dhsc, dg, gen_hm%a0_const)
      
      test10(iph)=dhab(1)
      !print*, iph, gen_hm%dens, gen_hm%a0
    enddo
  endif

! MRT functional response
  if (scantype=='mrt_func') then
    call system_clock(time_start)
    time_temp=time_start

    call construct_single_surf_ref()
    
    atm(:)%clw%dens  = 0.d0
    atm(:)%rain%dens = 0.d0
    atm(:)%ice%dens  = 0.d0
    atm(:)%snow%dens = 0.d0
    atm(:)%grpl%dens = 0.d0
    atm(7)%clw%dens  = gen_hm%dens
    atm(7)%rain%dens = gen_hm%dens
    atm(12)%ice%dens  = gen_hm%dens
    atm(12)%snow%dens = gen_hm%dens
    atm(12)%grpl%dens = gen_hm%dens

    call mrt( )
    !print*, Tbo_str_mat(1,:)

    name='mrt'
    time_start=time_temp
    call ex_time(1, name)
    test10(1)=testval(1)
    test10(2)=testval(2)
    test10(3)=testval(3)
    test10(4)=testval(4)
    test10(5)=testval(5)
    test10(6)=testval(6)
    test10(7)=testval(7)
    test10(8)=testval(8)
    !print*, time_seg(1), time_seg(2), time_seg(3)
    time_seg=0.d0
    num_call=0.d0
  endif

! MRT by channel
  if (scantype=='mrt_chan') then
    call calculate_atm_profile (atm, cur_lay, n_up_man, update_man)
    ivar=0
    do ichan=minchan, maxchan
      call extract_channel(ichan)
      call construct_single_surf_ref()

      call mrt( )

      ivar=ivar+1
      Tbo_sim(ivar,:)=Tbo_str_mat(1,:)
    enddo
    call cost_function ()
    ivar=1
    test=cost_tot
    test10(ivar)=cost_tot
    do ichan=minchan, maxchan
      ivar=ivar+1
      test10(ivar)=cost_chan(ivar-1)
    enddo

!    ivar=1
!    do ichan=minchan, maxchan
!     test10(ivar)=Tbo_sim(ivar,1)
!     ivar=ivar+1
!    enddo
  endif
!
! Planks function
  if (scantype=='plank') then
    call Plank_wavelength
    Btotal=Btotal+black
    print*, Btotal, wave
    call Plank_frequency
    test10(1)=black
  endif
!
! Planks function
  if (scantype=='ext_coef') call ABSORB_L93(Freq,tair,Psur,RH,kext)
!
! sunriuse and sunset
  if (scantype=='Sun') call SunriseSunset
!
! cosine zenith
  if (scantype=='cosz') then
    dec=23.5*sin(1.72e-2*(DOY-80))
    sind=sin(pi*dec/180.)
    cosd=cos(pi*dec/180.)
    call ZenithAngle (1, sind, cosd, lat, lon, tofday, cosz, test, test10)
  endif
!
! L-band soil moisture retrieval
   if(scantype=='L_band_VWC') then
    kroot(10)=5.5
    rho_om_max=140.
    poros_om=.9
    root_depth=0.7
    org_depth=0.18
    call dielectric_constant(depth, satfrac,sand,kroot(10),m_om,&
      root_depth, org_depth, rho_om_max, poros_om,test,test10)
   endif

  return
  end
!
