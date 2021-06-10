!
!=======================================================================
  subroutine scan_control(ix,iy)
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
  use netcdf_utilities_mod

integer, intent(in) :: ix
integer, intent(in) :: iy
!
!  IMPLICIT NONE

! local variables
  integer ichan ! channel index
  integer iph   ! phase index
  integer ivar  ! variable index index

  integer ichanidx, ncid, id
  real(8) :: starthydro, finishhydro, t_time, t_tot_chan
  real(8), allocatable :: vals(:,:,:)
  logical, parameter :: debug = .false.

  allocate(vals(maxi,maxj,nchannel))

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
      if (phase == 4) call h2o_mixed_dielectric(freq, tair, phase, epsil, depsil_dt, testvar1, testvar2, testvar3)
      if (phase == 5) call h2o_mixed_dielectric(freq, tair, phase, epsil, depsil_dt, testvar1, testvar2, testvar3)

      test10(iph)=dreal(depsil_dt)
    enddo
  endif

! hydrometeor scattering
! hm 5 phases
  if (scantype=='hm_ph5') then
    ncid = nc_open_file_readonly(file_index_table)

    ichanidx = 0

    atm(1)%clw%dens=gen_hm%dens
    atm(1)%rain%dens=gen_hm%dens
    atm(1)%ice%dens=gen_hm%dens
    atm(1)%snow%dens=gen_hm%dens
    atm(1)%grpl%dens=gen_hm%dens
    !print*, 'gen_hm%q   = ', gen_hm%q
    call calcprofile_d()
    do ichan=minchan, maxchan
      t_tot_chan = 0.0
      ichanidx = ichanidx+1
      call extract_channel(ichan)
      myfreq = channel%lo_freq
      do iph = 1, 5
        phase=iph
        if(phase == 1) gen_hm=atm(1)%clw
        if(phase == 2) gen_hm=atm(1)%rain
        if(phase == 3) gen_hm=atm(1)%ice
        if(phase == 4) gen_hm=atm(1)%snow
        if(phase == 5) gen_hm=atm(1)%grpl

        call cpu_time(starthydro)
        call hydrometeor_master_5ph_d(myfreq, phase, tair, gen_hm%p, gen_hm%q, gen_hm%k0, gen_hm%a0, &
                                     hab, hsc, g, dhab, dhsc, dg, gen_hm%a0_const, testvar1, testvar2, testvar3 )


        call cpu_time(finishhydro)

        t_time=finishhydro-starthydro

        write(444,*) '[', value(xvar), ',', value(yvar), ',', myfreq, ',', t_time, ',', phase, '],'
        print*,      '[', value(xvar), ',', value(yvar), ',', myfreq, ',', t_time, ',', phase, '],'
        t_tot_chan = t_tot_chan + t_time

        myphase=iph
        tarray(ix) = value(xvar)
        parray(iy) = value(yvar)
        myarray(iph,ichanidx,ix,iy,1) = hab
        myarray(iph,ichanidx,ix,iy,2) = hsc
        myarray(iph,ichanidx,ix,iy,3) = g
        myarray(iph,ichanidx,ix,iy,4) = dhab(1)
        myarray(iph,ichanidx,ix,iy,5) = dhsc(1)
        myarray(iph,ichanidx,ix,iy,6) = dg(1)
        myarray(iph,ichanidx,ix,iy,7) = dhab(2)
        myarray(iph,ichanidx,ix,iy,8) = dhsc(2)
        myarray(iph,ichanidx,ix,iy,9) = dg(2)
        myarray(iph,ichanidx,ix,iy,10) = dhab(3)
        myarray(iph,ichanidx,ix,iy,11) = dhsc(3)
        myarray(iph,ichanidx,ix,iy,12) = dg(3)

        !do id = 1,12
           !call get_hydro_d(ncid,myfreq,value(xvar),value(yvar),phase, ichan, id, vals)
        !enddo

        if(.false.) then
           print*, '-----------------------------'
           print*, '    hab = ', hab
           print*, '    hsc = ', hsc
           print*, '      g = ', g
           print*, 'dhab(1) = ', dhab(1)
           print*, 'dhsc(1) = ', dhsc(1)
           print*, '  dg(1) = ', dg(1)
           print*, 'dhab(2) = ', dhab(2)
           print*, 'dhsc(2) = ', dhsc(2)
           print*, '  dg(2) = ', dg(2)
           print*, 'dhab(3) = ', dhab(3)
           print*, 'dhsc(3) = ', dhsc(3)
           print*, '  dg(3) = ', dg(3)
           print*, '-----------------------------'
        endif

        myvals(1) = hab
        myvals(2) = hsc
        myvals(3) = g
        myvals(4) = dhab(1)
        myvals(5) = dhsc(1)
        myvals(6) = dg(1)
        myvals(7) = dhab(2)
        myvals(8) = dhsc(2)
        myvals(9) = dg(2)
        myvals(10) = dhab(3)
        myvals(11) = dhsc(3)
        myvals(12) = dg(3)

        !if (hab /= 0.) then
        !  test10(iph)=dlog10(hab)
        !else
        !  test10(iph)=-999d0
        !endif
        ! save as an array
!       ! absorbtion
        !test = hab
        ! scattering
        !test = hsc
        ! asymetry
        !test = g
      enddo
      !print*, 't_tot_chan = ', t_tot_chan, value(xvar), value(yvar)
    enddo
  call nc_close_file(ncid)
  endif

! MRT functional response
  if (scantype=='mrt_func') then
    call system_clock(time_start)
    time_temp=time_start

    call construct_single_surf_ref()

    call mrt( )
    !print*, 'Tbo_str_mat(1,:) = ', Tbo_str_mat(1,:)

    name='mrt'
    time_start=time_temp
    call ex_time(1, name)
    test10(1)=dTb_dT_obs_mat(5,1)
    test10(2)=dTb_dT_obs_mat(6,1)
    test10(3)=dTb_dT_obs_mat(7,1)
    test10(4)=dTb_dT_obs_mat(8,1)
    test10(5)=dTb_dT_obs_mat(9,1)
    test10(6)=dTb_dT_obs_mat(10,1)
    test10(7)=dTb_dT_obs_mat(11,1)
    !print*, time_seg(1), time_seg(2), time_seg(3)
    time_seg=0.d0
    num_call=0.d0
    !test10(2)=Tbo_str_mat(2,1)
    !test10(3)=Tbo_str_mat(3,1)
    !test10(4)=Tbo_str_mat(4,1)
    !test10(5)=Tbo_str_mat(5,1)
    !test10(6)=Tbo_str_mat(6,1)
    !test10(7)=Tbo_str_mat(7,1)
    !test10(8)=Tbo_str_mat(8,1)
  endif

! MRT by channel
  if (scantype=='mrt_chan') then
    ivar=0
    do ichan=minchan, maxchan
      call extract_channel(ichan)
      call construct_single_surf_ref()

      call mrt( )

      !print*, Tbo_str_mat(1,:)
      ivar=ivar+1
      Tbo_sim(ivar,:)=Tbo_str_mat(1,:)
      test10(ivar)=Tbo_str_mat(1,1)
    enddo
    call cost_function ()
    ivar=1
    test=cost_tot
    testvar1 = cloud%precip%top
    test10(ivar)=cost_tot
    do ichan=minchan, maxchan
      ivar=ivar+1
      test10(ivar)=cost_chan(ivar-1)
    enddo
  endif
!
! Planks function
  if (scantype=='plank') then
    call Plank_wavelength
    Btotal=Btotal+black
    !print*, '(Btotal, wave) = ', Btotal, wave
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
