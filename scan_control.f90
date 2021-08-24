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
  integer ichan ! (-) channel index
  integer iter  ! (-) iteration index
  integer iph   ! (-) phase index
  integer ivar  ! (-) variable index
  character(50) fmt    ! (-) format variable

  integer ichanidx, numchan
  real(8) :: starthydro, finishhydro, t_time, t_tot_chan
  real(8) :: z(12), tol

20      format ("[",4(f23.15,","),i2,"]")
21      format (5(f23.15,","))
22      format (A10, 5(f23.15,","))
!31      format (5(A23,","))
32      format (A10, 5(A23,","))

  tol = 0.1

! scasn MRT with generic cloud
  if (scantype=='cloud') then
     print*, "? BROKEN ?"
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
! hm 5 phases

  if (scantype=='hm_ph5') then

    ichanidx = 0

    atm(1)%clw%dens=gen_hm%dens
    atm(1)%rain%dens=gen_hm%dens
    atm(1)%ice%dens=gen_hm%dens
    atm(1)%snow%dens=gen_hm%dens
    atm(1)%grpl%dens=gen_hm%dens

    call calcprofile_d()

    do ichan=minchan, maxchan
        t_tot_chan = 0.0
        ichanidx = ichanidx+1
        call extract_channel(ichan)
        myfreq = channel%lo_freq
        ! do iph = 1, 5
        iph = 1
        phase=iph
        if(phase == 1) gen_hm=atm(1)%clw
        if(phase == 2) gen_hm=atm(1)%rain
        if(phase == 3) gen_hm=atm(1)%ice
        if(phase == 4) gen_hm=atm(1)%snow
        if(phase == 5) gen_hm=atm(1)%grpl

        phase=iph

        call cpu_time(starthydro)
        call hydrometeor_master_5ph_d(freq, &
                                      phase, &
                                      tair, &
                                      gen_hm%p, &
                                      gen_hm%q, &
                                      gen_hm%k0, &
                                      gen_hm%a0, &
                                      hab, &
                                      hsc, &
                                      g, &
                                      dhab, &
                                      dhsc, &
                                      dg, &
                                      gen_hm%a0_const)
        call cpu_time(finishhydro)
        t_time=finishhydro-starthydro
        if (dbg_index_table) write(*,'("TIME ORIG = ",f24.15)') t_time
        call cpu_time(starthydro)

        call get_hydro_d(value(xvar), &
                         value(yvar), &
                         phase, &
                         ichan, &
                         dens_lt, &
                         temp_lt, &
                         lookup_table, &
                         z)
        call cpu_time(finishhydro)
        t_time=finishhydro-starthydro
        if (dbg_index_table) write(*,'("TIME NEW  = ",f24.15)') t_time

        hab = z(1)
        hsc = z(2)
        g   = z(3)
        dhab(1) = z(4)
        dhsc(1) = z(5)
        dg  (1) = z(6)
        dhab(2) = z(6)
        dhsc(2) = z(8)
        dg  (2) = z(9)
        dhab(3) = z(10)
        dhsc(3) = z(11)
        dg  (3) = z(12)



        write(444,20)  value(xvar), value(yvar), myfreq, t_time, phase
        !if(dbg_index_table) then
           write(*,20) value(xvar), value(yvar), myfreq, t_time, phase
        !end if

        t_tot_chan = t_tot_chan + t_time

        myphase=iph
        if (gen_index_table) then
           parray(ix) = value(xvar)
           tarray(iy) = value(yvar)
           myarray(iph,ichanidx,ix,iy,1)  = hab
           myarray(iph,ichanidx,ix,iy,2)  = hsc
           myarray(iph,ichanidx,ix,iy,3)  = g
           myarray(iph,ichanidx,ix,iy,4)  = dhab(1)
           myarray(iph,ichanidx,ix,iy,5)  = dhsc(1)
           myarray(iph,ichanidx,ix,iy,6)  = dg(1)
           myarray(iph,ichanidx,ix,iy,7)  = dhab(2)
           myarray(iph,ichanidx,ix,iy,8)  = dhsc(2)
           myarray(iph,ichanidx,ix,iy,9)  = dg(2)
           myarray(iph,ichanidx,ix,iy,10) = dhab(3)
           myarray(iph,ichanidx,ix,iy,11) = dhsc(3)
           myarray(iph,ichanidx,ix,iy,12) = dg(3)
        endif

        if (dbg_index_table) then
           ! call get_hydro_d(value(xvar), &
           !                  value(yvar), &
           !                  phase, &
           !                  ichan, &
           !                  dens_lt, &
           !                  temp_lt, &
           !                  lookup_table, &
           !                  z)

           write(61,21) abs((z(1)-hab)/hab),        z(1), hab,     value(xvar), value(yvar)
           write(62,21) abs((z(2)-hsc)/hsc),        z(2), hsc,     value(xvar), value(yvar)
           write(63,21) abs((z(3)-g)/g),            z(3), g,       value(xvar), value(yvar)
           write(64,21) abs((z(4)-dhab(1))/dhab(1)),z(4), dhab(1), value(xvar), value(yvar)
           write(65,21) abs((z(5)-dhsc(1))/dhsc(1)),z(5), dhsc(1), value(xvar), value(yvar)
           write(66,21) abs((z(6)-dg(1))/dg(1)),    z(6), dg(1),   value(xvar), value(yvar)

           !print*, '-----------------------------'
           write(*,32) 'var' , 'percent_error', 'vals', 'real', 'dens', 'temp'
           write(*,22) 'hab' , abs((z(1)-hab)/hab),        z(1), hab,     value(xvar), value(yvar)
           write(*,22) 'hsc' , abs((z(2)-hsc)/hsc),        z(2), hsc,     value(xvar), value(yvar)
           write(*,22) 'g'   , abs((z(3)-g)/g),            z(3), g,       value(xvar), value(yvar)
           write(*,22) 'dhab', abs((z(4)-dhab(1))/dhab(1)),z(4), dhab(1), value(xvar), value(yvar)
           write(*,22) 'dhsc', abs((z(5)-dhsc(1))/dhsc(1)),z(5), dhsc(1), value(xvar), value(yvar)
           write(*,22) 'dg'  , abs((z(6)-dg(1))/dg(1)),    z(6), dg(1),   value(xvar), value(yvar)
        endif

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
      !enddo
      !print*, 't_tot_chan = ', t_tot_chan, value(xvar), value(yvar)
    enddo
  endif

    ! if (use_index_table) call nc_close_file(ncid)

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
    endif ! scantype=='mrt_func'

!   MRT by channel
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
      testvar1 = cloud%precip%top
      test10(ivar)=cost_tot
      do ichan=minchan, maxchan
        ivar=ivar+1
        test10(ivar)=cost_chan(ivar-1)
      enddo

!      ivar=1
!      do ichan=minchan, maxchan
!       test10(ivar)=Tbo_sim(ivar,1)
!       ivar=ivar+1
!      enddo
    endif ! scantype=='mrt_chan'

!   Optimize
    if (scantype=='mrt_opt') then
      call calc_init_precip( )
      update_man(7)%val1 = init_precip(1)
      atm = atm_init
      call calculate_atm_profile (atm, cur_lay, n_up_man, update_man)
      nits = 50

      fmt = '(a6,2x,a4,2x,5(a10,2x))'
      print(fmt), 'Case', 'iter', 'Precip', 'Cost', 'Residual', 'TB_sim', 'Tb_obs'
      fmt = '(a6,2x,i4,2x,5(f10.4,2x))'

      do iter = 1, nits
        precip_cur = update_man(7)%val1
        precip_del = 1.1d0*precip_cur
        if (precip_cur == 0.d0) precip_del = .1d0

!   base runs
        ivar=0
        do ichan=minchan, maxchan
          call extract_channel(ichan)
          call construct_single_surf_ref()
          call mrt( )
          ivar=ivar+1
          Tbo_sim(ivar,:)=Tbo_str_mat(1,:)
          res(ivar) = Tbo_obs(ivar,1)-Tbo_sim(ivar,1)
        enddo
        call cost_function ()
        if (iter==1) print(fmt), 'ok1 Opt',0, precip_cur, cost_tot, res(1), Tbo_sim(ivar,1), Tbo_obs(ivar,1)
        if(cost_tot<0.5d0 .or. dabs(res(1)) <.1 ) exit
        if(precip_cur < -12d0) exit

!   Calculate jacobian
        update_man(7)%val1=precip_del
        atm = atm_init
        call calculate_atm_profile (atm, cur_lay, n_up_man, update_man)
        ivar=0
        do ichan=minchan, maxchan
          call extract_channel(ichan)
          call construct_single_surf_ref()
          call mrt( )
          ivar=ivar+1
          Jacob(ivar)=(Tbo_str_mat(1,1)-Tbo_sim(ivar,1))/.1
        enddo

        call Optimizer_precip (atm, cur_lay)

        print(fmt), 'ok1 Opt',iter, precip_new, cost_tot, res(1), Tbo_sim(ivar,1), Tbo_obs(ivar,1)
        test10(1)=cost_tot
        test10(2)=testval(1)
        test10(3)=testval(2)
        test10(4)=testval(3)
        test10(5)=testval(4)
        test10(6)=testval(5)
        test10(7)=testval(6)
        test10(8)=testval(7)
        test10(9)=testval(8)
        atm = atm_init
        call calculate_atm_profile (atm, cur_lay, n_up_man, update_man)
      enddo
    endif ! scantype=='mrt_opt'

! Optimize
  if (scantype=='retrieval') then
    call basic_retrieval()
    call create_retrieval_output_file
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
