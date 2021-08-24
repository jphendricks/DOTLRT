!
!====================================================================
  subroutine basic_retrieval()
!====================================================================
! Does basic retrieval of FY data based on mrt_opt option
!
! History:
!  6/30/2021 Kevin Schaefer created routine
!--------------------------------------------------------------------
  use scan_Variables
  use dotlrt_variables
  use profiles
  use dotlrt_output
  implicit none
!
! internal variables
  integer ipt         ! (-) point index
  integer ivar        ! (-) variable index
  integer iter        ! (-) iteration index
  integer ichan       ! (-) channel index
  character(50) fmt   ! (-) format variable
  integer clock_rate  ! conversion between count and clock time
  integer clock_start ! time count at stat
  integer clock_stop  ! time count at end

! print message
  fmt = '(a8,2x,a4,2x,7(a10,2x))'
  print(fmt), 'point', 'iter', 'bad', 'ex', 'Precip', 'Cost', 'Residual', 'TB_sim', 'Tb_obs'
  fmt = '(i8,2x,i4,2x,7(f10.4,2x))'

! main swath point loop
  do ipt=1, 10 !swath_npts

! start measuring execution time per point
    call system_clock(clock_start) ! start counting

! get observed brightness temperatures
    swath_pt_num=ipt
    call extract_observed_brightness_temp (swath_pt_num)

! calculate first guess/initial precipitation
    call calc_init_precip( )
    update_man(7)%val1 = init_precip(1)

! set up initial atmospheric profile
    atm = atm_init
    call calculate_atm_profile (atm, cur_lay, n_up_man, update_man)

! Optimization loop to minimize cost function
    nits = 10
    do iter = 1, nits

! save current precip value
      precip_cur = update_man(7)%val1

! run forward model all channels with current value
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

! save optimization data
      swath_precip(ipt)=precip_cur
      swath_point(ipt)=ipt
      swath_iter(ipt)=iter
      swath_cost(ipt)=cost_tot
      swath_tb_sim(ipt)=Tbo_sim(1,1)
      swath_res(ipt)=res(1)
      call system_clock(clock_stop,clock_rate) ! stop counting
      swath_ex(ipt) = real(clock_stop-clock_start)/real(clock_rate)
      swath_bad(ipt) = 0.d0
      if(bad_point) swath_bad(ipt) = 1.d0

! print stuff to screen
      print(fmt), ipt,iter, swath_bad(ipt), swath_ex(ipt), precip_cur, cost_tot, res(1), Tbo_sim(1,1), Tbo_obs(1,1)

! check for bad data
      if(bad_point) exit

! check convergene
      if(cost_tot<0.1d0 .or. dabs(res(1)) <.1 ) exit
      if(precip_cur < -12d0) exit

! Calculate perturbed precip for finite difference Jacobian
      precip_del = 1.1d0*precip_cur
      if (precip_cur == 0.d0) precip_del = .1d0

! Run forward model for perturbed values to calculate finite difference jacobian
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

! Calculate new value
      call Optimizer_precip (atm, cur_lay)

! update the atmospheric profile
      atm = atm_init
      call calculate_atm_profile (atm, cur_lay, n_up_man, update_man)
    enddo
  enddo

  end subroutine basic_retrieval
!
!=======================================================================
  subroutine Optimizer_precip (atm_loc, layer_loc)
!=======================================================================
! estimates column total cloud and precip
!
! Modifications:
!  6/21/2021 Kevin Schaefer created subroutine
!----------------------------------------------------------------------
  use scan_Variables
  use dotlrt_variables
  use profiles
  use dotlrt_output

  implicit none

! input variables
  type(profile_type) atm_loc(max_nlev) ! (variable) local atmospheric profile
  type(hydro_layers) layer_loc     ! (varies) local hydrometeor layer variables

! local variables
  real(8) dden_dtot(nlev) ! (1/m) Jacobian hydro density wrt column total
  real(8) dden_dpk(nlev)  ! (g/m3/km) Jacobian hydro density wrt peak height
  real(8) dden_dstd(nlev) ! (g/m3/km) Jacobian hydro density wrt height standard deviation
  integer ilev  ! (-) level index
  integer ichan ! (-) channel index
  real(8) dtb_dw_loc(nlev,2) ! (K m3/g) local Jacobian brightness temperature wrt density
  real(8) jacobian        ! (varies) local Jacobian brightness temperature wrt hydrometeor layer parameters
  real(8) del             ! (varies) increment of hydrometeor layer parameters
  real(8) damp            ! (-) dampening coefficient

! Jacobian brightness temperature wrt cloud/precip density
! dTb_dw_str_mat(ilev,nstream,1,npol) = clw
! dTb_dw_str_mat(ilev,nstream,2,npol) = rain
! dTb_dw_str_mat(ilev,nstream,3,npol) = ice
! dTb_dw_str_mat(ilev,nstream,4,npol) = snow
! dTb_dw_str_mat(ilev,nstream,5,npol) = graupel
! dtb_dw_loc(ilev,1) = cloud
! dtb_dw_loc(ilev,2) = precip
  !print*, 'ilev, precip, rain, snow, graupel'
  do ilev =1, nlev
    dtb_dw_loc(ilev,1) = dTb_dw_str_mat(ilev,1,1,1) + dTb_dw_str_mat(ilev,1,3,1)
    dtb_dw_loc(ilev,2) = dTb_dw_str_mat(ilev,1,2,1) + dTb_dw_str_mat(ilev,1,4,1) + dTb_dw_str_mat(ilev,1,5,1)
    !print*, ilev, dtb_dw_loc(ilev,2), dTb_dw_str_mat(ilev,1,2,1), dTb_dw_str_mat(ilev,1,4,1), dTb_dw_str_mat(ilev,1,5,1)
  enddo

! Jacobian density wrt tot col precip
  call jacobian_per_layer_gaussian(atm_loc, layer_loc%precip, atm_loc(:)%cloud%ftot,dden_dtot, dden_dpk, dden_dstd)
  call dot_product(dtb_dw_loc(:,2), dden_dtot, nlev, jacobian)

! calculate increment
  damp=0.5d0
  del = 0.d0
  do ichan = 1, nchannel
    if (jacob(ichan) == 0.d0) then
      del = -damp*precip_cur*dabs(res(ichan))/res(ichan)
    else
      del = del + res(ichan)/jacob(ichan)
    endif
  enddo
  if(dabs(del) > 1.d0) then
    !print*, del
    del=1.d0*dabs(del)/del
  endif

  precip_new = precip_cur + del
  !print*, 'old: ', precip_cur, 'new: ',precip_new, 'Delta: ', del, 'jac: ', jacob(1)
  update_man(7)%val1 = precip_new
  !print*, man_loc(1)%val1, man_loc(1)%val2, man_loc(1)%val3
  !print*, man_loc(2)%val1, man_loc(2)%val2, man_loc(2)%val3

  testval(1)=precip_cur
  testval(2)=precip_new
  testval(3)=del
  testval(4)=jacobian
  testval(5)=jacob(1)
  testval(6)=res(1)
  testval(7)=Tbo_obs(1,1)
  testval(8)=Tbo_sim(1,1)

  return
  end subroutine Optimizer_precip
!
!=======================================================================
  subroutine Optimizer_col_tot (atm_loc, layer_loc, loc_n_man, man_loc)
!=======================================================================
! estimates column total cloud and precip
!
! Modifications:
!  6/21/2021 Kevin Schaefer created subroutine
!----------------------------------------------------------------------
  use scan_Variables
  use dotlrt_variables
  use profiles
  use dotlrt_output

  implicit none

! input variables
  type(profile_type) atm_loc(max_nlev) ! (variable) local atmospheric profile
  type(hydro_layers) layer_loc     ! (varies) local hydrometeor layer variables

! output variables
  integer loc_n_man                ! (varies) local number of manipulations
  type(manipulation) man_loc(1000) ! (varies) local manipulations for profile

! local variables
  real(8) dden_dtot(nlev) ! (1/m) Jacobian hydro density wrt column total
  real(8) dden_dpk(nlev)  ! (g/m3/km) Jacobian hydro density wrt peak height
  real(8) dden_dstd(nlev) ! (g/m3/km) Jacobian hydro density wrt height standard deviation
  integer ilev  ! (-) level index
  integer ichan ! (-) channel index
  integer ivar  ! (-) variable index
  real(8) param_cur(2)  ! (variable) vector of current parameters
  real(8) param_new(2)  ! (variable) vector of updated parameters
  real(8) dtb_dw_loc(nlev,2) ! (K m3/g) local Jacobian brightness temperature wrt density
  real(8) jacobian(2)        ! (varies) local Jacobian brightness temperature wrt hydrometeor layer parameters
  real(8) del(2)             ! (varies) lincrement of hydrometeor layer parameters

! transfer parameters to vector
  param_cur(1) = 10.d0**man_loc(1)%val1
  param_cur(2) = 10.d0**man_loc(2)%val1

! calculate residuals
  do ichan = 1, nchannel
    res(ichan) = Tbo_obs(ichan,1)-Tbo_sim(ichan,1)
!    print*, ichan, res(ichan), Tbo_obs(ichan,1), Tbo_sim(ichan,1)
  enddo

! Jacobian brightness temperature wrt cloud/precip density
! dTb_dw_str_mat(ilev,nstream,1,npol) = clw
! dTb_dw_str_mat(ilev,nstream,2,npol) = rain
! dTb_dw_str_mat(ilev,nstream,3,npol) = ice
! dTb_dw_str_mat(ilev,nstream,4,npol) = snow
! dTb_dw_str_mat(ilev,nstream,5,npol) = graupel
! dtb_dw_loc(ilev,1) = cloud
! dtb_dw_loc(ilev,2) = precip
  do ilev =1, nlev
    dtb_dw_loc(ilev,1) = dTb_dw_str_mat(ilev,1,1,1) + dTb_dw_str_mat(ilev,1,3,1)
    dtb_dw_loc(ilev,2) = dTb_dw_str_mat(ilev,1,2,1) + dTb_dw_str_mat(ilev,1,4,1) + dTb_dw_str_mat(ilev,1,5,1)
  enddo

! Jacobian density wrt hydrmeteor layer parameters
  call jacobian_per_layer_gaussian(atm_loc, layer_loc%cloud, atm_loc(:)%cloud%ftot,dden_dtot, dden_dpk, dden_dstd)
  call dot_product(dtb_dw_loc(:,1), dden_dtot, nlev, jacobian(1))

  call jacobian_per_layer_gaussian(atm_loc, layer_loc%precip, atm_loc(:)%precip%ftot,dden_dtot, dden_dpk, dden_dstd)
  call dot_product(dtb_dw_loc(:,2), dden_dtot, nlev, jacobian(2))

! calculate increment
  do ichan = 1, nchannel
    do ivar = 1, 2
      del(ivar) = res(ichan)/jacobian(ivar)
      if(dabs(del(ivar)) > 0.5*param_cur(ivar)) then
        del(ivar)=del(ivar)/dabs(del(ivar))*0.5*param_cur(ivar)
      endif
    enddo
  enddo
  param_new = param_cur + del
!  do ivar = 1, 2
!    print*, ivar, param_new(ivar),param_cur(ivar), del(ivar)
!  enddo
!print*, 'Total:', param_new(1), param_new(2)

! set number of update manipulations
  loc_n_man = 3

! load updated cloud parameters into update manipulation variable tree
  man_loc(1)%doit = .true.
  man_loc(1)%typ  = 'cloud'
  man_loc(1)%val1 = dlog10(param_new(1))
!  print*, man_loc(1)%val1, man_loc(1)%val2, man_loc(1)%val3

! load updated precipitation parameters into update manipulation variable tree
  man_loc(2)%doit = .true.
  man_loc(2)%typ  = 'precip'
  man_loc(2)%val1 = dlog10(param_new(2))
!  print*, man_loc(2)%val1, man_loc(2)%val2, man_loc(2)%val3

! print short profile
  man_loc(3)%doit = .false.
  man_loc(3)%typ  = 'print_prof'
  man_loc(3)%ind1 = 2

  return
  end subroutine Optimizer_col_tot
!
!=======================================================================
  subroutine Optimizer_full_vec (atm_loc, layer_loc, loc_n_man, man_loc)
!=======================================================================
! Constructs an atmospheric hydrometeor density profile from
! hydromet layer geometry
!
! Modifications:
!  6/8/2021 Kevin Schaefer created subroutine
!----------------------------------------------------------------------
  use scan_Variables
  use dotlrt_variables
  use profiles
  use dotlrt_output

  implicit none

! input variables
  type(profile_type) atm_loc(max_nlev) ! (variable) local atmospheric profile
  type(hydro_layers) layer_loc     ! (varies) local hydrometeor layer variables

! output variables
  integer loc_n_man                ! (varies) local number of manipulations
  type(manipulation) man_loc(1000) ! (varies) local manipulations for profile

! local variables
  real(8) dden_dtot(nlev) ! (1/m) Jacobian hydro density wrt column total
  real(8) dden_dpk(nlev)  ! (g/m3/km) Jacobian hydro density wrt peak height
  real(8) dden_dstd(nlev) ! (g/m3/km) Jacobian hydro density wrt height standard deviation
  integer ilev  ! (-) level index
  integer ichan ! (-) channel index
  integer ivar  ! (-) variable index
  real(8) param_cur(6)  ! (variable) vector of current parameters
  real(8) param_new(6)  ! (variable) vector of updated parameters
  real(8) dtb_dw_loc(nlev,2) ! (K m3/g) local Jacobian brightness temperature wrt density
  real(8) jacobian(6)        ! (varies) local Jacobian brightness temperature wrt hydrometeor layer parameters
  real(8) del(6)             ! (varies) lincrement of hydrometeor layer parameters

! transfer parameters to vector
  param_cur(1) = 10.d0**man_loc(1)%val1
  param_cur(2) = man_loc(1)%val2
  param_cur(3) = man_loc(1)%val3
  param_cur(4) = 10.d0**man_loc(2)%val1
  param_cur(5) = man_loc(2)%val2
  param_cur(6) = man_loc(2)%val3

! calculate residuals
  do ichan = 1, nchannel
    res(ichan) = Tbo_obs(ichan,1)-Tbo_sim(ichan,1)
!    print*, ichan, res(ichan), Tbo_obs(ichan,1), Tbo_sim(ichan,1)
  enddo

! Jacobian brightness temperature wrt cloud/precip density
! dTb_dw_str_mat(ilev,nstream,1,npol) = clw
! dTb_dw_str_mat(ilev,nstream,2,npol) = rain
! dTb_dw_str_mat(ilev,nstream,3,npol) = ice
! dTb_dw_str_mat(ilev,nstream,4,npol) = snow
! dTb_dw_str_mat(ilev,nstream,5,npol) = graupel
! dtb_dw_loc(ilev,1) = cloud
! dtb_dw_loc(ilev,2) = precip
  do ilev =1, nlev
    dtb_dw_loc(ilev,1) = dTb_dw_str_mat(ilev,1,1,1) + dTb_dw_str_mat(ilev,1,3,1)
    dtb_dw_loc(ilev,2) = dTb_dw_str_mat(ilev,1,2,1) + dTb_dw_str_mat(ilev,1,4,1) + dTb_dw_str_mat(ilev,1,5,1)
  enddo

! Jacobian density wrt hydrmeteor layer parameters
  call jacobian_per_layer_gaussian(atm_loc, layer_loc%cloud, atm_loc(:)%cloud%ftot,dden_dtot, dden_dpk, dden_dstd)
  call dot_product(dtb_dw_loc(:,1), dden_dtot, nlev, jacobian(1))
  call dot_product(dtb_dw_loc(:,1), dden_dpk, nlev,  jacobian(2))
  call dot_product(dtb_dw_loc(:,1), dden_dstd, nlev, jacobian(3))

  call jacobian_per_layer_gaussian(atm_loc, layer_loc%precip, atm_loc(:)%precip%ftot,dden_dtot, dden_dpk, dden_dstd)
  call dot_product(dtb_dw_loc(:,2), dden_dtot, nlev, jacobian(4))
  call dot_product(dtb_dw_loc(:,2), dden_dpk, nlev,  jacobian(5))
  call dot_product(dtb_dw_loc(:,2), dden_dstd, nlev, jacobian(6))

! calculate increment
  do ichan = 1, nchannel
    do ivar = 1, 6
      del(ivar) = res(ichan)/jacobian(ivar)
      if(dabs(del(ivar)) > 0.1*param_cur(ivar)) then
        del(ivar)=del(ivar)/dabs(del(ivar))*0.1*param_cur(ivar)
      endif
    enddo
  enddo
  param_new = param_cur + del
!  do ivar = 1, 6
!    print*, ivar, param_new(ivar),param_cur(ivar), del(ivar)
!  enddo
!print*, 'Total:', param_new(1), param_new(4)

! set number of update manipulations
  loc_n_man = 3

! load updated cloud parameters into update manipulation variable tree
  man_loc(1)%doit = .true.
  man_loc(1)%typ  = 'cloud'
  man_loc(1)%val1 = dlog10(param_new(1))
  man_loc(1)%val2 = param_new(2)
  man_loc(1)%val3 = param_new(3)

! load updated precipitation parameters into update manipulation variable tree
  man_loc(2)%doit = .true.
  man_loc(2)%typ  = 'precip'
  man_loc(2)%val1 = dlog10(param_new(4))
  man_loc(2)%val2 = param_new(5)
  man_loc(2)%val3 = param_new(6)

! print short profile
  man_loc(3)%doit = .false.
  man_loc(3)%typ  = 'print_prof'
  man_loc(3)%ind1 = 2

  return
  end subroutine Optimizer_full_vec
!
!=======================================================================
  subroutine calculate_atm_profile (atm_loc, layer_loc, loc_n_man, man_loc)
!=======================================================================
! Constructs an atmospheric hydrometeor density profile from
! hydromet layer geometry
!
! Modifications:
!  2/13/2021 Kevin Schaefer created subroutine
!  2/21/2021 Kevin Schaefer added all hydrometeor species
!  5/22/2021 Kevin Schaefer made generic for psuedo, init, and current profile
!----------------------------------------------------------------------
  use scan_Variables
  use dotlrt_variables
  use profiles
  use dotlrt_output

  implicit none

! input variables
  type(profile_type) atm_loc(max_nlev) ! (variable) local atmospheric profile
  type(hydro_layers) layer_loc ! local hydrometeor layer variables
  integer loc_n_man                  ! (-) local number of manipulations
  type(manipulation) man_loc(1000) ! local manipulations for profile

! local variables
  character*250 filename ! (-) local filename
  integer iman  ! (-) manipulation index
  integer ilev  ! (-) level index
  real(8) lapse ! (K/km) lapse rate

! loop through manipulations to construct atmospheric profile
  do iman=1, loc_n_man
    !print*, 'ok1 prof', iman, man_loc(iman)%doit, trim(man_loc(iman)%typ)

    if (man_loc(iman)%doit .and. man_loc(iman)%typ == 'clw') then
      atm_loc(:)%clw%dens = 0.d0
      call assign_hydromet_layer_variables (atm_loc, layer_loc%clw, man_loc(iman)%val1, man_loc(iman)%val2, man_loc(iman)%val3)
      call fraction_per_layer_gaussian(atm_loc, layer_loc%clw, atm_loc(:)%clw%ftot)
      call construct_hydromet_profile (atm_loc, layer_loc%clw%tot, atm_loc(:)%clw%ftot, atm_loc(:)%clw%dens)
      call calculate_hydromet_layer_variables (atm_loc, atm_loc(:)%clw%dens, layer_loc%clw)
    endif

    if (man_loc(iman)%doit .and. man_loc(iman)%typ == 'rain') then
      atm_loc(:)%rain%dens = 0.d0
      call assign_hydromet_layer_variables (atm_loc, layer_loc%rain, man_loc(iman)%val1, man_loc(iman)%val2, man_loc(iman)%val3)
      call fraction_per_layer_gaussian(atm_loc, layer_loc%rain, atm_loc(:)%rain%ftot)
      call construct_hydromet_profile (atm_loc, layer_loc%rain%tot, atm_loc(:)%rain%ftot, atm_loc(:)%rain%dens)
      call calculate_hydromet_layer_variables (atm_loc, atm_loc(:)%rain%dens, layer_loc%rain)
    endif

    if (man_loc(iman)%doit .and. man_loc(iman)%typ == 'ice') then
      atm_loc(:)%ice%dens = 0.d0
      call assign_hydromet_layer_variables (atm_loc, layer_loc%ice, man_loc(iman)%val1, man_loc(iman)%val2, man_loc(iman)%val3)
      call fraction_per_layer_gaussian(atm_loc, layer_loc%ice, atm_loc(:)%ice%ftot)
      call construct_hydromet_profile (atm_loc, layer_loc%ice%tot, atm_loc(:)%ice%ftot, atm_loc(:)%ice%dens)
      call calculate_hydromet_layer_variables (atm_loc, atm_loc(:)%ice%dens, layer_loc%ice)
    endif

    if (man_loc(iman)%doit .and. man_loc(iman)%typ == 'snow') then
      atm_loc(:)%snow%dens = 0.d0
      call assign_hydromet_layer_variables (atm_loc, layer_loc%snow, man_loc(iman)%val1, man_loc(iman)%val2, man_loc(iman)%val3)
      call fraction_per_layer_gaussian(atm_loc, layer_loc%snow, atm_loc(:)%snow%ftot)
      call construct_hydromet_profile (atm_loc, layer_loc%snow%tot, atm_loc(:)%snow%ftot, atm_loc(:)%snow%dens)
      call calculate_hydromet_layer_variables (atm_loc, atm_loc(:)%snow%dens, layer_loc%snow)
    endif

    if (man_loc(iman)%doit .and. man_loc(iman)%typ == 'grpl') then
      atm_loc(:)%grpl%dens = 0.d0
      call assign_hydromet_layer_variables (atm_loc, layer_loc%grpl, man_loc(iman)%val1, man_loc(iman)%val2, man_loc(iman)%val3)
      call fraction_per_layer_gaussian(atm_loc, layer_loc%grpl, atm_loc(:)%grpl%ftot)
      call construct_hydromet_profile (atm_loc, layer_loc%grpl%tot, atm_loc(:)%grpl%ftot, atm_loc(:)%grpl%dens)
      call calculate_hydromet_layer_variables (atm_loc, atm_loc(:)%grpl%dens, layer_loc%grpl)
    endif

    if (man_loc(iman)%doit .and. man_loc(iman)%typ == 'cloud') then
      atm_loc(:)%cloud%dens = 0.d0
      atm_loc(:)%clw%dens = 0.d0
      atm_loc(:)%ice%dens = 0.d0
      call assign_hydromet_layer_variables (atm_loc, layer_loc%cloud, man_loc(iman)%val1, man_loc(iman)%val2, man_loc(iman)%val3)
      call construct_cloud_density_profile (atm_loc, layer_loc)
      call calculate_hydromet_layer_variables (atm_loc, atm_loc(:)%clw%dens, layer_loc%clw)
      call calculate_hydromet_layer_variables (atm_loc, atm_loc(:)%ice%dens, layer_loc%ice)
      call calculate_hydromet_layer_variables (atm_loc, atm_loc(:)%hydro%dens, layer_loc%hydro)
    endif

    if (man_loc(iman)%doit .and. man_loc(iman)%typ == 'precip') then
      atm_loc(:)%precip%dens = 0.d0
      atm_loc(:)%rain%dens = 0.d0
      atm_loc(:)%snow%dens = 0.d0
      atm_loc(:)%grpl%dens = 0.d0
      call assign_hydromet_layer_variables (atm_loc, layer_loc%precip, man_loc(iman)%val1, man_loc(iman)%val2, man_loc(iman)%val3)
      call construct_precip_density_profile (atm_loc, layer_loc)
      call calculate_hydromet_layer_variables (atm_loc, atm_loc(:)%rain%dens, layer_loc%rain)
      call calculate_hydromet_layer_variables (atm_loc, atm_loc(:)%snow%dens, layer_loc%snow)
      call calculate_hydromet_layer_variables (atm_loc, atm_loc(:)%grpl%dens, layer_loc%grpl)
      call calculate_hydromet_layer_variables (atm_loc, atm_loc(:)%hydro%dens, layer_loc%hydro)
    endif

    if (man_loc(iman)%doit .and. man_loc(iman)%typ == 'temp') then
      lapse = (man_loc(iman)%val2 - man_loc(iman)%val1)/man_loc(iman)%val3
      do ilev = 1, nlev
        if (atm_loc(ilev)%hgt_mid < 15.d0) then
          atm_loc(ilev)%temp = man_loc(iman)%val1 + lapse*(atm_loc(ilev)%hgt_mid-atm_loc(1)%hgt_bot)
        else
          atm_loc(ilev)%temp = (man_loc(iman)%val1 + lapse*15.d0) - (atm_loc(ilev)%hgt_mid-15.d0)
        endif
      enddo
    endif

    if (man_loc(iman)%doit .and. man_loc(iman)%typ == 'read_sing') then
      filename = trim(files(man_loc(iman)%ind1)%path)
      call read_text_profile(atm_loc, filename)
      call calculate_hydromet_layer_variables (atm_loc, atm_loc(:)%clw%dens, layer_loc%clw)
      call calculate_hydromet_layer_variables (atm_loc, atm_loc(:)%rain%dens, layer_loc%rain)
      call calculate_hydromet_layer_variables (atm_loc, atm_loc(:)%ice%dens, layer_loc%ice)
      call calculate_hydromet_layer_variables (atm_loc, atm_loc(:)%snow%dens, layer_loc%snow)
      call calculate_hydromet_layer_variables (atm_loc, atm_loc(:)%grpl%dens, layer_loc%grpl)
      call calculate_hydromet_layer_variables (atm_loc, atm_loc(:)%cloud%dens, layer_loc%cloud)
      call calculate_hydromet_layer_variables (atm_loc, atm_loc(:)%precip%dens, layer_loc%precip)
      call calculate_hydromet_layer_variables (atm_loc, atm_loc(:)%hydro%dens, layer_loc%hydro)
    endif

    if (man_loc(iman)%doit .and. man_loc(iman)%typ == 'read_wrf') then
      call construct_atmospheric_profile(atm_loc, man_loc(iman)%ind1, man_loc(iman)%ind2)
      call calculate_hydromet_layer_variables (atm_loc, atm_loc(:)%clw%dens, layer_loc%clw)
      call calculate_hydromet_layer_variables (atm_loc, atm_loc(:)%rain%dens, layer_loc%rain)
      call calculate_hydromet_layer_variables (atm_loc, atm_loc(:)%ice%dens, layer_loc%ice)
      call calculate_hydromet_layer_variables (atm_loc, atm_loc(:)%snow%dens, layer_loc%snow)
      call calculate_hydromet_layer_variables (atm_loc, atm_loc(:)%grpl%dens, layer_loc%grpl)
      call calculate_hydromet_layer_variables (atm_loc, atm_loc(:)%cloud%dens, layer_loc%cloud)
      call calculate_hydromet_layer_variables (atm_loc, atm_loc(:)%precip%dens, layer_loc%precip)
      call calculate_hydromet_layer_variables (atm_loc, atm_loc(:)%hydro%dens, layer_loc%hydro)
    endif

    if (man_loc(iman)%doit .and. man_loc(iman)%typ == 'calc_tb') call calculate_brightness_temps(atm_loc)

    if (man_loc(iman)%doit .and. man_loc(iman)%typ == 'FY3_text') call read_observed_brightness_temp

    if (man_loc(iman)%doit .and. man_loc(iman)%typ == 'FY3_swath') call extract_observed_brightness_temp(man_loc(iman)%ind1)

    if (man_loc(iman)%doit .and. man_loc(iman)%typ == 'print_prof') then
      call print_reduced_profile(atm_loc, layer_loc, man_loc(iman)%ind1)
    endif

    if (man_loc(iman)%doit .and. man_loc(iman)%typ == 'write_prof') then
      filename = trim(files(man_loc(iman)%ind1)%path)
      call write_atmospheric_profile(atm_loc, filename)
    endif

    if (man_loc(iman)%doit .and. man_loc(iman)%typ == 'print_obs') call print_observed_brightness_temps

  enddo

  return
  end subroutine calculate_atm_profile
!
!====================================================================
  subroutine reconstruct_atm_layer_geometry(atm_loc)
!====================================================================
! recalculates atmospheric layer geometry in response to changes in
! cloud and precipitation layers
! done only when constructing an atmospheric profile from scratch
!
! History:
!  4/10/2021 Kevin Schaefer created routine
!  7/8/2021  Kevin Schaefer switched to local atm profile variable tree
!--------------------------------------------------------------------
  use scan_Variables
  use dotlrt_variables
  use profiles
  use dotlrt_output
  implicit none

! input variables
  type(profile_type) atm_loc(max_nlev) ! (variable) local atmospheric profile
!
! internal variables
  integer ilev         ! (-) level index
  integer ilev2        ! (-) level index
  real(8) delta_z_min  ! (km) minimum layer thickness
  real(8) delta_z_max  ! (km) maximum layer thickness
  real(8) scale_dz     ! (-) ratio between adjacent layer thicknesses
  real(8) hgt_max      ! (km) maximum height of atmospheric column
  real(8) del_z_scale  ! (km) scaler to normalize heights
  real(8) ratio        ! (-) scaling factor for interpolation
  character(50) fmt    ! (-) format variable
  logical split        ! (-) split layer flag
  logical print_flag   ! (-) print diagnostics flag

! set layer geometry parameters
  delta_z_min = 0.07d0
  delta_z_max = 0.7d0
  scale_dz = 1.2d0
  hgt_max= 25.d0

! clear out the atmospheric profile
! hydrometeors added later
  do ilev = 1, max_nlev
    atm_loc(ilev)%hgt_mid = 0.d0
    atm_loc(ilev)%hgt_top = 0.d0
    atm_loc(ilev)%hgt_bot = 0.d0
    atm_loc(ilev)%hgt_del = 0.d0
    atm_loc(ilev)%press   = 0.d0
    atm_loc(ilev)%temp    = 0.d0
    atm_loc(ilev)%humid   = 0.d0
  end do

! print cloud and precip geometry
  print_flag = .false.
  if (print_flag) then
    fmt = '(a8,2x,4(a9,2x))'
    print(fmt), 'Class', 'Bot (km)', 'Peak (km)', 'Top (km)', 'STD(km)'
    fmt = '(a8,2x,4(f9.3,2x))'
    print(fmt), 'precip', cur_lay%precip%bot, cur_lay%precip%pk, cur_lay%precip%top, cur_lay%precip%std
    print(fmt), 'cloud', cur_lay%cloud%bot, cur_lay%cloud%pk, cur_lay%cloud%top, cur_lay%cloud%std
    fmt = '(a3,2x,4(a8,2x))'
    print(fmt), 'Lay', 'Bot (km)', 'Mid (km)', 'Top (km)', 'Del(km)'
    fmt = '(i3,2x,4(f8.3,2x))'
  endif

! calculate new layer geometry
  nlev = 1
  atm_loc(1)%hgt_bot = 0.d0
  atm_loc(1)%hgt_del = delta_z_min
  atm_loc(1)%hgt_top = delta_z_min
  atm_loc(1)%hgt_mid = delta_z_min/2.d0
  do ilev = 2, max_nlev
    atm_loc(ilev)%hgt_del = scale_dz*atm_loc(ilev-1)%hgt_del
    atm_loc(ilev)%hgt_del = min(atm_loc(ilev)%hgt_del, delta_z_max)
    atm_loc(ilev)%hgt_bot = atm_loc(ilev-1)%hgt_top
    atm_loc(ilev)%hgt_top = atm_loc(ilev-1)%hgt_top + atm_loc(ilev)%hgt_del
    atm_loc(ilev)%hgt_mid = (atm_loc(ilev)%hgt_top + atm_loc(ilev)%hgt_bot)/2.d0
    nlev = nlev + 1

    if (print_flag) print(fmt), ilev, atm_loc(ilev)%hgt_bot, atm_loc(ilev)%hgt_mid, atm_loc(ilev)%hgt_top, atm_loc(ilev)%hgt_del
    if (atm_loc(ilev)%hgt_top > hgt_max) exit
  end do

! split layers to account for bottom of precip layer
  if (print_flag) then
    fmt = '(a3,2x,4(a8,2x))'
    print(fmt), 'Lay', 'Bot (km)', 'Mid (km)', 'Top (km)', 'Del(km)'
    fmt = '(i3,2x,4(f8.3,2x))'
  endif
  split = .true.
  if (cur_lay%precip%bot == 0.d0) split = .false.
  do ilev = 1, max_nlev
    if (atm_loc(ilev)%hgt_top > cur_lay%precip%bot .and. split) then
      split = .false.

! shift array up one level
      nlev = nlev+1
      do ilev2 = nlev, ilev+1, -1
        atm_loc(ilev2)%hgt_bot = atm_loc(ilev2-1)%hgt_bot
        atm_loc(ilev2)%hgt_mid = atm_loc(ilev2-1)%hgt_mid
        atm_loc(ilev2)%hgt_top = atm_loc(ilev2-1)%hgt_top
        atm_loc(ilev2)%hgt_del = atm_loc(ilev2-1)%hgt_del
      enddo

! adjust current level, which we split at precip%bot
      ! atm_loc(ilev+1)%hgt_bot unchanged
      atm_loc(ilev)%hgt_top = cur_lay%precip%bot
      atm_loc(ilev)%hgt_mid = (atm_loc(ilev)%hgt_top + atm_loc(ilev)%hgt_bot)/2.d0
      atm_loc(ilev)%hgt_del = atm_loc(ilev)%hgt_top - atm_loc(ilev)%hgt_bot

! adjust next level, which we split at precip%bot
      atm_loc(ilev+1)%hgt_bot = cur_lay%precip%bot
      !atm_loc(ilev+1)%hgt_top = unchanged
      atm_loc(ilev+1)%hgt_mid = (atm_loc(ilev+1)%hgt_top + atm_loc(ilev+1)%hgt_bot)/2.d0
      atm_loc(ilev+1)%hgt_del = atm_loc(ilev+1)%hgt_top - atm_loc(ilev+1)%hgt_bot
    endif
    if (print_flag) print(fmt), ilev, atm_loc(ilev)%hgt_bot, atm_loc(ilev)%hgt_mid, atm_loc(ilev)%hgt_top, atm_loc(ilev)%hgt_del
    if (atm_loc(ilev)%hgt_top > hgt_max) exit
  end do

! split layers to account for top of precip layer
  if (print_flag) then
    fmt = '(a3,2x,4(a8,2x))'
    print(fmt), 'Lay', 'Bot (km)', 'Mid (km)', 'Top (km)', 'Del(km)'
    fmt = '(i3,2x,4(f8.3,2x))'
  endif
  split = .true.
  do ilev = 1, max_nlev
    if (atm_loc(ilev)%hgt_top > cur_lay%precip%top .and. split) then
      split = .false.

! shift array up one level
      nlev = nlev+1
      do ilev2 = nlev, ilev+1, -1
        atm_loc(ilev2)%hgt_bot = atm_loc(ilev2-1)%hgt_bot
        atm_loc(ilev2)%hgt_mid = atm_loc(ilev2-1)%hgt_mid
        atm_loc(ilev2)%hgt_top = atm_loc(ilev2-1)%hgt_top
        atm_loc(ilev2)%hgt_del = atm_loc(ilev2-1)%hgt_del
      enddo

! adjust current level, which we split at precip%top
      ! atm_loc(ilev+1)%hgt_bot unchanged
      atm_loc(ilev)%hgt_top = cur_lay%precip%top
      atm_loc(ilev)%hgt_mid = (atm_loc(ilev)%hgt_top + atm_loc(ilev)%hgt_bot)/2.d0
      atm_loc(ilev)%hgt_del = atm_loc(ilev)%hgt_top - atm_loc(ilev)%hgt_bot

! adjust next level, which we split at precip%top
      atm_loc(ilev+1)%hgt_bot = cur_lay%precip%top
      !atm_loc(ilev+1)%hgt_top = unchanged
      atm_loc(ilev+1)%hgt_mid = (atm_loc(ilev+1)%hgt_top + atm_loc(ilev+1)%hgt_bot)/2.d0
      atm_loc(ilev+1)%hgt_del = atm_loc(ilev+1)%hgt_top - atm_loc(ilev+1)%hgt_bot
    endif
    if (print_flag) print(fmt), ilev, atm_loc(ilev)%hgt_bot, atm_loc(ilev)%hgt_mid, atm_loc(ilev)%hgt_top, atm_loc(ilev)%hgt_del
    if (atm_loc(ilev)%hgt_top > hgt_max) exit
  end do

! normalize the reference atmos profile to zero at bottom
  del_z_scale=atm_temp(1)%hgt_bot
  do ilev2 = 1, nlev_ref
    atm_temp(ilev2)%hgt_bot = atm_temp(ilev2)%hgt_bot - del_z_scale
    atm_temp(ilev2)%hgt_mid = atm_temp(ilev2)%hgt_mid - del_z_scale
    atm_temp(ilev2)%hgt_top = atm_temp(ilev2)%hgt_top - del_z_scale
    atm_temp(ilev2)%hgt_del = atm_temp(ilev2)%hgt_del - del_z_scale
  end do

! Interpolate temperature, pressure, and humidity to new center of each layer
! assume layer 1 matches
  atm_loc(1)%press = atm_temp(1)%press
  atm_loc(1)%temp = atm_temp(1)%temp
  atm_loc(1)%humid = atm_temp(1)%humid
  do ilev = 2, nlev
    do ilev2 = 1, nlev_ref
      if ( atm_temp(ilev2)%hgt_mid > atm_loc(ilev)%hgt_mid ) then
        !print*, ilev, ilev2, atm_loc(ilev)%hgt_mid, atm_temp(ilev2)%hgt_mid

        if(ilev2 == 1) then
          atm_loc(1)%press = atm_temp(1)%press
          atm_loc(1)%temp = atm_temp(1)%temp
          atm_loc(1)%humid = atm_temp(1)%humid
        else
          ratio = (atm_loc(ilev)%hgt_mid-atm_temp(ilev2-1)%hgt_mid) / (atm_temp(ilev2)%hgt_mid -atm_temp(ilev2-1)%hgt_mid )
          atm_loc(ilev)%press = ratio*atm_temp(ilev2)%press + (1.d0 - ratio)*scale*atm_temp(ilev2-1)%press
          atm_loc(ilev)%temp  = ratio*atm_temp(ilev2)%temp  + (1.d0 - ratio)*scale*atm_temp(ilev2-1)%temp
          atm_loc(ilev)%humid = ratio*atm_temp(ilev2)%humid + (1.d0 - ratio)*scale*atm_temp(ilev2-1)%humid
        endif

        exit
      endif
    enddo
  enddo

! Convert thickness to meters
  do ilev = 1, nlev
    atm_loc(ilev)%hgt_del = atm_loc(ilev)%hgt_del *1000.d0
  enddo

  end subroutine reconstruct_atm_layer_geometry
!
!====================================================================
  subroutine fraction_per_layer_gaussian(atm_loc,layer_loc, f_lay)
!====================================================================
! Calculates the hydrometeor fraction per layer using numerical integration
! assuming a Gaussian verticle distribution
!
! History:
!  5/8/2021 Kevin Schaefer created routine
!  7/8/2021 Kevin Schaefer switched to local atm profile variable tree
!--------------------------------------------------------------------
  use scan_Variables
  use dotlrt_variables
  use profiles
  use dotlrt_output
  implicit none

! input/output variable
  type(hydro_lay_geo) layer_loc ! input hydrometeor layer geometry
  type(profile_type) atm_loc(max_nlev) ! (variable) local atmospheric profile

! output variable
  real(8) f_lay(max_nlev) ! (-) hydrometeor fraction per layer of column total

! local variables
  integer ilev            ! (-) atm level index
  integer ngau            ! (-) number slices for Gaussian integral
  integer igau            ! (-) Gaussian integral index
  integer indx_top        ! (-) layer index top of hydrometeor layer
  integer indx_bot        ! (-) layer index bottom of hydrometeor layer
  real(8) del_h           ! (km) delta height for Gaussian integral
  real(8) loc_h           ! (km) current height for Gaussian integral
  real(8) f_gau(max_nlev) ! (-) unscaled gaussian fraction per layer
  real(8) hgt_var         ! (km) height variance
  real(8) total           ! (-) local total

! calculate height variance
  hgt_var = layer_loc%std*layer_loc%std

! find top and bottom of layer
  do ilev=1,nlev
    indx_bot=ilev
    if (atm_loc(ilev)%hgt_top > layer_loc%bot) exit
  enddo
  do ilev=nlev,1,-1
    indx_top=ilev
    if (atm_loc(ilev)%hgt_bot < layer_loc%top) exit
  enddo

! calculate raw, unscaled precip Gaussian fraction per layer
  f_gau(:) = 0.d0
  total=0.d0
  ngau=100
  do ilev=indx_bot,indx_top
    if (ilev == indx_bot) then ! adjust for partial of hydromet layer
      del_h= (atm_loc(ilev)%hgt_top - layer_loc%bot)/dble(ngau)
      loc_h = layer_loc%bot

    elseif (ilev == indx_top) then
      del_h= (layer_loc%top - atm_loc(ilev)%hgt_bot)/dble(ngau)
      loc_h = atm_loc(ilev)%hgt_bot

    else  ! adjust for partial of hydromet layer
      del_h= (atm_loc(ilev)%hgt_top - atm_loc(ilev)%hgt_bot)/dble(ngau)
      loc_h = atm_loc(ilev)%hgt_bot
    endif

    f_gau(ilev) = 0.d0
    do igau = 1, ngau
      f_gau(ilev) = f_gau(ilev)+ del_h*exp(-(loc_h-layer_loc%pk)**2.d0/hgt_var)
      loc_h = loc_h + del_h
    enddo
    if (f_gau(ilev) < 10.d-15) f_gau(ilev) = 0.d0

    total = total + f_gau(ilev)
    !print*, ilev, layer_loc%pk, layer_loc%std, f_gau(ilev)
  enddo

! Divide by the total integral to ensure all fractions sum to one
  do ilev=1,nlev
    f_lay(ilev) = f_gau(ilev)/total
    !if (ilev <20) print*, ilev, atm_loc(ilev)%hgt, f_gau(ilev), f_pre(ilev)
  enddo

  end subroutine fraction_per_layer_gaussian
!
!====================================================================
  subroutine fraction_frozen_hydro(atm_loc)
!====================================================================
! Calculates the liquid fraction of total hydrometeor based on freeze height
!
! History:
!  5/9/2021  Kevin Schaefer created routine
!  5/23/2021 Kevin Schaefer switched to local atmospheric profile
!  7/6/2021  Kevin Schaefer corrected to local atm profile for freeze height
!  7/8/2021  Kevin Schaefer corrected frz_height calculation
!--------------------------------------------------------------------
  use dotlrt_variables

  implicit none

! input variables
  type(profile_type) atm_loc(max_nlev) ! (variable) local atmospheric profile

! local variables
  integer ilev        ! (-) atm level index
  integer frz_lev     ! (-) freezing level index
  real(8) slope       ! (km/K) slope of height vs temp
  real(8) t_max       ! (K) maximum temp in profile
  real(8) t_min       ! (K) minimum temp in profile
  integer max_lev     ! (-) peak temp level index
  real(8) hgt_frz     ! (km) freezing height
  real(8) hgt_mix_top ! (km) top of mixed ice/liquid layer
  real(8) hgt_mix_bot ! (km) bottom of mixed ice/liquid layer

! calculate freeze height
  t_min = minval(atm_loc(:)%temp)
  t_max = maxval(atm_loc(:)%temp)
  if(t_max < t_frz) geom%h_frz = atm_loc(1)%hgt_bot
  if(t_min > t_frz) geom%h_frz = atm_loc(nlev)%hgt_top
  if(t_max >= t_frz .and. t_min <= t_frz) then
    do ilev = 1, nlev
      if (atm_loc(ilev)%temp == t_max) max_lev = ilev
    enddo
    do ilev = max_lev, nlev
      if (atm_loc(ilev)%temp < t_frz) exit
      frz_lev = ilev
    enddo
    slope = (atm_loc(frz_lev)%hgt_mid-atm_loc(frz_lev+1)%hgt_mid)/&
      (atm_loc(frz_lev)%temp-atm_loc(frz_lev+1)%temp)
    hgt_frz=atm_loc(frz_lev)%hgt_mid+slope*(t_frz-atm_loc(frz_lev)%temp)
  endif

! calculate mixed layer as +-1 km around freeze height
  hgt_mix_top = hgt_frz + 1.d0
  hgt_mix_bot = hgt_frz - 1.d0

! calculate frozen fractions of total hydrometeor per layer
! assume liquid transitions linearly from bottom to top of mixed layer
  do ilev=1,nlev
    if (atm_loc(ilev)%hgt_mid < hgt_mix_bot) then
      atm_loc(ilev)%f_froz = 0.d0
    elseif (atm_loc(ilev)%hgt_mid > hgt_mix_top) then
      atm_loc(ilev)%f_froz = 1.d0
    else
      atm_loc(ilev)%f_froz = 0.5d0*(atm_loc(ilev)%hgt_mid-hgt_frz)+0.5d0
    endif
  enddo

  end subroutine fraction_frozen_hydro

!=======================================================================
  subroutine construct_hydromet_profile (atm_loc, col_tot, frac_lay, prof_den)
!=======================================================================
!  Constructs an atmospheric profile of single hydrometeor
!
! Modifications:
!  2/28/2021 Kevin Schaefer created subroutine
!  4/1/2021  Kevin Schaefer added max density and vertical redistribution
!  5/23/2021 Kevin Schaefer removed max density and vertical redistribution
!  7/8/2021  Kevin Schaefer switched to local atm profile variable tree
!----------------------------------------------------------------------
  use scan_Variables
  use dotlrt_variables
  use profiles
  use dotlrt_output

  implicit none

! input variables
  type(profile_type) atm_loc(max_nlev) ! (variable) local atmospheric profile
  real(8) frac_lay(nlev) ! (-) fractions per layer of total colum hydromet mass
  real(8) col_tot        ! (g/m2) column total mass

! output variables
  real(8) prof_den(nlev) ! (g/m3) hydrometeor density profile

! local variables
  integer ilev      ! (-) atm level index

! Calculate densities
  do ilev=1,nlev
    prof_den(ilev) = frac_lay(ilev)*col_tot/atm_loc(ilev)%hgt_del
  enddo

  return
  end
!
!=======================================================================
  subroutine construct_precip_density_profile (atm_loc, layer_loc)
!=======================================================================
!  Calculates precipitation fraction of mass and density per layer
!
! Modifications:
!  3/16/2021 Kevin Schaefer created subroutine
!  5/8/2021  Kevin Schaefer moved Gaussian fraction to separate routine
!  5/9/2021  Kevin Schaefer moved frozen fraction fraction to separate routine
!  5/9/2021  Kevin Schaefer separated cloud and precip profiles into separate routines
!  5/23/2021 Kevin Schaefer corrected fraction calculations to atm_loc
!----------------------------------------------------------------------
  use scan_Variables
  use dotlrt_variables
  use profiles
  use dotlrt_output

  implicit none

! input/output variables
  type(profile_type) atm_loc(max_nlev) ! (variable) local atmospheric profile
  type(hydro_layers) layer_loc   ! (variable) local hydrometeor layer variables

! local variables
  integer ilev             ! (-) atm level index
  real(8) f_pre(max_nlev)  ! (-) precip fraction of total column precip per layer

! calculate fractions per layer
  call fraction_per_layer_gaussian(atm_loc, layer_loc%precip, f_pre)

! Calculate frozen fraction
  call fraction_frozen_hydro(atm_loc)

! calculate fraction of total column cloud water per species per layer
  do ilev=1,nlev
    atm_loc(ilev)%clw%fprecip  = 0.d0
    atm_loc(ilev)%rain%fprecip = f_pre(ilev)*(1.d0-atm_loc(ilev)%f_froz)
    atm_loc(ilev)%ice%fprecip  = 0.d0
    atm_loc(ilev)%snow%fprecip = f_pre(ilev)*0.5d0*atm_loc(ilev)%f_froz
    atm_loc(ilev)%grpl%fprecip = f_pre(ilev)*0.5d0*atm_loc(ilev)%f_froz
    atm_loc(ilev)%precip%ftot  = f_pre(ilev)
  enddo

! construct vertical density profiles
  call construct_hydromet_profile (atm_loc, layer_loc%precip%tot, atm_loc(:)%rain%fprecip, atm_loc(:)%rain%dens)
  call construct_hydromet_profile (atm_loc, layer_loc%precip%tot, atm_loc(:)%snow%fprecip, atm_loc(:)%snow%dens)
  call construct_hydromet_profile (atm_loc, layer_loc%precip%tot, atm_loc(:)%grpl%fprecip, atm_loc(:)%grpl%dens)
  call construct_hydromet_profile (atm_loc, layer_loc%precip%tot, atm_loc(:)%precip%ftot,  atm_loc(:)%precip%dens)

! calculate total hydrometeor density per layer
  do ilev=1,nlev
    atm_loc(ilev)%hydro%dens  = atm_loc(ilev)%precip%dens + atm_loc(ilev)%cloud%dens
  enddo

  return
  end subroutine construct_precip_density_profile

!=======================================================================
  subroutine construct_cloud_density_profile (atm_loc, layer_loc)
!=======================================================================
!  Calculates cloud fraction of mass and density per layer
!
! Modifications:
!  3/16/2021 Kevin Schaefer created subroutine
!  5/9/2021  Kevin Schaefer moved Gaussian & frozen frac to separate routine
!  5/9/2021  Kevin Schaefer separated cloud and precip profiles into separate routines
!  5/23/2021 Kevin Schaefer corrected fraction calculations to atm_loc
!----------------------------------------------------------------------
  use scan_Variables
  use dotlrt_variables
  use profiles
  use dotlrt_output

  implicit none

! input/output variables
  type(profile_type) atm_loc(max_nlev) ! (variable) local atmospheric profile
  type(hydro_layers) layer_loc   ! (variable) local hydrometeor layer variables

! local variables
  integer ilev             ! (-) atm level index
  real(8) f_cld(max_nlev)  ! (-) fraction of total column precip per layer

! calculate fractions per layer
  call fraction_per_layer_gaussian(atm_loc, layer_loc%cloud, f_cld)

! Calculate frozen fraction
  call fraction_frozen_hydro(atm_loc)

! calculate fraction of total column cloud water per species per layer
  do ilev=1,nlev
    atm_loc(ilev)%clw%fcloud  = f_cld(ilev)*(1.d0-atm_loc(ilev)%f_froz)
    atm_loc(ilev)%rain%fcloud = 0.d0
    atm_loc(ilev)%ice%fcloud  = f_cld(ilev)*atm_loc(ilev)%f_froz
    atm_loc(ilev)%snow%fcloud = 0.d0
    atm_loc(ilev)%grpl%fcloud = 0.d0
    atm_loc(ilev)%cloud%ftot  = f_cld(ilev)
  enddo

! construct vertical density profiles
  call construct_hydromet_profile (atm_loc, layer_loc%cloud%tot, atm_loc(:)%clw%fcloud, atm_loc(:)%clw%dens)
  call construct_hydromet_profile (atm_loc, layer_loc%cloud%tot, atm_loc(:)%ice%fcloud, atm_loc(:)%ice%dens)
  call construct_hydromet_profile (atm_loc, layer_loc%cloud%tot, atm_loc(:)%cloud%ftot, atm_loc(:)%cloud%dens)

! calculate total hydrometeor density per layer
  do ilev=1,nlev
    atm_loc(ilev)%hydro%dens  = atm_loc(ilev)%precip%dens + atm_loc(ilev)%cloud%dens
  enddo

  return
  end subroutine construct_cloud_density_profile
!
!=======================================================================
  subroutine calculate_hydromet_layer_variables (atm_loc, dens, hydro_geo_loc)
!=======================================================================
! assigns values to hydrometeor layer variables
!
! Modifications:
!  5/23/2021  Kevin Schaefer created subroutine from assign_hydromet_layer_variables
!----------------------------------------------------------------------
  use scan_Variables
  use dotlrt_variables
  use profiles
  use dotlrt_output

  implicit none

! input variables
  type(profile_type) atm_loc(max_nlev) ! (variable) local atmospheric profile
  real(8) dens(max_nlev) ! (g/m3) hydrometeor density profile
  type(hydro_lay_geo) hydro_geo_loc ! (-) local hydro layer geometry

! internal variables
  integer indx_top ! (-) level index
  integer indx_bot ! (-) level index
  integer ilev     ! (-) level index
  real(8) peak    ! (km) height of peak hydrometeor density

! calculate column total
  hydro_geo_loc%tot = 0.d0
  do ilev=1,nlev
   hydro_geo_loc%tot = hydro_geo_loc%tot + dens(ilev)*atm_loc(ilev)%hgt_del
  enddo

! check if the column has any hydrometeors at all
  if(hydro_geo_loc%tot == 0.d0) then
    hydro_geo_loc%bot = 0.d0
    hydro_geo_loc%top = 0.d0
    hydro_geo_loc%pk  = 0.d0
    hydro_geo_loc%std = 0.d0
    return
  endif

! find bottom of hydrometeor layer
  do ilev = 1, nlev
    indx_bot=ilev
    if(dens(ilev)/=0.d0) exit
  end do
  hydro_geo_loc%bot=atm_loc(indx_bot)%hgt_bot

! find top of hydrometeor layer
  do ilev = nlev,1,-1
    indx_top=ilev
    if(dens(ilev)/=0.d0) exit
  end do
  hydro_geo_loc%top=atm_loc(indx_top)%hgt_top

! calculate height of peak hydrometeor density
  peak=maxval(dens(indx_bot:indx_top))
  hydro_geo_loc%pk=atm_loc(indx_bot)%hgt_mid
  do ilev = indx_bot,indx_top
    if (dens(ilev) == peak) hydro_geo_loc%pk=atm_loc(ilev)%hgt_mid
  end do

! calculate std
  hydro_geo_loc%std = (hydro_geo_loc%top - hydro_geo_loc%bot)/4.d0

  return
  end subroutine calculate_hydromet_layer_variables
!
!=======================================================================
  subroutine assign_hydromet_layer_variables (atm_loc, hydro_geo_loc, col_tot, peak, std)
!=======================================================================
! assigns values to hydrometeor layer variables
!
! Modifications:
!  5/9/2021  Kevin Schaefer created subroutine
!  5/22/2021 Kevin Schaefer changed from psuedo dat only to generic
!  5/22/2021 Kevin Schaefer changed to a single hydrometeor
!  7/8/2021  Kevin Schaefer switched to local atm profile variable tree
!----------------------------------------------------------------------
  use scan_Variables
  use dotlrt_variables
  use profiles
  use dotlrt_output

  implicit none

! input variables
  type(profile_type) atm_loc(max_nlev) ! (variable) local atmospheric profile
  type(hydro_lay_geo) hydro_geo_loc ! (-) local hydro layer geometry
  real(8) col_tot ! (g/m2) column total hydrometeor
  real(8) peak    ! (km) height of peak hydrometeor density
  real(8) std     ! (km) height standard deviation of hydrometeor layer

  hydro_geo_loc%tot = 10.d0**col_tot
  hydro_geo_loc%pk  = peak
  hydro_geo_loc%std = std
  hydro_geo_loc%top = peak + 2.d0*std
  hydro_geo_loc%bot = peak - 2.d0*std
  if (hydro_geo_loc%bot < atm_loc(1)%hgt_bot) hydro_geo_loc%bot = atm_loc(1)%hgt_bot

  return
  end subroutine assign_hydromet_layer_variables
!
!=======================================================================
  subroutine calculate_brightness_temps (atm_loc)
!=======================================================================
! prints observed brightness temperatures
!
! Modifications:
!  5/23/2021 Kevin Schaefer created subroutine
!----------------------------------------------------------------------
  use scan_Variables
  use dotlrt_variables
  use profiles
  use dotlrt_output

  implicit none

! input variables
  type(profile_type) atm_loc(max_nlev) ! (variable) local atmospheric profile

! local variables
  integer ivar  ! (-) variable index
  integer ichan ! (-) channel index

  atm_temp = atm
  atm = atm_loc
  ivar=0
  do ichan=minchan, maxchan
    ivar=ivar+1
    call extract_channel(ichan)
    call construct_single_surf_ref()
    call mrt( )
    Tbo_obs(ivar,:) = Tbo_str_mat(1,:)
  enddo
  atm = atm_temp

  return
  end
!
!=======================================================================
  subroutine print_observed_brightness_temps ()
!=======================================================================
! prints observed brightness temperatures
!
! Modifications:
!  5/21/2021 Kevin Schaefer created subroutine
!----------------------------------------------------------------------
  use scan_Variables
  use dotlrt_variables
  use profiles
  use dotlrt_output

  implicit none

! local variables
  integer ivar  ! (-) variable index
  integer ichan ! (-) channel index
  character*100 fmt

  ivar=0
  fmt = '(a3,3x,a4,3x,3(a10,2x))'
  print(fmt), 'Num', 'chan', 'freq (GHz)', 'Tb H (K)', 'Tb V (K)'
  fmt = '(i3,3x,i4,3x,3(f10.3,2x))'
  do ichan=minchan, maxchan
    ivar = ivar+1
    call extract_channel(ichan)
    print(fmt), ivar, ichan, channel%lo_freq, Tbo_obs(ivar,:)
  enddo

  return
  end

!
!=======================================================================
  subroutine print_reduced_profile (atm_loc, layer_loc, print_flag)
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

! input variables
  type(profile_type) atm_loc(max_nlev) ! (variable) local atmospheric profile
  type(hydro_layers) layer_loc ! local hydrometeor layer variables
  integer print_flag ! (-) 1=full 2=short

! local variables
  integer ilev  ! (-) level index
  character*100 fmt

! recalculate reduced dimension stuff
  call construct_reduced_dim_atm_profile()

! print everything
  if (print_flag == 1) then
    fmt = '(a6,2x,4(a10,2x))'
    print(fmt),'phase', 'Col_tot', 'log10(tot)', 'h_peak', 'h_std'

    fmt = '(a6,2x,4(f10.3,2x))'
    print(fmt),'clw',  layer_loc%clw%tot,  log10(layer_loc%clw%tot), layer_loc%clw%pk, layer_loc%clw%std
    print(fmt),'rain', layer_loc%rain%tot, log10(layer_loc%rain%tot), layer_loc%rain%pk, layer_loc%rain%std
    print(fmt),'ice',  layer_loc%ice%tot,  log10(layer_loc%ice%tot), layer_loc%ice%pk, layer_loc%ice%std
    print(fmt),'snow', layer_loc%snow%tot, log10(layer_loc%snow%tot), layer_loc%snow%pk, layer_loc%snow%std
    print(fmt),'grpl', layer_loc%grpl%tot, log10(layer_loc%grpl%tot), layer_loc%grpl%pk, layer_loc%grpl%std
    print(fmt),'cloud', layer_loc%cloud%tot, log10(layer_loc%cloud%tot), layer_loc%cloud%pk, layer_loc%cloud%std
    print(fmt),'precip', layer_loc%precip%tot, log10(layer_loc%precip%tot), layer_loc%precip%pk, layer_loc%precip%std
    print(fmt),'hydro',layer_loc%hydro%tot,log10(layer_loc%hydro%tot), layer_loc%hydro%pk, layer_loc%hydro%std

    print*, 'Hydrometeor densities (g/m3)'
    fmt = '(a3,2x,10(a10,1x))'
    print(fmt), 'lev', 'ht (km)', 'temp', 'clw', 'rain', 'ice', 'snow', 'grpl', 'cloud', 'precip', 'hydro'
    fmt = '(i3,2x,10(f10.5,1x))'
    do ilev = 1, 25
      print(fmt), ilev, atm_loc(ilev)%hgt_mid,atm_loc(ilev)%temp, atm_loc(ilev)%clw%dens, atm_loc(ilev)%rain%dens, &
      atm_loc(ilev)%ice%dens,  atm_loc(ilev)%snow%dens,atm_loc(ilev)%grpl%dens, &
      atm_loc(ilev)%cloud%dens, atm_loc(ilev)%precip%dens, atm_loc(ilev)%hydro%dens
    enddo
  endif

! print summary
  if (print_flag == 2) then
    fmt = '(a6,2x,4(a10,2x))'
    print(fmt),'phase', 'Col_tot', 'log10(tot)', 'h_peak', 'h_std'

    fmt = '(a6,2x,4(f10.3,2x))'
    print(fmt),'cloud', layer_loc%cloud%tot, log10(layer_loc%cloud%tot), layer_loc%cloud%pk, layer_loc%cloud%std
    print(fmt),'precip', layer_loc%precip%tot, log10(layer_loc%precip%tot), layer_loc%precip%pk, layer_loc%precip%std
    print(fmt),'hydro',layer_loc%hydro%tot,log10(layer_loc%hydro%tot), layer_loc%hydro%pk, layer_loc%hydro%std
  endif

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
  real(8) num   ! (-) number channels, real version

! cost function
  ivar=0
  cost_tot=0.d0
  num=0.d0
  do ichan=minchan, maxchan
      ivar=ivar+1
      num=num+1.d0
      del1 = Tbo_obs(ivar,1)-Tbo_sim(ivar,1)
      del2 = Tbo_obs(ivar,2)-Tbo_sim(ivar,2)
      cost_chan(ivar) = (del1*del1 + del2*del2)/instr_spec(ichan)%dtrms/2.d0
      cost_tot = cost_tot + cost_chan(ivar)
  enddo

  cost_tot=cost_tot/num

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
      ZN=dcmplx(0.,0.)
      DO 10 I=1,44
       GAMMA=0.
       S=A(1,I)*P*V**3*EXP(A(2,I)*(1.-V))*1.E-6
       GAMMA=A(3,I)*(P*V**(0.8-A(4,I))+1.1*E*V)*1.E-3
       GAMMA=(GAMMA**2+(25*0.6E-4)**2)**0.5
       DELTA=(A(5,I)+A(6,I)*V)*(P+E)*(V**0.8)*1.E-3
       ZF=F/F0O2(I)*(dcmplx(1.,-DELTA)/dcmplx(F0O2(I)-F,-GAMMA)- &
                    dcmplx(1.,DELTA)/dcmplx(F0O2(I)+F,GAMMA))
       ZN=ZN+S*ZF
10    CONTINUE

! OXYGEN LINE ABSORPTION
! Cannot be less than 0.
      AT1=.182*F*AIMAG(ZN)
      IF(AT1.LT.0.) AT1=0.
!
! DRY AIR CONTINUUM
      ZN=dcmplx(0.,0.)
      So=6.14E-5*P*V**2
      GAMMAo=0.56E-3*(P+E)*V**.8
      ZFo=-F/dcmplx(F,GAMMAo)
      Sn=1.40E-12*p**2*V**3.5
      ZFn=dcmplx(0.,F/(1.93E-5*F**1.5+1.))
      ZN=So*ZFo+Sn*ZFn

! NONRESONAT DRY AIR ABSORPTION
      AT2=.182*F*AIMAG(ZN)
!
! WATER VAPOR
      ZN=dcmplx(0.,0.)
      DO 20 I=1,35
       GAMH=0.
       S=B(1,I)*E*V**3.5*EXP(B(2,I)*(1.-V))
! Doppler approximation.
       GAMH=B(3,I)*(P*V**B(5,I)+B(4,I)*E*V**B(6,I))*1.E-3
       GAMD2=1E-12/V*(1.46*F0H2O(I))**2
       GAMH=0.535*GAMH+(0.217*GAMH**2+GAMD2)**0.5
       DELH=0.
       ZF=F/F0H2O(I)*(dcmplx(1.,-DELH)/dcmplx(F0H2O(I)-F,-GAMH)- &
                   dcmplx(1.,DELH)/dcmplx(F0H2O(I)+F,GAMH))
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
        ZEp=Eps-f*((Eps-Epinf)/dcmplx(f,fD)+(Epinf-Eopt)/dcmplx(f,fS))
!
! ICE PERMITTIVITY [8]
      ELSE
        Ai=(62.*V-11.6)*1.E-4*EXP(-22.1*(V-1.))
        Bi=.542E-6*(-24.17+116.79/V+(V/(V-.9927))**2)
        Eps=3.15
! Complex Permittivity of Ice
        fice=f
        IF(f.LT..001)fice=.001
        ZEp=dcmplx(3.15,Ai/fice+Bi*fice)
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

!
!====================================================================
  subroutine fraction_per_layer_quad(atm_loc, layer_loc, f_lay)
!====================================================================
! Calculates the hydrometeor fraction per layer
! assuming a quadratic verticle distribution
!
! History:
!  5/26/2021 Kevin Schaefer created routine
!  7/8/2021  Kevin Schaefer switched to local atm profile variable tree
!--------------------------------------------------------------------
  use scan_Variables
  use dotlrt_variables
  use profiles
  use dotlrt_output
  implicit none

! input/output variables
  type(profile_type) atm_loc(max_nlev) ! (variable) local atmospheric profile
  type(hydro_lay_geo) layer_loc ! input hydrometeor layer geometry

! output variable
  real(8) f_lay(max_nlev) ! (-) hydrometeor fraction per layer of column total
  real(8) dden_dtot(max_nlev) ! (1/m) Jacobian density wrt total column density
  real(8) dden_dpk(max_nlev)  ! (g/m3/km) Jacobian density wrt peak height
  real(8) dden_dstd(max_nlev) ! (g/m3/km) Jacobian density wrt height standard deviation

! local variables
  integer ilev            ! (-) atm level index
  integer indx_top        ! (-) layer index top of hydrometeor layer
  integer indx_bot        ! (-) layer index bottom of hydrometeor layer
  real(8) z_top           ! (km) top of sublayer in physical coordinates
  real(8) z_bot           ! (km) bottom of sublayer in physical coordinates
  real(8) z_del(max_nlev) ! (m) thickness of sublayer in physical coordinates
  real(8) z_mid(max_nlev) ! (km) middle of sublayer in physical coordinates
  real(8) q_top(max_nlev) ! (-) top of sublayer in quadratic coordinates
  real(8) q_bot(max_nlev) ! (-) bottom of sublayer in quadratic coordinates
  real(8) fquad(max_nlev) ! (-) unscaled quadratic fraction per layer
  real(8) total           ! (-) local total
  real(8) dden_dq         ! (g/m3) Jacobian density wrt quadratic coordinate

! find top and bottom indicies of hydrometeor layer
  do ilev=1,nlev
    indx_bot=ilev
    if (atm_loc(ilev)%hgt_top > layer_loc%bot) exit
  enddo
  do ilev=nlev,1,-1
    indx_top=ilev
    if (atm_loc(ilev)%hgt_bot < layer_loc%top) exit
  enddo

! calculate unscaled fraction per layer
  fquad(:) = 0.d0
  z_del(:) = 0.d0
  total=0.d0
  do ilev=indx_bot,indx_top

! top and bottom of sublayer in physical coordinates
    if (ilev == indx_bot) then ! adjust for partial hydromet layer
      z_top = atm_loc(ilev)%hgt_top
      z_bot = layer_loc%bot

    elseif (ilev == indx_top) then ! adjust for partial hydromet layer
      z_top = layer_loc%top
      z_bot = atm_loc(ilev)%hgt_bot

    else
      z_top = atm_loc(ilev)%hgt_top
      z_bot = atm_loc(ilev)%hgt_bot
    endif
    z_del(ilev) = (atm_loc(ilev)%hgt_top - atm_loc(ilev)%hgt_bot)*1000
    z_mid(ilev) = (z_top - z_bot)/2.d0

! convert from physical to quadratic coordinates
    q_top(ilev)= (z_top-layer_loc%pk)/layer_loc%std
    q_bot(ilev)= (z_bot-layer_loc%pk)/layer_loc%std

! calculate fraction
    fquad(ilev) = 3.d0/8.d0*(q_top(ilev) - q_bot(ilev)) - (q_top(ilev)**3 - q_bot(ilev)**3)/12.d0
    if (fquad(ilev) < 10.d-15) fquad(ilev) = 0.d0

! calculate total
    total = total + fquad(ilev)
    !print*, ilev, layer_loc%pk, layer_loc%std, fquad(ilev)
   enddo

! Scale by dividing by the total integral to ensure all fractions sum to one
  do ilev=1,nlev
    f_lay(ilev) = fquad(ilev)/total
    !if (ilev <20) print*, ilev, atm_loc(ilev)%hgt, fquad(ilev), f_pre(ilev)
  enddo

! Calculate Jacobians
  dden_dtot(:) = 0.d0
  dden_dpk(:)  = 0.d0
  dden_dstd(:) = 0.d0
  do ilev=1,nlev
    dden_dtot(ilev) = layer_loc%tot/z_del(ilev)
    dden_dq = layer_loc%tot/z_del(ilev)*3.d0/32.d0*(q_bot(ilev)*q_bot(ilev) - q_top(ilev)*q_top(ilev))
    dden_dpk(ilev)  = -1.d0/layer_loc%std*layer_loc%tot/z_del(ilev)*dden_dq
    dden_dstd(ilev) = -(z_mid(ilev)- layer_loc%pk)/layer_loc%std/layer_loc%std*dden_dq
  enddo

  end subroutine fraction_per_layer_quad
!
!====================================================================
  subroutine dot_product(vec1, vec2, npts, dot_prod)
!====================================================================
! Calculates the jacobian per layer wrt hydrometeor parameters
! assuming a Gaussian verticle distribution
!
! History:
!  6/9/21/2021 Kevin Schaefer created routine
!--------------------------------------------------------------------
  use scan_Variables
  use dotlrt_variables
  use profiles
  use dotlrt_output
  implicit none

! input/output variable
  integer npts
  real(8) vec1(npts) ! (-) first vector
  real(8) vec2(npts) ! (-) second vector

! output variable
  real(8) dot_prod ! (-) dot product

! local variables
  integer ipt        ! (-) atm level index

! find top and bottom of hydrometeor layer
  dot_prod = 0.d0
!  print*, 'dot_product'
!  print*, 'ipt, dot_prod, vec1, vec2'
  do ipt=1,npts
    dot_prod = dot_prod + vec1(ipt)*vec2(ipt)
!    print*, ipt, dot_prod, vec1(ipt), vec2(ipt)
  enddo

  end subroutine dot_product
!
!====================================================================
  subroutine jacobian_per_layer_gaussian(atm_loc, layer_loc, f_lay, dden_dtot, dden_dpk, dden_dstd)
!====================================================================
! Calculates the jacobian per layer wrt hydrometeor parameters
! assuming a Gaussian verticle distribution
!
! History:
!  6/6/21/2021 Kevin Schaefer created routine
!--------------------------------------------------------------------
  use scan_Variables
  use dotlrt_variables
  use profiles
  use dotlrt_output
  implicit none

! input/output variable
  type(profile_type) atm_loc(max_nlev) ! (variable) local atmospheric profile
  type(hydro_lay_geo) layer_loc ! input hydrometeor layer geometry
  real(8) f_lay(nlev) ! (-) hydrometeor fraction per layer of column total

! output variable
  real(8) dden_dtot(nlev) ! (1/m) Jacobian hydro density wrt column total
  real(8) dden_dpk(nlev)  ! (g/m3/km) Jacobian hydro density wrt peak height
  real(8) dden_dstd(nlev) ! (g/m3/km) Jacobian hydro density wrt height standard deviation

! local variables
  integer ilev        ! (-) atm level index
  integer indx_top    ! (-) layer index top of hydrometeor layer
  integer indx_bot    ! (-) layer index bottom of hydrometeor layer
  real(8) hgt_var     ! (km) height variance
  real(8) loc_alpha   ! (-) local normalized vertical coordinate
  real(8) dalpha_dpk  ! (1/km) Jacobian normalized vertical coordinate wrt peak height
  real(8) dalpha_dstd ! (1/km) Jacobian normalized vertical coordinate wrt height standard deviation
  real(8) dden_dalpha ! (g/m3) Jacobian hydro density wrt normalized vertical coordinate

! calculate height variance
  hgt_var = layer_loc%std*layer_loc%std

! find top and bottom of hydrometeor layer
  do ilev=1,nlev
    indx_bot=ilev
    if (atm_loc(ilev)%hgt_top > layer_loc%bot) exit
  enddo
  do ilev=nlev,1,-1
    indx_top=ilev
    if (atm_loc(ilev)%hgt_bot < layer_loc%top) exit
  enddo

! clear out local Jacobian arrays
  dden_dtot(:) = 0.d0
  dden_dpk(:)  = 0.d0
  dden_dstd(:) = 0.d0

! calculate Jacobians per layer
!print*, 'jacobian_per_layer_gaussian'
!print*, 'ilev, dden_dtot, f_lay, hgt_del'

  do ilev=indx_bot,indx_top

! jacobian hydrometeor density wrt total colum hydrometeor mass
    dden_dtot(ilev) = f_lay(ilev)/atm_loc(ilev)%hgt_del
!    print*, ilev, dden_dtot(ilev), f_lay(ilev), atm_loc(ilev)%hgt_del

! Jacobian normalized vertical coordinate wrt peak height and height standard deviation
    dalpha_dpk = -1.d0/layer_loc%std
    dalpha_dstd = -(atm_loc(ilev)%hgt_mid - layer_loc%pk)/hgt_var

! Jacobian hydrometeor density wrt normalized vertical coordinate
    loc_alpha = (atm_loc(ilev)%hgt_mid - layer_loc%pk)/layer_loc%std
    dden_dalpha = -2.d0*loc_alpha*exp(-loc_alpha*loc_alpha)

! Jacobian hydrometeor density wrt peak height and height standard deviation
    dden_dpk(ilev)  = layer_loc%tot/atm_loc(ilev)%hgt_del*dden_dalpha*dalpha_dpk
    dden_dstd(ilev) = layer_loc%tot/atm_loc(ilev)%hgt_del*dden_dalpha*dalpha_dstd
  enddo

  end subroutine jacobian_per_layer_gaussian
