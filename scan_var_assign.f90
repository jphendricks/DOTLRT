!=======================================================================
  subroutine AssignVariableValues
!=======================================================================      
! Assigns values to all scanable variables.
!
! Modifications:
!  7/25/2001 Kevin Schaefer created routine
!  2/13/2021 Kevin Schaefer removed all variable not used for MRT
!----------------------------------------------------------------------
  use scan_Variables
  use dotlrt_variables
  use profiles
  use dotlrt_output
!
  implicit none

! local variables
  real(8) col_tot  ! (g/m2) local column total mass
  logical do_prof  ! (-) logical to construct vertical hydrometeor profile
  real(8),dimension(1) :: ppl
  real(8),dimension(1) :: ttl
  real(8),dimension(1) :: esst
  real(8),dimension(1) :: dtesst

! assign scanning values
  generic=value(1)
  Tplank=value(2)
  wave=value(3)*1.e-6
  freq=value(4)
  channel%lo_freq=value(4)
  RH=value(5)
  Tair=value(6)
  rhoair=value(7)
  psur=value(8)
  tofday=value(9)
  DOY=value(10)
  dec=value(11)
  cosz=value(12)
  Lat=value(13)
  Lon=value(14)
  obs_den=10.d0**value(15)
  gen_hm%dens=10.d0**value(16)

  col_tot=10.d0**value(17)
  if (Xvar==17.or.Yvar==17) then
    print*, 'var 17', col_tot, value(17)
    cloud%clw%tot=col_tot
    call construct_hydromet_profile (col_tot, atm(:)%clw%ftot, atm(:)%clw%dens)
  endif

  col_tot=10.d0**value(18)
  if (Xvar==18.or.Yvar==18) then
    print*, 'var 18', col_tot
    cloud%rain%tot=col_tot
    call construct_hydromet_profile (col_tot, atm(:)%rain%ftot, atm(:)%rain%dens)
  endif

  col_tot=10.d0**value(19)
  if (Xvar==19.or.Yvar==19) then
    print*, 'var 19 ', col_tot
    cloud%ice%tot=col_tot
    call construct_hydromet_profile (col_tot, atm(:)%ice%ftot, atm(:)%ice%dens)
  endif

  col_tot=10.d0**value(20)
  if (Xvar==20.or.Yvar==20) then
    print*, 'var 20 ', col_tot
    cloud%snow%tot=col_tot
    call construct_hydromet_profile (col_tot, atm(:)%snow%ftot, atm(:)%snow%dens)
  endif

  col_tot=10.d0**value(21)
  if (Xvar==21.or.Yvar==21) then
    print*, 'var 21 ', col_tot
    cloud%grpl%tot=col_tot
    call construct_hydromet_profile (col_tot, atm(:)%grpl%ftot, atm(:)%grpl%dens)
  endif
  
    col_tot=10.d0**value(22)
  if (Xvar==22.or.Yvar==22) then
    print*, 'var 22 ', col_tot
    cloud%cloud%tot=col_tot
    call construct_hydromet_profile (col_tot, atm(:)%clw%fcloud,  atm(:)%clw%dens)
    call construct_hydromet_profile (col_tot, atm(:)%ice%fcloud,  atm(:)%ice%dens)
  endif

  col_tot=10.d0**value(23)
  if (Xvar==23.or.Yvar==23) then
    print*, 'var 23 ', col_tot
    cloud%precip%tot=col_tot
    call construct_hydromet_profile (col_tot, atm(:)%rain%fhydro, atm(:)%rain%dens)
    call construct_hydromet_profile (col_tot, atm(:)%snow%fhydro, atm(:)%snow%dens)
    call construct_hydromet_profile (col_tot, atm(:)%grpl%fhydro, atm(:)%grpl%dens)
  endif


  col_tot=10.d0**value(24)
  if (Xvar==24.or.Yvar==24) then
    print*, 'var 24 ', col_tot
    cloud%hydro%tot=col_tot
    call construct_hydromet_profile (col_tot, atm(:)%clw%fhydro,  atm(:)%clw%dens)
    call construct_hydromet_profile (col_tot, atm(:)%rain%fhydro, atm(:)%rain%dens)
    call construct_hydromet_profile (col_tot, atm(:)%ice%fhydro,  atm(:)%ice%dens)
    call construct_hydromet_profile (col_tot, atm(:)%snow%fhydro, atm(:)%snow%dens)
    call construct_hydromet_profile (col_tot, atm(:)%grpl%fhydro, atm(:)%grpl%dens)
  endif

  cloud%cloud%top=value(25)
  cloud%cloud%bot=value(26)
  cloud%cloud%top= max(cloud%cloud%top,cloud%cloud%bot)
  cloud%precip%top=value(27)
  cloud%precip%bot=value(28)
  do_prof = .false.
  if(Xvar==25) do_prof = .true.
  if(yvar==25) do_prof = .true.
  if(Xvar==26) do_prof = .true.
  if(yvar==26) do_prof = .true.
  if(Xvar==27) do_prof = .true.
  if(yvar==27) do_prof = .true.
  if(Xvar==28) do_prof = .true.
  if(yvar==28) do_prof = .true.
  if (do_prof) then
    print*, cloud%cloud%top, cloud%cloud%bot, cloud%precip%bot
    call construct_precip_profile(cloud%cloud%top, cloud%cloud%bot, cloud%precip%bot)
  endif

  unused=value(29)
  unused=value(30)
  unused=value(31)
  unused=value(32)
  unused=value(33)
  unused=value(34)
  unused=value(35)
  unused=value(36)
  unused=value(37)
  unused=value(38)
  unused=value(39)
  unused=value(40)

! humidity
  call dtess_eau(1,ppl,ttl,esst,dtesst)

! generic data
  len=1

! root fraction
  kroot(1)=3.9
  kroot(2)=3.9
  kroot(3)=2.0
  kroot(4)=5.5
  kroot(5)=5.5
  kroot(6)=2.0
  kroot(7)=5.5
  kroot(8)=2.0
  kroot(9)=2.0
  kroot(10)=5.5
  kroot(11)=2.0
  kroot(12)=5.5

  return                                                                    
  end
