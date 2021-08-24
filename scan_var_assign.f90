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
  integer iman ! index of number of update manipulations
  real(8),dimension(1) :: ppl
  real(8),dimension(1) :: ttl
  real(8),dimension(1) :: esst
  real(8),dimension(1) :: dtesst

! initialize update manipulation count
  iman = 0
  call Set_up_update_man

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

  if (Xvar==17.or.Yvar==17) then ! clw tot
    update_man(1)%doit = .true.
    update_man(1)%val1=value(17)
  endif

  if (Xvar==18.or.Yvar==18) then ! rain tot
    update_man(2)%doit = .true.
    update_man(2)%val1=value(18)
  endif

  if (Xvar==19.or.Yvar==19) then ! ice tot
    update_man(3)%doit = .true.
    update_man(3)%val1=value(19)
  endif

  if (Xvar==20.or.Yvar==20) then ! snow tot
    update_man(4)%doit = .true.
    update_man(4)%val1=value(20)
  endif

  if (Xvar==21.or.Yvar==21) then ! grpl tot
    update_man(5)%doit = .true.
    update_man(5)%val1=value(21)
  endif
  
  if (Xvar==22.or.Yvar==22) then ! cloud tot
    update_man(6)%doit = .true.
    update_man(6)%val1=value(22)
  endif

  if (Xvar==23.or.Yvar==23) then ! precip tot
    update_man(7)%doit = .true.
    update_man(7)%val1=value(23)
  endif

  if (Xvar==24.or.Yvar==24) then ! hydro tot
    update_man(8)%doit = .true.
    update_man(8)%val1=value(24)
  endif

  if (Xvar==25.or.Yvar==25) then ! cloud peak
    update_man(6)%doit = .true.
    update_man(6)%val2=value(25)
  endif

  if (Xvar==26.or.Yvar==26) then ! cloud std
    update_man(6)%doit = .true.
    update_man(6)%val3=value(26)
  endif

  if (Xvar==27.or.Yvar==27) then ! precip peak
    update_man(7)%doit = .true.
    update_man(7)%val2=value(27)
  endif

  if (Xvar==28.or.Yvar==28) then ! precip std
    update_man(7)%doit = .true.
    update_man(7)%val3=value(28)
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

! some optional debugging code
!  iman = iman+1
!  n_up_man=iman
!  update_man(iman)%doit = .true.
!  update_man(iman)%typ = 'print_prof'
  
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
!=======================================================================
  subroutine Set_up_update_man
!=======================================================================      
! sets up update manipulation variable tree
!
! Modifications:
!  5/25/2021 Kevin Schaefer created routine
!----------------------------------------------------------------------
  use scan_Variables
  use dotlrt_variables
  use profiles
  use dotlrt_output
!
  implicit none

! local variables
  integer iman ! index of number of update manipulations

  n_up_man = 9

! assign values to update manipulation variable tree
  iman = 1
  update_man(iman)%doit = .false.
  update_man(iman)%typ = 'clw'
  update_man(iman)%val1=dlog10(cur_lay%clw%tot)
  update_man(iman)%val2=cur_lay%clw%pk
  update_man(iman)%val3=cur_lay%clw%std

  iman = 2
  update_man(iman)%doit = .false.
  update_man(iman)%typ = 'rain'
  update_man(iman)%val1=dlog10(cur_lay%rain%tot)
  update_man(iman)%val2=cur_lay%rain%pk
  update_man(iman)%val3=cur_lay%rain%std

  iman = 3
  update_man(iman)%doit = .false.
  update_man(iman)%typ = 'ice'
  update_man(iman)%val1=dlog10(cur_lay%ice%tot)
  update_man(iman)%val2=cur_lay%ice%pk
  update_man(iman)%val3=cur_lay%ice%std

  iman = 4
  update_man(iman)%doit = .false.
  update_man(iman)%typ = 'snow'
  update_man(iman)%val1=dlog10(cur_lay%snow%tot)
  update_man(iman)%val2=cur_lay%snow%pk
  update_man(iman)%val3=cur_lay%snow%std

  iman = 5
  update_man(iman)%doit = .false.
  update_man(iman)%typ = 'grpl'
  update_man(iman)%val1=dlog10(cur_lay%grpl%tot)
  update_man(iman)%val2=cur_lay%grpl%pk
  update_man(iman)%val3=cur_lay%grpl%std
  
  iman = 6
  update_man(iman)%doit = .false.
  update_man(iman)%typ = 'cloud'
  update_man(iman)%val1=typical(23)
  update_man(iman)%val2=typical(25)
  update_man(iman)%val3=typical(26)

  iman = 7
  update_man(iman)%doit= .false.
  update_man(iman)%typ = 'precip'
  update_man(iman)%val1=typical(24)
  update_man(iman)%val2=typical(27)
  update_man(iman)%val3=typical(28)

  iman = 8
  update_man(iman)%doit = .false.
  update_man(iman)%typ = 'hydro'
  update_man(iman)%val1=dlog10(cur_lay%hydro%tot)
  update_man(iman)%val2=cur_lay%hydro%pk
  update_man(iman)%val3=cur_lay%hydro%std

  iman = 9
  update_man(iman)%doit = .false.
  update_man(iman)%typ = 'print_prof'
  update_man(iman)%ind1=1

  return                                                                    
  end
