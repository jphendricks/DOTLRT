! RT_Jacobian.f90
! 28 March 2004    Bob Weber
subroutine RT_Jacobian()
  use variables
  implicit none
  integer ilr1, k, j, i
  real(8) al_abs, z
  do i0 = 1, nangover2
    do ilr1 = 1, nlr1
      ! abs_total1(ilr1) = absorption due to all atmospheric gases and hydrometeors
      ! dTbdTr(ilr1,i0,1)     =                              dtb_pl(ilr1,i0,1) /h1(ilr1) ! temperature up
      ! dTbdTr(ilr1,i0,2)     =                              dtb_mn(ilr1,i0,1) /h1(ilr1) ! temperature down
      ! dTbdKa(ilr1,i0,1)     =           abs_total1(ilr1) * dtb_pl(ilr1,i0,2) /h1(ilr1) ! absorption up
      ! dTbdKa(ilr1,i0,2)     =           abs_total1(ilr1) * dtb_mn(ilr1,i0,2) /h1(ilr1) ! absorption down
      ! dTbdKsliq(ilr1,i0,1)  = hydro_prof(ilr1,1)%cloudsc * dtb_pl(ilr1,i0,3) /h1(ilr1) ! liquid scatter up
      ! dTbdKsliq(ilr1,i0,2)  = hydro_prof(ilr1,1)%cloudsc * dtb_mn(ilr1,i0,3) /h1(ilr1) ! liquid scatter down
      ! dTbdgliq(ilr1,i0,1)   = hydro_prof(ilr1,1)%cloudg  * dtb_pl(ilr1,i0,4) /h1(ilr1) ! liquid asymmetry up
      ! dTbdgliq(ilr1,i0,2)   = hydro_prof(ilr1,1)%cloudg  * dtb_mn(ilr1,i0,4) /h1(ilr1) ! liquid asymmetry down
      ! dTbdKsrn(ilr1,i0,1)   = hydro_prof(ilr1,2)%cloudsc * dtb_pl(ilr1,i0,5) /h1(ilr1) ! rain scatter up
      ! dTbdKsrn(ilr1,i0,2)   = hydro_prof(ilr1,2)%cloudsc * dtb_mn(ilr1,i0,5) /h1(ilr1) ! rain scatter down
      ! dTbdgrn(ilr1,i0,1)    = hydro_prof(ilr1,2)%cloudg  * dtb_pl(ilr1,i0,6) /h1(ilr1) ! rain asymmetry up
      ! dTbdgrn(ilr1,i0,2)    = hydro_prof(ilr1,2)%cloudg  * dtb_mn(ilr1,i0,6) /h1(ilr1) ! rain asymmetry down
      ! dTbdKsice(ilr1,i0,1)  = hydro_prof(ilr1,3)%cloudsc * dtb_pl(ilr1,i0,7) /h1(ilr1) ! ice scatter up
      ! dTbdKsice(ilr1,i0,2)  = hydro_prof(ilr1,3)%cloudsc * dtb_mn(ilr1,i0,7) /h1(ilr1) ! ice scatter down
      ! dTbdgice(ilr1,i0,1)   = hydro_prof(ilr1,3)%cloudg  * dtb_pl(ilr1,i0,8) /h1(ilr1) ! ice asymmetry up
      ! dTbdgice(ilr1,i0,2)   = hydro_prof(ilr1,3)%cloudg  * dtb_mn(ilr1,i0,8) /h1(ilr1) ! ice asymmetry down
      ! dTbdKssnow(ilr1,i0,1) = hydro_prof(ilr1,4)%cloudsc * dtb_pl(ilr1,i0,9) /h1(ilr1) ! snow scatter up
      ! dTbdKssnow(ilr1,i0,2) = hydro_prof(ilr1,4)%cloudsc * dtb_mn(ilr1,i0,9) /h1(ilr1) ! snow scatter down
      ! dTbdgsnow(ilr1,i0,1)  = hydro_prof(ilr1,4)%cloudg  * dtb_pl(ilr1,i0,10)/h1(ilr1) ! snow asymmetry up
      ! dTbdgsnow(ilr1,i0,2)  = hydro_prof(ilr1,4)%cloudg  * dtb_mn(ilr1,i0,10)/h1(ilr1) ! snow asymmetry down
      ! dTbdKsgrpl(ilr1,i0,1) = hydro_prof(ilr1,5)%cloudsc * dtb_pl(ilr1,i0,11)/h1(ilr1) ! graupel scatter up
      ! dTbdKsgrpl(ilr1,i0,2) = hydro_prof(ilr1,5)%cloudsc * dtb_mn(ilr1,i0,11)/h1(ilr1) ! graupel scatter down
      ! dTbdggrpl(ilr1,i0,1)  = hydro_prof(ilr1,5)%cloudg  * dtb_pl(ilr1,i0,12)/h1(ilr1) ! graupel asymmetry up
      ! dTbdggrpl(ilr1,i0,2)  = hydro_prof(ilr1,5)%cloudg  * dtb_mn(ilr1,i0,12)/h1(ilr1) ! graupel asymmetry down
      dTbdTr(ilr1,i0,1)     = dtb_pl(ilr1,i0,1) ! temperature up
      dTbdTr(ilr1,i0,2)     = dtb_mn(ilr1,i0,1) ! temperature down
      dTbdKa(ilr1,i0,1)     = dtb_pl(ilr1,i0,2) ! absorption up
      dTbdKa(ilr1,i0,2)     = dtb_mn(ilr1,i0,2) ! absorption down
      dTbdKsliq(ilr1,i0,1)  = dtb_pl(ilr1,i0,3) ! liquid scatter up
      dTbdKsliq(ilr1,i0,2)  = dtb_mn(ilr1,i0,3) ! liquid scatter down
      dTbdgliq(ilr1,i0,1)   = dtb_pl(ilr1,i0,4) ! liquid asymmetry up
      dTbdgliq(ilr1,i0,2)   = dtb_mn(ilr1,i0,4) ! liquid asymmetry down
      dTbdKsrn(ilr1,i0,1)   = dtb_pl(ilr1,i0,5) ! rain scatter up
      dTbdKsrn(ilr1,i0,2)   = dtb_mn(ilr1,i0,5) ! rain scatter down
      dTbdgrn(ilr1,i0,1)    = dtb_pl(ilr1,i0,6) ! rain asymmetry up
      dTbdgrn(ilr1,i0,2)    = dtb_mn(ilr1,i0,6) ! rain asymmetry down
      dTbdKsice(ilr1,i0,1)  = dtb_pl(ilr1,i0,7) ! ice scatter up
      dTbdKsice(ilr1,i0,2)  = dtb_mn(ilr1,i0,7) ! ice scatter down
      dTbdgice(ilr1,i0,1)   = dtb_pl(ilr1,i0,8) ! ice asymmetry up
      dTbdgice(ilr1,i0,2)   = dtb_mn(ilr1,i0,8) ! ice asymmetry down
      dTbdKssnow(ilr1,i0,1) = dtb_pl(ilr1,i0,9) ! snow scatter up
      dTbdKssnow(ilr1,i0,2) = dtb_mn(ilr1,i0,9) ! snow scatter down
      dTbdgsnow(ilr1,i0,1)  = dtb_pl(ilr1,i0,10) ! snow asymmetry up
      dTbdgsnow(ilr1,i0,2)  = dtb_mn(ilr1,i0,10) ! snow asymmetry down
      dTbdKsgrpl(ilr1,i0,1) = dtb_pl(ilr1,i0,11) ! graupel scatter up
      dTbdKsgrpl(ilr1,i0,2) = dtb_mn(ilr1,i0,11) ! graupel scatter down
      dTbdggrpl(ilr1,i0,1)  = dtb_pl(ilr1,i0,12) ! graupel asymmetry up
      dTbdggrpl(ilr1,i0,2)  = dtb_mn(ilr1,i0,12) ! graupel asymmetry down
    end do ! ilr1
  end do ! i0
end subroutine RT_Jacobian
