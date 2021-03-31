!============================================================
subroutine RT_Jacobian()
!============================================================
! Calculates Jacobian
!
! History:
!  28 March 2004    Bob Weber createdroutine
!  9/26/20  Kevin Schaefer deleted unused variables
!------------------------------------------------------------

  use dotlrt_variables
  implicit none
  integer ilr1
  do i0 = 1, nangover2
    do ilr1 = 1, nlr1
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
