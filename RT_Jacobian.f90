!============================================================
subroutine RT_Jacobian()
!============================================================
! Calculates Jacobian wrt internal constants
!
! History:
!  3/28/2004  Bob Weber createdroutine
!  9/26/2020  Kevin Schaefer deleted unused variables
!  12/12/2020 Kevin Schaefer changed to standard indexing variables
!------------------------------------------------------------
use dotlrt_variables
implicit none

! internal variables
  integer ilev ! (-) vertical level index
  integer iang ! (-) stream angle index

  do iang = 1, nangover2
    do ilev = 1, nlev
      dTbdTr(ilev,iang,1)     = dtb_pl(ilev,iang,1) ! temperature up
      dTbdTr(ilev,iang,2)     = dtb_mn(ilev,iang,1) ! temperature down
      dTbdKa(ilev,iang,1)     = dtb_pl(ilev,iang,2) ! absorption up
      dTbdKa(ilev,iang,2)     = dtb_mn(ilev,iang,2) ! absorption down
      dTbdKsliq(ilev,iang,1)  = dtb_pl(ilev,iang,3) ! liquid scatter up
      dTbdKsliq(ilev,iang,2)  = dtb_mn(ilev,iang,3) ! liquid scatter down
      dTbdgliq(ilev,iang,1)   = dtb_pl(ilev,iang,4) ! liquid asymmetry up
      dTbdgliq(ilev,iang,2)   = dtb_mn(ilev,iang,4) ! liquid asymmetry down
      dTbdKsrn(ilev,iang,1)   = dtb_pl(ilev,iang,5) ! rain scatter up
      dTbdKsrn(ilev,iang,2)   = dtb_mn(ilev,iang,5) ! rain scatter down
      dTbdgrn(ilev,iang,1)    = dtb_pl(ilev,iang,6) ! rain asymmetry up
      dTbdgrn(ilev,iang,2)    = dtb_mn(ilev,iang,6) ! rain asymmetry down
      dTbdKsice(ilev,iang,1)  = dtb_pl(ilev,iang,7) ! ice scatter up
      dTbdKsice(ilev,iang,2)  = dtb_mn(ilev,iang,7) ! ice scatter down
      dTbdgice(ilev,iang,1)   = dtb_pl(ilev,iang,8) ! ice asymmetry up
      dTbdgice(ilev,iang,2)   = dtb_mn(ilev,iang,8) ! ice asymmetry down
      dTbdKssnow(ilev,iang,1) = dtb_pl(ilev,iang,9) ! snow scatter up
      dTbdKssnow(ilev,iang,2) = dtb_mn(ilev,iang,9) ! snow scatter down
      dTbdgsnow(ilev,iang,1)  = dtb_pl(ilev,iang,10) ! snow asymmetry up
      dTbdgsnow(ilev,iang,2)  = dtb_mn(ilev,iang,10) ! snow asymmetry down
      dTbdKsgrpl(ilev,iang,1) = dtb_pl(ilev,iang,11) ! graupel scatter up
      dTbdKsgrpl(ilev,iang,2) = dtb_mn(ilev,iang,11) ! graupel scatter down
      dTbdggrpl(ilev,iang,1)  = dtb_pl(ilev,iang,12) ! graupel asymmetry up
      dTbdggrpl(ilev,iang,2)  = dtb_mn(ilev,iang,12) ! graupel asymmetry down
    end do ! ilev
  end do ! iang

end subroutine RT_Jacobian
