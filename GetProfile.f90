subroutine getprofile( atminp )
  use variables
  implicit none
  integer i
  real(8) atminp(atm_inp%num_levels+1,9)

  character*120 debugout
  
  ! surface at first level: height = 0 (agl)
  surf_inp%surf_temp = atminp(1,3)

  do i = 1, atm_inp%num_levels
    atm_inp%prof(i)%height          = atminp(i+1,1) ! (km)
    atm_inp%prof(i)%pressure        = atminp(i+1,2) ! (mb)
    atm_inp%prof(i)%temperature     = atminp(i+1,3) ! (K)
    atm_inp%prof(i)%vapor_density   = atminp(i+1,4) ! (g/m^3)
    atm_inp%prof(i)%cloud_liq_dens  = atminp(i+1,5) ! (g/m^3)
    atm_inp%prof(i)%cloud_rn_dens   = atminp(i+1,6) ! (g/m^3)
    atm_inp%prof(i)%cloud_ice_dens  = atminp(i+1,7) ! (g/m^3)
    atm_inp%prof(i)%cloud_snow_dens = atminp(i+1,8) ! (g/m^3)
    atm_inp%prof(i)%cloud_grpl_dens = atminp(i+1,9) ! (g/m^3)
  end do ! i
  atm_inp%inf= .true.
  atm_inp%fiveph= .true.

 end subroutine getprofile
