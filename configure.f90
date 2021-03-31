!=================================================================================
subroutine setup_all_outputs( )
!=================================================================================
! This routine sets static parameters for DOTLRTv1.5
! and checks for bad inputs
!
! History:
!  5/27/2005 Bob Weber created routine
!  9/26/2020 Kevin Schaefer deleted unused variables
!  10/15/2020 Kevin Schaefer deleted all arguments duplicated in variables module 
!---------------------------------------------------------------------------------
use dotlrt_variables
use profiles
use dotlrt_output
implicit none

! Local variables
  logical bad

bad= .false.

!print message
  print*, 'Set up all Outputs'

! save single atmospheric profile
  if (save_sing_prof) then
    if(trim(prof_src)=='WRF') then
      call write_text_profile()
    else
      print*, 'save_sing_prof only valid if prof_src = WRF'
      stop
    endif
  endif

! Set up output radiation files
  call define_variables

! set up output write arrays
  if (save_rad_file) then
    allocate(Tbo_wrt(nlon,nlat,npol))
    allocate(Tbo_str_wrt(nlon,nlat,nstream_surf,npol))
    allocate(dTb_dp_str_wrt(nlat,nlon,nlev,nstream_surf,npol)) ! Jacobian brightness temperature wrt pressure 
    allocate(dTb_dT_str_wrt(nlat,nlon,nlev,nstream_surf,npol)) ! Jacobian brightness temperature wrt air temperature 
    allocate(dTb_dq_str_wrt(nlat,nlon,nlev,nstream_surf,npol)) ! Jacobian brightness temperature wrt absolute humidity 
    allocate(dTb_dc_str_wrt(nlat,nlon,nlev,nstream_surf,npol)) ! Jacobian brightness temperature wrt cloud liquid water 
    allocate(dTb_dr_str_wrt(nlat,nlon,nlev,nstream_surf,npol)) ! Jacobian brightness temperature wrt rain 
    allocate(dTb_di_str_wrt(nlat,nlon,nlev,nstream_surf,npol)) ! Jacobian brightness temperature wrt ice 
    allocate(dTb_ds_str_wrt(nlat,nlon,nlev,nstream_surf,npol)) ! Jacobian brightness temperature wrt snow 
    allocate(dTb_dg_str_wrt(nlat,nlon,nlev,nstream_surf,npol)) ! Jacobian brightness temperature wrt graupel 

    Tbo_wrt = missing
    Tbo_str_wrt = missing
    dTb_dp_str_wrt = missing
    dTb_dT_str_wrt = missing
    dTb_dq_str_wrt = missing
    dTb_dc_str_wrt = missing
    dTb_dr_str_wrt = missing
    dTb_di_str_wrt = missing
    dTb_ds_str_wrt = missing
    dTb_dg_str_wrt = missing
    
    call create_tb_file()
  endif

! write profile file
  if (save_prof_file) call create_profile_file()

end subroutine setup_all_outputs

!=================================================================================
subroutine setup_all_inputs( )
!=================================================================================
! This routine sets static parameters for DOTLRTv1.5
! and checks for bad inputs
!
! History:
!  5/27/2005 Bob Weber created routine
!  9/26/2020 Kevin Schaefer deleted unused variables
!  10/15/2020 Kevin Schaefer deleted all arguments duplicated in variables module 
!---------------------------------------------------------------------------------
use dotlrt_variables
use profiles
use dotlrt_output
implicit none

! Local variables
  integer kpts
  integer kind
  logical bad
  integer iang
  real(8) endpts(2)
  real(8) alpha
  real(8) beta
  real(8) b(32)
  real(8) quad(32)

!print message
  print*, 'Set up all Inputs'

! Read namel file
  call read_namel_dotlrt
! 
! check source
  bad=.true.
  if (trim(prof_src) == 'single' ) bad = .false.
  if (trim(prof_src) == 'WRF') bad = .false.
  if (bad) then
    print *, trim(prof_src)
    stop 'incorrect prof_src'
  endif

! get number of levels
  if (trim(prof_src) == 'single' ) then
    open(unit=20, file=trim(file_in), form='formatted', status='old')
    read(20,*) nlev
    close(unit=20)
    nlon=1
    nlat=1
  elseif (trim(prof_src) == 'WRF') then
    call read_WRF_netcdf_dimensions(file_in)
    if (lon_stop>nlon) lon_stop=nlon
    if (lat_stop>nlat) lat_stop=nlat
  endif

! nlev checks
  bad=.false.
  if (lon_strt>nlon) bad = .true.
  if (lon_stop>nlon) bad = .true.
  if (lat_strt>nlat) bad = .true.
  if (lat_stop>nlat) bad = .true.
  if (bad) then
    print*, 'bad grid parameters'
    print*, 'lon_strt: ', lon_strt, 'lon_stop: ', lon_stop
    print*, 'lat_strt: ', lat_strt, 'lat_stop: ', lat_stop
    stop
  endif

! assign nlev 
  if(flag_reduce_nvar) then
    nlev=new_nlev
  else
    nlev=nlev_max
  endif
  
! instrument specification parameters input file:
  call get_instr_spec( )

  nang = nstream ! number of stream angles, up and down (even)
  nangover2=nang/2
  if( 2*nangover2 /= nang ) then
    write(8,*) 'nang is odd integer, program STOPPED'
    stop
  end if

! stream angles and weights for Gaussian-Lobatto quadrature, including
! streams at 0 and 180 degrees
    kpts = 2
    kind = 1
    endpts(1) = -1.0D0
    endpts(2) = +1.0D0
    call gaussq(kind, nang, alpha, beta, kpts, endpts, b, quad, quad_wts)
    quad(1) = -1.0d0
    quad(nang) = +1.0d0
    do iang = 1, nang
      quad_ang(iang) = (180.0d0/pi)*acos(quad(nang-iang+1))
      cos_ang(iang)= dcos(quad_ang(iang)/180.d0*pi)
      sin_ang(iang)= dsin(quad_ang(iang)/180.d0*pi)
    enddo

! Henyey-Greenstein phase matrix
    call HG_phmat()

! read atmosphere profiles
  if (trim(prof_src) == 'single' ) then
    call read_text_profile()
  elseif (trim(prof_src) == 'WRF') then
    call get_wrf_data()
  endif

! create surface reflectance tables
  if (trim(prof_src) == 'WRF') then
    call construct_surf_ref_table()
  endif

! Allocate brightness temperature and jacobian output variables
  allocate(Tb_obs_mat(0:nlev,2))
  allocate(dTb_dT_obs_mat(nlev,2))
  allocate(dTb_dp_obs_mat(nlev,2))
  allocate(dTb_dq_obs_mat(nlev,2))
  allocate(dTb_dw_obs_mat(nlev,5,2))
  allocate(Tbo_str_mat(nstream/2, 2))
  allocate(dTb_dT_str_mat(nlev,nstream/2,2))
  allocate(dTb_dp_str_mat(nlev,nstream/2,2))
  allocate(dTb_dq_str_mat(nlev,nstream/2,2))
  allocate(dTb_dw_str_mat(nlev,nstream/2,5,2))

end subroutine setup_all_inputs
