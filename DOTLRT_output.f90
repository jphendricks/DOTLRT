!==============================================================
module dotlrt_output
!==============================================================
! Defines variables for all output from mrt used for data assimilation
!
! History:
!   10/15/2020 Kevin Schaefer created module
!   10/20/2020 Kevin Schaefer added netcdf variable specification tree
!   10/27/2020 Kevin Schaefer added all netcdf write variables
!--------------------------------------------------------------

implicit none

! Define ouptput variables for single grid cell
  real(8), allocatable :: Tb_obs_mat(:,:)         ! (K)      (nlev,npol)              Brightness temperature in sensor direction
  real(8), allocatable :: dTb_dT_obs_mat(:,:)     ! (K/K)    (nlev,npol)              Jacobian brightness temp wrt air temperature in sensor direction
  real(8), allocatable :: dTb_dp_obs_mat(:,:)     ! (K/mb)   (nlev,npol)              Jacobian brightness temp wrt pressure in sensor direction
  real(8), allocatable :: dTb_dq_obs_mat(:,:)     ! (K m3/g) (nlev,npol)              Jacobian brightness temp wrt specific humidity in sensor direction
  real(8), allocatable :: dTb_dw_obs_mat(:,:,:)   ! (K m3/g) (nvar,nlev,npol)         Jacobian brightness temp wrt hydrometeor in sensor direction
  real(8)  Tbo_mat(2)                             ! (K)      (npol)                   Brightness temperature at sensor nadir angle 
  real(8)  tau_mat(2)                             ! (-)      (npol)                   opacity at sensor
  real(8), allocatable :: Tbo_str_mat(:,:)        ! (K)      (nstream,npol)           Brightness temperature at top of atmosphere 
  real(8), allocatable :: dTb_dT_str_mat(:,:,:)   ! (K/K)    (nlev,nstream,npol)      Jacobian brightness temp wrt air temperature 
  real(8), allocatable :: dTb_dp_str_mat(:,:,:)   ! (K/mb)   (nlev,nstream,npol)      Jacobian brightness temp wrt pressure
  real(8), allocatable :: dTb_dq_str_mat(:,:,:)   ! (K m3/g) (nlev,nstream,npol)      Jacobian brightness temp wrt specific humidity
  real(8), allocatable :: dTb_dw_str_mat(:,:,:,:) ! (K m3/g) (nvar,nlev,nstream,npol) Jacobian brightness temp wrt hydrometeor 

! define netcdf write arrays for full domain
  real(8), allocatable :: Tbo_wrt(:,:,:)            ! (K)      (nlat,nlon,npol)              Brightness temperature at top of atmosphere in sensor direction
  real(8), allocatable :: Tbo_str_wrt(:,:,:,:)      ! (K)      (nlat,nlon,nstream,npol)      Brightness temperature at top of atmosphere per stream angle 
  real(8), allocatable :: dTb_dp_str_wrt(:,:,:,:,:) ! (K/mb)   (nlat,nlon,nlev,nstream,npol) Jacobian brightness temperature wrt pressure 
  real(8), allocatable :: dTb_dT_str_wrt(:,:,:,:,:) ! (K m3/g) (nlat,nlon,nlev,nstream,npol) Jacobian brightness temperature wrt air temperature 
  real(8), allocatable :: dTb_dq_str_wrt(:,:,:,:,:) ! (K m3/g) (nlat,nlon,nlev,nstream,npol) Jacobian brightness temperature wrt absolute humidity 
  real(8), allocatable :: dTb_dc_str_wrt(:,:,:,:,:) ! (K m3/g) (nlat,nlon,nlev,nstream,npol) Jacobian brightness temperature wrt cloud liquid water 
  real(8), allocatable :: dTb_dr_str_wrt(:,:,:,:,:) ! (K m3/g) (nlat,nlon,nlev,nstream,npol) Jacobian brightness temperature wrt rain 
  real(8), allocatable :: dTb_di_str_wrt(:,:,:,:,:) ! (K m3/g) (nlat,nlon,nlev,nstream,npol) Jacobian brightness temperature wrt ice 
  real(8), allocatable :: dTb_ds_str_wrt(:,:,:,:,:) ! (K m3/g) (nlat,nlon,nlev,nstream,npol) Jacobian brightness temperature wrt snow 
  real(8), allocatable :: dTb_dg_str_wrt(:,:,:,:,:) ! (K m3/g) (nlat,nlon,nlev,nstream,npol) Jacobian brightness temperature wrt graupel 

! output parameters
  integer, parameter :: max_out_var = 100  ! number output variables
  real, parameter :: missing = -999        ! (variable) standard missing value

! output control variables
  integer n_out_var ! number output variables

! Standard variable specification for netcdf
  type variable_spec
     character*20 name       ! (-) netcdf variable name
     character*100 long_name ! (-) long name or description
     character*20 units      ! (-) variable units
     integer ndim            ! (-) number of dimensions
     character*50 dim        ! (-) dimensions
     real(8) Missing         ! (-) missing value
     integer varid           ! (-) netcdf variable id
  end type variable_spec
  type (variable_spec) var(max_out_var)

end module dotlrt_output
