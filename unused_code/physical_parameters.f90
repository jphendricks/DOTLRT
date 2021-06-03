!=================================================================
      module physical_parameters
!=================================================================
! this modules specifices physical parameters and the units
! of those parameters (MKS is standard)

! Modifications
! 1/16/2021 Kevin Schaefer cleaned up code
! 1/16/2021 Kevin Schaefer got rid of kinds module
!-----------------------------------------------------------------

implicit none
save

! biophysical parameters
  real(8), parameter :: pi           = 3.14159265358979323846d0  ! (-) pi
  real(8), parameter :: pidaypy      = 0.0172142d0               ! (?) TBD
  real(8), parameter :: a            = 6.37122d+06               ! earth radius (m)
  real(8), parameter :: grav         = 9.8100d0                  ! (m/s) gravity constant
  real(8), parameter :: gravi        = 1.0d0/grav                ! inverse grav constant
  real(8), parameter :: omega        = 2.0d0*pi/86400.0d0        ! omega = earth angular velocity (1/s)
  real(8), parameter :: earth_area   = 5.100996990707616d+14     ! earth_area = surface area of earth (m2)
  real(8), parameter :: gas_const_R  = 287.000d0                 ! gas_const_R = gas constant for dry air (J/ (kg K))
  real(8), parameter :: spec_heat_cp = 1005.000d0                ! spec_heat_cp = specific heat at constant pressure (J/ (kg K))
  real(8), parameter :: kappa        = gas_const_R/spec_heat_cp  ! kappa = gas_const_R/spec_heat_cp
  real(8), parameter :: inv_kappa    = 1.0d0 / kappa             ! inv_kappa = inverse kappa
  real(8), parameter :: p0_sfc       = 1.0d+05                   ! p0_sfc = surface pressure (Pa)
  real(8), parameter :: inv_p0_sfc   = 1.0d0 / p0_sfc
  real(8), parameter :: tice         = 273.15d0                  ! freezing temperature of water at 1 atm (K)
  real(8), parameter :: hltm         = 2.25d+06                  ! latent heat of vaporization (J/kg)
  real(8), parameter :: gamfac       = hltm*5417.9827d0/spec_heat_cp ! a moist thermodynamic variable
  real(8), parameter :: delta        = 0.608d0                   ! molecular_weight_air/molecular_weight_water - 1
  real(8), parameter :: stefan       = 5.67d-08                  ! stefan boltzmann constant (W/m^2/K^4)
  real(8), parameter :: rv           = 4.61d+02                  ! gas constant for water vapor
  real(8), parameter :: pi2    = 2*pi !     pi2 = 2*pi
  real(8), parameter :: dtr    = pi/180.0d0                      ! Degrees To Radians conversion constant
  real(8), parameter :: rtd    = 180.0d0/pi                      ! Radians To Degrees conversion constant
  real(8), parameter :: alpha2 = 0.0                             !  alpha2 = rotation of axis of rotation from NP/SP
  real(8), parameter :: fPARmax=0.95d0                           ! (-) Max FPAR (98th percentile)
  real(8), parameter :: fPARmin=0.01d0                           ! (-) Min FPAR (2nd percentile)

end module physical_parameters

