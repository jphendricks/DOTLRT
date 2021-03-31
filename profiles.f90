!==============================================================
module profiles
!==============================================================
! Defines WRF output variables used as input to MRT
! the variable names match those in WRF output files
!
! History:
!  10/14/2020 Kevin Schaefer created module
!  10/15/2020 Kevin Schaefer added variables to calc temperature, pressure, height
!--------------------------------------------------------------

implicit none

! WRF grid variables
integer nlat     ! (-) number of latitude points
integer nlon     ! (-) number of longitude points
integer nlev_max ! (-) number of vertical levels in WRF output
integer ntime    ! (-) number of times

! temperature variables
real(8), allocatable :: T(:,:,:)      ! (K) perturbation potential temperature
real(8), allocatable :: temp(:,:,:)   ! (K) atmospheric temperature
real(8), allocatable :: TSK(:,:)      ! (K) surface skin temperature

! Variables to calculate height
real(8), allocatable :: PH(:,:,:)     ! (m2/s2) perturbation geopotential
real(8), allocatable :: PHB(:,:,:)    ! (m2/s2) base state geopotential
real(8), allocatable :: height(:,:,:) ! (km) height above sea level

! Variables to calculate pressure
real(8), allocatable :: P(:,:,:)      ! (Pa) perturbation atmospheric pressure
real(8), allocatable :: PB(:,:,:)     ! (Pa) base-state atmospheric pressure
real(8), allocatable :: press(:,:,:)  ! (Pa) total atmospheric pressure

! water related variables
real(8), allocatable :: QVAPOR(:,:,:) ! (g m-3) water vapor mixing ratio
real(8), allocatable :: QCLOUD(:,:,:) ! (g m-3) cloud liquid water mixing ratio
real(8), allocatable :: QICE(:,:,:)   ! (g m-3) ice mixing ratio
real(8), allocatable :: QSNOW(:,:,:)  ! (g m-3) snow mixing ratio
real(8), allocatable :: QGRAUP(:,:,:) ! (g m-3) graupel mixing ratio
real(8), allocatable :: QRAIN(:,:,:)  ! (g m-3) rain mixing ratio

! Land/Lake/Ocean mask variables
real(8), allocatable :: landmask(:,:) ! (-) land vs. ocean/water mask
real(8), allocatable :: lakemask(:,:) ! (-) lake vs. ocean mask

! optional wind speed variables when ocean_mod = 'Wilheit'
real(8), allocatable :: u10m(:,:)     ! (m s-1) u wind component
real(8), allocatable :: v10m(:,:)     ! (m s-1) v wind component
real(8), allocatable :: wind(:,:)     ! (m s-1) wind speed

! Misc variables we read in, but apparently do not use
real(8), allocatable :: PSFC(:,:)     ! (Pa) surface pressure
real(8), allocatable :: XLAT(:,:)     ! (deg) latitudes
real(8), allocatable :: XLONG(:,:)    ! (deg) longitudes
real(8), allocatable :: HGT(:,:)      ! (m) Terrain Height
real(8), allocatable :: q2(:,:)       ! (kg kg-1) water mixing ratio at 2 meters
real(8), allocatable :: t2(:,:)       ! (K) temperature at 2 meters

end module profiles
