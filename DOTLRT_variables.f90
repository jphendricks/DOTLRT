! ==============================================================
module dotlrt_variables
! ==============================================================
! Defines all input and internal variables for MRT
!
! History:
!  9/29/2020  Kevin Schaefer cleaned up code
!  10/5/2020  Kevin Schaefer cleaned up surf_inp_type
!  10/14/2020 Kevin Schaefer added execution control variables
!  10/15/2020 Kevin Schaefer added grav, gamma, cp; changed rgas
!  10/20/2020 Kevin Schaefer removed unused variables
!  11/3/2020  Kevin Schaefer removed unused hydromet variables
!  11/3/2020  Kevin Schaefer added hydromet cluster variables
!  12/1/2020  Kevin Schaefer added permmitivity of free space as constant
!  12/12/2020 Kevin Schaefer removed all but one nlev variable
!  1/24/2021  Kevin Schaefer restructured reduced dimension branch
!  5/16/2021  Kevin Schaefer changed d_albedo = 1.d-7 to albedo_diag=1.0d-12
! --------------------------------------------------------------

implicit none

! Execution control variables
integer chan_strt ! start index for channels
integer chan_stop ! stop index for channels
integer lon_strt  ! start index for longitude
integer lon_stop  ! stop index for longitude
integer lat_strt  ! start index for latitude
integer lat_stop  ! stop index for latitude
logical flag_print_full  ! (-) flag to print full messages
logical flag_read_anc    ! (-) flag to read in ancillary data from WRF
logical flag_reduce_nvar ! (-) flag to use only part of atmos profile from WRF

logical gen_index_table  ! (-) flag to use only part of atmos profile from WRF
logical use_index_table  ! (-) flag to use only part of atmos profile from WRF
logical dbg_index_table  ! (-) flag to use only part of atmos profile from WRF

integer new_nlev         ! (-) number of levels to use in reduced profile
character(20) prof_src   ! (-) source of atmospheric profile data
                         ! 'single' single column text file
                         ! 'WRF'    WRF output netcdf file
character(20) ocean_mod  ! (-) ocean model type
                         ! 'Fresnel' Fresnel reflectivity model
                         ! 'Wilheit' Wilheit reflectivity model
                         ! 'ITRA'    ITRA reflectivity model

! real(8), allocatable :: dens(:)
! real(8), allocatable :: temp(:)

! file names
character*250 file_in          ! (-) path to input atm profile
character*250 file_instr       ! (-) path to instrument spec file
character*250 file_var         ! (-) path to output variable definition file
character*250 file_index_table ! (-) path to output variable definition file
character*250 out_path         ! (-) path to output directory

! output control variables
logical save_rad_file  ! (-) save radiation output as netcdf file
logical save_jac_file  ! (-) save jacobian output as netcdf file
logical save_prof_file ! (-) save atmospheric profiles as netcdf file
logical save_sing_prof ! (-) save single atmospheric profile as text file
logical save_sing_rad  ! (-) save single atmospheric radiation output as text file
integer save_ilon      ! (-) longitude index of profile to save
integer save_ilat      ! (-) latitude index of profile to save

! execution time diagnostic variables
integer time_rate   ! conversion between count and clock time
integer time_start  ! time count at start
integer time_temp   ! time count at start
integer time_stop   ! time count at endreal(8)
logical print_ex    ! (-) flag to print execution times
integer nseg        ! (-) number of code segments
real(8) time_del          ! (s) delta time for individual call
real(8) time_seg(20)      ! (s) execution times for code segments
real(8) num_call(20)      ! (-) execution times for code segments
character*40 seg_name(20) ! (-) code segment names
character*40 name         ! (-) code segment name

! dimension parameters
integer, parameter :: max_nlev = 75           ! (-) Max number of input press levels to define atmos. state
integer, parameter :: max_nang = 16           ! (-) Max number of angles to define surface reflectivity
integer, parameter :: max_nchan = 64          ! (-) Maximum number of instrument channels
integer, parameter :: max_nfreq = 45          ! (-) Max num of frequencies for spectral passband integration
integer, parameter :: max_nvar = 16           ! (-) max number of variables to estimate (not same as number of phases)
integer, parameter :: max_nstream = 32        ! (-) max num quadrature angles
integer, parameter :: max_nphase = 5          ! (-) max number number hydrometeor phases
integer, parameter :: nphase = 5              ! (-) number of hydrometeor phases
integer, parameter :: nvar = 2 + 2 * nphase   ! (-) number of variables for radiation Jacobian
integer, parameter :: nvar_prof = 9           ! (-) number of variables in atmospheric profile
integer, parameter :: npol = 2                ! (-) number of polarizations
integer, parameter :: npseudo = 1             ! (-) number of pseudo-layers in each layer (fixed at 1: original layers)
integer, parameter :: nHG = 101               ! (-) number HG phase matrix frequencies
integer, parameter :: updown = 2              ! (-) number major directions (1 = up 2 = down)

! physical constant parameters
real(8), parameter :: grav = 9.80655d0        ! (m/s2) acceleration of gravity
real(8), parameter :: gravi        = 1.0d0/grav                ! (s/m) inverse grav constant
real(8), parameter :: hplank = 6.6252d-34     ! (J s) Plank's constant
real(8), parameter :: kboltz = 1.380662d-23   ! (J/K) Boltzmann's constant
real(8), parameter :: pi= 3.14159265358979d0  ! (-) value of pi
real(8), parameter :: rgas = 287.04d0         ! (J/K/kg) Ideal gas constant
real(8), parameter :: air_cp = 1004.d0        ! (J/K/kg) heat capacity of air (Note: not using Bolton's value of 1005.7)
real(8), parameter :: gamma=rgas/air_cp       ! (-) exponent to convert potential temperature to temperature
real(8), parameter :: press_ref=100000d0      ! (Pa) reference pressure for potential temperature
real(8), parameter :: base_pot=300.d0         ! (K) base state potential temperature
real(8), parameter :: speed_lt = 2.99792458d8 ! (m/s) speed of light
real(8), parameter :: e_free = 8.854e-12      ! (F/m) permmitivity of free spce
real(8), parameter :: max_freq = 1000.0d0     ! (GHz) Max frequency for absorption computations
real(8), parameter :: t_cosmic = 2.730d0      ! (K) cosmic background temperature
real(8), parameter :: land_refl = 0.05d0      ! (-) reference surface reflectivity for land
real(8), parameter :: water_refl = 0.5d0      ! (-) reference surface reflectivity for water (lake and ocean)
real(8), parameter :: sal_ocean = 0.035d0     ! (-) ocean reference salinity fraction (avg = 0.035)
real(8), parameter :: sal_lake = 0.d0         ! (-) freshwater lake salinity fraction
real(8), parameter :: d_albedo=1.0d-7         ! (-) lower limit scat/abs total for diagonalization layer inst parameters
real(8), parameter :: albedo_diag=1.0d-12     ! (-) lower limit on albedo to assume diagonal matrices
real(8), parameter :: t_frz=273.15d0          ! (K) freezing point of water

! quadrature parameters
integer nlev         ! (-) number of vertical levels
integer obs_lev
integer nstream      ! (-) total number of stream angles, including upward and downward propagating
integer nstream_surf ! (-) total number of stream angles for surface, upward only (=nstream/2)
integer nang         ! (-) equal to nstream
integer nangover2    ! (-) equal to nstream_surf
real(8) quad_ang(max_nstream) ! (-) quadrature angles
real(8) quad_wts(max_nstream) ! (-) quadrature weights
real(8) cos_ang(max_nstream)  ! (-) cosines of the quadrature angles
real(8) sin_ang(max_nstream)  ! (-) sines of the quadrature angles

real(8) obs_theta    ! (deg) zenith angle to observation
real(8) obs_height   ! (km) height of observation
real(8) ipol         ! (-) polarization number 1=vertical 0=horizontal
integer num_freqs    ! (-) number of freq points for passband quadrature
integer nsub_freq    ! (-) Number of frequency points in a sideband for passband quadrature
real(8) g_asymmetry(max_nlev)
real(8) passband_freq(max_nfreq)
real(8) testvar1 ! test diagnostic variable 1
real(8) testvar2 ! test diagnostic variable 2
real(8) testvar3 ! test diagnostic variable 3
real(8) testval(10) ! test diagnostic variable array

type gas_type
    real(8) absn2      ! (nepers/km) nitrogen absorption
    real(8) absh2o     ! (nepers/km) water absorption
    real(8) o2abs      ! (nepers/km) oxygen absorption
    real(8) dabsn2_dt  ! temperature derivative
    real(8) dabsn2_dp  ! pressure derivative
    real(8) dabsh2o_dt ! temperature derivative
    real(8) dabsh2o_dp ! pressure derivative
    real(8) dabsh2o_dw ! water vapor derivative
    real(8) do2abs_dt  ! temperature derivative
    real(8) do2abs_dp  ! pressure derivative
    real(8) do2abs_dw  ! water vapor derivative
end type gas_type
type (gas_type), dimension(max_nlev) :: gas_prof

type hydro_type
    real(8) cloudab         ! cloud particle absorption in nepers/km
    real(8) cloudsc         ! cloud particle scattering in nepers/km
    real(8) cloudg          ! cloud particle asymmetry factor (unitless)
    real(8) dcloudab_dt     ! temperature derivative
    real(8) dcloudab_dk0    ! k0 derivative
    real(8) dcloudab_da0    ! a0 derivative

    real(8) dcloudsc_dt     ! temperature derivative
    real(8) dcloudsc_dk0    ! k0 derivative
    real(8) dcloudsc_da0    ! a0 derivative

    real(8) dcloudg_dt      ! temperature derivative
    real(8) dcloudg_dk0     ! k0 derivative
    real(8) dcloudg_da0     ! a0 derivative ! not temperature or microwave frequency
end type hydro_type
type (hydro_type), dimension(max_nlev,max_nphase) :: hydro_prof

! brightness temperature and Jacobian in observation direction
real(8), allocatable :: Tb_obs_mat(:,:)         ! (K)      (nlev,npol)           Brightness temperature in sensor direction
real(8), allocatable :: dTb_dT_obs_mat(:,:)     ! (K/K)    (nlev,npol)           Jacobian brightness temp wrt air temperature in sensor direction
real(8), allocatable :: dTb_dp_obs_mat(:,:)     ! (K/mb)   (nlev,npol)           Jacobian brightness temp wrt pressure in sensor direction
real(8), allocatable :: dTb_dq_obs_mat(:,:)     ! (K m3/g) (nlev,npol)           Jacobian brightness temp wrt specific humidity in sensor direction
real(8), allocatable :: dTb_dw_obs_mat(:,:,:)   ! (K m3/g) (nvar,nlev,npol)      Jacobian brightness temp wrt hydrometeor in sensor direction
real(8)  Tbo_mat(2)                             ! (K)      (npol)                Brightness temperature at sensor nadir angle
real(8)  tau_mat(2)                             ! (-)      (npol)                opacity at sensor

! brightness temperature and Jacobian as function of stream angle
real(8), allocatable :: Tbo_str_mat(:,:)        ! (K)      (nstream,npol)        Brightness temperature at top of atmosphere
real(8), allocatable :: dTb_dT_str_mat(:,:,:)   ! (K/K)    (nlev,nstream,npol)   Jacobian brightness temp wrt air temperature
real(8), allocatable :: dTb_dp_str_mat(:,:,:)   ! (K/mb)   (nlev,nstream,npol)   Jacobian brightness temp wrt pressure
real(8), allocatable :: dTb_dq_str_mat(:,:,:)   ! (K m3/g) (nlev,nstream,npol)   Jacobian brightness temp wrt specific humidity
real(8), allocatable :: dTb_dw_str_mat(:,:,:,:) ! (K m3/g) (nlev,nstream,5,npol) Jacobian brightness temp wrt hydrometeor

! geophysical Jacobian
real(8), dimension(max_nlev) :: dKab_dT
real(8), dimension(max_nlev) :: dKab_dp
real(8), dimension(max_nlev) :: dKab_dq
real(8), dimension(max_nlev,max_nphase) :: dKsc_dT
real(8), dimension(max_nlev,max_nphase) :: dg_dT
real(8), dimension(max_nlev,max_nphase) :: dKab_dw
real(8), dimension(max_nlev,max_nphase) :: dKsc_dw
real(8), dimension(max_nlev,max_nphase) :: dg_dw

! total Jacobian
real(8), dimension(max_nlev,max_nstream,2) :: dTb_dT  ! (K K-1) Jacobian wrt air temperature
real(8), dimension(max_nlev,max_nstream,2) :: dTb_dp  ! (K mb-1) Jacobian wrt air pressure
real(8), dimension(max_nlev,max_nstream,2) :: dTb_dq  ! (K m3 g-1) Jacobian wrt air humidity
real(8), dimension(max_nlev,max_nstream,max_nphase,2) :: dTb_dw ! (K m3 g-1) Jacobian wrt hydrometeors

! radiative transfer Jacobian
real(8), dimension(max_nlev,max_nstream,2) :: dTbdTr      ! Jacobian temperature
real(8), dimension(max_nlev,max_nstream,2) :: dTbdKa      ! Jacobian absorption
real(8), dimension(max_nlev,max_nstream,2) :: dTbdKsliq   ! Jacobian scatter liquid
real(8), dimension(max_nlev,max_nstream,2) :: dTbdgliq    ! Jacobian asymmetry liquid
real(8), dimension(max_nlev,max_nstream,2) :: dTbdKsrn    ! Jacobian scatter rain
real(8), dimension(max_nlev,max_nstream,2) :: dTbdgrn     ! Jacobian asymmetry rain
real(8), dimension(max_nlev,max_nstream,2) :: dTbdKsice   ! Jacobian scatter ice
real(8), dimension(max_nlev,max_nstream,2) :: dTbdgice    ! Jacobian asymmetry ice
real(8), dimension(max_nlev,max_nstream,2) :: dTbdKssnow  ! Jacobian scatter snow
real(8), dimension(max_nlev,max_nstream,2) :: dTbdgsnow   ! Jacobian asymmetry snow
real(8), dimension(max_nlev,max_nstream,2) :: dTbdKsgrpl  ! Jacobian scatter graupel
real(8), dimension(max_nlev,max_nstream,2) :: dTbdggrpl   ! Jacobian asymmetry graupel

! Hydrometeor distribution state variable parameter branch
type hydrometeor_characteristics
    real(8) dens     ! (g/m3) hydrometeor density
    real(8) p        ! (?) hydrometeor particle size distribution exponent
    real(8) q        ! (?) hydrometeor particle size variance parameter
    real(8) k0       ! (num/m3/mm) reference hydrometeor number density in radius
    real(8) a0       ! (mm) mean hydrometeor radius
    real(8) dk0_dw   ! (?) k0 derivative wrt hydrometeor density
    real(8) da0_dw   ! (?) a0 derivative wrt hydrometeor density
    real(8) ftot     ! (-) fraction of total column value
    real(8) fprecip  ! (-) fraction of total precipitation per layer (rain + snow + graupel)
    real(8) fcloud   ! (-) fraction of total cloud per layer (CLW + ice)
    real(8) fhydro   ! (-) fraction of total hydrometeor per layer (CLW + ice + rain + snow + graupel)
    logical a0_const ! (-) true if density is zero
end type hydrometeor_characteristics
type(hydrometeor_characteristics) gen_hm

! atmospheric profile branch atm%
type profile_type
    ! Gaseous state variables
    real(8) press       ! (mb) pressure
    real(8) temp        ! (K) temperature
    real(8) humid       ! (g/m3) water vapor density
    real(8) f_froz      ! (-) frozen fraction of hydrometeors

    ! Derived radiative transfer quantities : frequency dependent
    real(8) abs_o2      ! oxygen and nitrogen absorption in nepers/km
    real(8) abs_h2o     ! water vapor absorption in nepers/km
    real(8) abs_cloud   ! cloud particle absorption in nepers/km
    real(8) scat_cloud  ! cloud particle scattering in nepers/km
    real(8) asymmetry   ! cloud particle asymmetry factor (unitless)
    real(8) ext_tot     ! total atmospheric extinction in nepers/km
    real(8) albedo      ! single scattering albedo
    real(8) bb_spec_int ! (W/m2/Hz/ster) black body intensity for either polarization at level in W/m**2/Hz/ster

    ! layer geometry variables
    real(8) hgt_mid     ! (km) height to middle of layer wrt to surface
    real(8) hgt_top     ! (km) height to top of layer wrt to surface
    real(8) hgt_bot     ! (km) height to bottom of layer wrt to surface
    real(8) hgt_del     ! (m) layer thickness (note the units are meters not km)

    ! Hydrometeor distribution state variable parameters
    type(hydrometeor_characteristics) clw
    type(hydrometeor_characteristics) rain
    type(hydrometeor_characteristics) ice
    type(hydrometeor_characteristics) snow
    type(hydrometeor_characteristics) grpl
    type(hydrometeor_characteristics) cloud
    type(hydrometeor_characteristics) precip
    type(hydrometeor_characteristics) hydro
end type profile_type
type(profile_type) atm(max_nlev)
type(profile_type) atm_temp(max_nlev)
type(profile_type) atm_psuedo(max_nlev)
type(profile_type) atm_init(max_nlev)

! hydrometeor characteristic and geometry branch
type hm_geometry_type
    real(8) tot     ! (g/m2) total column hydrometeor
    real(8) min     ! (g/m3) minimum hydrometeor value in cloud
    real(8) max     ! (g/m3) maximum hydrometeor value in cloud
    real(8) ave     ! (g/m3) average hydrometeor value in cloud
    real(8) std     ! (g/m3) standard deviation hydrometeor value in cloud
    real(8) frac    ! (-) fraction of total hydrometeor in cloud
    real(8) h_frz   ! (km) freezing height of cloud
    real(8) h_ave   ! (km) average height of cloud
    real(8) h_std   ! (km) height standard deviation of cloud
    real(8) h_top   ! (km) top of cloud
    real(8) h_bot   ! (km) bottom of cloud
    real(8) h_del   ! (km) thickness of cloud
    real(8) h_max   ! (km) height of max value in cloud
end type hm_geometry_type
type(hm_geometry_type) geom

! main cloud geometry branch
type atm_cloud_type
    type(hm_geometry_type) clw
    type(hm_geometry_type) rain
    type(hm_geometry_type) ice
    type(hm_geometry_type) snow
    type(hm_geometry_type) grpl
    type(hm_geometry_type) cloud
    type(hm_geometry_type) precip
    type(hm_geometry_type) hydro
end type atm_cloud_type
    type(atm_cloud_type) cld ! cloud geometry branch

! surface characteristics branch
type surface_charateristics
     integer(4) type         ! (-) surface type flag 2=land, 1=lake, 0=ocean
     character*20 ocean_mod  ! (-) ocean model: 'Fresnel','Wilheit',or 'ITRA'
     character*20 diel_mod   ! (-) water dielectric model: 'Salt','NaCl',or 'Fresh'
     real(8) temp            ! (K) surface skin temperature
     real(8) theta(max_nang) ! (deg) angle of incidence or reflection wrt surface normal
     real(8) vref(max_nang)  ! (-) vertically polarized surface reflectivity
     real(8) href(max_nang)  ! (-) horizontally polarized surface reflectivity
     real(8) sal             ! (-) water salinity ratio
     double complex diel     ! (-) dielectric constant
     real(8) wind            ! (m/s) surface wind speed
     integer nstream         ! (-) number of surface stream angles
end type surface_charateristics
type (surface_charateristics) :: surf
real(8), dimension(max_nang) :: surf_reflec ! (-) generic surface reflectivity for local calculations

! Instrument specifications branch
integer nchan  ! (-) number of channels
type chan_spec_type
     real(8) lo_freq    ! (GHz) LO frequency
     real(8) if1_freq   ! (GHz) IF1 frequency
     real(8) if2_freq   ! (GHz) IF2 frequency
     real(8) bandwidth  ! (GHz) bandwidth
     real(8) dtrms      ! (K) observation noise
     integer desig      ! (-) channel designation
     integer num        ! (-) channel number
     character(50) name ! (-) instrument or satellite name
end type chan_spec_type
type (chan_spec_type), dimension(max_nchan) :: instr_spec
type (chan_spec_type) :: channel

real(8), dimension(max_nlev) :: temperature
real(8), dimension(max_nlev) :: abs_O2
real(8), dimension(max_nlev) :: abs_H2O
real(8), dimension(max_nlev) :: abs_cloud
real(8), dimension(max_nlev) :: scat_cloud
real(8), dimension(max_nlev) :: temperature1
real(8), dimension(max_nlev) :: abs_O2_1
real(8), dimension(max_nlev) :: abs_H2O_1
real(8), dimension(max_nlev) :: abs_cloud1
real(8), dimension(max_nlev) :: scat_cloud1
real(8), dimension(max_nlev) :: al_gas1 ! (-) gas absorbtion (1/length) at a layer
real(8), dimension(max_nlev) :: dtemperature1
real(8), dimension(max_nlev) :: dabs_cloud1
real(8), dimension(max_nlev) :: dscat_cloud1
real(8), dimension(max_nlev) :: dal_gas1
real(8), dimension(max_nlev) :: h1
real(8), dimension(max_nlev) :: abs_total
real(8), dimension(max_nlev) :: abs_total1

real(8), dimension(nHG) :: HGg
real(8), dimension(nHG,max_nstream,max_nstream) :: HGph
real(8), dimension(nHG,max_nstream,max_nstream) :: dHGph

real(8), dimension(max_nlev+1) :: altitude ! (m) altitude above sea level
real(8), dimension(max_nlev+1) :: altitude1 ! (m) altitude above sea level

real(8), dimension(0:max_nlev,max_nstream) :: tb_pl ! (K) upwelling brightness temperatures at boundaries between layers
real(8), dimension(0:max_nlev,max_nstream) :: tb_mn ! (K) downwelling brightness temperatures at boundaries between layers
real(8), dimension(0:max_nlev,max_nstream,max_nvar):: dtb_pl ! (K) incremental btightness temperature
real(8), dimension(0:max_nlev,max_nstream,max_nvar):: dtb_mn ! (K/unit) derivatives of the Tb profiles with respect to a parameter

real(8), dimension(max_nlev,max_nstream,max_nstream) :: phase11
real(8), dimension(max_nlev,max_nstream,max_nstream) :: dphase11
real(8), dimension(max_nlev,max_nstream,max_nstream) :: phaseff  ! (?) scattering matrix in forward direction
real(8), dimension(max_nlev,max_nstream,max_nstream) :: phasefb  ! (?) scattering matrix in backward direction
real(8), dimension(max_nlev,max_nstream,max_nstream) :: dphaseff ! (?) derivatives of phase functions  forward scattering
real(8), dimension(max_nlev,max_nstream,max_nstream) :: dphasefb ! (?) derivatives of phase functions backward scattering
real(8), dimension(max_nlev,max_nstream,max_nstream) :: phaseff1
real(8), dimension(max_nlev,max_nstream,max_nstream) :: phasefb1
real(8), dimension(max_nlev,max_nstream,max_nstream) :: dphaseff1
real(8), dimension(max_nlev,max_nstream,max_nstream) :: dphasefb1
real(8), dimension(max_nlev,max_nstream,max_nstream) :: phase11_sc
real(8), dimension(max_nlev,max_nstream,max_nstream) :: phaseff_sc
real(8), dimension(max_nlev,max_nstream,max_nstream) :: phasefb_sc
real(8), dimension(max_nlev,max_nstream,max_nstream) :: phaseff1_sc
real(8), dimension(max_nlev,max_nstream,max_nstream) :: phasefb1_sc
real(8), dimension(max_nlev,max_nstream,max_nstream,max_nphase) :: dphaseff_g
real(8), dimension(max_nlev,max_nstream,max_nstream,max_nphase) :: dphasefb_g
real(8), dimension(max_nlev,max_nstream,max_nstream,max_nphase) :: dphaseff_sc
real(8), dimension(max_nlev,max_nstream,max_nstream,max_nphase) :: dphasefb_sc
real(8), dimension(max_nlev,max_nstream,max_nstream,max_nphase) :: d_phase11_g
real(8), dimension(max_nlev,max_nstream,max_nstream,max_nphase) :: d_phase11_sc
real(8), dimension(max_nlev,max_nstream,max_nstream,max_nphase) :: dphaseff1_g
real(8), dimension(max_nlev,max_nstream,max_nstream,max_nphase) :: dphasefb1_g
real(8), dimension(max_nlev,max_nstream,max_nstream,max_nphase) :: dphaseff1_sc
real(8), dimension(max_nlev,max_nstream,max_nstream,max_nphase) :: dphasefb1_sc

end module dotlrt_variables
