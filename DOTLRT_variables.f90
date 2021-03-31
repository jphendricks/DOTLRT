! ==============================================================
module dotlrt_variables
! ==============================================================
! Defines all input and internal variables for MRT
! 
! History:
!  9/29/2020 Kevin Schaefer cleaned up code
!  10/5/2020 Kevin Schaefer cleaned up surf_inp_type
!  10/14/2020 Kevin Schaefer added execution control variables
!  10/15/2020 Kevin Schaefer added grav, gamma, cp; changed rgas
!  10/20/2020 Kevin Schaefer removed unused variables
!  11/3/2020 Kevin Schaefer removed unused hydromet variables
!  11/3/2020 Kevin Schaefer added hydromet cluster variables
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
integer new_nlev         ! (-) number of levels to use in reduced profile
character(20) prof_src   ! (-) source of atmospheric profile data 
                         ! 'single' single column text file
                         ! 'WRF'    WRF output netcdf file
character(20) ocean_mod  ! (-) ocean model type 
                         ! 'Fresnel' Fresnel reflectivity model
                         ! 'Wilheit' Wilheit reflectivity model
                         ! 'ITRA'    ITRA reflectivity model

real(8) time_k, time_rt, time_io, times(2)

! file names
character*250 file_instr     ! (-) input instrument spec filename
character*250 file_single_in ! (-) input ascii single profile filename
character*250 file_wrf       ! (-) input wrf filename
character*250 file_var       ! (-) input variable definition filename
character*250 out_path       ! (-) output directory

! output control variables
logical save_rad_file  ! (-) save radiation output as netcdf file file
logical save_sing_prof ! (-) save single atmospheric profile as text file
logical save_sing_rad  ! (-) save single atmospheric radiation output as text file
integer save_ilon      ! (-) longitude index of profile to save
integer save_ilat      ! (-) latitude index of profile to save

! constant parameters
integer, parameter :: max_nlev = 75             ! Max # of input press levels to define atmos. state
integer, parameter :: max_num_angles = 16       ! Max # of angles to define surface reflectivity
real(8), parameter :: grav = 9.80655d0          ! (m/s2) acceleration of gravity
real(8), parameter :: hplank = 6.6252d-34       ! Plank's constant in J-s
real(8), parameter :: kboltz = 1.380662d-23     ! Boltzmann's constant in J/K
real(8), parameter :: pi= 3.14159265358979d0    ! not encircling any gravitational singularities
real(8), parameter :: rgas = 287.04d0           ! (J/K/kg) Ideal gas constant
real(8), parameter :: cp = 1004.d0              ! (J/K/kg) heat capacity of air (Note: not using Bolton's value of 1005.7)
real(8), parameter :: gamma=rgas/cp             ! (-) exponent to convert potential temperature to temperature
real(8), parameter :: press_ref=100000d0        ! (Pa) reference pressure for potential temperature
real(8), parameter :: base_pot_temp=300.d0      ! (K) base state potential temperature
real(8), parameter :: quad_weight_error =1.0d-5 ! Allowable unitarity error for all quadrature weights
real(8), parameter :: c = 2.99793d8             ! speed of light (m/s)
real(8), parameter :: std_temperature= 273.15d0 ! Standard temperature in K
real(8), parameter :: std_pressure = 1013.25d0  ! Standard pressure in mb
real(8), parameter :: max_freq = 1000.0d0       ! Max frequency (GHz) for absorption computations
real(8), parameter :: temp_incr = 0.01d0        ! absorption coeff, temperature inc weight functions
real(8), parameter :: wv_incr = 0.01d0          ! water vapor inc for calc deriv in absorption coeff
real(8), parameter :: btcb = 2.73d0             ! brightness temperature of cosmic background
real(8), parameter :: rad_earth = 6.3568d3      ! (km) Radius of the earth
real(8), parameter :: raythresh=0.1d0           ! Rayleigh scattering and extinction efficiency
real(8), parameter :: t_cb = 2.730d0            ! (K) cosmic background temperature
real(8), parameter :: land_refl = 0.05d0        ! (-) reference surface reflectivity for land
real(8), parameter :: water_refl = 0.5d0        ! (-) reference surface reflectivity for water (lake and ocean)
real(8), parameter :: salinity_ocean = 0.035d0  ! (-) ocean reference salinity fraction (avg = 0.035)
real(8), parameter :: salinity_lake = 0.d0      ! (-) freshwater lake salinity fraction
real(8), parameter :: d_albedo=1.0d-7           ! (-) lower limit scat/abs total for diagonalization layer inst parameters
integer, parameter :: max_nchannel = 64         ! (-) Maximum number of instrument channels
integer, parameter :: max_num_int_freqs = 45    ! Max num of spectral integration frequencies
integer, parameter :: max_nvar = 16
integer, parameter :: max_nstream = 32          ! (-) max num quadrature angles
integer, parameter :: max_nphase = 5            ! (-) max number number hydrometeor phases
integer, parameter :: nphase = 5                ! (-) number of hydrometeor phases
integer, parameter :: nvar = 2 + 2 * nphase     ! (-) number of variables for radiation Jacobian
integer, parameter :: nvar_prof = 9             ! (-) number of variables in atmospheric profile
integer, parameter :: npol = 2                  ! (-) number of polarizations

! grid related parameters
integer nlev         ! (-) number of vertical levels
integer nstream      ! (-) total number of stream angles, including upward and downward propagating
integer nstream_surf ! (-) total number of stream angles for surface, upward only (=nstream/2)
integer obs_lev
integer obs_lev1
integer ng
integer nlr
integer nlr1
integer m1
integer i0
integer nang
integer nangover2

real(8) obs_theta    ! (deg) zenith angle to observation
real(8) obs_height   ! (km) height of observation
real(8) ipol         ! (-) polarization number 1=vertical 0=horizontal
integer num_freqs    ! (-) number of freq points for passband quadrature
real(8) instrspec(5) ! (Ghz) instrument specification

integer nsub_freq ! Number of frequency points in a sideband for passband quadrature
integer nHG

integer :: num_quad_angles ! = max_nstream
real(8), dimension(max_num_int_freqs) :: passband_freq
real(8), dimension(max_nstream) :: quad_angle ! (-) quadrature angles
real(8), dimension(max_nlev) :: g_asymmetry
LOGICAL a0_is_constant(max_nlev,max_nphase)

type gas_type
    real(8) absn2      ! nitrogen absorption in nepers/km
    real(8) absh2o     ! water absorption in nepers/km
    real(8) o2abs      ! oxygen absorption in nepers/km
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

! geophysical Jacobian
real(8), dimension(max_nlev) :: dKab_dT, dKab_dp, dKab_dq
real(8), dimension(max_nlev,max_nphase) :: dKsc_dT, dg_dT, dKab_dw, dKsc_dw, dg_dw

! total Jacobian
real(8), dimension(max_nlev,max_nstream,2) :: dTb_dT, dTb_dp, dTb_dq
real(8), dimension(max_nlev,max_nstream,max_nphase,2) :: dTb_dw

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

type profile_type
    ! Gaseous state variables
    real(8) pressure    ! (mb) pressure
    real(8) temperature ! (K) temperature
    real(8) vapor_dens  ! (g/m3) water vapor density

    ! Hydrometeor distribution state variable parameters
    real(8) clw_p       ! cloud liquid water particle size distribution exponent
    real(8) clw_q       ! cloud liquid water particle size variance parameter
    real(8) clw_k0      ! cloud liquid water particle size normalization constant
    real(8) clw_a0      ! cloud liquid water particle size parameter in mm
    real(8) rain_p      ! cloud rain (precip-water) particle size dist. exponent
    real(8) rain_q      ! cloud rain particle size variance parameter
    real(8) rain_k0     ! cloud rain particle size normalization constant
    real(8) rain_a0     ! cloud rain  particle size parameter in mm
    real(8) ice_p       ! cloud frozen particle size distribution exponent
    real(8) ice_q       ! cloud frozen particle size variance parameter
    real(8) ice_k0      ! cloud frozen particle size normalization constant
    real(8) ice_a0      ! cloud frozen particle size parameter in mm
    real(8) snow_p      ! cloud snow particle size distribution exponent
    real(8) snow_q      ! cloud snow particle size variance parameter
    real(8) snow_k0     ! cloud snow particle size normalization constant
    real(8) snow_a0     ! cloud snow particle size parameter in mm
    real(8) grpl_p      ! cloud graupel particle size distribution exponent
    real(8) grpl_q      ! cloud graupel particle size variance parameter
    real(8) grpl_k0     ! cloud graupel particle size normalization constant
    real(8) grpl_a0     ! cloud graupel particle size parameter in mm

    real(8) hgt         ! (km) height to middle of layer wrt to surface
    real(8) hgt_top     ! (km) height to top of layer wrt to surface
    real(8) hgt_bot     ! (km) height to bottom of layer wrt to surface
    real(8) hgt_del     ! (m) layer thickness (not the units are meters not km)

    real(8) clw_dens    ! (g/m3) cloud liquid water density
    real(8) rain_dens   ! (g/m3) cloud rain density
    real(8) ice_dens    ! (g/m3) cloud ice density
    real(8) snow_dens   ! (g/m3) cloud snow density
    real(8) grpl_dens   ! (g/m3) cloud graupel density

    ! derivative of k0 and a0 with respect to water density
                        ! at each level for each phase of water
                        ! k0 and a0 depend only on water density
    real(8) dclw_k0_dw  ! k0 water density derivative
    real(8) dclw_a0_dw  ! a0 water density derivative
    real(8) drain_k0_dw ! k0 water density derivative
    real(8) drain_a0_dw ! a0 water density derivative
    real(8) dice_k0_dw  ! k0 water density derivative
    real(8) dice_a0_dw  ! a0 water density derivative
    real(8) dsnow_k0_dw ! k0 water density derivative
    real(8) dsnow_a0_dw ! a0 water density derivative
    real(8) dgrpl_k0_dw ! k0 water density derivative
    real(8) dgrpl_a0_dw ! a0 water density derivative

    ! Derived water vapor quantities
    real(8) h2o_v_sat   ! water vapor saturation pressure in mb
    real(8) rel_hum     ! relative humidity 0-100%

    ! Derived radiative transfer quantities : frequency dependent
    real(8) abs_o2      ! oxygen and nitrogen absorption in nepers/km
    real(8) abs_h2o     ! water vapor absorption in nepers/km
    real(8) abs_cloud   ! cloud particle absorption in nepers/km
    real(8) scat_cloud  ! cloud particle scattering in nepers/km
    real(8) asymmetry   ! cloud particle asymmetry factor (unitless)
    real(8) ext_tot     ! total atmospheric extinction in nepers/km
    real(8) albedo      ! single scattering albedo
    real(8) bb_spec_int ! brightness temperature or specific 
                        ! black body intensity for either polarization at level in W/m**2/Hz/ster
end type profile_type
type(profile_type) atm(max_nlev)

type reduced_profile_type
    real(8) clw_tot     ! (g/m2) total column cloud liquid water
    real(8) rain_tot    ! (g/m2) total column rain
    real(8) ice_tot     ! (g/m2) total column ice
    real(8) snow_tot    ! (g/m2) total column snow
    real(8) grpl_tot    ! (g/m2) total column graupel

    real(8) clw_h_ave   ! (km) average height cloud liquid water
    real(8) rain_h_ave  ! (km) average height rain
    real(8) ice_h_ave   ! (km) average height ice
    real(8) snow_h_ave  ! (km) average height snow
    real(8) grpl_h_ave  ! (km) average height graupel

    real(8) clw_h_std   ! (km) height standard deviation cloud liquid water
    real(8) rain_h_std  ! (km) height standard deviation rain
    real(8) ice_h_std   ! (km) height standard deviation ice
    real(8) snow_h_std  ! (km) height standard deviation snow
    real(8) grpl_h_std  ! (km) height standard deviation graupel
end type reduced_profile_type
type(reduced_profile_type) atm_reduced

! surface characteristics branch
type surf_inp_type   ! Assumes specular reflecting surface
     integer(4) land_flag          ! (-) land flag 2=land, 1=lake, 0=ocean
     character*20 ocean_mod        ! (-) ocean model: 'Fresnel','Wilheit',or 'ITRA'
     character*20 diel_mod         ! (-) water dielectric model: 'Salt','NaCl',or 'Fresh'
     real(8) surf_temp             ! (K) surface skin temperature
     real(8) theta(max_num_angles) ! (deg) angle of incidence or reflection wrt surface normal
     real(8) vr(max_num_angles)    ! (-) vertical surface reflectivity
     real(8) hr(max_num_angles)    ! (-) horizontal surface reflectivity
     real(8) salinity              ! (-) water salinity
     double complex dielectric     ! (-) dielectric constant
     real(8) windspeed             ! (m/s) surface wind speed
     integer nstream_surf          ! (-) number of surface stream angles
end type surf_inp_type
type (surf_inp_type) :: surf_inp   ! Surface information data

!------------------------------------------------------------------------------
! Instrument specifications branch
!------------------------------------------------------------------------------
! Frequency and noise specification for a single channel -flat passband model
integer nchannel  ! (-) number of channels
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
type (chan_spec_type), dimension(max_nchannel) :: instr_spec
type (chan_spec_type) :: channel

real(8), dimension(max_nstream) :: teta
real(8), dimension(max_nstream) :: cs
real(8), dimension(max_num_angles) :: surf_reflecv
real(8), dimension(max_num_angles) :: surf_reflech
real(8), dimension(max_num_angles) :: surf_reflec
real(8), dimension(max_num_angles) :: surfinp_reflecv
real(8), dimension(max_num_angles) :: surfinp_reflech
real(8), dimension(max_num_angles) :: reflec_angle_array
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
real(8), dimension(max_nlev) :: al_gas1
real(8), dimension(max_nlev) :: dtemperature1
real(8), dimension(max_nlev) :: dabs_cloud1
real(8), dimension(max_nlev) :: dscat_cloud1
real(8), dimension(max_nlev) :: dal_gas1
real(8), dimension(max_nlev) :: h1
real(8), dimension(max_nlev) :: abs_total
real(8), dimension(max_nlev) :: abs_total1

real(8), dimension(101) :: HGg
real(8), dimension(101,max_nstream,max_nstream) :: HGph
real(8), dimension(101,max_nstream,max_nstream) :: dHGph

real(8), dimension(max_nlev+1) :: altitude
real(8), dimension(max_nlev+1) :: altitude1

real(8), dimension(0:max_nlev,max_nstream) :: tb_pl
real(8), dimension(0:max_nlev,max_nstream) :: tb_mn
real(8), dimension(0:max_nlev,max_nstream,max_nvar):: dtb_pl
real(8), dimension(0:max_nlev,max_nstream,max_nvar):: dtb_mn

real(8), dimension(max_nlev,max_nstream,max_nstream) :: phase11
real(8), dimension(max_nlev,max_nstream,max_nstream) :: dphase11
real(8), dimension(max_nlev,max_nstream,max_nstream) :: phaseff
real(8), dimension(max_nlev,max_nstream,max_nstream) :: phasefb
real(8), dimension(max_nlev,max_nstream,max_nstream) :: dphaseff
real(8), dimension(max_nlev,max_nstream,max_nstream) :: dphasefb
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
real(8), dimension(max_nstream) :: cris_quad_wghts

end module dotlrt_variables
