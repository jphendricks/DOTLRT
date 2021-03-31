module variables
!    use type_kinds, ONLY: fp_kind
    implicit none
  integer instr_unit
  integer atm_unit
  integer rt_unit(4,32)
  integer Ksa_unit
  integer geo_unit
  integer tot_unit(32)
  integer jrec_Ksa
  integer jrec_tot(32)
  integer jrec_geo
  integer jrec_rt(44,32)
  integer dims(4)
    real(8) time_k, time_rt, time_io, times(2)

    integer, parameter :: max_num_levels = 128 !Max # of input press levels to define atmos. state
    integer nlev, ngrid

    integer, parameter :: max_num_angles = 16 ! 16  !Max # of angles to define surface reflectivity

    real(8), parameter  :: hplank = 6.6252d-34    !{Plank's constant in J-s
    real(8), parameter  :: kboltz = 1.380662d-23     !{Boltzmann's constant in J/K
    real(8), parameter  :: pi= 3.14159265358979d0 ! not encircling any gravitational singularities
    real(8), parameter  :: rgas = 8.31441d0        !Ideal gas constant J/mole-K
    real(8), parameter  :: quad_weight_error =1.0d-5  !Allowable unitarity error for all quadrature weights
    real(8), parameter  :: c = 2.99793d8            ! speed of light (m/s)
    real(8), parameter  :: std_temperature = 273.15d0 !Standard temperature in K
    real(8), parameter  :: std_pressure = 1013.25d0  !Standard pressure in mb
    real(8), parameter  :: max_freq = 1000.0d0      !Max frequency (GHz) for absorption computations
    real(8), parameter  :: temp_incr = 0.01d0   !absorption coeff, temperature inc weight functions
    real(8), parameter  :: wv_incr = 0.01d0       !water vapor inc for calc deriv in absorption coeff
    real(8), parameter  :: btcb = 2.73d0           !brightness temperature of cosmic background
    real(8), parameter  :: rad_earth = 6.3568d3  !Radius of the earth (km)
    real(8), parameter  :: raythresh=0.1d0    !Rayleigh scattering and extinction efficiency
    real(8), parameter  :: t_cb = 2.730d0 ! cosmic background temperature

    integer, parameter :: max_num_channels = 64 !Maximum number of instrument channels
    integer, parameter :: max_num_int_freqs = 45 !Max # of spectral integration frequencies
    integer, parameter :: max_nvar = 16
    integer, parameter :: max_num_quad_angles = 32 ! 16 !

    integer, parameter :: max_number_h2o_phases = 5 ! number of phases of water
    integer number_h2o_phases, jfreq(2)
    integer obs_lev, obs_lev1, ng, nlr, nlr1, m1, nvar, i0, nang, nangover2, int_method
    real(8) d_albedo

    real(8) :: inp_theta, inp_height, inp_pol
    integer :: num_surf_angles, nfreq, num_freqs
    ! Number of frequency points in a sideband for passband quadrature
    integer        :: num_sb_freqs, nHG
    !Instrument noise covariance matrix
    real(8), dimension(max_num_channels,max_num_channels) :: instr_noise_cov_matrix 

    integer   :: num_quad_angles ! = max_num_quad_angles

    real(8), dimension(max_num_int_freqs) :: passband_freq
    real(8), dimension(max_num_quad_angles) :: quad_angle_array 

    real(8), dimension(max_num_levels) :: g_asymmetry

    LOGICAL a0_is_constant(max_num_levels,max_number_h2o_phases)

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
    type (gas_type), dimension(max_num_levels) :: gas_prof 

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
        real(8) dcloudg_da0     ! a0 derivative                                        ! not temperature or microwave frequency
    end type hydro_type
    type (hydro_type), dimension(max_num_levels,max_number_h2o_phases) :: hydro_prof 

    ! geophysical Jacobian
    real(8), dimension(max_num_levels) :: dKab_dT, dKab_dp, dKab_dq
    real(8), dimension(max_num_levels,max_number_h2o_phases) :: dKsc_dT, dg_dT, dKab_dw, dKsc_dw, dg_dw

    ! total Jacobian
    real(8), dimension(max_num_levels,max_num_quad_angles,2) :: dTb_dT, dTb_dp, dTb_dq
    real(8), dimension(max_num_levels,max_num_quad_angles,max_number_h2o_phases,2) :: dTb_dw

    ! radiative transfer Jacobian
    real(8), dimension(max_num_levels,max_num_quad_angles,2) :: dTbdTr,      & ! temperature
                                           dTbdKa,      & ! absorption
                                           dTbdKsliq,   & ! scatter liquid
                                           dTbdgliq,    & ! asymmetry liquid
                                           dTbdKsrn,    & ! scatter rain
                                           dTbdgrn,     & ! asymmetry rain
                                           dTbdKsice,   & ! scatter ice
                                           dTbdgice,    & ! asymmetry ice
                                           dTbdKssnow,  & ! scatter snow
                                           dTbdgsnow,   & ! asymmetry snow
                                           dTbdKsgrpl,  & ! scatter graupel
                                           dTbdggrpl      ! asymmetry graupel

    type profile_type
               ! Gaseous state variables
                   real(8) pressure      !pressure of level in mb
                   real(8) temperature   !temperature of level in K
                   real(8) vapor_density !water vapor density of level in g/m**3
               ! Hydrometeor distribution state variable parameters
                   real(8) cloud_liq_p   !cloud liquid particle size distribution exponent
                   real(8) cloud_liq_q      !cloud liquid particle size variance parameter
                   real(8) cloud_liq_k0     !cloud liquid particle size normalization constant
                   real(8) cloud_liq_a0     !cloud liquid particle size parameter in mm
                   real(8) cloud_rn_p       !cloud rain (precip-water) particle size dist. exponent
                   real(8) cloud_rn_q       !cloud rain particle size variance parameter
                   real(8) cloud_rn_k0      !cloud rain particle size normalization constant
                   real(8) cloud_rn_a0      !cloud rain  particle size parameter in mm
                   real(8) cloud_ice_p      !cloud frozen particle size distribution exponent
                   real(8) cloud_ice_q      !cloud frozen particle size variance parameter
                   real(8) cloud_ice_k0     !cloud frozen particle size normalization constant
                   real(8) cloud_ice_a0     !cloud frozen particle size parameter in mm
                   real(8) cloud_snow_p     !cloud snow particle size distribution exponent
                   real(8) cloud_snow_q     !cloud snow particle size variance parameter
                   real(8) cloud_snow_k0    !cloud snow particle size normalization constant
                   real(8) cloud_snow_a0    !cloud snow particle size parameter in mm
                   real(8) cloud_grpl_p     !cloud graupel particle size distribution exponent
                   real(8) cloud_grpl_q     !cloud graupel particle size variance parameter
                   real(8) cloud_grpl_k0    !cloud graupel particle size normalization constant
                   real(8) cloud_grpl_a0    !cloud graupel particle size parameter in mm
                                           ! Derived height index
                   real(8) height           !height of level in km wrt to surface
                                           ! Derived hydrometeor distribution statistics
                   real(8) cloud_liq_dens   !cloud liquid water density in g/m**3
                   real(8) cloud_liq_rad    !cloud liquid particle average mode radius in mm
                   real(8) cloud_liq_std    !cloud liquid particle size standard deviation
                   real(8) cloud_liq_f      !cloud liquid fractional volume
                   real(8) cloud_rn_dens    !cloud rain  water density in g/m**3
                   real(8) cloud_rn_rad     !cloud rain particle average mode radius in mm
                   real(8) cloud_rn_std     !cloud rain particle size standard deviation
                   real(8) cloud_rn_f       !cloud rain fractional volume
                   real(8) cloud_ice_dens   !cloud frozen water density in g/m**3
                   real(8) cloud_ice_rad    !cloud frozen particle average mode radius in mm
                   real(8) cloud_ice_std    !cloud frozen particle size standard deviation
                   real(8) cloud_ice_f      !cloud frozen fractional volume
                   real(8) cloud_snow_dens  !cloud snow water density in g/m**3
                   real(8) cloud_snow_rad   !cloud snow particle average mode radius in mm
                   real(8) cloud_snow_std   !cloud snow particle size standard deviation
                   real(8) cloud_snow_f     !cloud snow fractional volume
                   real(8) cloud_grpl_dens  !cloud graupel water density in g/m**3
                   real(8) cloud_grpl_rad   !cloud graupel particle average mode radius in mm
                   real(8) cloud_grpl_std   !cloud graupel particle size standard deviation
                   real(8) cloud_grpl_f     !cloud graupel fractional volume
                                           ! derivative of k0 and a0 with respect to water density
                                           ! at each level for each phase of water
                                           ! k0 and a0 depend only on water density
                   real(8) dcloud_liq_k0_dw  ! k0 water density derivative
                   real(8) dcloud_liq_a0_dw  ! a0 water density derivative
                   real(8) dcloud_rn_k0_dw   ! k0 water density derivative
                   real(8) dcloud_rn_a0_dw   ! a0 water density derivative
                   real(8) dcloud_ice_k0_dw  ! k0 water density derivative
                   real(8) dcloud_ice_a0_dw  ! a0 water density derivative
                   real(8) dcloud_snow_k0_dw ! k0 water density derivative
                   real(8) dcloud_snow_a0_dw ! a0 water density derivative
                   real(8) dcloud_grpl_k0_dw ! k0 water density derivative
                   real(8) dcloud_grpl_a0_dw ! a0 water density derivative
                                           ! Derived water vapor quantities
                   real(8) h2o_v_sat        !water vapor saturation pressure in mb
                   real(8) rel_hum          !relative humidity 0-100%
                                           ! Derived radiative transfer quantities : frequency dependent
                   real(8) abs_o2           !oxygen and nitrogen absorption in nepers/km
                   real(8) abs_h2o          !water vapor absorption in nepers/km
                   real(8) abs_cloud        !cloud particle absorption in nepers/km
                   real(8) scat_cloud       !cloud particle scattering in nepers/km
                   real(8) asymmetry        !cloud particle asymmetry factor (unitless)
                   real(8) ext_tot          !total atmospheric extinction in nepers/km
                   real(8) albedo           !single scattering albedo
                   real(8) bb_spec_int      !brightness temperature or specific 
                        !black body intensity for either polarization at level in W/m**2/Hz/ster
    end type profile_type

    type atm_inp_type
        ! type (profile_type), dimension(max_num_levels) :: prof                      
        type(profile_type) prof(max_num_levels)
        integer   month
        integer   latitude     
        real(8)    h2o_v_scale_hgt 
        real(8)    surf_rel_hum    
        integer   num_levels      
        character*80 text          
        logical   inf             
        logical   fiveph
    end type atm_inp_type 
    type(atm_inp_type) :: atm_inp

    type surf_inp_type   !Assumes specular reflecting surface
         character*1 surf_type   !'L' for land, 'O' for ocean
         character*1 surf_subtype !ocean model type: 'F','W','I'
         real(8) freq             !frequency for which is surface model valid
         real(8) param1           !surface parameters, now used for: reflectivity
         real(8) param2           !salinity
         real(8) param3           !wind speed
         real(8) surf_temp        !surface kinetic temperature in K
         ! real(8), dimension(max_num_angles) :: theta !angle of incidence or reflection wrt surface normal
         ! real(8), dimension(max_num_angles) ::  vr  !vertical surface reflectivity
         ! real(8), dimension(max_num_angles) ::  hr !horizontal surface reflectivity
         real(8) theta(max_num_angles) !angle of incidence or reflection wrt surface normal
         real(8) vr(max_num_angles) !vertical surface reflectivity
         real(8) hr(max_num_angles) !horizontal surface reflectivity
         integer num_angles      !number of surface angles
         character text          !string descriptor
         logical inf           
    end type surf_inp_type
    type ( surf_inp_type)   :: surf_inp           !Surface information data

    ! Frequency and noise specification for a single channel -flat passband model
    type chan_spec_type 
         real(8)  :: lo_freq        !LO frequency     (GHz)
         real(8)  :: if1_freq       !IF1 frequency    (GHz)
         real(8)  :: if2_freq       !IF2 frequency    (GHz)
         real(8)  :: bandwidth      !filter bandwidth (GHz)
         real(8)  :: dtrms          !observation noise in K
         integer :: chan_desig     !channel designation
    end type chan_spec_type
    type(chan_spec_type)   :: chan_spec
    type(chan_spec_type)   :: chan

    !Instrument channel frequency and noise specifications - flat passband model
    type instr_spec_type
         type (chan_spec_type), DIMENSION(max_num_channels) :: chan
         integer                                            :: num_channels
         CHARACTER(80)                                       :: text
         LOGICAL                                            :: inf
    end type instr_spec_type
    type (instr_spec_type)  :: instr_spec

    real(8), dimension(max_num_quad_angles) :: teta, cs
    real(8), dimension(max_num_angles) :: surf_reflecv, surf_reflech, &
                                                   surf_reflec, surfinp_reflecv, &
                                                   surfinp_reflech, reflec_angle_array
    real(8), DIMENSION(max_num_levels) :: temperature, abs_O2, abs_H2O,     &
                                                  abs_cloud , scat_cloud, temperature1, &
                                                  abs_O2_1, abs_H2O_1, abs_cloud1,      &
                                                  scat_cloud1 , al_gas1, dtemperature1, &
                                                  dabs_cloud1, dscat_cloud1, dal_gas1, h1, &
                                                  abs_total, abs_total1

    real(8), DIMENSION(101) :: HGg
    real(8), dimension(101,max_num_quad_angles,max_num_quad_angles) :: HGph, dHGph

    real(8), DIMENSION(max_num_levels+1) :: altitude, altitude1

    real(8), DIMENSION(0:max_num_levels,max_num_quad_angles) :: tb_pl, tb_mn
    real(8), DIMENSION(0:max_num_levels,max_num_quad_angles,max_nvar):: dtb_pl, dtb_mn

    real(8), dimension(max_num_levels,max_num_quad_angles,max_num_quad_angles) :: phase11, dphase11, & 
                                                              phaseff, phasefb, dphaseff, dphasefb, &
                                                              phaseff1, phasefb1, dphaseff1, dphasefb1

    real(8), dimension(max_num_levels,max_num_quad_angles,max_num_quad_angles) :: phase11_sc, phaseff_sc, phasefb_sc, phaseff1_sc, phasefb1_sc
    real(8), dimension(max_num_levels,max_num_quad_angles,max_num_quad_angles,max_number_h2o_phases) :: dphaseff_g, dphasefb_g, &
                                                                                                                 dphaseff_sc, dphasefb_sc, &
                                                                                                                 d_phase11_g, d_phase11_sc, &
                                                                                                                 dphaseff1_g, dphasefb1_g, &
                                                                                                                 dphaseff1_sc, dphasefb1_sc
    real(8), dimension(max_num_quad_angles) :: cris_quad_wghts


      INTEGER, DIMENSION(4) :: dim_type 
      REAL(8), DIMENSION(59) :: data_array_type

      logical isfileopen

      TYPE prof_data_type 
         ! real(8), dimension(max_num_levels) :: u
         ! real(8), dimension(max_num_levels) :: v
         ! real(8), dimension(max_num_levels) :: pp
         ! real(8), dimension(max_num_levels) :: temperature
         ! real(8), dimension(max_num_levels) :: clw
         ! real(8), dimension(max_num_levels) :: rnw
         ! real(8), dimension(max_num_levels) :: ice
         ! real(8), dimension(max_num_levels) :: snow
         ! real(8), dimension(max_num_levels) :: grauple
         ! real(8), dimension(max_num_levels) :: altitude
         ! real(8), dimension(max_num_levels) :: q
         ! real(8), dimension(max_num_levels) :: sigmaH
         ! REAL(8)   :: terrain
         ! REAL(8)   :: ground_temp
         ! REAL(8)   :: pstarcrs
         ! REAL(8)   :: lat
         ! REAL(8)   :: lng
         ! REAL(8)   :: TopPress
         ! REAL(8)   :: Sealevel_press
         ! REAL(8)   :: Sealevel_temp
         ! REAL(8)   :: Lapse_rate
         ! INTEGER   :: X, Y, Z, T 
         real(8) u (max_num_levels)
         real(8) v (max_num_levels)
         real(8) pp (max_num_levels)
         real(8) temperature (max_num_levels)
         real(8) clw (max_num_levels)
         real(8) rnw (max_num_levels)
         real(8) ice (max_num_levels)
         real(8) snow (max_num_levels)
         real(8) grauple (max_num_levels)
         real(8) altitude (max_num_levels)
         real(8) q (max_num_levels)
         real(8) sigmaH (max_num_levels)
         REAL(8) terrain
         REAL(8) ground_temp
         REAL(8) pstarcrs
         REAL(8) lat
         REAL(8) lng
         REAL(8) TopPress
         REAL(8) Sealevel_press
         REAL(8) Sealevel_temp
         REAL(8) Lapse_rate
         INTEGER X, Y, Z, T 
      END TYPE prof_data_type

end module variables
