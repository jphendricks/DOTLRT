module sibtype

!----------------------------------------------------------------------
!
!   New SiB module: SiB returns to a single-point model, for 
!   reasons of compatability with BUGS and MPI offline code (SiBDRV)
!   as well as to adhere to the CCM4 coding standard
!
! Modifications
!  created Ian Baker, 06 Feb 2002
!  Kevin Schaefer changed tdew1/2 to sh1/2 (8/17/04)
!  Kevin Schaefer moved respfactor variables to param branch var tree (5/6/05)
!  Kevin Schaefer added auto, hetero, and total respiration (5/6/05)
!  Kevin Schaefer added CASA related variables (6/5/05)
!  Kevin Schaefer deleted dayflag (8/2/06)
!  Kevin Schaefer started new drvr variable branch (1/27/10)
!  Kevin Schaefer removed unused variables: zb,psb,tcc,vdsrad,bps(2) (1/25/10)
!  Kevin schaefer deleted all phosib potential rates (wegs, wags, whs, wci,
!    wsfht, wsflt, wsfws, ansqr, antemp, omepot, assimnp, assimci, assimpot) (5/28/13)
!----------------------------------------------------------------------

use kinds
use sib_const_module

implicit none

!jlc...These need to be public
public  sib_t
public  sib_local_vars

!jlc...Not sure what these need to be
public param_vars
public prognostic_vars
public diagnostic_vars
public sib_status

!------------------------------------------------------------------
!                   BOUNDARY CONDITION VARIABLES
!------------------------------------------------------------------
type param_vars

    !...boundary conditions--TIME INVARIANT
    real(kind=dbl_kind) :: biome        ! (-) biome type (see refs for description)
    real(kind=dbl_kind) :: lat          ! (deg) latitude
    real(kind=dbl_kind) :: lon          ! (deg) longitude
    real(kind=dbl_kind) :: chil         ! (-) leaf angle distribution factor
    real(kind=dbl_kind) :: lai_max      ! (m2/m2) max leaf area index per biome
    real(kind=dbl_kind) :: phc          ! (m) 1/2 crit leaf water pot limit
    real(kind=dbl_kind) :: z1           ! (m) canopy bottom
    real(kind=dbl_kind) :: z2           ! (m) canopy top
    real(kind=dbl_kind) :: poros(nsoil) ! (-) soil porosity (zero to one)
    real(kind=dbl_kind) :: poros_min    ! (-) mineral soil porosity (zero to one)
    real(kind=dbl_kind) :: poros_om     ! (-) organic soil porosity (zero to one)
    real(kind=dbl_kind) :: satco(nsoil) ! (m/s) hydraulic conductivity at saturation
    real(kind=dbl_kind) :: satco_min    ! (m/s) hydraulic conductivity at saturation of mineral soil
    real(kind=dbl_kind) :: satco_om     ! (m/s) hydraulic conductivity at saturation of organic soil
    real(kind=dbl_kind) :: bee(nsoil)   ! (-) Clapp & Hornberber 'b' exponent
    real(kind=dbl_kind) :: bee_min      ! (-) Clapp & Hornberber 'b' exponent for mineral soil
    real(kind=dbl_kind) :: bee_om       ! (-) Clapp & Hornberber 'b' exponent for organic soil
    real(kind=dbl_kind) :: phsat(nsoil) ! (m) soil tension at saturation
    real(kind=dbl_kind) :: phsat_min    ! (m) soil tension at saturation of mineral soil
    real(kind=dbl_kind) :: phsat_om     ! (m) soil tension at saturation of organic soil
    real(kind=dbl_kind) :: slope        ! (-) cosine of mean slope
    real(kind=dbl_kind) :: vcover       ! (-) fraction of vegetation cover (0-1)
    real(kind=dbl_kind) :: d_permc_min  ! (m) depth to top of permafrost carbon layer   
    real(kind=dbl_kind) :: d_permc_max  ! (m) depth to bottom of permafrost carbon layer   
    real(kind=dbl_kind) :: permc_den    ! (mole/m3) permafrost carbon density   
    real(kind=dbl_kind) :: permc_tot    ! (mole/m2) total permafrost carbon in soil column   
    real(kind=dbl_kind) :: vmax0(5)     ! (mol/m^2/sec) Rubisco vel of sun leaf
    !   (1) C3 photosynthesis
    !   (2) C4 photosynthesis
    !   (3) open
    !   (4) open
    !   (5) open
    real(kind=dbl_kind) :: trop(5)   ! temp coeff in GS-A model (K)
    !   (1) C3 photosynthesis
    !   (2) C4 photosynthesis
    !   (3) open
    !   (4) open
    !   (5) open
    real(kind=dbl_kind) :: trda(5)   ! temp coeff in GS-A model (K^-1)
    !   (1) C3 photosynthesis
    !   (2) C4 photosynthesis
    !   (3) open
    !   (4) open
    !   (5) open
    real(kind=dbl_kind) :: trdm(5)   ! temp coeff in GS-A model (K)
    !   (1) C3 photosynthesis
    !   (2) C4 photosynthesis
    !   (3) open
    !   (4) open
    !   (5) open
    real(kind=dbl_kind) :: respcp(5) ! respiration fraction of vmax0 (-)
    !   (1) C3 photosynthesis
    !   (2) C4 photosynthesis
    !   (3) open
    !   (4) open
    !   (5) open
    real(kind=dbl_kind) :: slti(5)   ! slope of lo-temp inhibition (K^-1)
    !   function (K^-1)
    !   (1) C3 photosynthesis
    !   (2) C4 photosynthesis
    !   (3) open
    !   (4) open
    !   (5) open
    real(kind=dbl_kind) :: shti(5)   ! slope of hi-temp inhibition (K^-1)
    !   function (K-1)
    !   (1) C3 photosynthesis
    !   (2) C4 photosynthesis
    !   (3) open
    !   (4) open
    !   (5) open
    real(kind=dbl_kind) :: hltii(5)  ! 1/2 point of lo-temp inhibition (K)
    !   function
    !   (1) C3 photosynthesis
    !   (2) C4 photosynthesis
    !   (3) open
    !   (4) open
    !   (5) open
    real(kind=dbl_kind) :: hhti(5)   ! 1/2 point of hi-temp inhibition (K)
    !   function
    !   (1) C3 photosynthesis
    !   (2) C4 photosynthesis
    !   (3) open
    !   (4) open
    !   (5) open
    real(kind=dbl_kind) :: effcon(5) ! quantum efficiency (mol/mol)
    !   (1) C3 photosynthesis
    !   (2) C4 photosynthesis
    !   (3) open
    !   (4) open
    !   (5) open
    real(kind=dbl_kind) :: binter(5) ! conductance-photosynthesis 
    !   intercept (mol m^-2 sec^-1)
    !   (1) C3 photosynthesis
    !   (2) C4 photosynthesis
    !   (3) open
    !   (4) open
    !   (5) open
    real(kind=dbl_kind) :: gradm(5)  ! conductance-photosynthesis slope 
    !   parameter (-)
    !   (1) C3 photosynthesis
    !   (2) C4 photosynthesis
    !   (3) open
    !   (4) open
    !   (5) open
    real(kind=dbl_kind) :: atheta(5) ! WC WE coupling parameter (-)
    !   (1) C3 photosynthesis
    !   (2) C4 photosynthesis
    !   (3) open
    !   (4) open
    !   (5) open
    real(kind=dbl_kind) :: btheta(5) ! WC&WE, WS coupling parameter (-)
    !   (1) C3 photosynthesis
    !   (2) C4 photosynthesis
    !   (3) open
    !   (4) open
    !   (5) open
    real(kind=dbl_kind) :: tfrost    ! temp below which frost reduces photosynthesis (K) 
    real(kind=dbl_kind) :: tcmin_rr  ! min canopy temp recovery rate (K/day) 
    real(kind=dbl_kind) :: soref(2)  ! soil reflectance (-)
    !   (1) shortwave
    !   (2) longwave
    real(kind=dbl_kind) :: sftf      ! (1/K) slope leaf/root mortality frost function
    real(kind=dbl_kind) :: hftf      ! (K) half point leaf/root mortality frost function
    real(kind=dbl_kind) :: tcftf     ! (K) critical underflow temperature for leaf/root mortality frost function
    real(kind=dbl_kind) :: sfrzi     ! (1/K) slope respiration freezing scaling factor
    real(kind=dbl_kind) :: hfrzi     ! (K) half point respiration freezing scaling factor
    real(kind=dbl_kind) :: tcfrzi    ! (K) critical underflow temperature for respiration freezing scaling factor
    real(kind=dbl_kind) :: zm        ! respiration parameter - exponent
    real(kind=dbl_kind) :: wopt      ! respiration parameter - optimum soil moisture (-)
    real(kind=dbl_kind) :: woptzm    ! wopt to the zm exponent
    real(kind=dbl_kind) :: wsat      ! respiration parameter (-)
    real(kind=dbl_kind) :: sandfrac  ! soil sand fraction
    real(kind=dbl_kind) :: clayfrac  ! soil clay fraction
    real(kind=dbl_kind) :: beta_clay  ! (-) exponent for liquid water fraction in frozen soil for clay
    real(kind=dbl_kind) :: beta_silt  ! (-) exponent for liquid water fraction in frozen soil for silt
    real(kind=dbl_kind) :: beta_sand  ! (-) exponent for liquid water fraction in frozen soil for sand
    real(kind=dbl_kind) :: beta_org   ! (-) exponent for liquid water fraction in frozen soil for organic matter
    real(kind=dbl_kind) :: tstar_clay ! (K) freezing point depression for clay
    real(kind=dbl_kind) :: tstar_silt ! (K) freezing point depression for silt
    real(kind=dbl_kind) :: tstar_sand ! (K) freezing point depression for sand
    real(kind=dbl_kind) :: tstar_org  ! (K) freezing point depression for organic matter
    
    real(kind=dbl_kind) :: physfrac(5)
    ! physiology fraction-5 elements must add up to 1.0
    !  (1) C3 vegetation
    !  (2) C4 vegetation
    !  (3) open
    !  (4) open
    !  (5) open
    real(kind=dbl_kind) :: physfrac1(5)
    real(kind=dbl_kind) :: physfrac2(5)
    real(kind=dbl_kind) :: physfrac3(5)
    ! physiology fraction-5 elements must add up to 1.0
    !  (1) C3 vegetation
    !  (2) C4 vegetation
    !  (3) open
    !  (4) open
    !  (5) open

    integer(kind=int_kind) :: phystype(5)
    ! physiology type-will be either 3 or 4
    !  (1) C3 vegetation - always 3
    !  (2) C4 vegetation - always 4
    !  (3) open
    !  (4) open
    !  (5) open
    real(kind=dbl_kind) :: kroot        ! (1/m) root density extinction coeficient
    real(kind=dbl_kind) :: rootf(nsoil) ! (-) root fraction per soil layer
    real(kind=dbl_kind) :: rootr(nsoil) ! (-) adjusted root fraction per soil layer
    real(kind=dbl_kind) :: fcarb(nsoil) ! (-) total soil carbon fraction per soil layer
    real(kind=dbl_kind) :: fcarb_lay(npoolmax,nsoil) ! (-) pool fraction per soil layer
    real(kind=dbl_kind) :: rootd        ! (m) maximum rooting depth
    integer(kind=int_kind) :: nroot     ! (-) lowest soil layer with roots
    integer(kind=int_kind) :: ncarb(npoolmax)     ! (-) lowest soil layer with carbon
    real(kind=dbl_kind) :: tran(2,2)    ! leaf transmittance (-)
    !  (1,1) - shortwave, green plants
    !  (1,2) - longwave, green plants
    !  (2,1) - shortwave, brown plants
    !  (2,2) - longwave, brown plants
    real(kind=dbl_kind) :: ref(2,2)  ! leaf reflectance (-)
    !  (1,1) - shortwave, green plants
    !  (1,2) - longwave, green plants
    !  (2,1) - shortwave, brown plants
    !  (2,2) - longwave, brown plants
    real(kind=dbl_kind) :: respfactor(nsoil) 
    ! factor for balancing soil respiration
    !  with annual assimilation (-)

    real(kind=dbl_kind) :: het_respfac(nsoil) 
    ! hetrotrophic respiration factor for balancing soil heterotrophic
    !  respiration with annual canopy net assimilation (-)

    real(kind=dbl_kind) :: auto_respfac 
    ! autotrophic respiration factor for balancing growth/maintenance
    !  respiration with annual canopy net assimilation (-)

    real(kind=dbl_kind) :: tot_an(13)       ! total canopy net assimilation
    real(kind=dbl_kind) :: tot_gpp(13)      ! total Gross Primary Production
    real(kind=dbl_kind) :: tot_rc(13)       ! total canopy autotrophic resp
    real(kind=dbl_kind) :: tot_fpar(13)     ! total fraction absorbed PAR
    real(kind=dbl_kind) :: tot_ss(13,nsoil) ! total soilscale
    real(kind=dbl_kind) :: tot_nee(13)      ! total NEE
    real(kind=dbl_kind) :: tot_het(13)      ! total heterotrophic resp
    real(kind=dbl_kind) :: tot_auto(13)     ! total autotrophic resp
    
    real(kind=dbl_kind) :: tcon_min         ! (W m^-1 K^-1) thermal conductivity, soil minerals 
    real(kind=dbl_kind) :: tcon_sat(nsoil)  ! (W m^-1 K^-1) thermal conductivity, saturated soil
    real(kind=dbl_kind) :: tcon_dry(nsoil)  ! (W m^-1 K^-1) thermal conductivity, dry soil     
    real(kind=dbl_kind) :: tcon_if(-nsnow+1:nsoil) ! (W m^-1 K^-1) ground/snow thermal conductivity at layer interface
    real(kind=dbl_kind) :: slamda(-nsnow+1:nsoil) ! CLM heat flux term (see begtem.F) (m^2 K W^-1)
    real(kind=dbl_kind) :: cpmin(nsoil)     ! heat capacity soil minerals (J/m^3/K)
    real(kind=dbl_kind) :: shcap(-nsnow+1:nsoil)  ! soil total heat capacity (J/m^2/K)
    real(kind=dbl_kind) :: satcap(2) ! saturation capacity depth for vegetation (1) and ground (2) (m)
    real(kind=dbl_kind) :: czc       ! canopy heat capacity (J m^-2 K-1)
    real(kind=dbl_kind) :: wilt(nsoil)      ! soil wilting point (volumetric)
    real(kind=dbl_kind) :: fieldcap(nsoil)  ! soil field capacity (volumetric)
    real(kind=dbl_kind) :: eta_wood         ! eta wood scaling factor

    !...boundary conditions-TIME VARYING
    real(kind=real_kind) :: NDVI      ! Normalized Difference Vegetation Index (-)
    real(kind=real_kind) :: NDVI1     ! Normalized Difference Vegetation Index (-)
    real(kind=real_kind) :: NDVI2     ! Normalized Difference Vegetation Index (-)
    real(kind=real_kind) :: NDVI3     ! Normalized Differenc eVegetation Index (-)
    real(kind=real_kind) :: NDVI_time1 ! NDVI time (day of year)
    real(kind=real_kind) :: NDVI_time2 ! NDVI time (day of year)
    real(kind=real_kind) :: NDVI_time3 ! NDVI time (day of year)
    integer(kind=int_kind) :: ndvi_period1  ! NDVI time (which period it's in)
    integer(kind=int_kind) :: ndvi_period2  ! NDVI time (which period it's in)
    integer(kind=int_kind) :: ndvi_period3  ! NDVI time (which period it's in)
    real(kind=dbl_kind) :: aparc      ! absorbed fraction of PAR (-)
    real(kind=dbl_kind) :: aparc1     ! absorbed fraction of PAR (-)
    real(kind=dbl_kind) :: aparc2     ! absorbed fraction of PAR (-)
    real(kind=dbl_kind) :: zlt        ! leaf area index (-)
    real(kind=dbl_kind) :: zlt1       ! leaf area index (-)
    real(kind=dbl_kind) :: zlt2       ! leaf area index (-)
    real(kind=dbl_kind) :: lai_delta  ! daily leaf area Index change rate (1/s)
    real(kind=dbl_kind) :: physfrac_del(5)   ! daily physfrac change rate (1/s)
    real(kind=dbl_kind) :: green      ! green fraction of LAI (-)
    real(kind=dbl_kind) :: green1     ! green fraction of LAI (-)
    real(kind=dbl_kind) :: green2     ! green fraction of LAI (-)
    real(kind=dbl_kind) :: z0d        ! roughness length (m)
    real(kind=dbl_kind) :: z0d1       ! roughness length (m)
    real(kind=dbl_kind) :: z0d2       ! roughness length (m)
    real(kind=dbl_kind) :: z0         ! roughness length adjust for 
    !   snow-covered canopy (m)
    real(kind=dbl_kind) :: z01        ! roughness length adjust for 
    !   snow-covered canopy (m)
    real(kind=dbl_kind) :: z02        ! roughness length adjust for 
    !   snow-covered canopy (m)
    real(kind=dbl_kind) :: zp_disp    ! zero-plane displacement (m)
    real(kind=dbl_kind) :: zp_disp1   ! zero-plane displacement (m)
    real(kind=dbl_kind) :: zp_disp2   ! zero-plane displacement (m)
    real(kind=dbl_kind) :: zpd_adj    ! zp_disp adjusted for snow on canopy (m)
    real(kind=dbl_kind) :: cc1        ! bulk pbl resistance coefficient (s/m)^0.5
    real(kind=dbl_kind) :: cc2        ! ground to CAS resistance (-)
    ! coefficient
    real(kind=dbl_kind) :: rbc        ! cc1 adjusted for snow (s/m)^0.5
    real(kind=dbl_kind) :: rbc1       ! cc1 adjusted for snow (s/m)^0.5
    real(kind=dbl_kind) :: rbc2       ! cc1 adjusted for snow (s/m)^0.5
    real(kind=dbl_kind) :: rdc        ! cc2 adjusted for snow (-)
    real(kind=dbl_kind) :: rdc1       ! cc2 adjusted for snow (-)
    real(kind=dbl_kind) :: rdc2       ! cc2 adjusted for snow (-)
    real(kind=dbl_kind) :: gmudmu     ! time-mean leaf projection (-)
    real(kind=dbl_kind) :: gmudmu1    ! time-mean leaf projection (-)
    real(kind=dbl_kind) :: gmudmu2    ! time-mean leaf projection (-)
    real(kind=dbl_kind) :: co2m_ppm   ! mixed layer CO2 concentration (ppmv)
    real(kind=dbl_kind) :: co2m_ppm1  ! mixed layer CO2 concentration (ppmv)
    real(kind=dbl_kind) :: co2m_ppm2  ! mixed layer CO2 concentration (ppmv)

    ! isotope stuff 
    real(kind=dbl_kind) :: d13cresp   ! del13C of respiration (per mil vs PDB)
    real(kind=dbl_kind) :: d13cresp1  ! del13C of respiration (per mil vs PDB)
    real(kind=dbl_kind) :: d13cresp2  ! del13C of respiration (per mil vs PDB)
    real(kind=dbl_kind) :: d13cresp3  ! del13C of respiration (per mil vs PDB)
!
! CASA related, biome dependent parameters
    real(kind=dbl_kind) q10          ! (-) base for respiration temp response function
    real(kind=dbl_kind) root_shoot   ! (-) ratio of root to leaf production
    real(kind=dbl_kind) spring_sum   ! (-) ratio of spring wood density to summer wood density
    real(kind=dbl_kind) metab        ! (-) metabolic fraction of biomass
    real(kind=dbl_kind) oxygen       ! (-) oxygen availability scaling factor
    real(kind=dbl_kind) lignin_str   ! (-) lignin fraction of surface and soil structural pools
    real(kind=dbl_kind) lignin       ! (-) lignin microbial efficiency
    real(kind=dbl_kind) lignin_cwd   ! (-) lignin fraction of coarse woody debris
    real(kind=dbl_kind) sla          ! (m2/mole C) specific leaf area
    real(kind=dbl_kind) sugar        ! (-) sugar fraction of storage pool
    real(kind=dbl_kind) fpar_max     ! (-) maximum observed fpar
    real(kind=dbl_kind) fpar_slp     ! (-) slope of fpar_max leaf growth scaling factor
!
! snow class parameters
    integer(kind=int_kind) snow_class  ! (-) snow class number
    real(kind=dbl_kind) d_min        ! (m) minimum depth hoar depth
    real(kind=dbl_kind) d_max        ! (m) maximum depth hoar depth
    real(kind=dbl_kind) f_bot_max    ! (-) maximum fraction of depth hoar
    real(kind=dbl_kind) den_ref_top  ! (kg/m3) minimum snow density
    real(kind=dbl_kind) den_ref_bot  ! (kg/m3) maximum snow density
    real(kind=dbl_kind) den_bulk_obs ! (kg/m3) observed snow bulk density
    real(kind=dbl_kind) den_bulk_std ! (kg/m3) standard dev of observed snow bulk density

end type param_vars

!------------------------------------------------------------------
!                   casa related variables
!------------------------------------------------------------------
type casa_vars
    real(kind=dbl_kind) carb_pool(npoolmax)           ! (mole/m2) pool sizes
    real(kind=dbl_kind) cpool_lay(npoolmax,nsoil)     ! (mole/m2) pool sizes per soil layer
    real(kind=dbl_kind) cpool_lay_frz(npoolmax,nsoil) ! (mole/m2) frozen pool sizes per soil layer
    real(kind=dbl_kind) carb_pool_old(npoolmax)       ! (mole/m2) pool sizes from previous time step
    real(kind=dbl_kind) pool_init(npoolmax)           ! (mole/m2) initial pool sizes
    real(kind=dbl_kind) pool_equib(npoolmax)          ! (mole/m2) estimated equilibrium pool sizes
    real(kind=dbl_kind) pool_loss(npoolmax)           ! (mole/m2) pool losses per time step
    real(kind=dbl_kind) pool_loss_lay(npoolmax,nsoil) ! (mole/m2) pool losses per time step per soil layer
    real(kind=dbl_kind) pool_gain(npoolmax)           ! (mole/m2) pool gains per time step
    real(kind=dbl_kind) pool_gain_lay(npoolmax,nsoil) ! (mole/m2) pool gains per time step per soil layer
    real(kind=dbl_kind) pool_delta(npoolmax)          ! (mole/m2) external pool gains/loss per time step
    real(kind=dbl_kind) pool_delta_lay(npoolmax,nsoil)! (mole/m2) external pool gains/loss per time step per soil layer
    real(kind=dbl_kind) resp_pool(npoolmax)           ! (mole/m2/s) respiration rate from each pool
    real(kind=dbl_kind) resp_pool_lay(npoolmax, nsoil)! (mole/m2/s) respiration rate from each pool per soil layer
    real(kind=dbl_kind) k_rate(npoolmax)              ! (1/s) unscaled pool decay rate constants
    real(kind=dbl_kind) scaled_k_rate(npoolmax)       ! (1/s) temperature/moisture scaled pool decay rate constants
    real(kind=dbl_kind) k_rate_lay(npoolmax,nsoil)    ! (1/s) total scaled pool decay rate constants per soil layer
    real(kind=dbl_kind) pool_frac(npoolmax)           ! (-) fraction of pool available for use
    real(kind=dbl_kind) biocon_eff(npoolmax,npoolmax) ! (-) biomass conversion efficiencies
    real(kind=dbl_kind) trans_frac(npoolmax,npoolmax) ! (-) transfer coefficients between pools
    real(kind=dbl_kind) resp_eff(npoolmax,npoolmax)   ! (-) respiration efficiencies
    real(kind=dbl_kind) tot_k_rate(npoolmax)  ! (1/s) sum of pool decay rate constants
    real(kind=dbl_kind) tot_assimn            ! (mole/m2) total assimn
    real(kind=dbl_kind) leaffrac              ! (-) average leaf fraction of store pool loses
    real(kind=dbl_kind) rootfrac              ! (-) average root fraction of store pool loses
    real(kind=dbl_kind) woodfrac              ! (-) average wood fraction of store pool loses
    real(kind=dbl_kind) sapwood               ! (mole/m2) the amount of sapwood in the wood pool
    real(kind=dbl_kind) lai                   ! (m2/m2) leaf area index derived from prognoistic leaf pool
    real(kind=dbl_kind) fpar                  ! (-) absorbed fraction of PAR based on prognostic LAI
    real(kind=dbl_kind) fpar_slp              ! (-) slope of fPAR based on prognostic LAI
    real(kind=dbl_kind) lai_chi_sqr           ! (-) chi squared stat between sim/obs LAI
    real(kind=dbl_kind) lai_err               ! (-) LAI uncertainty
!
! aggregate carbon pools
    real(kind=dbl_kind) carb_live    ! (mole/m2) live biomass
    real(kind=dbl_kind) carb_dead    ! (mole/m2) dead biomass
    real(kind=dbl_kind) carb_tot     ! (mole/m2) total carbon in ecosystem
    real(kind=dbl_kind) carb_litter  ! (mole/m2) litter (surface) carbon
    real(kind=dbl_kind) carb_soil    ! (mole/m2) total below surface dead carbon
    real(kind=dbl_kind) carb_awood   ! (mole/m2) above ground wood biomass
    real(kind=dbl_kind) carb_init    ! (mole/m2) total initial carbon all pools
    real(kind=dbl_kind) org_mat(nsoil)     ! (kg/m2) organic matter per soil layer
    real(kind=dbl_kind) org_mat_max(nsoil) ! (kg/m2) maximum possible organic matter per soil layer
    real(kind=dbl_kind) org_frac(nsoil)    ! (-) organic soil fraction per soil layer
    real(kind=dbl_kind) olt                ! (m) organic layer thickness
!
! disturbance or stand age related variables
    integer(kind=long_kind) dist_yr  ! (year) year in which disturbance occurs
    real(kind=dbl_kind) secd         ! (-) secondary land fraction
    logical(kind=log_kind) active    ! (-) flag indicating whether to execute SiB for this pixel
!
! derrived scaling factor data
    real (kind=dbl_kind) can_q10     ! (-) canopy temperature respiration scaling factor
    real (kind=dbl_kind) wood_grw_temp   ! (-) temperature scaling factor for wood growth
    real (kind=dbl_kind) wood_grw_moist  ! (-) soil moisture scaling factor for wood growth
    real (kind=dbl_kind) wood_grw_tot  ! (-) combined temp/soil moist scaling factor for wood growth
    real (kind=dbl_kind) soil_q10(nsoil)    ! (-) soil temperature respiration scaling factor per soil layer
    real (kind=dbl_kind) soil_frz(nsoil)    ! (-) soil respiration freezing inhibit function per soil layer
    real (kind=dbl_kind) soil_moist(nsoil)  ! (-) soil moisture respiration scaling factor per soil layer
    real (kind=dbl_kind) soil_scale(nsoil)  ! (-) combined resp scaling factor per soil layer
end type casa_vars
!
!------------------------------------------------------------------
! Seasonal diagnostics related variables
!------------------------------------------------------------------
type season_vars
    real(kind=dbl_kind) snow_melt     ! (day) DOY of snowmelt in spring
    real(kind=dbl_kind) first_snow    ! (day) DOY of first snowfall in autumn
    real(kind=dbl_kind) last_frost    ! (day) DOY of last frost in spring
    real(kind=dbl_kind) first_frost   ! (day) DOY of first frost in fall
    real(kind=dbl_kind) leaf_out      ! (day) DOY of leafout in spring
    real(kind=dbl_kind) leaf_fall     ! (day) DOY of leaf fall in autumn
    real(kind=dbl_kind) strt_thaw     ! (day) DOY of start of soil thaw in spring
    real(kind=dbl_kind) tot_thaw      ! (day) DOY of total thaw of soil column
    real(kind=dbl_kind) strt_frz      ! (day) DOY of start of soil freeze in autumn
    real(kind=dbl_kind) tot_frz       ! (day) DOY of total freezing of soil column
    real(kind=dbl_kind) max_act       ! (m) maximum active layer depth
    real(kind=dbl_kind) doy_max_act   ! (day) DOY of maximum active layer depth
    real(kind=dbl_kind) max_frz       ! (m) maximum frozen layer depth
    real(kind=dbl_kind) doy_max_frz   ! (day) DOY of maximum frozen layer depth
    real(kind=dbl_kind) deg_day       ! (degC day) cumulative degree days per month >5 degC
    real(kind=dbl_kind) cum_deg_day   ! (degC day) cumulative degree days per month
    real(kind=dbl_kind) chill_day     ! (day) cumulative chill days per month <5 degC
    real(kind=dbl_kind) cum_chill_day ! (day) cumulative degree days per month
    real(kind=dbl_kind) frosty_day    ! (day) frosty days per month
    real(kind=dbl_kind) snow_day      ! (day) snow days per month
    real(kind=dbl_kind) mean_ts       ! (K) daily mean temperature used for chill_day
    real(kind=dbl_kind) crit_deg_day  ! (degC day) critical thermal time for leaf-out
end type season_vars

!------------------------------------------------------------------
!                   PROGNOSTIC VARIABLES
!------------------------------------------------------------------
type prognostic_vars

    real(kind=dbl_kind) :: td(-nsnow+1:nsoil) ! soil temperature (K)
    real(kind=dbl_kind) :: tg        ! ground surface temp (K)
    real(kind=dbl_kind) :: ta        ! CAS temperature (K)
    real(kind=dbl_kind) :: tc        ! vegetation temperature (K)
    real(kind=dbl_kind) :: tcmin     ! min canopy temp for photosynthesis frost inhibition (K) 
    real(kind=dbl_kind) :: tha       ! CAS potential temperature (K)
    real(kind=dbl_kind) :: mra       ! CAS water vapor mixing ratio (kg/kg)
    real(kind=dbl_kind) :: ea        ! CAS water vapor pressure (hPa or mb)
    real(kind=dbl_kind) :: tke       ! turbulent kinetic energy (UNITS??)

    !...soil/snow
    real(kind=dbl_kind) :: www_liq(-nsnow+1:nsoil) ! (kg/m^2) soil liquid water
    real(kind=dbl_kind) :: www_ice(-nsnow+1:nsoil) ! (kg/m^2) soil ice
    real(kind=dbl_kind) :: vol_liq(-nsnow+1:nsoil) ! (-) soil liquid fraction of total layer volume
    real(kind=dbl_kind) :: vol_ice(-nsnow+1:nsoil) ! (-) soil ice fraction of total layer volume
    real(kind=dbl_kind) :: satfrac(-nsnow+1:nsoil) ! (-) saturation fraction of layer pore space (liq+ice)
    real(kind=dbl_kind) :: frozfrac(-nsnow+1:nsoil) ! (-) fraction of total water that is frozen
    real(kind=dbl_kind) :: liqfrac1(-nsnow+1:nsoil) ! (-) previous time step fraction of total soil volume that is liquid
    real(kind=dbl_kind) :: liqfrac2(-nsnow+1:nsoil) ! (-) current time step fraction of total soil volume that is liquid

    real(kind=dbl_kind) :: node_z(-nsnow+1:nsoil) ! (m) layer node depth
    real(kind=dbl_kind) :: z_top(-nsnow:nsoil)    ! (m) layer bottom depth
    real(kind=dbl_kind) :: z_bot(-nsnow:nsoil)    ! (m) layer top depth
    real(kind=dbl_kind) :: dz(-nsnow+1:nsoil)     ! (m) layer thickness
    real(kind=dbl_kind) :: dz_min(nsnow)          ! (m) min allowed snow layer thickness
    real(kind=dbl_kind) :: dz_max(nsnow)          ! (m) max allowed snow layer thickness
    real(kind=dbl_kind) :: dz_new(nsnow)          ! (m) min snow layer thickness to split into new layer
    real(kind=dbl_kind) :: snow_min_mass(nsnow)   ! (kg/m^2) min snow snow mass per layer
    real(kind=dbl_kind) :: snow_den(nsnow)        ! (kg/m3) snow density

    real(kind=dbl_kind) :: snow_veg   ! (kg/m^2) vegetation snow cover
    real(kind=dbl_kind) :: snow_mass  ! (kg/m^2) mass of snow on ground
    real(kind=dbl_kind) :: snow_depth ! (m) depth of snow on ground
    real(kind=dbl_kind) :: snow_bulk  ! (kg/m^3) bulk snow density of snowpack
    real(kind=dbl_kind) :: snow_fbot  ! (-) bottom density layer fraction of total snow depth
    real(kind=dbl_kind) :: snow_age   ! (tbd) non-dimensional snow age

    integer(kind=int_kind) :: nsl    ! number of (actual) snow layers
    !   THIS NUMBER IS NEGATIVE WHEN SNOW IS PRESENT
    real(kind=dbl_kind) :: capac(2)  ! vegetation and ground surface liquid
    !   water interception storage (kg m^-2)
    !   (1) canopy   (2) ground

    real(kind=dbl_kind) :: rst(6)    ! stomatal resistance (sec/m)
    !   (1) C3 photosynthesis
    !   (2) C4 photosynthesis
    !   (3) open
    !   (4) open
    !   (5) open
    !   (6) weighted sum over all phystypes

    real(kind=dbl_kind) :: pco2ap     ! (Pa) CAS CO2 partial pressure
    real(kind=dbl_kind) :: pco2ap_old ! (Pa) previous timestep pco2ap
    real(kind=dbl_kind) :: cas        ! (mole m^-2) CO2 store in Canopy air space
    real(kind=dbl_kind) :: cas_old    ! (mole m^-2) previous timestep CO2 store in Canopy air space
    real(kind=dbl_kind) :: expand     ! (mole m^-2) CAS expansion loss from previous timestep CO2 store

    !...isotope
    real(kind=dbl_kind) :: d13cca    ! del13C of canopy CO2 (per mil vs PDB)
    real(kind=dbl_kind) :: d13cm     ! del13C of ref level (per mil vs PDB)

end type prognostic_vars
!
!------------------------------------------------------------------
!                   Driver Variables
!------------------------------------------------------------------
type driver_vars
!
! input driver variables
    real(kind=dbl_kind) :: sw_dwn  ! surface incident shortwave radiation (W/m^2)
    real(kind=dbl_kind) :: sw_dwn1 ! surface incident shortwave radiation (W/m^2)
    real(kind=dbl_kind) :: sw_dwn2 ! surface incident shortwave radiation (W/m^2)

    real(kind=dbl_kind) :: dlwbot  ! surface incident longwave 
    real(kind=dbl_kind) :: dlwbot1 ! surface incident longwave 
    real(kind=dbl_kind) :: dlwbot2 ! surface incident longwave 

    real(kind=dbl_kind) :: tm      ! mixed layer temperature (K)
    real(kind=dbl_kind) :: tm1     ! mixed layer temperature (K)
    real(kind=dbl_kind) :: tm2     ! mixed layer temperature (K)

    real(kind=dbl_kind) :: mr      ! mixed layer water vapor mixing ratio (kg/kg)
    real(kind=dbl_kind) :: mr1     ! mixed layer water vapor mixing ratio (kg/kg)
    real(kind=dbl_kind) :: mr2     ! mixed layer water vapor mixing ratio (kg/kg)

    real(kind=dbl_kind) :: ps      ! surface pressure (hPa or mb)
    real(kind=dbl_kind) :: ps1     ! surface pressure (hPa or mb)
    real(kind=dbl_kind) :: ps2     ! surface pressure (hPa or mb)

    real(kind=dbl_kind) :: cupr    ! cumulus precipitation rate (mm/sec)
    real(kind=dbl_kind) :: cupr1   ! cumulus precipitation (mm)
    real(kind=dbl_kind) :: cupr2   ! cumulus precipitation (mm)

    real(kind=dbl_kind) :: lspr    ! stratiform precipitation rate (mm/sec) 
    real(kind=dbl_kind) :: lspr1   ! stratiform precipitation (mm)  
    real(kind=dbl_kind) :: lspr2   ! stratiform precipitation (mm) 

    real(kind=dbl_kind) :: spdm    ! wind speed (m/sec) 
    real(kind=dbl_kind) :: spdm1   ! wind speed (m/sec) 
    real(kind=dbl_kind) :: spdm2   ! wind speed (m/sec) 

    real(kind=dbl_kind) :: pco2m   ! mixed layer CO2 partial pressure (Pa)
    real(kind=dbl_kind) :: pco2m1  ! mixed layer CO2 partial pressure (Pa)
    real(kind=dbl_kind) :: pco2m2  ! mixed layer CO2 partial pressure (Pa)
!
! derived driver variables
    real(kind=dbl_kind) :: thm     ! mixed layer potential temperature (K)
    real(kind=dbl_kind) :: bps     ! (ps/1000)**kapa multiplying by bps turns a theta into a temperature
    real(kind=dbl_kind) :: ros     ! surface air density (kg/m^3)
    real(kind=dbl_kind) :: em      ! mixed layer water vapor pressure (hPa or mb)
    real(kind=dbl_kind) :: radvbc  ! visible beam radiation (W/m^2)
    real(kind=dbl_kind) :: radvdc  ! visible diffuse radiation (W/m^2)
    real(kind=dbl_kind) :: radnbc  ! nir beam radiation (W/m^2)
    real(kind=dbl_kind) :: radndc  ! nir diffuse radiation (W/m^2)
end type driver_vars

!------------------------------------------------------------------
!                   DIAGNOSTIC VARIABLES
!------------------------------------------------------------------
type diagnostic_vars

!itb_frost
    real(kind=dbl_kind) :: frost_stress  ! low-temp frost stress (-)
!itb_frost
    real(kind=dbl_kind) :: d_freeze  ! (m) freeze depth
    real(kind=dbl_kind) :: d_active  ! (m) active layer depth    
    real(kind=dbl_kind) :: thaw_dz   ! (m) thickness of unfrozen layer    
    real(kind=dbl_kind) :: frz_dz    ! (m) thickness of frozen layer    
    real(kind=dbl_kind) :: eastar    ! CAS saturation vapor pressure (hPa or mb)
    real(kind=dbl_kind) :: rha       ! CAS relative humidity (-)
    real(kind=dbl_kind) :: psy       ! psycrometric constant (gamma) (hPa K^-1)
    real(kind=dbl_kind) :: salb(2,2) ! total albedo
    !   (1,1) - visible, beam
    !   (1,2) - visible, diffuse
    !   (2,1) - nir, beam
    !   (2,2) - nir, diffuse
    real(kind=dbl_kind) :: salb_ave  ! average total albedo for all wavelengths (-)
    real(kind=dbl_kind) :: cas_cap_heat
    ! CAS heat capacity (J/m^2/K)
    real(kind=dbl_kind) :: cas_cap_vap
    ! CAS vapor capacity (J Pa^-1 m^-2)
    real(kind=dbl_kind) :: cas_cap_co2
    ! depth of 'canopy'  (m)

    real(kind=dbl_kind) :: cas_e_storage   ! CAS change in energy/timestep (W/m^-2)
    real(kind=dbl_kind) :: canex     ! snow depth on vegetation factor (-)
    real(kind=dbl_kind) :: wc        ! canopy wetness fraction (-)
    real(kind=dbl_kind) :: wg        ! ground wetness fraction (-)
    real(kind=dbl_kind) :: rstfac(4) ! stress factors (-)
    !  (1) leaf surface RH stress
    !  (2) rootzone water stress
    !  (3) temperature stress
    !  (4) product of factors 1-3
    real(kind=dbl_kind) :: wssp      ! water stress shape parameter. unitless.
                                     ! this variable controls the shape of the 
                                     ! water stress curve.

    real(kind=dbl_kind) :: paw_max     ! (kg m^-2) maximum plant available water (PAW) in soil column (liq only)
    real(kind=dbl_kind) :: paw(nsoil)  ! (kg m^-2) plant available water (PAW) per soil layer (liq only)
    real(kind=dbl_kind) :: paw_tot     ! (kg m^-2) plant available water (PAW) in soil column (liq only)
    real(kind=dbl_kind) :: pawfrac     ! (-) fraction of plant available water (PAW) available to plants (liq only)
    real(kind=dbl_kind) :: taw(nsoil)  ! (kg m^-2) total available water (TAW) per soil layer (liq + ice)
    real(kind=dbl_kind) :: taw_tot     ! (kg m^-2) total available water (TAW) in soil column (liq + ice)
    real(kind=dbl_kind) :: tawfrac     ! (-) fraction of total available water (PAW) available to plants (liq + ice)
    real(kind=dbl_kind) :: leaf_stress ! (-) leaf mortality drought stress


    real(kind=dbl_kind) :: areas     ! snow cover fraction (zero to one)
    real(kind=dbl_kind) :: a_areas   ! 'apparent' areas, for computation
    !   purposes (will have value 0.0 or 1.0)
    real(kind=dbl_kind) :: tsnow     ! snow surface temp (K)
    real(kind=dbl_kind) :: eff_poros(-nsnow+1:nsoil)
    ! effective porosity for liquid 
    !  (unitless)
    real(kind=dbl_kind) :: snowmelt  ! snow water converted to liquid
    !  and added to soil sfc (kg m^-2 sec^-1) 
    real(kind=dbl_kind) :: www_tot_soil  ! total soil water-all layers, water+ice (kg m^-2)
    real(kind=dbl_kind) :: roff      ! total subsurface runoff out of 
              !soil layers during whole sib timestep dtt (mm)
    real(kind=dbl_kind) :: roffo     ! overland runoff (mm)
    real(kind=dbl_kind) :: qqq       ! part of the subsurface runoff (mm)
    real(kind=dbl_kind) :: hr        ! soil surface relative humidity
    real(kind=dbl_kind) :: hrr       ! (copy) soil surface relative humidity
    real(kind=dbl_kind) :: soilscale(nsoil) 
    ! 'R-star' from Denning et al (1996) 
    !    (eqn 6) (UNITS?)

    real(kind=dbl_kind) :: soilq10(nsoil)
    ! soil temperature respiration dependence 
    !    function (UNITS?)
    real(kind=dbl_kind) :: www_inflow
    ! water inflow at ground surface 
    !    (kg m^-2 sec^-1)

    real(kind=dbl_kind) :: cu        ! momentum transfer coefficient (-)
    real(kind=dbl_kind) :: ct        ! thermal transfer coefficient (-)
    real(kind=dbl_kind) :: ustar     ! friction velocity (m sec^-1)
    real(kind=dbl_kind) :: drag(2)   ! drag (kg m^-2 sec^-1)
    real(kind=dbl_kind) :: ventmf    ! ventilation mass flux (kg m^-2 sec^-1)
    real(kind=dbl_kind) :: thvgm     ! sfc-reference height deficit of moisture
    !   UNITS ARE UNCLEAR

    real(kind=dbl_kind) :: ecmass    ! canopy evapotranspiration (kg m^-2 or 
    !                              mm water)
    real(kind=dbl_kind) :: egmass    ! ground evapotranspiration (kg m^-2 or 
    !                              mm water)
    real(kind=dbl_kind) :: chf       ! canopy heat storage flux (W m^-2)
    real(kind=dbl_kind) :: shf       ! soil heat storage flux (W m^-2)

    !...resistance
    real(kind=dbl_kind) :: ra        ! CAS-mixed layer resistance (sec/m) 
    real(kind=dbl_kind) :: rb        ! leaf-CAS resistance (sec/m)
    real(kind=dbl_kind) :: rc        ! canopy-CAS resistance (stomatal resistance + 2rb) (sec/m)
    real(kind=dbl_kind) :: rd        ! ground-CAS resistance (sec/m)
    real(kind=dbl_kind) :: rsoil     ! soil surface resistance (sec m^-1)
    real(kind=dbl_kind) :: rds       ! rd + rsoil (sec m^-1)

    real(kind=dbl_kind) :: ggl(6)    ! leaf conductance (sec/m)
    !   (1) C3 photosynthesis
    !   (2) C4 photosynthesis
    !   (3) open
    !   (4) open
    !   (5) open
    !   (6) weighted sum over all phystypes

    !...carbon

    real(kind=dbl_kind) :: pco2i(6)  ! leaf internal CO2 partial 
    !        pressure (Pa) 
    !   (1) C3 photosynthesis
    !   (2) C4 photosynthesis
    !   (3) open
    !   (4) open
    !   (5) open
    !   (6) weighted sum over all phystypes   
    real(kind=dbl_kind) :: pco2c(6)  ! chloroplast CO2 partial pressure (Pa)
    !   (1) C3 photosynthesis
    !   (2) C4 photosynthesis
    !   (3) open
    !   (4) open
    !   (5) open
    !   (6) weighted sum over all phystypes
    real(kind=dbl_kind) :: pco2s(6)  ! leaf surface CO2 partial
    !        pressure (Pa)
    !   (1) C3 photosynthesis
    !   (2) C4 photosynthesis
    !   (3) open
    !   (4) open
    !   (5) open
    !   (6) weighted sum over all phystypes

    real(kind=dbl_kind) :: thermk    ! canopy gap fraction for thermal IR 
    !                    radiation (-)
    real(kind=dbl_kind) :: tgeff     ! effective skin temp (K)
    ! (takes into effect vegetation, soil
    !       and snow)
    real(kind=dbl_kind) :: thgeff    ! effective skin potential temp (K)
    real(kind=dbl_kind) :: mrgeff    ! saturation mixing ratio of effective
    !   skin temperature (kg/kg) 
    real(kind=dbl_kind) :: radt(3)   ! net radiation (W/m^2)
    !   1) canopy leaves
    !   2) ground
    !   3) snow 
    real(kind=dbl_kind) :: swnet   ! net absorbed shortwave radiation (W/m^2)
    real(kind=dbl_kind) :: lwnet   ! net absorbed longwave radiation (W/m^2)
    real(kind=dbl_kind) :: p0        ! ground surface precip (after canopy
    !  interception, before snow/rain
    !  partition) (m/sec)
    real(kind=dbl_kind) :: pcpg_rain ! rain fraction of precip reaching ground
    !  (m/sec)
    real(kind=dbl_kind) :: pcpg_snow ! snow fraction of precip reaching ground
    !  (m/sec)
    real(kind=dbl_kind) :: cuprt     ! copy of cupr (m/sec)
    real(kind=dbl_kind) :: lsprt     ! copy of lspr (m/sec)

    real(kind=dbl_kind) :: radc3(2)  ! absorbed radiation (W/m^2)
    ! (1) - radiation absorbed by canopy
    ! (2) - radiation absorbed by ground 

    real(kind=dbl_kind) :: radfac(2,2,2)
    ! radiation absorption factors
    !   (1,1,1)
    !   (1,1,2)
    !   (1,2,1)
    !   (1,2,2)
    !   (2,1,1)
    !   (2,1,2)
    !   (2,2,1)
    !   (2,2,2)

    !...diagnostics/output
    real(kind=dbl_kind) :: hg        ! ground sfc sensible heat flux (W m^-2)
    real(kind=dbl_kind) :: hc        ! canopy sensible heat flux (W m^-2)
    real(kind=dbl_kind) :: hs        ! snow sensible heat flux (W m^-2)
    real(kind=dbl_kind) :: fss       ! CAS-BL sensible heat flux (W m^-2)
    real(kind=dbl_kind) :: fws       ! CAS-BL latent heat flux (W m^-2)
    !...ec,eg, and es are intermediate values used in delef.F and sibslv.F
    real(kind=dbl_kind) :: ec        ! canopy latent heat flux (W m^-2)
    real(kind=dbl_kind) :: eg        ! ground latent heat flux (W m^-2)
    real(kind=dbl_kind) :: es        ! snow   latent heat flux (W m^-2)


    real(kind=dbl_kind) :: egi       ! latent heat flux, ground interception
    !   (puddles) (W m^-2) 
    real(kind=dbl_kind) :: eci       ! latent heat flux, canopy interception
    !   (puddles) (W m^-2) 
    real(kind=dbl_kind) :: egs       ! latent heat flux, ground evaporation
    !    (W m^-2)
    real(kind=dbl_kind) :: ess       ! snow latent heat flux (W m^-2)
    real(kind=dbl_kind) :: ect       ! latent heat flux, canopy transpiration
    !                 (W m^-2)





    real(kind=dbl_kind) :: aparkk    ! canopy PAR use factor (-)
    real(kind=dbl_kind) :: resp_grnd ! ground respiration (mol m^-2 sec^-1)
    real(kind=dbl_kind) :: resp_tot  ! total respiration (mol m^-2 sec^-1)
    real(kind=dbl_kind) :: resp_het  ! heterotrophic resp (mol m^-2 sec^-1)
    real(kind=dbl_kind) :: resp_auto ! autotrophic respiration (mol m^-2 sec^-1)
    real(kind=dbl_kind) :: resp_can(6)  ! canopy auto resp (mol m^-2 sec^-1)
    !   (1) C3 photosynthesis
    !   (2) C4 photosynthesis
    !   (3) open
    !   (4) open
    !   (5) open
    !   (6) sum of physiological types 1-5
    real(kind=dbl_kind) :: pfd       ! incident PAR flux density 
    !           (moles m^-2 sec^-1)
    real(kind=dbl_kind) :: assim(6)  ! gross assimilation (mol m^-2 sec^-1)
    !   (1) C3 photosynthesis
    !   (2) C4 photosynthesis
    !   (3) open
    !   (4) open
    !   (5) open
    !   (6) weighted sum over all phystypes
    real(kind=dbl_kind) :: assimn(6) ! net assimilation (mol m^-2 sec^-1)
    !   (1) C3 photosynthesis
    !   (2) C4 photosynthesis
    !   (3) open
    !   (4) open
    !   (5) open
    !   (6) weighted sum over all phystypes

    real(kind=dbl_kind) :: cflux     ! carbon flux between CAS and reference 
    !   level (mol C  m^-2 sec^-1)
    real(kind=dbl_kind) :: nee     ! net ecosystem exchange (mol C  m^-2 sec^-1)
    real(kind=dbl_kind) :: npp     ! net primary productivity (mol C  m^-2 sec^-1)

!   ISOTOPE variables
    real(kind=dbl_kind) :: kiecps(6) ! Kinetic isotope effect during photosynthesis by phystype i 
    !   (1) C3 photosynthesis
    !   (2) C4 photosynthesis
    !   (3) open
    !   (4) open
    !   (5) open
    !   (6) weighted sum over all phystypes
    real(kind=dbl_kind) :: d13cassimn(6) ! del13C of CO2 assimilated by phystype i 
    !   (1) C3 photosynthesis
    !   (2) C4 photosynthesis
    !   (3) open
    !   (4) open
    !   (5) open
    !   (6) weighted sum over all phystypes
    real(kind=dbl_kind) :: c13assimn(6) ! total flux of C13 in CO2 assimilated 
    !  by phystype i 
    !   (1) C3 photosynthesis
    !   (2) C4 photosynthesis
    !   (3) open
    !   (4) open
    !   (5) open
    !   (6) weighted sum over all phystypes
    real(kind=dbl_kind) :: c12assimn(6) ! total flux of C12 in CO2 assimilated 
    !  by phystype i
    !   (1) C3 photosynthesis
    !   (2) C4 photosynthesis
    !   (3) open
    !   (4) open
    !   (5) open
    !   (6) weighted sum over all phystypes
    real(kind=dbl_kind) :: rcassimn(6) ! isotope ratio (13C/12C) of CO2 assimilated phystype i (unitless)
    !   (1) C3 photosynthesis
    !   (2) C4 photosynthesis
    !   (3) open
    !   (4) open
    !   (5) open
    !   (6) weighted sum over all phystypes

    real(kind=dbl_kind) :: flux13c   !
    real(kind=dbl_kind) :: flux12c   !
    real(kind=dbl_kind) :: flux_turb !
    real(kind=dbl_kind) :: testvar1  ! generic testing variable #1
    real(kind=dbl_kind) :: testvar2  ! generic testing variable #2
    real(kind=dbl_kind) :: testvar3  ! generic testing variable #3
    real(kind=dbl_kind) :: test2var1(nsoil)  ! generic testing variable #1
    real(kind=dbl_kind) :: test2var2(nsoil)  ! generic testing variable #2
    real(kind=dbl_kind) :: test2var3(nsoil)  ! generic testing variable #3

end type diagnostic_vars



!------------------------------------------------------------------
!                   STATUS VARIABLES
!------------------------------------------------------------------
type sib_status

    real(kind=dbl_kind) :: coszbar
    real(kind=dbl_kind) :: cosz
    real(kind=dbl_kind) :: julday
    integer(kind=int_kind) :: pt_num  ! (-) point number in nSiB vector

end type sib_status


!------------------------------------------------------------------
!                   TOP LEVEL ENCAPSULATING TYPE 'SIB_T'
!------------------------------------------------------------------
type sib_t

    type(param_vars) :: param          ! parameter variables
    type(prognostic_vars) :: prog      ! prognostic variables
    type(driver_vars) :: drvr          ! driver weather data variables
    type(casa_vars) :: casa            ! casa variables
    type(season_vars) :: seas          ! seasonal diagnostic variables
    type(diagnostic_vars) :: diag      ! diagnostic variables
    type(sib_status) :: stat           ! sib status variables

end type sib_t




!------------------------------------------------------------------
!                   LOCAL VARIABLES
!------------------------------------------------------------------
type sib_local_vars


    !...canopy air space (CAS) and canopy
    real(kind=dbl_kind) :: dtg       ! delta ground surface temp (K)
    real(kind=dbl_kind) :: dtd(-nsnow+1:nsoil)
    ! delta soil temperature (K)
    real(kind=dbl_kind) :: dtc       ! change in canopy temperature (K)
    real(kind=dbl_kind) :: dts       ! change in snow surface temperature (K)
    real(kind=dbl_kind) :: dth       ! change in ref level temperature (K)
    real(kind=dbl_kind) :: dqm       ! change in ref level moisture (Pa)
    real(kind=dbl_kind) :: dta       ! change in CAS temperature (K)
    real(kind=dbl_kind) :: dea       ! change in CAS moisture (Pa)

    real(kind=dbl_kind) :: etc       ! saturation vapor pressure at Tc (hPa)
    !   ('e-star' of Tc)
    real(kind=dbl_kind) :: getc      ! derivative of etc with respect to temp
    !   (d(etc)/dTc (hPa K^-1)
    real(kind=dbl_kind) :: etg       ! 'e-star' of ground surface (Pa)
    real(kind=dbl_kind) :: getg      ! d(etg)/dTg (hPa K^-1)
    real(kind=dbl_kind) :: ets       ! 'e-star' of snow surface (Pa)
    real(kind=dbl_kind) :: gets      ! d(ets)/dTs (hPa K^-1)

    real(kind=dbl_kind) :: fc        ! direction of vapor flux, canopy to CAS
    !  =1 when flux is canopy to CAS (etc>ea)
    !  =0 when flux is CAS to canopy (ea>etc)
    real(kind=dbl_kind) :: fg        ! direction of vapor flux, ground to CAS
    !  =1 when flux is ground to CAS (etg>ea)
    !  =0 when flux is CAS to ground (ea>etg)
    real(kind=dbl_kind) :: gect      ! dry fraction of veg / rc
    real(kind=dbl_kind) :: geci      ! wet fraction of veg / 2rb
    real(kind=dbl_kind) :: gegs      ! dry fraction of ground / rds
    real(kind=dbl_kind) :: gegi      ! wet fraction of ground /rd
    real(kind=dbl_kind) :: coc       ! gect + geci
    real(kind=dbl_kind) :: cog1      ! gegi + gegs*hrr
    real(kind=dbl_kind) :: cog2      ! gegi + gegs


    !...radiation


    integer(kind=int_kind) :: imelt(-nsnow+1:nsoil)
    ! flag for melting/freezing
    !   1=> melting, 2=> freezing
    real(kind=dbl_kind) :: frac_iceold(-nsnow+1:nsoil)
    ! start-of-timestep value for
    !   ice fraction of total liquid
    real(kind=dbl_kind) :: dtc4      ! d(canopy thermal em)/dT (W m^-2 K^-1)
    real(kind=dbl_kind) :: dtg4      ! d(ground thermal em)/dT (W m^-2 K^-1)
    real(kind=dbl_kind) :: dts4      ! d(snow thermal em)/dT   (W m^-2 K^-1)
    real(kind=dbl_kind) :: lcdtc     ! d(canopy thermal em)/dtc (W m^-2 K^-1)
    real(kind=dbl_kind) :: lcdtg     ! d(canopy thermal em)/dtg (W m^-2 K^-1)
    real(kind=dbl_kind) :: lcdts     ! d(canopy thermal em)/dts (W m^-2 K^-1)
    real(kind=dbl_kind) :: lgdtc     ! d(ground thermal em)/dtc (W m^-2 K^-1)
    real(kind=dbl_kind) :: lgdtg     ! d(ground thermal em)/dtg (W m^-2 K^-1)
    real(kind=dbl_kind) :: lsdts     ! d(snow thermal em)/dts   (W m^-2 K^-1)
    real(kind=dbl_kind) :: lsdtc     ! d(snow thermal em)/dtc   (W m^-2 K^-1)
    real(kind=dbl_kind) :: hcdtc     ! d(canopy H)/dtc  (W m^-2 K^-1)
    real(kind=dbl_kind) :: hcdta     ! d(canopy H)/dta  (W m^-2 K^-1)
    real(kind=dbl_kind) :: hgdta     ! d(ground H)/dta  (W m^-2 K^-1)
    real(kind=dbl_kind) :: hgdtg     ! d(ground H)/dtg  (W m^-2 K^-1)
    real(kind=dbl_kind) :: hsdta     ! d(snow H)/dta  (W m^-2 K^-1)
    real(kind=dbl_kind) :: hsdts     ! d(snow H)/dtsnow  (W m^-2 K^-1)
    real(kind=dbl_kind) :: hadta     ! d(CAS H)/dta  (W m^-2 K^-1)
    real(kind=dbl_kind) :: hadth     ! d(CAS H)/dtheta  (W m^-2 K^-1)
    real(kind=dbl_kind) :: ecdtc     ! d(canopy LE)/dtc (W m^-2 K^-1)
    real(kind=dbl_kind) :: ecdea     ! d(canopy LE)/dea (W m^-2 K^-1)
    real(kind=dbl_kind) :: egdtg     ! d(ground LE)/dtg (W m^-2 K^-1)
    real(kind=dbl_kind) :: egdea     ! d(ground LE)/dea (W m^-2 K^-1)
    real(kind=dbl_kind) :: esdts     ! d(snow LE)/dtsnow (W m^-2 Pa^-1)
    real(kind=dbl_kind) :: esdea     ! d(snow LE)/dea (W m^-2 K^-1)
    real(kind=dbl_kind) :: eadea     ! d(CAS LE)/dea (W m^-2 Pa^-1)
    real(kind=dbl_kind) :: eadem     ! d(CAS LE)/dem (W m^-2 Pa^-1)
    real(kind=dbl_kind) :: closs     ! canopy thermal loss     (W m^-2)
    real(kind=dbl_kind) :: gloss     ! ground thermal loss     (W m^-2)
    real(kind=dbl_kind) :: sloss     ! snow thermal loss       (W m^-2)
    real(kind=dbl_kind) :: fac1      ! effective ground cover for 
    !   thermal radiation (-)

    real(kind=dbl_kind) :: td_old(-nsnow+1:nsoil)
    ! prev timestep soil temperature (K)

end type sib_local_vars
!
! Casa pool configuration variables common to all biomes
type pool_config
  character*25 name        ! (-) name of pool configuration
  character*100 desc       ! (-) short description of pool configuration
  character*25 pool_name(npoolmax)       ! (-) pool names
  character*10 pool_name_short(npoolmax) ! (-) short pool names
  character*10 pool_type(npoolmax)       ! (-) pool types (live, dead, storage) for het/auto resp
  character*10 pool_loc(npoolmax)        ! (-) pool locations (soil, surface, canopy) for temp/moist scaling
  real(kind=dbl_kind) turnover(npoolmax,nbiomax)     ! (year) pool turnover times
  real(kind=dbl_kind) k_rate(npoolmax,nbiomax)       ! (1/s) instantaneous pool decay rate constants
  real(kind=dbl_kind) metab(nbiomax)                 ! (-) metabolic fraction of biomass
  real(kind=dbl_kind) lignin(nbiomax)                ! (-) lignin fraction
  real(kind=dbl_kind) q10(nbiomax)          ! (-) base for respiration temp response function
  real(kind=dbl_kind) root_shoot(nbiomax)   ! (-) ratio of root to leaf production
  real(kind=dbl_kind) spring_sum(nbiomax)    ! (-) ratio of spring wood density to summer wood density
  real(kind=dbl_kind) oxygen(nbiomax)       ! (-) oxygen availability scaling factor
  real(kind=dbl_kind) lignin_str(nbiomax)   ! (-) lignin fraction of surface and soil structural pools
  real(kind=dbl_kind) lignin_cwd(nbiomax)   ! (-) lignin fraction of coarse woody debris
  real(kind=dbl_kind) lma(nbiomax)          ! (g/m2) leaf mass per area
  real(kind=dbl_kind) sla(nbiomax)          ! (m2/mol C) specific leaf area
  real(kind=dbl_kind) sugar(nbiomax)        ! (-) sugar fraction of storage pool
  real(kind=dbl_kind) pool_frac(npoolmax,nbiomax)    ! (-) fraction of pool available for use
  real(kind=dbl_kind) biocon_eff(npoolmax,npoolmax,nbiomax) ! (-) biomass conversion efficiencies
  real(kind=dbl_kind) trans_frac(npoolmax,npoolmax,nbiomax) ! (-) transfer coefficients between pools
  real(kind=dbl_kind) resp_eff(npoolmax,npoolmax,nbiomax)   ! (-) respiration efficiencies 
  integer nlive_pool  ! (-) number of live biomass pools
  integer ndead_pool  ! (-) number of dead pools
  integer nsoil_pool  ! (-) number of pools in the soil
  integer nsurf_pool  ! (-) number of pools on the surface
  integer ncan_pool   ! (-) number of pools in the canopy
  integer ntrans      ! (-) total number of pool-to-pool transfers per time step
  integer indx_live(npoolmax)  ! (-) pool index numbers for all live pools
  integer indx_dead(npoolmax)  ! (-) pool index numbers for all dead pools
  integer indx_soil(npoolmax)  ! (-) pool index numbers for all soil pools
  integer indx_surf(npoolmax)  ! (-) pool index numbers for all surface pools
  integer indx_can(npoolmax)   ! (-) pool index numbers for all canopy pools
  integer indx_trans(npoolmax*npoolmax,2)   ! (-) pool-to-pool transfer indices
  integer indx_stor  ! (-) pool index number for storage pool
  integer indx_leaf  ! (-) pool index number for leaf pool
  integer indx_root  ! (-) pool index number for root pool
  integer indx_wood  ! (-) pool index number for wood pool
  integer indx_cwd   ! (-) pool index number for coarse woody debris pool
  integer indx_surfmic  ! (-) pool index number for surface microbial pool
  integer indx_slow  ! (-) pool index number for slow pool
end type pool_config
type(pool_config) config
!
! biome parameter table
type bio_config
  real(kind=dbl_kind) :: z2        ! (m) canopy top
  real(kind=dbl_kind) :: z1        ! (m) canopy bottom
  real(kind=dbl_kind) :: chil      ! (-) leaf angle distribution factor
  real(kind=dbl_kind) :: phc       ! (m) 1/2 crit leaf water pot limit
  real(kind=dbl_kind) :: tran(2,2) ! leaf transmittance (-)
                                   !  (1,1) - shortwave, green plants
                                   !  (1,2) - longwave, green plants
                                   !  (2,1) - shortwave, brown plants
                                   !  (2,2) - longwave, brown plants
  real(kind=dbl_kind) :: ref(2,2)  ! leaf reflectance (-)
                                   !  (1,1) - shortwave, green plants
                                   !  (1,2) - longwave, green plants
                                   !  (2,1) - shortwave, brown plants
                                   !  (2,2) - longwave, brown plants
  real(kind=dbl_kind) :: vmax0     ! (mol/m^2/sec) Rubisco vel of sun leaf
  real(kind=dbl_kind) :: effcon    ! quantum efficiency (mol/mol)
  real(kind=dbl_kind) :: gradm     ! conductance-photosynthesis slope parameter (-)
  real(kind=dbl_kind) :: binter    ! conductance-photosynthesis intercept (mol m^-2 sec^-1)
  real(kind=dbl_kind) :: atheta    ! WC WE coupling parameter (-)
  real(kind=dbl_kind) :: btheta    ! WC&WE, WS coupling parameter (-)
  real(kind=dbl_kind) :: trda      ! temp coeff in GS-A model (K^-1)
  real(kind=dbl_kind) :: trdm      ! temp coeff in GS-A model (K)
  real(kind=dbl_kind) :: trop      ! temp coeff in GS-A model (K)
  real(kind=dbl_kind) :: respcp    ! respiration fraction of vmax0 (-)
  real(kind=dbl_kind) :: slti      ! slope of lo-temp inhibition (K^-1) function (K^-1)
  real(kind=dbl_kind) :: hltii     ! 1/2 point of lo-temp inhibition (K) function
  real(kind=dbl_kind) :: shti      ! slop of hi-temp inhibition (K^-1) function (K-1)
  real(kind=dbl_kind) :: hhti      ! 1/2 point of hi-temp inhibition (K) function
  real(kind=dbl_kind) :: sftf      ! (1/K) slope leaf/root mortality frost function
  real(kind=dbl_kind) :: hftf      ! (K) half point leaf/root mortality frost function
  real(kind=dbl_kind) :: fpar_slp  ! (-) slope of fpar_max leaf growth scaling factor
  real(kind=dbl_kind) :: sfrzi     ! (1/K) slope respiration freezing scaling factor
  real(kind=dbl_kind) :: hfrzi     ! (K) half point respiration freezing scaling factor
  real(kind=dbl_kind) :: tfrost    ! (K) temp below which frost reduces photosynthesis 
  real(kind=dbl_kind) :: tcmin_rr  ! (K/day) min canopy temp recovery rate  
  real(kind=dbl_kind) :: kroot     ! (1/m) root density extinction coeficient
  real(kind=dbl_kind) :: rootd     ! (m) maximum rooting depth
end type bio_config
type(bio_config), allocatable :: C3_biotab(:)
type(bio_config), allocatable :: C4_biotab(:)

end module sibtype
