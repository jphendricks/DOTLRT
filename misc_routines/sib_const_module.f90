!...module to hold values for constants that are not contained in the
!...BUGS module physical_parameters.F

module sib_const_module

    use kinds
    use physical_parameters

    implicit none
    save

    real(kind=dbl_kind) :: dtt       ! model time step (seconds) 
    real(kind=dbl_kind) :: dti       ! inverse time step

    !--variables that retain a constant value throughout the simulation
    integer(kind=int_kind) ::     &
        nsib,           & !  number of SiB points in datasets
        subcount,       & !  actual number of SiB point in simulation
        snowl,          & !  number of (actual) snow layers
        ihr,            & !  global points in x-direction
        jhr,            & !  global points in y-direction
        nhr,            & !  ihr*jhr (total global points)
        nper              ! actual number of ndvi composite periods
integer(kind=int_kind) :: no_nsib ! number of sib points in northern hemisphere
integer(kind=int_kind) :: so_nsib ! number of sib points in southern hemisphere
        
    integer(kind=long_kind) :: &
        endtime,        & !  end time of integration -- units can vary
        starttime,      & !  start time of integration -- units can vary
        dtsib,          & !  timestep in seconds
        dtsibmetin,     & !  driver data input interval (seconds)
        dtsibres,       & !  restart interval (see namel_sibdrv for exp)
        dtsibout          !  output interval (see namel_sibdrv for exp)

     integer(kind=int_kind) :: &
        dtsibbcin,      & !  sib boundary condition input interval
        nsecond,        & !  simulation time, in seconds
        numsib,         & !  check of # of sib points--init_sibdrv
        nstepsib          !  number of timesteps integrated

    integer(kind=int_kind) :: endyear
    integer(kind=int_kind) :: ndtsibpbp

    !itb...to nsib points in init_sibdrv.F
    real(kind=real_kind), dimension(:), allocatable    ::    &
        latsib,       & !  SiB point latitude
        lonsib,       & !  SiB point longitude
        latitude,     &
        longitude,    &
        lonpbp,       &
        latpbp

real(kind=real_kind), parameter :: layer1    = 0.02  ! (m) thickness of first soil layer
real(kind=real_kind), parameter :: layer_rat = 1.229 ! (-) soil layer thickness ratio
integer(kind=int_kind), parameter :: nsoil   = 25  ! number of soil layers 
integer(kind=int_kind), parameter :: nsnow   = 5   ! max number of snow layers
integer(kind=int_kind), parameter :: physmax = 5   ! max number physiology types
                                                   ! (only C3 and C4 now, but capability for more)
integer(kind=int_kind), parameter :: npermax = 365 ! max number ndvi composite periods per year
integer(kind=int_kind), parameter :: npoolmax= 20  ! max number of carb pools
integer(kind=int_kind) :: npool                    ! number of carbon pools
integer(kind=int_kind), parameter :: nbiomax = 13  ! max number of biome types
integer(kind=int_kind) :: nbiome                   ! number of biome

!time variables
real(kind=real_kind) sin_dec  ! (-) sin solar declination
real(kind=real_kind) cos_dec  ! (-) cosine solar declination
real(kind=real_kind) lonearth ! (rad) Earth lon about Sun from vernal equinox
real(kind=real_kind) tau      ! (TBD) time
    integer(kind=int_kind) ::     & ! 
        startyear = 1,   & !  time manager stuff
        eqnx   = 80     !  day of vernal equinox

integer(kind=int_kind), dimension(:), allocatable :: subset   !  array of landpoint indices for subgrid
integer(kind=int_kind), dimension(:), allocatable :: latindex !  latitude index array of all landpoints
integer(kind=int_kind), dimension(:), allocatable :: lonindex !  longitude index array of all landpoints
integer(kind=int_kind), dimension(:), allocatable :: sublat   !  latitude index array of subset
integer(kind=int_kind), dimension(:), allocatable :: sublon   !  longitude index array of subset
integer(kind=int_kind), dimension(:), allocatable :: nohem_sib_indx  ! indx of sib points in northern hemisphere
integer(kind=int_kind), dimension(:), allocatable :: sohem_sib_indx  ! indx of sib points in southern hemisphere

    !itb...some SCALARS
real(kind=dbl_kind) :: c3day         !  timesteps per day
real(kind=dbl_kind) :: ztemp         !  height of temperature measurement (m)
real(kind=dbl_kind) :: zwind         !  height of wind measurement (m)

real(kind=dbl_kind), parameter :: tcon_ihoar = 0.183           !  (W/m/K) thermal conductivity of indurated depth hoar
real(kind=dbl_kind), parameter :: tcon_hoar = 0.072            !  (W/m/K) thermal conductivity of depth hoar
real(kind=dbl_kind), parameter :: den_snowfall_min = 50.       !  (kg/m3) minimum density of new snowfall
real(kind=dbl_kind), parameter :: version = 3.0                !  code version identifier
real(kind=dbl_kind), parameter :: snomel  = 3.705185e8         !  latent heat of fusion of ice (J m^-3) 
real(kind=dbl_kind), parameter :: hfus   = snomel/1000.0       !  latent heat of fusion of ice (J kg^-1) 
real(kind=dbl_kind), parameter :: cv     = 1952.0              !  (J deg^-1 kg^-1) specific heat of water vapor at constant pressure 
real(kind=dbl_kind), parameter :: cpice  = 2117.27             !  (J kg^-1 deg^-1) specific heat of ice 
real(kind=dbl_kind), parameter :: cpliq  = 4188.0              !  (J kg^-1 deg^-1) spec heat of liquid water 
real(kind=dbl_kind), parameter :: cpom   = 2.5e6               !  (J m^-3 deg^-1) spec heat of soil organic matter 
real(kind=dbl_kind), parameter :: clai   = 4.186*1000.0*0.2    !  (J m^-2 deg^-1) leaf heat capacity  
real(kind=dbl_kind), parameter :: cww    = 4.186*1000.0*1000.0 !  (J m^-3 deg^-1) water heat capacity 
real(kind=dbl_kind), parameter :: asnow  = 16.7                !  UNKNOWN
real(kind=dbl_kind), parameter :: rotper = 24.0                !  hours per day
real(kind=dbl_kind), parameter :: day    = rotper * 3600.0     !  seconds per day
real(kind=dbl_kind), parameter :: vkrmn  = 0.35                !  Von Karmann's constant (unitless)
real(kind=dbl_kind), parameter :: ribc   = 3.05                !  critical Richardson Number (unitless)
real(kind=dbl_kind), parameter :: pr0    = 0.74                !  turb Prandtl Number at neutral stblty
real(kind=dbl_kind), parameter :: tkemin = 0.01                !  minimum allowed value for tke
real(kind=dbl_kind), parameter :: rgfac  = 100.0/gas_const_r   !  
real(kind=dbl_kind), parameter :: cpdgrv = spec_heat_cp/grav   ! 
real(kind=dbl_kind), parameter :: po2m   = 20900.0             !  mixed layer O2 concentration
real(kind=dbl_kind), parameter :: perhl  = 102.7               !  (deg) longitude of perihelion     
real(kind=dbl_kind), parameter :: den_liq = 1000.0             !  (kg/m^3) density of liquid water
real(kind=dbl_kind), parameter :: den_ice = 917.0              !  (kg/m^3) density of ice 
real(kind=dbl_kind), parameter :: den_om_max = 140.0           !  (kg/m^3) maximum possible density of organic soil 
real(kind=dbl_kind), parameter :: den_min = 2700.0             !  (kg/m^3) density of mineral soil
real(kind=dbl_kind), parameter :: tcon_air  = 0.023            !  (W/m/K) thermal conductivity of air
real(kind=dbl_kind), parameter :: tcon_liq  = 0.6              !  (W/m/K) thermal conductivity of liquid water
real(kind=dbl_kind), parameter :: tcon_ice  = 2.29             !  (W/m/K) thermal conductivity of ice
real(kind=dbl_kind), parameter :: tcon_om   = 0.25             !  (W/m/K) thermal conductivity of soil organic matter
real(kind=dbl_kind), parameter :: tcon_om_dry   = 0.05         !  (W/m/K) thermal conductivity of dry organic soil
real(kind=dbl_kind), parameter :: snofac = hltm/(hltm + snomel * 1.E-3)  !  ratio of hltm to hltm+ht 
                                                               ! of fusion (see Sellers (1986) appendix B)
real(kind=dbl_kind), parameter :: wimp = 0.05                  !  water impermeable if effective porosity below this value
real(kind=dbl_kind), parameter :: phmin = -1.e8                !  minimum value for soil potential (mm)
real(kind=dbl_kind), parameter :: vwcmin = 0.1                 !  wilting point volumetric water content
real(kind=dbl_kind), parameter :: wtfact = 0.3                 !  fraction of area with high 
                                                               !  (HARDWIRE PATCH) water table
real(kind=dbl_kind), parameter :: ssi = 0.033                  !  irreducible water fraction of snow
real(kind=dbl_kind), parameter :: zlnd = 0.01                  !  roughness length for land (m)
real(kind=dbl_kind), parameter :: eccn   = 0.016715            !  eccentricity
real(kind=dbl_kind), parameter :: daypyr = 365.0               !  days per  year
real(kind=dbl_kind), parameter :: decmax = 23.441              !  max declination
real(kind=dbl_kind), parameter :: cn_fact = 0.5                !  Crank-Nicholson factor
real(kind=dbl_kind), parameter :: cosz_min = -0.1045           !  (-) min cosine of zenith angle
                                                               !  -0.1045 is 96 deg. which includes
                                                               !  civil twilight
real(kind=dbl_kind), parameter :: mwc = 12.                    ! (g mole^-1) molecular weight of carbon
real(kind=dbl_kind), parameter :: micro = 1.e6                 ! (microunits unit^-1) conversion to micro units
real(kind=dbl_kind), parameter :: rat_org_carb = 2.            ! (-) ratio of organic matter to carbon
real(kind=dbl_kind) :: k_sapwood = 1./20./365./24./3600.    ! (s^-1) sapwood rate constant (1/ 20 yr turnover time)
real(kind=dbl_kind) :: sapwood_frac = 0.7        ! (-) sapwood wood fraction
real(kind=dbl_kind) :: r_shadesun = 0.5          ! (-) ratio of shade leaf mass to sun leaf mass
real(kind=dbl_kind) :: tref_q10 = 298.15         ! (K) reference temperature for respiration q10 function
real(kind=dbl_kind) :: hfrzi = 271.              ! (K) half point respiration freeze temperature inhibition function
real(kind=dbl_kind) :: sfrzi = 1.                ! (1/K) slope respiration freeze temperature inhibition function
real(kind=dbl_kind), parameter :: rho_wood = 0.5               ! (g/cm3) wood density
real(kind=dbl_kind), parameter :: autofrac=0.5                 ! (-) fraction GPP to autotrophic respiration
real(kind=dbl_kind), parameter :: t_base=278.15         ! (K) base temperature for degree and chill days


!...CFRAX...CFRAX...CFRAX...CFRAX...CFRAX...CFRAX...CFRAX...CFRAX...CFRAX

!...Carbon isotopic fractionation constants (units = per mil)         
!...KIEC refers to Kinetic Isotope Effect (KIE) for Carbon (C), 
!...and can be converted to alpha notation by alpha = (1 - KIEC/1000). 
!...For a chemical reaction, alpha = Rreactant/Rproduct.  
!...KIEs are sometimes referred to as epsilon factors. 

real(kind=dbl_kind),parameter :: pdb = 0.0112372   ! 13C/12C ratio of Pee Dee Belemnite (no units)
real(kind=dbl_kind),parameter :: kieclfbl  = - 2.9 ! canopy air space to leaf boundary layer
real(kind=dbl_kind),parameter :: kiecstom  = - 4.4 ! leaf boundary layer to stomatal cavity
real(kind=dbl_kind),parameter :: kieclphas =  -0.7 ! liquid phase fractionation 
real(kind=dbl_kind),parameter :: kiecdis   =  -1.1 ! dissolution
real(kind=dbl_kind),parameter :: kiecrbsco = -28.2 ! C3 C-fixation enzyme rubisco
real(kind=dbl_kind),parameter :: tref = 298.16     ! standard temperature (K)
real(kind=dbl_kind),parameter :: pref = 101325.0   ! standard pressure (Pa)

end module sib_const_module
