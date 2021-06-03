!===================================================
  Module scan_Variables
!===================================================
! defines all variables used in scan
  use dotlrt_variables, only: pi

  IMPLICIT NONE

! program control variables
  character*25 scantype   ! type of scan
  logical calc_psuedo ! flag to calculate psuedo data
  logical calc_init  ! flag to calculate psuedo data
  Integer DevFlag    ! send plot to screen or postscript file
  Integer PlotFlag   ! plot data 
  Integer color_con  ! plot color contours 
  integer nline      ! number of plotted variables
  integer max_nline  ! max number of plotted variables
  integer LegFlag    ! flag to plot legends
  integer nLeg       ! number of legends
  integer flagDim    ! scan in 2-D or 3-D
  integer save_2d    ! save 2-dimensional data from scan
  integer save_3d    ! save 3-dimensional data from scan
  integer Yvar       ! flag for which Y variable to vary
  integer Xvar       ! flag for which X variable to vary
  integer numpts     ! total number of points in scan
  integer Yscale     ! auto or manual y value scaling
  character*8 barformat ! format for color bar labels
  Character *45 Junk ! junk variable for reading input descriptors
  Character *100 Plotfile ! file name for output plot

! channel variables
  integer nchannel  ! (-) total number of channels
  Integer minchan   ! (-) minimum channel number
  Integer maxchan   ! (-) maximum channel number

! scan output variables
  character*25 LabX                  ! Label for x-axis
  character*25 LabY                  ! Label for y-axis
  character*25 legend(100)           ! legends for plotted lines
  Character*100 Titlebase            ! title for individual plot
  Character*100 Title                ! title for individual plot
  Integer Scale                      ! Contour interval scaling factor for plotting
  real(8) base                       ! base for plotting contours
  integer NC                         ! Number of Contours in plot
  real(8), allocatable :: X(:,:)     ! generic x values
  real(8), allocatable :: Y(:,:)     ! generic Y values
  real(8), allocatable :: X1(:,:)    ! generic x values
  real(8), allocatable :: Y1(:,:)    ! generic Y values
  real(8), allocatable :: X3(:)      ! generic x values
  real(8), allocatable :: Y3(:)      ! generic y values
  real(8), allocatable :: zval(:,:,:)! generic Z values
  integer NumScanVar                 ! number of variables in table
  character (len=25), allocatable :: Label(:) ! labels for plotting
  real(8), allocatable :: Range(:)   ! range of scanned data values
  real(8), allocatable :: Start(:)   ! start or minimum scanned data value
  real(8), allocatable :: Stop(:)    ! stop or maximum scanned data value
  real(8), allocatable :: PlotMin(:) ! minimum value for plotting
  real(8), allocatable :: PlotMax(:) ! maximum value for plotting
  real(8), allocatable :: Typical(:) ! typical value when not scanning
  real(8), allocatable :: Value(:)   ! current value
  real(8) Dx
  real(8) Dy
  real(8) Dz
  real(8) Xmin, Xmax, Ymin, Ymax, Xval, Yval
  real(8) mincont, maxcont
  integer maxi, maxj
  real(8) test, test10(100) ! testing variables
  Character *40 file  ! Filename

! grid parameters
  real(8) lat    ! (deg) center latitude of grid cell
  real(8) lon    ! (deg) center Longitude of grid cell
  real(8) LatMin ! (deg) Minimum domain latitude (lower left domain corner)
  real(8) LonMin ! (deg) Minimum domain longitude (lower left domain corner)
  real(8) Dlat   ! (deg) Latitude Grid Spacing
  real(8) Dlon   ! (deg) Longitude grid spacing
  integer imin   ! (-) minimum longitude index for subset of domain grid
  integer imax   ! (-) maximum longitude index for subset of domain grid
  integer jmin   ! (-) minimum latitude index for subset of domain grid
  integer jmax   ! (-) maximum latitude index for subset of domain grid

! profile setup specifications
  integer n_psuedo      ! (-) number variables for for pseudo-data calculation
  integer n_init        ! (-) number variables for for pseudo-data calculation
  integer nout_man      ! (-) number of manipulations
  integer n_up_man      ! (-) number of manipulations in update
  type manipulation
    logical doit        ! perform the manipulation
    character*20 typ    ! type
    real(8) val1        ! 1st value
    real(8) val2        ! 2nd value
    real(8) val3        ! 3rd value
    real(8) val4        ! 4th value
    integer ind1        ! index value 1
    integer ind2        ! index value 2
    integer ind3        ! index value 3
    integer npath1      ! 1st path number
    integer npath2      ! 2nd path number
    integer npath3      ! 3rd path number
    Character*100 path1 ! 1st path
    Character*100 path2 ! 2st path
    Character*100 path3 ! 3st path
    integer null1       ! index of 1st null value
    integer null2       ! index of 2nd null value
    integer null3       ! index of 3rd null value
    character*20 txt1   ! text 1
    character*20 txt2   ! text 2
    character*20 txt3   ! text 3
    character*20 txt4   ! text 4
    logical saveit      ! save to file
  end type manipulation
  type(manipulation) out_man(1000)    ! manipulations for output
  type(manipulation) psuedo_man(1000) ! manipulations for psedo data
  type(manipulation) init_man(1000)   ! manipulations for initial profile
  type(manipulation) update_man(1000) ! manipulations for profile update

! biophysical parameters
  real(8), parameter :: Plank        = 6.63e-34                  ! (J s) plank's constant
  real(8), parameter :: Boltz        = 1.381e-23                 ! (J/k/molecule) Boltzman's constant
  real(8), parameter :: light        = 2.9973e8                  ! (m/s) speed of light
  real(8), parameter :: pidaypy      = 0.0172142d0               ! (?) TBD
  real(8), parameter :: a            = 6.37122d+06               ! (m) earth radius
  real(8), parameter :: grav         = 9.8100d0                  ! (m/s) gravity constant
  real(8), parameter :: gravi        = 1.0d0/grav                ! (s/m) inverse grav constant
  real(8), parameter :: omega        = 2.0d0*pi/86400.0d0        ! (1/s) earth angular velocity 
  real(8), parameter :: earth_area   = 5.100996990707616d+14     ! (m2) surface area of earth 
  real(8), parameter :: gas_const_R  = 287.000d0                 ! (J/kg/K))gas constant for dry air 
  real(8), parameter :: heat_cp      = 1005.000d0                ! (J/kg/K)) specific heat at constant pressure 
  real(8), parameter :: kappa        = gas_const_R/heat_cp       ! (?) gas_const_R/spec_heat_cp
  real(8), parameter :: inv_kappa    = 1.0d0 / kappa             ! (?) inverse kappa
  real(8), parameter :: p0_sfc       = 1.0d+05                   ! (Pa) surface pressure
  real(8), parameter :: inv_p0_sfc   = 1.0d0 / p0_sfc            ! (1/Pa) inverse surfae pressure
  real(8), parameter :: tice         = 273.15d0                  ! (K) freezing temperature of water at 1 atm
  real(8), parameter :: hltm         = 2.25d+06                  ! (J/kg) latent heat of vaporization
  real(8), parameter :: gamfac       = hltm*5417.9827d0/heat_cp  ! a moist thermodynamic variable
  real(8), parameter :: delta        = 0.608d0                   ! (-) molecular_weight_air/molecular_weight_water
  real(8), parameter :: stefan       = 5.67d-08                  ! (W/m^2/K^4) stefan boltzmann constant 
  real(8), parameter :: rv           = 4.61d+02                  ! gas constant for water vapor
  real(8), parameter :: pi2          = 2*pi                      ! (-) 2*pi
  real(8), parameter :: dtr          = pi/180.0d0                ! (deg/rad) Degrees To Radians conversion constant
  real(8), parameter :: rtd          = 180.0d0/pi                ! (rad/deg) Radians To Degrees conversion constant
  real(8), parameter :: alpha2       = 0.0                       ! alpha2 = rotation of axis of rotation from NP/SP
  real(8), parameter :: fPARmax      = 0.95d0                    ! (-) Max FPAR (98th percentile)
  real(8), parameter :: fPARmin      = 0.01d0                    ! (-) Min FPAR (2nd percentile)
  real(8), parameter :: rho_water    = 1000.                     ! (kg/m3) density of water
  real(8), parameter :: rho_ice      = 917.                      ! (kg/m3) density of ice

! Parameters turned off when linking to DOTLRT
!    real(8), parameter :: pi           = 3.14159265358979323846d0  ! (-) pi

! begin execution time parameters
  Integer Tstart    ! start execution time in internal units
  Integer Tend      ! stop execution time in internal units
  Integer clock2    ! Conversion factor betw internal and wall clock time
  real(8) Time         ! Execution time in seconds
  real(8) TotalTime         ! Execution time in seconds

! input/output files
  integer numfiles ! number of input files
  type in_files
    integer       num   ! file number
    character*20  type  ! file type
    character*250 path  ! file name
  end type in_files
  type(in_files), allocatable :: files(:) ! input file information
  character*250 obs_path  ! file name

! cost function variables
  integer nlev_ref   ! (-) number of layers for reference atmosphere profile
  real(8), allocatable :: Tbo_obs(:,:) ! (K) (nchannel,npol) observed brightness temperature at top of atmosphere
  real(8), allocatable :: Tbo_sim(:,:) ! (K) (nchannel,npol) simulated brightness temperature at top of atmosphere
  real(8), allocatable :: cost_chan(:) ! (K) cost function per channel
  real(8) cost_tot                     ! (K) total cost function
  real(8) obs_den                      ! (g/m2) observed hydrometeor column total

! hydrometeor_ph5 variables
  integer phase   ! (-) phase identifier (1=clw 2=rain 3=ice 4=snow 5 =grpl)
  real(8) hab     ! (-) cloud absorption coefficient
  real(8) hsc     ! (-) cloud scattering coefficient
  real(8) g       ! (-) cloud asymytry factor
  real(8) dhab(3) ! (?) derivative absorption wrt to temp, k0, a0
  real(8) dhsc(3) ! (?) derivative scattering wrt to temp, k0, a0
  real(8) dg(3)   ! (?) derivative asymytry factor wrt to temp, k0, a0
  real(8) tair    ! (K) air temperature
  double complex epsil     ! (-) ice dielectric constant
  double complex depsil_dt ! (1/K) derivative of ice dielectric constant wrt temperature

! hydrometeor profile variables
  type hydro_lay_geo
    real(8) bot  ! (km) height of bottom of hydrometeor layer
    real(8) top  ! (km) height of top of hydrometeor layer
    real(8) std  ! (km) height of peak value of hydrometeor layer
    real(8) pk   ! (km) height of peak value of hydrometeor layer
    real(8) tot  ! (g/m2) total column mass of hydrometeor
  end type hydro_lay_geo
  type hydro_layers
    type(hydro_lay_geo) clw    ! (-) cloud liquid water layer
    type(hydro_lay_geo) rain   ! (-) rain layer
    type(hydro_lay_geo) ice    ! (-) ice layer
    type(hydro_lay_geo) snow   ! (-) snow layer
    type(hydro_lay_geo) grpl   ! (-) graupel layer
    type(hydro_lay_geo) cloud  ! (-) cloud layer (clw + ice)
    type(hydro_lay_geo) precip ! (-) precip layer (rain + snow + grpl)
    type(hydro_lay_geo) hydro  ! (-) total hydromet layer (clw + rain + ice + snow + grpl)
  end type hydro_layers
  type(hydro_layers) cur_lay    ! hydrometeor layer variables for current profile
  type(hydro_layers) psuedo_lay ! hydrometeor layer variables for psuedo-data
  type(hydro_layers) init_lay   ! hydrometeor layer variables for initial guess

! planks law
  real(8) Tplank  ! (K) temperature used in plank's law
  real(8) Btotal  ! (W/m2) total flux over wavelength interval
  real(8) wave    ! (m) wavelength
  real(8) Dwave   ! (m) wavelength increment
  real(8) freq    ! (Ghz) frequency
  real(8) Dfreq   ! (1/s) frequency increment
  real(8) Black   ! (w/m2) blackbody emission
  real(8) RH      ! (%) relative humidity
  real(8) kext    ! (1/km) extinction (absorption) coefficient

! begin soil respiration parameters
  real(8) moistexp   ! moist exponent
  real(8) Wsat       ! parameter
  real(8) zm         ! skewness exponent
  real(8) clay       ! clay fraction
  real(8) sand       ! sand fraction
  real(8) soilQ10    ! soil respiration temperature factor
  real(8) Q10        ! soil respiration Q10 value
  real(8) respg      ! ground respiration (mole m-2 m-1)
  real(8) respT      ! total respiration (mole m-2 m-1)
  real(8) respfactor(7) 
  real(8) soilscale(7) 
  real(8) Moist(2)   ! respiration efficiency
  real(8) RconWat    ! (mol/m2/s) contribution of water to respg
  real(8) RconTem    ! (mol/m2/s) contribution of Temperature to respg
  real(8) xroot(7)   ! root fraction per soil layer
  real(8) kroot(12)  ! root density decay constant (1/m)
  logical forcerestore 

! Misc variables
  real(8) poros_om  ! (-) organic matter porosity
  complex(kind=8), allocatable ::  diel(:) ! dielectric constant 
  real(8), allocatable :: len_lay(:)  ! (m) correlation length scales of layer boundaries
  real(8), allocatable :: std_lay(:)  ! (m) standard deviation of layer boundaries
  real(8), allocatable :: dlayer(:)   ! (m) average z coordinate of layer boundaries
  real(8) satfrac       ! Fraction of saturation
  real(8) m_om          ! (kg/m2) mass of organic matter
  real(8) rho_om_max    ! (kg/m3) maximum density of organic matter based on peat
  real(8) dec     ! Solar declination
  real(8) DOY     ! day of year
  real(8) noon        ! GMT of local noon
  real(8) midnight    ! GMT of local midnight
  real(8) sunrise     ! GMT of local rise
  real(8) sunset      ! GMT of local set
  real(8) cosz        ! cosine of solar zenith angle
  real(8) tofday      ! time of day (GMT, in hours)
  real(8) sind        ! function of earth's orbit
  real(8) cosd        ! function of earth's orbit

! scan variables
  real(8) unused     ! (-) unused variable
  real(8) generic    ! generic X variable
  real(8) rhoair     ! (?) density of air
  real(8) Psur       ! (mb) surface total air pressure
  real(8) root_depth ! (m) maximum rooting depth
  real(8) org_depth  ! (m) organic layer depth
  real(8) depth      ! (m) soil depth
  integer len        ! size of array input

  end Module scan_Variables
!
