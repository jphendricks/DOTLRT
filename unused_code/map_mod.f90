!
!---------------------------------------------------
      Module Mapper_Variables
!---------------------------------------------------
! defines all variables used in mapper, input, and output routines
!
      IMPLICIT NONE
!
! begin Input program control variables
      Character *2 year  ! specified year for generating output
!
! Begin Input program option flags
      type Flags
        Integer Map      ! execute mapper subroutine
        Integer Mode     ! 'grid' or 'single' point mode
        Integer Sib      ! generate SiB BC input file
        Integer single   ! generate SiB BC input file ascii format single pt
        Integer EzPlot   ! generate generic ezplot output
        Integer Stats    ! generate output for statistics
        Integer monthly  ! generate monthly output files
        Integer GridKey  ! generate ascci file of lat/lon for each SiB point
        Integer Print    ! print statements to screen
        Integer SoRefTab ! use soil relectance look up table
        Integer SoilMap  ! use %clay/sand or soil type input maps
        Integer SoilProp ! soil properties from look up table or % sand/clay
        Integer fVCov    ! calculate fVCover or use input map
      end type Flags
      type(Flags) Flag
!
! Begin filenames for input lookup tables and maps
      type FileNames
        Character*100 BioMap   ! vegetation type map
        Character*100 BioTab   ! veg. type char. lookup table
        Character*100 MorphTab ! veg. morphilogical lookup table
        Character*100 SoilTab  ! soil type lookup table
        Character*100 SoRefVis ! soil ref. map visible
        Character*100 SoRefNIR ! soil ref map near IR
        Character*100 SoilMap  ! soil characteristic map
        Character*100 AeroVar  ! aerodynamic interpolation tables
        Character*100 PercClay ! percentage soil that is clay map
        Character*100 PercSand ! percentage soil that is sand map
        Character*100 fVCovMap ! fraction veg. cover map
        Character*100 sib_grid ! SiB gridmap file
        Character*100 sib_bc   ! SiB boundary condition file
        Character*100 sib_biom ! GCM biome file
        Character*100 stats    ! generic statistics file
        Character*100 Monthly  ! monthly files
        Character*100 GridKey  ! lat/lon grid key for SiB Points
        Character*100 asciiNDVI! ascii ndvi file for single point option
        Character*100 snowtab  ! snow class table
      end type FileNames
      type(FileNames) FileName    ! file of filenames for input tables and maps
      Character *70 NDVIfiles(100)  ! names for NDVI input files
!
! begin grid parameters
      real lat           ! center latitude of grid cell
      real lon           ! center Longitude of grid cell
      real LatMin        ! Minimum domain latitude (lower left domain corner)
      real LonMin        ! Minimum domain longitude (lower left domain corner)
      real Dlat          ! Latitude Grid Spacing
      real Dlon          ! Longitude grid spacing
      integer imin       ! minimum longitude index for subset of domain grid
      integer imax       ! maximum longitude index for subset of domain grid
      integer jmin       ! minimum latitude index for subset of domain grid
      integer jmax       ! maximum latitude index for subset of domain grid
!
! snow table variables
integer nsnowclass !number of snow classes
type snow_tab
    integer snow_class ! (-) snow class number
    character*20 name  ! (-) name of snow class
    real d_min         ! (m) minimum depth at which bottom layer can form
    real d_max         ! (m) depth at which f_bot_max reaches maximum value
    real f_bot_max     ! (-) maximum bottom layer fraction of snow pack
    real den_ref_top   ! (kg/m3) (kg m-3) reference density for top snow layer
    real den_ref_bot   ! (kg/m3) reference density for bottom snow layer
    real den_bulk_obs  ! (kg/m3) observed snow bulk density
    real den_bulk_std  ! (kg/m3) standard dev of observed snow bulk density
end type snow_tab
type(snow_tab), allocatable :: snowtab(:)
!
! begin input data variables
      integer BioNum     ! biome number read for point from FileName%BioMap
      real, allocatable :: NDVI(:)! array ofFASIR NDVI values for single pt
      real NDVIoffset
      real NDVIscale
!
! begin time dependant, output variables
      type time_dep_var
        real fPAR       ! Canopy absorbed fraction of PAR
        real LAI        ! Leaf-area index
        real Green      ! Canopy greeness fraction of LAI
        real zo         ! Canopy roughness coeff 
        real zp_disp    ! Zero plane displacement
        real RbC        ! RB Coefficient (c1)
        real RdC        ! RC Coefficient (c2)
        real gmudmu     ! Time-mean leaf projection
      end type time_dep_var
!
      type(time_dep_var), allocatable :: TimeVar(:)
!     TimeVar(time) ! time dependant variables
!
! begin aerodynamic interpolation tables
      real LAIgrid(50)   ! grid of LAI values for lookup table
      real fVCovergrid(50)! grid of fVCover values for interpolation table
!
      type aero_var
        real zo         ! Canopy roughness coeff 
        real zp_disp    ! Zero plane displacement
        real RbC        ! RB Coefficient
        real RdC        ! RC Coefficient
      end type aero_var
      type(aero_var), allocatable :: AeroVar(:,:,:)
!     AeroVar(Biome,LAI,fVCover) ! interpolation tables for aero variables
!
! begin Biome-dependent variables
! same structure applies to input lookup tables and output files
      type Biome_dep_var
        integer bioNum  ! biome or vegetation cover type
        real z2         ! Canopy top height (m)
        real z1         ! Canopy base height (m)
        real fVCover    ! Canopy cover fraction
        real ChiL       ! Leaf angle distribution factor
        real SoDep      ! Total depth of 3 soil layers (m)
        real RootD      ! Rooting depth (m)
        real Phi_half   ! 1/2 Critical leaf water potential limit (m)
        real LTran(2,2) ! Leaf transmittance for green/brown plants
        real LRef(2,2)  ! Leaf reflectance for green/brown plants
!                      For LTran and LRef:
!                      (1,1)=shortwave, green plants
!                      (2,1)=longwave, green plants
!                      (1,2)=shortwave, brown plants
!                      (2,2)=longwave, brown plants
        real vmax0      ! Rubisco velocity of sun-leaf (Mol m^-2 s^-1)
        real EffCon     ! Quantum efficiency (Mol Mol^-1)
        real gsSlope    ! Ball-Berry Conductance Slope Parameter
        real gsMin      ! Ball-Berry Conductance Intercept Parameter
        real Atheta     ! WC WE Coupling Parameter
        real Btheta     ! WC & WE, WS Coupling Parameter
        real TRDA       ! Temperature Coefficient in GS-A Model (K^-1)
        real TRDM       ! "" (K)
        real TROP       ! "" (K)
        real respcp     ! Respiration Fraction of Vmax
        real SLTI       ! Slope of low-temp inhibition function (K^-1)
        real HLTI       ! Slope of high-temp inhibition function (K^-1)
        real SHTI       ! 1/2 Point of low-temp inhibition function (K)
        real HHTI       ! 1/2 Point of high-temp inhibition function (K)
        real SoRef(2)   ! 2-stream soil and litter reflectivity
!                       Soref(1)=visible soil and litter reflectivity
!                       Soref(2)=Near IR soil and litter reflectivity
        real tfrost     ! temp below which frost reduces photosynthesis (K)
      end type Biome_dep_var
!
      type(Biome_dep_var) BioVar ! time ind., biome dependant variables
!
      type(Biome_dep_var), allocatable :: BioTab(:)
!     BioTab(Biome)   ! Lookup table of biome dependant variables
!
      type biome_morph_var
        real zc         ! Canopy inflection height (m)
        real LWidth     ! Leaf width
        real LLength    ! Leaf length
        real LAImax     ! Maximum LAI
        real stems      ! Stem area index
        real NDVImax    ! Maximum NDVI
        real NDVImin    ! Minimum NDVI
        real SRmax      ! Maximum simple ratio
        real SRmin      ! Minimum simple ratio
      end type biome_morph_var
!
      type(biome_morph_var), allocatable :: MorphTab(:)
!     MorphTab(Biome)   ! Lookup table of biome dependant morphology
!
! begin soil dependant variables
      integer SoilNum   ! soil type number
      type soil_Physical
        integer SoilNum ! soil type number
        real BEE     ! Soil wetness exponent
        real PhiSat  ! Soil tension at saturation
        real SatCo   ! Hydraulic conductivity at saturation
        real poros   ! Soil porosity
        real Slope   ! Cosine of mean slope
        real Wopt    ! optimal soil wetness for soil respiration
        real Skew    ! skewness exponent of soil respiration vs. wetness curve
        real RespSat ! assures soil respiration is 60-80% of max @ saturation
      end type soil_Physical
!
      type(soil_Physical) SoilVar  ! time ind., soil dependant variables
      type(soil_Physical), allocatable :: SoilTab(:)
!     SoilTab(nsoil)  ! Lookup table of Soil dependant variables
!
! begin Soil texture variables (Only clay and sand need to be specified)
      type soil_texture
        real clay       ! Percent clay content
        real silt       ! Percent silt content
        real sand       ! Percent sand content
        integer class   ! Soil texture class
      end type soil_texture
      type(soil_texture) text ! soil texture: percent sand, silt, and clay
!
! begin Input program control variables
      integer nm         ! total number of months in a year
      real DOYstart      ! Day of Year (DOY) first ndvi input map in time series
      real DDOY          ! Delta DOY between ndvi input maps in time series
      integer nv         ! total number of possible vegetation types
      integer ns         ! total number of possible soil types
      Integer minBiome   ! minimum biome number for biome subset option
      Integer maxBiome   ! maximum biome number for biome subset option
      integer imm        ! total number of longitude points
      integer jmm        ! total number of latitude points
!
      end Module Mapper_Variables
!
