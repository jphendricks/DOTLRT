!
!=======================================================================
      subroutine mapper( &
      lat, &
      doy, &
      prevNDVI, &
      curNDVI, &
      BioVar, &
      MorphTab, &
      TimeVar)
!=======================================================================
! calculates time dependant boundary condition variables for SiB.
!
IMPLICIT NONE
!
! begin input variables
      real lat         ! center latitude of grid cell
      real prevNDVI  ! curent NDVI values for a grid cell
      real curNDVI  ! curent NDVI values for a grid cell
!
! begin input Biome-dependent variables
      type Biome_dep_var
        real biome   ! biome or vegetation cover type
        real z2         ! Canopy base height (meters)
        real z1         ! Canopy top height (m)
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
        real gsSlope    ! Conductance-Photosynthesis Slope Parameter
        real gsMin      ! Conductance-Photosynthesis Intercept
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
        real SoRef(2)   ! 2-stream soil/litter ref
      end type Biome_dep_var
!
      type(Biome_dep_var) BioVar
!
! begin input biome dependant, physical morphology variables
      type biome_morph_var
        real zc        ! Canopy inflection height (m)
        real LWidth    ! Leaf width
        real LLength   ! Leaf length
        real LAImax    ! Maximum LAI
        real stems     ! Stem area index
        real NDVImax   ! Maximum NDVI
        real NDVImin   ! Minimum NDVI
        real SRmax     ! Maximum simple ratio
        real SRmin     ! Minimum simple ratio
      end type biome_morph_var
      type(biome_morph_var) MorphTab
!
! begin time dependant, output variables
      type time_dep_var
        real fPAR    ! Canopy absorbed fraction of PAR
        real LAI     ! Leaf-area index
        real Green   ! Canopy greeness fraction of LAI
        real zo      ! Canopy roughness coeff 
        real zp_disp ! Zero plane displacement
        real RbC     ! RB Coefficient (c1)
        real RdC     ! RC Coefficient (c2)
        real gmudmu  ! Time-mean leaf projection
      end type time_dep_var
      type(time_dep_var) TimeVar
!
! begin internal variables
      real DOY         ! Day of Year (DOY) of ndvi input map
      real prevfpar
      real, parameter :: fPARmax=0.95
!                      ! Maximum possible FPAR corresponding to 98th percentile
      real, parameter :: fPARmin=0.01
!                      ! Minimum possible FPAR corresponding to 2nd percentile
!     For more information on fPARmin and fPARmax, see
!     Sellers et al. (1994a, pg. 3532); Los (1998, pg. 29, 37-39)
!
!-----------------------------------------------------------------------
! Calculate time dependant variables
!-----------------------------------------------------------------------
! Calculate first guess fPAR (ave Simple Ratio (SR) and NDVI methods)

      call AverageAPAR (curNDVI,    &
         MorphTab%NDVImin, MorphTab%NDVImax,  &
         MorphTab%SRmin, MorphTab%SRmax, &
         fPARmax, fParmin, TimeVar%fPAR)

      call AverageAPAR (prevNDVI,    &
         MorphTab%NDVImin, MorphTab%NDVImax,  &
         MorphTab%SRmin, MorphTab%SRmax, &
         fPARmax, fParmin, prevfPAR)
!
! Calculate leaf area index (LAI) and greeness fraction (Green)
!   See S. Los et al 1998 section 4.2.
      call laigrn (TimeVar%fPAR, &
                   prevfPAR, &
                   fPARmax, &
                   BioVar%fVCover, &
                   MorphTab%stems, &
                   MorphTab%LAImax, &
                   TimeVar%Green, &
                   TimeVar%LAI)
!
! Interpolate to calculate aerodynamic, time varying variables
!      call AeroInterpolate (
!     &              TimeVar%LAI, 
!     &              BioVar%fVCover, 
!     &              LAIgrid,
!     &              fVCovergrid,
!     &              AeroVar,
!     &              TimeVar%zo,
!     &              TimeVar%zp_disp,
!     &              TimeVar%RbC,
!     &              TimeVar%RdC)
!
! Calculate mean leaf orientation to par flux (gmudmu)
      call Oldgmuder (lat, &
                    DOY, &
                    BioVar%ChiL, &
                    TimeVar%gmudmu)
!
! recalculate fPAR adjusting for Sun angle, vegetation cover fraction,
! and greeness fraction, and LAI
      call aparnew (TimeVar%LAI, &
                   TimeVar%Green, &
                   BioVar%LTran, &
                   BioVar%LRef, &
                   TimeVar%gmudmu, &
                   BioVar%fVCover, &
                   TimeVar%fPAR, &
                   fPARmax, &
                   fPARmin)
!
        return
        end
!
!=======================================================================
      subroutine laigrn (fPAR,fPARm,fPARmax,fVCover,stems, &
      LAImax,Green,LAI)
!=======================================================================
! calculate leaf area index (LAI) and greenness fraction (Green) from fPAR. 
! LAI is linear with vegetation fraction and exponential with fPAR.
! See Sellers et al (1994), Equations 7 through 13.
!                                                                               
      implicit none
!
! begin input variables
      real fPAR     ! fraction of PAR absorbed by plants at current time
      real fPARm    ! fraction of PAR absorbed by plants at previous time
      real fPARmax  ! maximum possible FPAR corresponding to 98th percentile
      real fVCover  ! vegetation cover fraction
      real stems    ! stem area index for the specific biome type
      real LAImax   ! maximum total leaf area index for specific biome type
!
! begin output variables
      real Green    ! greeness fraction of the total leaf area index
      real LAI      ! area average total leaf area index
!
! begin internal variables
      real LAIg     ! green leaf area index at current time
      real LAIgm    ! green leaf area index at previous time
      real LAId     ! dead leaf area index at current time
!
! Calculate current and previous green leaf area index (LAIg and LAIgm):
! LAIg is log-linear with fPAR.  Since measured fPAR is an area average, 
! divide by fVCover to get portion due to vegetation.  Since fVCover can
! be specified, check to assure that calculated fPAR does not exceed fPARMax.
!
      if(fPAR/fVCover.ge.fPARmax) then
        LAIg=LAImax
      else
        LAIg=alog(1.-fPAR/fVCover)*LAImax/alog(1-fPARmax)
      endif
!
      if(fPARm/fVCover.ge.fPARmax) then
        LAIgm=LAImax
      else
        LAIgm=alog(1.-fPARm/fVCover)*LAImax/alog(1-fPARmax)
      endif
!
! Calculate dead leaf area index (LAId):
! If LAIg is increasing or unchanged, the vegetation is in growth mode.
! LAId is then very small (very little dead matter).
! If LAIg is decreasing, the peak in vegetation growth has passed and
! leaves have begun to die off.  LAId is then half the change in LAIg,
! assuming half the dead leaves fall off.
!
!     Growth mode dead leaf area index:
      if (LAIg.ge.LAIgm) LAId=0.0001
!
!     die-off (post peak growth) dead leaf area index:
      if (LAIg.lt.LAIgm) LAId=0.5*(LAIgm-LAIg)
!
! Calculate area average, total leaf area index (LAI):
      LAI=(LAIg+LAId+stems)*fVCover
!
! Calculate greeness fraction (Green):
! Greeness fraction=(green leaf area index)/(total leaf area index)
      Green=LAIg/(LAIg+LAId+stems)
!      print*, 'fPARm=',fPARm,'fPAR=',fPAR,'LAId=',LAId,'green=', green
!
      return                                                                    
      end
!
!=======================================================================      
      subroutine Oldgmuder (Lat, DOY, ChiL, gmudmu)
!=======================================================================      
! calculates time mean leaf projection relative to the Sun.
!
      implicit none
!
! begin input variables
      real Lat      ! latitude in degrees
      real DOY      ! day-of-year (typically middle day of the month)
      real ChiL     ! leaf angle distribution factor
!
! begin output variables
      real gmudmu   ! time mean projected leaf area normal to Sun
!
! begin internal variables
      integer itime ! time counter
      real time     ! time from 0:00 Greenwhich Mean Time (GMT)
      real coshr    ! cosine of the Greenwhich Meridian (GM) Hour angle
      real cosz     ! cosine of the Solar zenith angle (mu)
      real chiv     ! dummy variable for leaf angle distribution factor
      real dec      ! declination of the Sun (Solar Declination)
      real sind     ! sine of the solar declination
      real cosd     ! cosine of the solar declination
      real pi       ! the constant pi
      real pi180    ! conversion factor from degrees to radians
      real aa       ! minimum possible LAI projection vs. cosine Sun angle
      real bb       ! slope leaf area projection vs. cosine Sun angle
      real cloud    ! Cloud cover fraction
      real coszbar  ! PAR weighted mean cosine of Sun Angle
      real swdown   ! downward shortwave flux onto canopy
      real pardif   ! Diffuse visible light onto canopy
      real pardir   ! Direct visible light onto canopy
      real difrat   ! diffuse fraction of SW light
      real vnrat    ! visible fraction of SW light
      real tor      ! (?) TBD
      real PARcoszTot ! Total cosz weighted PAR flux onto canopy during 24 hr
      real PARTot   ! total PAR flux onto canopy during 24 hour period
!
! Assign values to constants
      data pi /3.141592/
      pi180=pi/180.
!
! Assume a constant cloud fraction
      cloud=0.5
!
! Calculate solar declination in radians
      dec=pi180*23.5*sin(1.72e-2*(DOY-80.))
!
! Calculate sine and cosine of solar declination
      sind=sin(dec)                                                         
      cosd=cos(dec)
!
! PAR weighted mean cosz over 24 hour period.  Calculations do not
! match the G(mu)/mu described in Bonan (1996) and Sellers (1985).
      PARcoszTot=0.
      PARTot=0.
!
! Scan 24 hours in half hour increments
      do itime=1, 48, 1                                                     
!       time from zero Greenwhich Mean Time (GMT)
        time=0.5*real(itime) 
!
!       cosine of hour angle of Grenwhich Meridion (GM)
        coshr=cos(-pi+time/24.*2.*pi)
!
!       Cosine of the Sun angle (lon=GM=0 deg, latitude=Lat)
!       cosz=0.01 during night
        cosz=sin(Lat*pi180)*sind+cos(Lat*pi180)*cosd*coshr
        cosz=amax1(0.01, cosz) 
!
!       shortwave downward flux onto canopy based on solar constant
        tor=0.7**(1./cosz)
        swdown=1375.*cosz*(tor+0.271-0.294*tor)
!
!       diffuse fraction of shortwave light
        difrat=0.0604/(cosz-0.0223)+0.0683
        difrat=max(difrat,0.)
        difrat=min(difrat,1.)
        difrat=difrat+(1.-difrat)*cloud
!
!       visible fraction of shortwave flux
        vnrat=(580.-cloud*464.)/((580.-cloud*499.) &
       +(580.-cloud*464.))
!
!       direct and diffuse visible/PAR flux onto canopy
        pardir=(1.-difrat)*vnrat*swdown
        pardif=difrat*vnrat*swdown
!
!       Total cosz weighted PAR flux
        PARcoszTot=PARcoszTot+pardir*cosz+pardif*0.5
!
!       Total PAR flux
        PARTot=PARTot+pardir+pardif
      Enddo                                                                 
!
! PAR weighted mean cosz
      coszbar=PARcoszTot/PARTot
!
! LAI projection (G(mu)
      chiv=ChiL                                                               
      if (abs(chiv) .le. 0.01) chiv=0.01
!     minimum value of projected leaf area
      aa=0.5-0.633*chiv-0.33*chiv*chiv
!     slope of projected leaf area wrt cosine sun angle
      bb=0.877*(1.-2.*aa)                                             
!
! mean projected LAI in Sun direction using PAR weighted mean cosz
      gmudmu=(aa+bb*coszbar)/coszbar
!                                                                               
      return                                                                    
      end
!
!=======================================================================
      subroutine AeroInterpolate (LAI, fVCover,  &
      LAIgrid,fVCovergrid,AeroVar,zo,zp_disp,RbC,RdC)
!=======================================================================
! This subroutine calculates the aerodynamic parameters by bi-linear 
! interpolation from a lookup table of previously calculated values.  
! The interpolation table is a numpts x numpts LAI/fVCover grid with
! LAI ranging from 0.02 to 10 and fVCover ranging from 0.01 to 1.
!
      implicit none
!
! begin input variables
      real LAI            ! actual area averaged LAI for interpolation
      real fVCover        ! vegetation cover fraction for interpolation
      real LAIgrid(50)    ! grid of LAI values for lookup table 
      real fVCovergrid(50)! grid of fVCover values for interpolation table
      type aero_var
        real zo           ! Canopy roughness coeff 
        real zp_disp      ! Zero plane displacement
        real RbC          ! RB Coefficient
        real RdC          ! RC Coefficient
      end type aero_var
      type(aero_var) AeroVar(50,50) ! interpolation tables
!
! begin output variables
      real RbC            ! interpolated Rb coefficient
      real RdC            ! interpolated Rd coefficient
      real zo             ! interpolated roughness length
      real zp_disp        ! interpolated zero plane displacement
!
! begin internal variables
      integer i           ! index for LAI grid location
      integer j           ! index for fVCover grid location
      real LocLAI         ! local LAI var. to prevent changing main LAI value
      real LocfVCover     ! local fVCover var. to prevent changing fVCover value
      real DLAI           ! grid spacing between LAI values in tables
      real DfVCover       ! grid spacing between fVCover values in tables
!
! calculate grid spacing (assumed fixed)
      DLAI=LAIgrid(2)-LAIgrid(1)
      DfVCover=fVCovergrid(2)-fVCovergrid(1)
!
! Assign input LAI and fVCover to local variables and make sure
! they lie within the limits of the interpolation tables, assuring 
! the LAI and fVCover values returned from the subroutine are not modified.
      LocLAI=max(LAI,0.02)
      LocfVCover=max(fVCover,0.01)
!
! determine the nearest array location for the desired LAI and fVCover
      i=int(LocLAI/DLAI+1)
      j=int(LocfVCover/DfVCover+1)
      j=min(j,49)
!
! interpolate RbC variable
         call interpolate( &
      LAIgrid(i),  &
      LocLAI,  &
      DLAI, &
      fVCovergrid(j),  &
      LocfVCover,  &
      DfVCover, &
      AeroVar(i,j)%RbC,  &
      AeroVar(i+1,j)%RbC,  &
      AeroVar(i,j+1)%RbC,  &
      AeroVar(i+1,j+1)%RbC,  &
      RbC)
!
! interpolate RdC variable
         call interpolate( &
      LAIgrid(i),  &
      LocLAI,  &
      DLAI, &
      fVCovergrid(j),  &
      LocfVCover,  &
      DfVCover, &
      AeroVar(i,j)%RdC,  &
      AeroVar(i+1,j)%RdC,  &
      AeroVar(i,j+1)%RdC,  &
      AeroVar(i+1,j+1)%RdC,  &
      RdC)
!
! interpolate roughness length'
         call interpolate( &
      LAIgrid(i),  &
      LocLAI,  &
      DLAI, &
      fVCovergrid(j),  &
      LocfVCover,  &
      DfVCover, &
      AeroVar(i,j)%zo,  &
      AeroVar(i+1,j)%zo,  &
      AeroVar(i,j+1)%zo,  &
      AeroVar(i+1,j+1)%zo,  &
      zo)
!
! interpolate zero plane displacement
         call interpolate( &
      LAIgrid(i),  &
      LocLAI,  &
      DLAI, &
      fVCovergrid(j),  &
      LocfVCover,  &
      DfVCover, &
      AeroVar(i,j)%zp_disp,  &
      AeroVar(i+1,j)%zp_disp,  &
      AeroVar(i,j+1)%zp_disp,  &
      AeroVar(i+1,j+1)%zp_disp,  &
      zp_disp)
!
      return
      end
!
!=======================================================================
      subroutine NewAeroInterpolate (LAI, fVCover,  &
      LAIgrid,fVCovergrid,AeroTab,AeroVar)
!=======================================================================
! calculates the aerodynamic parameters by bi-linear 
! interpolation from a lookup table of previously calculated values.  
! The interpolation table is a numpts x numpts LAI/fVCover grid with
! LAI ranging from 0.02 to 10 and fVCover ranging from 0.01 to 1.
!
      implicit none 
!
! begin input variables
      real LAI            ! actual area averaged LAI for interpolation
      real fVCover        ! vegetation cover fraction for interpolation
      real LAIgrid(50)    ! grid of LAI values for lookup table 
      real fVCovergrid(50)! grid of fVCover values for interpolation table
      type aero_var
        real :: zo      ! Canopy roughness coeff 
        real :: zp_disp ! Zero plane displacement
        real :: RbC     ! RB Coefficient (c1 or cc1)
        real :: RdC     ! RC Coefficient (c2 or cc2)
        real :: G2      ! Ratio Ra (actual) to Ra (log-linear) for momentum
        real :: G3      ! Ratio Ra (actual) to Ra (log-linear) for heat transfer
        real :: CORB1   ! Non-nuetral correction for Ra between Ha and z2
        real :: CORB2   ! neutral value of RBB*U2^2 (RdC^2 for upper canopy)
        real :: HA      ! Canopy source height 
      end type aero_var
      type(aero_var) AeroTab(50,50)
!
! begin output variables
      type(aero_var) AeroVar
!
! begin internal variables
      integer i           ! index for LAI grid location
      integer j           ! index for fVCover grid location
      real LocLAI         ! local LAI var. to prevent changing main LAI value
      real LocfVCover     ! local fVCover var. to prevent changing fVCover value
      real DLAI           ! grid spacing between LAI values in tables
      real DfVCover       ! grid spacing between fVCover values in tables
!
! calculate grid spacing (assumed fixed)
      DLAI=LAIgrid(2)-LAIgrid(1)
      DfVCover=fVCovergrid(2)-fVCovergrid(1)
!
! Assign input LAI and fVCover to local variables and make sure
! they lie within the limits of the interpolation tables, assuring 
! the LAI and fVCover values returned from the subroutine are not modified.
      LocLAI=max(LAI,0.02)
      LocfVCover=max(fVCover,0.01)
!
! determine the nearest array location for the desired LAI and fVCover
      i=int(LocLAI/DLAI+1)
      j=int(LocfVCover/DfVCover+1)
      j=min(j,49)
!
! interpolate RbC variable
         call interpolate( &
      LAIgrid(i),  &
      LocLAI,  &
      DLAI, &
      fVCovergrid(j),  &
      LocfVCover,  &
      DfVCover, &
      AeroTab(i,j)%RbC,  &
      AeroTab(i+1,j)%RbC,  &
      AeroTab(i,j+1)%RbC,  &
      AeroTab(i+1,j+1)%RbC,  &
      AeroVar%RbC)

!
! interpolate RdC variable
         call interpolate( &
      LAIgrid(i),  &
      LocLAI,  &
      DLAI, &
      fVCovergrid(j),  &
      LocfVCover,  &
      DfVCover, &
      AeroTab(i,j)%RdC,  &
      AeroTab(i+1,j)%RdC,  &
      AeroTab(i,j+1)%RdC,  &
      AeroTab(i+1,j+1)%RdC,  &
      AeroVar%RdC)
!
! interpolate roughness length (z0)
         call interpolate( &
      LAIgrid(i),  &
      LocLAI,  &
      DLAI, &
      fVCovergrid(j),  &
      LocfVCover,  &
      DfVCover, &
      AeroTab(i,j)%zo,  &
      AeroTab(i+1,j)%zo,  &
      AeroTab(i,j+1)%zo,  &
      AeroTab(i+1,j+1)%zo,  &
      AeroVar%zo)
!
! interpolate zero plane displacement (zp_disp)
         call interpolate( &
      LAIgrid(i),  &
      LocLAI,  &
      DLAI, &
      fVCovergrid(j),  &
      LocfVCover,  &
      DfVCover, &
      AeroTab(i,j)%zp_disp,  &
      AeroTab(i+1,j)%zp_disp,  &
      AeroTab(i,j+1)%zp_disp,  &
      AeroTab(i+1,j+1)%zp_disp,  &
      AeroVar%zp_disp)
!
! interpolate Ra ratio for momentum (G2)
         call interpolate( &
      LAIgrid(i),  &
      LocLAI,  &
      DLAI, &
      fVCovergrid(j),  &
      LocfVCover,  &
      DfVCover, &
      AeroTab(i,j)%G2,  &
      AeroTab(i+1,j)%G2,  &
      AeroTab(i,j+1)%G2,  &
      AeroTab(i+1,j+1)%G2,  &
      AeroVar%G2)
!
! interpolate Ra ratio for Heat transfer (G3)
         call interpolate( &
      LAIgrid(i),  &
      LocLAI,  &
      DLAI, &
      fVCovergrid(j),  &
      LocfVCover,  &
      DfVCover, &
      AeroTab(i,j)%G3,  &
      AeroTab(i+1,j)%G3,  &
      AeroTab(i,j+1)%G3,  &
      AeroTab(i+1,j+1)%G3,  &
      AeroVar%G3)
!
! interpolate Non-nuetral correction for Ra between Ha and z2 (Corb1)
         call interpolate( &
      LAIgrid(i),  &
      LocLAI,  &
      DLAI, &
      fVCovergrid(j),  &
      LocfVCover,  &
      DfVCover, &
      AeroTab(i,j)%Corb1,  &
      AeroTab(i+1,j)%Corb1,  &
      AeroTab(i,j+1)%Corb1,  &
      AeroTab(i+1,j+1)%Corb1,  &
      AeroVar%Corb1)
!
! interpolate neutral value of RBB*U2^2 (Corb2)
         call interpolate( &
      LAIgrid(i),  &
      LocLAI,  &
      DLAI, &
      fVCovergrid(j),  &
      LocfVCover,  &
      DfVCover, &
      AeroTab(i,j)%Corb2,  &
      AeroTab(i+1,j)%Corb2,  &
      AeroTab(i,j+1)%Corb2,  &
      AeroTab(i+1,j+1)%Corb2,  &
      AeroVar%Corb2)
!
! interpolate Canopy source height (HA)
         call interpolate( &
      LAIgrid(i),  &
      LocLAI,  &
      DLAI, &
      fVCovergrid(j),  &
      LocfVCover,  &
      DfVCover, &
      AeroTab(i,j)%HA,  &
      AeroTab(i+1,j)%HA,  &
      AeroTab(i,j+1)%HA,  &
      AeroTab(i+1,j+1)%HA,  &
      AeroVar%HA)
!
      return
      end
!
!=======================================================================
      subroutine interpolate(x1, x, Dx, &
      y1, y, Dy, &
      z11, z21, z12, z22, z)
!=======================================================================
! calculates the value of z=f(x,y) by linearly interpolating
! between the 4 closest data points on a uniform grid.  The subroutine
! requires a grid point (x1, y1), the grid spacing (Dx and Dy), and the 
! 4 closest data points (z11, z21, z12, and z22).
!
! begin input variables
      real x1  ! the x grid location of z11
      real x   ! x-value at which you will interpolate z=f(x,y)
      real Dx  ! grid spacing in the x direction
      real y1  ! the y grid location of z11
      real y   ! y-value at which you will interpolate z=f(x,y)
      real Dy  ! grid spacing in the y direction
      real z11 ! f(x1, y1)
      real z21 ! f(x1+Dx, y1)
      real z12 ! f(x1, y1+Dy)
      real z22 ! f(x1+Dx, y1+Dy)
!
! begin output variables
      real z   ! f(x,y), the desired interpolated value
!
! begin internal variables
      real zp  ! z'=first interpolated value at (x, y1)
      real zpp ! z''=second interpolated value at (x, Y1+Dy)
!
! interpolate between z11 and z21 to calculate z' (zp) at (x, y1)
      zp=z11+(x-x1)*(z21-z11)/Dx
!
! interpolate between z12 and z22 to calculate z'' (zpp) at (x, Y1+Dy)
      zpp=z12+(x-x1)*(z22-z12)/Dx
!
! interpolate between zp and zpp to calculate z at (x,y)
      z=zp+(y-y1)*(zpp-zp)/Dy
!
      return
      end
!
!=======================================================================
      subroutine srapar (ndvi, SRmin, SRmax, fPAR, fPARmax, fParmin)
!=======================================================================      
! calculates Canopy absorbed fraction of Photosynthetically 
! Active Radiation (fPAR) using the Simple Ratio (sr) method 
! (Los et al. (1998), eqn 6). This empirical method assumes a linear
! relationship between fPAR and sr.
!----------------------------------------------------------------------
!
      implicit none
!
! begin input variables
      real ndvi    ! normalized difference vegetation index
      real SRmin   ! minimum simple ratio for vegetation type
      real SRmax   ! maximum simple ratio for vegetation type
      real fPARmax ! Maximum possible FPAR corresponding to 98th percentile
      real fPARmin ! Minimum possible FPAR corresponding to 2nd percentile
!
! begin output variables
      real fPAR    ! Canopy absorbed fraction of PAR
!
! begin internal variables
      real sr      ! simple ratio of near IR and visible radiances
!
! Calculate simple ratio (SR)
      sr=(1.+ndvi)/(1.-ndvi)
!
! Insure calculated SR value falls within physical limits for veg. type
      sr=max(sr,SRmin)
      sr=min(sr,SRmax)
!
! Calculate fPAR using SR method (Los et al. (1998), eqn 6)
      fPAR=(sr-SRmin)*(fPARmax-fPARmin)/(SRmax-SRmin)+fPARmin
!
      return                                                                    
      end
!
!=======================================================================
      subroutine NDVIapar (ndvi, NDVImin, NDVImax,  &
       fPAR, fPARmax, fParmin)
!=======================================================================      
! calculates Canopy absorbed fraction of Photosynthetically 
! Active Radiation (fPAR) using the NDVI method 
! (Los et al. (1998), eqn 7). This empirical method assumes a linear
! relationship between fPAR and NDVI.
!----------------------------------------------------------------------
!
      implicit none
!
! begin input variables
      real ndvi    ! normalized difference vegetation index
      real NDVImin ! minimum NDVI for vegetation type
      real NDVImax ! maximum NDVI for vegetation type
      real fPARmax ! Maximum possible FPAR corresponding to 98th percentile
      real fPARmin ! Minimum possible FPAR corresponding to 2nd percentile
!
! begin output variables
      real fPAR    ! Canopy absorbed fraction of PAR
!
! begin internal variables
      real LocNDVI ! local value of NDVI to prevent changes in input value
!
! switch to local value of ndvi to prevent any changes going back to main
      LocNDVI=NDVI
!
! Insure calculated NDVI value falls within physical limits for veg. type
      LocNDVI=max(LocNDVI,NDVImin)
      LocNDVI=min(LocNDVI,NDVImax)
!
! Calculate fPAR using NDVI method (Los et al. (1998), eqn 6)
      fPAR=(LocNDVI-NDVImin)*(fPARmax-fPARmin)/ &
          (NDVImax-NDVImin)+fPARmin
!
      return                                                                    
      end
!
!=======================================================================
      subroutine AverageAPAR (ndvi, NDVImin, NDVImax, SRmin, SRmax, &
       fPARmax, fParmin, fPAR)
!=======================================================================      
! calculates Canopy absorbed fraction of Photosynthetically 
! Active Radiation (fPAR) using an average of the Simple Ratio (sr) 
! and NDVI methods (Los et al. (1999), eqn 5-6).  The empirical
! SR method assumes a linear relationship between fPAR and SR.
! The NDVI method assumes a linear relationship between fPAR and NDVI.
!----------------------------------------------------------------------
!
      implicit none
!
! begin input variables
      real ndvi     ! normalized difference vegetation index
      real NDVImin  ! minimum NDVI for vegetation type
      real NDVImax  ! maximum NDVI for vegetation type
      real SRmin    ! minimum NDVI for vegetation type
      real SRmax    ! maximum NDVI for vegetation type
      real fPARmax  ! Maximum possible FPAR corresponding to 98th percentile
      real fPARmin  ! Minimum possible FPAR corresponding to 2nd percentile
!
! begin output variables
      real fPAR     ! Canopy absorbed fraction of PAR
!
! begin internal variables
      real LocNDVI  ! local value of NDVI to prevent changes in input value
      real sr       ! simple ratio of near IR and visible radiances
      real NDVIfPAR ! fPAR from NDVI method
      real SRfPAR   ! fPAR from SR method
!
! switch to local value of ndvi to prevent any changes going back to main
      LocNDVI=NDVI
!
! Insure calculated NDVI value falls within physical limits for veg. type
      LocNDVI=max(LocNDVI,NDVImin)
      LocNDVI=min(LocNDVI,NDVImax)
!
! Calculate simple ratio (SR)
      sr=(1.+LocNDVI)/(1.-LocNDVI)
!
! Calculate fPAR using SR method (Los et al. (1999), eqn 5)
      SRfPAR=(sr-SRmin)*(fPARmax-fPARmin)/(SRmax-SRmin)+fPARmin
!
! Calculate fPAR using NDVI method (Los et al. (1999), eqn 6)
      NDVIfPAR=(LocNDVI-NDVImin)*(fPARmax-fPARmin)/ &
          (NDVImax-NDVImin)+fPARmin
!
! take average of two methods
      fPAR=0.5*(SRfPAR+NDVIfPAR)
!
      return                                                                    
      end
!
!=======================================================================
      subroutine testfPAR (ndvi,NDVImin,NDVImax,SRmin,SRmax,fPARmax,fParmin,fPAR,test,test10)
!=======================================================================      
! calculates Canopy absorbed fraction of Photosynthetically 
! Active Radiation (fPAR) using an average of the Simple Ratio (sr) 
! and NDVI methods (Los et al. (1999), eqn 5-6).  The empirical
! SR method assumes a linear relationship between fPAR and SR.
! The NDVI method assumes a linear relationship between fPAR and NDVI.
!----------------------------------------------------------------------
!
      implicit none
!
! begin input variables
      real ndvi     ! normalized difference vegetation index
      real NDVImin  ! minimum NDVI for vegetation type
      real NDVImax  ! maximum NDVI for vegetation type
      real SRmin    ! minimum NDVI for vegetation type
      real SRmax    ! maximum NDVI for vegetation type
      real fPARmax  ! Maximum possible FPAR corresponding to 98th percentile
      real fPARmin  ! Minimum possible FPAR corresponding to 2nd percentile
!
! begin output variables
      real fPAR     ! Canopy absorbed fraction of PAR
      real test, test10(10)   ! test variables
!
! begin internal variables
      real LocNDVI  ! local value of NDVI to prevent changes in input value
      real sr       ! simple ratio of near IR and visible radiances
      real NDVIfPAR ! fPAR from NDVI method
      real SRfPAR   ! fPAR from SR method
!
! switch to local value of ndvi to prevent any changes going back to main
      LocNDVI=NDVI
!
! Insure calculated NDVI value falls within physical limits for veg. type
      LocNDVI=max(LocNDVI,NDVImin)
      LocNDVI=min(LocNDVI,NDVImax)
!
! Calculate simple ratio (SR)
      sr=(1.+LocNDVI)/(1.-LocNDVI)
!
! Calculate fPAR using SR method (Los et al. (1999), eqn 5)
      SRfPAR=(sr-SRmin)*(fPARmax-fPARmin)/(SRmax-SRmin)+fPARmin
!
! Calculate fPAR using NDVI method (Los et al. (1999), eqn 6)
      NDVIfPAR=(LocNDVI-NDVImin)*(fPARmax-fPARmin)/(NDVImax-NDVImin)+fPARmin
!
! take average of two methods
      fPAR=0.5*(SRfPAR+NDVIfPAR)
!
      test=fpar
      test10(1)=fpar
      return                                                                    
      end
!
!=======================================================================        
      subroutine aparnew (LAI,Green,LTran,LRef,gmudmu,fVCover, &
      fPAR, fPARmax, fPARmin)
!=======================================================================
! recomputes the Canopy absorbed fraction of Photosynthetically
! Active Radiation (fPAR), adjusting for solar zenith angle and the 
! vegetation cover fraction (fVCover) using a modified form of Beer's law.
! See Sellers et al. Part II (1996), eqns. 9-13.
!
      implicit none
!
! begin input variables
      real LAI       ! Leaf Area Index
      real Green     ! Greeness fraction of Leaf Area Index
      real LTran(2,2)! Leaf transmittance for green/brown plants
      real LRef(2,2) ! Leaf reflectance for green/brown plants
!                      For LTran and LRef:
!                        (1,1)=shortwave, green plants
!                        (2,1)=longwave, green plants
!                        (1,2)=shortwave, brown plants
!                        (2,2)=longwave, brown plants
      real gmudmu    ! daily Time-mean canopy optical depth
      real fVCover   ! Canopy cover fraction
      real fPARmax   ! Maximum possible FPAR corresponding to 98th percentile
      real fPARmin   ! Minimum possible FPAR corresponding to 2nd percentile
!
! begin output variables
      real fPAR      ! area average Canopy absorbed fraction of PAR
!
! begin internal variables
      real scatp     ! Canopy transmittance + reflectance coefficient wrt PAR
      real PARk      ! mean canopy absorption optical depth wrt PAR
!
! Calculate canopy transmittance + reflectance coefficient wrt PAR
! transmittance + reflectance coefficient=green plants + brown plants
      scatp=Green*(LTran(1,1)+LRef(1,1))+ &
      (1.-Green)*(LTran(1,2)+LRef(1,2))
!
! Calculate PAR absorption optical depth in canopy adjusting for 
! variance in projected leaf area wrt solar zenith angle
! (Sellers et al. Part II (1996), eqn. 13b)
! PAR absorption coefficient=(1-scatp)
      PARk=sqrt(1.-scatp)*gmudmu
!
! Calculate the new fPAR (Sellers et al. Part II (1996), eqn. 9)
      fPAR=(1.-exp(-PARk*LAI/fVCover))
!
! Ensure calculated fPAR falls within physical limits
      fPAR=amax1(fPARmin,fPAR)
      fPAR=amin1(fPARmax,fPAR)
!
! convert to area average
      fPAR=fVCover*fPAR
!
      return
      end
!
!******************************************************************
      subroutine textclass(text)
!******************************************************************
! Assigns soil texture classes based on the USDA texture triangle
! using subroutines developed by aris gerakis
!******************************************************************
!* +-----------------------------------------------------------------------
!* |                         T R I A N G L E
!* | Main program that calls WHAT_TEXTURE, a function that classifies soil
!* | in the USDA textural triangle using sand and clay %
!* +-----------------------------------------------------------------------
!* | Created by: aris gerakis, apr. 98 with help from brian baer
!* | Modified by: aris gerakis, july 99: now all borderline cases are valid
!* | Modified by: aris gerakis, 30 nov 99: moved polygon initialization to
!* |              main program
!* +-----------------------------------------------------------------------
!* | COMMENTS
!* | Supply a data file with two columns, in free format: 1st column sand,
!* |   2nd column clay %, no header.  The output is a file with the classes.
!* +-----------------------------------------------------------------------
!* | You may use, distribute and modify this code provided you maintain
!* ! this header and give appropriate credit.
!* +-----------------------------------------------------------------------
!
! Modifications:
!   Lara Prihodko customized triangle program for mapper (1/31/01)
!
      implicit none

      type text_type
        real :: clay  !  Percent clay content
        real :: silt  !  Percent silt content
        real :: sand  !  Percent sand content
        integer :: class  !  Soil texture class
      end type text_type

      type(text_type) :: text

      integer    :: what_texture
      real       :: sand, clay
      real       :: silty_loam(1:7,1:2), sandy(1:7,1:2), &
      silty_clay_loam(1:7,1:2), &
      loam(1:7,1:2), clay_loam(1:7,1:2), sandy_loam(1:7,1:2), &
      silty_clay(1:7,1:2), sandy_clay_loam(1:7,1:2), &
      loamy_sand(1:7,1:2), clayey(1:7,1:2), silt(1:7,1:2), &
      sandy_clay(1:7,1:2)

!Initalize polygon coordinates:

      data silty_loam/0, 0, 23, 50, 20, 8, 0, 12, 27, 27, 0, 0, 12, 0/
      data sandy/85, 90, 100, 0, 0, 0, 0, 0, 10, 0, 0, 0, 0, 0/
      data silty_clay_loam/0, 0, 20, 20, 0, 0, 0, 27, 40, 40, 27, 0, 0, 0/
      data loam/43, 23, 45, 52, 52, 0, 0, 7, 27, 27, 20, 7, 0, 0/
      data clay_loam/20, 20, 45, 45, 0, 0, 0, 27, 40, 40, 27, 0, 0, 0/
      data sandy_loam/50, 43, 52, 52, 80, 85, 70, 0, 7, 7, 20, 20, 15, 0/
      data silty_clay/0, 0, 20, 0, 0, 0, 0, 40, 60, 40, 0, 0, 0, 0/
      data sandy_clay_loam/52, 45, 45, 65, 80, 0, 0, 20, 27, 35, 35, 20, 0, 0/
      data loamy_sand/70, 85, 90, 85, 0, 0, 0, 0, 15, 10, 0, 0, 0, 0/
      data clayey/20, 0, 0, 45, 45, 0, 0, 40, 60, 100, 55, 40, 0, 0/
      data silt/0, 0, 8, 20, 0, 0, 0, 0, 12, 12, 0, 0, 0, 0/
      data sandy_clay/45, 45, 65, 0, 0, 0, 0, 35, 55, 35, 0, 0, 0, 0/


      sand = 0.0
      clay = 0.0

!Read input:

      sand = text%sand
      clay = text%clay

!Call function that estimates texture and put into structure:
      text%class = what_texture (sand, clay, silty_loam, sandy, &
      silty_clay_loam, &
      loam, clay_loam, sandy_loam, silty_clay, &
      sandy_clay_loam, loamy_sand, clayey, silt, &
      sandy_clay)

      return
      end
!
!******************************************************************
!* +-----------------------------------------------------------------------
!* | WHAT TEXTURE?
!* | Function to classify a soil in the triangle based on sand and clay %
!* +-----------------------------------------------------------------------
!* | Created by: aris gerakis, apr. 98
!* | Modified by: aris gerakis, june 99.  Now check all polygons instead of
!* | stopping when a right solution is found.  This to cover all borderline
!* | cases.
!* +-----------------------------------------------------------------------

      FUNCTION what_texture (sand, clay, silty_loam, sandy, &
      silty_clay_loam, loam, clay_loam, &
      sandy_loam, silty_clay, sandy_clay_loam, &
      loamy_sand, clayey, silt, sandy_clay)

      implicit none

!Declare arguments:

      real, intent(in) :: clay, sand, silty_loam(1:7,1:2), &
      sandy(1:7,1:2), &
      silty_clay_loam(1:7,1:2), loam(1:7,1:2), &
      clay_loam(1:7,1:2), sandy_loam(1:7,1:2), &
      silty_clay(1:7,1:2), sandy_clay_loam(1:7,1:2), &
      loamy_sand(1:7,1:2), clayey(1:7,1:2), silt(1:7,1:2), &
      sandy_clay(1:7,1:2)

!Declare local variables:

      logical :: inpoly
      integer :: texture, what_texture

!Find polygon(s) where the point is.

      texture = 13

      if (sand .gt. 0.0 .and. clay .gt. 0.0) then

      if (inpoly(silty_loam, 6, sand, clay)) then
            texture = 4
         endif
      if (inpoly(sandy, 3, sand, clay)) then
            texture = 1
         endif
      if (inpoly(silty_clay_loam, 4, sand, clay)) then
            texture = 10
         endif
      if (inpoly(loam, 5, sand, clay)) then
            texture = 6
         endif
      if (inpoly(clay_loam, 4, sand, clay)) then
            texture = 9
         endif
      if (inpoly(sandy_loam, 7, sand, clay)) then
            texture = 3
         endif
      if (inpoly(silty_clay, 3, sand, clay)) then
            texture = 11
         endif
      if (inpoly(sandy_clay_loam, 5, sand, clay)) then
            texture = 7
         endif
      if (inpoly(loamy_sand, 4, sand, clay)) then
            texture = 2
         endif
      if (inpoly(clayey, 5, sand, clay)) then
            texture = 12
         endif
      if (inpoly(silt, 4, sand, clay)) then
            texture = 5
         endif
      if (inpoly(sandy_clay, 3, sand, clay)) then
            texture = 8
         endif

      endif

      if (sand == 100) then
            texture = 1
      endif

      if (clay == 100) then
            texture = 12
      endif

      if (texture == 13 ) then
         texture = 13
      endif

      what_texture = texture


      END FUNCTION what_texture
!
!******************************************************************
!--------------------------------------------------------------------------
!                            INPOLY
!   Function to tell if a point is inside a polygon or not.
!--------------------------------------------------------------------------
!   Copyright (c) 1995-1996 Galacticomm, Inc.  Freeware source code.
!
!   Please feel free to use this source code for any purpose, commercial
!   or otherwise, as long as you don't restrict anyone else's use of
!   this source code.  Please give credit where credit is due.
!
!   Point-in-polygon algorithm, created especially for World-Wide Web
!   servers to process image maps with mouse-clickable regions.
!
!   Home for this file:  http://www.gcomm.com/develop/inpoly.c
!
!                                       6/19/95 - Bob Stein & Craig Yap
!                                       stein@gcomm.com
!                                       craig@cse.fau.edu
!--------------------------------------------------------------------------
!   Modified by:
!   Aris Gerakis, apr. 1998: 1.  translated to Fortran
!                            2.  made it work with real coordinates
!                            3.  now resolves the case where point falls
!                                on polygon border.
!   Aris Gerakis, nov. 1998: Fixed error caused by hardware arithmetic
!   Aris Gerakis, july 1999: Now all borderline cases are valid
!--------------------------------------------------------------------------
!   Glossary:
!   function inpoly: true=inside, false=outside (is target point inside
!                    a 2D polygon?)
!   poly(*,2):  polygon points, [0]=x, [1]=y
!   npoints: number of points in polygon
!   xt: x (horizontal) of target point
!   yt: y (vertical) of target point
!--------------------------------------------------------------------------

      FUNCTION inpoly (poly, npoints, xt, yt)

      implicit none

!Declare arguments:

      integer :: npoints
      real, intent(in)    :: poly(7, 2), xt, yt

!Declare local variables:

      real    :: xnew, ynew, xold, yold, x1, y1, x2, y2
      integer :: i
      logical :: inside, on_border, inpoly

      inside = .false.
      on_border = .false.

      if (npoints < 3)  then
        inpoly = .false.
        return
      end if

      xold = poly(npoints,1)
      yold = poly(npoints,2)

      do i = 1 , npoints
        xnew = poly(i,1)
        ynew = poly(i,2)

      if (xnew > xold)  then
          x1 = xold
          x2 = xnew
          y1 = yold
          y2 = ynew
      else
          x1 = xnew
          x2 = xold
          y1 = ynew
          y2 = yold
      end if

!The outer IF is the 'straddle' test and the 'vertical border' test.
!The inner IF is the 'non-vertical border' test and the 'north' test.

!The first statement checks whether a north pointing vector crosses
!(stradles) the straight segment.  There are two possibilities, depe-
!nding on whether xnew < xold or xnew > xold.  The '<' is because edge
!must be "open" at left, which is necessary to keep correct count when
!vector 'licks' a vertix of a polygon.

      if ((xnew < xt .and. xt <= xold) .or. (.not. xnew < xt .and. &
      .not. xt <= xold)) then
!The test point lies on a non-vertical border:
      if ((yt-y1)*(x2-x1) == (y2-y1)*(xt-x1)) then
        on_border = .true.
!Check if segment is north of test point.  If yes, reverse the
!value of INSIDE.  The +0.001 was necessary to avoid errors due
!arithmetic (e.g., when clay = 98.87 and sand = 1.13):
      elseif ((yt-y1)*(x2-x1) < (y2-y1)*(xt-x1) + 0.001) then
        inside = .not.inside ! cross a segment
      endif
!This is the rare case when test point falls on vertical border or
!left edge of non-vertical border. The left x-coordinate must be
!common.  The slope requirement must be met, but also point must be
!between the lower and upper y-coordinate of border segment.  There
!are two possibilities,  depending on whether ynew < yold or ynew >
!yold:
      elseif ((xnew == xt .or. xold == xt) .and. (yt-y1)*(x2-x1) == &
      (y2-y1)*(xt-x1) .and. ((ynew <= yt .and. yt <= yold) .or. &
      (.not. ynew < yt .and. .not. yt < yold))) then
        on_border = .true.
      endif

      xold = xnew
      yold = ynew

      enddo

!If test point is not on a border, the function result is the last state
!of INSIDE variable.  Otherwise, INSIDE doesn't matter.  The point is
!inside the polygon if it falls on any of its borders:

      if (.not. on_border) then
         inpoly = inside
      else
         inpoly = .true.
      endif
!
      END FUNCTION inpoly
!
!=======================================================================
      subroutine FractionVegCover(NDVI,fPARmax,fPARmin, &
      SRmax,SRmin,fVCover)
!=======================================================================
! calculates the vegetation cover fraction for a single pixel.
! The maximum fPAR for pixel during entire year determines vegetation
! cover fraction.  Maximum yearly NDVI corresponds to maximum fPAR.
! Calculate fPAR from Simple Ratio using an empirical linear relationship.
!
      implicit none
!
! begin input variables
      real NDVI      ! maximum FASIR NDVI value for a grid cell
      real fPARmax   ! Maximum possible FPAR corresponding to 98th percentile
      real fPARmin   ! Minimum possible FPAR corresponding to 2nd percentile
      real SRmax     ! Maximum simple ratio for biome type
      real SRmin     ! Minimum simple ratio for biome type
!
! begin output variables
      real fVCover   ! fractional vegetation cover
!
! begin internal variables
      real fPAR      ! maximum fPAR associated with maximum ndvi
!
! The maximum fPAR for pixel during entire year determines vegetation
! cover fraction.  Maximum yearly NDVI corresponds to maximum fPAR.
! Calculate fPAR from Simple Ratio using an empirical linear relationship.
!
       call srapar (NDVI,  &
                   SRmin, &
                   SRmax, &
                   fPAR, &
                   fPARmax, &
                   fPARmin)
!
! calculate fractional vegetation cover
      fVCover=fPAR/fPARmax
!
        return
        end
!
!=======================================================================
      subroutine SoilProperties(text,SoilVar)
!=======================================================================
! calculates soil physical properties given sand and clay content
!
! Modifications
!  Kevin Schaefer created subroutine for soil hydraulic properties (4/22/00)
!  Kevin Schaefer resp. variable curve fits from Raich et al., 1991 (6/19/00)
!  Kevin Schaefer combine code for hydraulic & respiration variables (3/30/01)
!
      implicit none
!
! begin Input variables
      type text_type
        real :: clay  !  Percent clay content
        real :: silt  !  Percent silt content
        real :: sand  !  Percent sand content
        integer :: class  !  Soil texture class
      end type text_type

      type(text_type) :: text
!
! begin Output soil property variables
      type soil_Physical
        integer SoilNum ! soil type number
        real BEE     ! Wetness exponent for soil conductivity (-)
        real PhiSat  ! Soil matrix potential (water tension) at Saturation (m)
        real SatCo   ! Soil Hydraulic Conductivity at Saturation (m/s)
        real poros   ! Soil Porosity or saturation water content (-)
        real Slope   ! Cosine of mean slope
        real Wopt    ! Optimal soil moisture for respiration (-)
        real Skew    ! skewness exponent of soil respiration vs. wetness curve
        real RespSat ! Parameter determining soil respiration at saturation (-)
      end type soil_Physical
!
      type(soil_Physical) SoilVar  ! time ind., soil dependant variables
!
! begin local variables
      real fclay   ! fraction of clay in soil
      real fsand   ! fraction of sand in soil
!
! calculate Soil hydraulic and thermal variables based on Klapp and Hornberger
      SoilVar%PhiSat=-10.*10**(1.88-0.0131*Text%Sand)/1000.
      SoilVar%poros=0.489-0.00126*Text%Sand
      SoilVar%SatCo=0.0070556*10**(-0.884+0.0153*Text%sand)/1000.
      SoilVar%bee=2.91+0.159*Text%Clay
!
! Calculate clay and sand fractions from percentages
      fclay=Text%clay/100.
      fsand=Text%sand/100.
!
! Calculate soil respiration variables based on curve fits to 
! data shown in Raich et al. (1991)
      SoilVar%Wopt=(-0.08*fclay**2+0.22*fclay+0.59)*100.
      SoilVar%Skew=-2*fclay**3-0.4491*fclay**2+0.2101*fclay+0.3478
      SoilVar%RespSat=0.25*fclay+0.5
!
! assign value for mean slope of terrain
      SoilVar%Slope=0.176
!
      return
      end
!
!=======================================================================      
      subroutine gmuder (Lat, DOY, ChiL, gmudmu)
!=======================================================================      
! calculates daily time mean optical depth of canopy relative to the Sun.
!
      implicit none
!
! begin input variables
      real Lat      ! latitude in degrees
      real DOY      ! day-of-year (typically middle day of the month)
      real ChiL     ! leaf angle distribution factor
!
! begin output variables
      real gmudmu   ! daily time mean canopy optical depth relative to Sun
!
! begin internal variables
      real mumax    ! max cosine of the Solar zenith angle (noon)
      real mumin    ! min cosine of the Solar zenith angle (rise/set)
      real dec      ! declination of the Sun (Solar Declination)
      real pi180    ! conversion factor from degrees to radians
      real aa       ! minimum possible LAI projection vs. cosine Sun angle
      real bb       ! slope leaf area projection vs. cosine Sun angle
!
! Calculate conversion factor from degrees to radians
      pi180=3.14159/180. 
!
! Calculate solar declination in degrees
      dec=23.5*sin(1.72e-2*(DOY-80.))
!
! Calculate maximum cosine of zenith angle corresponding to noon
      mumax=cos((dec-lat)*pi180)
      mumax=max(0.02, mumax)
!
! Assign min cosine zenith angle corresponding to start disc set (cos(89.4))
      mumin=0.01
!
! The projected leaf area relative to the Sun is G(mu)=aa+bb*mu
! Calculate minimum projected leaf area
      aa=0.5-0.633*ChiL-0.33*ChiL*ChiL
!
! Calculate slope of projected leaf area wrt cosine sun angle
      bb=0.877*(1.-2.*aa) 
!
! Calculate mean optical depth of canopy by integrating G(mu)/mu over
! all values of mu.  Since G(mu) has an analytical form, this comes to
      gmudmu=aa*alog(mumax/mumin)/(mumax-mumin)+bb
!
      return                                                                    
      end
!
