!
!=======================================================================
      SUBROUTINE respire(len,nsoil,nm, &
        assimn_M,soilscale_M,xroot,respFactor)
!=======================================================================
!  Calculates the annual respiration rate "respFactor" for each of nsoil
!   soil layer at each grid cell in the model, given monthly mean
!   maps of net carbon assimilation and "soilscale" at each level.
!
!  REFERENCES:
!  Denning et al., 1996a, Tellus 48B, 521-542
!  Denning at al., 1996b, Tellus 48B, 543-567 
!
! SoilScale is the product of a temperature and moisture response function,
! evaluated in each soil layer called (R*) 
! Denning et al (96a), eqn 6-9; (96b), eqn 8-9
!
! Modifications:
!  Kevin Schaefer made into SiB compatable subroutine
!
! Input variables
      INTEGER len      ! number of SiB points
      INTEGER nsoil 
      INTEGER nm              ! number of time points per year
      REAL assimn_M(len,nm) ! monthly mean net carbon assimilation (moles/m2/sec)
      REAL soilscale_M(len,nsoil,nm)  ! monthly mean Soilscale per soil layer
      REAL xroot(len,nsoil)     ! root Fraction per soil layer 
!
! Output variable
      REAL respFactor(len,nsoil) ! mean unstressed soil resp rate (mol/m2/sec)
!
! Local variables
      INTEGER DayPerMon(12)  ! Number of days in each month
      INTEGER month      ! Array index (month counter)
!
      REAL ANA(len)   ! (g C/m2)Total Annual net assimilation
      REAL SSA(len,nsoil) ! (-) Total Annual soilscale per soil layer
      REAL xAG(len)   ! (-) Fraction annual respiration from above-ground sources
      REAL xAGmin     ! (-) Min above-grnd resp fraction (low GPP ecosystems)
      REAL xAGmax     ! (-) Max above-grnd resp fraction (high GPP ecosystems)
      REAL ANAinflec  ! (gC/m2) ANA at inflection point for xAG function
      REAL kxAG       ! (m2/gC) Exponential constant for xAG function
      REAL xBG(len)   ! (-) Fraction annual respiration from below-ground sources
      REAL SecPerDay  ! (s/day) Number of seconds in a day
      REAL SecPerMon  ! (s/mon) Number of seconds in a given month
!
      DATA DayPerMon/31,28,31,30,31,30,31,31,30,31,30,31/
!
! assign local constants
      xAGmin=.10
      xAGmax=.75
      ANAinflec=1000.
      kxAG=5.e-3
      SecPerDay=24.*60.*60.
!
! Calculate annual total assimulation and soilscale
      ANA=0.
      SSA=0.
      do month=1, nm
        SecPerMon=DayPerMon(month)*SecPerDay
        ANA(:)=ANA(:)+assimn_M(:,month)*SecPerMon
        SSA(:,:)=SSA(:,:)+soilscale_M(:,:,month)*SecPerMon
      end do
!
! Calculate Above-ground fraction of respiration
! (assimilation converted from moles to grams C, factor of 12 below)
      xAG=xAGmin+(xAGmax-xAGmin)/(1.+exp(-kxAG*(ANA*12.-ANAinflec)))
!
! Calculate Below-ground fraction of respiration
      xBG=1.-xAG
!
! Calculate respfactors starting from the top and working down
      respFactor(:,7)=ANA(:)*(0.5*xAG(:))/SSA(:,7)
      respFactor(:,6)=ANA(:)*(0.5*xAG(:)+xBG(:)*xRoot(:,6))/SSA(:,6)
      respFactor(:,5)=ANA(:)*xBG(:)*xRoot(:,5)/SSA(:,5)
      respFactor(:,4)=ANA(:)*xBG(:)*xRoot(:,4)/SSA(:,4)
      respFactor(:,3)=ANA(:)*xBG(:)*xRoot(:,3)/SSA(:,3)
      respFactor(:,2)=0.
      respFactor(:,1)=0.
!
      RETURN
      END
!
