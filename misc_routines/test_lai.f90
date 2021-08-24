!
!=======================================================================
      subroutine TestLAI (fPAR,fPARm,fPARmax,fVCover,stems, &
       LAImax,Green,LAI, test, test10)
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
      real test, test10(10)     ! test variables
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

!
      if (LAIg.ge.LAIgm) LAId=0.0001
      if (LAIg.lt.LAIgm) LAId=LAIgm-LAIg
      
      test = lai
      test10(1) = lai
!
      return                                                                    
      end
