!
!=======================================================================
      subroutine ZenithAngle (len, sind, cosd, &
        latsib, lonsib, tofday, cosz, test, test10)
!=======================================================================      
! calculates the zenith angle for an array of lat/lon points.
! The cosH variable accounts for the difference between GMT and local time
!
! Modifications:
!   Kevin Schaefer created ZenithAngle subroutine (5/31/01)
!   Kevin Schaefer corrected hour angle calculation (5/31/01)
!   Kevin Schaefer removedminimum cosz limit (7/21/01)
!----------------------------------------------------------------------
!
      implicit none
!
! input variables
      integer len      ! number of points
      real sind        ! sine of solar declination
      real cosd        ! cosine of solar declination
      real latsib(len) ! SiB point latitude
      real lonsib(len) ! SiB point longitude
      real tofday      ! time of day (GMT, in hours)
      real test, test10(10)     ! test variables
!
! output variables
      real cosz(len)   ! cosine of Solar zenith angle
!
! internal variables
      real pid180      ! conversion from degrees to radians
      real sinlat      ! sine of latitude
      real coslat      ! cosine of latitude
      real HrAng       ! hour angle; longitude of Sun from Greenwhich meridian
      real cosH        ! cosine delta longitude between Sun & SiB point
      integer i        ! sib point index
!
! calculate conversion from degrees to radians
      pid180=3.14159/180.
!
! Calculate hour angle (longitude of Sun from Greenwhich meridian)
      HrAng=(12.-tofday)*360./24.
!
! Calculate cosine of solar zenith angle for each SiB point
      do i = 1,len
        cosH=cos(pid180*(HrAng-lonsib(i)))
        sinlat=sin(pid180*latsib(i))
        coslat=cos(pid180*latsib(i))
        cosz(i)=coslat*cosd*cosH+sinlat*sind
      enddo
      test = cosz(1)
      test10(1) = cosz(1)
!
      return                                                                    
      end
