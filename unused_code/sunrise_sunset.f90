!
!=======================================================================
      subroutine SunriseSunset
!=======================================================================      
! Calculates sunrise, sunset, noon, and midnight
!----------------------------------------------------------------------
!
      use Mapper_Variables
      use sibscan_Variables
!
      implicit none
!
! local variables
      real temp
      real nightlat ! latitude of 100% darkness
      real daylat   ! latitude of 100% sunlight
      real lenday   ! length of day
!
! calculate solar declination
      dec=23.5*sin(1.72e-2*(DOY-80))
      sind=sin(pi*dec/180.)
      cosd=cos(pi*dec/180.)
!
! Calculate noon and midnight
      noon=12.-24.*lon/360.
      midnight=12.-24.*(lon-180.)/360.
      if(midnight.gt.24.) midnight=midnight-24.
!
! Calculate sunrise and sunset
      nightlat=atan(-1./tan(dec*pi/180.))*180./pi
      daylat=atan(1./tan(dec*pi/180.))*180./pi
      if (dec*nightlat.lt.dec*lat.and.dec*lat.lt.dec*daylat) then
        temp=-tan(lat*pi/180.)*tan(dec*pi/180.)
        sunrise=12.-24./360.*(lon-180./pi*acos(temp))
        sunset=12.-24./360.*(lon+180./pi*acos(temp))
      else if (dec*nightlat.ge.dec*lat) then
        temp=1.
        sunrise=noon
        sunset=noon
      else if (dec*lat.ge.dec*daylat) then
        temp=-1.
        sunrise=midnight
        sunset=midnight
      endif
      if(sunrise.gt.24.) sunrise=sunrise-24.
      if(sunset.lt.0.) sunset=sunset+24.
      lenday=24./pi*acos(temp)
!
      return
      end
