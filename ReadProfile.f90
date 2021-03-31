!__________________________________________________
! unit ReadProfile
! Name:
!        ReadProfile
!
! Purpose:
!        sets the model up to read in mm5 prof%les
!
! Category:
!        ETL ET1 RTM
!
! Language:
!        Fortran 90
!
! Modules:
!        type_kinds for variable definitions
!
! Files Accessed:
!
!        ./prof/002708 - one profile converted mm5 output file
!
! Creation History:
!             Converted from ReadProfile.pas by Ronald Richter NOAA/ETL
!             Software Engineering Team, April 20, 2003
!
! Define or type a structure to contain an mm5 prof%le
!
!  module ReadProfile
!  use type_kinds

!      INTEGER, DIMENSION(4) :: dim_type 
!      REAL(fp_kind), DIMENSION(59) :: data_array_type

!      TYPE prof_data_type 
!         real(fp_kind), dimension(59) :: u
!         real(fp_kind), dimension(59) :: v
!         real(fp_kind), dimension(59) :: pp
!         real(fp_kind), dimension(59) :: temperature
!         real(fp_kind), dimension(59) :: clw
!         real(fp_kind), dimension(59) :: rnw
!         real(fp_kind), dimension(59) :: ice
!         real(fp_kind), dimension(59) :: snow
!         real(fp_kind), dimension(59) :: grauple
!         real(fp_kind), dimension(59) :: altitude
!         real(fp_kind), dimension(59) :: q
!         real(fp_kind), dimension(59) :: sigmaH
!         REAL(fp_kind)   :: terrain
!         REAL(fp_kind)   :: ground_temp
!         REAL(fp_kind)   :: pstarcrs
!         REAL(fp_kind)   :: lat
!         REAL(fp_kind)   :: lng
!         REAL(fp_kind)   :: TopPress
!         REAL(fp_kind)   :: Sealevel_press
!         REAL(fp_kind)   :: Sealevel_temp
!         REAL(fp_kind)   :: Lapse_rate
!         INTEGER         :: X, Y, Z, T 
!      END TYPE prof_data_type
!      logical isfileopen
!____________________________________________________________
!CONTAINS
!_____________________________________________________________

subroutine getprofile(prof, prof_file)
!  This subroutine inputs atmospheric data and surface data; then it computes and 
!  organizes that data for use in the rest of the program. 
use variables
! use type_kinds
implicit none
type (prof_data_type) :: prof
integer prof_file
! integer, dimension(4) :: dims
  REAL(8)    :: wind
integer i , ios, k1, k
  REAL(8)              :: r,g,ps0,phydro
  REAL(8)     :: Ryd,ma,mw,wa,xx  
ios = 1

 if ( ios > 0 ) then
!      read(UNIT=prof_file,FMT=110,IOSTAT=ios)  prof%X, prof%Y, prof%Z, prof%lat, &
      read(UNIT=prof_file,FMT=*,IOSTAT=ios)  prof%X, prof%Y, prof%Z, prof%lat, &
      prof%lng, prof%ground_temp, prof%terrain, prof%pstarcrs,    &
         prof%TopPress, prof%Sealevel_press, prof%Sealevel_temp,  &
         prof%Lapse_rate
! 110   format(1x,3i4,9f16.6)
      do i=1, dims(3)
!         read(UNIT=prof_file, FMT=120, IOSTAT=ios)  prof%u(i), prof%v(i), prof%temperature(i), prof%pp(i),  & 
         read(UNIT=prof_file, FMT=*, IOSTAT=ios)  prof%u(i), prof%v(i), prof%temperature(i), prof%pp(i),  &
         prof%q(i), prof%clw(i), prof%rnw(i), prof%ice(i), prof%snow(i),         &
         prof%grauple(i), prof%sigmaH(i)
      end do
! 120   format(1x,11f16.6)

     if (prof%terrain .LT. 0) then
         k1 = dims(3)
         surf_inp%surf_temp = prof%ground_temp
         wind = dsqrt(prof%u(k1)*prof%u(k1)+prof%v(k1)*prof%v(k1))
         call construct_surf('W','O',1.0d0, 1.0d0,0.035d0,wind)               
     else
         surf_inp%surf_temp = prof%ground_temp
         call construct_surf(' ','L',0.05d0,1.0d0,1.0d0,1.0d0)
     end if
     do i=1, atm_inp%num_levels
        k = atm_inp%num_levels + 1 - i
!       {altitude (m)}
        g = 9.81d0
        r = 287.04d0
        xx = prof%Sealevel_temp/prof%Lapse_rate
        ps0 = prof%Sealevel_press * dexp( -xx          &
            + ((xx*xx-2.d0*g*prof%terrain/(prof%Lapse_rate*r))**0.5d0) )
        ps0 = ps0-prof%TopPress
        phydro = ps0*prof%sigmaH(i)+prof%TopPress
        xx = phydro / prof%Sealevel_press
        xx = dlog(xx)
        prof%altitude(i) =-(r/g)*xx*( 0.5d0 * prof%Lapse_rate * xx   &
                           + prof%Sealevel_temp )
!       {if terrain height is above sea level,
!        then must be substracted from altitude}
        if (prof%terrain.GE. 0) then
            atm_inp%prof(k)%height=(prof%altitude(i)-prof%terrain)/1.0d3
        else
            atm_inp%prof(k)%height=prof%altitude(i)/1.0d3
        end if
!       {pressure (Pa)}
        atm_inp%prof(k)%pressure = prof%pstarcrs*prof%sigmaH(i)   &
                                  +prof%TopPress+prof%pp(i)
        atm_inp%prof(k)%temperature =prof%temperature(i)
        Ryd = 8.31436d0
        ma = 28.966d0
        mw = 18.016d0
        wa = mw/ma
        atm_inp%prof(k)%vapor_density = (atm_inp%prof(k)%pressure             &
                                      / (Ryd*prof%temperature(i)))           &
                                      * (prof%q(i)*ma/(1.0d0+prof%q(i)*wa))
        atm_inp%prof(k)%cloud_liq_dens = prof%clw(i)                         &
                                       * (atm_inp%prof(k)%pressure            &
                                       / (Ryd*prof%temperature(i)))          &
                                       * ma/(1.0d0+prof%q(i)*wa)               
        atm_inp%prof(k)%cloud_rn_dens = prof%rnw(i)                          &
                                      * (atm_inp%prof(k)%pressure             &
                                      / (Ryd*prof%temperature(i)))           &
                                      * ma/(1.0d0+prof%q(i)*wa)
        atm_inp%prof(k)%cloud_ice_dens = prof%ice(i)                         &
                                       * (atm_inp%prof(k)%pressure            &
                                       / (Ryd*prof%temperature(i)))          &
                                       * ma/(1.0d0+prof%q(i)*wa)
        atm_inp%prof(k)%cloud_snow_dens = prof%snow(i)                       &
                                        * (atm_inp%prof(k)%pressure           &
                                        / (Ryd*prof%temperature(i)))         &
                                        * ma/(1.0d0+prof%q(i)*wa)
        atm_inp%prof(k)%cloud_grpl_dens = prof%grauple(i)                    &
                                        * (atm_inp%prof(k)%pressure           &
                                        / (Ryd*prof%temperature(i)))         &
                                        * ma/(1.0d0+prof%q(i)*wa)
        ! pressure (mb) (1 mb = 100 Pa)
        atm_inp%prof(k)%pressure = atm_inp%prof(k)%pressure/100.0d0
     end do ! i
        atm_inp%inf= .true.
        atm_inp%fiveph= .true.
 else
      prof%X              = -1
      prof%Y              = -1
      prof%Z              = -1
      prof%lat            = -9999.0
      prof%lng            = -9999.0
      prof%ground_temp    = -9999.0
      prof%terrain        = -9999.0
      prof%pstarcrs       = -9999.0
      prof%TopPress       = -9999.0
      prof%Sealevel_press = -9999.0
      prof%Sealevel_temp  = -9999.0
      prof%Lapse_rate     = -9999.0

      do i = 1, dims(3) 
         prof%u(i)           = -9999.0
         prof%v(i)           = -9999.0
         prof%temperature(i) = -9999.0
         prof%pp(i)          = -9999.0
         prof%q(i)           = -9999.0
         prof%clw(i)         = -9999.0
         prof%rnw(i)         = -9999.0
         prof%ice(i)         = -9999.0
         prof%snow(i)        = -9999.0
         prof%grauple(i)     = -9999.0
         prof%sigmaH(i)      = -9999.0
      end do
        atm_inp%inf= .false.
        atm_inp%fiveph= .false.
 end if
!//   READ(prof_file)
 return
 end subroutine getprofile
!________________________________________________________________________
!SUBROUTINE closeprofilefile(prof_file)
!      integer prof_file

!      close(prof_file)
!      isfileopen = .false.

!end subroutine closeprofilefile


! end module readprofile
