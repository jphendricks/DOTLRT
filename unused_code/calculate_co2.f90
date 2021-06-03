!
!=======================================================================
      SUBROUTINE calculate_co2(year,lat,co2m,seas,test,test10)
!=======================================================================
! Calculate the active layer thickness given surface deformation assuming 
! exponential decrease in organic matter with depth
!
! Modifications:
!  Kevin Schaefer made subroutine (8/30/10)
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
! Input variables
real year     ! (yr) year
integer mon   ! (mon) month of year
real lat     ! (deg) year
!
! Output variable
real co2m  ! (ppm) mixed layer co2 concentration
!
! Local variables
integer i,j  ! indeces
real flask_lat(5)  ! flask latitudes
real flask_co2(5,12)  ! flask seasonal cycle
real year_per_month(12)  ! (year) fractional cumulative year at end of each month
real year_mid_month(12)  ! (year) fractional cumulative year at end of each month
real seas ! (ppm) seasonal variation adjustment
real temp ! temporary local variable
integer num1 ! southern flask number index for interpolation
integer num2 ! northern flask number index for interpolation
real del_time ! (day) delta time between trend values
real tpgf1    ! scaling factor between 1st flask value and current time
real tpgf2    ! scaling factor between current time and 2nd flask value
integer k1    ! 1st month index for interpolation
integer k2    ! 2nd month index for interpolation
real co2_1    ! scaling factor between 1st trend value and current time
real co2_2    ! scaling factor between current time and 2nd trend value
real test, test10(10)    ! test
!
! calculate global average co2 concentration
    if(year<2008.) then ! use curve fit of global observed co2
      co2m=280.+0.27*exp(0.019325329*(year-1700.))
    else ! use IPCC A1B projection to stabilize at 688 ppm by 2100
      co2m=3.3887*(year-2008.)+384.85
      co2m=min(co2m,700.) 
    endif
!
! determine month of year
    year_mid_month(:) = (/4.2465754E-02,0.1232877,0.2041096,0.2876712,0.3712329,&
      0.4547945,0.5383562,0.6232877,0.7068493,0.7904110,0.8739727,0.9575343 /)   
    year_per_month(:)=(/8.4931508E-02, 0.1616438, 0.2465753, 0.3287671, 0.4136986, &
      0.4958904, 0.5808219, 0.6657534, 0.7479451, 0.8328766, 0.9150684, 1. /)   
    j=year
    temp=year-real(j)
    mon=0
    do i=12,1,-1
      if(temp<=year_per_month(i))then
        mon=i
      endif
    enddo
!
! assign latitudes of flask data
    flask_lat(1)=90.     ! north pole
    flask_lat(2)=71.32   ! BRW (Barrow)
    flask_lat(3)=19.530  ! MLO (Mauna loa)
    flask_lat(4)=-89.980 ! SPO (South Pole Station)
    flask_lat(5)=-90.    ! south pole
!
! relative average seasonal cycle in co2 concentration (ppm)
! derived from flask observations trends and mean values removed
! BRW Season
    flask_co2(2,1)=0.32268710E+01
    flask_co2(2,2)=0.39178782E+01
    flask_co2(2,3)=0.42880478E+01
    flask_co2(2,4)=0.45407581E+01
    flask_co2(2,5)=0.48192682E+01
    flask_co2(2,6)=0.25147221E+01
    flask_co2(2,7)=-0.45384240E+01
    flask_co2(2,8)=-0.10298231E+02
    flask_co2(2,9)=-0.87144508E+01
    flask_co2(2,10)=-0.35581274E+01
    flask_co2(2,11)=0.55232322E+00
    flask_co2(2,12)=0.32477596E+01
!
! MLO Season
    flask_co2(3,1)=-0.21220398E+00
    flask_co2(3,2)=0.57764530E+00
    flask_co2(3,3)=0.16015472E+01
    flask_co2(3,4)=0.27818632E+01
    flask_co2(3,5)=0.33004551E+01
    flask_co2(3,6)=0.26449833E+01
    flask_co2(3,7)=0.18991274E+00
    flask_co2(3,8)=-0.19132940E+01
    flask_co2(3,9)=-0.33598199E+01
    flask_co2(3,10)=-0.31779656E+01
    flask_co2(3,11)=-0.17802587E+01
    flask_co2(3,12)=-0.33074304E+00
!
! SPO Season     
    flask_co2(4,1)=-0.40703124E+00
    flask_co2(4,2)=-0.10044072E+01
    flask_co2(4,3)=-0.10342996E+01
    flask_co2(4,4)=-0.86555952E+00
    flask_co2(4,5)=-0.61372006E+00
    flask_co2(4,6)=-0.31292197E+00
    flask_co2(4,7)=0.49972042E-01
    flask_co2(4,8)=0.52145582E+00
    flask_co2(4,9)=0.80552232E+00
    flask_co2(4,10)=0.91416389E+00
    flask_co2(4,11)=0.89274549E+00
    flask_co2(4,12)=0.79357713E+00
!
! Assume south pole has same values as SPO; north pole same as BRW
    do i=1,12
      flask_co2(1,i)=flask_co2(2,i)
      flask_co2(5,i)=flask_co2(4,i)
    enddo
!
! find bracketing flask stations
    lat=min(lat,90.)
    do i=5,1, -1
      if(flask_lat(i)<=lat)then
        num1=i-1
        num2=i
      endif
    enddo
    if(lat==90.) then
        num1=1
        num2=2
    endif
!
! Calculate scaling factors for time interpolation of flask values
    j=year
    temp=year-real(j)
    test=temp
    if (temp<year_mid_month(mon)) then
        k1=mon-1
        k2=mon 
        if(k1<1) then
        k1=12
        del_time=year_mid_month(k2)-year_mid_month(k1)+1.
      else
        del_time=year_mid_month(k2)-year_mid_month(k1)
      endif
        tpgf1=(year_mid_month(k2)-temp)/del_time
        tpgf2=1.0-tpgf1
    else
        k1=mon
        k2=mon+1
        if(k2>12) then
          k2=1
          del_time=year_mid_month(k2)-year_mid_month(k1)+1.
        else
          del_time=year_mid_month(k2)-year_mid_month(k1)
        endif
        tpgf2=(temp-year_mid_month(k1))/del_time
        tpgf1=1.0-tpgf2
    endif
!
! interpolate co2 values in time
    co2_1 = tpgf1*flask_co2(num1,k1) + tpgf2*flask_co2(num1,k2)
    co2_2 = tpgf1*flask_co2(num2,k1) + tpgf2*flask_co2(num2,k2)
!
! interpolate seasonal CO2 cycle in latitude between flask stations
    seas=(co2_1-co2_2)/(flask_lat(num1)-flask_lat(num2))*(lat-flask_lat(num2))+co2_2
    co2m=co2m+seas
    
    test = co2m
    test10(1) = co2m

    RETURN
    END
!
