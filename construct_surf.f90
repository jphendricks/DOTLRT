!======================================================================
subroutine construct_surf_characteristics(ilon, ilat) 
!======================================================================
! Constructs model of ocean surface from several options listed in the code below.
! Language: FORTRAN 90 Portland Group Compiler on Red Hat Linux
!
! History: 
!   Original from NOAA ETL/ET1 Pascal GMRT version by Bob Weber
!   April 2003 Ronald Richter Converted from original Pascal ron.richter@noaa.gov
!   9/26/2020 Kevin Schaefer deleted unused variables
!   10/5/2020 Kevin Schaefer commented the code
!   10/5/2020 Kevin Schaefer changed variable names and got rid of select case
!  10/16/2020 Kevin Schaefer switched to variable module and replaced duplicate variables
!  10/16/2020 Kevin Schaefer separated dielectric and reflectivity calculations
!----------------------------------------------------------------------
use dotlrt_variables
use profiles
implicit none
!
! input variables
!
! input variables
  integer ilon     ! (-) longitude index
  integer ilat     ! (-) latitude index

! local variables
  real(8) slope_var
  real(8) foam_frac ! (-) foam fraction of total area
  integer(4) iang   ! angle index
  double complex kappa             ! (-) complex dielectric constant of water
  double complex mu                ! (-) real component of complex number

! General surface characteristics
  surf_inp%surf_temp = TSK(ilon, ilat)
  surf_inp%land_flag = landmask(ilon, ilat)
  surf_inp%ocean_mod = ocean_mod
  surf_inp%nstream_surf = nstream_surf
  if (trim(surf_inp%ocean_mod) == 'Wilheit') surf_inp%windspeed = wind(ilon, ilat)

! common to all surface types
    do iang=1, surf_inp%nstream_surf
      surf_inp%theta(iang)=quad_angle(iang)  
    enddo 

!----------------------------------------------------------------------
! calculate dielectric constant
!----------------------------------------------------------------------
  if (surf_inp%land_flag==2) then ! land
    surf_inp%salinity=0.d0
    surf_inp%diel_mod='none'
    surf_inp%dielectric=dcmplx(0.0d0, 0.0d0)
  elseif (surf_inp%land_flag==1) then ! lake
    surf_inp%salinity=salinity_lake
    surf_inp%diel_mod='Fresh'
    call dielectric_const_water(channel%lo_freq)
  elseif (surf_inp%land_flag==0) then ! ocean
    surf_inp%salinity=salinity_ocean
    surf_inp%diel_mod='Salt'
    call dielectric_const_water(channel%lo_freq)
  endif

!----------------------------------------------------------------------
! Calculate reflectivities 
!----------------------------------------------------------------------
  if (surf_inp%land_flag==2) then ! land
    do iang=1, surf_inp%nstream_surf
      surf_inp%vr(iang)=land_refl 
      surf_inp%hr(iang)=land_refl 
    enddo 

  else ! lake or ocean
    if (trim(ocean_mod)=='Fresnel') then ! Fresnel reflectivity model
      do iang=1, surf_inp%nstream_surf
        call fresnel_refl(iang, kappa, mu)
      enddo 

    elseif (trim(ocean_mod)=='Wilheit') then ! Wilheit reflectivity model
      slope_var=0.003d0+0.0048d0*surf_inp%windspeed 
      if (channel%lo_freq < 35.0d0) slope_var=slope_var*(0.3d0+0.02d0*channel%lo_freq) 
      if (surf_inp%windspeed < 7.0d0) then
        foam_frac=0.0d0
      else
        foam_frac=(1.0d0-dexp(-channel%lo_freq/7.5d0))*(surf_inp%windspeed-7.0d0)*0.006d0
      endif
      do  iang=1, surf_inp%nstream_surf 
        ! Compute slope-averaged reflectivity using Gaussian distribution (Kirchoff approximation)
        call kirchoff_ocean_refl(iang, channel%lo_freq, slope_var)
        surf_inp%vr(iang)=surf_inp%vr(iang)*(1.0d0-foam_frac)
        surf_inp%hr(iang)=surf_inp%hr(iang)*(1.0d0-foam_frac)
      enddo

    elseif (trim(ocean_mod)=='ITRA') then ! ITRA reflectivity model
      do iang=1,  surf_inp%nstream_surf
        surf_inp%vr(iang)=0.638d0 - 0.00272d0*channel%lo_freq
        surf_inp%hr(iang)=0.638d0 - 0.00272d0*channel%lo_freq
      enddo 
    endif
  endif

return
end subroutine construct_surf_characteristics
