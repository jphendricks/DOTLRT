!
!======================================================================
subroutine construct_single_surf_ref() 
!======================================================================
! Constructs surface reflectivities for single point run
!
! History: 
!  12/9/2020 Kevin Schaefer created routine
!----------------------------------------------------------------------
use dotlrt_variables
use profiles
implicit none

! local variables
integer iang  ! (-) stream angle index
!
! stream angles
  do iang = 1, nstream_surf
    surf%theta(iang) = quad_ang(iang) ! (deg) stream angle
  enddo

! general surface characteristics
  surf%ocean_mod = ocean_mod
  surf%nstream = nstream_surf
  surf%temp = atm(1)%temp
  surf%type = 2
  if (trim(surf%ocean_mod) == 'Wilheit') surf%wind = 10.
  call construct_surf_characteristics(surf%href, surf%vref)

return
end subroutine construct_single_surf_ref
!
!======================================================================
subroutine construct_surf_ref_table() 
!======================================================================
! Constructs lookup table of surface reflectivities
!
! History: 
!  11/30/2020 Kevin Schaefer created routine
!  12/9/2020  Kevin Schaefer moved WRF stuff here from construct_surf_characteristics
!----------------------------------------------------------------------
use dotlrt_variables
use profiles
implicit none

! local variables
  integer ichan  ! (-) channel index
  integer ilon   ! (-) longitude index
  integer ilat   ! (-) latitude index

! print message
  print*, '    Create surface reflectance tables'

! general surface characteristics
  surf%ocean_mod = ocean_mod
  surf%nstream = nstream_surf

! loop through all channels and lociations
  do ichan = 1,nchan
    call extract_channel(ichan)
    do ilat = 1,nlat
      do ilon = 1,nlon
        surf%temp = TSK(ilon, ilat)
        surf%type = int(landmask(ilon, ilat))
        if (trim(surf%ocean_mod) == 'Wilheit') surf%wind = wind(ilon, ilat)
        call construct_surf_characteristics(sref_hor(ilon,ilat,ichan,:), sref_ver(ilon,ilat,ichan,:))
      enddo
    enddo
  enddo

return
end subroutine construct_surf_ref_table
!
!======================================================================
subroutine construct_surf_characteristics(surf_ref_hor, surf_ref_ver) 
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
!  11/30/2020 Kevin Schaefer changed reflectances to output variables in arg list
!----------------------------------------------------------------------
use dotlrt_variables
use profiles
implicit none
!
! output variables
  real(8) surf_ref_hor(max_nang) ! (-) horizontal surface reflectivity
  real(8) surf_ref_ver(max_nang) ! (-) vertical surface reflectivity

! local variables
  real(8) slope_var
  real(8) foam_frac ! (-) foam fraction of total area
  integer(4) iang   ! angle index

! common to all surface types
    do iang=1, surf%nstream
      surf%theta(iang)=quad_ang(iang)  
    enddo 

!----------------------------------------------------------------------
! calculate dielectric constant
!----------------------------------------------------------------------
  if (surf%type==2) then ! land
    surf%sal=0.d0
    surf%diel_mod='none'
    surf%diel=dcmplx(0.0d0, 0.0d0)
  elseif (surf%type==1) then ! lake
    surf%sal=sal_lake
    surf%diel_mod='Fresh'
    call dielectric_const_water(channel%lo_freq, surf%temp, surf%sal, surf%type)
  elseif (surf%type==0) then ! ocean
    surf%sal=sal_ocean
    surf%diel_mod='Salt'
    call dielectric_const_water(channel%lo_freq, surf%temp, surf%sal, surf%type)
  endif

!----------------------------------------------------------------------
! Calculate reflectivities 
!----------------------------------------------------------------------
  if (surf%type==2) then ! land
    do iang=1, surf%nstream
      surf_ref_ver(iang)=land_refl 
      surf_ref_hor(iang)=land_refl 
    enddo 

  else ! lake or ocean
    if (trim(ocean_mod)=='Fresnel') then ! Fresnel reflectivity model
      do iang=1, surf%nstream
        call fresnel_refl(iang, surf_ref_hor(iang), surf_ref_ver(iang))
      enddo 

    elseif (trim(ocean_mod)=='Wilheit') then ! Wilheit reflectivity model
      slope_var=0.003d0+0.0048d0*surf%wind 
      if (channel%lo_freq < 35.0d0) slope_var=slope_var*(0.3d0+0.02d0*channel%lo_freq) 
      if (surf%wind < 7.0d0) then
        foam_frac=0.0d0
      else
        foam_frac=(1.0d0-dexp(-channel%lo_freq/7.5d0))*(surf%wind-7.0d0)*0.006d0
      endif
      do  iang=1, surf%nstream 
        ! Compute slope-averaged reflectivity using Gaussian distribution (Kirchoff approximation)
        call kirchoff_ocean_refl(iang, channel%lo_freq, slope_var)
        surf_ref_ver(iang)=surf_ref_ver(iang)*(1.0d0-foam_frac)
        surf_ref_hor(iang)=surf_ref_hor(iang)*(1.0d0-foam_frac)
      enddo

    elseif (trim(ocean_mod)=='ITRA') then ! ITRA reflectivity model
      do iang=1,  surf%nstream
        surf_ref_ver(iang)=0.638d0 - 0.00272d0*channel%lo_freq
        surf_ref_hor(iang)=0.638d0 - 0.00272d0*channel%lo_freq
      enddo 
    endif
  endif

return

end subroutine construct_surf_characteristics
