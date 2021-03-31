!**********************************************************************
! NAME: construct_surf
!
! Purpose:Constructs model of ocean surface from several options listed in the code below.
!
! Language: FORTRAN 90 Portland Group Compiler on Red Hat Linux
!
! Modules:
! Uses variables.f90
!      types_kind
! Creation History: Original from NOAA ETL / ET1 Pascal GMRT version
!                   provided by Bob Weber
!
! Converted from original Pascal by:
!
! Ronald Richter NOAA/ETL SET April 2003  ron.richter@noaa.gov
!
!______________________________________________________________________
!
SUBROUTINE construct_surf(surf_type, surf_type1, refl,  freq, sal, wind_speed) 
  use variables
  implicit none
  character(1) surf_type, surf_type1
  integer(4) i, j
  real(8) refl, freq, sal, wind_speed, slope_var, foam_frac 

! Set up surf_inp structure, defined in variables.mod
!
          surf_inp%surf_type = surf_type1 
          surf_inp%surf_subtype = surf_type
          surf_inp%param1 = refl 
          surf_inp%freq = freq 
          surf_inp%param2 = sal 
          surf_inp%param3 = wind_speed 

!            {surface temperature in K}
       SELECT CASE (surf_type1)
              CASE ('L')
                 DO i=1, nangover2 ! surf_inp%num_angles
                     surf_inp%theta(i)=quad_angle_array(i) ! (i-1)*90.0d0/dble(surf_inp%num_angles-1) 
                     surf_inp%vr(i)=refl 
                     surf_inp%hr(i)=refl 
                 END DO 

              CASE ('O')
                SELECT CASE (surf_type) 
                   CASE ('F')
! {Fresnel reflectivity model}
                            DO i=1, nangover2 ! surf_inp%num_angles
                               surf_inp%theta(i)=quad_angle_array(i) ! (i-1)*90.0d0/dble(surf_inp%num_angles-1) 
                               call ocean_refl(freq,surf_inp%surf_temp,sal,surf_inp%theta(i),surf_inp%vr(i),surf_inp%hr(i)) 
                            END DO 
                            
! {Wilheit ocean reflectivity model
                  CASE ('W')           
                                    slope_var=0.003d0+0.0048d0*wind_speed 
                                    IF (freq .LT. 35.0d0) slope_var=slope_var*(0.3d0+0.02d0*freq) 
                                    IF (wind_speed .LT. 7.0d0) THEN
                                        foam_frac=0.0d0
                                    ELSE
                                        foam_frac=(1.0d0-dexp(-freq/7.5d0))*(wind_speed-7.0d0)*0.006d0
                                    END IF
                                    DO  i=1, nangover2 ! surf_inp%num_angles
                                            surf_inp%theta(i)=quad_angle_array(i) ! (i-1)*90.0d0/dble(surf_inp%num_angles-1)
!{Compute slope-averaged reflectivity using Gaussian distribution (Kirchoff approximation)}
                                            call kirchoff_ocean_refl(freq,surf_inp%surf_temp,sal,surf_inp%theta(i),slope_var,surf_inp%vr(i),surf_inp%hr(i))
                                            surf_inp%vr(i)=surf_inp%vr(i)*(1.0d0-foam_frac)
                                            surf_inp%hr(i)=surf_inp%hr(i)*(1.0d0-foam_frac)
                                    END DO
                  CASE ('I')
! {ITRA ocean reflectivity model}
                            refl=0.638d0 - 0.00272d0*freq 
                            DO i=1,  nangover2 ! surf_inp%num_angles
!                               BEGIN
                                    surf_inp%theta(i)=quad_angle_array(i) ! (i-1)*90.0d0/dble(surf_inp%num_angles-1) 
                                    surf_inp%vr(i)=refl
                                    surf_inp%hr(i)=refl
                            END DO 
!           
                  END SELECT ! subtype1
                END SELECT ! subtype 
            surf_inp%inf=.true. 
    return
END SUBROUTINE construct_surf