!=======================================================================
SUBROUTINE dtb94 ()
!=======================================================================
! This subroutine further prepares data for tb94 which calls CORE95.F90 which
! then calculates brightness temperature and the radiative Jacobian.
!  Below, that Jacobian is called the incremental profiles. 
!-----------------------------------------------------------------------
! Modification 04/24/03 : DTB9.F90 -> DTB93.F90
!   This subroutine excepts and transfers to TB94.F90 (former TB53.F90)
!   parameter _eps_err.
!  
!-----------------------------------------------------------------------
! Modification 12/17/02 : DTB9.F90 -> DTB93.F90
!   This subroutine has two additional input parameters (arrays):
!                     _dphaseff and _dphasefb
!  Those parameters are also interpolated within layers similarly to
!  the others.
!  Note, that incremental profiles are normalized by layer thickness
!   9/26/2020 Kevin Schaefer deleted unused variables
!-----------------------------------------------------------------------
!  This subroutine further prepares data for tb94.f90 which calls CORE95.F90 which
!  then calculates brightness temperature and (Tb) profiles
!  as well as incremental profiles with respect to variations of
!  gas absorption + hydrometeor absorption, scattering by cloud, and temperature
!  at different levels.
!  
!  Upwelling and downwelling brightness temperatures for different
!  altitudes are stored in file TB.RES. Incremental profiles are
!  stored in file DTB.RES. This file has the following structure:
!
!   altitude(ilr) , dTb_pl(obs_lev)/d(beta(ilr))
!                 , dTb_mn(obs_lev)/d(beta(ilr))
!
!  Here:
!
!         d/d(beta ) = scat_cloud(ilr)* d/d(scat_cloud(ilr))
!
!  where _ilr - is a number of original level at which absobtion 
!  by gases (which could be used also as a proxi for particles) or
!  scattering by cloud ( n0 in notation of paper) is varied.
!
!  Note, that there is a possibility to insert into original levels
!  _m1 "pseudo-levels". When m1=1 original levels are in fact
!  reproduced. Those pseudo-levels are uniformly spaced within original
!  layers, and temperature, gas and particle absorbtion, phase functions
!  are linearly iinterpolated. In the first layer no interpolation is
!  applied, and all parameters are kept constant corresponding to the
!  first level (the upper boundary of the first layer).
!
!
!                             INPUT
!   nlr       ( I) - total number of layers
!   m1        ( I) - number of pseudo-layers within each layer
!                    (when m1=1 original layers are reproduced)
!  obs_lev    ( I) - observation level at which Tb is observed. 
!                    Derivatives of this temperature with respect to
!                    absorbtion and scattering coefficients are stored
!                    in file DTB.RES. If _obs_lev > nlr it will be set
!                    obs_lev=nlr
!   nangover2      ( I) - number of angles
!   i0        ( I) - number of observation angle at which Tb is observed
!   nvar      ( I) - number of parameters subjected to variation
!  t_srf      (DP) - temperature of the surface
!  t_cb       (DP) - temperature of the cosmic background
!     altitude(nlr+1)  (DP) - an array with altitudes. Surface level
!                             corresponds _altitude(1) - this is 
!                             different from all other values below
!               
!   temperature (nlr)  (DP) - thermodynamic temperatures at the upper
!                             boundaries of layers. Thus, temperature(1)
!                             corresponds to _altitude(2) .
!        abs_O2 (nlr)  (DP) - absorbtions by O2 within layers
!       abs_H2O (nlr)  (DP) - absorbtions by H2O within layers
!     abs_cloud (nlr)  (DP) - absorbtion by cloud (particles) within 
!                             layers
!    scat_cloud (nlr)  (DP) - scattering by clouds within layers
!                             (parameter _n0 in theoretical paper)
!           cs (nangover2)  (DP) - array of cosines of incidence angles
!  surf_reflec (nangover2)  (DP) - array of surface reflectivities
!   phaseff(nlr,nangover2,nangover2)   (DP) - array of phase functions - 
!                                   forward scattering
!   phasefb(nlr,nangover2,nangover2)   (DP) - array of phase functions - 
!                                   backward scattering
!  dphaseff(nlr,nangover2,nangover2)   (DP) - array of derivatives of phase 
!                                   functions  forward scattering
!  dphasefb(nlr,nangover2,nangover2)   (DP) - array of derivatives of phase 
!                                   functions backward scattering
!
!                           OUTPUT
!  tb_pl1(0:nlr*m1,nangover2)  -   (DP)   array of upwelling brightness
!                                    temperatures at all levels
!                                    (including pseudo-levels)
!  tb_mn1(0:nlr*m1,nangover2)  -   (DP)   array of downwelling brightness
!                                    temperatures at all levels
!                                    (including pseudo-levels)
! dtb_pl1(0:nlr*m1,nangover2,nvar) - (DP) array of derivatives of upwelling
!                                    brightness temperatures at all
!                                    levels (including pseudo-levels)
!                                    ( see above ) for all variations
! dtb_mn1(0:nlr*m1,nangover2,nvar) - (DP) array of derivatives of downwelling
!                                    brightness temperatures at all
!                                    levels (including pseudo-levels)
!                                    for all variations
!-------------------------------------------------
!   ivar=1 : variation of scattering by clouds
!   ivar=2 : variation of (gas_abs+cloud_abs)
!   ivar=3 : variation of temperature
!   ivar=4 : variation of HG assymetry parameter
!------------------------------------------------- 
! eps_err                     - (DP) max allowed relative error in brightness
!                                    temperature allowed    
!-----------------------------------------------------------------------

use dotlrt_variables
IMPLICIT NONE

INTEGER :: i,j,k,k1,k2, ilr, iphase
real(8), DIMENSION(nlr+1)    :: h
real(8) k1m1, m1k1m1, k2m1, m1k2m1
EXTERNAL tb94

  DO ilr=1,nlr
    h(ilr)=altitude(ilr+1)-altitude(ilr)
  END DO

  !----  the first layer: no interpolation 
  altitude1(1)=altitude(1)
  k=0
  DO k1=1,m1
                  k = k+1

     k1m1 = dble(k1) / dble(m1)
     m1k1m1 = 1.0d0 - k1m1

     altitude1(k+1) = m1k1m1 * altitude(1) + k1m1 * altitude(2)

              h1(k) =  altitude1(k+1) - altitude1(k)
    temperature1(k) = temperature(1)
        abs_O2_1(k) =      abs_O2(1)
       abs_H2O_1(k) =     abs_H2O(1)
      abs_cloud1(k) =   abs_cloud(1)
     scat_cloud1(k) =  scat_cloud(1)
         al_gas1(k) =    abs_O2_1(1) + abs_H2O_1(1)
      abs_total1(k) =     al_gas1(k) + abs_cloud1(k)
    DO i=1,nangover2
      DO j=1,nangover2
        phaseff1(k,i,j) = phaseff(1,i,j)
        phasefb1(k,i,j) = phasefb(1,i,j)
        phaseff1_sc(k,i,j) = phaseff_sc(1,i,j)
        phasefb1_sc(k,i,j) = phasefb_sc(1,i,j)
        do iphase = 1, nphase
          dphasefb1_g(k,i,j,iphase)  = dphasefb_g(1,i,j,iphase)
          dphaseff1_g(k,i,j,iphase)  = dphaseff_g(1,i,j,iphase)
          dphasefb1_sc(k,i,j,iphase) = dphasefb_sc(1,i,j,iphase)
          dphaseff1_sc(k,i,j,iphase) = dphaseff_sc(1,i,j,iphase)
        end do
      END DO
    END DO

   END DO

  !----  higher layers: linear interpolation 
  DO ilr=2,nlr
    DO k1=1,m1
      k=k+1

      k1m1 = dble(k1) / dble(m1)
      m1k1m1 = 1.0d0 - k1m1

      altitude1(k+1) = m1k1m1 * altitude(ilr) + k1m1 * altitude(ilr+1)


      !----------------------------------------------------------------------
      !-----   linear iterpolation: if one sets k2=m1, then all parameters
      !-----   at pseudo-levels will be the same as at the basic levels. This
      !-----   can be used for checking purposes. The choise k2=k1 leads to
      !-----   linear interpolation.
      !
      !         k2=m1
      !----------------------------------------------------------------------
      k2=k1                          !   linear interpolation in effect
      h1(k) =   altitude1(k+1)-altitude1(k)

      k2m1 = dble(k2) / dble(m1)
      m1k2m1 = 1.0d0 - k2m1

      temperature1(k) = m1k2m1 * temperature(ilr-1) + k2m1 * temperature(ilr)
          abs_O2_1(k) = m1k2m1 *      abs_O2(ilr-1) + k2m1 *      abs_O2(ilr)
         abs_H2O_1(k) = m1k2m1 *     abs_H2O(ilr-1) + k2m1 *     abs_H2O(ilr)
        abs_cloud1(k) = m1k2m1 *   abs_cloud(ilr-1) + k2m1 *   abs_cloud(ilr)
       scat_cloud1(k) = m1k2m1 *  scat_cloud(ilr-1) + k2m1 *  scat_cloud(ilr)

           al_gas1(k) = abs_O2_1(k)+abs_H2O_1(k)
!!!! blw        abs_cloud1(k) = DMAX1( abs_cloud1(k), 1.0d-10 )  ! vacuum protection
        abs_total1(k) = al_gas1(k) + abs_cloud1(k)
      DO i=1,nangover2
        DO j=1,nangover2
          phaseff1(k,i,j)    = m1k2m1 * phaseff(ilr-1,i,j)    + k2m1 * phaseff(ilr,i,j)
          phasefb1(k,i,j)    = m1k2m1 * phasefb(ilr-1,i,j)    + k2m1 * phasefb(ilr,i,j)
          phaseff1_sc(k,i,j) = m1k2m1 * phaseff_sc(ilr-1,i,j) + k2m1 * phaseff_sc(ilr,i,j)
          phasefb1_sc(k,i,j) = m1k2m1 * phasefb_sc(ilr-1,i,j) + k2m1 * phasefb_sc(ilr,i,j)
          do iphase = 1, nphase
            dphasefb1_g(k,i,j,iphase)  = m1k2m1 * dphasefb_g(ilr-1,i,j,iphase)  &
                                                  +   k2m1 * dphasefb_g(ilr,i,j,iphase)
            dphaseff1_g(k,i,j,iphase)  = m1k2m1 * dphaseff_g(ilr-1,i,j,iphase)  &
                                                  +   k2m1 * dphaseff_g(ilr,i,j,iphase)
            dphasefb1_sc(k,i,j,iphase) = m1k2m1 * dphasefb_sc(ilr-1,i,j,iphase) &
                                                  +   k2m1 * dphasefb_sc(ilr,i,j,iphase)
            dphaseff1_sc(k,i,j,iphase) = m1k2m1 * dphaseff_sc(ilr-1,i,j,iphase) &
                                                  +   k2m1 * dphaseff_sc(ilr,i,j,iphase)
          end do ! iphase
        END DO
      END DO
    END DO ! k1
  END DO ! ilr

!-------  end of forming interpolated values

 CALL  tb94()

END SUBROUTINE dtb94
