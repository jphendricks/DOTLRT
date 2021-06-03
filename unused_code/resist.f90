! 
!============================================================
      subroutine Resistance(Tm,Ta,Tc,Tgs,rhoair,umSiB,hSiB, &
       Zlt,z0,zp_disp,g3,corb2,Rbc,Rdc,z2,zwind,zmet, &
       u2,Ra,Rb,Rd,drag,ustar,test,test10)
!============================================================
!   Calculates friction velocity, wind speed at canopy top,
!   aerodynamic resistance and drag.
! Modifications:
!   Tess Krebbs (6/21/00) rewrite original rasite code for SIB
!   Tess Krebbs (7/14/00) Clean up code, convert to fortran 90
!   Kevin Schaefer (8/15/00) more code cleanup, improve variable names
!   Kevin Schaefer (10/7/00) comment code, define variables
!   Kevin Schaefer (10/7/00) moved stability corrections into functions
!   Kevin Schaefer (10/7/00) added 0+-1 W/m**2 limit for stable case
!   Kevin Schaefer (10/7/00) Correct mistakes in RaZ2Zw calculations
!   Kevin Schaefer (2/14/00) Replaced iterative with exact U* solution
!   Kevin Schaefer (2/16/01) Replaced incorrect w/ exact RaHaZ2 solution
!   Kevin Schaefer (2/28/01) U2 calculation based on resistance not wind profile
!   Kevin Schaefer (2/28/01) Add TKE contribution to Um for low wind speed
!   Kevin Schaefer (3/1/01) bounded U2 to prevent ra<0 & Ra>>100 for low wind
! Assumptions:
!   1) rasite uses Monim-Obhukov similarity theory, which does not yeild
!      realistic results for very low wind speeds.
!   2) rasite assumes the input wind lies near the surface layer and above
!      the canopy top.
!   3) paulson psi-coefficients are constrained under unstable
!      conditions to prevent unrealistically low ra values.
! references:
!   Paulson C.A. (1970) 'Mathematical representation of wind
!   and temperature profiles in the unstable atmospheric surface
!   layer', J. Appl. Met., 9, 857-861.
!-------------------------------------------------------------
!
      implicit none
!
! input variables
      real Tm   ! air temperature at zmet (K)
      real Ta   ! Canopy air space temperature (K)
      real Tc   ! canopy temperature (K)
      real Tgs  ! ground surface temperature (K)
      real rhoair  ! air density (kg m^-3)
      real umSiB   ! wind speed at zwind (m s^-1)
      real hSiB ! sensible heat flux from surface (W m^-2)
      real Zlt  ! Leaf Area Index (m**2/m**2)
      real z0   ! roughness length (m)
      real zp_disp ! zero plane displacement (m)
      real g3   ! ratio ra_actual/ra_log-linear for heat between z2 & zmet
      real corb2  ! square of resistance coefficient for Ra (rac**2)
      real Rdc  ! resistance coefficient for Rd (-)
      real Rbc  ! resistance coeficient for Rb (s^.5/m^.5)
      real z2   ! height of canopy top (m)
      real zwind   ! height of momentum measurement (m)
      real zmet    ! height of temperature and humidity measurement (m)
!
! Output variables 
      real u2   ! velocity at top of canopy (m/s)
      real Ra   ! aero resistance between ha and zwind (s/m)
      real Rb      ! aero resistance between canopy and canopy air space (s/m)
      real Rd      ! aero resistance between ground and ha (s/m)
      real drag ! shear stress at canopy top (N/m**2)
      real ustar    ! friction velocity (m s^-1)
      real test, test10(10)   ! test variables
!
! Local variables
      real Um       ! local value of wind speed (UmSiB) (m/s)      
      real Hflux    ! local value sensible heat flux (hSiB) (W m^-2) 
      real zt       ! height of transition layer above canopy (m)
      real Obhukov  ! Obukhov length scale (m)
      real RaHaZ2   ! resistance for heat and momentum between ha and z2 (s/m)
      real RaZ2Zw   ! resistance for heat between z2 and zwind (s/m)
      real ustarHLT ! ustar at Hflux=HLT, stable case, for interpolation (m/s)
      real ustarHTL ! ustar at Hflux=HTL, stable case, for interpolation (m/s)
      real HLT      ! heat flux transition from laminar to turbulent flow (W m^-2)
      real HTL      ! heat flux transition from turbulent to laminar flow (W m^-2)
      real logZw    ! log thickness ratio from (zwind-z_disp)/z0
      real logZ2    ! log thickness ratio from (Z2-z_disp)/z0
      real logZtZ2  ! log thickness ratio from (zt-zp_disp)/(z2-zp_disp)
      real logZwZt  ! log thickness ratio from (zwind-zp_disp)/(zt-zp_disp)
      real logZwZ2  ! log thickness ratio from (zwind-zp_disp)/(z2-zp_disp)
      real Temp     ! dummy variable
      real Factor   ! turbulence scaling factor for resistance
      real UstarN   ! friction velocity for nuetral conditions (m/s)
      real Ustd     ! sqrt(Turb kinetic Energy); std deviation turb velocity (m/s)
      real U2N      ! nuetral value of U2 (m/s)
!
! Local constants
      real vkc      ! von Karman's coefficient
      real cpair    ! specific heat of air (J g^-1 K^-1)
      real pi       ! value of pi
      real grav     ! acceleration of gravity
      real ztz0     ! depth transition layer above canopy in roughness lengths
                    ! also G4: zt=z2+ztz0*z0 (Tuned to Shaw and Pereira model)
      real gamma    ! tuning coefficient from observations (Paulson (1970))
!
! Set local constants
      vkc = 0.41
      cpair = 1004.67
      pi = 3.14159
      grav = 9.8
      ztz0 = 11.875
      gamma = 16.0
!
! calculate transition height
      zt=z2+ztz0*z0
!
! Calculate logrithmic portions of log-linear wind profiles
      logZw=alog((zwind-zp_disp)/z0)
      logZ2=alog((Z2-zp_disp)/z0)
      logZtZ2=alog((zt-zp_disp)/(z2-zp_disp))
      logZwZt=alog((zwind-zp_disp)/(zt-zp_disp))
      logZwZ2=alog((zwind-zp_disp)/(z2-zp_disp))
!
! switch to local values of wind speed and sensible heat flux
      Um=UmSiB
      Hflux=hSiB
!
!-------------------------------------------------------------
! neutral Case (Hflux=0)
!-------------------------------------------------------------
! Calculate neutral friction velocity
      UstarN=vkc*um/logZw
!
!     Calc resistance from z2 to zwind: RaZ2Zw=Ra(z2-zt)+Ra(zt-zwind)
      if (zwind>zt) then
        RaZ2Zw=(G3*logZtZ2+logZwZt)/(vkc*ustarN)
      else
        RaZ2Zw=G3*logZwZ2/(vkc*ustarN)
      endif
!
! Calculate neutral wind speed at canopy top (U2N)
      U2N=Um-RaZ2Zw*ustarN**2+.1
      U2=U2N
!
!-------------------------------------------------------------
! Unstable case: Hflux > 0
!-------------------------------------------------------------
! psi is the Stability correction factor (Paulson, 1970, eqn 2)
!
      if (Hflux.gt.1.) then
!
!       calculate turbulent contribution to total velocity
        ustd=sqrt(2.*grav*Hflux*600./3./Tm/cpair/rhoair)
        Um=Um+Ustd
!
!       Calc friction velocity
        Obhukov=-UstarN**3.*rhoair*cpair*tm/(vkc*grav*Hflux)
        ustar=Um/(logZw-psiUnstabMom(zwind-zp_disp))*vkc
!
!       Calc resistance from z2 to zwind: RaZ2Zw=Ra(z2-zt)+Ra(zt-zwind)
        Obhukov=-ustar**3*rhoair*cpair*tm/(vkc*grav*Hflux)
        if (zwind>zt) then
          RaZ2Zw=G3*(logZtZ2 &
          -psiUnstabHeat(zt-zp_disp)+psiUnstabHeat(z2-zp_disp)) &
          +logZwZt &
          -psiUnstabHeat(zwind-zp_disp)+psiUnstabHeat(zt-zp_disp)

          RaZ2Zw=RaZ2Zw/(vkc*ustar)
        else
          RaZ2Zw=G3*(logZwZ2-psiUnstabHeat(zwind-zp_disp) &
          +psiUnstabHeat(z2-zp_disp))/(vkc*ustar)
        endif
!
!       Calculate wind speed at canopy top (U2)
        U2=Um-RaZ2Zw*ustar**2
        U2=max(u2,Ustd)
        U2=min(U2,Um)
!
!       Calculate resistance from ground to Ha
        Temp=MAX(0.,TGs-TA)                                  
        Factor=SQRT(1.+9.*grav*Temp*Z2/(TGS*U2*U2))         
        Rd=Rdc/(U2*Factor)
      Endif ! unstable case
!
!-------------------------------------------------------------
! Stable case: Hflux < 0
!-------------------------------------------------------------
! psi is the Stability correction factor (Paulson, 1970, eqn 2)
!
      if (Hflux.le.-1.) then  
!       The transition between laminar and turbulent depends on
!       whether the heat flux is increasing or decreasing.
!       This hysteresis results in two critical transition flux values:
!       1) hLT is the transition from laminar to turbulent flow
        hLT=-0.95*tm*rhoair*cpair/(2.0*4.7*grav*(zwind-zp_disp))* &
              (2.0*um/3.0)**3*(vkc/logZw)**2
!       2) hTL is the transition from turbulent to laminar flow (hTL<hLT always)
        hTL=5.*HLT
!
!       laminar flow (Hflux<hTL): Ra fixed at that for HTL to allow flux
        if (Hflux.le.HTL) then
!
!         Calculate friction velocity
!          ustar=vkc*um/(logZw+4.7)
          ustar=-Hflux*(zwind-zp_disp)*vkc*grav/rhoair/cpair/tm
          ustar=ustar**.3333
!
!         Calc resistance from z2 to zwind: RaZ2Zw=Ra(z2-zt)+Ra(zt-zwind)
          if (zmet.gt.zt) then    
            RaZ2Zw=(g3*logZtZ2+logZwZt)/(vkc*ustar)
          else
            RaZ2Zw=g3*logZwZ2/(vkc*ustar)
          endif
!
!         Calculate wind speed at canopy top (U2N)
          U2=Um-RaZ2Zw*ustar**2
          U2=max(u2,U2N)
          U2=min(U2,Um)
!          U2=ustar/vkc/(logz2+4.7)
        endif
!
!       transition flow (hTL<Hflux<hLT)
        if (Hflux.gt.HTL.and.Hflux.lt.HLT) then
!
!         Calculate friction velocity by interpolation between HTL and HLT
!          ustarHTL=vkc*um/(logZw+4.7)
          ustarHTL=-Hflux*(zwind-zp_disp)*vkc*grav/rhoair/cpair/tm
          ustarHTL=ustarHTL**.3333
          Obhukov=-UstarN**3.*rhoair*cpair*tm/(vkc*grav*HLT)
          ustarHLT=Um/(logZw-psiStabMom(zwind-zp_disp))*vkc
          ustar=(Hflux-hTL)/(hLT-hTL)*(ustarHLT-ustarHTL)+ustarHTL
!
!         Calc resistance from z2 to zwind: RaZ2Zw=Ra(z2-zt)+Ra(zt-zwind)
          Obhukov=-Ustar**3.*rhoair*cpair*tm/(vkc*grav*Hflux)
          if (zwind.gt.zt) then    
            RaZ2Zw=g3*(logZtZ2-psiStabHeat(zt-z2)) &
            +logZwZt-psiStabHeat(zwind-zt)
            RaZ2Zw=RaZ2Zw/(vkc*ustar)
          else
            RaZ2Zw=(g3*logZwZ2-psiStabHeat(zwind-z2))/(vkc*ustar)
          endif
!
!         calculate wind speed at canopy top
          U2=Um-RaZ2Zw*ustar**2
          U2=max(u2,U2N)
          U2=min(U2,Um)
        endif
!
!       turbulent flow (Hflux>hLT)
        if (Hflux.gt.HLT) then
!
!         Calculate friction velocity
          Obhukov=-UstarN**3.*rhoair*cpair*tm/(vkc*grav*Hflux)
          ustar=Um/(logZw-psiStabMom(zwind-zp_disp))*vkc
!
!         Calc resistance from z2 to zwind: RaZ2Zw=Ra(z2-zt)+Ra(zt-zwind)
          Obhukov=-ustar**3*rhoair*cpair*tm/(vkc*grav*Hflux)
          if (zwind.gt.zt) then    
            RaZ2Zw=g3*(logZtZ2-psiStabHeat(zt-z2)) &
            +logZwZt-psiStabHeat(zwind-zt)
            RaZ2Zw=RaZ2Zw/(vkc*ustar)
          else
            RaZ2Zw=(g3*logZwZ2-psiStabHeat(zwind-z2))/(vkc*ustar)
          endif
!
!         calculate wind speed at canopy top
          U2=Um-RaZ2Zw*ustar**2
        endif
!
!       Calculate resistance from ground to Ha        
        Rd=Rdc/U2
      endif ! end stable if statement
!
!-------------------------------------------------------------
! Calculate resistances (Ra, Rb, and Rd) and drag
!-------------------------------------------------------------
!
! Calculate resistance from Ha to z2
      RaHaZ2=sqrt(corb2)/u2
!
! Calculate total resistance from Ha to Zwind
      Ra=RaZ2Zw+RaHaZ2
!
! Calculate resistance from canopy to canopy air space
      Factor=ZLT/890.*(abs(TC-TA)*20.0)**0.25
!      Factor=ZLT/890.*(max(0.1,TC-TA)*20.0)**0.25
      Rb=1./(sqrt(U2)/Rbc+Factor)
!
! Calculate drag
      drag=rhoair*ustar*ustar
      test=drag
      test10(1)=drag
!
      return
!
      contains
!
!--------------------------------------------------------------------
      function psiUnstabMom(Dz) result(psi)
!--------------------------------------------------------------------
! Calculates the stability correction factor for momentum
! for the unstable case (Hflux>0) (Paulson, 1970)
!
      real psi ! Stability correction to wind profile (Paulson, 1970, eqn 2)
      real Dz  ! input delta z from zero plane displacment
      real x   ! dummy variable
!
      x=(1-gamma*Dz/Obhukov)**.25
      psi=2.*alog((1.+x)/2.)+alog((1.+x*x)/2.) &
                -2.*atan(x)+pi/2.
      psi=amin1(logZw*0.75,psi)
!
      end function psiUnstabMom
!
!--------------------------------------------------------------------
      function psiStabMom(Dz) result(psi)
!--------------------------------------------------------------------
! Calculates the stability correction factor for Momentum
! for the stable case (Hflux<=0)
!
      real psi ! Stability correction to wind profile (Paulson, 1970, eqn 2)
      real Dz  ! input delta z from zero plane displacment
!
      psi=-4.7*Dz/Obhukov
      psi=amax1(-4.7,psi)
!
      end function psiStabMom
!
!--------------------------------------------------------------------
      function psiUnstabHeat(Dz) result(psi)
!--------------------------------------------------------------------
! Calculates the stability correction factor for heat
! for the unstable case (Hflux>0) (Paulson, 1970)
!
      real psi ! Stability correction to wind profile (Paulson, 1970, eqn 2)
      real Dz  ! input delta z from zero plane displacment
      real x   ! dummy variable
!
      x=(1-gamma*Dz/Obhukov)**.25
      psi=2*alog((1+x*x)/2.)
!
      end function psiUnstabHeat
!
!--------------------------------------------------------------------
      function psiStabHeat(Dz) result(psi)
!--------------------------------------------------------------------
! Calculates the stability correction factor for heat
! for the stable case (Hflux<=0)
!
      real psi ! Stability correction to wind profile (Paulson, 1970, eqn 2)
      real Dz  ! input delta z from zero plane displacment
!
      psi=-4.7*Dz/Obhukov
      psi=amax1(-4.7,psi)
!
      end function psiStabHeat
!
      end
