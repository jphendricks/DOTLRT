!======================================================================
SUBROUTINE hydrometeor_master_5ph_d( freq, phase, temp, p, q, k0, a0,       &
                                     hab, hsc, g,                           &
                                     dhab, dhsc, dg,                        &
                                     a0_const, testvar1, testvar2, testvar3 )
!======================================================================
! This file contains subroutine HYDROMETEOR_MASTER_5PH_d.F90
! and all slave subroutines called by it
! NOTE:  it is used here, that (p=0, q=0) or (p=0, q=1) only
!----------------------------------------------------------------------
! HYDROMETEOR_MASTER_5PH_d.F90 was created by modifying subroutine HYDROMETEOR_EXT_5PH.F90
! From CALCPR0FILE the size distributions only let either k0 or a0 vary, but not both.
! For phases 1 (cloud liquid), 3 (ice),    k0 varies with water density, but a0 is fixed.
! For phases 2(rain), 4(snow), 5(graupel), a0 varies with water density, but k0 is fixed.
! This master routine is for a general size distribution for which both k0 and a0 
! vary with water density.
! Following this subroutine and in this file, is a subroutine specialized for phase = 1 & 3
! as well as another subroutine specialized for phase = 2,4,& 5 (which is not presently used, but is saved)
! Differentiation was redone by Reg Hill 8/12/03; testing and corrections were finished in Sept. 03
! The logical variable "is_a0_constant" was inserted by Reg Hill 8/28/03 to identify cases for
! which a0 is constant because a great savings of computer time is obtained for that case. 
!-----------------------------------------------------------------------
!   This is a version of the original subroutine which calculates along
!   with parameters _hab, _hsc, _g vectors _dhab(1-3), _dhsc(1-3), _dg(1-3)
!   which are derivatives of the appropriate values with respect to
!   parameters _temp, _k0, _a0, appropriately
!   [so that, e.g.   _dhsc(1)=d(_hsc)/d(_temp) 
!                  , _dhsc(2)=d(_hsc)/d(_k0  )
!                  , _dhsc(3)=d(_hsc)/d(_a0  ) ].
! Warning file with reference number _outfile should be opened outside
! of the subroutine; in the present version it is set _outfile=5.
!-----------------------------------------------------------------------
! History:
!  9/26/2020 Kevin Schaefer deleted unused variables
!  2/2/2021  Kevin Schaefer cleaned up code, deleted unused code
!-----------------------------------------------------------------------
IMPLICIT NONE

! Inputs
integer,intent(in) :: phase
real(8),intent(in) :: freq
real(8),intent(in) :: temp
real(8),intent(in) :: p
real(8),intent(in) :: q
real(8),intent(in) :: k0
real(8),intent(in) :: a0
logical, intent(in) :: a0_const ! this single value corresponds to one atmospheric level and one phase

! output
 real(8), intent(out) :: hab     ! (-) cloud absorption coefficient
 real(8), intent(out) :: hsc     ! (-) cloud scattering coefficient
 real(8), intent(out) :: g       ! (-) cloud asymytry factor
 real(8), intent(out) :: dhab(3) ! (?) derivative absorption wrt to temp, k0, a0
 real(8), intent(out) :: dhsc(3) ! (?) derivative scattering wrt to temp, k0, a0
 real(8), intent(out) :: dg(3)   ! (?) derivative asymytry factor wrt to temp, k0, a0
! nvars x nhydro x nchan x dx = 12 x 5 x {15,24}{FY,GEMS}
! local variables
 integer, parameter :: m1=20
 integer, parameter :: m2=10
 integer, parameter :: m3=8
 real(8), parameter :: raythresh=0.1d0
 integer i
 integer numpanels
 logical fixnumpanels
 real(8) pi
 real(8) c
 real(8) speclength
 real(8) f
 real(8) meanradrel
 real(8) sigmarel
 real(8) ka
 real(8) hex
 real(8) etasc 
 real(8) etaex
 real(8) assymetry 
 real(8) dka_da0
 real(8) dhsc_dt 
 real(8) dhex_dt 
 real(8) dg_dt
 real(8) dhsc_da0
 real(8) dhex_da0
 real(8) dg_da0 
 real(8) detasc_dt 
 real(8) detaex_dt 
 real(8) dassymetry_dt
 real(8) detasc_da0
 real(8) detaex_da0
 real(8) dassymetry_da0
 real(8) weight 
 real(8) alph
 real(8) alph1
 real(8) dalph
 real(8) dalph1
 real(8) alphu
 real(8) w
 real(8) dw_dt
 real(8) x
 real(8) dx
 real(8) ddalph_da0
 real(8) dalph1_da0
 real(8) dalph_da0
 real(8) dweight_da0
 real(8) dc_da0
 real(8) dw_da0
 real(8) df_da0
 real(8) dw_dk0
 real(8) df_dk0
 double complex epsil
 double complex depsil_dt
 double complex index
 double complex dindex_dt
 real(8), parameter :: t_frz=273.15d0 ! (K) freezing point of water
  real(8) testvar1 ! test diagnostic variable 1
  real(8) testvar2 ! test diagnostic variable 2
  real(8) testvar3 ! test diagnostic variable 3

! initialize some variables 
  numpanels = 0
  fixnumpanels = .false.
  pi=4.0d0*DATAN(1.d0)
  hab=0.d0
  hsc=0.d0
  g=0.d0
  dhab=0.d0
  dhsc=0.d0
  dg=0.d0

! This subroutine avoids calculating derivatives wrt a0 when a0 is a constant.
  if (a0_const) then
    CALL hydrometeor_a0_is_constant_d(freq,phase,temp,p,q,k0,a0, &
                               hab, hsc, g, dhab,dhsc,dg,numpanels,fixnumpanels)
    return
  end if

! Main part
  IF(k0 > 0.d0) THEN

! calculate dielectric constant
    if (phase == 1) then ! clw
      call h2o_liquid_dielectric(freq, temp, epsil, depsil_dt)
      speclength=1.d0

    elseif (phase == 2) then ! rain
      call h2o_liquid_dielectric(freq, temp, epsil, depsil_dt)
      speclength=1.d0

    elseif (phase == 3) then ! ice
      call h2o_ice_dielectric(freq, dmin1(t_frz,temp), epsil, depsil_dt) 
      speclength=1.028d0   

    elseif (phase == 4) then ! snow
      call h2o_mixed_dielectric(freq, temp, phase, epsil, depsil_dt, testvar1, testvar2, testvar3) 
      speclength=1.028d0   

    elseif (phase == 5) then ! graupel
      call h2o_mixed_dielectric(freq, temp, phase, epsil, depsil_dt, testvar1, testvar2, testvar3) 
      speclength=1.028d0   
    endif

    IF( DIMAG(epsil) > 0.0d0 .AND. DREAL(epsil) /= 1.0d0 ) THEN  
      IF( p >= 0.d0 .AND. q > 0.d0 ) THEN
        meanradrel=1.d0
        sigmarel=1.d0
      ELSE
        meanradrel=1.d0
        sigmarel=0.d0
      END IF
      dka_da0=2.d0*pi*freq/3.0d2*meanradrel*speclength
      ka=dka_da0*a0
 
! absorption & scattering 
      IF( ka*(1.0d0+DBLE(m2)*sigmarel) > raythresh/10.0d0 ) THEN  
        index=DCONJG(SQRT(epsil))
        dindex_dt=DCONJG(depsil_dt)/index/2.0d0

! averaging over distribution 
        IF( p >= 0.d0 .AND. q > 0.0d0 ) THEN
          dalph = 1.d0/dble(m1)
          ddalph_da0 = 0.d0 
          IF( ka > 1.d0 ) THEN
            IF( ka > dble(m3) ) THEN
              dalph= dalph/dble(m3)
              ddalph_da0 = 0.d0  
            ELSE
              dalph= dalph/ka
              ddalph_da0 = -dalph/ka*dka_da0 !other derivatives are zero
            END IF
          END IF
          alph1 = 1.0d0-m2*sigmarel
          dalph1 = 0.d0
          IF( alph1 <= 0.d0 ) THEN
            alph1    =  dalph
            dalph1_da0 = ddalph_da0  !other derivatives are zero
          END IF
          alphu=1.0d0+m2*sigmarel
          IF (fixnumpanels .eqv. .false.) THEN
            numpanels=INT2((alphu-alph1)/dalph+1.5d0)
          END IF 
          dalph= (alphu -  alph1)/numpanels
          ddalph_da0 = (- dalph1_da0)/numpanels !other derivatives are zero
          hex=0.d0
          hsc=0.d0
          g=0.d0

          dhex_dt=0.d0
          dhsc_dt=0.d0
          dg_dt=0.d0

          dhex_da0=0.d0
          dhsc_da0=0.d0
          dg_da0=0.d0

          alph    =  alph1
          dalph_da0 = dalph1_da0 

! trapezoidal integration
          DO i=0,numpanels ! trapezoidal integration
            x= ka*alph
            dx =  dka_da0*alph+ka*dalph_da0
            CALL ksph_d( x, index, dx , dindex_dt, &
                  etasc, etaex, assymetry, &
                  detasc_dt , detaex_dt , dassymetry_dt, &
                  detasc_da0, detaex_da0, dassymetry_da0 )

            weight = DEXP((p+2)*DLOG(alph) - DEXP(q*DLOG(meanradrel*alph)))
            dweight_da0 = weight*((p+2)*dalph_da0/alph - DEXP(q*DLOG(meanradrel*alph)) *q*dalph_da0/alph  ) 

            if( i==0 .OR. i==numpanels ) then
              etasc = etasc/2.0d0
              etaex = etaex/2.0d0
              assymetry = assymetry/2.0d0
              detasc_dt = detasc_dt/2.0d0
              detaex_dt = detaex_dt/2.0d0
              dassymetry_dt = dassymetry_dt/2.0d0 
              detasc_da0 = detasc_da0/2.0d0
              detaex_da0 = detaex_da0/2.0d0
              dassymetry_da0 = dassymetry_da0/2.0d0
            end if
            hsc=hsc + etasc*weight
            hex=hex + etaex*weight
            g= g+assymetry*etasc*weight
            dhsc_dt = dhsc_dt + detasc_dt*weight
            dhex_dt = dhex_dt + detaex_dt*weight 
            dg_dt = dg_dt+ (dassymetry_dt*etasc + assymetry*detasc_dt)*weight
            dhsc_da0= dhsc_da0 + detasc_da0*weight + etasc*dweight_da0 
            dhex_da0= dhex_da0 + detaex_da0*weight + etaex*dweight_da0 
            dg_da0= dg_da0 + etasc*weight*dassymetry_da0 + assymetry*weight*detasc_da0 + assymetry*etasc*dweight_da0 
            alph =  alph + dalph
            dalph_da0 = dalph_da0 + ddalph_da0
          END DO

          IF( hsc /= 0.d0 ) THEN
            g = g/hsc
            dg_dt = dg_dt /hsc - g/hsc*dhsc_dt ! note that g here is the previous g divided by hsc
            dg_da0= dg_da0/hsc - g/hsc*dhsc_da0  ! ditto
          ELSE
          IF( g /= 0.d0 ) THEN
            g = 0.d0
            dg_dt = 0.d0
            dg_da0 = 0.d0      
          END IF
        END IF

        c = dalph*DEXP((p+3)*DLOG(meanradrel))
        dc_da0 = ddalph_da0*c/dalph   
        dhsc_dt =c*dhsc_dt
        dhex_dt =c*dhex_dt
        dhsc_da0 = c*dhsc_da0 + dc_da0*hsc 
        dhex_da0 = c*dhex_da0 + dc_da0*hex 
        hsc = c*hsc  
        hex = c*hex
   
        ELSE    ! averaging over distribution
          CALL ksph_d( ka, index, dka_da0, dindex_dt, hsc, hex, g, dhsc_dt, dhex_dt, dg_dt, dhsc_da0, dhex_da0, dg_da0)
        END IF  ! absorption & scattering

        dg(1) = dg_dt    ! from either part of IF, this is where to define dg
        dg(3) = dg_da0
        dg(2) = 0.d0  ! the k0 derivative 
        c=pi*1.0d-3 * speclength**2
        w =  c*k0*a0**3
        dw_dk0 =  w/k0      
        dw_da0 =  w/a0*3.d0 

! calculate derivatives first, since _hsc=w*hsc, _hab=w*hex-hsc will be re-defined
        dhsc(1)=w*dhsc_dt
        dhab(1)=w*dhex_dt - dhsc(1) ! why introduce arrays here?
        dhsc(2)= dw_dk0 *hsc 
        dhab(2)= dw_dk0 *hex - dhsc(2)
        dhsc(3)= dw_da0*hsc+w*dhsc_da0 
        dhab(3)= dw_da0*hex + w*dhex_da0 - dhsc(3)    
        hsc=w*hsc     
        hab=w*hex-hsc 
        
! Rayleigh absorption without scattering
      ELSE   
        IF( a0 > 0.d0 ) THEN  
          IF( p >= 0.d0 .AND. q > 0.d0 ) THEN
            df_dk0= 8.0d0*pi*1.0d-9 *a0**4   ! so far p=0, q=1
          ELSE
            df_dk0= 4.0d0*pi/3.0d0*1.0d-9*a0**4 
          END IF
          f=df_dk0*k0 ! f = 8.0d0*pi*1.0d-9 * k0 *a0**4 ELSE f= 4.0d0*pi/3.0d0*1.0d-9 * k0 *a0**4
          df_da0=4.d0*f/a0 
        ELSE                  ! a0 < 0
          f=0.d0 !all values of array set to zero
          df_dk0= 0.d0 ! ditto
          df_da0= 0.d0 ! ditto
        END IF

        epsil=epsil + 2.0d0 !!!!!why?
        hsc=0.d0
        g=0.d0
        dhsc=0.d0 ! arrays set to zero 
        dg=0.d0 ! ditto
        c = 0.1885d6*freq*speclength**3
        w =  c*DIMAG(-1.0d0/epsil)
        dw_dt =  c*DIMAG(depsil_dt/epsil**2)
        hab=w*f
        dhab(1)= dw_dt*f 
        dhab(2)= w*df_dk0
        dhab(3)= w*df_da0
      END IF ! end of abs & scattering / Rayleigh absorption without scattering
   END IF ! DIMAG(epsil) > 0.d0 .AND. DREAL(epsil) /= 1.d0
 END IF ! end of the main part

! Check calculated quantities against constraints
  IF(hab < 0.d0) hab=0.d0
  IF(hsc < 0.d0) hsc=0.d0
  IF(g < -1.d0) g=-1.d0
  IF(g >  1.d0) g=1.d0

END SUBROUTINE hydrometeor_master_5ph_d

!======================================================================
SUBROUTINE hydrometeor_a0_is_constant_d(freq,phase,temp,p,q,k0,a0, &
                               hab, hsc, g,   &
                               dhab,dhsc,dg,numpanels,fixnumpanels)
!======================================================================
! This file contains subroutine HYDROMETEOR_a0_is_constant_d.F90
! and all slave subroutines called by it
!     NOTE:  it is used here, that (p=0, q=0) or (p=0, q=1) only
!From CALCPR0FILE the size distributions only let either k0 or a0 vary, but not both.
!For phases 1 (cloud liquid), 3 (ice),    k0 varies with water density, but a0 is fixed.
!For phases 2(rain), 4(snow), 5(graupel), a0 varies with water density, but k0 is fixed.
!This cloud liquid and ice subroutine therefore contains no derivatives wrt a0 
!  Differentiation was completed by Reg Hill on 8/12/03.   
!-----------------------------------------------------------------------
!   This is a version of the original subroutine which calculates along
!   with parameters _hab, _hsc, _g vectors _dhab(1-3), _dhsc(1-3), _dg(1-3)
!   which are derivatives of the appropriate values with respect to
!   parameters _temp, _k0, _a0, appropriately
!   [so that, e.g.   _dhsc(1)=d(_hsc)/d(_temp) 
!                  , _dhsc(2)=d(_hsc)/d(_k0  )
!                  , _dhsc(3)=d(_hsc)/d(_a0  ) ].
!NOTE: dg_dk0 =0 , 
! Warning file with reference number _outfile should be opened outside
! of the subroutine; in the present version it is set _outfile=5.
!
! Modifications
! 2/2/2021 Kevin Schaefer cleaned up code, deleted unused code
!-----------------------------------------------------------------------

IMPLICIT NONE

real(8) :: freq, temp,p,q,k0,a0
real(8),               INTENT(OUT) ::  hab,  hsc,  g
real(8), DIMENSION(3), INTENT(OUT) :: dhab, dhsc, dg
INTEGER         ,               INTENT( IN) :: phase
INTEGER         , PARAMETER :: m1=20, m2=10, m3=8
real(8), PARAMETER :: raythresh=0.1d0
INTEGER :: i, numpanels
real(8) ::  pi, c , speclength, f, meanradrel, sigmarel, ka   &
                   , hex, etasc , etaex, assymetry                     &
                   ,  dhsc_dt , dhex_dt , dg_dt                        &
                   , detasc_dt , detaex_dt , dassymetry_dt             &
                   , detasc_da0, detaex_da0, dassymetry_da0            &
                   , weight , alph, alph1, dalph, alphu, w , dw_dt
real(8) :: dw_dk0, dhsc_da0, dhex_da0, dg_da0
DOUBLE COMPLEX   ::   epsil, depsil_dt, index, dindex_dt
  real(8) testvar1 ! test diagnostic variable 1
  real(8) testvar2 ! test diagnostic variable 2
  real(8) testvar3 ! test diagnostic variable 3
LOGICAL fixnumpanels
real(8) ::   df_dk0
   pi=4.0d0*DATAN(1.d0)
   hab=0.d0
   hsc=0.d0
     g=0.d0
  dhab=0.d0
  dhsc=0.d0
    dg=0.d0

!------------------------------------------------------
 IF( k0 > 0.d0 ) THEN            !        the main part
!------------------------------------------------------

   IF( phase < 3) THEN     ! liquid case

    CALL h2o_liquid_dielectric(freq,temp,epsil,depsil_dt)

    speclength=1.d0

   ELSE 

    IF( phase == 3)  THEN                                                
     CALL h2o_ice_dielectric( freq, DMIN1(273.15d0,temp), epsil, depsil_dt) 
    ELSE    ! ice/liquid mixture
    CALL h2o_mixed_dielectric(freq, temp, phase, epsil, depsil_dt, testvar1, testvar2, testvar3)
    END IF

    speclength=1.028d0   
                                        
   END IF                  ! end of dielectric constant calculation

!-----------------------------------------------------------------------
   IF( DIMAG(epsil) > 0.0d0 .AND. DREAL(epsil) /= 1.0d0 ) THEN  
!-----------------------------------------------------------------------

     IF( p >= 0.d0 .AND. q > 0.d0 ) THEN 
       meanradrel=1.d0  
         sigmarel=1.d0   
     ELSE
       meanradrel=1.d0
         sigmarel=0.d0
     END IF

      ka = 2.d0*pi*freq/3.0d2*meanradrel*speclength*a0

!--------------------------------------------------------------------------
   IF( ka*(1.0d0+DBLE(m2)*sigmarel) > raythresh/10.0d0 ) THEN  ! absorption & scattering
!--------------------------------------------------------------------------

     index   =DCONJG(SQRT(epsil))

    dindex_dt=DCONJG(depsil_dt)/index/2.0d0

   IF( p >= 0.d0 .AND. q > 0.0d0 ) THEN     ! averaging over distribution

     dalph      = 1.d0/dble(m1)

     IF( ka > 1.d0 ) THEN
      IF( ka > dble(m3) ) THEN
       dalph     = dalph/dble(m3)
      ELSE
       dalph     =  dalph/ka
      END IF
     END IF

     alph1=1.0d0-m2*sigmarel

     IF( alph1 <= 0.d0 ) THEN
        alph1    =  dalph
        END IF
        alphu=1.0d0+m2*sigmarel
IF (fixnumpanels .eqv. .false.) THEN
     numpanels=INT2((alphu-alph1)/dalph+1.5d0) ! IIDINT((alphu-alph1)/dalph+1.5d0) !?? differentiate an integer?? NO
END IF 
    
        dalph     = (alphu -  alph1    )/numpanels

           hex=0.d0
           hsc=0.d0
             g=0.d0

       dhex_dt=0.d0
       dhsc_dt=0.d0
         dg_dt=0.d0

        alph    =  alph1
     DO i=0,numpanels         ! trapezoidal integration

      CALL ksph_d( ka*alph, index, 0.d0 , dindex_dt              &
                  , etasc, etaex, assymetry                                 &
                  , detasc_dt , detaex_dt , dassymetry_dt                   &
                  , detasc_da0, detaex_da0, dassymetry_da0 )

       weight    =    DEXP(  (p+2)*DLOG(alph)    - DEXP(q*DLOG(meanradrel*alph)))


       if( i==0 .OR. i==numpanels ) then
          etasc = etasc/2.0d0
          etaex = etaex/2.0d0
          assymetry = assymetry/2.0d0

          detasc_dt = detasc_dt/2.0d0
          detaex_dt = detaex_dt/2.0d0
          dassymetry_dt = dassymetry_dt/2.0d0 
       end if

       hsc=hsc+          etasc*weight
       hex=hex+          etaex*weight
         g=  g+assymetry*etasc*weight

      dhsc_dt = dhsc_dt  + detasc_dt                    *weight
     
      dhex_dt = dhex_dt  + detaex_dt                    *weight 
   
        dg_dt = dg_dt    + ( dassymetry_dt * etasc                          &
                         +    assymetry    *detasc_dt  )*weight

         alph     =  alph     +  dalph
     END DO

     IF( hsc /= 0.d0 ) THEN
       g   =  g    /hsc
      dg_dt = dg_dt /hsc - g/hsc*dhsc_dt ! note that g here is the previous g divided by hsc
     ELSE
      IF( g /= 0.d0 ) THEN
        g    = 0.d0
      dg_dt = 0.d0
      END IF
     END IF

      c     =  dalph*DEXP((p+3)*DLOG(meanradrel))

    dhsc_dt =c*dhsc_dt
    dhex_dt =c*dhex_dt

   hsc = c * hsc  
   hex = c * hex

!---------------------------------------------------------------------------
    ELSE    ! p < 0.d0 .OR. q <= 0 ! IF( p >= 0.d0 .AND. q > 0.0d0 ) THEN     ! averaging over distribution
!---------------------------------------------------------------------------

      CALL ksph_d( ka      , index, 0.d0     , dindex_dt              &
                  , hsc, hex, g, dhsc_dt , dhex_dt , dg_dt                  &
                               , dhsc_da0, dhex_da0, dg_da0    )
 
!---------------------------------------------------------------------------
    END IF  ! IF( ka*(1.0d0+DBLE(m2)*sigmarel) > raythresh/10.0d0 ) THEN  ! absorption & scattering
!---------------------------------------------------------------------------
    dg(1) = dg_dt    !!from either part of IF, this is where to define dg
    dg(3) = 0.d0  !dg_da0
    dg(2) = 0.d0   

    c=pi*1.0d-3 * speclength**2
    w     =  c*k0*a0**3
   dw_dk0 =  w/k0

    dhsc(1)=w*dhsc_dt
    dhab(1)=w*dhex_dt - dhsc(1) 

    dhsc(2)= dw_dk0 *hsc 
    dhab(2)= dw_dk0 *hex - dhsc(2) 

    dhab(3)= 0.d0 !dw_da0*hex + w*dhex_da0 - dhsc_da0    
    dhsc(3)= 0.d0 !dw_da0*hsc+w*dhsc_da0
     hsc=w*hsc     
     hab=w*hex-hsc
      
!--------------------------------------------------------------------------
   ELSE   ! Rayleigh absorption without scattering
!--------------------------------------------------------------------------

    IF( a0 > 0.d0 ) THEN  

     IF( p >= 0.d0 .AND. q > 0.d0 ) THEN

      df_dk0= 8.0d0*pi*1.0d-9 *a0**4   ! so far p=0, q=1

     ELSE
      df_dk0= 4.0d0*pi/3.0d0*1.0d-9*a0**4 
      
     END IF
     f=df_dk0*k0 ! f = 8.0d0*pi*1.0d-9 * k0 *a0**4 ELSE f= 4.0d0*pi/3.0d0*1.0d-9 * k0 *a0**4

    ELSE                  ! a0 < 0
      f=0.d0 !all values of array set to zero
     df_dk0= 0.d0 ! ditto
    END IF                ! a0 > 0, a0 < 0

    epsil=epsil + 2.0d0 !!!!!why?

     hsc=0.d0
       g=0.d0
    dhsc=0.d0 ! arrays set to zero 
      dg=0.d0 ! ditto

     c = 0.1885d6*freq*speclength**3
     w    =  c*DIMAG(-1.0d0/epsil)
    dw_dt =  c*DIMAG(depsil_dt/epsil**2)

     hab=w*f

    dhab(1)= dw_dt*f !+ w * df_dt , but df_dt is zero here
    dhab(2)= w*df_dk0
    dhab(3)= 0.d0 ! w*df_da0

    END IF ! end of abs & scattering / Rayleigh absorption without scattering
   END IF ! DIMAG(epsil) > 0.d0 .AND. DREAL(epsil) /= 1.d0
 END IF                          ! end of the main part

! Check calculated quantities against constraints
    IF( hab < 0.d0 ) hab=0.d0
    IF( hsc < 0.d0 ) hsc=0.d0
    IF( g < -1.d0 ) g=-1.d0
    IF( g >  1.d0 ) g=1.d0

END SUBROUTINE hydrometeor_a0_is_constant_d


SUBROUTINE ksph_d( x , m , dx, dm,  ksc   ,  kex   ,  g                 &               
                                 , dksc_dt, dkex_dt, dg_dt              &
                                 , dksc_da, dkex_da, dg_da )

!-----------------------------------------------------------------------
!    SPHERICAL PARTICLE SCATTERING AND EXTINCTION EFFICIENCY 
!    A. GASIEWSKI 8/19/88, 7/15/92
!    RAYLEIGH LIMIT FROM WISCOMBE, APPL.OPT. V.19 p1505
!    MIE COEFFICIENT UPWARD RECURSION FROM DEIRMENDJIAN, 1969
!    "EM SCATTERING ON SPHERICAL POLYDISPERSIONS"
!    X = 2*PI*RADIUS/WAVELENGTH
!    M = COMPLEX INDEX OF REFRACTION
!    KSC= EFFICIENCY FACTOR FOR SCATTERING
!    KEX = EFFICIENCY FACTOR FOR EXTINCTION
!    G = SCATTERING ASYMMETRY FACTOR
!    IMAG(M)<0 FOR LOSSY MEDIA
!-----------------------------------------------------------------------
!
!   Rewritten for FORTRAN 90 by A. Voronovich May 25, 2003
!
!   Index of refraction _m depends on a real parameter _t 
!   (which could be, e.g., temperature T). Then:
!
!     dksc_dt=d(ksc)/dt=d(ksc)/dm*dm/dt, where dm/dt=_dm
!
!     and similarly for _kex, _g.
!
!   Similarly, _x depends on a real parameter _a, and
!
!     dksc_da=d(ksc)/d(_a)=d(ksc)/dx*dx/d(_a), where dx/d(_a)=_dx
!
!     and similarly for _kex, _g.
!
! Warning file with reference number _outfile should be opened outside
! of the subroutine; in the present version it is set _outfile=5.
!
!----------------------------------------------------------------------

IMPLICIT NONE

real(8), INTENT( IN) :: x, dx
DOUBLE COMPLEX  , INTENT( IN) :: m, dm

real(8), INTENT(OUT) :: ksc,kex,g, dksc_dt ,dkex_dt ,dg_dt    &
                                          , dksc_da ,dkex_da ,dg_da


INTEGER          :: n, nmax
real(8) :: x2, rem, f1, delksc, c, dx2, dc_dt, dc_da
DOUBLE COMPLEX   :: y, m2, m2m1, a1, a2, b1, d, wn, wnm1, wnm2          &
                   , an, smalan, smalbn, smlan0, smlbn0, dm2            &
                   , dy_dt, da1_dt, da2_dt, db1_dt, dd_dt               &
                   , dy_da, da1_da, da2_da, db1_da, dd_da               &
                   ,  dan_dt, dsmalan_dt, dsmlan0_dt                    &
                   ,  dan_da, dsmalan_da, dsmlan0_da                    &
                   ,          dsmalbn_dt, dsmlbn0_dt                    &
                   ,          dsmalbn_da, dsmlbn0_da                    &
                   ,  dwn, dwnm1, dwnm2, an_new, dan_new_dt, dan_new_da &
                   , ddelksc_dt , ddelksc_da             
 
DOUBLE COMPLEX  , PARAMETER :: iz=(0.d0,1.d0)

!---------------------------------

  y =m*x
  x2=x*x

  dy_dt= dm* x
  dy_da=  m*dx
  dx2=2*x*dx

!-----------------------------------------------------
  IF( CDABS(y) < 0.25d0 ) THEN   ! Rayleigh limit
!-----------------------------------------------------

   m2  =m*m
   m2m1=m2-1.0d0

   dm2=2*m*dm

   d=m2+2.0d0+(1.0d0-0.7d0*m2)*x2-(8.0d0*m2*m2-385.0d0*m2+350.0d0)/1400.0d0*x2*x2              &
       +2.0d0*iz*m2m1/3.0d0*x*x2*(1.0d0-x2/10.0d0)

   dd_dt=dm2-0.7d0*dm2*x2-(16.0d0*m2-385.0d0)*dm2/1400.0d0*x2*x2                   &
       +2.0d0*iz*dm2 /3.0d0*x*x2*(1.0d0-x2/10.0d0)

   dd_da=(1.0d0-0.7d0*m2)*dx2-(8.0d0*m2*m2-385.0d0*m2+350.0d0)/700.0d0*x2*dx2       &
       +2.0d0*iz*m2m1/3.0d0*((dx*x2+x*dx2)*(1.0d0-x2/10.0d0)-x*x2*dx2/10.0d0)

   a1=2.0d0*iz/3.0d0 *m2m1*   (1.0d0-x2/10.0d0+ (4.0d0*m2+5.0d0)/1400.0d0*x2*x2)/d
   b1=  iz/45.0d0*m2m1*x2*(1.0d0+(2.0d0*m2-5.0d0)/70.0d0*x2)/(1.0d0-(2.0d0*m2-5.0d0)/30.0d0*x2)
   a2=  iz/15.0d0*m2m1*x2*(1.0d0-x2/14.0d0)/(2.0d0*m2+3.0d0-(2.0d0*m2-7.0d0)/14.0d0*x2)

    c=CDABS(a1)**2+CDABS(b1)**2+1.66667d0*CDABS(a2)**2

   da1_dt=2.0d0*iz/3.0d0 * (1.0d0-x2/10.0d0+ (4.0d0*m2+5.0d0)/1400.0d0*x2*x2)/d*(dm2-m2m1*dd_dt/d)
   da1_da=2.0d0*iz/3.0d0 *m2m1/d*(-(1.0d0-x2/10.0d0+(4.0d0*m2+5.0d0)/1400.0d0*x2*x2)*dd_da/d       &
                          +(-dx2/10.0d0+(4.0d0*m2+5.0d0)/700.0d0*x2*dx2 )          )

   db1_dt=  iz/45.0d0*  x2*( dm2*(1.0d0+(2.0d0*m2-5.0d0)/70.0d0*x2)/(1.0d0-(2.0d0*m2-5.0d0)/30.0d0*x2)     &
                      +     m2m1*dm2*x2*2.0d0/(1.0d0-(2.0d0*m2-5.0d0)/30.0d0*x2)**2/21.0d0 )

   db1_da=  iz/45.0d0 * m2m1*( dx2*(1.0d0+(2.0d0*m2-5.0d0)/70.0d0*x2)/(1.0d0-(2.0d0*m2-5.0d0)/30.0d0*x2)     &
                        +x2*dx2*(2.0d0*m2-5.0d0) /(1.0d0-(2.0d0*m2-5.0d0)/30.0d0*x2)**2/21.0d0 )

   da2_dt=  iz/3.0d0 *x2*(1.0d0-x2/14.0d0)*(1.0d0+x2/14.0d0)/(2.0d0*m2+3.0d0-(2.0d0*m2-7.0d0)/14.0d0*x2)**2*dm2

   da2_da=  iz/15.0d0 *m2m1*(2.0d0*m2+3.0d0-x2/7.0d0*((2.0d0-x2/14.0d0)*m2+x2/4.0d0+3.0d0))*dx2         &
                                        /(2.0d0*m2+3.0d0-(2.0d0*m2-7.0d0)/14.0d0*x2)**2

    dc_dt=2.0d0*DREAL(a1*DCONJG(da1_dt)+          b1*DCONJG(db1_dt)        &
                                   +1.66667d0*a2*DCONJG(da2_dt) )

    dc_da=2*DREAL(a1*DCONJG(da1_da)+          b1*DCONJG(db1_da)        &
                                   +1.66667d0*a2*DCONJG(da2_da) )
   ksc=6.0d0*x2*x2*c
   kex=6.0d0*x*DREAL(a1+b1+1.66667d0*a2)
     g=DREAL(a1*DCONJG(b1+a2))/c

   dksc_dt=6.0d0*x2* x2*dc_dt
   dksc_da=6.0d0*x2*(x2*dc_da+2*dx2*c)

   dkex_dt=6.0d0* x*DREAL(da1_dt+db1_dt+1.66667d0*da2_dt)

   dkex_da= 6.0d0*dx*DREAL( a1   + b1   +1.66667d0*a2   )                  &
           +6.0d0* x*DREAL(da1_da+db1_da+1.66667d0*da2_da)

     dg_dt=-dc_dt/c*              DREAL( a1*DCONJG( b1   + a2   ))/c   &
           + DREAL(da1_dt*DCONJG(b1+a2) +a1*DCONJG(db1_dt+da2_dt))/c

     dg_da=-dc_da/c*              DREAL( a1*DCONJG( b1   + a2   ))/c   &
           + DREAL(da1_da*DCONJG(b1+a2) +a1*DCONJG(db1_da+da2_da))/c

!-----------------------------------------------------
  ELSE                           ! Mie case
!-----------------------------------------------------

   rem=DREAL(m)
   f1 =-8.0d0+(26.22d0+(-0.4474d0+(2.04d-3 -1.75d-4*rem)*rem**3)*rem)*rem**2

   IF( -DIMAG(m)*x > f1/2.0d0 ) THEN
    ! WRITE(*,'(/x,a)') ' KSPH_d : '
    ! WRITE(*,'(/x,a,4(x,e10.3))') ' Upward recursion unstable:    &
    !                x, Re(m), Im(m), f1= ', x, DREAL(m), DIMAG(m), f1
   END IF

!                    Mie recursion
!-------------------------------------------------
   n   =1
   nmax=INT2(x+4.05d0*x**0.333333d0+2)

   wnm1=DCMPLX(DSIN(x),DCOS(x))
   wnm2=-iz*wnm1
   wn  =wnm1/x-wnm2

  dwnm1=DCMPLX(DCOS(x),-DSIN(x))*dx
  dwnm2=-iz*dwnm1
  dwn  =(dwnm1-wnm1/x*dx)/x-dwnm2

   a1=EXP(-2*iz*y)
   an=-1.0d0/y+1.0d0/(1.0d0/y-iz*(1.0d0+a1)/(1.0d0-a1))

  da1_dt=-2*iz*a1*dy_dt
  da1_da=-2*iz*a1*dy_da

  dan_dt=dy_dt/y**2 + (dy_dt/y**2+2*iz*da1_dt/(1.0d0-a1)**2)                 &
                     /(1.0d0/y-iz*(1.0d0+a1)/(1.0d0-a1))**2

  dan_da=dy_da/y**2 + (dy_da/y**2+2*iz*da1_da/(1.0d0-a1)**2)                 &
                     /(1.0d0/y-iz*(1.0d0+a1)/(1.0d0-a1))**2

   d=an/m+1.0d0/x
   smalan=(d*DREAL(wn)-DREAL(wnm1))/(d*wn-wnm1)

   dd_dt=(dan_dt-an*dm/m)/m
   dd_da= dan_da/m-dx/x/x

   dsmalan_dt= (wn*DREAL(wnm1)-wnm1*DREAL(wn))*dd_dt/(d*wn-wnm1)**2
   dsmalan_da= (wn*DREAL(wnm1)-wnm1*DREAL(wn))*dd_da/(d*wn-wnm1)**2      &
              +( (d*DREAL(dwn)-DREAL(dwnm1))*(d* wn- wnm1)               &
                -(d*DREAL( wn)-DREAL( wnm1))*(d*dwn-dwnm1))              &
                                                    /(d*wn-wnm1)**2

   d=an*m+1.0d0/x
   smalbn=(d*DREAL(wn)-DREAL(wnm1))/(d*wn-wnm1)

   dd_dt= dan_dt*m+an*dm
   dd_da= dan_da*m-dx/x/x

   dsmalbn_dt= (wn*DREAL(wnm1)-wnm1*DREAL(wn))*dd_dt/(d*wn-wnm1)**2
   dsmalbn_da= (wn*DREAL(wnm1)-wnm1*DREAL(wn))*dd_da/(d*wn-wnm1)**2      &
              +( (d*DREAL(dwn)-DREAL(dwnm1))*(d* wn- wnm1)               &
                -(d*DREAL( wn)-DREAL( wnm1))*(d*dwn-dwnm1))              &
                                                    /(d*wn-wnm1)**2

   ksc=3.0d0*(CDABS(smalan)**2+CDABS(smalbn)**2)
   kex=3.0d0* DREAL(smalan+smalbn)
     g=0.d0

  dksc_dt=6.0d0*DREAL(smalan*DCONJG(dsmalan_dt)+smalbn*DCONJG(dsmalbn_dt))
  dksc_da=6.0d0*DREAL(smalan*DCONJG(dsmalan_da)+smalbn*DCONJG(dsmalbn_da))

  dkex_dt=3.0d0* DREAL(dsmalan_dt+dsmalbn_dt)
  dkex_da=3.0d0* DREAL(dsmalan_da+dsmalbn_da)

    dg_dt=0.d0
    dg_da=0.d0

    DO
      n=n+1
      wnm2=wnm1
      wnm1=wn
      smlan0=smalan
      smlbn0=smalbn

     dwnm2=dwnm1
     dwnm1=dwn
     dsmlan0_dt=dsmalan_dt
     dsmlan0_da=dsmalan_da
     dsmlbn0_dt=dsmalbn_dt
     dsmlbn0_da=dsmalbn_da


     wn=(2.0d0*n-1.0d0)/x*wnm1-wnm2
    dwn=(dwnm1-wnm1*dx/x)*(2.0d0*n-1.0d0)/x-dwnm2

      an_new   =-n/y+1.0d0/(n/y-an)
     dan_new_dt=n/y/y*dy_dt+(n/y/y*dy_dt+dan_dt)/(n/y-an)**2
     dan_new_da=n/y/y*dy_da+(n/y/y*dy_da+dan_da)/(n/y-an)**2

       an   = an_new
      dan_dt=dan_new_dt
      dan_da=dan_new_da

      d=an/m+n/x
     smalan=(d*DREAL(wn)-DREAL(wnm1))/(d*wn-wnm1)

     dd_dt=(dan_dt-an*dm/m)/m
     dd_da=dan_da/m-n*dx/x/x


   dsmalan_dt= (wn*DREAL(wnm1)-wnm1*DREAL(wn))*dd_dt/(d*wn-wnm1)**2
   dsmalan_da= (wn*DREAL(wnm1)-wnm1*DREAL(wn))*dd_da/(d*wn-wnm1)**2      &
              +( (d*DREAL(dwn)-DREAL(dwnm1))*(d* wn- wnm1)               &
                -(d*DREAL( wn)-DREAL( wnm1))*(d*dwn-dwnm1))              &
                                                    /(d*wn-wnm1)**2

      d=an*m+n/x
     smalbn=(d*DREAL(wn)-DREAL(wnm1))/(d*wn-wnm1)

     dd_dt=dan_dt*m+an*dm
     dd_da=dan_da*m-n*dx/x/x

   dsmalbn_dt= (wn*DREAL(wnm1)-wnm1*DREAL(wn))*dd_dt/(d*wn-wnm1)**2
   dsmalbn_da= (wn*DREAL(wnm1)-wnm1*DREAL(wn))*dd_da/(d*wn-wnm1)**2      &
              +( (d*DREAL(dwn)-DREAL(dwnm1))*(d* wn- wnm1)               &
                -(d*DREAL( wn)-DREAL( wnm1))*(d*dwn-dwnm1))              &
                                                    /(d*wn-wnm1)**2
     delksc=    dble(2*n+1)*(CDABS(smalan)**2+CDABS(smalbn)**2)

    ddelksc_dt=2.0d0*dble(2*n+1)*DREAL( smalan*DCONJG(dsmalan_dt)                &
                               +smalbn*DCONJG(dsmalbn_dt)  )
    ddelksc_da=2.0d0*dble(2*n+1)*DREAL( smalan*DCONJG(dsmalan_da)                &
                               +smalbn*DCONJG(dsmalbn_da)  )  

      ksc   = ksc   + delksc

     dksc_dt=dksc_dt+ddelksc_dt
     dksc_da=dksc_da+ddelksc_da

       kex   = kex   +dble(2*n+1)*DREAL( smalan   +smalbn   )

      dkex_dt=dkex_dt+dble(2*n+1)*DREAL(dsmalan_dt+dsmalbn_dt)
      dkex_da=dkex_da+dble(2*n+1)*DREAL(dsmalan_da+dsmalbn_da)

        g=g+(dble(n)*dble(n)-1.0d0)/dble(n)* DREAL( smlan0*DCONJG(smalan)               &
                                +smlbn0*DCONJG(smalbn))              &
           +(2.0d0*dble(n)-1.d0)/(dble(n)-1.0d0)/dble(n)*DREAL(smlan0*DCONJG(smlbn0))

       dg_dt=dg_dt+(dble(n)*dble(n)-1.d0)/dble(n)*DREAL(  smlan0   *DCONJG(dsmalan_dt)  &
                                      +dsmlan0_dt*DCONJG( smalan   )  &
                                      + smlbn0   *DCONJG(dsmalbn_dt)  &
                                      +dsmlbn0_dt*DCONJG( smalbn   )) &
         +(2.0d0*dble(n)-1.d0)/(dble(n)-1.d0)/dble(n)*DREAL(  smlan0   *DCONJG(dsmlbn0_dt)  &
                                      +dsmlan0_dt*DCONJG( smlbn0   ))

       dg_da=dg_da+(dble(n)*dble(n)-1.d0)/dble(n)*DREAL(  smlan0   *DCONJG(dsmalan_da)  &
                                      +dsmlan0_da*DCONJG( smalan   )  &
                                      + smlbn0   *DCONJG(dsmalbn_da)  &
                                      +dsmlbn0_da*DCONJG( smalbn   )) &
         +(2.0d0*dble(n)-1.d0)/(dble(n)-1.d0)/dble(n)*DREAL(  smlan0   *DCONJG(dsmlbn0_da)  &
                                      +dsmlan0_da*DCONJG( smlbn0   ))
      IF( n > nmax .OR. DABS(delksc) <= (1.0d-5*DABS(ksc)) ) EXIT
!!!!blw     IF( n > nmax .and. DABS(delksc) < (1.0d-5*DABS(ksc)) ) EXIT
    END DO

   ksc   =2.0d0/x2* ksc
  dksc_dt=2.0d0/x2*dksc_dt
  dksc_da=2.0d0/x2*dksc_da-ksc*dx2/x2


   kex   =2.0d0/x2* kex
  dkex_dt=2.0d0/x2*dkex_dt
  dkex_da=2.0d0/x2*dkex_da-kex*dx2/x2

     g=4.0d0/x2/ksc*g

    dg_dt=4.0d0/x2/ksc*dg_dt-g* dksc_dt/ksc
    dg_da=4.0d0/x2/ksc*dg_da-g*(dksc_da/ksc+dx2/x2)

!-----------------------------------------------------
  END IF                         ! Rayleigh/Mie cases
!-----------------------------------------------------

!--------------------------------------------------------------------
END SUBROUTINE ksph_d

SUBROUTINE REALPART(j,tt,nREAL,dnREAL_dt) ! written and tested by R. Hill Dec. 2003
IMPLICIT NONE
! the following two data statements contain polynomial coefficients for fits to the real part
! of the refractive index(in the data statement for index_data_transp in 
! SUBROUTINE h2o_ice_dielectric) for the first 9 wavelengths, next 6, and the last 20 wavelengths.
!R. Hill 12/5/03
real(8), INTENT(OUT) :: nREAL, dnREAL_dt
real(8), INTENT(IN)  :: tt ! temperature C
INTEGER ,         INTENT(IN)  :: j !index created in h2o_ice_dielectric to identify wavelength
INTEGER :: k 
real(8), DIMENSION(9)   :: first9
data first9 / &
   1.7948000000000d+000, &
   1.7921000000000d+000, &
   1.7884000000000d+000, &
   1.7860000000000d+000, &
   1.7843000000000d+000, &
   1.7832000000000d+000, &
   1.7825000000000d+000, &
   1.7820000000000d+000, &
   1.7816250000000d+000  /

real(8), DIMENSION(3,6)   :: next6
data next6 /           & ! 10 to 15
  -1.1249506638294d-007, 8.6699086830078d-006, 1.7816246293210d+000, &
  -3.6128402055896d-007, 8.7961603783943d-006, 1.7819277498148d+000, &
  -6.6391517773279d-007, 3.0003298458813d-005, 1.7830883015433d+000, &
  -1.1028402227018d-006, 4.2666779551484d-005, 1.7844273444445d+000, &
  -1.4720571837293d-006, 4.5725084706441d-005, 1.7853397148766d+000, &
  -7.3235697645008d-007, 1.1054287910918d-004, 1.7864637891976d+000  /

real(8), DIMENSION(2,20)   :: last20 ! 16 to 35
data last20 / &
 1.6217271474507d-004, 1.7870117133670d+000, &
 1.6481396417089d-004, 1.7871435002297d+000, &
 1.6651355075793d-004, 1.7872300413413d+000, &
 1.6745521359669d-004, 1.7872752870923d+000, &
 1.6745521359669d-004, 1.7872752870923d+000, &
 1.6832797427654d-004, 1.7872940514469d+000, &
 1.6825907211761d-004, 1.7872675700505d+000, &
 1.6913183279745d-004, 1.7872863344051d+000, &
 1.6913183279745d-004, 1.7872863344051d+000, &
 1.6913183279745d-004, 1.7872863344051d+000, &
 1.6913183279745d-004, 1.7872863344051d+000, &
 1.6913183279745d-004, 1.7872863344051d+000, &
 1.6913183279745d-004, 1.7872863344051d+000, &
 1.6913183279745d-004, 1.7872863344051d+000, &
 1.6913183279745d-004, 1.7872863344051d+000, &
 1.6913183279745d-004, 1.7872863344051d+000, &
 1.6913183279745d-004, 1.7872863344051d+000, &
 1.6913183279745d-004, 1.7872863344051d+000, &
 1.6913183279745d-004, 1.7872863344051d+000, &
 1.6913183279745d-004, 1.7872863344051d+000  /
! Note: All real part of refractive index data are equal below wl(22), so fit coefficients
! above are likewise equal. R. Hill 12/5/03

! j identifies the wavelength wl(j)
IF ( j < 10 ) THEN ! use first9 for j 
 nREAL      = first9(j)
dnREAL_dt   = 0.d0
ELSE IF (j<16) THEN ! use next6
k = j-9
 nREAL    = tt*( tt*next6(1,k)   + next6(2,k)   ) + next6(3,k)
! to test derivatives, it is useful to comment out next6(3,k), like this:!+ next6(3,k)
! otherwise,the finite difference derivatives become inaccurate. R. Hill 12/5/03 
dnREAL_dt = 2.d0*tt*next6(1,k)   + next6(2,k)
ELSE
k = j-15
 nREAL   =  tt*last20(1,k)   + last20(2,k)
dnREAL_dt   =  last20(1,k)
END IF

END SUBROUTINE REALPART


SUBROUTINE INVPOL(j, tt, nIMAG, dnIMAG_dt)! written and tested by R. Hill Dec. 2003
! calculate imaginary part of the refractive index for designated wavelengths
IMPLICIT NONE
INTEGER,          INTENT (IN) :: j ! designates wavelength 'wl' at which to evaluate curve fit
real(8), INTENT (IN) :: tt ! Celsius temperature
real(8):: nIMAG, dnIMAG_dt !  imaginary part of refractive index,
                                    ! and its derivative w.r.t. temperature
real(8) g ! internal variable for efficiency

! The following data statement contains coefficients of the fit
! to the imaginary part of the refractive index (in the data statement for
! index_data_transp in SUBROUTINE h2o_ice_dielectric).
!  Note from the code below that  tt  is the Celsius temperature.
! The first 2 coefficients (say, a,b) are for a term 1/(tt*a + b) and the last
! 3 coefficients (say, c,d,e) are for a term c*tt**2 + d*tt + e.
! Thus the fit is: [1/(tt*a + b)] +  [c*tt**2 + d*tt + e]. 
! The coefficients are to be used for wavelengths j = 1 through 14 and 25 through 28
! R. Hill 12/5/03
real(8), DIMENSION(5,35)   :: inv1pol2 
data inv1pol2 /                                  &
  -3.6646030650334d-001, 5.8627814653732d+001, -5.6289115820330d-007, -3.6437671234587d-005, -2.1850687789092d-004, &
  -5.3564166271146d-001, 5.9695753886372d+001, -1.0319722462800d-006, -6.7328211582158d-005, -4.1561749613748d-004, &
  -8.2322045951163d-001, 6.1704865204199d+001, -1.4191129833009d-006, -9.3502209112711d-005, -6.0179659075622d-004, &
  -1.1723024142989d+000, 6.3719389174892d+001, -2.0606619336682d-006, -1.3690062472861d-004, -9.1280069498642d-004, &
  -1.4473257183158d+000, 6.6753708467432d+001, -2.1358722751881d-006, -1.4261436670996d-004, -9.7590638335297d-004, &
  -1.7634344732996d+000, 6.9578518818719d+001, -2.1965218589651d-006, -1.4741712506459d-004, -1.0367300726803d-003, &
  -2.0540588761093d+000, 7.3165400934154d+001, -2.1177213896409d-006, -1.4266257303249d-004, -1.0242901197260d-003, &
  -2.3363121335135d+000, 7.7160560304587d+001, -1.9891761734163d-006, -1.3440970890144d-004, -9.8220769972350d-004, &
  -3.0261762777990d+000, 8.5521035325192d+001, -1.7480099700531d-006, -1.1889767796159d-004, -9.0190301559228d-004, &
  -3.9076692933771d+000, 9.4302974782422d+001, -1.5817743859706d-006, -1.0833524045946d-004, -8.5380783824600d-004, &
  -4.8353827989517d+000, 1.0516073297629d+002, -1.3015136443794d-006, -8.9637479750714d-005, -7.2885287370630d-004, &
  -8.1481383434313d+000, 1.5037890929986d+002, -4.4590690632704d-007, -3.1133588023336d-005, -2.7209008704293d-004, &
  -1.3680499628054d+001, 2.1649942168884d+002, -4.3662690185494d-009, -4.0220126349026d-007, -6.1982044288857d-006, &
  -2.2867936944882d+001, 3.1207897896411d+002, 1.8978626424265d-007, 1.3412210812321d-005, 1.2757248555890d-004, &
  -3.7946707866107d+001, 4.5153390966321d+002, 2.6425769702522d-007, 1.8902935491902d-005, 1.8967674083130d-004, &
  -5.3369424932138d+001, 5.7976041222038d+002, 2.7705355411137d-007, 1.9954012634026d-005, 2.0636564127610d-004, &
  -5.6607748706596d+001, 6.2576627621266d+002, 2.5257566611353d-007, 1.8174420393184d-005, 1.8721015034086d-004, &
  -6.1524348915724d+001, 7.1798150108606d+002, 2.0676962242144d-007, 1.4845107528360d-005, 1.5132879590947d-004, &
  -6.6190964495118d+001, 8.1601622873824d+002, 1.8131342498278d-007, 1.2989980063212d-005, 1.3114605693184d-004, &
  -7.0527890249583d+001, 9.1379366031428d+002, 1.3840616153874d-007, 9.8791214053461d-006, 9.8065106596828d-005, &
  -7.7524900128393d+001, 1.1735033726439d+003, 8.4445687564476d-008, 5.9804422310025d-006, 5.7163455847967d-005, &
  -8.3818293369649d+001, 1.4325413793441d+003, 2.6525199940735d-008, 1.8532266156849d-006, 1.6658273319943d-005, &
  -9.2265040122672d+001, 2.0949255109051d+003, 2.4066893258298d-009, 1.6315844938011d-007, 1.2742293727007d-006, &
  -9.9951829128780d+001, 2.8225590495227d+003, -8.9234216015160d-009, -6.1224314737768d-007, -4.9191379321349d-006, &
  -1.0504365518827d+002, 3.0798364263161d+003, -1.3264457891021d-008, -9.0719064417845d-007, -7.1483028938154d-006, &
  -1.1282512224493d+002, 3.2860020466099d+003, -9.2201693849082d-009, -6.3153168086588d-007, -5.0034822248015d-006, &
  -1.2140900905999d+002, 3.4270972883971d+003, -8.6972952888480d-009, -5.9627273373904d-007, -4.7748819475824d-006, &
  -1.3143348350302d+002, 3.5291405915935d+003, -7.9670524995063d-009, -5.4734524504591d-007, -4.4374428025952d-006, &
  -1.3850265759660d+002, 3.6491138049557d+003, -3.7007245291592d-010, -2.5594452374589d-008, -2.1461637574783d-007, &
  -1.5034753893115d+002, 3.6874838143447d+003, -1.0143796578698d-009, -6.7974551357539d-008, -4.8919549622996d-007, &
  -1.5996055442541d+002, 3.7334671274726d+003, -2.6939732226019d-009, -1.8237250159238d-007, -1.3823682023493d-006, &
  -1.8155681042950d+002, 3.7490822365838d+003, -2.6412081217632d-009, -1.7840934677766d-007, -1.3374178108536d-006, &
  -2.3459704621644d+002, 3.5808997154793d+003, -5.5089817846779d-009, -3.8078506224346d-007, -3.2312031777942d-006, &
  -2.8001390051706d+002, 3.2742395630101d+003, -1.7513788188842d-008, -1.2408415579980d-006, -1.1846364735637d-005, &
  -3.1392066762651d+002, 2.9491559498164d+003, -3.5847407796915d-008, -2.5683353071739d-006, -2.5772355068430d-005  /

 g = ( tt*inv1pol2(1,j)   + inv1pol2(2,j) )**(-1) 
nIMAG       =       g             +  tt*(tt*inv1pol2(3,j)  + inv1pol2(4,j) ) +inv1pol2(5,j)
dnIMAG_dt   = -inv1pol2(1,j)*g**2 + 2.d0*tt*inv1pol2(3,j)  + inv1pol2(4,j)
END SUBROUTINE INVPOL



SUBROUTINE INV(j, tt, nIMAG, dnIMAG_dt) ! written and tested by R. Hill Nov-Dec 2003
! calculate imaginary part of the refractive index for designated wavelengths
IMPLICIT NONE
INTEGER,          INTENT (IN) :: j ! designates wavelength 'wl' at which to evaluate curve fit
real(8), INTENT (IN) :: tt ! Celsius temperature
real(8), INTENT (OUT):: nIMAG, dnIMAG_dt ! imaginary part of refractive index,
                                                  ! and its derivative w.r.t. temperature
! the following data statement contains coefficients of the fit of the inverse of a cubic polynomial
! to the imaginary part of the refractive index (in the data statement for index_data_transp in SUBROUTINE h2o_ice_dielectric)
! they are to be used for wavelengths j = 15 through 24 and 29 through 35
real(8), DIMENSION(4,35)   :: inv3
data inv3 /                                  &
   2.6652569025449d-004, 2.4034948209791d-002, 1.3594739286085d-001, 6.0172239030401d+001, &
   8.0925301816504d-005, 1.1345083692264d-002, -1.2733941677063d-001, 6.1589791486567d+001, &
  -1.8185522550068d-004, -7.0338805002987d-003, -5.7420017961343d-001, 6.3948780877919d+001, &
  -5.0599123599133d-004, -2.7835328168715d-002, -1.0016663767510d+000, 6.7052873844535d+001, &
  -6.6383769740496d-004, -3.7802676435433d-002, -1.3045950448505d+000, 7.0674989837053d+001, &
  -8.4185914881998d-004, -4.9075460873528d-002, -1.6541644442663d+000, 7.4151644915034d+001, &
  -1.0716461778711d-003, -6.5369231915819d-002, -2.0981247947635d+000, 7.7966172790974d+001, &
  -1.1274601959110d-003, -6.8017169663409d-002, -2.3348592220893d+000, 8.2477793199243d+001, &
  -1.6355513674774d-003, -1.0495337810578d-001, -3.4211738266478d+000, 9.1021766641600d+001, &
  -2.2996261286780d-003, -1.5285518892590d-001, -4.7972074301033d+000, 1.0017515106770d+002, &
  -2.8893559553520d-003, -1.9716004937659d-001, -6.2756782458059d+000, 1.1074102235416d+002, &
  -5.1262094434226d-003, -3.8293424144816d-001, -1.2797540947652d+001, 1.4861086289756d+002, &
  -8.4037844366686d-003, -6.6767086421874d-001, -2.3823796320229d+001, 1.9955261998004d+002, &
  -1.2617007918334d-002, -1.0544928507280d+000, -4.1387948644352d+001, 2.6829590250710d+002, &
  -1.8331877526127d-002, -1.6145983335876d+000, -7.0066929859937d+001, 3.5887976394655d+002, &
  -2.2151994571666d-002, -2.0412144928653d+000, -9.7880097335250d+001, 4.3605385878007d+002, &
  -2.4440569684843d-002, -2.2382645773284d+000, -1.0484282548797d+002, 4.7208364219783d+002, &
  -2.2627817539964d-002, -2.1176556298154d+000, -1.0905088493726d+002, 5.5971080954168d+002, &
  -2.3795593639468d-002, -2.2475316815648d+000, -1.1747123435145d+002, 6.4232825931223d+002, &
  -2.7655125200579d-002, -2.5341398291676d+000, -1.2520104589463d+002, 7.3937440432658d+002, &
  -2.1480308873378d-002, -2.0362987190448d+000, -1.2428639047368d+002, 1.0140920643001d+003, &
  -3.7446905846765d-002, -3.1277666796579d+000, -1.3866274655010d+002, 1.3032764940870d+003, &
  -2.3318418147914d-002, -1.8916349683074d+000, -1.2287550206565d+002, 2.0341652282776d+003, &
  -2.2232202294734d-002, -1.6338421737541d+000, -1.1841433171892d+002, 2.8243737488408d+003, &
  -3.9328158681127d-002, -2.8987243330348d+000, -1.3824057810470d+002, 3.0800529016709d+003, &
  -5.0156889519076d-002, -3.8225360480873d+000, -1.6310370872250d+002, 3.2420292146538d+003, &
  -2.4177262289528d-002, -1.7235813999299d+000, -1.3812480953278d+002, 3.4478039852859d+003, &
  -1.7864365808577d-002, -1.2158186676432d+000, -1.4013722463817d+002, 3.5647644333674d+003, &
  -2.6204766467031d-003, -1.9926238664454d-001, -1.4110116431325d+002, 3.6469742654755d+003, &
   2.6146950028394d-002, 2.1288023974753d+000, -1.1553672107419d+002, 3.7583296157264d+003, &
   4.4619123566491d-002, 3.6663628768891d+000, -9.8432634616846d+001, 3.8661995898838d+003, &
   5.5206512028210d-002, 4.5116474894132d+000, -1.0698908904300d+002, 3.9046187270077d+003, &
   6.3931363815253d-002, 5.3025719389708d+000, -1.4332187954792d+002, 3.7884473538927d+003, &
   4.9176439163727d-002, 4.5622161482880d+000, -1.7925520847355d+002, 3.6041105396961d+003, &
   2.3392135372844d-002, 3.1201672780260d+000, -2.0595902976795d+002, 3.4405792315857d+003  /

 nIMAG   = ( tt*( tt*( tt*inv3(1,j)      +      inv3(2,j) )  + inv3(3,j)   ) + inv3(4,j))**(-1)
dnIMAG_dt=     -( tt*( tt*3.d0*inv3(1,j) + 2.d0*inv3(2,j) )  + inv3(3,j))*nIMAG**2
END SUBROUTINE INV

!======================================================================
SUBROUTINE h2o_mixed_dielectric(freq, temp, phase, epsil, depsil_dt, testvar1, testvar2, testvar3)
!======================================================================
! Complex dielectric constant for snow or graupel particle composed of ice, water and air
! Uses dielectric mixing theory [Matthew and Sadiku, 1985] 
! Rx=(epsilon_x-1)/(epsilon_x+F), epsilon=kp-j*kpp
! "subscripts" i=ice, w=water, s=snow, a=air, m=mixed_ice-water-air_particle 
!
! References
!  Sadiku, Matthew N O (1985), Refractive index of snow at microwave
!    frequencies, Applied Optics, 24(4), pp.572-575
!  Nishitsuji, A, and M Hirayama (1971), Anomalous attenuation of radio waves due to snowfall,
!    Electronics and Communications in Japan, Vol 54-B, No 11, pp. 27-33
!
! Modifications:
! 2/1/1994   GMSJ Developed routine 
! 5/28/2003  A. Voronovich converted original subroutine to FORTRAN 90
! 7/10/2003  R. Hill correctly coded analytic derivatives depsil_dt = d(epsil)/d(Temperature)
! 11/26/2003 R. Hill changed water fraction to function with continuous derivatives
! 12/5/2003  R. Hill changed form factor to function with continuous derivatives
! 2/6/2021   Kevin Schaefer cleaned up code and deleted unused code
! 2/6/2021   Kevin Schaefer changed name to h2o_mixed_dielectric
!----------------------------------------------------------------------
  implicit none

! inputs
  real(8), intent(in) :: freq  ! (GHz) frequency range of [1..1000] GHz
  real(8), intent(in) :: temp  ! (K) air temperature
  integer, intent(in) :: phase ! (-) phase

! outputs
  double complex, intent(out) :: epsil     ! (-) mixed phase dielectric constant
  double complex, intent(out) :: depsil_dt ! (1/K) derivative of mixed phase dielectric constant wrt temperature

! Internal variables
  real(8) w  ! (-) fraction liquid water in mixed phase, not total water mass
  real(8) f  ! (-) form factor
  real(8) ci ! (-) ratio of density of graupel/snow to density of ice
  real(8) cw ! (-) ratio of density of graupel/snow to density of water
  real(8) a
  real(8) arg
  real(8) dw_dt
  real(8) df_dx
  real(8) df_dt
  real(8) dy3n2_dx
  real(8) dpolyfit_dx
  real(8) dE_dx
  real(8) dx_dw
  real(8) polyfit
  real(8) E
  real(8) x
  real(8) y3n2
  real(8) xzeta 
  double complex epsilw     ! (-) dielectric water
  double complex epsili     ! (-) dielectric ice
  double complex depsilw_dt ! (1/K) derivative dielectric water wrt temperature
  double complex depsili_dt ! (1/K) derivative dielectric ice wrt temperature
  double complex rw         ! (-) refractive index water
  double complex ri         ! (-) refractive index ice
  double complex drw_dt     ! (1/K) derivative refractive index water wrt temperature
  double complex dri_dt     ! (1/K) derivative refractive index ice wrt temperature
  double complex rm         ! (-) refractive index mixed phase
  double complex drm_dt     ! (1/K) derivative refractive index mixed phase wrt temperature
  real(8) testvar1 ! test diagnostic variable 1
  real(8) testvar2 ! test diagnostic variable 2
  real(8) testvar3 ! test diagnostic variable 3
  real(8) wmax
  real(8) w_liq_trans

! parameters
  real(8), parameter :: b = 0.12959822356838d0
  real(8), parameter :: point = 16.d0
  real(8), parameter :: A3 = 1.4324101899475d0
  real(8), parameter :: zeta3 = -1.7495065207990d0
  real(8), parameter :: coefs(3) = (/ -2.0079700293637d-1, 6.8397684390001d-1, -4.2709198240675d-1 /)
     
! assign ci and cw according to phase type
  if ( phase == 4 ) then ! snow
    ci = 0.109d0 ! density of snow/density of ice 
    cw = 0.100d0 ! density of snow/density of water
  elseif ( phase == 5) then ! graupel  
    ci = 0.436d0 ! density of graupel/density of ice   
    cw = 0.400d0 ! density of graupel/density of water  
  endif

! test code to tansition to 100% water
  a = b*(320.d0 - 265.65d0)
  wmax = (1.d0 + dtanh(a))/(90.d0*b)
  w_liq_trans = (1.d0-wmax)/(1.d0 + dexp(2.d0*(273.15d0-temp)))

! fraction of liquid water content in the mixed phase
  a = b*(temp - 265.65d0)
  w = (1.d0 + dtanh(a))/(90.d0*b)
  dw_dt = dcosh(a)**(-2)/(90.d0)  

! determine form factor [Nishitsuji, A, and M Hirayama, 1971, fig 1]
  x = -DLOG10(w)
  IF (x < 3.999d0) THEN ! despite this IF  statement, this function and 
    ! its derivative are continuous at x = 4.d0
    arg = 1.d0/(x**2-point)
    E = DEXP(arg) ! E is a 'testing function', such that
    ! at x = sqrt(point) = 4, E and all of its derivatives are zero. 
    polyfit = x*(x*(x*coefs(1) + coefs(2)) + coefs(3)) + 1.d0
    xzeta = A3*x**zeta3
    y3n2 = DLOG10(2.d0) + xzeta*E*polyfit ! this is log10(form factor)
    f = 10.d0**y3n2
    dx_dw = -DLOG10(DEXP(1.d0))/w
    dE_dx = -2.d0*x*E*arg**2
    dpolyfit_dx = (x*(x*3.d0*coefs(1) + 2.d0*coefs(2)) + coefs(3))
    dy3n2_dx = xzeta*((zeta3/x)*E*polyfit + dE_dx*polyfit + E*dpolyfit_dx)
    df_dx = DLOG(10.d0)*f*dy3n2_dx
    df_dt = df_dx*dx_dw*dw_dt   
  ELSE
    f = 2.d0  ! i.e., 10.d0**DLOG10(2.d0)
    df_dt = 0.d0
  END IF
  if(f < 2.d0) then
    f = 2.d0
    df_dt = 0.d0
  endif

! dielectric of pure water and pure ice
  CALL h2o_liquid_dielectric(freq, 273.15d0 , epsilw, depsilw_dt) 
  CALL h2o_ice_dielectric(freq,DMIN1(temp,273.15d0), epsili, depsili_dt)

! complex indeces of refraction
  rw=(epsilw-1.0d0)/(epsilw+f)
  ri=(epsili-1.0d0)/(epsili+f)
  rm=(1.0d0-w)*ci*ri + w*cw*rw

! derivatives complex indeces of refraction
  drw_dt = -rw*df_dt/(epsilw+f)
  dri_dt =  (  depsili_dt - ri*(depsili_dt+df_dt)  )/(epsili+f)
  drm_dt = (1.d0-w)*ci*dri_dt + w*cw*drw_dt + dw_dt*(-ci*ri + cw*rw)

! dielectric constant
  epsil    = (1.0d0+f*rm)/(1.0d0-rm)
  depsil_dt = (df_dt*rm + drm_dt*(f + epsil))/(1.0d0-rm)
  testvar1=dreal(depsil_dt)

 END SUBROUTINE h2o_mixed_dielectric

!======================================================================
SUBROUTINE h2o_liquid_dielectric(freq, temp, epsil, depsil_dt)
!======================================================================
! Liquid water complex dielectric constant from Liebe [1983]
! Also see Wilheit and Chang [1979]
! DIMAG(epsil) assumed positive for lossy media.
!
! Modifications:
! 5/27/2003  A Voronovich converted routine to FORTRAN 90
! 2/6/2021   Kevin Schaefer cleaned up code
!----------------------------------------------------------------------
  implicit none

! inputs
  real(8), intent(in) :: temp ! (K) air temperature
  real(8), intent(in) :: freq ! (GHz) frequency

! outputs
  double complex, intent(out) :: epsil     ! (-) liquid water dielectric constant
  double complex, intent(out) :: depsil_dt ! (1/K) derivative of liquid water dielectric constant wrt temperature

! Internal Variables
  real(8) theta     ! (-) ratio of temperature to reference temperature
  real(8) dtheta    ! (1/K) derivative of theta wrt temperature
  real(8) ksmkinf
  real(8) dksmkinf  ! (?) derivative of ksmkinf
  real(8) ftt
  real(8) dftt      ! (?) derivative of ftt wrt temperature
  real(8) mag2pole
  real(8) dmag2pole ! (?) derivative of mag2pole wrt temperature
  real(8) kpp
  real(8) dkpp      ! (?) derivative of kpp wrt temperature

! ratio to referencetemperature
  theta = 300.0d0/temp
  dtheta = -theta/temp

! exponential scaling factor
  ftt=freq*4.17d-5*theta*dexp(7.13d0*theta)
  dftt=freq*4.17d-5*(1.0d0+7.13d0*theta)*dexp(7.13d0*theta)*dtheta

  ksmkinf =185.0d0-113.0d0/theta
  dksmkinf = (113.0d0/(theta*theta))*dtheta

  mag2pole = 1.0d0 + ftt*ftt
  dmag2pole = 2.0d0*ftt*dftt

  kpp = ksmkinf/mag2pole
  dkpp = dksmkinf/mag2pole-(kpp/mag2pole)*dmag2pole

  epsil = DCMPLX(kpp+4.9d0, kpp*ftt)
  depsil_dt = DCMPLX(dkpp,dkpp*ftt+kpp*dftt)

END SUBROUTINE h2o_liquid_dielectric

!======================================================================
SUBROUTINE h2o_ice_dielectric(freq, temp, epsil, depsil_dt)
!======================================================================
! Frozen water complex dielectric constant from tabulations by 
! valid temperature range: -80 to 0 deg C
! Valid frequency range: 1 to 1000 GHz
! DIMAG(epsil) assumed positive for lossy media
!
! Reference:
!  Warren, SG (1984), Applied Optics, Vol 23, No 8, pp 1206-1225
!
! Modifications:
!  5/28/2003   A Voronovich converted original subroutine into FORTRAN 90
!  11/17/2003  R. Hill extrapolated Warren [1984] -80 to -60 deg C
!  2/6/2021    Kevin Schaefer cleaned up code and deleted unused code
!  2/6/2021    Kevin Schaefer changed name to h2o_ice_dielectric
!----------------------------------------------------------------------
  implicit none

! inouts
  real(8), intent(in) :: freq ! (GHz) frequency
  real(8), intent(in) :: temp ! (K) air temperature

! outputs
  double complex, intent(out) :: epsil     ! (-) ice dielectric constant
  double complex, intent(out) :: depsil_dt ! (1/K) derivative of ice dielectric constant wrt temperature

! internal variables
  integer i       ! (-) wavelength table index 1
  integer j       ! (-) wavelength table index 2
  real(8) lambda  ! (mm) wavelength
  real(8) tt      ! (deg C) air temperature in centigrade
  real(8) nr 
  real(8) ni 
  real(8) dnr_dt 
  real(8) dni_dt
  real(8) ratio 
  real(8) prod 
  real(8) slope 
  real(8) intercept 
  real(8) dslope_dt 
  real(8) dintercept_dt 
  real(8) nREAL(35) 
  real(8) nIMAG(35) 
  real(8) dnREAL_dt(35) 
  real(8) dnIMAG_dt(35)
  double complex n
  double complex dn_dt
 
! Parameters
  real(8), PARAMETER :: c=2.99793d8  ! (m/s) speed of light
  real(8) wl(35) ! (mm) wavelength

! Assign wavelength table values
  data wl / &
   .29850d+00,.31620d+00,.35480d+00,.39810d+00,.44670d+00,.50120d+00,.56230d+00, &
   .63100d+00,.79430d+00,.10000d+01,.12590d+01,.25000d+01,.50000d+01,.10000d+02, &
   .20000d+02,.32000d+02,.35000d+02,.40000d+02,.45000d+02,.50000d+02,.60000d+02, &
   .70000d+02,.90000d+02,.11100d+03,.12000d+03,.13000d+03,.14000d+03,.15000d+03, &
   .16000d+03,.17000d+03,.18000d+03,.20000d+03,.25000d+03,.29000d+03,.32000d+03  /

! Convert temperature to centigrade
  tt=temp-273.15d0
  
! calculate wavelength (mm) from frequency
  lambda=c/freq *1.d-6

! Check if wavelength and temperature fall within table limits
  if(lambda >= wl(1) .and. lambda <= wl(35)) then ! ok
    if(tt <= 0.d0 .and. tt >= -80.d0) then ! ok

      ! wavelength table index
      j=2 
      do 
        if(lambda <= wl(j) .OR. j == 35) exit
        j=j+1
      enddo  

      ! Calculate imaginary part to refractive indeces at wl(j) and wl(j-1)
      ! Choose the best type of curve fit.
      if((j <= 14) .or. (j >= 25 .and. j <= 28) ) then
        call INVPOL(j, tt, nIMAG(j), dnIMAG_dt(j))
      else
        call INV(j, tt, nIMAG(j), dnIMAG_dt(j))
      endif

      i = j-1
      if( (i <= 14) .or. (i >= 25 .and. i <= 28)) then
        call INVPOL(i, tt, nIMAG(i), dnIMAG_dt(i))
      else
        call INV(i, tt, nIMAG(i), dnIMAG_dt(i))
      endif

      ! Linearly nterpolate imaginary part of refractive index in natural log space [Warren, 1984, table 2]
      ! ln(nIMAG) = a + b*ln(lambda)
      ratio = DLOG(wl(j)/wl(i))
      prod = DLOG(wl(j)*wl(i))
      slope = DLOG(nIMAG(j)/nIMAG(i))/ratio
      dslope_dt = (dnIMAG_dt(j)/nIMAG(j) - dnIMAG_dt(i)/nIMAG(i))/ratio
      intercept = (DLOG(nIMAG(j)*nIMAG(i)) - slope*prod)/2.d0
      dintercept_dt = (dnIMAG_dt(j)/nIMAG(j) + dnIMAG_dt(i)/nIMAG(i) - dslope_dt*prod )/2.d0
      ni = DEXP(intercept + slope*DLOG(lambda))
      dni_dt =  ni*(dintercept_dt + dslope_dt*DLOG(lambda))

      ! Calculate real part of refractive index
      call REALPART(j,tt,nREAL(j),dnREAL_dt(j))
      
      ! Check if need to inearly nterpolate real part of refractive index in natural log space [Warren, 1984, table 2]
      if(j < 23) then
        i=j-1
        call REALPART(i,tt,nREAL(i),dnREAL_dt(i))
        slope = (nREAL(j) - nREAL(i))/ratio
        dslope_dt = (dnREAL_dt(j) - dnREAL_dt(i))/ratio
        intercept = (nREAL(j) +  nREAL(i) - slope*prod)/2.d0
        dintercept_dt = (dnREAL_dt(j) + dnREAL_dt(i) - dslope_dt*prod )/2.d0
        nr =  intercept + slope*DLOG(lambda)
        dnr_dt = dintercept_dt + dslope_dt*DLOG(lambda)
      else
        nr = nREAL(j)
        dnr_dt = dnREAL_dt(j)
      endif

    else ! temperature out of range (less than -80,0)
      nr = 1.0d0
      ni = 0.0d0
      dnr_dt = 0.0d0
      dni_dt = 0.0d0
    endif ! temperature in range (-80,0)
  else !  wavelength out of range
    nr = 1.0d0
    ni = 0.0d0
    dnr_dt = 0.0d0
    dni_dt = 0.0d0
  endif ! wavelength inf range

! index of refraction
  n = DCMPLX(nr, ni)
  dn_dt = DCMPLX(dnr_dt,dni_dt)

! dielectric constant
  epsil    =  n*n
  depsil_dt = 2.0d0*n*dn_dt

END SUBROUTINE h2o_ice_dielectric
