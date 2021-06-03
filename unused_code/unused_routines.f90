
!=======================================================================
  subroutine construct_hydromet_profile (col_tot, frac_lay, prof_den)
!=======================================================================      
!  Constructs an atmospheric profile of single hydrometeor  
!
! Modifications:
!  2/28/2021 Kevin Schaefer created subroutine
!  4/1/2021  Kevin Schaefer added vertical redistribution
!----------------------------------------------------------------------
  use scan_Variables
  use dotlrt_variables
  use profiles
  use dotlrt_output

  implicit none

! input variables
  real(8) frac_lay(nlev) ! (-) fractions per layer of total colum hydromet mass
  real(8) col_tot        ! (g/m2) column total mass  

! output variables
  real(8) prof_den(nlev) ! (g/m3) hydrometeor density profile

! local variables
  integer ilev      ! (-) atm level index
  integer ibad      ! (-) bad value index
  integer nbad      ! (-) total number of bad layers
  real(8) mass_max(nlev) ! (g/m2) maximum mass per level
  real(8) mass_lay(nlev) ! (g/m2) mass per level
  real(8) mass_ext(nlev) ! (g/m2) extra mass per level
  real(8) x_mass    ! (g/m2) extra hydromet mass
  real(8) max_den   ! (g/m3) maximum possible hydromet density
  real(8) frac_tot  ! (g/m3) maximum possible hydromet density

! maximum allowed density
  max_den = 1d4

! calculate initial mass per layer
  do ilev=1,nlev
    mass_lay(ilev) = frac_lay(ilev)*col_tot
    mass_max(ilev) = max_den*atm(ilev)%hgt_del
  enddo
  stop

! redistribute mass if required
  do ibad = 1,10  ! Iterate a max of 10 times
    nbad=0
    x_mass = 0.d0
    frac_tot=0.d0
    do ilev=1,nlev
      mass_ext(ilev) = 0.
      if(mass_lay(ilev) > mass_max(ilev)) then
        mass_ext(ilev) = mass_lay(ilev) - mass_max(ilev)
        mass_lay(ilev) = mass_max(ilev)
        x_mass = x_mass + mass_ext(ilev)
        nbad = nbad+1
        frac_tot = frac_tot + frac_lay(ilev)
      endif
      !if (mass_lay(ilev) > 0.d0) print*, ilev, mass_lay(ilev), mass_max(ilev)
    enddo

    do ilev=1,nlev
      if(mass_lay(ilev) > mass_max(ilev)) then
        mass_lay(ilev) = mass_lay(ilev) + frac_lay(ilev)/frac_tot*x_mass
      endif
    enddo
    
    if (nbad == 0) exit
  enddo

! calculate density per layer
  do ilev=1,nlev
    prof_den(ilev) = mass_lay(ilev)/atm(ilev)%hgt_del
    !if (prof_den(ilev) > 0.d0) print*, ilev, prof_den(ilev), max_den
  enddo

  return                                                                    
  end


!-----------------------------------------------------------------------
!          This file contains subroutine HYDROMETEOR_a0_is_not_constant_d.F90
!          and all slave subroutines called by it
!
!     NOTE:  it is used here, that (p=0, q=0) or (p=0, q=1) only
!
!-----------------------------------------------------------------------
SUBROUTINE hydrometeor_a0_is_not_constant_d(freq,phase,temp,p,q,k0,a0, &
                               hab, hsc, g,   &
                               dhab,dhsc,dg,numpanels)

!From CALCPR0FILE the size distributions only let either k0 or a0 vary, but not both.
!For phases 1 (cloud liquid), 3 (ice),    k0 varies with water density, but a0 is fixed.
!For phases 2(rain), 4(snow), 5(graupel), a0 varies with water density, but k0 is fixed.
! Therefore there are no derivatives with respect to k0 in this subroutine
!  Differentiation was completed by Reg Hill on 8/12/03.   
!-----------------------------------------------------------------------
!   This is a version of the original subroutine which calculates along
!   with parameters _hab, _hsc, _g vectors _dhab(1-3), _dhsc(1-3), _dg(1-3)
!   which are derivatives of the appropriate values with respect to
!   parameters _temp, _k0, _a0, appropriately
!   [so that, e.g.   _dhsc(1)=d(_hsc)/d(_temp) 
!                  , _dhsc(2)=d(_hsc)/d(_k0  )
!                  , _dhsc(3)=d(_hsc)/d(_a0  ) ].
!
! Warning file with reference number _outfile should be opened outside
! of the subroutine; in the present version it is set _outfile=5.
!
!-----------------------------------------------------------------------

IMPLICIT NONE
                 
real(8),               INTENT( IN) :: freq, temp,p,q,k0,a0
real(8),               INTENT(OUT) ::  hab,  hsc,  g
real(8), DIMENSION(3), INTENT(OUT) :: dhab, dhsc, dg
INTEGER         ,               INTENT( IN) :: phase

INTEGER         , PARAMETER :: m1=20, m2=10, m3=8
real(8), PARAMETER :: raythresh=0.1d0

INTEGER :: i, numpanels


real(8) ::  pi, c , speclength, f, meanradrel, sigmarel, ka   &
                   , hex, etasc , etaex, assymetry                     &
                   , dka_da0, dhsc_dt , dhex_dt , dg_dt                &
                   ,          dhsc_da0, dhex_da0, dg_da0               &
                   , detasc_dt , detaex_dt , dassymetry_dt             &
                   , detasc_da0, detaex_da0, dassymetry_da0            &
                   , weight , alph, alph1, dalph, alphu, w, dw_dt 

real(8) :: ddalph_da0, dalph1_da0, dalph_da0, dweight_da0, dc_da0, dw_da0, df_da0

! real(8) :: dw_dk0, df_dk0

DOUBLE COMPLEX   ::   epsil, depsil_dt, index, dindex_dt
 
                   ! epsil=kp+i*kpp

EXTERNAL  h2o_liq_diel_d                        &
         ,h2o_ice_diel_d                        &
         ,h2omixeddiel_d                        &
         ,        ksph_d                        


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

    CALL h2o_liq_diel_d(freq,temp,epsil,depsil_dt)

    speclength=1.d0

   ELSE

    IF( phase == 3)  THEN                                                
     CALL h2o_ice_diel_d( freq, DMIN1(273.15d0,temp), epsil, depsil_dt) 
    ELSE    ! ice/liquid mixture
     CALL h2omixeddiel_d( freq, temp, phase         , epsil, depsil_dt)
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

     dka_da0=2.d0*pi*freq/3.0d2*meanradrel*speclength

      ka    =dka_da0*a0
 
! note  dka_dk0 = 0.d0 !!!!!!!!!!!!!!!!!!!!!!!  

!--------------------------------------------------------------------------
   IF( ka*(1.0d0+DBLE(m2)*sigmarel) > raythresh/10.0d0 ) THEN  ! absorption & scattering
!--------------------------------------------------------------------------

     index   =DCONJG(SQRT(epsil))

    dindex_dt=DCONJG(depsil_dt)/index/2.0d0
! note dindex_k0 or _a0 = 0.d0 !!!!!!!!!!!!!!!!!!!!1

   IF( p >= 0.d0 .AND. q > 0.0d0 ) THEN     ! averaging over distribution

     dalph      = 1.d0/dble(m1)
 ddalph_da0 = 0.d0 !!!!!!!!!!!!!!!!!!!!!

     IF( ka > 1.d0 ) THEN
      IF( ka > dble(m3) ) THEN
       dalph     = dalph/dble(m3)
  ddalph_da0 = 0.d0   !!!!!!!!!!!!!!!!!!!!
      ELSE
       dalph     =  dalph/ka
  ddalph_da0 = -dalph/ka*dka_da0 !!!!!!!other derivatives are zero!!!!!!!!!!!!
      END IF
     END IF

      alph1     = 1.0d0-m2*sigmarel
     dalph1_da0 = 0.d0
     IF( alph1 <= 0.d0 ) THEN
  alph1    =  dalph
     dalph1_da0 = ddalph_da0  !!!!!!!!!!!!!other derivatives are zero!!!!!!!!!!!!!!!
  END IF
        alphu=1.0d0+m2*sigmarel
     numpanels=INT2((alphu-alph1)/dalph+1.5d0) !?? differentiate an integer?? NO
        dalph     = (alphu -  alph1    )/numpanels
   ddalph_da0 = (      - dalph1_da0)/numpanels !!!!!!!!!!other derivatives are zero!!!!!!!!!!!1


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
       dalph_da0 = dalph1_da0 !!!!!!!!!!!!!!other derivatives are zero!!!!!!!!!!!!!!!!!!!!!!!!
! test a derivative point


     DO i=0,numpanels         ! trapezoidal integration

      CALL ksph_d( ka*alph, index, dka_da0*alph , dindex_dt              &
                  , etasc, etaex, assymetry                                 &
                  , detasc_dt , detaex_dt , dassymetry_dt                   &
                  , detasc_da0, detaex_da0, dassymetry_da0 )

       weight    =    DEXP(  (p+2)*DLOG(alph)    - DEXP(q*DLOG(meanradrel*alph)))
      dweight_da0 = weight*(  (p+2)*dalph_da0/alph - &       !!!!!other derivatives are zero
       DEXP(q*DLOG(meanradrel*alph)) *q*dalph_da0/alph  )  !!!!!other derivatives are zero!!!!!!!!!!!!!!!!!!
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

       hsc=hsc+          etasc*weight
       hex=hex+          etaex*weight
         g=  g+assymetry*etasc*weight

      dhsc_dt = dhsc_dt  + detasc_dt                    *weight
    
      dhex_dt = dhex_dt  + detaex_dt                    *weight 
   
        dg_dt = dg_dt    + ( dassymetry_dt * etasc                          &
                         +    assymetry    *detasc_dt  )*weight
 
      dhsc_da0= dhsc_da0 + detasc_da0                   *weight &
                         +  etasc                       *dweight_da0 !!!!!!!!!!!!!!!!!!!

      dhex_da0= dhex_da0 + detaex_da0                   *weight &
                    +  etaex                       *dweight_da0 !!!!!!!!!!!!!!!!! 
        dg_da0= dg_da0   + ( dassymetry_da0 * etasc                         &
                         +    assymetry     *detasc_da0)*weight &
                         +  (assymetry*etasc)           *dweight_da0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!

        alph     =  alph     +  dalph
        dalph_da0 = dalph_da0 + ddalph_da0 !!!!!!!!!!!!!!!!!!!!!!!!!!!! 
     END DO

     IF( hsc /= 0.d0 ) THEN
       g   =  g    /hsc
      dg_dt = dg_dt /hsc - g/hsc*dhsc_dt   ! note that g here is the previous g divided by hsc
      dg_da0= dg_da0/hsc - g/hsc*dhsc_da0  ! ditto
     ELSE
      IF( g /= 0.d0 ) THEN
       ! WRITE(outfile,'(30x,a )') ' HYDROMETEOR_EXT_5PH_d: '
       ! WRITE(outfile,'( 5x,a/)') ' INCONSISTENT SCATTERING AND ASSYMETRY '
       g    = 0.d0
      dg_dt = 0.d0
      dg_da0= 0.d0      
      END IF
     END IF

      c     =  dalph*DEXP((p+3)*DLOG(meanradrel))
     dc_da0 = ddalph_da0*c/dalph   !!!!!!!!!!!!!!!!!!!!!!!!!
 !   hsc = c * hsc  !moved to below derivatives
 !   hex = c * hex  !

    dhsc_dt =c*dhsc_dt
    dhex_dt =c*dhex_dt
          
    dhsc_da0 =  c    * dhsc_da0 &
             + dc_da0 *  hsc !!!the old value of hsc must be used!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

    dhex_da0 =  c     *dhex_da0 &
             + dc_da0 * hex !!!!!!!!the old value of hex must be used!!!!!!!!!!!!!!!!!!!!!!!! 
   hsc = c * hsc  
   hex = c * hex  
!---------------------------------------------------------------------------
    ELSE    ! p < 0.d0 .OR. q <= 0 ! IF( p >= 0.d0 .AND. q > 0.0d0 ) THEN     ! averaging over distribution
!---------------------------------------------------------------------------

      CALL ksph_d( ka      , index, dka_da0       , dindex_dt              &
                  , hsc, hex, g, dhsc_dt , dhex_dt , dg_dt                  &
                               , dhsc_da0, dhex_da0, dg_da0    )

!   dg(1) = dg_dt    !!!!!!!!!!!!!!!! why put it into array now?
!   dg(3) = dg_da0
!dg(2) = 0.d0   !!!!!!!!!!!!!!!! why put it into array now?

!---------------------------------------------------------------------------
    END IF  ! IF( ka*(1.0d0+DBLE(m2)*sigmarel) > raythresh/10.0d0 ) THEN  ! absorption & scattering
!---------------------------------------------------------------------------
    dg(1) = dg_dt    !!from either part of IF, this is where to define dg
    dg(3) = dg_da0
dg(2) = 0.d0  ! no k0 derivative 

    c=pi*1.0d-3 * speclength**2
    w     =  c*k0*a0**3
!  dw_dk0 =  w/k0      !!!!!!!!!!!!!!!!!!!!!!!!!!!!
   dw_da0 =  w/a0*3.d0 !!!!!!!!!!!!!!!!!!!!!!!!!!! also defined below 

! calculate derivatives first, since _hsc=w*hsc, _hab=w*hex-hsc will be re-defined

 !   hsc=w*hsc        ! moved to below derivates for
 !  hab=w*hex-hsc     ! numerical efficiency and clarity

    dhsc(1)=w*dhsc_dt
    dhab(1)=w*dhex_dt - dhsc(1) ! why introduce arrays here?


    dhab(2)= 0.d0 ! dw_dk0 *hex !-dhsc_dk0 I think that this second terms is zero
    dhsc(2)= 0.d0 ! dw_dk0 *hsc 

    dhsc(3)= dw_da0*hsc+w*dhsc_da0
    dhab(3)= dw_da0*hex + w*dhex_da0 - dhsc(3)   

     hsc=w*hsc     
     hab=w*hex-hsc       
!--------------------------------------------------------------------------
   ELSE   ! Rayleigh absorption without scattering
!--------------------------------------------------------------------------

    IF( a0 > 0.d0 ) THEN  

     IF( p >= 0.d0 .AND. q > 0.d0 ) THEN

      f= 8.0d0*pi*1.0d-9 *a0**4 *k0   ! so far p=0, q=1

     ELSE
      f= 4.0d0*pi/3.0d0*1.0d-9*a0**4 *k0   
      
     END IF

     df_da0=4.d0*f/a0 
 !   df_dt = 0.d0

    ELSE                  ! a0 < 0
      f=0.d0 !all values of array set to zero
!     df_dk0= 0.d0 ! ditto
     df_da0= 0.d0 ! ditto
 !   df_dt = 0.d0
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
    dhab(2)= 0.d0 !w*df_dk0
    dhab(3)= w*df_da0

!----------------------------------------------------------------------------
    END IF ! end of abs & scattering / Rayleigh absorption without scattering
!----------------------------------------------------------------------------
   END IF ! DIMAG(epsil) > 0.d0 .AND. DREAL(epsil) /= 1.d0
!-----------------------------------------------------------------------
!------------------------------------------------------
 END IF                          ! end of the main part
!------------------------------------------------------

! Check calculated quantities against constraints

    IF( hab < 0.d0 )THEN
     ! WRITE(outfile,'(50x,a )') ' HYDROMETEOR_EXT_5PH_d: ' 
     ! WRITE(outfile,'(x,a   )') ' WARNING : Absorption calculation negative; &
     !                            using zero instead'
     hab=0.d0
    END IF

    IF( hsc < 0.d0 )THEN
     ! WRITE(outfile,'(50x,a )') ' HYDROMETEOR_EXT_5PH_d: ' 
     ! WRITE(outfile,'(x,a   )') ' WARNING : Scattering calculation negative; &
     !                            using zero instead'
     hsc=0.d0
    END IF

    IF( g < -1.d0 )THEN
     ! WRITE(outfile,'(50x,a )') ' HYDROMETEOR_EXT_5PH_d: ' 
     ! WRITE(outfile,'(x,a   )') ' WARNING : Asymmetry factor _g less than -1; &
     !                            using _g=-1 instead'
     g=-1.d0
    END IF

    IF( g >  1.d0 )THEN
     ! WRITE(outfile,'(50x,a )') ' HYDROMETEOR_EXT_5PH_d: '
     ! WRITE(outfile,'(x,a   )') ' WARNING : Asymmetry factor _g exceeds    1; &
     !                           using _g= 1 instead'
     g=1.d0
    END IF

!-----------------------------------------------------------------------
END SUBROUTINE hydrometeor_a0_is_not_constant_d
!
!======================================================================
subroutine prognostic_fpar(sib,generic)
!======================================================================
! calculates fpar based on prognostic LAI
!
! Modifications:
!   Kevin Schaefer created routine (6/3/13)
!----------------------------------------------------------------------

use kinds
use sibtype
use sib_const_module

implicit none    

type(sib_t), dimension(subcount) :: sib ! Sib variable tree
real(kind=real_kind) :: generic        ! Area average canopy absorbed fraction  of PAR

! LOCAL VARIABLES
integer i  ! index
! For LTran and LRef:
!   (1,1) : shortwave, green plants
!   (2,1) : longwave, green plants
!   (1,2) : shortwave, brown plants
!   (2,2) : longwave, brown plants
real(kind=dbl_kind) :: fpar_slp    ! slope of fPAR
real(kind=dbl_kind) :: scatp       ! Canopy transmittance + reflectance wrt PAR
real(kind=dbl_kind) :: park        ! Mean canopy absorption optical depth wrt PAR
real(kind=dbl_kind) :: leaf_grw_max ! (-) scaling factor to keep leaf_pool near max observed value
real(kind=dbl_kind) :: leaf_inhib ! (-) scaling factor to keep leaf_pool near max observed value
real(kind=dbl_kind) :: leaf_exp ! (-) scaling factor to keep leaf_pool near max observed value
real(kind=dbl_kind) :: templ ! (-) scaling factor to keep leaf_pool near max observed value
!
! loop through points
    do i = 1, subcount
!
! canopy extinction coeficient
      scatp = sib(i)%param%green * (sib(i)%param%tran(1,1) + sib(i)%param%ref(1,1)) + &
        (1. - sib(i)%param%green) * (sib(i)%param%tran(1,2) + sib(i)%param%ref(1,2))
      park = sqrt(1. - scatp) * sib(i)%param%gmudmu
!
! fpar obtained by integrating Beer's Law from canopy top to bottom
      sib(i)%casa%fpar = sib(i)%param%vcover * (1. - exp(-park * sib(i)%casa%lai / sib(i)%param%vcover))
      sib(i)%casa%fpar = amax1(fparmin, sib(i)%casa%fpar)
      sib(i)%casa%fpar = amin1(fparmax, sib(i)%casa%fpar)
!
! fpar slope scaling factor: fparslope divided by max fparslope
      sib(i)%casa%fpar_slp = park*exp(-park * sib(i)%casa%lai / sib(i)%param%vcover)
!
! exponential growth limits
      fpar_slp=1.
      leaf_grw_max=exp(-fpar_slp*sib(i)%casa%fpar/fparmax)
!
! growth inhibition function based on fpar
      fpar_slp=10.
      templ = 1. + exp(fpar_slp*(sib(i)%casa%fpar-generic))
      leaf_inhib = 1.0/templ
!
! growth inhibition function based on lai
      fpar_slp=5.
      templ=exp(fpar_slp*(sib(i)%casa%lai-sib(i)%param%lai_max))
      templ=max(templ,1.e-7)
      leaf_exp=1.0/(1.+templ)

      sib(i)%diag%testvar1=leaf_grw_max
      sib(i)%diag%testvar2=leaf_inhib
      sib(i)%diag%testvar3=leaf_exp

    enddo
!
end subroutine prognostic_fpar
!
!=======================================================================
      subroutine soil_organic_layer(num_bio)
!=========================================================================
! soil thermodynamic properties
!
! Modifications:
!  Kevin Schaefer created routine (5/28/13)
!--------------------------------------------------------------------------
!
    use Mapper_Variables
    use sibscan_Variables
!
    IMPLICIT NONE
!
! define internal variables
integer num_bio
integer j
real(kind=dbl_kind)carb_lay(nsoil)     ! (kg/m2) Carbon per layer
real(kind=dbl_kind)carb_lay_max(nsoil) ! (kg/m2) maximum Carbon per layer
real(kind=dbl_kind) :: totalroot       ! (?) total root density in soil column
real(kind=dbl_kind) :: carbon          ! (mol C m-2) total carbon in soil column
real(kind=dbl_kind) :: Del_carb        ! (mol C m-2) carbon moved to next lower layer in soil column
real(kind=dbl_kind) :: olt             ! (m) organic layer thickness
!
! pure constituent densities
    rho_min=2700.
    rho_om_max=140*0.5*1000./12.
    sib(1)%param%rootd    = biotab(num_bio)%rootd
!
! root fraction per soil layer
    sib(1)%param%kroot=kroot(num_bio)
    totalroot = (1.0 - exp(-sib(1)%param%kroot*sib(1)%param%rootd))/sib(1)%param%kroot
    ztop = 0.0
    do j=1,nsoil
      zbot = ztop + sib(1)%prog%dz(j)
      if(sib(1)%param%rootd>=zbot) then
        sib(1)%param%rootf(j) = (exp(-sib(1)%param%kroot*ztop) &
          - exp(-sib(1)%param%kroot*zbot))/(sib(1)%param%kroot * totalroot)
      elseif(sib(1)%param%rootd>ztop.and.sib(1)%param%rootd<=zbot) then
        sib(1)%param%rootf(j) = (exp(-sib(1)%param%kroot*ztop) &
          - exp(-sib(1)%param%kroot*sib(1)%param%rootd))/(sib(1)%param%kroot * totalroot)
      else
        sib(1)%param%rootf(j) = 0.
      endif
      ztop = zbot
    enddo
!
! total carbon in soil column
    carbon=m_om*1000./12.
!
! theoretical maximum Organic layer thickness (OLT)
    olt=carbon/rho_om_max
    test10(1)=olt
!
! initial vertical distribution
    carb_lay(:)=0.
    do j=1,nsoil
      carb_lay(j)=sib(1)%param%rootf(j)*carbon
      carb_lay_max(j)=rho_om_max*sib(1)%prog%dz(j)
    enddo
!
! initial Organic layer thickness (OLT)
    olt=0.
    do j=1,nsoil
      f_om=carb_lay(j)/carb_lay_max(j)
      if(f_om>1.) olt=olt+sib(1)%prog%dz(j)
    enddo
    test10(2)=olt
!
! test code to vertically redistribute organic material
    do j=1,nsoil-1
      if(carb_lay(j)>carb_lay_max(j)) then
        del_carb=carb_lay(j)-carb_lay_max(j)
        carb_lay(j)=carb_lay(j)-del_carb
        carb_lay(j+1)=carb_lay(j+1)+del_carb
      endif
    enddo

!    test10(1)=carb_lay(1)/carb_lay_max(1)
!    test10(2)=carb_lay(2)/carb_lay_max(2)
!    test10(3)=carb_lay(3)/carb_lay_max(3)
!    test10(4)=carb_lay(4)/carb_lay_max(4)
!    test10(5)=carb_lay(5)/carb_lay_max(5)
!    test10(6)=carb_lay(6)/carb_lay_max(6)
!    test10(7)=carb_lay(7)/carb_lay_max(7)
!    test10(8)=carb_lay(8)/carb_lay_max(8)
!    test10(9)=carb_lay(9)/carb_lay_max(9)
!    test10(10)=carb_lay(10)/carb_lay_max(10)
!
! simulated OLT
    olt=0.
    do j=1,nsoil
      f_om=carb_lay(j)/carb_lay_max(j)
      if(f_om==1.) olt=olt+sib(1)%prog%dz(j)
    enddo
    test10(3)=olt
!
! soil porosity
    poros_min=0.489-0.00126*Sand
    poros_om=.9
    poros=(1.-f_om)*poros_min+f_om*poros_om

    return
    end
!
!=======================================================================
      subroutine soil_thermo_properties(ix,iy)
!=========================================================================
! soil thermodynamic properties
!
! Modifications:
!  Kevin Schaefer split off from main program (5/27/13)
!--------------------------------------------------------------------------
!
    use Mapper_Variables
    use sibscan_Variables
!
    IMPLICIT NONE
!
integer ix,iy
real junkvar, var1, var2
!
! define internal variables
!
! pure constituent densities
    rho_min=2700.
    rho_om_max=140.
    if(ix==1.and.iy==1) then
      deform=0.
      var1=0.
      var2=0.
      sum=0.
    endif
!
! other constants
    Lfus=0.3336e6
    deltaT=0.5
    lay_dz=0.02
!
! Constituent thermal properties
    c_min=(2.128*sand+2.385*clay)/(sand+clay)*1.e6 ! J/m3/K
    c_liq=4188.    ! J/kg/K
    c_ice=2117.27  ! J/kg/K
    c_om=2.5e6     ! J/m3/K
    tcon_min=(8.80*sand+2.92*clay)/(sand+clay)
    tcon_dry_om=0.05

! root fraction
    kroot(1)=3.9
    kroot(2)=3.9
    kroot(3)=2.0
    kroot(4)=5.5
    kroot(5)=5.5
    kroot(6)=2.0
    kroot(7)=5.5
    kroot(8)=2.0
    kroot(9)=2.0
    kroot(10)=5.5
    kroot(11)=2.0
    kroot(12)=5.5
!
! traditional root fraction from SiB
! equation valid as long as generic is depth such that dx = dz
    f_root=(exp(-kRoot(10)*generic)-exp(-kRoot(10)*(generic+dx)))/(1.-exp(-kRoot(10)*1.))
!
! organic soil fraction
    f_om=f_root*m_om/rho_om_max/dx
    if(f_om>1.) f_om=1.

! simpler organic soil fraction for remote sensing ALT
    f_om=kRoot(10)*m_om/(1.-exp(-kRoot(10)*1.))*exp(-kRoot(10)*generic)/rho_om_max
    if(f_om>1.) f_om=1.
    if(generic>1.) f_om=0.
    sum=sum+f_om*dx*rho_om_max
!
! soil porosity
    poros_min=0.489-0.00126*Sand
    poros_om=.9
    poros=(1.-f_om)*poros_min+f_om*poros_om

    test10(1)=f_om
    test10(2)=deform
    test10(3)=var2
    test10(4)=deform
    test10(5)=tcon_ice**((1.-v_water)*poros)
    test10(6)=v_ice+v_liq+v_min+v_om+v_air
    test=deform
!
! surface deformation
    if(generic>0.) then
      var1=var1+(rho_water-rho_ice)/rho_ice*dx*100.            ! pure water
      deform=deform+poros*(rho_water-rho_ice)/rho_ice*dx*100.  ! organic mix
      var2=var2+poros_min*(rho_water-rho_ice)/rho_ice*dx*100.  ! pure mineral soil
    endif
!
! water content
    satfrac=(m_liq/rho_water+m_ice/rho_ice)/poros/lay_dz
    if(satfrac>1.) satfrac=1.
              Theta=WWW(1)*poros_min
!
! volume fractions
    v_water=0.
    if(satfrac>0.) v_water=m_liq/rho_water/(m_liq/rho_water+m_ice/rho_ice)
!    satfrac=www(1)
    v_ice=(1.-v_water)*satfrac*poros
    v_liq=v_water*satfrac*poros
    v_air=(1.-satfrac)*poros
    v_om=f_om*(1.-poros_om)
    v_min=(1.-f_om)*(1.-poros_min)
    v_sol=1.-poros
!
! thermal conductivity
    rho_dry_bulk=rho_min*(1.-poros)
    tcon_dry_min=(0.135*rho_dry_bulk+64.7)/(2700-0.947*rho_dry_bulk)
    tcon_dry=(1.-f_om)*tcon_dry_min+f_om*tcon_dry_om
    tcon_sat=tcon_min**v_min*tcon_om**v_om*tcon_liq**(v_water*poros)*tcon_ice**((1.-v_water)*poros)
             ! if(v_water==0.) then
             !   kersten=satfrac
             ! else
    if(satfrac<=0.1) then
      kersten=0.
    else
      kersten=1.+log10(satfrac)
    endif
    thermcon=kersten*tcon_sat+(1.-kersten)*tcon_dry
!
! heat capacity
    heat=(1.-poros)*(1.-f_om)*c_min+(1.-poros)*f_om*c_om+m_liq*c_liq/lay_dz+m_ice*c_ice/lay_dz
    junkvar=v_min*c_min+v_om*c_om+m_liq*c_liq/lay_dz+m_ice*c_ice/lay_dz
!
!hydraulic conductivity
    satco_min=0.0070556*10.**(-.884+.0153*sand)/1000.
    satco_om=1.e-4
    satco=(1.-f_om)*satco_min+f_om*satco_om
    frzscale=10.**(-(1.25*(satco*1000.-3.)**(2.)*v_ice))
!
! Matric Potential
    phsat_min=-10.*10**(1.88-0.0131*sand)/1000.
    phsat_om=-10.3/1000.
    phsat=(1.-f_om)*phsat_min+f_om*phsat_om
    bee_min = 2.91+0.159*clay
    bee_om=2.7
    bee=(1.-f_om)*bee_min+f_om*bee_om

    satfrac=(m_liq/rho_water+m_ice/rho_ice)/poros/lay_dz
    if(satfrac>1.) then
      heat=0.
      junkvar=0.
      thermcon=0.
      tcon_sat=0.
    endif
    if(f_om>1.) then
      heat=0.
      junkvar=0.
      thermcon=0.
      satco=0.
      phsat=0.
      bee=0.
    endif
    if(sand+clay>100.) then
      tcon_min=0.
      tcon_dry_min=0.
      tcon_dry=0.
      thermcon=0.
      tcon_sat=0.
    endif
!
! liquid water fraction for simple soil model
    if(td(1)>(tf)) then ! liquid water
      fw=1.
      fw_slope=0.
    else if (td(1)<(tf-1.)) then ! ice
      fw=0.
      fw_slope=0.
    else ! mix between ice and water
      fw=td(1)-tf+1.
      fw_slope=1.
                
      fw=1.581976707*exp(-(td(1)-tf)**2.)-0.581976707
      fw_slope=1.581976707*exp(-(td(1)-tf)**2.)*(-2.)*(td(1)-tf)

      fw=1./(1.+EXP(20.*(tf-td(1)-.5)))
      fw=-2.*(td(1)-tf)**3.-3.*(td(1)-tf)**2.+1.
      fw_slope=-6.*(td(1)-tf)**2.-6.*(td(1)-tf)
    end if
!
! simple soil model heat capacity
    c_liq=4.188e6  ! J/m3/K
    c_ice=2.094e6  ! J/m3/K
    csoilthaw=fw*theta*c_liq
    csoilfroz=(1.-fw)*theta*c_ice
    heat=c_min+csoilfroz+csoilthaw+Lfus*rho_water*Theta*fw_slope
! 
! simple soil model thermal conductivity
    tcon_thaw=0.15+(tcon_min*tcon_liq**Theta-0.15)*satfrac
    tcon_froz=0.15+(tcon_min*tcon_ice**Theta-0.15)*satfrac
    ThermCon=(1.-fw)*tcon_froz+fw*tcon_thaw
!
! damping depth
    damp=sqrt(tau*365.*24.*3600./pi*ThermCon/heat)

    return
    end
!
!=======================================================================
      subroutine Respiration_scale_factors(biome_num)
!=========================================================================
! respiration temperature and moisture scaling factors
!
! Modifications:
!  Kevin Schaefer split off from main program (1/10/13)
!--------------------------------------------------------------------------
!
    use Mapper_Variables
    use sibscan_Variables
!
    IMPLICIT NONE
!
    integer biome_num
!
! define internal variables
    real var1, var2
    real temp_frz   ! (-) respiration freezing scaling factor
    real root_die_frost       ! (K) local temperature
    real casa_q10   ! (-) original casa Q10 scaling factor
    real casa_moist ! (-) original casa moisture scaling factor
!
! Q10 scaling factor
    soilQ10=exp(log(Q10)/10.*(td(1)-298.15))
!
! CASA temperature respiration scaling factor
    casa_q10=exp(log(2.)/10.*(td(1)-303.15))
!
! frozen soil respiration inhinbition function
! sfrzi and hfrzi defined in sib_const_module.f90
    sfrzi=generic
    sfrzi=1.
    hfrzi=271.
    temp_frz=1./(1.+exp(sfrzi*(hfrzi-td(1))))
!
! frozen soil exp decrease scaling factor based on Mikan et al. [2002]
    sib(1)%param%sfrzi=2.
    sib(1)%param%hfrzi=273.15
    sib(1)%param%tcfrzi=-5.*log(10.)/sib(1)%param%sfrzi+sib(1)%param%hfrzi
    if(td(1)>sib(1)%param%tcfrzi) then
      temp_frz=exp(sib(1)%param%sfrzi*(td(1)-sib(1)%param%hfrzi))
      temp_frz=min(temp_frz,1.)
    else
      temp_frz=0.
    endif
    test10(1)=temp_frz
!
! Freeze factor formulation exactly based on Mikan et al [2002] Incubation curve fits
    sib(1)%param%sfrzi=.5
    sib(1)%param%hfrzi=273.15
    sib(1)%param%tcfrzi=-5.*log(10.)/sib(1)%param%sfrzi+sib(1)%param%hfrzi
    if(td(1)>sib(1)%param%tcfrzi) then
      temp_frz=exp(sib(1)%param%sfrzi*(td(1)-sib(1)%param%hfrzi))
      temp_frz=min(temp_frz,1.)
    else
      temp_frz=0.
    endif
    test10(2)=temp_frz
!
! Romanovsky Osterkamp [2000] liquid water fraction
    sib(1)%param%hfrzi=273.15
    temp_frz=0.12*(sib(1)%param%hfrzi-td(1))**(-0.5)
    temp_frz=min(temp_frz,1.)
    test10(3)=temp_frz
!
! Monson et al. [2006]
    sib(1)%param%hfrzi=273.15
    temp_frz=exp(log(105./10.)*(td(1)-sib(1)%param%hfrzi))
    temp_frz=min(temp_frz,1.)
    test10(4)=temp_frz

    test=temp_frz

!
! liquid water fraction
!
! combined temperature freeze scaling factor
    var1=temp_frz*soilQ10
    if(var1<1.e-4) var1=0.
    var1=max(var1,1.e-4)

!    test10(3)=soilQ10
!    test10(4)=var1
!    test=temp_frz*soilQ10
!
! respiration soil moisture scaling function
    var1=clay/100.
    Wsat=0.25*var1+0.5
    zm=-2.*var1**3-0.4491*var1**2+.2101*var1+.3478
    Wopt=-0.08*var1*var1+0.22*var1+0.59
    MoistExp=((WWW(1)**zm-Wopt**zm)/(1-Wopt**zm))**2
    MoistExp=min(MoistExp,10.)
    moist(1)=0.8*Wsat**MoistExp+0.2
!
! leaf and fine root mortality increases under frozen soil conditions
    sib(1)%param%sftf=1.
    sib(1)%param%hftf=273.15
    sib(1)%param%tcftf=log(4.9e6)/sib(1)%param%sftf+sib(1)%param%hftf
    if(sib(1)%prog%td(1)<sib(1)%param%tcftf) then
      root_die_frost=1. + EXP(sib(1)%param%sftf*(sib(1)%prog%td(1)-sib(1)%param%hftf))
      root_die_frost=1. + 49./root_die_frost
    else
      root_die_frost=1.
    endif
!    test10(1)=root_die_frost
!
! not sure, might be wilting point
    var2=((1-Wopt**zm)*sqrt(log(0.375)/log(wsat))+ Wopt**zm)**(-zm)
    var2=(phsat/0.02/BioTab(biome_num)%phi_half)**(1./bee)
!
! CASA soil moisture scaling factor
    casa_moist=1./(1.+exp((.52-WWW(1))*35.))
!
!
    return
    end
!
!=======================================================================
      subroutine phosib_temp_factor(biome_num)
!=========================================================================
! photosynthesis temperature scaling factor
!
! Modifications:
!  Kevin Schaefer split off from main program (8/10/09)
!--------------------------------------------------------------------------
!
      use Mapper_Variables
      use sibscan_Variables
!
      IMPLICIT NONE
!
    integer biome_num
!
! define internal variables
    real templ   ! low temperature factor
    real temph   ! high temperature factor
    real tempf   ! freezing temperature factor
    real rstfac3 ! overall temperature scaling factor  
!
! original temperature response
    sib(1)%param%slti(1) =BioTab(biome_num)%SLTI
    sib(1)%param%hltii(1)=BioTab(biome_num)%HLTI
    sib(1)%param%shti(1) =BioTab(biome_num)%SHTI
    sib(1)%param%hhti(1) =BioTab(biome_num)%hhTI
    temph = 1. + EXP(sib(1)%param%shti(1) * (sib(1)%prog%tc - sib(1)%param%hhti(1) ))
    templ = 1. + EXP(sib(1)%param%slti(1) * (sib(1)%param%hltii(1) - sib(1)%prog%tc))
    rstfac3=1./templ/temph
    test10(1)=rstfac3
!
! 10C linear shutdown
    hfti = 283.15
    if(sib(1)%prog%tc < sib(1)%prog%tcmin) sib(1)%prog%tcmin = sib(1)%prog%tc
    sib(1)%prog%tcmin = sib(1)%prog%tcmin + (2.0_dbl_kind *generic)
    sib(1)%diag%frost_stress=(sib(1)%prog%tcmin-tice)/(10.)
    sib(1)%diag%frost_stress=min(sib(1)%diag%frost_stress,1._dbl_kind)
    sib(1)%diag%frost_stress=max(sib(1)%diag%frost_stress,0._dbl_kind)
    test10(2)=rstfac3*sib(1)%diag%frost_stress
    test=rstfac3*sib(1)%diag%frost_stress
!
! 5C linear shutdown
    hfti = 278.15
    if(sib(1)%prog%tc < sib(1)%prog%tcmin) sib(1)%prog%tcmin = sib(1)%prog%tc
    sib(1)%prog%tcmin = sib(1)%prog%tcmin + (2.0_dbl_kind *generic)
    sib(1)%diag%frost_stress=(sib(1)%prog%tcmin-tice)/(5.)
    sib(1)%diag%frost_stress=min(sib(1)%diag%frost_stress,1._dbl_kind)
    sib(1)%diag%frost_stress=max(sib(1)%diag%frost_stress,0._dbl_kind)
    test10(3)=rstfac3*sib(1)%diag%frost_stress
!
! Frost inhibition function
    sfti = 0.6
    hfti = 269.15
    if(sib(1)%prog%tc < sib(1)%prog%tcmin) sib(1)%prog%tcmin = sib(1)%prog%tc
    sib(1)%prog%tcmin = sib(1)%prog%tcmin + (2.0_dbl_kind *generic)
    tempf = 1. + EXP(sfti * (hfti - sib(1)%prog%tcmin))
    sib(1)%diag%frost_stress = 1.0/tempf
    test10(4)=rstfac3*sib(1)%diag%frost_stress
!
! simple frost on/off switch
    if(sib(1)%prog%tc < sib(1)%prog%tcmin) sib(1)%prog%tcmin = sib(1)%prog%tc
    sib(1)%prog%tcmin = sib(1)%prog%tcmin + (2.0_dbl_kind *generic)
    if(sib(1)%prog%tcmin<273.15) then
      sib(1)%diag%frost_stress=0.
    else
      sib(1)%diag%frost_stress=1.
    endif
    test10(5)=rstfac3*sib(1)%diag%frost_stress
!
! high temperature response
    sib(1)%param%slti(1) =BioTab(biome_num)%SLTI
    sib(1)%param%hltii(1)=BioTab(biome_num)%HLTI
    sib(1)%param%shti(1) =BioTab(biome_num)%SHTI
    sib(1)%param%hhti(1) =BioTab(biome_num)%hhTI+10.
    temph = 1. + EXP(sib(1)%param%shti(1) * (sib(1)%prog%tc - sib(1)%param%hhti(1) ))
    templ = 1. + EXP(sib(1)%param%slti(1) * (sib(1)%param%hltii(1) - sib(1)%prog%tc))
    rstfac3=1./templ/temph
!    test10(2)=rstfac3
!
      return
      end
!
!=======================================================================
      subroutine leaf_temp_factor(biome_num)
!=========================================================================
! leaf growth temperature scaling factor
!
! Modifications:
!  Kevin Schaefer split off from main program (8/10/09)
!--------------------------------------------------------------------------
!
      use Mapper_Variables
      use sibscan_Variables
!
      IMPLICIT NONE
!
    integer biome_num
!
! define internal variables
    real temp    ! temp variable
    real templ   ! low temperature factor
    real temph   ! high temperature factor
    real tempf   ! freezing temperature factor
    real rstfac3 ! overall temperature scaling factor 
    real lma(12)        ! leaf mass per area
    real leaf_grw_max    ! (-) scaling factor to restrict leaf pool to < LAImax
    real leaf_grw_fpar   ! (-) scaling factor to keep leaf_pool consistent with fPAR
    real leaf_grw_test   ! (-) scaling factor to test
    real leaf_grw_frz   ! (-) scaling factor to slow leaf_pool at low temp
    real leaf_grw_tot   ! (-) total leaf pool growth scaling factor
    real root_grw_tot   ! (-) root pool growth scaling factor
    real leaf_die_frz  ! (-) leaf pool freezing scaling factor
    real wood_grw_nsc  ! (-) wood pool nsc scaling factor
    real lma_scale    ! (-) leaf biomass scaling factor
    real max_leaf     ! (mol/m2) maximum leaf pool size
    real root_shoot   ! (-) root to shoot ratio
    real can_q10      ! (-) canopy temperature scaling factor
    real del_fpar  ! (-) fpar change rate

    lma(1)=40.0
    lma(2)=33.0
    lma(3)=33.0
    lma(4)=100.0
    lma(5)=33.0
    lma(6)=25.0
    lma(7)=25.0
    lma(8)=25.0
    lma(9)=28.0
    lma(10)=50.0
    lma(11)=50.0
    lma(12)=16.7
!
! fpar leaf growth scaling factor
    leaf_grw_fpar=(fpar-fPARmin/(fPARmax-fPARmin))
    leaf_grw_fpar=max(0.0,leaf_grw_fpar)
    leaf_grw_fpar=min(1.0,leaf_grw_fpar)
    leaf_grw_fpar=exp(-4.0*fpar/fPARmax)
!
! maximum leaf pool size
    lma_scale=0.5
    max_leaf=lma_scale*MorphTab(biome_num)%LAImax*lma(biome_num)/mwc
!
! max LAI leaf growth scaling factor
    leaf_grw_max=1.-generic
    leaf_grw_max=max(0.0,leaf_grw_max)
    leaf_grw_max=min(1.0,leaf_grw_max)
    leaf_grw_max=1.
!
! growth rate (generic is delta(fpar)
    del_fpar=generic/0.9
    leaf_grw_max=1./(1.+exp(-50.*del_fpar))
!
! scaling factor to keep near obs lai
! generic is Leaf_pool/storemin
    leaf_grw_test=4./(1.+exp(10.*(generic-.5)))
!
! leaf growth frost/freezing scaling factor
    leaf_grw_frz=1./(1.+EXP(1.3*(273.-tc)))
    leaf_grw_frz=max(0.1,leaf_grw_frz)
!
! total leaf growth scaling factor
    leaf_grw_tot=leaf_grw_fpar*leaf_grw_max*leaf_grw_frz

    test10(1)=leaf_grw_tot
    test10(2)=exp(-4.0*fpar/fPARmax)
    test10(3)=exp(-5.0*fpar/fPARmax)
    test=leaf_grw_tot
!
! root growth scaling factor
    root_shoot=1.
    root_grw_tot=root_shoot*leaf_grw_tot
!
! LMA scaling factor
    lma_scale=0.5*(1.+(fPARmax-fpar)/(fPARmax-fPARmin))
!
! assign constants
    sfti = 0.6
    hfti = 269.15
    sib(1)%param%slti(1) =BioTab(biome_num)%SLTI
    sib(1)%param%hltii(1)=BioTab(biome_num)%HLTI
    sib(1)%param%shti(1) =BioTab(biome_num)%SHTI
    sib(1)%param%hhti(1) =BioTab(biome_num)%hhTI
    sib(1)%param%trop(1) =BioTab(biome_num)%trop
    sib(1)%param%trda(1) =BioTab(biome_num)%trda
    sib(1)%param%trdm(1) =BioTab(biome_num)%trdm
!
! leaf decay temperature scaling factor
    can_q10 = 2.0**(0.1*(tc-sib(1)%param%trop(1)))
    if(sib(1)%prog%tc >= sib(1)%param%trdm(1) )then
      can_q10 = can_q10/(1.+EXP(sib(1)%param%trda(1)*(tc-sib(1)%param%trdm(1))))
    endif
!
! photosynthesis freezing inhibition function
    tempf = 1. + EXP(sfti * (hfti - tc))
    sib(1)%diag%frost_stress = 1.0/tempf
!
! photosynthesis high/low temperature inhibition functions
    temph = 1. + EXP(sib(1)%param%shti(1) * (tc - sib(1)%param%hhti(1) ))
    templ = 1. + EXP(sib(1)%param%slti(1) * (sib(1)%param%hltii(1) - tc))
!
! total photosynthesis temperature scaling factor
    rstfac3 = 1./(temph*templ*tempf)
!
! leaf mortality freezing function
    sfti = generic
    hfti = 273.15
    temp = 1. + EXP(sfti * (tc-hfti))
    leaf_die_frz=1.+49./temp
!
! wood growth nsc scaling factor
    wood_grw_nsc=1./(1.+EXP(10.*(1.-generic)))
!
      return
      end
!
!=======================================================================
      subroutine LAI_from_leaf_pool()
!=========================================================================
! conversion factor from leaf pool size to LAI
!
! Modifications:
!  Kevin Schaefer split off from main program (1/19/10)
!--------------------------------------------------------------------------
!
      use Mapper_Variables
      use sibscan_Variables
!
      IMPLICIT NONE
!
! define internal variables
    integer i,j,n  ! indeces
    real lma(12)    ! (g/m2) leaf mass per area
    real sla(12)    ! (m2/g) specific leaf area (leaf area per mass)
    integer indx(7) ! indeces of biomes
    real lai_max    ! maximum LAI
    real mass_max   ! leaf mass for laimax
    real mass_one   ! leaf mass for lai=1
    real mass_half  ! leaf mass half point
    real mass       ! leaf mass
    real dmass   ! mass increment
    real fpar_min ! sla slope
    real fpar_max   ! sla y-intercept
    real sla_slope ! sla slope
    real sla_int   ! sla y-intercept
    real sla_loc   ! sla value
    integer nsteps ! number integration steps
    integer istart  ! biome start
    integer istop   ! biome stop
!
! set lma values
    lma(1)=40.0
    lma(2)=33.0
    lma(3)=33.0
    lma(4)=100.0
    lma(5)=33.0
    lma(6)=25.0
    lma(7)=25.0
    lma(8)=25.0
    lma(9)=28.0
    lma(10)=50.0
    lma(11)=50.0
    lma(12)=16.7
!
! calculate sla values
    do j=1,12
      sla(j)=mwc/lma(j)
!      print*, j,MorphTab(j)%LAImax
    enddo
!
! pick indeces from lowest to highest lma
    indx(1)=12
    indx(2)=6
    indx(3)=9
    indx(4)=5
    indx(5)=1
    indx(6)=10
    indx(7)=4

    istart=1
    istop=7
    do i=istart,istop
      j=indx(i)
!
! constant sla approach
      lai=sla(j)*generic
!      test10(i)=lai
!      test10(1)=lai
!
! 2-segment integral approach
      lai=0.
      dmass=0.01
      nsteps=generic/dmass
      mass_max=MorphTab(j)%LAImax/sla(j)*2./3.
      mass=0.
      do n=1,nsteps
        mass=mass+dmass
        if(mass<=mass_max) then
          sla_slope=sla(j)/mass_max
          sla_int=sla(j)
          sla_loc=sla_slope*mass+sla_int
        else
          sla_loc=2.*sla(j)
        endif
        lai=lai+dmass*sla_loc
      enddo
!      test10(i)=lai
!      test10(2)=lai
!
! 2-segment analytical approach
      lai=0.
      dmass=0.01
      mass_max=MorphTab(j)%LAImax/sla(j)*2./3.
      if(generic<=mass_max) then
        sla_slope=sla(j)/mass_max
        sla_int=sla(j)
        lai=sla_slope/2.*generic*generic+sla_int*generic
      else
        lai=MorphTab(j)%LAImax+2*sla(j)*(generic-mass_max)
      endif
!      test10(i)=lai
!      test10(3)=lai
!
! 3-segment integral approach
      lai=0.
      dmass=0.01
      nsteps=generic/dmass
      mass_one=1./sla(j)
      mass_max=(MorphTab(j)%LAImax/sla(j)+mass_one/4.)*2./3.
      mass=0.
      do n=1,nsteps
        mass=mass+dmass
        if(mass<=mass_one) then
          sla_loc=sla(j)
        elseif(mass>mass_one.and.mass<=mass_max) then
          sla_slope=sla(j)/(mass_max-mass_one)
          sla_int=sla(j)*(mass_max-2.*mass_one)/(mass_max-mass_one)
          sla_loc=sla_slope*mass+sla_int
        else
          sla_loc=2.*sla(j)
        endif
        lai=lai+dmass*sla_loc
      enddo
!      test10(i)=lai
!      test10(4)=lai
!
! 3-segment analytical approach
      mass_one=1./sla(j)
      mass_max=(MorphTab(j)%LAImax/sla(j)+mass_one/4.)*2./3.
      sla_slope=sla(j)/(mass_max-mass_one)
      sla_int=sla(j)*(mass_max-2.*mass_one)/(mass_max-mass_one)
      if(generic<=mass_one) then
        lai=sla(j)*generic
        sla_loc=sla(j)
      elseif(generic>mass_one.and.generic<=mass_max) then
        lai=1.+sla_slope/2.*(generic*generic-mass_one*mass_one)+sla_int*(generic-mass_one)
        sla_loc=sla_slope*generic+sla_int
      else
        lai=1.+sla_slope/2.*(mass_max*mass_max-mass_one*mass_one)+sla_int*(mass_max-mass_one)
        lai=lai+2.*sla(j)*(generic-mass_max)
        sla_loc=2*sla(j)
      endif
      test10(i)=lai
!      test10(5)=lai
!
! 3-segment analytical approach with constant lai_max
      lai_max=8.
      mass_one=1./sla(j)
      mass_max=(lai_max/sla(j)+mass_one/4.)*2./3.
      sla_slope=sla(j)/(mass_max-mass_one)
      sla_int=sla(j)*(mass_max-2.*mass_one)/(mass_max-mass_one)
      if(generic<=mass_one) then
        lai=sla(j)*generic
        sla_loc=sla(j)
      elseif(generic>mass_one.and.generic<=mass_max) then
        lai=1.+sla_slope/2.*(generic*generic-mass_one*mass_one)+sla_int*(generic-mass_one)
        sla_loc=sla_slope*generic+sla_int
      else
        lai=1.+sla_slope/2.*(mass_max*mass_max-mass_one*mass_one)+sla_int*(mass_max-mass_one)
        lai=lai+2.*sla(j)*(generic-mass_max)
        sla_loc=2*sla(j)
      endif
!      test10(i)=lai
!      test10(6)=lai
!
! exponential analytical approach
      mass_half=6.
      lai=sla(j)*(2.*generic+mass_half*exp(-generic/mass_half)-mass_half)
!      test10(i)=lai
!      test10(6)=lai

      mass_half=1.
      sla_loc=sla(j)*(2.-exp(-generic/mass_half))
!      test10(1)=sla_loc
      mass_half=2.
      sla_loc=sla(j)*(2.-exp(-generic/mass_half))
!      test10(2)=sla_loc
      mass_half=3.
      sla_loc=sla(j)*(2.-exp(-generic/mass_half))
!      test10(3)=sla_loc
      mass_half=4.
      sla_loc=sla(j)*(2.-exp(-generic/mass_half))
!      test10(4)=sla_loc
      mass_half=5.
      sla_loc=sla(j)*(2.-exp(-generic/mass_half))
!      test10(5)=sla_loc
      mass_half=6.
      sla_loc=sla(j)*(2.-exp(-generic/mass_half))
!      test10(6)=sla_loc
!
! original exponential approach in SiBCASA
      mass_max=0.7213*lma(j)/mwc*MorphTab(j)%LAImax
      if(generic<=mass_max) then
        lai=-0.6931*generic*mwc/lma(j)/MorphTab(j)%LAImax
        lai=MorphTab(j)%LAImax/(-0.6931)*alog(lai+1)
      else
        lai=MorphTab(j)%LAImax+4.*(generic-mass_max)*mwc/lma(j)
      endif
!      test10(i)=lai
!      test10(7)=lai
!
! calculate fpar with this lai
      fPAR_min=0.01
      fPAR_max=0.95
      call aparnew (LAI, &
                    Green, &
                    BioTab(j)%LTran, &
                    BioTab(j)%LRef, &
                    gmudmu, &
                    fVCover, &
                    fPAR, &
                    fPAR_max, &
                    fPAR_min)
!      test10(i)=fpar

    enddo
!    test=lai
!
      return
      end
!
!
!=======================================================================
      subroutine sla_scaling()
!=========================================================================
! conversion factor from leaf pool size to LAI
!
! Modifications:
!  Kevin Schaefer split off from main program (1/19/10)
!--------------------------------------------------------------------------
!
      use Mapper_Variables
      use sibscan_Variables
!
      IMPLICIT NONE
!
! define internal variables
    integer i,j,n  ! indeces
    real lma(12)    ! (g/m2) leaf mass per area
    real sla(12)    ! (m2/g) specific leaf area (leaf area per mass)
    integer indx(7) ! indeces of biomes
    real leaf_min   ! (mole) min leaf pool size
    real leaf_max   ! (mole) max leaf pool size
    real lai_slope  ! slope of LAI vs. leaf pool curve
    real mass_max   ! maximum leafmass
    real dmass   ! mass increment
    real sla_slope ! sla slope
    real sla_int   ! sla y-intercept
    real sla_loc   ! sla value
    integer nsteps ! number integration steps
    integer istart  ! biome start
    integer istop   ! biome stop
!
! set lma values
    lma(1)=40.0
    lma(2)=33.0
    lma(3)=33.0
    lma(4)=100.0
    lma(5)=33.0
    lma(6)=25.0
    lma(7)=25.0
    lma(8)=25.0
    lma(9)=28.0
    lma(10)=50.0
    lma(11)=50.0
    lma(12)=16.7
!
! calculate sla values
    do i=1,12
      sla(i)=mwc/lma(i)
    enddo
!
! pick indeces from lowest to highest lma
    indx(1)=12
    indx(2)=6
    indx(3)=9
    indx(4)=5
    indx(5)=1
    indx(6)=10
    indx(7)=4

    istart=1
    istop=7
    do i=istart,istop
      j=indx(i)
!
! 2-segment integral approach
      lai=0.
      dmass=0.01
      nsteps=generic/dmass
      do n=1,nsteps
        if(LAI<=MorphTab(j)%LAImax) then
          sla_slope=-sla(j)*.5/MorphTab(j)%LAImax
          sla_int=sla(j)
          sla_loc=sla_slope*LAI+sla_int
        endif
        if(LAI>MorphTab(j)%LAImax) then
          sla_loc=0.5*sla(j)
        endif
        lai=lai+dmass*sla_loc
      enddo
!      test10(i)=lai
!      test10(1)=lai
!
! 3-segment integral approach
      lai=0.
      dmass=0.01
      nsteps=generic/dmass
      do n=1,nsteps
        if(LAI<=1.) then
          sla_loc=sla(j)
        endif
        if(LAI>1..and.lai<=MorphTab(j)%LAImax) then
          sla_slope=sla(j)*.5/(1.0-MorphTab(j)%LAImax)
          sla_int=sla(j)*(1.0-2.0*MorphTab(j)%LAImax)*.5/(1.0-MorphTab(j)%LAImax)
          sla_loc=sla_slope*LAI+sla_int
        endif
        if(LAI>MorphTab(j)%LAImax) then
          sla_loc=0.5*sla(j)
        endif
        lai=lai+dmass*sla_loc
      enddo
!      test10(i)=lai
!      test10(2)=lai
!
! constant lma approach
      lai_slope=mwc/lma(j)
      lai=lai_slope*generic
!      test10(i)=lai
!      test10(3)=lai
!
! segmented approach
      leaf_min=lma(j)/mwc
      leaf_max=MorphTab(j)%LAImax*lma(j)/mwc*0.5
      if(generic<=leaf_min) then
        lai_slope=mwc/lma(j)
      elseif(generic>leaf_min.and.generic<=leaf_max) then
        lai_slope=(generic+leaf_max-2*leaf_min)/(leaf_max-leaf_min)*mwc/lma(j)
      elseif(generic>leaf_max) then
        lai_slope=2.*mwc/lma(j)
      endif
      lai=lai_slope*generic
!      test10(i)=lai
!      test10(4)=lai
!
! continuous approach linear LMA
      leaf_max=MorphTab(j)%LAImax*lma(j)/mwc*0.5
      lai_slope=(generic/leaf_max+1.)*mwc/lma(j)
      lai_slope=min(2.,lai_slope)
      lai=lai_slope*generic
!      test10(i)=lai
!      test10(5)=lai
!
! exponential approach
      mass_max=0.7213*lma(j)/mwc*MorphTab(j)%LAImax
      if(generic<=mass_max) then
        lai=-0.6931*generic*mwc/lma(j)/MorphTab(j)%LAImax
        lai=MorphTab(j)%LAImax/(-0.6931)*alog(lai+1)
      else
        lai=MorphTab(j)%LAImax+4.*(generic-mass_max)*mwc/lma(j)
      endif
!      test10(i)=lai
!      test10(6)=lai
!
! corrected exponential approach
      sla_slope=log(0.5)/MorphTab(j)%LAImax
      mass_max=sla(j)/sla_slope*(exp(sla_slope*MorphTab(j)%LAImax)-1.)
      if(generic<=mass_max) then
!        lai=log(1.+sla_slope*generic/sla(j))/sla_slope
      else
        lai=MorphTab(j)%LAImax+(generic-mass_max)*sla(j)
      endif
      test10(i)=lai
!      test10(6)=lai
!
! exponential integral approach
      lai=0.
      dmass=0.01
      nsteps=generic/dmass
      do n=1,nsteps
        if(LAI<=MorphTab(j)%LAImax) then
          sla_slope=-sla(j)*.5/MorphTab(j)%LAImax
          sla_int=sla(j)
          sla_loc=sla(j)*exp(-0.6931/MorphTab(j)%LAImax*lai)
        endif
        if(LAI>MorphTab(j)%LAImax) then
          sla_loc=0.5*sla(j)
        endif
        lai=lai+dmass*sla_loc
      enddo
!      test10(i)=lai
!      test10(1)=lai

    enddo
    test=lai
!
      return
      end
!
!=======================================================================
      subroutine liquid_water_fraction(sib,test,test10)
!=========================================================================
! calculates the liquid water fraction in soils below freezing as a
! weighted average of clay, silt, sand, and organic matter liquid fractions
! output is liquid water fraction relative to fully saturated soil
!
! Modifications:
!  Kevin Schaefer created subroutine (1/17/14)
!--------------------------------------------------------------------------
!
  use kinds
  use sibtype
  use cfrax
  use sib_const_module
  use physical_parameters
!
  IMPLICIT NONE
!
! input/output variables
type(sib_t), intent(inout) :: sib
!
! define internal variables
  integer(kind=int_kind) :: j  ! index 
  real(kind=dbl_kind) :: sand_loc   ! (-) sand fraction of soil
  real(kind=dbl_kind) :: silt_loc   ! (-) silt fraction of soil
  real(kind=dbl_kind) :: clay_loc   ! (-) clay fraction of soil
  real(kind=dbl_kind) :: min_loc    ! (-) mineral fraction of soil 
  real(kind=dbl_kind) :: water_clay ! (m3/m3) silt fraction of soil 
  real(kind=dbl_kind) :: water_silt ! (m3/m3) silt fraction of soil 
  real(kind=dbl_kind) :: water_sand ! (m3/m3) silt fraction of soil 
  real(kind=dbl_kind) :: water_org  ! (m3/m3) silt fraction of soil
  real(kind=dbl_kind) :: liqfrac_crit_clay        ! (-) critical liquid fraction for transition from bulk to film freezing for clay  
  real(kind=dbl_kind) :: liqfrac_crit_silt        ! (-) critical liquid fraction for transition from bulk to film freezing for silt  
  real(kind=dbl_kind) :: liqfrac_crit_sand        ! (-) critical liquid fraction for transition from bulk to film freezing for sand  
  real(kind=dbl_kind) :: liqfrac_crit_org         ! (-) critical liquid fraction for transition from bulk to film freezing for org  
  real(kind=dbl_kind) :: min_clay        ! (-) minimum liquid fraction at -8 C for clay  
  real(kind=dbl_kind) :: min_silt        ! (-) minimum liquid fraction at -8 C for silt  
  real(kind=dbl_kind) :: min_sand        ! (-) minimum liquid fraction at -8 C for sand  
  real(kind=dbl_kind) :: min_org         ! (-) minimum liquid fraction at -8 C for org  
  real test, test10(10)     ! test variables
!
! switch current and previous liquid water fractions
    sib%prog%liqfrac1(:)=sib%prog%liqfrac2(:)
!
! loop through soil layers and calculate liquid water fraction
    do j=1,nsoil
!
! clay, silt, sand fractions of total soil volume
      min_loc=1._dbl_kind-sib%casa%org_frac(j)
      sand_loc=min_loc*sib%param%sandfrac/100._dbl_kind
      clay_loc=min_loc*sib%param%clayfrac/100._dbl_kind
      silt_loc=min_loc*(100._dbl_kind-sib%param%sandfrac-sib%param%clayfrac)/100._dbl_kind
!
! constituent liquid water fractions
      if(sib%prog%td(j)>tice) then
        water_clay=1._dbl_kind
        water_silt=1._dbl_kind
        water_sand=1._dbl_kind
        water_org=1._dbl_kind
      else
        water_clay=((sib%param%hfrzi-sib%prog%td(j))/sib%param%tstar_clay)**sib%param%beta_clay
        water_silt=((sib%param%hfrzi-sib%prog%td(j))/sib%param%tstar_silt)**sib%param%beta_silt
        water_sand=((sib%param%hfrzi-sib%prog%td(j))/sib%param%tstar_sand)**sib%param%beta_sand
        water_org=((sib%param%hfrzi-sib%prog%td(j))/sib%param%tstar_org)**sib%param%beta_org
      endif
!
! weighted average
      sib%prog%liqfrac2(j)=sand_loc*water_sand+silt_loc*water_silt+clay_loc*water_clay+sib%casa%org_frac(j)*water_org
      sib%prog%liqfrac2(j)=min(1.0_dbl_kind,sib%prog%liqfrac2(j))

      sib%diag%test2var1(j)=sib%prog%liqfrac2(j)
      sib%diag%test2var2(j)=silt_loc
      sib%diag%test2var2(j)=clay_loc
      test10(1)=water_clay
      test10(2)=water_silt
      test10(3)=water_sand
      test10(4)=water_org
      test10(5)=sib%prog%liqfrac2(j)
    enddo
!
! critical liquid water fractions
    liqfrac_crit_clay=((sib%param%hfrzi-tice)/sib%param%tstar_clay)**sib%param%beta_clay
    liqfrac_crit_silt=((sib%param%hfrzi-tice)/sib%param%tstar_silt)**sib%param%beta_silt
    liqfrac_crit_sand=((sib%param%hfrzi-tice)/sib%param%tstar_sand)**sib%param%beta_sand
    liqfrac_crit_org=((sib%param%hfrzi-tice)/sib%param%tstar_org)**sib%param%beta_org
!
! minimum liquid water
    min_clay=((sib%param%hfrzi-263.15)/sib%param%tstar_clay)**sib%param%beta_clay
    min_silt=((sib%param%hfrzi-263.15)/sib%param%tstar_silt)**sib%param%beta_silt
    min_sand=((sib%param%hfrzi-263.15)/sib%param%tstar_sand)**sib%param%beta_sand
    min_org=((sib%param%hfrzi-263.15)/sib%param%tstar_org)**sib%param%beta_org

    test10(1)=water_clay
    test10(2)=water_silt
    test10(3)=water_sand
    test10(4)=water_org
    test=min_org

  end subroutine liquid_water_fraction
!
!=======================================================================
      subroutine liquid_water_fraction_test
!=========================================================================
! liquid water fraction in soils below freezing
!
! Modifications:
!  Kevin Schaefer created subroutine (1/15/14)
!--------------------------------------------------------------------------
!
    use Mapper_Variables
    use sibscan_Variables
use kinds
use sibtype
use cfrax
use sib_const_module
use physical_parameters
!
    IMPLICIT NONE
!
! input/output variables
!type(sib_t), intent(inout) :: sib
!type(sib_local_vars), intent(inout) :: sib_loc ! variables local to SiB
!
! define internal variables
    real var1  ! temp variable
    real temp_frz   ! (-) respiration freezing scaling factor
    real casa_q10   ! (-) original casa Q10 scaling factor
    real num     ! (1/m3) number soil particles per volume soil
    real vol     ! (m3/m3) volume soil particles per volume soil
    real surf    ! (m2/g) specific surface area soil particles per volume soil
    real water   ! (m3/m3) volume liquid water per total soil volume
    real alpha_loc  ! (-) amplitude
    real beta_loc   ! (-) liq water exponent
    real Tstar   ! (deg C) temperature below Oc where freezing begins 
    real sand_loc,silt_loc, clay_loc, min_loc   ! (%) silt fraction of soil 
    real tstar_clay, beta_clay, water_clay   ! (%) silt fraction of soil 
    real tstar_silt, beta_silt, water_silt   ! (%) silt fraction of soil 
    real tstar_sand, beta_sand, water_sand   ! (%) silt fraction of soil 
    real tstar_org,  beta_org,  water_org   ! (%) silt fraction of soil 
    real tstar_combo,  beta_combo,  water_combo   ! (%) silt fraction of soil 
    real water_ave   ! (%) silt fraction of soil 
    real delta_r  ! (m) thin film radius
!
! number soil particles per volume (square packing)
    var1=radius*1.e-6
    num=1./(2.*var1)**3.
    vol=4./3.*3.14159*var1**3.*num
    poros=1.-vol
    surf=4.*3.14159*var1*var1/den_min/1000.*num
!
! clay, silt, sand formulation
    sib(1)%param%hfrzi=273.15
    sand_loc=sand/100.
    clay_loc=clay/100.
    silt_loc=1.-sand_loc-clay_loc
    min_loc=1.-f_om
    tstar_clay=0.3
    tstar_silt=0.1
    tstar_sand=0.01
    tstar_org =0.01
    beta_clay=-0.3
    beta_silt=-0.5
    beta_sand=-0.9
    beta_org =-1.0
! clay
    tstar=tstar_clay
    beta_loc=beta_clay   
    if(td(1)<sib(1)%param%hfrzi-tstar) then
      water=((sib(1)%param%hfrzi-td(1))/tstar)**beta_loc
    else
      water=1.
    endif
    water_clay=water
! silt
    tstar=tstar_silt
    beta_loc=beta_silt   
    if(td(1)<sib(1)%param%hfrzi-tstar) then
      water=((sib(1)%param%hfrzi-td(1))/tstar)**beta_loc
    else
      water=1.
    endif
    water_silt=water
! sand
    tstar=tstar_sand
    beta_loc=beta_sand   
    if(td(1)<sib(1)%param%hfrzi-tstar) then
      water=((sib(1)%param%hfrzi-td(1))/tstar)**beta_loc
    else
      water=1.
    endif
    water_sand=water
! organic
    tstar=tstar_org
    beta_loc=beta_org   
    if(td(1)<sib(1)%param%hfrzi-tstar) then
      water=((sib(1)%param%hfrzi-td(1))/tstar)**beta_loc
    else
      water=1.
    endif
    water_org=water
!
! weighted average
    water_ave=min_loc*(sand_loc*water_sand+ silt_loc*water_silt+clay_loc*water_clay)+f_om*water_org
!
! thin film radius
    delta_r=(water_ave*poros+vol)/vol
    delta_r=(delta_r**.333-1.)*radius
!
! combined alpha/beta
    beta_combo =min_loc*(sand_loc*beta_sand+ silt_loc*beta_silt+clay_loc*beta_clay)+f_om*beta_org
    tstar_combo= min_loc*(sand_loc*tstar_sand+silt_loc*tstar_silt+clay_loc*tstar_clay)+f_om*tstar_org

    beta_loc=beta_combo
    tstar=tstar_combo
    if(td(1)<sib(1)%param%hfrzi-tstar) then
      water=poros*((sib(1)%param%hfrzi-td(1))/tstar)**beta_loc
    else
      water=poros
    endif
    water_combo=water
!
! liquid water fraction parameters
    alpha_loc=den_min*vol/den_liq*exp(0.5519*log(surf)+0.2168)/100.
    beta_loc =-exp(-0.264*log(surf)+.3711)
!
! liquid water content kurylyk Watanabi [2013]
    sib(1)%param%hfrzi=273.15
    if(td(1)<sib(1)%param%hfrzi) then
      water=alpha_loc*(sib(1)%param%hfrzi-td(1))**beta_loc
      water=min(water,poros)
    else
      water=poros
    endif
!
! Romanovsky Osterkamp [2000] liquid water fraction for barrow
    sib(1)%param%hfrzi=273.15
    water=alpha*(sib(1)%param%hfrzi-td(1))**beta
    water=((sib(1)%param%hfrzi-td(1))/0.1)**beta
    water=min(water,poros)
    test=water
!
! Romanovsky Osterkamp [2000] liquid water fraction for barrow
    sib(1)%param%hfrzi=273.15
!
! peat 0.0-0.2 m
    water=0.026*(sib(1)%param%hfrzi-td(1))**(-0.38)
    water=min(water,poros)
!
! silt 0.2-0.3 m
    water=0.12*(sib(1)%param%hfrzi-td(1))**(-0.5)
    water=min(water,poros)
!
! silt 0.3-1.0 m
    water=0.064*(sib(1)%param%hfrzi-td(1))**(-0.38)
    water=min(water,poros)
!
! liquid water fraction for simple soil model
    if(td(1)>(tf)) then ! liquid water
      fw=1.
      fw_slope=0.
    else if (td(1)<(tf-1.)) then ! ice
      fw=0.
      fw_slope=0.
    else ! mix between ice and water
      fw=td(1)-tf+1.
      fw_slope=1.
                
      fw=1.581976707*exp(-(td(1)-tf)**2.)-0.581976707
      fw_slope=1.581976707*exp(-(td(1)-tf)**2.)*(-2.)*(td(1)-tf)

      fw=1./(1.+EXP(20.*(tf-td(1)-.5)))
      fw=-2.*(td(1)-tf)**3.-3.*(td(1)-tf)**2.+1.
      fw_slope=-6.*(td(1)-tf)**2.-6.*(td(1)-tf)
    end if
!
! Q10 scaling factor
    soilQ10=exp(log(Q10)/10.*(td(1)-298.15))
!
! CASA temperature respiration scaling factor
    casa_q10=exp(log(2.)/10.*(td(1)-303.15))
!
! frozen soil respiration inhinbition function
! sfrzi and hfrzi defined in sib_const_module.f90
    sfrzi=generic
    sfrzi=1.
    hfrzi=271.
    temp_frz=1./(1.+exp(sfrzi*(hfrzi-td(1))))
!
! frozen soil exp decrease scaling factor based on Mikan et al. [2002]
    sib(1)%param%sfrzi=2.
    sib(1)%param%hfrzi=273.15
    sib(1)%param%tcfrzi=-5.*log(10.)/sib(1)%param%sfrzi+sib(1)%param%hfrzi
    if(td(1)>sib(1)%param%tcfrzi) then
      temp_frz=exp(sib(1)%param%sfrzi*(td(1)-sib(1)%param%hfrzi))
      temp_frz=min(temp_frz,1.)
    else
      temp_frz=0.
    endif
!    test10(1)=temp_frz
!
! Freeze factor formulation exactly based on Mikan et al [2002] Incubation curve fits
    sib(1)%param%sfrzi=.5
    sib(1)%param%hfrzi=273.15
    sib(1)%param%tcfrzi=-5.*log(10.)/sib(1)%param%sfrzi+sib(1)%param%hfrzi
    if(td(1)>sib(1)%param%tcfrzi) then
      temp_frz=exp(sib(1)%param%sfrzi*(td(1)-sib(1)%param%hfrzi))
      temp_frz=min(temp_frz,1.)
    else
      temp_frz=0.
    endif
!
! Romanovsky Osterkamp [2000] liquid water fraction
    sib(1)%param%hfrzi=273.15
    temp_frz=0.12*(sib(1)%param%hfrzi-td(1))**(-0.5)
    temp_frz=min(temp_frz,1.)
!    test10(3)=temp_frz
!
! Monson et al. [2006]
    sib(1)%param%hfrzi=273.15
    temp_frz=exp(log(105./10.)*(td(1)-sib(1)%param%hfrzi))
    temp_frz=min(temp_frz,1.)
!    test10(4)=temp_frz
!
    return
    end
!
!=======================================================================
      subroutine external_carb_in_out(sib,test,test10)
!=========================================================================
! calculates internal inputs and outputs of carbon to the system of pools
! representing redistribution of carbon due to an external process
! or disturbance.  The external i/o currently includes:
!  1) input to the storage pool due to photosynthesis
! this is not the final version that went into SiBCASA

! liquid water fraction of total soil volume in soils below freezing
!
! Modifications:
!  Kevin Schaefer created subroutine (1/16/14)
!--------------------------------------------------------------------------
!
  use kinds
  use sibtype
  use cfrax
  use sib_const_module
  use physical_parameters
!
  IMPLICIT NONE
!
! input/output variables
type(sib_t), intent(inout) :: sib
!
! define internal variables
  integer(kind=int_kind) :: j  ! index 
  real(kind=dbl_kind) :: Tstar      ! (K) freezing point depression 
  real(kind=dbl_kind) :: beta       ! (-) liq water exponent
  real(kind=dbl_kind) :: sand_loc   ! (-) sand fraction of soil
  real(kind=dbl_kind) :: silt_loc   ! (-) silt fraction of soil
  real(kind=dbl_kind) :: clay_loc   ! (-) clay fraction of soil
  real(kind=dbl_kind) :: min_loc    ! (-) mineral fraction of soil 
  real(kind=dbl_kind) :: water_clay ! (m3/m3) silt fraction of soil 
  real(kind=dbl_kind) :: water_silt ! (m3/m3) silt fraction of soil 
  real(kind=dbl_kind) :: water_sand ! (m3/m3) silt fraction of soil 
  real(kind=dbl_kind) :: water_org  ! (m3/m3) silt fraction of soil 
  real(kind=dbl_kind) :: water_ave  ! (m3/m3) silt fraction of soil 
  real test, test10(10)     ! test variables
!
! loop through soil layers
  do j=1,nsoil
!
! clay, silt, sand fractions of total soil volume
    min_loc=1.-sib%casa%org_frac(j)
    sand_loc=min_loc*sib%param%sandfrac/100.
    clay_loc=min_loc*sib%param%clayfrac/100.
    silt_loc=min_loc*(100.-sib%param%sandfrac-sib%param%clayfrac)/100.
!
! clay
    tstar=sib%param%tstar_clay
    beta=sib%param%beta_clay   
    if(sib%prog%td(j)<sib%param%hfrzi-tstar) then
      water_clay=sib%param%poros(j)*((sib%param%hfrzi-sib%prog%td(j))/tstar)**beta
    else
      water_clay=sib%param%poros(j)
    endif
!
! silt
    tstar=sib%param%tstar_silt
    beta=sib%param%beta_silt   
    if(sib%prog%td(j)<sib%param%hfrzi-tstar) then
      water_silt=sib%param%poros(j)*((sib%param%hfrzi-sib%prog%td(j))/tstar)**beta
    else
      water_silt=sib%param%poros(j)
    endif
!
! sand
    tstar=sib%param%tstar_sand
    beta=sib%param%beta_sand   
    if(sib%prog%td(j)<sib%param%hfrzi-tstar) then
      water_sand=sib%param%poros(j)*((sib%param%hfrzi-sib%prog%td(j))/tstar)**beta
    else
      water_sand=sib%param%poros(j)
    endif
!
! organic
    tstar=sib%param%tstar_org
    beta=sib%param%beta_org   
    if(sib%prog%td(j)<sib%param%hfrzi-tstar) then
      water_org=sib%param%poros(j)*((sib%param%hfrzi-sib%prog%td(j))/tstar)**beta
    else
      water_org=sib%param%poros(j)
    endif
!
! weighted average
    water_ave=sand_loc*water_sand+ silt_loc*water_silt+clay_loc*water_clay+sib%casa%org_frac(j)*water_org

    test = water_clay
    test10(1)=water_clay
    test10(2)=water_silt
    test10(3)=water_sand
    test10(4)=water_org
    test10(5)=water_ave

  enddo

  return
  end
!
!
!
!=======================================================================
      subroutine switch_biome_table
!=========================================================================
! reads input parameters for the scan program
!
! Modifications:
!  Kevin Schaefer created routine (10/23/13)
!--------------------------------------------------------------------------
!
      use Mapper_Variables
      use sibscan_Variables
!
      IMPLICIT NONE
!
integer i,k  ! indeces
character*20 phystype
!
! determine phys type
   phystype='c3'
   if(bionum==6)  phystype='c4'
   if(bionum==7)  phystype='c4'
   if(bionum==8)  phystype='c4'
   if(bionum==11) phystype='c4'
!
! set physfrac
   if(phystype=='c3') then
     sib(1)%param%physfrac(1) = 1.
     sib(1)%param%physfrac(2) = 0.
     sib%param%phystype(1)=3.
     sib%param%phystype(2)=4.
   elseif(phystype=='c4') then
     sib(1)%param%physfrac(1) = 1.
     sib(1)%param%physfrac(2) = 0.
     sib%param%phystype(1)=4.
     sib%param%phystype(2)=3.
   endif
!
! read look-up table of biome dependant variables
   if(phystype=='c3') then
     k=0
     do i=1,numfiles
       if(files(i)%type=='c3_biome_tab') k=i
     enddo
     if(k==0) then
       print*, 'Error: no c3_biome_tab specified'
       stop
     endif
     Filename%BioTab=files(k)%path
   elseif(phystype=='c4') then
     k=0
     do i=1,numfiles
       if(files(i)%type=='c4_biome_tab') k=i
     enddo
     if(k==0) then
       print*, 'Error: no c4_biome_tab specified'
       stop
     endif
     Filename%BioTab=files(k)%path
   endif
   call ReadBioTable

   return
   end
!

