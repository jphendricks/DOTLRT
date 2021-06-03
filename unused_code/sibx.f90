!
!=======================================================================
      subroutine aero (zlt,fv,chil, &
       z2,z1,zc,lw,ll, &
       cc1,cc2,z0,d)                                      
!=======================================================================       
!     Two options for calculation of aerodynamic parameters:
!     Modified version of SIBX, or simple equations.
!     Source of simplified equations (Niall Hanan):
!     r(b)-Choudhury and Monteith, 1988: Q.J.R.Meteorological
!     Soc.114,373-398.
!     r(a_incanopy)-derived from: Shuttleworth and Wallace, 1985:
!     Q.J.R.Meteorological Soc. 111, 839-855.
!     Dolman, 1993: Ag. For. Met. 65, 21-45.
!     These eqns assume that K-theory holds within the canopy,
!     which is not quite true but is a good first approximation (sibx
!     makes similar assumptions).
!
!     Hanan has substituted r(d) for ra_low== in-canopy aerodynamic 
!     resistance for a 1-layer model with bare soil beneath. r(d) has
!     essentially the same definition.  The original ra_low differs in
!     that it allows for a second vegetation layer under the first
!     layer. This is a special case (LAI and roughness of
!     the lower layer approach zero i.e. no second layer only soil).      
!-----------------------------------------------------------------------
!
      implicit none
!
      integer, parameter :: nv=13
      real*8 dcc1, dcc2, dz0, dd
      real cc1, cc2, z0, d
      real zwind     ! Reference height for momentum, meters
      		     ! In practice, zwind is the height of measurement
      real zmet      ! Reference height for weather, meters
      		     ! In practice, the height of measurement
      real z2	     ! Canopy height
      real z1	     ! Canopy base
      real zc        ! Inflection point for leaf density
      real zs        ! Soil roughness length, meters
      real lw        ! Leaf width
      real ll        ! Leaf length
      real zlt       ! Leaf-area index
      real fv        ! Fractional vegetation cover
      real chil      ! Leaf angle distribution factor, or Xl
      real*8 dzwind  ! Reference height for momentum, meters
      		     ! In practice, zwind is the height of measurement
      real*8 dzmet   ! Reference height for weather, meters
      		     ! In practice, the height of measurement
      real*8 dz2     ! Canopy height
      real*8 dz1     ! Canopy base
      real*8 dzc     ! Inflection point for leaf density
      real*8 dzs     ! Soil roughness length, meters
      real*8 dlw     ! Leaf width
      real*8 dll     ! Leaf length
      real*8 dzlt    ! Leaf-area index
      real*8 dfv     ! Fractional vegetation cover
      real*8 dchil   ! Leaf angle distribution factor, or Xl
      real dummy(13,17), dlook(13,17),zlook(13,17) 
      real rin,prop
      integer inp, in, iden, iv, ivtype
      integer switch
!
      switch=1
!
      if (switch.eq.1) then
!
!----------------------------------------------------------------------
! Sibx
!----------------------------------------------------------------------
! Zwind and zmet values are the measurement heights for wind 
! weather.  They must be higher than the canopy top (Z2).  Set them
! for local domains, or use values of 80-100 for GCM.
! Zs value from Sellers et al Part II pp. 723.
! G1 and Zt/Z0 are prescribed in SIBX, 
! values from Sellers et al. Part I
! Also in SIBX, set IGCM=0 or 1.  This affects the 
! resistance values.
!
      zwind=100.
      zmet=100.
      zs=0.05
      dzwind=dble(zwind)
      dzmet=dble(zmet)
      dzs=dble(zs)
      dz2=dble(z2)
      dz1=dble(z1)
      dzc=dble(zc)
      dlw=dble(lw)
      dll=dble(ll)
      dzlt=dble(zlt)
      dfv=dble(fv)
      dchil=dble(chil)
      dcc1=0.d0
      dcc2=0.d0
      dz0=0.d0
      dd=0.d0
!
      call sibx(dzwind, dzmet, dz2, dzc, dz1, dzs, dzlt, dfv, &
       dll, dlw, dchil, dcc1, dcc2, dz0, dd)
!
      cc1=dcc1
      cc2=dcc2
      z0=dz0
      d=dd
!
      else if (switch.ne.1) then
!
!----------------------------------------------------------------------
! Hanan's equations
!----------------------------------------------------------------------
! the eddy decay coefficient is taken to be ~ 2.5
! von karman's constant is 0.41
!                                                                               
      open(unit=13, file='sibaero.data', form='formatted')
!
      do iv=1, nv                                                        
      read(13, *) ivtype                                                      
      read(13, *) (zlook (iv,iden), iden=1, 17)
      read(13, *) (dlook (iv,iden), iden=1, 17)
      read(13, *) (dummy (iv,iden), iden=1, 17)
      read(13, *) (dummy (iv,iden), iden=1, 17)         
      enddo
!
      rin= (zlt+0.01)*2.+1.                                      
      inp=int(rin+1.)                                                         
      inp=min0 (17, inp)                                                    
      in=inp-1                                                             
      prop=rin-in*1.0
!
! Interpolate d, since we can't diagnose it without sibx.
      d=dlook(iv,in)+prop*(dlook(iv,inp)-dlook(iv,in))
!
! Hanan's equations:
      cc1=(100*2.5/(2*zlt))*(sqrt(lw)/(1-exp(-1.25)))
      cc2=z2*exp(2.5)*(1-exp(-2.5*(zs+d)/z2))
      cc2=cc2*alog((z2-d)/zs)/(2.5*0.41*0.41)
!
! z0 is from Sellers et al Part II Equation 6
! only valid when 0.5 < zlt < 5.0
! for zlt > 5.0, use sibx
!      z0=z2*(1-0.91*exp(-0.0075*zlt))
! Or, interpolate z0 using sibx values
       z0=zlook(iv,in)+prop*(zlook(iv,inp)-zlook(iv,in))
      endif
!
      return                                                                    
      end
!======================================================================
	subroutine sibx(ZWIND, ZMET, Z2, ZC, Z1, ZS, ZLT, VCOV1, &
                    ZLEN, ZLW, CHIL, RBC, RDC, Z0, D )
!======================================================================
!
!     PROGRAM CALCULATES AERODYNAMIC PROPERTIES OF CANOPIES USING
!     TRIANGULAR DISTRIBUTION OF LEAF AREA DENSITY.
!
!======================================================================
!
!     (1)  FURTHER DOCUMENTATION IN SUBROUTINE MOMOPT.
!
!     (2)  VERSION CALCULATES G3 IN THE EVENT THAT ZMET IS DIFFERENT
!          FROM ZWIND.
!
!     (3)  GCM-COMPATIBLE PARAMETERS PRODUCED BY SETTING VARIABLE
!          GCM=1. OTHERWISE SET GCM=0, E.G. WHEN USING RASITE IN 1-D
!          STUDIES.
!
!----------------------------------------------------------------------
!
!     REFERENCES
!     ----------
!
!         SELLERS P.J. , Y. MINTZ, Y.S. SUD, A. DALCHER (1986) 'A SIMPLE
!         BIOSPHERE MODEL (SIB) FOR USE WITHIN GENERAL CIRCULATION
!         MODELS', ATMOS. SCI., 43, 6, 505-531.
!
!         SELLERS P.J. , W.J. SHUTTLEWORTH, J.L. DORMAN, A. DALCHER,
!         J.M. ROBERTS (1989) ' CALIBRATING THE SIMPLE BIOSPHERE
!         MODEL (SIB) FOR AMAZONIAN TROPICAL FOREST USING FIELD AND
!         REMOTE SENSING DATA. PART 1 , AVERAGE CALIBRATION WITH FIELD
!         DATA ', ( SEE APPENDIX I ), J. APPL. MET.
!
!         SHAW R.H. , A.R. PEREIRA (1982) 'AERODYNAMIC ROUGHNESS OF A
!         CANOPY: A NUMERICAL EXPERIMENT', AGRIC MET., 26, 51-65.
!
!-----------------------------------------------------------------------
!
!     INPUT VARIABLES passed from mapper
!
!     In sibx   In mapper           
!
!     Z2        veg%z2d        Canopy top height
!     Z1        veg%z1d        Canopy base height
!     ZC        veg%zcd        Canopy inflection  height
!     ZLT       fields%zltf    Leaf area index LAI
!     VCOV1     fields%fvf     Fractional vegetation cover
!     ZLEN      veg%lld        Leaf length
!     ZLW       veg%lwd        Leaf width
!     CHIL      veg%child      Leaf angle distribution factor

!     INPUT CONSTANTS set in subroutine aero
!
!     ZWIND   Measurement height for wind
!     ZMET    Measurement height for weather
!     ZS      Soil roughness length

!     CONSTANTS set in sibx below
!
!     G1      Ratio of eddy diffusivity extrapolated from a log-linear
!             profile (Km) to augmented eddy diffusivity at the 
!             canopy top.
!
!     ZTZ0    Also G4.  Ratio of ZT to Z0.  ZT is the transition
!             height between the perturbed surface layer and the
!             inertial boundary layer. Z) is the canopy roughness
!             length.

!     OUTPUT VARIABLES
!
!     RBC     Resistance to heat transfer between leafs and the 
!             canopy airspace.
!
!     RDC     Resistance to heat transfer between the ground and 
!             the canopy airspace.
!
!     Z0      Canopy roughness length
!
!     D       Zero plane displacement

!    INTERMEDIATE VARIABLES
!
!    PS       Leaf shelter factor
!    REYNO    Reynolds number
!    CDD      Leaf drag coefficient
!    DL       Mean leaf-area density
!    CS       Heat-mass transfer coefficient
!    HA       Canopy source height


       IMPLICIT REAL*8 (A-H, O-Z)

	G1 = 1.449
	ZTZ0 = 11.785
	IGCM = 1

      FLOWB = 1./3.141592 * ( 1.-CHIL)
      FLOWB = DMAX1( FLOWB, 0.00001D0)
      REYNO = ( ZLW+ZLEN )/2. * 10000. / 0.15
      CDD= 1.328*2. / DSQRT(REYNO) + 0.45*FLOWB**1.6
      DL = (ZLT / ( Z2-Z1 ))*VCOV1
       PS = 1 + DL**0.6
      CS = 90. * DSQRT(dble(ZLW))
!
!-----------------------------------------------------------------------
!
      
      CALL MOMOPT ( Z2, ZC, Z1, ZS, ZWIND, ZMET, ZLT, G1, ZTZ0, &
                    CDD, PS, CS, &
                    Z0, D, RBC, RDC, &
                    HA, G2, G3, CORB1, CORB2, U2FAC )
!
!-----------------------------------------------------------------------
!
      IF ( IGCM .EQ. 0 ) U2FAC = 1.
      RBC = RBC / DSQRT(U2FAC)
      RDC = RDC / U2FAC
      Z0Z2 = Z0/Z2
      DZ2 = D/Z2
      HAZ2 = HA/Z2
      USTARU = 0.41 / DLOG ( ( ZWIND-D ) / Z0 )
 
      RETURN
      END

!=======================================================================
      SUBROUTINE MOMOPT ( Z2, ZC, Z1, ZS, ZWIND, ZMET, ZLT, G1, ZTZ0, &
                          CDD, PS, CS, &
                          Z0, D, RBC, RDC, &
                          HA, G2, G3, CORB1, CORB2, U2FAC )
!======================================================================
!
!     SUBROUTINE TO CALCULATE AERODYNAMIC PARAMETERS FOR A VEGETATED
!     SURFACE. THEORETICAL TREATMENT FOLLOWS FIRST-ORDER CLOSURE MODEL
!     DESCRIBED IN SELLERS ET AL. (1986) AND UPDATED IN SELLERS ET AL.
!     (1989). LOGLINEAR PROFILE IS ASSUMED ABOVE CANOPY EXCEPT BELOW
!     TRANSITION HEIGHT, ZL, WHERE DIFFUSION COEFFICIENT IS ENHANCED.
!     VALUES OF G1 AND ZTZ0 WERE OBTAINED BY FORCING THE RESULTS
!     TO FIT THOSE OF SHAW AND PEREIRA ( 1982): G1=1.449, ZTZ0=11.785
!
!-----------------------------------------------------------------------
!
!    (1) INPUT PARAMETERS ARE DEFINED IN SELLERS ET AL. (1989) AND IN
!        PROGRAM SIBX, ALSO OUTPUT PARAMETERS EXCEPT FOR IN (2) BELOW.
!
!    (2) DISTORTION OF PROFILE NEAR SURFACE NECESSITATES CALCULATION OF
!        THE FOLLOWING PARAMETERS.
!
!      G1, G2, G3, ZTZ0, CORB1, CORB2, HA
!
!      G1     : RATIO OF KM(ACTUAL) TO KM(LOG-LINEAR) AT Z = Z2
!      G2     : RATIO OF RA(ACTUAL) TO RA(LOG-LINEAR) FOR MOMENTUM
!               BETWEEN: Z = Z2 AND Z = ZX, WHERE ZX = MIN(ZL,ZWIND)
!      G3     : RATIO OF RA(ACTUAL) TO RA(LOG-LINEAR) FOR HEAT
!               BETWEEN: Z = Z2 AND Z = ZX, WHERE ZX = MIN(ZL,ZMET)
!      ZTZ0   : PARAMETER TO DETERMINE DEPTH OF TRANSITION LAYER ABOVE
!               CANOPY, ZL. ZL = Z2 + ZTZ0 * Z0
!      CORB1  : NON-NEUTRAL CORRECTION FOR CALCULATION OF AERODYNAMIC
!               RESISTANCE BETWEEN HA AND Z2. WHEN MULTIPLIED BY
!               H*RBB/TM GIVES BULK ESTIMATE OF LOCAL RICHARDSON NUMBER.
!               RBB = RA FOR HEAT BETWEEN HA AND Z2.
!               CORB2 = 9*G/( RHOAIR*CPAIR* (DU/DZ)**2 )
!      CORB2  : NEUTRAL VALUE OF RBB*U2 ( SQUARED ), EQUIVALENT TO
!               RDC**2 FOR UPPER CANOPY
!      HA     : CANOPY SOURCE HEIGHT FOR HEAT
!      U2FAC  : U2(ACTUAL) / U2(LOG-LINEAR ESTIMATE). THIS IS ONLY USED
!               FOR GCM APPLICATIONS.
!
!    (3) OUTPUT IS COMPATIBLE WITH SUBROUTINE RASITE
!
!-----------------------------------------------------------------------
!
!     REFERENCES
!     ----------
!
!         SELLERS P.J. , Y. MINTZ, Y.S. SUD, A. DALCHER (1986) 'A SIMPLE
!         BIOSPHERE MODEL (SIB) FOR USE WITHIN GENERAL CIRCULATION
!         MODELS', ATMOS. SCI., 43, 6, 505-531.
!
!         SELLERS P.J. , W.J. SHUTTLEWORTH, J.L. DORMAN, A. DALCHER,
!         J.M. ROBERTS (1989) ' CALIBRATING THE SIMPLE BIOSPHERE
!         MODEL (SIB) FOR AMAZONIAN TROPICAL FOREST USING FIELD AND
!         REMOTE SENSING DATA. PART 1 , AVERAGE CALIBRATION WITH FIELD
!         DATA ', ( SEE APPENDIX I ), J. APPL. MET.
!
!         SHAW R.H. , A.R. PEREIRA (1982) 'AERODYNAMIC ROUGHNESS OF A
!         CANOPY: A NUMERICAL EXPERIMENT', AGRIC MET., 26, 51-65.
!
!
!-----------------------------------------------------------------------

      IMPLICIT REAL*8 (A-H, O-Z)
      
      DIMENSION ZMAT(4,5), WORK(4,5), COEFS(4)
      DATA VKC/0.41/

      PZT = 1.
      ZT = Z1
      CDG= (VKC * PZT / DLOG(ZT/ZS)) ** 2 / (2.0 - PZT)

      CALL DENCAL(Z2, Z1, ZC, ZLT, A0L, B0L, A0U, B0U)

!     Set up simple bisection.
!     rl is the low guess (dependent variable is negative)
!     rh is the high guess (dependent variable is positive)

!     Make sure root is bracketed by rl and rh

      rl = 1.0d0*Z2
      rh = 0.001d0*Z2
      SIGMA = 0.01d0*Z2
      dxold = rh - rl
      dx=dxold
      
      MAXIT = 200
      xacc = 0.0001d0*Z2
      
!    SIGMA is effectively a mixing length, which
!    must be tuned iteratively below.  In general,
!    values for SIGMA, D, and Z0 should be of order Z2.
!    See Shaw and Pereira 1982 for help with scaling.
!
!    SIGMA ITERATION LOOP STARTS HERE
!
      do j=1,MAXIT
      
!     These are the coefficients of the Bessel
!     equation Ai and Bi.     

      AT = 2.*A0U * CDD / (PS*SIGMA)
      BT = 2.*B0U * CDD / (PS*SIGMA)
      AL = 2.*A0L * CDD / (PS*SIGMA)
      BL = 2.*B0L * CDD / (PS*SIGMA)

!     WINVAL solves the Bessel equation for Ya and Yb                 

      CALL WINVAL(Z2, AT, BT, ALPHA, BETA)

      ZMAT(1,1) = ALPHA
      ZMAT(1,2) = BETA
      
      ZMAT(1,5) = 1.

      CALL WINVAL(ZC, AT, BT, ALPHA, BETA)

      ZMAT(2,1) = ALPHA
      ZMAT(2,2) = BETA

      CALL WINVAL(ZC, AL, BL, ALPHA, BETA)

      ZMAT(2,3) = -ALPHA
      ZMAT(2,4) = -BETA
      
!     GRADVA solves the Bessel equation for 
!     dYa/dz and dYb/dz

      CALL GRADVA(ZC, AT, BT, ALPHA, BETA)

      ZMAT(3,1) = ALPHA
      ZMAT(3,2) = BETA
      
      CALL GRADVA(ZC, AL, BL, ALPHA, BETA)

      ZMAT(3,3) = -ALPHA
      ZMAT(3,4) = -BETA
      
      CALL WINVAL(Z1, AL, BL, ALPHA, BETA)
      CALL GRADVA(Z1, AL, BL, ALP, BET)

      ZMAT(4,3) = SIGMA*ALP - CDG*ALPHA
      ZMAT(4,4) = SIGMA*BET - CDG*BETA
      ZMAT(4,5) = 0.
      
      N   = 4
      NP1 = 5
      
!     GAUSSD solves a matrix of velocity and the vertical
!     gradient of velocity with boundary conditions.
!     Outputs the coefficients alpha and beta.

      CALL GAUSSD(ZMAT, N, NP1, COEFS, WORK)
      
!     These two expressions of the zero plane displacement
!     must converge at Z=Z2 for a root value which is SIGMA.
      
      CALL DEVAL1(COEFS, AT, BT, Z2, G1, SIGMA, VKC, D1)

      CALL DEVAL2(COEFS, AT, BT, AL, BL, A0U, B0U, A0L, B0L, ZS, Z2, &
                 ZC, Z1, PS, CDD, SIGMA, D2 )

      Y = (D1 - D2)
            
      if (dabs(Y).lt.xacc) goto 100
      
!     Keep root bracketed 
      
      if (Y.lt.0.) then
        rl = SIGMA
      else
        rh = SIGMA
      endif      

!     Take a bisection step to find the root of Y.

      dxold=dx
      dx=0.5*(rh-rl)
      SIGMA=rl+dx
      
      SIGMA = DMAX1( SIGMA, 0.0000001D0 )
           
      enddo

100   D = (D1+D2) / 2.
!
!-----------------------------------------------------------------------
!     ITERATIVE CALCULATION OF ZO AS A FUNCTION OF G2
!     CALCULATION OF G2 ASSUMES LINEAR VARIATION OF KM BETWEEN Z2 AND
!     ZL SO THAT :
!                   KM(ACTUAL) = G1* KM(LOG-LINEAR) AT Z = Z2
!                   KM(ACTUAL) =     KM(LOG-LINEAR) AT Z = ZL
!
!     INTEGRATION FOLLOWS EQUATIONS (A.10) AND (A.11) IN APPENDIX A OF
!     SELLERS ET AL. ( 198X) ' CALIBRATING THE SIB MODEL FOR TROPICAL
!     FOREST...'.
!-----------------------------------------------------------------------
!     Now solve for Z0 using bisection method
                
      rl = 0.50d0*Z2
      rh = 0.01d0*Z2
      Z0 = 0.05d0*Z2
      dxold = rl - rh
      dx=dxold
      
      MAXIT = 200
      xacc = 0.00001d0
      
      do j=1,MAXIT
      
      ZL = Z2 + ZTZ0*Z0
      AA = 1.-G1
      BB = G1*ZL - Z2 + D*(G1-1.)
      CC = -D* ( G1*ZL - Z2 )
      BIGH = CC/AA - BB*BB/(4.*AA*AA)
      SQHH = DSQRT( DABS(BIGH) )
      V2 = Z2 + BB/(2.*AA)
      VX = ZWIND + BB/(2.*AA)
      IF(ZL .LT. ZWIND) VX = ZL + BB/(2.*AA)
      IF ( BIGH .LT. 0. ) GO TO 200
      ARG1 = DATAN ( VX/SQHH ) / SQHH
      ARG2 = DATAN ( V2/SQHH ) / SQHH
      ISOL = 1
      GO TO 300
200   ARG1 = DLOG ( DABS(VX-SQHH) / (VX+SQHH) ) / ( 2.*SQHH )
      ARG2 = DLOG ( DABS(V2-SQHH) / (V2+SQHH) ) / ( 2.*SQHH )
      ISOL = 2
300   ARG3 = 0.
      IF(ZL .LT. ZWIND) ARG3 = DLOG ( (ZWIND-D)/(ZL-D) )
      RASTUB = ( ZL-Z2 )/AA* ( ARG1 - ARG2 ) + ARG3

      ARGEX = RASTUB + G1*VKC*VKC*(Z2-D)/SIGMA
      Z0F = (ZWIND-D) * DEXP(-ARGEX)
      Z0F = DMAX1( Z0F, 0.00001D0 )
!
      Y = Z0F - Z0
      
      if (dabs(Y).lt.xacc) goto 330 
      
      if (Y.lt.0.) then
        rl = Z0F
      else
        rh = Z0F
      endif

      dxold=dx
      dx=0.5*(rh-rl)
      Z0=rl+dx      
      
      Z0 = DMAX1( Z0, 0.00001D0 )
      
      enddo
      
330   Z0 = Z0F

      TOPG2 = RASTUB - ARG3
      ZUP = ZWIND
      IF ( ZL .LT. ZWIND ) ZUP = ZL
      BOTG2 = DLOG( (ZUP-D)/(Z2-D) )
      G2 = TOPG2/BOTG2
!
!-----------------------------------------------------------------------
!     CALCULATION OF G3 ASSUMES LINEAR VARIATION OF KH BETWEEN Z2 AND
!     ZL IN EXACTLY THE SAME WAY AS FOR G2. CALCULATION MUST BE REPEATED
!     AS ZMET MAY NOT BE EQUAL TO ZWIND.
!       N.B.:
!            IF ZMET = ZWIND ,        G2 = G3.
!            IF ZMET, ZWIND .GT. ZL , G2 = G3
!
!     ITERATION UNNECESSARY AS ZL HAS BEEN DETERMINED IN G2 CALCULATION
!-----------------------------------------------------------------------
!
      AA = 1.-G1
      BB = G1*ZL - Z2 + D*(G1-1.)
      CC = -D* ( G1*ZL - Z2 )
      BIGH = CC/AA - BB*BB/(4.*AA*AA)
      SQHH = DSQRT( DABS(BIGH) )
      V2 = Z2 + BB/(2.*AA)
      VX = ZMET + BB/(2.*AA)
      IF(ZL .LT. ZMET) VX = ZL + BB/(2.*AA)
      IF ( BIGH .LT. 0. ) GO TO 400
      ARG1 = DATAN ( VX/SQHH ) / SQHH
      ARG2 = DATAN ( V2/SQHH ) / SQHH
      ISOL = 1
      GO TO 500
400   ARG1 = DLOG ( DABS(VX-SQHH) / (VX+SQHH) ) / ( 2.*SQHH )
      ARG2 = DLOG ( DABS(V2-SQHH) / (V2+SQHH) ) / ( 2.*SQHH )
      ISOL = 2
500   ARG3 = 0.
      IF(ZL .LT. ZMET) ARG3 = DLOG ( (ZMET-D)/(ZL-D) )
      RASTUB = ( ZL-Z2 )/AA* ( ARG1 - ARG2 ) + ARG3
!
      TOPG3 = RASTUB - ARG3
      ZUP = ZMET
      IF ( ZL .LT. ZMET ) ZUP = ZL
      BOTG3 = DLOG( (ZUP-D)/(Z2-D) )
      G3 = TOPG3/BOTG3
!
!-----------------------------------------------------------------------
!
      CALL CALFAC ( Z2, ZC, Z1, CDG, COEFS, AT, BT, AL, BL, &
                    A0U, B0U, A0L, B0L, SIGMA, CS, PS, &
                    RBC, HA, RDC, CORB1, CORB2 )
!
      ARG1 = DLOG((ZL-D)/Z0)
      ARG2 = DLOG((ZL-D)/(Z2-D))
      U2FAC = ( ARG1 - G2*ARG2 ) / ( ARG1 - ARG2 )
 
      RETURN
      END

!=======================================================================
      SUBROUTINE DENCAL( Z2, Z1, ZC, ZLT, A0L, B0L, A0U, B0U)
!======================================================================
      IMPLICIT REAL*8 (A-H, O-Z)
!
      ZLMAX = 2.*( ZLT/( Z2-Z1 ) - 0.001 ) + 0.001
!
      B0L = (ZLMAX-0.001)/(ZC-Z1)
      A0L =  ZLMAX-B0L*ZC
      B0U = (ZLMAX-0.001)/(ZC-Z2)
      A0U =  ZLMAX-B0U*ZC
!
      RETURN
      END
!
!=======================================================================
      SUBROUTINE WINVAL( Z, A, B, ALPHA, BETA)
!======================================================================
      IMPLICIT REAL*8 (A-H, O-Z)
      REAL*8 A,B,ZNU,ARG
!
      ZNU = 1./DABS(B)**(2./3.) * (A + B*Z)
      ARG = 2./3. * ZNU**(3./2.)
      
      ORDER = 1./3.

      call bessik(ARG,ORDER,BI,BK,BIP,BKP)
      
! The XMIN parameter in bessik is 2.0d0
! If argument > 2 then use exponential scaling

      if (ARG.lt.(2.0d0)) then

      ALPHA = DSQRT(ZNU) * BI
      BETA  = DSQRT(ZNU) * BK
      
      else if (ARG.ge.(2.0d0)) then
      
      ALPHA = DSQRT(ZNU) * BI*DEXP(ARG)
      BETA  = DSQRT(ZNU) * BK*DEXP(-ARG)
      
      endif
      
      RETURN
      END
!
!=======================================================================
      SUBROUTINE GRADVA(Z, A, B, ALPHA, BETA)
!======================================================================
      IMPLICIT REAL*8 (A-H, O-Z)
!
      ZNU = 1./DABS(B)**(2./3.) * (A + B*Z )
      ARG = 2./3. * ZNU**(3./2.)
      
!      print*, 'ZNU in GRADVAL is', ZNU
!      print*, 'Z in gradval is', Z
!
      ORDER = 1./3.
      
      call bessik(ARG,ORDER,BI13,BK13,BIP,BKP)
      
      ORDER = 4./3.
      
      call bessik(ARG,ORDER,BI43,BK43,BIP,BKP)

      ORDER= 2./3.
      
      call bessik(ARG,ORDER,BI23,BK23,BIP,BKP) 
      
!      BIm23 = BI23 + (2./3.141592)*sin((2./3.)*3.141592)*BK23
      
      if (ARG.lt.(2.0d0)) then
      
      ALPHA = (1./DSQRT(ZNU))*BI13 + ZNU*BI43

!      ALPHA = ZNU * BIm23

      ALPHA = ALPHA/2. * (DABS(B)**(1./3.))*DSIGN(1.D0,B)

!      BETA = BK13 - ZNU*BK43 + 1./(2.*DSQRT(ZNU)) * BK13
!      BETA = - ZNU/2.*(BK43 + BK23) + 1./(2.*DSQRT(ZNU)) * BK13
      BETA = -ZNU*BK23
      
      BETA = BETA/2. *(DABS(B)**(1./3.))*DSIGN(1.D0,B)
      
      else if (ARG.ge.(2.0d0)) then
      
      ALPHA=(1./DSQRT(ZNU))*BI13*DEXP(ARG) + ZNU*BI43*DEXP(ARG)

!      ALPHA = ZNU * BIm23 * DEXP(ARG)     
      ALPHA = ALPHA/2. * (DABS(B)**(1./3.))*DSIGN(1.D0,B)

!      BETA = BK13 * DEXP(-ARG) - ZNU * BK43 * DEXP(-ARG)
!     & + 1./(2.*DSQRT(ZNU)) * BK13 * DEXP(-ARG)

!      BETA = - ZNU/2. * (BK43 + BK23) * DEXP(-ARG)
!     & + 1./(2.*DSQRT(ZNU)) * BK13 * DEXP(-ARG)

      BETA = -ZNU*BK23*DEXP(-ARG)
      
      BETA = BETA/2. *(DABS(B)**(1./3.))*DSIGN(1.D0,B)
      
      endif

      RETURN
      END
      
!
!=======================================================================
      SUBROUTINE DEVAL1(COEFS, AT, BT, Z2, G1, SIGMA, VKC, D1)
!=======================================================================
      IMPLICIT REAL*8 (A-H, O-Z)
      DIMENSION COEFS(4)
!
      CALL WINVAL(Z2, AT, BT, A1, B1)
      Y = COEFS(1)*A1 + COEFS(2)*B1
!
      CALL GRADVA(Z2,AT,BT,A2,B2)
      
      YDASH = COEFS(1)*A2 + COEFS(2)*B2

      D1 = Z2 - 1./(G1*VKC)*DSQRT( Y*SIGMA/YDASH )
!
      RETURN
      END
!
!=======================================================================
      SUBROUTINE DEVAL2(COEFS, AT, BT, AL, BL, A0U, B0U, A0L, B0L, &
                      ZS, Z2, ZC, Z1, G1, VKC, PS, CDD, SIGMA, D2 )
!=======================================================================
      IMPLICIT REAL*8 (A-H, O-Z)
      DIMENSION COEFS(4)
!
      TOP   = 0.
      BOT   = 0.
      A     = AT
      B     = BT
      ZUP   = Z2
      ZLO   = ZC
      A0    = A0U
      B0    = B0U
      ALPHA = COEFS(1)
      BETA  = COEFS(2)
!
      DO 1000 ILEVEL = 1,2
!
      CALL YINTNU (A, B, A0, B0,  ZUP, ZLO, ALPHA, BETA, ZN1, ZN2)
!
      TOP = TOP + ZN2
      BOT = BOT + ZN1
!
      A     = AL
      B     = BL
      ZUP   = ZC
      ZLO   = Z1
      A0 = A0L
      B0 = B0L
      ALPHA = COEFS(3)
      BETA  = COEFS(4)
!
 1000 CONTINUE
!
      CALL WINVAL(Z1, AL, BL, ALP, BET)
!
      VAL = COEFS(3)*ALP + COEFS(4)*BET
      TOR1 = ((0.41*VAL)/alog(Z1/ZS))**2.d0
      BOTTOR = TOR1*PS/CDD*SIGMA
!
      D2 = TOP/(BOT+BOTTOR)
!
      RETURN
      END
!
!=======================================================================
      SUBROUTINE YINTNU(A, B, A0, B0, ZUP, ZLO, ALPHA, BETA, ZN1, ZN2 )
!=======================================================================
      IMPLICIT REAL*8 (A-H, O-Z)
!
!    CALCULATE INTEGRALS OF U*U*LD AND U*U*LD*Z NUMERICALLY
!
      ZN1 = 0.
      ZN2 = 0.
!
      DZ = (ZUP-ZLO)/20.
!
      DO 1000 I =1,20
!
      ZIN = I*DZ - DZ*0.5 + ZLO
      
      ZLTZ = A0 + B0 * ZIN
!
      CALL WINVAL ( ZIN, A, B, A1, B1 )
      VAL = ALPHA * A1 + BETA * B1
!
      ZN1 = ZN1 + VAL * ZLTZ * DZ
      ZN2 = ZN2 + VAL * ZLTZ * ZIN * DZ
!
 1000 CONTINUE
!
      RETURN
      END
!
!=======================================================================
      SUBROUTINE CALFAC ( Z2, ZC, Z1, CDG, COEFS, AT, BT, AL, BL, &
                          A0U, B0U, A0L, B0L, SIGMA, CS, PS, &
                          RBC, HA, RDC, CORB1, CORB2 )
!=======================================================================
!
      IMPLICIT REAL*8 ( A-H, O-Z )
      DIMENSION COEFS(4), HEIGHT(40), UU(40), GRADU(40), RBADD(40)
!
      RBC = 0.
      A = AL
      B = BL
      ZUP = ZC
      ZLO = Z1
      A0 = A0L
      B0 = B0L
      ALPHA = COEFS(3)
      BETA = COEFS(4)
!
      DO 1000 ILEVEL = 1, 2
      DZ = ( ZUP - ZLO ) / 20.
      DO 2000 IZ = 1, 20
!
      INDEX = (ILEVEL-1)* 20 + IZ
      ZIN = ZLO + DZ* (IZ-0.5)
      HEIGHT(INDEX) = ZIN
!
      CALL WINVAL ( ZIN, A, B, A1, B1 )
      UU(INDEX) = ALPHA*A1 + BETA *B1
!
      CALL GRADVA ( ZIN, A, B, A2, B2 )
      GRADU(INDEX) = ( ALPHA*A2 + BETA*B2 ) / DSQRT( UU(INDEX) )
!
      ZLTZ = A0 + B0 * ZIN
      RBC = RBC + ZLTZ * DSQRT( DSQRT ( UU(INDEX ) ) ) * DZ
      RBADD(INDEX) = RBC
!
!
!      WIND PROFILE AND LOCAL GRADIENT OF WIND PROFILE MAY BE
!      OUTPUT BY ACTIVATING WRITE STATEMENT BELOW.
!
!
!      WRITE ( 6, 900 ) HEIGHT(INDEX), UU(INDEX), GRADU(INDEX)
!900   FORMAT ( 3X,'HEIGHT ',F6.3,' UU ', F7.3,' DUDZ ', F10.8 )
!901   FORMAT ( 3X,3F10.4)
!
!
2000  CONTINUE
!
      A = AT
      B = BT
      ZUP = Z2
      ZLO = ZC
      A0 = A0U
      B0 = B0U
      ALPHA = COEFS(1)
      BETA = COEFS(2)
!
1000  CONTINUE
!
      TARGET = RBC/2.
      RBC = 1./RBC * CS * PS
!
!     CALCULATION OF SOURCE HEIGHT , HA
!
      DO 3000 IZ = 1, 40
      IF ( TARGET .GT. RBADD(IZ) ) GO TO 100
      WEIGHT = 1.- ( RBADD(IZ) - TARGET ) / ( RBADD(IZ) - RBADD(IZ-1) )
      HA = HEIGHT(IZ-1) + ( HEIGHT(IZ) - HEIGHT(IZ-1) ) * WEIGHT
      GO TO 200
100   CONTINUE
3000  CONTINUE
200   CONTINUE
!
!     CALCULATION OF RDC
!
      CALL WINVAL ( Z1, AL, BL, A1, B1 )
      U1 = DSQRT( COEFS(3)*A1 + COEFS(4)*B1 )
      RDC = 1./ ( U1*CDG )
!
      DO 4000 IZ = 1, 40
      DZ = HEIGHT(2) - HEIGHT(1)
      IF ( IZ .GT. 21 ) DZ = HEIGHT(40) - HEIGHT(39)
      ZINK = 1./( SIGMA*DSQRT(UU(IZ)) ) * DZ
      IF ( HEIGHT(IZ+1) .GT. HA ) GO TO 300
      RDC = RDC + ZINK
!
4000  CONTINUE
!
300   RDC = RDC + WEIGHT*ZINK
      DZSAVE = DZ
!
      CORB2 = ZINK * ( 1.-WEIGHT )
      ISTART = IZ + 1
      DO 5000 IZ = ISTART, 39
      DZ = HEIGHT(2) - HEIGHT(1)
      IF ( IZ .GT. 21 ) DZ = HEIGHT(40) - HEIGHT(39)
      CORB2 = CORB2 + 1. / ( SIGMA * DSQRT(UU(IZ)) ) * DZ
5000  CONTINUE
      CORB2 = CORB2**2
!
      CORB1 = GRADU(ISTART-1)**2 *( 1.-WEIGHT)*DZSAVE
      DO 6000 IZ = ISTART, 39
      DZ = HEIGHT(2) - HEIGHT(1)
      IF ( IZ .GT. 21 ) DZ = HEIGHT(40) - HEIGHT(39)
      CORB1 = CORB1 + GRADU(IZ)**2 * DZ
6000  CONTINUE
      CORB1 = CORB1 / ( Z2 - HA )
      CORB1 = 9. * 9.81 / ( 1010. * 1.2 * CORB1 )
!
      RETURN
      END
!
!=======================================================================
      SUBROUTINE GAUSSD(A,N,NP1,X,WORK)
!=======================================================================
!
!     SOLVE A LINEAR SYSTEM BY GAUSSIAN ELIMINATION.  DEVELOPED BY
!     DR. CHIN-HOH MOENG.  A IS THE MATRIX OF COEFFICIENTS, WITH THE
!     VECTOR OF CONSTANTS APPENDED AS AN EXTRA COLUMN.  X IS THE VECTOR
!     CONTAINING THE RESULTS.  THE INPUT MATRIX IS NOT DESTROYED.
!
      IMPLICIT REAL*8 ( A- H , O - Z )
      DIMENSION A(4,5),WORK(4,5),X(4)
!
      DO 1000 I=1,N
      DO 1000 J=1,NP1
1000  WORK(I,J)=A(I,J)
!
      DO 20 I=2,N
      DO 20 J=I,N
!
      R=WORK(J,I-1)/WORK(I-1,I-1)
!
      DO 20 K=1,NP1
20    WORK(J,K)=WORK(J,K)-R*WORK(I-1,K)
!
      DO 30 I=2,N
      K=N-I+2
      R=WORK(K,NP1)/WORK(K,K)
!
      DO 30 J=I,N
      L=N-J+1
30    WORK(L,NP1)=WORK(L,NP1)-R*WORK(L,K)
!
      DO 40 I=1,N
40    X(I)=WORK(I,NP1)/WORK(I,I)
!
      RETURN
      END
      
      	subroutine bessik(x,xnu,ri,rk,rip,rkp)
	
	integer MAXIT
	DOUBLE PRECISION ri,rip,rk,rkp,x,xnu,XMIN
	DOUBLE PRECISION EPS,FPMIN,PI
	PARAMETER (EPS=1.e-16,FPMIN=1.e-30,MAXIT=10000,XMIN=2., &
        PI=3.141592653589793d0)
       
       INTEGER i,l,nl
       DOUBLE PRECISION a,a1,b,c,d,del,del1,delh,dels,e,f,fact, &
       fact2,ff,gam1,gam2,gammi,gampl,h,p,pimu,q,q1,q2, &
       qnew,ril,ril1,rimu,rip1,ripl,ritemp,rk1,rkmu,rkmup, &
       rktemp,s,sum,sum1,x2,xi,xi2,xmu,xmu2
      
       if (x.le.0..or.xnu.lt.0) pause 'bad arguments in bessik'
       
       nl=int(xnu+.5d0)
       xmu=xnu-nl
       xmu2=xmu*xmu
       xi=1.d0/x
       xi2=2.d0*xi
       h=xnu*xi
       if(h.lt.FPMIN)h=FPMIN
       b=xi2*xnu
       d=0.d0
       c=h
       do i=1,MAXIT
       	b=b+xi2
       	d=1.d0/(b+d)
       	c=b+1.d0/c
       	del=c*d
       	h=del*h
       	if(abs(del-1.d0).lt.EPS)goto 1
       enddo
       
       pause 'x too large in bessik; try asymptotic expansion'
       
1      continue

       ril=FPMIN
       ripl=h*ril
       ril1=ril
       rip1=ripl
       fact=xnu*xi
       
       do l=nl,1,-1
       	ritemp=fact*ril+ripl
       	fact=fact-xi
       	ripl=fact*ritemp+ril
       	ril=ritemp
       enddo
       
      
       f=ripl/ril
       
       
       if(x.lt.XMIN) then
       	x2=.5d0*x
       	pimu=PI*xmu
       	if(abs(pimu).lt.EPS) then
          fact=1.0d0
        else
          fact=pimu/sin(pimu)
        endif
        d=-log(x2)
        e=xmu*d
        if(abs(e).lt.EPS) then
           fact2=1.d0
        else
           fact2=sinh(e)/e
        endif
        call beschb(xmu,gam1,gam2,gampl,gammi)
        ff=fact*(gam1*cosh(e)+gam2*fact2*d)
        sum=ff
        e=exp(e)
        p=0.5d0*e/gampl
        q=0.5d0/(e*gammi)
        c=1.d0
        d=x2*x2
        sum1=p
        
        do i=1,MAXIT
          ff=(i*ff+p+q)/(i*i-xmu2)
          c=c*d/i
          p=p/(i-xmu)
          q=q/(i+xmu)
          del=c*ff
          sum=sum+del
          del1=c*(p-i*ff)
          sum1=sum1+del1
          if(abs(del).lt.abs(sum)*EPS) goto 2
        enddo
        
        pause 'bessk series failed to converge'
        
2       continue

        rkmu=sum
        rk1=sum1*xi2
        
       else
        
        b=2.d0*(1.d0+x)
        d=1.d0/b
        delh=d
        h=delh
        q1=0.d0
        q2=1.d0
        a1=.25d0-xmu2
        c=a1
        q=c
        a=-a1
        s=1.d0+q*delh
                
        do i=2,MAXIT
         a=a-2*(i-1)
         c=-a*c/i
         qnew=(q1-b*q2)/a
         q1=q2
         q2=qnew
         q=q+c*qnew
         b=b+2.d0
         d=1.d0/(b+a*d)
         delh=(b*d-1.d0)*delh
         h=h+delh
         dels=q*delh
         s=s+dels
         if(abs(dels/s).lt.EPS) goto 3
        enddo
        
        pause 'bessik: failure to converge in CF2'
3       continue

        h=a1*h
        
! Exponentially scaled for all four functions

        rkmu=sqrt(PI/(2.d0*x))/s
        rk1=rkmu*(xmu+x+.5d0-h)*xi
        
       endif
       
       rkmup=xmu*xi*rkmu-rk1
       rimu=xi/(f*rkmu-rkmup)
       ri=(rimu*ril1)/ril
       rip=(rimu*rip1)/ril

       do i=1,nl
         rktemp=(xmu+i)*xi2*rk1+rkmu
         rkmu=rk1
         rk1=rktemp
       enddo
       
       rk=rkmu
       rkp=xnu*xi*rkmu-rk1

       return
       
       END
       
!=======================================================================
       subroutine beschb(x,gam1,gam2,gampl,gammi)
!=======================================================================
       INTEGER NUSE1,NUSE2
       DOUBLE PRECISION gam1,gam2,gammi,gampl,x
       PARAMETER (NUSE1=7,NUSE2=8)
       
       DOUBLE PRECISION xx,c1(7),c2(8),chebev
       SAVE c1,c2
       DATA c1/ -1.142022680371168d0, 6.5165112670737d-3, &
       3.087090173086d-4, -3.4706269649d-6, 6.9437664d-9, &
       3.67795d-11, -1.356d-13/
      
       DATA c2/ 1.843740587300905d0, -7.68528408447867d-2, &
       1.2719271366546d-3, -4.9717367042d-6, -3.31261198d-8, &
       2.423096d-10, -1.702d-13, -1.49d-15 /
      
       xx=8.d0*x*x-1.d0
       gam1=chebev(-1.0d0,1.0d0,c1,NUSE1,xx)
       gam2=chebev(-1.0d0,1.0d0,c2,NUSE2,xx)
       gampl=gam2-x*gam1
       gammi=gam2+x*gam1
       
       return
       END 
       
       function chebev(a,b,c,m,x)
       INTEGER m
       DOUBLE PRECISION chebev,a,b,x,c(m)
       INTEGER j
       DOUBLE PRECISION d,dd,sv,y,y2
       
       if ((x-a)*(x-b).gt.0) pause 'x not in range in chebev'
       d=0.
       dd=0.
       y=(2.*x-a-b)/(b-a)
       y2=2.*y
       
       do j=m,2,-1
       	sv=d
       	d=y2*d-dd+c(j)
       	dd=sv
       enddo
       
       chebev = y2*d-dd+0.5*c(1)
       
       return 
       END
