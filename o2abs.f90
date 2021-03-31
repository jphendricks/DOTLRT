!=============================================================================
real(8) function o2abs(temperature, pressure, vapor_density, frequency, do2abs_t, do2abs_v,do2abs_p)
!=============================================================================
! Real function o2abs  calculates the oxygen (o2) resonant and nonresonant absorption.
!
! temperature (K)
! pressure (mb)
! vapor_density (g/m**3) (water vapor density - enters linewidth calculation
! due to greater broadening efficiency of H2O)
! frequency (GHz)
! o2abs = air absorption coefficient due to oxygen , in (Np/km)
! frequencies, strength, linewidths from:
!    Liebe, Radio Science v.20, p. 1069 (1985)
! except for n=1 from:
!    Setzer&Pickett, J. Chem. Phys. v.67, p.340 (1977)
! interference from SY3, with LAGM=.16
! reference for equations:
!    P.W. Rosenkranz, JQSRT v.39, p.287 (1988)
! line are arranged 1-,1+,3-,3+,etc.
! m.k.: frequency of 773.8387 o2 abs line was changed to 773.8397GHz
! as a result of comparison between JPL online
! catalog and Rosenkranz abs line centers
!----------------------------------------------------------------------------
! History:
!   10/28/1988 pwr created routine
!   3/22/1991 AJG allow for adjustment of:
!            1, O2 linewidth temperature exponent XX
!            2, 118 line strength
!            3, O2-N2 low temperature exponent
!            4, O2-N2 temperature exponent breakpoint
!   4/3/1991 AJG include non-overlapping SMMW lines
!   6/1/2003 Ron Richter converted from Pascal to Fortran
!   9/26/2020 Kevin Schaefer deleted unused variables, tabs, unused code
!----------------------------------------------------------------------------
  implicit none
  real(8) temperature, pressure, vapor_density, frequency
  real(8)  B,    PRESDA,    PRESWV,    TH,    DFNR,    SUM1,    DF
  real(8) dB_t, dPRESDA_t, dPRESWV_t, dTH_t, dDFNR_t, dSUM1_t, dDF_t
  real(8)       dPRESDA_v, dPRESWV_v,        dDFNR_v, dSUM1_v, dDF_v
  real(8)                                    dDFNR_p, dSUM1_p, dDF_p

  real(8)  SF1,    SF2,    DEN,    G,    BFAC,    Y,    STR, f1, f2, f1sq, f2sq
  real(8) dSF1_t, dSF2_t, dDEN_t, dG_t, dBFAC_t, dY_t, dSTR_t
  real(8) dSF1_v, dSF2_v, dDEN_v,                dY_v, dSTR_v 
  real(8) dSF1_p, dSF2_p, dDEN_p,                dY_p, dSTR_p   
  real(8)   denomSUM1, ddenomSUM1_t, ddenomSUM1_v, ddenomSUM1_p 
  real(8) do2abs_t, do2abs_v, do2abs_p
  real(8) wb300, XX, TN2, YY, AA, BB, CC, Fp04 ! these could be PARAMETERs

  real(8)  DF2,    AAA,    TH3,    THp2,    TH1,    TH16
  real(8) dDF2_t, dAAA_t, dTH3_t, dTHp2_t, dTH1_t, dTH16_t
  real(8) dDF2_v,         dTH3_v, dTHp2_v
  real(8) dDF2_p, dAAA_p, dTH3_p, dTHp2_p
  real(8) SUMAND,  dSUMAND_t, dSUMAND_v, dSUMAND_p 
   
  integer       K, K2
! const
  real(8), dimension(34) :: F  =                                                   &
             (/118.7503d0, 56.2648d0, 62.4863d0, 58.4466d0, 60.3061d0, 59.5910d0, &
                59.1642d0, 60.4348d0, 58.3239d0, 61.1506d0, 57.6125d0, 61.8002d0, &
                56.9682d0, 62.4112d0, 56.3634d0, 62.9980d0, 55.7838d0, 63.5685d0, &
                55.2214d0, 64.1278d0, 54.6712d0, 64.6789d0, 54.1300d0, 65.2241d0, &
                53.5957d0, 65.7648d0, 53.0669d0, 66.3021d0, 52.5424d0, 66.8368d0, &
                52.0214d0, 67.3696d0, 51.5034d0, 67.9009d0/)
! width in GHz/bar
  real(8), dimension(34) :: w300 =                                                 &
                  (/1.630d0, 1.646d0, 1.468d0, 1.449d0, 1.382d0, 1.360d0,         &
                    1.319d0, 1.297d0, 1.266d0, 1.248d0, 1.221d0, 1.207d0,         &
                    1.181d0, 1.171d0, 1.144d0, 1.139d0, 1.110d0, 1.108d0,         &
                    1.079d0, 1.078d0, 1.050d0, 1.050d0, 1.020d0, 1.020d0,         &
                    1.000d0, 1.000d0, 0.970d0, 0.970d0, 0.940d0, 0.940d0,         &
                    0.920d0, 0.920d0, 0.890d0, 0.890d0/)

  real(8), dimension(34) :: S300  =                                                &
                  (/0.2936d-14, 0.8079d-15, 0.2480d-14, 0.2228d-14,               &
                    0.3351d-14, 0.3292d-14, 0.3721d-14, 0.3891d-14,               &
                    0.3640d-14, 0.4005d-14, 0.3227d-14, 0.3715d-14,               &
                    0.2627d-14, 0.3156d-14, 0.1982d-14, 0.2477d-14,               &
                    0.1391d-14, 0.1808d-14, 0.9124d-15, 0.1230d-14,               &
                    0.5603d-15, 0.7842d-15, 0.3228d-15, 0.4689d-15,               &
                    0.1748d-15, 0.2632d-15, 0.8898d-16, 0.1389d-15,               &
                    0.4264d-16, 0.6899d-16, 0.1924d-16, 0.3229d-16,               &
                    0.8191d-17, 0.1423d-16/)
! interference in 1/bar
  real(8), dimension(34) :: U  =                                                   & 
                  (/-0.0247d0,    0.2881d0,    -0.4290d0,     0.6848d0,           &
                    -0.7170d0,    0.8266d0,    -0.6032d0,     0.5664d0,           &
                    -0.2635d0,    0.1731d0,    -0.2414d0,     0.1738d0,           &
                    -0.0556d0,   -0.0048d0,    -0.0596d0,     0.0134d0,           &
                    -0.0920d0,    0.0541d0,    -0.1151d0,     0.0814d0,           &
                    -0.0706d0,    0.0415d0,    -0.0314d0,     0.0069d0,           &
                    -0.0066d0,   -0.0143d0,     0.0252d0,    -0.0428d0,           &
                     0.0579d0,   -0.0726d0,     0.0883d0,    -0.1002d0,           &
                     0.1165d0,   -0.1255d0/)
  real(8), dimension(34) :: V  =                                                   &
                  (/ 0.0003d0,   -0.0069d0,     0.0238d0,    -0.0647d0,           &
                     0.0916d0,   -0.1413d0,     0.1858d0,    -0.2323d0,           &
                     0.2686d0,   -0.3039d0,     0.3536d0,    -0.3797d0,           &
                     0.4104d0,   -0.4277d0,     0.4750d0,    -0.4860d0,           &
                     0.5025d0,   -0.5079d0,     0.5514d0,    -0.5525d0,           &
                     0.5520d0,   -0.5520d0,     0.5520d0,    -0.5520d0,           &
                     0.5520d0,   -0.5520d0,     0.5520d0,    -0.5520d0,           &
                     0.5520d0,   -0.5520d0,     0.5520d0,    -0.5520d0,           &
                     0.5520d0,   -0.5520d0/)
! {data for SMMW lines (from Liebe)
! center frequency in GHz}
  real(8), dimension(6) :: F0 = (/368.4983d0, 424.7631d0, 487.2494d0, 715.3931d0, 773.8397d0, 834.1453d0/)
! A1 IN KHZ/MB}
  real(8), dimension(6) :: A1 = (/6.79d-6, 6.38d-5, 2.35d-5, 9.96d-6, 6.71d-5, 1.8d-5/)
! A2 (UNITLESS)}
  real(8), dimension(6) :: A2 = (/0.020d0, 0.011d0, 0.011d0, 0.089d0, 0.079d0, 0.079d0/)
! A3 IN GHZ/MB}
  real(8), dimension(6) :: A3 = (/1.92d-3, 1.916d-3,1.92d-3, 1.81d-3, 1.81d-3, 1.81d-3/)

  wb300 = 0.48d0
  XX    = 0.8d0
  TN2   = 210.0d0
  YY    = 0.8d0
  AA    = 1.0d0
  BB    = 0.0d0
  CC    = 0.0d0

   PRESWV   = vapor_density * temperature / 217.0d0
  dPRESWV_t = vapor_density/217.0d0
  dPRESWV_v = temperature/217.0d0

   PRESDA   = pressure - PRESWV
  dPRESDA_t = - dPRESWV_t
  dPRESDA_v = - dPRESWV_v

   TH   = 300.0d0 / temperature
  dTH_t = -TH/temperature

   B   =   TH**XX
  dB_t = XX*B/TH*dTH_t  ! =XX*B/TH*( -TH/temperature)

   AAA   = 1.d0
  dAAA_t = AA * (  (pressure/1013.0d0)**BB) * ( CC*TH**CC/TH*dTH_t ) ! note that dAAA_t=0.d0 because CC=0
  dAAA_p = AA *  (BB*(pressure/1013.0d0)**BB )/(pressure) * (TH**CC)

   TH1   = 1.0d0 - TH
  dTH1_t =       -dTH_t

   TH16   = 6.89526d-3 *  TH1
  dTH16_t = 6.89526d-3 * dTH1_t

   THp2   =  PRESDA  *(TH**0.2d0)                                     + 1.10d0 *   PRESWV   * TH
  dTHp2_t = dPRESDA_t*(TH**0.2d0) + PRESDA*(0.2d0*TH**0.2d0/TH*dTH_t) + 1.10d0 * (dPRESWV_t * TH + PRESWV * dTH_t )
  dTHp2_v = dPRESDA_v*(TH**0.2d0)                                     + 1.10d0 *  dPRESWV_v * TH
  dTHp2_p =           (TH**0.2d0) ! dPRESDA_p = 1. , dTH_p = 0                    dPRESWV_p = 0

   TH3   =  PRESDA   * (TH**3)
  dTH3_t = dPRESDA_t * (TH**3) +  PRESDA * (3.d0*TH**3/TH*dTH_t)
  dTH3_v = dPRESDA_v * (TH**3)
  dTH3_p =             (TH**3) ! dPRESDA_p = 1.
  
  Fp04 = 0.04191d0*frequency
  if (temperature < TN2) THEN
     G   = ((300.0d0/TN2)**XX) * ((TN2/temperature)**YY)
     dG_t = G*(-YY/temperature)
  else
     G   =  B
     dG_t = dB_t

  end if
 DEN =   0.001d0 * (0.22d0 *   PRESDA   * B                  + 0.78d0 *   PRESDA   * G                  + 1.1d0 *   PRESWV   * TH)
dDEN_t = 0.001d0 * (0.22d0 * (dPRESDA_t * B + PRESDA * dB_t) + 0.78d0 * (dPRESDA_t * G + PRESDA * dG_t) &
                                                             +1.1d0*(dPRESWV_t*TH+PRESWV*dTH_t))
dDEN_v = 0.001d0 * (0.22d0 * (dPRESDA_v * B )                + 0.78d0 * (dPRESDA_v * G )                &
                                                             + 1.1d0 * (dPRESWV_v * TH ) )
dDEN_p = 0.001d0 * (0.22d0 *              B                  + 0.78d0 * G )! dPRESDA_p = 1, ! dPRESWV_p = 0

 DFNR   = wb300 * DEN
dDFNR_t = wb300 * dDEN_t
dDFNR_v = wb300 * dDEN_v
dDFNR_p = wb300 * dDEN_p

  denomSUM1   =  TH    * (frequency*frequency + DFNR*DFNR)
 ddenomSUM1_t = dTH_t  * (frequency*frequency + DFNR*DFNR) + TH * (2.d0*DFNR*dDFNR_t)
 ddenomSUM1_v =                                            + TH * (2.d0*DFNR*dDFNR_v)
 ddenomSUM1_p =                                            + TH * (2.d0*DFNR*dDFNR_p)

 SUM1   = 1.6d-17 * (frequency*frequency) *   DFNR  / denomSUM1
dSUM1_t = 1.6d-17 * (frequency*frequency) *( dDFNR_t/ denomSUM1 &
                                            - DFNR * ddenomSUM1_t  /denomSUM1**2 )
dSUM1_v = 1.6d-17 * (frequency*frequency) *( dDFNR_v/ denomSUM1 &
                                            - DFNR * ddenomSUM1_v  /denomSUM1**2 )
dSUM1_p = 1.6d-17 * (frequency*frequency) *( dDFNR_p/ denomSUM1&
                                            - DFNR * ddenomSUM1_p  /denomSUM1**2 )
! outputs from here show that derivatives that contribute to SUM1, dSUM1, which are
!dDFNR, ddenomSUM1, dTH_t, dDEN, dPRESDA, dB_t, dG_t, dPRESWV, are all correct. 
 
  do K = 1, 34 
    K2 = K / 2
    if ( (K2*2) /= K ) then
       BFAC   = exp( K * (K+1) * TH16 )
       dBFAC_t = exp( K * (K+1) * TH16 )*(K * (K+1) * dTH16_t)
    end if

     DF   = w300(K) *  DEN
     dDF_t = w300(K) * dDEN_t
     dDF_v = w300(K) * dDEN_v
     dDF_p = w300(K) * dDEN_p

     DF2   =      DF *  DF
     dDF2_t = 2.d0*DF * dDF_t
     dDF2_v = 2.d0*DF * dDF_v
     dDF2_p = 2.d0*DF * dDF_p


     Y    = DEN   * (U(K) + V(K) * TH)
     dY_t = dDEN_t * (U(K) + V(K) * TH) + DEN * (V(K) * dTH_t)
    dY_v = dDEN_v * (U(K) + V(K) * TH)
    dY_p = dDEN_p * (U(K) + V(K) * TH)

     STR   = S300(K) *  BFAC !these 3 lines would be more efficient within an ELSE of the IF below
     dSTR_t = S300(K) * dBFAC_t
     dSTR_p = 0.d0

    if (K == 1) then
       STR   =  STR   * AAA
       dSTR_t = dSTR_t * AAA + STR * dAAA_t
       dSTR_p =                STR * dAAA_p
    end if

    f1 = frequency - F(K)
    f2 = frequency + F(K)
     SF1   =    (DF   + f1* Y  )/(f1*f1 + DF2)
    dSF1_t = ( (dDF_t + f1*dY_t)*(f1*f1 + DF2) - (DF + f1*Y)*dDF2_t ) / ( f1*f1 + DF2 )**2
    dSF1_v = ( (dDF_v + f1*dY_v)*(f1*f1 + DF2) - (DF + f1*Y)*dDF2_v ) / ( f1*f1 + DF2 )**2
    dSF1_p = ( (dDF_p + f1*dY_p)*(f1*f1 + DF2) - (DF + f1*Y)*dDF2_p ) / ( f1*f1 + DF2 )**2

     SF2   =    (DF   - f2* Y  )/(f2*f2 + DF2)
    dSF2_t = ( (dDF_t - f2*dY_t)*(f2*f2 + DF2) - (DF + f2*Y)*dDF2_t ) / ( (f2*f2 + DF2) )**2
    dSF2_v = ( (dDF_v - f2*dY_v)*(f2*f2 + DF2) - (DF + f2*Y)*dDF2_v ) / ( (f2*f2 + DF2) )**2
    dSF2_p = ( (dDF_p - f2*dY_p)*(f2*f2 + DF2) - (DF + f2*Y)*dDF2_p ) / ( (f2*f2 + DF2) )**2

  SUM1   =  SUM1   +    STR * (SF1 + SF2)                           * (frequency*frequency)/ ( F(K)*F(K))
 dSUM1_t = dSUM1_t + ( dSTR_t*(SF1 + SF2) + STR*(dSF1_t + dSF2_t) ) * (frequency*frequency)/ ( F(K)*F(K))
 dSUM1_v = dSUM1_v + (                      STR*(dSF1_v + dSF2_v) ) * (frequency*frequency)/ ( F(K)*F(K))
 dSUM1_p = dSUM1_p + ( dSTR_p*(SF1 + SF2) + STR*(dSF1_p + dSF2_p) ) * (frequency*frequency)/ ( F(K)*F(K))

  end do !ON K

! output from a RETURN here shows that the dSUM1 are correct to here 

  dSUM1_t = (0.5034d12/ 3.14159d0) *( dSUM1_t * TH3 + SUM1 * dTH3_t )
  dSUM1_v = (0.5034d12/ 3.14159d0) *( dSUM1_v * TH3 + SUM1 * dTH3_v  )
  dSUM1_p = (0.5034d12/ 3.14159d0) *( dSUM1_p * TH3 + SUM1 * dTH3_p )
   SUM1   = (0.5034d12/ 3.14159d0) *   SUM1   * TH3
! output from a RETURN here shows that the dSUM1 are correct to here, BEWARE:
! SUM1 = ...  statement MUST be after the derivatives are calculated

! ADD NON-OVERLAPPING SMMW O2 LINES}

  do K = 1,  6
!    DF  =   A3(K)*(PRESDA*(TH**0.2d0) + 1.10d0 * PRESWV * TH)
     DF    = A3(K)* THp2
    dDF_t  = A3(K)*dTHp2_t
    dDF_v  = A3(K)*dTHp2_v
    dDF_p  = A3(K)*dTHp2_p
! the above DF derivatives have been verified

     DF2   =      DF * DF
    dDF2_t = 2.d0*DF * dDF_t
    dDF2_v = 2.d0*DF * dDF_v
    dDF2_p = 2.d0*DF * dDF_p
     ! the above DF derivatives have been verified


     STR =    frequency/F0(K)*A1(K) *    DF  *TH3*EXP(A2(K)*TH1)
    dSTR_t =  frequency/F0(K)*A1(K) * ( dDF_t*TH3*EXP(A2(K)*TH1) + DF*dTH3_t*EXP(A2(K)*TH1) +DF*TH3*EXP(A2(K)*TH1)*A2(K)*dTH1_t )
        dSTR_v =  frequency/F0(K)*A1(K) * ( dDF_v*TH3*EXP(A2(K)*TH1) + DF*dTH3_v*EXP(A2(K)*TH1) )
     dSTR_p =  frequency/F0(K)*A1(K) * ( dDF_p*TH3*EXP(A2(K)*TH1) + DF*dTH3_p*EXP(A2(K)*TH1) )
! the above STR derivatives have been verified

    f1sq = (frequency - F0(K))**2
    f2sq = (frequency + F0(K))**2
     SF1   = 1.0d0   /(f1sq + DF2)
    dSF1_t = -dDF2_t /(f1sq + DF2)**2
    dSF1_v = -dDF2_v /(f1sq + DF2)**2
    dSF1_p = -dDF2_p /(f1sq + DF2)**2

     SF2   = 1.0d0   /(f2sq + DF2)
    dSF2_t = -dDF2_t /(f2sq + DF2)**2
    dSF2_v = -dDF2_v /(f2sq + DF2)**2
    dSF2_p = -dDF2_p /(f2sq + DF2)**2
! the above SF2 derivatives have been verified

      SUMAND   = Fp04 *   STR  *(SF1 + SF2)
     dSUMAND_t = Fp04 *( dSTR_t*(SF1 + SF2) + STR*(dSF1_t + dSF2_t) )
     dSUMAND_v = Fp04 *( dSTR_v*(SF1 + SF2) + STR*(dSF1_v + dSF2_v) )
     dSUMAND_p = Fp04 *( dSTR_p*(SF1 + SF2) + STR*(dSF1_p + dSF2_p) )
! the summand derivatives have been verified

     SUM1   =  SUM1   +  SUMAND
    dSUM1_t = dSUM1_t + dSUMAND_t
    dSUM1_v = dSUM1_v + dSUMAND_v
    dSUM1_p = dSUM1_p + dSUMAND_p
 
  end do 

   o2abs   =  SUM1
  do2abs_t = dSUM1_t
  do2abs_v = dSUM1_v
  do2abs_p = dSUM1_p

  return

end function o2abs
!
