!=================================================================================
real(8) function abh2o (temperature, pressure, RHO, frequency, dabh2o_t, dabh2o_v, dabh2o_p)
!=================================================================================
! calculates water vapor absorption coefficient
! assumes vvw line shape
!
! history:
!   2/17/1991 AJG version 10.1
!   3/1/1997 rewritten for Pascal by MK
!   8/1/2003 adapted from version 10, pwr differentiated and tested by R. Hill
!   9/26/2020 Kevin Schaefer deleted tabs

!outputs:
! abh2o = absorption coeff. in the atmosphere due to water vapor (Nepers/km)
! derivative of abh2o with respect to temperature: dabh2o_t  (Nepers/km/K)
! derivative of abh2o with respect to water vapor density: dabh2o_v  (Nepers/km/(grams/cubic meter))
! derivative of abh2o with respect to pressure dabh2o_p    (Nepers/km/mb)

!inputs:
! temperature (K)
! pressure (mb)
! RHO = water vapor density (g/m3), named humid in o2abs
! frequency (GHz)
!
! line parameters from H. J. Liebe, Radio Science V.20(5), pp.1069-1089 (1985)
! updated in FREQUENZ V.41, pp. 31-36 (1987)
! includes 30 MPM H2O lines and appropriate continuum term for this line base
! 
implicit none
  real(8) temperature, pressure, RHO, frequency
  real(8) f1, f2, TI1, TIP3, W2
  real(8) PVAP, PDA, TI, TI3, WFAC, SUM1, WIDTH, S
  real(8) dABH2O_t, dABH2O_v, dABH2O_p, dPVAP_t , dPVAP_v , dPVAP_p, dPDA_t, dPDA_v, dPDA_p,&
         dTIP3_t, dTIP3_v, dTIP3_p, dWFAC_t, dWFAC_v, dWFAC_p, dSUM1_t, dSUM1_v, dSUM1_p,&
       dWIDTH_t, dWIDTH_v, dWIDTH_p, dW2_t, dW2_v, dW2_p, dS_t, dS_v, dS_p, dTI1_t 
  real(8) dTI_t, dTI3_t
  integer I 
!const
  real(8), dimension(30) :: S1 =                                                                  &
               (/ 0.0109d0,  0.00011d0, 0.00007d0, 0.23000d0,  0.00464d0, 0.1540d0, 0.0001d0,    &
                  1.1900d0,  0.00044d0, 0.00637d0, 0.09210d0,  0.01940d0, 1.0600d0, 0.0330d0,    &
                  0.1280d0,  0.02530d0, 0.00374d0, 0.00125d0, 51.00000d0, 0.5090d0, 0.0274d0,    &
                 25.0000d0,  0.00130d0, 0.01330d0, 0.00550d0,  0.00380d0, 0.0183d0, 0.8560d0,    &
                  0.9160d0, 13.80000d0/)

  real(8), dimension(30) :: B2 =                                                                  &
               (/ 2.143d0, 8.730d0, 8.347d0, 0.653d0, 6.156d0, 1.515d0, 9.802d0,                 &
                  1.018d0, 7.318d0, 5.015d0, 3.561d0, 5.015d0, 1.370d0, 3.561d0,                 &
                  2.342d0, 2.814d0, 6.693d0, 6.693d0, 0.114d0, 2.150d0, 7.767d0,                 &
                  0.336d0, 8.113d0, 7.989d0, 7.845d0, 8.360d0, 5.039d0, 1.369d0,                 &
                  1.842d0, 0.178d0/)

  real(8), dimension(30) :: W3 =                                                                 &
               (/ 2.784d-3, 2.760d-3, 2.700d-3, 3.164d-3, 2.140d-3, 2.970d-3, 2.650d-3,         &
                  3.036d-3, 1.900d-3, 1.370d-3, 1.640d-3, 1.440d-3, 2.380d-3, 1.820d-3,         &    
                  1.980d-3, 2.490d-3, 1.150d-3, 1.190d-3, 3.000d-3, 2.230d-3, 3.000d-3,         &
                  2.860d-3, 1.410d-3, 2.860d-3, 2.860d-3, 2.640d-3, 2.340d-3, 2.530d-3,         &
                  2.400d-3, 2.860d-3/)

  real(8), dimension(30) :: FL =                                                                                   &
              (/ 22.235080d0,  67.813960d0, 119.995940d0, 183.310117d0, 321.225644d0, 325.152919d0, 336.187000d0, &
                380.197372d0, 390.134508d0, 437.346667d0, 439.150812d0, 443.018295d0, 448.001075d0, 470.888947d0, &
                474.689127d0, 488.491133d0, 503.568532d0, 504.482692d0, 556.936002d0, 620.700807d0, 658.006500d0, &
                752.033227d0, 841.073593d0, 859.865000d0, 899.407000d0, 902.555000d0, 906.205524d0, 916.171582d0, &
                970.315022d0, 987.926764d0/)

  if (RHO <= 0) then
      abh2o = 0.0d0
     dabh2o_t = 0.0d0
     dabh2o_v = 0.0d0
     dabh2o_p = 0.0d0
           return
  end if
  PVAP = RHO * temperature / 217.0d0
  dPVAP_t = RHO  / 217.0d0
  dPVAP_v = temperature / 217.0d0
  dPVAP_p = 0.d0

  PDA = pressure - PVAP
  dPDA_t = - dPVAP_t
  dPDA_v = - dPVAP_v
  dPDA_p = 1.d0

  TI = 300.0d0 / temperature
  dTI_t = -TI/temperature

  TI3 = TI**3.5d0
  dTI3_t = 3.5d0*TI3/TI*dTI_t

  TI1 = 1.0d0 - TI
  dTI1_t = - dTI_t

  TIP3 = PVAP * TI3
  dTIP3_t = dPVAP_t * TI3 + PVAP * dTI3_t
  dTIP3_v = dPVAP_v * TI3
  dTIP3_p = dPVAP_p * TI3

   WFAC   =  PDA   * (TI** 0.8d0)               + 4.8d0 *  PVAP   * TI
  dWFAC_t = dPDA_t * (TI** 0.8d0)               + 4.8d0 * dPVAP_t * TI &
          +  PDA   * (0.8d0*TI**(-0.2d0)*dTI_t) + 4.8d0 *  PVAP   * dTI_t
  dWFAC_v = dPDA_v * (TI** 0.8d0) + 4.8d0 * dPVAP_v * TI
  dWFAC_p = dPDA_p * (TI** 0.8d0) + 4.8d0 * dPVAP_p * TI

!  SUM1 = (1.13d-8 * PDA * TI**3.0d0 + 3.57d-7 * PVAP * TI**10.5d0 ) * PVAP * frequency
   SUM1 =   frequency*( PVAP   * (1.13d-8 *  PDA   * TI**3           + 3.57d-7 *  PVAP   * TI**10.5d0 )  )
  dSUM1_t = frequency*(dPVAP_t * (1.13d-8 *  PDA   * TI**3           + 3.57d-7 *  PVAP   * TI**10.5d0 ) &
                     +  PVAP   * (1.13d-8 * dPDA_t * TI**3           + 3.57d-7 * dPVAP_t * TI**10.5d0 &
                        + 1.13d-8 *  PDA * 3.0d0*TI**2*dTI_t + 3.57d-7 *  PVAP * 10.5d0*TI**9.5d0*dTI_t  ) )
  dSUM1_v = frequency*(  dPVAP_v *  (1.13d-8 * PDA    * TI**3    + 3.57d-7 *  PVAP   * TI**10.5d0 ) &
                     +    PVAP   *  (1.13d-8 * dPDA_v * TI**3    + 3.57d-7 * dPVAP_v * TI**10.5d0 )  )
  dSUM1_p = frequency*(  dPVAP_p *  (1.13d-8 * PDA    * TI**3    + 3.57d-7 *  PVAP   * TI**10.5d0 ) &
                     +    PVAP   *  (1.13d-8 * dPDA_p * TI**3    + 3.57d-7 * dPVAP_p * TI**10.5d0 )  )



  do I = 1,  30 
     WIDTH   = W3(I) * WFAC 
    dWIDTH_t = W3(I) * dWFAC_t
    dWIDTH_v = W3(I) * dWFAC_v
    dWIDTH_p = W3(I) * dWFAC_p

     W2   = WIDTH * WIDTH
    dW2_t = 2.d0*WIDTH * dWIDTH_t
    dW2_v = 2.d0*WIDTH * dWIDTH_v
    dW2_p = 2.d0*WIDTH * dWIDTH_p

    f1 = (frequency - FL(I))
    f2 = (frequency + FL(I))
!   S = S1(I) * PVAP * TI3 * Exp(B2(I) * (1.0d0 - TI))
     S   = S1(I) * TIP3 * Exp(B2(I) * TI1)
   dS_t = S1(I)* Exp(B2(I) * TI1) *( dTIP3_t  + TIP3 *(B2(I) * dTI1_t) )
   dS_v = S1(I)* Exp(B2(I) * TI1) *( dTIP3_v   )
   dS_p = S1(I)* Exp(B2(I) * TI1) *( dTIP3_p   )  


!   SUM1 = SUM1 + S * (frequency /FL(I)) * (WIDTH /(f1*f1) + (WIDTH*WIDTH))  &
!        + WIDTH/((f2*f2)+ (WIDTH*WIDTH))
SUM1    =  SUM1   +  S   * ( frequency /FL(I)) * ( WIDTH   /(f1*f1 + W2)&
                                                    + WIDTH   /(f2*f2+ W2)  )
dSUM1_t = dSUM1_t + dS_t * ( frequency /FL(I)) *  (WIDTH   /(f1*f1 + W2) + WIDTH/(f2*f2+ W2)  ) &
                 +  S   * ( frequency /FL(I)) * (((f1*f1 + W2)*dWIDTH_t - WIDTH*dW2_t )/(f1*f1 + W2)**2 &
                                              + ((f2*f2 + W2)*dWIDTH_t - WIDTH*dW2_t )/(f2*f2 + W2)**2  )
dSUM1_v = dSUM1_v + dS_v * ( frequency /FL(I)) * (WIDTH /(f1*f1 + W2) + WIDTH/(f2*f2+ W2)  ) &
                 +  S   * ( frequency /FL(I)) * (((f1*f1 + W2)*dWIDTH_v - WIDTH*dW2_v )/(f1*f1 + W2)**2 &
                                             + ((f2*f2 + W2)*dWIDTH_v - WIDTH*dW2_v )/(f2*f2 + W2)**2  )
dSUM1_p = dSUM1_p + dS_p * ( frequency /FL(I)) * (WIDTH /(f1*f1 + W2) + WIDTH/(f2*f2+ W2)  ) &
                 +  S   * ( frequency /FL(I)) * (((f1*f1 + W2)*dWIDTH_p - WIDTH*dW2_p )/(f1*f1 + W2)**2 &
                                            + ((f2*f2 + W2)*dWIDTH_p - WIDTH*dW2_p )/(f2*f2 + W2)**2  )

  end do

  abh2o    = 0.0419d0 * frequency *  SUM1
  dabh2o_t = 0.0419d0 * frequency * dSUM1_t
  dabh2o_v = 0.0419d0 * frequency * dSUM1_v
  dabh2o_p = 0.0419d0 * frequency * dSUM1_p

end function abh2o
