!======================================================================
SUBROUTINE InterpSR(QN, SRN, QAngleA, SRAngleA, SRVGiven, SRHGiven, SRVInt, SRHInt)
!======================================================================
! This subroutine interpolates horizontal and vertical reflectivities between
! array values at fixed angles to obtain
! horizontal and vertical reflectivities at the desired angle.
!
! History:
!   9/26/20 Kevin Schaefer deleted tabs
!----------------------------------------------------------------------

IMPLICIT NONE
INTEGER, INTENT(IN) :: QN        ! Number of quadrature angles
INTEGER, INTENT(IN) :: SRN     ! Number of reflectivity angles
DOUBLE PRECISION, DIMENSION(QN), INTENT(IN) :: QAngleA     
!  Array of Quad angles (X1in)
DOUBLE PRECISION, DIMENSION(SRN), INTENT(IN) :: SRAngleA
!  Array of reflec angles (X2in)
DOUBLE PRECISION, DIMENSION(SRN), INTENT(IN) :: SRVGiven
!  Array of input vertical reflectivities (Y1in)
DOUBLE PRECISION, DIMENSION(SRN), INTENT(IN) :: SRHGiven   
!  Array of input horizontal reflectivities (Y2in)
DOUBLE PRECISION, DIMENSION(QN + 2), INTENT(OUT) :: SRVInt
!  Array of interpolated vert reflectivities (Y1out)
DOUBLE PRECISION, DIMENSION(QN + 2), INTENT(OUT) :: SRHInt    
!  Array of interpolated horiz reflectivities (Y2out)

! Position 1 of interpolated arrays already has SRGiven(1)
! Position SRN of interpolated arrays already has SRGiven(SRN)

INTEGER :: i
INTEGER :: lowP, middleP, highP
DOUBLE PRECISION :: currentA
DOUBLE PRECISION :: deltaA, Aup
DOUBLE PRECISION :: deltaSRV, SRVup
DOUBLE PRECISION :: deltaSRH, SRHup
REAL(8), parameter :: tol = 0.00001

DO i = 1, QN
!
! Bracket the angle
  currentA = QAngleA(i)
  lowP = 1
  highP = SRN

  DO
    middleP = (lowP + highP) / 2
    IF(currentA >= SRAngleA(middleP)) THEN
      lowP = middleP
    ELSE
      highP = middleP
    END IF

    IF (highP - lowP <= 1) EXIT
  END DO

  IF ( abs(currentA - SRAngleA(1)) < tol) THEN
         lowP = 1
  ELSE IF ( abs(currentA - SRAngleA(QN)) < tol) THEN
     lowP = SRN - 1
  END IF

  highP = lowP + 1

!  The angle is now bracketed between lowP and highP
! Interpolate the surface reflectivities
  deltaA = SRAngleA(highP) - SRAngleA(lowP)
  deltaSRV = SRVGiven(highP) - SRVGiven(lowP)
  deltaSRH = SRHGiven(highP) - SRHGiven(lowP)
  Aup = currentA - SRAngleA(lowP)
  SRVup = deltaSRV * (Aup / deltaA)
  SRHup = deltaSRH * (Aup / deltaA)
  SRVInt(i + 1) = SRVGiven(lowP) + SRVup
  SRHInt(i + 1) = SRHGiven(lowP) + SRHup

END DO

RETURN
END SUBROUTINE
