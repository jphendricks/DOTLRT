!=======================================================================
      subroutine gaussq(kind, n, alpha, beta, kpts, endpts, b, t, w)
!=======================================================================
! 
!           this set of routines computes the nodes t(j) and weights
!        w(j) for gaussian-type quadrature rules with pre-assigned
!        nodes.  these are used when one wishes to approximate
!
!                  integral (from a to b)  f(x) w(x) dx
!
!                              n
!        by                   sum w  f(t )
!                             j=1  j    j
!
!        (note w(x) and w(j) have no connection with each other.)
!        here w(x) is one of six possible non-negative weight
!        functions (listed below), and f(x) is the
!        function to be integrated.  gaussian quadrature is particularly
!        useful on infinite intervals (with appropriate weight
!        functions), since then other techniques often fail.
!
!           associated with each weight function w(x) is a set of
!        orthogonal polynomials.  the nodes t(j) are just the zeroes
!        of the proper n-th degree polynomial.
!
!     input parameters (all real numbers are in double precision)
!
!        kind     an integer between 1 and 6 giving the type of
!                 quadrature rule:
!
!        kind = 1:  legendre quadrature, w(x) = 1 on (-1, 1)
!        kind = 2:  chebyshev quadrature of the first kind
!                  w(x) = 1/sqrt(1 - x*x) on (-1, +1)
!        kind = 3:  chebyshev quadrature of the second kind
!                   w(x) = sqrt(1 - x*x) on (-1, 1)
!        kind = 4:  hermite quadrature, w(x) = exp(-x*x) on
!                   (-infinity, +infinity)
!        kind = 5:  jacobi quadrature, w(x) = (1-x)**alpha * (1+x)**
!                   beta on (-1, 1), alpha, beta .gt. -1.
!                   note: kind=2 and 3 are a special case of this.
!        kind = 6:  generalized laguerre quadrature, w(x) = exp(-x)*
!                   x**alpha on (0, +infinity), alpha .gt. -1
!
!        n        the number of points used for the quadrature rule
!        alpha    real parameter used only for gauss-jacobi and gauss-
!                 laguerre quadrature (otherwise use 0.d0).
!        beta     real parameter used only for gauss-jacobi quadrature--
!                 (otherwise use 0.d0)
!        kpts     (integer) normally 0, unless the left or right end-
!                 point (or both) of the interval is required to be a
!                 node (this is called gauss-radau or gauss-lobatto
!                 quadrature).  then kpts is the number of fixed
!                 endpoints (1 or 2).
!        endpts   real array of length 2.  contains the values of
!                 any fixed endpoints, if kpts = 1 or 2.
!        b        real scratch array of length n
!
!     output parameters (both double precision arrays of length n)
!
!        t        will contain the desired nodes.
!        w        will contain the desired weights w(j).
!
!     underflow may sometimes occur, but is harmless.
!
!     references
!        1.  golub, g. h., and welsch, j. h., "calculation of gaussian
!            quadrature rules," mathematics of computation 23 (april,
!            1969), pp. 221-230.
!        2.  golub, g. h., "some modified matrix eigenvalue problems,"
!            siam review 15 (april, 1973), pp. 318-334 (section 7).
!        3.  stroud and secrest, gaussian quadrature formulas, prentice-
!            hall, englewood cliffs, n.j., 1966.
!
! History
!   1/20/1975 original version from stanford
!   12/21/1983 eric grosse changes:
!          imtql2 => gausq2
!          hex constant => d1mach (from core library)
!          compute pi using datan
!          removed accuracy claims, description of method
!          added single precision version
!   9/26/20 Kevin Schaefer deleted tabs
!--------------------------------------------------------------------------
!
      double precision b(n), t(n), w(n), endpts(2), muzero, t1, &
        gam, solve, dsqrt, alpha, beta
!
      call class (kind, n, alpha, beta, b, t, muzero)
!
!           the matrix of coefficients is assumed to be symmetric.
!           the array t contains the diagonal elements, the array
!           b the off-diagonal elements.
!           make appropriate changes in the lower right 2 by 2
!           submatrix.
!
      if (kpts.eq.0)  go to 100
      if (kpts.eq.2)  go to  50
!
!           if kpts=1, only t(n) must be changed
!
      t(n) = solve(endpts(1), n, t, b)*b(n-1)**2 + endpts(1)
      go to 100
!
!           if kpts=2, t(n) and b(n-1) must be recomputed
!
   50 gam = solve(endpts(1), n, t, b)
      t1 = ((endpts(1) - endpts(2))/(solve(endpts(2), n, t, b) - gam))
      b(n-1) = dsqrt(t1)
      t(n) = endpts(1) + gam*t1
!
!           note that the indices of the elements of b run from 1 to n-1
!           and thus the value of b(n) is arbitrary.
!           now compute the eigenvalues of the symmetric tridiagonal
!           matrix, which has been modified as necessary.
!           the method used is a ql-type method with origin shifting
!
  100 w(1) = 1.0d0
      do 105 i = 2, n
  105    w(i) = 0.0d0
!
      call gausq2 (n, t, b, w, ierr)
      do 110 i = 1, n
  110    w(i) = muzero * w(i) * w(i)
!
      return
      end

!**********************************************************************
!**********************************************************************
      double precision function solve(shift, n, a, b)
!
!       this procedure performs elimination to solve for the
!       n-th component of the solution delta to the equation
!
!             (jn - shift*identity) * delta  = en,
!
!       where en is the vector of all zeroes except for 1 in
!       the n-th position.
!
!       the matrix jn is symmetric tridiagonal, with diagonal
!       elements a(i), off-diagonal elements b(i).  this equation
!       must be solved to obtain the appropriate changes in the lower
!       2 by 2 submatrix of coefficients for orthogonal polynomials.
!
!
      double precision shift, a(n), b(n), alpha
!
      alpha = a(1) - shift
      nm1 = n - 1
      do 10 i = 2, nm1
   10    alpha = a(i) - shift - b(i-1)**2/alpha
      solve = 1.0d0/alpha
      return
      end

!**********************************************************************
!**********************************************************************
      subroutine class(kind, n, alpha, beta, b, a, muzero)
!
!           this procedure supplies the coefficients a(j), b(j) of the
!        recurrence relation
!
!             b p (x) = (x - a ) p   (x) - b   p   (x)
!              j j            j   j-1       j-1 j-2
!
!        for the various classical (normalized) orthogonal polynomials,
!        and the zero-th moment
!
!             muzero = integral w(x) dx
!
!        of the given polynomial's weight function w(x).  since the
!        polynomials are orthonormalized, the tridiagonal matrix is
!        guaranteed to be symmetric.
!
!        the input parameter alpha is used only for laguerre and
!        jacobi polynomials, and the parameter beta is used only for
!        jacobi polynomials.  the laguerre and jacobi polynomials
!        require the gamma function.
!
      double precision a(n), b(n), muzero, alpha, beta
      double precision abi, a2b2, dgamma, pi, dsqrt, ab
!
      pi = 4.0d0 * datan(1.0d0)
      nm1 = n - 1
      go to (10, 20, 30, 40, 50, 60), kind
!
!              kind = 1:  legendre polynomials p(x)
!              on (-1, +1), w(x) = 1.
!
   10 muzero = 2.0d0
      do 11 i = 1, nm1
         a(i) = 0.0d0
         abi = i
   11    b(i) = abi/dsqrt(4*abi*abi - 1.0d0)
      a(n) = 0.0d0
      return
!
!              kind = 2:  chebyshev polynomials of the first kind t(x)
!              on (-1, +1), w(x) = 1 / sqrt(1 - x*x)
!
  20 muzero = pi
      do 21 i = 1, nm1
         a(i) = 0.0d0
   21    b(i) = 0.5d0
      b(1) = dsqrt(0.5d0)
      a(n) = 0.0d0
      return
!
!              kind = 3:  chebyshev polynomials of the second kind u(x)
!              on (-1, +1), w(x) = sqrt(1 - x*x)
!
   30 muzero = pi/2.0d0
      do 31 i = 1, nm1
         a(i) = 0.0d0
   31    b(i) = 0.5d0
      a(n) = 0.0d0
      return
!
!              kind = 4:  hermite polynomials h(x) on (-infinity,
!              +infinity), w(x) = exp(-x**2)
!
   40 muzero = dsqrt(pi)
      do 41 i = 1, nm1
         a(i) = 0.0d0
   41    b(i) = dsqrt(i/2.0d0)
      a(n) = 0.0d0
      return
!
!              kind = 5:  jacobi polynomials p(alpha, beta)(x) on
!              (-1, +1), w(x) = (1-x)**alpha + (1+x)**beta, alpha and
!              beta greater than -1
!
   50 ab = alpha + beta
      abi = 2.0d0 + ab
      muzero = 2.0d0 ** (ab + 1.0d0) * dgamma(alpha + 1.0d0) * &
      dgamma(beta + 1.0d0) / dgamma(abi)
      a(1) = (beta - alpha)/abi
      b(1) = dsqrt(4.0d0*(1.0d0 + alpha)*(1.0d0 + beta)/ &
      ((abi + 1.0d0)* abi*abi))
      a2b2 = beta*beta - alpha*alpha
      do 51 i = 2, nm1
         abi = 2.0d0*i + ab
         a(i) = a2b2/((abi - 2.0d0)*abi)
   51    b(i) = dsqrt (4.0d0*i*(i + alpha)*(i + beta)*(i + ab)/ &
         ((abi*abi - 1)*abi*abi))
      abi = 2.0d0*n + ab
      a(n) = a2b2/((abi - 2.0d0)*abi)
      return
!
!              kind = 6:  laguerre polynomials l(alpha)(x) on
!              (0, +infinity), w(x) = exp(-x) * x**alpha, alpha greater
!              than -1.
!
   60 muzero = dgamma(alpha + 1.0d0)
      do 61 i = 1, nm1
         a(i) = 2.0d0*i - 1.0d0 + alpha
   61    b(i) = dsqrt(i*(i + alpha))
      a(n) = 2.0d0*n - 1 + alpha
      return
      end

!**********************************************************************
!**********************************************************************
      subroutine gausq2(n, d, e, z, ierr)
!
!     this subroutine is a translation of an algol procedure,
!     num. math. 12, 377-383(1968) by martin and wilkinson,
!     as modified in num. math. 15, 450(1970) by dubrulle.
!     handbook for auto. comp., vol.ii-linear algebra, 241-248(1971).
!     this is a modified version of the 'eispack' routine imtql2.
!
!     this subroutine finds the eigenvalues and first components of the
!     eigenvectors of a symmetric tridiagonal matrix by the implicit ql
!     method.
!
!     on input:
!
!        n is the order of the matrix;
!
!        d contains the diagonal elements of the input matrix;
!
!        e contains the subdiagonal elements of the input matrix
!          in its first n-1 positions.  e(n) is arbitrary;
!
!        z contains the first row of the identity matrix.
!
!      on output:
!
!        d contains the eigenvalues in ascending order.  if an
!          error exit is made, the eigenvalues are correct but
!          unordered for indices 1, 2, ..., ierr-1;
!
!        e has been destroyed;
!
!        z contains the first components of the orthonormal eigenvectors
!          of the symmetric tridiagonal matrix.  if an error exit is
!          made, z contains the eigenvectors associated with the stored
!          eigenvalues;
!
!        ierr is set to
!          zero       for normal return,
!          j          if the j-th eigenvalue has not been
!                     determined after 30 iterations.
!
!     ------------------------------------------------------------------
!
      integer i, j, k, l, m, n, ii, mml, ierr
      real*8 d(n), e(n), z(n), b, c, f, g, p, r, s, machep
      real*8 dsqrt, dabs, dsign, d1mach
!
      machep=d1mach(4)
!
      ierr = 0
      if (n .eq. 1) go to 1001
!
      e(n) = 0.0d0
      do 240 l = 1, n
         j = 0
!     :::::::::: look for small sub-diagonal element ::::::::::
  105    do 110 m = l, n
            if (m .eq. n) go to 120
            if (dabs(e(m)) .le. machep * (dabs(d(m)) + dabs(d(m+1)))) &
              go to 120
  110    continue
!
  120    p = d(l)
         if (m .eq. l) go to 240
         if (j .eq. 30) go to 1000
         j = j + 1
!     :::::::::: form shift ::::::::::
         g = (d(l+1) - p) / (2.0d0 * e(l))
         r = dsqrt(g*g+1.0d0)
         g = d(m) - p + e(l) / (g + dsign(r, g))
         s = 1.0d0
         c = 1.0d0
         p = 0.0d0
         mml = m - l
!
!     :::::::::: for i=m-1 step -1 until l do -- ::::::::::
         do 200 ii = 1, mml
            i = m - ii
            f = s * e(i)
            b = c * e(i)
            if (dabs(f) .lt. dabs(g)) go to 150
            c = g / f
            r = dsqrt(c*c+1.0d0)
            e(i+1) = f * r
            s = 1.0d0 / r
            c = c * s
            go to 160
  150       s = f / g
            r = dsqrt(s*s+1.0d0)
            e(i+1) = g * r
            c = 1.0d0 / r
            s = s * c
  160       g = d(i+1) - p
            r = (d(i) - g) * s + 2.0d0 * c * b
            p = s * r
            d(i+1) = g + p
            g = c * r - b
!     :::::::::: form first component of vector ::::::::::
            f = z(i+1)
            z(i+1) = s * z(i) + c * f
  200       z(i) = c * z(i) - s * f
!
         d(l) = d(l) - p
         e(l) = g
         e(m) = 0.0d0
         go to 105
  240 continue
!
!     :::::::::: order eigenvalues and eigenvectors ::::::::::
      do 300 ii = 2, n
         i = ii - 1
         k = i
         p = d(i)
!
         do 260 j = ii, n
            if (d(j) .ge. p) go to 260
            k = j
            p = d(j)
  260    continue
!
         if (k .eq. i) go to 300
         d(k) = d(i)
         d(i) = p
         p = z(i)
         z(i) = z(k)
         z(k) = p
  300 continue
!
      go to 1001
!     :::::::::: set error -- no convergence to an
!                eigenvalue after 30 iterations ::::::::::
 1000 ierr = l
 1001 return
!     :::::::::: last card of gausq2 ::::::::::
      end

!**********************************************************************
!**********************************************************************
      DOUBLE PRECISION FUNCTION D1MACH(I)
      INTEGER I
!
!  DOUBLE-PRECISION MACHINE CONSTANTS
!  D1MACH( 1) = B**(EMIN-1), THE SMALLEST POSITIVE MAGNITUDE.
!  D1MACH( 2) = B**EMAX*(1 - B**(-T)), THE LARGEST MAGNITUDE.
!  D1MACH( 3) = B**(-T), THE SMALLEST RELATIVE SPACING.
!  D1MACH( 4) = B**(1-T), THE LARGEST RELATIVE SPACING.
!  D1MACH( 5) = LOG10(B)
!
      INTEGER SMALL(2)
      INTEGER LARGE(2)
      INTEGER RIGHT(2)
      INTEGER DIVER(2)
      INTEGER LOG10(2)
      INTEGER SC, CRAY1(38), J
      COMMON /D9MACH/ CRAY1
      SAVE SMALL, LARGE, RIGHT, DIVER, LOG10, SC
      DOUBLE PRECISION DMACH(5)
      EQUIVALENCE (DMACH(1),SMALL(1))
      EQUIVALENCE (DMACH(2),LARGE(1))
      EQUIVALENCE (DMACH(3),RIGHT(1))
      EQUIVALENCE (DMACH(4),DIVER(1))
      EQUIVALENCE (DMACH(5),LOG10(1))
!  THIS VERSION ADAPTS AUTOMATICALLY TO MOST CURRENT MACHINES.
!  R1MACH CAN HANDLE AUTO-DOUBLE COMPILING, BUT THIS VERSION OF
!  D1MACH DOES NOT, BECAUSE WE DO NOT HAVE QUAD CONSTANTS FOR
!  MANY MACHINES YET.
!  TO COMPILE ON OLDER MACHINES, ADD A C IN COLUMN 1
!  ON THE NEXT LINE
      DATA SC/0/
!  AND REMOVE THE C FROM COLUMN 1 IN ONE OF THE SECTIONS BELOW.
!  CONSTANTS FOR EVEN OLDER MACHINES CAN BE OBTAINED BY
!          mail netlib@research.bell-labs.com
!          send old1mach from blas
!  PLEASE SEND CORRECTIONS TO dmg OR ehg@bell-labs.com.
!
!     MACHINE CONSTANTS FOR THE HONEYWELL DPS 8/70 SERIES.
!      DATA SMALL(1),SMALL(2) / O402400000000, O000000000000 /
!      DATA LARGE(1),LARGE(2) / O376777777777, O777777777777 /
!      DATA RIGHT(1),RIGHT(2) / O604400000000, O000000000000 /
!      DATA DIVER(1),DIVER(2) / O606400000000, O000000000000 /
!      DATA LOG10(1),LOG10(2) / O776464202324, O117571775714 /, SC/987/
!
!     MACHINE CONSTANTS FOR PDP-11 FORTRANS SUPPORTING
!     32-BIT INTEGERS.
!      DATA SMALL(1),SMALL(2) /    8388608,           0 /
!      DATA LARGE(1),LARGE(2) / 2147483647,          -1 /
!      DATA RIGHT(1),RIGHT(2) /  612368384,           0 /
!      DATA DIVER(1),DIVER(2) /  620756992,           0 /
!      DATA LOG10(1),LOG10(2) / 1067065498, -2063872008 /, SC/987/
!
!     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES.
!      DATA SMALL(1),SMALL(2) / O000040000000, O000000000000 /
!      DATA LARGE(1),LARGE(2) / O377777777777, O777777777777 /
!      DATA RIGHT(1),RIGHT(2) / O170540000000, O000000000000 /
!      DATA DIVER(1),DIVER(2) / O170640000000, O000000000000 /
!      DATA LOG10(1),LOG10(2) / O177746420232, O411757177572 /, SC/987/
!
!     ON FIRST CALL, IF NO DATA UNCOMMENTED, TEST MACHINE TYPES.
      IF (SC .NE. 987) THEN
         DMACH(1) = 1.D13
         IF (      SMALL(1) .EQ. 1117925532 &
            .AND. SMALL(2) .EQ. -448790528) THEN
!           *** IEEE BIG ENDIAN ***
            SMALL(1) = 1048576
            SMALL(2) = 0
            LARGE(1) = 2146435071
            LARGE(2) = -1
            RIGHT(1) = 1017118720
            RIGHT(2) = 0
            DIVER(1) = 1018167296
            DIVER(2) = 0
            LOG10(1) = 1070810131
            LOG10(2) = 1352628735
         ELSE IF ( SMALL(2) .EQ. 1117925532 &
            .AND. SMALL(1) .EQ. -448790528) THEN
!           *** IEEE LITTLE ENDIAN ***
            SMALL(2) = 1048576
            SMALL(1) = 0
            LARGE(2) = 2146435071
            LARGE(1) = -1
            RIGHT(2) = 1017118720
            RIGHT(1) = 0
            DIVER(2) = 1018167296
            DIVER(1) = 0
            LOG10(2) = 1070810131
            LOG10(1) = 1352628735
         ELSE IF ( SMALL(1) .EQ. -2065213935 &
            .AND. SMALL(2) .EQ. 10752) THEN
!               *** VAX WITH D_FLOATING ***
            SMALL(1) = 128
            SMALL(2) = 0
            LARGE(1) = -32769
            LARGE(2) = -1
            RIGHT(1) = 9344
            RIGHT(2) = 0
            DIVER(1) = 9472
            DIVER(2) = 0
            LOG10(1) = 546979738
            LOG10(2) = -805796613
         ELSE IF ( SMALL(1) .EQ. 1267827943 &
            .AND. SMALL(2) .EQ. 704643072) THEN
!               *** IBM MAINFRAME ***
            SMALL(1) = 1048576
            SMALL(2) = 0
            LARGE(1) = 2147483647
            LARGE(2) = -1
            RIGHT(1) = 856686592
            RIGHT(2) = 0
            DIVER(1) = 873463808
            DIVER(2) = 0
            LOG10(1) = 1091781651
            LOG10(2) = 1352628735
         ELSE IF ( SMALL(1) .EQ. 1120022684 &
            .AND. SMALL(2) .EQ. -448790528) THEN
!           *** CONVEX C-1 ***
            SMALL(1) = 1048576
            SMALL(2) = 0
            LARGE(1) = 2147483647
            LARGE(2) = -1
            RIGHT(1) = 1019215872
            RIGHT(2) = 0
            DIVER(1) = 1020264448
            DIVER(2) = 0
            LOG10(1) = 1072907283
            LOG10(2) = 1352628735
         ELSE IF ( SMALL(1) .EQ. 815547074 &
            .AND. SMALL(2) .EQ. 58688) THEN
!           *** VAX G-FLOATING ***
            SMALL(1) = 16
            SMALL(2) = 0
            LARGE(1) = -32769
            LARGE(2) = -1
            RIGHT(1) = 15552
            RIGHT(2) = 0
            DIVER(1) = 15568
            DIVER(2) = 0
            LOG10(1) = 1142112243
            LOG10(2) = 2046775455
         ELSE
            DMACH(2) = 1.D27 + 1
            DMACH(3) = 1.D27
            LARGE(2) = LARGE(2) - RIGHT(2)
            IF (LARGE(2) .EQ. 64 .AND. SMALL(2) .EQ. 0) THEN
               CRAY1(1) = 67291416
               DO 10 J = 1, 20
                  CRAY1(J+1) = CRAY1(J) + CRAY1(J)
 10               CONTINUE
               CRAY1(22) = CRAY1(21) + 321322
               DO 20 J = 22, 37
                  CRAY1(J+1) = CRAY1(J) + CRAY1(J)
 20               CONTINUE
               IF (CRAY1(38) .EQ. SMALL(1)) THEN
!                  *** CRAY ***
                  CALL I1MCRY(SMALL(1), J, 8285, 8388608, 0)
                  SMALL(2) = 0
                  CALL I1MCRY(LARGE(1), J, 24574, 16777215, 16777215)
                  CALL I1MCRY(LARGE(2), J, 0, 16777215, 16777214)
                  CALL I1MCRY(RIGHT(1), J, 16291, 8388608, 0)
                  RIGHT(2) = 0
                  CALL I1MCRY(DIVER(1), J, 16292, 8388608, 0)
                  DIVER(2) = 0
                  CALL I1MCRY(LOG10(1), J, 16383, 10100890, 8715215)
                  CALL I1MCRY(LOG10(2), J, 0, 16226447, 9001388)
               ELSE
                  WRITE(*,9000)
                  STOP 779
                  END IF
            ELSE
               WRITE(*,9000)
               STOP 779
               END IF
            END IF
         SC = 987
         END IF
!    SANITY CHECK
      IF (DMACH(4) .GE. 1.0D0) STOP 778
      IF (I .LT. 1 .OR. I .GT. 5) THEN
         WRITE(*,*) 'D1MACH(I): I =',I,' is out of bounds.'
         STOP
         END IF
      D1MACH = DMACH(I)
      RETURN
 9000 FORMAT(/' Adjust D1MACH by uncommenting data statements'/ &
      ' appropriate for your machine.')
! /* Standard C source for D1MACH -- remove the * in column 1 */
!#include <stdio.h>
!#include <float.h>
!#include <math.h>
!double d1mach_(long *i)
!{
!       switch(*i){
!         case 1: return DBL_MIN;
!         case 2: return DBL_MAX;
!         case 3: return DBL_EPSILON/FLT_RADIX;
!         case 4: return DBL_EPSILON;
!         case 5: return log10((double)FLT_RADIX);
!         }
!       fprintf(stderr, "invalid argument: d1mach(%ld)\n", *i);
!       exit(1); return 0; /* some compilers demand return values */
!}
      END

!**********************************************************************
!**********************************************************************
      SUBROUTINE I1MCRY(A, A1, B, C, D)
!***  SPECIAL COMPUTATION FOR OLD CRAY MACHINES ****
      INTEGER A, A1, B, C, D
      A1 = 16777216*B + C
      A = 16777216*A1 + D
      END
      INTEGER FUNCTION I1MACH(I)
      INTEGER I
!
!    I1MACH( 1) = THE STANDARD INPUT UNIT.
!    I1MACH( 2) = THE STANDARD OUTPUT UNIT.
!    I1MACH( 3) = THE STANDARD PUNCH UNIT.
!    I1MACH( 4) = THE STANDARD ERROR MESSAGE UNIT.
!    I1MACH( 5) = THE NUMBER OF BITS PER INTEGER STORAGE UNIT.
!    I1MACH( 6) = THE NUMBER OF CHARACTERS PER CHARACTER STORAGE UNIT.
!    INTEGERS HAVE FORM SIGN ( X(S-1)*A**(S-1) + ... + X(1)*A + X(0) )
!    I1MACH( 7) = A, THE BASE.
!    I1MACH( 8) = S, THE NUMBER OF BASE-A DIGITS.
!    I1MACH( 9) = A**S - 1, THE LARGEST MAGNITUDE.
!    FLOATS HAVE FORM  SIGN (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
!               WHERE  EMIN .LE. E .LE. EMAX.
!    I1MACH(10) = B, THE BASE.
!  SINGLE-PRECISION
!    I1MACH(11) = T, THE NUMBER OF BASE-B DIGITS.
!    I1MACH(12) = EMIN, THE SMALLEST EXPONENT E.
!    I1MACH(13) = EMAX, THE LARGEST EXPONENT E.
!  DOUBLE-PRECISION
!    I1MACH(14) = T, THE NUMBER OF BASE-B DIGITS.
!    I1MACH(15) = EMIN, THE SMALLEST EXPONENT E.
!    I1MACH(16) = EMAX, THE LARGEST EXPONENT E.
!
      INTEGER IMACH(16), OUTPUT, SC, SMALL(2)
      SAVE IMACH, SC
      REAL RMACH
      EQUIVALENCE (IMACH(4),OUTPUT), (RMACH,SMALL(1))
      INTEGER I3, J, K, T3E(3)
      DATA T3E(1) / 9777664 /
      DATA T3E(2) / 5323660 /
      DATA T3E(3) / 46980 /
!  THIS VERSION ADAPTS AUTOMATICALLY TO MOST CURRENT MACHINES,
!  INCLUDING AUTO-DOUBLE COMPILERS.
!  TO COMPILE ON OLDER MACHINES, ADD A C IN COLUMN 1
!  ON THE NEXT LINE
      DATA SC/0/
!  AND REMOVE THE C FROM COLUMN 1 IN ONE OF THE SECTIONS BELOW.
!  CONSTANTS FOR EVEN OLDER MACHINES CAN BE OBTAINED BY
!          mail netlib@research.bell-labs.com
!          send old1mach from blas
!  PLEASE SEND CORRECTIONS TO dmg OR ehg@bell-labs.com.
!
!     MACHINE CONSTANTS FOR THE HONEYWELL DPS 8/70 SERIES.
!
!      DATA IMACH( 1) /    5 /
!      DATA IMACH( 2) /    6 /
!      DATA IMACH( 3) /   43 /
!      DATA IMACH( 4) /    6 /
!      DATA IMACH( 5) /   36 /
!      DATA IMACH( 6) /    4 /
!      DATA IMACH( 7) /    2 /
!      DATA IMACH( 8) /   35 /
!      DATA IMACH( 9) / O377777777777 /
!      DATA IMACH(10) /    2 /
!      DATA IMACH(11) /   27 /
!      DATA IMACH(12) / -127 /
!      DATA IMACH(13) /  127 /
!      DATA IMACH(14) /   63 /
!      DATA IMACH(15) / -127 /
!      DATA IMACH(16) /  127 /, SC/987/
!
!     MACHINE CONSTANTS FOR PDP-11 FORTRANS SUPPORTING
!     32-BIT INTEGER ARITHMETIC.
!
!      DATA IMACH( 1) /    5 /
!      DATA IMACH( 2) /    6 /
!      DATA IMACH( 3) /    7 /
!      DATA IMACH( 4) /    6 /
!      DATA IMACH( 5) /   32 /
!      DATA IMACH( 6) /    4 /
!      DATA IMACH( 7) /    2 /
!      DATA IMACH( 8) /   31 /
!      DATA IMACH( 9) / 2147483647 /
!      DATA IMACH(10) /    2 /
!      DATA IMACH(11) /   24 /
!      DATA IMACH(12) / -127 /
!      DATA IMACH(13) /  127 /
!      DATA IMACH(14) /   56 /
!      DATA IMACH(15) / -127 /
!      DATA IMACH(16) /  127 /, SC/987/
!
!     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES.
!
!     NOTE THAT THE PUNCH UNIT, I1MACH(3), HAS BEEN SET TO 7
!     WHICH IS APPROPRIATE FOR THE UNIVAC-FOR SYSTEM.
!     IF YOU HAVE THE UNIVAC-FTN SYSTEM, SET IT TO 1.
!
!      DATA IMACH( 1) /    5 /
!      DATA IMACH( 2) /    6 /
!      DATA IMACH( 3) /    7 /
!      DATA IMACH( 4) /    6 /
!      DATA IMACH( 5) /   36 /
!      DATA IMACH( 6) /    6 /
!      DATA IMACH( 7) /    2 /
!      DATA IMACH( 8) /   35 /
!      DATA IMACH( 9) / O377777777777 /
!      DATA IMACH(10) /    2 /
!      DATA IMACH(11) /   27 /
!      DATA IMACH(12) / -128 /
!      DATA IMACH(13) /  127 /
!      DATA IMACH(14) /   60 /
!      DATA IMACH(15) /-1024 /
!      DATA IMACH(16) / 1023 /, SC/987/
!
      IF (SC .NE. 987) THEN
!        *** CHECK FOR AUTODOUBLE ***
         SMALL(2) = 0
         RMACH = 1E13
         IF (SMALL(2) .NE. 0) THEN
!           *** AUTODOUBLED ***
            IF (      (SMALL(1) .EQ. 1117925532 &
                .AND. SMALL(2) .EQ. -448790528) &
            .OR.     (SMALL(2) .EQ. 1117925532 &
                .AND. SMALL(1) .EQ. -448790528)) THEN
!               *** IEEE ***
               IMACH(10) = 2
               IMACH(14) = 53
               IMACH(15) = -1021
               IMACH(16) = 1024
            ELSE IF ( SMALL(1) .EQ. -2065213935 &
               .AND. SMALL(2) .EQ. 10752) THEN
!               *** VAX WITH D_FLOATING ***
               IMACH(10) = 2
               IMACH(14) = 56
               IMACH(15) = -127
               IMACH(16) = 127
            ELSE IF ( SMALL(1) .EQ. 1267827943 &
               .AND. SMALL(2) .EQ. 704643072) THEN
!               *** IBM MAINFRAME ***
               IMACH(10) = 16
               IMACH(14) = 14
               IMACH(15) = -64
               IMACH(16) = 63
            ELSE
               WRITE(*,9010)
               STOP 777
               END IF
            IMACH(11) = IMACH(14)
            IMACH(12) = IMACH(15)
            IMACH(13) = IMACH(16)
         ELSE
            RMACH = 1234567.
            IF (SMALL(1) .EQ. 1234613304) THEN 
!               *** IEEE ***
               IMACH(10) = 2
               IMACH(11) = 24
               IMACH(12) = -125
               IMACH(13) = 128
               IMACH(14) = 53
               IMACH(15) = -1021
               IMACH(16) = 1024
               SC = 987
            ELSE IF (SMALL(1) .EQ. -1271379306) THEN
!               *** VAX ***
               IMACH(10) = 2
               IMACH(11) = 24
               IMACH(12) = -127
               IMACH(13) = 127
               IMACH(14) = 56
               IMACH(15) = -127
               IMACH(16) = 127
               SC = 987
            ELSE IF (SMALL(1) .EQ. 1175639687) THEN
!               *** IBM MAINFRAME ***
               IMACH(10) = 16
               IMACH(11) = 6
               IMACH(12) = -64
               IMACH(13) = 63
               IMACH(14) = 14
               IMACH(15) = -64
               IMACH(16) = 63
               SC = 987
            ELSE IF (SMALL(1) .EQ. 1251390520) THEN
!              *** CONVEX C-1 ***
               IMACH(10) = 2
               IMACH(11) = 24
               IMACH(12) = -128
               IMACH(13) = 127
               IMACH(14) = 53
               IMACH(15) = -1024
               IMACH(16) = 1023
            ELSE
               DO 10 I3 = 1, 3
                  J = SMALL(1) / 10000000
                  K = SMALL(1) - 10000000*J
                  IF (K .NE. T3E(I3)) GO TO 20
                  SMALL(1) = J
 10               CONTINUE
!              *** CRAY T3E ***
               IMACH( 1) = 5
               IMACH( 2) = 6
               IMACH( 3) = 0
               IMACH( 4) = 0
               IMACH( 5) = 64
               IMACH( 6) = 8
               IMACH( 7) = 2
               IMACH( 8) = 63
               CALL I1MCR1(IMACH(9), K, 32767, 16777215, 16777215)
               IMACH(10) = 2
               IMACH(11) = 53
               IMACH(12) = -1021
               IMACH(13) = 1024
               IMACH(14) = 53
               IMACH(15) = -1021
               IMACH(16) = 1024
               GO TO 35
 20            CALL I1MCR1(J, K, 16405, 9876536, 0)
               IF (SMALL(1) .NE. J) THEN
                  WRITE(*,9020)
                  STOP 777
                  END IF
!              *** CRAY 1, XMP, 2, AND 3 ***
               IMACH(1) = 5
               IMACH(2) = 6
               IMACH(3) = 102
               IMACH(4) = 6
               IMACH(5) = 46
               IMACH(6) = 8
               IMACH(7) = 2
               IMACH(8) = 45
               CALL I1MCR1(IMACH(9), K, 0, 4194303, 16777215)
               IMACH(10) = 2
               IMACH(11) = 47
               IMACH(12) = -8188
               IMACH(13) = 8189
               IMACH(14) = 94
               IMACH(15) = -8141
               IMACH(16) = 8189
               GO TO 35
               END IF
            END IF
         IMACH( 1) = 5
         IMACH( 2) = 6
         IMACH( 3) = 7
         IMACH( 4) = 6
         IMACH( 5) = 32
         IMACH( 6) = 4
         IMACH( 7) = 2
         IMACH( 8) = 31
         IMACH( 9) = 2147483647
 35      SC = 987
         END IF
 9010 FORMAT(/' Adjust autodoubled I1MACH by uncommenting data'/ &
       ' statements appropriate for your machine and setting'/ &
       ' IMACH(I) = IMACH(I+3) for I = 11, 12, and 13.')
 9020 FORMAT(/' Adjust I1MACH by uncommenting data statements'/ &
       ' appropriate for your machine.')
      IF (I .LT. 1  .OR.  I .GT. 16) GO TO 40
      I1MACH = IMACH(I)
      RETURN
 40   WRITE(*,*) 'I1MACH(I): I =',I,' is out of bounds.'
      STOP
! /* C source for I1MACH -- remove the * in column 1 */
! /* Note that some values may need changing. */
!#include <stdio.h>
!#include <float.h>
!#include <limits.h>
!#include <math.h>
!
!long i1mach_(long *i)
!{
!       switch(*i){
!         case 1:  return 5;    /* standard input */
!         case 2:  return 6;    /* standard output */
!         case 3:  return 7;    /* standard punch */
!         case 4:  return 0;    /* standard error */
!         case 5:  return 32;   /* bits per integer */
!         case 6:  return sizeof(int);
!         case 7:  return 2;    /* base for integers */
!         case 8:  return 31;   /* digits of integer base */
!         case 9:  return LONG_MAX;
!         case 10: return FLT_RADIX;
!         case 11: return FLT_MANT_DIG;
!         case 12: return FLT_MIN_EXP;
!         case 13: return FLT_MAX_EXP;
!         case 14: return DBL_MANT_DIG;
!         case 15: return DBL_MIN_EXP;
!         case 16: return DBL_MAX_EXP;
!         }
!       fprintf(stderr, "invalid argument: i1mach(%ld)\n", *i);
!       exit(1);return 0; /* some compilers demand return values */
!}
      END

!**********************************************************************
!**********************************************************************
      SUBROUTINE I1MCR1(A, A1, B, C, D)
!***  SPECIAL COMPUTATION FOR OLD CRAY MACHINES ****
      INTEGER A, A1, B, C, D
      A1 = 16777216*B + C
      A = 16777216*A1 + D
      END
      REAL FUNCTION R1MACH(I)
      INTEGER I
!
!  SINGLE-PRECISION MACHINE CONSTANTS
!  R1MACH(1) = B**(EMIN-1), THE SMALLEST POSITIVE MAGNITUDE.
!  R1MACH(2) = B**EMAX*(1 - B**(-T)), THE LARGEST MAGNITUDE.
!  R1MACH(3) = B**(-T), THE SMALLEST RELATIVE SPACING.
!  R1MACH(4) = B**(1-T), THE LARGEST RELATIVE SPACING.
!  R1MACH(5) = LOG10(B)
!
      INTEGER SMALL(2)
      INTEGER LARGE(2)
      INTEGER RIGHT(2)
      INTEGER DIVER(2)
      INTEGER LOG10(2)
!     needs to be (2) for AUTODOUBLE, HARRIS SLASH 6, ...
      INTEGER SC
      SAVE SMALL, LARGE, RIGHT, DIVER, LOG10, SC
      REAL RMACH(5)
      EQUIVALENCE (RMACH(1),SMALL(1))
      EQUIVALENCE (RMACH(2),LARGE(1))
      EQUIVALENCE (RMACH(3),RIGHT(1))
      EQUIVALENCE (RMACH(4),DIVER(1))
      EQUIVALENCE (RMACH(5),LOG10(1))
      INTEGER J, K, L, T3E(3)
      DATA T3E(1) / 9777664 /
      DATA T3E(2) / 5323660 /
      DATA T3E(3) / 46980 /
!  THIS VERSION ADAPTS AUTOMATICALLY TO MOST CURRENT MACHINES,
!  INCLUDING AUTO-DOUBLE COMPILERS.
!  TO COMPILE ON OLDER MACHINES, ADD A C IN COLUMN 1
!  ON THE NEXT LINE
      DATA SC/0/
!  AND REMOVE THE C FROM COLUMN 1 IN ONE OF THE SECTIONS BELOW.
!  CONSTANTS FOR EVEN OLDER MACHINES CAN BE OBTAINED BY
!          mail netlib@research.bell-labs.com
!          send old1mach from blas
!  PLEASE SEND CORRECTIONS TO dmg OR ehg@bell-labs.com.
!
!     MACHINE CONSTANTS FOR THE HONEYWELL DPS 8/70 SERIES.
!      DATA RMACH(1) / O402400000000 /
!      DATA RMACH(2) / O376777777777 /
!      DATA RMACH(3) / O714400000000 /
!      DATA RMACH(4) / O716400000000 /
!      DATA RMACH(5) / O776464202324 /, SC/987/
!
!     MACHINE CONSTANTS FOR PDP-11 FORTRANS SUPPORTING
!     32-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL).
!      DATA SMALL(1) /    8388608 /
!      DATA LARGE(1) / 2147483647 /
!      DATA RIGHT(1) /  880803840 /
!      DATA DIVER(1) /  889192448 /
!      DATA LOG10(1) / 1067065499 /, SC/987/
!      DATA RMACH(1) / O00040000000 /
!      DATA RMACH(2) / O17777777777 /
!      DATA RMACH(3) / O06440000000 /
!      DATA RMACH(4) / O06500000000 /
!      DATA RMACH(5) / O07746420233 /, SC/987/
!
!     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES.
!      DATA RMACH(1) / O000400000000 /
!      DATA RMACH(2) / O377777777777 /
!      DATA RMACH(3) / O146400000000 /
!      DATA RMACH(4) / O147400000000 /
!      DATA RMACH(5) / O177464202324 /, SC/987/
!
      IF (SC .NE. 987) THEN 
!        *** CHECK FOR AUTODOUBLE ***
         SMALL(2) = 0
         RMACH(1) = 1E13
         IF (SMALL(2) .NE. 0) THEN
!           *** AUTODOUBLED ***
            IF (      SMALL(1) .EQ. 1117925532 &
               .AND. SMALL(2) .EQ. -448790528) THEN
!              *** IEEE BIG ENDIAN ***
               SMALL(1) = 1048576
               SMALL(2) = 0
               LARGE(1) = 2146435071
               LARGE(2) = -1
               RIGHT(1) = 1017118720
               RIGHT(2) = 0
               DIVER(1) = 1018167296
               DIVER(2) = 0
               LOG10(1) = 1070810131
               LOG10(2) = 1352628735
            ELSE IF ( SMALL(2) .EQ. 1117925532 &
               .AND. SMALL(1) .EQ. -448790528) THEN
!              *** IEEE LITTLE ENDIAN ***
               SMALL(2) = 1048576
               SMALL(1) = 0
               LARGE(2) = 2146435071
               LARGE(1) = -1
               RIGHT(2) = 1017118720
               RIGHT(1) = 0
               DIVER(2) = 1018167296
               DIVER(1) = 0
               LOG10(2) = 1070810131
               LOG10(1) = 1352628735
            ELSE IF ( SMALL(1) .EQ. -2065213935 &
               .AND. SMALL(2) .EQ. 10752) THEN
!              *** VAX WITH D_FLOATING ***
               SMALL(1) = 128
               SMALL(2) = 0
               LARGE(1) = -32769
               LARGE(2) = -1
               RIGHT(1) = 9344
               RIGHT(2) = 0
               DIVER(1) = 9472
               DIVER(2) = 0
               LOG10(1) = 546979738
               LOG10(2) = -805796613
            ELSE IF ( SMALL(1) .EQ. 1267827943 &
               .AND. SMALL(2) .EQ. 704643072) THEN
!              *** IBM MAINFRAME ***
               SMALL(1) = 1048576
               SMALL(2) = 0
               LARGE(1) = 2147483647
               LARGE(2) = -1
               RIGHT(1) = 856686592
               RIGHT(2) = 0
               DIVER(1) = 873463808
               DIVER(2) = 0
               LOG10(1) = 1091781651
               LOG10(2) = 1352628735
            ELSE
               WRITE(*,9010)
               STOP 777
               END IF
         ELSE
            RMACH(1) = 1234567.
            IF (SMALL(1) .EQ. 1234613304) THEN
!              *** IEEE ***
               SMALL(1) = 8388608
               LARGE(1) = 2139095039
               RIGHT(1) = 864026624
               DIVER(1) = 872415232
               LOG10(1) = 1050288283
            ELSE IF (SMALL(1) .EQ. -1271379306) THEN
!              *** VAX ***
               SMALL(1) = 128
               LARGE(1) = -32769
               RIGHT(1) = 13440
               DIVER(1) = 13568
               LOG10(1) = 547045274
            ELSE IF (SMALL(1) .EQ. 1175639687) THEN
!              *** IBM MAINFRAME ***
               SMALL(1) = 1048576
               LARGE(1) = 2147483647
               RIGHT(1) = 990904320
               DIVER(1) = 1007681536
               LOG10(1) = 1091781651
            ELSE IF (SMALL(1) .EQ. 1251390520) THEN
!              *** CONVEX C-1 ***
               SMALL(1) = 8388608
               LARGE(1) = 2147483647
               RIGHT(1) = 880803840
               DIVER(1) = 889192448
               LOG10(1) = 1067065499
            ELSE
               DO 10 L = 1, 3
                  J = SMALL(1) / 10000000
                  K = SMALL(1) - 10000000*J
                  IF (K .NE. T3E(L)) GO TO 20
                  SMALL(1) = J
 10               CONTINUE
!              *** CRAY T3E ***
               CALL I1MCRA(SMALL(1), K, 16, 0, 0)
               CALL I1MCRA(LARGE(1), K, 32751, 16777215, 16777215)
               CALL I1MCRA(RIGHT(1), K, 15520, 0, 0)
               CALL I1MCRA(DIVER(1), K, 15536, 0, 0)
               CALL I1MCRA(LOG10(1), K, 16339, 4461392, 10451455)
               GO TO 30
 20            CALL I1MCRA(J, K, 16405, 9876536, 0)
               IF (SMALL(1) .NE. J) THEN
                  WRITE(*,9020)
                  STOP 777
                  END IF
!              *** CRAY 1, XMP, 2, AND 3 ***
               CALL I1MCRA(SMALL(1), K, 8195, 8388608, 1)
               CALL I1MCRA(LARGE(1), K, 24574, 16777215, 16777214)
               CALL I1MCRA(RIGHT(1), K, 16338, 8388608, 0)
               CALL I1MCRA(DIVER(1), K, 16339, 8388608, 0)
               CALL I1MCRA(LOG10(1), K, 16383, 10100890, 8715216)
               END IF
            END IF
 30      SC = 987
         END IF
!     SANITY CHECK
      IF (RMACH(4) .GE. 1.0) STOP 776
      IF (I .LT. 1 .OR. I .GT. 5) THEN
         WRITE(*,*) 'R1MACH(I): I =',I,' is out of bounds.'
         STOP
         END IF
      R1MACH = RMACH(I)
      RETURN
 9010 FORMAT(/' Adjust autodoubled R1MACH by getting data'/ &
       ' appropriate for your machine from D1MACH.')
 9020 FORMAT(/' Adjust R1MACH by uncommenting data statements'/ &
       ' appropriate for your machine.')
! /* C source for R1MACH -- remove the * in column 1 */
!#include <stdio.h>
!#include <float.h>
!#include <math.h>
!float r1mach_(long *i)
!{
!       switch(*i){
!         case 1: return FLT_MIN;
!         case 2: return FLT_MAX;
!         case 3: return FLT_EPSILON/FLT_RADIX;
!         case 4: return FLT_EPSILON;
!         case 5: return log10((double)FLT_RADIX);
!         }
!       fprintf(stderr, "invalid argument: r1mach(%ld)\n", *i);
!       exit(1); return 0; /* else complaint of missing return value */
!}
      END

!**********************************************************************
!**********************************************************************
      SUBROUTINE I1MCRA(A, A1, B, C, D)
!**** SPECIAL COMPUTATION FOR CRAY MACHINES ****
      INTEGER A, A1, B, C, D
      A1 = 16777216*B + C
      A = 16777216*A1 + D
      END SUBROUTINE I1MCRA
      
!======================================================================
      function alog (x)
!======================================================================
! series for aln on the interval 0. to 3.46021d-03
!                                        with weighted error   1.50e-16
!                                         log weighted error  15.82
!                               significant figures required  15.65
!                                    decimal places required  16.21
! History:
!   june 1977 edition.   w. fullerton, c3, los alamos scientific lab.
!----------------------------------------------------------------------
      dimension alncs(6), center(4), alncen(5)
      external csevl, inits, r1mach
!
      data alncs( 1) /   1.3347199877973882e0   /
      data alncs( 2) /    .000693756283284112e0 /
      data alncs( 3) /    .000000429340390204e0 /
      data alncs( 4) /    .000000000289338477e0 /
      data alncs( 5) /    .000000000000205125e0 /
      data alncs( 6) /    .000000000000000150e0 /
!
      data center(1) / 1.0 /
      data center(2) / 1.25 /
      data center(3) / 1.50 /
      data center(4) / 1.75 /
!
      data alncen(  1) / 0.0e0                                    /
      data alncen(  2) /  .223143551314209755e0                   /
      data alncen(  3) /  .405465108108164381e0                   /
      data alncen(  4) /  .559615787935422686e0                   /
      data alncen(  5) /  .693147180559945309e0                   /
!
      data aln2 / 0.068147180559945309e0 /
      data nterms / 0 /
!
      if (nterms.eq.0) nterms = inits (alncs, 6, 28.9*r1mach(3))
!
      call r9upak (x, y, n)
!
      xn = n - 1
      y = 2.0*y
      ntrval = 4.0*y - 2.5
      if (ntrval.eq.5) t = ((y-1.0)-1.0) / (y+2.0)
      if (ntrval.lt.5) t = (y-center(ntrval))/(y+center(ntrval))
      t2 = t*t
!
      alog = 0.625*xn + (aln2*xn + alncen(ntrval) + 2.0*t + &
        t*t2*csevl(578.0*t2-1.0, alncs, nterms) )
!
      return
      end

!**********************************************************************
!**********************************************************************
      function csevl (x, cs, n)
! april 1977 version.  w. fullerton, c3, los alamos scientific lab.
!
! evaluate the n-term chebyshev series cs at x.  adapted from
! r. broucke, algorithm 446, c.a.c.m., 16, 254 (1973).  also see fox
! and parker, chebyshev polys in numerical analysis, oxford press, p.56.
!
!             input arguments --
! x      value at which the series is to be evaluated.
! cs     array of n terms of a chebyshev series.  in eval-
!        uating cs, only half the first coef is summed.
! n      number of terms in array cs.
!
      dimension cs(1)
!
      b1 = 0.
      b0 = 0.
      twox = 2.*x
      do 10 i=1,n
        b2 = b1
        b1 = b0
        ni = n + 1 - i
        b0 = twox*b1 - b2 + cs(ni)
 10   continue
!
      csevl = 0.5 * (b0-b2)
!
      return
      end

!**********************************************************************
!**********************************************************************
      subroutine d9gaml (xmin, xmax)
! june 1977 edition.   w. fullerton, c3, los alamos scientific lab.
!
! calculate the minimum and maximum legal bounds for x in gamma(x).
! xmin and xmax are not the only bounds, but they are the only non-
! trivial ones to calculate.
!
!             output arguments --
! xmin   dble prec minimum legal value of x in gamma(x).  any smaller
!        value of x might result in underflow.
! xmax   dble prec maximum legal value of x in gamma(x).  any larger
!        value of x might cause overflow.
!
      double precision xmin, xmax, alnbig, alnsml, xln, xold, d1mach, &
        dlog
      external d1mach, dlog
!
      alnsml = dlog(d1mach(1))
      xmin = -alnsml
      do 10 i=1,10
        xold = xmin
        xln = dlog(xmin)
        xmin = xmin - xmin*((xmin+0.5d0)*xln - xmin - 0.2258d0 + &
          alnsml) / (xmin*xln+0.5d0)
        if (dabs(xmin-xold).lt.0.005d0) go to 20
 10   continue
 20   xmin = -xmin + 0.01d0
!
      alnbig = dlog (d1mach(2))
      xmax = alnbig
      do 30 i=1,10
        xold = xmax
        xln = dlog(xmax)
        xmax = xmax - xmax*((xmax-0.5d0)*xln - xmax + 0.9189d0 - &
        alnbig) / (xmax*xln-0.5d0)
        if (dabs(xmax-xold).lt.0.005d0) go to 40
 30   continue
 40   xmax = xmax - 0.01d0
      xmin = dmax1 (xmin, -xmax+1.d0)
!
      return
      end

!**********************************************************************
!**********************************************************************
      double precision function d9lgmc (x)
! august 1977 edition.  w. fullerton, c3, los alamos scientific lab.
!
! compute the log gamma correction factor for x .ge. 10. so that
! dlog (dgamma(x)) = dlog(dsqrt(2*pi)) + (x-.5)*dlog(x) - x + d9lgmc(x)
!
      double precision x, algmcs(15), xbig, xmax, dcsevl, d1mach, &
        dexp, dlog, dsqrt
      external d1mach, dcsevl, dexp, dlog, dsqrt, initds
!
! series for algm       on the interval  0.          to  1.00000e-02
!                                        with weighted error   1.28e-31
!                                         log weighted error  30.89
!                               significant figures required  29.81
!                                    decimal places required  31.48
!
      data algmcs(  1) /  .1666389480451863247205729650822d+0      /
      data algmcs(  2) / -.1384948176067563840732986059135d-4      /
      data algmcs(  3) /  .9810825646924729426157171547487d-8      /
      data algmcs(  4) / -.1809129475572494194263306266719d-10     /
      data algmcs(  5) /  .6221098041892605227126015543416d-13     /
      data algmcs(  6) / -.3399615005417721944303330599666d-15     /
      data algmcs(  7) /  .2683181998482698748957538846666d-17     /
      data algmcs(  8) / -.2868042435334643284144622399999d-19     /
      data algmcs(  9) /  .3962837061046434803679306666666d-21     /
      data algmcs( 10) / -.6831888753985766870111999999999d-23     /
      data algmcs( 11) /  .1429227355942498147573333333333d-24     /
      data algmcs( 12) / -.3547598158101070547199999999999d-26     /
      data algmcs( 13) /  .1025680058010470912000000000000d-27     /
      data algmcs( 14) / -.3401102254316748799999999999999d-29     /
      data algmcs( 15) /  .1276642195630062933333333333333d-30     /
!
      data nalgm, xbig, xmax / 0, 2*0.d0 /
!
      if (nalgm.ne.0) go to 10
      nalgm = initds (algmcs, 15, sngl(d1mach(3)) )
      xbig = 1.0d0/dsqrt(d1mach(3))
      xmax = dexp (dmin1(dlog(d1mach(2)/12.d0), -dlog(12.d0*d1mach(1))))
!
 10   continue
      if (x.ge.xmax) go to 20
!
      d9lgmc = 1.d0/(12.d0*x)
      if (x.lt.xbig) d9lgmc = dcsevl (2.0d0*(10.d0/x)**2-1.d0, algmcs, &
        nalgm) / x
      return
!
 20   d9lgmc = 0.d0
      return
!
      end

!**********************************************************************
!**********************************************************************
      double precision function d9pak (y, n)
! december 1979 edition. w. fullerton, c3, los alamos scientific lab.
!
! pack a base 2 exponent into floating point number x.  this routine is
! almost the inverse of d9upak.  it is not exactly the inverse, because
! dabs(x) need not be between 0.5 and 1.0.  if both d9pak and 2.d0**n
! were known to be in range we could compute
!                d9pak = x * 2.0d0**n
!
      double precision y, aln2b, aln210, d1mach
      external d1mach, i1mach
      data nmin, nmax / 2*0 /
      data aln210 / 3.321928094887362347870319429489d0 /
!
      if (nmin.ne.0) go to 10
      aln2b = 1.0d0
      if (i1mach(10).ne.2) aln2b = d1mach(5)*aln210
      nmin = aln2b*dble(float(i1mach(15)))
      nmax = aln2b*dble(float(i1mach(16)))
!
 10   call d9upak (y, d9pak, ny)
!
      nsum = n + ny
      if (nsum.lt.nmin) go to 40
!
      if (nsum.eq.0) return
      if (nsum.gt.0) go to 30
!
 20   d9pak = 0.5d0*d9pak
      nsum = nsum + 1
      if (nsum.ne.0) go to 20
      return
!
 30   d9pak = 2.0d0*d9pak
      nsum = nsum - 1
      if (nsum.ne.0) go to 30
      return
!
 40   continue 
      d9pak = 0.0d0
      return
!
      end

!**********************************************************************
!**********************************************************************
      subroutine d9upak (x, y, n)
! august 1980 portable edition.  w. fullerton, los alamos scientific lab
!
! unpack floating point number x so that x = y * 2.0**n, where
! 0.5 .le. abs(y) .lt. 1.0 .
!
      double precision x, y, absx
!
      absx = dabs(x)
      n = 0
      y = 0.0d0
      if (x.eq.0.0d0) return
!
 10   if (absx.ge.0.5d0) go to 20
      n = n - 1
      absx = absx*2.0d0
      go to 10
!
 20   if (absx.lt.1.0d0) go to 30
      n = n + 1
      absx = absx*0.5d0
      go to 20
!
 30   y = dsign (absx, x)
      return
!
      end

!**********************************************************************
!**********************************************************************
      double precision function dcsevl (x, a, n)
!
! evaluate the n-term chebyshev series a at x.  adapted from
! r. broucke, algorithm 446, c.a.c.m., 16, 254 (1973).
!
!             input arguments --
! x      dble prec value at which the series is to be evaluated.
! a      dble prec array of n terms of a chebyshev series.  in eval-
!        uating a, only half the first coef is summed.
! n      number of terms in array a.
!
      double precision a(n), x, twox, b0, b1, b2
!
      twox = 2.0d0*x
      b1 = 0.d0
      b0 = 0.d0
      do 10 i=1,n
        b2 = b1
        b1 = b0
        ni = n - i + 1
        b0 = twox*b1 - b2 + a(ni)
 10   continue
!
      dcsevl = 0.5d0 * (b0-b2)
!
      return
      end

!**********************************************************************
!**********************************************************************
      double precision function dexp (x)
! may 1980 edition.   w. fullerton, c3, los alamos scientific lab.
      double precision x, expcs(14), twon16(17), aln216, f, xint, xmax,&
        xmin, y,  d1mach, dint, d9pak, dcsevl, dlog
      external d1mach, d9pak, dcsevl, dint, dlog, initds
!
! series for exp        on the interval -1.00000e+00 to  1.00000e+00
!                                        with weighted error   2.30e-34
!                                         log weighted error  33.64
!                               significant figures required  32.28
!                                    decimal places required  34.21
!
      data expcs(  1) / .866569493314985712733404647266231d-1    /
      data expcs(  2) / .938494869299839561896336579701203d-3    /
      data expcs(  3) / .677603970998168264074353014653601d-5    /
      data expcs(  4) / .366931200393805927801891250687610d-7    /
      data expcs(  5) / .158959053617461844641928517821508d-9    /
      data expcs(  6) / .573859878630206601252990815262106d-12   /
      data expcs(  7) / .177574448591421511802306980226000d-14   /
      data expcs(  8) / .480799166842372422675950244533333d-17   /
      data expcs(  9) / .115716376881828572809260000000000d-19   /
      data expcs( 10) / .250650610255497719932458666666666d-22   /
      data expcs( 11) / .493571708140495828480000000000000d-25   /
      data expcs( 12) / .890929572740634240000000000000000d-28   /
      data expcs( 13) / .148448062907997866666666666666666d-30   /
      data expcs( 14) / .229678916630186666666666666666666d-33   /
!
! twon16(i) is 2.0**((i-1)/16) - 1.0
      data twon16(  1) / 0.0d0                                   /
      data twon16(  2) / .44273782427413840321966478739929d-1    /
      data twon16(  3) / .90507732665257659207010655760707d-1    /
      data twon16(  4) / .13878863475669165370383028384151d0     /
      data twon16(  5) / .18920711500272106671749997056047d0     /
      data twon16(  6) / .24185781207348404859367746872659d0     /
      data twon16(  7) / .29683955465100966593375411779245d0     /
      data twon16(  8) / .35425554693689272829801474014070d0     /
      data twon16(  9) / .41421356237309504880168872420969d0     /
      data twon16( 10) / .47682614593949931138690748037404d0     /
      data twon16( 11) / .54221082540794082361229186209073d0     /
      data twon16( 12) / .61049033194925430817952066735740d0     /
      data twon16( 13) / .68179283050742908606225095246642d0     /
      data twon16( 14) / .75625216037329948311216061937531d0     /
      data twon16( 15) / .83400808640934246348708318958828d0     /
      data twon16( 16) / .91520656139714729387261127029583d0     /
      data twon16( 17) / 1.d0                                    /
!
! aln216 is 16.0/alog(2.) - 23.0
      data aln216 / .83120654223414517758794896030274d-1    /
      data nterms, xmin, xmax /0, 2*0.d0 /
!
      if (nterms.ne.0) go to 10
      nterms = initds (expcs, 14, 0.1*sngl(d1mach(3)))
      xmin = dlog (d1mach(1)) + .01d0
      xmax = dlog (d1mach(2)) - 0.001d0
!
 10   if (x.lt.xmin) go to 20
!
      xint = dint (x)
      y = x - xint
!
      y = 23.d0*y + x*aln216
      n = y
      f = y - dble(float(n))
      n = 23.d0*xint + dble(float(n))
      n16 = n/16
      if (n.lt.0) n16 = n16 - 1
      ndx = n - 16*n16 + 1
!
      dexp = 1.0d0 + (twon16(ndx) + f*(1.0d0 + twon16(ndx)) * &
        dcsevl (f, expcs, nterms) )
!
      dexp = d9pak (dexp, n16)
      return
!
 20   continue 
      dexp = 0.d0
      return
!
      end

!**********************************************************************
!**********************************************************************
      double precision function dgamma (x)
! jan 1984 edition.  w. fullerton, c3, los alamos scientific lab.
! jan 1994 wpp@ips.id.ethz.ch, ehg@research.att.com   declare xsml
      double precision x, gamcs(42), dxrel, pi, sinpiy, sq2pil, xmax, &
        xmin, y, d9lgmc, dcsevl, d1mach, dexp, dint, dlog, &
        dsin, dsqrt, xsml
      external d1mach, d9lgmc, dcsevl, dexp, dint, dlog, dsin, dsqrt, &
        initds
!
! series for gam        on the interval  0.          to  1.00000e+00
!                                        with weighted error   5.79e-32
!                                         log weighted error  31.24
!                               significant figures required  30.00
!                                    decimal places required  32.05
!
      data gamcs(  1) / .8571195590989331421920062399942d-2      /
      data gamcs(  2) / .4415381324841006757191315771652d-2      /
      data gamcs(  3) / .5685043681599363378632664588789d-1      /
      data gamcs(  4) / -.4219835396418560501012500186624d-2     /
      data gamcs(  5) / .1326808181212460220584006796352d-2      /
      data gamcs(  6) / -.1893024529798880432523947023886d-3     /
      data gamcs(  7) / .3606925327441245256578082217225d-4      /
      data gamcs(  8) / -.6056761904460864218485548290365d-5     /
      data gamcs(  9) / .1055829546302283344731823509093d-5      /
      data gamcs( 10) / -.1811967365542384048291855891166d-6     /
      data gamcs( 11) / .3117724964715322277790254593169d-7      /
      data gamcs( 12) / -.5354219639019687140874081024347d-8     /
      data gamcs( 13) / .9193275519859588946887786825940d-9      /
      data gamcs( 14) / -.1577941280288339761767423273953d-9     /
      data gamcs( 15) / .2707980622934954543266540433089d-10     /
      data gamcs( 16) / -.4646818653825730144081661058933d-11    /
      data gamcs( 17) / .7973350192007419656460767175359d-12     /
      data gamcs( 18) / -.1368078209830916025799499172309d-12    /
      data gamcs( 19) / .2347319486563800657233471771688d-13     /
      data gamcs( 20) / -.4027432614949066932766570534699d-14    /
      data gamcs( 21) / .6910051747372100912138336975257d-15     /
      data gamcs( 22) / -.1185584500221992907052387126192d-15    /
      data gamcs( 23) / .2034148542496373955201026051932d-16     /
      data gamcs( 24) / -.3490054341717405849274012949108d-17    /
      data gamcs( 25) / .5987993856485305567135051066026d-18     /
      data gamcs( 26) / -.1027378057872228074490069778431d-18    /
      data gamcs( 27) / .1762702816060529824942759660748d-19     /
      data gamcs( 28) / -.3024320653735306260958772112042d-20    /
      data gamcs( 29) / .5188914660218397839717833550506d-21     /
      data gamcs( 30) / -.8902770842456576692449251601066d-22    /
      data gamcs( 31) / .1527474068493342602274596891306d-22     /
      data gamcs( 32) / -.2620731256187362900257328332799d-23    /
      data gamcs( 33) / .4496464047830538670331046570666d-24     /
      data gamcs( 34) / -.7714712731336877911703901525333d-25    /
      data gamcs( 35) / .1323635453126044036486572714666d-25     /
      data gamcs( 36) / -.2270999412942928816702313813333d-26    /
      data gamcs( 37) / .3896418998003991449320816639999d-27     /
      data gamcs( 38) / -.6685198115125953327792127999999d-28    /
      data gamcs( 39) / .1146998663140024384347613866666d-28     /
      data gamcs( 40) / -.1967938586345134677295103999999d-29    /
      data gamcs( 41) / .3376448816585338090334890666666d-30     /
      data gamcs( 42) / -.5793070335782135784625493333333d-31    /
!
      data pi / 3.14159265358979323846264338327950d0 /
! sq2pil is 0.5*alog(2*pi) = alog(sqrt(2*pi))
      data sq2pil / 0.91893853320467274178032973640562d0 /
      data ngam, xmin, xmax, xsml, dxrel / 0, 4*0.d0 /
!
      if (ngam.ne.0) go to 10
      ngam = initds (gamcs, 42, 0.1*sngl(d1mach(3)) )
!
      call d9gaml (xmin, xmax)
      xsml = dexp (dmax1 (dlog(d1mach(1)), -dlog(d1mach(2)))+0.01d0)
      dxrel = dsqrt (d1mach(4))
!
 10   y = dabs(x)
      if (y.gt.10.d0) go to 50
!
! compute gamma(x) for -xbnd .le. x .le. xbnd.  reduce interval and find
! gamma(1+y) for 0.0 .le. y .lt. 1.0 first of all.
!
      n = x
      if (x.lt.0.d0) n = n - 1
      y = x - dble(float(n))
      n = n - 1
      dgamma = 0.9375d0 + dcsevl (2.d0*y-1.d0, gamcs, ngam)
      if (n.eq.0) return
!
      if (n.gt.0) go to 30
!
! compute gamma(x) for x .lt. 1.0
!
      n = -n
!
      do 20 i=1,n
        dgamma = dgamma/(x+dble(float(i-1)) )
 20   continue
      return
!
! gamma(x) for x .ge. 2.0 and x .le. 10.0
!
 30   do 40 i=1,n
        dgamma = (y+dble(float(i))) * dgamma
 40   continue
      return
!
! gamma(x) for dabs(x) .gt. 10.0.  recall y = dabs(x).
!
 50   continue 
!
      dgamma = 0.d0
      if (x.lt.xmin) return
!
      dgamma = dexp ((y-0.5d0)*dlog(y) - y + sq2pil + d9lgmc(y) )
      if (x.gt.0.d0) return
!
      sinpiy = dsin (pi*y)
!
      dgamma = -pi/(y*sinpiy*dgamma)
!
      return
      end

!**********************************************************************
!**********************************************************************
      double precision function dint (x)
! december 1983 edition. w. fullerton, c3, los alamos scientific lab.
!
! dint is the double precision equivalent of aint.  this portable
! version is quite efficient when the argument is reasonably small (a
! common case), and so no faster machine-dependent version is needed.
!
      double precision x, xscl, scale, xbig, xmax, part, d1mach, &
        dlog
      external d1mach, dlog, i1mach, r1mach
      data npart, scale, xbig, xmax / 0, 3*0.0d0 /
!
      if (npart.ne.0) go to 10
      ibase = i1mach(10)
      xmax = 1.0d0/d1mach (4)
      xbig = amin1 (float (i1mach(9)), 1.0/r1mach(4))
      scale = ibase**int(dlog(xbig)/dlog(dble(float(ibase)))-0.5d0)
      npart = dlog(xmax)/dlog(scale) + 1.0d0
!
 10   if (x.lt.(-xbig) .or. x.gt.xbig) go to 20
!
      dint = int(sngl(x))
      return
!
 20   xscl = dabs(x)
      if (xscl.gt.xmax) go to 50
!
      do 30 i=1,npart
        xscl = xscl/scale
 30   continue
!
      dint = 0.0d0
      do 40 i=1,npart
        xscl = xscl*scale
        ipart = xscl
        part = ipart
        xscl = xscl - part
        dint = dint*scale + part
 40   continue
!
      if (x.lt.0.0d0) dint = -dint
      return
!
 50   continue 
      dint = x
      return
!
      end

!**********************************************************************
!**********************************************************************
      double precision function dlog (x)
! june 1977 edition.   w. fullerton, c3, los alamos scientific lab.
      double precision x, alncs(11), center(4), alncen(5), aln2, y, t, &
        t2, xn,  dcsevl, d1mach
      external d1mach, dcsevl, initds
!
! series for aln        on the interval  0.          to  3.46021e-03
!                                        with weighted error   4.15e-32
!                                         log weighted error  31.38
!                               significant figures required  31.21
!                                    decimal places required  31.90
!
      data alncs(  1) / .13347199877973881561689386047187d+1     /
      data alncs(  2) / .69375628328411286281372438354225d-3     /
      data alncs(  3) / .42934039020450834506559210803662d-6     /
      data alncs(  4) / .28933847795432594580466440387587d-9     /
      data alncs(  5) / .20512517530340580901741813447726d-12    /
      data alncs(  6) / .15039717055497386574615153319999d-15    /
      data alncs(  7) / .11294540695636464284521613333333d-18    /
      data alncs(  8) / .86355788671171868881946666666666d-22    /
      data alncs(  9) / .66952990534350370613333333333333d-25    /
      data alncs( 10) / .52491557448151466666666666666666d-28    /
      data alncs( 11) / .41530540680362666666666666666666d-31    /
!
      data center(1) / 1.0d0 /
      data center(2) / 1.25d0 /
      data center(3) / 1.50d0 /
      data center(4) / 1.75d0 /
!
      data alncen(  1) / 0.0d0                                         /
      data alncen(  2) / .22314355131420975576629509030983d0     /
      data alncen(  3) / .40546510810816438197801311546434d0     /
      data alncen(  4) / .55961578793542268627088850052682d0     /
      data alncen(  5) / .69314718055994530941723212145817d0     /
!
! aln2 = alog(2.0) - 0.625
      data aln2 / 0.06814718055994530941723212145818d0 /
      data nterms / 0 /
!
      if (nterms.eq.0) nterms = initds (alncs, 11, 28.9*sngl(d1mach(3)))
!
      call d9upak (x, y, n)
!
      xn = n - 1
      y = 2.0d0*y
      ntrval = 4.0d0*y - 2.5d0
!
      if (ntrval.eq.5) t = ((y-1.0d0)-1.0d0) / (y+2.0d0)
      if (ntrval.lt.5) t = (y-center(ntrval)) / (y+center(ntrval))
      t2 = t*t
      dlog = 0.625d0*xn + (aln2*xn + alncen(ntrval) + 2.0d0*t + &
        t*t2*dcsevl(578.d0*t2-1.0d0, alncs, nterms) )
!
      return
      end

!**********************************************************************
!**********************************************************************
      double precision function dsin (x)
! august 1980 edition.  w. fullerton, los alamos scientific lab.
!
! this routine is based on the algorithm of cody and waite in
! argonne tm-321, software manual working note number 1
!
      double precision x, sincs(15), pihi, pilo, pirec, pi2rec, xsml, &
        xwarn, xmax, y, xn, sgn, f, dint, dcsevl, d1mach, dsqrt
      external d1mach, dcsevl, dint, dsqrt, initds
!
! series for sin    on the interval  0.00000e+00 to  2.46740e+00
!                                        with weighted error   2.56e-34
!                                         log weighted error  33.59
!                               significant figures required  33.01
!                                    decimal places required  34.18
!
      data sincs(  1) / -0.374991154955873175839919279977323464d0/
      data sincs(  2) / -0.181603155237250201863830316158004754d0/
      data sincs(  3) /  0.005804709274598633559427341722857921d0/
      data sincs(  4) / -0.000086954311779340757113212316353178d0/
      data sincs(  5) /  0.000000754370148088851481006839927030d0/
      data sincs(  6) / -0.000000004267129665055961107126829906d0/
      data sincs(  7) /  0.000000000016980422945488168181824792d0/
      data sincs(  8) / -0.000000000000050120578889961870929524d0/
      data sincs(  9) /  0.000000000000000114101026680010675628d0/
      data sincs( 10) / -0.000000000000000000206437504424783134d0/
      data sincs( 11) /  0.000000000000000000000303969595918706d0/
      data sincs( 12) / -0.000000000000000000000000371357734157d0/
      data sincs( 13) /  0.000000000000000000000000000382486123d0/
      data sincs( 14) / -0.000000000000000000000000000000336623d0/
      data sincs( 15) /  0.000000000000000000000000000000000256d0/
!
! pihi + pilo = pi.  pihi is exactly representable on all machines
! with at least 8 bits of precision.  whether it is exactly
! represented depends on the compiler.  this routine is more
! accurate if it is exactly represented.
      data pihi   / 3.140625d0 /
      data pilo   / 9.6765358979323846264338327950288d-4/
      data pirec  / 0.31830988618379067153776752674503d0 /
      data pi2rec / 0.63661977236758134307553505349006d0 /
      data ntsn, xsml, xwarn, xmax / 0, 3*0.0d0 /
!
      if (ntsn.ne.0) go to 10
      ntsn = initds (sincs, 15, 0.1*sngl(d1mach(3)))
!
      xsml = dsqrt (2.0d0*d1mach(3))
      xmax = 1.0d0/d1mach(4)
      xwarn = dsqrt (xmax)
!
 10   y = dabs (x)
!
      dsin = x
      if (y.lt.xsml) return
!
      xn = dint (y*pirec+0.5d0)
      n2 = dmod (xn, 2.0d0) + 0.5d0
      sgn = x
      if (n2.ne.0) sgn = -sgn
      f = (y-xn*pihi) - xn*pilo
!
      dsin = f + f*dcsevl(2.0d0*(f*pi2rec)**2-1.0d0, sincs, ntsn)
      if (sgn.lt.0.0d0) dsin = -dsin
      if (dabs(dsin).gt.1.0d0) dsin = dsign (1.0d0, dsin)
!
      return
      end

!**********************************************************************
!**********************************************************************
      double precision function dsqrt (x)
! june 1977 edition.   w. fullerton, c3, los alamos scientific lab.
      double precision x, sqrt2(3), y,  d9pak, d1mach
      external alog, d1mach, d9pak
      data sqrt2(1) / 0.70710678118654752440084436210485d0 /
      data sqrt2(2) / 1.0d0 /
      data sqrt2(3) / 1.41421356237309504880168872420970d0 /
!
      data niter / 0 /
!
      if (niter.eq.0) niter = 1.443*alog(-0.104*alog(0.1*sngl(d1mach(3)))) + 1.0
!
      if (x.le.0.d0) go to 20
!
      call d9upak (x, y, n)
      ixpnt = n/2
      irem = n - 2*ixpnt + 2
!
! the approximation below has accuracy of 4.16 digits.
      z = y
      dsqrt = .261599e0 + z*(1.114292e0 + z*(-.516888e0 + z*.141067e0))
!
      do 10 iter=1,niter
        dsqrt = dsqrt + 0.5d0*(y - dsqrt*dsqrt) / dsqrt
 10   continue
!
      dsqrt = d9pak (sqrt2(irem)*dsqrt, ixpnt)
      return
!
 20   continue
      dsqrt = 0.0d0
      return
!
      end

!**********************************************************************
!**********************************************************************
      integer function i8save(isw,ivalue,set)
!
!  if (isw = 1) i8save returns the current error number and
!               sets it to ivalue if set = .true. .
!
!  if (isw = 2) i8save returns the current recovery switch and
!               sets it to ivalue if set = .true. .
!
      logical set
!
      integer iparam(2)
!  iparam(1) is the error number and iparam(2) is the recovery switch.
!
!  start execution error free and with recovery turned off.
!
      data iparam(1) /0/,  iparam(2) /2/
!
      i8save=iparam(isw)
      if (set) iparam(isw)=ivalue
!
      return
!
      end

!**********************************************************************
!**********************************************************************
      function initds (dos, nos, eta)
! june 1977 edition.   w. fullerton, c3, los alamos scientific lab.
!
! initialize the double precision orthogonal series dos so that initds
! is the number of terms needed to insure the error is no larger than
! eta.  ordinarily eta will be chosen to be one-tenth machine precision.
!
!             input arguments --
! dos    dble prec array of nos coefficients in an orthogonal series.
! nos    number of coefficients in dos.
! eta    requested accuracy of series.
!
      double precision dos(nos)
!
      err = 0.
      do 10 ii=1,nos
        i = nos + 1 - ii
        err = err + abs(sngl(dos(i)))
        if (err.gt.eta) go to 20
 10   continue
!
 20   continue 
      initds = i
!
      return
      end

!**********************************************************************
!**********************************************************************
      function inits (os, nos, eta)
! april 1977 version.  w. fullerton, c3, los alamos scientific lab.
!
! initialize the orthogonal series so that inits is the number of terms
! needed to insure the error is no larger than eta.  ordinarily, eta
! will be chosen to be one-tenth machine precision.
!
!             input arguments --
! os     array of nos coefficients in an orthogonal series.
! nos    number of coefficients in os.
! eta    requested accuracy of series.
!
      dimension os(nos)
!
      err = 0.
      do 10 ii=1,nos
        i = nos + 1 - ii
        err = err + abs(os(i))
        if (err.gt.eta) go to 20
 10   continue
!
 20   continue
      inits = i
!
      return
      end

!**********************************************************************
!**********************************************************************
      subroutine r9upak (x, y, n)
! august 1980 portable edition.  w. fullerton, los alamos scientific lab
!
! unpack floating point number x so that x = y * 2.0**n, where
! 0.5 .le. abs(y) .lt. 1.0 .
!
      absx = abs(x)
      n = 0
      y = 0.0
      if (x.eq.0.0) return
!
 10   if (absx.ge.0.5) go to 20
      n = n - 1
      absx = absx*2.0
      go to 10
!
 20   if (absx.lt.1.0) go to 30
      n = n + 1
      absx = absx*0.5
      go to 20
!
 30   y = sign (absx, x)
      return
!
      end
