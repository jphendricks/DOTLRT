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

      END

!**********************************************************************
!**********************************************************************
      SUBROUTINE I1MCRA(A, A1, B, C, D)
!**** SPECIAL COMPUTATION FOR CRAY MACHINES ****
      INTEGER A, A1, B, C, D
      A1 = 16777216*B + C
      A = 16777216*A1 + D
      END SUBROUTINE I1MCRA
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
      real(8), parameter :: tol = 0.000001

      absx = abs(x)
      n = 0
      y = 0.0
      if (abs(x) > tol) return
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

