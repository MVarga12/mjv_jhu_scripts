
      real, allocatable :: X(:), F(:), D(:), Wk(:), Fe(:), Xe(:), Fd(:)
      logical spline

      parameter (spline=.false.)

      open (9,file='productregion.dat')
      N=0
      do
        read (9,*,iostat=ios) xxx,yyy; if (ios /= 0) exit; N=N+1
      enddo 
      rewind (9)
      Nf  = 2*N
      Nwk = 2*N
      allocate( X(N), F(N), D(N), Wk(Nwk), Fe(Nf), Xe(Nf), Fd(Nf) )

      do i=1,N
         read (9,*) X(i), F(i)
      enddo

      call pchez (N,X,F,D,spline,Wk,Nwk,ierr)
      dx = ( X(N)-X(1) ) / real(Nf-1)
      do i=1,Nf
        Xe(i) = X(1) + (i-1) * dx
      enddo

      call pchev (N,X,F,D,Nf,Xe,Fe,Fd,ierr)     ! this does the fitting

      open (10,file='fit.dat')
      do i=1,Nf
        write (10,*) Xe(i), Fe(i)               ! fitted function
      enddo

      a=X(1)
      b=X(N)
      rintegral = pchqa(N,X,F,D,a,b,ierr)       ! integration
      print *, rintegral

      end

c **********************************************************************

      SUBROUTINE PCHEZ(N,X,F,D,SPLINE,WK,LWK,IERR)
C***BEGIN PROLOGUE  PCHEZ
C***DATE WRITTEN   870821   (YYMMDD)
C***REVISION DATE  870908   (YYMMDD)
C***CATEGORY NO.  E1B
C***KEYWORDS  CUBIC HERMITE MONOTONE INTERPOLATION, SPLINE
C             INTERPOLATION, EASY TO USE PIECEWISE CUBIC INTERPOLATION
C***AUTHOR  KAHANER, D.K., (NBS)
C             SCIENTIFIC COMPUTING DIVISION
C             NATIONAL BUREAU OF STANDARDS
C             GAITHERSBURG, MARYLAND 20899
C             (301) 975-3808
C***PURPOSE  Easy to use spline or cubic Hermite interpolation.
C***DESCRIPTION
C
C          PCHEZ:  Piecewise Cubic Interpolation, Easy to Use.
C
C     From the book "Numerical Methods and Software"
C          by  D. Kahaner, C. Moler, S. Nash
C               Prentice Hall 1988
C
C     Sets derivatives for spline (two continuous derivatives) or
C     Hermite cubic (one continuous derivative) interpolation.
C     Spline interpolation is smoother, but may not "look" right if the
C     data contains both "steep" and "flat" sections.  Hermite cubics
C     can produce a "visually pleasing" and monotone interpolant to
C     monotone data. This is an easy to use driver for the routines
C     by F. N. Fritsch in reference (4) below. Various boundary
C     conditions are set to default values by PCHEZ. Many other choices
C     are available in the subroutines PCHIC, PCHIM and PCHSP.
C
C     Use PCHEV to evaluate the resulting function and its derivative.
C
C ----------------------------------------------------------------------
C
C  Calling sequence:   CALL  PCHEZ (N, X, F, D, SPLINE, WK, LWK, IERR)
C
C     INTEGER  N, IERR,  LWK
C     REAL  X(N), F(N), D(N), WK(*)
C     LOGICAL SPLINE
C
C   Parameters:
C
C     N -- (input) number of data points.  (Error return if N.LT.2 .)
C           If N=2, simply does linear interpolation.
C
C     X -- (input) real array of independent variable values.  The
C           elements of X must be strictly increasing:
C                X(I-1) .LT. X(I),  I = 2(1)N.
C           (Error return if not.)
C
C     F -- (input) real array of dependent variable values to be inter-
C           polated.  F(I) is value corresponding to X(I).
C
C     D -- (output) real array of derivative values at the data points.
C
C     SPLINE -- (input) logical variable to specify if the interpolant
C           is to be a spline with two continuous derivaties
C           (set SPLINE=.TRUE.) or a Hermite cubic interpolant with one
C           continuous derivative (set SPLINE=.FALSE.).
C        Note: If SPLINE=.TRUE. the interpolating spline satisfies the
C           default "not-a-knot" boundary condition, with a continuous
C           third derivative at X(2) and X(N-1). See reference (3).
C              If SPLINE=.FALSE. the interpolating Hermite cubic will be
C           monotone if the input data is monotone. Boundary conditions are
C           computed from the derivative of a local quadratic unless this
C           alters monotonicity.
C
C     WK -- (scratch) real work array, which must be declared by the calling
C           program to be at least 2*N if SPLINE is .TRUE. and not used
C           otherwise.
C
C     LWK -- (input) length of work array WK. (Error return if
C           LWK.LT.2*N and SPLINE is .TRUE., not checked otherwise.)
C
C     IERR -- (output) error flag.
C           Normal return:
C              IERR = 0  (no errors).
C           Warning error:
C              IERR.GT.0  (can only occur when SPLINE=.FALSE.) means that
C                 IERR switches in the direction of monotonicity were detected.
C                 When SPLINE=.FALSE.,  PCHEZ guarantees that if the input
C                 data is monotone, the interpolant will be too. This warning
C                 is to alert you to the fact that the input data was not
C                 monotone.
C           "Recoverable" errors:
C              IERR = -1  if N.LT.2 .
C              IERR = -3  if the X-array is not strictly increasing.
C              IERR = -7  if LWK is less than 2*N and SPLINE is .TRUE.
C             (The D-array has not been changed in any of these cases.)
C               NOTE:  The above errors are checked in the order listed,
C                   and following arguments have **NOT** been validated.
C
C ----------------------------------------------------------------------
C***REFERENCES  1. F.N.FRITSCH AND R.E.CARLSON, 'MONOTONE PIECEWISE
C                 CUBIC INTERPOLATION,' SIAM J.NUMER.ANAL. 17, 2 (APRIL
C                 1980), 238-246.
C               2. F.N.FRITSCH AND J.BUTLAND, 'A METHOD FOR CONSTRUCTING
C                 LOCAL MONOTONE PIECEWISE CUBIC INTERPOLANTS,' LLNL
C                 PREPRINT UCRL-87559 (APRIL 1982).
C               3. CARL DE BOOR, A PRACTICAL GUIDE TO SPLINES, SPRINGER-
C                 VERLAG (NEW YORK, 1978).  (ESP. CHAPTER IV, PP.49-62.)
C               4. F.N.FRITSCH, 'PIECEWISE CUBIC HERMITE INTERPOLATION
C                 PACKAGE, FINAL SPECIFICATIONS', LAWRENCE LIVERMORE
C                 NATIONAL LABORATORY, COMPUTER DOCUMENTATION UCID-30194,
C                 AUGUST 1982.
C***ROUTINES CALLED  PCHIM,PCHSP
C***END PROLOGUE  PCHEZ
      INTEGER  N, LWK, IERR
      REAL  X(N), F(N), D(N), WK(LWK)
      LOGICAL SPLINE
C
C  DECLARE LOCAL VARIABLES.
C
      INTEGER IC(2), INCFD
      REAL  VC(2)
      DATA IC(1) /0/
      DATA IC(2) /0/
      DATA INCFD /1/
C
C
C***FIRST EXECUTABLE STATEMENT  PCHEZ
C
      IF ( SPLINE ) THEN
        CALL  PCHSP (IC, VC, N, X, F, D, INCFD, WK, LWK, IERR)
      ELSE
        CALL  PCHIM (N, X, F, D, INCFD, IERR)
      ENDIF
C
C  ERROR CONDITIONS ALREADY CHECKED IN PCHSP OR PCHIM

      RETURN
C------------- LAST LINE OF PCHEZ FOLLOWS ------------------------------
      END

      SUBROUTINE PCHEV(N,X,F,D,NVAL,XVAL,FVAL,DVAL,IERR)
C***BEGIN PROLOGUE  PCHEV
C***DATE WRITTEN   870828   (YYMMDD)
C***REVISION DATE  870828   (YYMMDD)
C***CATEGORY NO.  E3,H1
C***KEYWORDS  CUBIC HERMITE OR SPLINE DIFFERENTIATION,CUBIC HERMITE
C             EVALUATION,EASY TO USE SPLINE OR CUBIC HERMITE EVALUATOR
C***AUTHOR  KAHANER, D.K., (NBS)
C             SCIENTIFIC COMPUTING DIVISION
C             NATIONAL BUREAU OF STANDARDS
C             ROOM A161, TECHNOLOGY BUILDING
C             GAITHERSBURG, MARYLAND 20899
C             (301) 975-3808
C***PURPOSE  Evaluates the function and first derivative of a piecewise
C            cubic Hermite or spline function at an array of points XVAL,
C            easy to use.
C***DESCRIPTION
C
C          PCHEV:  Piecewise Cubic Hermite or Spline Derivative Evaluator,
C                  Easy to Use.
C
C     From the book "Numerical Methods and Software"
C          by  D. Kahaner, C. Moler, S. Nash
C                 Prentice Hall 1988
C
C     Evaluates the function and first derivative of the cubic Hermite
C     or spline function defined by  N, X, F, D, at the array of points XVAL.
C
C     This is an easy to use driver for the routines by F.N. Fritsch
C     described in reference (2) below. Those also have other capabilities.
C
C ----------------------------------------------------------------------
C
C  Calling sequence: CALL  PCHEV (N, X, F, D, NVAL, XVAL, FVAL, DVAL, IERR)
C
C     INTEGER  N, NVAL, IERR
C     REAL  X(N), F(N), D(N), XVAL(NVAL), FVAL(NVAL), DVAL(NVAL)
C
C   Parameters:
C
C     N -- (input) number of data points.  (Error return if N.LT.2 .)
C
C     X -- (input) real array of independent variable values.  The
C           elements of X must be strictly increasing:
C             X(I-1) .LT. X(I),  I = 2(1)N. (Error return if not.)
C
C     F -- (input) real array of function values.  F(I) is
C           the value corresponding to X(I).
C
C     D -- (input) real array of derivative values.  D(I) is
C           the value corresponding to X(I).
C
C  NVAL -- (input) number of points at which the functions are to be
C           evaluated. ( Error return if NVAL.LT.1 )
C
C  XVAL -- (input) real array of points at which the functions are to
C           be evaluated.
C
C          NOTES:
C           1. The evaluation will be most efficient if the elements
C              of XVAL are increasing relative to X;
C              that is,   XVAL(J) .GE. X(I)
C              implies    XVAL(K) .GE. X(I),  all K.GE.J .
C           2. If any of the XVAL are outside the interval [X(1),X(N)],
C              values are extrapolated from the nearest extreme cubic,
C              and a warning error is returned.
C
C  FVAL -- (output) real array of values of the cubic Hermite function
C           defined by  N, X, F, D  at the points  XVAL.
C
C  DVAL -- (output) real array of values of the first derivative of
C           the same function at the points  XVAL.
C
C  IERR -- (output) error flag.
C           Normal return:
C              IERR = 0  (no errors).
C           Warning error:
C              IERR.GT.0  means that extrapolation was performed at
C                 IERR points.
C           "Recoverable" errors:
C              IERR = -1  if N.LT.2 .
C              IERR = -3  if the X-array is not strictly increasing.
C              IERR = -4  if NVAL.LT.1 .
C           (Output arrays have not been changed in any of these cases.)
C               NOTE:  The above errors are checked in the order listed,
C                   and following arguments have **NOT** been validated.
C              IERR = -5  if an error has occurred in the lower-level
C                         routine CHFDV.  NB: this should never happen.
C                         Notify the author **IMMEDIATELY** if it does.
C
C ----------------------------------------------------------------------
C***REFERENCES  1. F.N.FRITSCH AND R.E.CARLSON, 'MONOTONE PIECEWISE
C                 CUBIC INTERPOLATION,' SIAM J.NUMER.ANAL. 17, 2 (APRIL
C                 1980), 238-246.
C               2. F.N.FRITSCH, 'PIECEWISE CUBIC HERMITE INTERPOLATION
C                 PACKAGE, FINAL SPECIFICATIONS', LAWRENCE LIVERMORE
C                 NATIONAL LABORATORY, COMPUTER DOCUMENTATION UCID-30194,
C                 AUGUST 1982.
C***ROUTINES CALLED  PCHFD
C***END PROLOGUE  PCHEV
      INTEGER  N, NVAL, IERR
      REAL  X(N), F(N), D(N), XVAL(NVAL), FVAL(NVAL), DVAL(NVAL)
C
C  DECLARE LOCAL VARIABLES.
C
      INTEGER INCFD
      LOGICAL SKIP
      DATA SKIP /.TRUE./
      DATA INCFD /1/

C
C
C***FIRST EXECUTABLE STATEMENT  PCHEV
C
      CALL PCHFD(N,X,F,D,INCFD,SKIP,NVAL,XVAL,FVAL,DVAL,IERR)
C
C
 5000 CONTINUE
      RETURN
C
C------------- LAST LINE OF PCHEV FOLLOWS ------------------------------
      END

      REAL FUNCTION PCHQA(N,X,F,D,A,B,IERR)
C***BEGIN PROLOGUE  PCHQA
C***DATE WRITTEN   870829   (YYMMDD)
C***REVISION DATE  870829   (YYMMDD)
C***CATEGORY NO.  E3,H2A2
C***KEYWORDS  EASY TO USE CUBIC HERMITE OR SPLINE INTEGRATION
C             NUMERICAL INTEGRATION, QUADRATURE
C***AUTHOR  KAHANER, D.K., (NBS)
C             SCIENTIFIC COMPUTING DIVISION
C             NATIONAL BUREAU OF STANDARDS
C             ROOM A161, TECHNOLOGY BUILDING
C             GAITHERSBURG, MARYLAND 20899
C             (301) 975-3808
C***PURPOSE  Evaluates the definite integral of a piecewise cubic Hermite
C            or spline function over an arbitrary interval, easy to use.
C***DESCRIPTION
C
C          PCHQA:  Piecewise Cubic Hermite or Spline Integrator,
C                  Arbitrary limits, Easy to Use.
C
C          From the book "Numerical Methods and Software"
C                  by  D. Kahaner, C. Moler, S. Nash
C                          Prentice Hall 1988
C
C     Evaluates the definite integral of the cubic Hermite or spline
C     function defined by  N, X, F, D  over the interval [A, B].  This
C     is an easy to use driver for the routine PCHIA by F.N. Fritsch
C     described in reference (2) below. That routine also has other
C     capabilities.
C ----------------------------------------------------------------------
C
C  Calling sequence:
C
C           VALUE = PCHQA (N, X, F, D, A, B, IERR)
C
C     INTEGER  N, IERR
C     REAL  X(N), F(N), D(N), A, B
C
C   Parameters:
C
C     VALUE -- (output) VALUE of the requested integral.
C
C     N -- (input) number of data points.  (Error return if N.LT.2 .)
C
C     X -- (input) real array of independent variable values.  The
C           elements of X must be strictly increasing:
C                X(I-1) .LT. X(I),  I = 2(1)N.
C           (Error return if not.)
C
C     F -- (input) real array of function values.  F(I) is
C           the value corresponding to X(I).
C
C     D -- (input) real array of derivative values.  D(I) is
C           the value corresponding to X(I).
C
C     A,B -- (input) the limits of integration.
C           NOTE:  There is no requirement that [A,B] be contained in
C                  [X(1),X(N)].  However, the resulting integral value
C                  will be highly suspect, if not.
C
C     IERR -- (output) error flag.
C           Normal return:
C              IERR = 0  (no errors).
C           Warning errors:
C              IERR = 1  if  A  is outside the interval [X(1),X(N)].
C              IERR = 2  if  B  is outside the interval [X(1),X(N)].
C              IERR = 3  if both of the above are true.  (Note that this
C                        means that either [A,B] contains data interval
C                        or the intervals do not intersect at all.)
C           "Recoverable" errors:
C              IERR = -1  if N.LT.2 .
C              IERR = -3  if the X-array is not strictly increasing.
C                (Value has not been computed in any of these cases.)
C               NOTE:  The above errors are checked in the order listed,
C                   and following arguments have **NOT** been validated.
C
C***REFERENCES  1. F.N.FRITSCH AND R.E.CARLSON, 'MONOTONE PIECEWISE
C                 CUBIC INTERPOLATION,' SIAM J.NUMER.ANAL. 17, 2 (APRIL
C                 1980), 238-246.
C               2. F.N.FRITSCH, 'PIECEWISE CUBIC HERMITE INTERPOLATION
C                 PACKAGE, FINAL SPECIFICATIONS', LAWRENCE LIVERMORE
C                 NATIONAL LABORATORY, COMPUTER DOCUMENTATION UCID-30194,
C                 AUGUST 1982.
C***ROUTINES CALLED  PCHIA
C***END PROLOGUE  PCHQA
      INTEGER  N, IERR
      REAL  X(N), F(N), D(N), A, B
C
C  DECLARE LOCAL VARIABLES.
C
      INTEGER  INCFD
      REAL  PCHIA
      LOGICAL SKIP
C
C  INITIALIZE.
C
      DATA  INCFD /1/
      DATA  SKIP /.TRUE./
C
C
C***FIRST EXECUTABLE STATEMENT  PCHQA

      PCHQA  =  PCHIA( N, X, F, D, INCFD, SKIP, A, B, IERR )
C
C ERROR MESSAGES ARE FROM LOWER LEVEL ROUTINES
      RETURN
C
C------------- LAST LINE OF PCHQA FOLLOWS ------------------------------
      END
c ***********************************************************************

      SUBROUTINE PCHIM(N,X,F,D,INCFD,IERR)
C***BEGIN PROLOGUE  PCHIM
C***DATE WRITTEN   811103   (YYMMDD)
C***REVISION DATE  870707   (YYMMDD)
C***CATEGORY NO.  E1B
C***KEYWORDS  LIBRARY=SLATEC(PCHIP),
C             TYPE=SINGLE PRECISION(PCHIM-S DPCHIM-D),
C             CUBIC HERMITE INTERPOLATION,MONOTONE INTERPOLATION,
C             PIECEWISE CUBIC INTERPOLATION
C***AUTHOR  FRITSCH, F. N., (LLNL)
C             MATHEMATICS AND STATISTICS DIVISION
C             LAWRENCE LIVERMORE NATIONAL LABORATORY
C             P.O. BOX 808  (L-316)
C             LIVERMORE, CA  94550
C             FTS 532-4275, (415) 422-4275
C***PURPOSE  Set derivatives needed to determine a monotone piecewise
C            cubic Hermite interpolant to given data.  Boundary values
C            are provided which are compatible with monotonicity.  The
C            interpolant will have an extremum at each point where mono-
C            tonicity switches direction.  (See PCHIC if user control is
C            desired over boundary or switch conditions.)
C***DESCRIPTION
C
C          PCHIM:  Piecewise Cubic Hermite Interpolation to
C                  Monotone data.
C
C     Sets derivatives needed to determine a monotone piecewise cubic
C     Hermite interpolant to the data given in X and F.
C
C     Default boundary conditions are provided which are compatible
C     with monotonicity.  (See PCHIC if user control of boundary con-
C     ditions is desired.)
C
C     If the data are only piecewise monotonic, the interpolant will
C     have an extremum at each point where monotonicity switches direc-
C     tion.  (See PCHIC if user control is desired in such cases.)
C
C     To facilitate two-dimensional applications, includes an increment
C     between successive values of the F- and D-arrays.
C
C     The resulting piecewise cubic Hermite function may be evaluated
C     by PCHFE or PCHFD.
C
C ----------------------------------------------------------------------
C
C  Calling sequence:
C
C        PARAMETER  (INCFD = ...)
C        INTEGER  N, IERR
C        REAL  X(N), F(INCFD,N), D(INCFD,N)
C
C        CALL  PCHIM (N, X, F, D, INCFD, IERR)
C
C   Parameters:
C
C     N -- (input) number of data points.  (Error return if N.LT.2 .)
C           If N=2, simply does linear interpolation.
C
C     X -- (input) real array of independent variable values.  The
C           elements of X must be strictly increasing:
C                X(I-1) .LT. X(I),  I = 2(1)N.
C           (Error return if not.)
C
C     F -- (input) real array of dependent variable values to be inter-
C           polated.  F(1+(I-1)*INCFD) is value corresponding to X(I).
C           PCHIM is designed for monotonic data, but it will work for
C           any F-array.  It will force extrema at points where mono-
C           tonicity switches direction.  If some other treatment of
C           switch points is desired, PCHIC should be used instead.
C                                     -----
C     D -- (output) real array of derivative values at the data points.
C           If the data are monotonic, these values will determine a
C           a monotone cubic Hermite function.
C           The value corresponding to X(I) is stored in
C                D(1+(I-1)*INCFD),  I=1(1)N.
C           No other entries in D are changed.
C
C     INCFD -- (input) increment between successive values in F and D.
C           This argument is provided primarily for 2-D applications.
C           (Error return if  INCFD.LT.1 .)
C
C     IERR -- (output) error flag.
C           Normal return:
C              IERR = 0  (no errors).
C           Warning error:
C              IERR.GT.0  means that IERR switches in the direction
C                 of monotonicity were detected.
C           "Recoverable" errors:
C              IERR = -1  if N.LT.2 .
C              IERR = -2  if INCFD.LT.1 .
C              IERR = -3  if the X-array is not strictly increasing.
C             (The D-array has not been changed in any of these cases.)
C               NOTE:  The above errors are checked in the order listed,
C                   and following arguments have **NOT** been validated.
C
C***REFERENCES  1. F.N.FRITSCH AND R.E.CARLSON, 'MONOTONE PIECEWISE
C                 CUBIC INTERPOLATION,' SIAM J.NUMER.ANAL. 17, 2 (APRIL
C                 1980), 238-246.
C               2. F.N.FRITSCH AND J.BUTLAND, 'A METHOD FOR CONSTRUCTING
C                 LOCAL MONOTONE PIECEWISE CUBIC INTERPOLANTS,' LLNL
C                 PREPRINT UCRL-87559 (APRIL 1982).
C***ROUTINES CALLED  PCHST,XERROR
C***END PROLOGUE  PCHIM
C
C ----------------------------------------------------------------------
C
C  Change record:
C     82-02-01   1. Introduced  PCHST  to reduce possible over/under-
C                   flow problems.
C                2. Rearranged derivative formula for same reason.
C     82-06-02   1. Modified end conditions to be continuous functions
C                   of data when monotonicity switches in next interval.
C                2. Modified formulas so end conditions are less prone
C                   of over/underflow problems.
C     82-08-03   Minor cosmetic changes for release 1.
C
C ----------------------------------------------------------------------
C
C  Programming notes:
C
C     1. The function  PCHST(ARG1,ARG2)  is assumed to return zero if
C        either argument is zero, +1 if they are of the same sign, and
C        -1 if they are of opposite sign.
C     2. To produce a double precision version, simply:
C        a. Change PCHIM to DPCHIM wherever it occurs,
C        b. Change PCHST to DPCHST wherever it occurs,
C        c. Change all references to the Fortran intrinsics to their
C           double precision equivalents,
C        d. Change the real declarations to double precision, and
C        e. Change the constants ZERO and THREE to double precision.
C
C  DECLARE ARGUMENTS.
C
      INTEGER  N, INCFD, IERR
      REAL  X(N), F(INCFD,N), D(INCFD,N)
C
C  DECLARE LOCAL VARIABLES.
C
      INTEGER  I, NLESS1
      REAL  DEL1, DEL2, DMAX, DMIN, DRAT1, DRAT2, DSAVE,
     *      H1, H2, HSUM, HSUMT3, THREE, W1, W2, ZERO
      REAL  PCHST
      DATA  ZERO /0./,  THREE /3./
C
C  VALIDITY-CHECK ARGUMENTS.
C
C***FIRST EXECUTABLE STATEMENT  PCHIM
      IF ( N.LT.2 )  GO TO 5001
      IF ( INCFD.LT.1 )  GO TO 5002
      DO 1  I = 2, N
         IF ( X(I).LE.X(I-1) )  GO TO 5003
    1 CONTINUE
C
C  FUNCTION DEFINITION IS OK, GO ON.
C
      IERR = 0
      NLESS1 = N - 1
      H1 = X(2) - X(1)
      DEL1 = (F(1,2) - F(1,1))/H1
      DSAVE = DEL1
C
C  SPECIAL CASE N=2 -- USE LINEAR INTERPOLATION.
C
      IF (NLESS1 .GT. 1)  GO TO 10
      D(1,1) = DEL1
      D(1,N) = DEL1
      GO TO 5000
C
C  NORMAL CASE  (N .GE. 3).
C
   10 CONTINUE
      H2 = X(3) - X(2)
      DEL2 = (F(1,3) - F(1,2))/H2
C
C  SET D(1) VIA NON-CENTERED THREE-POINT FORMULA, ADJUSTED TO BE
C     SHAPE-PRESERVING.
C
      HSUM = H1 + H2
      W1 = (H1 + HSUM)/HSUM
      W2 = -H1/HSUM
      D(1,1) = W1*DEL1 + W2*DEL2
      IF ( PCHST(D(1,1),DEL1) .LE. ZERO)  THEN
         D(1,1) = ZERO
      ELSE IF ( PCHST(DEL1,DEL2) .LT. ZERO)  THEN
C        NEED DO THIS CHECK ONLY IF MONOTONICITY SWITCHES.
         DMAX = THREE*DEL1
         IF (ABS(D(1,1)) .GT. ABS(DMAX))  D(1,1) = DMAX
      ENDIF
C
C  LOOP THROUGH INTERIOR POINTS.
C
      DO 50  I = 2, NLESS1
         IF (I .EQ. 2)  GO TO 40
C
         H1 = H2
         H2 = X(I+1) - X(I)
         HSUM = H1 + H2
         DEL1 = DEL2
         DEL2 = (F(1,I+1) - F(1,I))/H2
   40    CONTINUE
C
C        SET D(I)=0 UNLESS DATA ARE STRICTLY MONOTONIC.
C
         D(1,I) = ZERO
         IF ( PCHST(DEL1,DEL2) )  42, 41, 45
C
C        COUNT NUMBER OF CHANGES IN DIRECTION OF MONOTONICITY.
C
   41    CONTINUE
         IF (DEL2 .EQ. ZERO)  GO TO 50
         IF ( PCHST(DSAVE,DEL2) .LT. ZERO)  IERR = IERR + 1
         DSAVE = DEL2
         GO TO 50
C
   42    CONTINUE
         IERR = IERR + 1
         DSAVE = DEL2
         GO TO 50
C
C        USE BRODLIE MODIFICATION OF BUTLAND FORMULA.
C
   45    CONTINUE
         HSUMT3 = HSUM+HSUM+HSUM
         W1 = (HSUM + H1)/HSUMT3
         W2 = (HSUM + H2)/HSUMT3
         DMAX = AMAX1( ABS(DEL1), ABS(DEL2) )
         DMIN = AMIN1( ABS(DEL1), ABS(DEL2) )
         DRAT1 = DEL1/DMAX
         DRAT2 = DEL2/DMAX
         D(1,I) = DMIN/(W1*DRAT1 + W2*DRAT2)
C
   50 CONTINUE
C
C  SET D(N) VIA NON-CENTERED THREE-POINT FORMULA, ADJUSTED TO BE
C     SHAPE-PRESERVING.
C
      W1 = -H2/HSUM
      W2 = (H2 + HSUM)/HSUM
      D(1,N) = W1*DEL1 + W2*DEL2
      IF ( PCHST(D(1,N),DEL2) .LE. ZERO)  THEN
         D(1,N) = ZERO
      ELSE IF ( PCHST(DEL1,DEL2) .LT. ZERO)  THEN
C        NEED DO THIS CHECK ONLY IF MONOTONICITY SWITCHES.
         DMAX = THREE*DEL2
         IF (ABS(D(1,N)) .GT. ABS(DMAX))  D(1,N) = DMAX
      ENDIF
C
C  NORMAL RETURN.
C
 5000 CONTINUE
      RETURN
C
C  ERROR RETURNS.
C
 5001 CONTINUE
C     N.LT.2 RETURN.
      IERR = -1
      CALL XERROR ('PCHIM -- NUMBER OF DATA POINTS LESS THAN TWO'
     *           , 44, IERR, 1)
      RETURN
C
 5002 CONTINUE
C     INCFD.LT.1 RETURN.
      IERR = -2
      CALL XERROR ('PCHIM -- INCREMENT LESS THAN ONE'
     *           , 32, IERR, 1)
      RETURN
C
 5003 CONTINUE
C     X-ARRAY NOT STRICTLY INCREASING.
      IERR = -3
      CALL XERROR ('PCHIM -- X-ARRAY NOT STRICTLY INCREASING'
     *           , 40, IERR, 1)
      RETURN
C------------- LAST LINE OF PCHIM FOLLOWS ------------------------------
      END

      SUBROUTINE PCHSP(IC,VC,N,X,F,D,INCFD,WK,NWK,IERR)
C***BEGIN PROLOGUE  PCHSP
C***DATE WRITTEN   820503   (YYMMDD)
C***REVISION DATE  870707   (YYMMDD)
C***CATEGORY NO.  E1B
C***KEYWORDS  LIBRARY=SLATEC(PCHIP),
C             TYPE=SINGLE PRECISION(PCHSP-S DPCHSP-D),
C             CUBIC HERMITE INTERPOLATION,PIECEWISE CUBIC INTERPOLATION,
C             SPLINE INTERPOLATION
C***AUTHOR  FRITSCH, F. N., (LLNL)
C             MATHEMATICS AND STATISTICS DIVISION
C             LAWRENCE LIVERMORE NATIONAL LABORATORY
C             P.O. BOX 808  (L-316)
C             LIVERMORE, CA  94550
C             FTS 532-4275, (415) 422-4275
C***PURPOSE  Set derivatives needed to determine the Hermite represen-
C            tation of the cubic spline interpolant to given data, with
C            specified boundary conditions.
C***DESCRIPTION
C
C          PCHSP:   Piecewise Cubic Hermite Spline
C
C     Computes the Hermite representation of the cubic spline inter-
C     polant to the data given in X and F satisfying the boundary
C     conditions specified by IC and VC.
C
C     To facilitate two-dimensional applications, includes an increment
C     between successive values of the F- and D-arrays.
C
C     The resulting piecewise cubic Hermite function may be evaluated
C     by PCHFE or PCHFD.
C
C     NOTE:  This is a modified version of C. de Boor'S cubic spline
C            routine CUBSPL.
C
C ----------------------------------------------------------------------
C
C  Calling sequence:
C
C        PARAMETER  (INCFD = ...)
C        INTEGER  IC(2), N, NWK, IERR
C        REAL  VC(2), X(N), F(INCFD,N), D(INCFD,N), WK(NWK)
C
C        CALL  PCHSP (IC, VC, N, X, F, D, INCFD, WK, NWK, IERR)
C
C   Parameters:
C
C     IC -- (input) integer array of length 2 specifying desired
C           boundary conditions:
C           IC(1) = IBEG, desired condition at beginning of data.
C           IC(2) = IEND, desired condition at end of data.
C
C           IBEG = 0  to set D(1) so that the third derivative is con-
C              tinuous at X(2).  This is the "not a knot" condition
C              provided by de Boor'S cubic spline routine CUBSPL.
C              < This is the default boundary condition. >
C           IBEG = 1  if first derivative at X(1) is given in VC(1).
C           IBEG = 2  if second derivative at X(1) is given in VC(1).
C           IBEG = 3  to use the 3-point difference formula for D(1).
C                     (Reverts to the default b.c. if N.LT.3 .)
C           IBEG = 4  to use the 4-point difference formula for D(1).
C                     (Reverts to the default b.c. if N.LT.4 .)
C          NOTES:
C           1. An error return is taken if IBEG is out of range.
C           2. For the "natural" boundary condition, use IBEG=2 and
C              VC(1)=0.
C
C           IEND may take on the same values as IBEG, but applied to
C           derivative at X(N).  In case IEND = 1 or 2, the value is
C           given in VC(2).
C
C          NOTES:
C           1. An error return is taken if IEND is out of range.
C           2. For the "natural" boundary condition, use IEND=2 and
C              VC(2)=0.
C
C     VC -- (input) real array of length 2 specifying desired boundary
C           values, as indicated above.
C           VC(1) need be set only if IC(1) = 1 or 2 .
C           VC(2) need be set only if IC(2) = 1 or 2 .
C
C     N -- (input) number of data points.  (Error return if N.LT.2 .)
C
C     X -- (input) real array of independent variable values.  The
C           elements of X must be strictly increasing:
C                X(I-1) .LT. X(I),  I = 2(1)N.
C           (Error return if not.)
C
C     F -- (input) real array of dependent variable values to be inter-
C           polated.  F(1+(I-1)*INCFD) is value corresponding to X(I).
C
C     D -- (output) real array of derivative values at the data points.
C           These values will determine the cubic spline interpolant
C           with the requested boundary conditions.
C           The value corresponding to X(I) is stored in
C                D(1+(I-1)*INCFD),  I=1(1)N.
C           No other entries in D are changed.
C
C     INCFD -- (input) increment between successive values in F and D.
C           This argument is provided primarily for 2-D applications.
C           (Error return if  INCFD.LT.1 .)
C
C     WK -- (scratch) real array of working storage.
C
C     NWK -- (input) length of work array.
C           (Error return if NWK.LT.2*N .)
C
C     IERR -- (output) error flag.
C           Normal return:
C              IERR = 0  (no errors).
C           "Recoverable" errors:
C              IERR = -1  if N.LT.2 .
C              IERR = -2  if INCFD.LT.1 .
C              IERR = -3  if the X-array is not strictly increasing.
C              IERR = -4  if IBEG.LT.0 or IBEG.GT.4 .
C              IERR = -5  if IEND.LT.0 of IEND.GT.4 .
C              IERR = -6  if both of the above are true.
C              IERR = -7  if NWK is too small.
C               NOTE:  The above errors are checked in the order listed,
C                   and following arguments have **NOT** been validated.
C             (The D-array has not been changed in any of these cases.)
C              IERR = -8  in case of trouble solving the linear system
C                         for the interior derivative values.
C             (The D-array may have been changed in this case.)
C             (             Do **NOT** use it!                )
C
C***REFERENCES  CARL DE BOOR, A PRACTICAL GUIDE TO SPLINES, SPRINGER-
C                 VERLAG (NEW YORK, 1978), PP. 53-59.
C***ROUTINES CALLED  PCHDF,XERROR
C***END PROLOGUE  PCHSP
C
C ----------------------------------------------------------------------
C
C  Change record:
C     82-08-04   Converted to SLATEC library version.
C     87-07-07   Minor cosmetic changes to prologue.
C
C ----------------------------------------------------------------------
C
C  Programming notes:
C
C     To produce a double precision version, simply:
C        a. Change PCHSP to DPCHSP wherever it occurs,
C        b. Change the real declarations to double precision, and
C        c. Change the constants ZERO, HALF, ... to double precision.
C
C  DECLARE ARGUMENTS.
C
      INTEGER  IC(2), N, INCFD, NWK, IERR
      REAL  VC(2), X(N), F(INCFD,N), D(INCFD,N), WK(2,N)
C
C  DECLARE LOCAL VARIABLES.
C
      INTEGER  IBEG, IEND, INDEX, J, NM1
      REAL  G, HALF, ONE, STEMP(3), THREE, TWO, XTEMP(4), ZERO
      REAL  PCHDF
C
      DATA  ZERO /0./,  HALF /0.5/,  ONE /1./,  TWO /2./,  THREE /3./
C
C  VALIDITY-CHECK ARGUMENTS.
C
C***FIRST EXECUTABLE STATEMENT  PCHSP
      IF ( N.LT.2 )  GO TO 5001
      IF ( INCFD.LT.1 )  GO TO 5002
      DO 1  J = 2, N
         IF ( X(J).LE.X(J-1) )  GO TO 5003
    1 CONTINUE
C
      IBEG = IC(1)
      IEND = IC(2)
      IERR = 0
      IF ( (IBEG.LT.0).OR.(IBEG.GT.4) )  IERR = IERR - 1
      IF ( (IEND.LT.0).OR.(IEND.GT.4) )  IERR = IERR - 2
      IF ( IERR.LT.0 )  GO TO 5004
C
C  FUNCTION DEFINITION IS OK -- GO ON.
C
      IF ( NWK .LT. 2*N )  GO TO 5007
C
C  COMPUTE FIRST DIFFERENCES OF X SEQUENCE AND STORE IN WK(1,.). ALSO,
C  COMPUTE FIRST DIVIDED DIFFERENCE OF DATA AND STORE IN WK(2,.).
      DO 5  J=2,N
         WK(1,J) = X(J) - X(J-1)
         WK(2,J) = (F(1,J) - F(1,J-1))/WK(1,J)
    5 CONTINUE
C
C  SET TO DEFAULT BOUNDARY CONDITIONS IF N IS TOO SMALL.
C
      IF ( IBEG.GT.N )  IBEG = 0
      IF ( IEND.GT.N )  IEND = 0
C
C  SET UP FOR BOUNDARY CONDITIONS.
C
      IF ( (IBEG.EQ.1).OR.(IBEG.EQ.2) )  THEN
         D(1,1) = VC(1)
      ELSE IF (IBEG .GT. 2)  THEN
C        PICK UP FIRST IBEG POINTS, IN REVERSE ORDER.
         DO 10  J = 1, IBEG
            INDEX = IBEG-J+1
C           INDEX RUNS FROM IBEG DOWN TO 1.
            XTEMP(J) = X(INDEX)
            IF (J .LT. IBEG)  STEMP(J) = WK(2,INDEX)
   10    CONTINUE
C                 --------------------------------
         D(1,1) = PCHDF (IBEG, XTEMP, STEMP, IERR)
C                 --------------------------------
         IF (IERR .NE. 0)  GO TO 5009
         IBEG = 1
      ENDIF
C
      IF ( (IEND.EQ.1).OR.(IEND.EQ.2) )  THEN
         D(1,N) = VC(2)
      ELSE IF (IEND .GT. 2)  THEN
C        PICK UP LAST IEND POINTS.
         DO 15  J = 1, IEND
            INDEX = N-IEND+J
C           INDEX RUNS FROM N+1-IEND UP TO N.
            XTEMP(J) = X(INDEX)
            IF (J .LT. IEND)  STEMP(J) = WK(2,INDEX+1)
   15    CONTINUE
C                 --------------------------------
         D(1,N) = PCHDF (IEND, XTEMP, STEMP, IERR)
C                 --------------------------------
         IF (IERR .NE. 0)  GO TO 5009
         IEND = 1
      ENDIF
C
C --------------------( BEGIN CODING FROM CUBSPL )--------------------
C
C  **** A TRIDIAGONAL LINEAR SYSTEM FOR THE UNKNOWN SLOPES S(J) OF
C  F  AT X(J), J=1,...,N, IS GENERATED AND THEN SOLVED BY GAUSS ELIM-
C  INATION, WITH S(J) ENDING UP IN D(1,J), ALL J.
C     WK(1,.) AND WK(2,.) ARE USED FOR TEMPORARY STORAGE.
C
C  CONSTRUCT FIRST EQUATION FROM FIRST BOUNDARY CONDITION, OF THE FORM
C             WK(2,1)*S(1) + WK(1,1)*S(2) = D(1,1)
C
      IF (IBEG .EQ. 0)  THEN
         IF (N .EQ. 2)  THEN
C           NO CONDITION AT LEFT END AND N = 2.
            WK(2,1) = ONE
            WK(1,1) = ONE
            D(1,1) = TWO*WK(2,2)
         ELSE
C           NOT-A-KNOT CONDITION AT LEFT END AND N .GT. 2.
            WK(2,1) = WK(1,3)
            WK(1,1) = WK(1,2) + WK(1,3)
            D(1,1) =((WK(1,2) + TWO*WK(1,1))*WK(2,2)*WK(1,3)
     *                        + WK(1,2)**2*WK(2,3)) / WK(1,1)
         ENDIF
      ELSE IF (IBEG .EQ. 1)  THEN
C        SLOPE PRESCRIBED AT LEFT END.
         WK(2,1) = ONE
         WK(1,1) = ZERO
      ELSE
C        SECOND DERIVATIVE PRESCRIBED AT LEFT END.
         WK(2,1) = TWO
         WK(1,1) = ONE
         D(1,1) = THREE*WK(2,2) - HALF*WK(1,2)*D(1,1)
      ENDIF
C
C  IF THERE ARE INTERIOR KNOTS, GENERATE THE CORRESPONDING EQUATIONS AND
C  CARRY OUT THE FORWARD PASS OF GAUSS ELIMINATION, AFTER WHICH THE J-TH
C  EQUATION READS    WK(2,J)*S(J) + WK(1,J)*S(J+1) = D(1,J).
C
      NM1 = N-1
      IF (NM1 .GT. 1)  THEN
         DO 20 J=2,NM1
            IF (WK(2,J-1) .EQ. ZERO)  GO TO 5008
            G = -WK(1,J+1)/WK(2,J-1)
            D(1,J) = G*D(1,J-1)
     *                  + THREE*(WK(1,J)*WK(2,J+1) + WK(1,J+1)*WK(2,J))
            WK(2,J) = G*WK(1,J-1) + TWO*(WK(1,J) + WK(1,J+1))
   20    CONTINUE
      ENDIF
C
C  CONSTRUCT LAST EQUATION FROM SECOND BOUNDARY CONDITION, OF THE FORM
C           (-G*WK(2,N-1))*S(N-1) + WK(2,N)*S(N) = D(1,N)
C
C     IF SLOPE IS PRESCRIBED AT RIGHT END, ONE CAN GO DIRECTLY TO BACK-
C     SUBSTITUTION, SINCE ARRAYS HAPPEN TO BE SET UP JUST RIGHT FOR IT
C     AT THIS POINT.
      IF (IEND .EQ. 1)  GO TO 30
C
      IF (IEND .EQ. 0)  THEN
         IF (N.EQ.2 .AND. IBEG.EQ.0)  THEN
C           NOT-A-KNOT AT RIGHT ENDPOINT AND AT LEFT ENDPOINT AND N = 2.
            D(1,2) = WK(2,2)
            GO TO 30
         ELSE IF ((N.EQ.2) .OR. (N.EQ.3 .AND. IBEG.EQ.0))  THEN
C           EITHER (N=3 AND NOT-A-KNOT ALSO AT LEFT) OR (N=2 AND *NOT*
C           NOT-A-KNOT AT LEFT END POINT).
            D(1,N) = TWO*WK(2,N)
            WK(2,N) = ONE
            IF (WK(2,N-1) .EQ. ZERO)  GO TO 5008
            G = -ONE/WK(2,N-1)
         ELSE
C           NOT-A-KNOT AND N .GE. 3, AND EITHER N.GT.3 OR  ALSO NOT-A-
C           KNOT AT LEFT END POINT.
            G = WK(1,N-1) + WK(1,N)
C           DO NOT NEED TO CHECK FOLLOWING DENOMINATORS (X-DIFFERENCES).
            D(1,N) = ((WK(1,N)+TWO*G)*WK(2,N)*WK(1,N-1)
     *                  + WK(1,N)**2*(F(1,N-1)-F(1,N-2))/WK(1,N-1))/G
            IF (WK(2,N-1) .EQ. ZERO)  GO TO 5008
            G = -G/WK(2,N-1)
            WK(2,N) = WK(1,N-1)
         ENDIF
      ELSE
C        SECOND DERIVATIVE PRESCRIBED AT RIGHT ENDPOINT.
         D(1,N) = THREE*WK(2,N) + HALF*WK(1,N)*D(1,N)
         WK(2,N) = TWO
         IF (WK(2,N-1) .EQ. ZERO)  GO TO 5008
         G = -ONE/WK(2,N-1)
      ENDIF
C
C  COMPLETE FORWARD PASS OF GAUSS ELIMINATION.
C
      WK(2,N) = G*WK(1,N-1) + WK(2,N)
      IF (WK(2,N) .EQ. ZERO)   GO TO 5008
      D(1,N) = (G*D(1,N-1) + D(1,N))/WK(2,N)
C
C  CARRY OUT BACK SUBSTITUTION
C
   30 CONTINUE
      DO 40 J=NM1,1,-1
         IF (WK(2,J) .EQ. ZERO)  GO TO 5008
         D(1,J) = (D(1,J) - WK(1,J)*D(1,J+1))/WK(2,J)
   40 CONTINUE
C --------------------(  END  CODING FROM CUBSPL )--------------------
C
C  NORMAL RETURN.
C
      RETURN
C
C  ERROR RETURNS.
C
 5001 CONTINUE
C     N.LT.2 RETURN.
      IERR = -1
      CALL XERROR ('PCHSP -- NUMBER OF DATA POINTS LESS THAN TWO'
     *           , 44, IERR, 1)
      RETURN
C
 5002 CONTINUE
C     INCFD.LT.1 RETURN.
      IERR = -2
      CALL XERROR ('PCHSP -- INCREMENT LESS THAN ONE'
     *           , 32, IERR, 1)
      RETURN
C
 5003 CONTINUE
C     X-ARRAY NOT STRICTLY INCREASING.
      IERR = -3
      CALL XERROR ('PCHSP -- X-ARRAY NOT STRICTLY INCREASING'
     *           , 40, IERR, 1)
      RETURN
C
 5004 CONTINUE
C     IC OUT OF RANGE RETURN.
      IERR = IERR - 3
      CALL XERROR ('PCHSP -- IC OUT OF RANGE'
     *           , 24, IERR, 1)
      RETURN
C
 5007 CONTINUE
C     NWK TOO SMALL RETURN.
      IERR = -7
      CALL XERROR ('PCHSP -- WORK ARRAY TOO SMALL'
     *           , 29, IERR, 1)
      RETURN
C
 5008 CONTINUE
C     SINGULAR SYSTEM.
C   *** THEORETICALLY, THIS CAN ONLY OCCUR IF SUCCESSIVE X-VALUES   ***
C   *** ARE EQUAL, WHICH SHOULD ALREADY HAVE BEEN CAUGHT (IERR=-3). ***
      IERR = -8
      CALL XERROR ('PCHSP -- SINGULAR LINEAR SYSTEM'
     *           , 31, IERR, 1)
      RETURN
C
 5009 CONTINUE
C     ERROR RETURN FROM PCHDF.
C   *** THIS CASE SHOULD NEVER OCCUR ***
      IERR = -9
      CALL XERROR ('PCHSP -- ERROR RETURN FROM PCHDF'
     *           , 32, IERR, 1)
      RETURN
C------------- LAST LINE OF PCHSP FOLLOWS ------------------------------
      END

      REAL FUNCTION PCHST(ARG1,ARG2)
C***BEGIN PROLOGUE  PCHST
C***REFER TO  PCHCE,PCHCI,PCHCS,PCHIM
C***ROUTINES CALLED  (NONE)
C***DESCRIPTION
C
C         PCHST:  PCHIP Sign-Testing Routine.
C
C
C     Returns:
C        -1. if ARG1 and ARG2 are of opposite sign.
C         0. if either argument is zero.
C        +1. if ARG1 and ARG2 are of the same sign.
C
C     The object is to do this without multiplying ARG1*ARG2, to avoid
C     possible over/underflow problems.
C
C  Fortran intrinsics used:  SIGN.
C
C ----------------------------------------------------------------------
C
C  Programmed by:  Fred N. Fritsch,  FTS 532-4275, (415) 422-4275,
C                  Mathematics and Statistics Division,
C                  Lawrence Livermore National Laboratory.
C
C  Change record:
C     82-08-05   Converted to SLATEC library version.
C
C ----------------------------------------------------------------------
C
C  Programming notes:
C
C     To produce a double precision version, simply:
C        a. Change PCHST to DPCHST wherever it occurs,
C        b. Change all references to the Fortran intrinsics to their
C           double presision equivalents,
C        c. Change the real declarations to double precision, and
C        d. Change the constants  ZERO  and  ONE  to double precision.
C***END PROLOGUE  PCHST
      REAL  ARG1, ARG2
C
C  DECLARE LOCAL VARIABLES.
C
      REAL  ONE, ZERO
      DATA  ZERO /0./,  ONE /1./
C
C  PERFORM THE TEST.
C
C***FIRST EXECUTABLE STATEMENT  PCHST
      PCHST = SIGN(ONE,ARG1) * SIGN(ONE,ARG2)
      IF ((ARG1.EQ.ZERO) .OR. (ARG2.EQ.ZERO))  PCHST = ZERO
C
      RETURN
C------------- LAST LINE OF PCHST FOLLOWS ------------------------------
      END

      REAL FUNCTION PCHDF(K,X,S,IERR)
C***BEGIN PROLOGUE  PCHDF
C***REFER TO  PCHCE,PCHSP
C***ROUTINES CALLED  XERROR
C***DESCRIPTION
C
C          PCHDF:   PCHIP Finite Difference Formula
C
C     Uses a divided difference formulation to compute a K-point approx-
C     imation to the derivative at X(K) based on the data in X and S.
C
C     Called by  PCHCE  and  PCHSP  to compute 3- and 4-point boundary
C     derivative approximations.
C
C ----------------------------------------------------------------------
C
C     On input:
C        K      is the order of the desired derivative approximation.
C               K must be at least 3 (error return if not).
C        X      contains the K values of the independent variable.
C               X need not be ordered, but the values **MUST** be
C               distinct.  (Not checked here.)
C        S      contains the associated slope values:
C                  S(I) = (F(I+1)-F(I))/(X(I+1)-X(I)), I=1(1)K-1.
C               (Note that S need only be of length K-1.)
C
C     On return:
C        S      will be destroyed.
C        IERR   will be set to -1 if K.LT.2 .
C        PCHDF  will be set to the desired derivative approximation if
C               IERR=0 or to zero if IERR=-1.
C
C ----------------------------------------------------------------------
C
C  Reference:  Carl de Boor, A Practical Guide to Splines, Springer-
C              Verlag (New York, 1978), pp. 10-16.
C
C ----------------------------------------------------------------------
C
C  Programmed by:  Fred N. Fritsch,  FTS 532-4275, (415) 422-4275,
C                  Mathematics and Statistics Division,
C                  Lawrence Livermore National Laboratory.
C
C  Change record:
C     82-08-05   Converted to SLATEC library version.
C
C ----------------------------------------------------------------------
C
C  Programming notes:
C
C     To produce a double precision version, simply:
C        a. Change PCHDF to DPCHDF wherever it occurs,
C        b. Change the real declarations to double precision, and
C        c. Change the constant ZERO to double precision.
C***END PROLOGUE  PCHDF
      INTEGER  K, IERR
      REAL  X(K), S(K)
C
C  DECLARE LOCAL VARIABLES.
C
      INTEGER  I, J
      REAL  VALUE, ZERO
      DATA  ZERO /0./
C
C  CHECK FOR LEGAL VALUE OF K.
C
C***FIRST EXECUTABLE STATEMENT  PCHDF
      IF (K .LT. 3)  GO TO 5001
C
C  COMPUTE COEFFICIENTS OF INTERPOLATING POLYNOMIAL.
C
      DO 10  J = 2, K-1
         DO 9  I = 1, K-J
            S(I) = (S(I+1)-S(I))/(X(I+J)-X(I))
    9    CONTINUE
   10 CONTINUE
C
C  EVALUATE DERIVATIVE AT X(K).
C
      VALUE = S(1)
      DO 20  I = 2, K-1
         VALUE = S(I) + VALUE*(X(K)-X(I))
   20 CONTINUE
C
C  NORMAL RETURN.
C
      IERR = 0
      PCHDF = VALUE
      RETURN
C
C  ERROR RETURN.
C
 5001 CONTINUE
C     K.LT.3 RETURN.
      IERR = -1
      CALL XERROR ('PCHDF -- K LESS THAN THREE'
     *           , 26, IERR, 1)
      PCHDF = ZERO
      RETURN
C------------- LAST LINE OF PCHDF FOLLOWS ------------------------------
      END

      SUBROUTINE PCHFD(N,X,F,D,INCFD,SKIP,NE,XE,FE,DE,IERR)
C***BEGIN PROLOGUE  PCHFD
C***DATE WRITTEN   811020   (YYMMDD)
C***REVISION DATE  870707   (YYMMDD)
C***CATEGORY NO.  E3,H1
C***KEYWORDS  LIBRARY=SLATEC(PCHIP),
C             TYPE=SINGLE PRECISION(PCHFD-S DPCHFD-D),
C             CUBIC HERMITE DIFFERENTIATION,CUBIC HERMITE EVALUATION,
C             HERMITE INTERPOLATION,PIECEWISE CUBIC EVALUATION
C***AUTHOR  FRITSCH, F. N., (LLNL)
C             MATHEMATICS AND STATISTICS DIVISION
C             LAWRENCE LIVERMORE NATIONAL LABORATORY
C             P.O. BOX 808  (L-316)
C             LIVERMORE, CA  94550
C             FTS 532-4275, (415) 422-4275
C***PURPOSE  Evaluate a piecewise cubic hermite function and its first
C            derivative at an array of points.  May be used by itself
C            for Hermite interpolation, or as an evaluator for PCHIM
C            or PCHIC.  If only function values are required, use
C            PCHFE instead.
C***DESCRIPTION
C
C          PCHFD:  Piecewise Cubic Hermite Function and Derivative
C                  evaluator
C
C     Evaluates the cubic Hermite function defined by  N, X, F, D,  to-
C     gether with its first derivative, at the points  XE(J), J=1(1)NE.
C
C     If only function values are required, use PCHFE, instead.
C
C     To provide compatibility with PCHIM and PCHIC, includes an
C     increment between successive values of the F- and D-arrays.
C
C ----------------------------------------------------------------------
C
C  Calling sequence:
C
C        PARAMETER  (INCFD = ...)
C        INTEGER  N, NE, IERR
C        REAL  X(N), F(INCFD,N), D(INCFD,N), XE(NE), FE(NE), DE(NE)
C        LOGICAL  SKIP
C
C        CALL  PCHFD (N, X, F, D, INCFD, SKIP, NE, XE, FE, DE, IERR)
C
C   Parameters:
C
C     N -- (input) number of data points.  (Error return if N.LT.2 .)
C
C     X -- (input) real array of independent variable values.  The
C           elements of X must be strictly increasing:
C                X(I-1) .LT. X(I),  I = 2(1)N.
C           (Error return if not.)
C
C     F -- (input) real array of function values.  F(1+(I-1)*INCFD) is
C           the value corresponding to X(I).
C
C     D -- (input) real array of derivative values.  D(1+(I-1)*INCFD) is
C           the value corresponding to X(I).
C
C     INCFD -- (input) increment between successive values in F and D.
C           (Error return if  INCFD.LT.1 .)
C
C     SKIP -- (input/output) logical variable which should be set to
C           .TRUE. if the user wishes to skip checks for validity of
C           preceding parameters, or to .FALSE. otherwise.
C           This will save time in case these checks have already
C           been performed (say, in PCHIM or PCHIC).
C           SKIP will be set to .TRUE. on normal return.
C
C     NE -- (input) number of evaluation points.  (Error return if
C           NE.LT.1 .)
C
C     XE -- (input) real array of points at which the functions are to
C           be evaluated.
C
C
C          NOTES:
C           1. The evaluation will be most efficient if the elements
C              of XE are increasing relative to X;
C              that is,   XE(J) .GE. X(I)
C              implies    XE(K) .GE. X(I),  all K.GE.J .
C           2. If any of the XE are outside the interval [X(1),X(N)],
C              values are extrapolated from the nearest extreme cubic,
C              and a warning error is returned.
C
C     FE -- (output) real array of values of the cubic Hermite function
C           defined by  N, X, F, D  at the points  XE.
C
C     DE -- (output) real array of values of the first derivative of
C           the same function at the points  XE.
C
C     IERR -- (output) error flag.
C           Normal return:
C              IERR = 0  (no errors).
C           Warning error:
C              IERR.GT.0  means that extrapolation was performed at
C                 IERR points.
C           "Recoverable" errors:
C              IERR = -1  if N.LT.2 .
C              IERR = -2  if INCFD.LT.1 .
C              IERR = -3  if the X-array is not strictly increasing.
C              IERR = -4  if NE.LT.1 .
C           (Output arrays have not been changed in any of these cases.)
C               NOTE:  The above errors are checked in the order listed,
C                   and following arguments have **NOT** been validated.
C              IERR = -5  if an error has occurred in the lower-level
C                         routine CHFDV.  NB: this should never happen.
C                         Notify the author **IMMEDIATELY** if it does.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  CHFDV,XERROR
C***END PROLOGUE  PCHFD
C
C ----------------------------------------------------------------------
C
C  Change record:
C     82-08-03   Minor cosmetic changes for release 1.
C     87-07-07   Minor cosmetic changes to prologue.
C
C ----------------------------------------------------------------------
C
C  Programming notes:
C
C     1. To produce a double precision version, simply:
C        a. Change PCHFD to DPCHFD, and CHFDV to DCHFDV, wherever they
C           occur,
C        b. Change the real declaration to double precision,
C
C     2. Most of the coding between the call to CHFDV and the end of
C        the IR-loop could be eliminated if it were permissible to
C        assume that XE is ordered relative to X.
C
C     3. CHFDV does not assume that X1 is less than X2.  thus, it would
C        be possible to write a version of PCHFD that assumes a strict-
C        ly decreasing X-array by simply running the IR-loop backwards
C        (and reversing the order of appropriate tests).
C
C     4. The present code has a minor bug, which I have decided is not
C        worth the effort that would be required to fix it.
C        If XE contains points in [X(N-1),X(N)], followed by points .LT.
C        X(N-1), followed by points .GT.X(N), the extrapolation points
C        will be counted (at least) twice in the total returned in IERR.
C
C  DECLARE ARGUMENTS.
C
      INTEGER  N, INCFD, NE, IERR
      REAL  X(N), F(INCFD,N), D(INCFD,N), XE(NE), FE(NE), DE(NE)
      LOGICAL  SKIP
C
C  DECLARE LOCAL VARIABLES.
C
      INTEGER  I, IERC, IR, J, JFIRST, NEXT(2), NJ
C
C  VALIDITY-CHECK ARGUMENTS.
C
C***FIRST EXECUTABLE STATEMENT  PCHFD
      IF (SKIP)  GO TO 5
C
      IF ( N.LT.2 )  GO TO 5001
      IF ( INCFD.LT.1 )  GO TO 5002
      DO 1  I = 2, N
         IF ( X(I).LE.X(I-1) )  GO TO 5003
    1 CONTINUE
C
C  FUNCTION DEFINITION IS OK, GO ON.
C
    5 CONTINUE
      IF ( NE.LT.1 )  GO TO 5004
      IERR = 0
      SKIP = .TRUE.
C
C  LOOP OVER INTERVALS.        (   INTERVAL INDEX IS  IL = IR-1  . )
C                              ( INTERVAL IS X(IL).LE.X.LT.X(IR) . )
      JFIRST = 1
      IR = 2
   10 CONTINUE
C
C     SKIP OUT OF LOOP IF HAVE PROCESSED ALL EVALUATION POINTS.
C
         IF (JFIRST .GT. NE)  GO TO 5000
C
C     LOCATE ALL POINTS IN INTERVAL.
C
         DO 20  J = JFIRST, NE
            IF (XE(J) .GE. X(IR))  GO TO 30
   20    CONTINUE
         J = NE + 1
         GO TO 40
C
C     HAVE LOCATED FIRST POINT BEYOND INTERVAL.
C
   30    CONTINUE
         IF (IR .EQ. N)  J = NE + 1
C
   40    CONTINUE
         NJ = J - JFIRST
C
C     SKIP EVALUATION IF NO POINTS IN INTERVAL.
C
         IF (NJ .EQ. 0)  GO TO 50
C
C     EVALUATE CUBIC AT XE(I),  I = JFIRST (1) J-1 .
C
C       ----------------------------------------------------------------
        CALL CHFDV (X(IR-1),X(IR), F(1,IR-1),F(1,IR), D(1,IR-1),D(1,IR),
     *              NJ, XE(JFIRST), FE(JFIRST), DE(JFIRST), NEXT, IERC)
C       ----------------------------------------------------------------
         IF (IERC .LT. 0)  GO TO 5005
C
         IF (NEXT(2) .EQ. 0)  GO TO 42
C        IF (NEXT(2) .GT. 0)  THEN
C           IN THE CURRENT SET OF XE-POINTS, THERE ARE NEXT(2) TO THE
C           RIGHT OF X(IR).
C
            IF (IR .LT. N)  GO TO 41
C           IF (IR .EQ. N)  THEN
C              THESE ARE ACTUALLY EXTRAPOLATION POINTS.
               IERR = IERR + NEXT(2)
               GO TO 42
   41       CONTINUE
C           ELSE
C              WE SHOULD NEVER HAVE GOTTEN HERE.
               GO TO 5005
C           ENDIF
C        ENDIF
   42    CONTINUE
C
         IF (NEXT(1) .EQ. 0)  GO TO 49
C        IF (NEXT(1) .GT. 0)  THEN
C           IN THE CURRENT SET OF XE-POINTS, THERE ARE NEXT(1) TO THE
C           LEFT OF X(IR-1).
C
            IF (IR .GT. 2)  GO TO 43
C           IF (IR .EQ. 2)  THEN
C              THESE ARE ACTUALLY EXTRAPOLATION POINTS.
               IERR = IERR + NEXT(1)
               GO TO 49
   43       CONTINUE
C           ELSE
C              XE IS NOT ORDERED RELATIVE TO X, SO MUST ADJUST
C              EVALUATION INTERVAL.
C
C              FIRST, LOCATE FIRST POINT TO LEFT OF X(IR-1).
               DO 44  I = JFIRST, J-1
                  IF (XE(I) .LT. X(IR-1))  GO TO 45
   44          CONTINUE
C              NOTE-- CANNOT DROP THROUGH HERE UNLESS THERE IS AN ERROR
C                     IN CHFDV.
               GO TO 5005
C
   45          CONTINUE
C              RESET J.  (THIS WILL BE THE NEW JFIRST.)
               J = I
C
C              NOW FIND OUT HOW FAR TO BACK UP IN THE X-ARRAY.
               DO 46  I = 1, IR-1
                  IF (XE(J) .LT. X(I)) GO TO 47
   46          CONTINUE
C              NB-- CAN NEVER DROP THROUGH HERE, SINCE XE(J).LT.X(IR-1).
C
   47          CONTINUE
C              AT THIS POINT, EITHER  XE(J) .LT. X(1)
C                 OR      X(I-1) .LE. XE(J) .LT. X(I) .
C              RESET IR, RECOGNIZING THAT IT WILL BE INCREMENTED BEFORE
C              CYCLING.
               IR = MAX0(1, I-1)
C           ENDIF
C        ENDIF
   49    CONTINUE
C
         JFIRST = J
C
C     END OF IR-LOOP.
C
   50 CONTINUE
      IR = IR + 1
      IF (IR .LE. N)  GO TO 10
C
C  NORMAL RETURN.
C
 5000 CONTINUE
      RETURN
C
C  ERROR RETURNS.
C
 5001 CONTINUE
C     N.LT.2 RETURN.
      IERR = -1
      CALL XERROR ('PCHFD -- NUMBER OF DATA POINTS LESS THAN TWO'
     *           , 44, IERR, 1)
      RETURN
C
 5002 CONTINUE
C     INCFD.LT.1 RETURN.
      IERR = -2
      CALL XERROR ('PCHFD -- INCREMENT LESS THAN ONE'
     *           , 32, IERR, 1)
      RETURN
C
 5003 CONTINUE
C     X-ARRAY NOT STRICTLY INCREASING.
      IERR = -3
      CALL XERROR ('PCHFD -- X-ARRAY NOT STRICTLY INCREASING'
     *           , 40, IERR, 1)
      RETURN
C
 5004 CONTINUE
C     NE.LT.1 RETURN.
      IERR = -4
      CALL XERROR ('PCHFD -- NUMBER OF EVALUATION POINTS LESS THAN ONE'
     *           , 50, IERR, 1)
      RETURN
C
 5005 CONTINUE
C     ERROR RETURN FROM CHFDV.
C   *** THIS CASE SHOULD NEVER OCCUR ***
      IERR = -5
      CALL XERROR ('PCHFD -- ERROR RETURN FROM CHFDV -- FATAL'
     *           , 41, IERR, 2)
      RETURN
C------------- LAST LINE OF PCHFD FOLLOWS ------------------------------
      END
      SUBROUTINE CHFDV(X1,X2,F1,F2,D1,D2,NE,XE,FE,DE,NEXT,IERR)
C***BEGIN PROLOGUE  CHFDV
C***DATE WRITTEN   811019   (YYMMDD)
C***REVISION DATE  870707   (YYMMDD)
C***CATEGORY NO.  E3,H1
C***KEYWORDS  LIBRARY=SLATEC(PCHIP),
C             TYPE=SINGLE PRECISION(CHFDV-S DCHFDV-D),
C             CUBIC HERMITE DIFFERENTIATION,CUBIC HERMITE EVALUATION,
C             CUBIC POLYNOMIAL EVALUATION
C***AUTHOR  FRITSCH, F. N., (LLNL)
C             MATHEMATICS AND STATISTICS DIVISION
C             LAWRENCE LIVERMORE NATIONAL LABORATORY
C             P.O. BOX 808  (L-316)
C             LIVERMORE, CA  94550
C             FTS 532-4275, (415) 422-4275
C***PURPOSE  Evaluate a cubic polynomial given in Hermite form and its
C            first derivative at an array of points.  While designed for
C            use by PCHFD, it may be useful directly as an evaluator for
C            a piecewise cubic Hermite function in applications, such as
C            graphing, where the interval is known in advance.
C            If only function values are required, use CHFEV instead.
C***DESCRIPTION
C
C        CHFDV:  Cubic Hermite Function and Derivative Evaluator
C
C     Evaluates the cubic polynomial determined by function values
C     F1,F2 and derivatives D1,D2 on interval (X1,X2), together with
C     its first derivative, at the points  XE(J), J=1(1)NE.
C
C     If only function values are required, use CHFEV, instead.
C
C ----------------------------------------------------------------------
C
C  Calling sequence:
C
C        INTEGER  NE, NEXT(2), IERR
C        REAL  X1, X2, F1, F2, D1, D2, XE(NE), FE(NE), DE(NE)
C
C        CALL  CHFDV (X1,X2, F1,F2, D1,D2, NE, XE, FE, DE, NEXT, IERR)
C
C   Parameters:
C
C     X1,X2 -- (input) endpoints of interval of definition of cubic.
C           (Error return if  X1.EQ.X2 .)
C
C     F1,F2 -- (input) values of function at X1 and X2, respectively.
C
C     D1,D2 -- (input) values of derivative at X1 and X2, respectively.
C
C     NE -- (input) number of evaluation points.  (Error return if
C           NE.LT.1 .)
C
C     XE -- (input) real array of points at which the functions are to
C           be evaluated.  If any of the XE are outside the interval
C           [X1,X2], a warning error is returned in NEXT.
C
C     FE -- (output) real array of values of the cubic function defined
C           by  X1,X2, F1,F2, D1,D2  at the points  XE.
C
C     DE -- (output) real array of values of the first derivative of
C           the same function at the points  XE.
C
C     NEXT -- (output) integer array indicating number of extrapolation
C           points:
C            NEXT(1) = number of evaluation points to left of interval.
C            NEXT(2) = number of evaluation points to right of interval.
C
C     IERR -- (output) error flag.
C           Normal return:
C              IERR = 0  (no errors).
C           "Recoverable" errors:
C              IERR = -1  if NE.LT.1 .
C              IERR = -2  if X1.EQ.X2 .
C                (Output arrays have not been changed in either case.)
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  XERROR
C***END PROLOGUE  CHFDV
C
C ----------------------------------------------------------------------
C
C  Change record:
C     82-08-03   Minor cosmetic changes for release 1.
C
C ----------------------------------------------------------------------
C
C  Programming notes:
C
C     To produce a double precision version, simply:
C        a. Change CHFDV to DCHFDV wherever it occurs,
C        b. Change the real declaration to double precision,
C        c. Change the constant ZERO to double precision, and
C        d. Change the names of the Fortran functions:  AMAX1, AMIN1.
C
C  DECLARE ARGUMENTS.
C
      INTEGER  NE, NEXT(2), IERR
      REAL  X1, X2, F1, F2, D1, D2, XE(NE), FE(NE), DE(NE)
C
C  DECLARE LOCAL VARIABLES.
C
      INTEGER  I
      REAL  C2, C2T2, C3, C3T3, DEL1, DEL2, DELTA, H, X, XMI, XMA, ZERO
      DATA  ZERO /0./
C
C  VALIDITY-CHECK ARGUMENTS.
C
C***FIRST EXECUTABLE STATEMENT  CHFDV
      IF (NE .LT. 1)  GO TO 5001
      H = X2 - X1
      IF (H .EQ. ZERO)  GO TO 5002
C
C  INITIALIZE.
C
      IERR = 0
      NEXT(1) = 0
      NEXT(2) = 0
      XMI = AMIN1(ZERO, H)
      XMA = AMAX1(ZERO, H)
C
C  COMPUTE CUBIC COEFFICIENTS (EXPANDED ABOUT X1).
C
      DELTA = (F2 - F1)/H
      DEL1 = (D1 - DELTA)/H
      DEL2 = (D2 - DELTA)/H
C                                           (DELTA IS NO LONGER NEEDED.)
      C2 = -(DEL1+DEL1 + DEL2)
      C2T2 = C2 + C2
      C3 = (DEL1 + DEL2)/H
C                               (H, DEL1 AND DEL2 ARE NO LONGER NEEDED.)
      C3T3 = C3+C3+C3
C
C  EVALUATION LOOP.
C
      DO 500  I = 1, NE
         X = XE(I) - X1
         FE(I) = F1 + X*(D1 + X*(C2 + X*C3))
         DE(I) = D1 + X*(C2T2 + X*C3T3)
C          COUNT EXTRAPOLATION POINTS.
         IF ( X.LT.XMI )  NEXT(1) = NEXT(1) + 1
         IF ( X.GT.XMA )  NEXT(2) = NEXT(2) + 1
C        (NOTE REDUNDANCY--IF EITHER CONDITION IS TRUE, OTHER IS FALSE.)
  500 CONTINUE
C
C  NORMAL RETURN.
C
      RETURN
C
C  ERROR RETURNS.
C
 5001 CONTINUE
C     NE.LT.1 RETURN.
      IERR = -1
      CALL XERROR ('CHFDV -- NUMBER OF EVALUATION POINTS LESS THAN ONE'
     *           , 50, IERR, 1)
      RETURN
C
 5002 CONTINUE
C     X1.EQ.X2 RETURN.
      IERR = -2
      CALL XERROR ('CHFDV -- INTERVAL ENDPOINTS EQUAL'
     *           , 33, IERR, 1)
      RETURN
C------------- LAST LINE OF CHFDV FOLLOWS ------------------------------
      END

      REAL FUNCTION PCHIA(N,X,F,D,INCFD,SKIP,A,B,IERR)
C***BEGIN PROLOGUE  PCHIA
C***DATE WRITTEN   820730   (YYMMDD)
C***REVISION DATE  870707   (YYMMDD)
C***CATEGORY NO.  E3,H2A2
C***KEYWORDS  LIBRARY=SLATEC(PCHIP),
C             TYPE=SINGLE PRECISION(PCHIA-S DPCHIA-D),
C             CUBIC HERMITE INTERPOLATION,NUMERICAL INTEGRATION,
C             QUADRATURE
C***AUTHOR  FRITSCH, F. N., (LLNL)
C             MATHEMATICS AND STATISTICS DIVISION
C             LAWRENCE LIVERMORE NATIONAL LABORATORY
C             P.O. BOX 808  (L-316)
C             LIVERMORE, CA  94550
C             FTS 532-4275, (415) 422-4275
C***PURPOSE  Evaluate the definite integral of a piecewise cubic
C            Hermite function over an arbitrary interval.
C***DESCRIPTION
C
C          PCHIA:  Piecewise Cubic Hermite Integrator, Arbitrary limits
C
C     Evaluates the definite integral of the cubic Hermite function
C     defined by  N, X, F, D  over the interval [A, B].
C
C     To provide compatibility with PCHIM and PCHIC, includes an
C     increment between successive values of the F- and D-arrays.
C
C ----------------------------------------------------------------------
C
C  Calling sequence:
C
C        PARAMETER  (INCFD = ...)
C        INTEGER  N, IERR
C        REAL  X(N), F(INCFD,N), D(INCFD,N), A, B
C        LOGICAL  SKIP
C
C        VALUE = PCHIA (N, X, F, D, INCFD, SKIP, A, B, IERR)
C
C   Parameters:
C
C     VALUE -- (output) VALUE of the requested integral.
C
C     N -- (input) number of data points.  (Error return if N.LT.2 .)
C
C     X -- (input) real array of independent variable values.  The
C           elements of X must be strictly increasing:
C                X(I-1) .LT. X(I),  I = 2(1)N.
C           (Error return if not.)
C
C     F -- (input) real array of function values.  F(1+(I-1)*INCFD) is
C           the value corresponding to X(I).
C
C     D -- (input) real array of derivative values.  D(1+(I-1)*INCFD) is
C           the value corresponding to X(I).
C
C     INCFD -- (input) increment between successive values in F and D.
C           (Error return if  INCFD.LT.1 .)
C
C     SKIP -- (input/output) logical variable which should be set to
C           .TRUE. if the user wishes to skip checks for validity of
C           preceding parameters, or to .FALSE. otherwise.
C           This will save time in case these checks have already
C           been performed (say, in PCHIM or PCHIC).
C           SKIP will be set to .TRUE. on return with IERR.GE.0 .
C
C     A,B -- (input) the limits of integration.
C           NOTE:  There is no requirement that [A,B] be contained in
C                  [X(1),X(N)].  However, the resulting integral value
C                  will be highly suspect, if not.
C
C     IERR -- (output) error flag.
C           Normal return:
C              IERR = 0  (no errors).
C           Warning errors:
C              IERR = 1  if  A  is outside the interval [X(1),X(N)].
C              IERR = 2  if  B  is outside the interval [X(1),X(N)].
C              IERR = 3  if both of the above are true.  (Note that this
C                        means that either [A,B] contains data interval
C                        or the intervals do not intersect at all.)
C           "Recoverable" errors:
C              IERR = -1  if N.LT.2 .
C              IERR = -2  if INCFD.LT.1 .
C              IERR = -3  if the X-array is not strictly increasing.
C                (Value has not been computed in any of these cases.)
C               NOTE:  The above errors are checked in the order listed,
C                   and following arguments have **NOT** been validated.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  CHFIV,PCHID,XERROR
C***END PROLOGUE  PCHIA
C
C ----------------------------------------------------------------------
C
C
C  Change record:
C     82-08-04   Converted to SLATEC library version.
C     87-07-07   Corrected double precision conversion instructions.
C
C ----------------------------------------------------------------------
C
C  Programming notes:
C
C     To produce a double precision version, simply:
C        a. Change PCHIA to DPCHIA wherever it occurs,
C        b. Change PCHID to DPCHID wherever it occurs,
C        c. Change CHFIV to DCHFIV wherever it occurs,
C        d. Change the real declarations to double precision,  and
C        e. Change the constant  ZERO  to double precision.
C
C  DECLARE ARGUMENTS.
C
      INTEGER  N, INCFD, IERR
      REAL  X(N), F(INCFD,N), D(INCFD,N), A, B
      LOGICAL  SKIP
C
C  DECLARE LOCAL VARIABLES.
C
      INTEGER  I, IA, IB, IERD, IERV, IL, IR
      REAL  VALUE, XA, XB, ZERO
      REAL  CHFIV, PCHID
C
C  INITIALIZE.
C
      DATA  ZERO /0./
C
C  VALIDITY-CHECK ARGUMENTS.
C
C***FIRST EXECUTABLE STATEMENT  PCHIA
      IF (SKIP)  GO TO 5
C
      IF ( N.LT.2 )  GO TO 5001
      IF ( INCFD.LT.1 )  GO TO 5002
      DO 1  I = 2, N
         IF ( X(I).LE.X(I-1) )  GO TO 5003
    1 CONTINUE
C
C  FUNCTION DEFINITION IS OK, GO ON.
C
    5 CONTINUE
      SKIP = .TRUE.
      IERR = 0
      IF ( (A.LT.X(1)) .OR. (A.GT.X(N)) )  IERR = IERR + 1
      IF ( (B.LT.X(1)) .OR. (B.GT.X(N)) )  IERR = IERR + 2
C
C  COMPUTE INTEGRAL VALUE.
C
      IF (A .EQ. B)  THEN
         VALUE = ZERO
      ELSE
         XA = AMIN1 (A, B)
         XB = AMAX1 (A, B)
         IF (XB .LE. X(2))  THEN
C           INTERVAL IS TO LEFT OF X(2), SO USE FIRST CUBIC.
C                   --------------------------------------------
            VALUE = CHFIV (X(1),X(2), F(1,1),F(1,2),
     *                                D(1,1),D(1,2), A, B, IERV)
C                   --------------------------------------------
            IF (IERV .LT. 0)  GO TO 5004
         ELSE IF (XA .GE. X(N-1))  THEN
C           INTERVAL IS TO RIGHT OF X(N-1), SO USE LAST CUBIC.
C                   -----------------------------------------------
            VALUE = CHFIV(X(N-1),X(N), F(1,N-1),F(1,N),
     *                                 D(1,N-1),D(1,N), A, B, IERV)
C                   -----------------------------------------------
            IF (IERV .LT. 0)  GO TO 5004
         ELSE
C           'NORMAL' CASE -- XA.LT.XB, XA.LT.X(N-1), XB.GT.X(2).
C      ......LOCATE IA AND IB SUCH THAT
C               X(IA-1).LT.XA.LE.X(IA).LE.X(IB).LE.XB.LE.X(IB+1)
            IA = 1
            DO 10  I = 1, N-1
               IF (XA .GT. X(I))  IA = I + 1
   10       CONTINUE
C             IA = 1 IMPLIES XA.LT.X(1) .  OTHERWISE,
C             IA IS LARGEST INDEX SUCH THAT X(IA-1).LT.XA,.
C
            IB = N
            DO 20  I = N, IA, -1
               IF (XB .LT. X(I))  IB = I - 1
   20       CONTINUE
C             IB = N IMPLIES XB.GT.X(N) .  OTHERWISE,
C             IB IS SMALLEST INDEX SUCH THAT XB.LT.X(IB+1) .
C
C     ......COMPUTE THE INTEGRAL.
            IERV = 0
            IF (IB .LT. IA)  THEN
C              THIS MEANS IB = IA-1 AND
C                 (A,B) IS A SUBSET OF (X(IB),X(IA)).
C                      ------------------------------------------------
               VALUE = CHFIV (X(IB),X(IA), F(1,IB),F(1,IA),
     *                                     D(1,IB),D(1,IA), A, B, IERV)
C                      ------------------------------------------------
               IF (IERV .LT. 0)  GO TO 5004
            ELSE
C
C              FIRST COMPUTE INTEGRAL OVER (X(IA),X(IB)).
               IF (IB .EQ. IA)  THEN
                  VALUE = ZERO
               ELSE
C                         ---------------------------------------------
                  VALUE = PCHID (N, X, F, D, INCFD, SKIP, IA, IB, IERD)
C                         ---------------------------------------------
                  IF (IERD .LT. 0)  GO TO 5005
               ENDIF
C
C              THEN ADD ON INTEGRAL OVER (XA,X(IA)).
               IF (XA .LT. X(IA))  THEN
                  IL = MAX0 (1, IA-1)
                  IR = IL + 1
C                                 -------------------------------------
                  VALUE = VALUE + CHFIV (X(IL),X(IR), F(1,IL),F(1,IR),
     *                                D(1,IL),D(1,IR), XA, X(IA), IERV)
C                                 -------------------------------------
                  IF (IERV .LT. 0)  GO TO 5004
               ENDIF
C
C              THEN ADD ON INTEGRAL OVER (X(IB),XB).
               IF (XB .GT. X(IB))  THEN
                  IR = MIN0 (IB+1, N)
                  IL = IR - 1
C                                 -------------------------------------
                  VALUE = VALUE + CHFIV (X(IL),X(IR), F(1,IL),F(1,IR),
     *                                D(1,IL),D(1,IR), X(IB), XB, IERV)
C                                 -------------------------------------
                  IF (IERV .LT. 0)  GO TO 5004
               ENDIF
C
C              FINALLY, ADJUST SIGN IF NECESSARY.
               IF (A .GT. B)  VALUE = -VALUE
            ENDIF
         ENDIF
      ENDIF
C
C  NORMAL RETURN.
C
      PCHIA = VALUE
      RETURN
C
C  ERROR RETURNS.
C
 5001 CONTINUE
C     N.LT.2 RETURN.
      IERR = -1
      CALL XERROR ('PCHIA -- NUMBER OF DATA POINTS LESS THAN TWO'
     *           , 44, IERR, 1)
      RETURN
C
 5002 CONTINUE
C     INCFD.LT.1 RETURN.
      IERR = -2
      CALL XERROR ('PCHIA -- INCREMENT LESS THAN ONE'
     *           , 32, IERR, 1)
      RETURN
C
 5003 CONTINUE
C     X-ARRAY NOT STRICTLY INCREASING.
      IERR = -3
      CALL XERROR ('PCHIA -- X-ARRAY NOT STRICTLY INCREASING'
     *           , 40, IERR, 1)
      RETURN
C
 5004 CONTINUE
C     TROUBLE IN CHFIV.  (SHOULD NEVER OCCUR.)
      IERR = -4
      CALL XERROR ('PCHIA -- TROUBLE IN CHFIV'
     *           , 25, IERR, 1)
      RETURN
C
 5005 CONTINUE
C     TROUBLE IN PCHID.  (SHOULD NEVER OCCUR.)
      IERR = -5
      CALL XERROR ('PCHIA -- TROUBLE IN PCHID'
     *           , 25, IERR, 1)
      RETURN
C------------- LAST LINE OF PCHIA FOLLOWS ------------------------------
      END
      REAL FUNCTION PCHID(N,X,F,D,INCFD,SKIP,IA,IB,IERR)
C***BEGIN PROLOGUE  PCHID
C***DATE WRITTEN   820723   (YYMMDD)
C***REVISION DATE  870707   (YYMMDD)
C***CATEGORY NO.  E1B,H2A2
C***KEYWORDS  LIBRARY=SLATEC(PCHIP),
C             TYPE=SINGLE PRECISION(PCHID-S DPCHID-D),
C             CUBIC HERMITE INTERPOLATION,NUMERICAL INTEGRATION,
C             QUADRATURE
C***AUTHOR  FRITSCH, F. N., (LLNL)
C             MATHEMATICS AND STATISTICS DIVISION
C             LAWRENCE LIVERMORE NATIONAL LABORATORY
C             P.O. BOX 808  (L-316)
C             LIVERMORE, CA  94550
C             FTS 532-4275, (415) 422-4275
C***PURPOSE  Evaluate the definite integral of a piecewise cubic
C            Hermite function over an interval whose endpoints are
C            data points.
C***DESCRIPTION
C
C          PCHID:  Piecewise Cubic Hermite Integrator, Data limits
C
C     Evaluates the definite integral of the cubic Hermite function
C     defined by  N, X, F, D  over the interval [X(IA), X(IB)].
C
C     To provide compatibility with PCHIM and PCHIC, includes an
C     increment between successive values of the F- and D-arrays.
C
C ----------------------------------------------------------------------
C
C  Calling sequence:
C
C        PARAMETER  (INCFD = ...)
C        INTEGER  N, IA, IB, IERR
C        REAL  X(N), F(INCFD,N), D(INCFD,N)
C        LOGICAL  SKIP
C
C        VALUE = PCHID (N, X, F, D, INCFD, SKIP, IA, IB, IERR)
C
C   Parameters:
C
C     VALUE -- (output) VALUE of the requested integral.
C
C     N -- (input) number of data points.  (Error return if N.LT.2 .)
C
C     X -- (input) real array of independent variable values.  The
C           elements of X must be strictly increasing:
C                X(I-1) .LT. X(I),  I = 2(1)N.
C           (Error return if not.)
C
C     F -- (input) real array of function values.  F(1+(I-1)*INCFD) is
C           the value corresponding to X(I).
C
C     D -- (input) real array of derivative values.  D(1+(I-1)*INCFD) is
C           the value corresponding to X(I).
C
C     INCFD -- (input) increment between successive values in F and D.
C           (Error return if  INCFD.LT.1 .)
C
C     SKIP -- (input/output) logical variable which should be set to
C           .TRUE. if the user wishes to skip checks for validity of
C           preceding parameters, or to .FALSE. otherwise.
C           This will save time in case these checks have already
C           been performed (say, in PCHIM or PCHIC).
C           SKIP will be set to .TRUE. on return with IERR = 0 or -4.
C
C     IA,IB -- (input) indices in X-array for the limits of integration.
C           both must be in the range [1,N].  (Error return if not.)
C           No restrictions on their relative values.
C
C     IERR -- (output) error flag.
C           Normal return:
C              IERR = 0  (no errors).
C           "Recoverable" errors:
C              IERR = -1  if N.LT.2 .
C              IERR = -2  if INCFD.LT.1 .
C              IERR = -3  if the X-array is not strictly increasing.
C              IERR = -4  if IA or IB is out of range.
C                (Value has not been computed in any of these cases.)
C               NOTE:  The above errors are checked in the order listed,
C                   and following arguments have **NOT** been validated.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  XERROR
C***END PROLOGUE  PCHID
C
C ----------------------------------------------------------------------
C
C  Change record:
C     82-08-04   Converted to SLATEC library version.
C
C ----------------------------------------------------------------------
C
C  Programming notes:
C
C     To produce a double precision version, simply:
C        a. Change PCHID to DPCHID wherever it occurs,
C        b. Change the real declarations to double precision,  and
C        c. Change the constants ZERO, HALF, SIX to double precision.
C
C  DECLARE ARGUMENTS.
C
      INTEGER  N, INCFD, IA, IB, IERR
      REAL  X(N), F(INCFD,N), D(INCFD,N)
      LOGICAL  SKIP
C
C  DECLARE LOCAL VARIABLES.
C
      INTEGER  I, IUP, LOW
      REAL  H, HALF, SIX, SUM, VALUE, ZERO
C
C  INITIALIZE.
C
      DATA  ZERO /0./,  HALF /0.5/,  SIX /6./
C
C  VALIDITY-CHECK ARGUMENTS.
C
C***FIRST EXECUTABLE STATEMENT  PCHID
      IF (SKIP)  GO TO 5
C
      IF ( N.LT.2 )  GO TO 5001
      IF ( INCFD.LT.1 )  GO TO 5002
      DO 1  I = 2, N
         IF ( X(I).LE.X(I-1) )  GO TO 5003
    1 CONTINUE
C
C  FUNCTION DEFINITION IS OK, GO ON.
C
    5 CONTINUE
      SKIP = .TRUE.
      IF ((IA.LT.1) .OR. (IA.GT.N))  GO TO 5004
      IF ((IB.LT.1) .OR. (IB.GT.N))  GO TO 5004
      IERR = 0
C
C  COMPUTE INTEGRAL VALUE.
C
      IF (IA .EQ. IB)  THEN
         VALUE = ZERO
      ELSE
         LOW = MIN0(IA, IB)
         IUP = MAX0(IA, IB) - 1
         SUM = ZERO
         DO 10  I = LOW, IUP
            H = X(I+1) - X(I)
            SUM = SUM + H*( (F(1,I) + F(1,I+1)) +
     *                      (D(1,I) - D(1,I+1))*(H/SIX) )
   10    CONTINUE
         VALUE = HALF * SUM
         IF (IA .GT. IB)  VALUE = -VALUE
      ENDIF
C
C  NORMAL RETURN.
C
      PCHID = VALUE
      RETURN
C
C  ERROR RETURNS.
C
 5001 CONTINUE
C     N.LT.2 RETURN.
      IERR = -1
      CALL XERROR ('PCHID -- NUMBER OF DATA POINTS LESS THAN TWO'
     *           , 44, IERR, 1)
      RETURN
C
 5002 CONTINUE
C     INCFD.LT.1 RETURN.
      IERR = -2
      CALL XERROR ('PCHID -- INCREMENT LESS THAN ONE'
     *           , 32, IERR, 1)
      RETURN
C
 5003 CONTINUE
C     X-ARRAY NOT STRICTLY INCREASING.
      IERR = -3
      CALL XERROR ('PCHID -- X-ARRAY NOT STRICTLY INCREASING'
     *           , 40, IERR, 1)
      RETURN
C
 5004 CONTINUE
C     IA OR IB OUT OF RANGE RETURN.
      IERR = -4
      CALL XERROR ('PCHID -- IA OR IB OUT OF RANGE'
     *           , 30, IERR, 1)
      RETURN
C------------- LAST LINE OF PCHID FOLLOWS ------------------------------
      END
      REAL FUNCTION CHFIV(X1,X2,F1,F2,D1,D2,A,B,IERR)
C***BEGIN PROLOGUE  CHFIV
C***REFER TO  PCHIA
C***ROUTINES CALLED  XERROR
C***REVISION DATE  870707   (YYMMDD)
C***DESCRIPTION
C
C          CHFIV:  Cubic Hermite Function Integral Evaluator.
C
C     Called by  PCHIA  to evaluate the integral of a single cubic (in
C     Hermite form) over an arbitrary interval (A,B).
C
C ----------------------------------------------------------------------
C
C  Calling sequence:
C
C        INTEGER  IERR
C        REAL  X1, X2, F1, F2, D1, D2, A, B
C        REAL  VALUE, CHFIV
C
C        VALUE = CHFIV (X1, X2, F1, F2, D1, D2, A, B, IERR)
C
C   Parameters:
C
C     VALUE -- (output) VALUE of the requested integral.
C
C     X1,X2 -- (input) endpoints if interval of definition of cubic.
C           (Must be distinct.  Error return if not.)
C
C     F1,F2 -- (input) function values at the ends of the interval.
C
C     D1,D2 -- (input) derivative values at the ends of the interval.
C
C     A,B -- (input) endpoints of interval of integration.
C
C     IERR -- (output) error flag.
C           Normal return:
C              IERR = 0 (no errors).
C           "Recoverable errors":
C              IERR = -1  if X1.EQ.X2 .
C                (VALUE has not been set in this case.)
C
C***END PROLOGUE  CHFIV
C
C ----------------------------------------------------------------------
C
C  Programmed by:  Fred N. Fritsch,  FTS 532-4275, (415) 422-4275,
C                  Mathematics and Statistics Division,
C                  Lawrence Livermore National Laboratory.
C
C  Change record:
C     82-08-05   Converted to SLATEC library version.
C
C ----------------------------------------------------------------------
C
C  Programming notes:
C
C     To produce a double precision version, simply:
C        a. Change CHFIV to DCHFIV wherever it occurs,
C        b. Change the real declarations to double precision, and
C        c. Change the constants HALF, TWO, ... to double precision.
C
C  DECLARE ARGUMENTS.
C
      INTEGER  IERR
      REAL  X1, X2, F1, F2, D1, D2, A, B
C
C  DECLARE LOCAL VARIABLES.
C
      REAL  DTERM, FOUR, FTERM, H, HALF, PHIA1, PHIA2, PHIB1, PHIB2,
     *      PSIA1, PSIA2, PSIB1, PSIB2, TA1, TA2, TB1, TB2, THREE, TWO,
     *      UA1, UA2, UB1, UB2
C
C  INITIALIZE.
C
      DATA  HALF /0.5/,  TWO /2./,  THREE /3./,  FOUR /4./,  SIX /6./
C
C  VALIDITY CHECK INPUT.
C
C***FIRST EXECUTABLE STATEMENT  CHFIV
      IF (X1 .EQ. X2)  GO TO 5001
      IERR = 0
C
C  COMPUTE INTEGRAL.
C
      H = X2 - X1
      TA1 = (A - X1) / H
      TA2 = (X2 - A) / H
      TB1 = (B - X1) / H
      TB2 = (X2 - B) / H
C
      UA1 = TA1**3
      PHIA1 = UA1 * (TWO - TA1)
      PSIA1 = UA1 * (THREE*TA1 - FOUR)
      UA2 = TA2**3
      PHIA2 =  UA2 * (TWO - TA2)
      PSIA2 = -UA2 * (THREE*TA2 - FOUR)
C
      UB1 = TB1**3
      PHIB1 = UB1 * (TWO - TB1)
      PSIB1 = UB1 * (THREE*TB1 - FOUR)
      UB2 = TB2**3
      PHIB2 =  UB2 * (TWO - TB2)
      PSIB2 = -UB2 * (THREE*TB2 - FOUR)
C
      FTERM =   F1*(PHIA2 - PHIB2) + F2*(PHIB1 - PHIA1)
      DTERM = ( D1*(PSIA2 - PSIB2) + D2*(PSIB1 - PSIA1) )*(H/SIX)
C
C  RETURN VALUE.
C
      CHFIV = (HALF*H) * (FTERM + DTERM)
      RETURN
C
C  ERROR RETURN.
C
 5001 CONTINUE
      IERR = -1
      CALL XERROR ('CHFIV -- X1 EQUAL TO X2'
     *           , 23, IERR, 1)
      RETURN
C------------- LAST LINE OF CHFIV FOLLOWS ------------------------------
      END
      SUBROUTINE XERROR(MESSG,NMESSG,NERR,LEVEL)
C***BEGIN PROLOGUE  XERROR
C***DATE WRITTEN   790801   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  R3C
C***KEYWORDS  ERROR,XERROR PACKAGE
C***AUTHOR  JONES, R. E., (SNLA)
C***PURPOSE  Processes an error (diagnostic) message.
C***DESCRIPTION
C    From the book "Numerical Methods and Software"
C       by  D. Kahaner, C. Moler, S. Nash
C           Prentice Hall 1988
C     Abstract
C        XERROR processes a diagnostic message, in a manner
C        determined by the value of LEVEL and the current value
C        of the library error control flag, KONTRL.
C        (See subroutine XSETF for details.)
C
C     Description of Parameters
C      --Input--
C        MESSG - the Hollerith message to be processed, containing
C                no more than 72 characters.
C        NMESSG- the actual number of characters in MESSG.
C        NERR  - the error number associated with this message.
C                NERR must not be zero.
C        LEVEL - error category.
C                =2 means this is an unconditionally fatal error.
C                =1 means this is a recoverable error.  (I.e., it is
C                   non-fatal if XSETF has been appropriately called.)
C                =0 means this is a warning message only.
C                =-1 means this is a warning message which is to be
C                   printed at most once, regardless of how many
C                   times this call is executed.
C
C     Examples
C        CALL XERROR('SMOOTH -- NUM WAS ZERO.',23,1,2)
C        CALL XERROR('INTEG  -- LESS THAN FULL ACCURACY ACHIEVED.',
C                    43,2,1)
C        CALL XERROR('ROOTER -- ACTUAL ZERO OF F FOUND BEFORE INTERVAL F
C    1ULLY COLLAPSED.',65,3,0)
C        CALL XERROR('EXP    -- UNDERFLOWS BEING SET TO ZERO.',39,1,-1)
C
C     Latest revision ---  19 MAR 1980
C     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
C***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
C                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
C                 1982.
C***ROUTINES CALLED  XERRWV
C***END PROLOGUE  XERROR
      CHARACTER*(*) MESSG
C***FIRST EXECUTABLE STATEMENT  XERROR
      CALL XERRWV(MESSG,NMESSG,NERR,LEVEL,0,0,0,0,0.,0.)
      RETURN
      END
      SUBROUTINE XERRWV(MESSG,NMESSG,NERR,LEVEL,NI,I1,I2,NR,R1,R2)
C***BEGIN PROLOGUE  XERRWV
C***DATE WRITTEN   800319   (YYMMDD)
C***REVISION DATE  870916   (YYMMDD)
C***CATEGORY NO.  R3C
C***KEYWORDS  ERROR,XERROR PACKAGE
C***AUTHOR  JONES, R. E., (SNLA)
C***PURPOSE  Processes error message allowing 2 integer and two real
C            values to be included in the message.
C***DESCRIPTION
C    From the book "Numerical Methods and Software"
C       by  D. Kahaner, C. Moler, S. Nash
C           Prentice Hall 1988
C     Abstract
C        XERRWV processes a diagnostic message, in a manner
C        determined by the value of LEVEL and the current value
C        of the library error control flag, KONTRL.
C        (See subroutine XSETF for details.)
C        In addition, up to two integer values and two real
C        values may be printed along with the message.
C
C     Description of Parameters
C      --Input--
C        MESSG - the Hollerith message to be processed.
C        NMESSG- the actual number of characters in MESSG.
C        NERR  - the error number associated with this message.
C                NERR must not be zero.
C        LEVEL - error category.
C                =2 means this is an unconditionally fatal error.
C                =1 means this is a recoverable error.  (I.e., it is
C                   non-fatal if XSETF has been appropriately called.)
C                =0 means this is a warning message only.
C                =-1 means this is a warning message which is to be
C                   printed at most once, regardless of how many
C                   times this call is executed.
C        NI    - number of integer values to be printed. (0 to 2)
C        I1    - first integer value.
C        I2    - second integer value.
C        NR    - number of real values to be printed. (0 to 2)
C        R1    - first real value.
C        R2    - second real value.
C
C     Examples
C        CALL XERRWV('SMOOTH -- NUM (=I1) WAS ZERO.',29,1,2,
C    1   1,NUM,0,0,0.,0.)
C        CALL XERRWV('QUADXY -- REQUESTED ERROR (R1) LESS THAN MINIMUM (
C    1R2).,54,77,1,0,0,0,2,ERRREQ,ERRMIN)
C
C     Latest revision ---  16 SEPT 1987
C     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
C***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
C                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
C                 1982.
C***ROUTINES CALLED  FDUMP,I1MACH,J4SAVE,XERABT,XERCTL,XERPRT,XERSAV,
C                    XGETUA
C***END PROLOGUE  XERRWV
      CHARACTER*(*) MESSG
      CHARACTER*20 LFIRST
      CHARACTER*37 FORM
      DIMENSION LUN(5)
C     GET FLAGS
C***FIRST EXECUTABLE STATEMENT  XERRWV
      LKNTRL = J4SAVE(2,0,.FALSE.)
      MAXMES = J4SAVE(4,0,.FALSE.)
C     CHECK FOR VALID INPUT
      IF ((NMESSG.GT.0).AND.(NERR.NE.0).AND.
     1    (LEVEL.GE.(-1)).AND.(LEVEL.LE.2)) GO TO 10
         IF (LKNTRL.GT.0) CALL XERPRT('FATAL ERROR IN...',17)
         CALL XERPRT('XERROR -- INVALID INPUT',23)
         IF (LKNTRL.GT.0) CALL FDUMP
         IF (LKNTRL.GT.0) CALL XERPRT('JOB ABORT DUE TO FATAL ERROR.',
     1  29)
         IF (LKNTRL.GT.0) CALL XERSAV(' ',0,0,0,KDUMMY)
         CALL XERABT('XERROR -- INVALID INPUT',23)
         RETURN
   10 CONTINUE
C     RECORD MESSAGE
      JUNK = J4SAVE(1,NERR,.TRUE.)
      CALL XERSAV(MESSG,NMESSG,NERR,LEVEL,KOUNT)
C     LET USER OVERRIDE
      LFIRST = MESSG
      LMESSG = NMESSG
      LERR = NERR
      LLEVEL = LEVEL
      CALL XERCTL(LFIRST,LMESSG,LERR,LLEVEL,LKNTRL)
C     RESET TO ORIGINAL VALUES
      LMESSG = NMESSG
      LERR = NERR
      LLEVEL = LEVEL
      LKNTRL = MAX0(-2,MIN0(2,LKNTRL))
      MKNTRL = IABS(LKNTRL)
C     DECIDE WHETHER TO PRINT MESSAGE
      IF ((LLEVEL.LT.2).AND.(LKNTRL.EQ.0)) GO TO 100
      IF (((LLEVEL.EQ.(-1)).AND.(KOUNT.GT.MIN0(1,MAXMES)))
     1.OR.((LLEVEL.EQ.0)   .AND.(KOUNT.GT.MAXMES))
     2.OR.((LLEVEL.EQ.1)   .AND.(KOUNT.GT.MAXMES).AND.(MKNTRL.EQ.1))
     3.OR.((LLEVEL.EQ.2)   .AND.(KOUNT.GT.MAX0(1,MAXMES)))) GO TO 100
         IF (LKNTRL.LE.0) GO TO 20
            CALL XERPRT(' ',1)
C           INTRODUCTION
            IF (LLEVEL.EQ.(-1)) CALL XERPRT
     1('WARNING MESSAGE...THIS MESSAGE WILL ONLY BE PRINTED ONCE.',57)
            IF (LLEVEL.EQ.0) CALL XERPRT('WARNING IN...',13)
            IF (LLEVEL.EQ.1) CALL XERPRT
     1      ('RECOVERABLE ERROR IN...',23)
            IF (LLEVEL.EQ.2) CALL XERPRT('FATAL ERROR IN...',17)
   20    CONTINUE
C        MESSAGE
         CALL XERPRT(MESSG,LMESSG)
         CALL XGETUA(LUN,NUNIT)
         ISIZEI = LOG10(FLOAT(I1MACH(9))) + 1.0
         ISIZEF = LOG10(FLOAT(I1MACH(10))**I1MACH(11)) + 1.0
         DO 50 KUNIT=1,NUNIT
            IUNIT = LUN(KUNIT)
C            IF (IUNIT.EQ.0) IUNIT = I1MACH(4)
            DO 22 I=1,MIN(NI,2)
               WRITE (FORM,21) I,ISIZEI
   21          FORMAT ('(11X,21HIN ABOVE MESSAGE, I',I1,'=,I',I2,')   ')
               IF(IUNIT.EQ.0)THEN
                 IF (I.EQ.1) WRITE (*,FORM) I1
                 IF (I.EQ.2) WRITE (*,FORM) I2
               ELSE
                 IF (I.EQ.1) WRITE (IUNIT,FORM) I1
                 IF (I.EQ.2) WRITE (IUNIT,FORM) I2
               ENDIF
   22       CONTINUE
            DO 24 I=1,MIN(NR,2)
               WRITE (FORM,23) I,ISIZEF+10,ISIZEF
   23          FORMAT ('(11X,21HIN ABOVE MESSAGE, R',I1,'=,E',
     1         I2,'.',I2,')')
               IF(IUNIT.EQ.0)THEN
                 IF (I.EQ.1) WRITE (*,FORM) R1
                 IF (I.EQ.2) WRITE (*,FORM) R2
               ELSE
                 IF (I.EQ.1) WRITE (IUNIT,FORM) R1
                 IF (I.EQ.2) WRITE (IUNIT,FORM) R2
               ENDIF
   24       CONTINUE
            IF (LKNTRL.LE.0) GO TO 40
C              ERROR NUMBER
               IF(IUINT.EQ.0)THEN
                 WRITE(*,30) LERR
               ELSE
                 WRITE (IUNIT,30) LERR
               ENDIF
   30          FORMAT (15H ERROR NUMBER =,I10)
   40       CONTINUE
   50    CONTINUE
C        TRACE-BACK
         IF (LKNTRL.GT.0) CALL FDUMP
  100 CONTINUE
      IFATAL = 0
      IF ((LLEVEL.EQ.2).OR.((LLEVEL.EQ.1).AND.(MKNTRL.EQ.2)))
     1IFATAL = 1
C     QUIT HERE IF MESSAGE IS NOT FATAL
      IF (IFATAL.LE.0) RETURN
      IF ((LKNTRL.LE.0).OR.(KOUNT.GT.MAX0(1,MAXMES))) GO TO 120
C        PRINT REASON FOR ABORT
         IF (LLEVEL.EQ.1) CALL XERPRT
     1   ('JOB ABORT DUE TO UNRECOVERED ERROR.',35)
         IF (LLEVEL.EQ.2) CALL XERPRT
     1   ('JOB ABORT DUE TO FATAL ERROR.',29)
C        PRINT ERROR SUMMARY
         CALL XERSAV(' ',-1,0,0,KDUMMY)
  120 CONTINUE
C     ABORT
      IF ((LLEVEL.EQ.2).AND.(KOUNT.GT.MAX0(1,MAXMES))) LMESSG = 0
      CALL XERABT(MESSG,LMESSG)
      RETURN
      END
      SUBROUTINE XERSAV(MESSG,NMESSG,NERR,LEVEL,ICOUNT)
C***BEGIN PROLOGUE  XERSAV
C***DATE WRITTEN   800319   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  Z
C***KEYWORDS  ERROR,XERROR PACKAGE
C***AUTHOR  JONES, R. E., (SNLA)
C***PURPOSE  Records that an error occurred.
C***DESCRIPTION
C    From the book "Numerical Methods and Software"
C       by  D. Kahaner, C. Moler, S. Nash
C           Prentice Hall 1988
C     Abstract
C        Record that this error occurred.
C
C     Description of Parameters
C     --Input--
C       MESSG, NMESSG, NERR, LEVEL are as in XERROR,
C       except that when NMESSG=0 the tables will be
C       dumped and cleared, and when NMESSG is less than zero the
C       tables will be dumped and not cleared.
C     --Output--
C       ICOUNT will be the number of times this message has
C       been seen, or zero if the table has overflowed and
C       does not contain this message specifically.
C       When NMESSG=0, ICOUNT will not be altered.
C
C     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
C     Latest revision ---  19 Mar 1980
C***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
C                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
C                 1982.
C***ROUTINES CALLED  I1MACH,S88FMT,XGETUA
C***END PROLOGUE  XERSAV
      INTEGER LUN(5)
      CHARACTER*(*) MESSG
      CHARACTER*20 MESTAB(10),MES
      DIMENSION NERTAB(10),LEVTAB(10),KOUNT(10)
      SAVE MESTAB,NERTAB,LEVTAB,KOUNT,KOUNTX
C     NEXT TWO DATA STATEMENTS ARE NECESSARY TO PROVIDE A BLANK
C     ERROR TABLE INITIALLY
      DATA KOUNT(1),KOUNT(2),KOUNT(3),KOUNT(4),KOUNT(5),
     1     KOUNT(6),KOUNT(7),KOUNT(8),KOUNT(9),KOUNT(10)
     2     /0,0,0,0,0,0,0,0,0,0/
      DATA KOUNTX/0/
C***FIRST EXECUTABLE STATEMENT  XERSAV
      IF (NMESSG.GT.0) GO TO 80
C     DUMP THE TABLE
         IF (KOUNT(1).EQ.0) RETURN
C        PRINT TO EACH UNIT
         CALL XGETUA(LUN,NUNIT)
         DO 60 KUNIT=1,NUNIT
            IUNIT = LUN(KUNIT)
            IF (IUNIT.EQ.0) IUNIT = I1MACH(4)
C           PRINT TABLE HEADER
            WRITE (IUNIT,10)
   10       FORMAT (32H0          ERROR MESSAGE SUMMARY/
     1      51H MESSAGE START             NERR     LEVEL     COUNT)
C           PRINT BODY OF TABLE
            DO 20 I=1,10
               IF (KOUNT(I).EQ.0) GO TO 30
               WRITE (IUNIT,15) MESTAB(I),NERTAB(I),LEVTAB(I),KOUNT(I)
   15          FORMAT (1X,A20,3I10)
   20       CONTINUE
   30       CONTINUE
C           PRINT NUMBER OF OTHER ERRORS
            IF (KOUNTX.NE.0) WRITE (IUNIT,40) KOUNTX
   40       FORMAT (41H0OTHER ERRORS NOT INDIVIDUALLY TABULATED=,I10)
            WRITE (IUNIT,50)
   50       FORMAT (1X)
   60    CONTINUE
         IF (NMESSG.LT.0) RETURN
C        CLEAR THE ERROR TABLES
         DO 70 I=1,10
   70       KOUNT(I) = 0
         KOUNTX = 0
         RETURN
   80 CONTINUE
C     PROCESS A MESSAGE...
C     SEARCH FOR THIS MESSG, OR ELSE AN EMPTY SLOT FOR THIS MESSG,
C     OR ELSE DETERMINE THAT THE ERROR TABLE IS FULL.
      MES = MESSG
      DO 90 I=1,10
         II = I
         IF (KOUNT(I).EQ.0) GO TO 110
         IF (MES.NE.MESTAB(I)) GO TO 90
         IF (NERR.NE.NERTAB(I)) GO TO 90
         IF (LEVEL.NE.LEVTAB(I)) GO TO 90
         GO TO 100
   90 CONTINUE
C     THREE POSSIBLE CASES...
C     TABLE IS FULL
         KOUNTX = KOUNTX+1
         ICOUNT = 1
         RETURN
C     MESSAGE FOUND IN TABLE
  100    KOUNT(II) = KOUNT(II) + 1
         ICOUNT = KOUNT(II)
         RETURN
C     EMPTY SLOT FOUND FOR NEW MESSAGE
  110    MESTAB(II) = MES
         NERTAB(II) = NERR
         LEVTAB(II) = LEVEL
         KOUNT(II)  = 1
         ICOUNT = 1
         RETURN
      END
      SUBROUTINE FDUMP
C***BEGIN PROLOGUE  FDUMP
C***DATE WRITTEN   790801   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  Z
C***KEYWORDS  ERROR,XERROR PACKAGE
C***AUTHOR  JONES, R. E., (SNLA)
C***PURPOSE  Symbolic dump (should be locally written).
C***DESCRIPTION
C    From the book "Numerical Methods and Software"
C       by  D. Kahaner, C. Moler, S. Nash
C           Prentice Hall 1988
C        ***Note*** Machine Dependent Routine
C        FDUMP is intended to be replaced by a locally written
C        version which produces a symbolic dump.  Failing this,
C        it should be replaced by a version which prints the
C        subprogram nesting list.  Note that this dump must be
C        printed on each of up to five files, as indicated by the
C        XGETUA routine.  See XSETUA and XGETUA for details.
C
C     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
C     Latest revision ---  23 May 1979
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  FDUMP
C***FIRST EXECUTABLE STATEMENT  FDUMP
      RETURN
      END
      SUBROUTINE XGETUA(IUNITA,N)
C***BEGIN PROLOGUE  XGETUA
C***DATE WRITTEN   790801   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  R3C
C***KEYWORDS  ERROR,XERROR PACKAGE
C***AUTHOR  JONES, R. E., (SNLA)
C***PURPOSE  Returns unit number(s) to which error messages are being
C            sent.
C***DESCRIPTION
C    From the book "Numerical Methods and Software"
C       by  D. Kahaner, C. Moler, S. Nash
C           Prentice Hall 1988
C     Abstract
C        XGETUA may be called to determine the unit number or numbers
C        to which error messages are being sent.
C        These unit numbers may have been set by a call to XSETUN,
C        or a call to XSETUA, or may be a default value.
C
C     Description of Parameters
C      --Output--
C        IUNIT - an array of one to five unit numbers, depending
C                on the value of N.  A value of zero refers to the
C                default unit, as defined by the I1MACH machine
C                constant routine.  Only IUNIT(1),...,IUNIT(N) are
C                defined by XGETUA.  The values of IUNIT(N+1),...,
C                IUNIT(5) are not defined (for N .LT. 5) or altered
C                in any way by XGETUA.
C        N     - the number of units to which copies of the
C                error messages are being sent.  N will be in the
C                range from 1 to 5.
C
C     Latest revision ---  19 MAR 1980
C     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
C***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
C                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
C                 1982.
C***ROUTINES CALLED  J4SAVE
C***END PROLOGUE  XGETUA
      DIMENSION IUNITA(5)
C***FIRST EXECUTABLE STATEMENT  XGETUA
      N = J4SAVE(5,0,.FALSE.)
      DO 30 I=1,N
         INDEX = I+4
         IF (I.EQ.1) INDEX = 3
         IUNITA(I) = J4SAVE(INDEX,0,.FALSE.)
   30 CONTINUE
      RETURN
      END
      INTEGER FUNCTION I1MACH(I)
C***BEGIN PROLOGUE  I1MACH
C***DATE WRITTEN   750101   (YYMMDD)
C***REVISION DATE  840405   (YYMMDD)
C***CATEGORY NO.  R1
C***KEYWORDS  MACHINE CONSTANTS
C***AUTHOR  FOX, P. A., (BELL LABS)
C           HALL, A. D., (BELL LABS)
C           SCHRYER, N. L., (BELL LABS)
C***PURPOSE  Returns integer machine dependent constants
C***DESCRIPTION
C
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C   These machine constant routines must be activated for
C   a particular environment.
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C     I1MACH can be used to obtain machine-dependent parameters
C     for the local machine environment.  It is a function
C     subroutine with one (input) argument, and can be called
C     as follows, for example
C
C          K = I1MACH(I)
C
C     where I=1,...,16.  The (output) value of K above is
C     determined by the (input) value of I.  The results for
C     various values of I are discussed below.
C
C  I/O unit numbers.
C    I1MACH( 1) = the standard input unit.
C    I1MACH( 2) = the standard output unit.
C    I1MACH( 3) = the standard punch unit.
C    I1MACH( 4) = the standard error message unit.
C
C  Words.
C    I1MACH( 5) = the number of bits per integer storage unit.
C    I1MACH( 6) = the number of characters per integer storage unit.
C
C  Integers.
C    assume integers are represented in the S-digit, base-A form
C
C               sign ( X(S-1)*A**(S-1) + ... + X(1)*A + X(0) )
C
C               where 0 .LE. X(I) .LT. A for I=0,...,S-1.
C    I1MACH( 7) = A, the base.
C    I1MACH( 8) = S, the number of base-A digits.
C    I1MACH( 9) = A**S - 1, the largest magnitude.
C
C  Floating-Point Numbers.
C    Assume floating-point numbers are represented in the T-digit,
C    base-B form
C               sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
C
C               where 0 .LE. X(I) .LT. B for I=1,...,T,
C               0 .LT. X(1), and EMIN .LE. E .LE. EMAX.
C    I1MACH(10) = B, the base.
C
C  Single-Precision
C    I1MACH(11) = T, the number of base-B digits.
C    I1MACH(12) = EMIN, the smallest exponent E.
C    I1MACH(13) = EMAX, the largest exponent E.
C
C  Double-Precision
C    I1MACH(14) = T, the number of base-B digits.
C    I1MACH(15) = EMIN, the smallest exponent E.
C    I1MACH(16) = EMAX, the largest exponent E.
C
C  To alter this function for a particular environment,
C  the desired set of DATA statements should be activated by
C  removing the C from column 1.  Also, the values of
C  I1MACH(1) - I1MACH(4) should be checked for consistency
C  with the local operating system.
C***REFERENCES  FOX P.A., HALL A.D., SCHRYER N.L.,*FRAMEWORK FOR A
C                 PORTABLE LIBRARY*, ACM TRANSACTIONS ON MATHEMATICAL
C                 SOFTWARE, VOL. 4, NO. 2, JUNE 1978, PP. 177-188.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  I1MACH
C
      INTEGER IMACH(16),OUTPUT
      EQUIVALENCE (IMACH(4),OUTPUT)
C
C
C     MACHINE CONSTANTS FOR THE CDC CYBER 170 SERIES (FTN5).
C
C      DATA IMACH( 1) /    5 /
C      DATA IMACH( 2) /    6 /
C      DATA IMACH( 3) /    7 /
C      DATA IMACH( 4) /    6 /
C      DATA IMACH( 5) /   60 /
C      DATA IMACH( 6) /   10 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   48 /
C      DATA IMACH( 9) / O"00007777777777777777" /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   48 /
C      DATA IMACH(12) / -974 /
C      DATA IMACH(13) / 1070 /
C      DATA IMACH(14) /   96 /
C      DATA IMACH(15) / -927 /
C      DATA IMACH(16) / 1070 /
C
C     MACHINE CONSTANTS FOR THE CDC CYBER 200 SERIES
C
C     DATA IMACH( 1) /      5 /
C     DATA IMACH( 2) /      6 /
C     DATA IMACH( 3) /      7 /
C     DATA IMACH( 4) /      6 /
C     DATA IMACH( 5) /     64 /
C     DATA IMACH( 6) /      8 /
C     DATA IMACH( 7) /      2 /
C     DATA IMACH( 8) /     47 /
C     DATA IMACH( 9) / X'00007FFFFFFFFFFF' /
C     DATA IMACH(10) /      2 /
C     DATA IMACH(11) /     47 /
C     DATA IMACH(12) / -28625 /
C     DATA IMACH(13) /  28718 /
C     DATA IMACH(14) /     94 /
C     DATA IMACH(15) / -28625 /
C     DATA IMACH(16) /  28718 /
C
C
C     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES.
C
C     DATA IMACH( 1) /    5 /
C     DATA IMACH( 2) /    6 /
C     DATA IMACH( 3) /    7 /
C     DATA IMACH( 4) /6LOUTPUT/
C     DATA IMACH( 5) /   60 /
C     DATA IMACH( 6) /   10 /
C     DATA IMACH( 7) /    2 /
C     DATA IMACH( 8) /   48 /
C     DATA IMACH( 9) / 00007777777777777777B /
C     DATA IMACH(10) /    2 /
C     DATA IMACH(11) /   47 /
C     DATA IMACH(12) / -929 /
C     DATA IMACH(13) / 1070 /
C     DATA IMACH(14) /   94 /
C     DATA IMACH(15) / -929 /
C     DATA IMACH(16) / 1069 /
C
C     MACHINE CONSTANTS FOR THE CRAY 1
C
C     DATA IMACH( 1) /   100 /
C     DATA IMACH( 2) /   101 /
C     DATA IMACH( 3) /   102 /
C     DATA IMACH( 4) /   101 /
C     DATA IMACH( 5) /    64 /
C     DATA IMACH( 6) /     8 /
C     DATA IMACH( 7) /     2 /
C     DATA IMACH( 8) /    63 /
C     DATA IMACH( 9) /  777777777777777777777B /
C     DATA IMACH(10) /     2 /
C     DATA IMACH(11) /    47 /
C     DATA IMACH(12) / -8189 /
C     DATA IMACH(13) /  8190 /
C     DATA IMACH(14) /    94 /
C     DATA IMACH(15) / -8099 /
C     DATA IMACH(16) /  8190 /
C
C
C     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
C     THE XEROX SIGMA 5/7/9, THE SEL SYSTEMS 85/86, AND
C     THE PERKIN ELMER (INTERDATA) 7/32.
C
C     DATA IMACH( 1) /   5 /
C     DATA IMACH( 2) /   6 /
C     DATA IMACH( 3) /   7 /
C     DATA IMACH( 4) /   6 /
C     DATA IMACH( 5) /  32 /
C     DATA IMACH( 6) /   4 /
C     DATA IMACH( 7) /  16 /
C     DATA IMACH( 8) /  31 /
C     DATA IMACH( 9) / Z7FFFFFFF /
C     DATA IMACH(10) /  16 /
C     DATA IMACH(11) /   6 /
C     DATA IMACH(12) / -64 /
C     DATA IMACH(13) /  63 /
C     DATA IMACH(14) /  14 /
C     DATA IMACH(15) / -64 /
C     DATA IMACH(16) /  63 /
C
C     MACHINE CONSTANTS FOR THE IBM PC FAMILY (D. KAHANER NBS)
C
      DATA IMACH/5,6,0,6,32,4,2,31,2147483647,2,24,
     * -125,127,53,-1021,1023/
C               NOTE! I1MACH(3) IS NOT WELL DEFINED AND IS SET TO ZERO.
C
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR).
C
C     DATA IMACH( 1) /    5 /
C     DATA IMACH( 2) /    6 /
C     DATA IMACH( 3) /    5 /
C     DATA IMACH( 4) /    6 /
C     DATA IMACH( 5) /   36 /
C     DATA IMACH( 6) /    5 /
C     DATA IMACH( 7) /    2 /
C     DATA IMACH( 8) /   35 /
C     DATA IMACH( 9) / "377777777777 /
C     DATA IMACH(10) /    2 /
C     DATA IMACH(11) /   27 /
C     DATA IMACH(12) / -128 /
C     DATA IMACH(13) /  127 /
C     DATA IMACH(14) /   54 /
C     DATA IMACH(15) / -101 /
C     DATA IMACH(16) /  127 /
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR).
C
C     DATA IMACH( 1) /    5 /
C     DATA IMACH( 2) /    6 /
C     DATA IMACH( 3) /    5 /
C     DATA IMACH( 4) /    6 /
C     DATA IMACH( 5) /   36 /
C     DATA IMACH( 6) /    5 /
C     DATA IMACH( 7) /    2 /
C     DATA IMACH( 8) /   35 /
C     DATA IMACH( 9) / "377777777777 /
C     DATA IMACH(10) /    2 /
C     DATA IMACH(11) /   27 /
C     DATA IMACH(12) / -128 /
C     DATA IMACH(13) /  127 /
C     DATA IMACH(14) /   62 /
C     DATA IMACH(15) / -128 /
C     DATA IMACH(16) /  127 /
C
C
C     MACHINE CONSTANTS FOR THE SUN-3 (INCLUDES THOSE WITH 68881 CHIP,
C       OR WITH FPA BOARD. ALSO INCLUDES SUN-2 WITH SKY BOARD. MAY ALSO
C       WORK WITH SOFTWARE FLOATING POINT ON EITHER SYSTEM.)
C
C     DATA IMACH( 1) /    5 /
C      DATA IMACH( 2) /    6 /
C      DATA IMACH( 3) /    6 /
C      DATA IMACH( 4) /    0 /
C      DATA IMACH( 5) /   32 /
C      DATA IMACH( 6) /    4 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   31 /
C      DATA IMACH( 9) / 2147483647 /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   24 /
C      DATA IMACH(12) / -125 /
C      DATA IMACH(13) /  128 /
C      DATA IMACH(14) /   53 /
C      DATA IMACH(15) / -1021 /
C      DATA IMACH(16) /  1024 /
C
C
C     MACHINE CONSTANTS FOR THE VAX 11/780
C
C     DATA IMACH(1) /    5 /
C     DATA IMACH(2) /    6 /
C     DATA IMACH(3) /    5 /
C     DATA IMACH(4) /    6 /
C     DATA IMACH(5) /   32 /
C     DATA IMACH(6) /    4 /
C     DATA IMACH(7) /    2 /
C     DATA IMACH(8) /   31 /
C     DATA IMACH(9) /2147483647 /
C     DATA IMACH(10)/    2 /
C     DATA IMACH(11)/   24 /
C     DATA IMACH(12)/ -127 /
C     DATA IMACH(13)/  127 /
C     DATA IMACH(14)/   56 /
C     DATA IMACH(15)/ -127 /
C     DATA IMACH(16)/  127 /
C
C***FIRST EXECUTABLE STATEMENT  I1MACH
      IF (I .LT. 1  .OR.  I .GT. 16)
     1   CALL XERROR ( 'I1MACH -- I OUT OF BOUNDS',25,1,2)
C
      I1MACH=IMACH(I)
      RETURN
C
      END
      FUNCTION J4SAVE(IWHICH,IVALUE,ISET)
C***BEGIN PROLOGUE  J4SAVE
C***REFER TO  XERROR
C    From the book "Numerical Methods and Software"
C       by  D. Kahaner, C. Moler, S. Nash
C           Prentice Hall 1988
C     Abstract
C        J4SAVE saves and recalls several global variables needed
C        by the library error handling routines.
C
C     Description of Parameters
C      --Input--
C        IWHICH - Index of item desired.
C                = 1 Refers to current error number.
C                = 2 Refers to current error control flag.
C                 = 3 Refers to current unit number to which error
C                    messages are to be sent.  (0 means use standard.)
C                 = 4 Refers to the maximum number of times any
C                     message is to be printed (as set by XERMAX).
C                 = 5 Refers to the total number of units to which
C                     each error message is to be written.
C                 = 6 Refers to the 2nd unit for error messages
C                 = 7 Refers to the 3rd unit for error messages
C                 = 8 Refers to the 4th unit for error messages
C                 = 9 Refers to the 5th unit for error messages
C        IVALUE - The value to be set for the IWHICH-th parameter,
C                 if ISET is .TRUE. .
C        ISET   - If ISET=.TRUE., the IWHICH-th parameter will BE
C                 given the value, IVALUE.  If ISET=.FALSE., the
C                 IWHICH-th parameter will be unchanged, and IVALUE
C                 is a dummy parameter.
C      --Output--
C        The (old) value of the IWHICH-th parameter will be returned
C        in the function value, J4SAVE.
C
C     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
C    Adapted from Bell Laboratories PORT Library Error Handler
C     Latest revision ---  23 MAY 1979
C***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
C                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
C                 1982.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  J4SAVE
      LOGICAL ISET
      INTEGER IPARAM(9)
      SAVE IPARAM
      DATA IPARAM(1),IPARAM(2),IPARAM(3),IPARAM(4)/0,2,0,10/
      DATA IPARAM(5)/1/
      DATA IPARAM(6),IPARAM(7),IPARAM(8),IPARAM(9)/0,0,0,0/
C***FIRST EXECUTABLE STATEMENT  J4SAVE
      J4SAVE = IPARAM(IWHICH)
      IF (ISET) IPARAM(IWHICH) = IVALUE
      RETURN
      END
      SUBROUTINE XERABT(MESSG,NMESSG)
C***BEGIN PROLOGUE  XERABT
C***DATE WRITTEN   790801   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  R3C
C***KEYWORDS  ERROR,XERROR PACKAGE
C***AUTHOR  JONES, R. E., (SNLA)
C***PURPOSE  Aborts program execution and prints error message.
C***DESCRIPTION
C    From the book "Numerical Methods and Software"
C       by  D. Kahaner, C. Moler, S. Nash
C           Prentice Hall 1988
C     Abstract
C        ***Note*** machine dependent routine
C        XERABT aborts the execution of the program.
C        The error message causing the abort is given in the calling
C        sequence, in case one needs it for printing on a dayfile,
C        for example.
C
C     Description of Parameters
C        MESSG and NMESSG are as in XERROR, except that NMESSG may
C        be zero, in which case no message is being supplied.
C
C     Written by Ron Jones, with SLATEC Common Math Library Subcommittee
C     Latest revision ---  19 MAR 1980
C***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
C                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
C                 1982.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  XERABT
      CHARACTER*(*) MESSG
C***FIRST EXECUTABLE STATEMENT  XERABT
      STOP
      END
      SUBROUTINE XERCTL(MESSG1,NMESSG,NERR,LEVEL,KONTRL)
C***BEGIN PROLOGUE  XERCTL
C***DATE WRITTEN   790801   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  R3C
C***KEYWORDS  ERROR,XERROR PACKAGE
C***AUTHOR  JONES, R. E., (SNLA)
C***PURPOSE  Allows user control over handling of individual errors.
C***DESCRIPTION
C    From the book "Numerical Methods and Software"
C       by  D. Kahaner, C. Moler, S. Nash
C           Prentice Hall 1988
C     Abstract
C        Allows user control over handling of individual errors.
C        Just after each message is recorded, but before it is
C        processed any further (i.e., before it is printed or
C        a decision to abort is made), a call is made to XERCTL.
C        If the user has provided his own version of XERCTL, he
C        can then override the value of KONTROL used in processing
C        this message by redefining its value.
C        KONTRL may be set to any value from -2 to 2.
C        The meanings for KONTRL are the same as in XSETF, except
C        that the value of KONTRL changes only for this message.
C        If KONTRL is set to a value outside the range from -2 to 2,
C        it will be moved back into that range.
C
C     Description of Parameters
C
C      --Input--
C        MESSG1 - the first word (only) of the error message.
C        NMESSG - same as in the call to XERROR or XERRWV.
C        NERR   - same as in the call to XERROR or XERRWV.
C        LEVEL  - same as in the call to XERROR or XERRWV.
C        KONTRL - the current value of the control flag as set
C                 by a call to XSETF.
C
C      --Output--
C        KONTRL - the new value of KONTRL.  If KONTRL is not
C                 defined, it will remain at its original value.
C                 This changed value of control affects only
C                 the current occurrence of the current message.
C***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
C                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
C                 1982.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  XERCTL
      CHARACTER*20 MESSG1
C***FIRST EXECUTABLE STATEMENT  XERCTL
      RETURN
      END
      SUBROUTINE XERPRT(MESSG,NMESSG)
C***BEGIN PROLOGUE  XERPRT
C***DATE WRITTEN   790801   (YYMMDD)
C***REVISION DATE  870916   (YYMMDD)
C***CATEGORY NO.  Z
C***KEYWORDS  ERROR,XERROR PACKAGE
C***AUTHOR  JONES, R. E., (SNLA)
C***PURPOSE  Prints error messages.
C***DESCRIPTION
C    From the book "Numerical Methods and Software"
C       by  D. Kahaner, C. Moler, S. Nash
C           Prentice Hall 1988
C     Abstract
C        Print the Hollerith message in MESSG, of length NMESSG,
C        on each file indicated by XGETUA.
C     Latest revision ---  16 SEPT 1987
C***REFERENCES  JONES R.E., KAHANER D.K., "XERROR, THE SLATEC ERROR-
C                 HANDLING PACKAGE", SAND82-0800, SANDIA LABORATORIES,
C                 1982.
C***ROUTINES CALLED  I1MACH,S88FMT,XGETUA
C***END PROLOGUE  XERPRT
      INTEGER LUN(5)
      CHARACTER*(*) MESSG
C     OBTAIN UNIT NUMBERS AND WRITE LINE TO EACH UNIT
C***FIRST EXECUTABLE STATEMENT  XERPRT
      CALL XGETUA(LUN,NUNIT)
      LENMES = LEN(MESSG)
      DO 20 KUNIT=1,NUNIT
         IUNIT = LUN(KUNIT)
C         IF (IUNIT.EQ.0) IUNIT = I1MACH(4)
         DO 10 ICHAR=1,LENMES,72
            LAST = MIN0(ICHAR+71 , LENMES)
            IF(IUNIT.EQ.0)THEN
              WRITE (*,'(1X,A)') MESSG(ICHAR:LAST)
            ELSE
              WRITE (IUNIT,'(1X,A)') MESSG(ICHAR:LAST)
            ENDIF
   10    CONTINUE
   20 CONTINUE
      RETURN
      END
