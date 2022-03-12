MODULE FortranFunctions
  PRIVATE
  PUBLIC :: PI, ABSO, COSINE, SINE, DIV!, SQR

  REAL*8, PARAMETER :: PI = 3.14159265358979D0

  CONTAINS

  REAL*8 FUNCTION ABSO(B)
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: B

    IF (B .LT. 0.0D0) THEN
      ABSO = -B
    ELSE
      ABSO = B
    END IF
  END FUNCTION ABSO

  REAL*8 FUNCTION COSINE(num)
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: num

    COSINE = SINE(PI/2.0D0 - num)
  END FUNCTION COSINE

! Bhaskara approximation:   temp = (PI-x)*x
!                           sin(x) = (16*temp) / (5*PI*PI - 4*temp)
! Second approximation:     temp = (x/PI)*(x/PI - 1)
!                           sin(x) = (temp/10)*(36*temp - 31)
! Weight average: Bhaskara -> 0.385  Second -> 0.615
!
! sin(x) approximation with weight average: temp = (x/PI)*(x/PI - 1)
!   sin(x) = temp(2.21652(temp - 31/36) - 1.5372/(1.25 + temp))
  REAL*8 FUNCTION SINE(num)
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: num
    REAL*8 :: x
    REAL*8 :: temp

    x = MOD(ABSO(num),PI)

    temp = (x/PI)*(x/PI - 1.0D0)
    SINE = temp*(2.21652D0*(temp - 31.0D0/36.0D0) - 1.5372D0/(1.25D0 + temp))

    !Adjust for negative angles and shift the graph down by 2.6D-5
    IF (num > 0.0D0) THEN
      IF (MOD(num,2.0D0*PI) > PI) THEN
        SINE = -SINE - 2.6D-5
      ELSE
        SINE = SINE + 2.6D-5
      END IF
    ELSE
      IF (MOD(num,2.0D0*PI) > -PI) THEN
        SINE = -SINE - 2.6D-5
      ELSE
        SINE = SINE + 2.6D-5
      END IF
    END IF
  END FUNCTION SINE

!Bhaskara approximation:   temp = (PI-x)*x
!                          sin(x) = (16*temp) / (5*PI*PI - 4*temp)
  REAL*8 FUNCTION SINB(num)
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: num
    REAL*8 :: x
    REAL*8 :: temp

    x = MOD(ABSO(num),PI)

    temp = (PI-x)*x
    SINB = (16.0D0*temp) / (5.0D0*PI*PI - 4.0D0*temp)

    !Adjust for negative angles
    IF (num > 0.0D0) THEN
      IF (MOD(num,2.0D0*PI) > PI) SINB = -SINB
    ELSE
      IF (MOD(num,2.0D0*PI) > -PI) SINB = -SINB
    END IF
  END FUNCTION SINB

  REAL*8 FUNCTION SQR(B) !This function finds the square root of B
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: B
    REAL*8 :: SQA, SQB, SQC, EPS
    INTEGER :: i

    SQR = 1.0D0     !SQR needs a first value to start the iterations
    SQA = -0.1D0    !SQA, SQB & SQC are dummy variables to keep track
    SQB = -0.2D0    !of how SQR evolves through the iterations.
    SQC = -0.3D0

!      S =~(B)^0.5
!  E + S = (B)^0.5
! (E + S)^2 = B
!  E^2 + 2SE + S^2 = B, But E is very small, so E^2=0:
!        E = DIV(B - S^2,2S), this way:
!  S =~ B^0.5 = S + E -> S =~ S + E, this updates S to a new value.
! So by repeating this multiple times we will get B^0.5

    DO i=1, 100
      EPS = DIV(B - SQR*SQR, 2.0D0*SQR) !Get the value of E according to the current S
      SQR = SQR + EPS !Update S to S+E
      SQA = SQB       !Update the values of our dummys
      SQB = SQC
      SQC = SQR
      IF(SQC .EQ. SQB) THEN           !If the last two values are the same then the code stops.
        EXIT
      ELSE IF (SQC .EQ. SQA) THEN     !If the ith and ith-2 value are the same, then the
        SQR = DIV(SQA + SQB, 2.0D0) !code is cycling, so we stop it and average between this
        EXIT                        !values
      END IF
    END DO
  END FUNCTION SQR

  REAL*8 FUNCTION DIV(AA,BB) !This function returns A/B
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: AA, BB

    DIV = AA*RECIDUALIPROCAL(BB)
  END FUNCTION DIV

  REAL*8 FUNCTION RECIDUALIPROCAL(ZZ)
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: ZZ
    REAL*8 :: FAC, RECID, Z

    Z = ZZ

    !stop if 1/0
    IF(Z .EQ. 0.0D0) THEN
      RECID = 1.0D0/Z
      IF (ISNAN(RECID)) STOP 'HI'
      WRITE(*,*) RECID
    END IF

    !ALLOW NEGATIVE NUMBERS
    FAC = 1.0D0
    IF(Z .LT. 0.0D0)THEN
      FAC = -1.0D0
      Z = ABSO(Z)
    END IF

    IF (Z .GT. 1.0D0)THEN                           !IF Z > 1
      CALL GT1DIVIDE(Z - 1.0D0,RECID)
    ELSE IF(Z .GT. 0.0D0 .AND. Z .LT. 1.0D0) THEN   !IF 0 < Z < 1
      CALL LT1DIVIDE(1.0D0 - Z,RECID)
    ELSE                                            !IF Z = 1
      RECID = 1.0D0
    END IF

    RECIDUALIPROCAL = FAC*RECID
  END FUNCTION RECIDUALIPROCAL

  !DIVIDER (1/(1+x))   0<x<1
  REAL*8 FUNCTION DIVIDER(X)
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: X
    REAL*8 :: TWOTHIRDS
    REAL*8 :: Y

    IF (X .LT. 0.0D0) STOP 'X LT 0'
    IF (X .GT. 1.0D0) STOP 'X GT 1'

    TWOTHIRDS = 0.666666666666667D0

    IF (X .GT. 0.5D0) THEN
      Y = TWOTHIRDS*(X - 0.5D0)
      DIVIDER = DIVIDE(Y)*TWOTHIRDS
    ELSE
      DIVIDER = DIVIDE(X)
    END IF
  END FUNCTION DIVIDER

  !DIVIDE (1/(1+x))   0<x<0.5
  REAL*8 FUNCTION DIVIDE(X)
    IMPLICIT NONE

    REAL*8, INTENT(IN) :: X
    INTEGER :: i, N = 50
    REAL*8 :: P, RECID

    IF (X .LT. 0.0D0) STOP 'X LT 0   DIVIDE'
    IF (X .GT. 0.5D0) STOP 'X GT 0.5 DIVIDE'

    RECID = 0.0D0
    P = 1.0D0

    DO i=1, N
      RECID = RECID + P
      P = -P*X
    END DO

    DIVIDE = RECID
  END FUNCTION DIVIDE

  SUBROUTINE GT1DIVIDE(X,RECID)
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: X
    REAL*8, INTENT(OUT) :: RECID
    REAL*8 :: ONED2N(0:64), ARE(0:64), TWON(0:64)
    REAL*8 :: Y
    INTEGER :: i

    ONED2N(0) = 1.0D0
    DO i=1, 64
      ONED2N(i) = ONED2N(i-1)*0.5D0
    END DO

    ARE(0) = 0.5D0
    DO i=1, 64
      ARE(i) = ONED2N(i)*DIVIDE(ONED2N(i))
    END DO

    TWON(0) = 1.0D0
    DO i=1, 64
      TWON(i) = TWON(i-1)*2.0D0
    END DO

    IF (X .LT. 1.0D0) THEN
      RECID = DIVIDER(X)
    ELSE
      DO i=0, 63
        IF (X .GE. TWON(i) .AND. X .LT. TWON(i+1)) THEN
          Y = (X - TWON(i))*ARE(i)
          RECID = DIVIDER(Y)
          RECID = ARE(i)*RECID
        END IF
      END DO
    END IF
  END SUBROUTINE GT1DIVIDE

  !new subroutine for 1/Z, where 0<Z<1
  !1/(1-x)=(1+x+x**2+....x**N)
  SUBROUTINE LT1DIVIDE(X, RECID)  !Upgrade this part
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: X
    REAL*8, INTENT(OUT) :: RECID
    REAL*8 :: P
    INTEGER :: i, N

    P = 1.0D0
    N = 1000

    IF (X .LT. 0.0D0) STOP 'X LT 0 LT1DIVIDE'
    IF (X .GT. 1.0D0) STOP 'X GT 1 LT1DIVIDE'
    RECID = 0.0D0

    DO i=1, N
      RECID = RECID + P
      P=P*X
    END DO
  END SUBROUTINE LT1DIVIDE

  !1/Z = 1/(0.5**N + (Z-0.5**N)
  !    = (1./0.5**N) 1/(1+ (Z-0.5**N)/(0.5**N))
  !    = 2**N [1/(1 + (2^N(Z-0.5^N)))]
  SUBROUTINE TNYDIVIDE(Z,RECIDUAL)
    IMPLICIT REAL*8 (A-H,O-Z)
    REAL*8, INTENT(IN) :: Z
    REAL*8, INTENT(INOUT) :: RECIDUAL
    REAL*8 :: PP, PR
    INTEGER :: i

    PP = 1.0D0
    PR = 1.0D0
    DO i=0, 64
      IF(Z .GE. PR .and. Z .LT. PR*2.0D0)THEN
        X = (Z-PR)*PP
        RECIDUAL = PP*DIVIDER(X)
      END IF
      PR = PR*0.5D0
      PP = PP*2.0D0
    END DO
  END SUBROUTINE TNYDIVIDE
END MODULE FortranFunctions
