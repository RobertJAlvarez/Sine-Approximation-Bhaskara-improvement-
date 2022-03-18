! Programmed by Robert Alvarez
! Last modified: March 18th 2022
!
! Compute the SIN(theta), where theta is in radians, using a better approximation derived from Bhaskara's approximation.
MODULE SINE_COSINE
  PRIVATE
  PUBLIC :: PI, COSINE, SINE

  REAL*8, PARAMETER :: PI = 3.14159265358979D0

  CONTAINS

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

    x = MOD(ABS(num),PI)

    temp = (x/PI)*(x/PI - 1.0D0)
    SINE = temp*(2.21652D0*(temp - 31.0D0/36.0D0) - 1.5372D0/(1.25D0 + temp))

    !Adjust for negative angles and shift the graph down by 2.6D-5 to minimize error
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
END MODULE SINE_COSINE
