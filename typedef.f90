MODULE typedef

IMPLICIT NONE

!Symbolic names for kind types of 4-, 2-, and 1-byte integers: 
INTEGER, PARAMETER :: I4B = SELECTED_INT_KIND(9) 
INTEGER, PARAMETER :: I2B = SELECTED_INT_KIND(4) 
INTEGER, PARAMETER :: I1B = SELECTED_INT_KIND(2)

! Not available in Intel ?: 
! INTEGER, PARAMETER :: QP = SELECTED_REAL_KIND( P=16,  R=275 ) 

!Symbolic names for kind types of single- and DOUBLE-PRECISION reals: 
INTEGER, PARAMETER :: sp = KIND(1.0) 
INTEGER, PARAMETER :: dp = KIND(1.0D0) 

!Symbolic names for kind types of single- and DOUBLE-PRECISION COMPLEX: 
INTEGER, PARAMETER :: spc = KIND((1.0,1.0)) 
INTEGER, PARAMETER :: dpc = KIND((1.0D0,1.0D0)) 

!Symbolic name for kind TYPE of default LOGICAL: 
INTEGER, PARAMETER :: LGT = KIND(.TRUE.) 

!imaginary unit
COMPLEX(KIND=dpc), PARAMETER :: iunit = (0.d0, 1.d0)           
!imagianry zero
COMPLEX(KIND=dpc), PARAMETER :: czero = (0.d0, 0.d0)           
!pi
real(KIND=dp), parameter :: pi=3.14159265358979323846264338327950d0

 TYPE igrid
    !specifies a rectangular 2d grid for plotting
    INTEGER :: xs
    INTEGER :: ys
    REAL(KIND=dp), POINTER, DIMENSION(:) :: x1a
    REAL(KIND=dp), POINTER, DIMENSION(:) :: x2a
 END TYPE igrid

CONTAINS

END MODULE typedef
