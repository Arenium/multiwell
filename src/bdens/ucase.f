      SUBROUTINE ucase ( string )
!
!  Shifts an arbitrary character string to UPPER CASE on any processor.
!
! From Stephen J. Chapman, "Fortran 90/95 for Scientist and Engineers, 
! Second Edition", McGraw Hill Higher Education, Boston, 2004), 
! pp. 434-435.
!
      CHARACTER(len=*), INTENT(INOUT) :: string
      INTEGER :: i
      INTEGER :: length

      length = LEN ( string )
   
      DO i = 1 , length
         IF ( LGE(string(i:i),'a') .AND. LLE(string(i:i),'z') ) THEN
            string(i:i) = ACHAR ( IACHAR ( string(i:i) ) - 32 )
         ENDIF
      END DO
      
      END SUBROUTINE ucase

