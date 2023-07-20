      SUBROUTINE qtrans( TT, AMU, Sunits, Ctrans, Strans, Htrans,
     &           Qtranslat )
! Computes heat capacity, enthalpy, and enthalpy function for classical translations
!
      IMPLICIT NONE
      REAL(8), INTENT(IN) ::            TT             ! Temperature (K)
      REAL(8), INTENT(IN) ::            AMU            ! Mass (atomic mass units)
      CHARACTER(LEN=3), INTENT(IN) ::   Sunits         ! Standard state (atm, bar, molecule/cm3, NOT)
      REAL(8), INTENT(OUT) ::           Ctrans         ! Heat capacity (units of Rgas)
      REAL(8), INTENT(OUT) ::           Strans         ! Entropy (units of Rgas)
      REAL(8), INTENT(OUT) ::           Htrans         ! Enthalpy function (units of Rgas)
      REAL(8), INTENT(OUT) ::           Qtranslat      ! Translational partition fxn
      
      IF ( Sunits.EQ.'ATM' ) THEN                      ! TRANSLATION PARTITION FUNCTION (q/V)
         Qtranslat = ( (AMU*TT)**1.5 )*0.02560786*TT   ! std. state= 1 atm
      ELSEIF ( Sunits.EQ.'BAR' ) THEN
         Qtranslat = ( (AMU*TT)**1.5 )*0.02594716*TT   ! std. state= 1 bar
      ELSEIF ( Sunits.EQ.'MCC' ) THEN
         Qtranslat = ( (AMU*TT)**1.5 )*1.879333E+20    ! std. state= 1 molecule/cc
      ELSEIF ( Sunits.EQ.'NOT' ) THEN
         Qtranslat = 1.0d+00                           ! Translation ignored
      ELSE
         WRITE(*,*) 
         WRITE(*,*) '***FATAL: Sunits not recognized: ', Sunits, ' ****'
         WRITE(*,*)          
      ENDIF
      
      IF ( Sunits.NE.'NOT' ) THEN
         Ctrans = 2.5
         Strans = 2.5 + log( Qtranslat )
         Htrans = 2.5*TT             ! H(T)-H(0)
      ELSE
         Ctrans = 0.0
         Strans = 0.0 
         Htrans = 0.0             ! H(T)-H(0)  
      ENDIF
          
      RETURN
      END SUBROUTINE qtrans
      
